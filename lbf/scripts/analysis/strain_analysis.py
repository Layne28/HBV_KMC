#!/usr/bin/env python3
"""Elastic strain in a disassembling capsid, from HBV_KMC restart snapshots.

Strain is the actual value minus the preferred value the simulation uses
(System.cpp):
  bond length    l - l0[type], or l - l0h[type] for a holo dimer       (one per dimer)
  binding angle  phi - phi0[get_angle_type(type, next type)]           (one per contact)
                 phi = angle at the shared vertex between a dimer and the
                 dimer it is bonded to (dimer_bend_energy)
  dihedral       theta - theta0[get_dihedral_type(...)]                (one per dimer with
                 complete triangles on both sides)
                 theta = angle between the two triangles' normals, taken negative
                 when the fold is concave (bend_energy). For dihedral types with
                 theta0 <= 0 the code uses 1 - cos(theta), so the strain is |theta|.
Angles are reported in degrees, lengths in simulation units.

Outputs (in --outdir):
  strain_summary.csv          per snapshot: sweep, dimers, holo fraction, mean/RMS strains
  strain_vs_time.png          dimers, holo fraction and RMS strains vs sweep
  strain_distributions.png    strain histograms at several disassembly stages
  holo_correlation.png        P(holo at r | holo at 0) and the length-strain correlation C(r)
  strain_snapshot_<sweep>.png bonds colored by length / binding-angle / dihedral strain

Usage:
  python strain_analysis.py RUN_DIR INPUT.in [--outdir DIR] [--stages 5] [--frames S1 S2 ...]
RUN_DIR is a seed directory holding restart_lammps_<sweep>.dat files, e.g.
data/disassembly_crosslinked_extreme_new/seed=1.
"""

import argparse
import csv
import glob
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["text.usetex"] = False  # labels use Unicode symbols, which LaTeX rejects
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from mpl_toolkits.mplot3d.art3d import Line3DCollection

# Colors: diverging blue <-> red with a gray midpoint for signed strain, and a
# one-hue blue ramp (light -> dark = early -> late) for disassembly stages.
DIVERGING = LinearSegmentedColormap.from_list("strain", ["#2a78d6", "#f0efec", "#e34948"])
STAGE_RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]
NO_VALUE = "#c3c2b7"  # dimers without that kind of strain (no contact / no two faces)
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e4e3df"


# ----------------------------------------------------------------------------
# Parameters: same defaults and override rules as System::System and
# System::do_paramdict_assign
# ----------------------------------------------------------------------------
def read_params(path):
    params = {}
    with open(path) as f:
        for line in f:
            line = "".join(line.split())
            if ":" not in line:
                continue
            key, value = line.split(":", 1)
            params[key] = value  # later entries win, as in ParamDict
    return params


def preferred_values(params):
    l0 = [1.05, 0.95, 0.95, 1.05]
    l0h = [1.05, 0.95, 0.95, 1.05]
    phi0 = [1.05, 1.17, 0.98, 1.05]
    theta0 = [-1.0, -1.0, -1.0, -1.0]
    get = lambda k: float(params[k])
    for i in range(4):
        if f"phi{i}" in params:
            phi0[i] = get(f"phi{i}")
    if "theta0" in params:
        theta0[0] = get("theta0")
    if "theta1" in params:
        theta0[1] = get("theta1")
    if "l0_CD" in params:
        l0[0] = l0[3] = get("l0_CD")
    if "l0_AB" in params:
        l0[1] = l0[2] = get("l0_AB")
    if "l0_CD_holo" in params:
        l0h[0] = l0h[3] = get("l0_CD_holo")
    if "l0_AB_holo" in params:
        l0h[1] = l0h[2] = get("l0_AB_holo")
    for i in range(4):
        if f"l{i}" in params:
            l0[i] = get(f"l{i}")
    for i in range(4):
        if f"l{i}h" in params:
            l0h[i] = get(f"l{i}h")
        elif ("l0_CD_holo" if i in (0, 3) else "l0_AB_holo") not in params:
            l0h[i] = l0[i]
    theta0[2] = theta0[1]
    theta0[3] = 0.0
    return np.array(l0), np.array(l0h), np.array(phi0), np.array(theta0)


def check_against_run(run_dir, l0, l0h):
    """Warn if the bond lengths differ from those the run printed."""
    path = os.path.join(run_dir, "prod", "parameters_run.out")
    if not os.path.exists(path):
        return
    with open(path) as f:
        for line in f:
            m = re.match(r"bond type (\d) l0_apo ([\d.]+) l0_holo ([\d.]+)", line)
            if m:
                i, apo, holo = int(m.group(1)), float(m.group(2)), float(m.group(3))
                if abs(apo - l0[i]) > 1e-3 or abs(holo - l0h[i]) > 1e-3:
                    print(f"WARNING: bond type {i}: input gives l0={l0[i]}, l0h={l0h[i]} but the run used "
                          f"l0={apo}, l0h={holo}. Is this the right input file?", file=sys.stderr)


def angle_type(etype, nexttype):
    if etype in (0, 3) and nexttype in (0, 3):
        return 0
    if etype == 1 and nexttype == 2:
        return 1
    if etype in (0, 3) and nexttype == 1:
        return 2
    if etype == 2 and nexttype in (0, 3):
        return 2
    return 0


def dihedral_type(e, n, p, oe, on, op):
    CD = (0, 3)
    if e in CD and n == 1 and p == 2 and oe in CD and on == 1 and op == 2:
        return 0
    if e == 1 and n == 2 and p in CD and oe == 2 and on in CD and op == 1:
        return 1
    if e == 2 and n in CD and p == 1 and oe == 1 and on == 2 and op in CD:
        return 1
    if (e, n, p, oe, on, op) in ((0, 1, 2, 3, 3, 3), (3, 3, 3, 0, 1, 2), (3, 1, 2, 0, 0, 0), (0, 0, 0, 3, 1, 2)):
        return 0
    return 3


# ----------------------------------------------------------------------------
# Snapshots
# ----------------------------------------------------------------------------
def load_snapshot(path):
    """Restart file -> vertex positions and halfedge arrays (all 0-based)."""
    section = None
    X, btype, vin, vout, holo = [], [], [], [], []
    nxt, prv = {}, {}
    with open(path) as f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if t[0] in ("Atoms", "Bonds", "Angles", "Impropers"):
                section = t[0]
                continue
            if section == "Atoms" and len(t) >= 7:
                X.append([float(t[4]), float(t[5]), float(t[6])])
            elif section == "Bonds" and len(t) >= 4:
                btype.append(int(t[1]) - 1)
                vin.append(int(t[2]) - 1)
                vout.append(int(t[3]) - 1)
                holo.append(int(t[4]) if len(t) > 4 else 0)
            elif section == "Angles" and len(t) >= 5:
                h = int(t[2]) - 1
                nxt[h] = int(t[3]) - 1
                prv[h] = int(t[4]) - 1
    n = len(btype)
    s = dict(X=np.array(X), type=np.array(btype), vin=np.array(vin), vout=np.array(vout),
             holo=np.array(holo, dtype=bool),
             nxt=np.array([nxt.get(h, -1) for h in range(n)]),
             prv=np.array([prv.get(h, -1) for h in range(n)]))
    by_ends = {(a, b): h for h, (a, b) in enumerate(zip(vin, vout))}
    s["op"] = np.array([by_ends[(b, a)] for a, b in zip(vin, vout)])
    return s


def compute_strains(s, l0, l0h, phi0, theta0):
    X, typ, vin, vout = s["X"], s["type"], s["vin"], s["vout"]
    nxt, prv, op = s["nxt"], s["prv"], s["op"]
    vec = X[vout] - X[vin]  # hevec
    length = np.linalg.norm(vec, axis=1)

    # one entry per dimer: the halfedge with the lower index
    dimers = np.where(np.arange(len(typ)) < op)[0]
    pref = np.where(s["holo"], l0h[typ], l0[typ])
    he_lstrain = length - pref
    d_lstrain = 0.5 * (he_lstrain[dimers] + he_lstrain[op[dimers]])
    d_mid = 0.5 * (X[vin[dimers]] + X[vout[dimers]])
    d_holo = s["holo"][dimers]
    dimer_of = np.empty(len(typ), dtype=int)
    dimer_of[dimers] = np.arange(len(dimers))
    dimer_of[op[dimers]] = np.arange(len(dimers))

    # binding angles: one per contact h -> next
    a_strain, a_dimers = [], []
    for h in np.where(nxt >= 0)[0]:
        n = nxt[h]
        a, b = -vec[h], vec[n]  # both leave the shared vertex vout[h]
        c = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
        phi = np.arccos(min(1.0, c))  # as dimer_bend_energy: only clamped above
        a_strain.append(phi - phi0[angle_type(typ[h], typ[n])])
        a_dimers.append((dimer_of[h], dimer_of[n]))
    a_strain = np.degrees(np.array(a_strain))

    # dihedrals: halfedges with complete triangles on both sides
    def normal(h):
        nv = np.cross(vec[h], vec[nxt[h]])
        return nv / np.linalg.norm(nv)

    he_dstrain = np.full(len(typ), np.nan)
    for h in range(len(typ)):
        o = op[h]
        if min(nxt[h], prv[h], nxt[o], prv[o]) < 0:
            continue
        nh, no = normal(h), normal(o)
        ndot = np.dot(no, nh)
        theta = 0.0 if ndot >= 1 else np.arccos(ndot)
        t = dihedral_type(typ[h], typ[nxt[h]], typ[prv[h]], typ[o], typ[nxt[o]], typ[prv[o]])
        if theta0[t] > 0:
            concave = np.dot(X[vout[nxt[o]]] - X[vout[nxt[h]]], nh) > 0
            he_dstrain[h] = (-theta if concave else theta) - theta0[t]
        else:
            he_dstrain[h] = theta
    d_dstrain = np.degrees(0.5 * (he_dstrain[dimers] + he_dstrain[op[dimers]]))

    # per-dimer mean binding-angle strain over the contacts it takes part in
    d_astrain = np.full(len(dimers), np.nan)
    if len(a_strain):
        sums = np.zeros(len(dimers))
        counts = np.zeros(len(dimers))
        for value, (i, j) in zip(a_strain, a_dimers):
            for k in (i, j):
                sums[k] += value
                counts[k] += 1
        has = counts > 0
        d_astrain[has] = sums[has] / counts[has]

    return dict(length=d_lstrain, angle=a_strain, dihedral=d_dstrain[~np.isnan(d_dstrain)],
                d_angle=d_astrain, d_dihedral=d_dstrain, mid=d_mid, holo=d_holo,
                ends=np.stack([X[vin[dimers]], X[vout[dimers]]], axis=1))


def find_snapshots(run_dir):
    files = glob.glob(os.path.join(run_dir, "restart_lammps_*.dat"))
    snaps = []
    for f in files:
        m = re.search(r"restart_lammps_(\d+)\.dat$", f)
        if m:
            snaps.append((int(m.group(1)), f))
    return sorted(snaps)


# ----------------------------------------------------------------------------
# Spatial correlations
# ----------------------------------------------------------------------------
def pair_counts(st, edges):
    """Histogram dimer pairs (i != j) by center distance: all pairs, holo-anything
    pairs, holo-holo pairs, and sum of delta_i delta_j of the length strain."""
    mid, holo = st["mid"], st["holo"]
    d = np.linalg.norm(mid[:, None, :] - mid[None, :, :], axis=2)
    off = ~np.eye(len(mid), dtype=bool)
    delta = st["length"] - st["length"].mean()
    prod = delta[:, None] * delta[None, :]
    r = d[off]
    h_any = np.repeat(holo, len(mid)).reshape(len(mid), len(mid))[off]
    h_h = (holo[:, None] & holo[None, :])[off]
    return (np.histogram(r, edges)[0],
            np.histogram(r[h_any], edges)[0],
            np.histogram(r[h_h], edges)[0],
            np.histogram(r, edges, weights=prod[off])[0],
            (delta ** 2).sum(), len(delta))


# ----------------------------------------------------------------------------
# Plot helpers
# ----------------------------------------------------------------------------
def style_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK_SECONDARY)
    ax.tick_params(colors=INK_SECONDARY, labelsize=9)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def stage_colors(k):
    idx = np.linspace(0, len(STAGE_RAMP) - 1, k).round().astype(int)
    return [STAGE_RAMP[i] for i in idx]


def active_part(results):
    """Snapshots up to the first one that reaches the final size: after that the
    run only holds the end state (e.g. the last triangle)."""
    final = results[-1]["n"]
    last = len(results) - 1
    while last > 0 and results[last - 1]["n"] == final:
        last -= 1
    return results[:last + 1]


def make_stages(results, k):
    """Group snapshots into k stages of equal width in the number of dimers left."""
    ns = np.array([r["n"] for r in results])
    edges = np.linspace(ns.max() + 1e-9, ns.min() - 1e-9, k + 1)
    stages = []
    for a, b in zip(edges[:-1], edges[1:]):
        members = [r for r in results if b < r["n"] <= a]
        if members:
            lo, hi = min(r["n"] for r in members), max(r["n"] for r in members)
            stages.append((f"{hi}–{lo} dimers" if hi != lo else f"{hi} dimers", members))
    return stages


def plot_distributions(stages, outdir):
    """Small multiples: one row per disassembly stage, one column per strain kind."""
    kinds = [("length", "Bond length strain  l − l₀", "simulation units"),
             ("angle", "Binding angle strain  φ − φ₀", "degrees"),
             ("dihedral", "Dihedral strain  θ − θ₀", "degrees")]
    fig, axes = plt.subplots(len(stages), 3, figsize=(13, 1.9 * len(stages) + 0.8),
                             sharex="col", squeeze=False, constrained_layout=True)
    for col, (key, title, unit) in enumerate(kinds):
        allv = np.concatenate([np.concatenate([r["strain"][key] for r in m]) for _, m in stages])
        bins = np.linspace(*np.percentile(allv, [0.5, 99.5]), 30) if len(allv) else 10
        for row, (label, members) in enumerate(stages):
            ax = axes[row, col]
            v = np.concatenate([r["strain"][key] for r in members])
            if len(v):
                ax.hist(v, bins=bins, density=True, color=STAGE_RAMP[2], edgecolor="white", linewidth=0.5)
                ax.axvline(np.median(v), color=INK, linewidth=1.5)
            ax.axvline(0, color=INK_SECONDARY, linewidth=1, linestyle=":")
            ax.text(0.98, 0.92, f"n = {len(v)}" + (f"\nmedian {np.median(v):.3g}" if len(v) else ""),
                    transform=ax.transAxes, ha="right", va="top", fontsize=8, color=INK_SECONDARY)
            style_axes(ax)
            if row == 0:
                ax.set_title(title, color=INK, fontsize=11, loc="left")
            if col == 0:
                ax.set_ylabel(label + "\n\ndensity", color=INK_SECONDARY, fontsize=9)
            if row == len(stages) - 1:
                ax.set_xlabel(unit, color=INK_SECONDARY, fontsize=9)
    fig.suptitle("Strain distributions by disassembly stage (black line: median, dotted: zero strain)",
                 color=INK, fontsize=12, x=0.01, ha="left")
    fig.savefig(os.path.join(outdir, "strain_distributions.png"), dpi=150)
    plt.close(fig)


def plot_vs_time(results, outdir):
    sweeps = np.array([r["sweep"] for r in results])
    rows = [("Dimers remaining", [r["n"] for r in results], ""),
            ("Holo fraction", [r["holo_fraction"] for r in results], ""),
            ("RMS bond length strain", [r["rms_length"] for r in results], "simulation units"),
            ("RMS binding angle strain", [r["rms_angle"] for r in results], "degrees"),
            ("RMS dihedral strain", [r["rms_dihedral"] for r in results], "degrees")]
    fig, axes = plt.subplots(len(rows), 1, figsize=(8, 10), sharex=True, constrained_layout=True)
    for ax, (title, y, unit) in zip(axes, rows):
        ax.plot(sweeps, y, color=STAGE_RAMP[2], linewidth=2)
        ax.set_title(title, color=INK, fontsize=10, loc="left")
        ax.set_ylabel(unit, color=INK_SECONDARY, fontsize=9)
        style_axes(ax)
    axes[-1].set_xlabel("sweep", color=INK_SECONDARY, fontsize=9)
    fig.savefig(os.path.join(outdir, "strain_vs_time.png"), dpi=150)
    plt.close(fig)


MIN_PAIRS = 30  # distance bins with fewer dimer pairs than this are not plotted


def plot_correlations(stages, edges, outdir):
    centers = 0.5 * (edges[:-1] + edges[1:])
    colors = stage_colors(len(stages))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    for (label, members), c in zip(stages, colors):
        n_all = sum(m["pairs"][0] for m in members)
        n_h = sum(m["pairs"][1] for m in members)
        n_hh = sum(m["pairs"][2] for m in members)
        cov = sum(m["pairs"][3] for m in members)
        var = sum(m["pairs"][4] for m in members) / max(1, sum(m["pairs"][5] for m in members))
        holo_frac = (sum(np.sum(m["strain"]["holo"]) for m in members) /
                     max(1, sum(len(m["strain"]["holo"]) for m in members)))
        with np.errstate(invalid="ignore", divide="ignore"):
            p = np.where(n_h >= MIN_PAIRS, n_hh / n_h, np.nan)
            corr = np.where(n_all >= MIN_PAIRS, cov / n_all / var, np.nan) if var > 0 else np.full_like(centers, np.nan)
        ax1.plot(centers, p, color=c, linewidth=2, marker="o", markersize=4, label=label)
        ax1.axhline(holo_frac, color=c, linewidth=1, linestyle="--")
        ax2.plot(centers, corr, color=c, linewidth=2, marker="o", markersize=4, label=label)
    ax1.set_title("P(holo at r | holo at 0)   (dashed: overall holo fraction)", color=INK, fontsize=11, loc="left")
    fig.suptitle(f"Distance bins with fewer than {MIN_PAIRS} dimer pairs are omitted", color=INK_SECONDARY,
                 fontsize=9, x=0.01, ha="left")
    ax1.set_ylim(-0.02, 1.02)
    ax2.set_title("Bond length strain correlation  C(r) = ⟨δlᵢ δlⱼ⟩ / ⟨δl²⟩", color=INK, fontsize=11, loc="left")
    ax2.axhline(0, color=INK_SECONDARY, linewidth=1, linestyle=":")
    for ax in (ax1, ax2):
        ax.set_xlabel("distance r between dimer centers (simulation units)", color=INK_SECONDARY, fontsize=9)
        style_axes(ax)
    ax1.set_ylabel("probability", color=INK_SECONDARY, fontsize=9)
    ax2.set_ylabel("correlation", color=INK_SECONDARY, fontsize=9)
    ax1.legend(title="stage", fontsize=8, title_fontsize=8, frameon=False, labelcolor=INK)
    fig.savefig(os.path.join(outdir, "holo_correlation.png"), dpi=150)
    plt.close(fig)


def plot_snapshot(r, limits, outdir):
    st = r["strain"]
    panels = [("Bond length strain", st["length"], limits["length"], "simulation units"),
              ("Binding angle strain (mean over the dimer's contacts)", st["d_angle"], limits["angle"], "degrees"),
              ("Dihedral strain", st["d_dihedral"], limits["dihedral"], "degrees")]
    fig = plt.figure(figsize=(16, 5.6), constrained_layout=True)
    ends = st["ends"]
    center = ends.reshape(-1, 3).mean(axis=0)
    span = np.abs(ends.reshape(-1, 3) - center).max()
    widths = np.where(st["holo"], 3.5, 1.5)
    for k, (title, values, lim, unit) in enumerate(panels):
        ax = fig.add_subplot(1, 3, k + 1, projection="3d")
        norm = TwoSlopeNorm(vcenter=0.0, vmin=-lim, vmax=lim)
        has = ~np.isnan(values)
        if (~has).any():
            ax.add_collection3d(Line3DCollection(ends[~has], colors=NO_VALUE, linewidths=widths[~has],
                                                 linestyles="dotted"))
        lc = Line3DCollection(ends[has], cmap=DIVERGING, norm=norm, linewidths=widths[has])
        lc.set_array(values[has])
        ax.add_collection3d(lc)
        for axis in "xyz":
            getattr(ax, f"set_{axis}lim")(center["xyz".index(axis)] - span, center["xyz".index(axis)] + span)
        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=20, azim=35)
        ax.set_axis_off()
        ax.set_title(title, color=INK, fontsize=11, loc="left")
        cb = fig.colorbar(lc, ax=ax, shrink=0.7, pad=0.02, extend="both")
        cb.set_label(unit, color=INK_SECONDARY, fontsize=9)
        cb.ax.tick_params(labelsize=8, colors=INK_SECONDARY)
    fig.suptitle(f"Sweep {r['sweep']}: {r['n']} dimers, {int(st['holo'].sum())} holo   "
                 "(thick = holo, thin = apo, dotted gray = no value)", color=INK, fontsize=12, x=0.01, ha="left")
    fig.savefig(os.path.join(outdir, f"strain_snapshot_{r['sweep']}.png"), dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir", help="seed directory with restart_lammps_<sweep>.dat files")
    ap.add_argument("input", help="the .in file the run used (for the preferred values)")
    ap.add_argument("--outdir", default=None, help="output directory (default: RUN_DIR/strain)")
    ap.add_argument("--stages", type=int, default=4, help="number of disassembly stages for the distributions")
    ap.add_argument("--include-final", action="store_true",
                    help="keep the snapshots after the structure reaches its final size in the stage plots")
    ap.add_argument("--frames", type=int, nargs="*", default=None,
                    help="sweeps to draw (default: 4 spread over the run)")
    ap.add_argument("--every", type=int, default=1, help="use every N-th snapshot")
    ap.add_argument("--dr", type=float, default=0.25, help="distance bin width for the correlations")
    args = ap.parse_args()

    outdir = args.outdir or os.path.join(args.run_dir, "strain")
    os.makedirs(outdir, exist_ok=True)
    l0, l0h, phi0, theta0 = preferred_values(read_params(args.input))
    check_against_run(args.run_dir, l0, l0h)
    print(f"l0 {l0}  l0h {l0h}  phi0 {phi0}  theta0 {theta0}")

    snaps = find_snapshots(args.run_dir)[::args.every]
    if not snaps:
        sys.exit(f"no restart_lammps_*.dat files in {args.run_dir}")
    results = []
    for sweep, path in snaps:
        st = compute_strains(load_snapshot(path), l0, l0h, phi0, theta0)
        rms = lambda v: float(np.sqrt(np.mean(v ** 2))) if len(v) else np.nan
        mean = lambda v: float(np.mean(v)) if len(v) else np.nan
        results.append(dict(sweep=sweep, n=len(st["length"]), strain=st,
                            holo_fraction=float(st["holo"].mean()),
                            mean_length=mean(st["length"]), rms_length=rms(st["length"]),
                            mean_angle=mean(st["angle"]), rms_angle=rms(st["angle"]),
                            mean_dihedral=mean(st["dihedral"]), rms_dihedral=rms(st["dihedral"])))
    print(f"{len(results)} snapshots, sweeps {results[0]['sweep']}–{results[-1]['sweep']}, "
          f"{results[0]['n']} → {results[-1]['n']} dimers")

    with open(os.path.join(outdir, "strain_summary.csv"), "w", newline="") as f:
        cols = ["sweep", "n", "holo_fraction", "mean_length", "rms_length",
                "mean_angle", "rms_angle", "mean_dihedral", "rms_dihedral"]
        w = csv.writer(f)
        w.writerow(["sweep", "dimers", "holo_fraction", "mean_length_strain", "rms_length_strain",
                    "mean_angle_strain_deg", "rms_angle_strain_deg", "mean_dihedral_strain_deg",
                    "rms_dihedral_strain_deg"])
        for r in results:
            w.writerow([r[c] for c in cols])

    rmax = max(np.ptp(r["strain"]["mid"], axis=0).max() for r in results if r["n"] > 1) * 1.05
    edges = np.arange(0, rmax + args.dr, args.dr)
    for r in results:
        r["pairs"] = pair_counts(r["strain"], edges)

    staged = results if args.include_final else active_part(results)
    if len(staged) < len(results):
        print(f"stages use sweeps {staged[0]['sweep']}–{staged[-1]['sweep']}; the {len(results) - len(staged)} later "
              f"snapshots only hold the final {results[-1]['n']}-dimer structure (--include-final to keep them)")
    stages = make_stages(staged, args.stages)
    plot_vs_time(results, outdir)
    plot_distributions(stages, outdir)
    plot_correlations(stages, edges, outdir)

    # shared, symmetric color limits so snapshots are comparable
    limits = {}
    for key, per_dimer in (("length", "length"), ("angle", "d_angle"), ("dihedral", "d_dihedral")):
        v = np.concatenate([r["strain"][per_dimer] for r in results])
        v = np.abs(v[~np.isnan(v)])
        limits[key] = float(np.percentile(v, 98)) if len(v) else 1.0
    if args.frames:
        chosen = [r for r in results if r["sweep"] in set(args.frames)]
        missing = set(args.frames) - {r["sweep"] for r in chosen}
        if missing:
            print(f"WARNING: no snapshot for sweeps {sorted(missing)}", file=sys.stderr)
    else:
        # spread over the disassembly by dimers remaining, not by time
        targets = np.linspace(staged[0]["n"], staged[-1]["n"], 4)
        chosen = []
        for t in targets:
            r = next(r for r in staged if r["n"] <= t + 1e-9)
            if r not in chosen:
                chosen.append(r)
    for r in chosen:
        plot_snapshot(r, limits, outdir)
    print(f"wrote figures and strain_summary.csv to {outdir}")


if __name__ == "__main__":
    main()
