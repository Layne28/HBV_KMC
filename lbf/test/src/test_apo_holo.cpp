// Standalone test for the apo/holo conformational switch (MC::attempt_switch_apo_holo).
//
// Builds the smallest possible system -- two vertices joined by a single AB
// edge (one halfedge pair, no triangle/mesh machinery needed since the
// apo/holo move only touches g.he[].holo, g.dG_apoholo and g.T) -- and
// repeatedly attempts the apo/holo switch on that one edge. The holo state
// is parameterized to have an equilibrium bond length 25% longer than the
// apo state, and dG_apoholo is set to 2 (in units of kT, since T=1 in this
// code). Detailed balance predicts a time-averaged population ratio
// N_apo/N_holo = exp(dG_apoholo/T) = exp(2); this test measures that ratio
// from the simulation and checks it against the analytical prediction.

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>

#include "../../src/System.hpp"
#include "../../src/MC.hpp"
#include "../../src/ParamDict.hpp"
#include "../../src/CustomRandom.hpp"

int main(int argc, char *argv[])
{
    long unsigned int seed = 1;
    if (argc > 1) seed = std::strtoul(argv[1], nullptr, 10);

    //Physical setup for this test
    const double l0_apo = 1.0;
    const double l0_holo = 1.25 * l0_apo; // holo state is 25% longer
    const double dG_apoholo = 2.0;        // G_holo - G_apo, in units of kT (T=1)
    const long int nsteps = 5000000;
    const long int burnin = 10000;
    const double tol = 0.05; // allowed relative error on the measured ratio

    //Minimal parameter set. Most of these aren't exercised by the apo/holo
    //move itself, but are set to sane values so System() constructs cleanly.
    ParamDict params;
    params.add_entry("epsilon0", "4200");
    params.add_entry("kappa0", "40");
    params.add_entry("kappaPhi0", "800");
    params.add_entry("theta0", "0.240");
    params.add_entry("theta1", "0.480");
    params.add_entry("gb0", "-9.85");
    params.add_entry("muCD", "-10.8");
    params.add_entry("dmu", "-4.5");
    params.add_entry("mudrug", "-1000.0");
    params.add_entry("gdrug0", "-6.0");
    params.add_entry("l0_AB", std::to_string(l0_apo));
    params.add_entry("l0_AB_holo", std::to_string(l0_holo));
    params.add_entry("dG_apoholo", std::to_string(dG_apoholo));

    gsl_rng *rg = CustomRandom::init_rng(seed);

    System sys(params, rg);
    MC solver(sys, params, rg);

    //Build a single AB edge: 2 vertices, 1 edge (= a pair of opposite halfedges)
    double v0[3] = {0.0, 0.0, 0.0};
    double v1[3] = {l0_apo, 0.0, 0.0};
    sys.add_vertex(v0);
    sys.add_vertex(v1);
    sys.add_edge_type(0, 1, 2); // type 2 = AB; also creates the opposite (type 1 = BA) halfedge

    int heid0 = 0; //id of the AB halfedge just created
    int heindex0 = sys.heidtoindex[heid0];

    //Sanity check: stretch_energy should use l0 (apo) or l0h (holo) depending
    //on the edge's current state, with the edge length l fixed at l0_apo.
    sys.he[heindex0].holo = false;
    double e_apo = sys.stretch_energy(heindex0);
    sys.he[heindex0].holo = true;
    double e_holo = sys.stretch_energy(heindex0);
    double e_apo_expected = 0.0; // l == l0_apo exactly, by construction
    double e_holo_expected = 0.5 * sys.epsilon[2] * (l0_apo - l0_holo) * (l0_apo - l0_holo);
    std::cout << "Stretch energy check: apo = " << e_apo << " (expected " << e_apo_expected
               << "), holo = " << e_holo << " (expected " << e_holo_expected << ")" << std::endl;
    bool stretch_ok = (std::fabs(e_apo - e_apo_expected) < 1e-8) && (std::fabs(e_holo - e_holo_expected) < 1e-8);

    //Reset to the apo state before running the switch dynamics
    sys.he[heindex0].holo = false;

    long int n_apo = 0, n_holo = 0;
    for (long int i = 0; i < nsteps; i++)
    {
        solver.attempt_switch_apo_holo(sys, heid0);
        if (i >= burnin)
        {
            if (sys.he[heindex0].holo)
                n_holo++;
            else
                n_apo++;
        }
    }

    double ratio_apo_holo = (double)n_apo / (double)n_holo;
    double ratio_holo_apo = (double)n_holo / (double)n_apo;
    double expected_apo_holo = std::exp(dG_apoholo);
    double expected_holo_apo = std::exp(-dG_apoholo);

    double rel_err = std::fabs(ratio_apo_holo - expected_apo_holo) / expected_apo_holo;

    std::cout << "Samples: " << (n_apo + n_holo) << " (apo = " << n_apo << ", holo = " << n_holo << ")" << std::endl;
    std::cout << "Measured apo/holo ratio:  " << ratio_apo_holo << " (expected exp(+dG_apoholo) = " << expected_apo_holo << ")" << std::endl;
    std::cout << "Measured holo/apo ratio:  " << ratio_holo_apo << " (expected exp(-dG_apoholo) = " << expected_holo_apo << ")" << std::endl;
    std::cout << "Relative error on apo/holo ratio: " << rel_err * 100.0 << "%" << std::endl;

    bool ratio_ok = rel_err < tol;

    if (stretch_ok && ratio_ok)
    {
        std::cout << "PASS: apo/holo population ratio matches exp(dG_apoholo) within " << tol * 100.0 << "%." << std::endl;
        return 0;
    }
    else
    {
        if (!stretch_ok)
            std::cout << "FAIL: stretch_energy did not use l0/l0h as expected." << std::endl;
        if (!ratio_ok)
            std::cout << "FAIL: measured apo/holo ratio deviates from exp(dG_apoholo) by more than " << tol * 100.0 << "%." << std::endl;
        return 1;
    }
}
