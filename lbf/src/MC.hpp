//MC contains functions for doing various Monte Carlo moves, i.e. running
//the "dynamics"

#ifndef MC_HPP
#define MC_HPP

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include "MC.hpp"
#include "CustomRandom.hpp"
#include "System.hpp"

class System;

class MC
{
private:
    friend class System;
public:

    /*** Variables ***/
    double ks0; //rate of subunit binding relative to elastic relaxation
    double kd0; //rate of drug binding relative to elastic relaxation
    double kl0 = 1.0; //apo/holo switch attempts per edge per sweep (1.0 = one attempt per edge on average)
    int do_vertex_only = 0; //Only do vertex relaxation moves (no adding dimers or drug)
    int debug_sheet = 0; //Flag to print out probabilities of CD addition moves
    int debug_dimer_drug_removal = 0;
    int test_monomer_removal = 0;
    int allow_trimer_moves = 0; //Enable attempt_add/remove_trimer_dimer in sweep() (off by default: with typical mu values these can dissolve the seed triangle entirely, which existing configs don't expect)
    int allow_bridge_moves = 0; //Enable attempt_add/remove_monomer_bridge in sweep() (off by default, same reasoning as allow_trimer_moves)
    int allow_bridge_add_moves = 0; //Additionally enable attempt_add_monomer_bridge specifically (needs allow_bridge_moves too). Off by default: validated extensively for attempt_remove_monomer_bridge, but attempt_add_monomer_bridge's boundary-chain bookkeeping still has known unresolved edge cases on large, highly-doubleboundary structures (many split/merge cycles) -- see MC.cpp. Safe to enable removal alone without this.
    int allow_interior_moves = 0; //Enable attempt_remove/add_interior_dimer in sweep(): remove a dimer from between two complete faces (e.g. from a closed shell) and its reverse. Off by default so existing configs are unaffected.

    int frame = 0;
    int sweep_count = 0;
    int monomeradded = 0;
    int dimeradded = 0;
    int monomerremoved = 0;
    int dimerremoved = 0;
    int drugadded = 0;
    int drugremoved = 0;
    int typechanged = 0;
    int apoholochanged = 0;
    int interiorremoved = 0;
    int interioradded = 0;
    int fusion = 0;
    int fission = 0;
    int wedgefusion = 0;
    int wedgefission = 0;
    int boundtri = 0;
    int unboundtri = 0;
    int trimeradded = 0;
    int trimerremoved = 0;
    int bridgeadded = 0;
    int bridgeremoved = 0;
    int binding = 0;
    int unbinding = 0;
    int minhe_fission = 50;

    int del_vert_counter = 0;
    int no_del_vertex_counter = 0;

    /*** Random Number Generator ***/    
    gsl_rng *rg;

    /*** Methods ***/

    //constructor
    //TODO: change from "g" to something more descriptive throughout
    MC(System &g, ParamDict &theParams, gsl_rng *&the_rg);

    //destructor
    ~MC();

    //MC moves
    void sweep(System &g);
    void run_relax(System &g, int nsteps);
    void run_production(System &g, int nsteps);

    void move_vertices(System &g);
    int move_one_vertex(System &g, int vid0);
    int attempt_add_monomer(System &g, int heid0);
    int attempt_add_dimer(System &g, int heid0);
    int attempt_add_monomer_dimer(System &g, int heid0);
    int attempt_remove_monomer(System &g, int heid0);
    int attempt_remove_dimer(System &g, int heid0);
    int attempt_remove_monomer_dimer(System &g, int heid0);
    int old_attempt_vertex_fusion(System &g, int heid0);
    int attempt_vertex_fusion(System &g);

    int attempt_add_monomer_dimer_drug(System &g, int heid0);
    int attempt_remove_monomer_dimer_drug(System &g, int heid0);

    int attempt_wedge_fusion(System &g);
    int attempt_wedge_fission(System &g);

    int attempt_fusion(System &g);
    int attempt_fission(System &g);

    int old_attempt_vertex_fission(System &g, int heid0);
    int attempt_vertex_fission(System &g);
    int attempt_change_edge_type(System &g, int heid0);
    int attempt_change_edge_type_tri(System &g, int heid0);
    int attempt_switch_apo_holo(System &g, int heid0);
    double propose_new_edge_states(System &g, const std::vector<int> &heids, double &proposal_factor);
    double conf_free_energy(System &g, const std::vector<int> &heids);
    double solution_conf_free_energy(System &g);
    double local_elastic_energy(System &g, const std::vector<int> &heids);
    bool interior_dimer_closed_ok(System &g, int heid0);
    int attempt_remove_interior_dimer(System &g);
    int attempt_add_interior_dimer(System &g);
    int old_attempt_bind_wedge_dimer(System &g, int heid0);
    int attempt_bind_wedge_dimer(System &g, int heid0);

    int old_attempt_unbind_wedge_dimer(System &g, int heid0);
    int attempt_unbind_wedge_dimer(System &g, int vid0);

    int attempt_bind_triangle(System &g, int heid0);
    int attempt_unbind_triangle(System &g, int heid0);

    int attempt_add_trimer_dimer(System &g);
    int attempt_remove_trimer_dimer(System &g);
    bool is_pendant_boundary_loop(System &g, int heid0);

    int attempt_add_monomer_bridge(System &g);
    int attempt_remove_monomer_bridge(System &g);
    int fresh_boundary_index(System &g);

    int attempt_add_drug(System &g, int heid0);
    int attempt_remove_drug(System &g, int heid0);
    void make_seed(System &g);
    void make_seed_T3(System &g);
    void get_dimer_etypes(int etypeheid0, int etypenew1, int etypenew2);
    int force_add_monomer_with_next(System &g, int heid0, int xid);

    int check_vout_in_vid(System &g);
    int check_vid_in_vout(System &g, int vid);

};

#endif
