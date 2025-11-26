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
    int do_vertex_only = 0; //Only do vertex relaxation moves (no adding dimers or drug)
    int debug_sheet = 0; //Flag to print out probabilities of CD addition moves

    int frame = 0;
    int sweep_count = 0;
    int monomeradded = 0;
    int dimeradded = 0;
    int monomerremoved = 0;
    int dimerremoved = 0;
    int drugadded = 0;
    int drugremoved = 0;
    int typechanged = 0;
    int fusion = 0;
    int fission = 0;
    int wedgefusion = 0;
    int wedgefission = 0;
    int boundtri = 0;
    int binding = 0;
    int unbinding = 0;
    int minhe_fission = 50;

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
    int old_attempt_bind_wedge_dimer(System &g, int heid0);
    int attempt_bind_wedge_dimer(System &g, int heid0);

    int old_attempt_unbind_wedge_dimer(System &g, int heid0);
    int attempt_unbind_wedge_dimer(System &g, int vid0);

    int attempt_bind_triangle(System &g, int heid0);
    int attempt_unbind_triangle(System &g, int heid0);

    int attempt_add_drug(System &g, int heid0);
    int attempt_remove_drug(System &g, int heid0);
    void make_seed(System &g);
    void make_seed_T3(System &g);
    void get_dimer_etypes(int etypeheid0, int etypenew1, int etypenew2);
    int force_add_monomer_with_next(System &g, int heid0, int xid);


};

#endif