//A LabBench consists of (1) a System, (2) an MC (solver), (3) an Observer, and (4) a ParamDict.
//The LabBench uses the ParamDict to read in and store parameters.
//It creates a System from the parameters, and uses an MC to advance the System in time.
//The Observer writes data from the System to file over time.
//A LabBench optionally takes in a filename and a seed.
//If the filename is not empty, then the ParamDict reads in parameters from the corresponding file.
//If the seed is specified, then LabBench uses that number to seed the RNG. Otherwise, it uses a seed of 1.

#ifndef LABBENCH_HPP
#define LABBENCH_HPP

#include <string>
#include <iostream>

#include "System.hpp"
#include "MC.hpp"
#include "Observer.hpp"
#include "ParamDict.hpp"
#include "IO.hpp"

class LabBench
{
public:
    ParamDict params;
    System sys;
    Observer obs;
    MC solver;
    gsl_rng *rg;
    int seed;
    int stop_early = 0;
    int stop_slow = 1;
    int max_size = -1;

    std::string simulation = "standard";
    std::string initial_config = "triangle";
    int equil_steps = 0;
    int production_steps = 0;
    int info_freq = 1;

    /*** Methods ***/

    //constructor
    LabBench(ParamDict& theParams, gsl_rng*& theGen);

    //destructor
    ~LabBench();

    //run trajectory
    void run_equil(int nstps=-1);
    void run(int nstps=-1, std::string subdir="/", int config_freq=-1, int therm_freq=-1);

    //ffs functions
    auto run_ffs_stage1(int N0, double l0, double la, double lb);

    //simulations
    void do_simulation(std::string expt); //TODO: keep this Public, move rest to private
    void run_standard_simulation(); //equilibrate and then observe an unperturbed system
    void run_ffs_simulation(); //forward flux sampling
    void run_us_simulation(); //umbrella sampling
    void test_dimer_drug_removal(); 
    void test_dimer_removal();
    void test_monomer_removal();

};

#endif
