//An Observer writes out information about a System to file

#ifndef OBSERVER_HPP
#define OBSERVER_HPP

#include <string>
#include <iomanip>
#include <filesystem>
#include "MC.hpp"
#include "System.hpp"
//#include <H5Cpp.h>
//#include "hdf5"
//#include <highfive/H5File.hpp>
//#include <highfive/H5DataSet.hpp>
//#include <highfive/H5DataSpace.hpp>

namespace fs = std::filesystem;

class System;
class MC;

class Observer
{
private:
    friend class System;
    friend class MC;

public:
    std::string output_dir = "./"; //where to dump output

    int particles_freq=10;
    int thermo_freq=10;
    int freq_log=10;
    int print_freq=1000; //how often to print to console
    //int do_h5md=1;

    /*** Methods ***/

    //constructor
    Observer(ParamDict &theParams);

    //destructor
    ~Observer();

    //Output
    //void open_h5md(System &theSys, std::string subdir);
    //void dump_h5md(System &theSys, std::string subdir); //write all particle data to hdf5 file
    void dump_parameters(System &theSys, std::string subdir);
    void dump_lammps_traj(System &g, int time0);
    void dump_lammps_traj_restart(System &g, int time0);
    void dump_lammps_data_file(System &g, int time0);
    void dump_lammps_traj_dimers(System &g, int time0);
    void dump_lammps_data_dimers(System &g, int time0);

    void dump_restart_lammps_data_file(System &g, int time0);
    void dump_data_frame(System &g, FILE *f, int time);
    void dump_analysis(System &g, FILE *ofile, int sweep = -1, int seed = -1, int seconds = -1);

};


#endif
