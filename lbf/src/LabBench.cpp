#include "LabBench.hpp"

LabBench::LabBench(ParamDict& theParams, gsl_rng*& theGen) : sys(theParams, theGen), obs(theParams), solver(sys, theParams, theGen)
{

    params = theParams;
    rg = theGen;
    sys.set_obs(obs);

    if(theParams.is_key("equil_steps")) equil_steps = std::stoi(theParams.get_value("equil_steps"));
    if(theParams.is_key("production_steps")) production_steps = std::stoi(theParams.get_value("production_steps"));
    if(theParams.is_key("info_freq")) info_freq = std::stoi(theParams.get_value("info_freq"));
    if(theParams.is_key("simulation")) simulation = theParams.get_value("simulation");
    if(theParams.is_key("initial_config")) initial_config = theParams.get_value("initial_config");
    if(theParams.is_key("seed")) seed = std::stoi(theParams.get_value("seed"));
    if(theParams.is_key("stop_early")) stop_early = std::stoi(theParams.get_value("stop_early"));
}

LabBench::~LabBench() {}

void LabBench::run_equil(int nstps)
{
    if (nstps==-1) nstps = this->equil_steps;
    std::cout << "Running equilibration for " << nstps << " steps." << std::endl;
    solver.run_relax(sys, nstps);
    std::cout << "Equilibration complete." << std::endl;
}

void LabBench::run(int nstps, std::string subdir, int config_freq, int therm_freq)
{
    if (nstps==-1) nstps = this->production_steps;
    if (config_freq==-1) config_freq = this->obs.particles_freq;
    if (therm_freq==-1) therm_freq = this->obs.thermo_freq;

    time_t timer1, timer2;
    int seconds;
    time(&timer1);
    int frame = 0; //track how many configurations have been written
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
    int deletednorate = 0;

    int minHE_update_neigh = 150;
    int lastNhe = 0;
    int lastNheGrowth=0;
    int npace=0;
    double avgpace=0;
    int avgAddInterval=2000;//10000;

    //Write out parameters to file
    obs.dump_parameters(sys, subdir);

    // set up an output file
    FILE *ofile, *finalfile, *fi, *paramfile, *anglefile;

    ofile = fopen((obs.output_dir + "/" + subdir + "/energy.dat").c_str(), "a");
    anglefile = fopen((obs.output_dir + "/" + subdir + "/angle_bonds.dat").c_str(), "a");

    for (int i=0; i<nstps; i++) {
        //Record data
        if (i%info_freq==0) std::cout << "step " << i << std::endl;

        //Do some checks and output data
        double ee = 0;
        //if (i % (obs.freq_log) == 0 ) 
	if (i % (therm_freq) == 0 )
        {
            
            if (sys.Nhe==6) recenter(sys); //set com to zero
            sys.update_boundary();
            sys.check_odd_neigh();
            ee = sys.compute_energy();

            time(&timer2);
            seconds = difftime(timer2, timer1);

            dump_analysis(sys, ofile, i, seed, seconds);
            //dump_angle_bonds(sys, anglefile, i, seed, seconds);
            if(solver.do_vertex_only == 0) {
                dump_angle_bonds(sys, anglefile, i, seed, seconds);
            }
            dump_lammps_data_file(sys, 22222222);
            //dump_lammps_traj_restart(g, sweep);
            dump_lammps_data_dimers(sys, 11111111);

        }
        if (i % (obs.freq_restart) == 0 ) 
        {
            dump_restart_lammps_data_file(sys, i);
        }

        if (i % obs.particles_freq == 0)
        {
            obs.dump_lammps_traj_dimers(sys, i);
            //obs.dump_lammps_traj_angles(sys, i);
        }

        if (i % obs.print_freq == 0 )
        {
            time(&timer2);
            seconds = difftime(timer2, timer1);
            cout << "###################################################################" << endl;
            cout << " ################ RUN TIME " << seconds << " SECONDS ###############" << endl;
            cout << " ############# SWEEP " << i << "##############" << endl;
            cout << "######### FRAME " << frame << " ##############" << endl;
            cout << "######### ENERGY " << ee << " ##############" << endl;
            cout << "######### ENERGY PER DIMER " << 2 * ee / sys.Nhe << " ##############" << endl;
            cout << "#########  NHE " << sys.Nhe << " #######################" << endl;
            cout << "#########  NHESURF " << sys.boundary.size() << " ##############" << endl;
            cout << "#########  NVSURF " << sys.boundaryv.size() << " ##############" << endl;
            cout << "#########  NV_BONDSURF " << sys.boundaryvbond.size() << " ##############" << endl;
            cout << "#########  NV5 " << sys.Nv5 << " ##############" << endl;
            cout << "#########  MONOMER ADDED " << monomeradded << " ##############" << endl;
            cout << "#########  MONOMER REMOVED " << monomerremoved << " ##############" << endl;
            cout << "#########  DIMER ADDED " << dimeradded << " ##############" << endl;
            cout << "#########  DIMER REMOVED " << dimerremoved << " ##############" << endl;
            cout << "#########  no rate REMOVED " << deletednorate << " ##############" << endl;
            cout << "#########  Surface bound " << solver.binding << " ##############" << endl;
            cout << "#########  Surface Unbound " << solver.unbinding << " ##############" << endl;
            cout << "#########  DrugAdded " << drugadded << " ##############" << endl;
            cout << "#########  DrugRemoved " << drugremoved << " ##############" << endl;
            cout << "#########  ND " << sys.Nd << " ##############" << endl;
            cout << "#########  TYPE CHANGED " << typechanged << " ##############" << endl;
            cout << "#########  WEDGE FUSION " << wedgefusion << " ##############" << endl;
            cout << "#########  WEDGE FISSION " << wedgefission << " ##############" << endl;
            cout << "#########  FUSION " << fusion << " ##############" << endl;
            cout << "#########  FISSION " << fission << " ##############" << endl;
            cout << "#########  FUSION HALFEDGES " << sys.fusionhe.size() << " ##############" << endl;
            cout << "#########  WEDGE FUSION HALFEDGES " << sys.fusionwedgehe.size() << " ##############" << endl;
            cout << "#########  ALL NEIGH " << sys.all_neigh << " ##############" << endl;
            cout << "#########  Nboundary " << sys.Nboundary << " ##############" << endl;
            cout << "#########  Bound Triangle " << boundtri << " ##############" << endl;
            cout << "#########  Acceptance vmove " << (1.0 * sys.accepted_vmove) / (1.0 * (sys.accepted_vmove + sys.rejected_vmove)) << "#################" << endl;
            cout << "#########  Nvlast " << sys.Nvlast << " Nhelast " << sys.Nhelast << " ################" << endl;
            cout << "#########  T4 " << sys.NCD_T4_in << " T3 " << sys.NCD_T3_in << " ################" << endl;
            cout << "#########  NCD_Hex" << sys.NCD_Hex << "#################"<<endl;
            cout << "#########  avgAddInterval "<<avgAddInterval<<" ###############"<<endl;
        }

        int cc = check_bind_triangle(sys);
        if (cc > 0)
        {
            cout << "bound triangle" << endl;
            sys.update_boundary();
            boundtri += cc;
        }

        sys.check_odd_neigh();

        if (sys.Nhe>70 && sys.Nhe < minHE_update_neigh && i % 10000 == 0)
        {
            double thispace=float(sys.Nhe-lastNhe)/10000.0; //pace of adding Nhe per sweep_count
            cout << "thispace "<< thispace <<endl;

            avgpace=(npace*avgpace+thispace)/(npace+1.0); //average pace of adding 

            avgAddInterval=int(pow(10,(-1 * int(floor(log10(avgpace))) ) ) );
            cout << "avgpace" << avgpace << " avgAddInterval " <<avgAddInterval << endl;
            
            npace++;
            lastNhe=sys.Nhe;
        }

        if (sys.Nhe > minHE_update_neigh && i % 10000 == 0)
        {
            sys.update_neigh();
            if (sys.find_overlap_all() < 0)
            {
                cout << "error overlap" << endl;
                dump_lammps_data_dimers(sys, 5555555);
                exit(-1);
            }
        }
        //LBF check if capsid growth is stalled
        if (sys.Nhe > minHE_update_neigh && i % (100*avgAddInterval) == 0 && stop_early==1)
        {
            if (abs( sys.Nhe - lastNheGrowth)<=4){
                std::cout << "Assembly is >35 edges and growth has slowed to <2 edges per" << 20*avgAddInterval << ". Stopping." << std::endl;
                sys.update_boundary();
                dump_restart_lammps_data_file(sys, i);
                time(&timer2);
                seconds = difftime(timer2, timer1);
                dump_analysis(sys, ofile, i, seed, seconds);
                exit(-1);
            }
            lastNheGrowth = sys.Nhe;
            lastNhe = sys.Nhe;
        }

        //see if capsid is growing or it is stalled in mixed morphology
        /*
        if (sys.Nhe > minHE_update_neigh && i % (10*avgAddInterval) == 0)
        {    
            std::cout << "Stalled? Last Nhe: " << lastNhe << " current Nhe: " << sys.Nhe << std::endl;
            sys.update_geometry_parameters();
            if ( sys.Nhe - lastNhe<=2 ){
                if ((sys.Nhe >= 220 && sys.NCD_T4_in >=26 && sys.NCD_T3_in >= 3  && sys.Nsurf > 10 ) || 
                    (sys.Nhe >= 160 && sys.NCD_T4_in >= 3 && sys.NCD_T3_in >=16  && sys.Nsurf > 10) || 
                    (sys.Nhe >= 200 && sys.NCD_T4_in >= 5 && sys.NCD_T3_in >=5  && sys.Nsurf > 10) ||
                    (sys.Nhe >= 200 && sys.NCD_T4_in >= 5 && sys.NCD_Hex>=1) || (sys.Nhe >= 200 && sys.NCD_T4_in >= 50 && sys.NCD_T3_in >=1 ))
                    {
                        cout << "STOP for now - mixed morph" << endl;
                        sys.update_boundary();
                        //dump_lammps_traj_dimers(sys, int(i));
                        //dump_lammps_data_dimers(sys, 44444444);
                        //dump_lammps_data_dimers(sys, 11111111);
                        time(&timer2);
                        seconds = difftime(timer2, timer1);
                        dump_analysis(sys, ofile, i, seed, seconds);
                       exit(-1);
                    }
            }
            
            
            if (i % (100*avgAddInterval) == 0)
            {

                //if (sys.NCD_T4_in>0 && sys.NCD_T3_in>0 && abs( sys.Nhe - lastNheGrowth)<=4 ){
                if ((sys.NCD_T4_in>0 && abs( sys.Nhe - lastNheGrowth)<=4) || (sys.NCD_T3_in>0 && abs( sys.Nhe - lastNheGrowth)<=4)){
                    fprintf(stderr, "STOP for now - not growing\n");
                    sys.update_boundary();
                    //dump_lammps_traj_dimers(g, int(sweep_count));
                    //dump_lammps_data_dimers(sys, 333333333);
                    //dump_lammps_data_dimers(sys, 11111111);
                    dump_restart_lammps_data_file(sys, i);
                    time(&timer2);
                    seconds = difftime(timer2, timer1);
                    dump_analysis(sys, ofile, i, seed, seconds);
                  exit(-1);
                }
                lastNheGrowth = sys.Nhe;
            }
            lastNhe = sys.Nhe;
        }
        */
        /*
        if (i == 200000000)
        {

            fprintf(stderr, "STOP for now - too long\n");
            sys.update_boundary();
            obs.dump_lammps_traj_dimers(sys, int(i));
            //obs.dump_lammps_traj_angles(sys, int(i));
            //dump_lammps_traj_restart(sys, int(sweep_count));
            dump_lammps_data_dimers(sys, 77777777);
            dump_restart_lammps_data_file(sys, i);
            time(&timer2);
            seconds = difftime(timer2, timer1);
            dump_analysis(sys, ofile, i, seed, seconds);
            exit(-1);
        }
        */

        //Stop early if you reach T3 or T4
        if (sys.Nhe==240 && sys.NCD_T4_in==60 && stop_early==1){
            std::cout << "Assembled T4 capsid! Stopping now." << std::endl;
            sys.update_boundary();
            obs.dump_lammps_traj_dimers(sys, int(i));
            dump_lammps_data_dimers(sys, 9999999);
            dump_restart_lammps_data_file(sys, i);
            time(&timer2);
            seconds = difftime(timer2, timer1);
            dump_analysis(sys, ofile, i, seed, seconds);
            exit(-1);
        }
        if (sys.Nhe==180 && sys.NCD_T3_in==30 && stop_early==1){
            std::cout << "Assembled T3 capsid! Stopping now." << std::endl;
            sys.update_boundary();
            obs.dump_lammps_traj_dimers(sys, int(i));
            dump_lammps_data_dimers(sys, 9999999);
            dump_restart_lammps_data_file(sys, i);
            time(&timer2);
            seconds = difftime(timer2, timer1);
            dump_analysis(sys, ofile, i, seed, seconds);
            exit(-1);
        }

        // if (sys.Nhe >= 310 || sys.Nv >= 65)
        // {

        // // fprintf(stderr, "STOP for now - too large\n");
        //     sys.update_boundary();
        //     //dump_lammps_traj_dimers(sys, int(i));
        //     //dump_lammps_data_dimers(sys, 88888888);
        //     //dump_lammps_data_dimers(sys, 11111111);
        //     dump_restart_lammps_data_file(sys, i);
        //     time(&timer2);
        //     seconds = difftime(timer2, timer1);
        //     dump_analysis(sys, ofile, i, seed, seconds);
        // // exit(-1);
        // }
        //sweep_count++;

        //Dump configuration
        //dump_lammps_traj_dimers(sys, i);
        //dump_restart_lammps_data_file(sys, i);

        /***************************** */
        //Advance dynamics
        solver.sweep(sys);
        /***************************** */
    }

    obs.dump_lammps_traj_dimers(sys, frame++);
    //obs.dump_lammps_traj_angles(sys, frame++);
    dump_lammps_data_file(sys, 22222222);
    dump_lammps_data_dimers(sys, 11111111);
    dump_restart_lammps_data_file(sys, nstps);

    time(&timer2);
    seconds = difftime(timer2, timer1);
    cout << " ################  FULL RUN TIME " << seconds << " SECONDS ###############" << endl;
    cout << " ############# SWEEP " << nstps << "##############" << endl;

    double ee = sys.compute_energy();
    cout << "######### ENERGY " << ee << " ##############" << endl;
    cout << "######### ENERGY PER DIMER " << 2 * ee / sys.Nhe << " ##############" << endl;
    cout << "#########  NHE " << sys.Nhe << " #######################" << endl;
    cout << "#########  NHESURF " << sys.boundary.size() << " ##############" << endl;
    cout << "#########  NVSURF " << sys.boundaryv.size() << " ##############" << endl;
    cout << "#########  NV_BONDSURF " << sys.boundaryvbond.size() << " ##############" << endl;
    cout << "#########  NV5 " << sys.Nv5 << " ##############" << endl;
    cout << "#########  MONOMER ADDED " << monomeradded << " ##############" << endl;
    cout << "#########  MONOMER REMOVED " << monomerremoved << " ##############" << endl;
    cout << "#########  DIMER ADDED " << dimeradded << " ##############" << endl;
    cout << "#########  DIMER REMOVED " << dimerremoved << " ##############" << endl;
    cout << "#########  no rate REMOVED " << deletednorate << " ##############" << endl;
    cout << "#########  Surface bound " << solver.binding << " ##############" << endl;
    cout << "#########  Surface Unbound " << solver.unbinding << " ##############" << endl;
    cout << "#########  DrugAdded " << drugadded << " ##############" << endl;
    cout << "#########  DrugRemoved " << drugremoved << " ##############" << endl;
    cout << "#########  ND " << sys.Nd << " ##############" << endl;
    cout << "#########  TYPE CHANGED " << typechanged << " ##############" << endl;
    cout << "#########  WEDGE FUSION " << wedgefusion << " ##############" << endl;
    cout << "#########  WEDGE FISSION " << wedgefission << " ##############" << endl;
    cout << "#########  FUSION " << fusion << " ##############" << endl;
    cout << "#########  FISSION " << fission << " ##############" << endl;
    cout << "#########  FUSION HALFEDGES " << sys.fusionhe.size() << " ##############" << endl;
    cout << "#########  WEDGE FUSION HALFEDGES " << sys.fusionwedgehe.size() << " ##############" << endl;
    cout << "#########  ALL NEIGH " << sys.all_neigh << " ##############" << endl;
    cout << "#########  Nboundary " << sys.Nboundary << " ##############" << endl;
    cout << "#########  Bound Triangle " << boundtri << " ##############" << endl;
    cout << "#########  NCD_Hex" << sys.NCD_Hex << "#################"<<endl;

    time(&timer2);
    seconds = difftime(timer2, timer1);

    dump_analysis(sys, ofile, nstps, seed, seconds);
    finalfile = fopen((obs.output_dir + "/last.dat").c_str(), "w");
    dump_analysis(sys, finalfile, nstps, seed, seconds);

    fclose(ofile);
    fclose(finalfile);
}

void LabBench::do_simulation(std::string expt)
{
    if (expt=="standard") {
        std::cout << "Running standard MC simulation." << std::endl;
        this->run_standard_simulation();
    } 
    else if (expt=="ffs") {
        std::cout << "Doing forward flux sampling." << std::endl;
        this->run_ffs_simulation();
    } 
    else if (expt=="us") {
        std::cout << "Doing umbrella sampling." << std::endl;
        this->run_us_simulation();
    } 
    else {
        std::cout << "This simulation has not been designed yet.\n" << std::endl;
        exit(0);
    }
}

void LabBench::run_standard_simulation()
{

    std::cout << "Creating initial configuration..." << std::endl;
    //TODO: add a flag in .in file to specify how to create the initial configuration
    if(initial_config=="triangle"){
        make_initial_triangle(sys);
    }
    else if(initial_config=="diamond_AB"){
        make_initial_diamond_AB(sys);
    }
    else if(initial_config=="diamond_CD"){
        //std::cout << "TEST" << std::endl;
        make_initial_diamond_CD(sys);
    }
    else if(initial_config=="diamond_DC"){
        make_initial_diamond_DC(sys);
    }
    else if(initial_config=="pentamer"){
        make_initial_pentamer(sys);
    }
    else if(ends_with(initial_config,".dat")){
        std::cout << "Attempting to read data from " << initial_config << std::endl;
        read_restart_lammps_data_file(sys, initial_config.data());
        //  FILE *f;
        //  f = fopen(initial_config.data(), "r");
        //  read_restart_lammps_data_traj(sys, f, -1);
        //make_initial_from_file(sys, initial_config);
    }
    else{
        std::cout << "Error: initial configuration type not recognized." << std::endl;
        exit(-1);
    }

    std::cout << "Equilibrating..." << std::endl;
    this->run_equil(this->equil_steps);

    std::cout << "Doing production run..." << std::endl;
    this->run(this->production_steps, "/prod", this->obs.particles_freq, this->obs.thermo_freq);
}

/******************************/
/*** Forward Flux Sampling ****/
/******************************/

auto LabBench::run_ffs_stage1(int N0, double l0, double la, double lb)
{
//     //This function has to be placed before "run_ffs_simulation"
//     //because of the use of auto return type
//     //TODO: handle case where system reaches state B (rare but possible)

//     //Define struct for output of stage 1
//     struct result {
//         double time; //total time to collect N0 configurations
//         std::vector<std::vector<Particle>> configs; //vector of configurations crossing la
//     };

//     //Define variables
//     std::vector<std::vector<Particle>> configs; //vector of configurations

//     int config_counter = 0; //no. of configs reaching la (out of N0)
//     double time = 0;        //total time to get N0 configs

//     int was_in_a = 1;       //keep track of whether system was in state A
//                             //before crossing l0

//     //Dynamics loop
//     int nsteps = 0;
//     int outfreq = 100000;
//     while (config_counter<N0){

//         if (nsteps % outfreq == 0){
//             std::cout << "time: " << time << " timesteps: " << nsteps << std::endl;
//         }

//         //Run one step
//         solver.update(sys, sys.dt);
//         time += sys.dt;
//         nsteps++;

//         if (sys.get_order_parameter() <= la){
//             was_in_a = 1;
//         }
//         //Append configuration and increment counter if 
//         //system has crossed lambda_a
//         if (sys.get_order_parameter() > l0 && was_in_a==1){
//             std::cout << "particle crossed lambda_0 at time " << time << std::endl;
//             configs.push_back(sys.particles);
//             config_counter++;
//             was_in_a = 0;
//         }
//     }
    

//     return result {time, configs};
}

void LabBench::run_ffs_simulation()
{
//     //Do forward flux sampling (ffs) using original "direct" algorithm
//     //INPUT:
//     //  -N0 (number of points at first interface)
//     //  -M0 (number of "firing runs" at first interface)
//     //  -nint (number of interfaces)
//     //  -la (op value demarking edge of state A)
//     //  -lb (op value demarking edge of state B)
//     //  -op (order parameter)

//     int N0 = 5;
//     int M0 = 5;
//     int nint = 5;
//     double la = 0.0;
//     double lb = 1.0;
//     std::string op = "single_particle_x";
//     int savefreq = 10;

//     if(params.is_key("N0")) N0 = std::stoi(params.get_value("N0"));
//     if(params.is_key("M0")) M0 = std::stoi(params.get_value("M0"));
//     if(params.is_key("nint")) nint = std::stoi(params.get_value("nint"));
//     if(params.is_key("la")) la = std::stod(params.get_value("la"));
//     if(params.is_key("lb")) lb = std::stod(params.get_value("lb"));
//     if(params.is_key("order_parameter")) op = params.get_value("order_parameter");

//     sys.order_parameter = op;

//     //Define interfaces
//     //For now, use simplest method: nint evenly spaced
//     //interfaces between la and lb
//     std::vector<double> lambdas(nint, 0.0);
//     for(int i=0; i<nint; i++){
//         lambdas[i] = la + (i+1)*(lb-la)/nint;
//         std::cout << lambdas[i] << std::endl;
//     }

//     /*****************************/
//     //Stage 1: sampling the A basin
//     /*****************************/

//     //TODO: reset system to make sure it starts in A basin
//     auto [Tstage1, stage1_configs] = run_ffs_stage1(N0, lambdas[0], la, lb);

//     //Compute flux from A to 0
//     double phi_A0 = N0/Tstage1;
//     std::cout << "Phi_A,0: " << phi_A0 << std::endl;

//     /*****************************/
//     //Stage 2: crossing interfaces
//     /*****************************/
    
//     //TODO: save configurations of successful trajectories 
//     //Create "tree" structure to save branched trajectories
//     struct segment
//     {
//         std::vector<std::vector<Particle>> data; //Trajectory
//         int parent;                              //Index of parent trajectory
//         int depth;                               //Which interface
//     };

//     //Create vector of vectors of trajectory segments 
//     //(for each interface and each successful trial)
//     std::vector<std::vector<segment>> trial_trajs;

//     //Vector to store transition probabilities at each interface
//     std::vector<double> transition_probs(nint-1, 0);

//     //Loop through interfaces
//     for(int i=1; i<nint; i++){

//         std::cout << "transition from lambda_" << (i-1) << " to lambda_" << i << std::endl;

//         //Hold successful trial trajectories for current interface
//         std::vector<segment> trials_curr;

//         int Ni=0; //number of trajectories crossing lambda_(i+1)

//         //Launch M_i (=M0 for now) trial trajectories
//         //randomly sampled from the N_i (=N0 for now) stored configurations
//         for(int j=0; j<M0; j++){

//             //std::cout << "num to select from: " << N0 << std::endl;

//             int index = gsl_rng_uniform_int(rg, N0);
//             //std::cout << index << std::endl;
//             std::vector<Particle> config_curr;
//             if (i==1){
//                 config_curr = stage1_configs[index];
//             }
//             else {
//                 config_curr = (trial_trajs[i-2][index].data).back(); //might need to change if we append A->0 segments
//             }
            
//             //Assign system the selected configuration
//             //TODO: if pbc, may also need to store and select
//             //periodic images, update cell list, etc
//             sys.particles = config_curr;

//             //Create current trajectory
//             std::vector<std::vector<Particle>> traj_curr;
//             //TODO: might need to change below if we end up 
//             //saving pre-lambda_0 segments
//             if(i==1){
//                 traj_curr.push_back(config_curr);
//             }

//             //Dynamics loop
//             int done = 0;
//             while (done==0){
//                 solver.update(sys, sys.dt);         //Advance time
//                 traj_curr.push_back(sys.particles); //Append config

//                 if (sys.get_order_parameter() <= la){
//                     //std::cout << "system returned to state A" << std::endl;
//                     done = 1;
//                 }
                
//                 if (sys.get_order_parameter() > lambdas[i]){
//                     std::cout << "particle crossed lambda_" << i << std::endl;
//                     std::cout << "(value " << sys.get_order_parameter() << ")" << std::endl;
//                     done = 1;
//                     Ni++;
//                     //Save trajectory segment
//                     segment mySeg;
//                     mySeg.data = traj_curr;
//                     mySeg.parent = index;
//                     mySeg.depth = (i-1);
//                     trials_curr.push_back(mySeg);
//                 }
//             }
//         }
//         trial_trajs.push_back(trials_curr);
//         N0 = trials_curr.size();
//         transition_probs[i-1] = (1.0*Ni)/M0;
//     }

//     std::cout << "Transition probabilities:" << std::endl;
//     for(int i=0; i<(nint-1); i++){
//         std::cout << i << "-->" << (i+1) << ": " << transition_probs[i] << std::endl;
//     }

//     double kAB = phi_A0;
//     for (int i=0; i<(nint-1); i++){
//         kAB *= transition_probs[i];
//     }
//     std::cout << "Rate constant: " << kAB << std::endl;

//     int n_successful = (trial_trajs.back()).size();

//     std::cout << "No. of successful trajectories: " << n_successful << std::endl;

//     //Gather transition paths
//     std::cout << "Gathering transition paths..." << std::endl;
//     std::vector<std::vector<std::vector<Particle>>> transition_paths;
//     for(int i=0; i<n_successful; i++){
//         transition_paths.push_back((trial_trajs.back())[i].data);
//         int curr_index = i;
//         for(int j=nint-2; j>0; j--){
//             int myindex = trial_trajs[j][curr_index].parent;
//             curr_index = myindex;
//             //std::cout << j << " " << myindex << std::endl;
//             transition_paths[i].insert((transition_paths[i]).begin(), (trial_trajs[j-1][myindex]).data.begin(), (trial_trajs[j-1][myindex]).data.end());
//         }
//     }

//     //Print interface info to file
//     std::ofstream myfile;
//     myfile.open(obs.output_dir + "/ffs/interfaces.txt");
//     myfile << "order parameter: " << sys.order_parameter << std::endl;
//     myfile << "no. of interfaces: " << nint << std::endl;
//     myfile << la << std::endl;
//     for(int i=0; i<lambdas.size(); i++){
//         myfile << lambdas[i] << std::endl;
//     }
//     myfile.close();

//     //Print rate to file
//     std::ofstream ratefile;
//     ratefile.open(obs.output_dir + "/ffs/rate.txt");
//     ratefile << "Format: (1) rate constant (2) initial flux (3+) transition probabilities" << std::endl;
//     ratefile << kAB << std::endl;
//     ratefile << phi_A0 << std::endl;
//     for(int i=0; i<(nint-1); i++){
//         ratefile << transition_probs[i] << std::endl;
//     }
//     ratefile.close();

//     //Print transition paths to file
//     int nevery = 10;
//     for(int i=0; i<n_successful; i++){
//         if (i % nevery==0){
//             std::string subdir = "ffs/path=" + std::to_string(i);
//             obs.open_h5md(sys, subdir);
//             for(int j=0; j<transition_paths[i].size(); j++){
//                 sys.particles = transition_paths[i][j];
//                 sys.time = j*sys.dt;
//                 obs.dump_h5md(sys, subdir);
//             }
//         }        
//     }
}

/******************************/
/*** Umbrella Sampling ********/
/******************************/
void LabBench::run_us_simulation(){}
