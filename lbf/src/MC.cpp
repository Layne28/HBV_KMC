#include "MC.hpp"

MC::MC(System &g, ParamDict &theParams, gsl_rng *&the_rg)
{

    if(theParams.is_key("ks0")) ks0 = std::stod(theParams.get_value("ks0"));
    if(theParams.is_key("kd0")) kd0 = std::stod(theParams.get_value("kd0"));
    if(theParams.is_key("do_vertex_only")) do_vertex_only = std::stod(theParams.get_value("do_vertex_only"));
    if(theParams.is_key("debug_sheet")) debug_sheet = std::stoi(theParams.get_value("debug_sheet"));

    //***Set RNG***
    rg = the_rg;

}

MC::~MC() 
{
}

void MC::run_relax(System &g, int nsteps)
{
    for (int n=0; n<nsteps; n++){
        move_vertices(g);
        g.update_boundary();
    }
}

// void MC::run_production(System &g, int nsteps)
// {

//     for (int n=0; n<nsteps; n++){
//         sweep(g);
//         double ee = 0;
//         if (sweep_count % (g.get_obs().freq_log) == 0 ) 
//         {
            
//             if (g.Nhe==6) recenter(g);
//             g.update_boundary();
//             g.check_odd_neigh();
//             ee = g.compute_energy();

//             dump_lammps_traj_dimers(g, sweep_count);

//             time(&timer2);
//             seconds = difftime(timer2, timer1);

//             dump_analysis(g, ofile, sweep_count, seed, seconds);
//             dump_lammps_data_file(g, 22222222);
//             //dump_lammps_traj_restart(g, sweep);
//             dump_lammps_data_dimers(g, 11111111);

//             dump_restart_lammps_data_file(g, sweep_count);
//         }

//         if (sweep_count % freq_out == 0 )
//         {
//             time(&timer2);
//             seconds = difftime(timer2, timer1);
//             cout << "###################################################################" << endl;
//             cout << " ################ RUN TIME " << seconds << " SECONDS ###############" << endl;
//             cout << " ############# SWEEP " << sweep_count << "##############" << endl;
//             cout << "######### FRAME " << frame << " ##############" << endl;
//             cout << "######### ENERGY " << ee << " ##############" << endl;
//             cout << "######### ENERGY PER DIMER " << 2 * ee / g.Nhe << " ##############" << endl;
//             cout << "#########  NHE " << g.Nhe << " #######################" << endl;
//             cout << "#########  NHESURF " << g.boundary.size() << " ##############" << endl;
//             cout << "#########  NVSURF " << g.boundaryv.size() << " ##############" << endl;
//             cout << "#########  NV_BONDSURF " << g.boundaryvbond.size() << " ##############" << endl;
//             cout << "#########  NV5 " << g.Nv5 << " ##############" << endl;
//             cout << "#########  MONOMER ADDED " << monomeradded << " ##############" << endl;
//             cout << "#########  MONOMER REMOVED " << monomerremoved << " ##############" << endl;
//             cout << "#########  DIMER ADDED " << dimeradded << " ##############" << endl;
//             cout << "#########  DIMER REMOVED " << dimerremoved << " ##############" << endl;
//             cout << "#########  no rate REMOVED " << deletednorate << " ##############" << endl;
//             cout << "#########  Surface bound " << binding << " ##############" << endl;
//             cout << "#########  Surface Unbound " << unbinding << " ##############" << endl;
//             cout << "#########  DrugAdded " << drugadded << " ##############" << endl;
//             cout << "#########  DrugRemoved " << drugremoved << " ##############" << endl;
//             cout << "#########  ND " << g.Nd << " ##############" << endl;
//             cout << "#########  TYPE CHANGED " << typechanged << " ##############" << endl;
//             cout << "#########  WEDGE FUSION " << wedgefusion << " ##############" << endl;
//             cout << "#########  WEDGE FISSION " << wedgefission << " ##############" << endl;
//             cout << "#########  FUSION " << fusion << " ##############" << endl;
//             cout << "#########  FISSION " << fission << " ##############" << endl;
//             cout << "#########  FUSION HALFEDGES " << g.fusionhe.size() << " ##############" << endl;
//             cout << "#########  WEDGE FUSION HALFEDGES " << g.fusionwedgehe.size() << " ##############" << endl;
//             cout << "#########  ALL NEIGH " << g.all_neigh << " ##############" << endl;
//             cout << "#########  Nboundary " << g.Nboundary << " ##############" << endl;
//             cout << "#########  Bound Triangle " << boundtri << " ##############" << endl;
//             cout << "#########  Acceptance vmove " << (1.0 * g.accepted_vmove) / (1.0 * (g.accepted_vmove + g.rejected_vmove)) << "#################" << endl;
//             cout << "#########  Nvlast " << g.Nvlast << " Nhelast " << g.Nhelast << " ################" << endl;
//             cout << "#########  T4 " << g.NCD_T4_in << " T3 " << g.NCD_T3_in << " ################" << endl;
//         cout << "#########  NCD_Hex" << g.NCD_Hex << "#################"<<endl;
//             cout << "#########  avgAddInterval "<<avgAddInterval<<" ###############"<<endl;
//         }

//         int cc = check_bind_triangle(g);
//         if (cc > 0)
//         {
//             cout << "bound triangle" << endl;
//             g.update_boundary();
//             boundtri += cc;
//         }

//         g.check_odd_neigh();
        
//         if (g.Nhe>70 && g.Nhe < minHE_update_neigh && sweep_count % 10000 == 0)
//         {
//             double thispace=float(g.Nhe-lastNhe)/10000.0; //pace of adding Nhe per sweep_count
//             cout << "thispace "<< thispace <<endl; 
            
//             avgpace=(npace*avgpace+thispace)/(npace+1.0); //average pace of adding 

//             avgAddInterval=int(pow(10,(-1 * int(floor(log10(avgpace))) ) ) );
//             cout << "avgpace" << avgpace << " avgAddInterval " <<avgAddInterval << endl;
            
//             npace++;
//             lastNhe=g.Nhe;
//         }

//         if (g.Nhe > minHE_update_neigh && sweep_count % 10000 == 0)
//         {
//             g.update_neigh();
//             if (g.find_overlap_all() < 0)
//             {
//                 cout << "error overlap" << endl;
//                 dump_lammps_data_dimers(g, 5555555);
//                 exit(-1);
//             }
//         }
//         //see if capsid is growing or it is stalled in mixed morphology
//         if (g.Nhe > minHE_update_neigh && sweep_count % (10*avgAddInterval) == 0)
//         {    
//             update_geometry_parameters(g);
//             if ( g.Nhe - lastNhe<=2 ){
//                 if ((g.Nhe >= 220 && g.NCD_T4_in >=26 && g.NCD_T3_in >= 3  && g.Nsurf > 10 ) || 
//                     (g.Nhe >= 160 && g.NCD_T4_in >= 3 && g.NCD_T3_in >=16  && g.Nsurf > 10) || 
//                     (g.Nhe >= 200 && g.NCD_T4_in >= 5 && g.NCD_T3_in >=5  && g.Nsurf > 10) )
//                     {
//             //         cout << "STOP for now - mixed morph" << endl;
//                         g.update_boundary();
//                         dump_lammps_traj_dimers(g, int(sweep_count));
//                         dump_lammps_data_dimers(g, 44444444);
//                         dump_lammps_data_dimers(g, 11111111);
//                         time(&timer2);
//                         seconds = difftime(timer2, timer1);
//                         dump_analysis(g, ofile, sweep_count, seed, seconds);
//             //           exit(-1);
//                     }
//             }
            
            
//             if (sweep_count % (100*avgAddInterval) == 0)
//             {
                
//                 if (g.NCD_T4_in>0 && g.NCD_T3_in>0 && abs( g.Nhe - lastNheGrowth)<=4 ){
//                     fprintf(stderr, "STOP for now - not growing\n");
//                     g.update_boundary();
//                     //dump_lammps_traj_dimers(g, int(sweep_count));
//                     dump_lammps_data_dimers(g, 333333333);
//                     dump_lammps_data_dimers(g, 11111111);
//                     dump_restart_lammps_data_file(g, sweep_count);
//                     time(&timer2);
//                     seconds = difftime(timer2, timer1);
//                     dump_analysis(g, ofile, sweep_count, seed, seconds);
//                 //  exit(-1);
//                 }
//                 lastNheGrowth = g.Nhe;
//             }
//             lastNhe = g.Nhe;            
//         }

//         if (sweep_count == 200000000)
//         {

//             fprintf(stderr, "STOP for now - too long\n");
//             g.update_boundary();
//             dump_lammps_traj_dimers(g, int(sweep_count));
//             //dump_lammps_traj_restart(g, int(sweep_count));
//             dump_lammps_data_dimers(g, 77777777);
//             dump_restart_lammps_data_file(g, sweep_count);
//             time(&timer2);
//             seconds = difftime(timer2, timer1);
//             dump_analysis(g, ofile, sweep_count, seed, seconds);
//             exit(-1);
//         }

//         if (g.Nhe >= 310 || g.Nv >= 65)
//         {

//         // fprintf(stderr, "STOP for now - too large\n");
//             g.update_boundary();
//             dump_lammps_traj_dimers(g, int(sweep_count));
//             dump_lammps_data_dimers(g, 88888888);
//             dump_lammps_data_dimers(g, 11111111);
//             dump_restart_lammps_data_file(g, sweep_count);
//             time(&timer2);
//             seconds = difftime(timer2, timer1);
//             dump_analysis(g, ofile, sweep_count, seed, seconds);
//         // exit(-1);
//         }
//         sweep_count++;
//     }

//     //dump_lammps_traj_restart(g, int(sweep_count));
//     dump_lammps_traj_dimers(g, int(sweep_count));

//     //if (g.Nboundary == 1)
//     dump_restart_lammps_data_file(g, sweep_count);
//     //equilibrating final structure
//     for (int rstep = 0; rstep < (10 * freq_log); rstep++)
//     {
//         move_vertices(g, r);
//         g.update_boundary();
//         int ind1 = gsl_rng_uniform_int(r, g.Nhe);
//         int e1 = g.he[ind1].id;
//         int x = -1;

//         //change edge type
//         x = attempt_change_edge_type(g, e1, r);
//         if (x >= 0)
//             typechanged++;

//         if (sweep_count % freq_vis == 0)
//         {

//             dump_lammps_traj_dimers(g, int(sweep_count));
//         }
//         if (sweep_count % freq_out == 0)
//         {
//             time(&timer2);
//             seconds = difftime(timer2, timer1);
//             cout << "###################################################################" << endl;
//             cout << " ################  RUN TIME " << seconds << " SECONDS ###############" << endl;
//             cout << " ############# SWEEP " << sweep_count << "##############" << endl;

//             double ee = g.compute_energy();
//             cout << "######### ENERGY " << ee << " ##############" << endl;
//             cout << "######### ENERGY PER DIMER " << 2 * ee / g.Nhe << " ##############" << endl;
//         }
//         if (sweep_count % freq_log == 0)
//         {

//             time(&timer2);
//             seconds = difftime(timer2, timer1);

//             dump_analysis(g, ofile, sweep_count, seed, seconds);
//         }
//         sweep_count++;
//     }
// }

void MC::sweep(System &g)
{
    double ps_attempt = 0;
    int ss = 0;
    int ssadd = 0;

    /**************************************************/
    /*** Attempt MC moves at prescribed frequencies ***/
    /**************************************************/

    /*** Vertex relaxation ***/
    move_vertices(g);

    if (do_vertex_only==1) return;

    /*** Conformational change ***/
    if (g.Nhe==6){
        int ind = gsl_rng_uniform_int(rg, g.boundary.size());
            //cout <<"ind id" << ind<<endl;
        int hh = g.boundary[ind];
        int x = attempt_change_edge_type_tri(g, hh);
        if (x >= 0)
            typechanged++; 
    }
    else{
        for (int nc = 0; nc < g.Nhe/2; nc++)
        {
            int ind1 = gsl_rng_uniform_int(rg, g.Nhe);
            int e1 = g.he[ind1].id;
            int x = -1;
            x = attempt_change_edge_type(g, e1);
            if (x >= 0)
                typechanged++;
        }
    }

    /*** Add monomers/dimers ***/
    ps_attempt = ks0 * g.Nsurf;
    if (gsl_rng_uniform(rg) < ps_attempt)
    {
        //std::cout << "Attempt to add monomer/dimer" << std::endl;
        //std::cout << ps_attempt << std::endl;
        int ind = gsl_rng_uniform_int(rg, g.boundary.size());
        int e = g.boundary[ind];
        if (g.check_inside_overlap(e) > 0)
        {
            if(debug_sheet==1){
                std::cout << "Attempt to add monomer/dimer" << std::endl;
                //std::cout << ps_attempt << std::endl;
            }

            //TESTING
            //ssadd = -1;
            ssadd = attempt_add_monomer_dimer(g, e);
            //std::cout << ssadd << std::endl;
            if (ssadd > 1){
                dimeradded++;
                if(debug_sheet==1){
                    std::cout << "Dimer added" << std::endl;
                }
            }
            else if (ssadd > 0){
                monomeradded++;
                if(debug_sheet==1){
                    std::cout << "Monomer added" << std::endl;
                }
            }
            ssadd = -1;
            g.update_boundary();
        }
    }

    //Attempt to add drug-bound dimer
    if (gsl_rng_uniform(rg) < ps_attempt)
    {
        int ind = gsl_rng_uniform_int(rg, g.boundary.size());
        int e = g.boundary[ind];
        if (g.check_inside_overlap(e) > 0)
        {
            ssadd = attempt_add_monomer_dimer_drug(g, e);
            if (ssadd > 1){
                dimeradded++;
                drugadded++;
            }
            ssadd = -1;
            g.update_boundary();
        }
    }


    if (g.Nhe>6) {

        /*** Remove monomers/dimers ***/
        ps_attempt = ks0 * g.Nsurf;
        if (gsl_rng_uniform(rg) < ps_attempt)
        {
            //std::cout << "Attempt to remove monomer/dimer" << std::endl;
            //std::cout << ps_attempt << std::endl;
            //if (sweep > 0 && g.Nhe > 6)
            if (g.Nhe > 6)
            {
                //std::cout << "Attempt to remove monomer/dimer 2" << std::endl;
                int ind = gsl_rng_uniform_int(rg, g.boundary.size());
                int e = g.boundary[ind];
                if (g.no_bond_boundary(e) > 0)
                {
                    ss = attempt_remove_monomer_dimer(g, e);
                    if (ss > 1)
                        dimerremoved++;
                    else if (ss > 0)
                        monomerremoved++;
                    ss = 0;
                    g.update_boundary();
                }
            }
        }

        //Attempt to remove drug-bound dimer
        if (gsl_rng_uniform(rg) < ps_attempt)
        {
            if (g.Nhe > 6)
            {
                int ind = gsl_rng_uniform_int(rg, g.boundary.size());
                int e = g.boundary[ind];
                if (g.no_bond_boundary(e) > 0)
                {
                    ss = attempt_remove_monomer_dimer_drug(g, e);
                    ss = 0;
                    if (ss > 1)
                    {
                        dimerremoved++;
                        drugremoved++;
                    }
                    //else if (ss > 0)
                        //  monomerremoved++;
                    ss = 0;
                    g.update_boundary();
                }
            }
        }

        if (g.Nhe > 15) 
        {
            int cc = check_bind_triangle(g);
            if (cc > 0)
            {
                cout << "bound triangle" << endl;
                g.update_boundary();
                boundtri += cc;
            }
            else 
            {
                /*** Bind wedge ***/
                //if (gsl_rng_uniform(r) < pb_attempt){
                if(g.boundary.size()>0) {

                    int ind = gsl_rng_uniform_int(rg, g.boundary.size());
                    int hh = g.boundary[ind];
                    int tt = -1;
                    if (g.no_bond_boundary(hh) > 0)
                    {
                        tt = attempt_bind_wedge_dimer(g, hh);
                        if (tt > 0)
                            binding++;
                        tt = -1;
                    }
                    g.update_boundary();
                }

                /*** Unind wedge ***/
                if (g.boundaryvbond.size() > 0) // &&  gsl_rng_uniform(r) < pb_attempt)
                {
                    int ind = gsl_rng_uniform_int(rg, g.boundary.size());
                    int hh = g.boundary[ind];
                    if ((g.is_bond_in_boundary(hh) > 0) || (g.is_bond_out_boundary(hh) > 0))
                    {
                        //cout <<" trying unbind vv is " << vv <<endl;
                        int tt = attempt_unbind_wedge_dimer(g, hh);
                        if (tt > 0)
                            unbinding++;
                        tt = -1;
                        g.update_boundary();
                    }
                }

                if (g.Nhe > minhe_fission && g.Nsurf > 3){ // dont try if only last triangle is open

                    int movetype = gsl_rng_uniform_int(rg, 4);
                    switch (movetype)
                    {
            
                    case 0:
                        /*** Wedge fusion ***/
                        if (g.all_neigh > 0 )
                        {
                        //double kf0=1;
                        //double pf_attempt=kf0*g.Nsurf;
                        //if (gsl_rng_uniform(r) < pf_attempt)
                        //{
                        int ff = attempt_wedge_fusion(g);
                        g.update_boundary();
                        if (ff > 0)
                            wedgefusion++;
                        ff = 0;
                        //}
                        }
                        break;

                    case 1:
                        /*** Wedge fission ***/
                        {
                            //pf_attempt=kf0; //g.Nsurf
                            //if (gsl_rng_uniform(r) < pf_attempt)
                            //{
                            int ff = attempt_wedge_fission(g);
                            g.update_boundary();
                            if (ff > 0)
                                wedgefission++;
                            ff = 0;
                            //}
                        }
                        break;

                    case 2:
                        /*** Fusion ***/
                        if (g.all_neigh > 0 )
                        {
                            //double kf0=1;
                            //double pf_attempt=kf0;//g.Nsurf;
                            //if (gsl_rng_uniform(r) < pf_attempt)
                            //{
                            int ff = attempt_fusion(g);
                            g.update_boundary();
                            if (ff > 0)
                                fusion++;
                            ff = 0;
                        }
                        break;

                    case 3:
                        /*** Fission ***/
                        
                        {
                        //pf_attempt=kf0*.25; //g.Nsurf
                        //if (gsl_rng_uniform(r) < pf_attempt)
                        //{
                        int ff = attempt_fission(g);
                        //std::cout << "done attempting fission. updating boundary..." << std::endl;
                        g.update_boundary();
                        /*for (vector<VTX>::iterator it2 = g.v.begin(); it2 != g.v.end(); ++it2){
                                std::cout << "it2->vid: " << it2->vid << std::endl;
                        }*/
                        for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
                        {
                            //std::cout << "it->vout: " << it->vout << std::endl;
                            int vout_in_v = 0;
                            for (vector<VTX>::iterator it2 = g.v.begin(); it2 != g.v.end(); ++it2){
                                //if(it->vout==776732) std::cout << "check: it->vout: " << it->vout << " it2->vid: " << it2->vid << std::endl;
                                if(it->vout==it2->vid) vout_in_v = 1;
                            }
                            if (vout_in_v==0){
                                std::cout << "WARNING: after update boundary: vout (" << it->vout << ") is not in the vextex id list!" << std::endl;
                            }
                        }
                        if (ff > 0)
                            fission++;
                        ff = 0;
                        //}
                        }
                        break;

                    }
                }


            }

        }
    }
    
    /*** Drug moves ***/
    
    double d_attempt = kd0 * g.Nhe;

    /*** Drug binding ***/
    if (gsl_rng_uniform(rg) < d_attempt)
    {

        int ind = gsl_rng_uniform_int(rg, g.Nhe);
        //fprintf(stderr, "Attempt add drug, Nd=%d index=%d\n", g.Nd, ind );
        int e = g.he[ind].id;
        if ((g.he[ind].type == 0 || g.he[ind].type == 3))
        { //no_bond_boundary(e)>0) {
            //fprintf(stderr, "Attempt delete monomer\n" );
            ss = attempt_add_drug(g, e);
            if (ss > 0)
            { //cout << "drug added " << endl;
                drugadded++;
                //g.Nd++;
            }

            //g.update_index();

            ss = 0;
        }
    }
    g.update_boundary();
    

    /*** Drug unbinding ***/
    if (gsl_rng_uniform(rg) < d_attempt && g.Nd > 0)
    {

        int ind = gsl_rng_uniform_int(rg, g.Nhe);
        //fprintf(stderr, "Attempt delete drug, Nd=%d index=%d\n", g.Nd, ind );
        int e = g.he[ind].id;
        if (g.he[ind].type == 0 || g.he[ind].type == 3)
        { //no_bond_boundary(e)>0) {
            //fprintf(stderr, "Attempt delete monomer\n" );
            ss = attempt_remove_drug(g, e);
            if (ss > 0)
            { //cout << "drug removed " << endl;
                drugremoved++;
                //g.Nd--;
            }
            
            ss = 0;
        }
    }
    
    g.update_boundary();
                
}

void MC::move_vertices(System &g)
{
    /*** Attempt to move all vertices ***/
    //std::cout << "inin move_vertices"<<endl;
    if(no_del_vertex_counter>0) std::cout << "moving vertices after not deleting vertex" << std::endl;
    double *newv = new double[3];
    double *oldv = new double[3];

    int overlapflag = -1;
    double e1 = 0, e2 = 0, de = 0;
    //double e11=0;

    for (unsigned int i = 0; i < g.v.size(); i++)
    {
        //g.update_normals();
        //g.validate_surface();
        /* to test overlap */
        overlapflag = -1;
        // for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
        //{
        //    g.update_half_edge(it->id);

        //}
        //g.update_normals();
        int ind = gsl_rng_uniform_int(rg, g.v.size());
        /* test !!!! */
        //int ind = i;
        //std::cout << " move " << ind <<endl;
        /* save initial energy */
        e1 = g.vertex_energy(g.v[ind].vid);
        /* test !!!!*/
        //e11=g.compute_energy();
        /* save old vertex coordinates */
        oldv[0] = g.v[ind].co[0];
        oldv[1] = g.v[ind].co[1];
        oldv[2] = g.v[ind].co[2];

        /* move the vertex to coordinates newv */
        g.move_v_epsilon(g.xi, oldv, newv, rg);
        if(no_del_vertex_counter>0) std::cout << "moved to new coords" << std::endl;

        /* update the coordinates of the vertex  */
        g.v[ind].co[0] = newv[0];
        g.v[ind].co[1] = newv[1];
        g.v[ind].co[2] = newv[2];

        /* update vertex geometry , normals , excluders */
        if(no_del_vertex_counter>0) std::cout << "about to update" << std::endl;
        g.update_geometry_vertex(ind);
        g.update_normals_vertex(ind);
        g.update_excluder_top_vertex(ind);
        if(no_del_vertex_counter>0) std::cout << "done" << std::endl;

        //for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
        //{
        //    g.update_half_edge(it->id);

        //}
        //g.update_normals();
        /* save new energy */
        e2 = g.vertex_energy(g.v[ind].vid);

        /* check overlap */
        if(no_del_vertex_counter>0) std::cout << "checking overlap" << std::endl;
        if (g.Nhe > 20)
        {
            if (g.check_overlap_g(g.v[ind].vid) < 0)
            {
                overlapflag = 1;
            }
        }
        if(no_del_vertex_counter>0) std::cout << "done" << std::endl;

        de = e2 - e1;
        /* test !!!!*/
        //std::cout << "de                          "<< de  << " e2 "<< e2  << " e1 "<< e1 <<endl;
        //std::cout << "g.compute_energy()-e11      "<< g.compute_energy()-e11  << " g.compute_energy() "<< g.compute_energy() << " e11 "<< e11 <<endl;
        //std::cout << "de-(g.compute_energy()-e11) " << de-(g.compute_energy()-e11) << endl;

        /*if (abs(de-(g.compute_energy()-e11)>0.0000000001)) { 
                
                std::cout << "de                          "<< de  << " e2 "<< e2  << " e1 "<< e1 <<endl; 
                std::cout << "g.compute_energy()-e11      "<< g.compute_energy()-e11  << " g.compute_energy() "<< g.compute_energy() << " e11 "<< e11 <<endl;
                std::cout << "de-(g.compute_energy()-e11) " << de-(g.compute_energy()-e11) << endl; 

                std::cout << "vindex is " <<  ind <<endl;
                std::cout << "hein.size() " << g.v[ind].hein.size() <<endl;
                for (vector<int>::iterator ithe = g.v[ind].hein.begin(); ithe != g.v[ind].hein.end(); ithe++)
                {
                    int heindex0=g.heidtoindex[*ithe];
                    std::cout <<"heindex " <<heindex0<<endl;
                    std::cout << g.he[heindex0].n[0] <<" "<<g.he[heindex0].n[1] <<" "<<g.he[heindex0].n[2] <<" "<<endl;
                    std::cout << "next " << g.he[heindex0].nextid <<endl;
                    std::cout << "prev " << g.he[heindex0].previd <<endl;
                    if (g.is_surface(g.he[heindex0].opid)>0) 
                    {
                        int opindex0= g.heidtoindex[g.he[heindex0].opid];
                        std::cout <<"opindex " <<opindex0<<endl;
                        std::cout << g.he[opindex0].n[0] <<" "<<g.he[opindex0].n[1] <<" "<<g.he[opindex0].n[2] <<" "<<endl;
                        }
                }
                //g.validate_surface();
                std::exit(-1);
                
            }*/

        double crit = exp((-de) / g.T);

        /* reject the move if crit is not met or overlap */
        if (gsl_rng_uniform(rg) > crit || overlapflag == 1)
        {
            /* update vertex coordinates to oldv */
            g.v[ind].co[0] = oldv[0];
            g.v[ind].co[1] = oldv[1];
            g.v[ind].co[2] = oldv[2];

            /* update geometry normals excluders */
            g.update_geometry_vertex(ind);
            g.update_normals_vertex(ind);
            g.update_excluder_top_vertex(ind);
            //std::cout << " not accepted " << endl;
            g.rejected_vmove++;
        }
        else {
            g.accepted_vmove++;
            //############ move dist
            //ofstream myfile;
            //myfile.open ("movedist.txt", ios::out | ios::app);
            //myfile <<dis<<endl;
            //myfile.close();
            //############
        }
    
    }

    /* after all moves update surface neighbors */

    if(no_del_vertex_counter>0) std::cout << "updating surface neighbors" << std::endl;
    for (vector<int>::iterator it = g.boundaryv.begin(); it != g.boundaryv.end(); ++it)
    {
        g.update_neigh_vertex_and_neigh(*it);
    }
    if(no_del_vertex_counter>0) std::cout << "done" << std::endl;

    if(no_del_vertex_counter>0){ std::cout << "deleting newv, oldv" << std::endl;
    std::cout << "newv: " << newv[0] << " " << newv[1] << " " << newv[2] << std::endl;
    std::cout << "oldv: " << oldv[0] << " " << oldv[1] << " " << oldv[2] << std::endl;
}
    //delete[] newv;
    //delete[] oldv;
    if(no_del_vertex_counter>0) std::cout << "done" << std::endl;
}


int MC::move_one_vertex(System &g, int vid0)
{
    //std::cout << "inin move_vertices"<<endl;
    double *newv = new double[3];
    double *oldv = new double[3];

    int overlapflag = -1;
    double e1 = 0, e2 = 0, de = 0;
    
    overlapflag = -1;

    int ind = g.vidtoindex[vid0];
    /* test !!!! */
    //int ind = i;
    //std::cout << " move " << ind <<endl;
    /* save initial energy */
    e1 = g.vertex_energy(g.v[ind].vid);
    /* test !!!!*/
    //e11=g.compute_energy();
    /* save old vertex coordinates */
    oldv[0] = g.v[ind].co[0];
    oldv[1] = g.v[ind].co[1];
    oldv[2] = g.v[ind].co[2];

    /* move the vertex to coordinates newv */
    g.move_v_epsilon(g.xi, oldv, newv, rg);

    /* update the coordinates of the vertex  */
    g.v[ind].co[0] = newv[0];
    g.v[ind].co[1] = newv[1];
    g.v[ind].co[2] = newv[2];

    /* update vertex geometry , normals , excluders */
    g.update_geometry_vertex(ind);
    g.update_normals_vertex(ind);
    g.update_excluder_top_vertex(ind);

        
    e2 = g.vertex_energy(g.v[ind].vid);

    /* check overlap */
    if (g.Nhe > 20)
    {
        if (g.check_overlap_g(g.v[ind].vid) < 0)
        {
            overlapflag = 1;
        }
    }

    de = e2 - e1;
    

    double crit = exp((-de) / g.T);

        /* reject the move if crit is not met or overlap */
    if (gsl_rng_uniform(rg) <crit && overlapflag != 1)
    {
        //############ move dist
        //ofstream myfile;
        //myfile.open ("movedist.txt", ios::out | ios::app);
        //myfile <<dis<<endl;
        //myfile.close();
        //############
        g.update_neigh_vertex_and_neigh(vid0);
        return(1);
    }
    else{
        /* update vertex coordinates to oldv */
        g.v[ind].co[0] = oldv[0];
        g.v[ind].co[1] = oldv[1];
        g.v[ind].co[2] = oldv[2];

        /* update geometry normals excluders */
        g.update_geometry_vertex(ind);
        g.update_normals_vertex(ind);
        g.update_excluder_top_vertex(ind);
        //std::cout << " not accepted " << endl;
        
    }
    g.update_neigh_vertex_and_neigh(vid0);
    return(-1);
    
}


int MC::attempt_add_monomer_dimer(System &g, int heid0) //!!! Should update with pre_oipen wedge
{
    if(debug_sheet==1){
        std::cout << "in attempt_add_monomer-dimer" <<endl;
    }
    g.update_normals();

    if (g.is_boundary(heid0) < 0)
    {
        std::cout << " cannot add not on the !" << endl;
        std::exit(-1);
    }

    int etypenew1 = -1;
    int etypenew2 = -1;

    int heindex0 = g.heidtoindex[heid0];
    int etypeheid0 = g.he[heindex0].type;

    int bi = g.he[heindex0].boundary_index;
    //temp test
    if (bi == -1)
    {
        std::cout << "wrong boundary index in add monomer dimer" << endl;
        std::exit(-1);
    }

    //should not add if wedge fusion is possible
    g.update_fusion_pairs_he();

    if (g.he[heindex0].next_wedge_fusion_heid != -1 || g.he[heindex0].prev_wedge_fusion_heid != -1)
    {
        return -1;
    }
    //monomer with next

    //if (g.he[heindex0].nextid!=-1 && g.he[heindex0].previd!=-1) {std::cout <<" wrong not on boundary" <<endl; std::exit(-1); }
    int xid = -1;

    int vid0 = g.is_bond_out_boundary(heid0);
    int vin0 = g.is_bond_in_boundary(heid0);

    if (vid0 >= 0 && g.v[g.vidtoindex[g.he[heindex0].vin]].hein.size() < 6 && g.v[g.vidtoindex[g.he[g.heidtoindex[g.he[heindex0].nextid]].vin]].hein.size() < 6)
    {
        // no nead for this check bind triangle will take care of this
        /*if (g.he[g.heidtoindex[g.he[heindex0].nextid]].next_wedge_fusion_heid != -1 || g.he[g.heidtoindex[g.he[heindex0].nextid]].prev_wedge_fusion_heid != -1)
        {
            //std::cout << "g.he[heindex0].next_wedge_fusion_heid "<<g.he[heindex0].next_wedge_fusion_heid<<endl;
            //std::cout << "g.he[g.heidtoindex[g.he[heindex0].previd]].next_wedge_fusion_heid " << g.he[g.heidtoindex[g.he[heindex0].previd]].next_wedge_fusion_heid<<endl;
            //std::cout << "g.he[g.heidtoindex[g.he[heindex0].nextid]].next_wedge_fusion_heid" <<g.he[g.heidtoindex[g.he[heindex0].nextid]].next_wedge_fusion_heid <<endl;
            //std::exit(-1);
            return -1;
        }*/

        //std::cout << "00 add monomer with next"<<endl;
        /* ToDo unify add with next and add with prev */
        //int thisprev_heid=heid0;
        //int thisnext_heid=g.he[heindex0].nextid;
        xid = g.he[heindex0].nextid;
        int xidindex = g.heidtoindex[xid];

        if (g.is_bond_in_boundary(xid) != vid0 || g.is_bond_out_boundary(xid) != -1)
        {
            std::cout << "wrong bound on boundary" << endl;
            std::exit(-1);
        }
        else
        {
            //std::cout << "try add monomer" <<endl;
            int etypenew = -1;

            int etypexid = g.he[g.heidtoindex[xid]].type;
            // SELECTED TYPES

           /* if ((etypeheid0 == 2) && ((etypexid == 0) || (etypexid == 3)))
            {
                etypenew = 1;
            }
            else if (((etypeheid0 == 3) || (etypeheid0 == 0)) && (etypexid == 1))
            {
                etypenew = 2;
            }
            else if (((etypeheid0 == 3) && (etypexid == 3)))
            {
                etypenew = 3;
            }
            else if (((etypeheid0 == 0) && (etypexid == 0)) || ((etypeheid0 == 1) && (etypexid == 2)))
            {
                etypenew = 0;
            }
            else if ((etypeheid0 == 1) && (etypexid == 2))
            {
                if (gsl_rng_uniform(r) < 0.5)
                {
                    etypenew = 0;
                }
                else
                {
                    etypenew = 3;
                }
            }
            else
            {
                etypenew = gsl_rng_uniform_int(r, 4);
            }*/
	        etypenew=gsl_rng_uniform_int(rg, 4);
            double e1 = g.bend_energy(heindex0) + g.bend_energy(xidindex);
            // TYPES BASED ON Concentration
            /*if (gsl_rng_uniform(r) < g.cdProb) {
                    if (gsl_rng_uniform(r) < 0.5) {etypenew=0;}
                    else {etypenew=3; }
                } 
                else{
                    if (gsl_rng_uniform(r) < 0.5) {etypenew=1;}
                    else {etypenew=2; } 

                }  */
            // etypenew=gsl_rng_uniform_int(r,4);

            bool drug0 = 0;
            //bool drug1=0;
            /*if ((etypenew==3) ||(etypenew==0) ) {
                if (gsl_rng_uniform(r) < g.drugProb) {drug0=1; }
                if (gsl_rng_uniform(r) < g.drugProb) {drug1=1;  }
                }*/
            //std::cout << "g.he[heindex0].din" << g.he[heindex0].din<<endl;
            double gbb = g.find_dg(etypenew, etypeheid0, g.he[heindex0].din);
            gbb += g.find_dg(etypexid, etypenew, drug0);
            //std::cout << "add_monomer with prev" <<endl;

            //double e1=g.bend_energy(heindex0) + g.bend_energy(xidindex);
            int x = g.add_monomer(heid0, xid, etypenew);
            //int Nhelast2index = g.heidtoindex[g.Nhelast - 2];
            if (x > 0)
            {
                //double de=g.monomer_energy(g.Nhelast-1);
                int voutid = g.he[heindex0].vout;
                int vinid = g.he[xidindex].vin;
                vector<int> vecupdate;
                vecupdate.push_back(vinid);
                vecupdate.push_back(voutid);
                g.update_half_edge(heid0);
                g.update_half_edge(g.he[heindex0].opid);
                g.update_half_edge(xid);
                g.update_half_edge(g.he[xidindex].opid);
                g.update_half_edge(g.Nhelast - 1);
                //g.he[g.heidtoindex(g.Nhelast-1)].din=drug1;
                g.update_half_edge(g.Nhelast - 2);
                //g.he[g.heidtoindex[g.Nhelast-2]].din=drug0;
                g.he[g.heidtoindex[g.Nhelast - 1]].boundary_index = bi; // update boundary index
                g.update_boundary();

                //g.update_neigh_vertex(vinid);
                //g.update_neigh_vertex(voutid);
                g.update_geometry_vertex(g.vidtoindex[vinid]);
                g.update_normals_vertex(g.vidtoindex[vinid]);
                g.update_excluder_top_vertex(g.vidtoindex[vinid]);
                g.update_geometry_vertex(g.vidtoindex[voutid]);
                g.update_normals_vertex(g.vidtoindex[voutid]);
                g.update_excluder_top_vertex(g.vidtoindex[voutid]);

                double de = g.stretch_energy(g.heidtoindex[g.Nhelast - 2]);

                de += g.dimer_bend_energy(g.heidtoindex[g.Nhelast - 2]) + g.dimer_bend_energy(xidindex);
                de += g.bend_energy(heindex0) + g.bend_energy(xidindex)-e1; // g.dimer_bend_energy(heindex0);
                //de+=  ; //g.monomer_energy(heid0);
                //std::cout << " de is " << de <<endl;
                //double e2=g.compute_energy();
                //de += g.gb * 3 - g.mu;
                //gbb=gb0next+gb0prev;
                de += gbb - g.mu[etypenew];
                if(debug_sheet==1){
                    std::cout << "Binding free energy: " << gbb << ", chemical potential: " << g.mu[etypenew] << ", total: " << gbb - g.mu[etypenew] << std::endl;
                    if(de<0){
                        std::cout << "Favorable free energy! Found critical nucleus! Press any button to continue..." << std::endl;
                        std::cin.get();
                    }
                }
                /*int vinid=g.he[g.heidtoindex[g.Nhelast - 1]].vin;
                    int vindex0=g.vidtoindex[vinid];
                    g.update_neigh_vertex(vinid);
                    if (g.v[vindex0].vneigh.size() > 0)
                    {
                        for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                        {	
                            g.update_neigh_vertex(*it);
                        }
                    }
                    int voutid=g.he[g.heidtoindex[g.Nhelast - 1]].vout;
                    vindex0=g.vidtoindex[voutid];
                    if (g.v[vindex0].vneigh.size() > 0)
                    {
                        for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                        {	
                            g.update_neigh_vertex(*it);
                        }
                    }*/

                //double crit = 2.*g.z*g.K*g.K*g.K*exp((-de)/g.T);

                double crit = exp((-de) / g.T) / 2;
                int overlapflag = -1;

                /* check overlap */
                for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
                {
                    if (g.check_overlap_g(*it) < 0)
                        overlapflag = 1;
                }

                if (gsl_rng_uniform(rg) < crit && overlapflag == -1)
                {
                    //ToDo
                    // update next_prevoius surface

                    g.he[g.heidtoindex[g.Nhelast - 1]].boundary_index = bi;
                    g.v[g.vidtoindex[vid0]].doubleboundary = -1;
                    //no nead for this -- it should be closed by check bind triangle
                    /*if (g.he[g.heidtoindex[g.he[xidindex].nextid_boundary]].vout == g.he[g.heidtoindex[g.he[heindex0].previd_boundary]].vin)
                    {
                        g.set_prev_next(g.he[heindex0].previd_boundary, g.he[xidindex].nextid_boundary, g.Nhelast - 1);
                        g.set_prev_next(g.he[xidindex].nextid_boundary, g.Nhelast - 1, g.he[heindex0].previd_boundary);
                        g.set_prev_next(g.Nhelast - 1, g.he[heindex0].previd_boundary, g.he[xidindex].nextid_boundary);
                        g.Nboundary--;

                    }*/

                    //std::cout << "00 added monomer with next accepted"<<endl;
                    vecupdate.clear();
                    return 1;
                }
                /*deleting the added monomer */
                else
                {
                    /* ToDo update this with delete monomer */

                    int nextid0 = g.he[g.heidtoindex[g.Nhelast - 2]].nextid;
                    int previd0 = g.he[g.heidtoindex[g.Nhelast - 2]].previd;
                    int nextidboundary0 = g.he[g.heidtoindex[g.Nhelast - 1]].nextid_boundary;
                    int previdboundary0 = g.he[g.heidtoindex[g.Nhelast - 1]].previd_boundary;
                    //int xidindex=g.heidtoindex[xid];
                    if (nextid0 == -1 || previd0 == -1)
                    {
                        std::cout << "not accepted in add monomer!" << endl;
                        std::exit(-1);
                    }
                    if (g.delete_edge(g.Nhelast - 1) > 0)
                    {

                        g.he[heindex0].previd = -1;
                        g.he[xidindex].nextid = -1;

                        g.update_half_edge(heid0);
                        g.update_half_edge(g.he[heindex0].opid);
                        g.update_half_edge(xid);
                        g.update_half_edge(g.he[xidindex].opid);

                        // update boundary_index
                        g.he[heindex0].boundary_index = bi;
                        g.he[xidindex].boundary_index = bi;
                        g.set_prev_next_boundary(xid, nextidboundary0);
                        g.set_prev_next_boundary(previdboundary0, heid0);
                        g.set_prev_next_boundary(heid0, xid);
                        //std::cout << "updating index after deleting monomer (l 1077)" << std::endl;
                        g.update_index();
                        //std::cout << " MONOMER REMOVED AFTER addition xid1" <<endl;

                        for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
                        {
                            g.update_neigh_vertex(*it);
                        }

                        // no need to update nex previous surface
                        //std::cout << "00 added monomer with next removed"<<endl;
                        vecupdate.clear();
                        return -1;
                    }
                    else
                    {
                        std::cout << " could not delete_monomer HERE 4444";
                        std::exit(-1);
                    }
                }
            }
            else
            {

                //std::cout << "could not add monomer" <<endl;
                return -1;
            }
        }
    }

    //adding monomer  with previous
    else if (vin0 >= 0 && g.v[g.vidtoindex[g.he[heindex0].vout]].hein.size() < 6 && g.v[g.vidtoindex[g.he[g.heidtoindex[g.he[heindex0].previd]].vin]].hein.size() < 6)
    {
        //std::cout << "00 add monomer with previous"<<endl;
        xid = g.he[heindex0].previd;
        int xidindex = g.heidtoindex[xid];
        double e1 = g.bend_energy(heindex0) + g.bend_energy(xidindex);
        if (g.he[g.heidtoindex[g.he[heindex0].previd]].next_wedge_fusion_heid != -1 || g.he[g.heidtoindex[g.he[heindex0].previd]].prev_wedge_fusion_heid != -1)
        {
            return -1;
        }

        if (g.is_bond_out_boundary(xid) != vin0 || g.is_bond_in_boundary(xid) != -1)
        {
            std::cout << "wrong bound on boundary" << endl;
            std::exit(-1);
        }
        else
        {

            int etypexid = g.he[g.heidtoindex[xid]].type;
            int etypenew = -1;
            //std::cout  << "add_monomer with previous" <<endl;
            //BASED on types

            //get_monomer_etype(etypexid,etypeheid0,etypenew)
           /* if ((etypexid == 2) && ((etypeheid0 == 0) || (etypeheid0 == 3)))
            {
                etypenew = 1;
            }
            else if (((etypexid == 3) || (etypexid == 0)) && (etypeheid0 == 1))
            {
                etypenew = 2;
            }
            else if (((etypexid == 3) && (etypeheid0 == 3)))
            {
                etypenew = 3;
            }
            else if (((etypexid == 0) && (etypeheid0 == 0)) || ((etypexid == 1) && (etypeheid0 == 2)))
            {
                etypenew = 0;
            }
            else if ((etypexid == 1) && (etypeheid0 == 2))
            {
                if (gsl_rng_uniform(r) < .5)
                {
                    etypenew = 0;
                }
                else
                {
                    etypenew = 3;
                }
            }
            else
            {
                etypenew = gsl_rng_uniform_int(r, 4);
            }*/
	        etypenew=gsl_rng_uniform_int(rg, 4);
            /* if (gsl_rng_uniform(r) < g.cdProb) {
                    if (gsl_rng_uniform(r) < 0.5) {etypenew=0;}
                    else {etypenew=3; }
                } 
                else{
                    if (gsl_rng_uniform(r) < 0.5) {etypenew=1;}
                    else {etypenew=2; } 
                }*/
            //CHOOSING RANDOM
            //etypenew=gsl_rng_uniform_int(r,4);

            //else { std::cout << }
            bool drug0 = 0;
            //bool drug1=0;
            /*if ((etypenew==3) ||(etypenew==0) ) {
                    if (gsl_rng_uniform(r) < g.drugProb) {drug0=1;}
                    if (gsl_rng_uniform(r) < g.drugProb) {drug1=1;}
                }*/
            double gbb = g.find_dg(etypenew, etypexid, g.he[g.heidtoindex[xid]].din);
            gbb += g.find_dg(etypeheid0, etypenew, drug0);

            //double e1=g.bend_energy(heindex0) + g.bend_energy(xidindex);
            int x = g.add_monomer(xid, heid0, etypenew);
            if (x > 0)
            {
                int vinid = g.he[heindex0].vin;
                int voutid = g.he[xidindex].vout;
                vector<int> vecupdate;
                vecupdate.push_back(vinid);
                vecupdate.push_back(voutid);
                //double de = (g.stretch_energy(g.heidtoindex[g.Nhelast-2]) + g.bend_energy(g.heidtoindex[g.Nhelast-2]));
                g.update_half_edge(heid0);
                g.update_half_edge(g.he[g.heidtoindex[heid0]].opid);
                g.update_half_edge(xid);
                g.update_half_edge(g.he[g.heidtoindex[xid]].opid);
                g.update_half_edge(g.Nhelast - 1);
                //g.he[g.heidtoindex(g.Nhelast-1)].din=drug1;
                g.update_half_edge(g.Nhelast - 2);
                //g.he[g.heidtoindex[g.Nhelast-2]].din=drug0;
                g.update_boundary();

                
                //g.update_neigh_vertex(vinid);
                //g.update_neigh_vertex(voutid);
                g.update_geometry_vertex(g.vidtoindex[vinid]);
                g.update_normals_vertex(g.vidtoindex[vinid]);
                g.update_excluder_top_vertex(g.vidtoindex[vinid]);
                g.update_geometry_vertex(g.vidtoindex[voutid]);
                g.update_normals_vertex(g.vidtoindex[voutid]);
                g.update_excluder_top_vertex(g.vidtoindex[voutid]);



                double de = g.stretch_energy(g.heidtoindex[g.Nhelast - 2]);
                de += g.dimer_bend_energy(g.heidtoindex[g.Nhelast - 2]) + g.dimer_bend_energy(heindex0);
                de += g.bend_energy(heindex0) + g.bend_energy(xidindex)-e1;

                //std::cout << " crit is " << crit << endl;
                de += gbb - g.mu[etypenew];
                if(debug_sheet==1){
                    std::cout << "Binding free energy: " << gbb << ", chemical potential: " << g.mu[etypenew] << ", total: " << gbb - g.mu[etypenew] << std::endl;
                    if(de<0){
                        std::cout << "Favorable free energy! Found critical nucleus! Press any button to continue..." << std::endl;
                        std::cin.get();
                    }
                }

                /*int vinid=g.he[g.heidtoindex[g.Nhelast - 1]].vin;
                    int vindex0=g.vidtoindex[vinid];
                    
                    if (g.v[vindex0].vneigh.size() > 0)
                    {
                        for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                        {	
                            g.update_neigh_vertex(*it);
                        }
                    }
                    int voutid=g.he[g.heidtoindex[g.Nhelast - 1]].vout;
                    vindex0=g.vidtoindex[voutid];
                    if (g.v[vindex0].vneigh.size() > 0)
                    {
                        for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                        {	
                            g.update_neigh_vertex(*it);
                        }
                    }*/

                double crit = exp((-de) / g.T) / 2;
                int overlapflag = -1;

                /* check overlap */
                for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
                {
                    if (g.check_overlap_g(*it) < 0)
                        overlapflag = 1;
                }

                if (gsl_rng_uniform(rg) < crit && overlapflag == -1)
                {
                    // ToDo
                    // update next_previous boundary
                    //std::cout << "00 monomer added"<<endl;
                    g.he[g.heidtoindex[g.Nhelast - 1]].boundary_index = bi;
                    g.v[g.vidtoindex[vin0]].doubleboundary = -1;

                    // if triangle otherside close it
                    if (g.he[g.heidtoindex[g.he[xidindex].previd_boundary]].vin == g.he[g.heidtoindex[g.he[heindex0].nextid_boundary]].vout)
                    {
                        g.set_prev_next(g.he[xidindex].previd_boundary, g.he[heindex0].nextid_boundary, g.Nhelast - 1);
                        g.set_prev_next(g.he[heindex0].nextid_boundary, g.Nhelast - 1, g.he[xidindex].previd_boundary);
                        g.set_prev_next(g.Nhelast - 1, g.he[xidindex].previd_boundary, g.he[heindex0].nextid_boundary);
                        g.Nboundary--;
                    }
                    //std::cout <<"00 added with previous accepted" <<endl;
                    vecupdate.clear();
                    return 1;
                }
                /*deleting the added monomer */
                else
                {

                    int nextid0 = g.he[g.heidtoindex[g.Nhelast - 2]].nextid;
                    int previd0 = g.he[g.heidtoindex[g.Nhelast - 2]].previd;
                    int nextidboundary0 = g.he[g.heidtoindex[g.Nhelast - 1]].nextid_boundary;
                    int previdboundary0 = g.he[g.heidtoindex[g.Nhelast - 1]].previd_boundary;

                    if (nextid0 == -1 || previd0 == -1)
                    {
                        std::cout << "not accepted in add monomer!" << endl;
                        std::exit(-1);
                    }
                    //std::cout << "!!!!!!!" <<endl;
                    int xidindex = g.heidtoindex[xid];
                    if (g.delete_edge(g.he[g.Nhe - 1].id) > 0)
                    {

                        g.he[heindex0].nextid = -1;
                        g.he[xidindex].previd = -1;

                        g.update_half_edge(heid0);
                        g.update_half_edge(g.he[heindex0].opid);
                        g.update_half_edge(xid);
                        g.update_half_edge(g.he[xidindex].opid);

                        // update // update boundary_index

                        g.he[heindex0].boundary_index = bi;
                        g.he[xidindex].boundary_index = bi;
                        g.set_prev_next_boundary(heid0, nextidboundary0);
                        g.set_prev_next_boundary(previdboundary0, xid);
                        g.set_prev_next_boundary(xid, heid0);
//                        std::cout << "updating index after deleting monomer (l 1316)" << std::endl;
                        g.update_index();

                        //std::cout << " MONOMER REMOVED AFTER addition xid2" <<endl;

                        for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
                        {
                            g.update_neigh_vertex(*it);
                        }

                        //ToDo update next_previous boundary

                        //std::cout <<"00 added with previous not accepted" <<endl;
                        vecupdate.clear();
                        return -1;
                    }
                    else
                    {
                        std::cout << " could not delete_monomer HERE 555";
                        std::exit(-1);
                    }
                }
            }
            else
            {
                //std::cout << "could not add monomer" << endl;
                return -1;
            }
        }
    }
    else if (g.v[g.vidtoindex[g.he[heindex0].vin]].hein.size() < 6 && (g.v[g.vidtoindex[g.he[heindex0].vout]].hein.size()) < 6) //adding dimer
    {
        if(debug_sheet==1){
            std::cout << "Attempting to add dimer" << std::endl;
        }
        if (vin0 >= 0 || vid0 >= 0)
        {
            return (-1);
        }

        //get_dimer_etypes(*etypeheid0,*etypenew1,etypenew2,r);
        if (etypeheid0 == 0)
        {
            if (gsl_rng_uniform(rg) < .5)
            {
                etypenew1 = 1;
                etypenew2 = 2;
            } //ba
            else
            {
                etypenew1 = 0;
                etypenew2 = 0;
            }
        }

        else if (etypeheid0 == 3)
        {
            if (gsl_rng_uniform(rg) < .5)
            {
                etypenew1 = 3;
                etypenew2 = 3;
            } //ba
            else
            {
                etypenew1 = 1;
                etypenew2 = 2;
            }
        }

        else if (etypeheid0 == 1)
        {
            if (gsl_rng_uniform(rg) < .5)
            {
                etypenew1 = 2;
                etypenew2 = 0;
            } //ba
            else
            {
                etypenew1 = 2;
                etypenew2 = 3;
            }
        }

        else if (etypeheid0 == 2)
        {
            if (gsl_rng_uniform(rg) < .5)
            {
                etypenew1 = 0;
                etypenew2 = 1;
            }
            else
            {
                etypenew1 = 3;
                etypenew2 = 1;
            }
        }
        //std::cout << "adding dimer"<<endl;
        // Based on concentration
        /*    if (gsl_rng_uniform(r) < g.cdProb) {
                    if (gsl_rng_uniform(r) < 0.5) {etypenew1=0;}
                    else {etypenew1=3; }
                } 
                else{
                    if (gsl_rng_uniform(r) < 0.5) {etypenew1=1;}
                    else {etypenew1=2; } 

                } 

                if (gsl_rng_uniform(r) < g.cdProb) {
                    if (gsl_rng_uniform(r) < 0.5) {etypenew2=0;}
                    else {etypenew2=3; }
                } 
                else{
                    if (gsl_rng_uniform(r) < 0.5) {etypenew2=1;}
                    else {etypenew2=2; } 

                }*/
        //std::cout << "new types are etypenew1= " <<etypenew1 << " etypenew2 "  << etypenew2<<endl;
        //etypenew1=gsl_rng_uniform_int(r,4);
        //etypenew2=gsl_rng_uniform_int(r,4);

        bool drug1 = 0;
        bool drug2 = 0;
        //bool drug10=0;
        //bool drug20=0;
        /*if ((etypenew1==3) ||(etypenew1==0) ) {
                if (gsl_rng_uniform(r) < g.drugProb) {drug1=1;}
                if (gsl_rng_uniform(r) < g.drugProb) {drug10=1;}
            }
            if ((etypenew2==3) ||(etypenew2==0) ) {
                if (gsl_rng_uniform(r) < g.drugProb) {drug2=1;}
                if (gsl_rng_uniform(r) < g.drugProb) {drug20=1;}
            }*/

        //std::cout << "types " << etype1 <<" " << etype2<<endl;
        //double gbb=g.find_gbb(etypeheid0,etypenew1,etypenew2);
        double gbb = g.find_dg(etypeheid0, etypenew1, drug1);
        gbb += g.find_dg(etypenew1, etypenew2, drug2);
        gbb += g.find_dg(etypenew2, etypeheid0, g.he[heindex0].din);
        if(debug_sheet==1){
            std::cout << "Binding interface type: " << etypeheid0 << "-" << etypenew1 << ", energy 1: " << g.find_dg(etypeheid0, etypenew1, drug1) << std::endl;
            std::cout << "Binding interface type: " << etypenew1 << "-" << etypenew2 << ", energy 2: " << g.find_dg(etypenew1, etypenew2, drug2) << std::endl;
            std::cout << "Binding interface type: " << etypenew2 << "-" << etypeheid0 << ", energy 3: " << g.find_dg(etypenew2, etypeheid0, drug2) << std::endl;
        }
        //if
        //std::cout << "gbb is  " << gbb <<endl;
        double *dis_vector;
        dis_vector= new double[3];
        double dis_new = g.add_dimer(heid0, rg, etypenew1, etypenew2,dis_vector);


        if (dis_new < 0)
        {
            //std::cout << "could not add dimer " << endl;
            return (-1);
        }
        ///else {
        //std::cout << " dimer added , now check if it stays"  <<endl;
        // }


        g.update_half_edge(g.Nhelast - 2);
        //g.he[g.heidtoindex[g.Nhelast-2]].din=drug2; //double check in add_dimer
        g.update_half_edge(g.Nhelast - 1);
        //g.he[g.get_heindex(g.Nhelast - 1)].din=drug20;
        g.update_half_edge(g.Nhelast - 4);
        //g.he[g.get_heindex(g.Nhelast - 4)].din=drug1;
        g.update_half_edge(g.Nhelast - 3);
        //g.he[g.get_heindex(g.Nhelast - 3)].din=drug10;
        g.update_boundary();

        ////std::cout<< "op of edge" << g.Nhelast-2 << " is " << g.he[g.heidtoindex[g.Nhelast-2]].opid<<endl;
        int index2 = g.heidtoindex[g.Nhelast - 2];
        int index4 = g.heidtoindex[g.Nhelast - 4];
        double de = 0;
        if(debug_sheet==1){
            std::cout << "Total elastic energy: " << g.stretch_energy(index2) + g.dimer_bend_energy(index2) + g.stretch_energy(index4) + g.dimer_bend_energy(index4) + g.bend_energy(heindex0) + g.dimer_bend_energy(heindex0)<< std::endl;
        }
        else{
            de = g.stretch_energy(index2) + g.dimer_bend_energy(index2);
            de += g.stretch_energy(index4) + g.dimer_bend_energy(index4);
            de += g.bend_energy(heindex0) + g.dimer_bend_energy(heindex0);
        }

        /* test !!!!
            if (abs(de-(g.compute_energy()-e11)>0.0000000001))  
            {//std::cout << "in add dimer des don't match"<<endl;
            std::cout<<"de-(g.compute_energy()-e11)" << de-(g.compute_energy()-e11) << endl; 
            std::cout << "e11 "<< e11 <<endl; std::cout << "de "<< de <<endl; 
            std::cout << "g.compute_energy() "<< g.compute_energy() <<endl;  
            std::cout << "g.compute_energy() -e11 "<< g.compute_energy()-e11 <<endl; std::exit(-1);} */

        //**************************
        // gaussian correction
        //**************************
        
        //std::cout << "dis_vector" << dis_vector[0] << " " << dis_vector[1] << " " << dis_vector[2] <<endl;
        
        //std::cout << "pdfs " << gsl_ran_gaussian_pdf(dis_vector[0], g.l_thermal_sigma) << " "<<gsl_ran_gaussian_pdf(dis_vector[1], g.l_thermal_sigma) << " "<< gsl_ran_gaussian_pdf(dis_vector[2], g.l_thermal_kappa) <<endl;
        //std::cout<<"pdfs calac "<<endl;
        
        
       
        //double vp = 1/(gsl_ran_gaussian_pdf(dis_new, g.gaussian_sigma));

        double vp = pow((sqrt(2*M_PI)*g.gaussian_sigma),3)/( exp(-((dis_new*dis_new)/(2*g.gaussian_sigma*g.gaussian_sigma))) );

        de += gbb - (g.mu[etypenew1] + g.mu[etypenew2]);
        if(debug_sheet==1){
            std::cout << "Total Binding free energy: " << gbb << ", total chemical potential: " << g.mu[etypenew1] + g.mu[etypenew2] << ", total: " << de << std::endl;
            if(de<0){
                std::cout << "Favorable free energy! Found critical nucleus! Press any button to continue..." << std::endl;
                std::cin.get();
            }
        }
        delete[] dis_vector;    

        double crit = 0;
        if (debug_sheet==1){
            crit = exp(-de / g.T);
            std::cout << "vp: " << vp << std::endl;
            std::cout << "log(2vp): " << log(2*vp) << std::endl;
        }
        else{
            crit = 2 * vp * exp(-de / g.T);
        }
        //std::cout << " crit is " << crit << endl;
        int overlapflag = -1;
        //g.update_index();
        g.update_neigh_vertex(g.Nvlast - 1);
        g.update_geometry_vertex(g.vidtoindex[g.Nvlast - 1]);
        g.update_normals_vertex(g.vidtoindex[g.Nvlast - 1]);
        g.update_excluder_top_vertex(g.vidtoindex[g.Nvlast - 1]);
        //if (g.Nhe>190) std::cout<<" in add dimer updating neigh of added vetex"<<endl;
        /* save vecupdate */
        int vindex0 = g.vidtoindex[g.Nvlast - 1];
        vector<int> vecupdate;
        if (g.v[vindex0].vneigh.size() > 0)
        {
            for (vector<int>::iterator it = g.v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
            {
                //if (g.Nhe>190) std::cout << "neigh of new vertex are" << *it <<endl;
                vecupdate.push_back(*it);
                g.update_neigh_vertex(*it);
            }
        }

        if (g.check_overlap_g(g.Nvlast - 1) < 0)
        {
            overlapflag = 1;
        }

        if (gsl_rng_uniform(rg) < crit && overlapflag == -1) //dimer added
        {
            // ToDo
            // update boundary index
            // update next_previous boundary
            g.he[g.heidtoindex[g.Nhelast - 3]].boundary_index = bi;
            g.he[g.heidtoindex[g.Nhelast - 1]].boundary_index = bi;
            //std::cout << "00 added dimer accepeted" <<endl;
            vecupdate.clear();

            //##############################
            //############# dist add accepted
            /*ofstream myfile;
            myfile.open ("adddist.txt", ios::out | ios::app);
            myfile <<dis_new<<endl;
            myfile.close();*/

            return 2;
        }
        /* removing the added dimer */
        else
        {
            //##############################
            //############# dist add accepted
            /*ofstream myfile;
            myfile.open ("reject_adddist.txt", ios::out | ios::app);
            myfile <<dis_new<<endl;
            myfile.close();*/

            int nextidboundary0 = g.he[g.heidtoindex[g.he[g.heidtoindex[g.he[g.heidtoindex[heid0]].nextid]].opid]].nextid_boundary;
            int previdboundary0 = g.he[g.heidtoindex[g.he[g.heidtoindex[g.he[g.heidtoindex[heid0]].previd]].opid]].previd_boundary;

            if (g.delete_edge(g.Nhelast - 1) < 0)
            {
                std::cout << " could not delete_added_dimer HERE 777";
                std::exit(-1);
            }

//            std::cout << "updating index after deleting dimer (l 1606)" << std::endl;
            g.update_index();
            if (g.delete_edge(g.Nhelast - 3) < 0)
            {
                std::cout << " could not delete_dimer HERE 888";
                std::exit(-1);
            }

//            std::cout << "updating index after deleting dimer pt 2 (l 1614)" << std::endl;
            g.update_index();

            if (g.delete_vertex(g.Nvlast - 1) > 0)
            {
                g.set_prev_next(heid0, -1, -1);

//                std::cout << "updating index after deleting vertex (l 1621)" << std::endl;
                g.update_index();
                /* update neigh of vecupdate*/
                for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
                {
                    g.update_neigh_vertex(*it);
                }
                //ToDo
                //update next_ previous boundary
                //update boundary index
                g.he[g.heidtoindex[heid0]].boundary_index = bi;
                g.set_prev_next_boundary(heid0, nextidboundary0);
                g.set_prev_next_boundary(previdboundary0, heid0);

                vecupdate.clear();

                //std::cout << "00 added dimer edges and vertices are deleted" <<endl;
                //heindex0=g.heidtoindex[heid0];
                // MAYBE NOT NEEDED XXX

                //g.update_boundary();
                return 0;
            }
            else
            {
                std::cout << "HERE 1111";
                std::exit(-1);
            }
            //vecupdate.clear();
            //return -1;
        }
    }
    return -1;
}

int MC::attempt_add_monomer_dimer_drug(System &g, int heid0) //!!! Should update with pre_oipen wedge
{
    //std::cout << "in attempt_add_monomer-dimer" <<endl;
    g.update_normals();

    if (g.is_boundary(heid0) < 0)
    {
        std::cout << " cannot add not on the !" << endl;
        std::exit(-1);
    }
    //double gbb=0;

    int etypenew1 = -1;
    int etypenew2 = -1;

    int heindex0 = g.heidtoindex[heid0];
    int etypeheid0 = g.he[heindex0].type;

    int bi = g.he[heindex0].boundary_index;
    //temp test
    if (bi == -1)
    {
        std::cout << "wrong boundary index in add monomer dimer" << endl;
        std::exit(-1);
    }

    //should not add if wedge fusion is possible
    g.update_fusion_pairs_he();

    if (g.he[heindex0].next_wedge_fusion_heid != -1 || g.he[heindex0].prev_wedge_fusion_heid != -1)
    {
        return -1;
    }
    //monomer with next

    //if (g.he[heindex0].nextid!=-1 && g.he[heindex0].previd!=-1) {std::cout <<" wrong not on boundary" <<endl; std::exit(-1); }
    //int xid = -1;

    int vid0 = g.is_bond_out_boundary(heid0);
    int vin0 = g.is_bond_in_boundary(heid0);
    //adding monomer  with previous
    
    if (g.v[g.vidtoindex[g.he[heindex0].vin]].hein.size() < 6 && (g.v[g.vidtoindex[g.he[heindex0].vout]].hein.size()) < 6) //adding dimer
    {
        if (vin0 >= 0 || vid0 >= 0)
        {
            return (-1);
        }

        //get_dimer_etypes(*etypeheid0,*etypenew1,etypenew2,r);
        

        if (etypeheid0 == 3 )
            {
                etypenew1 = 3;
                etypenew2 = 3;
            } 
        
        if (etypeheid0 == 0)
            {
                etypenew1 = 0;
                etypenew2 = 0;
            }
        else{
            return(-1);
        }
        

        

        
        //std::cout << "adding dimer"<<endl;
        // Based on concentration
        /*    if (gsl_rng_uniform(r) < g.cdProb) {
                    if (gsl_rng_uniform(r) < 0.5) {etypenew1=0;}
                    else {etypenew1=3; }
                } 
                else{
                    if (gsl_rng_uniform(r) < 0.5) {etypenew1=1;}
                    else {etypenew1=2; } 

                } 

                if (gsl_rng_uniform(r) < g.cdProb) {
                    if (gsl_rng_uniform(r) < 0.5) {etypenew2=0;}
                    else {etypenew2=3; }
                } 
                else{
                    if (gsl_rng_uniform(r) < 0.5) {etypenew2=1;}
                    else {etypenew2=2; } 

                }*/
        //std::cout << "new types are etypenew1= " <<etypenew1 << " etypenew2 "  << etypenew2<<endl;
        //etypenew1=gsl_rng_uniform_int(r,4);
        //etypenew2=gsl_rng_uniform_int(r,4);

        bool drug1 = 0;
        bool drug2 = 0;
        //bool drug10=0;
        //bool drug20=0; 
        /*if ((etypenew1==3) ||(etypenew1==0) ) {
                if (gsl_rng_uniform(r) < g.drugProb) {drug1=1;}
                if (gsl_rng_uniform(r) < g.drugProb) {drug10=1;}
            }
            if ((etypenew2==3) ||(etypenew2==0) ) {
                if (gsl_rng_uniform(r) < g.drugProb) {drug2=1;}
                if (gsl_rng_uniform(r) < g.drugProb) {drug20=1;}
            }*/

        //std::cout << "types " << etype1 <<" " << etype2<<endl;
        //double gbb=g.find_gbb(etypeheid0,etypenew1,etypenew2);
        drug2 = 1;
        double gbb = g.find_dg(etypeheid0, etypenew1, drug1);
        gbb += g.find_dg(etypenew1, etypenew2, drug2);
        gbb += g.find_dg(etypenew2, etypeheid0, g.he[heindex0].din);
        //if
        //std::cout << "gbb is  " << gbb <<endl;
        double *dis_vector;
        dis_vector= new double[3];
        double dis_new = g.add_dimer(heid0, rg, etypenew1, etypenew2,dis_vector);


        if (dis_new < 0)
        {
            //std::cout << "could not add dimer " << endl;
            return (-1);
        }
        ///else {
        //std::cout << " dimer added , now check if it stays"  <<endl;
        // }


        g.update_half_edge(g.Nhelast - 2);
        g.he[g.heidtoindex[g.Nhelast-2]].din=drug2; //double check in add_dimer
        g.update_half_edge(g.Nhelast - 1);
        //g.he[g.get_heindex(g.Nhelast - 1)].din=drug20;
        g.update_half_edge(g.Nhelast - 4);
        g.he[g.heidtoindex[g.Nhelast - 4]].din=drug1;
        g.update_half_edge(g.Nhelast - 3);
        //g.he[g.get_heindex(g.Nhelast - 3)].din=drug10;
        g.update_boundary();

        ////std::cout<< "op of edge" << g.Nhelast-2 << " is " << g.he[g.heidtoindex[g.Nhelast-2]].opid<<endl;
        int index2 = g.heidtoindex[g.Nhelast - 2];
        int index4 = g.heidtoindex[g.Nhelast - 4];
        double de = 0;
        if(debug_sheet!=1){
            de = g.stretch_energy(index2) + g.dimer_bend_energy(index2);
            de += g.stretch_energy(index4) + g.dimer_bend_energy(index4);
            de += g.bend_energy(heindex0) + g.dimer_bend_energy(heindex0);
        }

        /* test !!!!
            if (abs(de-(g.compute_energy()-e11)>0.0000000001))  
            {//std::cout << "in add dimer des don't match"<<endl;
            std::cout<<"de-(g.compute_energy()-e11)" << de-(g.compute_energy()-e11) << endl; 
            std::cout << "e11 "<< e11 <<endl; std::cout << "de "<< de <<endl; 
            std::cout << "g.compute_energy() "<< g.compute_energy() <<endl;  
            std::cout << "g.compute_energy() -e11 "<< g.compute_energy()-e11 <<endl; std::exit(-1);} */

        //**************************
        // gaussian correction
        //**************************
        
        //std::cout << "dis_vector" << dis_vector[0] << " " << dis_vector[1] << " " << dis_vector[2] <<endl;
        
        //std::cout << "pdfs " << gsl_ran_gaussian_pdf(dis_vector[0], g.l_thermal_sigma) << " "<<gsl_ran_gaussian_pdf(dis_vector[1], g.l_thermal_sigma) << " "<< gsl_ran_gaussian_pdf(dis_vector[2], g.l_thermal_kappa) <<endl;
        //std::cout<<"pdfs calac "<<endl;
        
        
       
        //double vp = 1/(gsl_ran_gaussian_pdf(dis_new, g.gaussian_sigma));

        double vp = pow((sqrt(2*M_PI)*g.gaussian_sigma),3)/( exp(-((dis_new*dis_new)/(2*g.gaussian_sigma*g.gaussian_sigma))) );

        de += gbb - (g.mu[etypenew1] + g.mu[etypenew2]);
        delete[] dis_vector;    
        double crit = 2 * vp * exp(-de / g.T);
        //std::cout << " crit is " << crit << endl;
        int overlapflag = -1;
        //g.update_index();
        g.update_neigh_vertex(g.Nvlast - 1);
        g.update_geometry_vertex(g.vidtoindex[g.Nvlast - 1]);
        g.update_normals_vertex(g.vidtoindex[g.Nvlast - 1]);
        g.update_excluder_top_vertex(g.vidtoindex[g.Nvlast - 1]);
        //if (g.Nhe>190) std::cout<<" in add dimer updating neigh of added vetex"<<endl;
        /* save vecupdate */
        int vindex0 = g.vidtoindex[g.Nvlast - 1];
        vector<int> vecupdate;
        if (g.v[vindex0].vneigh.size() > 0)
        {
            for (vector<int>::iterator it = g.v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
            {
                //if (g.Nhe>190) std::cout << "neigh of new vertex are" << *it <<endl;
                vecupdate.push_back(*it);
                g.update_neigh_vertex(*it);
            }
        }

        if (g.check_overlap_g(g.Nvlast - 1) < 0)
        {
            overlapflag = 1;
        }

        if (gsl_rng_uniform(rg) < crit && overlapflag == -1) //dimer added
        {
            //std::cout << "Added drug-bound dimer" << endl;
            // ToDo
            // update boundary index
            // update next_previous boundary
            g.he[g.heidtoindex[g.Nhelast - 3]].boundary_index = bi;
            g.he[g.heidtoindex[g.Nhelast - 1]].boundary_index = bi;
            g.Nd++;
            //if(drug2==1){g.Nd++;}
            //std::cout << "00 added dimer accepeted" <<endl;
            vecupdate.clear();

            //##############################
            //############# dist add accepted
            /*ofstream myfile;
            myfile.open ("adddist.txt", ios::out | ios::app);
            myfile <<dis_new<<endl;
            myfile.close();*/

            return 2;
        }
        /* removing the added dimer */
        else
        {
            //##############################
            //############# dist add accepted
            /*ofstream myfile;
            myfile.open ("reject_adddist.txt", ios::out | ios::app);
            myfile <<dis_new<<endl;
            myfile.close();*/

            int nextidboundary0 = g.he[g.heidtoindex[g.he[g.heidtoindex[g.he[g.heidtoindex[heid0]].nextid]].opid]].nextid_boundary;
            int previdboundary0 = g.he[g.heidtoindex[g.he[g.heidtoindex[g.he[g.heidtoindex[heid0]].previd]].opid]].previd_boundary;
            g.he[g.heidtoindex[g.Nhelast-2]].din=0;
            //g.Nd--;
            if (g.delete_edge(g.Nhelast - 1) < 0)
            {
                std::cout << " could not delete_added_dimer HERE 777";
                std::exit(-1);
            }

//            std::cout << "updating index after deleting dimer with drug (l 1902)" << std::endl;
            g.update_index();
            if (g.delete_edge(g.Nhelast - 3) < 0)
            {
                std::cout << " could not delete_dimer HERE 888";
                std::exit(-1);
            }

//            std::cout << "updating index after deleting dimer with drug, pt2 (l 1910)" << std::endl;
            g.update_index();

            if (g.delete_vertex(g.Nvlast - 1) > 0)
            {
                g.set_prev_next(heid0, -1, -1);

//                std::cout << "updating index after deleting vertex in dimer with drug (l 1917)" << std::endl;
                g.update_index();
                /* update neigh of vecupdate*/
                for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
                {
                    g.update_neigh_vertex(*it);
                }
                //ToDo
                //update next_ previous boundary
                //update boundary index
                g.he[g.heidtoindex[heid0]].boundary_index = bi;
                g.set_prev_next_boundary(heid0, nextidboundary0);
                g.set_prev_next_boundary(previdboundary0, heid0);

                vecupdate.clear();

                //std::cout << "00 added dimer edges and vertices are deleted" <<endl;
                //heindex0=g.heidtoindex[heid0];
                // MAYBE NOT NEEDED XXX

                //g.update_boundary();
                return 0;
            }
            else
            {
                std::cout << "HERE 1111";
                std::exit(-1);
            }
            //vecupdate.clear();
            //return -1;
        }
    }
    return -1;
}

int MC::attempt_remove_monomer_dimer(System &g, int heid0) /* 102220 THIS NEEDS UPDATE -> MOVE GEOMETRY TO geometry! */
{
    //std::cout << "in attempt_remove_monomer_dimer" << endl;
    //std::cout << "Nd is " <<g.Nd <<endl;

    if (!(g.is_boundary(heid0) > 0))
    {
        std::cout << "not on boundary cannot remove" << endl;
        std::exit(-1);
    }

    g.update_half_edge(heid0);
    int heindex0 = g.heidtoindex[heid0]; // index of this edge
    //int heid0type=g.he[heindex0].type;
    if ((g.he[heindex0].nextid != -1) || (g.he[heindex0].previd != -1))
    {
        //std::cout << " has next or previous cannot remove" << endl;
        return -1;
    }

    int bi = g.he[heindex0].boundary_index;
    if (bi == -1)
    {
        std::cout << "Error in boundary index in remove_monomer_dimer" << endl;
        std::exit(-1);
    }

    g.update_half_edge(g.he[heindex0].opid);

    int heopindex0 = g.heidtoindex[g.he[heindex0].opid]; // indexd of opposite edge
    int optype = g.he[heopindex0].type;
    int nextopid0 = g.he[heopindex0].nextid; // id of next of opposite edge
    int prevopid0 = g.he[heopindex0].previd;
    //std::cout << " nextopid0 " << nextopid0 << " prevopid0 " <<prevopid0 <<endl;
    if ((nextopid0 == -1) || (prevopid0 == -1))
    {
        std::cout << "wrong geometry" << endl;
        std::exit(-1);
    }
    int nextopindex0 = g.heidtoindex[nextopid0]; // id of prev of opposite edge
    int prevopindex0 = g.heidtoindex[prevopid0];
    int opnexttype = g.he[nextopindex0].type;
    int opprevtype = g.he[prevopindex0].type;

    if ((g.he[heindex0].din == 1) || g.he[heopindex0].din == 1)
    {
        //std::cout << " din=1 "<<endl;
        return (-1);
    }                                    // WITH DRUG NO REMOVAL
    //LBF 12/13/25: Is there a problem here?
    g.update_boundary();
    int heid_prev_boundary = g.he[nextopindex0].opid; // now back to this side // ToDo this should be previd_boundary
    int heid_next_boundary = g.he[prevopindex0].opid; // after vertex // ToDo this should be nextid_boundary
    int heid_prev_boundary2 = g.he[heindex0].previd_boundary;
    int heid_next_boundary2 = g.he[heindex0].nextid_boundary;
    //std::cout << "compare: " << heid_prev_boundary << " vs " << heid_prev_boundary2 << "; " << heid_next_boundary << " vs " << heid_next_boundary2 << std::endl;
    //the following lines cause an error in new_vertex_edge
    //heid_prev_boundary = g.he[heindex0].previd_boundary;
    //heid_next_boundary = g.he[heindex0].nextid_boundary;

    double gbb = 0;

    // remove monomer
    // first try remove monomer, if the edge is not a wedge, remove dimer
    if (g.is_boundary(heid_prev_boundary) < 0 && g.is_boundary(heid_next_boundary) < 0)
    {
     if ((g.is_vboundary(g.he[nextopindex0].vout) < 0 || g.Nboundary != 1) && \
       (g.v[g.vidtoindex[g.he[heindex0].vout]].doubleboundary==-1 && g.v[g.vidtoindex[g.he[heindex0].vin]].doubleboundary==-1 ) )// if Nboundary>0 allow for double boundary;
    {                                                                                                                               //delete monomer
        //std::cout << "in deleting monomer" <<endl;
        //if (g.is_vboundary(g.he[g.heidtoindex[nextopid0]].vout) > 0) { std::cout << "v on boundary wrong geometry!"<<endl; std::exit(-1);}
        if ((g.he[nextopindex0].din == 1) || g.he[prevopindex0].din == 1)
        {
            return (-1);
        }

        double de = -(g.stretch_energy(heindex0));
        int nextboundary0 = g.he[heindex0].nextid_boundary;
        int prevboundary0 = g.he[heindex0].previd_boundary;
        //if (g.he[heindex0].previd!=-1) { de-=g.dimer_bend_energy(g.get_heindex(g.he[heindex0].previd)); }
        //if (g.he[heindex0].nextid!=-1) {de-=g.dimer_bend_energy(heindex0); }
        de -= (g.dimer_bend_energy(heopindex0) + g.dimer_bend_energy(prevopindex0));
        //de -= g.bend_energy(nextopindex0) + g.bend_energy(prevopindex0); //g.monomer_energy(heid0);

        gbb += g.find_dg(opprevtype, optype, g.he[heopindex0].din);
        gbb += g.find_dg(optype, opnexttype, g.he[nextopindex0].din);

        de -= (gbb - g.mu[g.he[heindex0].type]);
        double crit = 2 * exp(-de / g.T); ///(2.0*g.z*g.K*g.K*g.K);

        if (g.Test_assembly == 1)
        {
            std::cout << "crit is " << crit << endl;
            crit = 1;
        }
        if (gsl_rng_uniform(rg) < crit)
        {
            //g.Nd-=g.he[heopindex0].din;
            //g.Nd-=g.he[heindex0].din;
            //std::cout << "REMOVING MONOMER two drugs removed " << g.he[heopindex0].din + g.he[heindex0].din <<endl;

            //std::cout << "de is" <<de <<endl;
            //std::cout << "crit is  " << crit <<endl;
            //if ( g.he[heindex0].previd!=-1) {g.he[g.get_heindex(g.he[heindex0].previd)].nextid=-1;}
            //if (g.he[heindex0].nextid!=-1) {g.he[g.get_heindex(g.he[heindex0].nextid)].previd=-1; }
            //std::cout << "g.he[heopindex0].previ " << g.he[heopindex0].previd <<endl;
            //std::cout << "g.he[prevopindex0].nextid " << g.he[prevopindex0].nextid <<endl;
            //std::cout << "g.he[heopindex0].nextid " << g.he[heopindex0].nextid <<endl;
            //std::cout << "g.he[nextopindex0].previd" << g.he[nextopindex0].previd <<endl;
            if (g.he[heopindex0].previd != -1)
            {
                g.he[prevopindex0].nextid = -1;
            }
            if (g.he[heopindex0].nextid != -1)
            {
                g.he[nextopindex0].previd = -1;
            }

            //ToDo -> Done
            //update new edges boundary index
            //update new edges nex_boundary
            g.he[nextopindex0].boundary_index = bi;
            g.he[prevopindex0].boundary_index = bi;

            g.set_prev_next_boundary(nextopid0, prevopid0);
            g.set_prev_next_boundary(prevboundary0, nextopid0);
            g.set_prev_next_boundary(prevopid0, nextboundary0);

            if (g.is_vboundary(g.he[nextopindex0].vout) > 0)
                g.v[g.vidtoindex[g.he[nextopindex0].vout]].doubleboundary = 1;

            int vidin = g.he[heindex0].vin;
            int vidout = g.he[heindex0].vout;

            int x = g.delete_edge(heid0);
            //g.update_edge();_

            if (x < 0)
            {
                std::cout << "could not delete " << endl;
                std::exit(-1);
            }
            //

            //g.set_prev_next(nextopid0, -1, prevopid0);
            //g.set_prev_next(prevopid0, nextopid0, -1);
            //std::cout << "00 monomer removed" <<endl;
            //std::cout << " 004 g.Nd is " <<g.Nd<<endl;
//            std::cout << "updating index after removing monomer (l 2092)" << std::endl;
            g.update_index();
            g.update_neigh_vertex(vidin);
            g.update_neigh_vertex(vidout);

            //std::cout << " 005 g.Nd is " <<g.Nd<<endl;

            return 1;
        }
        else
        {
            //std::cout << " remove monomer not accepted  " << endl;
            //g.add_monomer(nextopid0, prevopid0,optype);

            return -1;
        }
    }
    }
            else // edge is not wedge so remove dimer //with previous or next
    {

        int heindex_prev_boundary = g.heidtoindex[heid_prev_boundary];
        int heindex_next_boundary = g.heidtoindex[heid_next_boundary];
        //double *vco=new double[3];
        //remove dimer
        if (g.is_boundary(heid_prev_boundary) > 0 && g.he[heindex_prev_boundary].previd == -1 && g.v[g.vidtoindex[g.he[heindex0].vin]].doubleboundary==-1) //delete dimer this and next (inside)(nextopindex) / this and prev (on boundary)  there should be nno bonds between this and previous
        {
            //std::cout << "00 remove dimer with next (previd_boundary) "<<endl;  //remove this and previd_boundary=heid_prev_boundary
            if (g.he[heindex_prev_boundary].din == 1 || g.he[nextopindex0].din == 1)
            {
                return (-1);
            } //NOREMOVAL DRUG

            // int yid=g.he[heindex0].nextid;
            int vi = g.he[heopindex0].vout;

            //Reject moves that attempt to remove dimers with more than one vertex
            int still_in_vout = check_vid_in_vout(g, vi);
            if(still_in_vout>2){
                std::cout << "Rejecting attempt to remove dimer where vertex is shared by others" << std::endl;
                return (-1);
            }

            //gbb+= g.find_gbb(optype,opnexttype,opprevtype); //whole triangle
            gbb += g.find_dg(optype, opnexttype, g.he[g.heidtoindex[nextopid0]].din);
            gbb += g.find_dg(opnexttype, opprevtype, g.he[prevopindex0].din);
            gbb += g.find_dg(opprevtype, optype, g.he[heopindex0].din);

            double de = 0;
            if(debug_sheet!=1){
                de = -(g.stretch_energy(heindex0) + g.stretch_energy(g.heidtoindex[heid_prev_boundary]));
                de -= g.bend_energy(g.heidtoindex[heid_next_boundary]);
                de -= (g.dimer_bend_energy(heopindex0));
                de -= (g.dimer_bend_energy(nextopindex0) + g.dimer_bend_energy(prevopindex0));
            }
            //if (yid!=-1) {  de-=g.dimer_bend_energy(heindex0);}

            //if (xid!=-1 && g.he[xindex].nextid==heid_prev_boundary) {
            //    gbb+=g.find_dg( g.he[xindex].type,g.he[g.heidtoindex[heid_prev_boundary]].type);
            //     de-=g.dimer_bend_energy(xindex);
            //xidconnect=1;
            //}

            
            /* test !!!!
                if (abs(de-(g.compute_energy()-e11)>0.0000000001)) 
                {//std::cout << "in remove dimer 111 des don't match" << de-(g.compute_energy()-e11) << endl; 
                std::cout << "e11 "<< e11 <<endl; std::cout << "de "<< de <<endl; std::cout << "g.compute_energy() "<< g.compute_energy() <<endl;  
                std::cout << "g.compute_energy() -e11 "<< g.compute_energy()-e11 <<endl; std::exit(-1);} */

            de -= (gbb - (g.mu[g.he[heindex0].type] + g.mu[g.he[heindex_prev_boundary].type]));
            //
            //**************************
            // gaussian correction
            //**************************
            double *tempv1 = new double[3];
            //double *dis_vector=new double[3];
            int heidtemp=g.he[nextopindex0].nextid;
            g.new_vertex_edge(g.heidtoindex[heidtemp], tempv1,g.he[heopindex0].type); // heopindex0 = g.heidtoindex[g.he[g.heidtoindex[heidtemp]].nextid]
            //std::cout<<"tmpv1 "<< tempv1[0] <<" "<<tempv1[1] <<" "<<tempv1[2] <<" "<<endl;
            double dis_new=veclen(tempv1,g.v[g.vidtoindex[vi]].co);
            //std::cout<<dis_new << " " << isnan(dis_new)<<endl;
            //if (isnan(dis_new)) std::exit(-1);
            //subvec(tempv1,g.v[g.vidtoindex[vi]].co,dis_vector);
            //double dis_new=g.find_project_dist_axes(g.heidtoindex[heidtemp], tempv1, g.v[g.vidtoindex[vi]].co , dis_vector );
            //############# dist remove
            /*ofstream myfile;
            myfile.open ("removedist.txt", ios::out | ios::app);
            myfile <<dis_new<<endl;
            myfile.close();*/
            //std::cout<<x<<endl;
            delete[] tempv1;
            //double vp = 1/(gsl_ran_gaussian_pdf(dis_new, g.gaussian_sigma));
            double crit = 0;
            if (debug_sheet==1){
                crit = exp(-de / g.T);
            }
            else{
                double vp = pow((sqrt(2*M_PI)*g.gaussian_sigma),3)/( exp(-((dis_new*dis_new)/(2*g.gaussian_sigma*g.gaussian_sigma))) );
                crit = exp(-de / g.T) / (2 * vp); 
            }
            //std::cout << " crit is " << crit << endl;
            if (gsl_rng_uniform(rg) < crit) //delete dimer this and next (inside)(nextopindex) / this and prev (on boundary)
            {
                int nextidboundary0 = g.he[g.heidtoindex[heid0]].nextid_boundary;
                int previdboundary0 = g.he[g.heidtoindex[heid_prev_boundary]].previd_boundary;

                vector<int> vecupdate;
                if (g.v[g.vidtoindex[vi]].vneigh.size() > 0)
                {
                    for (vector<int>::iterator it = g.v[g.vidtoindex[vi]].vneigh.begin(); it != g.v[g.vidtoindex[vi]].vneigh.end(); ++it)
                    {
                        vecupdate.push_back(*it);
                    }
                }
                //std::cout << " 006 g.Nd is " <<g.Nd<<endl;

                //WITHDRUG No removal
                //g.Nd-=g.he[heindex0].din;
                //g.Nd-=g.he[heopindex0].din;
                //g.Nd-=g.he[heindex_prev_boundary].din;
                //g.Nd-=g.he[nextopindex0].din;

                //std::cout << "REMOVING DIMER four drugs removed " << g.he[heopindex0].din + g.he[heindex0].din +  g.he[heindex_prev_boundary].din + g.he[nextopindex0].din<<endl;
                //std::cout << "removing heid0 and heid_prev_boundary " << heid0 <<" " << heid_prev_boundary <<endl;
                int success = g.delete_edge(heid0);
                if (success < 0)
                {
                    std::cout << "!!!" << endl;
                    std::exit(-1);
                }
//                std::cout << "updating index after removing dimer (l 2216)" << std::endl;
                g.update_index();
                success = g.delete_edge(heid_prev_boundary); // next of op
                if (success < 0)
                {
                    std::cout << "!!!" << endl;
                    std::exit(-1);
                }
//                std::cout << "updating index after removing dimer, pt 2 (l 2224)" << std::endl;
                g.update_index();
                int x = g.delete_vertex(vi);
                if (x < 0)
                {
                    std::cout << "Error deleting vertex" << std::endl;
                    std::exit(-1);
                }

                g.set_prev_next(prevopid0, -1, -1); //prev of op
                g.he[g.heidtoindex[prevopid0]].boundary_index = bi;
                g.set_prev_next_boundary(prevopid0, nextidboundary0);
                g.set_prev_next_boundary(previdboundary0, prevopid0);
                //update new edges nex_boundary

//                std::cout << "updating index after removing vertex in remove dimer (l 2238)" << std::endl;
                g.update_index();
                for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
                {
                    g.update_neigh_vertex(*it);
                }

                //ToDo
                //update new edges boundary index

                //std::cout << "00 remove dimer with next accepted "<<endl;
                vecupdate.clear();
                return 2;
            }
            // no need for anything if removal cannot happen
            else
            {
                //std::cout << "00 remove dimer with previous not accepted "<<endl;
                return 0;
            }
        }
        else if (g.is_boundary(heid_next_boundary) > 0 && g.he[heindex_next_boundary].nextid == -1) //delete dimer this and previous (inside)(prevopindex) / this and next (on boundary) there should be no bonds between next_boudary and its next
        {                                                                 //delete dimer this and prev

            //std::cout << "00 remove dimer with previous (heid_next_boundary) "<<endl; //with opid_boundary=heid_next_boundary
            //if (g.he[g.heidtoindex[heid_next_boundary]].din==1) return -1;
            //if (g.he[heindex_next_boundary].din==1 || g.he[g.heidtoindex[prevopid0]].din==1 ) { return(-1);} //WITHDRUG no removal
            if (g.he[heindex_next_boundary].din == 1 || g.he[prevopindex0].din == 1)
            {
                return (-1);
            } //NOREMOVAL DRUG
            
            int vi = g.he[heopindex0].vin; //this should be the vout of heid0
            //std::cout << "Check: should be equal: " << vi << " " << g.he[heindex0].vout << std::endl;
            //Reject moves that attempt to remove dimers with more than one vertex
            int still_in_vout = check_vid_in_vout(g, vi);
            if(still_in_vout>2){
                std::cout << "Rejecting attempt to remove dimer where vertex is shared by others" << std::endl;
                return (-1);
            }

            //gbb+= g.find_gbb(optype,opnexttype,opprevtype); //whole triangle
            gbb += g.find_dg(optype, opnexttype, g.he[g.heidtoindex[nextopid0]].din);
            gbb += g.find_dg(opnexttype, opprevtype, g.he[g.heidtoindex[prevopid0]].din);
            gbb += g.find_dg(opprevtype, optype, g.he[heopindex0].din);
            double de = -(g.stretch_energy(heindex0) + g.stretch_energy(g.heidtoindex[heid_next_boundary]));
            de -= g.bend_energy(g.heidtoindex[heid_prev_boundary]);

            de -= (g.dimer_bend_energy(heopindex0));
            de -= (g.dimer_bend_energy(g.heidtoindex[nextopid0]) + g.dimer_bend_energy(g.heidtoindex[prevopid0]));

            de -= (gbb - (g.mu[g.he[heindex0].type] + g.mu[g.he[heindex_next_boundary].type]));

            //**************************
            // gaussian correction
            //**************************
            double *tempv1 = new double[3];
            //double *dis_vector=new double[3];
            int heidtemp=g.he[heopindex0].nextid;
            g.new_vertex_edge(g.heidtoindex[heidtemp], tempv1,g.he[prevopindex0].type); // prevopindex = ? g.heidtoindex[g.he[g.heidtoindex[heidtemp]].nextid]
            double dis_new=veclen(tempv1,g.v[g.vidtoindex[vi]].co);
            //std::cout<<"tmpv1 "<< tempv1[0] <<" "<<tempv1[1] <<" "<<tempv1[2] <<" "<<endl;
            //std::cout<<dis_new << " " << isnan(dis_new)<<endl;
            if (isnan(dis_new)){
                std::cout << "error: nan" << std::endl;
                std::exit(-1);
            }
            //subvec(tempv1,g.v[g.vidtoindex[vi]].co,dis_vector);
            // dis_new=g.find_project_dist_axes(g.heidtoindex[heidtemp], tempv1, g.v[g.vidtoindex[vi]].co , dis_vector );
            //############# dist remove
            /*ofstream myfile;
            myfile.open ("removedist.txt", ios::out | ios::app);
            myfile <<dis_new<<endl;
            myfile.close();*/
            //std::cout<<x<<endl;
            delete[] tempv1;

            //double vp = 1/(gsl_ran_gaussian_pdf(dis_new, g.gaussian_sigma));

            double vp = pow((sqrt(2*M_PI)*g.gaussian_sigma),3)/( exp(-((dis_new*dis_new)/(2*g.gaussian_sigma*g.gaussian_sigma))) );

            double crit = exp(-de / g.T) / (2 * vp);
            //delete[] dis_vector;
            //std::cout << " crit is " << crit << endl;
            if (gsl_rng_uniform(rg) < crit)
            {
                //std::cout << " 008 g.Nd is " <<g.Nd<<endl;
                int nextidboundary0 = g.he[g.heidtoindex[heid_next_boundary]].nextid_boundary;
                int previdboundary0 = g.he[g.heidtoindex[heid0]].previd_boundary;
                //g.Nd-=g.he[heindex0].din;
                //g.Nd-=g.he[heopindex0].din;
                //g.Nd-=g.he[heindex_next_boundary].din;
                //g.Nd-=g.he[prevopindex0].din;
                //std::cout << "REMOVING DIMER four drugs removed " << g.he[heopindex0].din + g.he[heindex0].din +  g.he[heindex_next_boundary].din + g.he[prevopindex0].din<<endl;
                vector<int> vecupdate;
                if (g.v[g.vidtoindex[vi]].vneigh.size() > 0)
                {
                    for (vector<int>::iterator it = g.v[g.vidtoindex[vi]].vneigh.begin(); it != g.v[g.vidtoindex[vi]].vneigh.end(); ++it)
                    {
                        vecupdate.push_back(*it);
                    }
                }
                /*
                int still_in_vout = check_vid_in_vout(g, vi);
                if(still_in_vout>2){
                    std::cout << "WARNING: about to delete edges with a vertex that has other edges!!" << std::endl;
                    std::cout << still_in_vout << std::endl;
                    g.get_obs().dump_lammps_traj_dimers(g, 1000000);
                }
                */
                int success = g.delete_edge(heid0);
                if (success < 0)
                {
                    std::cout << "!!!heid_next_boundary-0" << endl;
                    std::exit(-1);
                }
//                std::cout << "updating index after removing dimer (l 2335)" << std::endl;
                g.update_index();
                success = g.delete_edge(heid_next_boundary); // next of op
                if (success < 0)
                {
                    std::cout << "!!!heid_next_boundary-2" << endl;
                    std::exit(-1);
                }
                /*
                still_in_vout = check_vid_in_vout(g, vi);
                if(still_in_vout>0){
                    std::cout << "WARNING: about to delete vertex that is still a vout!!" << std::endl;
                    std::cout << still_in_vout << std::endl;
                    g.get_obs().dump_lammps_traj_dimers(g, 1000001);
                }
                */
                //std::cout << "updating index after removing dimer (l 2343)" << std::endl;
                g.update_index();
                //std::cout << "checking verts before deleting vertex" << std::endl;
                check_vout_in_vid(g);
		//241125 between here and next g.update_index() something is going wrong
		//FIXED: check whether the vertex is still associated with other half edges, if so don't delete it
		        if(check_vid_in_vout(g,vi)==0){
                    //std::cout << "deleting vertex " << vi << std::endl;
                    int x = g.delete_vertex(vi);
                    del_vert_counter++;
                    //std::cout << "Deleted "<< del_vert_counter << " times" << std::endl;
                    //std::cout << "Checking verts on l 2346" << std::endl;
                    check_vout_in_vid(g);
                    if (x < 0)
                    {
                        std::cout << "error deleting vertex" << std::endl;
                        std::exit(-1);
                    }
                }
                else{
                    //Make sure vertex neighbors get updated after deleting intervening edges
                    //std::cout << "updating vertex neighbors after deleting edges" << std::endl;
                    no_del_vertex_counter++;

                    g.update_geometry_vertex(g.vidtoindex[vi]);
                    g.update_normals_vertex(g.vidtoindex[vi]);
                    g.update_excluder_top_vertex(g.vidtoindex[vi]);
                    g.update_neigh_vertex_and_neigh(vi);
                }
                g.set_prev_next(nextopid0, -1, -1); //prev of op
                //std::cout << "Checking verts on l 2352" << std::endl;
                //check_vout_in_vid(g);
                g.he[g.heidtoindex[nextopid0]].boundary_index = bi;
                //std::cout << "Checking verts on l 2356" << std::endl;
                //check_vout_in_vid(g);
                g.set_prev_next_boundary(previdboundary0, nextopid0);
                //std::cout << "Checking verts on l 2359" << std::endl;
                //check_vout_in_vid(g);
                g.set_prev_next_boundary(nextopid0, nextidboundary0);
                //std::cout << "Checking verts on l 2362" << std::endl;
                //check_vout_in_vid(g);
                //update new edges nex_boundary

                //g.update_half_edge(nextopid0);
                //g.update_half_edge(heid_prev_boundary);

                //multvec(g.v[g.vidtoindex(vi)].co,1,vco);

                //double de= g.compute_energy()-e1;
                //gbb=gb0next+gbnextprev+gb0prev;

                //delete[] vco;
                //std::cout << "updating index after deleting vertex in remove dimer (l 2366)" << std::endl;
                g.update_index(); //241125 getting segfault here
                for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
                {
                    g.update_neigh_vertex(*it);
                }

                //ToDo
                //update new edges boundary index

                //std::cout << "00 remove dimer with previous accepted"<<endl;
                vecupdate.clear();
                return 2;
            }
            else
            { //removel not succesful
                //std::cout << "00 remove dimer with previous not accepted"<<endl;
                return 0;
            }
        }
        //delete[] vco;
    }

    return -1;
}

int MC::attempt_remove_monomer_dimer_drug(System &g, int heid0) /* 102220 THIS NEEDS UPDATE -> MOVE GEOMETRY TO geometry! */
{
    //std::cout << "in attempt_remove_monomer_dimer" << endl;
    //std::cout << "Nd is " <<g.Nd <<endl;

    if (!(g.is_boundary(heid0) > 0))
    {
        std::cout << "not on boundary cannot remove" << endl;
        std::exit(-1);
    }

    g.update_half_edge(heid0);
    int heindex0 = g.heidtoindex[heid0]; // index of this edge
    //int heid0type=g.he[heindex0].type;
    if ((g.he[heindex0].nextid != -1) || (g.he[heindex0].previd != -1))
    {
        //std::cout << " has next or previous cannot remove" << endl;
        return -1;
    }

    int bi = g.he[heindex0].boundary_index;
    if (bi == -1)
    {
        std::cout << "Error in boundary index in remove_monomer_dimer" << endl;
        std::exit(-1);
    }

    g.update_half_edge(g.he[heindex0].opid);

    int heopindex0 = g.heidtoindex[g.he[heindex0].opid]; // indexd of opposite edge
    int optype = g.he[heopindex0].type;
    int nextopid0 = g.he[heopindex0].nextid; // id of next of opposite edge
    int prevopid0 = g.he[heopindex0].previd;
    //std::cout << " nextopid0 " << nextopid0 << " prevopid0 " <<prevopid0 <<endl;
    if ((nextopid0 == -1) || (prevopid0 == -1))
    {
        std::cout << "wrong geometry" << endl;
        std::exit(-1);
    }
    int nextopindex0 = g.heidtoindex[nextopid0]; // id of prev of opposite edge
    int prevopindex0 = g.heidtoindex[prevopid0];
    int opnexttype = g.he[nextopindex0].type;
    int opprevtype = g.he[prevopindex0].type;

    if ((g.he[heindex0].din == 1))
    {
        //std::cout << " din=1 "<<endl;
        return (-1);
    }                                    // WITH DRUG NO REMOVAL
    int heid_prev_boundary = g.he[nextopindex0].opid; // now back to this side // ToDo this should be previd_boundary
    int heid_next_boundary = g.he[prevopindex0].opid; // after vertex // ToDo this should be nextid_boundary

    double gbb = 0;

    // remove monomer
    // first try remove monomer, if the edge is not a wedge, remove dimer
    if (g.is_boundary(heid_prev_boundary) < 0 && g.is_boundary(heid_next_boundary) < 0)
    {
     if ((g.is_vboundary(g.he[nextopindex0].vout) < 0 || g.Nboundary != 1) && \
       (g.v[g.vidtoindex[g.he[heindex0].vout]].doubleboundary==-1 && g.v[g.vidtoindex[g.he[heindex0].vin]].doubleboundary==-1 ) )// if Nboundary>0 allow for double boundary;
    {                                                                                                                               //delete monomer
        return -1;                                                                                                                               //delete monomer
        //std::cout << "in deleting monomer" <<endl;
        //if (g.is_vboundary(g.he[g.heidtoindex[nextopid0]].vout) > 0) { std::cout << "v on boundary wrong geometry!"<<endl; std::exit(-1);}
        /*if ((g.he[nextopindex0].din == 1) || g.he[prevopindex0].din == 1)
        {
            return (-1);
        }

        double de = -(g.stretch_energy(heindex0));
        int nextboundary0 = g.he[heindex0].nextid_boundary;
        int prevboundary0 = g.he[heindex0].previd_boundary;
        //if (g.he[heindex0].previd!=-1) { de-=g.dimer_bend_energy(g.get_heindex(g.he[heindex0].previd)); }
        //if (g.he[heindex0].nextid!=-1) {de-=g.dimer_bend_energy(heindex0); }
        de -= (g.dimer_bend_energy(heopindex0) + g.dimer_bend_energy(prevopindex0));
        //de -= g.bend_energy(nextopindex0) + g.bend_energy(prevopindex0); //g.monomer_energy(heid0);

        gbb += g.find_dg(opprevtype, optype, g.he[heopindex0].din);
        gbb += g.find_dg(optype, opnexttype, g.he[nextopindex0].din);

        de -= (gbb - g.mu[g.he[heindex0].type]);
        double crit = 2 * exp(-de / g.T); ///(2.0*g.z*g.K*g.K*g.K);

        if (g.Test_assembly == 1)
        {
            std::cout << "crit is " << crit << endl;
            crit = 1;
        }
        if (gsl_rng_uniform(r) < crit)
        {
            //g.Nd-=g.he[heopindex0].din;
            //g.Nd-=g.he[heindex0].din;
            //std::cout << "REMOVING MONOMER two drugs removed " << g.he[heopindex0].din + g.he[heindex0].din <<endl;

            //std::cout << "de is" <<de <<endl;
            //std::cout << "crit is  " << crit <<endl;
            //if ( g.he[heindex0].previd!=-1) {g.he[g.get_heindex(g.he[heindex0].previd)].nextid=-1;}
            //if (g.he[heindex0].nextid!=-1) {g.he[g.get_heindex(g.he[heindex0].nextid)].previd=-1; }
            //std::cout << "g.he[heopindex0].previ " << g.he[heopindex0].previd <<endl;
            //std::cout << "g.he[prevopindex0].nextid " << g.he[prevopindex0].nextid <<endl;
            //std::cout << "g.he[heopindex0].nextid " << g.he[heopindex0].nextid <<endl;
            //std::cout << "g.he[nextopindex0].previd" << g.he[nextopindex0].previd <<endl;
            if (g.he[heopindex0].previd != -1)
            {
                g.he[prevopindex0].nextid = -1;
            }
            if (g.he[heopindex0].nextid != -1)
            {
                g.he[nextopindex0].previd = -1;
            }

            //ToDo -> Done
            //update new edges boundary index
            //update new edges nex_boundary
            g.he[nextopindex0].boundary_index = bi;
            g.he[prevopindex0].boundary_index = bi;

            g.set_prev_next_boundary(nextopid0, prevopid0);
            g.set_prev_next_boundary(prevboundary0, nextopid0);
            g.set_prev_next_boundary(prevopid0, nextboundary0);

            if (g.is_vboundary(g.he[nextopindex0].vout) > 0)
                g.v[g.vidtoindex[g.he[nextopindex0].vout]].doubleboundary = 1;

            int vidin = g.he[heindex0].vin;
            int vidout = g.he[heindex0].vout;

            int x = g.delete_edge(heid0);
            //g.update_edge();_

            if (x < 0)
            {
                std::cout << "could not delete " << endl;
                std::exit(-1);
            }
            //

            //g.set_prev_next(nextopid0, -1, prevopid0);
            //g.set_prev_next(prevopid0, nextopid0, -1);
            //std::cout << "00 monomer removed" <<endl;
            //std::cout << " 004 g.Nd is " <<g.Nd<<endl;
            g.update_index();
            g.update_neigh_vertex(vidin);
            g.update_neigh_vertex(vidout);

            //std::cout << " 005 g.Nd is " <<g.Nd<<endl;

            return 1;
        }
        else
        {
            //std::cout << " remove monomer not accepted  " << endl;
            //g.add_monomer(nextopid0, prevopid0,optype);

            return -1;
        }*/
    }
    }
            else // edge is not wedge so remove dimer //with previous or next
    {

        int heindex_prev_boundary = g.heidtoindex[heid_prev_boundary];
        int heindex_next_boundary = g.heidtoindex[heid_next_boundary];
        //double *vco=new double[3];
        //remove dimer
        if (g.is_boundary(heid_prev_boundary) > 0 && g.he[heindex_prev_boundary].previd == -1 && g.v[g.vidtoindex[g.he[heindex0].vin]].doubleboundary==-1) //delete dimer this and next (inside)(nextopindex) / this and prev (on boundary)  there should be nno bonds between this and previous
        {
            //std::cout << "00 remove dimer with next (previd_boundary) "<<endl;  //remove this and previd_boundary=heid_prev_boundary
            if (g.he[g.heidtoindex[nextopid0]].din == 1 )
            {
            if (g.he[heindex0].type == 0 || g.he[heindex0].type==3)
            {
            if(g.he[heindex_prev_boundary].type==0 || g.he[heindex_prev_boundary].type==3)
            {
            

            // int yid=g.he[heindex0].nextid;
            int vi = g.he[heopindex0].vout;

            //Reject moves that attempt to remove dimers with more than one vertex
            int still_in_vout = check_vid_in_vout(g, vi);
            if(still_in_vout>2){
                std::cout << "Rejecting attempt to remove dimer where vertex is shared by others" << std::endl;
                return (-1);
            }

            //gbb+= g.find_gbb(optype,opnexttype,opprevtype); //whole triangle
            gbb += g.find_dg(optype, opnexttype, g.he[g.heidtoindex[nextopid0]].din);
            gbb += g.find_dg(opnexttype, opprevtype, g.he[prevopindex0].din);
            gbb += g.find_dg(opprevtype, optype, g.he[heopindex0].din);
            double de = -(g.stretch_energy(heindex0) + g.stretch_energy(g.heidtoindex[heid_prev_boundary]));
            de -= g.bend_energy(g.heidtoindex[heid_next_boundary]);

            de -= (g.dimer_bend_energy(heopindex0));
            de -= (g.dimer_bend_energy(nextopindex0) + g.dimer_bend_energy(prevopindex0));
            //if (yid!=-1) {  de-=g.dimer_bend_energy(heindex0);}

            //if (xid!=-1 && g.he[xindex].nextid==heid_prev_boundary) {
            //    gbb+=g.find_dg( g.he[xindex].type,g.he[g.heidtoindex[heid_prev_boundary]].type);
            //     de-=g.dimer_bend_energy(xindex);
            //xidconnect=1;
            //}

            
            /* test !!!!
                if (abs(de-(g.compute_energy()-e11)>0.0000000001)) 
                {//std::cout << "in remove dimer 111 des don't match" << de-(g.compute_energy()-e11) << endl; 
                std::cout << "e11 "<< e11 <<endl; std::cout << "de "<< de <<endl; std::cout << "g.compute_energy() "<< g.compute_energy() <<endl;  
                std::cout << "g.compute_energy() -e11 "<< g.compute_energy()-e11 <<endl; std::exit(-1);} */

            de -= (gbb - (g.mu[g.he[heindex0].type] + g.mu[g.he[heindex_prev_boundary].type]));
            //
            //**************************
            // gaussian correction
            //**************************
            double *tempv1 = new double[3];
            //double *dis_vector=new double[3];
            int heidtemp=g.he[nextopindex0].nextid;
            g.new_vertex_edge(g.heidtoindex[heidtemp], tempv1,g.he[heopindex0].type); // heopindex0 = g.heidtoindex[g.he[g.heidtoindex[heidtemp]].nextid]
            //std::cout<<"tmpv1 "<< tempv1[0] <<" "<<tempv1[1] <<" "<<tempv1[2] <<" "<<endl;
            double dis_new=veclen(tempv1,g.v[g.vidtoindex[vi]].co);
            //std::cout<<dis_new << " " << isnan(dis_new)<<endl;
            //if (isnan(dis_new)) std::exit(-1);
            //subvec(tempv1,g.v[g.vidtoindex[vi]].co,dis_vector);
            //double dis_new=g.find_project_dist_axes(g.heidtoindex[heidtemp], tempv1, g.v[g.vidtoindex[vi]].co , dis_vector );
            //############# dist remove
            /*ofstream myfile;
            myfile.open ("removedist.txt", ios::out | ios::app);
            myfile <<dis_new<<endl;
            myfile.close();*/
            //std::cout<<x<<endl;
            delete[] tempv1;
            //double vp = 1/(gsl_ran_gaussian_pdf(dis_new, g.gaussian_sigma));
            double vp = pow((sqrt(2*M_PI)*g.gaussian_sigma),3)/( exp(-((dis_new*dis_new)/(2*g.gaussian_sigma*g.gaussian_sigma))) );
            double crit = exp(-de / g.T) / (2 * vp); 
            //std::cout << " crit is " << crit << endl;
            if (gsl_rng_uniform(rg) < crit) //delete dimer this and next (inside)(nextopindex) / this and prev (on boundary)
            {
                int nextidboundary0 = g.he[g.heidtoindex[heid0]].nextid_boundary;
                int previdboundary0 = g.he[g.heidtoindex[heid_prev_boundary]].previd_boundary;

                vector<int> vecupdate;
                if (g.v[g.vidtoindex[vi]].vneigh.size() > 0)
                {
                    for (vector<int>::iterator it = g.v[g.vidtoindex[vi]].vneigh.begin(); it != g.v[g.vidtoindex[vi]].vneigh.end(); ++it)
                    {
                        vecupdate.push_back(*it);
                    }
                }
                //std::cout << " 006 g.Nd is " <<g.Nd<<endl;

                //WITHDRUG No removal
                g.Nd-=g.he[heindex0].din;
                g.Nd-=g.he[heopindex0].din;
                g.Nd-=g.he[heindex_prev_boundary].din;
                g.Nd-=g.he[nextopindex0].din;

                //std::cout << "REMOVING DIMER four drugs removed " << g.he[heopindex0].din + g.he[heindex0].din +  g.he[heindex_prev_boundary].din + g.he[nextopindex0].din<<endl;
                //std::cout << "removing heid0 and heid_prev_boundary " << heid0 <<" " << heid_prev_boundary <<endl;
                int success = g.delete_edge(heid0);
                if (success < 0)
                {
                    std::cout << "!!!" << endl;
                    std::exit(-1);
                }
//                std::cout << "updating index after removing dimer w drug (l 2650)" << std::endl;
                g.update_index();
                success = g.delete_edge(heid_prev_boundary); // next of op
                if (success < 0)
                {
                    std::cout << "!!!" << endl;
                    std::exit(-1);
                }
//                std::cout << "updating index after removing dimer w drug (l 2658)" << std::endl;
                g.update_index();
                int x = g.delete_vertex(vi);
                if (x < 0)
                {
                    std::cout << "error deleting vertex" << std::endl;
                    std::exit(-1);
                }

                g.set_prev_next(prevopid0, -1, -1); //prev of op
                g.he[g.heidtoindex[prevopid0]].boundary_index = bi;
                g.set_prev_next_boundary(prevopid0, nextidboundary0);
                g.set_prev_next_boundary(previdboundary0, prevopid0);
                //update new edges nex_boundary

//                std::cout << "updating index after removing vertex in dimer w drug (l 2672)" << std::endl;
                g.update_index();
                for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
                {
                    g.update_neigh_vertex(*it);
                }

                //ToDo
                //update new edges boundary index

                //std::cout << "00 remove dimer with next accepted "<<endl;
                vecupdate.clear();
                return 2;
            }
            // no need for anything if removal cannot happen
            else
            {
                //std::cout << "00 remove dimer with previous not accepted "<<endl;
                return 0;
            }
        }
        }
        }
        }
        else if (g.is_boundary(heid_next_boundary) > 0 && g.he[heindex_next_boundary].nextid == -1) //delete dimer this and previous (inside)(prevopindex) / this and next (on boundary) there should be no bonds between next_boudary and its next
        {                                                                 //delete dimer this and prev

            //std::cout << "00 remove dimer with previous (heid_next_boundary) "<<endl; //with opid_boundary=heid_next_boundary
            //if (g.he[g.heidtoindex[heid_next_boundary]].din==1) return -1;
            //if (g.he[heindex_next_boundary].din==1 || g.he[g.heidtoindex[prevopid0]].din==1 ) { return(-1);} //WITHDRUG no removal
            if (g.he[heopindex0].din == 1)
            {
            if (g.he[heindex0].type == 0 || g.he[heindex0].type==3)
            {
            if(g.he[heindex_next_boundary].type==0 || g.he[heindex_next_boundary].type==3)
            {
            
             //NOREMOVAL DRUG
            int vi = g.he[heopindex0].vin;

            //Reject moves that attempt to remove dimers with more than one vertex
            int still_in_vout = check_vid_in_vout(g, vi);
            if(still_in_vout>2){
                std::cout << "Rejecting attempt to remove dimer where vertex is shared by others" << std::endl;
                return (-1);
            }

            //gbb+= g.find_gbb(optype,opnexttype,opprevtype); //whole triangle
            gbb += g.find_dg(optype, opnexttype, g.he[g.heidtoindex[nextopid0]].din);
            gbb += g.find_dg(opnexttype, opprevtype, g.he[g.heidtoindex[prevopid0]].din);
            gbb += g.find_dg(opprevtype, optype, g.he[heopindex0].din);
            double de = -(g.stretch_energy(heindex0) + g.stretch_energy(g.heidtoindex[heid_next_boundary]));
            de -= g.bend_energy(g.heidtoindex[heid_prev_boundary]);

            de -= (g.dimer_bend_energy(heopindex0));
            de -= (g.dimer_bend_energy(g.heidtoindex[nextopid0]) + g.dimer_bend_energy(g.heidtoindex[prevopid0]));

            de -= (gbb - (g.mu[g.he[heindex0].type] + g.mu[g.he[heindex_next_boundary].type]));

            //**************************
            // gaussian correction
            //**************************
            double *tempv1 = new double[3];
            //double *dis_vector=new double[3];
            int heidtemp=g.he[heopindex0].nextid;
            g.new_vertex_edge(g.heidtoindex[heidtemp], tempv1,g.he[prevopindex0].type); // prevopindex = ? g.heidtoindex[g.he[g.heidtoindex[heidtemp]].nextid]
            double dis_new=veclen(tempv1,g.v[g.vidtoindex[vi]].co);
            //std::cout<<"tmpv1 "<< tempv1[0] <<" "<<tempv1[1] <<" "<<tempv1[2] <<" "<<endl;
            //std::cout<<dis_new << " " << isnan(dis_new)<<endl;
            if (isnan(dis_new)){
                std::cout << "error: nan" << std::endl;
                std::exit(-1);
            }
            //subvec(tempv1,g.v[g.vidtoindex[vi]].co,dis_vector);
            // dis_new=g.find_project_dist_axes(g.heidtoindex[heidtemp], tempv1, g.v[g.vidtoindex[vi]].co , dis_vector );
            //############# dist remove
            /*ofstream myfile;
            myfile.open ("removedist.txt", ios::out | ios::app);
            myfile <<dis_new<<endl;
            myfile.close();*/
            //std::cout<<x<<endl;
            delete[] tempv1;

            //double vp = 1/(gsl_ran_gaussian_pdf(dis_new, g.gaussian_sigma));

            double vp = pow((sqrt(2*M_PI)*g.gaussian_sigma),3)/( exp(-((dis_new*dis_new)/(2*g.gaussian_sigma*g.gaussian_sigma))) );

            double crit = exp(-de / g.T) / (2 * vp);
            //delete[] dis_vector;
            //std::cout << " crit is " << crit << endl;
            if (gsl_rng_uniform(rg) < crit)
            {
                //std::cout << " 008 g.Nd is " <<g.Nd<<endl;
                int nextidboundary0 = g.he[g.heidtoindex[heid_next_boundary]].nextid_boundary;
                int previdboundary0 = g.he[g.heidtoindex[heid0]].previd_boundary;
                g.Nd-=g.he[heindex0].din;
                g.Nd-=g.he[heopindex0].din;
                g.Nd-=g.he[heindex_next_boundary].din;
                g.Nd-=g.he[prevopindex0].din;
                //std::cout << "REMOVING DIMER four drugs removed " << g.he[heopindex0].din + g.he[heindex0].din +  g.he[heindex_next_boundary].din + g.he[prevopindex0].din<<endl;
                vector<int> vecupdate;
                if (g.v[g.vidtoindex[vi]].vneigh.size() > 0)
                {
                    for (vector<int>::iterator it = g.v[g.vidtoindex[vi]].vneigh.begin(); it != g.v[g.vidtoindex[vi]].vneigh.end(); ++it)
                    {
                        vecupdate.push_back(*it);
                    }
                }
                int success = g.delete_edge(heid0);
                if (success < 0)
                {
                    std::cout << "!!!heid_next_boundary-0" << endl;
                    std::exit(-1);
                }
//                std::cout << "updating index after removing dimer (l 2776)" << std::endl;
                g.update_index();
                success = g.delete_edge(heid_next_boundary); // next of op
                if (success < 0)
                {
                    std::cout << "!!!heid_next_boundary-2" << endl;
                    std::exit(-1);
                }
//                std::cout << "updating index after removing dimer (l 2784)" << std::endl;
                g.update_index();
                int x = g.delete_vertex(vi);
                if (x < 0)
                {
                    std::cout << "error deleting vertex" << std::endl;
                    std::exit(-1);
                }
                g.set_prev_next(nextopid0, -1, -1); //prev of op
                g.he[g.heidtoindex[nextopid0]].boundary_index = bi;
                g.set_prev_next_boundary(previdboundary0, nextopid0);
                g.set_prev_next_boundary(nextopid0, nextidboundary0);
                //update new edges nex_boundary

                //g.update_half_edge(nextopid0);
                //g.update_half_edge(heid_prev_boundary);

                //multvec(g.v[g.vidtoindex(vi)].co,1,vco);

                //double de= g.compute_energy()-e1;
                //gbb=gb0next+gbnextprev+gb0prev;

                //delete[] vco;
//                std::cout << "updating index after removing vertex in dimer (l 2806)" << std::endl;
                g.update_index();
                for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
                {
                    g.update_neigh_vertex(*it);
                }

                //ToDo
                //update new edges boundary index

                //std::cout << "00 remove dimer with previous accepted"<<endl;
                vecupdate.clear();
                return 2;
            }
            else
            { //removel not succesful
                //std::cout << "00 remove dimer with previous not accepted"<<endl;
                return 0;
            }
        }
        }
        }
        }
        //delete[] vco;
    }

    return -1;
}

int MC::attempt_wedge_fusion(System &g)
{
    //ToDo : update it to read pairs from a vector of pairs
    //std::cout << "in attempt wedge Fusion  " << endl;
    
//    std::cout << "updating index in attempt wedge fusion (l 2839)" << std::endl;
    g.update_index();

    g.update_fusion_pairs_he();

    /*for (vector<int>::iterator it = g.boundary.begin(); it != g.boundary.end(); ++it)
        {
            int heindex00=g.heidtoindex[*it];
            std::cout<< " he is " << *it << " he_prev_wedge_fusion_heid " << g.he[heindex00].prev_wedge_fusion_heid << endl;
            std::cout<< " he is " << *it << " he_next_wedge_fusion_heid " << g.he[heindex00].next_wedge_fusion_heid << endl;
            std::cout<< " he is " << *it << " he_prev_fusion_heid " << g.he[heindex00].prev_fusion_heid << endl;
            std::cout<< " he is " << *it << " he_next_fusion_heid " << g.he[heindex00].next_fusion_heid << endl;
        }*/

    if (g.fusionwedgehe.size() == 0)
    {
        return -1;
    }
    /*for (vector<int>::iterator it = g.fusionwedgehe.begin(); it != g.fusionwedgehe.end(); ++it)
        {
            int heindex00=g.heidtoindex[*it];
            if (g.he[heindex00].prev_wedge_fusion_heid!=-1) std::cout<< " he is " << *it << " he_prev_wedge_fusion_heid " << g.he[heindex00].prev_wedge_fusion_heid << endl;
            if (g.he[heindex00].next_wedge_fusion_heid!=-1) std::cout<< " he is " << *it << " he_next_wedge_fusion_heid " << g.he[heindex00].next_wedge_fusion_heid << endl;
            if (g.he[heindex00].prev_fusion_heid!=-1) std::cout<< " he is " << *it << " he_prev_fusion_heid " << g.he[heindex00].prev_fusion_heid << endl;
            if (g.he[heindex00].next_fusion_heid!=-1) std::cout<< " he is " << *it << " he_next_fusion_heid " << g.he[heindex00].next_fusion_heid << endl;
            //dump_restart_lammps_data_file(g,11111111);
            //std::cout << "exit for now"<<endl;
            //std::exit(-1);
        }
    for (vector<int>::iterator it = g.fusionwedgehe.begin(); it != g.fusionwedgehe.end(); ++it)
        {
            int heindex00=g.heidtoindex[*it];
            if (g.he[heindex00].prev_wedge_fusion_heid!=-1) std::cout<< " In fusion he is " << *it << " he_prev_wedge_fusion_heid " << g.he[heindex00].prev_wedge_fusion_heid << endl;
            if (g.he[heindex00].next_wedge_fusion_heid!=-1) std::cout<< " In fusion he is " << *it << " he_next_wedge_fusion_heid " << g.he[heindex00].next_wedge_fusion_heid << endl;
            if (g.he[heindex00].prev_fusion_heid!=-1) std::cout<< " In fusion he is " << *it << " he_prev_fusion_heid " << g.he[heindex00].prev_fusion_heid << endl;
            if (g.he[heindex00].next_fusion_heid!=-1) std::cout<< " In fusion he is " << *it << " he_next_fusion_heid " << g.he[heindex00].next_fusion_heid << endl;
            //dump_restart_lammps_data_file(g,11111111);
            //std::cout << "exit for now"<<endl;
            //std::exit(-1);
        }*/

    int ind = gsl_rng_uniform_int(rg, g.fusionwedgehe.size());

    int heid0 = g.fusionwedgehe[ind];
    //std::cout << "00 heid for fusion is " << heid0 <<endl;
    int heindex0 = g.heidtoindex[heid0];
    int heindex_next = -1;
    int heindex_prev = -1;
    int vidi = -1;
    int vidj = -1;
    if ((g.he[heindex0].prev_wedge_fusion_heid == -1) && (g.he[heindex0].next_wedge_fusion_heid == -1))
    {
        std::cout << "why wrong choice in wege fusion" << endl;
        std::exit(-1);
        return -1;
    }
    else if ((g.he[heindex0].prev_wedge_fusion_heid != -1) && (g.he[heindex0].next_wedge_fusion_heid != -1))
    { // if has both choose randomly
        if (gsl_rng_uniform(rg) < .5)
        { // with prev
            heindex_next = heindex0;
            heindex_prev = g.heidtoindex[g.he[heindex_next].prev_wedge_fusion_heid];
            vidi = g.he[heindex_next].vin;
            vidj = g.he[heindex_prev].vout;
        }
        else
        { // with next
            heindex_prev = heindex0;
            heindex_next = g.heidtoindex[g.he[heindex_prev].next_wedge_fusion_heid];
            vidi = g.he[heindex_next].vin;
            vidj = g.he[heindex_prev].vout;
        }
    }
    else if (g.he[heindex0].prev_wedge_fusion_heid != -1)
    {
        heindex_next = heindex0;
        heindex_prev = g.heidtoindex[g.he[heindex_next].prev_wedge_fusion_heid];
        vidi = g.he[heindex_next].vin;
        vidj = g.he[heindex_prev].vout;
    }
    else if (g.he[heindex0].next_wedge_fusion_heid != -1)
    {
        heindex_prev = heindex0;
        heindex_next = g.heidtoindex[g.he[heindex_prev].next_wedge_fusion_heid];
        vidi = g.he[heindex_next].vin;
        vidj = g.he[heindex_prev].vout;
    }
    int heid_next = g.he[heindex_next].id;
    int heid_prev = g.he[heindex_prev].id;

    int nextidboundary0 = g.he[heindex_prev].nextid_boundary;
    int previdboundary0 = g.he[heindex_next].previd_boundary;

    int bi = g.he[heindex_next].boundary_index;
    //std::cout << "wedge fusion pair chosen" <<endl;

    int heidm = g.he[heindex_next].nextid;
    int bondnextm = 1;
    int bondmprev = -1;
    if (heidm == -1)
    {
        heidm = g.he[heindex_prev].previd;
        bondnextm = -1;
        bondmprev = 1;
    }
    if (heidm == -1)
    {
        std::cout << "WRONG Geometry in wedge fusion" << endl;
        std::exit(-1);
    }
    // now save th status
    //std::cout << "heid_prev " << heid_prev << " heid_next " << heid_next<< "heidm" << heidm <<endl;

    double e1 = 0;

    //e1 += g.vertex_energy(vidi);
    //e1 += g.vertex_energy(vidj);
    e1 = g.compute_energy();

    //save vecupdate old neighbors except vidi , vidj
    int vindexi = g.vidtoindex[vidi];
    int vindexj = g.vidtoindex[vidj];

    vector<int> vecupdate;
    for (vector<int>::iterator it = g.v[vindexj].vneigh.begin(); it != g.v[vindexj].vneigh.end(); ++it)
    {
        if (*it != vidi)
            vecupdate.push_back(*it);
    }

    for (vector<int>::iterator it = g.v[vindexi].vneigh.begin(); it != g.v[vindexi].vneigh.end(); ++it)
    {
        if (*it != vidj)
            vecupdate.push_back(*it);
    }

    //save old vertex i
    VTX *vtxi;
    vtxi = new VTX;
    g.save_vtx(vidi, vtxi);

    //save old vertex j
    VTX *vtxj;
    vtxj = new VTX;
    g.save_vtx(vidj, vtxj);

    double *tempv = new double[3];
    double *newv = new double[3];
    centvec(g.v[vindexi].co, g.v[vindexj].co, newv); //tempv if random
    //no random position
    //g.move_p(tempv, newv, r);

    //add_new vertex  vid is newvid
    g.add_vertex(newv);
    int newvid = g.Nvlast - 1;
    //update index
//    std::cout << "updating index after adding vertex in wedge fusion (l 2996)" << std::endl;
    g.update_index(); // work ids after this
    //std::cout << " newvid is " << newvid << " vidtoindex[newvid] is " << g.vidtoindex[newvid] << " Nv is " << g.Nv << endl;
    delete[] tempv;
    delete[] newv;

    //std::cout << "new vertex added" << endl;

    //merging vertices
    //std::cout<< " After addition g.Nv " << g.Nv <<endl;
    int newvindex = g.vidtoindex[newvid];
    vindexi = g.vidtoindex[vidi];
    vindexj = g.vidtoindex[vidj];
    for (vector<int>::iterator ithe = g.v[vindexi].hein.begin(); ithe != g.v[vindexi].hein.end(); ++ithe)
    {

        g.v[newvindex].hein.push_back(*ithe);
        int heindex0 = g.heidtoindex[*ithe];
        g.he[heindex0].vout = newvid;
        g.he[g.heidtoindex[g.he[heindex0].opid]].vin = newvid;
        g.update_half_edge(g.he[heindex0].id);
        g.update_half_edge(g.he[heindex0].opid);
    }

    for (vector<int>::iterator ithe = g.v[vindexj].hein.begin(); ithe != g.v[vindexj].hein.end(); ++ithe)
    {

        g.v[newvindex].hein.push_back(*ithe);
        int heindex0 = g.heidtoindex[*ithe];
        g.he[heindex0].vout = newvid;
        g.he[g.heidtoindex[g.he[heindex0].opid]].vin = newvid;
        g.update_half_edge(g.he[heindex0].id);
        g.update_half_edge(g.he[heindex0].opid);
    }

    if (vindexi > vindexj)
    {
        g.delete_vertex(vidi);
        g.delete_vertex(vidj);
    }
    else
    {
        g.delete_vertex(vidj);
        g.delete_vertex(vidi);
    }

    //std::cout<< " After delete g.Nv " << g.Nv <<endl;
//    std::cout << "updating index after deleting vertex in wedge fusion (l 3043)" << std::endl;
    g.update_index();

    //can do only for those not connected
    //g.he[g.heidtoindex[heid_next]].previd=heid_prev;
    //g.he[g.heidtoindex[heid_prev]].nextid=heid_next;
    //g.set_prev_next(hei_, heidprev, heidnext);
    //g.set_prev_next(heidprev, heidnext, heidm);
    //g.set_prev_next(heidnext, heidm, heidprev);
    //std::cout << "in wedge fusion set next prev " << endl;
    //std::cout << " heid_next, heid_prev, heidm "<< heid_next<< " "<< heid_prev << " " << heidm << endl;
    g.set_prev_next(heid_next, heid_prev, heidm);
    //std::cout<<"heid_next next_prev set"<<endl;
    g.set_prev_next(heidm, heid_next, heid_prev);
    //std::cout<<"heidm next_prev set"<<endl;
    g.set_prev_next(heid_prev, heidm, heid_next);
    //std::cout<<"heid_prev next_prev set"<<endl;

    //std::cout << "in wedge fusion next prev is set " << endl;

    //g.update_normals();
    // //std::cout << "in fusion normald updated " << endl;
    //final configuration  // binding energies!
    double e2 = g.find_dg(g.he[g.heidtoindex[heid_prev]].type, g.he[g.heidtoindex[heid_next]].type, g.he[heindex_next].din);

    if (bondnextm < 0)
    {
        e2 += g.find_dg(g.he[g.heidtoindex[heid_next]].type, g.he[g.heidtoindex[heidm]].type, g.he[g.heidtoindex[heidm]].din);
    }
    if (bondmprev < 0)
    {
        e2 += g.find_dg(g.he[g.heidtoindex[heidm]].type, g.he[g.heidtoindex[heid_prev]].type, g.he[g.heidtoindex[heid_prev]].din);
    }
    //std::cout << " now new structure " << endl;
    g.update_normals();
    g.update_excluder_top_vertex(g.vidtoindex[newvid]);

    //std::cout << " now new structure " << endl;
    //int overlap = -1;

    e2 += g.compute_energy(); //g.vertex_energy(newvid);

    // //std::cout << "in fusion normald updated " << endl;
    //final configuration  // binding energies!

    //double vp = 1; //4 / 3. * 4 * M_PI * g.xi * g.xi * g.xi;
    double vp = g.l_thermal_kappa * g.l_thermal_kappa * g.l_thermal_kappa;
    double de = e2 - e1;
    //std::cout <<"de is " <<de <<endl;
    double crit = exp(-de / g.T) * 2 / (vp); //double check
    if (g.Test_assembly == 1)
    {
        crit = 1;
    }
    //std::cout << "crit is "<<crit <<endl;
    //std::cout << " for now neglecting overlap " <<endl;
    int overlapflag = -1;
    //
    /* update neighbors of "neighbors of vtxi and vtxj"  */
    for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
    {
        g.update_neigh_vertex(*it);
    }

    g.update_neigh_vertex(newvid);
    int newvindex0 = g.vidtoindex[newvid];
    if (g.v[newvindex0].vneigh.size() > 0)
    {
        for (vector<int>::iterator it = g.v[newvindex0].vneigh.begin(); it != g.v[newvindex0].vneigh.end(); ++it)
        {
            g.update_neigh_vertex(*it);
            //vecupdate.push_back(*it);
        }
    }

    if (g.check_overlap_g(newvid) < 0)
    {
        overlapflag = 1;
    }
    else
    {
        for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
        {

            if (g.check_overlap_g(*it) < 0)
            {
                overlapflag = 1;
            }
        }
    }
    if (gsl_rng_uniform(rg) < crit && overlapflag == -1)
    {
        std::cout << "wedge fusion  crit is met" << endl;

        //g.check_odd_neigh();
        // update boundary nextid
        g.set_prev_next_boundary(previdboundary0, nextidboundary0);
        vecupdate.clear();
        delete vtxi;
        delete vtxj;
        return 1;
    }
    else
    {
        //std::cout << "fusion not accepted  now move things back to what it was " << endl; //
        //std::exit(-1);
        g.he[g.heidtoindex[heid_prev]].nextid = -1;
        g.he[g.heidtoindex[heid_next]].previd = -1;

        if (bondnextm < 0)
        {

            g.he[g.heidtoindex[heid_next]].nextid = -1;
            g.he[g.heidtoindex[heidm]].previd = -1;
        }
        if (bondmprev < 0)
        {
            g.he[g.heidtoindex[heid_prev]].previd = -1;
            g.he[g.heidtoindex[heidm]].nextid = -1;
        }
        //std::cout<<"attempt_vertex_fusion not accepted" <<endl;
        //std::cout <<"update boundary_index"<<endl;
        g.he[g.heidtoindex[heidm]].boundary_index = bi;
        g.he[g.heidtoindex[heid_prev]].boundary_index = bi;
        g.he[g.heidtoindex[heid_next]].boundary_index = bi;

        //

        //update boundary nextid
        //std::cout <<"update boundary_nextid 0"<<endl;
        g.set_prev_next_boundary(previdboundary0, heid_next);
        //std::cout <<"update boundary_nextid 1"<<endl;
        g.set_prev_next_boundary(heid_next, heidm);

        g.set_prev_next_boundary(heidm, heid_prev);
        //std::cout <<"update boundary_nextid "<<endl;

        g.set_prev_next_boundary(heid_prev, nextidboundary0);

        //std::cout << "fusion not accepted deleting  newvid "<< newvid << endl;
        //vector<int> vecupdate;
        int vindextemp = g.vidtoindex[newvid];
        //std::cout << "newvid  vindex is "<< vindextemp << endl;
        //std::cout << " now delete its neighbors " << endl;
        if (g.v[vindextemp].vneigh.size() > 0)
        {
            for (vector<int>::iterator it = g.v[vindextemp].vneigh.begin(); it != g.v[vindextemp].vneigh.end(); ++it)
            {
                vecupdate.push_back(*it);
            }
        }

        //std::cout << "g.Nv " << g.Nv <<endl;
        //std::cout << "fusion not accepted adding vtxi back" << endl;
        g.add_vertex(vtxi->co);

        int lastindex = g.vidtoindex[g.Nvlast - 1];
        for (vector<int>::iterator ithe = vtxi->hein.begin(); ithe != vtxi->hein.end(); ++ithe)
        {

            int heindex = g.heidtoindex[*ithe];
            g.v[lastindex].hein.push_back(*ithe);
            g.he[heindex].vout = g.v[lastindex].vid;
            g.he[g.heidtoindex[g.he[heindex].opid]].vin = g.v[lastindex].vid;
            g.update_half_edge(g.he[heindex].id);
            g.update_half_edge(g.he[heindex].opid);
        }
        //std::cout << "g.Nv " << g.Nv <<endl;
        //std::cout << "fusion not accepted adding vtxj back" << endl;
        g.add_vertex(vtxj->co);

        lastindex = g.vidtoindex[g.Nvlast - 1];
        for (vector<int>::iterator ithe = vtxj->hein.begin(); ithe != vtxj->hein.end(); ++ithe)
        {
            int heindex = g.heidtoindex[*ithe];
            g.v[lastindex].hein.push_back(*ithe);
            ////std::cout << "updating vout of he " << g.he[heindex].id << "index " << heindex << " *ithe " << *ithe <<endl;
            g.he[heindex].vout = g.v[lastindex].vid;
            g.he[g.heidtoindex[g.he[heindex].opid]].vin = g.v[lastindex].vid;
            g.update_half_edge(g.he[heindex].id);
            g.update_half_edge(g.he[heindex].opid);
        }

        //std::cout <<" deleting ne vertex" <<endl;
        g.delete_vertex(newvid);

        //std::cout <<" deleted , updating index" <<endl;
        //g.update_index();
        //std::cout << "g.Nv " << g.Nv <<endl;
        // std::cout << "fusion not accepted updating index here !!!!!!!!!!!! "<< endl;
//        std::cout << "updating index after deleting vertex in fusion (l 3233)" << std::endl;
        g.update_index();

        for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
        {

            //std::cout <<" *it"<< *it <<endl;
            if (*it != newvid)
            {
                g.update_neigh_vertex(*it);
            }
        }
        g.update_neigh_vertex(g.Nvlast - 1);
        int vindex0 = g.vidtoindex[g.Nvlast - 1];
        if (g.v[vindex0].vneigh.size() > 0)
        {
            for (vector<int>::iterator it = g.v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
            {
                g.update_neigh_vertex(*it);
            }
        }
        g.update_neigh_vertex(g.Nvlast - 2);
        vindex0 = g.vidtoindex[g.Nvlast - 2];
        if (g.v[vindex0].vneigh.size() > 0)
        {
            for (vector<int>::iterator it = g.v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
            {
                g.update_neigh_vertex(*it);
            }
        }

        //g.check_odd_neigh();
        delete vtxi;
        delete vtxj;
        vecupdate.clear();
        return 0;
    }
    return -1;
}

int MC::attempt_wedge_fission(System &g)
{
    //std::cout << "in attempt wedge fission  " << endl;
    int vid0 = -1;

    if (g.boundary.size() == 0)
    {
        std::cout << "what ! attempt_wedge_fission!" << endl;
        std::exit(-1);
    }

    //std::cout << "000000000000000000000000000000000000"<<endl;
    //std::cout << "00 in wedge fission g.boundary.size()" << g.boundary.size() <<endl;
    int ind = gsl_rng_uniform_int(rg, g.boundary.size());

    int heid0 = g.boundary[ind];
    int heindex0 = g.heidtoindex[heid0];

    int bi = g.he[heindex0].boundary_index;

    //std::cout << "00 Starting wedge Fission  heid0 " << heid0 <<" boundary index is " << g.he[heindex0].boundary_index << endl;

    if (gsl_rng_uniform(rg) < .5)
    {
        vid0 = g.he[heindex0].vin;
        //std::cout << " wedge Fission  vin " << vid0 << endl;
        if (g.is_bond_in_boundary(heid0) > 0)
        { //std::cout << "isbound_in"<<endl;
            return (-1);
        }
    }
    else
    {
        vid0 = g.he[heindex0].vout;
        //std::cout << " wedge Fission  vin " << vid0 <<  endl;
        if (g.is_bond_out_boundary(heid0) > 0)
        { //std::cout << "isbound_out"<<endl;
            return (-1);
        }
    }

    if (g.v[g.vidtoindex[vid0]].doubleboundary != -1)
    { //std::cout <<"doubleboundary"<<endl;
        return -1;
    } // vertex should not be double boundary
    // if  other side is also bound retuen -1;
    //if (g.Test_assembly==1) {
    //    vid0=6;
    //}

    int vindex0 = g.vidtoindex[vid0];
    if (g.is_bond_vboundary(vid0) > 0)
    {
        //std::cout << "00 no wedge fission of vertex Bound!  " << vid0 <<endl;
        return (-1);
    }

    int nextidboundary0 = g.v[vindex0].heboundaryoutid;
    int previdboundary0 = g.he[g.heidtoindex[nextidboundary0]].previd_boundary;
    //std::cout << "00 wedge fission of vertex  " << vid0 <<endl;
    //std::cout << "in wedge fission Nv is " << g.Nv << " Nvlast is " << g.Nvlast <<endl;

    if (g.v[vindex0].hein.size() < 4)
    { //std::cout << "hein<4" <<endl;
        return (-1);
    }

    // choose where to break ->xid
    int indhe = gsl_rng_uniform_int(rg, g.v[vindex0].hein.size());

    int xid = g.v[vindex0].hein[indhe];

    int xindex = g.heidtoindex[xid];
    int xidopid = g.he[xindex].opid;
    int xidnextid = g.he[xindex].nextid;
    int xidprevid = g.he[xindex].previd;

    //std::cout << "00 in wedge fission xid is " << xid << "xidopid is" << xidopid  <<endl;
    if (g.is_boundary(xid) > 0 || g.is_boundary(xidopid) > 0)
    {
        //std::cout << "fission return -1 111"<<endl;
        return -1;
    }

    if (xidnextid == -1)
    {
         std::cout << " why xidnextid  is -1" << endl;
        std::exit(-1);
    }
    int xidopidnext = g.he[g.heidtoindex[xidnextid]].opid;
    if (g.is_boundary(xidopidnext) > 0)
    {
        //std::cout << "fission return -1 222"<<endl;
        return -1;
    }

    if (g.is_vboundary(g.he[xindex].vin) > 0 || g.is_vboundary(g.he[g.heidtoindex[xidnextid]].vout) > 0)
    {
        //std::cout << "fission return -1 333"<<endl;
        return -1;
    }

    /* starting fission */

    //std::cout << "fission of vertex  " << vid0 <<endl;

    // saving neighbors of this
    vector<int> vecupdate;
    int vindextemp = g.vidtoindex[vid0];
    for (vector<int>::iterator it = g.v[vindextemp].vneigh.begin(); it != g.v[vindextemp].vneigh.end(); ++it)
    {
        vecupdate.push_back(*it);
    }

    double e1 = 0;
    //energy of initial configuration
    //e1 += g.vertex_energy(vid0);
    e1 += g.compute_energy();

    e1 += g.find_dg(g.he[g.heidtoindex[xid]].type, g.he[g.heidtoindex[xidnextid]].type, g.he[g.heidtoindex[xidnextid]].din);
    //int breakprevxid=-1;
    //int breaknextprev=-1;
    if (gsl_rng_uniform(rg) < .5)
    { //BREAK prev- this
        // breakprevxid=g.he[g.heidtoindex[xidprevid]].vout;
        e1 += g.find_dg(g.he[g.heidtoindex[xidprevid]].type, g.he[g.heidtoindex[xid]].type, g.he[g.heidtoindex[xid]].din);
        g.he[g.heidtoindex[xidprevid]].nextid = -1;
        g.he[g.heidtoindex[xid]].previd = -1;
        //std::cout << "00 Fission Breaking bond between xidprevid " << xidprevid << " and xid " << xid << endl;
        //add a flag here
    }
    else
    {
        //breaknextprev=g.he[g.heidtoindex[xidnextid]].vout;
        //BREAK next-prev
        e1 += g.find_dg(g.he[g.heidtoindex[xidnextid]].type, g.he[g.heidtoindex[xidprevid]].type, g.he[g.heidtoindex[xidprevid]].din);
        g.he[g.heidtoindex[xidnextid]].nextid = -1;
        g.he[g.heidtoindex[xidprevid]].previd = -1;
        //std::cout << "00 Fission Breaking bond between xidnext " << xidnextid << " and xidprevid " << xidprevid << endl;
    }

    // BREAK xid-next
    //std::cout << "00 Fission Breaking bond between xid " << xid << " and xidnextid " << xidnextid << endl;
    g.he[g.heidtoindex[xid]].nextid = -1;
    g.he[g.heidtoindex[xidnextid]].previd = -1;

    int newvid1 = -1;
    int newvid2 = -1;

    double *newv1 = new double[3];
    double *newv2 = new double[3];
    double *tempv = new double[3]; // this is vector from this to new_moved

    g.move_p(g.v[vindex0].co, newv1, rg); // move one in one direction //

    //std::cout << "veclen(g.v[vindex0].co, tempv) " << veclen(g.v[vindex0].co, newv1) <<endl<<endl;

    subvec(g.v[vindex0].co, newv1, tempv); //find connection vector

    //std::cout << "norm(tempv) " << norm(tempv) <<endl<<endl;

    subvec(tempv, g.v[vindex0].co, newv2);

    //double dis1=veclen(g.v[vindex0].co,newv2);
    //std::cout<< " distance of g.v[vindex0].co,newv2 " << dis1 <<endl<<endl;

    //double dis=veclen(newv1,newv2);
    //std::cout<< " distance of two new vertex " << dis <<endl<<endl;

    int xidnextopindex = g.heidtoindex[g.he[g.heidtoindex[xidopid]].nextid];
    double len1 = veclen(g.he[xidnextopindex].hecent, newv1);

    double len2 = veclen(g.he[xidnextopindex].hecent, newv2);

    //std::cout << "len1 is "<< len1 << " len2 is " << len2 <<endl;

    /* store old vtx */
    VTX *vtxi;
    vtxi = new VTX;
    g.save_vtx(vid0, vtxi);

    //std::cout << "vtxi updated hein vtxi->hein.size()" <<vtxi->hein.size() <<endl;

    g.add_vertex(newv1);
    g.add_vertex(newv2);

    //std::cout << "in wedge fission Nv is " << g.Nv << " Nvlast is " << g.Nvlast <<endl;
    if (len1 <= len2)
    {
        newvid1 = g.Nvlast - 2; //xid->yid.vout
        newvid2 = g.Nvlast - 1; //xidnextid->zid.vin
    }
    else
    {
        newvid1 = g.Nvlast - 1;
        newvid2 = g.Nvlast - 2;
    }

    delete[] newv1;
    delete[] newv2;
    delete[] tempv;

    /* transfer edges to new vertices newvid1*/
    int yid = xid;
    int yidopid = g.he[g.heidtoindex[yid]].opid;
    while (true)
    {
        if (g.is_boundary(yidopid) > 0)
        { //std::cout << "yid  shoul be -1" << yid <<endl;
            break;
        };
        g.v[g.vidtoindex[newvid1]].hein.push_back(yid);
        yidopid = g.he[g.heidtoindex[yid]].opid;

        g.he[g.heidtoindex[yid]].vout = newvid1;
        g.he[g.heidtoindex[yidopid]].vin = newvid1;
        g.update_half_edge(yid);
        g.update_half_edge(yidopid);

        yid = g.he[g.heidtoindex[yidopid]].previd;
    }
    //std::cout << " newvid1 updated" <<endl;

    /* transfer edges to new vertices newvid2 */
    int zid = xidnextid;
    int zidopid = -1; //g.he[g.heidtoindex[yid]].opid;
    while (true)
    {

        zidopid = g.he[g.heidtoindex[zid]].opid;
        g.v[g.vidtoindex[newvid2]].hein.push_back(zidopid);

        g.he[g.heidtoindex[zid]].vin = newvid2;
        g.he[g.heidtoindex[zidopid]].vout = newvid2;
        g.update_half_edge(zid);
        g.update_half_edge(zidopid);

        zid = g.he[g.heidtoindex[zidopid]].nextid;
        if (zid == -1)
        { //std::cout << "zid  is -1" << zid <<endl;
            break;
        }
    }

    //std::cout << "00 now delete the vertex wedge fission" <<endl;

    g.delete_vertex(vid0);

    //g.update_boundary(); /test !!!!!!!!!!!!!!!!!!!!!!!!!!!!! 0425
    // calculate energy

    g.update_normals();

    //std::cout << "00 normals calculated wedge fission" <<endl;
    double e2 = 0;
    //energy of final

    //g.update_excluder_top_vertex(newvid2);
    //g.update_excluder_top_vertex(newvid1);
    //e2 += g.vertex_energy(newvid1);
    //e2 += g.vertex_energy(newvid2);
    e2 += g.compute_energy();

    //double vp = 1; //4 / 3. * 4 * M_PI * g.xi * g.xi * g.xi;
    double vp = g.l_thermal_kappa * g.l_thermal_kappa * g.l_thermal_kappa;
    double de = e2 - e1;
    double crit = exp(-de / g.T) * (vp * 2);
    if (g.Test_assembly == 1)
    {
        crit = 1;
    }

    int overlapflag = -1;

//    std::cout << "updating index in wedge fission (l 3547)" << std::endl;
    g.update_index();
    //std::cout << "00 indices updated wedge fission" <<endl;
    // update old neighbors
    for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
    {

        if (*it != vid0)
            g.update_neigh_vertex(*it);
    }

    //update new and neighbors
    g.update_neigh_vertex_and_neigh(g.Nvlast - 2);
    g.update_neigh_vertex_and_neigh(g.Nvlast - 1);

    if (g.check_overlap_g(g.Nvlast - 1) < 0)
    {
        overlapflag = 1;
    }
    else if (g.check_overlap_g(g.Nvlast - 2) < 0)
    {
        overlapflag = 1;
    }
    else
    {
        for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
        {

            if (g.check_overlap_g(*it) < 0)
            {
                overlapflag = 1;
            }
        }
    }

    if (gsl_rng_uniform(rg) < crit && overlapflag == -1)
    {
        std::cout << "wedge fission accepted added vid " << newvid1 << " " << newvid2 << endl;
        //std::cout << "#############" << endl;

        //update boundary index
        g.he[g.heidtoindex[xid]].boundary_index = bi; //this is hein of vid0
        g.he[g.heidtoindex[xidprevid]].boundary_index = bi;
        g.he[g.heidtoindex[xidnextid]].boundary_index = bi;

        //update nextid_boundary
        g.he[g.heidtoindex[nextidboundary0]].previd = -1;
        g.he[g.heidtoindex[previdboundary0]].nextid = -1;
        std::cout << "nextidboundary0 " << nextidboundary0 << " previdboundary0 " << previdboundary0 << endl;
        std::cout << " xid " << xid << " xidnextid " << xidnextid << " xidprevid " << xidprevid << endl;
        g.set_prev_next_boundary(xid, nextidboundary0);
        g.set_prev_next_boundary(previdboundary0, xidnextid);
        g.set_prev_next_boundary(xidprevid, xid);
        g.set_prev_next_boundary(xidnextid, xidprevid);

        delete vtxi;
        //g.check_odd_neigh();
        return 1;
    }
    else // if not accepted put things back
    {

        //vector<int> vecupdate;
        int vindextemp = g.vidtoindex[g.Nvlast - 1];
        if (g.v[vindextemp].vneigh.size() > 0)
        {
            for (vector<int>::iterator it = g.v[vindextemp].vneigh.begin(); it != g.v[vindextemp].vneigh.end(); ++it)
            {
                vecupdate.push_back(*it);
            }
        }
        vindextemp = g.vidtoindex[g.Nvlast - 2];
        if (g.v[vindextemp].vneigh.size() > 0)
        {
            for (vector<int>::iterator it = g.v[vindextemp].vneigh.begin(); it != g.v[vindextemp].vneigh.end(); ++it)
            {
                vecupdate.push_back(*it);
            }
        }
        //std::cout << "wedge fission not accepted " << endl; // now move things back to what it was

        g.add_vertex(vtxi->co);
        //update added neighbors

        //std::cout << "vtxi is added back -- updating indices g.Nvlast " << g.Nvlast  << " g.Nv "<< g.Nv <<endl;

        int lastindex = g.vidtoindex[g.Nvlast - 1];

        for (vector<int>::iterator ithe = vtxi->hein.begin(); ithe != vtxi->hein.end(); ++ithe)
        {
            g.v[lastindex].hein.push_back(*ithe);
            int heindex = g.heidtoindex[*ithe];
            g.he[heindex].vout = g.v[lastindex].vid;
            g.he[g.heidtoindex[g.he[heindex].opid]].vin = g.v[lastindex].vid;

            g.update_half_edge(g.he[heindex].id);
            g.update_half_edge(g.he[heindex].opid);
            //std::cout << " heid is " << g.he[heindex].id <<endl;
            //std::cout << " heid opid is " << g.he[heindex].opid <<endl;
        }

        g.set_prev_next(xid, xidprevid, xidnextid);
        g.set_prev_next(xidprevid, xidnextid, xid);
        g.set_prev_next(xidnextid, xid, xidprevid);

        int newvindex1 = g.vidtoindex[newvid1];
        int newvindex2 = g.vidtoindex[newvid2];
        //std::cout << " try deleting BEFORE!!!!! putting back the xtzi" <<endl;
        if (newvindex1 > newvindex2)
        {
            g.delete_vertex(g.v[newvindex1].vid);
            g.delete_vertex(g.v[newvindex2].vid);
        }
        else
        {
            g.delete_vertex(g.v[newvindex2].vid);
            g.delete_vertex(g.v[newvindex1].vid);
        }
        //std::cout << "00 no fission updating neighbor list"<<endl;
        //UPDATE INDICES

        //std::cout << "0000 fission" <<endl;
//        std::cout << "updating index in wedge fission (l 3669" << std::endl;
        g.update_index();
        for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
        {

            //std::cout <<" *it"<< *it <<endl;
            if ((*it != newvid1) && (*it != newvid2))
            {
                g.update_neigh_vertex(*it);
            }
        }

        //std::cout << "00 in wedge fission now update geometry" <<endl;

        g.update_geometry_vertex(g.vidtoindex[g.Nvlast - 1]);
        g.update_normals_vertex(g.vidtoindex[g.Nvlast - 1]);

        g.update_neigh_vertex_and_neigh(g.Nvlast - 1);

        /*
            g.check_odd_neigh();
            
            for (unsigned int i = 0; i < g.v.size(); i++)
            {    
                if (g.check_overlap_g(g.v[i].vid) < 0) { 
                    
                    dump_lammps_data_file(g, 5555555);
                    std::cout << "after move vertex" <<endl;
                    std::cout <<"OVERLAP vid " << g.v[i].vid << "vindex is " << g.vidtoindex[g.v[i].vid] << endl;
                    //g.find_overlap_g(g.v[ind].vid);
                    std::exit(-1); 
                } 
            } */
        vecupdate.clear();
        delete vtxi;
        return 0;
    }

    return -1;
}


int MC::attempt_change_edge_type_tri(System &g, int heid0)
{
    //std::cout << "inin attempt change type"<<endl;
    if (g.Nhe!=6) {
        std::cout <<"error not triangle in attempt_change_edge_type_tri"<<endl;
        std::exit(-1);
    }

    //we are dealing with inner triangle so look at the opid
    int heindex0 = g.heidtoindex[g.he[g.heidtoindex[heid0]].opid];
    int nextid0 = g.he[heindex0].nextid; // indexd of next of opposite edge
    int previd0 = g.he[heindex0].previd; // indexd of prev of opposite edge

    if (nextid0 == -1 || previd0==-1)
    {
        
        
        std::cout<< "error wrong triangle"<<endl;
        std::exit(-1);
    }
  

    //g.update_half_edge(heid0);
    //g.update_half_edge(g.he[heindex0].opid);

    //in this function heindex is the opposite of the heid 
    //inner triangle is set then oppsit hakf edges are set

    int newtype = -1;
    int next_newtype=-1;
    int prev_newtype=-1;

    //type selection:
    int newtypeselector = gsl_rng_uniform_int(rg, 6);

    switch (newtypeselector)
    {
    case 0: // 0 can be selected for drug now keep it to 3 
        {
        newtype=3;
        next_newtype=3;
        prev_newtype=3;
        }
        break;
    case 1: 
        {
        newtype=1;
        next_newtype=2;
        prev_newtype=0;
        }
        break;
    case 2: 
        {
        newtype=2;
        next_newtype=0;
        prev_newtype=1;
        }
        break;
    case 3: 
        {
        newtype=3;
        next_newtype=1;
        prev_newtype=2;
        }
        break;
    case 4: 
        {
        newtype=3;
        next_newtype=3;
        prev_newtype=3;
        }
        break;
    case 5: 
        {
        newtype=3;
        next_newtype=3;
        prev_newtype=3;
        }
        break;

    }

    //now save energy
    //3 binding energies and 3 dimer_bending 
    
    



    int nextindex0 = g.heidtoindex[nextid0];
    int nexttype = g.he[nextindex0].type;

    int previndex0 = g.heidtoindex[previd0];
    int prevtype = g.he[previndex0].type;

    int heidtype = g.he[heindex0].type;

    
    double e1 = 0;

    
    e1 += g.find_dg(heidtype, nexttype, g.he[nextindex0].din);
    e1 += g.dimer_bend_energy(heindex0);
    
    e1 += g.find_dg(prevtype, heidtype, g.he[heindex0].din);
    e1 += g.dimer_bend_energy(previndex0);
    
    e1 += g.find_dg(nexttype, prevtype,g.he[previndex0].din);
    e1 +=  g.dimer_bend_energy(nextindex0);

   
    //change type
    g.he[heindex0].type = newtype;
    g.he[nextindex0].type = next_newtype;
    g.he[previndex0].type = prev_newtype;

    double e2 = 0;
    e1 += g.find_dg(newtype, next_newtype, g.he[nextindex0].din);
    e1 += g.dimer_bend_energy(heindex0);
    
    e1 += g.find_dg(prev_newtype, newtype, g.he[heindex0].din);
    e1 += g.dimer_bend_energy(previndex0);
    
    e1 += g.find_dg(next_newtype, prev_newtype,g.he[previndex0].din);
    e1 +=  g.dimer_bend_energy(nextindex0);


    

    double de = e2 - (g.mu[newtype] + g.mu[next_newtype] +g.mu[prev_newtype])  - (e1 - g.mu[heidtype] +g.mu[nexttype] + g.mu[prevtype]) ;
    double crit = exp(-de / g.T);
    if (gsl_rng_uniform(rg) < crit)
    {
        //std::cout <<" change accepted";
        //update_opposite
        for (vector<HE>::iterator ithe = g.he.begin(); ithe != g.he.end(); ++ithe)
        {
            if (ithe->nextid!=-1 && ithe->previd!=-1)
                {
                    switch (ithe->type)
                    {
                    case 0:
                        g.he[g.heidtoindex[ithe->opid]].type=3;
                        break;
                    case 1:
                        g.he[g.heidtoindex[ithe->opid]].type=2;
                        break;
                    case 2:
                        g.he[g.heidtoindex[ithe->opid]].type=1;
                        break;
                    case 3:
                        g.he[g.heidtoindex[ithe->opid]].type=0;
                        break;
                    
                    default:
                        break;
                    }
                }
        }


        return newtype;
    }
    else
    {
        // back to old type
        g.he[heindex0].type = heidtype; 
        g.he[nextindex0].type = nexttype; 
        g.he[previndex0].type = prevtype;
        g.update_half_edge(heid0);
        g.update_half_edge(g.he[heindex0].opid);
        g.update_normals_vertex(g.vidtoindex[g.he[heindex0].vin]);
        g.update_normals_vertex(g.vidtoindex[g.he[nextindex0].vin]);
        //std::cout << " change_type not accepted"<<endl;
        return -1;
    }
}



int MC::attempt_change_edge_type(System &g, int heid0)
{
    //std::cout << "inin attempt change type"<<endl;
    int newtype = -1;
    int newoptype = -1;
    newtype = gsl_rng_uniform_int(rg, 4);

    int heindex0 = g.heidtoindex[heid0];
    int oldtype = g.he[heindex0].type;

    int heopindex0 = g.heidtoindex[g.he[heindex0].opid]; // indexd of opposite edge
    int optype = g.he[heopindex0].type;

    if (oldtype == 3 || oldtype == 0)
    {
        if (g.he[heindex0].din == 1 || g.he[heopindex0].din == 1)
        {
            return -1;
        }
        if (g.he[heindex0].nextid != -1)
        {
            if (g.he[g.heidtoindex[g.he[heindex0].nextid]].din == 1)
                return -1;
        }
        if (g.he[heopindex0].nextid != -1)
        {
            if (g.he[g.heidtoindex[g.he[heopindex0].nextid]].din == 1)
                return -1;
        }
    }

    if (oldtype == newtype)
        return -1;
    if (newtype == 0)
    {
        newoptype = 3;
    }
    else if (newtype == 3)
    {
        newoptype = 0;
    }
    else if (newtype == 2)
    {
        newoptype = 1;
    }
    else if (newtype == 1)
    {
        newoptype = 2;
    }

    //std::cout << "edge heid0  "  << heid0 <<"heindex0 "<< heindex0 << "opid " << g.he[heindex0].opid <<endl;
    double e1 = 0;

    g.update_half_edge(heid0);
    g.update_half_edge(g.he[heindex0].opid);

    int nextid0 = g.he[heindex0].nextid; // indexd of next  edge

    if (nextid0 != -1)
    {
        int nextindex0 = g.heidtoindex[nextid0];
        int nexttype = g.he[nextindex0].type;
        e1 += g.find_dg(oldtype, nexttype, g.he[nextindex0].din);
        e1 += g.bend_energy(nextindex0);
    }
    //std::cout << "one  " <<endl;
    int previd0 = g.he[heindex0].previd; // indexd of prev  edge
    if (previd0 != -1)
    {
        int previndex0 = g.heidtoindex[previd0];
        int prevtype = g.he[previndex0].type;
        e1 += g.find_dg(prevtype, oldtype, g.he[heindex0].din);
        e1 += g.bend_energy(previndex0) + g.dimer_bend_energy(previndex0);
    }

    //std::cout << "222  " <<endl;

    int nextopid0 = g.he[heopindex0].nextid; // id of next of opposite edge
    if (nextopid0 != -1)
    {
        int nextopindex0 = g.heidtoindex[nextopid0]; // id of prev of opposite edge
        int opnexttype = g.he[nextopindex0].type;
        e1 += g.find_dg(optype, opnexttype, g.he[nextopindex0].din);
        e1 += g.bend_energy(nextopindex0);
    }
    //std::cout << "333  " <<endl;
    int prevopid0 = g.he[heopindex0].previd;
    if (prevopid0 != -1)
    {
        int prevopindex0 = g.heidtoindex[prevopid0];
        int opprevtype = g.he[prevopindex0].type;
        e1 += g.find_dg(opprevtype, optype, g.he[heopindex0].din);
        e1 += g.bend_energy(prevopindex0) + g.dimer_bend_energy(prevopindex0);
    }

    e1 += g.bend_energy(heindex0);
    e1 += g.dimer_bend_energy(heindex0) + g.dimer_bend_energy(heopindex0);

    //std::cout << "333  " <<endl;

    /*if (oldtype == 1 || oldtype == 2)
        {
            if (gsl_rng_uniform(r) < 0.5)
            {
                newtype = 0;
                newoptype = 3;
            }
            else
            {
                newtype = 3;
                newoptype = 0;
            }
        }
        else if (oldtype == 3 || oldtype == 0)
        {
            if (gsl_rng_uniform(r) < 0.5)
            {
                newtype = 1;
                newoptype = 2;
            }
            else
            {
                newtype = 2;
                newoptype = 1;
            }
        }*/

    g.he[heindex0].type = newtype;
    g.he[heopindex0].type = newoptype;

    double e2 = g.bend_energy(heindex0);
    e2 += g.dimer_bend_energy(heindex0) + g.dimer_bend_energy(heopindex0);
    if (nextid0 != -1)
    {
        int nextindex0 = g.heidtoindex[nextid0];
        int nexttype = g.he[nextindex0].type;
        e2 += g.find_dg(newtype, nexttype, g.he[nextindex0].din);
        e2 += g.bend_energy(nextindex0);
    }
    //e2+=g.find_dg(newtype,nexttype,g.he[nextindex0].din);
    if (previd0 != -1)
    {
        int previndex0 = g.heidtoindex[previd0];
        int prevtype = g.he[previndex0].type;
        e2 += g.find_dg(prevtype, newtype, g.he[heindex0].din);
        e2 += g.bend_energy(previndex0) + g.dimer_bend_energy(previndex0);
    }
    //e2+=g.find_dg(prevtype,newtype,g.he[heindex0].din);
    if (nextopid0 != -1)
    {
        int nextopindex0 = g.heidtoindex[nextopid0]; // id of prev of opposite edge
        int opnexttype = g.he[nextopindex0].type;
        e2 += g.find_dg(newoptype, opnexttype, g.he[nextopindex0].din);
        e2 += g.bend_energy(nextopindex0);
    }

    if (prevopid0 != -1)
    {
        int prevopindex0 = g.heidtoindex[prevopid0];
        int opprevtype = g.he[prevopindex0].type;
        e2 += g.find_dg(opprevtype, newoptype, g.he[heopindex0].din);
        e2 += g.bend_energy(prevopindex0) + g.dimer_bend_energy(prevopindex0);
    }
    //e2+=g.find_dg(newoptype,opnexttype,g.he[nextopindex0].din);
    //e2+=g.find_dg(opprevtype,newoptype,g.he[heopindex0].din);

    double de = e2 - g.mu[newtype] - (e1 - g.mu[oldtype]);
    double crit = exp(-de / g.T);
    if (gsl_rng_uniform(rg) < crit)
    {
        //std::cout <<" change accepted";
        return newtype;
    }
    else
    {
        g.he[heindex0].type = oldtype; // back to old type
        g.he[heopindex0].type = optype;
        g.update_half_edge(heid0);
        g.update_half_edge(g.he[heindex0].opid);
        //std::cout << " change_type not accepted"<<endl;
        return -1;
    }
}

/***** This function : ***********************/
/* two halfedges will bind and a doubleboundary vertex will be created. */
/* boundary edges should be updated */
/*************************************/
int MC::attempt_fusion(System &g)
{

    //std::cout << "in attempt_fusion " << endl;
    //if (is_bond_vboundary(vid0)>0) return -1;

//    std::cout << "updating index in attemtp fusion (l 4085)" << std::endl;
    g.update_index();
    //std::cout << " Starting Fusion g.Nv " << g.Nv << endl;
    g.update_fusion_pairs_he();

    if (g.fusionhe.size() == 0)
    {
        return -1;
    }

    /*for (vector<int>::iterator it = g.fusionwedgehe.begin(); it != g.fusionwedgehe.end(); ++it)
        {
            int heindex00=g.heidtoindex[*it];
            std::cout<< " he is " << *it << " he_prev_wedge_fusion_heid " << g.he[heindex00].prev_wedge_fusion_heid << endl;
            std::cout<< " he is " << *it << " he_next_wedge_fusion_heid " << g.he[heindex00].next_wedge_fusion_heid << endl;
            std::cout<< " he is " << *it << " he_prev_fusion_heid " << g.he[heindex00].prev_fusion_heid << endl;
            std::cout<< " he is " << *it << " he_next_fusion_heid " << g.he[heindex00].next_fusion_heid << endl;
        }

    for (vector<int>::iterator it = g.fusionhe.begin(); it != g.fusionhe.end(); ++it)
        {
            int heindex00=g.heidtoindex[*it];
            if (g.he[heindex00].prev_wedge_fusion_heid!=-1) std::cout<< " In fusion he is " << *it << " he_prev_wedge_fusion_heid " << g.he[heindex00].prev_wedge_fusion_heid << endl;
            if (g.he[heindex00].next_wedge_fusion_heid!=-1) std::cout<< " In fusion he is " << *it << " he_next_wedge_fusion_heid " << g.he[heindex00].next_wedge_fusion_heid << endl;
            if (g.he[heindex00].prev_fusion_heid!=-1) std::cout<< " In fusion he is " << *it << " he_prev_fusion_heid " << g.he[heindex00].prev_fusion_heid << endl;
            if (g.he[heindex00].next_fusion_heid!=-1) std::cout<< " In fusion he is " << *it << " he_next_fusion_heid " << g.he[heindex00].next_fusion_heid << endl;
            //dump_restart_lammps_data_file(g,11111111);
            //std::cout << "exit for now"<<endl;
            //std::exit(-1);
        }*/

    int ind = gsl_rng_uniform_int(rg, g.fusionhe.size());
    int heid0 = g.fusionhe[ind];
    //std::cout << "heid for fusion " << heid0 <<endl;
    //if (g.Test_assembly==1) { heid0=}

    int heindex0 = g.heidtoindex[heid0];
    //if v[vindex0].vneigh.size()==0) return -1;
    int vidi = -1;
    int vidj = -1;
    int heid_next = -1;
    int heid_prev = -1;
    int heindex_prev = -1;
    int heindex_next = -1;
    if ((g.he[heindex0].prev_fusion_heid == -1) && (g.he[heindex0].next_fusion_heid == -1))
    {
        std::cout << "WRONG choice in fusions" << endl;
        std::exit(-1);
        return -1;
    }
    else if ((g.he[heindex0].prev_fusion_heid != -1) && (g.he[heindex0].next_fusion_heid != -1))
    { // if has both choose randomly
        if (gsl_rng_uniform(rg) < .5)
        { // with prev
            heindex_next = heindex0;
            heindex_prev = g.heidtoindex[g.he[heindex_next].prev_fusion_heid];
            vidi = g.he[heindex_next].vin;
            vidj = g.he[heindex_prev].vout;
        }
        else
        {
            heindex_prev = heindex0;
            heindex_next = g.heidtoindex[g.he[heindex_prev].next_fusion_heid];
            vidi = g.he[heindex_next].vin;
            vidj = g.he[heindex_prev].vout;
        }
    }
    else if (g.he[heindex0].prev_fusion_heid != -1)
    {
        heindex_next = heindex0;
        heindex_prev = g.heidtoindex[g.he[heindex_next].prev_fusion_heid];
        vidi = g.he[heindex_next].vin;
        vidj = g.he[heindex_prev].vout;
    }
    else if (g.he[heindex0].next_fusion_heid != -1)
    {
        heindex_prev = heindex0;
        heindex_next = g.heidtoindex[g.he[heindex_prev].next_fusion_heid];
        vidi = g.he[heindex_next].vin;
        vidj = g.he[heindex_prev].vout;
    }
    heid_next = g.he[heindex_next].id;
    heid_prev = g.he[heindex_prev].id;

    //double check if it is not wedge
    if (g.he[heindex_next].nextid_boundary == g.he[heindex_prev].previd_boundary)
    {
        std::cout << "error in fusion, this is a wedge!! " << endl;
        std::exit(-1);
    }
    // now save the status

    /* save the other side to update boundary index  heidbegin and heidend are begin and end of the loop for rupdating the bound_index*/
    int heidbegin = g.he[heindex_prev].nextid_boundary;
    int heidend = g.he[heindex_next].previd_boundary;

    if (heidbegin == -1 || heidend == -1)
    {
        std::cout << "error in fusion wrong heidbegin heidend" << endl;
        std::exit(-1);
    }
    //std::cout << "00 in fusion heidend " <<heidend<<endl;
    //std::cout << "00 in fusion heidbegin " <<heidbegin<<endl;

    double e1 = 0;
    e1 += g.compute_energy();
    //e1 += g.vertex_energy(vidi);
    //e1 += g.vertex_energy(vidj);
    //std::cout << " e1 is " << e1 << endl;

    //double e1test=0;

    //save vecupdate old neighbors except vidi , vidj
    int vindexi = g.vidtoindex[vidi];
    int vindexj = g.vidtoindex[vidj];

    vector<int> vecupdate;
    for (vector<int>::iterator it = g.v[vindexj].vneigh.begin(); it != g.v[vindexj].vneigh.end(); ++it)
    {
        if (*it != vidi)
            vecupdate.push_back(*it);
    }

    for (vector<int>::iterator it = g.v[vindexi].vneigh.begin(); it != g.v[vindexi].vneigh.end(); ++it)
    {
        if (*it != vidj)
            vecupdate.push_back(*it);
    }

    //save old vertex i
    VTX *vtxi;
    vtxi = new VTX;
    g.save_vtx(vidi, vtxi);

    //save old vertex j
    VTX *vtxj;
    vtxj = new VTX;
    g.save_vtx(vidj, vtxj);

    double *tempv = new double[3];
    double *newv = new double[3];
    centvec(g.v[vindexi].co, g.v[vindexj].co, newv); //tempv if random
    //no random position
    //g.move_p(tempv, newv, r);

    //add_new vertex  vid is newvid
    g.add_vertex(newv);
    int newvid = g.Nvlast - 1;
    //update index
//    std::cout << "updating index after add vertex (l 4234)" << std::endl;
    g.update_index(); // work ids after this
    //std::cout << " newvid is " << newvid << " vidtoindex[newvid] is " << g.vidtoindex[newvid] << " Nv is " << g.Nv << endl;
    delete[] tempv;
    delete[] newv;

    //std::cout << "new vertex added" << endl;

    //merging vertices
    //std::cout<< " After addition g.Nv " << g.Nv <<endl;
    int newvindex = g.vidtoindex[newvid];
    vindexi = g.vidtoindex[vidi];
    vindexj = g.vidtoindex[vidj];
    for (vector<int>::iterator ithe = g.v[vindexi].hein.begin(); ithe != g.v[vindexi].hein.end(); ++ithe)
    {

        g.v[newvindex].hein.push_back(*ithe);
        int heindex0 = g.heidtoindex[*ithe];
        g.he[heindex0].vout = newvid;
        g.he[g.heidtoindex[g.he[heindex0].opid]].vin = newvid;
        g.update_half_edge(g.he[heindex0].id);
        g.update_half_edge(g.he[heindex0].opid);
    }

    for (vector<int>::iterator ithe = g.v[vindexj].hein.begin(); ithe != g.v[vindexj].hein.end(); ++ithe)
    {

        g.v[newvindex].hein.push_back(*ithe);
        int heindex0 = g.heidtoindex[*ithe];
        g.he[heindex0].vout = newvid;
        g.he[g.heidtoindex[g.he[heindex0].opid]].vin = newvid;
        g.update_half_edge(g.he[heindex0].id);
        g.update_half_edge(g.he[heindex0].opid);
    }

    if (vindexi > vindexj)
    {
        g.delete_vertex(vidi);
        g.delete_vertex(vidj);
    }
    else
    {
        g.delete_vertex(vidj);
        g.delete_vertex(vidi);
    }

    //std::cout<< " After delete g.Nv " << g.Nv <<endl;
//    std::cout << "updating index after delete vertex (l 4281)" << std::endl;
    g.update_index();

    //can do only for those not connected
    g.he[g.heidtoindex[heid_next]].previd = heid_prev;
    g.he[g.heidtoindex[heid_prev]].nextid = heid_next;
    g.set_prev_next_boundary(heid_prev, heid_next);
    //g.set_prev_next(hei_, heidprev, heidnext);
    //g.set_prev_next(heidprev, heidnext, heidm);
    //g.set_prev_next(heidnext, heidm, heidprev);
    //std::cout << "in fusion next prev is set " << endl;

    g.update_boundary();
    g.update_normals();
    g.update_excluder_top_vertex(g.vidtoindex[newvid]);

    double e2 = g.find_dg(g.he[g.heidtoindex[heid_prev]].type, g.he[g.heidtoindex[heid_next]].type, g.he[heindex_next].din);

    //std::cout << " e2 bind is " << e2 << endl;
    //std::cout << " now new structure " << endl;
    //int overlap = -1;

    if (g.Test_assembly == 1)
    {
        //double e2test=g.compute_energy();
        //std::cout << "de test total is " << e2test-e1test<<endl;

        int vid0 = newvid;
        int vindex0 = g.vidtoindex[vid0];
        double v_eng = 0;
        for (vector<int>::iterator ithe = g.v[vindex0].hein.begin(); ithe != g.v[vindex0].hein.end(); ++ithe)
        {
            int heindex = g.heidtoindex[*ithe];
            std::cout << " hein is " << *ithe << " stretch_energy is " << g.stretch_energy(heindex) << " bend_energy is" << g.bend_energy(heindex) << " dimer_bend is " << g.dimer_bend_energy(heindex) << endl;
            v_eng += g.stretch_energy(heindex);
            v_eng += g.bend_energy(heindex);
            v_eng += g.dimer_bend_energy(heindex);

            if (g.he[heindex].previd != -1)
            {

                int previndex = g.heidtoindex[g.he[heindex].previd];
                v_eng += g.dimer_bend_energy(previndex);
                std::cout << "adding other side previd is " << g.he[heindex].previd << " dimer_bend is " << g.dimer_bend_energy(previndex) << endl;
                //if (is_boundary(he[heindex].previd)<0) {
                v_eng += g.bend_energy(previndex);
                //}
            }
            if (g.he[heindex].nextid != -1)
            {

                int nextindex = g.heidtoindex[g.he[heindex].nextid];
                std::cout << "adding other side nextid is" << g.he[heindex].nextid << " dimer_bend is " << g.dimer_bend_energy(nextindex) << endl;
                v_eng += g.dimer_bend_energy(nextindex);
            }
        }
        std::cout << "total is " << v_eng << endl;
    }
    e2 += g.compute_energy();

    //std::cout << " e2 is " << e2 << endl;
    //final configuration  // binding energies!

    //double vp = 1; //4 / 3. * 4 * M_PI * g.xi * g.xi * g.xi;
    double vp = g.l_thermal_kappa * g.l_thermal_kappa * g.l_thermal_kappa;
    double de = e2 - e1;
    //std::cout <<"de is " <<de <<endl;
    double crit = exp(-de / g.T) * 2 / (vp); //double check
    if (g.Test_assembly == 1)
    {
        std::cout << "crit is " << crit << endl;
        crit = 1;
    }
    //std::cout << "crit is "<<crit <<endl;
    //std::cout << " for now neglecting overlap " <<endl;
    int overlapflag = -1;
    //
    /* update neighbors of "neighbors of vtxi and vtxj"  */
    for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
    {
        g.update_neigh_vertex(*it);
    }

    g.update_neigh_vertex(newvid);
    int newvindex0 = g.vidtoindex[newvid];
    if (g.v[newvindex0].vneigh.size() > 0)
    {
        for (vector<int>::iterator it = g.v[newvindex0].vneigh.begin(); it != g.v[newvindex0].vneigh.end(); ++it)
        {
            g.update_neigh_vertex(*it);
            //vecupdate.push_back(*it);
        }
    }

    if (g.check_overlap_g(newvid) < 0)
    {
        overlapflag = 1;
    }
    else
    {
        for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
        {

            if (g.check_overlap_g(*it) < 0)
            {
                overlapflag = 1;
            }
        }
    }
    if (gsl_rng_uniform(rg) < crit && overlapflag == -1)
    {

        //std::cout << "00 in fusion fusion crit is met"<<endl;
        //std::cout << "00 in fusion g.Nboundary "<< g.Nboundary << " g.Nboundarylast " << g.Nboundarylast<< endl;
        //g.check_odd_neigh();
        g.v[newvindex0].doubleboundary = 1;

        /* update the boundary */ //already set
        //g.set_prev_next_boundary(heid_prev,heid_next);

        /* other side , change the boundary index */

        //std::cout << "00 in fusion other side , add Nboundary , change the boundary index "<<endl;

        g.set_prev_next_boundary(heidend, heidbegin);

        int checkendheid = heidbegin;
        int bi = g.Nboundarylast;
        while (checkendheid != heidend)
        {

            //std::cout << "00 in fusion updating boundary index for heid checkendheid= "<< checkendheid << endl;
            //std::cout << "00 before update boundary_index is "<< g.he[g.heidtoindex[checkendheid]].boundary_index <<endl;
            g.he[g.heidtoindex[checkendheid]].boundary_index = bi;
            //std::cout << "00 after update boundary_index is "<< g.he[g.heidtoindex[checkendheid]].boundary_index <<endl;
            checkendheid = g.he[g.heidtoindex[checkendheid]].nextid_boundary;
            if (checkendheid == -1)
            {
                std::cout << "00 error in fusion updating boundary index " << endl;
                std::exit(-1);
            }
        }
        g.he[g.heidtoindex[heidend]].boundary_index = bi;

        g.Nboundary++;
        g.Nboundarylast++;
        //std::cout << "00 in fusion g.Nboundary "<< g.Nboundary << " g.Nboundarylast " << g.Nboundarylast<< endl;
        //temp test
        checkendheid = heid_next;
        while (checkendheid != heid_prev)
        {

            //std::cout << "00 in fusion other side heid checkendheid= "<< checkendheid << endl;
            //std::cout << "00 other side g.he[g.heidtoindex[checkendheid]].boundary_index is "<< g.he[g.heidtoindex[checkendheid]].boundary_index <<endl;
            checkendheid = g.he[g.heidtoindex[checkendheid]].nextid_boundary;
            if (checkendheid == -1)
            {
                std::cout << "00 error in other side boundary index " << endl;
                std::exit(-1);
            }
        }
        // remove it after test

        vecupdate.clear();
        delete vtxi;
        delete vtxj;
        return 1;
    }
    else
    {
        //std::cout << "fusion not accepted  now move things back to what it was " << endl; //

        g.he[g.heidtoindex[heid_prev]].nextid = -1;
        g.he[g.heidtoindex[heid_next]].previd = -1;

        g.set_prev_next_boundary(heid_prev, heidbegin);
        g.set_prev_next_boundary(heidend, heid_next);

        //std::cout << "fusion not accepted deleting  newvid "<< newvid << endl;
        //vector<int> vecupdate;
        int vindextemp = g.vidtoindex[newvid];
        //std::cout << "newvid  vindex is "<< vindextemp << endl;
        //std::cout << " now delete its neighbors " << endl;
        if (g.v[vindextemp].vneigh.size() > 0)
        {
            for (vector<int>::iterator it = g.v[vindextemp].vneigh.begin(); it != g.v[vindextemp].vneigh.end(); ++it)
            {
                vecupdate.push_back(*it);
            }
        }

        //std::cout << "g.Nv " << g.Nv <<endl;
        //std::cout << "fusion not accepted adding vtxi back" << endl;
        g.add_vertex(vtxi->co);

        int lastindex = g.vidtoindex[g.Nvlast - 1];
        for (vector<int>::iterator ithe = vtxi->hein.begin(); ithe != vtxi->hein.end(); ++ithe)
        {

            int heindex = g.heidtoindex[*ithe];
            g.v[lastindex].hein.push_back(*ithe);
            g.he[heindex].vout = g.v[lastindex].vid;
            g.he[g.heidtoindex[g.he[heindex].opid]].vin = g.v[lastindex].vid;
            g.update_half_edge(g.he[heindex].id);
            g.update_half_edge(g.he[heindex].opid);
        }
        //std::cout << "g.Nv " << g.Nv <<endl;
        //std::cout << "fusion not accepted adding vtxj back" << endl;
        g.add_vertex(vtxj->co);

        lastindex = g.vidtoindex[g.Nvlast - 1];
        for (vector<int>::iterator ithe = vtxj->hein.begin(); ithe != vtxj->hein.end(); ++ithe)
        {
            int heindex = g.heidtoindex[*ithe];
            g.v[lastindex].hein.push_back(*ithe);
            ////std::cout << "updating vout of he " << g.he[heindex].id << "index " << heindex << " *ithe " << *ithe <<endl;
            g.he[heindex].vout = g.v[lastindex].vid;
            g.he[g.heidtoindex[g.he[heindex].opid]].vin = g.v[lastindex].vid;
            g.update_half_edge(g.he[heindex].id);
            g.update_half_edge(g.he[heindex].opid);
        }

        //std::cout <<" deleting ne vertex" <<endl;
        g.delete_vertex(newvid);

        //std::cout <<" deleted , updating index" <<endl;
        //g.update_index();
        //std::cout << "g.Nv " << g.Nv <<endl;
        // std::cout << "fusion not accepted updating index here !!!!!!!!!!!! "<< endl;
//        std::cout << "updating index after delete vertex (l 4510)" << std::endl;
        g.update_index();

        for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
        {

            //std::cout <<" *it"<< *it <<endl;
            if (*it != newvid)
            {
                g.update_neigh_vertex(*it);
            }
        }
        g.update_neigh_vertex(g.Nvlast - 1);
        int vindex0 = g.vidtoindex[g.Nvlast - 1];
        if (g.v[vindex0].vneigh.size() > 0)
        {
            for (vector<int>::iterator it = g.v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
            {
                g.update_neigh_vertex(*it);
            }
        }
        g.update_neigh_vertex(g.Nvlast - 2);
        vindex0 = g.vidtoindex[g.Nvlast - 2];
        if (g.v[vindex0].vneigh.size() > 0)
        {
            for (vector<int>::iterator it = g.v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
            {
                g.update_neigh_vertex(*it);
            }
        }
        //std::cout<<"attempt_vertex_fusion not accepted" <<endl;
        //g.check_odd_neigh();
        delete vtxi;
        delete vtxj;
        vecupdate.clear();
        return 0;
    }
    return -1;
}

/*****************Fission *********************/
/* Nboundary should be>1 for this to happen */

int MC::attempt_fission(System &g)
{
    //std::cout << "in attempt Fission  " << endl;
    if (g.Nboundary < 2)
        return -1;

    int ind = gsl_rng_uniform_int(rg, g.boundary.size()); //choose a boundary halfedge
    int heid0 = g.boundary[ind];
    int heindex0 = g.heidtoindex[heid0];
    int vid0 = -1;
    int heid_prev = -1;
    int heid_next = -1;

    //randomly choose vin or vout vertex for opening
    //choose previuos and next heid
    if (g.v[g.vidtoindex[g.he[heindex0].vin]].doubleboundary != -1)
    {
        vid0 = g.he[heindex0].vin;
        heid_prev = g.he[heindex0].previd;
        heid_next = heid0;
    }
    else if (g.v[g.vidtoindex[g.he[heindex0].vout]].doubleboundary != -1)
    {
        vid0 = g.he[heindex0].vout;
        heid_prev = heid0;
        heid_next = g.he[heindex0].nextid;
    }

    else
    {
        return -1;
    }

    if (heid_prev == -1 || heid_next == -1)
    {
        return -1;
    } // vertex should be bound

    if (g.he[g.heidtoindex[heid_prev]].boundary_index != g.he[g.heidtoindex[heid_next]].boundary_index)
    {
        std::cout << "error in fission wrong boundary index " << g.he[g.heidtoindex[heid_prev]].boundary_index << "and " << g.he[g.heidtoindex[heid_next]].boundary_index << endl;
        std::exit(-1);
    }

    // vertex should be double boundary
    // if  other side is also bound return -1; ?

    //save the vertex
    vector<int> vecupdate;
    int vindex0 = g.vidtoindex[vid0];
    for (vector<int>::iterator it = g.v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
    {
        vecupdate.push_back(*it);
    }

    double e1 = 0;
    //energy of initial configuration
    e1 += g.vertex_energy(vid0);
    e1 += g.find_dg(g.he[g.heidtoindex[heid_prev]].type, g.he[g.heidtoindex[heid_next]].type, g.he[g.heidtoindex[heid_next]].din);

    //ToDo make it through function
    //now do the unbind_fission
    //break the bond
    g.he[g.heidtoindex[heid_prev]].nextid = -1;
    g.he[g.heidtoindex[heid_next]].previd = -1;

    //Find new next previuos boundary
    /*int nextid_boundary = -1;
    int newprevid_boundary = -1;
    
    for (vector<int>::iterator it = g.v[vindex0].hein.begin(); it != g.v[vindex0].hein.end(); ++it)
    {
        if (*it != heid_prev && g.is_boundary(*it) > 0)
            newprevid_boundary = *it;
        if (*it != heid_next && g.is_boundary(g.he[g.heidtoindex[*it]].opid) > 0)
            newnextid_boundary = *it;
    }*/

    /*temp test */
    int nextidboundary0 = -1;
    int previdboundary0 = -1;
    if (g.v[g.vidtoindex[vid0]].heboundaryoutid == heid_next)
    {
        nextidboundary0 = g.v[g.vidtoindex[vid0]].heboundaryoutid2;
        previdboundary0 = g.he[g.heidtoindex[nextidboundary0]].previd_boundary;
    }
    else if (g.v[g.vidtoindex[vid0]].heboundaryoutid2 == heid_next)
    {
        nextidboundary0 = g.v[g.vidtoindex[vid0]].heboundaryoutid;
        previdboundary0 = g.he[g.heidtoindex[nextidboundary0]].previd_boundary;
    }
    else
    {
        std::cout << "error in fission heboundaryoutid" << endl;
    }

    if (nextidboundary0 == -1 || previdboundary0 == -1)
    {
        std::cout << "error in fission" << endl;
        std::exit(-1);
    }

    int newvid_prev = -1;
    int newvid_next = -1;

    double *newv1 = new double[3];
    double *newv2 = new double[3];
    double *tempv = new double[3];
    g.move_p( g.v[vindex0].co, tempv, rg); // move one in one direction //
    multvec(tempv, .5, newv1);

    // move the other one in opposit direction
    subvec(newv1, g.v[vindex0].co, tempv);
    addvec(g.v[vindex0].co, tempv, newv2);

    //int xidnextopindex = g.heidtoindex[g.he[g.heidtoindex[xidopid]].nextid];
    double len_prev_1 = veclen(g.he[g.heidtoindex[heid_prev]].hecent, newv1);
    double len_prev_2 = veclen(g.he[g.heidtoindex[heid_prev]].hecent, newv2);

    //std::cout << "len1 is "<< len1 << " len2 is " << len2 <<endl;

    /* store old vtx */
    VTX *vtxi;
    vtxi = new VTX;
    g.save_vtx(vid0, vtxi);

    //std::cout << "vtxi updated hein vtxi->hein.size()" <<vtxi->hein.size() <<endl;

    g.add_vertex(newv1);
    g.add_vertex(newv2);

    //std::cout << "in fission Nv is " << g.Nv << " Nvlast is " << g.Nvlast <<endl;
    if (len_prev_1 <= len_prev_2)
    {
        newvid_prev = g.Nvlast - 2; // first added vertex
        newvid_next = g.Nvlast - 1; // later added vertex
    }
    else
    {
        newvid_prev = g.Nvlast - 1; // later added vertex
        newvid_next = g.Nvlast - 2; // first added vertex
    }

    delete[] newv1;
    delete[] newv2;
    delete[] tempv;

    // now transfer edges of the heid_prev side
    std::cout << endl;
    //std::cout << "in fission now transfer edges of the heid_prev side" <<endl;
    int yid = heid_prev;
    int yidopid = g.he[g.heidtoindex[yid]].opid;
    if (g.is_boundary(yidopid) > 0)
    {
        std::cout << "error in fssion wrong geometry" << endl;
        std::exit(-1);
    } //temp test
    while (true)
    {
        //std::cout << "updating vout of yid " << yid << " to " << newvid_prev << endl;
        //std::cout << "in fission boundary index:  g.he[g.heidtoindex[yid]].boundaryindex " << g.he[g.heidtoindex[yid]].boundary_index<<endl;
        //std::cout << "in fission boundary index:  g.he[g.heidtoindex[yidopid]].boundaryindex " << g.he[g.heidtoindex[yidopid]].boundary_index<<endl;
        g.v[g.vidtoindex[newvid_prev]].hein.push_back(yid);
        //std::cout << "cleared vidtoindex" << std::endl;
        yidopid = g.he[g.heidtoindex[yid]].opid;
        //std::cout << "cleared heidtoindex" << std::endl;
        g.he[g.heidtoindex[yid]].vout = newvid_prev;
        g.he[g.heidtoindex[yidopid]].vin = newvid_prev;
        g.update_half_edge(yid);
        g.update_half_edge(yidopid);
        if (g.is_boundary(yidopid) > 0)
        { //std::cout << "yid  shoul be -1" << yid <<endl;
            break;
        }
        yid = g.he[g.heidtoindex[yidopid]].previd;
    }
    std::cout << endl;
    //now transfer edges of the heid_next side
    //std::cout << "in fission now transfer edges of the heid_next side" <<endl;

    int zid = heid_next;
    int zidopid = -1; //g.he[g.heidtoindex[yid]].opid;
    while (true)
    {

        zidopid = g.he[g.heidtoindex[zid]].opid;
        //std::cout << "updating vout of zidopid " << zidopid << " to " << newvid_next << endl;
        //std::cout << "in fission boundary index:  g.he[g.heidtoindex[zid]].boundaryindex " << g.he[g.heidtoindex[zid]].boundary_index<<endl;
        //std::cout << "in fission boundary index:  g.he[g.heidtoindex[zidopid]].boundaryindex " << g.he[g.heidtoindex[zidopid]].boundary_index<<endl;
        g.v[g.vidtoindex[newvid_next]].hein.push_back(zidopid);
        //std::cout << "cleared vidtoindex" << std::endl;

        g.he[g.heidtoindex[zid]].vin = newvid_next;
        //std::cout << "cleared heidtoindex" << std::endl;
        g.he[g.heidtoindex[zidopid]].vout = newvid_next;
        g.update_half_edge(zid);
        g.update_half_edge(zidopid);

        zid = g.he[g.heidtoindex[zidopid]].nextid;
        if (zid == -1)
        { //std::cout << "zid  is -1" << zid <<endl;
            break;
        }
    }

    std::cout << "in fission now delete the vertex" <<endl;

    g.delete_vertex(vid0);
    std::cout << "done" << std::endl;
    //std::cout << "in fission now update boundary" <<endl;
    //g.update_boundary(); should not update the boundary before
    // calculate energy

    std::cout << "in fission now update index " <<endl;
    ////std::cout << "in fission test vid outs "<<endl;
    /*
    for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
    {
    	std::cout << "heid " << it->id << " it->vout " << it->vout <<endl;

     }
    */

//    std::cout << "updating index in wedge fission (l 4774)" << std::endl;
    g.update_index();

    //std::cout << "done updating index" << std::endl;
    //std::cout << "in fission now update normals " <<endl;
    g.update_normals();

    double e2 = 0;
    //energy of final
    //std::cout << "in fission now calculate energy" <<endl;
    e2 += g.vertex_energy(newvid_next);
    //std::cout << "in fission now calculate energy" <<endl;
    e2 += g.vertex_energy(newvid_prev);

    //double vp = 1; //4 / 3. * 4 * M_PI * g.xi * g.xi * g.xi;
    double vp = g.l_thermal_kappa * g.l_thermal_kappa * g.l_thermal_kappa;
    double de = e2 - e1;
    double crit = exp(-de / g.T) * (vp * 2);
    //if (g.Test_assembly==1) { crit=1;}

    int overlapflag = -1;
    //std::cout << "in fission now update index" <<endl;
//    std::cout << "updating index in wedge fission (l 4796)" << std::endl;
    g.update_index();

    // update old neighbors
    //std::cout << "in fission now add neighbors" <<endl;

    /*if vecupdate.size()for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
    {
        
        if (*it != vid0)
            g.update_neigh_vertex(*it);
    }*/

    //update new and neighbors
    //std::cout << "in fission update new and neighbors" <<endl;
    g.update_neigh_vertex_and_neigh(newvid_next);
    g.update_neigh_vertex_and_neigh(newvid_prev);

    if (g.check_overlap_g(newvid_prev) < 0)
    {
        overlapflag = 1;
    }
    else if (g.check_overlap_g(newvid_next) < 0)
    {
        overlapflag = 1;
    }
    else
    {
        for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
        {

            if (g.check_overlap_g(*it) < 0)
            {
                overlapflag = 1;
            }
        }
    }

    if ((gsl_rng_uniform(rg) < crit && overlapflag == -1) || g.Test_assembly == 1)
    {
        //std::cout << "fission accepted added vid " << newvid_prev << " " << newvid_next << endl;
        std::cout << "fission accepted opeened heid_prev" << heid_prev << " and heid_next" << heid_next << endl;

        //std::cout << "#############" << endl;
        /* read the boundary index from one side */
        /* find other side begin */
        /* start a loop from that begin */
        /* until arrive at the new vertex */
        /* update the boundary index array */

        int nb = g.he[g.heidtoindex[heid_prev]].boundary_index;
        int heidupdate = nextidboundary0;

        while (heidupdate != previdboundary0)
        {
            g.he[g.heidtoindex[heidupdate]].boundary_index = nb;
            heidupdate = g.he[g.heidtoindex[heidupdate]].nextid_boundary;
        }
        g.he[g.heidtoindex[heidupdate]].boundary_index = nb;

        //set boundary next previous
        //std::cout << "in fission update boundary heid_prev" << heid_prev << " and , nextidboundary0 " << nextidboundary0<<endl;

        g.set_prev_next_boundary(heid_prev, nextidboundary0);
        //std::cout << "in fission update boundary previdboundary0" << previdboundary0 << " and heid_next " << heid_next <<endl;
        g.set_prev_next_boundary(previdboundary0, heid_next);
        g.Nboundary--;

        delete vtxi;
        //g.check_odd_neigh();
        return 1;
    }
    else
    {
        //std::cout << "in fission not acepted " << endl;
        //vector<int> vecupdate;
        int vindextemp = g.vidtoindex[g.Nvlast - 1];
        if (g.v[vindextemp].vneigh.size() > 0)
        {
            for (vector<int>::iterator it = g.v[vindextemp].vneigh.begin(); it != g.v[vindextemp].vneigh.end(); ++it)
            {
                vecupdate.push_back(*it);
            }
        }
        vindextemp = g.vidtoindex[g.Nvlast - 2];
        if (g.v[vindextemp].vneigh.size() > 0)
        {
            for (vector<int>::iterator it = g.v[vindextemp].vneigh.begin(); it != g.v[vindextemp].vneigh.end(); ++it)
            {
                vecupdate.push_back(*it);
            }
        }
        std::cout << "fission not accepted " << endl; // now move things back to what it was

        g.add_vertex(vtxi->co);
        //update added neighbors

        std::cout << "vtxi is added back -- updating indices g.Nvlast " << g.Nvlast  << " g.Nv "<< g.Nv <<endl;

        int lastindex = g.vidtoindex[g.Nvlast - 1];

        for (vector<int>::iterator ithe = vtxi->hein.begin(); ithe != vtxi->hein.end(); ++ithe)
        {
            g.v[lastindex].hein.push_back(*ithe);
            int heindex = g.heidtoindex[*ithe];
            g.he[heindex].vout = g.v[lastindex].vid;
            g.he[g.heidtoindex[g.he[heindex].opid]].vin = g.v[lastindex].vid;

            g.update_half_edge(g.he[heindex].id);
            g.update_half_edge(g.he[heindex].opid);
            //std::cout << " heid is " << g.he[heindex].id <<endl;
            //std::cout << " heid opid is " << g.he[heindex].opid <<endl;
        }

        g.he[g.heidtoindex[heid_prev]].nextid = heid_next;
        g.he[g.heidtoindex[heid_next]].previd = heid_prev;
        //g.set_prev_next(xid, xidprevid, xidnextid);
        //g.set_prev_next(xidprevid, xidnextid, xid);
        //g.set_prev_next(xidnextid, xid, xidprevid);

        int newvindex1 = g.vidtoindex[newvid_next];
        int newvindex2 = g.vidtoindex[newvid_prev];
        //std::cout << " try deleting BEFORE!!!!! putting back the xtzi" <<endl;
        if (newvindex1 > newvindex2)
        {
            g.delete_vertex(g.v[newvindex1].vid);
            g.delete_vertex(g.v[newvindex2].vid);
        }
        else
        {
            g.delete_vertex(g.v[newvindex2].vid);
            g.delete_vertex(g.v[newvindex1].vid);
        }
        ///std::cout <<"no fission updating neighbor list"<<endl;
        //UPDATE INDICES

        //std::cout << "0000 fission" <<endl;
//        std::cout << "updating index in fission rejection (l 4932)..." << std::endl;
        g.update_index();
        for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
        {

            //std::cout <<" *it"<< *it <<endl;
            if ((*it != newvid_next) && (*it != newvid_prev))
            {
                g.update_neigh_vertex(*it);
            }
        }

        //std::cout << "2222 fission" <<endl;

        g.update_neigh_vertex_and_neigh(g.Nvlast - 1);

        std::cout << "done updating index and neigh vertices" << std::endl;

        
        /*
            g.check_odd_neigh();
            
            for (unsigned int i = 0; i < g.v.size(); i++)
            {    
                if (g.check_overlap_g(g.v[i].vid) < 0) { 
                    
                    dump_lammps_data_file(g, 5555555);
                    std::cout << "after move vertex" <<endl;
                    std::cout <<"OVERLAP vid " << g.v[i].vid << "vindex is " << g.vidtoindex[g.v[i].vid] << endl;
                    //g.find_overlap_g(g.v[ind].vid);
                    std::exit(-1); 
                } 
            } */
        vecupdate.clear();
        delete vtxi;

        //Check that vouts are in vertex id list
        for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
        {
            //std::cout << "it->vout: " << it->vout << std::endl;
            int vout_in_v = 0;
            for (vector<VTX>::iterator it2 = g.v.begin(); it2 != g.v.end(); ++it2){
                //std::cout << "it2->vid: " << it2->vid << std::endl;
                if(it->vout==it2->vid) vout_in_v = 1;
            }
            if (vout_in_v==0){
                std::cout << "WARNING: in fission rejection: vout (" << it->vout << ") is not in the vextex id list!" << std::endl;
            }
        }


        return 0;
    }

    return -1;
}

/************************** unbind wedge **************************/
/******************************************************************/

int MC::attempt_unbind_wedge_dimer(System &g, int heid0)
{
    ////std::cout << "in attempt unbind_wedge "<< endl;

    int heindex0 = g.heidtoindex[heid0];   // this edge can bcome next or prev
    int heid_next = g.he[heindex0].nextid; // next bound on boundary
    int heid_prev = g.he[heindex0].previd; // prev bound on boundary
    //std::cout << "heid0 "<< heid0 << " heindex0 " << heindex0 << " heid_next_surf " << g.he[heindex0].nextid_boundary << " heid_prev_surf  " << g.he[heindex0].previd_boundary <<endl;
    //std::cout << "in unbind wedge heid0 "<< heid0 << " heindex0 " << heindex0 << " heid_next " << g.he[heindex0].nextid << " heid_prev " << g.he[heindex0].previd <<endl;

    if ((heid_next != -1) && (heid_prev != -1))
    {
        std::cout << heid0 << " heindex0 " << heindex0 << " heid_next_surf " << g.he[heindex0].nextid_boundary << " heid_prev_surf  " << g.he[heindex0].previd_boundary << endl;
        std::cout << " both next and previus on boundary !  wrong geometry!?" << endl;
        std::exit(-1);
    }

    int previndex = -1;
    //g.heidtoindex[heid_prev];
    int nextindex = -1;
    //g.heidtoindex[heid_next];

    if (heid_next != -1)
    {
        //std::cout << " unbind this and next" <<endl;
        nextindex = g.heidtoindex[heid_next];

        previndex = heindex0;
    }
    else if (heid_prev != -1)
    {
        //std::cout << " unbind this and prev" <<endl;
        previndex = g.heidtoindex[heid_prev];

        nextindex = heindex0;
    }
    else
    {
        std::cout << " no next no previus why in unbind!?" << endl;
        std::exit(-1);
    }
    int vid0 = g.he[nextindex].vin;
    /* if it is double boundary, other side should be bonded */
    if (g.v[g.vidtoindex[vid0]].doubleboundary == 1)
    {
        if (g.he[nextindex].id == g.v[g.vidtoindex[vid0]].heboundaryoutid)
        {
            if (g.he[g.heidtoindex[g.v[g.vidtoindex[vid0]].heboundaryoutid2]].previd == -1)
            {
                return (-1);
            }
        }
        else if (g.he[nextindex].id == g.v[g.vidtoindex[vid0]].heboundaryoutid2)
        {
            if (g.he[g.heidtoindex[g.v[g.vidtoindex[vid0]].heboundaryoutid]].previd == -1)
            {
                return (-1);
            }
        }
    }

    //std::cout << heid0 << heid0 << " nextindex " << nextindex << " nextid " << g.he[nextindex].id << " previndex " << previndex << " previd " << g.he[previndex].id <<endl;

    //now next and prev are obvoius
    double de = -g.find_dg(g.he[previndex].type, g.he[nextindex].type, g.he[nextindex].din);
    de -= g.dimer_bend_energy(previndex); //+g.bend_energy(nextindex)+g.bend_energy(previndex);
    double crit = exp((-de) / g.T);

    if (g.Test_assembly == 1)
    {
        std::cout << " unbind wedge de is " << de << " crit is " << crit << endl;
        crit = 1;
    }
    if (gsl_rng_uniform(rg) < crit)
    {
        g.he[previndex].nextid = -1;
        g.he[nextindex].previd = -1;
        //std::cout <<" unbind done " << g.he[previndex].id << " and " << g.he[nextindex].id <<endl;
        return 1;
    }
    //std::cout << " could not unbind crit not met" <<endl;
    return -1;
}

/* new test

int attempt_bind_wedge_dimer(System &g, int heid0, gsl_rng *r)
{
}*/

int MC::attempt_bind_wedge_dimer(System &g, int heid0)
{
    ////std::cout << "in attempt_bind_wedge_dimer heid0 " << heid0 << endl;
    int heindex0 = g.heidtoindex[heid0]; // this edge on boundary
    if ((g.he[heindex0].nextid != -1) || (g.he[heindex0].previd != -1))
    { // it should not be bound from any end
        std::cout << "why bond edge in attempt bind wedge" << endl;
        std::exit(-1);
    }

    int heidnext = -1;
    int heidprev = -1;
    int heindexnext = -1;
    int heindexprev = -1;

    //if (gsl_rng_uniform(r) <.5){
    //out << "gsl_rng_uniform(r) <.5" <<endl;
    heidnext = g.next_open_wedge(heid0); // also tests if next is not bound
    if (heidnext != -1)
    {
        if (heidnext != g.he[heindex0].nextid_boundary)
        {
            std::cout << "00 in bind wedge dimer error nextid_boundary is " << g.he[heindex0].nextid_boundary << " heidnext is " << heidnext << endl;
        }
    } //}
      //else {
    //std::cout << "gsl_rng_uniform(r) >.5" <<endl;
    heidprev = g.pre_open_wedge(heid0);
    //}

    if (heidprev != -1)
    {
        if (heidprev != g.he[heindex0].previd_boundary)
        {
            std::cout << "00 in bind wedge dimer error previd_boundary is " << g.he[heindex0].previd_boundary << " heidprev is " << heidprev << endl;
        }
    }

    if ((heidnext != -1) && (heidprev != -1))
    {
        if (gsl_rng_uniform(rg) < .5)
        {
            heindexnext = g.heidtoindex[heidnext];
            heindexprev = heindex0;
        }
        else
        {
            heindexprev = g.heidtoindex[heidprev];
            heindexnext = heindex0;
        }
    }
    else if (heidnext != -1)
    {
        heindexnext = g.heidtoindex[heidnext];
        heindexprev = heindex0;
    }
    else if (heidprev != -1)
    {
        heindexprev = g.heidtoindex[heidprev];
        heindexnext = heindex0;
    }
    else
    {

        return -1;
    }

    // double check ??
    //if (g.he[heindexprev].previd!=-1) { //std::cout<< " other end is bound in open wedge" <<endl;
    //if (g.no_bond_boundary(heidnext) < 0) { //std::cout<< " other end is bound in open wedge" <<endl;

    if (g.he[heindexnext].nextid != -1)
    {
        std::cout << " wrong geometry in bind-wedge next has next " << endl;
        std::exit(-1);
    } //
    if (g.he[heindexprev].previd != -1)
    {
        std::cout << " wrong geometry in bind-wedge prev has prev " << endl;
        std::exit(-1);
    } //
    //int heindexnext=g.heidtoindex[heidnext];
    //std::cout <<"next open wedge is  " << heidnext << endl;

    //if (g.he[heindexprev].nextid!=-1 || g.he[heindexnext].previd!=-1) { std::cout << " wrong geometry in bind-wedge " <<endl; std::exit(-1);}

    heidnext = g.he[heindexnext].id;
    heidprev = g.he[heindexprev].id;
    //std::cout << "heidprev " << heidprev <<  " heindexnext " << heindexnext << endl;
    g.he[heindexprev].nextid = g.he[heindexnext].id;
    g.he[heindexnext].previd = g.he[heindexprev].id;

    //std::cout << "heindnext"  << heidnext << " vin " << g.he[heindexnext].vin << endl;
    //std::cout << "heindprev"  << heidprev <<"vout " << g.he[heindexnext].vin << endl;
    if (g.he[heindexnext].vin != g.he[heindexnext].vin)
    {
        std::cout << " WHAT in bind wedge!" << endl;
        std::exit(-1);
    }

    //g.get_normal(heidnext);
    //g.get_normal(heidprev);

    //std::cout << "calculate de " <<endl;
    double de = g.dimer_bend_energy(heindexprev); //+ g.bend_energy(g.heidtoindex[heidnext])+g.bend_energy(g.heidtoindex[heidprev]);
    //std::cout << "1 de is "<< de<<endl;
    de += g.find_dg(g.he[heindexprev].type, g.he[heindexnext].type, g.he[heindexnext].din);
    //std::cout << "2 de is "<< de<<endl;
    double crit = exp((-de) / g.T);
    if (g.Test_assembly == 1)
    {
        crit = 1;
    }
    if (gsl_rng_uniform(rg) < crit)
    {
        //std::cout << "bound! heidpre to heid next"  << heidprev  << heidnext << endl;
        return (1);
    }
    else
    {
        //std::cout << "Did not bound! crit is " << crit <<endl;
        g.he[heindexprev].nextid = -1;
        g.he[heindexnext].previd = -1;
        //g.update_normals();
    }
    return (-1);
}

int MC::attempt_bind_triangle(System &g, int heid0) //
{
    //std::cout << "in attempt_bind_triangle heid0 " << heid0 << endl;

    int heindex0 = g.heidtoindex[heid0]; // this edge on boundary
    /* if triangle */
    int bi = g.he[heindex0].boundary_index;

    //std::cout << "in attempt_bind_triangle heindex0 " << heindex0 <<  " boundary_index " <<bi << endl;

    if (g.is_bond_in_boundary(heid0) > 0 || g.is_bond_out_boundary(heid0) > 0)
        return -1; // should be unbound
    int nextboundaryid0 = g.he[heindex0].nextid_boundary;
    int prevboundaryid0 = g.he[heindex0].previd_boundary;

    if (nextboundaryid0 == -1 || prevboundaryid0 == -1)
    {
        std::cout << "error in attempt_bind_triangle heindex0 " << endl;
        std::exit(-1);
    }

    int nextboundaryindex0 = g.heidtoindex[nextboundaryid0]; // next of heid0
    int prevboundaryindex0 = g.heidtoindex[prevboundaryid0]; // prev of heid0

    if (g.he[prevboundaryindex0].previd != g.he[heindex0].nextid_boundary)
        return -1; // should be infront of a wedge

    //bind the other two junctions of triangle
    //std::cout << "00 set_prev_next heid0 " << heid0 << " prev: g.he[heindex0].previd_boundary " <<g.he[heindex0].previd_boundary <<" next: g.he[heindex0].nextid_boundary "<<  g.he[heindex0].nextid_boundary <<endl;
    g.set_prev_next(heid0, prevboundaryid0, nextboundaryid0);
    g.set_prev_next(prevboundaryid0, nextboundaryid0, heid0);
    g.set_prev_next(nextboundaryid0, heid0, prevboundaryid0);

    std::cout << "00 id of g.he[heindex0].nextid_boundary " << g.he[heindex0].nextid_boundary << " its previd " << g.he[g.heidtoindex[g.he[heindex0].nextid_boundary]].previd << " its nextid " << g.he[g.heidtoindex[g.he[heindex0].nextid_boundary]].nextid << endl;
    //std::cout << "00 in attempt_bind_triangle boundary_index is " << bi << endl;

    g.update_boundary();

    g.update_normals_vertex(g.vidtoindex[g.he[heindex0].vout]);
    g.update_normals_vertex(g.vidtoindex[g.he[nextboundaryindex0].vout]);
    g.update_normals_vertex(g.vidtoindex[g.he[prevboundaryindex0].vout]);

    //std::cout << "00 in attempt_bind_triangle heid0 " << heid0 << endl;
    //std::cout << "00 in attempt_bind_triangle g.he[heindex0].previd_boundary " << g.he[heindex0].previd_boundary << endl;
    //std::cout << "00 in attempt_bind_triangle g.he[heindex0].previd " << g.he[heindex0].previd << endl;
    //std::cout << "00 in attempt_bind_triangle g.he[heindex0].nextid_boundary " << g.he[heindex0].nextid_boundary << endl;
    //std::cout << "00 in attempt_bind_triangle g.he[heindex0].nextid " << g.he[heindex0].nextid << endl;

    //std::cout << "calculate de " <<endl;
    double de = g.dimer_bend_energy(prevboundaryindex0) + g.dimer_bend_energy(heindex0);
    de += g.bend_energy(prevboundaryindex0) + g.bend_energy(heindex0) + g.bend_energy(nextboundaryindex0);
    //std::cout << "1 de is "<< de<<endl;
    de += g.find_dg(g.he[prevboundaryindex0].type, g.he[heindex0].type, g.he[heindex0].din);
    de += g.find_dg(g.he[heindex0].type, g.he[nextboundaryindex0].type, g.he[nextboundaryindex0].din);

    //std::cout << "2 de is "<< de<<endl;
    double crit = exp((-de) / g.T);
    if (g.Test_assembly == 1)
    {
        crit = 1;
    }
    if (gsl_rng_uniform(rg) < crit || g.Test_assembly == 1)
    {
        g.Nboundary--;
        std::cout << "bound! triangle " << heid0 << endl;
        return (1);
    }
    else
    {
        std::cout << "Did not bind triangle! crit is " << crit << endl;
        g.he[heindex0].boundary_index = bi;
        g.he[prevboundaryindex0].boundary_index = bi;
        g.he[nextboundaryindex0].boundary_index = bi;

        g.set_prev_next_boundary(g.he[prevboundaryindex0].id, g.he[heindex0].id);
        g.set_prev_next_boundary(g.he[heindex0].id, g.he[nextboundaryindex0].id);
        g.set_prev_next_boundary(g.he[nextboundaryindex0].id, g.he[prevboundaryindex0].id);

        g.he[heindex0].nextid = -1;
        g.he[heindex0].previd = -1;
        g.he[nextboundaryindex0].previd = -1;
        g.he[prevboundaryindex0].nextid = -1;
        return (-1);
    }
    return (-1);
}

int MC::attempt_unbind_triangle(System &g, int heid0)
{

    //std::cout << "in attempt_unbind_triangle heid0 " << heid0 << endl;
    int heindex0 = g.heidtoindex[heid0]; // this edge on boundary
    /* if triangle */

    if (g.is_boundary(heid0) > 0 || g.is_boundary(g.he[heindex0].opid) > 0)
        return -1; // should not be on surface
    int nextindex0 = g.heidtoindex[g.he[heindex0].nextid];
    int previndex0 = g.heidtoindex[g.he[heindex0].previd];

    if (g.is_boundary(g.he[nextindex0].opid) > 0 || g.is_boundary(g.he[previndex0].opid > 0))
        return -1; // sound not be the last triangle

    //if (g.he[previndex0].previd!=g.he[heindex0].nextid_boundary) return -1; // should be infront of a wedge

    //std::cout << "calculate de " <<endl;
    double de = -g.dimer_bend_energy(previndex0) + g.dimer_bend_energy(heindex0);
    de -= g.bend_energy(previndex0) + g.bend_energy(heindex0) + g.bend_energy(nextindex0);
    //std::cout << "1 de is "<< de<<endl;
    de -= g.find_dg(g.he[previndex0].type, g.he[heindex0].type, g.he[heindex0].din);
    de -= g.find_dg(g.he[heindex0].type, g.he[nextindex0].type, g.he[nextindex0].din);

    double crit = exp((-de) / g.T) / 2;
    if (gsl_rng_uniform(rg) < crit || g.Test_assembly == 1)
    {
        int bi = g.Nboundarylast;
        g.he[heindex0].boundary_index = bi;
        g.he[previndex0].boundary_index = bi;
        g.he[nextindex0].boundary_index = bi;

        g.set_prev_next_boundary(g.he[previndex0].id, g.he[heindex0].id);
        g.set_prev_next_boundary(g.he[heindex0].id, g.he[nextindex0].id);
        g.set_prev_next_boundary(g.he[nextindex0].id, g.he[previndex0].id);

        g.he[heindex0].nextid = -1;
        g.he[heindex0].previd = -1;
        g.he[nextindex0].previd = -1;
        g.he[previndex0].nextid = -1;

        g.Nboundary++;
        g.Nboundarylast++;
        std::cout << "ubound! triangle" << heid0 << endl;
        return (1);
    }
    else
    {
        std::cout << "Did not unbind triangle! crit is " << crit << endl;
        return (-1);
    }
    return (-1);
}

int MC::attempt_add_drug(System &g, int heid0)
{

    //std::cout << "in attempt adding drug " << heid0 <<endl;
    //std::cout << " g.Nd is " <<g.Nd<<endl;
    int heindex0 = g.heidtoindex[heid0];
    if (g.is_boundary(heid0) > 0)
    {
        return (-1);
    }
    if (g.he[heindex0].din == 1)
        return 0; //already has drug
    int hetype = g.he[heindex0].type;
    //int nexttype = -1;
    int previd0 = g.he[heindex0].previd;
    if (previd0 == -1)
    {
        return -1;
    } //
    int previndex0 = g.heidtoindex[previd0];

    //if (hetype == 1 || hetype == 2 || g.he[previndex0].type==1) //nexttype == 1 || nexttype == 2 ) //
    //{
    //    return -1;
    //}
    double e1 = 0; //(g.gdrug-g.mudrug);
                   // int previd0 = g.he[heindex0].previd;

    // if (previd0 != -1)
    // {

    e1 += g.dimer_bend_energy(previndex0);
    e1 += g.find_dg(g.he[previndex0].type, hetype, 0);
    //}

    g.he[heindex0].din = 1;

    double e2 = 0;
    //if (previd0 != -1)
    //{
    //int previndex0 = g.heidtoindex[previd0];
    e2 += g.dimer_bend_energy(previndex0);
    e2 += g.find_dg(g.he[previndex0].type, hetype, 1);
    //}
    //else
    //{
    //  e2 += (g.gdrug-g.mudrug);
    //}
    //std::cout << "Energy to add drug is " << e2 - e1 << " and should be -4.9" << endl;
    double crit = exp((-(e2 - e1)) / g.T);
    if (gsl_rng_uniform(rg) < crit)
    {
        //std::cout << "drug added on edge " << heindex0 <<endl;
        //std::cout << "prev is " << g.he[heindex0].previd <<endl;
        //std::cout << "prev is " << g.he[heindex0].previd <<endl;
        //std::cout<<"edgetype is " << hetype <<endl;

        //std::cout << " de is is " <<e2-e1<<endl;
        //std::cout << "crit is" << crit<<endl;
        g.update_half_edge(heid0);
        g.Nd++;
        //std::cout << " 010 g.Nd is " <<g.Nd<<endl;
        return 1;
    }
    else
    {
        g.he[heindex0].din = 0;
        return -1;
    }
}
int MC::attempt_remove_drug(System &g, int heid0)
{
    //std::cout << "in attempt remove drug " << heid0 <<endl;
    //std::cout << " g.Nd is " <<g.Nd<<endl;
    /* no drug on boundary*/
    int heindex0 = g.heidtoindex[heid0];
    int previd0 = g.he[heindex0].previd;
    if (previd0 == -1)
    {
        return -1;
    }
    //if (previd0==-1 && g.he[heindex0].din!=0) { std::cout <<"WRONG DRUG POSITION N REMOVE" <<endl; std::exit(-1);}

    if (g.he[heindex0].din == 0)
        return 0; //no drug
    int hetype = g.he[heindex0].type;
    //int nexttype=-1;
    //if (g.he[heindex0].nextid!=-1) { nexttype=g.he[g.heidtoindex[g.he[heindex0].nextid]].type; }
    //if (hetype == 1 || hetype == 2) //|| nexttype == 1 || nexttype == 2 ) //

    //{
    //    std::cout << " why drug on other than 0 3 g.he[heindex0].type" << g.he[heindex0].type << endl;
    //std::exit(-1);
    //}
    double e1 = 0; //(g.gdrug-g.mudrug);

    //if (previd0 != -1)
    //{
    int previndex0 = g.heidtoindex[previd0];
    e1 += g.dimer_bend_energy(previndex0);
    e1 += g.find_dg(g.he[previndex0].type, hetype, 1);
    //std::cout<< "g.dimer_bend_energy(previndex0); " << g.dimer_bend_energy(previndex0)<<endl;
    /*}
        else
        {
            e1 += (g.gdrug[][]-g.mudrug);
        }*/

    //std::cout << "g.he[heindex0].din "  <<  g.he[heindex0].din <<endl;
    g.he[heindex0].din = 0;
    //std::cout << "try removal g.he[heindex0].din "  <<  g.he[heindex0].din <<endl;
    double e2 = 0;
    //if (previd0 != -1)
    //{
    //int previndex0 = g.heidtoindex[previd0];

    e2 += g.dimer_bend_energy(previndex0);
    e2 += g.find_dg(g.he[previndex0].type, hetype, 0);
    //std::cout<< "g.dimer_bend_energy(previndex0); " << g.dimer_bend_energy(previndex0)<<endl;
    //}
    //std::cout << "Energy to remove drug is " << e2 - e1 << " and should be 4.9" << endl;
    double crit = exp((-(e2 - e1)) / g.T);
    //std::cout << " de removal drug is " <<e2-e1<<endl;
    //std::cout << " 011 g.Nd is " <<g.Nd<<endl;
    if (gsl_rng_uniform(rg) < crit)
    {
        //std::cout << "drug removed from edge " << heindex0 <<endl;
        // std::cout<<"edgetype is " << hetype <<endl;

        //std::cout << "crit is" << crit<<endl;
        //std::cout << "drug removed" <<endl;
        g.Nd--;
        //std::cout << " 012 g.Nd is " <<g.Nd<<endl;
        return 1;
    }
    else
    {
        g.he[heindex0].din = 1;
        //std::cout << " did not remove drug " <<g.Nd<<endl;
        //std::cout << " 013 g.Nd is " <<g.Nd<<endl;

        //std::cout <<"g.he[g.heidtoindex[heid0)].din"<<g.he[g.heidtoindex[heid0)].din<<endl;
        return 0;
    }
}
void MC::get_dimer_etypes(int etypeheid0, int etypenew1, int etypenew2)
{

    if (etypeheid0 == 0)
    {
        if (gsl_rng_uniform(rg) < .5)
        {
            etypenew1 = 1;
            etypenew2 = 2;
        } //ba
        else
        {
            etypenew1 = 0;
            etypenew2 = 0;
        }
    }

    else if (etypeheid0 == 3)
    {
        if (gsl_rng_uniform(rg) < .5)
        {
            etypenew1 = 3;
            etypenew2 = 3;
        } //ba
        else
        {
            etypenew1 = 1;
            etypenew2 = 2;
        }
    }

    else if (etypeheid0 == 1)
    {
        if (gsl_rng_uniform(rg) < .5)
        {
            etypenew1 = 2;
            etypenew2 = 0;
        } //ba
        else
        {
            etypenew1 = 2;
            etypenew2 = 3;
        }
    }

    else if (etypeheid0 == 2)
    {
        if (gsl_rng_uniform(rg) < .5)
        {
            etypenew1 = 0;
            etypenew2 = 1;
        }
        else
        {
            etypenew1 = 3;
            etypenew2 = 1;
        }
    }
    //std::cout << "adding dimer"<<endl;
    // Based on concentration
    /*    if (gsl_rng_uniform(r) < g.cdProb) {
                    if (gsl_rng_uniform(r) < 0.5) {etypenew1=0;}
                    else {etypenew1=3; }
                } 
                else{
                    if (gsl_rng_uniform(r) < 0.5) {etypenew1=1;}
                    else {etypenew1=2; } 

                } 

                if (gsl_rng_uniform(r) < g.cdProb) {
                    if (gsl_rng_uniform(r) < 0.5) {etypenew2=0;}
                    else {etypenew2=3; }
                } 
                else{
                    if (gsl_rng_uniform(r) < 0.5) {etypenew2=1;}
                    else {etypenew2=2; } 

                }*/
    //std::cout << "new types are etypenew1= " <<etypenew1 << " etypenew2 "  << etypenew2<<endl;
    //etypenew1=gsl_rng_uniform_int(r,4);
    //etypenew2=gsl_rng_uniform_int(r,4);
}

int MC::force_add_monomer_with_next(System &g, int heid0, int xid)
{

    //std::cout << "in force_add_monomer_with_next" <<endl;
    g.update_normals();

    if (g.is_boundary(heid0) < 0)
    {
        std::cout << " cannot add not on the !" << endl;
        std::exit(-1);
    }
    //double gbb=0;

    int heindex0 = g.heidtoindex[heid0];
    int etypeheid0 = g.he[heindex0].type;

    int bi = g.he[heindex0].boundary_index;
    //temp test
    if (bi == -1)
    {
        std::cout << "wrong boundary index in add monomer dimer" << endl;
        std::exit(-1);
    }

    //should not add if wedge fusion is possible
    g.update_fusion_pairs_he();

    if (g.he[heindex0].next_wedge_fusion_heid != -1 || g.he[heindex0].prev_wedge_fusion_heid != -1)
    {
        return -1;
    }
    //monomer with next

    //if (g.he[heindex0].nextid!=-1 && g.he[heindex0].previd!=-1) {std::cout <<" wrong not on boundary" <<endl; std::exit(-1); }
    //int xid = -1;

    int vid0 = g.is_bond_out_boundary(heid0);

    if (vid0 >= 0 && g.v[g.vidtoindex[g.he[heindex0].vin]].hein.size() < 6 && g.v[g.vidtoindex[g.he[g.heidtoindex[g.he[heindex0].nextid]].vin]].hein.size() < 6)
    {

        if (g.he[g.heidtoindex[g.he[heindex0].nextid]].next_wedge_fusion_heid != -1 || g.he[g.heidtoindex[g.he[heindex0].nextid]].prev_wedge_fusion_heid != -1)
        {
            //std::cout << "g.he[heindex0].next_wedge_fusion_heid "<<g.he[heindex0].next_wedge_fusion_heid<<endl;
            //std::cout << "g.he[g.heidtoindex[g.he[heindex0].previd]].next_wedge_fusion_heid " << g.he[g.heidtoindex[g.he[heindex0].previd]].next_wedge_fusion_heid<<endl;
            //std::cout << "g.he[g.heidtoindex[g.he[heindex0].nextid]].next_wedge_fusion_heid" <<g.he[g.heidtoindex[g.he[heindex0].nextid]].next_wedge_fusion_heid <<endl;
            //std::exit(-1);
            return -1;
        }

        //std::cout << "00 add monomer with next"<<endl;
        /* ToDo unify add with next and add with prev */
        //int thisprev_heid=heid0;
        //int thisnext_heid=g.he[heindex0].nextid;
        //xid = g.he[heindex0].nextid;
        int xidindex = g.heidtoindex[xid];

        if (g.is_bond_in_boundary(xid) != vid0 || g.is_bond_out_boundary(xid) != -1)
        {
            std::cout << "wrong bound on boundary" << endl;
            std::exit(-1);
        }
        else
        {
            //std::cout << "try add monomer" <<endl;
            int etypenew = -1;

            int etypexid = g.he[g.heidtoindex[xid]].type;
            // SELECTED TYPES

            /*if ((etypeheid0 == 2) && ((etypexid == 0) || (etypexid == 3)))
            {
                etypenew = 1;
            }
            else if (((etypeheid0 == 3) || (etypeheid0 == 0)) && (etypexid == 1))
            {
                etypenew = 2;
            }
            else if (((etypeheid0 == 3) && (etypexid == 3)))
            {
                etypenew = 3;
            }
            else if (((etypeheid0 == 0) && (etypexid == 0)) || ((etypeheid0 == 1) && (etypexid == 2)))
            {
                etypenew = 0;
            }
            else if ((etypeheid0 == 1) && (etypexid == 2))
            {
                if (gsl_rng_uniform(r) < 0.5)
                {
                    etypenew = 0;
                }
                else
                {
                    etypenew = 3;
                }
            }
            else
            {
                etypenew = gsl_rng_uniform_int(r, 4);
            }*/
		etypenew=gsl_rng_uniform_int(rg, 4);
            // TYPES BASED ON Concentration
            /*if (gsl_rng_uniform(r) < g.cdProb) {
                    if (gsl_rng_uniform(r) < 0.5) {etypenew=0;}
                    else {etypenew=3; }
                } 
                else{
                    if (gsl_rng_uniform(r) < 0.5) {etypenew=1;}
                    else {etypenew=2; } 

                }  */
            // etypenew=gsl_rng_uniform_int(r,4);

            bool drug0 = 0;
            //bool drug1=0;
            /*if ((etypenew==3) ||(etypenew==0) ) {
                if (gsl_rng_uniform(r) < g.drugProb) {drug0=1; }
                if (gsl_rng_uniform(r) < g.drugProb) {drug1=1;  }
                }*/
            //std::cout << "g.he[heindex0].din" << g.he[heindex0].din<<endl;
            double gbb = g.find_dg(etypenew, etypeheid0, g.he[heindex0].din);
            gbb += g.find_dg(etypexid, etypenew, drug0);
            //std::cout << "add_monomer with prev" <<endl;

            //double e1=g.bend_energy(heindex0) + g.bend_energy(xidindex);
            int x = g.add_monomer(heid0, xid, etypenew);
            //int Nhelast2index = g.heidtoindex[g.Nhelast - 2];
            if (x > 0)
            {
                //double de=g.monomer_energy(g.Nhelast-1);
                int voutid = g.he[heindex0].vout;
                int vinid = g.he[xidindex].vin;
                vector<int> vecupdate;
                vecupdate.push_back(vinid);
                vecupdate.push_back(voutid);
                g.update_half_edge(heid0);
                g.update_half_edge(g.he[heindex0].opid);
                g.update_half_edge(xid);
                g.update_half_edge(g.he[xidindex].opid);
                g.update_half_edge(g.Nhelast - 1);
                //g.he[g.heidtoindex(g.Nhelast-1)].din=drug1;
                g.update_half_edge(g.Nhelast - 2);
                //g.he[g.heidtoindex[g.Nhelast-2]].din=drug0;
                g.he[g.heidtoindex[g.Nhelast - 1]].boundary_index = bi; // update boundary index
                g.update_boundary();
                double de = g.stretch_energy(g.heidtoindex[g.Nhelast - 2]);

                de += g.dimer_bend_energy(g.heidtoindex[g.Nhelast - 2]) + g.dimer_bend_energy(xidindex);
                de += g.bend_energy(heindex0) + g.bend_energy(xidindex); //-e1; // g.dimer_bend_energy(heindex0);
                //de+=  ; //g.monomer_energy(heid0);
                //std::cout << " de is " << de <<endl;
                //double e2=g.compute_energy();
                //de += g.gb * 3 - g.mu;
                //gbb=gb0next+gb0prev;
                de += gbb - g.mu[etypenew];
                /*int vinid=g.he[g.heidtoindex[g.Nhelast - 1]].vin;
                    int vindex0=g.vidtoindex[vinid];
                    g.update_neigh_vertex(vinid);
                    if (g.v[vindex0].vneigh.size() > 0)
                    {
                        for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                        {	
                            g.update_neigh_vertex(*it);
                        }
                    }
                    int voutid=g.he[g.heidtoindex[g.Nhelast - 1]].vout;
                    vindex0=g.vidtoindex[voutid];
                    if (g.v[vindex0].vneigh.size() > 0)
                    {
                        for (vector<int>::iterator it =g. v[vindex0].vneigh.begin(); it != g.v[vindex0].vneigh.end(); ++it)
                        {	
                            g.update_neigh_vertex(*it);
                        }
                    }*/

                //double crit = 2.*g.z*g.K*g.K*g.K*exp((-de)/g.T);

                double crit = 1; //exp((-de) / g.T) / 2;
                int overlapflag = -1;

                /* check overlap */
                for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
                {
                    if (g.check_overlap_g(*it) < 0)
                        overlapflag = 1;
                }

                if (gsl_rng_uniform(rg) < crit && overlapflag == -1)
                {
                    //ToDo
                    // update next_prevoius surface

                    g.he[g.heidtoindex[g.Nhelast - 1]].boundary_index = bi;
                    g.v[g.vidtoindex[vid0]].doubleboundary = -1;
                    //if triangle otherside close it
                    if (g.he[g.heidtoindex[g.he[xidindex].nextid_boundary]].vout == g.he[g.heidtoindex[g.he[heindex0].previd_boundary]].vin)
                    {
                        g.set_prev_next(g.he[heindex0].previd_boundary, g.he[xidindex].nextid_boundary, g.Nhelast - 1);
                        g.set_prev_next(g.he[xidindex].nextid_boundary, g.Nhelast - 1, g.he[heindex0].previd_boundary);
                        g.set_prev_next(g.Nhelast - 1, g.he[heindex0].previd_boundary, g.he[xidindex].nextid_boundary);
                        g.Nboundary--;
                    }

                    //std::cout << "00 added monomer with next accepted"<<endl;
                    vecupdate.clear();
                    return 1;
                }
                /*deleting the added monomer */
                else
                {
                    /* ToDo update this with delete monomer */

                    int nextid0 = g.he[g.heidtoindex[g.Nhelast - 2]].nextid;
                    int previd0 = g.he[g.heidtoindex[g.Nhelast - 2]].previd;
                    int nextidboundary0 = g.he[g.heidtoindex[g.Nhelast - 1]].nextid_boundary;
                    int previdboundary0 = g.he[g.heidtoindex[g.Nhelast - 1]].previd_boundary;
                    //int xidindex=g.heidtoindex[xid];
                    if (nextid0 == -1 || previd0 == -1)
                    {
                        std::cout << "not accepted in add monomer!" << endl;
                        std::exit(-1);
                    }
                    if (g.delete_edge(g.Nhelast - 1) > 0)
                    {

                        g.he[heindex0].previd = -1;
                        g.he[xidindex].nextid = -1;

                        g.update_half_edge(heid0);
                        g.update_half_edge(g.he[heindex0].opid);
                        g.update_half_edge(xid);
                        g.update_half_edge(g.he[xidindex].opid);

                        // update boundary_index
                        g.he[heindex0].boundary_index = bi;
                        g.he[xidindex].boundary_index = bi;
                        g.set_prev_next_boundary(xid, nextidboundary0);
                        g.set_prev_next_boundary(previdboundary0, heid0);
                        g.set_prev_next_boundary(heid0, xid);
//                        std::cout << "updating index (l 5821)" << std::endl;
                        g.update_index();
                        //std::cout << " MONOMER REMOVED AFTER addition xid1" <<endl;

                        for (vector<int>::iterator it = vecupdate.begin(); it != vecupdate.end(); ++it)
                        {
                            g.update_neigh_vertex(*it);
                        }

                        // no need to update nex previous surface
                        //std::cout << "00 added monomer with next removed"<<endl;
                        vecupdate.clear();
                        return -1;
                    }
                    else
                    {
                        std::cout << " could not delete_monomer HERE 4444";
                        std::exit(-1);
                    }
                }
            }
            else
            {

                //std::cout << "could not add monomer" <<endl;
                return -1;
            }
        }
    }
    return 0;
}

int MC::check_vout_in_vid(System &g){
    for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
    {
        //std::cout << "it->vout: " << it->vout << std::endl;
        int vout_in_v = 0;
        for (vector<VTX>::iterator it2 = g.v.begin(); it2 != g.v.end(); ++it2){
            //std::cout << "it2->vid: " << it2->vid << std::endl;
            if(it->vout==it2->vid) vout_in_v = 1;
        }
        if (vout_in_v==0){
            std::cout << "WARNING: vout (" << it->vout << ") is not in the vextex id list!" << std::endl;
            //exit(-1);
        }
    }
}

int MC::check_vid_in_vout(System &g, int vid){
    int vid_in_vout = 0;
    for (vector<HE>::iterator it = g.he.begin(); it != g.he.end(); ++it)
    {
        if(it->vout == vid) vid_in_vout++;
    }
    if (vid_in_vout>0){
        return vid_in_vout;
    }
    else{
        return 0;
    }
}
