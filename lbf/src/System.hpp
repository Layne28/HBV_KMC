//A System consists of a set of Particles along with boundary conditions and control parameters (temperature, external field, etc.)

#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <vector>
#include <iostream>
#include <cmath>
#include <cstring>
#include "ParamDict.hpp"
#include "Observer.hpp"
#include "CustomRandom.hpp"
#include "tri_tri_intersect.hpp"

using namespace std;

class Observer;

struct HE
{
	int vin;  //beginning vertex id
	int vout; //end vertex id
	int id;	  // id of halfedge
	int type; //type of halfedge
	//type 0: CD
	//type 1: BA
	//type 2: AB
	//type 3: DC

	int nextid; // id of next halfedge
	int previd; // id of previous halfedge
	int opid;	//id of opposite halfedge

	int nextid_boundary; // id of next halfedge on boundary for edge on boundary , -1 if edge not on boundary
	int previd_boundary; // id of previous halfedge on boundary for edge on boundary , -1 if edge not on boundary

	double hevec[3];  // 3 component vector of this halfedge from vin to vout
	double hecent[3]; // xyz coordinates of the center of the halfedge (vin+vout)/2
	double n[3];	  // 3 component vector of the normal of this halfedge
	double hetop[3];  // xyz coordinates of the excluder calculated from normal and center of edge

	double l; // length of halfedge
	//bool dout;
	bool din; // drug boolian at the vin

	int prev_fusion_heid;		//previous fusion pair halfedge id
	int next_fusion_heid;		//next fusion pair halfedge id
	int prev_wedge_fusion_heid; //previous wedge_fusion pair halfedge id
	int next_wedge_fusion_heid; //next wedge_fusion pair halfedge id
	//double nextangle;
	int boundary_index; //id of boundary that this halfedge belong to,-1 if not on boundary
};

struct VTX
{
	long long int vid;		  // id of this vertex
	double co[3];	  // xyz coordinates of this vertex
	vector<int> hein; // vector of halfedges entering this vertex
	//int fusion_vid;
	int heboundaryoutid;
	int heboundaryoutid2;
	//int hesurfinid;
	int doubleboundary; //id of halfedge that is between two boundaries

	vector<int> vneigh; // vector of neighbor vertices
};
typedef struct HE HE;
typedef struct VTX VTX;

class System
{
private:
    Observer *obs; //use to write output
	gsl_rng *rg;

public:
    //Variables
	int seed = 0;
    bool Test_assembly; // for trial runs =1
	int Nvlast;
	int Nhelast;
	int Ntype;
	int Nv;
	int Nv5;
	int Nv6;
	int Nhe;
	int Nsurf;
	int Nd;
	int NAB;
	int NAB_in;
	int NCD;
	int NCD_T4;
	int NCD_T3;
	int NCD_T4_in;
	int NCD_T3_in;
	int NCD_Hex;
	int NCD_other;
	int Nv_in;
	int Nhe_in;
	int Nboundary;
	int Nboundarylast;
	int accepted_vmove;
	int rejected_vmove;

	int lenpoints;

	int all_neigh;
	double gb0;
	double dg;
	double mudimer;
	double *mu;
	double dmu;
	double mudrug;
	double ks0;
	double kd0;
	double dg12;
	double dg01;
	double dg20;
	double dg33;
	double dg00;
	double dgother;
	double gdrug0;

	double gaussian_sigma;
	double l_thermal_sigma;
	double l_thermal_kappa;
	double theta_thermal_kappa;

	double drugProb;
	double T;

	double xi;

	double *epsilon;
	double *kappa;
	double *phi0;
	double *kappaPhi;
	double *theta0;
	double *l0;
	double **gb;
	double **gdrug;
	
	double **dist_points;

	vector<VTX> v; /**< Current positions of the vertices. */
	vector<HE> he; /**< Pairs i,j of vertex indices that are connected. */
	//vector<VTX> tricent; // center of triangles
	
	int *vidtoindex;
	int *heidtoindex;
	vector<int> boundary; // all boundary he
	vector<int> boundaryv;
	vector<int> boundaryvbond;
	vector<int> fusionv;
	vector<int> fusionhe;
	vector<int> fusionwedgehe;
    
    //Variables for specific experiments
    std::string order_parameter; //order parameter for ffs, umbrella sampling, ...
                                 //TODO: add support for multiple order parameters


    /*** Methods ***/

    //constructor
    //TODO: add default ParamDict
    System(ParamDict &theParams, gsl_rng *&the_rg);

    //destructor
    ~System();

    //Make changes to System state
    void do_paramdict_assign(ParamDict &theParams);
    //void do_particle_init();
    //void do_spring_init();
    void set_obs(Observer &anObs);
    void initialize(int Ntyp0);
	void update_index();
	void update_neigh();
	void update_neigh_vertex(int vid0);
	void update_neigh_vertex_and_neigh(int vid0);
	void update_boundary();
	void update_excluder_top();
	void update_excluder_top_he(int heid);
	void update_normals_vertex(int vindex0);
	void update_geometry_vertex(int vindex0);
	void update_excluder_top_vertex(int vindex0);
	void update_geometry_parameters();

    //Miscellaneous
    int is_boundary(int heid0);
	int is_bond_in_boundary(int heid0);
	int is_bond_out_boundary(int heid0);
	int no_bond_boundary(int heid0);
	int is_vboundary(int vid);
	int is_bond_vboundary(int vid);
	int open_wedge(int heid0, int *flag);
	int pre_open_wedge(int heid0);
	int next_open_wedge(int heid0);
	int connected(int vid0, int vid1);
	int connectedH(int heid0, int heid1);
	int next_connected_boundary(int vid0, int vid1);
	int not_cross_edge(int heid0, int heid1);
	void add_vertex(double *xyz);
	int add_edge_type(int vin, int vout, int etype);
	void add_half_edge_type(int vin0, int vout0, int etype, int b_index);
	int check_overlap_centerh(double *tempcenter);
	int check_overlap_centerv(double *tempcenter);
	int do_intersect(int heid1, int heid2);
	int remove_neigh(int vid0, int removevid);
	int check_overlap_g(int vid0);
	int find_overlap_g(int vid0);
	int find_overlap_all();
	int check_inside_overlap(int heid0);
	void update_fusion_heid(int heidsurf0);
	int is_same_triangle(int heid0, int heid1);
	int opposite_edge(int heid0);
	void set_prev_next(int heid0, int previd0, int nextid);
	int get_prev_next(int heid0, int *previd0, int *nextid0);
	double add_dimer(int heid0, gsl_rng *r, int typenext, int typeprev, double *distance_vector);
	int force_add_dimer(int heid0, double *newv, int tyepnext, int typeprev);
	int remove_dimer(int heindex0, int heindexhext0);
	void make_hexamer();
	void he_initialize(int heindex, int heid0, int vin0, int vout0, int etype, int b_index);
	int add_monomer(int nextofnewid, int prevofnewid, int etype);
	double find_gbb(int etypenew, int etypenextofnew, int etypeprevifnew);
	double find_dg(int type, int typenext, bool drug);
	int force_add_monomer(int nextofnewid, int prevofnewid, int etype);
	int add_monomer_dimer(int heid0);
	int remove_monomer_dimer(int heid0, gsl_rng *r);
	void new_vertex(int heindex0, double *newv);
	void new_vertex_edge(int heindex0, double *newv,int et);
	double new_vertex_edge_and_move(int heindex0, double *newv, int etnew, gsl_rng *r);
	void new_vertex_points(int heindex0, int ind_point, double *newv);
	double move_p(double *pi, double *pf, gsl_rng *r);
	double move_p_gaussian(double len_v, double *pi, double *pf, gsl_rng *r);
	double move_p_gaussian_axes(int heindex0,double len_v, double *pi, double *pf, gsl_rng *r , double *dis_vector);
	double move_p_rotate_axes(int heindex0,double len_v, double *pi, double *pf, gsl_rng *r, double *dis_vector);
	double find_project_dist_axes(int heindex0, double *newv, double *oldv , double *dis_vector_project );
	void move_v_epsilon(double len, double *pi, double *pf, gsl_rng *r);
	void hecenter(int heindex, double *vcenter);
	void helen(int heindex, double *helen);
	int delete_vertex(int vid0);
	int delete_edge(int heid0);
	int get_normal(int heid0);
	void update_normals();
	void update_edge(int heid0);
	void update_half_edge(int heid0);
	double stretch_energy(int heindex0);
	int check_bind_wedge(int heid0);
	double bend_energy(int heindex0);
	double dimer_bend_energy(int heindex0);
	double dimer_energy(int heid0, int heid1);
	double monomer_energy(int heid0);
	double compute_bind_energy();
	double compute_energy();
	double vertex_energy(int vid0);
	void dump_parameters();
	void get_next_fusion_heid(int heidsurf0);
	void update_fusion_pairs_he();
	void save_vtx(int vid0, VTX *tempvtx);
	void check_odd_neigh();
	void set_prev_next_boundary(int previd0, int nextid0);

    //Get (some of these could be made static)
    Observer get_obs();

    //Methods for specific simulations
    double get_order_parameter();
};

/*** Miscellaneous functions. TODO: move these to other appropriate classes/.cpp files ***/
//e.g. all of the "dump" functions could be moved to Observer
void read_points(System &g);
void subvec(double *vinit, double *vfin, double *vec);
void addvec(double *vinit, double *vfin, double *vec);
void multvec(double *vinit, double scalar, double *vec);
void centvec(double *vinit, double *vfin, double *vec);
double veclen(double *vinit, double *vfin);
double norm(double *v);
double dot(double *v1, double *v2);
void cross(double *v1, double *v2, double *res);
void randvec(double *v, gsl_rng *r);
void dump_lammps_traj(System &g, int time0);
void dump_lammps_traj_dimers(System &g, int time0);
void dump_lammps_traj_restart(System &g, int time0);
void dump_lammps_data_file(System &g, int time0);
void dump_lammps_data_dimers(System &g, int time0);
void move_vertex(System &g, gsl_rng *r);
void rotatevec(double *vec, double *axis, double angle, double *vec2);
void read_lammps_data(System &g, char filename[]);
int read_restart_lammps_data_file(System &g, char filename[]);
int read_restart_lammps_data_traj(System &g, FILE *trajfile, int step);
void dump_restart_lammps_data_file(System &g, int time0);
void dump_data_frame(System &g, FILE *f, int time);
void recenter(System &g);
int surfclosev(System &g);
void make_initial_triangle(System &g);
void make_initial_diamond_T4(System &g);
void make_initial_diamond_CD(System &g);
void make_initial_pentamer(System &g);
int check_bind_triangle(System &g);
void dump_analysis(System &g, FILE *ofile, int sweep, int seed, int seconds);

#endif