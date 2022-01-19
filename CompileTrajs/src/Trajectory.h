#ifndef TRAJECTORY_H
#define TRAJECTORY_H
//#include "Statistics.h"
#include "Input_Class.h"

class Trajectory
{
 public:
	int iNode;
	int iProc;
	int iTraj;
	int conv;
	int recomb_check;
	
	int iPES;           
	
	int parity = 0;    // 1:Odd j only, 2: Even j only; 0: both odd/even allowed
	
	std :: string Dir;        // Directroy
	std :: string PaQ_fname;  // File name of the PaQEvo file

	int NCoords;     // Number of coordinates of P/Q
	int NAtoms;
	
	double b1;           // Impact Parameter of A-B
	double b2;        // Impact Parameter of AB-C
	
	double b1_max;        // Maximum Impact Parameter of A-B
	double b2_max;     // Maximum Impact Parameter of AB-C
	
	double E1;       // Initial Relative translational energy of A-B
	double E2;    // Initial Relative translational energy of AB-C
	
	double E_int;     // Internal energy of the final molecule
	double E_kin_in;  // Total Kinetic energy of A + B + C
	double E_rel_fin; // Relative kinetic energy of O2+O
	
	double d1;         // Initial Distance between A and B 
	double d2;      // Initial Distance between C and COM of A,B
	
	double gamma;     // Gamma Paramter
	double alpha;     // alpha parameter
	double omega;     // omega parameter
	
	double tau;       // Bunker's Lifetime
	
	int nsteps;           // Number of timesteps in the trajectory    
	
	double v_ini;      // Initial Vibrational q.n
	double j_ini;      // Initial Rotational q.n
	double arr_ini;    // Initial arrangement q.n

	double v;          // Final Vibrational q.n
	double j;          // Final Rotational q.n
	double arr;        // Final Arrangement q.n
	
	double*  t;        // Time
	double** P;        // Velocities
	double** Q;        // Positions
	double** R;        // Internuclear distances
	double** Acc;      // Accelerations
	
	double* Pini;      // Initial Velocities
	double* Qini;      // Initial Positions

	double* Pfin;		// Final Velocities 
	double* Qfin;           // Final Positions
	double  H_fin;		// Final Hamiltonian	
		
	double mu1;      // Reduced mass of A and B 
	double mu2;     // Reduced mass of AB and C
	
	vector<vector<double>> r_min;   // Array storing the internuclear distance when each pair reaches its minimum distance for the first time
	vector<vector<double>> t_min;   // Array storing thhe time when each pair reaches its minimum distance for the first time
	vector<vector<int>> i_min;   // Array storing the index
	
	double r_direct;   // Radius for direct 3B Collision  (Defined in the constructor)
	double t_cmplx;    // Time difference for direct 3B collision (Defined in the constructor)
	double t_skip;     // Time skip parameter for direct 3B collision
	double t_del;
	int cmplx_id;      // Index of pair that forms a ""complex""
	
	bool assigned;     // A flag to keep track of wheter the pathway has already been assigned
	int path_idx;      // Index of recombination pathway; 0 : No recombination; 1: Lindemann; 2: Chaperon; 3: Direct
	
	//public:
	
	Trajectory();
	~Trajectory();

	void Initialize_Trajectory(int i, Input_Class* Input, double** Traj_tot, double** PaQSol, double** arr_matrix);
    
    //void Initialize_Trajectory(int iN, int iP, int iT, const std::string& Traj_dir, int NTraj_Global, double** Traj_tot, Statistics* Stat, vector<vector<double>> PaQSol);
   
	void Determine_pathway(std :: string Proc_Dir, Input_Class* Input);
	void Read_PaQEvo(std :: string fname, Input_Class* Input);
	void Calc_inter_distance();
	void Calc_acceleration(int method);
	int  Check_Direct();
	
    void Find_minima(int num, double t_beg);
    //int Check_Direct(Statistics* Stat);
    void Determine_cmplx();
    void Eval_tau(double e_LJ, double sigma_LJ);
    void Calc_Traj_Params();
    
    void Adjust_QN();
    
    // New functions to check if determining actual gamma works better in differentiating the mechanism
    void Calc_Jacobi();
    //int Check_Direct_New(Statistics* Stat);
    
};

#endif /* TRAJECTORY_H */
