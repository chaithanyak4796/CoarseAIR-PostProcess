#ifndef STATISTICS_H
#define STATISTICS_H
#include <iostream>
#include <vector>
#include "Input_Class.h"
#include "Trajectory.h"

//using namespace std;

class Statistics
{
	int num_path ;     // Number of reaction pathways
	int* count_path;  // Count of the # of reactions in each pathway
	
	double** arr_matrix;  // Matrix for the arrangement q.n (Add STS later)

	int NTraj;            // Number of trajectories
	double** Traj_tot;    // Info of trajectories.out
	double** PaQSol;      // Info of PaQSol.out
	
 	vector<vector<vector<double> >> Path_Traj ;
	vector<vector<double> > Traj_En;

	std :: string Proc_Dir;

	FILE * ftraj;
	std :: string fname_traj;
	const char* ftraj_header;
	const char* ftraj_format;
	
	public:
	
	Statistics();
	~Statistics();
	
	void Initialize_Statistics(int iP, Input_Class* Input);
	void Read_Traj(const std :: string fname);
	void Read_PaQSol(const std :: string fname);

	void Process_Trajs(int iP, Input_Class* Input, int i_Debug_Loc);
	
	double Get_arr_matrix(int i, int j);
	double Get_count_path(int path);
	void Update_count(int path);
	void WriteOutput(Trajectory* Traj);
	void InitializeOutput_files(Input_Class* Input, int iP);
	
	//void WriteOutput_Old(const std :: string& Dir);
	void UpdatePath_Traj(int iPath, int iPES, int iP, int iT, double E1, double E2, double b1, double b2, double b1_max, double b2_max, double d1, double d2, double v, double j, double arr, double omega, double E_int, double T_OP);
	void MoveOutput(const std :: string& FromDir, const std :: string& ToDir );
	
	void UpdateTraj_En(int iPath, int iP, int iT, double E1, double E2, double E_kin_in, double E_int);
};

#endif /* STATISTICS_H */
