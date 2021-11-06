#ifndef STATISTICS_H
#define STATISTICS_H
#include "Input_Class.h"

using namespace std;

class Statistics
{
  int n_paths;
  int NTraj;
  int NRec;

  double Temp;

  std :: string Stat_Dir;
  std :: string Traj_tot_fname;

  int nrings1;       // Number of rings for b1
  vector<double> b1_max;    // Array for b1_max
  double* Area_1;     // Array of rings area

  int nrings2;     // Number of rings for b2
  vector<double> b2_max;    // Array for b2_max
  double* Area_2;     // Array of rings area
  
  double** Traj;

  double omega_min;
  double omega_max;
  
  double** TrajsperBox;    // Number of trajs per (b1,b2) box. Third dimension is recombined or not
  double** TrajsRecperBox;
  double**** TrajsRec;     // Number of trajs that recombined into (b1,b2) box with (path_idx,istate)
  double**** Tau_Prob;     // Stores tau_OP*Prob(b1,b2,path_idx,istate). Required for rate constant
  double**** Tau_Prob2;    // Stores \sigma_{i=1}^{N_r} (tau_i^2)/N. Required for rate constant std. 
  
  double** Opacity_1;       // Opacity functions
  double** Opacity_2;       // Opacity functions
  
  double eps_LJ;
  double sig_LJ;

  int Poisson_treat;
  int resolve_path;

  vector<int> FinState;     // Vector containing final states
  vector<vector<int>> StateCount;   // Vecotr containing number of trajectories in each state

  int nstates;              // Total number of states
  
  int QNmax  = 1000;        // Max QN. This allows us to convert (v,j)-> i
  int istate;
  int vib,rot;
  double tau;
  
  int nlevels;
  double** levels;
  int missing;

  double** k_state;          // Recombination Rate Constant (State Dependent)
  double** k_state_var;      // Std Deviation in State Dependent Rate Constant
  double* krec;              // Overall Recombination Rate Constant
  double* krec_var;          // Std Deviation in Overall recombination rate constant
  
  // Private functions
  void SetRings();
  void Merge_Traj_files(Input_Class* Input);
  void Read_Traj_tot();
  void ReadLevels(Input_Class* Input);
  int  vjtoi(int v, int j);
  void itovj(int& v, int& j, int& i);

  void Compute_Probabilities();
  void Compute_RateConstant(Input_Class* Input);
  
 public:
   Statistics();
   ~Statistics();
  void Initialize_Statistics (Input_Class* Input);
  void Process_Statistics(Input_Class* Input);
  void WriteOutput(Input_Class* Input);
 
};

#endif /* STATISTICS_J */
