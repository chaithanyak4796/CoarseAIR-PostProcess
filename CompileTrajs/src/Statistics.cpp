#include<iostream>
#include <cstdlib>
#include<fstream>
#include<math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include "Main.h"
#include "Input_Class.h"
#include "Global.h"
#include "Statistics.h"
#include "Logger.h"
#include "Trajectory.h"
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>

using namespace std;

Statistics :: Statistics() // Constructor
{
	//num_path = 1;
	
}

Statistics :: ~Statistics() // Destructor
{
	delete []count_path;
	
	for(int i=0; i<3; i++)
	{
		delete[] arr_matrix[i];
	}
	delete[] arr_matrix;
	
	Path_Traj.clear();
}

/**---------------------------------------------------------------------------------------------------------------------------------------**/

void Statistics :: Initialize_Statistics(int iP, Input_Class* Input)
{
	// This function initializes the variable members of this class.
	int i_Debug_Loc = 1;
	std :: string Debug = " [Initialize Statistics] :";
	if(i_Debug_Loc) Write(" ");
	if(i_Debug_Loc) Write(Debug,"Entering");
	if(i_Debug_Loc) Write(Debug,"Initializing Statistics for iProc = ",iP);
	
	num_path = 7;	
	// path = 0  ==> Lindemann
	// path = 1  ==> Chaperon
	// path = 2  ==> Direct
	// path = 3  ==> Other
	// path = 4  ==> No recombination
	// path = 5  ==> Unconverged Trajectories
	// path = 6  ==> System began with molecules
	
	if(i_Debug_Loc) Write(Debug,"Number of recombination pathways being considered, num_path =",num_path);
	
	Path_Traj.resize(num_path);
	for(int i=0; i<num_path; i++)
	{
		Path_Traj[i].resize(17);
	}
	
	if(i_Debug_Loc) Write(Debug,"The vector Path_traj has been resized with indices of :",num_path,"x 15"); 
	
	Traj_En.resize(7);
	if(i_Debug_Loc) Write(Debug," The Traj_En vector has been resized. ");

	
	count_path = new int [num_path];
	
	for(int i=0; i<num_path; i++)
	{
		count_path[i] = 0;
	}
	
	if(i_Debug_Loc) Write(Debug,"Initialized all elements of count_path to 0");
	
	arr_matrix = new double *[3];
	for (int i=0; i<3; i++)
	{
		arr_matrix[i] = new double [3];
	}
	
	// Initialize the arr_matrix
	arr_matrix[0][0] = 16.0;
	arr_matrix[0][1] = 17.0;
	arr_matrix[0][2] = 19.0;
	arr_matrix[1][0] = 32.0;
	arr_matrix[1][1] = 33.0;
	arr_matrix[1][2] = 35.0;
	arr_matrix[2][0] = 48.0;
	arr_matrix[2][1] = 49.0;
	arr_matrix[2][2] = 51.0;
	
	if(i_Debug_Loc) Write(Debug,"Initialized the arr_matrix with appropriate values");

	std :: string Traj_fname, PaQSol_fname;
	Proc_Dir = Input->Main_Dir + Input->Bins_Dir + "Node_1/Proc_" + to_string(iP) + "/";

	Traj_fname = Proc_Dir + "trajectories.out";
	if(i_Debug_Loc) Write(Debug, "Reading trajectories file : ", Traj_fname);
	Read_Traj(Traj_fname);
	if(i_Debug_Loc) Write(Debug, "Done reading trajectories file");
	if(i_Debug_Loc) Write(Debug, "Number of trajs in file = ", NTraj);
	if(i_Debug_Loc) Write(Debug, "Traj_tot[10][3] = ",Traj_tot[10][3]);

	PaQSol_fname = Proc_Dir + "PaQSol.out";
	if(i_Debug_Loc) Write(Debug, "Reading PaQSol file : ", PaQSol_fname);
	Read_PaQSol(PaQSol_fname);
	if(i_Debug_Loc) Write(Debug, "Done reading PaQSol file");
	if(i_Debug_Loc)Write(Debug, "PaQSol[10][3] = ",PaQSol[10][3]);

	
	if(i_Debug_Loc) Write(Debug,"Exiting");
}

/**--------------------------------------------------------------------------------------------------------------------------------------------------**/

double Statistics :: Get_arr_matrix(int i, int j)
{
	return arr_matrix[i][j];
}

void Statistics :: Update_count(int path)
{
	count_path[path] ++;
}

double Statistics :: Get_count_path(int path)
{
	return count_path[path];
}

/**-------------------------------------------------------------------------------------------------------------------------------------------------------**/
void Statistics :: Read_Traj(const std :: string fname)
{
 
  NTraj    = Find_length(fname,1);
  Traj_tot = new double* [NTraj];
  int n_cols = 14;
  for (int i=0; i<NTraj; i++)
    {
      Traj_tot[i] = new double [n_cols];
      for (int j=0; j<n_cols; j++)
	{
	  Traj_tot[i][j] = 0;
	}
    }

  ifstream f_traj;
  f_traj.open(fname.c_str());
  if(!f_traj)
    {
      Write("Error opening file : ",fname);
      exit(0);
    }

  int r_skip  = 1;
  std :: string line;
  for (int i=0; i<r_skip; i++)
    std :: getline(f_traj,line);

  for (int i=0; i<NTraj; i++)
    {
      for (int j=0; j<n_cols; j++)
	{
	  f_traj>>Traj_tot[i][j];
	}
    }
  
  f_traj.close();
}

/**-------------------------------------------------------------------------------------------------------------------------------------------------------**/
void Statistics :: Read_PaQSol(const std :: string fname)
{
  int len = Find_length(fname,1);
  if(len != NTraj)
    {
      Write(" WARNING :  The number of trajectories in PaQSol and trajectories.out do not match !!!!!!!");
    }
  PaQSol = new double* [NTraj];
  int n_cols = 28;
  for (int i=0; i<NTraj; i++)
    {
      PaQSol[i] = new double [n_cols];
      for (int j=0; j<n_cols; j++)
	{
	  PaQSol[i][j] = 0;
	}
    }

  ifstream f_sol;
  f_sol.open(fname.c_str());
  if(!f_sol)
    {
      Write("Error opening file : ",fname);
      exit(0);
    }

  int r_skip  = 1;
  std :: string line;
  for (int i=0; i<r_skip; i++)
    std :: getline(f_sol,line);

  int iT = -1;
  for (int i=0; i<len; i++)
    {
      iT++;
      for (int j=0; j<n_cols; j++)
	{
	  f_sol>>PaQSol[iT][j];
	}
      if( PaQSol[iT][0] != Traj_tot[iT][0] )
	{
	  iT--;
	}
    }
  
  f_sol.close();
}

/**-------------------------------------------------------------------------------------------------------------------------------------------------------**/
void Statistics :: Process_Trajs(int iP, Input_Class* Input, int i_Debug_Loc)
{
  
  std :: string Debug = "  [Process_Trajs] :";
  if(i_Debug_Loc) Write(Debug, "Entering");
  
  if(i_Debug_Loc) Write(Debug, "Initializing Output files");
  InitializeOutput_files(Input, iP);
  if(i_Debug_Loc) Write(Debug,"Done Initializing Output files");
  
  int iT, path_idx;
  
  for (int i=0; i<NTraj; i++)
    {
      if(i_Debug_Loc) Write(" ");
      iT = int(Traj_tot[i][0]);
      Trajectory* Traj = new Trajectory;
      Traj->Initialize_Trajectory(i, Input, Traj_tot, PaQSol, arr_matrix);

      // Determine reaction pathway
      if(Input->determine_pathway == 1 and Traj->recomb_check == 1)
	{
	  path_idx = Traj->Determine_pathway(Proc_Dir, Input, arr_matrix);
	  Update_count(path_idx-1);
	  UpdatePath_Traj( path_idx, Traj->iNode, Traj->iProc, Traj->E1*Hartree_eV, Traj->E2*Hartree_eV, Traj->b1, Traj->b2, Traj->b1_max, Traj->b2_max,
			  Traj->d1, Traj->d2, Traj->v, Traj->j, Traj->arr, Traj->omega, Traj->E_int*Hartree_eV, Traj->tau, Traj->ratio );

	  if(Input->write_misc == 1)  // Compute and write SAI, ratio of time scales
	    {
	      Traj->Evaluate_SAI();
	    }
	}
      WriteOutput(Traj, Input);
      delete Traj;
    }

  fclose(ftraj);
  if(Input->write_path_out && Input->determine_pathway)
    {
      fclose(fL);
      fclose(fC);
      fclose(fD);
      fclose(fN);
    }

  if (Input->determine_pathway)
    Write_Count(Input);

  if(Input->determine_pathway && Input->write_misc)
    fclose(f_misc);
  
  if(i_Debug_Loc) Write(Debug, "Exiting");
}
/**-------------------------------------------------------------------------------------------------------------------------------------------------------**/
void Statistics :: WriteOutput(Trajectory* Traj, Input_Class* Input)
{
  // This function writes the output of the statistics class
  fprintf(ftraj, ftraj_format, Traj->iPES, Traj->iProc, Traj->iTraj, Traj->E1*Hartree_eV, Traj->E2*Hartree_eV, Traj->b1_max, Traj->b1, Traj->b2_max, Traj->b2,
	  Traj->v, Traj->j, Traj->arr, Traj->omega, Traj->path_idx, Traj->E_int*Hartree_eV, Traj->tau);

  if(Input->write_path_out && Input->determine_pathway)
    {
      if(Traj->path_idx == 0)
	{
	  fprintf(fN, ftraj_format, Traj->iPES, Traj->iProc, Traj->iTraj, Traj->E1*Hartree_eV, Traj->E2*Hartree_eV, Traj->b1_max, Traj->b1, Traj->b2_max, Traj->b2,
		  Traj->v, Traj->j, Traj->arr, Traj->omega, Traj->path_idx, Traj->E_int*Hartree_eV, Traj->tau);
	}
      else if(Traj->path_idx == 1)
	{
	  fprintf(fL, ftraj_format, Traj->iPES, Traj->iProc, Traj->iTraj, Traj->E1*Hartree_eV, Traj->E2*Hartree_eV, Traj->b1_max, Traj->b1, Traj->b2_max, Traj->b2,
		  Traj->v, Traj->j, Traj->arr, Traj->omega, Traj->path_idx, Traj->E_int*Hartree_eV, Traj->tau);
	}
      else if(Traj->path_idx == 2)
	{
	  fprintf(fC, ftraj_format, Traj->iPES, Traj->iProc, Traj->iTraj, Traj->E1*Hartree_eV, Traj->E2*Hartree_eV, Traj->b1_max, Traj->b1, Traj->b2_max, Traj->b2,
		  Traj->v, Traj->j, Traj->arr, Traj->omega, Traj->path_idx, Traj->E_int*Hartree_eV, Traj->tau);
	}
      else if (Traj->path_idx == 3)
	{
	  fprintf(fD, ftraj_format, Traj->iPES, Traj->iProc, Traj->iTraj, Traj->E1*Hartree_eV, Traj->E2*Hartree_eV, Traj->b1_max, Traj->b1, Traj->b2_max, Traj->b2,
		  Traj->v, Traj->j, Traj->arr, Traj->omega, Traj->path_idx, Traj->E_int*Hartree_eV, Traj->tau);
	}
    }

  if(Input->determine_pathway && Input->write_misc && Traj->recomb_check)
    {
      fprintf(f_misc, fmisc_format, Traj->iProc, Traj->iTraj, Traj->SAI, Traj->ratio, Traj->cmplx_id, Traj->omega);
    }
}
/**-------------------------------------------------------------------------------------------------------------------------------------------------------**/
void Statistics :: InitializeOutput_files(Input_Class* Input, int iP)
{
  // Overall Traj file
  Output_Dir = Input->Main_Dir + Input->Bins_Dir + "Node_1/Proc_" + to_string(iP) + "/";
  
  fname_traj = Output_Dir + "Trajectories-3B.out";
  ftraj      = fopen(fname_traj.c_str(),"w");
  if(!ftraj)
    {
      Write("Error opening file :", fname_traj);
      exit(0);
    }

   ftraj_header = " iPES iProc  iTraj   E1[eV]    E2[eV]  b1_max[Bo]  b1[Bo]  b2_max[Bo]  b2[Bo]       v          j        arr        omega   Path_idx E_int[eV]     tau_OP[s]\n";
   ftraj_format = "%4d  %4d  %6d  %5.2E  %5.2E  %5.2E  %5.2E  %5.2E  %5.2E  %+5.2E  %+5.2E  %5.2E  %+6.4E  %2d    %+6.3E   %6.8E\n";

  fprintf(ftraj,ftraj_header);

  // Initializing the path specific trajectoru files
  if(Input->write_path_out && Input->determine_pathway)
    {
      fname_L = Output_Dir + "Lindemann_Traj.out";
      fname_C = Output_Dir + "Chaperon_Traj.out";
      fname_D = Output_Dir + "Direct_Traj.out";
      fname_N = Output_Dir + "No_Recomb_Traj.out";

      fL = fopen(fname_L.c_str(), "w");
      fC = fopen(fname_C.c_str(), "w");
      fD = fopen(fname_D.c_str(), "w");
      fN = fopen(fname_N.c_str(), "w");

      if(!fC || !fL || !fD || !fN)
	{
	  Write(" Error opening path specific files ");
	  exit(0);
	}

      fprintf(fL, ftraj_header);
      fprintf(fC, ftraj_header);
      fprintf(fD, ftraj_header);
      fprintf(fN, ftraj_header);
    }

  if(Input->determine_pathway && Input->write_misc)
    {
      fname_misc = Output_Dir + "Misc.out";
      f_misc     = fopen(fname_misc.c_str(), "w");

      if(!f_misc)
	{
	  Write(" Error opening file : ", fname_misc);
	  exit(0);
	}

      fmisc_header = " iProc  iTraj   SAI[rad]     ratio      cmplx_id    omega\n";
      fmisc_format = "%4d  %6d  %8.4E  %8.4E     %1d  %8.4E\n";

      fprintf(f_misc, fmisc_header);
    }

}
/**-------------------------------------------------------------------------------------------------------------------------------------------------------**/
void Statistics :: Write_Count(Input_Class* Input)
{
  // Count file
  std :: string fname_count = Output_Dir + "Count.out";
  FILE * f_count;
  f_count      = fopen(fname_count.c_str(),"w");
  if(!f_count)
    {
      Write("Error opening file :", fname_count);
      exit(0);
    }

  int f = 100;
  int tot = 0;
  for (int i=0; i<4; i++)
    tot+= count_path[i];

  // Write the count information
  fprintf(f_count,"L : %5d  %6.4f\n", count_path[0],double(f*count_path[0])/tot);
  fprintf(f_count,"C : %5d  %6.4f\n", count_path[1],double(f*count_path[1])/tot);
  fprintf(f_count,"D : %5d  %6.4f\n", count_path[2],double(f*count_path[2])/tot);
  fprintf(f_count,"N : %5d   ----\n", count_path[4]);

  fclose(f_count);
}

/**--------------------------------------------------------------------------------------------------------------------------------------------------------**/
void Statistics :: UpdatePath_Traj(int iPath, int iPES, int iP, int iT, double E1, double E2, double b1, double b2, double b1_max, double b2_max, double d1, double d2, double v, double j, double arr, double omega, double E_int, double T_OP)
{
	/*omp_lock_t writelock;
	omp_init_lock(&writelock);
	omp_set_lock(&writelock);*/
	{
	Path_Traj[iPath][0].push_back(iPES);
	Path_Traj[iPath][1].push_back(iP);
	Path_Traj[iPath][2].push_back(iT);
	Path_Traj[iPath][3].push_back(E1);
	Path_Traj[iPath][4].push_back(E2);
	Path_Traj[iPath][5].push_back(b1_max);
	Path_Traj[iPath][6].push_back(b1);
	Path_Traj[iPath][7].push_back(b2_max);
	Path_Traj[iPath][8].push_back(b2);	
	Path_Traj[iPath][9].push_back(d1);
	Path_Traj[iPath][10].push_back(d2);
	Path_Traj[iPath][11].push_back(v);
	Path_Traj[iPath][12].push_back(j);
	Path_Traj[iPath][13].push_back(arr);
	Path_Traj[iPath][14].push_back(omega);
	Path_Traj[iPath][15].push_back(E_int);
	Path_Traj[iPath][16].push_back(T_OP);
}
	//omp_unset_lock(&writelock);
	
}
/**--------------------------------------------------------------------------------------------------------------------------------------------------------**/
void Statistics :: UpdateTraj_En(int iPath, int iP, int iT, double E1, double E2, double E_kin_in, double E_int)
{
     Traj_En[0].push_back(iPath);
     Traj_En[1].push_back(iP);
     Traj_En[2].push_back(iT);
     Traj_En[3].push_back(E1);
     Traj_En[4].push_back(E2);
     Traj_En[5].push_back(E_kin_in);
     Traj_En[6].push_back(E_int);
}
/**---------------------------------------------------------------------------------------------------------------------------------------------------------**/
void Statistics :: MoveOutput(const std :: string& FromDir, const std :: string& ToDir )
{
	int i_Debug_Loc = 1;
	std :: string Debug = "  [MoveOutput]";
	
	if(i_Debug_Loc) Write();
	if(i_Debug_Loc) Write(Debug,"Entering");
	
	const std :: string Dir = ToDir + "Stat_files_2/";
	struct stat buffer;
	int Dir_exist = stat ((Dir).c_str(), &buffer);
	
	if(Dir_exist != 0)
	{
		if(i_Debug_Loc) Write(Debug,"Directory does not exist. Attempting to create.");
		int Dir_cr = mkdir((Dir).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		
		if(!Dir_cr)
		{
			if(i_Debug_Loc) Write(Debug,"Directory successfully created.");
		}
		else
		{
			if(i_Debug_Loc) Write(Debug,"Unable to create Directory :",Dir);
			exit(0);
		}
	}
	
	else
	{
		if(i_Debug_Loc) Write(Debug,"Directory Exists");
	}
	
	if(i_Debug_Loc) Write(Debug,"Moving Files now");
	int status [9] = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
	
	status[0] = rename((FromDir+"Count.out").c_str(), (Dir+"Count.out").c_str());
	status[1] = rename((FromDir+"Lindemann_Traj.out").c_str(), (Dir+"Lindemann_Traj.out").c_str());
	status[2] = rename((FromDir+"Chaperon_Traj.out").c_str(), (Dir+"Chaperon_Traj.out").c_str());
	status[3] = rename((FromDir+"Direct_Traj.out").c_str(), (Dir+"Direct_Traj.out").c_str());
	status[4] = rename((FromDir+"Others_Traj.out").c_str(), (Dir+"Others_Traj.out").c_str());
	status[5] = rename((FromDir+"Unconverged_Traj.out").c_str(), (Dir+"Unconverged_Traj.out").c_str());
	status[6] = rename((FromDir+"Bad_Initial_Cond_Traj.out").c_str(), (Dir+"Bad_Initial_Cond_Traj.out").c_str());
	status[7] = rename((FromDir+"No_Recomb_Traj.out").c_str(), (Dir+"No_Recomb_Traj.out").c_str());
	status[8] = rename((FromDir+"Recomb_Energy.out").c_str(), (Dir+"Recomb_Energy.out").c_str());
	//std :: filestream :: copy (FromDir.c_str(), Dir.c_str());
	
	for(int i=0; i<8; i++)
	{
		if(status[i] != 0)
		{
			if(i_Debug_Loc) Write(Debug,"Error in moving file #", i+1);
			if(i_Debug_Loc) Write(Debug,"Target Directory :", FromDir);
			if(i_Debug_Loc) Write(Debug,"Target Directory :", Dir);
			exit(0);
		}
	}
	
	if(i_Debug_Loc) Write(Debug,"All files moved successfully");
	
	if(i_Debug_Loc) Write(Debug,"Exiting");
}
