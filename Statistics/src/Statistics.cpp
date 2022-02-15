//--------- Performing all calculations in Hartree units ---------------//
//-------- Final values may be converted in the output section --------//
#include<iostream>
#include <cstdlib>
#include<fstream>
#include<math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include "Main.h"
#include "Global.h"
#include "Logger.h"
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "Statistics.h"
#include <algorithm>
#include "Input_Class.h"
#include <unistd.h>
#include "Reactions.h"

using namespace std;

Statistics :: Statistics() // Constructor
{
  NTraj   = 1;
  n_paths = 1;

  nrings1 = 1;
  nrings2 = 2;
  
}

Statistics :: ~Statistics() // Destructor
{
}


/**-------------------------------------------------------------------------------------------------------**/
void Statistics :: Initialize_Statistics (Input_Class* Input)
{
  // This function initializes the variable members of this class.
  int i_Debug_Loc = 1;
  std :: string Debug = " [Initialize Statistics] :";
  if(i_Debug_Loc) Write(Debug,"Entering");
  
  Temp      = Input->Temp;
  omega_min = Input->omega_min;
  omega_max = Input->omega_max;

  resolve_path = Input->resolve_path;
  n_paths      = 3;
  if(resolve_path == 0)
    n_paths = 1;
  Poisson_treat = Input->Poisson_treat;

  // Initialize the reactions object
  React = new Reactions;
  React->Initialize_Reactions(Input);
  num_reactions = React->num_reactions;

  // Merge the trajectory files
  Merge_Traj_files(Input);
  NTraj = Find_length(Traj_tot_fname,1);
  if(i_Debug_Loc) Write(Debug, "Total number of trajectories = ", NTraj);

  // Read the merged file
  if(i_Debug_Loc) Write(Debug, "Reading the merged trajectories file");
  Read_Traj_tot();
  if(i_Debug_Loc) Write(Debug, "Done reading the merged trajectories file");
  if(i_Debug_Loc) Write(Debug, "Traj[187][5] = ", Traj[187][5]);
  if(i_Debug_Loc) Write(Debug, "Number of recombined trajectories = ",NRec);

  // Setting the impact parameter rings
  SetRings();

  // Initialize the TrajsPerbox array
  TrajsperBox    = new double*  [nrings1];
  TrajsRecperBox = new double** [nrings1];
  
  for(int i=0; i<nrings1; i++)
    {
      TrajsperBox[i]    = new double [nrings2];
      TrajsRecperBox[i] = new double* [nrings2];
      
      for (int j=0; j<nrings2; j++)
	{
	  TrajsperBox[i][j]    = 0;
	  TrajsRecperBox[i][j] = new double[React->num_reactions];
	  for (int k=0; k<React->num_reactions; k++)
	    TrajsRecperBox[i][j][k] = 0;
	}
    }
  if(i_Debug_Loc) Write(Debug, "Done initializing Trajsper box");

  // Initialize the Opacity functions
  Opacity_1 = new double** [nrings1];   
  for (int i=0; i<nrings1; i++)
    {
      Opacity_1[i] = new double* [2];
      for (int j=0; j<2; j++)
	{
	  Opacity_1[i][j] = new double [num_reactions];
	  for (int k=0; k<num_reactions; k++)
	    Opacity_1[i][j][k] = 0;
	}
    }

  Opacity_2 = new double** [nrings2];
  for (int i=0; i<nrings2; i++)
    {
      Opacity_2[i] = new double* [2];
      for (int j=0; j<2; j++)
	{
	  Opacity_2[i][j]= new double [num_reactions];
	  for (int k=0; k<num_reactions; k++)
	    Opacity_2[i][j][k] = 0;
	}
    }
  if(i_Debug_Loc) Write(Debug, "Done initializing Opacity functions.");

  // Read the levels file
  ReadLevels(Input);
  
  if(i_Debug_Loc) Write(Debug,"Exiting");
}

/**-------------------------------------------------------------------------------------------------------**/
void Statistics :: Process_Statistics(Input_Class* Input)
{
  int i_Debug_Loc = 1;
  std :: string Debug = " [Process Statistics] :";

  if(i_Debug_Loc) Write();
  if(i_Debug_Loc) Write(Debug,"Entering");

  /* ======= Computing Probabilities and Cross Sections ================================*/
  Compute_Probabilities();

  Compute_RateConstant(Input);

  /* ======= Writing Outputs ================================*/
  WriteOutput(Input);

  if(i_Debug_Loc) Write(Debug,"Exiting");
}
/**-------------------------------------------------------------------------------------------------------**/
void Statistics :: Compute_Probabilities()
{
  /* This function:
   * Calculates the total number of trajectories in a box (b1,b2,omega)
   * Determines whether the final state already exists, updates the count if it exists and adds the final state if it does not exist
   * Then counts the trajectories and computes probabilities 
   */

  int i_Debug_Loc = 1;
  int i_Debug_St  = 0;  // For the states identification part

  std :: string Debug = "  [Compute_Probabilites] :";

  if(i_Debug_Loc) Write(Debug,"Entering");

  if(i_Debug_Loc) Write(Debug,"Identifying the final states");

  std :: vector<int> :: iterator itr;
  double b1,b2,om;
  int pid, s, count, reac_id;

  FinState.resize(num_reactions);

  for(int iTraj = 0; iTraj<NTraj; iTraj++)
    {
      //========================================= Identify all the final states =============================================================//
      if(i_Debug_St) Write(Debug, "iTraj = ", iTraj);

      vib = int(Traj[iTraj][9]);
      rot = int(Traj[iTraj][10]);
      istate = vjtoi(vib,rot);
      arr = Traj[iTraj][11];

      if(istate >= 0)    // Only bound states (non-recombined products have been assigned v=j=-1)
	{
	  reac_id = React->Identify_reaction(arr);
	  if(i_Debug_St) Write(Debug, " reac_id = ", reac_id);

	  if(FinState[reac_id].size() == 0)
	    {
	      FinState[reac_id].push_back(istate);
	      if(i_Debug_St) Write(Debug,"This state does not exist. Adding new state, istate =",istate);
	    }
	  else
	    {
	      // Identify if this state already exists in FinState
	      itr = std :: find(FinState[reac_id].begin(), FinState[reac_id].end() , istate);
	      if(itr != FinState[reac_id].end() )
		{
		  if(i_Debug_St) Write(Debug,"This state already exists. istate =",istate);
		}
	      else
		{
		  if(i_Debug_St) Write(Debug,"This state is new. Adding new state, istate =",istate);
		  FinState[reac_id].push_back(istate);
		}
	    }

	}
    }

  if(i_Debug_Loc) 
    {
      for(int i=0; i<num_reactions; i++)
	Write(Debug,"Total number of states found for reaction ", i, " = ", FinState[i].size());
    }

  //nstates = FinState.size();

  // Initializing the StateCount vector;
  StateCount.resize(num_reactions);
  for (int i=0; i<num_reactions; i++)
    {
      StateCount[i].resize(FinState[i].size());
      for(int j=0; j<FinState[i].size(); j++)
	{
	  StateCount[i][j].resize(5);
	  StateCount[i][j][0] = 0;    // Count
	  StateCount[i][j][1] = 0;    // Whether the state is B(0)/QB(1)/UQB(3)
	  StateCount[i][j][2] = 0;    // Number of Lind trajs
	  StateCount[i][j][3] = 0;    // Number of Chaperon trajs
	  StateCount[i][j][4] = 0;    // Number of Direct trajs
	}
      if(i_Debug_Loc) Write(Debug,"Initialized the StateCount vector for reaction ",i);
    }

  // Fully initialzing the TrajsRec, Tau_Prob and Tau_Prob2 arrays
  if(i_Debug_Loc) Write(Debug, "Initializing TrajsRec, Tau_Prob and Tau_Prob2 arrays.");

  TrajsRec.resize(num_reactions);
  Tau_Prob.resize(num_reactions);
  Tau_Prob2.resize(num_reactions);

  for(int i=0; i<num_reactions; i++)
    {
      TrajsRec[i].resize(nrings1);
      Tau_Prob[i].resize(nrings1);
      Tau_Prob2[i].resize(nrings1);

      for (int j=0; j<nrings1; j++)
	{
	  TrajsRec[i][j].resize(nrings2);
	  Tau_Prob[i][j].resize(nrings2);
	  Tau_Prob2[i][j].resize(nrings2);

	  for(int k=0; k<nrings2; k++)
	    {
	      TrajsRec[i][j][k].resize(n_paths);
	      Tau_Prob[i][j][k].resize(n_paths);
	      Tau_Prob2[i][j][k].resize(n_paths);

	      for(int p=0; p<n_paths; p++)
		{
		  TrajsRec[i][j][k][p].resize(FinState[i].size());
		  Tau_Prob[i][j][k][p].resize(FinState[i].size());
		  Tau_Prob2[i][j][k][p].resize(FinState[i].size());
		  
		  for(int l=0; l<FinState[i].size(); l++)
		    {
		      TrajsRec[i][j][k][p][l]  = 0;
		      Tau_Prob[i][j][k][p][l]  = 0.0;
		      Tau_Prob2[i][j][k][p][l] = 0.0;
		      }
		    }

		}
	 }

      if(i_Debug_Loc) Write(Debug, " Shape of TrajsRec  for reaction #", i, ": ", TrajsRec[i].size(), TrajsRec[i][0].size(), TrajsRec[i][0][0].size(), TrajsRec[i][0][0][0].size());
      if(i_Debug_Loc) Write(Debug, " Shape of Tau_Prob  for reaction #", i, ": ", Tau_Prob[i].size(), Tau_Prob[i][0].size(), Tau_Prob[i][0][0].size(), Tau_Prob[i][0][0][0].size());
      if(i_Debug_Loc) Write(Debug, " Shape of Tau_Prob2 for reaction #", i, ": ", Tau_Prob2[i].size(), Tau_Prob2[i][0].size(), Tau_Prob2[i][0][0].size(), Tau_Prob2[i][0][0][0].size());
    }

  if(i_Debug_Loc) Write(Debug, "Initialized TrajsRec, Tau_Prob and Tau_Prob2 arrays.");

  count = 0;
  int iden = 0; // Identified or not
  if(i_Debug_Loc) Write(Debug, "Collecting histogram info for", NTraj, " trajectories");
  int ntraj_count = 0;

  for(int iTraj=0; iTraj<NTraj; iTraj++)
    {
      iden = 0;

      //================================= Identify the boxes to which this trajectory belongs ===========================================================//
      b1 = Traj[iTraj][6];
      b2 = Traj[iTraj][8];
      om = Traj[iTraj][12];
      pid = int(Traj[iTraj][13]);

      vib = int(Traj[iTraj][9]);
      rot = int(Traj[iTraj][10]);
      arr = Traj[iTraj][11];

      tau = Traj[iTraj][15]/au_s;   // converting from s to atomic units

      istate = vjtoi(vib,rot);    // If you want to only resolve vibrationally, pass rot = 0

      for(int i1=0; i1<nrings1; i1++)
	{
	  
	  if( (i1<nrings1-1 && b1_max[i1]<=b1 && b1_max[i1+1]>b1) || (i1==nrings1-1 && b1_max[i1]<=b1 && b1_max[i1+1]>=b1) )
	    {
	      for(int j1=0; j1<num_reactions; j1++)
		Opacity_1[i1][0][j1]++;

	      for (int i2=0; i2<nrings2; i2++)
		{
		  if( (i2<nrings2-1 && b2_max[i2]<=b2 && b2_max[i2+1]>b2) || (i2==nrings2-1 && b2_max[i2]<=b2 && b2_max[i2+1]>=b2) )
		    {
		      for (int j2=0; j2<num_reactions; j2++)
			Opacity_2[i2][0][j2]++;

		      ntraj_count++;
		      TrajsperBox[i1][i2]++;
		      
		      if(istate >=0)
			{
			  reac_id = React->Identify_reaction(arr);
			  Opacity_1[i1][1][reac_id]++;
			  Opacity_2[i2][1][reac_id]++;
			  TrajsRecperBox[i1][i2][reac_id]++;
			  
			  iden = 1;
			  for(int i3=0; i3<n_paths; i3++)
			    {
			      if (i3+1 == pid) // pid = 1 (L); 2(C); 3(D)
				{
				  itr = find(FinState[reac_id].begin(),FinState[reac_id].end(),istate);

				  if(itr == FinState[reac_id].end())
				    {
				      if(i_Debug_Loc) Write(Debug,"Error: This state not identified earlier. Stopping. istate=",istate);
				      if(i_Debug_Loc) Write(Debug,"iProc = ",Traj[iTraj][1],"iTraj = ",Traj[iTraj][2]);
				      exit(0);
				    }
				  else
				    {
				      s = std :: distance(FinState[reac_id].begin(),itr);
				      TrajsRec[reac_id][i1][i2][i3][s]++;
				      if(Poisson_treat == 0)
					{
					  Tau_Prob[reac_id][i1][i2][i3][s]  += tau*fabs(om);
					  Tau_Prob2[reac_id][i1][i2][i3][s] += pow(tau*fabs(om),2);
					}
				      else
					{
					  Tau_Prob[reac_id][i1][i2][i3][s]  += pow(tau,2)*fabs(om);
					  Tau_Prob2[reac_id][i1][i2][i3][s] += pow(tau*tau*fabs(om),2);
					}

				      count++;
				      StateCount[reac_id][s][0]++;
				      StateCount[reac_id][s][1] = int(Traj[iTraj][11]) % 16; 
				      StateCount[reac_id][s][i3+2]++;
				    }
				  
				}
			    }
			}
		    }
		}
	    }
	}
      if(istate >=0 && iden ==0 )
	{
	  if(i_Debug_Loc) Write(Debug," This recombined traj has not been counted. iTraj =",iTraj,"b1 =",b1,"Loc_ID =",Traj[iTraj][1],",",Traj[iTraj][2], "istate = ",istate,"iden = ",iden);
	}
    }
  
  if(i_Debug_Loc) Write(Debug,"Done counting the trajectories. Calculating probabilities.");
  if(i_Debug_Loc) Write(Debug," ntraj_count - NTraj =",ntraj_count-NTraj);
  if(i_Debug_Loc) Write(Debug,"Nrec-count = ", NRec-count);
  if(NRec != count)
    {
      if(i_Debug_Loc) Write(Debug," WARNING !!!!!!!!!!  The number of recombination trajectories counted don't add up. Check limits on the variables.....");
      if(i_Debug_Loc) Write(Debug,".... Problem likely in limits of omega");
      //exit(0);
    }

  if(i_Debug_Loc)
    {
      double sum1, sum2;
      sum1=0; sum2=0;
      cout << Debug << " Opacity_1 : ";
      for(int i=0; i<nrings1; i++)
	{
	  sum1 += Opacity_1[i][0][0];
	  cout << Opacity_1[i][0][0]<< " ";
	}
      cout << "Sum = "<< sum1<<endl;
      cout << Debug << " Opacity_2 : ";
      for(int i=0; i<nrings2; i++)
	{
	  sum2 += Opacity_2[i][0][0];
	  cout << Opacity_2[i][0][0]<< " ";
	}
      cout << "Sum = "<< sum2<<endl;
      }


  // ==================================================== Calculate the probabilities ==================================================================//
  for (int r=0; r<num_reactions; r++)
    {
      for (int i1=0; i1<nrings1; i1++)
	{
	  for(int i2=0; i2<nrings2; i2++)
	    {
	      for(int i3=0; i3<n_paths; i3++)
		{
		  for(int i4=0; i4<FinState[r].size(); i4++)
		    {
		      if(TrajsperBox[i1][i2] > 0)
			{
			  TrajsRec[r][i1][i2][i3][i4]  /= TrajsperBox[i1][i2];
			  Tau_Prob[r][i1][i2][i3][i4]  /= TrajsperBox[i1][i2];
			  Tau_Prob2[r][i1][i2][i3][i4] /= TrajsperBox[i1][i2];
			}
		      else
			{
			  TrajsRec[r][i1][i2][i3][i4]  = 0;
			  Tau_Prob[r][i1][i2][i3][i4]  = 0;
			  Tau_Prob2[r][i1][i2][i3][i4] = 0;
			}
		    }
		}
	    }
	}
    }
  if(i_Debug_Loc) Write(Debug,"Done calculating Probabilities.");
  
  if(i_Debug_Loc) Write(Debug,"Exiting");
}
/**-------------------------------------------------------------------------------------------------------**/
void Statistics :: Compute_RateConstant(Input_Class* Input)
{
  // This function uses the probabilities and computes the recombination rate constants (State Dependent and Overall)
  int i_Debug_Loc = 1;
  std :: string Debug = "  [Compute_RateConstant] :";

  if(i_Debug_Loc) Write();
  if(i_Debug_Loc) Write(Debug,"Entering");

  if(i_Debug_Loc) Write(Debug,"Initializing the k_state and k_state_var arrays");

  if(i_Debug_Loc) Write(Debug, "num_reactions = ", num_reactions);


  k_state     = new double** [num_reactions];
  k_state_var = new double** [num_reactions];
  krec        = new double*  [num_reactions];
  krec_var    = new double*  [num_reactions];

  for (int r=0; r<num_reactions; r++)
    {
      k_state[r]     = new double* [n_paths+1];    // n_paths+1 because the last index is used to store the overall rates
      k_state_var[r] = new double* [n_paths+1];
      krec[r]        = new double  [n_paths+1];
      krec_var[r]    = new double  [n_paths+1];

      for (int i=0; i<n_paths+1; i++)
	{
	  krec[r][i] = 0;
	  krec_var[r][i] = 0;

	  k_state[r][i]     = new double [FinState[r].size()];
	  k_state_var[r][i] = new double [FinState[r].size()];

	  for (int j=0; j<FinState[r].size(); j++)
	    {
	      k_state[r][i][j]     = 0;
	      k_state_var[r][i][j] = 0;
	    }
	}
    }

  if(i_Debug_Loc) Write(Debug,"Done initializing the k_state and k_state_var arrays");

  if(i_Debug_Loc) Write(Debug,"Calculating the rate constants");

  for (int r = 0; r<num_reactions; r++)
    {
      for (int iring1=0; iring1<nrings1; iring1++)
	{
	  for(int iring2=0; iring2<nrings2; iring2++)
	    {
	      for(int s=0; s<FinState[r].size(); s++)
		{
	   
		  for(int p=0; p<n_paths; p++)
		    {
		      k_state[r][p][s] += Area_1[iring1]*Area_2[iring2]*Tau_Prob[r][iring1][iring2][p][s];
		      krec[r][p]       += Area_1[iring1]*Area_2[iring2]*Tau_Prob[r][iring1][iring2][p][s];

		      k_state_var[r][p][s] += pow(Area_1[iring1],2)*pow(Area_2[iring2],2)*(Tau_Prob2[r][iring1][iring2][p][s] - pow(Tau_Prob[r][iring1][iring2][p][s],2))/TrajsperBox[iring1][iring2];
		      krec_var[r][p]       += pow(Area_1[iring1],2)*pow(Area_2[iring2],2)*(Tau_Prob2[r][iring1][iring2][p][s] - pow(Tau_Prob[r][iring1][iring2][p][s],2))/TrajsperBox[iring1][iring2];
		    }
		
		}

	    }
	}
      
    }

  double m1, m2, m3, mu1, mu2, k_const;

  for (int r=0; r<num_reactions; r++)
    {
      krec[r][n_paths]     = 0;
      krec_var[r][n_paths] = 0;
  
      // Calculating the sum of rate constants across all paths
      for (int p=0; p<n_paths; p++)
	{
	  for(int s=0; s<FinState[r].size(); s++)
	    {
	      k_state[r][n_paths][s]     += k_state[r][p][s];
	      k_state_var[r][n_paths][s] += k_state_var[r][p][s];
	    }

	  krec[r][n_paths]     += krec[r][p];
	  krec_var[r][n_paths] += krec_var[r][p];
	}

      // Calculating the std deviation from the variance
      for (int p=0; p<n_paths+1; p++)
	{
	  for (int s=0; s<FinState[r].size(); s++)
	    {
	      k_state_var[r][p][s] = sqrt(k_state_var[r][p][s]);  // sqrt
	    }
	  krec_var[r][p] = sqrt(krec_var[r][p]); // sqrt
	}

      m1  = Input->Atom_Masses[0];
      m2  = Input->Atom_Masses[1];
      m3  = Input->Atom_Masses[2];
      mu1 = (m1*m2)/(m1+m2);
      mu2 = (m3*(m1+m2))/(m1+m2+m3);

      k_const = (React->symm_fac)*8*Kb_au*Temp*(omega_max-omega_min)/M_PI/sqrt(mu1*mu2);
      //double k_const = 4*Kb_au*Temp*(1)/M_PI/sqrt(mu1*mu2);

      for (int p=0; p<n_paths+1; p++)
	{
	  for (int s=0; s<FinState[r].size(); s++)
	    {
	      k_state[r][p][s]     *= k_const;
	      k_state_var[r][p][s] *= k_const;
	    }

	  krec[r][p]     *= k_const;
	  krec_var[r][p] *= k_const;
	}

      if(i_Debug_Loc) Write(Debug,"Done calculating the rate constants for reaction #",r);
      if(i_Debug_Loc) Write(Debug," krec [atomic units] =",krec[r][n_paths]);
    }

  if(i_Debug_Loc) Write(Debug,"Exiting");
  
}
/**-------------------------------------------------------------------------------------------------------**/
void Statistics :: WriteOutput(Input_Class* Input)
{
  int i_Debug_Loc = 1;
  std :: string Debug = " [Write Output] :";

  if(i_Debug_Loc) Write();
  if(i_Debug_Loc) Write(Debug,"Entering");

  for (int r=0; r<num_reactions; r++)
    {
      if(i_Debug_Loc) Write(Debug, "Writing the output for reaction #",r);

      // TrajsperBox file
      if(i_Debug_Loc) Write(Debug, "Writing the TrajsperBox file.");
      FILE * fbox;
      std :: string fname_box = Stat_Dir + "TrajsPerImpBox_" + to_string(r+1) + ".out";

      fbox = fopen(fname_box.c_str(),"w");
      if(!fbox)
	{
	  if(i_Debug_Loc) Write(Debug,"Error opening file:",fname_box);
	  exit(0);
	}

      fprintf(fbox,"    b1        b2        P(b1,b2)      N(b1)     N(b2)\n");
      for(int i=0; i<nrings1; i++)
	{
	  for(int j=0; j<nrings2; j++)
	    {
	      fprintf(fbox,"%6.4E  %6.4E  %6.4E  %5d  %5d\n",b1_max[i],b2_max[j],TrajsRecperBox[i][j][r]/TrajsperBox[i][j],int(TrajsRecperBox[i][j][r]),int(TrajsperBox[i][j]));
	    }
	}
      fclose(fbox);

      // Opacity functions
      if(i_Debug_Loc) Write(Debug,"Writing the opacity functions files");
      FILE * fimp1, *fimp2;
      std :: string fname_op1 = Stat_Dir + "Opacity_Func1_" + to_string(r+1) + ".out";
      std :: string fname_op2 = Stat_Dir + "Opacity_Func2_" + to_string(r+1) + ".out";

      fimp1 = fopen(fname_op1.c_str(),"w");
      fimp2 = fopen(fname_op2.c_str(),"w");

      if(!fimp1)
	{
	  if(i_Debug_Loc) Write(Debug,"Error opening file:",fname_op1);
	  exit(0);
	}

      if(!fimp2)
	{
	  if(i_Debug_Loc) Write(Debug,"Error opening file:",fname_op2);
	  exit(0);
	}

      fprintf(fimp1,"    b1        N_tot(b1)    N_rec(b1)    ratio\n");
      fprintf(fimp2,"    b2        N_tot(b2)    N_rec(b2)    ratio\n");


      for(int i=0; i<nrings1; i++)
	{
	  fprintf(fimp1,"%6.4E   %6.4E   %6.4E   %6.4E\n",(b1_max[i]+b1_max[i+1])/2,Opacity_1[i][0][r],Opacity_1[i][1][r],Opacity_1[i][1][r]/Opacity_1[i][0][r]);
	}
      for(int i=0; i<nrings2; i++)
	{
	  fprintf(fimp2,"%6.4E   %6.4E   %6.4E   %6.4E\n",(b2_max[i]+b2_max[i+1])/2,Opacity_2[i][0][r],Opacity_2[i][1][r],Opacity_2[i][1][r]/Opacity_2[i][0][r]);
	}
      fclose(fimp1);
      fclose(fimp2);
      if(i_Debug_Loc) Write(Debug,"Done writing the opacity functions files.");
  
      // FinStates.out
      FILE * fst;

      std :: string fname_state = Stat_Dir + "FinStates_" + to_string(r+1) + ".out";
      fst     = fopen(fname_state.c_str(),"w");

      std :: vector<int> :: iterator itr;
      int s;

      if(!fst)
	{
	  if(i_Debug_Loc) Write(Debug,"Error opening file:",fname_state);
	  exit(0);
	}

      // Writing the final states
      fprintf(fst,"  i        v    j      N_i    B/QB    Energy [eV]    k_rec[m^6/s]    k_rec_std[m^6/s]   N_L    N_C    N_D\n");
      int v,j,i;

      double en;
      int nl;

      double k_conv_fac = pow(bo_Angs*Angs_m,6)/au_s;
      missing = 0;
      for(int istate=0; istate<FinState[r].size(); istate++)
	{
	  i = FinState[r][istate];
	  itovj(v,j,i);

	  nl = -1;
	  for(int l=0; l<levels[r].size(); l++)
	    {
	      if(int(levels[r][l][0]) == i)
		{
		  nl = l;
		  break;
		}
	    }
	  if(nl == -1)
	    {
	      //if(i_Debug_Loc) Write(Debug,"WARNING. Level not found in levels file !!!!!!!!!!!!!!");
	      missing++;
	      en = 0;
	    }
	  else
	    {
	      en = levels[r][nl][1] * Hartree_eV ;
	    }

	  fprintf(fst,"%6d  %4d  %4d  %6d  %4d     %10.6E     %10.6E    %10.6E   %4d   %4d   %4d\n",i,v,j,StateCount[r][istate][0], StateCount[r][istate][1], en, k_state[r][n_paths][istate]*k_conv_fac, k_state_var[r][n_paths][istate]*k_conv_fac, StateCount[r][istate][2], StateCount[r][istate][3], StateCount[r][istate][4]);

	  //fprintf(fst,"%6d  %4d  %4d  %6d  %4d     %10.6E  %4d   %4d   %4d\n",i,v,j,StateCount[r][istate][0], StateCount[r][istate][1], en, StateCount[r][istate][2], StateCount[r][istate][3], StateCount[r][istate][4]);

	}

      if(i_Debug_Loc) Write(Debug,"Number of levels found by QCT, but not in levels file =",missing);

      fclose(fst);
      if(i_Debug_Loc) Write(Debug,"Written the Final States file");

      // The Rate constants
      if(i_Debug_Loc) Write(Debug,"Writing the overall recombination rate constant.");
      FILE * fks;
      std :: string fname_kstate = Stat_Dir + "Total_Rate_Constant_" + to_string(r+1) + ".out";

      fks = fopen(fname_kstate.c_str(),"w");
      if(!fks)
	{
	  if(i_Debug_Loc) Write(Debug," Error opening file:",fname_kstate);
	  exit(0);
	}

      fprintf(fks,"   Temp       path_idx   k_rec[m^6/s]   k_rec_std[m^6/s]\n");
      
      int beg;
      if (resolve_path) beg = 0;
      else              beg = n_paths;
       
      for (int p=beg; p<n_paths+1; p++)
	{
	  fprintf(fks,"%.6f      %2d       %10.6E    %10.6E\n", Temp, p, krec[r][p]*k_conv_fac, krec_var[r][p]*k_conv_fac);
	}

      fclose(fks);

      if(i_Debug_Loc) Write(Debug, "Done writing the overall recombination rate constant.");
      if(i_Debug_Loc) Write(Debug, "Done writing the output files for reaction #",r);
      if(i_Debug_Loc) Write(Debug, " ");
    }

  if(i_Debug_Loc) Write(Debug,"Exiting");
}
/**-------------------------------------------------------------------------------------------------------**/
/**-------------------------------------------------------------------------------------------------------**/

/**-------------------------------------------------------------------------------------------------------**/
void Statistics :: Merge_Traj_files(Input_Class* Input)
{
  int i_Debug_Loc = 1;
  std :: string Debug = " [Merge_Traj_files] :";
  if(i_Debug_Loc) Write(Debug,"Entering");

  Stat_Dir = Input->Source_Dir + "Statistics/";
  if(i_Debug_Loc) Write (Debug, "Stat_Dir = ", Stat_Dir);

  std :: string command;

  if(i_Debug_Loc) Write(Debug,"Checking if directory exists");
  struct stat buffer;
  int exist = stat((Stat_Dir).c_str(), &buffer);
  if(exist == 0)
    {
      if(i_Debug_Loc) Write(Debug, "Directory exists. Deleting it now !!");
      command = "rm -r " + Stat_Dir;
      system(command.c_str());
    }

  command = "mkdir -p " + Stat_Dir;

  if(i_Debug_Loc) Write(Debug," Creating Stat directory : ", Stat_Dir);
  system(command.c_str());

  Traj_tot_fname = Stat_Dir + "Trajectories-3B-Tot.out";
  if(i_Debug_Loc) Write(Debug, "Traj_tot_fname = ",Traj_tot_fname);

  // Merging Trajectories-3B.out
  if(i_Debug_Loc) Write(Debug, "Merging Trajectories-3B.out files.");
  int icase = Input->case_beg;
  std :: string case_file = Input->Source_Dir + Input->Temp_pref + to_string(icase) + "/Stat_files/Trajectories-3B.out";
  command = "scp " + case_file + " " + Traj_tot_fname;
  std :: system(command.c_str());

  for (int i=2; i<=Input->ncases; i++)
    {
      icase = Input->case_beg + i - 1;
      case_file = Input->Source_Dir +Input->Temp_pref + to_string(icase) +"/Stat_files/Trajectories-3B.out";
      command = "tail -n+2 " + case_file + " >> " + Traj_tot_fname;
      std :: system(command.c_str());
    }
  if(i_Debug_Loc) Write(Debug, "Done Merging Trajectories-3B.out files.");
  
  if(i_Debug_Loc) Write(Debug,"Exiting");
}
/**-------------------------------------------------------------------------------------------------------**/
void Statistics :: Read_Traj_tot()
{
  Traj = new double* [NTraj];
  
  for (int i=0; i<NTraj; i++)
    {
      Traj[i] = new double [ncols];
      for(int j=0; j<ncols; j++)
	{
	  Traj[i][j] = 0;
	}
    }

  ifstream ftraj;
  ftraj.open(Traj_tot_fname.c_str());
  if(!ftraj)
    {
      Write("Error opening file : ",Traj_tot_fname);
      exit(0);
    }

  int r_skip = 1;
  std :: string line;
  for (int i=0; i<r_skip; i++)
    std :: getline(ftraj, line);
  
  NRec = 0; 
  for (int i=0; i<NTraj; i++)
    {
      for (int j=0; j<ncols; j++)
	{
	  ftraj>>Traj[i][j];
	}
      if(Traj[i][13] > 0)  // Path_id
	NRec++;
    }
  ftraj.close();
}

/**-------------------------------------------------------------------------------------------------------**/
void Statistics :: SetRings()
{
  // This function sets the information regarding the impact parameter rings
  int i_Debug_Loc = 1;
  std :: string Debug = "  [SetRings] :";

  if(i_Debug_Loc) Write(Debug,"Entering");

  // First determine the number of rings
  b1_max.push_back(0);
  b2_max.push_back(0);
  b1_max.push_back(Traj[0][5]);
  b2_max.push_back(Traj[0][7]);

  std :: vector<double> :: iterator it1, it2;

  for(int i=0; i<NTraj; i++)
    {
      it1 = find(b1_max.begin(),b1_max.end(), Traj[i][5]);
      if(it1 == b1_max.end())
	{
	  b1_max.push_back(Traj[i][5]);
	  //if(i_Debug_Loc) Write(Debug,"New ring found. Adding ring with b1_max =",Traj[i][5]);
	}

      it2 = find(b2_max.begin(),b2_max.end(), Traj[i][7]);
      if(it2 == b2_max.end())
	{
	  b2_max.push_back(Traj[i][7]);
	  //if(i_Debug_Loc) Write(Debug,"New ring found. Adding ring with b2_max =",Traj[i][7]);
	}

      if(Traj[i][5] > 20 )
	{
	  Write(Debug,"bmax_1 > 20 for iTraj =",i,"with bmax_1 =",Traj[i][5],"Loc_Traj_ID =",Traj[i][2]);
	  exit(0);
	}
    }

  nrings1 = b1_max.size() - 1;
  nrings2 = b2_max.size() - 1;

  /*===========================Sort the rings =========================================*/
  sort(b1_max.begin(), b1_max.end());
  sort(b2_max.begin(), b2_max.end());

  /*========================== Calculate the area of the rings ======================*/
  Area_1 = new double [nrings1];
  Area_2 = new double [nrings2];
  
  for(int i=0; i<nrings1;i++)
    {
      Area_1[i] = M_PI * (pow(b1_max[i+1],2) - pow(b1_max[i],2)) ;
    }

  for(int i=0; i<nrings2;i++)
    {
      Area_2[i] = M_PI * (pow(b2_max[i+1],2) - pow(b2_max[i],2)) ;
    }

  if(i_Debug_Loc)
    {
      Write(Debug,"Number of rings for b1, nrings_1 =",nrings1);
      Write(Debug,"Number of rings for b2, nrings_2 =",nrings2);

      cout<<Debug<<" b1_max [] = ";
      for (int i=0; i<b1_max.size(); i++) cout<<b1_max[i]<<" ";
      cout<<endl<<Debug<<" b2_max [] = ";
      for (int i=0; i<b2_max.size(); i++) cout<<b2_max[i]<<" ";
      cout<<endl;

      cout<<Debug<<" Area_1 [] = ";
      for (int i=0; i<nrings1; i++) cout<<Area_1[i]<<" ";
      cout<<endl<<Debug<<" Area_2 [] = ";
      for (int i=0; i<nrings2; i++) cout<<Area_2[i]<<" ";
      cout<<endl;

    }
  
  if(i_Debug_Loc) Write(Debug,"Exiting");
}
/**-------------------------------------------------------------------------------------------------------**/
void Statistics :: ReadLevels(Input_Class* Input)
{
  // This function reads the levels file
  int i_Debug_Loc = 1;
  std :: string Debug = "  [Read_Levels] :";

  if(i_Debug_Loc) Write(Debug,"Entering");
  
  std :: string fname, line;
  int r_skip = 15;
  ifstream fin;
  double temp[11];

  // Initializing the levels
  levels.resize(React->num_reactions);

  for (int m=0; m<React->num_reactions; m++)
    {
      
      if(i_Debug_Loc) Write(Debug, "Reading the levels file for molecule #", m);
      if(i_Debug_Loc) Write(Debug, "Levels file name : ",React->levels_fname[m]);

      nlevels = Find_length(React->levels_fname[m],15);
      if(i_Debug_Loc) Write(Debug,"Number of levels in levels file nlevels =",nlevels);

      fin.open(React->levels_fname[m].c_str());
      if(!fin)
	{
	  Write(Debug,"Error opening file:",React->levels_fname[m]);
	  exit(0);
	}

      for(int i=0; i<r_skip; i++)
	{
	  std :: getline(fin,line);        // Comments
	}

      levels[m].resize(nlevels);

      for(int i=0; i<nlevels; i++)
	{
	  levels[m][i].resize(2);

	  for(int j=0; j<11; j++)
	    {
	      fin>>temp[j];
	    }

	  vib = int(temp[0]);
	  rot = int(temp[1]);

	  istate = vjtoi(vib,rot);

	  levels[m][i][0] = istate;
	  levels[m][i][1] = temp[2];
	} 
      fin.close();
      if(i_Debug_Loc) Write(Debug, " Done reading the levels file for molecule #", m);
    }
  if(i_Debug_Loc) Write(Debug,"Exiting");
}
/**-------------------------------------------------------------------------------------------------------**/
int Statistics :: vjtoi(int v, int j)
{
  return QNmax*v+j;
}

void Statistics :: itovj(int& v, int& j, int& i)
{
  v = i/QNmax;
  j = i%QNmax;
}
/**-------------------------------------------------------------------------------------------------------**/
/**-------------------------------------------------------------------------------------------------------**/
/**-------------------------------------------------------------------------------------------------------**/
/**-------------------------------------------------------------------------------------------------------**/
/**-------------------------------------------------------------------------------------------------------**/
