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
#include "Trajectory.h"
#include "Logger.h"
#include <algorithm>
#include "Statistics.h"
#include <list>
#include "Input_Class.h"
using namespace std;

Trajectory:: Trajectory()  // Constructor
{
  iNode = 0;
  iProc = 0;
  iTraj = 0;
  conv  = 0;

  iPES = 1;

  b1 = 0;
  b2 = 0;

  b1_max = 0;
  b2_max = 0;

  E1  = 0;
  E2 = 0;

  E_kin_in  = 0;
  E_rel_fin = 0;

  H_fin = 0;

  gamma = 0;

  d1 = 0;
  d2 = 0;

  tau = 0;

  nsteps = 0;

  v_ini = 0;
  j_ini = 0;

  v   = 0;
  j   = 0;
  arr = 0;

  cmplx_id = -1;

  path_idx = 0;
  ratio    = 0.0;

}

Trajectory :: ~Trajectory() // Destructor
{
  /*if(recomb_check)
    {
      for(int i=0; i<nsteps; i++)
	{
	  delete[] P[i];
	  delete[] Q[i];
	  delete[] R[i];
	  delete[] Acc[i];
	}
      delete[] t;
      delete[] P;
      delete[] Q;
      delete[] R;
      delete[] Acc;
      }*/
  
}

//----------------------------------------------------------------------------------//

void Trajectory :: Initialize_Trajectory(int i, Input_Class* Input, double** Traj_tot, double** PaQSol, double** arr_matrix)
{
  int i_Debug_Loc = i_Debug;
  //int i_Debug_Loc = 1;
  std :: string Debug = "   [Initialize_Trajectory] :";
  if(i_Debug_Loc) Write(Debug, "Entering");

  iTraj  = int(Traj_tot[i][0]);
  iPES   = int(Traj_tot[i][1]);
  b1_max = Traj_tot[i][2];
  b1     = Traj_tot[i][3];
  b2_max = Traj_tot[i][4];
  b2     = Traj_tot[i][5];
  j_ini  = Traj_tot[i][6];
  v_ini  = Traj_tot[i][7];
  arr_ini= Traj_tot[i][8];
  j      = Traj_tot[i][9];
  v      = Traj_tot[i][10];
  arr    = Traj_tot[i][11];
  iNode  = int(Traj_tot[i][12]);
  iProc  = int(Traj_tot[i][13]);
  NCoords= 3*Input->NAtoms;
  NAtoms = Input->NAtoms;;
  
  double m1 = Input->Atom_Masses[0];
  double m2 = Input->Atom_Masses[1];
  double m3 = Input->Atom_Masses[2];

  mu1 = m1*m2/(m1+m2);
  mu2 = (m1+m2)*m3/(m1+m2+m3);
  M = m1 + m2;

  if(i_Debug_Loc) Write(Debug,"Initializing Traj # ",i,",iTraj = ",iTraj);
  if(i_Debug_Loc) Write(Debug,"Before adjusting QNs: v =",v,"j =",j);
  Adjust_QN();
  if(i_Debug_Loc) Write(Debug,"After adjusting QNs: v =",v,"j =",j);

  // Check if the trajectory leads to a recombination event
  recomb_check=0;
  for (int p=0; p<3; p++)
    {
      for(int q=0; q<3; q++)
	{
	  if(arr == arr_matrix[p][q])
	    recomb_check = 1;
	}
    }
  if(i_Debug_Loc) Write(Debug, "recomb_check = ", recomb_check);

  // Reading the initial and final PaQ
  Pini = new double [NCoords];
  Qini = new double [NCoords];
  Pfin = new double [NCoords];
  Qfin = new double [NCoords];

  P1 = new double [NCoords];
  P2 = new double [NCoords];
  Q1 = new double [NCoords];
  Q2 = new double [NCoords];

  for(int p=0; p<NCoords; p++)
    {
      Pini[p] = 0;
      Qini[p] = 0;
      Pfin[p] = 0;
      Qfin[p] = 0;
    }

  for (int p=0; p<NCoords-3; p++)
    {  // Indices assume NAtoms = 3
      Pini[p] = PaQSol[i][p+3];
      Qini[p] = PaQSol[i][p+9];
      Pfin[p] = PaQSol[i][p+16];
      Qfin[p] = PaQSol[i][p+22];
    }

  // Compute initial PaQ of the 3rd atom
  for (int p=0; p<3; p++)
    {
      Pini[p+6] = -1*(Input->Atom_Masses[0]*Pini[p] + Input->Atom_Masses[1]*Pini[p+3])/Input->Atom_Masses[2];
      Qini[p+6] = -1*(Input->Atom_Masses[0]*Qini[p] + Input->Atom_Masses[1]*Qini[p+3])/Input->Atom_Masses[2];
      Pfin[p+6] = -1*(Input->Atom_Masses[0]*Pfin[p] + Input->Atom_Masses[1]*Pfin[p+3])/Input->Atom_Masses[2];
      Qfin[p+6] = -1*(Input->Atom_Masses[0]*Qfin[p] + Input->Atom_Masses[1]*Qfin[p+3])/Input->Atom_Masses[2];
    }

  // Compute the initial P1,P2 and Q1,Q2
  for(int p=0; p<3; p++)
    {
      P1[p] = ( m1*Pini[p] + m2*Pini[p+3] )/M;
      Q1[p] = ( m1*Qini[p] + m2*Qini[p+3] )/M;

      P2[p] = Pini[p+6] - P1[p];
      Q2[p] = Qini[p+6] - Q1[p];
    }

  H_fin = PaQSol[i][15];
  
  if(i_Debug_Loc)
    {
      cout << Debug << " Pini : ";
      for (int p=0; p<NCoords; p++)
	{
	  cout << Pini[p] <<"\t";
	}
      cout<<endl;
    }

  // Calculate the relative translational energies
  E1 = pow((Pini[0] - Pini[3]),2) + pow((Pini[1] - Pini[4]),2) + pow((Pini[2] - Pini[5]),2);
  E2 = pow(P2[0],2) + pow(P2[1],2) + pow(P2[2],2);

  E1 *= 0.5 * mu1;
  E2 *= 0.5 * mu2;

  if(i_Debug_Loc) Write(Debug,"E1 [eV] = ",E1*Hartree_eV);
  if(i_Debug_Loc) Write(Debug,"E2 [eV] = ",E2*Hartree_eV);

  // Verify impact params
  double b1_calc, b2_calc, ct1, ct2;
  double R1[3], R2[3],V1[3], V2[3];
  for (int p=0; p<3; p++)
    {
      R1[p] = Qini[p] - Qini[p+3];
      V1[p] = Pini[p] - Pini[p+3];
      R2[p] = Q2[p];
      V2[p] = P2[p];
    }

  ct1 = pow(dot(R1,V1),2)/(dot(R1,R1)*dot(V1,V1));
  ct2 = pow(dot(R2,V2),2)/(dot(R2,R2)*dot(V2,V2));

  if(ct1 > 1.0) ct1 = 1.0;
  if(ct2 > 1.0) ct2 = 1.0;

  b1_calc = pow(dot(R1,R1),0.5) * pow(1-ct1,0.5);
  b2_calc = pow(dot(R2,R2),0.5) * pow(1-ct2,0.5);

  if(i_Debug_Loc) Write(Debug," Calc b1 = ",b1_calc);
  if(i_Debug_Loc) Write(Debug," Calc b2 = ",b2_calc);

  if( fabs(b1-b1_calc) > 1E-4 || fabs(b2-b2_calc) > 1E-4)
    {
      if(i_Debug_Loc) Write(Debug," ERROR : Mismatch between impact params calcualted from PaQSol and trajectories.out file");
      cout << " ERROR : Mismatch between impact params calcualted from PaQSol and trajectories.out file";
      exit(0);
    }

  // Calculate the internal energy of the molecule
  if(recomb_check)
    {
      path_idx = 1; // In case we check the recombination pathways, this will be updated again.
      
      int mol1,mol2,mol3;
      double m_int_1, m_int_2, m_int_3;

      mol1 = -1; mol2 = -1; mol3 = -1;
      if( arr == arr_matrix[0][0] || arr == arr_matrix[0][1] || arr == arr_matrix[0][2] )
	{
	  mol1 = 0; m_int_1 = m1;
	  mol2 = 3; m_int_2 = m2;
	  mol3 = 6; m_int_3 = m3;
	}
      if( arr == arr_matrix[1][0] || arr == arr_matrix[1][1] || arr == arr_matrix[1][2] )
	{
	  mol1 = 0; m_int_1 = m1;
	  mol2 = 6; m_int_2 = m3;
	  mol3 = 3; m_int_3 = m2;
	}
      if( arr == arr_matrix[2][0] || arr == arr_matrix[2][1] || arr == arr_matrix[2][2] )
	{
	  mol1 = 3; m_int_1 = m2;
	  mol2 = 6; m_int_2 = m3;
	  mol3 = 0; m_int_3 = m1;
	}

      if(i_Debug_Loc) Write(Debug,"mol1 = ",mol1,", mol2 = ",mol2,", mol3 = ",mol3);

      E_int = 0;
      double E_CM  = 0;   // Translational energy of the molecule
      double E_3   = 0;   // Translational energy of the atom
      double m_CM  = m_int_1 + m_int_2; // Mass of the molecule
      double V_CM[3]  = {0,0,0};           // COM  of the molecule

      for (int p=0; p<3; p++)
	{
	  E_3+= 0.5 * m_int_3 * pow(Pfin[p+mol3],2);
	  V_CM[p] = (m_int_1 * Pfin[p+mol1] + m_int_2 * Pfin[p+mol2])/(m_CM);
	}

      for(int p=0; p<3; p++)
	{
	  E_CM += 0.5 * m_CM * pow(V_CM[p],2);
	}

      E_int = H_fin - E_3 - E_CM;
    }
  else
    {
      E_int = 0;
    }

  if(i_Debug_Loc) Write(Debug,"E_int [eV] = ",E_int*Hartree_eV);  // This is negative in some cases. Not a problem. Has to do with the reference values.

  // Calculate the initial distances
  d1 = pow(Qini[0] - Qini[3],2) + pow(Qini[1] - Qini[4],2) + pow(Qini[2] - Qini[5],2);
  d2 = pow(Qini[6] - Q1[0],2)  + pow(Qini[7] - Q1[1],2)  + pow(Qini[8] - Q1[2],2);

  d1 = pow(d1,0.5);
  d2 = pow(d2,0.5);
  if(i_Debug_Loc) Write(Debug,"d1 [Bo] = ",d1);
  if(i_Debug_Loc) Write(Debug,"d2 [Bo] = ",d2);

  // Evaluating Bunker's lifetime and omega
  Eval_tau(Input->eps_LJ, Input->sig_LJ);
  Calc_Traj_Params();
  if(i_Debug_Loc) Write(Debug, "omega = ",omega);

  if(recomb_check == 0)
    {
      v = -1;
      j = -1;
    }
  
  if(i_Debug_Loc) Write(Debug, "Exiting");
}

/*--------------------------------------------------------------------------------------------------------------------------------------------------*/
void Trajectory :: Eval_tau(double e_LJ, double sigma_LJ)
{
  int i_Debug_Loc = i_Debug;
  //int i_Debug_Loc = 0;
  std :: string Debug = "   [Eval_tau] :";

  long double mass = mu1*me;
  long double Etr  = E1*Hartree_eV*eV2J;

  if(i_Debug_Loc) Write(Debug,"mass =",mass);
  if(i_Debug_Loc) Write(Debug,"Etr  =",Etr);

  tau = 1.5*sigma_LJ*pow(mass,0.5);
  tau = tau*pow(e_LJ,1.0/6.0);
  tau = tau*pow(Etr,-2.0/3.0);

  if(i_Debug_Loc) Write(Debug,"tau =",tau/au_s);
  if(i_Debug_Loc) Write(Debug,"Exiting");
}

/*--------------------------------------------------------------------------------------------------------------------------------------------------*/
void Trajectory :: Calc_Traj_Params()
{
  int i_Debug_Loc = i_Debug;
  //int i_Debug_Loc = 0;
  std :: string Debug = "  [Calc_Traj_Params]";

  if(i_Debug_Loc) Write(Debug,"Entering");

  // Calculating omega
  double t1,t2;
  double r1[3]  = {Qini[3]-Qini[0], Qini[4]-Qini[1], Qini[5]-Qini[2]};
  double rcm[3] = {Q1[0], Q1[1], Q1[2]};
  double r2[3]  = {Qini[6]-rcm[0], Qini[7]-rcm[1], Qini[8]-rcm[2]};

  double v1[3]  = {Pini[3]-Pini[0], Pini[4]-Pini[1], Pini[5]-Pini[2]};
  double vcm[3] = {P1[0], P1[1], P1[2]};
  double v2[3]  = {Pini[6]-vcm[0], Pini[7]-vcm[1], Pini[8]-vcm[2]};

  t1 = -dot(r1,v1)/(dot(v1,v1));
  t2 = -dot(r2,v2)/(dot(v2,v2));

  omega = 1*(t2-t1)/(tau/au_s);

  if(i_Debug_Loc) Write(Debug,"t1 [s]  =",t1*au_s);
  if(i_Debug_Loc) Write(Debug,"tau [s] =",tau);
  if(i_Debug_Loc) Write(Debug,"Omega = ",omega);
  
  if(i_Debug_Loc) Write(Debug,"Exiting");
}

/*--------------------------------------------------------------------------------------------------------------------------------------------------*/
void Trajectory :: Adjust_QN()
{
	// This function adjusts the quantum numbers for the trajectories.
        int i_Debug_Loc = i_Debug; 
	//int i_Debug_Loc = 0;
	std :: string Debug = "  [Adjust QN] :";
	
	if(i_Debug_Loc) Write(Debug,"Entering");
	
	v     = floor(v);
	arr   = floor(arr);
	
	double j1,j2,jf;
	
	if(parity == 0)  // Both Odd/Even allowed
	{
		j = floor(j);
	} 
	else if(parity == 1) // Only odd allowed
	{
		j  = j - 0.5;
		j1 = floor(j);
		
		if(int(j1) % 2 == 1 )
		{
			j1 = j1;
			j2 = j1+2;
		}
		else
		{
			j1 = j1 - 1;
			j2 = j1 + 2;
		}
		
		if(j-j1 <= 1)
		{
			j = j1;
		}	
		else
		{
			j = j2;
		}
		
		if(j < 0)
		{
			j = 1;   // To avoid negative j becuase of rounding down.
		}
		
	}
	else if(parity == 2) // Only even allowed
	{
		j  = j - 0.5;
		j1 = floor(j);
		
		if(int(j1) % 2 == 0 )
		{
			j1 = j1;
			j2 = j1+2;
		}
		else
		{
			j1 = j1 - 1;
			j2 = j1 + 2;
		}
		
		if(j-j1 <= 1)
		{
			j = j1;
		}	
		else
		{
			j = j2;
		}
		
		if(j < 0)
		{
			j = 0;   // To avoid negative j becuase of rounding down.
		}
	}

	if(i_Debug_Loc) Write(Debug,"Vibrational QN v   =",v);
	if(i_Debug_Loc) Write(Debug,"Rotational  QN j   =",j);
	if(i_Debug_Loc) Write(Debug,"Arrangement QN arr =",arr);


	if(i_Debug_Loc) Write(Debug,"Exiting");
}
	
/*--------------------------------------------------------------------------------------------------------------------------------------------------*/
int Trajectory :: Determine_pathway(std :: string Proc_Dir, Input_Class* Input, double** arr_matrix)
{
  int i_Debug_Loc = i_Debug;
  //int i_Debug_Loc = 1;
  std :: string Debug = " [Determine pathway] :";
  if(i_Debug_Loc) Write(Debug,"Entering");

  // First read the PaQEvo file
  std :: string fname = Proc_Dir + "PaQEvo-" + to_string(iTraj) + ".out";
  if(i_Debug_Loc) Write(Debug, "PaQEvo file name = ",fname);

  ifstream fcheck;
  fcheck.open(fname.c_str());
  if(!fcheck)
    {
      Write(Debug,"Error opening file :",fname);
      exit(0);
    }
  fcheck.close();

  nsteps = Find_length(fname,1);
  if(i_Debug_Loc) Write(Debug,"Number of timesteps in this trajectory =",nsteps);

  t = new double  [nsteps];
  P = new double* [nsteps];
  Q = new double* [nsteps];
  R = new double* [nsteps];
  Acc = new double* [nsteps];

  for (int i=0; i<nsteps; i++)
    {
      P[i]   = new double [NCoords];
      Q[i]   = new double [NCoords];
      R[i]   = new double [3];   // Has to be adjusted for NAtoms>3
      Acc[i] = new double [3];
    }

  if(i_Debug_Loc) Write(Debug,"Reading the PaQ Evo file and calculating the PaQ for the 3rd atom");
  Read_PaQEvo(fname, Input);
  //if(i_Debug_Loc) Write(Debug,"P[17][2] = ",P[17][2]);

  if(i_Debug_Loc) Write(Debug,"Calculating the inter atomic distances.");
  Calc_inter_distance();
  if(i_Debug_Loc) Write(Debug, "Initial Distances = ", R[0][0], R[0][1], R[0][2]);

  if(i_Debug_Loc) Write(Debug, "Calculating the acclerations");
  Calc_acceleration(0);

  if(i_Debug_Loc) Write(Debug, "Checking if the pathway is Direct");
  int direct_check = Check_Direct();

  if(!direct_check)
    {
      if(i_Debug_Loc) Write(Debug,"The current trajectory corresponds to a Direct 3B recombination !!!!!!");
      path_idx = 3;
    }
  else
    {
      if( arr == arr_matrix[cmplx_id][0] || arr == arr_matrix[cmplx_id][1] || arr == arr_matrix[cmplx_id][2] )
	{
	  if(i_Debug_Loc) Write(Debug,"This trajectory corresponds to a Lindemann Mechanism");
	  path_idx = 1;
	}
      else
	{
	  if(i_Debug_Loc) Write(Debug,"This trajectory corresponds to a Chaperon Mechanism");
	  path_idx = 2;
	}
    }

  return path_idx;
  if(i_Debug_Loc) Write(Debug,"Exiting");
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------*/
int Trajectory :: Check_Direct()
{
  int i_Debug_Loc = i_Debug;
  //int i_Debug_Loc = 1;
  int direct_check = 0 ;

  std :: string Debug = "  [Check_Direct_New]  :";
  if(i_Debug_Loc) Write(Debug,"Entering");

  // Getting t_2B and t_3B
  double t_first [3];   // Array containing the first timestamps when the atoms feel a force
  double t_first_idx [3];   // Array containing the first timestamps (index) when the atoms feel a force
  for (int k=0; k<3; k++)
    {
      t_first[k] = -1;
      for(int i = 1; i<nsteps-1; i++)
	{
	  //if(i_Debug_Loc) Write(Debug,t[i]*au_ps,Acc[i][j]);
	  if(Acc[i][k] > 1.0E-4 )  // threshold of acceleration [Bo/ps^2]
	    {
	      if(i_Debug_Loc) Write(Debug,"Acc[i][k] = ", Acc[i][k]);  // There is some error here. Start from here.
	      t_first[k] = t[i]*au_ps;
	      t_first_idx[k] = i;
	      break;
	    }
	}
    }

  if(i_Debug_Loc) Write(Debug,"The time-stamps (in ps) of first interactions are ",t_first[0],t_first[1],t_first[2]);

  double t_2B, t_3B, t_OP;
  t_2B = 0; t_3B = 0; t_OP = 0;

  t_3B = t_first[get_max(t_first[0],t_first[1],t_first[2])];
  t_2B = t_first[get_min(t_first,3)];

  // Determine the pair that forms a "complex"
  Determine_complex(t_first_idx);
  if(i_Debug_Loc) Write(Debug,"The pair that forms a complex is iPair =",cmplx_id);
  t_OP  = Determine_t_min(cmplx_id);
  t_OP *= au_ps;
  
  if(i_Debug_Loc) Write(Debug,"t_2B = ",t_2B);
  if(i_Debug_Loc) Write(Debug,"t_3B = ",t_3B);

  if(i_Debug_Loc) Write(Debug,"t_OP = ",t_OP);

  ratio = (t_OP - t_3B)/(t_OP - t_2B);

  if(i_Debug_Loc) Write(Debug,"ratio = ",ratio);

  if(ratio < 0.95) // Another treshold
    {
      direct_check++;   // direct_check == 0  ====> Direct recombination
    }
  
  return direct_check;
  if(i_Debug_Loc) Write(Debug,"Exiting");
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------*/
void Trajectory :: Determine_complex(double t_first_idx[3])
{
  // This function determines the id of the OP
  // It looks for the first time stamp when each atoms face a force and identifies the pair that feels the same force at that time.
  
  int i_Debug_Loc = i_Debug;
  //int i_Debug_Loc = 0;
  std :: string Debug = "   [Determine_cmplx_new]  :";
  if(i_Debug_Loc) Write(Debug,"Entering");

  int t2B;
  t2B = get_min(t_first_idx,3);
  t2B = t_first_idx[t2B];
  if(i_Debug_Loc) Write(Debug,"t_first_idx = ",t_first_idx[0],t_first_idx[1],t_first_idx[2]);
  if(i_Debug_Loc) Write(Debug,"t2B = ",t2B);
  if(i_Debug_Loc) Write(Debug,"Acc[t2B] = ",Acc[t2B][0],Acc[t2B][1],Acc[t2B][2]);

  double eps = 1E-5;

  double da12, da13, da23;

  da12 = fabs(Acc[t2B][0]-Acc[t2B][1]);
  da13 = fabs(Acc[t2B][0]-Acc[t2B][2]);
  da23 = fabs(Acc[t2B][1]-Acc[t2B][2]);

  double da[3] = {da12, da13, da23};

  if(i_Debug_Loc)
    {
      Write(Debug,"da12 = ",da12);
      Write(Debug,"da13 = ",da13);
      Write(Debug,"da23 = ",da23);
    }

  cmplx_id = get_min(da,3);

  if(i_Debug_Loc) Write(Debug,"cmplx_id = ",cmplx_id);
  if(i_Debug_Loc) Write(Debug,"Exiting");
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------*/
double Trajectory :: Determine_t_min(int pair)
{
  // This function determines the first instance of local minima for a pair of atoms
  // Currently, this is primarily used for determining the t_OP timestamp

  double t_min = 0;
  
  for (int i=1; i<nsteps-1; i++)
    {
      if(R[i][pair] < R[i-1][pair] && R[i][pair] < R[i+1][pair])
	{
	  t_min = t[i];
	  break;
	}
    }

  return t_min;
}
 
/*--------------------------------------------------------------------------------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------------------------------------*/
void Trajectory :: Calc_inter_distance()
{
  for(int i=0; i<nsteps; i++)
    {
      R[i][0] = pow((Q[i][0] - Q[i][3]),2) + pow((Q[i][1] - Q[i][4]),2) + pow((Q[i][2] - Q[i][5]),2);   // R12
      R[i][1] = pow((Q[i][0] - Q[i][6]),2) + pow((Q[i][1] - Q[i][7]),2) + pow((Q[i][2] - Q[i][8]),2);   // R13
      R[i][2] = pow((Q[i][3] - Q[i][6]),2) + pow((Q[i][4] - Q[i][7]),2) + pow((Q[i][5] - Q[i][8]),2);   // R23

      R[i][0] = sqrt(R[i][0]);
      R[i][1] = sqrt(R[i][1]);
      R[i][2] = sqrt(R[i][2]);
    }
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------*/
void Trajectory :: Calc_acceleration(int method)
{
  // This function calculates the acceleration of all 3 particles from velocities
  // method = 0 : Central Differencing
  // method = 1 : Backward Differencing

  // Initialize the array (Move this to Initialize() later if required)
  double temp[9], arr[3];

  for (int i = 1; i< nsteps-1; i++)
    {
      for (int k=0; k<9; k++)
	{
	  if(method == 0)
	    {
	      temp[k] = (P[i+1][k] - P[i-1][k])/(t[i+1] - t[i-1]);
	    }
	  else if(method == 1)
	    {
	      temp[k] = (P[i][k] - P[i-1][k])/(t[i] - t[i-1]);
	    }

	  temp[k] = temp[k]/pow((au_ps),2);  // Converting to Bo/ps^2
	}

      for(int k=0; k<3; k++)
	{
	  arr[0] = temp[3*k];
	  arr[1] = temp[3*k+1];
	  arr[2] = temp[3*k+2];

	  Acc[i][k] = sqrt(dot(arr,arr));
	}

    }
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------*/
void Trajectory :: Read_PaQEvo(std :: string fname, Input_Class* Input)
{
  int i_Debug_Loc = i_Debug;
  //int i_Debug_Loc = 1;
  std :: string Debug = " [Read_PaQEvo] :";
  if(i_Debug_Loc) Write(Debug,"Entering");
  
  double m1 = Input->Atom_Masses[0];
  double m2 = Input->Atom_Masses[1];
  double m3 = Input->Atom_Masses[2];
  
  ifstream fr;
  fr.open(fname.c_str());
  if(!fr)
    {
      Write(" Error opening file : ",fname);
      exit(0);
    }
    
  int r_skip = 1;
  std :: string line;
  for (int i=0; i<r_skip; i++)
    std :: getline(fr,line);

  int nsteps_ac = nsteps;
  int idx = -1;
  
  double** arr;
  arr = new double* [nsteps];
  for (int i=0; i<nsteps; i++)
    {
      arr[i] = new double[15];
      for(int k=0; k<15; k++)
	{
	  fr >> arr[i][k];
	}

      if(idx >=0 && arr[i][0] == t[idx])
	{
	  // Time steps are repeated
	  nsteps_ac--;
	}
      
      else
	{
	  idx ++;
	  // Copy the data into appropriate arrays
	  t[idx]    = arr[i][0];
	  P[idx][0] = arr[i][2];
	  P[idx][1] = arr[i][3];
	  P[idx][2] = arr[i][4];
	  P[idx][3] = arr[i][5];
	  P[idx][4] = arr[i][6];
	  P[idx][5] = arr[i][7];
	  Q[idx][0] = arr[i][8];
	  Q[idx][1] = arr[i][9];
	  Q[idx][2] = arr[i][10];
	  Q[idx][3] = arr[i][11];
	  Q[idx][4] = arr[i][12];
	  Q[idx][5] = arr[i][13];

	  // Calculate the P and Q for the third atom.

	  P[idx][6] = -(m1/m3) * P[idx][0] - (m2/m3) * P[idx][3];
	  P[idx][7] = -(m1/m3) * P[idx][1] - (m2/m3) * P[idx][4];
	  P[idx][8] = -(m1/m3) * P[idx][2] - (m2/m3) * P[idx][5];

	  Q[idx][6] = -(m1/m3) * Q[idx][0] - (m2/m3) * Q[idx][3];
	  Q[idx][7] = -(m1/m3) * Q[idx][1] - (m2/m3) * Q[idx][4];
	  Q[idx][8] = -(m1/m3) * Q[idx][2] - (m2/m3) * Q[idx][5];
	}
    }

  fr.close();

  for(int i=0; i<nsteps; i++)
    delete[] arr[i];
  delete[] arr;

  if(i_Debug_Loc) Write(Debug," Number of repititions = ", nsteps - nsteps_ac);
  
  nsteps = nsteps_ac;
  
  if(i_Debug_Loc) Write(Debug,"Exiting");
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------*/
void Trajectory :: Evaluate_SAI()
{
  int i_Debug_Loc = i_Debug;
  //int i_Debug_Loc = 1;
  std :: string Debug = " [Evaluate_SAI] :";
  if(i_Debug_Loc) Write(Debug,"Entering");

  n_rev = 0;

  int p1, p2;
  
  if(cmplx_id == 0)
    {
      p1 = 0; p2 = 3;
    }
  else if (cmplx_id == 1)
    {
      p1 = 0; p2 = 6;
    }
  else
    {
      p1 = 3; p2 = 6;
    }

  double vel_ini1[3] = {P[0][p1], P[0][p1+1], P[0][p1+2]};
  double vel_ini2[3] = {P[0][p2], P[0][p2+1], P[0][p2+2]};

  double vel_fin1[3], vel_fin2[3];
  double rel_ini[3], rel_fin[3];

  double temp;
  
  for (int j=0; j<3; j++)
    rel_ini[j] = vel_ini2[j] - vel_ini1[j];
  
  for (int i=1; i<nsteps-1; i++)
    {
      for (int j=0; j<3; j++)
	{
	  vel_fin1[j] = P[i][p1+j];
	  vel_fin2[j] = P[i][p2+j];

	  rel_fin[j] = vel_fin2[j] - vel_fin1[j];
	}
      temp = comp_angle(rel_ini, rel_fin);

      if(SAI >= 1.5*M_PI && SAI < 2*M_PI)  // If SAI is in the 4th quadrant
	{
	  if (temp < 0.5*M_PI) // and completes on revolution
	    {
	      n_rev++;
	    }
	}

      SAI = temp;

      if(t[i] >= t_3B/au_ps)   // Impact
	break;
    }

  SAI += n_rev * 2 * M_PI;

  if(i_Debug_Loc) Write("SAI = ", SAI);
  if(i_Debug_Loc) Write(Debug,"Exiting");
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------*/

double Trajectory :: comp_angle(double v1[3], double v2[3])
{
  // Computes the angle betweeen two vectors in the same plane
  // Angle will be in the range of [0,2*pi)
  
  double cross_prod[3];
  double dot_prod  = 0;
  double cross_mag = 0;

  cross_prod[0] = v1[1]*v2[2] - v1[2]*v2[1];
  cross_prod[1] = v1[2]*v2[0] - v1[0]*v2[2];
  cross_prod[2] = v1[0]*v2[1] - v1[1]*v2[0];

  for (int i=0; i<3; i++)
    {
      cross_mag += pow(cross_prod[i],2);
      dot_prod  += v1[i] * v2[i];
    }

  cross_mag = sqrt(cross_mag);

  double theta = atan2(cross_mag, dot_prod);
  if(theta < 0)
    theta += 2*M_PI;

  return theta;
}
