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

  n_t = 0;

  v_ini = 0;
  j_ini = 0;

  v   = 0;
  j   = 0;
  arr = 0;

  //mu1_2  = m1*m2/(m1+m2);
  //mu12_3 = (m1+m2)*m3/(m1+m2+m3);


  r_direct = 10;           // Currently using 8 bohrs. Update later
  //t_cmplx  = 1E-13/au_s;  // Currently using this value (in a.u)
  t_del   = 1.0*t_14;  // Currently using this value (in a.u)
  // The above three values can be updated as per our wish. If I want to give different trehsholds for different energies,
  // just initialize them in the Initialize_Trajectory function()

  cmplx_id = -1;



}

Trajectory :: ~Trajectory() // Destructor
{
  
}

//----------------------------------------------------------------------------------//

void Trajectory :: Initialize_Trajectory(int i, Input_Class* Input, double** Traj_tot, double** PaQSol, double** arr_matrix)
{
  int i_Debug_Loc = 0;
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

  double m1 = Input->Atom_Masses[0];
  double m2 = Input->Atom_Masses[1];
  double m3 = Input->Atom_Masses[2];

  mu1 = m1*m2/(m1+m2);
  mu2 = (m1+m2)*m3/(m1+m2+m3);

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

  for (int p=0; p<3; p++)
    {
      Pini[p+6] = -1*(Input->Atom_Masses[0]*Pini[p] + Input->Atom_Masses[1]*Pini[p+3])/Input->Atom_Masses[2];
      Qini[p+6] = -1*(Input->Atom_Masses[0]*Qini[p] + Input->Atom_Masses[1]*Qini[p+3])/Input->Atom_Masses[2];
      Pfin[p+6] = -1*(Input->Atom_Masses[0]*Pfin[p] + Input->Atom_Masses[1]*Pfin[p+3])/Input->Atom_Masses[2];
      Qfin[p+6] = -1*(Input->Atom_Masses[0]*Qfin[p] + Input->Atom_Masses[1]*Qfin[p+3])/Input->Atom_Masses[2];
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
  E2 = pow((0.5*(Pini[0]+Pini[3]) - Pini[6]),2) + pow((0.5*(Pini[1]+Pini[4]) - Pini[7]),2) + pow((0.5*(Pini[2]+Pini[5]) - Pini[8]),2);

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
      R2[p] = Qini[p+6] - 0.5*(Qini[p] + Qini[p+3]);
      V2[p] = Pini[p+6] - 0.5*(Pini[p] + Pini[p+3]);
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
      int mol1,mol2,mol3;
      mol1 = -1; mol2 = -1; mol3 = -1;
      if( arr == arr_matrix[0][0] || arr == arr_matrix[0][1] || arr == arr_matrix[0][2] )
	{
	  mol1 = 0;
	  mol2 = 3;
	  mol3 = 6;
	}
      if( arr == arr_matrix[1][0] || arr == arr_matrix[1][1] || arr == arr_matrix[1][2] )
	{
	  mol1 = 0;
	  mol2 = 6;
	  mol3 = 3;
	}
      if( arr == arr_matrix[2][0] || arr == arr_matrix[2][1] || arr == arr_matrix[2][2] )
	{
	  mol1 = 3;
	  mol2 = 6;
	  mol3 = 0;
	}

      if(i_Debug_Loc) Write(Debug,"mol1 = ",mol1,", mol2 = ",mol2,", mol3 = ",mol3);

      E_int = 0;
      double V2_CM = 0;

      double E_3 = 0;

      for (int p=0; p<3; p++)
	{
	  E_3+= pow(Pfin[p+mol3],2);
	  V2_CM+= pow(0.5*(Pfin[p+mol1]+Pfin[p+mol2]),2);
	}

      E_int = H_fin - 0.5*m1*E_3 - 0.5*(2*m1)*V2_CM;
    }
  else
    {
      E_int = 0;
    }

  if(i_Debug_Loc) Write(Debug,"E_int [eV] = ",E_int*Hartree_eV);  // This is negative in some cases. Not a problem. Has to do with the reference values.

  // Calculate the initial distances
  d1 = pow(Qini[0] - Qini[3],2) + pow(Qini[1] - Qini[4],2) + pow(Qini[2] - Qini[5],2);
  d2 = pow(Qini[6] - 0.5*(Qini[0]+Qini[3]),2) + pow(Qini[7] - 0.5*(Qini[1]+Qini[4]),2) + pow(Qini[8] - 0.5*(Qini[2]+Qini[5]),2);

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
  int i_Debug_Loc = 0;
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
  int i_Debug_Loc = 0;
  std :: string Debug = "  [Calc_Traj_Params]";

  if(i_Debug_Loc) Write(Debug,"Entering");

  // Calculating omega
  double t1,t2;
  double r1[3]  = {Qini[3]-Qini[0], Qini[4]-Qini[1], Qini[5]-Qini[2]};
  double rcm[3] = {0.5*(Qini[3]+Qini[0]), 0.5*(Qini[4]+Qini[1]), 0.5*(Qini[5]+Qini[2])};
  double r2[3]  = {Qini[6]-rcm[0], Qini[7]-rcm[1], Qini[8]-rcm[2]};

  double v1[3]  = {Pini[3]-Pini[0], Pini[4]-Pini[1], Pini[5]-Pini[2]};
  double vcm[3] = {0.5*(Pini[3]+Pini[0]), 0.5*(Pini[4]+Pini[1]), 0.5*(Pini[5]+Pini[2])};
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
	int i_Debug_Loc = 0;
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
	
