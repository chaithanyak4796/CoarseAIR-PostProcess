#include<iostream>
#include<cstdlib>
#include<fstream>
#include<math.h>
#include<sstream>
#include <stdlib.h>
#include <vector>
//#include "Main.h"
#include "Global.h"
#include "Logger.h"
#include "Input_Class.h"
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

Input_Class :: Input_Class() // Constructor
{
  NProcs = 20;
  NAtoms = 3;
    
}

void Input_Class :: Read_Input(const std::string& Inp_fname)
{
  int i_Debug_Loc = 1;
  std :: string Debug = "  [Read_Input] : ";
  if(i_Debug_Loc) Write();
  if(i_Debug_Loc) Write(Debug,"Entering");

  if(i_Debug_Loc) Write(Debug,"Reading ", Inp_fname);

  ifstream finp;
  finp.open(Inp_fname.c_str());
  if(!finp)
    {
      Write(Debug,"Error opening file : ", Inp_fname);
      exit(0);
    }

  std :: string line;
  std :: getline(finp,line);      // Header
  std :: getline(finp,line);      // Blank line
  std :: getline(finp,line);      // Comment
  std :: getline(finp,Main_Dir);  // Main_Dir
  if(i_Debug_Loc) Write(Debug,"Main_Dir =", Main_Dir);

  std :: getline(finp,line);     // Blank line
  std :: getline(finp,line);     // Comment
  std :: getline(finp,Bins_Dir); // Bins_Dir
  if(i_Debug_Loc) Write(Debug,"Bins_Dir =", Bins_Dir);
  
  std :: getline(finp,line);     // Blank line
  std :: getline(finp,line);     // Comment
  std :: getline(finp,line);     // NProcs
  NProcs = stoi(line);
  if(i_Debug_Loc) Write(Debug,"Number of procs = ",NProcs);

  std :: getline(finp,line);     // Blank line
  std :: getline(finp,line);     // Comment
  std :: getline(finp,line);     // determine_pathway
  determine_pathway = stoi(line);
  if(i_Debug_Loc) Write(Debug,"Determine Pathway? = ",determine_pathway);

  std :: getline(finp,line);     // Blank line
  std :: getline(finp,line);     // Comment
  std :: getline(finp,line);     // NAtoms
  NAtoms = stoi(line);
  if(i_Debug_Loc) Write(Debug,"Number of atoms = ",NAtoms);

  Atom_Names  = new char   [NAtoms];
  Atom_Masses = new double [NAtoms];
  
  std :: getline(finp,line);     // Blank line
  std :: getline(finp,line);     // Comment
  for (int i=0; i<NAtoms; i++)
    {
      std :: getline(finp,line);  // Atom_Names
      Atom_Names[i] = line[0];     
      if(i_Debug_Loc)  Write(Debug,"Atom_Names[",i,"]   = ",Atom_Names[i]);

      if     (Atom_Names[i] == 'O') Atom_Masses[i] = mass_O;
      else if(Atom_Names[i] == 'N') Atom_Masses[i] = mass_N;
      else
	{
	  Write(Debug, "No atom mass found for given atom name : ",Atom_Names[i]);
	  exit(0);
	}
      if(i_Debug_Loc) Write(Debug,"Atoms_Masses[",i,"] = ",Atom_Masses[i]);
    }

  std :: getline(finp,line);    // Blank line
  std :: getline(finp,line);    // Comment
  std :: getline(finp,line);    // eps_LJ
  eps_LJ = stod(line);
  std :: getline(finp,line);    // sig_LJ
  sig_LJ = stod(line);

  if(i_Debug_Loc) Write(Debug, "eps_LJ [K] = ",eps_LJ);
  if(i_Debug_Loc) Write(Debug, "sig_LJ [A] = ",sig_LJ);

  eps_LJ *= Kb;
  sig_LJ *= Angs_m;

  if(i_Debug_Loc) Write(Debug, "eps_LJ [J] = ",eps_LJ);
  if(i_Debug_Loc) Write(Debug, "sig_LJ [m] = ",sig_LJ);

  
  if(i_Debug_Loc) Write(Debug,"Exiting");
}
