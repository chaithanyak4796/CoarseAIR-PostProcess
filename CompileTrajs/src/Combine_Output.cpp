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
#include <bits/stdc++.h>
#include <iostream>

using namespace std;

void Combine_Output(Input_Class* Input)
{
  int i_Debug_Loc = 1;
  int del_Proc_files = 1;  // delete the processor level output files
  
  std :: string Debug = "[Combine_Output] :";
  if(i_Debug_Loc) Write(Debug,"Entering");

  std :: string Output_Dir = Input->Main_Dir + "Stat_files/";
  std :: string command;
  
  if(i_Debug_Loc) Write(Debug,"Checking if directory exists");
  struct stat buffer;
  int exist = stat((Output_Dir).c_str(), &buffer);
  if(exist == 0)
    {
      if(i_Debug_Loc) Write(Debug, "Directory exists. Deleting it now !!");
      command = "rm -r " + Output_Dir;
      system(command.c_str());
    }
  
  command = "mkdir -p " + Output_Dir;
  
  if(i_Debug_Loc) Write(Debug," Creating output directory : ", Output_Dir);
  system(command.c_str());

  std :: string Proc_pref   = Input->Main_Dir + Input->Bins_Dir + "Node_1/Proc_";
  
  // Merging Trajectories-3B.out
  Merge_files(Proc_pref, Input->NProcs, Output_Dir, "Trajectories-3B.out", del_Proc_files);
  if(i_Debug_Loc) Write(Debug, "Done combining Trajectories-3B.out");

  // Merging Path-specific trajectory files
  if(Input->determine_pathway && Input->write_path_out)
    {
      Merge_files(Proc_pref, Input->NProcs, Output_Dir, "Lindemann_Traj.out", del_Proc_files);
      Merge_files(Proc_pref, Input->NProcs, Output_Dir, "Chaperon_Traj.out", del_Proc_files);
      Merge_files(Proc_pref, Input->NProcs, Output_Dir, "Direct_Traj.out", del_Proc_files);
      Merge_files(Proc_pref, Input->NProcs, Output_Dir, "No_Recomb_Traj.out", del_Proc_files);
    }
  if(i_Debug_Loc) Write(Debug, "Done combining path-specific trajectory files");

  Merge_Count_files(Proc_pref, Input->NProcs, Output_Dir, del_Proc_files);
  
  if(i_Debug_Loc) Write(Debug,"Exiting");
}


//-----------------------------------------------------------------------------------------------------------------------------//

void Delete_Proc_files(const std :: string& file_name, int del)
{
  if(del)
    {
      std :: string cmd = "rm " + file_name;
      std :: system(cmd.c_str());
    }
}

//-----------------------------------------------------------------------------------------------------------------------------//

void Merge_files(const std :: string& Proc_pref, int NProcs, const std :: string& Output_Dir, const std :: string& file_name, int del_Proc_files)
{
  int iProc = 1;
  
  std :: string Merged_file = Output_Dir + file_name;
  std :: string Proc_file, command;

  Proc_file = Proc_pref + to_string(iProc) + "/" + file_name;
  command = "scp " + Proc_file + " " +  Merged_file;
  std :: system(command.c_str());
  Delete_Proc_files(Proc_file, del_Proc_files);

  for (iProc=2; iProc <= NProcs; iProc++)
    {
      Proc_file = Proc_pref + to_string(iProc) + "/" + file_name;
      command ="tail -n+2 " + Proc_file + " >> " + Merged_file;
      std :: system(command.c_str());
      Delete_Proc_files(Proc_file, del_Proc_files);
    }
  
}

//-----------------------------------------------------------------------------------------------------------------------------// 

void Merge_Count_files(const std :: string& Proc_pref, int NProcs, const std :: string& Output_Dir, int del_Proc_files)
{

  std :: string Merged_file = Output_Dir + "Count.out";
  std :: string Proc_file, command;

  int count[3];
  double count_p[3];
  string c1, c2;
  int count_tot[3] = {0,0,0};
  int tot = 0;
  int f = 100;

  for (int iProc = 1; iProc <= NProcs; iProc++)
    {
        Proc_file = Proc_pref + to_string(iProc) + "/Count.out";
	std :: ifstream fread(Proc_file);
	if(!fread.fail())
	  {
	    for (int i=0; i<3; i++)
	      {
		fread >> c1 >> c2 >> count[i] >> count_p[i];
		count_tot[i]+= count[i];
		tot+= count[i];
	      }
	  }
	fread.close();
    }

  // Write the overall Count info
  FILE * fc;
  fc = fopen(Merged_file.c_str(), "w");
  fprintf(fc, "L : %5d  %6.4f\n", count_tot[0], double(f*count_tot[0])/tot);
  fprintf(fc, "C : %5d  %6.4f\n", count_tot[1], double(f*count_tot[1])/tot);
  fprintf(fc, "D : %5d  %6.4f\n", count_tot[2], double(f*count_tot[2])/tot);
  fprintf(fc, "T : %5d\n",        tot);

  fclose(fc);

}

//-----------------------------------------------------------------------------------------------------------------------------// 
