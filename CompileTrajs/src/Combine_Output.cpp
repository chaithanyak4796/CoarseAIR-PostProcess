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

using namespace std;

void Combine_Output(Input_Class* Input)
{
  int i_Debug_Loc = 1;
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

  // Merging Trajectories-3B.out
  int iProc = 1;
  std :: string Proc_pref   = Input->Main_Dir + Input->Bins_Dir + "Node_1/Proc_";
  std :: string Merged_file = Output_Dir + "Trajectories-3B.out";
  std :: string Proc_file, line;
  command = "scp " + Proc_pref + to_string(iProc) + "/Trajectories-3B.out  " + Merged_file;
  std :: system(command.c_str());

  for (iProc=2; iProc <= Input->NProcs; iProc++)
    {
      command = "tail -n+2 " + Proc_pref + to_string(iProc) + "/Trajectories-3B.out  >> " + Merged_file;
      std :: system(command.c_str());
    }
  
  if(i_Debug_Loc) Write(Debug, "Done comibining Trajectories-3B.out");
  
  if(i_Debug_Loc) Write(Debug,"Exiting");
}
