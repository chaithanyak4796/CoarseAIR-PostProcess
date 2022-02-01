#include<iostream>
#include <cstdlib>
#include<fstream>
#include<math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include "Input_Class.h"
#include "Main.h"
#include "Global.h"
#include "Statistics.h"
#include "Logger.h"
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
using namespace std;

int main(int argc, char *argv[])
{
  time_t tbeg, tend;
  time(&tbeg);
  
  int i_Debug_Loc = 1;
  std :: string Debug = "[Main] : ";
  std :: string Inp_fname = "./Input_files/Files.inp";

  Input_Class* Input = new Input_Class;
  if(i_Debug_Loc) Write(Debug, "Calling ReadInput() ");
  Input->Read_Input(Inp_fname);
  if(i_Debug_Loc) Write(Debug, "Done Reading Inputs");

  if(i_Debug_Loc) Write(Debug, "Number of procs = ",Input->NProcs);

  

  #pragma omp parallel
  {
    #pragma omp for
    for (int iP=1; iP<=Input->NProcs; iP++)
      {
	Statistics* Stat = new Statistics;
	//sleep(1);
	Stat->Initialize_Statistics(iP, Input);
	Stat->Process_Trajs(iP, Input, 0);
	delete Stat;
    }
  }

  // Combining the output files
  Combine_Output(Input);

  time(&tend);
  int tot_time = difftime(tend,tbeg);
  int hours    = tot_time/3600;
  int minutes  = (tot_time-hours*3600)/60;
  int seconds  =  (tot_time)%60;

  cout<<endl<<"Total Time Taken = "<<hours<<" hours "<<minutes<<" minutes "<<seconds<<" seconds"<<endl;

  return 0;
}
