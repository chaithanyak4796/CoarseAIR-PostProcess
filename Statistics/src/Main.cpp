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
#include "Input_Class.h"

using namespace std;

int main()
{
  time_t tbeg, tend;
  time(&tbeg);

  int i_Debug_Loc = 1;
  std :: string Debug = "[Main] : ";
  std :: string Inp_fname = "./Input_files/Stat_Input.inp";

  Input_Class* Input = new Input_Class;
  if(i_Debug_Loc) Write(Debug, "Calling ReadInput() ");
  Input->Read_Input(Inp_fname);
  if(i_Debug_Loc) Write(Debug, "Done Reading Inputs");

  if(i_Debug_Loc) Write(Debug, "Initializing Statistics class");
  Statistics* Stat = new Statistics;
  Stat->Initialize_Statistics(Input);
  
  Stat->Process_Statistics(Input);
  
  time(&tend);
  int tot_time = difftime(tend,tbeg);
  int hours    = tot_time/3600;
  int minutes  = (tot_time-hours*3600)/60;
  int seconds  =  (tot_time)%60;

  cout<<endl<<"Total Time Taken = "<<hours<<" hours "<<minutes<<" minutes "<<seconds<<" seconds"<<endl;

  return 0;
}
