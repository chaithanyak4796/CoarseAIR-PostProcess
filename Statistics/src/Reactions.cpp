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
#include "Reactions.h"
#include <algorithm>
#include "Input_Class.h"
#include <unistd.h>

using namespace std;

Reactions :: Reactions() // Constructor
{
  num_reactions = 1;
}

Reactions :: ~Reactions() // Destructor
{
}

/**-------------------------------------------------------------------------------------------------------**/
void Reactions :: Initialize_Reactions(Input_Class* Input)
{
  int i_Debug_Loc = 1;
  std :: string Debug = " [Initialize Reactions] :";
  if(i_Debug_Loc) Write(Debug,"Entering");

    // Initializing the arr_matrix
  arr_matrix = new double *[3];
  for (int i=0; i<3; i++)
    {
      arr_matrix[i] = new double [3];
    }

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


  // Idenitfy the number of recombination channels                                                                                                                                
  if (Input->Atom_Names [0] == 'O' && Input->Atom_Names [1] == 'O' && Input->Atom_Names [2] == 'O')
    {
      num_reactions = 1; 
      reactions = new std :: string [num_reactions];
      reactions[0] = "O + O + O -> O2 + O";
      symm_fac = 0.5;

      levels_fname = new std :: string [num_reactions];
      levels_fname[0] = "./Input_files/levels_O2.inp";

      reaction_arr.resize(1);
      reaction_arr[0].resize(9);
      for(int i=0; i<3; i++)
	{
	  for(int j=0; j<3; j++)
	    reaction_arr[0].push_back(arr_matrix[i][j]);
	}
      
    }
  else if (Input->Atom_Names [0] == 'N' && Input->Atom_Names [1] == 'N' && Input->Atom_Names [2] == 'N')
    {
      num_reactions = 1;
      reactions= new std :: string [num_reactions];
      reactions[0] = "N + N + N -> N2 + N";
      symm_fac = 0.5;

      levels_fname = new std :: string [num_reactions];
      levels_fname[0] = "./Input_files/levels_N2.inp";

      reaction_arr.resize(1);
      reaction_arr[0].resize(9);
      for(int i=0; i<3;i++)
	{
          for(int j=0; j<3; j++)
            reaction_arr[0].push_back(arr_matrix[i][j]);
        }

    }
  else if (Input->Atom_Names [0] == 'N' && Input->Atom_Names [1] == 'N' && Input->Atom_Names [2] == 'O')
    {
      num_reactions = 2;
      reactions= new std :: string [num_reactions];
      reactions[0] = "N + N + O -> N2 + O";
      reactions[1] = "N + N + O -> NO + N";
      symm_fac = 0.5;
      
      levels_fname = new std :: string [num_reactions];
      levels_fname[0] = "./Input_files/levels_N2_for_N2O.inp";
      levels_fname[1] = "./Input_files/levels_NO_for_N2O.inp";

      reaction_arr.resize(2);
      reaction_arr[0].resize(3);
      reaction_arr[1].resize(6);

      for (int i=0; i<3; i++)
	{
	  reaction_arr[0].push_back(arr_matrix[0][i]);
	  reaction_arr[1].push_back(arr_matrix[1][i]);
	  reaction_arr[1].push_back(arr_matrix[2][i]);
	}
    }

  if(i_Debug_Loc)
    {
      Write(Debug, "Number of recombination reactions = ", num_reactions);
      for (int i=0; i< num_reactions; i++)
        {
          Write(Debug, "Reaction [", i, "] : ", reactions[i]);
        }

    }



  if(i_Debug_Loc) Write(Debug,"Exiting");
}

/**-------------------------------------------------------------------------------------------------------**/
int Reactions :: Identify_reaction(double arr)
{
  if(num_reactions == 1)
    return 0;
  else if(num_reactions == 2)
    {
      for (int i=0; i<2; i++)
	{
	  for (int j=0; j<reaction_arr[i].size(); j++)
	    {
	      if (arr == reaction_arr[i][j])
		return i;
	    }
	}
    }

      
}
/**-------------------------------------------------------------------------------------------------------**/
/**-------------------------------------------------------------------------------------------------------**/
/**-------------------------------------------------------------------------------------------------------**/
