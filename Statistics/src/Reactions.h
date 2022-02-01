#ifndef REACTIONS_H
#define REACTIONS_H
#include "Input_Class.h"

using namespace std;

class Reactions
{
   public:
      int num_reactions;            // Number of distinct recombination reactions;  
      std :: string* reactions;     // Strings containing reactions
      vector<vector<double>> reaction_arr;        // The arrangement q.n corresponding to each reaction

      double** arr_matrix;  // Matrix for the arrangement q.n
      double symm_fac;

      std :: string* levels_fname;   // String of file names for levels

      Reactions();
      ~Reactions();

      void Initialize_Reactions(Input_Class* Input);
      int Identify_reaction(double arr);
  
};

#endif /* STATISTICS_J */
