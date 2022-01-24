#ifndef INPUT_CLASS_H
#define INPUT_CLASS_H
#include <vector>
using namespace std;

class Input_Class
{
  public:
  
  std :: string Source_Dir;
  std :: string Temp_pref;
  std :: string levels_fname;
  
  double Temp;       // Temperature
  double omega_min;  // Range of omega (Min)
  double omega_max;  // Range of omega (Max)
  double case_beg;   // Case_idx (beg)
  double case_end;   // Case_idx (end)
  int ncases;        // Number of cases
  
  int NProcs;   // Number of processors
  int NAtoms;  // Number of atoms

  char* Atom_Names;     // Name of atoms
  double* Atom_Masses;  // Mass of the atoms

  double eps_LJ;
  double sig_LJ;

  int resolve_path;           // Resolve recombination pathways?

  int Poisson_treat;    // How to treat the Poisson distribution?

  //public:
  Input_Class();
  void Read_Input(const std::string& Inp_fname);

};

#endif /* INPUT_H */
