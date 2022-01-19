#ifndef INPUT_CLASS_H
#define INPUT_CLASS_H

class Input_Class
{
  public:
  
  std :: string Main_Dir;
  std :: string Bins_Dir;

  int NProcs;   // Number of processors
  int NAtoms;  // Number of atoms

  int determine_pathway;  // Flag to determine pathways?

  char* Atom_Names;     // Name of atoms
  double* Atom_Masses;  // Mass of the atoms

  double eps_LJ;
  double sig_LJ;

  //public:
  Input_Class();
  void Read_Input(const std::string& Inp_fname);

};

#endif /* INPUT_H */
