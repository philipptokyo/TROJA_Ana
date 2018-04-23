#ifndef __COMPOUND_HH
#define __COMPOUND_HH

#include <iostream>
#include <math.h>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <fstream>

#include "Nucleus.hh"

class Compound{
 public:
  Compound(char*);
  Compound(Nucleus*);
  ~Compound();
  //Compound(int nofelements, Nucleus* nuclei, double* fracs);
  Compound(Nucleus* n1, double f1, Nucleus* n2, double f2);
  void SetNofElements(int nofelements){
    fNofElements = nofelements;
  };
  int GetNofElements(){
    return fNofElements;
  };
  double GetMass(){
    return fMass;
  };
  char* GetSymbol(){
    return fSymbol;
  };
  Nucleus* GetNucleus(int);
  double GetFrac(int);
 private:
  Nucleus** fNuclei;
  double* fFrac;
  int fNofElements;
  double fMass;
  char* fSymbol;
};
#endif
