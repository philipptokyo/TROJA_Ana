#ifndef __RECONSTRUCTION_HH
#define __RECONSTRUCTION_HH

#include <iostream>
#include <fstream>
#include <string>

#include "TMath.h"
#include "TSpline.h"
#include "TGraph.h"

#include "Nucleus.hh"
#include "Compound.hh"


#ifndef PI
#define PI                       (TMath::Pi())
#endif

using namespace std;

class Reconstruction {
public:
  Reconstruction();
  Reconstruction(Nucleus* projectile, Compound* target);
  Reconstruction(Nucleus* projectile, Compound* target, double thickness);
  void SetTargetThickness(double thickness){
    fTargetThickness = thickness;
  }
  void SetProj(Nucleus* projectile){
    fProj = projectile;
  }
  void SetTarget(Compound* target){
    fTarget = target;
  }
  double StoppingPower(double energy, bool gaseous);
  double StoppingPower(Nucleus* target, double energy, bool gaseous);
  double CompoundRange(double energy, int limit);
  double EnergyAfter(double energy, int limit);
  double EnergyLoss(double energy, int limit){
    return energy - EnergyAfter(energy, limit);
  };
  
  double EnergyStraggling(double dE_dx, double dE_dx_after_target, double energy_loss);
  double AngularStraggling(double energy);
  double ChargeState(double energy);
  void Print(double energy);
  void Print(double from, double to, int steps);
  TSpline3* Energy2Range(double emax, double size);
  TSpline3* Range2Energy(double emax, double size);
  TSpline3* Energy2EnergyLoss(double emax, double size);
  TSpline3* EnergyAfter2EnergyLoss(double emax, double size);
  TSpline3* EnergyLoss2Energy(double emax, double size);
  TSpline3* Energy2EnergyAfter(double emax, double size);
  TGraph* EnergyAfter2Energy(double emax, double size);
private:
  double fTargetThickness;
  Nucleus* fProj;
  Compound* fTarget;
  double a_h(int index, int z);
  double a_he(int index, int z);
  double b_he(int index, int z);
  double shell_correction(int z);
  //ClassDef(Reconstruction, 1);
};

#endif
