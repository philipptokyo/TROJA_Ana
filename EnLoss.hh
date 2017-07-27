#ifndef EnLoss_hh
#define EnLoss_hh

#include "LibPerso.h"

#include "InputInfo.hh"
#include "DetectorInfo.hh"

#define maxEnLossEntriesE 50
#define maxEnLossEntriesX 50

typedef struct _dEdX
{
  
  string materialname;
  string comment;
  
  Int_t ebins, xbins;
  Double_t e[maxEnLossEntriesE];
  Double_t x[maxEnLossEntriesX];
  Double_t enLoss[maxEnLossEntriesE][maxEnLossEntriesX]; // E bin, X bin

} dEdX;



class EnLoss
{
  
  public:

    EnLoss();
    EnLoss(InputInfo *i, DetectorInfo* d);
    ~EnLoss();

    void CollectData(string filename, Int_t index);
    Double_t CalcEnLoss(Double_t en, Double_t dist, Int_t index);
    Double_t CalcEnLoss(Double_t en, Double_t dist, Int_t index, Bool_t warnings);

    Double_t CalcParticleEnergy(Double_t en, Double_t dist, Int_t index);


  private:
    
    void InvertData(Int_t index, Double_t prec);

    InputInfo* fInfo;
    DetectorInfo* fDetInfo;

    dEdX fEnLoss[1];
    dEdX fEnLossInv[1];


};

#endif

