#ifndef Analysis_hh
#define Analysis_hh

#include "LibPerso.h"

#include "InputInfo.hh"
#include "DetectorInfo.hh"


class Analysis
{
  
  public:
    
    Bool_t Init();
    void ResetVariables();
    void Analysis1();
    Bool_t GetCuts();
    void Analysis2();


    Analysis(InputInfo *i, DetectorInfo* d);
    
    ~Analysis();


  private:
    
    InputInfo* info;
    DetectorInfo* detInfo;

    TFile* fileBeam;
    TFile* infile;
    TFile* fileAnalysis;

    //TTree* treeHeader;
    TTree* tree;
    TTree* treeBeam;
    TTree* treeAnalysis1;
    TTree* treeAnalysis2;
    
    Int_t projA, projZ, targA, targZ;
    Float_t massProj, massTarget, massLight, massHeavy, qValue;

    Float_t energyKinProj;  
    Float_t beamX, beamY, beamZ;
    Float_t beamTheta, beamPhi;
    Double_t genLightEnergy, genLightTheta, genLightPhi;
    Float_t genExcEn;


    Int_t           eventNumber;
    Int_t           detHit[maxDetectors]; // bool: 0 if no hit, 1 if hit in detector []
    Double_t        detEnergyLoss[maxDetectors];
    Double_t        detEnergyLossNotSmeared[maxDetectors];
    Int_t           detStripX[maxDetectors];
    Int_t           detStripY[maxDetectors];
    Double_t        recoPosX[maxDetectors];
    Double_t        recoPosY[maxDetectors];
    Double_t        recoPosZ[maxDetectors];
    // for bug fixing:
    Int_t           firstDetID; // find out which detector fired first
    Double_t        FIx; // first interaction point
    Double_t        FIy;
    Double_t        FIz;
    Int_t           FIdetID; //from sim file
    //Int_t           detHitID[maxDetectors]={-1};

    // define output of analysis
    Int_t detHitMul; // detector hit multiplicity
    Int_t detHitMul3; // aux, count events with more than 2 detectors fired
    Double_t gammaProj, gammaLight;
    Double_t energyTotProj, energyTotLight;
    //  Double_t energyKinHeavy=0.0, energyTotHeavy=0.0, gammaHeavy=0.0, momentumHeavy=0.0;
    Double_t energyKinLight; // is sum of all energy losses
    Double_t momentumProj, momentumLight;
    Double_t simDetectorHitPos[3]; // x, y, z; position used for analysis
    Double_t thetaLightLab, thetaLightCM, phiLight;
    Double_t miss;



};

#endif

