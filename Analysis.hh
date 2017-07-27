#ifndef Analysis_hh
#define Analysis_hh

#include "LibPerso.h"

#include "InputInfo.hh"
#include "DetectorInfo.hh"

#include "Nucleus.hh"
#include "Kinematics.hh"
#include "EnLoss.hh"

// maxCutType defined in InputInfo.hh
#define maxCuts 10 // per cut file 
#define maxGammaMulGen 2

class Analysis
{
  
  public:
    
    Bool_t Init();
    void ResetVariables();
    void CreateHeader();
    void Analysis1();
    Bool_t GetCuts();
    void Analysis2();
    void MissingMass(Int_t channel);
    void AnalysisGamma();

    Analysis(InputInfo *i, DetectorInfo* d);
    
    ~Analysis();


  private:

    void ProcessDetectorHits();
    
    InputInfo* info;
    DetectorInfo* detInfo;

    EnLoss* fEnLoss;

    TFile* fileBeam;
    TFile* infile;
    TFile* fileAnalysis;

    TTree* tree;
    TTree* treeBeam;
    TTree* treeAnaHeader[maxCutType];
    TTree* treeAnalysis1;
    TTree* treeAnalysis2[maxCutType];
    
    Int_t projA, projZ, targA, targZ;
    Int_t recoA, recoZ, ejecA, ejecZ;
    Float_t massProj, massTarget, massLight, massHeavy, qValue;

    Nucleus* nucProj;
    Nucleus* nucTarg;
    Nucleus* nucEjec[maxCutType];
    Nucleus* nucReco[maxCutType];

    Float_t energyKinProj;  
    Float_t beamX, beamY, beamZ;
    Float_t beamA, beamB;
    Float_t beamTheta, beamPhi;
    Double_t genLightEnergy, genLightTheta, genLightThetaCM, genLightPhi;
    Float_t genExcEn;
    Int_t genState;


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

    Int_t type; // reaction type; 0 = elastic scattering, 1 = (d,p) 

    Double_t gammaProj, gammaLight;
    Double_t energyTotProj, energyTotLight;
    //  Double_t energyKinHeavy=0.0, energyTotHeavy=0.0, gammaHeavy=0.0, momentumHeavy=0.0;
    Double_t energyKinLight, energyKinLightUncorr; // is sum of all energy losses
    Double_t momentumProj, momentumLight;
    Double_t simDetectorHitPos[3]; // x, y, z; position used for analysis
    Double_t thetaLightLab, thetaLightCM, phiLight;
    Double_t miss;
    
    //Double_t energyGamma;
    //Double_t energyGammaDC;
    //Int_t grapeDet;

    Double_t targetEnergyLoss;
    
    Int_t genGammaMul;
    Float_t genGammaERest[maxGammaMulGen], genGammaELab[maxGammaMulGen];

    Int_t grapeDetMul;
    Int_t grapeCryMul[grapeMaxDet];
    Int_t grapeSegMul[grapeMaxDet][grapeMaxCry];
    
    Double_t grapeDetEnergy[grapeMaxDet];
    Double_t grapeCryEnergy[grapeMaxDet][grapeMaxCry];
    Double_t grapeSegEnergy[grapeMaxDet][grapeMaxCry][grapeMaxSeg];
    
    Double_t grapeSumEnergyDC;
    Double_t grapeDetEnergyDC[grapeMaxDet];
    Double_t grapeCryEnergyDC[grapeMaxDet][grapeMaxCry];
    Double_t grapeSegEnergyDC[grapeMaxDet][grapeMaxCry][grapeMaxSeg];






    
    Bool_t cutExists[maxCutType][maxCuts];
    TCutG* cut[maxCutType][maxCuts];

    TRandom3* randomizer; 

    // define histograms
    TH1F* hMiss;
    //TH2F* hMissTheta;
  
    TH2F* hdEE_A;
    TH2F* hdEE_B;
    TH2F* hdEE_F;
    TH2F* hEth;
  
    TH1F* hThetaLab;
    TH1F* hThetaCM;

};

#endif

