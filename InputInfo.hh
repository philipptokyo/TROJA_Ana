#ifndef InputInfo_hh
#define InputInfo_hh

#include "LibPerso.h"

#define maxNumberOfStates 20
#define maxCutType 3

using namespace std;

typedef struct _tar
{
  
  string material;
  Double_t density;
  Double_t densityOffset;
  
  Double_t center[3];
  Double_t size[3];

} tar;



class InputInfo
{
  public:
  
  //todo: private variables plus getter and setter
  
  // file names
  char fOutFileNameReaction[500];
  char fOutFileNameFresco[500];
  char fOutFileNamePace[500];
  char fOutFileNameMakeEvents[500];
  char fOutFileNameTroja[500];
  char fInFileNameGeometry[500];
  char fOutFileNameAnalysis[500];
  char fOedoSimFileName[500];
  char fSourceTypeName[50]; // alpha, beta-, beta+, gamma
  Int_t fSourceTypeId; //   0      1      2      3
  Int_t fSourceLines;
  Float_t fSourceEnergies[maxNumberOfStates];
  
  char fFileNameCuts[maxCutType][500];
  Int_t fNumberOfCutFiles;
  
  // general, for all programs
  Int_t fNumberEvents;
  
  // for reaction and event generator
  Int_t fProjA, fProjZ, fTargetA, fTargetZ, fLightA, fLightZ;
  Int_t fNumberOfStates;
  //Float_t fMaxExEnergy;
  Float_t fStateEnergy[maxNumberOfStates+1]; // index 0 is elastic
  Int_t   fStateGammaMul[maxNumberOfStates+1]; // [state][energy]
  Float_t fStateGammaEnergies[maxNumberOfStates+1][maxNumberOfStates+1]; // [state][energy]
  Float_t fStateSpecFact[maxNumberOfStates+1]; // spectroscopic factors
  Float_t fBeamEnergy;
  Float_t fAngleMin, fAngleMax; // minimum and maximum angle in CM
  Float_t fElasticDownscale;
  
  // event generator related
  Bool_t HaveOedoSimFileName(){return fHaveOedoSimFileName;};
  Bool_t HaveFrescoFileName(){return fHaveFrescoFileName;};
  void UnsetFrescoFileName(){fHaveFrescoFileName=false;};
  Bool_t FrescoHeaderOnly(){return fFrescoHeaderOnly;};
  Bool_t HavePaceFileName(){return fHavePaceFileName;};
  Bool_t PaceOnly(){return fPaceOnly;};
  Bool_t ProfileBeamE(){return fProfileE;};
  Bool_t ProfileBeamX(){return fProfileX;};
  Bool_t ProfileBeamY(){return fProfileY;};
  Bool_t ProfileBeamA(){return fProfileA;};
  Bool_t ProfileBeamB(){return fProfileB;};
  Bool_t NoBeamTracking(){return fNoBeamTracking;};
  Bool_t IncludeElastic(){return fIncludeElastic;};
  Bool_t Source(){return fSource;};
  Bool_t AddGammas(){return fAddGammas;}
  
  //analysis related
  Float_t fResTargetX, fResTargetY, fResTargetZ; // should be obtained from beam tracking detector(s)
  Float_t fResTargetA, fResTargetB; // should be obtained from beam tracking detector(s)
  //Float_t fResDet1X, fResDet1Y, fResDet1Z; // this should be given in detector coordinates
  //Float_t fResDet1E; // energy (loss) resolution in first detector
  //Float_t fResDet2E; // energy (loss) resolution in second detector
  Float_t fResBeamE;


  void SetTargetPosition(Double_t x, Double_t y, Double_t z) {fTarget.center[0]=x; fTarget.center[1]=y; fTarget.center[2]=z; }
   void SetTargetSize(Double_t x, Double_t y, Double_t z) {fTarget.size[0]=x; fTarget.size[1]=y; fTarget.size[2]=z; }
   void SetTargetMaterial(string mat) {fTarget.material=mat; }
   void SetTargetDensity(Double_t den) {fTarget.density=den;}
   void SetTargetDensityOffset(Double_t deno) {fTarget.densityOffset=deno;}

   Double_t GetTargetPosition(Int_t i) {return fTarget.center[i]; }
       Double_t GetTargetSize(Int_t i) {return fTarget.size[i]; }
           Double_t GetTargetDensity() {return fTarget.density; }
     Double_t GetTargetDensityOffset() {return fTarget.densityOffset; }
            string GetTargetMaterial() {return fTarget.material; }

  
  //
  void parse(char filename[100]);
  
  InputInfo();
  ~InputInfo();
  
  
  
  private:
  
  
  Bool_t fHaveOedoSimFileName; 
  Bool_t fHaveFrescoFileName;
  Bool_t fFrescoHeaderOnly;
  Bool_t fHavePaceFileName;
  Bool_t fPaceOnly;
  Bool_t fProfileE, fProfileX, fProfileY, fProfileA, fProfileB;
  Bool_t fNoBeamTracking;
  Bool_t fIncludeElastic; 
  Bool_t fSource;
  Bool_t fAddGammas;

  tar fTarget;
  
  	
  
};

#endif
