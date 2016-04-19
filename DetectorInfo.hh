#ifndef DetectorInfo_h
#define DetectorInfo_h 1

//#include "globals.hh"

//#include "/home/apollo/geant4.10.00.p04/source/externals/clhep/include/CLHEP/Units/SystemOfUnits.h"
#include "LibPerso.h"

//#define maxDetectors 10
#include "/home/philipp/sim/troja/include/DetectorGlobals.hh"

using namespace std;

typedef struct _geo
{

  string name;
  string type;

  Double_t center[3];
  Double_t rotation[3]; // rotation angle 
  Double_t size[5]; 
  // size for box: x, y, z direction (detector coordinate system)
  // size for tube: rmin, rmax, dz, phistart, phiD

  Int_t noStrips[2]; // in x and y

  Int_t resType; // type of detector resolution: 0=no smearing; 1=gaus
  //const Int_t resNoPars[2] = { 0, 1 };
  static const Int_t resNoPars[2]; // type 1: 1 parameter // initialized in DetectorInfo.cc // index is 'resType'
  Double_t resPar[1]; // type 1: sigma of the gausian in MeV // index is from 0 to 'resNoPars'

} geo;



typedef struct _dat
{
  Int_t eventNumber;

  Double_t fIX, fIY, fIZ; // first interaction point, x, y, z

  //Int_t noOfDet[maxDetectors];
  Int_t fIDetID; // ID of detector with first interaction point
  
  Int_t hitMul; // number of hits, i.e. number of detectors with energy deposition
  Int_t haveHit[maxDetectors]; // flag if detector was hit by particle
  Int_t haveHitID[maxDetectors]; // aux, for quicker debugging... 
  
  Double_t energy[maxDetectors];
  Double_t energyNotSmeared[maxDetectors];
  Int_t stripX[maxDetectors];
  Int_t stripY[maxDetectors];

  Double_t hitPositionX[maxDetectors]; // x, y, z in cartesian coordinates, origin: center of target
  Double_t hitPositionY[maxDetectors]; // x, y, z in cartesian coordinates, origin: center of target
  Double_t hitPositionZ[maxDetectors]; // x, y, z in cartesian coordinates, origin: center of target

} dat;




class DetectorInfo
{
  public:
    DetectorInfo();
    virtual ~DetectorInfo();

    void Parse(string filename);    
    void CheckInput();

    void ResetData();
    void ClearGeometry();
    
    void CalcStripNumbers(Int_t detID, Double_t hx, Double_t hy, Double_t hz, Int_t &stripx, Int_t &stripy);
    void CalcHitPosition(Int_t detID, Int_t stripx, Int_t stripy);
    void CalcHitPosition(Int_t detID, Int_t stripx, Int_t stripy, Double_t &hx, Double_t &hy, Double_t &hz);

    Bool_t IsInFront(Int_t id1, Int_t id2); // checks if detector 1 is in front of detector 2

    Double_t SmearEnergy(Int_t detID, Double_t energy);
        

    // getter
   Int_t GetMaxNoDetectors() const { return (Int_t)maxDetectors; }
    Int_t GetNoOfDetectors() const { return fNoOfDet; }
    
      Double_t GetCenterX(Int_t d) const { return detGeo[d].center[0] ;}
      Double_t GetCenterY(Int_t d) const { return detGeo[d].center[1] ;}
      Double_t GetCenterZ(Int_t d) const { return detGeo[d].center[2] ;}

    Double_t GetRotationX(Int_t d) const { return detGeo[d].rotation[0] ;}
    Double_t GetRotationY(Int_t d) const { return detGeo[d].rotation[1] ;}
    Double_t GetRotationZ(Int_t d) const { return detGeo[d].rotation[2] ;}

        Double_t GetSize0(Int_t d) const { return detGeo[d].size[0] ;}
        Double_t GetSize1(Int_t d) const { return detGeo[d].size[1] ;}
        Double_t GetSize2(Int_t d) const { return detGeo[d].size[2] ;}
        Double_t GetSize3(Int_t d) const { return detGeo[d].size[3] ;}
        Double_t GetSize4(Int_t d) const { return detGeo[d].size[4] ;}

       Int_t GetNoStripsX(Int_t d) const { return detGeo[d].noStrips[0]; }
       Int_t GetNoStripsY(Int_t d) const { return detGeo[d].noStrips[1]; }
       
//         G4bool HaveHit(Int_t d) const { return detGeo[d].haveHit; }

                 string GetName(Int_t d) { return detGeo[d].name; }
                 string GetType(Int_t d) { return detGeo[d].type; }

         Int_t GetResType(Int_t d) const { return detGeo[d].resType; }
         Int_t GetResNoPars(Int_t d) const { return detGeo[d].resNoPars[detGeo[d].resType]; }
Double_t GetResPar(Int_t d, Int_t p) const { return detGeo[d].resPar[p] ;}
    

    TRandom3* GetRandomizer() { return fRandomizer; }

    TRotation* GetRotationMatrix(Int_t d);
    TRotation GetInverseRotationMatrix(Int_t d);


    // setter
    
    void SetNoOfDetectors(Int_t n) { fNoOfDet=n; }
    
      void SetCenterX(Int_t d, Double_t c) { detGeo[d].center[0]=c ;}
      void SetCenterY(Int_t d, Double_t c) { detGeo[d].center[1]=c ;}
      void SetCenterZ(Int_t d, Double_t c) { detGeo[d].center[2]=c ;}

    void SetRotationX(Int_t d, Double_t r) { detGeo[d].rotation[0]=r;}
    void SetRotationY(Int_t d, Double_t r) { detGeo[d].rotation[1]=r;}
    void SetRotationZ(Int_t d, Double_t r) { detGeo[d].rotation[2]=r;}

        void SetSize0(Int_t d, Double_t s) { detGeo[d].size[0]=s ;}
        void SetSize1(Int_t d, Double_t s) { detGeo[d].size[1]=s ;}
        void SetSize2(Int_t d, Double_t s) { detGeo[d].size[2]=s ;}
        void SetSize3(Int_t d, Double_t s) { detGeo[d].size[3]=s ;}
        void SetSize4(Int_t d, Double_t s) { detGeo[d].size[4]=s ;}

    void SetNoStripsX(Int_t d, Int_t n)    { detGeo[d].noStrips[0]=n; }
    void SetNoStripsY(Int_t d, Int_t n)    { detGeo[d].noStrips[1]=n; }

//    void SetHaveHit(Int_t d, G4bool h)   { detGeo[d].haveHit=h; }

           void SetName(Int_t d, string n) { detGeo[d].name=n; }
           void SetType(Int_t d, string n) { detGeo[d].type=n; }
   
         void SetResType(Int_t d, Int_t t) { detGeo[d].resType=t; }
void SetResPar(Int_t d, Int_t p, Double_t v) { detGeo[d].resPar[p]=v ;}
   
   
   
   
   
   
   
    
    dat detData; // todo: this should be private, ne


  private:

    TRandom3* fRandomizer;

    geo detGeo[maxDetectors];

    Int_t fNoOfDet; // number of detectors


};



#endif
