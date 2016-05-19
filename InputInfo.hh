#ifndef InputInfo_hh
#define InputInfo_hh

#include "LibPerso.h"

#define maxNumberOfStates 10

class InputInfo
{
	public:
	
        //todo: private variables plus getter and setter

        // file names
        char fOutFileNameReaction[500];
        char fOutFileNameFresco[500];
	char fOutFileNameMakeEvents[500];
	char fOutFileNameTroja[500];
	char fInFileNameGeometry[500];
        char fOutFileNameAnalysis[500];
        char fOedoSimFileName[500];
        
        // general, for all programs
        Int_t fNumberEvents;

        // for reaction and event generator
        Int_t fProjA, fProjZ, fTargetA, fTargetZ, fLightA, fLightZ;
        Int_t fNumberOfStates;
        //Float_t fMaxExEnergy;
        Float_t fStateEnergy[maxNumberOfStates+1]; // index 0 is elastic
	Float_t fBeamEnergy;
        Float_t fAngleMin, fAngleMax; // minimum and maximum angle in CM
        Float_t fElasticDownscale;

        // event generator related
        Bool_t HaveOedoSimFileName(){return fHaveOedoSimFileName;};
        Bool_t HaveFrescoFileName(){return fHaveFrescoFileName;};
        void UnsetFrescoFileName(){fHaveFrescoFileName=false;};
        Bool_t FrescoHeaderOnly(){return fFrescoHeaderOnly;};
        Bool_t ProfileBeamE(){return fProfileE;};
        Bool_t ProfileBeamX(){return fProfileX;};
        Bool_t ProfileBeamY(){return fProfileY;};
        Bool_t ProfileBeamA(){return fProfileA;};
        Bool_t ProfileBeamB(){return fProfileB;};
        Bool_t NoBeamTracking(){return fNoBeamTracking;};
        Bool_t IncludeElastic(){return fIncludeElastic;}
	
        //analysis related
        Float_t fResTargetX, fResTargetY, fResTargetZ; // should be obtained from beam tracking detector(s)
        Float_t fResTargetA, fResTargetB; // should be obtained from beam tracking detector(s)
        //Float_t fResDet1X, fResDet1Y, fResDet1Z; // this should be given in detector coordinates
        //Float_t fResDet1E; // energy (loss) resolution in first detector
        //Float_t fResDet2E; // energy (loss) resolution in second detector
        Float_t fResBeamE;
	

        //
        void parse(char filename[100]);
	
	InputInfo();
	~InputInfo();
	
        
        
        private:
       
        
        Bool_t fHaveOedoSimFileName; 
        Bool_t fHaveFrescoFileName;
        Bool_t fFrescoHeaderOnly;
        Bool_t fProfileE, fProfileX, fProfileY, fProfileA, fProfileB;
        Bool_t fNoBeamTracking;
        Bool_t fIncludeElastic; 
        	
	
};

#endif
