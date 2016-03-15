#ifndef InputInfo_hh
#define InputInfo_hh

#include "LibPerso.h"

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
        Float_t fMaxExEnergy;
	Float_t fBeamEnergy;

        // event generator related
        Bool_t HaveOedoSimFileName(){return fHaveOedoSimFileName;};
        Bool_t ProfileBeamE(){return fProfileE;};
        Bool_t ProfileBeamX(){return fProfileX;};
        Bool_t ProfileBeamY(){return fProfileY;};
        Bool_t ProfileBeamA(){return fProfileA;};
        Bool_t ProfileBeamB(){return fProfileB;};
	
        //analysis related
        Float_t fResTargetX, fResTargetY, fResTargetZ; // should be obtained from beam tracking detector(s)
        Float_t fResTargetA, fResTargetB; // should be obtained from beam tracking detector(s)
        Float_t fResDet1X, fResDet1Y, fResDet1Z; // this should be given in detector coordinates
        Float_t fResDet1E; // energy (loss) resolution in first detector
        Float_t fResDet2E; // energy (loss) resolution in second detector
        Float_t fResBeamE;
	

        //
        void parse(char filename[100]);
	
	InputInfo();
	~InputInfo();
	
        
        
        private:
       
        
        Bool_t fHaveOedoSimFileName; 
        Bool_t fProfileE, fProfileX, fProfileY, fProfileA, fProfileB; 
        	
	
};

#endif
