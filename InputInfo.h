#ifndef InputInfo_h
#define InputInfo_h

#include "LibPerso.h"

class InputInfo
{
	public:
	
        //todo: make proper getter and setter

        Int_t fNumberEvents;
	Float_t fBeamEnergy;
        Float_t fResDet1X, fResDet1Y, fResDet1Z; // this should be given in detector coordinates
        Float_t fResDet1E; // energy (loss) resolution in first detector
        Float_t fResDet2E; // energy (loss) resolution in second detector
        Float_t fResTargetX, fResTargetY, fResTargetZ; // should be obtained from beam tracking detector(s)
        Float_t fResTargetA, fResTargetB; // should be obtained from beam tracking detector(s)
	
        char fOutfilenameReaction[300];
	
	char fOutfilenameMakeEvents[300];
	
	char fOutfilenameTroja[300];
	
        char fOedoSimFileName[300];
        
        //
        Bool_t HaveOedoSimFileName(){return fHaveOedoSimFileName;};
        Bool_t ProfileBeamE(){return fProfileE;};
        Bool_t ProfileBeamX(){return fProfileX;};
        Bool_t ProfileBeamY(){return fProfileY;};
        Bool_t ProfileBeamA(){return fProfileA;};
        Bool_t ProfileBeamB(){return fProfileB;};
	
        void parse(char filename[100]);
	
	InputInfo();
	~InputInfo();
	
        
        
        private:
       
        
        Bool_t fHaveOedoSimFileName; 
        Bool_t fProfileE, fProfileX, fProfileY, fProfileA, fProfileB; 
        	
	
};

#endif
