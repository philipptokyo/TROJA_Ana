#ifndef InputInfo_h
#define InputInfo_h

#include "LibPerso.h"

class InputInfo
{
	public:
	
        //todo: make proper getter and setter

        Int_t fNumberEvents;
	Float_t fBeamEnergy;
	
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
