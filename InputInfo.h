#ifndef InputInfo_h
#define InputInfo_h

#include "LibPerso.h"

class InputInfo
{
	public:
	
	Float_t fBeamEnergy;
	
	char fInfilenameFromReaction[200];
	
	char fRootfilename[200];
	
	Int_t fNumberEvents;
        
        char fOedoSimFileName[300];
        
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
