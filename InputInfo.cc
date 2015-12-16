
#include "InputInfo.h"

using namespace std;

InputInfo::InputInfo(){
	
	fNumberEvents=100;
	fBeamEnergy=10.0; // in MeV/u

        fHaveOedoSimFileName=false;
        
        fProfileE=false;	
        fProfileX=false;	
        fProfileY=false;	
        fProfileA=false;	
        fProfileB=false;	
}

InputInfo::~InputInfo(){
	
}

void InputInfo::parse(char filename[100]){
	
	cout << endl;
	cout << "Parsing of input file" << endl;
	
	
	//check if in/out file exists
	ifstream fin;
	ifstream FileTest(filename);
	if(!FileTest){
		cout << "Can not open file '" << filename << "'! Does it exist?" << endl;
		abort();
	}else{
		cout << "Opening file '" << filename << "'" << endl;
		fin.open(filename);
	}
	
	
	//parse input file line-wise
	string line;
	Int_t counter=0;
	const Int_t stopper=1000; 
	
	const Int_t maxArg=2;
	char temp[maxArg][500];
	
	
	
	while(fin.good())  {
		
		
		//reset values
		for(Int_t i=0; i< maxArg; i++){temp[i][0]='\0';}
		
		
		getline(fin,line);
		//printf("Parsing line '%s'\n",line.c_str());
		
		
		
		//skip comments
		line.copy(temp[0],2,0); //first two characters
		temp[0][2]='\0';
		//printf("temp1=%s\n",temp1);
		if(strcmp(temp[0],"//")==0){continue;}
		
		//parse the line
		//write each word to a char
		//delimiter is whitespace
		std::istringstream iss(line);
		for(Int_t i=0; i<maxArg; i++){iss >> temp[i];}
		
		//skip empty lines
		if(strcmp(temp[0],"")==0){continue;}
		
		
		
		//get input parameter:
		
		
		if(strcmp(temp[0],"output_rootfile_reaction")==0)  {
			strcpy(fOutfilenameReaction,temp[1]);
			cout << "Output file of reactions is '" << fOutfilenameReaction << "'" << endl;
		}
		else if(strcmp(temp[0],"output_rootfile_makeEvents")==0){
			strcpy(fOutfilenameMakeEvents,temp[1]);
			cout << "Output file of makeEvents is '" << fOutfilenameMakeEvents << "'" << endl;
		}
		else if(strcmp(temp[0],"beam_energy")==0){
			fBeamEnergy=atof(temp[1]);
			cout << "Beam energy is set to '" << fBeamEnergy << "' MeV/u" << endl;
		}
		else if(strcmp(temp[0],"number_events")==0){
			fNumberEvents=atoi(temp[1]);
			cout << "Number of events to generate is set to " << fNumberEvents << endl;
		}
		else if(strcmp(temp[0],"beam_profile_file_oedo")==0){
			strcpy(fOedoSimFileName,temp[1]);
                        fHaveOedoSimFileName=true;
			cout << "Root file name with OEDO beam profile is set to " << fOedoSimFileName << endl;
		}
		else if(strcmp(temp[0],"beam_profile_energy")==0){
                        fProfileE=true;
			cout << "Beam profile - energy - requested " << endl;
		}
		else if(strcmp(temp[0],"beam_profile_position_x")==0){
                        fProfileX=true;
			cout << "Beam profile - x position - requested " << endl;
		}
		else if(strcmp(temp[0],"beam_profile_position_y")==0){
                        fProfileY=true;
			cout << "Beam profile - y position - requested " << endl;
		}
		else if(strcmp(temp[0],"beam_profile_angular_a")==0){
                        fProfileA=true;
			cout << "Beam profile - a angle (dx) - requested " << endl;
		}
		else if(strcmp(temp[0],"beam_profile_angular_b")==0){
                        fProfileB=true;
			cout << "Beam profile - b angle (dy) - requested " << endl;
		}
		else if(strcmp(temp[0],"output_rootfile_troja")==0){
			strcpy(fOutfilenameTroja,temp[1]);
			cout << "Output file of troja is '" << fOutfilenameTroja << "'" << endl;
		}
		else if(strcmp(temp[0],"resolution_target_x")==0){
			fResTargetX=atof(temp[1]);
			cout << "Resolution of target X position is set to '" << fResTargetX << "' mm (sigma of a Gaussian)" << endl;
		}
		else if(strcmp(temp[0],"resolution_target_y")==0){
			fResTargetY=atof(temp[1]);
			cout << "Resolution of target Y position is set to '" << fResTargetY << "' mm (sigma of a Gaussian)" << endl;
		}
		else if(strcmp(temp[0],"resolution_target_z")==0){
			fResTargetZ=atof(temp[1]);
			cout << "Resolution of target Z position is set to '" << fResTargetZ << "' mm (sigma of a Gaussian)" << endl;
		}
		else if(strcmp(temp[0],"resolution_target_a")==0){
			fResTargetA=atof(temp[1]);
			cout << "Resolution of target A angle is set to '" << fResTargetA << "' mrad (sigma of a Gaussian)" << endl;
		}
		else if(strcmp(temp[0],"resolution_target_b")==0){
			fResTargetB=atof(temp[1]);
			cout << "Resolution of target B angle is set to '" << fResTargetB << "' mrad (sigma of a Gaussian)" << endl;
		}
		else if(strcmp(temp[0],"resolution_detector1_x")==0){
			fResDet1X=atof(temp[1]);
			cout << "Resolution of detector 1 X-position is set to '" << fResDet1X << "' mm (sigma of a Gaussian)" << endl;
		}
		else if(strcmp(temp[0],"resolution_detector1_y")==0){
			fResDet1Y=atof(temp[1]);
			cout << "Resolution of detector 1 Y-position is set to '" << fResDet1Y << "' mm (sigma of a Gaussian)" << endl;
		}
		else if(strcmp(temp[0],"resolution_detector1_z")==0){
			fResDet1Z=atof(temp[1]);
			cout << "Resolution of detector 1 Z-position is set to '" << fResDet1Z << "' mm (sigma of a Gaussian)" << endl;
		}
		else if(strcmp(temp[0],"resolution_detector1_e")==0){
			fResDet1E=atof(temp[1]);
			cout << "Resolution of detector 1 energy is set to '" << fResDet1E << "' MeV (sigma of a Gaussian)" << endl;
		}
		else if(strcmp(temp[0],"resolution_detector2_e")==0){
			fResDet2E=atof(temp[1]);
			cout << "Resolution of detector 2 energy is set to '" << fResDet2E << "' MeV (sigma of a Gaussian)" << endl;
		}
		else {
			cout<<"Could not read your input keyword '" << temp[0] << "'. Aborting program."<<endl; 
			abort();
		}
		
		
		//count lines
		counter++;
		if(counter>stopper){
			cout << "Reached " << counter << " lines in file '" << filename << "'! Stopping!" << endl;
			abort();
		}
	} // end of reading input file
	
	cout << endl;
	
	
}
