/**
 * 
 * main.cc for analysis of simulated transfer reactions
 * 
 * 
 **/





#include "main.h"

#include "LibPerso.h"

using namespace std;



Int_t main(Int_t argc, char **argv){
  
  
  
  
  
  
  // if root shall stop the program before it finished, comment this in
//  TApplication *theApp=new TApplication("theApp",  0, 0);
  
  
  cout << "Welcome to troja analyzer" << endl;
  
  
  cout << "Input arguments: " << endl;
  for (int i = 1; i < argc; i++) {
  	cout << argv[i] << endl;
  }
  
  if(argc!=2){
  	cout << "Please give 1 input argument: root file with simulated data, e.g. './analysis sim.root' " << endl; 
  	return 0;
  }
  
  
  
  TFile* infile = TFile::Open(argv[1],"read");
  
  if(!infile){
    
   cout << "Rootfile not found!" << endl;
   return 0;
  
  }
 
 

  
 
 
  
  // if histograms shall be plotted, run theApp
  // otherwise program closes
//  theApp->Run();
  
  
  return 0;
}


