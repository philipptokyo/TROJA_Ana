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
  TApplication *theApp=new TApplication("theApp",  0, 0);
  
  
  cout << "Welcome to troja analyzer" << endl;
  
  
  cout << "Input arguments: " << endl;
  for (int i = 1; i < argc; i++) {
  	cout << argv[i] << endl;
  }
  
  if(argc!=2){
  	cout << "Please give 1 input argument: root file with simulated data, e.g. './analysis sim.root' " << endl; 
  	return 0;
  }
  
  
  
  // define some constants
  
  // nuclear masses in MeV/u
  /// masses should be defined in an external file!
  // const Float_t massTarget = 1875.628; // deuteron
  const Float_t massProj = 122855.922;   // 132Sn
  const Float_t massLight = 938.279;     // light ejectile, proton
  const Float_t massHeavy = 123793.125;  // heavy ejectile, 133Sn
  
  // beam energy
  // should be written from root file
  Float_t energyKinProj = 10.0*132.0;  // 10 MeV/u
 
  
 
 
 
  
  TFile* infile = TFile::Open(argv[1],"read");
  
  if(!infile){
   cout << "Rootfile not found!" << endl;
   return 0;
  }
 
  TTree* tree=(TTree*)infile->Get("troja");
  //TTree* tree=(TTree*)infile->Get("events"); //simulation input
  if(!tree){
    cout << "TTree 'troja' not found!" << endl;
    //cout << "TTree 'events' not found!" << endl;
    return 0;
  }

  
  
  // at the moment:
  //  only 132Sn (d,p) is implemented
  // 132Sn momentum/energy is fix
  Float_t gammaProj = (energyKinProj)/massProj + 1.0;
  //Float_t betaProj = TMath::Sqrt(1.0-(1.0/(gammaProj*gammaProj))); //just for cross checking
  Float_t momentumProj = massProj*TMath::Sqrt(gammaProj*gammaProj-1.0); 
  Float_t energyTotProj = massProj*gammaProj; //total energy
  
  //cout << "Incoming momentum " << momentumProj << ", kin energy " << energyKinProj << ", total energy " << energyTotProj << ", beta " << betaProj << ", gamma " << gammaProj << endl;
   
  TVector3 vProj(0.0, 0.0, momentumProj); 
  TLorentzVector lProj;
  lProj.SetPxPyPzE(0.0, 0.0, momentumProj, energyTotProj);
 
 


  
  // Declaration of leaf types
//  Int_t        eventNumber; //simulation input
//  Double_t        energy;   //simulation input
  Double_t        eventNumber;
  Double_t        energyLoss;
  Double_t        energyTotal;
  Double_t        theta;
  Double_t        phi;
  
//  tree->SetBranchAddress("energy", &energy); //simulation input
  tree->SetBranchAddress("eventNumber", &eventNumber);
  tree->SetBranchAddress("energyLoss", &energyLoss);
  tree->SetBranchAddress("energyTotal", &energyTotal);
  tree->SetBranchAddress("theta", &theta);
  tree->SetBranchAddress("phi", &phi); 



  //TH1F* hMiss=new TH1F("hMiss", "Missing Mass", 1000, -20.0, 20.0);
  TH1F* hMiss=new TH1F("hMiss", "Missing Mass", 2000, -10.0, 10.0);
  //TH2F* hMissTheta=new TH2F("hMissTheta", "Missing mass vs. theta proton", 360,0,180,1000,-20,20);

  
  Int_t nevents=tree->GetEntries();
  cout << "Number of entries found in tree: " << nevents << endl;
  
  for(Int_t e=0; e<nevents; e++){
    tree->GetEvent(e);

    // take only events in backward direction
    // in very forward direction are usually punch-troughs of protons through the detector

    if(theta<90.){
      continue;
    }

    theta/=180.0/TMath::Pi(); //this is important if simulation output is used

    // get total energy and momentum of the light ejectile
    Float_t energyKinLight = energyLoss + energyTotal; //kinetic energy of proton
    //Float_t energyKinLight = energy; //kinetic energy of proton //simulation input
    Float_t gammaLight = energyKinLight/massLight+1.0;
    Float_t energyTotLight = gammaLight*massLight;
    Float_t momentumLight = massLight*TMath::Sqrt(gammaLight*gammaLight-1.0);

    //cout << "Light ejectile momentum " << momentumLight << ", kin energy " << energyKinLight << ", gamma " << gammaLight << ", theta " << theta*180.0/TMath::Pi() << ", phi " << phi << endl;

    TVector3 vLight(0.0, 0.0, 0.0); //momentum of proton
    vLight.SetMagThetaPhi(momentumLight, theta, phi);

    TLorentzVector lLight;
    lLight.SetPxPyPzE(vLight.X(), vLight.Y(), vLight.Z(), energyTotLight);
    //cout << "Light ejectile momentum " << lLight.Pt() << ", kin energy " << lLight.E()-massLight << ", total energy " << lLight.E() << endl;



    //calculate momentum vector of heavy ejectile
    TVector3 vHeavy=vProj;
    vHeavy-=vLight;

    Float_t momentumHeavy = vHeavy.Mag();
    
    Float_t gammaHeavy = TMath::Sqrt((momentumHeavy*momentumHeavy)/(massHeavy*massHeavy)+1.0);
    //Float_t energyKinHeavy = (gammaHeavy-1.0)*massHeavy;     

    //generate the Lorentz vector of the outgoing heavy particle
    Float_t energyTotHeavy = gammaHeavy*massHeavy;

    TLorentzVector lHeavy;
    lHeavy.SetPxPyPzE(vHeavy.X(), vHeavy.Y(), vHeavy.Z(), energyTotHeavy);

    //cout << "Heavy ejectile momentum " << momentumHeavy << ", kin energy " << energyTotHeavy-massHeavy << ", total energ " << energyTotHeavy << endl;
    //cout << endl; 
    


    //get the excitation energy from the missing mass
    //TLorentzVector lMiss = lHeavy - lLight - lProj;
    
    //cout << "Missing mass " << lMiss.E() << endl;
    
    Float_t miss = (lHeavy.E()-massHeavy) + (lLight.E()-massLight) - (lProj.E()-massProj) ;

    //cout << "Missing mass " << miss << endl;

   
    hMiss->Fill(miss);
    //hMissTheta->Fill(vLight.Theta()*180.0/TMath::Pi(), miss);
  
  }

  cout << nevents << " analyzed" << endl;
 
  hMiss->Draw();
  //hMissTheta->Draw();
  
   
  // if histograms shall be plotted, run theApp
  // otherwise program closes
  theApp->Run();
  //delete theApp;  
  
  return 0;
}


