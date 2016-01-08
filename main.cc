/**
 * 
 * main.cc for analysis of simulated transfer reactions
 * 
 * 
 **/


#include "main.hh"

#include "LibPerso.h"

#include "InputInfo.hh"

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
  	cout << "Please give 1 input argument: text file with input information, e.g. './analysis input.txt' " << endl; 
  	return 0;
  }

  InputInfo* info=new InputInfo();
  info->parse(argv[1]);
  
  TRandom3* randomizer = new TRandom3();
  randomizer->SetSeed(0); 
  
  // define some constants
  
  // nuclear masses in MeV/u
  /// masses should be defined in an external file!
  const Float_t massProj = 122855.922;   // 132Sn
  const Float_t massLight = 938.279;     // light ejectile, proton
//  const Float_t massLight = 1875.628;     // light ejectile, deuteron
  const Float_t massHeavy = 123793.125;  // heavy ejectile, 133Sn
  
 
 
 
  
  TFile* fileBeam = TFile::Open(info->fOutFileNameMakeEvents,"read");
  
  if(!fileBeam){
   cout << "1st root file not found!" << endl;
   return 0;
  }
 
  TTree* treeBeam=(TTree*)fileBeam->Get("events");
  //TTree* tree=(TTree*)infile->Get("events"); //simulation input
  if(!treeBeam){
    cout << "TTree 'events' not found in 1st root file!" << endl;
    //cout << "TTree 'events' not found!" << endl;
    return 0;
  }

  
  
  TFile* infile = TFile::Open(info->fOutFileNameTroja,"read");
  
  if(!infile){
   cout << "2nd root file not found!" << endl;
   return 0;
  }
 
  TTree* tree=(TTree*)infile->Get("troja");
  //TTree* tree=(TTree*)infile->Get("events"); //simulation input
  if(!tree){
    cout << "TTree 'troja' not found in 2nd root file!" << endl;
    //cout << "TTree 'events' not found!" << endl;
    return 0;
  }



  // Get tree with light particle simulation

  // Declaration of leaf types
//  Int_t        eventNumber; //simulation input
//  Double_t        energy;   //simulation input
  Double_t        eventNumber=0.0;
  Double_t        energyLoss=0.0;
  Double_t        energyTotal=0.0;
  Double_t        x=0.0, y=0.0, z=0.0;
  Double_t        theta=0.0; //should not be used, doesn't include vertex
  Double_t        phi=0.0;   //should not be used, doesn't include vertex
  
//  tree->SetBranchAddress("energy", &energy); //simulation input
  tree->SetBranchAddress("eventNumber", &eventNumber);
  tree->SetBranchAddress("energyLoss", &energyLoss);
  tree->SetBranchAddress("energyTotal", &energyTotal);
  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("z", &z);
  tree->SetBranchAddress("theta", &theta);
  tree->SetBranchAddress("phi", &phi); 



  // Initialize tree with projectile information
  //Int_t           beamMassNumber;
  //Int_t           beamChargeNumber;
  //Float_t         beamEnergy;
  Float_t         energyKinProj = 10.0*132.0;  // 10 MeV/u
  Float_t         beamX;
  Float_t         beamY;
  Float_t         beamZ;
  Float_t         beamTheta;
  Float_t         beamPhi;

  //treeBeam->SetBranchAddress("beamMassNumber", &beamMassNumber);     //needs proper implementation
  //treeBeam->SetBranchAddress("beamChargeNumber", &beamChargeNumber); //needs proper implementation
  //treeBeam->SetBranchAddress("beamEnergy", &beamEnergy);             //needs proper implementation
  treeBeam->SetBranchAddress("beamEnergy", &energyKinProj); // is in MeV/u
  treeBeam->SetBranchAddress("beamX", &beamX);
  treeBeam->SetBranchAddress("beamY", &beamY);
  treeBeam->SetBranchAddress("beamZ", &beamZ);
  treeBeam->SetBranchAddress("beamTheta", &beamTheta);
  treeBeam->SetBranchAddress("beamPhi", &beamPhi);








  //TH1F* hMiss=new TH1F("hMiss", "Missing Mass", 1000, -20.0, 20.0);
  //TH1F* hMiss=new TH1F("hMiss", "Missing Mass", 2000, -10.0, 10.0);
  TH1F* hMiss=new TH1F("hMiss", "Missing Mass", 2000, -6.0, 1.0);
  //TH2F* hMissTheta=new TH2F("hMissTheta", "Missing mass vs. theta proton", 360,0,180,1000,-20,20);






  
  Int_t nevents=tree->GetEntries();
  cout << "Number of entries found in tree: " << nevents << endl;

  Int_t neventsBeam=treeBeam->GetEntries();
  if(nevents!=neventsBeam){
    if(nevents<neventsBeam){
      cout << "INFO: " << nevents << " simulated and " << neventsBeam << " found in simulation input. Analyzing only " << nevents << "." << endl;
    }
    if(nevents>neventsBeam){
      cout << "ERROR: more events simulated (" << nevents << ") than found in simulation input file (" << neventsBeam << ")! Check this!" << endl;
      return 0;
    }
  }
  
  for(Int_t e=0; e<nevents; e++){
    tree->GetEvent(e);
    treeBeam->GetEvent(e); // todo: make tree friend instead
    

    // smear out data with detector position resolutions
    beamX=randomizer->Gaus(beamX, info->fResTargetX);
    beamY=randomizer->Gaus(beamY, info->fResTargetY);
    beamZ=randomizer->Gaus(beamZ, info->fResTargetZ);

    x=randomizer->Gaus(x, info->fResDet1X); 
    y=randomizer->Gaus(y, info->fResDet1Y); 
    z=randomizer->Gaus(z, info->fResDet1Z);
    
    // smear out data with detector energy resolutions
    energyLoss=randomizer->Gaus(energyLoss, info->fResDet1E);  
    energyTotal=randomizer->Gaus(energyTotal, info->fResDet2E);  

    // smearing with angles of the incoming beam is done below






    // projectile kinematics

    energyKinProj*=132.0; // todo: should be massProj
    // Projectile data
    // at the moment from simulation input
    // todo: separate simulation including incoming tracking

    energyKinProj = randomizer->Gaus(energyKinProj, info->fResBeamE);
      
    //  only 132Sn (d,p) is implemented
    // 132Sn momentum/energy is fix
    Float_t gammaProj = (energyKinProj)/massProj + 1.0;
    //Float_t betaProj = TMath::Sqrt(1.0-(1.0/(gammaProj*gammaProj))); //just for cross checking
    Float_t momentumProj = massProj*TMath::Sqrt(gammaProj*gammaProj-1.0); 
    Float_t energyTotProj = massProj*gammaProj; //total energy
    
     
    TVector3 vProj(0.0, 0.0, momentumProj); 
    vProj.SetMagThetaPhi(momentumProj, beamTheta, beamPhi); // comment out this line to see the effect of no beam profile correction
    
    // rotate by beam angular resolution
    vProj.RotateY(randomizer->Gaus(0.0, (info->fResTargetA)/1000.0));
    vProj.RotateX(randomizer->Gaus(0.0, (info->fResTargetB)/1000.0));

    TLorentzVector lProj;
    lProj.SetPxPyPzE(vProj.X(), vProj.Y(), vProj.Z(), energyTotProj);





    // light ejectile kinematics

    // take only events in backward direction
    // in very forward direction are usually punch-troughs of protons through the detector
    if(theta<90.0){
      continue;
    }

    theta/=180.0/TMath::Pi(); //this is important if simulation output is used

    // get total energy and momentum of the light ejectile
    Float_t energyKinLight = energyLoss + energyTotal; //kinetic energy of proton
    //Float_t energyKinLight = energy; //kinetic energy of proton //simulation input
    Float_t gammaLight = energyKinLight/massLight+1.0;
    Float_t energyTotLight = gammaLight*massLight;
    Float_t momentumLight = massLight*TMath::Sqrt(gammaLight*gammaLight-1.0);


    TVector3 vLight(x-beamX, y-beamY, z-beamZ); //momentum of proton
    vLight.SetMag(momentumLight);
    //vLight.SetMagThetaPhi(momentumLight, theta, phi);

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


