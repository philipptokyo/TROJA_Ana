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
//  const Float_t massProj = 122855.922;   // 132Sn
//  const Float_t massLight = 938.279;     // light ejectile, proton
////  const Float_t massLight = 1875.628;     // light ejectile, deuteron
//  const Float_t massHeavy = 123793.125;  // heavy ejectile, 133Sn
  
 
 
 
  
  TFile* fileBeam = TFile::Open(info->fOutFileNameMakeEvents,"read");
  
  if(!fileBeam){
   cout << "makeEvents root file not found!" << endl;
   return 0;
  }

  TTree* treeHeader=(TTree*)fileBeam->Get("header");
  if(!treeHeader){
    cout << "TTree 'header' not found in makeEvents root file!" << endl;
    return 0;
  }


  Int_t projA=132;
  //Float_t projMass=0.0, targetMass=0.0, lightMass=0.0, heavyMass=0.0, qValue=0.0;
  Float_t massProj=0.0, massTarget=0.0, massLight=0.0, massHeavy=0.0, qValue=0.0;

  treeHeader->SetBranchAddress("projA", &projA);
  treeHeader->SetBranchAddress("projMass", &massProj);
  treeHeader->SetBranchAddress("targetMass", &massTarget);
  treeHeader->SetBranchAddress("lightMass", &massLight);
  treeHeader->SetBranchAddress("heavyMass", &massHeavy);
  treeHeader->SetBranchAddress("qValue", &qValue);

  treeHeader->GetEntry(0);
  //printf("Obtained masses: projectile %f, target %f, light ejectile %f, heavy ejectile %f; Q-value %f\n", massProj, massTarget, massLight, massHeavy, qValue);





 
  TTree* treeBeam=(TTree*)fileBeam->Get("events");
  //TTree* tree=(TTree*)infile->Get("events"); //simulation input
  if(!treeBeam){
    cout << "TTree 'events' not found in makeEvents root file!" << endl;
    //cout << "TTree 'events' not found!" << endl;
    return 0;
  }

  
  
  TFile* infile = TFile::Open(info->fOutFileNameTroja,"read");
  
  if(!infile){
   cout << "geant (troja) root file not found!" << endl;
   return 0;
  }
 
  TTree* tree=(TTree*)infile->Get("troja");
  //TTree* tree=(TTree*)infile->Get("events"); //simulation input
  if(!tree){
    cout << "TTree 'troja' not found in geant root file!" << endl;
    //cout << "TTree 'events' not found!" << endl;
    return 0;
  }



  // Get tree with light particle simulation

  // Declaration of leaf types
//  Int_t        eventNumber; //simulation input
//  Double_t        energy;   //simulation input
  Double_t        eventNumber=0.0;
  Double_t        energyLoss=0.0;
  Double_t        energySum=0.0;
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
  Float_t         energyKinProj = 10.0*132.0;  // 10 MeV/u
  Float_t         beamX, beamY, beamZ;
  Float_t         beamTheta, beamPhi;
  Double_t        genLightEnergy, genLightTheta, genLightPhi;
  Float_t         genExcEn=0.0;

  treeBeam->SetBranchAddress("lightEnergy", &genLightEnergy);
  treeBeam->SetBranchAddress("lightTheta", &genLightTheta);
  treeBeam->SetBranchAddress("lightPhi", &genLightPhi);
  treeBeam->SetBranchAddress("beamEnergy", &energyKinProj); // is in MeV/u
  treeBeam->SetBranchAddress("beamX", &beamX);
  treeBeam->SetBranchAddress("beamY", &beamY);
  treeBeam->SetBranchAddress("beamZ", &beamZ);
  treeBeam->SetBranchAddress("beamTheta", &beamTheta);
  treeBeam->SetBranchAddress("beamPhi", &beamPhi);
  treeBeam->SetBranchAddress("excitationEnergy", &genExcEn);






  // define output of analysis
  Float_t gammaProj, gammaLight, gammaHeavy;
  Float_t energyTotProj, energyTotLight, energyTotHeavy;
  Float_t momentumProj, momentumLight, momentumHeavy; 
  Float_t miss=0.0;
  
  TFile* fileAnalysis = new TFile(info->fOutFileNameAnalysis, "recreate");
  fileAnalysis->cd();

  // define histograms
  //TH1F* hMiss=new TH1F("hMiss", "Missing Mass", 1000, -20.0, 20.0);
  //TH1F* hMiss=new TH1F("hMiss", "Missing Mass", 2000, -10.0, 10.0);
  //TH1F* hMiss=new TH1F("hMiss", "Missing Mass", 2000, -6.0, 1.0);
  TH1F* hMiss=new TH1F("hMiss", "Missing Mass", 2000, -5.0, 1.0);
  //TH2F* hMissTheta=new TH2F("hMissTheta", "Missing mass vs. theta proton", 360,0,180,1000,-20,20);


  // define tree
  TTree* treeAnalysis = new TTree();

  // write generated data to tree
  // these values are with resolutions!!!!!!!!!!!!!
  treeAnalysis->Branch("genLightEnergy", &genLightEnergy, "genLightEnergy/D");
  treeAnalysis->Branch("genLightTheta", &genLightTheta, "genLightTheta/D");
  treeAnalysis->Branch("genLightPhi", &genLightPhi, "genLightPhi/D");
  treeAnalysis->Branch("genBeamX", &beamX, "genBeamX/F");
  treeAnalysis->Branch("genBeamY", &beamY, "genBeamY/F");
  treeAnalysis->Branch("genBeamZ", &beamZ, "genBeamZ/F");
  treeAnalysis->Branch("genBeamEnergy", &energyKinProj, "genBeamEnergy/F");
  treeAnalysis->Branch("genBeamTheta", &beamTheta, "genBeamTheta/F");
  treeAnalysis->Branch("genBeamPhi", &beamPhi, "genBeamPhi/F");
  treeAnalysis->Branch("genExcEn", &genExcEn, "genExcEn/F");

  // write simulated data to tree
  // these values are with resolutions!!!!!!!!!!!!!
  treeAnalysis->Branch("simLightEnergy1", &energyLoss, "simLightEnergy1/D");
  treeAnalysis->Branch("simLightEnergy2", &energyTotal, "simLightEnergy2/D");
  treeAnalysis->Branch("simLightEnergySum", &energySum, "simLightEnergySum/D");
  treeAnalysis->Branch("simLightTheta", &theta, "simLightTheta/D");
  treeAnalysis->Branch("simLightPhi", &phi, "simLightPhi/D");

  // new analysis data
  treeAnalysis->Branch("anaMissingMass", &miss, "anaMissingMass/F");
  treeAnalysis->Branch("anaProjGamma", &gammaProj, "anaProjGamma/F");
  treeAnalysis->Branch("anaProjTotalEnergy", &energyTotProj, "anaProjTotalEnergy/F");
  treeAnalysis->Branch("anaProjMomentum", &momentumProj, "anaProjMomentum/F");
  treeAnalysis->Branch("anaLightGamma", &gammaLight, "anaLightGamma/F");
  treeAnalysis->Branch("anaLightTotalEnergy", &energyTotLight, "anaLightTotalEnergy/F");
  treeAnalysis->Branch("anaLightMomentum", &momentumLight, "anaLightMomentum/F");
  treeAnalysis->Branch("anaHeavyGamma", &gammaHeavy, "anaHeavyGamma/F");
  treeAnalysis->Branch("anaHeavyTotalEnergy", &energyTotHeavy, "anaHeavyTotalEnergy/F");
  treeAnalysis->Branch("anaHeavyMomentum", &momentumHeavy, "anaHeavyMomentum/F");

  














  
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

    //todo: reset vairables!
    miss=1000.0;



    tree->GetEvent(e);
    treeBeam->GetEvent(e); // todo: make tree friend instead
    


    // take only events in backward direction
    // in very forward direction are usually punch-troughs of protons through the detector
    if(theta<90.0){
      continue;
    }
    
    theta/=180.0/TMath::Pi(); // this is important if simulation output is used
    
    // smear out data with detector position resolutions
    beamX=randomizer->Gaus(beamX, info->fResTargetX); // in mm
    beamY=randomizer->Gaus(beamY, info->fResTargetY);
    beamZ=randomizer->Gaus(beamZ, info->fResTargetZ);

    x=randomizer->Gaus(x, info->fResDet1X); 
    y=randomizer->Gaus(y, info->fResDet1Y); 
    z=randomizer->Gaus(z, info->fResDet1Z);
    
    // smear out data with detector energy resolutions
    energyLoss=randomizer->Gaus(energyLoss, info->fResDet1E); // in MeV 
    energyTotal=randomizer->Gaus(energyTotal, info->fResDet2E);  
    energySum=energyLoss+energyTotal;



    // projectile kinematics

    energyKinProj*=(Float_t)projA; 
    // Projectile data
    // at the moment from simulation input
    // todo: separate simulation including incoming tracking

    energyKinProj = randomizer->Gaus(energyKinProj, info->fResBeamE);
      
    //  only 132Sn (d,p) is implemented
    // 132Sn momentum/energy is fix
    gammaProj = (energyKinProj)/massProj + 1.0;
    //Float_t betaProj = TMath::Sqrt(1.0-(1.0/(gammaProj*gammaProj))); //just for cross checking
    momentumProj = massProj*TMath::Sqrt(gammaProj*gammaProj-1.0); 
    energyTotProj = massProj*gammaProj; //total energy
    
     
    TVector3 vProj(0.0, 0.0, momentumProj); 
    vProj.SetMagThetaPhi(momentumProj, beamTheta, beamPhi); // comment out this line to see the effect of no beam profile correction
//TVector3 vProj(0.0, 0.0, 1.0); 
//vProj.SetTheta(beamTheta); 
    
    // rotate by beam angular resolution
    vProj.RotateY(randomizer->Gaus(0.0, (info->fResTargetA))); // resolutions in mrad
    vProj.RotateX(randomizer->Gaus(0.0, (info->fResTargetB)));

    TLorentzVector lProj;
    lProj.SetPxPyPzE(vProj.X(), vProj.Y(), vProj.Z(), energyTotProj);
    ////lProj.Boost(0.0, 0.0, -betaProj);
    ////printf("Proj mass %f, from Lorentz %f\n", massProj, lProj.E());
//TLorentzVector lProj(vProj, energyKinProj+massProj);
//lProj.SetRho(TMath::Sqrt( (energyKinProj+massProj) * (energyKinProj+massProj) - (massProj*massProj) ));

//printf("Proj Ekin %f, Lor E %f\n", energyKinProj, lProj.E());








    // light ejectile kinematics


theta=genLightTheta;

    // get total energy and momentum of the light ejectile
    
    energySum=energyLoss+energyTotal;
    //gammaLight = energySum/massLight+1.0;      // simulated
gammaLight = genLightEnergy/massLight+1.0; // generated
   
   
   
    energyTotLight = gammaLight*massLight;
    //energyTotLight = energyKinLight+massLight;
    momentumLight = massLight*TMath::Sqrt(gammaLight*gammaLight-1.0);


    TVector3 vLight(x-beamX, y-beamY, z-beamZ); //momentum direction of proton
    vLight.SetMag(momentumLight);
vLight.SetMagThetaPhi(momentumLight, theta, phi);
    ////vLight.SetMagThetaPhi(momentumLight, genLightTheta, genLightPhi);

    TLorentzVector lLight;
    lLight.SetPxPyPzE(vLight.X(), vLight.Y(), vLight.Z(), energyTotLight);
    ////cout << "Light ejectile momentum " << lLight.Pt() << ", kin energy " << lLight.E()-massLight << ", total energy " << lLight.E() << endl;
    ////lLight.Boost(0.0, 0.0, -betaProj);
    ////printf("Light mass %f, from Lorentz %f\n", massLight, lLight.E());
//TVector3 vLight(0.0, 0.0, 1.0);
//vLight.SetTheta(theta/TMath::Pi()*180.0);
//TLorentzVector lLight(vLight, energyKinLight+massLight);
//lLight.SetRho(TMath::Sqrt( (energyKinLight+massLight) * (energyKinLight+massLight) - (massLight*massLight) ));











    //calculate momentum vector of heavy ejectile
    TVector3 vHeavy=vProj;
    vHeavy-=vLight;
    

    momentumHeavy = vHeavy.Mag();
    gammaHeavy = TMath::Sqrt((momentumHeavy*momentumHeavy)/(massHeavy*massHeavy)+1.0);
    //Float_t energyKinHeavy = (gammaHeavy-1.0)*massHeavy;     

    //generate the Lorentz vector of the outgoing heavy particle
    energyTotHeavy = gammaHeavy*massHeavy;
    //energyTotHeavy = energyKinHeavy+massHeavy;

    TLorentzVector lHeavy;
    lHeavy.SetPxPyPzE(vHeavy.X(), vHeavy.Y(), vHeavy.Z(), energyTotHeavy);
//TLorentzVector lHeavy(vHeavy, energyTotHeavy);
//lLight.SetRho(TMath::Sqrt( (energyTotHeavy) * (energyTotHeavy) - (massHeavy*massHeavy) ));



    //lHeavy.Boost(0.0, 0.0, -betaProj);
    //printf("Heavy mass %f, from Lorentz %f\n", massHeavy, lHeavy.E());



    //get the excitation energy from the missing mass
    TLorentzVector lMiss = lHeavy + lLight - lProj;
    //printf("lMiss: E %f, Mag %f, P %f, T %f, X %f, Y %f, Z %f\n", lMiss.E(), lMiss.Mag(), lMiss.P(), lMiss.T(), lMiss.X(), lMiss.Y(), lMiss.Z());
    miss = lMiss.E() - (massHeavy+massLight-massProj) - qValue;
    
    //miss = (lHeavy.E()-massHeavy) + (lLight.E()-massLight) - (lProj.E()-massProj) - qValue;
    //miss = (lHeavy.E()-massHeavy) + (lLight.E()-massLight) - qValue;

    //cout << "Missing mass " << miss << endl;

    ////if(miss>1){
    //  printf("miss mass %f, light theta %f, lightEnergy %f\n", miss, theta, energyKinLight);
    //  printf("light 3 vector: %f %f %f\n", vLight.X(), vLight.Y(), vLight.Z());
    //  printf("heavy 3 vector: %f %f %f\n", vHeavy.X(), vHeavy.Y(), vHeavy.Z());
    //  printf("light gamma*mass %f,  mass+ekin %f\n",gammaLight*massLight, massLight+energySum);

    ////}

   
    hMiss->Fill(miss);
    //hMissTheta->Fill(vLight.Theta()*180.0/TMath::Pi(), miss);


    energyKinProj/=(Float_t)projA; // MeV/u`

    treeAnalysis->Fill();
  
  }

  cout << nevents << " analyzed" << endl;
  
  fileAnalysis->cd();
  treeAnalysis->Write("analysis");
  //fileAnalysis->Close();

  hMiss->Draw();
  //hMissTheta->Draw();
  
   
  // if histograms shall be plotted, run theApp
  // otherwise program closes
  theApp->Run();
  delete theApp;  
  
  return 0;
}


