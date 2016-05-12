/**
 * 
 * main.cc for analysis of simulated transfer reactions
 * 
 * 
 **/


#include "main.hh"

#include "LibPerso.h"

#include "InputInfo.hh"
#include "DetectorInfo.hh"

#include "/home/philipp/sim/troja/include/DetectorGlobals.hh"

using namespace std;







//void FitTwoGaus(TH1F* hist)
//{
//
//  TF1* fit1 = new TF1("gaus1", "gaus(0)", -5, 1);
//  TF1* fit2 = new TF1("gaus2", "gaus(0)", -5, 1);
//
//  fit1->SetParameters(100, -4, 0.1);
//  fit2->SetParameters(100, 0, 0.1);
//
//  hist->Fit(fit1 , "","", -4.5, -3.5);
//  hist->Fit(fit2 , "","", -0.5, 0.5);
//
//}








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

  InputInfo* info = new InputInfo();
  info->parse(argv[1]);

  DetectorInfo* detInfo = new DetectorInfo();
  detInfo->Parse(info->fInFileNameGeometry);
  
  TRandom3* randomizer = new TRandom3();
  randomizer->SetSeed(0); 
  
  // define some constants
  
  // nuclear masses in MeV/u
  // masses should be defined in an external file!
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


  
  // tree with generated events / projectile information
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





  // Get tree with light particle simulation

  // Declaration of leaf types
  
  Int_t           eventNumber=0;
  Int_t           detHit[maxDetectors]={0}; // bool: 0 if no hit, 1 if hit in detector []
  Double_t        detEnergyLoss[maxDetectors]={0.0};
  Double_t        detEnergyLossNotSmeared[maxDetectors]={0.0};
  Int_t           detStripX[maxDetectors]={-1};
  Int_t           detStripY[maxDetectors]={-1};
  Double_t        recoPosX[maxDetectors]={NAN};
  Double_t        recoPosY[maxDetectors]={NAN};
  Double_t        recoPosZ[maxDetectors]={NAN};
  // for bug fixing:
  Int_t           firstDetID = -1; // find out which detector fired first
  Double_t        FIx=0.0; // first interaction point
  Double_t        FIy=0.0;
  Double_t        FIz=0.0;
  Int_t           FIdetID=-1; //from sim file
  //Int_t           detHitID[maxDetectors]={-1};
  
//  tree->SetBranchAddress("energy", &energy); //simulation input
  tree->SetBranchAddress("eventNumber", &eventNumber);
  tree->SetBranchAddress("detHit", detHit);
  tree->SetBranchAddress("energy", detEnergyLoss);
  tree->SetBranchAddress("energyNotSmeared", detEnergyLossNotSmeared);
  tree->SetBranchAddress("stripX", detStripX);
  tree->SetBranchAddress("stripY", detStripY);
  tree->SetBranchAddress("FIx", &FIx);
  tree->SetBranchAddress("FIy", &FIy);
  tree->SetBranchAddress("FIz", &FIz);
  tree->SetBranchAddress("FIdetID", &FIdetID);
  //tree->SetBranchAddress("recoPosX", recoPosX);
  //tree->SetBranchAddress("recoPosY", recoPosY);
  //tree->SetBranchAddress("recoPosZ", recoPosZ);
  //tree->SetBranchAddress("detHitID", detHitID);





  // define output of analysis
  Int_t detHitMul=0; // detector hit multiplicity
  Int_t detHitMul3 = 0; // aux, count events with more than 2 detectors fired
  Double_t gammaProj=0.0, gammaLight=0.0;
  Double_t energyTotProj=0.0, energyTotLight=0.0;
//  Double_t energyKinHeavy=0.0, energyTotHeavy=0.0, gammaHeavy=0.0, momentumHeavy=0.0;
  Double_t energyKinLight=0.0; // is sum of all energy losses
  Double_t momentumProj=0.0, momentumLight=0.0; 
  Double_t simDetectorHitPos[3]={0.0}; // x, y, z; position used for analysis
  Double_t thetaLightLab=0.0, thetaLightCM=0.0, phiLight=0.0;
  Double_t miss=0.0;
  
  TFile* fileAnalysis = new TFile(info->fOutFileNameAnalysis, "recreate");
  fileAnalysis->cd();

  char tmpName[50];

  // define histograms
  //TH1F* hMiss=new TH1F("hMiss", "Missing Mass", 1000, -20.0, 20.0);
  //TH1F* hMiss=new TH1F("hMiss", "Missing Mass", 2000, -10.0, 10.0);
  //TH1F* hMiss=new TH1F("hMiss", "Missing Mass", 2000, -6.0, 1.0);
  TH1F* hMiss=new TH1F("hMiss", "Missing Mass", 600, -5.0, 1.0);
  //TH2F* hMissTheta=new TH2F("hMissTheta", "Missing mass vs. theta proton", 360,0,180,1000,-20,20);
  
  TH2F* hdEE=new TH2F("hdEE", "delta E vs. E proton", 1000,0,20,100,0,8);
  TH2F* hEth=new TH2F("hEth", "E proton vs. theta lab", 1800,0,180,500,0,50);

  TH1F* hThetaLab = new TH1F("hThetaLab","Theta Lab",1800,0,180);
  TH1F* hThetaCM = new TH1F("hThetaCM","Theta CM",1800,0,180);


  // define tree
  TTree* treeAnalysis = new TTree();

  treeAnalysis->Branch("eventNumber", &eventNumber, "eventNumber/I");
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
  sprintf(tmpName, "simDetectorEnergy[%d]/D", maxDetectors);
  treeAnalysis->Branch("simDetectorEnergyLoss", detEnergyLoss, tmpName);
  treeAnalysis->Branch("simDetectorEnergyLossNotSmeared", detEnergyLossNotSmeared, tmpName);
  treeAnalysis->Branch("simLightKinEnergy", &energyKinLight, "simLightEnergy/D"); // sum of detector energy losses

  treeAnalysis->Branch("simFIx", &FIx, "simFIx/D");
  treeAnalysis->Branch("simFIy", &FIy, "simFIy/D");
  treeAnalysis->Branch("simFIz", &FIz, "simFIz/D");
  treeAnalysis->Branch("simFIdetID", &FIdetID, "simFIdetID/I");
  treeAnalysis->Branch("anaFIdetID", &firstDetID, "anaFIdetID/I");

  sprintf(tmpName, "detStripX[%d]/I", maxDetectors);
  treeAnalysis->Branch("detStripX", detStripX, tmpName);
  sprintf(tmpName, "detStripY[%d]/I", maxDetectors);
  treeAnalysis->Branch("detStripY", detStripY, tmpName);


  
  treeAnalysis->Branch("simLightThetaLab", &thetaLightLab, "simLightThetaLab/D");
  treeAnalysis->Branch("simLightThetaCM", &thetaLightCM, "simLightThetaCM/D");
  treeAnalysis->Branch("simLightPhi", &phiLight, "simLightPhi/D");


  // new analysis data
  
  treeAnalysis->Branch("anaDetectorHitMul", &detHitMul, "anaDetectorHitMul/I");
  sprintf(tmpName, "anaDetectorHitPos[3]/D");
  treeAnalysis->Branch("anaDetectorHitPos", simDetectorHitPos, tmpName);
  
  treeAnalysis->Branch("anaMissingMass", &miss, "anaMissingMass/D");
  //treeAnalysis->Branch("anaProjGamma", &gammaProj, "anaProjGamma/F");
  //treeAnalysis->Branch("anaProjTotalEnergy", &energyTotProj, "anaProjTotalEnergy/F");
  //treeAnalysis->Branch("anaProjMomentum", &momentumProj, "anaProjMomentum/F");
  //treeAnalysis->Branch("anaLightGamma", &gammaLight, "anaLightGamma/F");
  //treeAnalysis->Branch("anaLightTotalEnergy", &energyTotLight, "anaLightTotalEnergy/F");
  //treeAnalysis->Branch("anaLightMomentum", &momentumLight, "anaLightMomentum/F");
  //treeAnalysis->Branch("anaHeavyGamma", &gammaHeavy, "anaHeavyGamma/F");
  //treeAnalysis->Branch("anaHeavyTotalEnergy", &energyTotHeavy, "anaHeavyTotalEnergy/F");
  //treeAnalysis->Branch("anaHeavyMomentum", &momentumHeavy, "anaHeavyMomentum/F");

  



  TStopwatch* watch = new TStopwatch();
  watch->Start();

  
  Int_t nevents=tree->GetEntries();
  cout << "Number of entries found in tree: " << nevents << endl;
  cout << "Starting analysis ..." << endl;

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

  Int_t goodEvents=0;
  

// **********************************************************************************************************************************************************
// ************************************************** start event loop **************************************************************************************
// **********************************************************************************************************************************************************


  for(Int_t e=0; e<nevents; e++){

    //todo: reset vairables!
    miss=1000.0;
    energyKinLight=0.0;
    FIdetID=-1;
    firstDetID=-1;

    for(Int_t p=0;p<3;p++){
      simDetectorHitPos[p]=0.0;
    }

    for(Int_t d=0; d<maxDetectors; d++){
      detEnergyLoss[d] = 0.0;
      detHit[d] = 0;
      recoPosX[d] = NAN;
      recoPosY[d] = NAN;
      recoPosZ[d] = NAN;
    }



    tree->GetEvent(e);
    treeBeam->GetEvent(e); // todo: make tree friend instead
    
    detHitMul = 0; // detector hit multiplicity
    //Int_t firstDetID = -1; // find out which detector fired first
    //Int_t seconDetID = -1; // find out which detector fired second 

    // sum up all energy losses
    for(Int_t d=0; d<maxDetectors; d++){
      detHitMul += detHit[d];
      if(detHit[d]==1){

        energyKinLight+=detEnergyLoss[d];
        
        // reconstruct position from strip number
        detInfo->CalcHitPosition(d, detStripX[d], detStripY[d], recoPosX[d], recoPosY[d], recoPosZ[d]);

        if(detInfo->IsInFront(d, firstDetID)){
          firstDetID = d;
        }
        //else{
        //  seconDetID = d;
        //}


      }
    }
    
    if(detHitMul==0){
      continue;
    }
    else if(detHitMul>2){
      //printf("Warning: more than 2 detectors fired in event %d! Such events are not handled properly!\n", e);
      detHitMul3++;

      // todo: a more sophisticated routine to find the first hit
      // but second order problem, it happens in about 700 events out of 10 mio! 
    }

    goodEvents++;
    
    //printf("event %d, have %d detector hits with sum energy %f\n", e, detHitSum, energyKinLight);

    //simDetectorHitPos[0] = FIx;
    //simDetectorHitPos[1] = FIy;
    //simDetectorHitPos[2] = FIz;
    //simDetectorHitPos[0] = recoPosX[FIdetID];
    //simDetectorHitPos[1] = recoPosY[FIdetID];
    //simDetectorHitPos[2] = recoPosZ[FIdetID];
    simDetectorHitPos[0] = recoPosX[firstDetID];
    simDetectorHitPos[1] = recoPosY[firstDetID];
    simDetectorHitPos[2] = recoPosZ[firstDetID];

    // aux:
    //Double_t x = FIx;
    //Double_t y = FIy;
    //Double_t z = FIz;

    // take only events in backward direction
    // in very forward direction are usually punch-troughs of protons through the detector
//    if(theta<90.0){
//    if(theta==0.0){
//      continue;
//    }
    
//    theta/=180.0/TMath::Pi(); // is actually not used anymore
    
    // smear out data with detector position resolutions
    beamX=randomizer->Gaus(beamX, info->fResTargetX); // in mm
    beamY=randomizer->Gaus(beamY, info->fResTargetY);
    beamZ=randomizer->Gaus(beamZ, info->fResTargetZ);
    
    // obsolete: 
//    x=randomizer->Gaus(x, info->fResDet1X); 
//    y=randomizer->Gaus(y, info->fResDet1Y); 
//    z=randomizer->Gaus(z, info->fResDet1Z);
    
    // smear out data with detector energy resolutions
//    energyLoss=randomizer->Gaus(energyLoss, info->fResDet1E); // in MeV 
//    energyTotal=randomizer->Gaus(energyTotal, info->fResDet2E);  
//    energySum=energyLoss+energyTotal;



    // projectile kinematics

    energyKinProj*=(Float_t)projA; 
    // Projectile data
    // at the moment from simulation input
    // todo: separate simulation including incoming tracking

    energyKinProj = randomizer->Gaus(energyKinProj, info->fResBeamE);
      
    gammaProj = (energyKinProj)/massProj + 1.0;
    //Float_t betaProj = TMath::Sqrt(1.0-(1.0/(gammaProj*gammaProj))); //just for cross checking
    momentumProj = massProj*TMath::Sqrt(gammaProj*gammaProj-1.0); 
    energyTotProj = massProj*gammaProj; //total energy
    
     
    TVector3 vProj(0.0, 0.0, momentumProj); 
    vProj.SetMagThetaPhi(momentumProj, beamTheta, beamPhi); // comment out this line to see the effect of no beam profile correction
    
    // rotate by beam angular resolution
    vProj.RotateY(randomizer->Gaus(0.0, (info->fResTargetA))); // resolutions in mrad
    vProj.RotateX(randomizer->Gaus(0.0, (info->fResTargetB)));



    // center of mass kinematic values
    Float_t energyCm = TMath::Sqrt(massProj*massProj + massTarget*massTarget + 2.0*energyTotProj*massTarget);
    Float_t betaCm = momentumProj/(energyTotProj+massTarget);
    TVector3 vCm(vProj); // direction of projectile including beam profile
    vCm.SetMag(betaCm);





    // light ejectile kinematics

    // get total energy and momentum of the light ejectile
    
    //energySum=energyLoss+energyTotal;
    //gammaLight = energySum/massLight+1.0;   
    gammaLight = energyKinLight/massLight+1.0;      // simulated
    //gammaLight = genLightEnergy/massLight+1.0; // generated
    //theta=genLightTheta; // generated theta
   
   
    energyTotLight = gammaLight*massLight;
    //energyTotLight = energyKinLight+massLight;
    momentumLight = massLight*TMath::Sqrt(gammaLight*gammaLight-1.0);

    //printf("lightTheta %f, lightEnergy %f \n", theta, energySum);
    //printf("x %f, y %f, z %f, beamX %f, beamY %f, beamZ %f\n",x ,y, z, beamX, beamY, beamZ);
    //TVector3 vLight(x-beamX, y-beamY, z-beamZ); //momentum direction of proton
    TVector3 vLight(simDetectorHitPos[0]-beamX, simDetectorHitPos[1]-beamY, simDetectorHitPos[2]-beamZ); //momentum direction of proton
    //printf("vLight.Mag %f\n", vLight.Mag());
    vLight.SetMag(momentumLight);
    //vLight.SetMagThetaPhi(momentumLight, theta, phi); // without beam position spread
    
    // for the root tree
    thetaLightLab = vLight.Theta(); 
    phiLight = vLight.Phi(); 

    TLorentzVector lLight;
    lLight.SetVect(vLight);
    lLight.SetE(energyTotLight);

    TLorentzVector lL=lLight;

    lLight.Boost(-vCm);

    //thetaLightCM=lLight.Theta();

    lL.Boost(vCm);
    thetaLightCM=lL.Theta();

    hThetaLab->Fill(thetaLightLab*180.0/TMath::Pi());
    hThetaCM->Fill(thetaLightCM*180.0/TMath::Pi());




    // heavy ejectile kinematics in center of mass system
    TLorentzVector lHeavy;
    lHeavy.SetVect(-lLight.Vect());
    lHeavy.SetE(energyCm-lLight.E());

    miss = -lHeavy.M()+massHeavy;

    //printf("miss %f\n", miss);
    

    // fill histograms

    hMiss->Fill(miss);
    //hMissTheta->Fill(vLight.Theta()*180.0/TMath::Pi(), miss);

    hdEE->Fill(energyKinLight,detEnergyLoss[firstDetID]);
    hEth->Fill(thetaLightLab*180.0/TMath::Pi(), energyKinLight);


    energyKinProj/=(Float_t)projA; // MeV/u`

    treeAnalysis->Fill();
  
  }

  cout << nevents << " events processed! " << goodEvents << " events used in analysis! Ratio: " << (Float_t)goodEvents/(Float_t)nevents*100.0 << "%" << endl;
  cout << detHitMul3 << " events with detector hit multiplicity 3 or larger (angle reconstruction might not be correct for these events)" << endl;

  watch->Stop();
  cout << "Took: real time " << watch->RealTime() << "sec., CPU time " << watch->CpuTime() << " sec." << endl;
  cout << endl;
  
  fileAnalysis->cd();
  treeAnalysis->Write("analysis");

  hMiss->Write("missingMass");
  hdEE->Write("dEE");
  hEth->Write("Eth");
  hThetaLab->Write("hThetaLab");
  hThetaCM->Write("hThetaCM");

  printf("Analyzed events writen to file '%s'\n", info->fOutFileNameAnalysis);
  




  // plot histograms
  
  TCanvas* candEE = new TCanvas();
  candEE->cd();
  hdEE->Draw("colz");


  TCanvas* canTheta = new TCanvas();
  canTheta->Divide(1,2);
  canTheta->cd(1);
  hThetaLab->Draw();
  canTheta->cd(2);
  hThetaCM->Draw();

  

  TCanvas* canMissMass = new TCanvas();
  canMissMass->cd();
  hMiss->GetXaxis()->SetTitle("E_{miss} / MeV");
  hMiss->GetYaxis()->SetTitle("cts");
  gStyle->SetOptStat(0);
  hMiss->Draw();
  //hMissTheta->Draw();

  //FitTwoGaus(hMiss);  
   
  TF1* fit1 = new TF1("gaus1", "gaus(0)", -5, 1);
  TF1* fit2 = new TF1("gaus2", "gaus(0)", -5, 1);

  fit1->SetParameters(100, -2, 0.1);
  fit2->SetParameters(100, 0, 0.1);

  fit1->SetParLimits(0,0,1000000);
  fit1->SetParLimits(1,-4.5,-1.5);
  fit1->SetParLimits(2,0.0,1.0  );
  fit2->SetParLimits(0,0,1000000);
  fit2->SetParLimits(1,-0.5,0.5);
  fit2->SetParLimits(2,0.0,1.0  );

  hMiss->Fit(fit1 , "","", -4.5, -1.5);
  hMiss->Fit(fit2 , "","", -0.5, 0.5);

  fit1->Draw("same"); 
  

  TCanvas* canEth = new TCanvas();
  canEth->cd();
  hEth->Draw("colz");



  //printf("Integrals: Fit1 = %f, Fit2 = %f\n", fit1->Integral(), fit2->Integral()); 
  

  
  
  
  
  
  // if histograms shall be plotted, run theApp
  // otherwise program closes
  theApp->Run();
//  fileAnalysis->Close();
//  delete theApp;  
  
  return 0;
}


