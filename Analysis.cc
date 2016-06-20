#include "Analysis.hh"




Analysis::Analysis(InputInfo *i, DetectorInfo* d){
  info = i;
  detInfo = d;

  Bool_t success = Init();
  if(success){
    cout << "Analysis: initalization successfull! " << endl;
  }else{
    cout << "Analysis: Error during initialization! aborting..." << endl;
    abort();
  }
}

Analysis::~Analysis(){

}


void Analysis::ResetVariables(){

  massProj=0.0; 
  massTarget=0.0; 
  massLight=0.0; 
  massHeavy=0.0;
  qValue=0.0;

  energyKinProj = 10.0*132.0;  // 10 MeV/u
  beamX = 0.0; 
  beamY = 0.0;  
  beamZ = 0.0;    
  beamTheta = 0.0;
  beamPhi = 0.0;
  genLightEnergy = 0.0;
  genLightTheta = 0.0;
  genLightPhi = 0.0;
  genExcEn = 0.0;

  eventNumber=0; // read from tree
  for(Int_t i=0; i<maxDetectors; i++){
    detHit[i]=0; // bool: 0 if no hit, 1 if hit in detector []
    detEnergyLoss[i]=0.0;
    detEnergyLossNotSmeared[i]=0.0;
    detStripX[i]=-1;
    detStripY[i]=-1;
    recoPosX[i]=NAN;
    recoPosY[i]=NAN;
    recoPosZ[i]=NAN;
  }
  // for bug fixing:
  firstDetID = -1; // find out which detector fired first
  FIx=0.0; // first interaction point
  FIy=0.0;
  FIz=0.0;
  FIdetID=-1; //from sim file
  //detHitID[maxDetectors]={-1};


  // define output of analysis
  detHitMul=0; // detector hit multiplicity
  //detHitMul3 = 0; // aux, count events with more than 2 detectors fired
  gammaProj=0.0;
  gammaLight=0.0;
  energyTotProj=0.0;
  energyTotLight=0.0;
  energyKinLight=0.0; // is sum of all energy losses
  momentumProj=0.0;
  momentumLight=0.0;
  for(Int_t i=0; i<3; i++){simDetectorHitPos[i]=0.0;} // x, y, z; position used for analysis
  thetaLightLab=0.0;
  thetaLightCM=0.0;
  phiLight=0.0;
  miss=0.0;



}

Bool_t Analysis::Init(){

  ResetVariables();
   
  fileBeam = TFile::Open(info->fOutFileNameMakeEvents,"read");

  if(!fileBeam){
   cout << "makeEvents root file not found!" << endl;
   return false;
  }


  projA=info->fProjA;
  projZ=info->fProjZ;
  targA=info->fTargetA;
  targZ=info->fTargetZ;

  
  treeBeam=(TTree*)fileBeam->Get("events");
  //TTree* tree=(TTree*)infile->Get("events"); //simulation input
  if(!treeBeam){
    cout << "TTree 'events' not found in makeEvents root file!" << endl;
    //cout << "TTree 'events' not found!" << endl;
    return false;
  }



  infile = TFile::Open(info->fOutFileNameTroja,"read");

  if(!infile){
   cout << "geant (troja) root file not found!" << endl;
   return false;
  }

  tree=(TTree*)infile->Get("troja");
  //TTree* tree=(TTree*)infile->Get("events"); //simulation input
  if(!tree){
    cout << "TTree 'troja' not found in geant root file!" << endl;
    //cout << "TTree 'events' not found!" << endl;
    return false;
  }



  // tree with generated events / projectile information

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


  fileAnalysis = new TFile(info->fOutFileNameAnalysis, "recreate");

  char tmpName[50];


  // define tree
  treeAnalysis1 = new TTree();

  sprintf(tmpName, "detectorEnergy[%d]/D", maxDetectors);
  treeAnalysis1->Branch("detectorEnergyLoss", detEnergyLoss, tmpName);
  treeAnalysis1->Branch("detectorEnergyLossNotSmeared", detEnergyLossNotSmeared, tmpName);
  sprintf(tmpName, "detectorStripX[%d]/I", maxDetectors);
  treeAnalysis1->Branch("detectorStripX", detStripX, tmpName);
  sprintf(tmpName, "detectorStripY[%d]/I", maxDetectors);
  treeAnalysis1->Branch("detectorStripY", detStripY, tmpName);

  treeAnalysis1->Branch("simLightKinEnergy", &energyKinLight, "simLightEnergy/D"); // sum of detector energy losses

  treeAnalysis1->Branch("simFIx", &FIx, "simFIx/D");
  treeAnalysis1->Branch("simFIy", &FIy, "simFIy/D");
  treeAnalysis1->Branch("simFIz", &FIz, "simFIz/D");
  treeAnalysis1->Branch("simFIdetID", &FIdetID, "simFIdetID/I");
  treeAnalysis1->Branch("anaFIdetID", &firstDetID, "anaFIdetID/I");




  treeAnalysis2 = new TTree();

  treeAnalysis2->Branch("eventNumber", &eventNumber, "eventNumber/I");
  // write generated data to tree
  // these values are with resolutions!!!!!!!!!!!!!
  treeAnalysis2->Branch("genLightEnergy", &genLightEnergy, "genLightEnergy/D");
  treeAnalysis2->Branch("genLightTheta", &genLightTheta, "genLightTheta/D");
  treeAnalysis2->Branch("genLightPhi", &genLightPhi, "genLightPhi/D");
  treeAnalysis2->Branch("genBeamX", &beamX, "genBeamX/F");
  treeAnalysis2->Branch("genBeamY", &beamY, "genBeamY/F");
  treeAnalysis2->Branch("genBeamZ", &beamZ, "genBeamZ/F");
  treeAnalysis2->Branch("genBeamEnergy", &energyKinProj, "genBeamEnergy/F");
  treeAnalysis2->Branch("genBeamTheta", &beamTheta, "genBeamTheta/F");
  treeAnalysis2->Branch("genBeamPhi", &beamPhi, "genBeamPhi/F");
  treeAnalysis2->Branch("genExcEn", &genExcEn, "genExcEn/F");

  // write simulated data to tree
  // these values are with resolutions!!!!!!!!!!!!!
  sprintf(tmpName, "detectorEnergy[%d]/D", maxDetectors);
  treeAnalysis2->Branch("detectorEnergyLoss", detEnergyLoss, tmpName);
  treeAnalysis2->Branch("detectorEnergyLossNotSmeared", detEnergyLossNotSmeared, tmpName);
  sprintf(tmpName, "detectorStripX[%d]/I", maxDetectors);
  treeAnalysis2->Branch("detectorStripX", detStripX, tmpName);
  sprintf(tmpName, "detectorStripY[%d]/I", maxDetectors);
  treeAnalysis2->Branch("detectorStripY", detStripY, tmpName);



  treeAnalysis2->Branch("simLightKinEnergy", &energyKinLight, "simLightEnergy/D"); // sum of detector energy losses

  treeAnalysis2->Branch("simFIx", &FIx, "simFIx/D");
  treeAnalysis2->Branch("simFIy", &FIy, "simFIy/D");
  treeAnalysis2->Branch("simFIz", &FIz, "simFIz/D");
  treeAnalysis2->Branch("simFIdetID", &FIdetID, "simFIdetID/I");
  treeAnalysis2->Branch("anaFIdetID", &firstDetID, "anaFIdetID/I");



  treeAnalysis2->Branch("simLightThetaLab", &thetaLightLab, "simLightThetaLab/D");
  treeAnalysis2->Branch("simLightThetaCM", &thetaLightCM, "simLightThetaCM/D");
  treeAnalysis2->Branch("simLightPhi", &phiLight, "simLightPhi/D");


  // new analysis data

  treeAnalysis2->Branch("anaDetectorHitMul", &detHitMul, "anaDetectorHitMul/I");
  sprintf(tmpName, "anaDetectorHitPos[3]/D");
  treeAnalysis2->Branch("anaDetectorHitPos", simDetectorHitPos, tmpName);

  treeAnalysis2->Branch("anaMissingMass", &miss, "anaMissingMass/D");
  //treeAnalysis2->Branch("anaProjGamma", &gammaProj, "anaProjGamma/F");
  //treeAnalysis2->Branch("anaProjTotalEnergy", &energyTotProj, "anaProjTotalEnergy/F");
  //treeAnalysis2->Branch("anaProjMomentum", &momentumProj, "anaProjMomentum/F");
  //treeAnalysis2->Branch("anaLightGamma", &gammaLight, "anaLightGamma/F");
  //treeAnalysis2->Branch("anaLightTotalEnergy", &energyTotLight, "anaLightTotalEnergy/F");
  //treeAnalysis2->Branch("anaLightMomentum", &momentumLight, "anaLightMomentum/F");
  //treeAnalysis2->Branch("anaHeavyGamma", &gammaHeavy, "anaHeavyGamma/F");
  //treeAnalysis2->Branch("anaHeavyTotalEnergy", &energyTotHeavy, "anaHeavyTotalEnergy/F");
  //treeAnalysis2->Branch("anaHeavyMomentum", &momentumHeavy, "anaHeavyMomentum/F");


  Int_t nevents=tree->GetEntries();
  cout << "Number of entries found in tree: " << nevents << endl;
  //cout << "Starting analysis ..." << endl;

  Int_t neventsBeam=treeBeam->GetEntries();
  if(nevents!=neventsBeam){
    if(nevents<neventsBeam){
      cout << "INFO: " << nevents << " simulated and " << neventsBeam << " found in simulation input. Analyzing only " << nevents << "." << endl;
    }
    if(nevents>neventsBeam){
      cout << "ERROR: more events simulated (" << nevents << ") than found in simulation input file (" << neventsBeam << ")! Check this!" << endl;
      return false;
    }
  }


  
  return true;

} // Init






void Analysis::ProcessDetectorHits(){  // private

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

  //simDetectorHitPos[0] = FIx;
  //simDetectorHitPos[1] = FIy;
  //simDetectorHitPos[2] = FIz;
  simDetectorHitPos[0] = recoPosX[firstDetID];
  simDetectorHitPos[1] = recoPosY[firstDetID];
  simDetectorHitPos[2] = recoPosZ[firstDetID];


} // ProcessDetectorHits






void Analysis::Analysis1(){
  
  // define histograms
  TH2F* hdEE=new TH2F("hdEE", "delta E vs. E proton", 1000,0,20,100,0,8);


  Int_t goodEvents=0;
  Int_t nevents=tree->GetEntries();

  for(Int_t e=0; e<nevents; e++){

    ResetVariables();


    tree->GetEvent(e);
    //treeBeam->GetEvent(e); // todo: make tree friend instead

    ProcessDetectorHits(); // private function, manipulates private variables only

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


    hdEE->Fill(energyKinLight,detEnergyLoss[firstDetID]);

    treeAnalysis1->Fill();


  } // end of event loop


  cout << nevents << " events processed! " << goodEvents << " events used in analysis! Ratio: " << (Float_t)goodEvents/(Float_t)nevents*100.0 << "%" << endl;
  cout << detHitMul3 << " events with detector hit multiplicity 3 or larger (angle reconstruction might not be correct for these events)" << endl;

  fileAnalysis->cd();
  treeAnalysis1->Write("analysis1");

//  hdEE->Write("dEE");

  printf("Analysis1: Analyzed events writen to file '%s'\n\n", info->fOutFileNameAnalysis);



} // Analysis1


void Analysis::Analysis2(){
 
  TRandom3* randomizer = new TRandom3();
  randomizer->SetSeed(0);
 



  TTree* treeHeader=(TTree*)fileBeam->Get("header");
  if(!treeHeader){
    cout << "TTree 'header' not found in makeEvents root file!" << endl;
    //return 0;
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
  printf("Obtained masses: projectile %f, target %f, light ejectile %f, heavy ejectile %f; Q-value %f\n", massProj, massTarget, massLight, massHeavy, qValue);












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


  Int_t goodEvents=0;
  Int_t nevents=tree->GetEntries();
  
  // reset the aux vaiable 
  detHitMul3=0;
  
  for(Int_t e=0; e<nevents; e++){
    
    ResetVariables();

    treeHeader->GetEntry(0); // todo: take these values from Kinematics/Nucleus/mass file

    tree->GetEvent(e);
    treeBeam->GetEvent(e); // todo: make tree friend instead
    
    ProcessDetectorHits();

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
    if(info->NoBeamTracking()){
      beamX=0.0;
      beamY=0.0;
      beamZ=0.0;
    }else{
      beamX=randomizer->Gaus(beamX, info->fResTargetX); // in mm
      beamY=randomizer->Gaus(beamY, info->fResTargetY);
      beamZ=randomizer->Gaus(beamZ, info->fResTargetZ);
    }
    


    // projectile kinematics

    energyKinProj*=(Float_t)projA; 
    // Projectile data
    // at the moment from simulation input
    // todo: separate simulation including incoming tracking

    if(!info->NoBeamTracking()){
      energyKinProj = randomizer->Gaus(energyKinProj, info->fResBeamE);
    }
      
    gammaProj = (energyKinProj)/massProj + 1.0;
    //Float_t betaProj = TMath::Sqrt(1.0-(1.0/(gammaProj*gammaProj))); //just for cross checking
    momentumProj = massProj*TMath::Sqrt(gammaProj*gammaProj-1.0); 
    energyTotProj = massProj*gammaProj; //total energy
    
     
    TVector3 vProj(0.0, 0.0, 1.0); 
    if(!info->NoBeamTracking()){
      vProj.SetMagThetaPhi(1.0, beamTheta, beamPhi);                                         // comment out this line to see the effect of no beam profile correction
    
      // rotate by beam angular resolution
      vProj.RotateY(randomizer->Gaus(0.0, (info->fResTargetA)/1000.0)); // resolutions in mrad
      vProj.RotateX(randomizer->Gaus(0.0, (info->fResTargetB)/1000.0));
    }
    vProj.SetMag(momentumProj);

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
    //hEth->Fill(thetaLightLab*180.0/TMath::Pi(), energyKinLight);


    energyKinProj/=(Float_t)projA; // MeV/u`

    treeAnalysis2->Fill();
  
  }

  cout << nevents << " events processed! " << goodEvents << " events used in analysis! Ratio: " << (Float_t)goodEvents/(Float_t)nevents*100.0 << "%" << endl;
  cout << detHitMul3 << " events with detector hit multiplicity 3 or larger (angle reconstruction might not be correct for these events)" << endl;



  fileAnalysis->cd();
  treeAnalysis2->Write("analysis2");

  hdEE->Write("dEE");
  hMiss->Write("missingMass");
  hEth->Write("Eth");
  hThetaLab->Write("hThetaLab");
  hThetaCM->Write("hThetaCM");

  printf("Analysis2: Analyzed events writen to file '%s'\n\n", info->fOutFileNameAnalysis);



} // Analysis2

