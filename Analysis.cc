#include "Analysis.hh"

#include "Nucleus.hh"

#define warningType maxCutType+1


Analysis::Analysis(InputInfo *i, DetectorInfo* d){
  
  info = i;
  detInfo = d;
  
  for(Int_t f=0; f<maxCutType; f++){
  for(Int_t i=0; i<maxCuts; i++){
    cutExists[f][i]=false; // -1: doesn't exist, >-1: type of cut/reaction (0 = elastic, 1 = (d,p), ...)
  }
  }

  Bool_t success = Init();
  if(success){
    cout << "Analysis: initalization successfull! " << endl;
    CreateHeader();
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

  //energyKinProj = 10.0*132.0;  // 10 MeV/u
  energyKinProj = 0.0;  
  beamX = 0.0; 
  beamY = 0.0;  
  beamZ = 0.0;    
  beamA = 0.0;  
  beamB = 0.0;    
  beamTheta = 0.0;
  beamPhi = 0.0;
  genLightEnergy = 0.0;
  genLightTheta = 0.0;
  genLightThetaCM = 0.0;
  genLightPhi = 0.0;
  genExcEn = 0.0;
  genState = -1;

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

  type = -1;

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

  energyGamma=0.0;
  grapeDet=-1;



}



Bool_t Analysis::Init(){

  ResetVariables();


  randomizer = new TRandom3();
  randomizer->SetSeed(0);

  // define histograms
  //hMiss=new TH1F("hMiss", "Missing Mass", 1000, -20.0, 20.0);
  //hMiss=new TH1F("hMiss", "Missing Mass", 2000, -10.0, 10.0);
  //hMiss=new TH1F("hMiss", "Missing Mass", 2000, -6.0, 1.0);
  hMiss=new TH1F("hMiss", "Missing Mass", 600, -5.0, 1.0);
  //hMissTheta=new TH2F("hMissTheta", "Missing mass vs. theta proton", 360,0,180,1000,-20,20);

  //hdEE=new TH2F("hdEE_analysis2", "delta E vs. E, Analysis2", 1000,0,50,100,0,8);
  hdEE_A=new TH2F("hdEE_A_analysis2", "delta E vs. E, all directions, Analysis2", 1000,0,25,100,0,8);
  hdEE_B=new TH2F("hdEE_B_analysis2", "delta E vs. E, backward direction, Analysis2", 1000,0,25,100,0,8);
  hdEE_F=new TH2F("hdEE_F_analysis2", "delta E vs. E, forward direction, Analysis2", 1000,0,25,100,0,8);
  
  hEth=new TH2F("hEth", "E proton vs. theta lab", 180,0,180,600,0,60);

  hThetaLab = new TH1F("hThetaLab","Theta Lab",180,0,180);
  hThetaCM = new TH1F("hThetaCM","Theta CM",180,0,180);
   
  fileBeam = TFile::Open(info->fOutFileNameMakeEvents,"read");

  if(!fileBeam){
   cout << "makeEvents root file not found!" << endl;
   return false;
  }

  projA = info->fProjA;
  projZ = info->fProjZ;
  targA = info->fTargetA;
  targZ = info->fTargetZ;
  
  treeBeam = (TTree*)fileBeam->Get("events");
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

  tree = (TTree*)infile->Get("troja");
  //TTree* tree=(TTree*)infile->Get("events"); //simulation input
  if(!tree){
    cout << "TTree 'troja' not found in geant root file!" << endl;
    //cout << "TTree 'events' not found!" << endl;
    return false;
  }

  // tree with generated events / projectile information

  treeBeam->SetBranchAddress("lightEnergy", &genLightEnergy);
  treeBeam->SetBranchAddress("lightTheta", &genLightTheta);
  treeBeam->SetBranchAddress("lightThetaCM", &genLightThetaCM);
  treeBeam->SetBranchAddress("lightPhi", &genLightPhi);
  treeBeam->SetBranchAddress("beamEnergy", &energyKinProj); // is in AMeV
  treeBeam->SetBranchAddress("beamX", &beamX);
  treeBeam->SetBranchAddress("beamY", &beamY);
  treeBeam->SetBranchAddress("beamZ", &beamZ);
  treeBeam->SetBranchAddress("beamA", &beamA);
  treeBeam->SetBranchAddress("beamB", &beamB);
  treeBeam->SetBranchAddress("beamTheta", &beamTheta);
  treeBeam->SetBranchAddress("beamPhi", &beamPhi);
  treeBeam->SetBranchAddress("state", &genState);
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
  
  tree->SetBranchAddress("grapeEnergy", &energyGamma);
  tree->SetBranchAddress("grapeDetector", &grapeDet);


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
  //treeAnalysis1->Branch("simLightThetaLab", &thetaLightLab, "simLightThetaLab/D");
  treeAnalysis1->Branch("anaFIdetID", &firstDetID, "anaFIdetID/I");
  
  treeAnalysis1->Branch("energyGamma", &energyGamma, "energyGamma/D");
  treeAnalysis1->Branch("energyGammaDC", &energyGammaDC, "energyGammaDC/D");



  for(Int_t f=0; f<maxCutType; f++){
    
    treeAnalysis2[f] = new TTree(Form("analysis2_%d", f), Form("Analysis2 (missing mass) tree, reaction channel/type %d", f));

    treeAnalysis2[f]->Branch("eventNumber", &eventNumber, "eventNumber/I");
    // write generated data to tree
    // these values are with resolutions!!!!!!!!!!!!!
    treeAnalysis2[f]->Branch("genLightEnergy", &genLightEnergy, "genLightEnergy/D");
    treeAnalysis2[f]->Branch("genLightThetaLab", &genLightTheta, "genLightThetaLab/D");
    treeAnalysis2[f]->Branch("genLightThetaCM", &genLightThetaCM, "genLightThetaCM/D");
    treeAnalysis2[f]->Branch("genLightPhi", &genLightPhi, "genLightPhi/D");
    treeAnalysis2[f]->Branch("genBeamX", &beamX, "genBeamX/F");
    treeAnalysis2[f]->Branch("genBeamY", &beamY, "genBeamY/F");
    treeAnalysis2[f]->Branch("genBeamZ", &beamZ, "genBeamZ/F");
    treeAnalysis2[f]->Branch("genBeamA", &beamA, "genBeamA/F");
    treeAnalysis2[f]->Branch("genBeamB", &beamB, "genBeamB/F");
    treeAnalysis2[f]->Branch("genBeamEnergy", &energyKinProj, "genBeamEnergy/F");
    treeAnalysis2[f]->Branch("genBeamTheta", &beamTheta, "genBeamTheta/F");
    treeAnalysis2[f]->Branch("genBeamPhi", &beamPhi, "genBeamPhi/F");
    treeAnalysis2[f]->Branch("genExcEn", &genExcEn, "genExcEn/F");
    treeAnalysis2[f]->Branch("genState", &genState, "genState/I");

    // write simulated data to tree
    // these values are with resolutions!!!!!!!!!!!!!
    sprintf(tmpName, "detectorEnergy[%d]/D", maxDetectors);
    treeAnalysis2[f]->Branch("detectorEnergyLoss", detEnergyLoss, tmpName);
    treeAnalysis2[f]->Branch("detectorEnergyLossNotSmeared", detEnergyLossNotSmeared, tmpName);
    sprintf(tmpName, "detectorStripX[%d]/I", maxDetectors);
    treeAnalysis2[f]->Branch("detectorStripX", detStripX, tmpName);
    sprintf(tmpName, "detectorStripY[%d]/I", maxDetectors);
    treeAnalysis2[f]->Branch("detectorStripY", detStripY, tmpName);



    treeAnalysis2[f]->Branch("simLightKinEnergy", &energyKinLight, "simLightEnergy/D"); // sum of detector energy losses

    treeAnalysis2[f]->Branch("simFIx", &FIx, "simFIx/D");
    treeAnalysis2[f]->Branch("simFIy", &FIy, "simFIy/D");
    treeAnalysis2[f]->Branch("simFIz", &FIz, "simFIz/D");
    treeAnalysis2[f]->Branch("simFIdetID", &FIdetID, "simFIdetID/I");
    treeAnalysis2[f]->Branch("anaFIdetID", &firstDetID, "anaFIdetID/I");

    treeAnalysis2[f]->Branch("simLightThetaLab", &thetaLightLab, "simLightThetaLab/D");
    treeAnalysis2[f]->Branch("simLightThetaCM", &thetaLightCM, "simLightThetaCM/D");                           // !!!!!!!!!!!!!!!!!!!!!!!!!!!!! needs bug fixing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    treeAnalysis2[f]->Branch("simLightPhi", &phiLight, "simLightPhi/D");


    // new analysis data

    treeAnalysis2[f]->Branch("anaDetectorHitMul", &detHitMul, "anaDetectorHitMul/I");
    sprintf(tmpName, "anaDetectorHitPos[3]/D");
    treeAnalysis2[f]->Branch("anaDetectorHitPos", simDetectorHitPos, tmpName);

    treeAnalysis2[f]->Branch("anaType", &type, "anaType/I");
    
    treeAnalysis2[f]->Branch("anaMissingMass", &miss, "anaMissingMass/D");
    //treeAnalysis2[f]->Branch("anaProjGamma", &gammaProj, "anaProjGamma/F");
    //treeAnalysis2[f]->Branch("anaProjTotalEnergy", &energyTotProj, "anaProjTotalEnergy/F");
    //treeAnalysis2[f]->Branch("anaProjMomentum", &momentumProj, "anaProjMomentum/F");
    //treeAnalysis2[f]->Branch("anaLightGamma", &gammaLight, "anaLightGamma/F");
    //treeAnalysis2[f]->Branch("anaLightTotalEnergy", &energyTotLight, "anaLightTotalEnergy/F");
    //treeAnalysis2[f]->Branch("anaLightMomentum", &momentumLight, "anaLightMomentum/F");
    //treeAnalysis2[f]->Branch("anaHeavyGamma", &gammaHeavy, "anaHeavyGamma/F");
    //treeAnalysis2[f]->Branch("anaHeavyTotalEnergy", &energyTotHeavy, "anaHeavyTotalEnergy/F");
    //treeAnalysis2[f]->Branch("anaHeavyMomentum", &momentumHeavy, "anaHeavyMomentum/F");
    
    treeAnalysis2[f]->Branch("energyGamma", &energyGamma, "energyGamma/D");
    treeAnalysis2[f]->Branch("energyGammaDC", &energyGammaDC, "energyGammaDC/D");

  } // treeAnalysis2


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





void Analysis::CreateHeader(){
  
  cout << "Creating header with following information about the nuclei:" << endl;

  fileAnalysis->cd(); 
  
  
  projA = info->fProjA;
  projZ = info->fProjZ;
  targA = info->fTargetA;
  targZ = info->fTargetZ;
 
  
  char* massFile = (char*)"/home/philipp/programme/makeEvents/mass.dat";
  
  nucProj = new Nucleus(projZ, projA-projZ, massFile);
  nucTarg = new Nucleus(targZ, targA-targZ, massFile);
  
  massProj = nucProj->GetMass();
  massTarget = nucTarg->GetMass(); 
  
  for(Int_t f=0; f<maxCutType; f++){
    switch(f) {
      case 0: { // elastic scattering
        recoA = targA;
        recoZ = targZ;
        ejecA = projA;
        ejecZ = projZ;
        break;
      }
      case 1: { // (d,p) one neutron transfer
        recoA = targA-1;
        recoZ = targZ;
        ejecA = projA+1;
        ejecZ = projZ;
        break;
      }
      case 2: { // (t,p) one neutron transfer
        recoA = targA-2;
        recoZ = targZ;
        ejecA = projA+2;
        ejecZ = projZ;
        break;
      }
    }
    
    nucReco[f] = new Nucleus(recoZ, recoA-recoZ, massFile);
    nucEjec[f] = new Nucleus(ejecZ, ejecA-ejecZ, massFile);

    massLight = nucReco[f]->GetMass();
    massHeavy = nucEjec[f]->GetMass();

    //qValue = nucProj->GetMass() + nucTarg->GetMass() - nucReco[f]->GetMass() - nucEjec[f]->GetMass() ;
    qValue = massProj + massTarget - massLight - massHeavy; 
    
    treeAnaHeader[f] = new TTree(Form("treeAnaHeader_channel%d", f), Form("Header, channel %d", f));

    treeAnaHeader[f]->Branch("projA", &projA, "projA/I");
    treeAnaHeader[f]->Branch("projZ", &projZ, "projZ/I");
    treeAnaHeader[f]->Branch("projMass", &massProj, "projMass/F");
    treeAnaHeader[f]->Branch("targA", &targA, "targA/I");
    treeAnaHeader[f]->Branch("targZ", &targZ, "targZ/I");
    treeAnaHeader[f]->Branch("targetMass", &massTarget, "targetMass/F");
    treeAnaHeader[f]->Branch("recoA", &recoA, "recoA/I");
    treeAnaHeader[f]->Branch("recoZ", &recoZ, "recoZ/I");
    treeAnaHeader[f]->Branch("recoMass", &massLight, "recoMass/F");
    treeAnaHeader[f]->Branch("ejecA", &ejecA, "ejecA/I");
    treeAnaHeader[f]->Branch("ejecZ", &ejecZ, "ejecZ/I");
    treeAnaHeader[f]->Branch("ejecMass", &massHeavy, "ejecMass/F");
    treeAnaHeader[f]->Branch("qValue", &qValue, "qValue/F");

    treeAnaHeader[f]->Fill();
    treeAnaHeader[f]->Write(Form("anaHeader%d", f));
    
    printf("Wrote for reaction channel/type %d: projectile A/Z %d/%d, target A/Z %d/%d, recoiled A/Z %d/%d, ejectile A/Z %d/%d\n", f, projA,projZ, targA,targZ, recoA,recoZ, ejecA,ejecZ);
    printf("Wrote for reaction channel/type %d: mass projectile %f, mass target %f, mass recoiled %f, mass ejectile %f, q-value %f\n", f, massProj, massTarget, massLight, massHeavy, qValue);


  }

}





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
  TH2F* hdEE_A=new TH2F("hdEE_A_analysis1", "delta E vs. E, all directions, Analysis1 ", 1000,0,25,100,0,8);
  TH2F* hdEE_B=new TH2F("hdEE_B_analysis1", "delta E vs. E, backward direction, Analysis1 ", 1000,0,25,100,0,8);
  TH2F* hdEE_F=new TH2F("hdEE_F_analysis1", "delta E vs. E, forward direction, Analysis1 ", 1000,0,25,100,0,8);


  Int_t goodEvents=0;
  Int_t nevents=tree->GetEntries();

  cout << "Analysis1: start event loop (" << nevents << " events)" << endl; 
  for(Int_t e=0; e<nevents; e++){
    
    if(e%100000==0){
      cout << "Processing event " << e << endl;
    }

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


    hdEE_A->Fill(energyKinLight,detEnergyLoss[firstDetID]);
    if(recoPosZ[firstDetID] < 0.0){ //backward
      hdEE_B->Fill(energyKinLight,detEnergyLoss[firstDetID]);
    }
    if(recoPosZ[firstDetID] > 0.0){ //forward
      hdEE_F->Fill(energyKinLight,detEnergyLoss[firstDetID]);
    }

    treeAnalysis1->Fill();


  } // end of event loop


  cout << nevents << " events processed! " << goodEvents << " events used in analysis! Ratio: " << (Float_t)goodEvents/(Float_t)nevents*100.0 << "%" << endl;
  cout << detHitMul3 << " events with detector hit multiplicity 3 or larger (angle reconstruction might not be correct for these events)" << endl;

  fileAnalysis->cd();
  treeAnalysis1->Write("analysis1");

  hdEE_A->Write("dEE1_A");
  hdEE_B->Write("dEE1_B");
  hdEE_F->Write("dEE1_F");

  printf("Analysis1: Analyzed events writen to file '%s'\n\n", info->fOutFileNameAnalysis);



} // Analysis1















void Analysis::Analysis2(){
 

  Bool_t doAnalysisType[maxCutType] = {false};
  //for(Int_t f=0; f<maxCutType; f++){
  //  doAnalysisType[f] = false;
  //}
  Int_t goodEvents=0;
  Int_t multipleTypeEvents=0;
  Int_t nevents=tree->GetEntries();
  
  // reset the aux vaiable 
  detHitMul3=0;
  
  cout << "Analysis2: start event loop (" << nevents << " events)" << endl; 
  for(Int_t e=0; e<nevents; e++){
    
    if(e%100000==0){
      cout << "Processing event " << e << endl;
    }
     
    ResetVariables();
    for(Int_t f=0; f<maxCutType; f++){
      doAnalysisType[f] = false;
    }


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
    
    
    // select reaction channel from grahical cuts
    for(Int_t f=0; f<maxCutType; f++){
      for(Int_t i=0; i<maxCuts; i++){
        if(cutExists[f][i]){
          if(cut[f][i]->IsInside(energyKinLight, detEnergyLoss[firstDetID])){
            doAnalysisType[f] = true;
            type = f;
          }
        }
      }

      for(Int_t ff=0; ff<f; ff++){
        if(doAnalysisType[f] && doAnalysisType[ff]){
          multipleTypeEvents++;
          type = warningType; 
        }
      }
    }
    
    
    for(Int_t f=0; f<maxCutType; f++){
      if(doAnalysisType[f]){
        MissingMass(f);
        if(type!=warningType){
          goodEvents++;
        }
      }
    }

    
  
  } // event loop

  cout << nevents << " events processed! " << goodEvents << " good events (unambuguious dE-E cut)! Ratio: " << (Float_t)goodEvents/(Float_t)nevents*100.0 << "%" << endl;
  cout << multipleTypeEvents << " multiple type events (within dE-E cuts of two reaction channels). Both channels were analyzed!" << endl;
  cout << detHitMul3 << " events with detector hit multiplicity 3 or larger (angle reconstruction might not be correct for these events)" << endl;



  fileAnalysis->cd();
  for(Int_t f=0; f<maxCutType; f++){
    treeAnalysis2[f]->Write(Form("analysis2_%d", f));
  }

  hdEE_A->Write("dEE2_A");
  hdEE_B->Write("dEE2_B");
  hdEE_F->Write("dEE2_F");
  
  hMiss->Write("missingMass");
  hEth->Write("Eth");
  hThetaLab->Write("hThetaLab");
  hThetaCM->Write("hThetaCM");

  printf("Analysis2: Analyzed events writen to file '%s'\n\n", info->fOutFileNameAnalysis);



} // Analysis2






  // with Kinematics class 
void Analysis::MissingMass(Int_t channel){  

  
  // load masses from header
  treeAnaHeader[channel]->GetEvent(0);

  
  // projectile kinematics
  energyKinProj*=(Float_t)projA; 

  // smear out data with detector position resolutions
  if(info->NoBeamTracking()){
    beamX=0.0;
    beamY=0.0;
    beamZ=0.0;
    beamA=0.0;
    beamB=0.0;
    energyKinProj=info->fBeamEnergy*(Float_t)projA;
  }else{
    beamX=randomizer->Gaus(beamX, info->fResTargetX); // in mm
    beamY=randomizer->Gaus(beamY, info->fResTargetY);
    beamZ=randomizer->Gaus(beamZ, info->fResTargetZ);
    beamA=randomizer->Gaus(beamA, info->fResTargetA);
    beamB=randomizer->Gaus(beamB, info->fResTargetB);
    energyKinProj = randomizer->Gaus(energyKinProj, info->fResBeamE);
  }


   
//  TVector3 vProj(0.0, 0.0, 1.0); 
//  if(!info->NoBeamTracking()){
//    vProj.SetMagThetaPhi(1.0, beamTheta, beamPhi);                                         // comment out this line to see the effect of no beam profile correction
//  
//    // rotate by beam angular resolution
//    vProj.RotateY(randomizer->Gaus(0.0, (info->fResTargetA)/1000.0)); // resolutions in mrad
//    vProj.RotateX(randomizer->Gaus(0.0, (info->fResTargetB)/1000.0));
//  }
 
  energyTotLight = energyKinLight+massLight;

  TVector3 vLight(simDetectorHitPos[0]-beamX, simDetectorHitPos[1]-beamY, simDetectorHitPos[2]-beamZ); // take beam position into account
  
  vLight.RotateY(-beamA/1000.0);
  vLight.RotateX(-beamB/1000.0);
  
  vLight.SetMag(1.0);

  //vLight.SetMagThetaPhi(momentumLight, theta, phi); // without beam position spread


  // energy loss in target, correction
  // roughly
  //if(vLight.Theta()<(110.0*TMath::Pi()/180.0)){
  //  energyKinLight+=0.14;
  //}else{
  //  energyKinLight+=0.14;
  //}
  //energyKinLight+=0.45;

  TLorentzVector lLight(vLight, energyTotLight*1000.0);
  if(lLight.Mag()>0){
    lLight.SetRho( TMath::Sqrt( (energyKinLight+massLight)*(energyKinLight+massLight) - massLight*massLight )*1000 );
  }

  // for the root tree
  thetaLightLab = vLight.Theta(); 
  phiLight = vLight.Phi(); 

  //Kinematics* kine = new Kinematics(nucProj, nucTarg, energyKinProj);
  Kinematics* kine = new Kinematics(nucProj, nucTarg, nucReco[channel], nucEjec[channel], energyKinProj, 0.0);
  //printf("proj %s, targ %s, reco %s, ejec %s, ", nucProj->GetSymbol(), nucTarg->GetSymbol(), nucReco[channel]->GetSymbol(), nucEjec[channel]->GetSymbol());
  //printf("energy kin proj %f ", energyKinProj);
  
  //kine->SetAngles(thetaLightLab, 2, 0);
  kine->Final(thetaLightLab, 2);
  
  miss = -kine->GetExcEnergy(lLight)/1000.0;
  
  thetaLightCM = -kine->GetThetacm(2) + TMath::Pi();


  hThetaLab->Fill(thetaLightLab*180.0/TMath::Pi());
  hThetaCM->Fill(thetaLightCM*180.0/TMath::Pi());


  // fill histograms

  hMiss->Fill(miss);
  //hMissTheta->Fill(vLight.Theta()*180.0/TMath::Pi(), miss);

  hdEE_A->Fill(energyKinLight,detEnergyLoss[firstDetID]);
  if(simDetectorHitPos[2] < 0.0){
    hdEE_B->Fill(energyKinLight,detEnergyLoss[firstDetID]);
  }
  if(simDetectorHitPos[2] > 0.0){
    hdEE_F->Fill(energyKinLight,detEnergyLoss[firstDetID]);
  }

  hEth->Fill(thetaLightLab*180.0/TMath::Pi(), energyKinLight);
  

  
  // Doppler correction
  if(grapeDet>-1){
    Float_t projGamma = (energyKinProj / massProj) + 1.0;
    Float_t projBeta = TMath::Sqrt(1.0 - 1.0/(projGamma*projGamma));

    Float_t th=0.0;
    if(grapeDet<6) th=70.0;
    if(grapeDet>5 && grapeDet<12) th=90.0;
    if(grapeDet>11) th=110.0;

    energyGammaDC = energyGamma * projGamma * (1.0 - (projBeta * TMath::Cos((th/180.0)*TMath::Pi())));



  }else{
    energyGammaDC=NAN;
  }








  energyKinProj/=(Float_t)projA; // AMeV

  treeAnalysis2[channel]->Fill();

  delete kine;

} // MissingMass
























//  // with Lorentz boost by my own
//void Analysis::MissingMass(Int_t channel){  
//
////  TTree* treeHeader=(TTree*)fileBeam->Get("header");
////  if(!treeHeader){
////    cout << "TTree 'header' not found in makeEvents root file!" << endl;
////    //return 0;
////  }
////
////
////  Int_t projA=132;
////  //Float_t projMass=0.0, targetMass=0.0, lightMass=0.0, heavyMass=0.0, qValue=0.0;
////  Float_t massProj=0.0, massTarget=0.0, massLight=0.0, massHeavy=0.0, qValue=0.0;
////
////  treeHeader->SetBranchAddress("projA", &projA);
////  treeHeader->SetBranchAddress("projMass", &massProj);
////  treeHeader->SetBranchAddress("targetMass", &massTarget);
////  treeHeader->SetBranchAddress("lightMass", &massLight);
////  treeHeader->SetBranchAddress("heavyMass", &massHeavy);
////  treeHeader->SetBranchAddress("qValue", &qValue);
////
////  treeHeader->GetEntry(0); // todo: take these values from Kinematics/Nucleus/mass file
////  //printf("Obtained masses: projectile %f, target %f, light ejectile %f, heavy ejectile %f; Q-value %f\n", massProj, massTarget, massLight, massHeavy, qValue);
//
//  
//  // load masses from header
//  treeAnaHeader[channel]->GetEvent(0);
//
//
//    // aux:
//    //Double_t x = FIx;
//    //Double_t y = FIy;
//    //Double_t z = FIz;
//
//    // take only events in backward direction
//    // in very forward direction are usually punch-troughs of protons through the detector
//  //    if(theta<90.0){
//  //    if(theta==0.0){
//  //      continue;
//  //    }
//      
//  //    theta/=180.0/TMath::Pi(); // is actually not used anymore
//    
//    // smear out data with detector position resolutions
//    if(info->NoBeamTracking()){
//      beamX=0.0;
//      beamY=0.0;
//      beamZ=0.0;
//    }else{
//      beamX=randomizer->Gaus(beamX, info->fResTargetX); // in mm
//      beamY=randomizer->Gaus(beamY, info->fResTargetY);
//      beamZ=randomizer->Gaus(beamZ, info->fResTargetZ);
//    }
//    
//
//
//    // projectile kinematics
//
//    energyKinProj*=(Float_t)projA; 
//    // Projectile data
//    // at the moment from simulation input
//    // todo: separate simulation including incoming tracking
//
//    if(!info->NoBeamTracking()){
//      energyKinProj = randomizer->Gaus(energyKinProj, info->fResBeamE);
//    }
//      
//    gammaProj = (energyKinProj)/massProj + 1.0;
//    //Float_t betaProj = TMath::Sqrt(1.0-(1.0/(gammaProj*gammaProj))); //just for cross checking
//    momentumProj = massProj*TMath::Sqrt(gammaProj*gammaProj-1.0); 
//    energyTotProj = massProj*gammaProj; //total energy
//    
//     
//    TVector3 vProj(0.0, 0.0, 1.0); 
//    if(!info->NoBeamTracking()){
//      vProj.SetMagThetaPhi(1.0, beamTheta, beamPhi);                                         // comment out this line to see the effect of no beam profile correction
//    
//      // rotate by beam angular resolution
//      vProj.RotateY(randomizer->Gaus(0.0, (info->fResTargetA)/1000.0)); // resolutions in mrad
//      vProj.RotateX(randomizer->Gaus(0.0, (info->fResTargetB)/1000.0));
//    }
//    vProj.SetMag(momentumProj);
//
//    // center of mass kinematic values
//    Float_t energyCm = TMath::Sqrt(massProj*massProj + massTarget*massTarget + 2.0*energyTotProj*massTarget);
//    Float_t betaCm = momentumProj/(energyTotProj+massTarget);
//    TVector3 vCm(vProj); // direction of projectile including beam profile
//    vCm.SetMag(betaCm);
//
//
//
//
//    // light ejectile kinematics
//
//    // get total energy and momentum of the light ejectile
//    
//    //energySum=energyLoss+energyTotal;
//    //gammaLight = energySum/massLight+1.0;   
//    gammaLight = energyKinLight/massLight+1.0;      // simulated
//    //gammaLight = genLightEnergy/massLight+1.0; // generated
//    //theta=genLightTheta; // generated theta
//   
//   
//    energyTotLight = gammaLight*massLight;
//    //energyTotLight = energyKinLight+massLight;
//    momentumLight = massLight*TMath::Sqrt(gammaLight*gammaLight-1.0);
//
//    //printf("lightTheta %f, lightEnergy %f \n", theta, energySum);
//    //printf("x %f, y %f, z %f, beamX %f, beamY %f, beamZ %f\n",x ,y, z, beamX, beamY, beamZ);
//    //TVector3 vLight(x-beamX, y-beamY, z-beamZ); //momentum direction of proton
//    TVector3 vLight(simDetectorHitPos[0]-beamX, simDetectorHitPos[1]-beamY, simDetectorHitPos[2]-beamZ); //momentum direction of proton
//    //printf("vLight.Mag %f\n", vLight.Mag());
//    vLight.SetMag(momentumLight);
//    //vLight.SetMagThetaPhi(momentumLight, theta, phi); // without beam position spread
//
//
//
//    // for the root tree
//    thetaLightLab = vLight.Theta(); 
//    phiLight = vLight.Phi(); 
//
//
//    TLorentzVector lLight;
//    lLight.SetVect(vLight);
//    lLight.SetE(energyTotLight);
//
////    TLorentzVector lL=lLight;
//
//    lLight.Boost(-vCm);
//    thetaLightCM = TMath::Pi() - lLight.Theta();
//
//
//
////vLight.SetTheta(-vLight.Theta() + TMath::Pi());
////lL.SetTheta(vLight.Theta() - TMath::Pi());
//
////    lL.Boost(-vCm);
////    thetaLightCM=lL.Theta();
//
//
//    hThetaLab->Fill(thetaLightLab*180.0/TMath::Pi());
//    hThetaCM->Fill(thetaLightCM*180.0/TMath::Pi());
//
//
//
//    // heavy ejectile kinematics in center of mass system
//    TLorentzVector lHeavy;
//    lHeavy.SetVect(-lLight.Vect());
//    lHeavy.SetE(energyCm-lLight.E());
//
//    miss = -lHeavy.M()+massHeavy;
//
//    //printf("miss %f\n", miss);
//    
//
//    // fill histograms
//
//    hMiss->Fill(miss);
//    //hMissTheta->Fill(vLight.Theta()*180.0/TMath::Pi(), miss);
//
//    hdEE->Fill(energyKinLight,detEnergyLoss[firstDetID]);
//    //hEth->Fill(thetaLightLab*180.0/TMath::Pi(), energyKinLight);
//
//
//    energyKinProj/=(Float_t)projA; // MeV/u`
//
//    treeAnalysis2[channel]->Fill();
//
//} // MissingMass


























Bool_t Analysis::GetCuts(){
  
  Int_t filesFound=0;
  TFile* file[maxCutType];
  Bool_t haveACut = false;
  char tmpName[100];

  for(Int_t f=0; f<maxCutType; f++){
    
    if(strcmp(info->fFileNameCuts[f],"")==0){
      continue;
    }else{
      cout << "Opening file '" << info->fFileNameCuts[f] << "'" << endl;
      file[f] = TFile::Open(info->fFileNameCuts[f],"read");
      if(!file[f]){
        cout << "Cuts file not found: " << info->fFileNameCuts[f] << endl;
        return false;
      }
    }

    filesFound++;

    TIter next(file[f]->GetListOfKeys());
    TKey *key;
    Int_t nkeys = file[f]->GetListOfKeys()->GetSize();
    printf("Found %d keys\n", nkeys);
    //while ((key = (TKey*)next())) {
    for(Int_t i=0; i<nkeys; i++){
      key = (TKey*)next();
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TCutG")) {continue;}
      cut[f][i] = (TCutG*)key->ReadObj();
      sprintf(tmpName, "%s", cut[f][i]->GetName());
      printf("Found cut named '%s', gave index %d\n", tmpName, i);

      cut[f][i]->SetTitle(Form("dE-E cut, type %d, number %d, original name '%s'", f, i, tmpName));
      cutExists[f][i] = true;
      if(!haveACut){haveACut = true;}

      fileAnalysis->cd();
      cut[f][i]->Write(Form("cut_%d_%d", f, i));

    } // loop over keys
  } // loop over files



  return haveACut;


}




