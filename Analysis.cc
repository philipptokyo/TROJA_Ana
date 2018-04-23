#include "Analysis.hh"

#include "Nucleus.hh"
#include "Compound.hh"
#include "Reconstruction.hh"

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

  fEnLoss = new EnLoss();
  fEnLoss->CollectData("/home/philipp/analysis/troja/EnLossData/enLoss_tables_p_in_CD2.csv", 0);

  if(detInfo->IncludeDali()){
    for(Int_t i=0; i<NUMBEROFDALI2CRYSTALS; i++){
    for(Int_t j=0; j<NUMBEROFDALI2CRYSTALS; j++){
      daliAddbackTable[i][j]=0;
    }
    }
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
  betaProj=0.0;
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
  miss=NAN;

  //energyGamma=0.0;
  //grapeDet=-1;
  
  grapeDetMul=0;
  grapeSumEnergyDC=0.0;
  for(Int_t d=0; d<grapeMaxDet; d++){
    grapeDetEnergy[d]=NAN;
    grapeDetEnergyDC[d]=NAN;
    grapeCryMul[d]=0;
    for(Int_t c=0; c<grapeMaxCry; c++){
      grapeCryEnergy[d][c]=NAN;
      grapeCryEnergyDC[d][c]=NAN;
      grapeSegMul[d][c]=0;
      for(Int_t s=0; s<grapeMaxSeg; s++){
        grapeSegEnergy[d][c][s]=NAN;
        grapeSegEnergyDC[d][c][s]=NAN;
      }
    }
  }

  daliCrystalMult=0;
  daliAddbackMult=0;
  daliEnergyDopplerSum=NAN;
  for(Int_t d=0; d<NUMBEROFDALI2CRYSTALS; d++){

    daliCrystalFlag[d]=false;
    daliCrystalEnergy[d]=NAN;
    daliCrystalTime[d]=NAN;
    daliFITime[d]=NAN;
    daliFIX[d]=NAN;
    daliFIY[d]=NAN;
    daliFIZ[d]=NAN;

    daliEnergy[d]=NAN; // with resolution
    daliEnergyDoppler[d]=NAN; // with resolution, doppler corrected

  }

  for(Int_t d=0; d<NUMBEROFDALI2ADDBACKCRYSTALS; d++){
    daliAddbackCrystalMult[d] = 0;
    daliAddbackFlag[d]=false;
    daliAddbackEnergy[d]=NAN; // addbacked, ith resolution
    daliAddbackEnergyDoppler[d]=NAN; // addbacked, with resolution, doppler corrected
    daliCrystalUsedForAddback[d]=false;
  }


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

  MakeSplineEnAfter2EnLoss();
  
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

  treeHeader = (TTree*)infile->Get("header");
  //TTree* tree=(TTree*)infile->Get("events"); //simulation input
  if(!treeHeader){
    cout << "TTree 'header' not found in geant root file! That can be ignored, but some anaysis steps will not work (e.g. Dali Doppler correction)" << endl;
    //return false;
  }else{
    treeHeader->SetBranchAddress("daliPosX", (detInfo->daliHead.fDaliPosX));
    treeHeader->SetBranchAddress("daliPosY", (detInfo->daliHead.fDaliPosY));
    treeHeader->SetBranchAddress("daliPosZ", (detInfo->daliHead.fDaliPosZ));
    treeHeader->SetBranchAddress("daliTheta", (detInfo->daliHead.fDaliTheta));
    treeHeader->SetBranchAddress("daliPhi", (detInfo->daliHead.fDaliPhi));
    treeHeader->SetBranchAddress("daliDistance", (detInfo->daliHead.fDaliDistance));
    treeHeader->SetBranchAddress("daliTimeResolution", (detInfo->daliHead.fDaliTimeResolution));
    treeHeader->SetBranchAddress("daliEnergyResolutionInd", (detInfo->daliHead.fDaliEnergyResolutionInd));

    treeHeader->GetEvent(0);
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
  
  treeBeam->SetBranchAddress("gammaMul", &genGammaMul);
  treeBeam->SetBranchAddress("gammaERest", genGammaERest);
  treeBeam->SetBranchAddress("gammaE", genGammaELab);

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
  
  //tree->SetBranchAddress("grapeEnergy", &energyGamma);
  //tree->SetBranchAddress("grapeDetector", &grapeDet);
  
  tree->SetBranchAddress("targetEnergyLoss", &targetEnergyLoss);
  
  if(detInfo->IncludeGrape()){
    tree->SetBranchAddress("grapeDetMul", &grapeDetMul); // ID of the detector with the first interaction point

    tree->SetBranchAddress("grapeCryMul", grapeCryMul); // ID of the detector with the first interaction point
    tree->SetBranchAddress("grapeSegMul", grapeSegMul); // ID of the detector with the first interaction point

    tree->SetBranchAddress("grapeDetEnergy", grapeDetEnergy); // ID of the detector with the first interaction point
    tree->SetBranchAddress("grapeCryEnergy", grapeCryEnergy); // ID of the detector with the first interaction point
    tree->SetBranchAddress("grapeSegEnergy", grapeSegEnergy); // ID of the detector with the first interaction point
  }

  if(detInfo->IncludeDali()){
    
    tree->SetBranchAddress("DALI2Flag",          daliCrystalFlag);
    tree->SetBranchAddress("DALI2EnergyNotCor",  daliCrystalEnergy);
    tree->SetBranchAddress("DALI2Mult",          &daliCrystalMult);
    tree->SetBranchAddress("DALI2Time",          daliCrystalTime);
    tree->SetBranchAddress("DALI2FITime",        daliFITime);
    tree->SetBranchAddress("DALI2FIX",           daliFIX);
    tree->SetBranchAddress("DALI2FIY",           daliFIY);
    tree->SetBranchAddress("DALI2FIZ",           daliFIZ);

  }


  fileAnalysis = new TFile(info->fOutFileNameAnalysis, "recreate");
  //fEnAfter2EnLoss->Write("EnergyAfter2EnergyLoss");
  fileAnalysis->mkdir("EnergyAfter2EnergyLoss");
  fileAnalysis->cd("EnergyAfter2EnergyLoss");
  for(Int_t s=0; s<maxEnLossSplines; s++){
    fEnAfter2EnLoss[s]->Write(Form("spline%03dmum", s));
  }
  fileAnalysis->cd();

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
  
  //treeAnalysis1->Branch("energyGamma", &energyGamma, "energyGamma/D");
  //treeAnalysis1->Branch("energyGammaDC", &energyGammaDC, "energyGammaDC/D");



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

    treeAnalysis2[f]->Branch("genGammaMul", &genGammaMul, "genGammaMul/I");
    sprintf(tmpName, "genGammaERest[%d]/F", maxGammaMulGen);
    treeAnalysis2[f]->Branch("genGammaERest", genGammaERest, tmpName);
    sprintf(tmpName, "genGammaELab[%d]/F", maxGammaMulGen);
    treeAnalysis2[f]->Branch("genGammaELab", genGammaELab, tmpName);
    
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
    treeAnalysis2[f]->Branch("simLightKinEnergyUncorr", &energyKinLightUncorr, "simLightEnergyUncorr/D"); // sum of detector energy losses
    treeAnalysis2[f]->Branch("simTargetEnergyLoss", &targetEnergyLoss, "simTargetEnergyLoss/D"); // sum of detector energy losses

    treeAnalysis2[f]->Branch("simFIx", &FIx, "simFIx/D");
    treeAnalysis2[f]->Branch("simFIy", &FIy, "simFIy/D");
    treeAnalysis2[f]->Branch("simFIz", &FIz, "simFIz/D");
    treeAnalysis2[f]->Branch("simFIdetID", &FIdetID, "simFIdetID/I");
    treeAnalysis2[f]->Branch("anaFIdetID", &firstDetID, "anaFIdetID/I");

    treeAnalysis2[f]->Branch("simLightThetaLab", &thetaLightLab, "simLightThetaLab/D");
    treeAnalysis2[f]->Branch("simLightThetaCM", &thetaLightCM, "simLightThetaCM/D");     
    treeAnalysis2[f]->Branch("simLightPhi", &phiLight, "simLightPhi/D");
    
    
    treeAnalysis2[f]->Branch("anaProjBeta", &betaProj, "anaProjBeta/D");
    treeAnalysis2[f]->Branch("anaProjGamma", &gammaProj, "anaProjGamma/D");


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
    
    //treeAnalysis2[f]->Branch("energyGamma", &energyGamma, "energyGamma/D");
    //treeAnalysis2[f]->Branch("energyGammaDC", &energyGammaDC, "energyGammaDC/D");
   
    if(detInfo->IncludeGrape()){
      treeAnalysis2[f]->Branch("grapeDetMul", &grapeDetMul, "grapeDetMul/I"); // ID of the detector with the first interaction point

      sprintf(tmpName, "grapeCryMul[%d]/I", grapeMaxDet);
      treeAnalysis2[f]->Branch("grapeCryMul", grapeCryMul, tmpName); // ID of the detector with the first interaction point
      sprintf(tmpName, "grapeSegMul[%d][%d]/I", grapeMaxDet, grapeMaxCry);
      treeAnalysis2[f]->Branch("grapeSegMul", grapeSegMul, tmpName); // ID of the detector with the first interaction point
  
      sprintf(tmpName, "grapeDetEnergy[%d]/D", grapeMaxDet);
      treeAnalysis2[f]->Branch("grapeDetEnergy", grapeDetEnergy, tmpName); // ID of the detector with the first interaction point
      sprintf(tmpName, "grapeCryEnergy[%d][%d]/D", grapeMaxDet, grapeMaxCry);
      treeAnalysis2[f]->Branch("grapeCryEnergy", grapeCryEnergy, tmpName); // ID of the detector with the first interaction point
      sprintf(tmpName, "grapeSegEnergy[%d][%d][%d]/D", grapeMaxDet, grapeMaxCry, grapeMaxSeg);
      treeAnalysis2[f]->Branch("grapeSegEnergy", grapeSegEnergy, tmpName); // ID of the detector with the first interaction point
  
      treeAnalysis2[f]->Branch("grapeSumEnergyDC", &grapeSumEnergyDC, "grapeSumEnergyDC/D"); // ID of the detector with the first interaction point
      sprintf(tmpName, "grapeDetEnergyDC[%d]/D", grapeMaxDet);
      treeAnalysis2[f]->Branch("grapeDetEnergyDC", grapeDetEnergyDC, tmpName); // ID of the detector with the first interaction point
      sprintf(tmpName, "grapeCryEnergyDC[%d][%d]/D", grapeMaxDet, grapeMaxCry);
      treeAnalysis2[f]->Branch("grapeCryEnergyDC", grapeCryEnergyDC, tmpName); // ID of the detector with the first interaction point
      sprintf(tmpName, "grapeSegEnergyDC[%d][%d][%d]/D", grapeMaxDet, grapeMaxCry, grapeMaxSeg);
      treeAnalysis2[f]->Branch("grapeSegEnergyDC", grapeSegEnergyDC, tmpName); // ID of the detector with the first interaction point
    }


    if(detInfo->IncludeDali()){
        
      sprintf(tmpName, "DALI2Mult/I");
      treeAnalysis2[f]->Branch("DALI2Mult", &(daliCrystalMult), tmpName);
      sprintf(tmpName, "DALI2Flag[%d]/O", NUMBEROFDALI2CRYSTALS);
      treeAnalysis2[f]->Branch("DALI2Flag", (daliCrystalFlag), tmpName);
      sprintf(tmpName, "DALI2EnergyNotCor[%d]/F", NUMBEROFDALI2CRYSTALS);
      treeAnalysis2[f]->Branch("DALI2EnergyNotCor", (daliCrystalEnergy), tmpName);
      sprintf(tmpName, "DALI2Time[%d]/F", NUMBEROFDALI2CRYSTALS);
      treeAnalysis2[f]->Branch("DALI2Time", (daliCrystalTime), tmpName);
      sprintf(tmpName, "DALI2FITime[%d]/D", NUMBEROFDALI2CRYSTALS);
      treeAnalysis2[f]->Branch("DALI2FITime", (daliFITime), tmpName);
      sprintf(tmpName, "DALI2FIX[%d]/F", NUMBEROFDALI2CRYSTALS);
      treeAnalysis2[f]->Branch("DALI2FIX", (daliFIX), tmpName);
      sprintf(tmpName, "DALI2FIY[%d]/F", NUMBEROFDALI2CRYSTALS);
      treeAnalysis2[f]->Branch("DALI2FIY", (daliFIY), tmpName);
      sprintf(tmpName, "DALI2FIZ[%d]/F", NUMBEROFDALI2CRYSTALS);
      treeAnalysis2[f]->Branch("DALI2FIZ", (daliFIZ), tmpName);


      sprintf(tmpName, "DALI2Energy[%d]/F", NUMBEROFDALI2CRYSTALS);
      treeAnalysis2[f]->Branch("DALI2Energy", (daliEnergy), tmpName); // with resolution
      sprintf(tmpName, "DALI2EnergyDoppler[%d]/F", NUMBEROFDALI2CRYSTALS);
      treeAnalysis2[f]->Branch("DALI2EnergyDoppler", (daliEnergyDoppler), tmpName); // with resolution
      treeAnalysis2[f]->Branch("DALI2EnergyDopplerSum", &(daliEnergyDopplerSum), "DALI2EnergyDopplerSum/F"); 
      
      
      sprintf(tmpName, "DALI2AddbackMult/I");
      treeAnalysis2[f]->Branch("DALI2AddbackMult", &(daliAddbackMult), tmpName);
      sprintf(tmpName, "DALI2AddbackCrystalMult[%d]/I", NUMBEROFDALI2ADDBACKCRYSTALS);
      treeAnalysis2[f]->Branch("DALI2AddbackCrystalMult", (daliAddbackCrystalMult), tmpName);
      sprintf(tmpName, "DALI2AddbackFlag[%d]/O", NUMBEROFDALI2CRYSTALS);
      treeAnalysis2[f]->Branch("DALI2AddbackFlag", (daliAddbackFlag), tmpName);
      sprintf(tmpName, "DALI2AddbackEnergy[%d]/F", NUMBEROFDALI2CRYSTALS);
      treeAnalysis2[f]->Branch("DALI2AddbackEnergy", (daliAddbackEnergy), tmpName); 
      sprintf(tmpName, "DALI2AddbackEnergyDoppler[%d]/F", NUMBEROFDALI2CRYSTALS);
      treeAnalysis2[f]->Branch("DALI2AddbackEnergyDoppler", (daliAddbackEnergyDoppler), tmpName); 

    }

  } // treeAnalysis2





  if(detInfo->IncludeDali()){
    

    printf("Reading Dali detector positions\n");
    // get detector positions
    FILE *fFileIn  = fopen("/home/philipp/sim/troja/dali2_geometry_in.txt","r");
    for(Int_t i=0; !feof(fFileIn) && i<NUMBEROFDALI2CRYSTALS; i++)  {
      //fscanf(fFileIn,"%f %f %f %f %f %f %f %i",&x,&y,&z,&psi,&theta,&phi,&rotSign,&detType);
      Float_t trash1[4]={0.0};
      Int_t   trash2=0;
      fscanf(fFileIn,"%lf %lf %lf %f %f %f %f %i", &detInfo->daliHead.fDaliPosX[i], 
                                                &detInfo->daliHead.fDaliPosY[i], 
                                                &detInfo->daliHead.fDaliPosZ[i], 
                                                &trash1[0], &trash1[1], &trash1[2], &trash1[3],
                                                &trash2
                                                );
    }
    fclose(fFileIn);
    


    // create an addback table
    printf("Creating Dali addback table\n");
    FILE *fAddbackTableIn  = fopen("AddbackTable.out","w");
    Float_t dummy[3][2];
    Bool_t inTable;
    Int_t counter = 0;
    for(Int_t i=0;i<NUMBEROFDALI2CRYSTALS;i++)  {
      fprintf(fAddbackTableIn," %i",i);
      for(Int_t j=i+1;j<NUMBEROFDALI2CRYSTALS;j++)  {
        
        dummy[0][0] = detInfo->daliHead.fDaliPosX[i];
        dummy[1][0] = detInfo->daliHead.fDaliPosY[i];
        dummy[2][0] = detInfo->daliHead.fDaliPosZ[i];
        dummy[0][1] = detInfo->daliHead.fDaliPosX[j];
        dummy[1][1] = detInfo->daliHead.fDaliPosY[j];
        dummy[2][1] = detInfo->daliHead.fDaliPosZ[j];
        
        //printf("Distance between %d and %d is ", i, j);
        inTable = IncludeAddbackTable(dummy);  
        if(inTable && counter< NUMBEROFDALI2ADDBACKCRYSTALS) {
          fprintf(fAddbackTableIn," %i",j);
          daliAddbackTable[i][counter] = j;
          counter++;
        }
        if(counter == NUMBEROFDALI2ADDBACKCRYSTALS)  { //Too many detectors 
          cout<<"You have to increase the variable NUMBEROFDALI2ADDBACKCRYSTALS!!!"<<endl;
          abort();
        }
      }
      counter = 0;
      fprintf(fAddbackTableIn," \n");
    }
    fclose(fAddbackTableIn);
    printf("Addback table created\n");
  }




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


void Analysis::MakeSplineEnAfter2EnLoss(){
  printf("Creating ELoss spline\n");
  //Nucleus* projectile, Compound* target, double thickness
  
  //printf("Creating Compound target and proton\n");
  Compound* comTarg = new Compound((char*)"DPE");
  //printf("Compound mass %lf, symbol %s\n", comTarg->GetMass(), comTarg->GetSymbol());
  char* massFile = (char*)"/home/philipp/programme/makeEvents/mass.dat";
  Nucleus* prot = new Nucleus(1, 0, massFile); 
  //printf("Proton: mass %lf, symbol %s\n", prot->GetMass(), prot->GetSymbol());
  
  //printf("Creating Reconstruction\n");
  Double_t thknss = 0.5*0.819*100.0;
  //printf("Target thickness %lf mg/cm2\n", thknss);
  Reconstruction* recons = new Reconstruction(prot, comTarg, thknss*2.0); // mg/cm2, Recon takes half
  //recons->Print(30);
  for(Int_t s=0; s<maxEnLossSplines; s++){
    thknss = s/1000.0*0.819*100.0;
    //thknss = s/1000.0*1.00*100.0;
    recons->SetTargetThickness(thknss*2.0);
    fEnAfter2EnLoss[s] = recons->EnergyAfter2EnergyLoss(30.0, 1.0);
  }
  
  //printf("Getting spline\n");
  //fEnAfter2EnLoss = recons->EnergyAfter2EnergyLoss(30.0, 1.0);
  
  printf("Splines created\n");
}


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

//printf("ProcessDetectorHits\n");

// sum up all energy losses
  for(Int_t d=0; d<maxDetectors; d++){
    detHitMul += detHit[d];
//printf("d %d, detHit[d] %d, detStripX[d] %d, detStripY[d] %d\n", d, detHit[d], detStripX[d], detStripY[d]);
    //if(detHit[d]==1){
    if( (detHit[d]==1) ){

      energyKinLight+=detEnergyLoss[d];

      if(detInfo->IsPosDet(d) && (detStripX[d]>-1) && (detStripY[d]>-1)){
//printf("Detectr %d is pos det\n", d);
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
  }


  // if no position sensitive detector was hit, the event will be skipped
  if(firstDetID<0){
    detHitMul=0;
  }


//if(detHitMul>0){
//  printf("used detector is %d, hitmul is %d\n", firstDetID, detHitMul);
//}


  //simDetectorHitPos[0] = FIx;
  //simDetectorHitPos[1] = FIy;
  //simDetectorHitPos[2] = FIz;
  simDetectorHitPos[0] = recoPosX[firstDetID];
  simDetectorHitPos[1] = recoPosY[firstDetID];
  simDetectorHitPos[2] = recoPosZ[firstDetID];

//printf("position: %f %f %f\n", recoPosX[firstDetID], recoPosY[firstDetID], recoPosZ[firstDetID]);


} // ProcessDetectorHits






void Analysis::Analysis1(){
  // this function is only processing the silicon detector hits
  // for delta E - E plot 
  // after selection in the PID plot, rerun the program with 
  // cuts for analysis2 (missing mass and ggf. gamma processing)
  
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
  //Int_t nevents=tree->GetEntries();
  Int_t nevents=info->fNumberEvents;
  
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

//printf("after ProcessDetectorHits:\ndetHitMul %d\n", detHitMul);
    if(detHitMul==0){
    //if(detHitMul==0 && grapeDetMul==0){ // no hits in silicons and no hits in Gamma Det
      
      continue;
    }

//printf("\n\nEvent %d\n", e);
//printf("detHitMul %d, grapeDetMul %d\n", detHitMul, grapeDetMul);


    if(detHitMul>2){
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
    
    // do the actual missing mass analysis
    for(Int_t f=0; f<maxCutType; f++){
      if(doAnalysisType[f]){
        MissingMass(f);
        //printf("missing mass in loop %f\n", miss);
        if(type!=warningType){
          goodEvents++;
        }
      
        // now, after MissingMass, all information for gamma analysis are available
        if((detInfo->IncludeGrape()) && (grapeDetMul>0)){
          AnalysisGrape();
          //printf("grapeSumEnergyDC %f\n", grapeSumEnergyDC);
        }
        
        if(detInfo->IncludeDali() && (daliCrystalMult>0)){
          //printf("daliCrystalMult %d\n", daliCrystalMult);
          AnalysisDali();
        }

        
        // finally, fill the tree
        //if(doAnalysisType[f] || grapeDetMul>0){
        treeAnalysis2[f]->Fill();
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
  
//printf("missing mass\nposition: %f %f %f\n\n", simDetectorHitPos[0], simDetectorHitPos[1], simDetectorHitPos[2]);

  vLight.RotateY(-beamA/1000.0);
  vLight.RotateX(-beamB/1000.0);
  
  vLight.SetMag(1.0);

  //vLight.SetMagThetaPhi(momentumLight, theta, phi); // without beam position spread

  // for the root tree
  thetaLightLab = vLight.Theta(); 
  phiLight = vLight.Phi(); 


  // energy loss in target, correction
  energyKinLightUncorr=energyKinLight;
  
  //Double_t enloca=fEnLoss->CalcEnLoss(energyKinLight, (detInfo->GetTargetSize(2)*1000.0/2.0)/TMath::Abs(TMath::Cos(vLight.Theta())), 0);
  //energyKinLight+=enloca;
  
  //energyKinLight=fEnLoss->CalcParticleEnergy(energyKinLightUncorr, (detInfo->GetTargetSize(2)*1000.0/2.0)/TMath::Abs(TMath::Cos(vLight.Theta())), 0);

  Int_t pathlength = (Int_t)(detInfo->GetTargetSize(2)*1000.0/2.0)/TMath::Abs(TMath::Cos(vLight.Theta()));
  //printf("pathlength %d\n", pathlength);
  if(pathlength >= maxEnLossSplines){
    printf("No spline for pathlength %d mum available! Increase 'maxEnLossSplines'!\n", pathlength);
    abort();
  }
  energyKinLight += fEnAfter2EnLoss[pathlength]->Eval(energyKinLightUncorr);

  //energyKinLight = energyKinLightUncorr+targetEnergyLoss; // hack for testing the code
  

  TLorentzVector lLight(vLight, energyTotLight*1000.0);
  //printf("lLight.Mag() %f, ", lLight.Mag());
  if(lLight.Mag()>0){
    lLight.SetRho( TMath::Sqrt( (energyKinLight+massLight)*(energyKinLight+massLight) - massLight*massLight )*1000 ); // keV
  }

  //Kinematics* kine = new Kinematics(nucProj, nucTarg, energyKinProj);
  Kinematics* kine = new Kinematics(nucProj, nucTarg, nucReco[channel], nucEjec[channel], energyKinProj, 0.0);
  //printf("proj %s, targ %s, reco %s, ejec %s, ", nucProj->GetSymbol(), nucTarg->GetSymbol(), nucReco[channel]->GetSymbol(), nucEjec[channel]->GetSymbol());
  //printf("energy kin proj %f ", energyKinProj);
  
  //kine->SetAngles(thetaLightLab, 2, 0);
  kine->Final(thetaLightLab, 2);
  
  miss = -kine->GetExcEnergy(lLight)/1000.0; // MeV

  thetaLightCM = -kine->GetThetacm(2) + TMath::Pi();
  
  //printf("reco theta lab %f, kin en %f, missing mass %f\n", thetaLightLab, energyKinLight, miss);

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
  

  

  energyKinProj/=(Float_t)projA; // AMeV

  gammaProj = ((energyKinProj*(Float_t)projA) / massProj) + 1.0;
  betaProj = TMath::Sqrt(1.0 - 1.0/(gammaProj*gammaProj));

  delete kine;

} // MissingMass






void Analysis::AnalysisGrape(){
  
  treeAnaHeader[0]->GetEvent(0); // only projectile information used; is same for all channels
  
  // Doppler correction
    
  //energyKinProj*=(Float_t)projA; 
  
//  Double_t projGamma = ((energyKinProj*(Float_t)projA) / massProj) + 1.0;
//  Double_t projBeta = TMath::Sqrt(1.0 - 1.0/(projGamma*projGamma));

  Double_t projGamma = gammaProj; // todo: remove this
  Double_t projBeta = betaProj;

  Double_t grapeTheta[grapeMaxDet]={
                                    125.0, 125.0, 125.0, 125.0, 125.0, 125.0,
                                     90.0,  90.0,  90.0,  90.0,  90.0,  90.0, 
                                     55.0,  55.0,  55.0,  55.0,  55.0,  55.0 
                                    };
  Double_t df = 0.0;
  //todo: get theta of detector from header, or the GeConfig.txt

  ////proposal:
  //if(grapeDet<6) th=70.0;
  //if(grapeDet>5 && grapeDet<12) th=90.0;
  //if(grapeDet>11) th=110.0;
  
  //// Ota-san's code:
  //if(grapeDet<6) th=55.0;
  //if(grapeDet>5 && grapeDet<12) th=90.0;
  //if(grapeDet>11) th=125.0;

  //energyGammaDC = energyGamma * projGamma * (1.0 - (projBeta * TMath::Cos((th/180.0)*TMath::Pi())));
  
  for(Int_t d=0; d<grapeMaxDet; d++){
    df = projGamma * (1.0 - (projBeta * TMath::Cos((grapeTheta[d]/180.0)*TMath::Pi())));
    grapeDetEnergyDC[d] = grapeDetEnergy[d] * df;

    if(grapeDetEnergyDC[d]>0.0){
      grapeSumEnergyDC += grapeDetEnergyDC[d];
    }

    //if(!TMath::IsNaN(grapeDetEnergyDC[d])){
    //  printf("Proj kin En %f, mass %f, gamma %f, proj beta %f, Gamma E %f, df %f, gamma E dc %f\n", energyKinProj, massProj, projGamma, projBeta, grapeDetEnergy[d], df, grapeDetEnergyDC[d]);
    //  printf("Det %d, theta %f Gamma E %f, df %f, gamma E dc %f\n", d, grapeTheta[d], grapeDetEnergy[d], df, grapeDetEnergyDC[d]);
    //}
    for(Int_t c=0; c<grapeMaxCry; c++){
      grapeCryEnergyDC[d][c] = grapeCryEnergy[d][c] * df;
      for(Int_t s=0; s<grapeMaxSeg; s++){
        grapeSegEnergyDC[d][c][s] = grapeSegEnergy[d][c][s] * df;
      }
    }
  }
  
  //energyKinProj/=(Float_t)projA; 

} // AnalysisGrape



void Analysis::AnalysisDali(){
  
  if(TMath::IsNaN(daliEnergyDopplerSum)){daliEnergyDopplerSum=0.0;}

  Int_t crystalMultDali2=0;
  Int_t dali2_crystalFired[NUMBEROFDALI2ADDBACKCRYSTALS] = {-1};


  //cout << endl;
  for(Int_t d=0; d<NUMBEROFDALI2CRYSTALS; d++){
    
    // printf("Dali crystal %d, flag %d, energy %lf\n", d, daliCrystalFlag[d], daliEnergy[d]);

    if(daliCrystalFlag[d]){
      
      // smear with energy resolution; todo: do this in sim already
      daliEnergy[d] = daliCrystalEnergy[d];
      
      // todo: dali threshold here
      if(daliEnergy[d]>0.001){
        dali2_crystalFired[crystalMultDali2] = d;
        crystalMultDali2++;
      }



      // for doppler correction:
      // get theta of crystal
      TVector3 vdt(
                  detInfo->daliHead.fDaliPosX[d] - beamX,
                  detInfo->daliHead.fDaliPosY[d] - beamY,
                  detInfo->daliHead.fDaliPosZ[d] - beamZ
                  );
      Float_t daliTheta = vdt.Theta();
      //printf("Dali theta from positon %lf,  from header %lf\n", daliTheta/TMath::Pi()*180.0, detInfo->daliHead.fDaliTheta[d]);


      //Float_t df = gammaProj * (1.0 - (betaProj * TMath::Cos((daliTheta/180.0)*TMath::Pi())));
      Float_t df = gammaProj * (1.0 - (betaProj * TMath::Cos(daliTheta)));
      daliEnergyDoppler[d] = daliEnergy[d] * df;

      daliEnergyDopplerSum+=daliEnergyDoppler[d];

    }
  }
  daliCrystalMult = crystalMultDali2;

  // do the addback
  if(daliCrystalMult>=1)  {
    
    //daliAddbackMult=1; // there is at least one addback gamma

    for(Int_t i = 0; i<daliCrystalMult; i++)  {  
      if(daliCrystalUsedForAddback[dali2_crystalFired[i]]==true) {continue;}
      
      Float_t dummyEnergy = daliEnergyDoppler[dali2_crystalFired[i]];
      //float dummyEnergyWithThreshold[11];
      //for(int ppp=0;ppp<11;ppp++)  {
      //  dummyEnergyWithThreshold[ppp]= dali2_dopplerEnergy[i];
      //}

      daliCrystalUsedForAddback[dali2_crystalFired[i]]=true; 
      daliAddbackCrystalMult[daliAddbackMult]++;
      
       
      for(Int_t j = i; j<crystalMultDali2; j++){
        //if(daliCrystalUsedForAddback[dali2_crystalFired[j]]==false && daliEnergyDoppler[dali2_crystalFired[i]]>0.0){
        if(daliCrystalUsedForAddback[dali2_crystalFired[j]]==false){
          for(Int_t k = 0; k<NUMBEROFDALI2ADDBACKCRYSTALS; k++) {
            if(dali2_crystalFired[j] == daliAddbackTable[dali2_crystalFired[i]][k]){
      
              daliCrystalUsedForAddback[dali2_crystalFired[j]]=true;
              dummyEnergy = dummyEnergy + daliEnergyDoppler[dali2_crystalFired[j]];
              
              daliAddbackCrystalMult[daliAddbackMult]++;
              
              //for(int ppp=0;ppp<11;ppp++)  {
              //  if(dali2_dopplerEnergy[i]>=ppp*50&& dali2_dopplerEnergy[j]>=ppp*50) 
              //    dummyEnergyWithThreshold[ppp] = dummyEnergyWithThreshold[ppp] +  dali2_dopplerEnergy[j];
              //}
            }
          }
        }
      }
      //h_dali2_doppler_addback[0]->Fill(dummyEnergy);
      
      
      daliAddbackEnergyDoppler[daliAddbackMult]=dummyEnergy;
      daliAddbackMult++;
      
      ////for(int ppp=0;ppp<11;ppp++)  {
      ////  h_doppler_addback_with_threshold[0][ppp]->Fill(dummyEnergyWithThreshold[ppp]);
      ////}
      //if(posThatsItDali2[2]<0) {
      //  h_dali2_doppler_addback[1]->Fill(dummyEnergy);
      //}
      //else {
      //  h_dali2_doppler_addback[2]->Fill(dummyEnergy);
      //}   
    }
  } // end of addback

}


Bool_t Analysis::IncludeAddbackTable(Float_t det[3][2])  {

  Float_t distance = TMath::Sqrt(TMath::Power(det[0][0]-det[0][1],2) +
                                 TMath::Power(det[1][0]-det[1][1],2) + 
                                 TMath::Power(det[2][0]-det[2][1],2));

  //cout<<"Distance: "<<distance<<endl;
  if( (distance > maxAddbackDistance) || (distance < 0.001) ) return false; // 'maxAddbackDistance' defined in /home/philipp/sim/troja/include/DaliGlobals.hh
  else return true;
} //end of IncludeAddbackTable











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




