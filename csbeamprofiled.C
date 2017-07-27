{
  
  char tmp1[200], tmp2[200];
  //Int_t thBins=100;
  //Float_t thMin=1.5, thMax=2.8;
  Int_t thBins=180;
  //Float_t thMin=0.0, thMax=3.1;
  Float_t thMin=0.0, thMax=180.0;
  // missing mass cut
  //Float_t mmcMin=-0.2, mmcMax=0.2; // ground state (1)
  Float_t mmcMin=-0.1, mmcMax=0.2; // ground state (1)
  //Float_t mmcMin=-1.0, mmcMax=-0.7; // state 2
  //Float_t mmcMin=-1.45, mmcMax=-1.0; // state 3
  //Float_t mmcMin=-1.8, mmcMax=-1.45; // state 4
  //Float_t mmcMin=-2.2, mmcMax=-1.8; // state 5
  Float_t scaling[200]={1.0}; 
  
  Int_t state=4;
  Float_t massNr=56.0;
  Int_t rebinnr=4;
  
  //Int_t maxEvents=500;
  //Float_t scalefactor=0.15;
  //Int_t maxEvents=1000;
  //Float_t scalefactor=0.2;
  //Int_t maxEvents=1200;
  //Float_t scalefactor=0.07;
  Int_t maxEvents=1500;
  Float_t scalefactor=0.047;
  Float_t specFact=0.36;

  TFile* fileFre;
  TFile* filePhy;
  TFile* fileAcc;
  TFile* fileNor;
  
  TH2F* h2FreCM; // fresco histogram
  
  TH2F* h2GenLab;// generated events (from fresco)
  TH2F* h2GenCM;// generated events (from fresco)
  
  TH2F* h2PhyCM;
  TH2F* h2PhyLab;
  
  TH2F* h2NorCM;
  TH2F* h2NorLab;
  
  TH2F* h2AccCM;
  TH2F* h2AccLab;
  
  TH2F* h2AccLabNor;
  TH2F* h2AccCMNor;
  
  TH2F* h2PhyCMCorr; 
  TH2F* h2PhyLabCorr;
  
  //TH2F* h2LabE_CM_Mat; 
  //TH3F* h3LabE_CM_Mat; 
  
//  TH2F* h2Phy1 = new TH2F("h2Phy1", "Simulation", 22, 600.0, 2700.0, thBins, thMin, thMax);
//  TH2F* h2Phy2 = new TH2F("h2Phy2", "Simulation, with beam profile", 21, 600.0, 2600.0, thBins, thMin, thMax);
//  TH2F* h2PhyCorr1 = new TH2F("h2PhyCorr1", "Simulation, Acceptance corrected", 22, 600.0, 2700.0, thBins, thMin, thMax);
//  TH2F* h2PhyCorr2 = new TH2F("h2PhyCorr2", "Simulation, Acceptance corrected, with beam profile", 21, 600.0, 2600.0, thBins, thMin, thMax);
//  TH2F* h2Fre1 = new TH2F("h2Fre1", "Fresco", 22, 600.0, 2700.0, thBins, thMin, thMax);
//  TH2F* h2Fre2 = new TH2F("h2Fre2", "Fresco, with beam profile", 21, 600.0, 2600.0, thBins, thMin, thMax);
  
  
  // include the beam profile
//  TFile* fileBeam = TFile::Open("/mnt/raid/OEDO/OEDO_Matsushita/132Sn/3rd_order/oedo_132sn_10MeV_short.root", "read");
//  fileBeam->cd();
//  events->Draw("energy>>beamProfiletmp(100,0,35)");
//  TH1F* hBeamProfile = (TH1F*)beamProfiletmp->Clone();

  
  //scaling[f]=hBeamProfile->GetBinContent(b);
//  scaling[f]=hBeamProfile->GetBinContent(b)/hBeamProfile->Integral();
//  printf("Run %d, Beam energy %f, intensity %f\n", 600 + 100*f, e, scaling[f]);

  //sprintf(tmp2, "132Sn_jones");
  sprintf(tmp2, "56Cr_exp_dL");

  sprintf(tmp1, "/mnt/raid/OEDO/philipp/events/%s_phy.root", tmp2);
  //printf("Opening file '%s'\n", tmp);
  fileFre = TFile::Open(tmp1 , "read");
  
  sprintf(tmp1, "/mnt/raid/OEDO/philipp/analysis/%s_phy.root", tmp2);
  //printf("Opening file '%s'\n", tmp);
  filePhy = TFile::Open(tmp1, "read");
  
  sprintf(tmp1, "/mnt/raid/OEDO/philipp/events/%s_acc.root", tmp2);
  //printf("Opening file '%s'\n", tmp);
  fileNor = TFile::Open(tmp1 , "read");
  
  sprintf(tmp1, "/mnt/raid/OEDO/philipp/analysis/%s_acc.root", tmp2);
  //printf("Opening file '%s'\n", tmp);
  fileAcc = TFile::Open(tmp1 , "read");
  
  


  // *****************************************************************************
  // ************************ Lab System *****************************************
  // *****************************************************************************
  
  // fresco events
  fileFre->cd();
  h2FreCM = (TH2F*)fileFre->Get("histCScmFresco2d_01"); // state 1, ground state

  Int_t eBins = h2FreCM->GetYaxis()->GetNbins();
  Float_t eMin, eMax;
  eMin = h2FreCM->GetYaxis()->GetBinLowEdge(0);
  eMax = h2FreCM->GetYaxis()->GetBinUpEdge(eBins);
  printf("Found %d beam energy bins from %f to %f MeV\n", eBins, eMin, eMax);

  thBins = h2FreCM->GetXaxis()->GetNbins();
  thMin = h2FreCM->GetXaxis()->GetBinLowEdge(0);
  thMax = h2FreCM->GetXaxis()->GetBinUpEdge(thBins);
  
  thBins = 90;
  thMin = 0;
  thMax = 180.0;
  printf("Found %d theta CM bins from %f to %f\n", thBins, thMin, thMax);
  
  
  sprintf(tmp1, "beamEnergy*%f:lightTheta/TMath::Pi()*180.0>>hfretmp(%f, %f, %f, %f, %f, %f)",massNr, thBins, thMin, thMax, eBins, eMin, eMax);
  //sprintf(tmp2, "missingMass>%f && missingMass<%f", mmcMin, mmcMax);
  sprintf(tmp2, "state==%d", state);
  events->Draw(tmp1, tmp2,"");
  h2GenLab = (TH2F*)hfretmp->Clone();
  h2GenLab->SetTitle( Form("Generated events, #vartheta_{lab}, E_{beam} vs. #vartheta_{lab}") );

  //sprintf(tmp1, "lightThetaCM:lightTheta>>hmat(%f, %f, %f, %f, %f, %f)", thBins, thMin, thMax, thBins, thMin, thMax);
  //sprintf(tmp2, "state==1");
  //events->Draw(tmp1, tmp2,"");
  //h2LabE_CM_Mat = (TH2F*)hmat->Clone();

  // physics/fresco simulation
  filePhy->cd();
  sprintf(tmp1, "genBeamEnergy*%f:simLightThetaLab/TMath::Pi()*180.0>>hphytmp(%f, %f, %f, %f, %f, %f)", massNr, thBins, thMin, thMax, eBins, eMin, eMax);
  //sprintf(tmp2, "anaMissingMass>%f && anaMissingMass<%f", mmcMin, mmcMax);
  sprintf(tmp2, "genState==%d", state);
  analysis2_1->Draw(tmp1, tmp2,"",maxEvents);
  h2PhyLab = (TH2F*)hphytmp->Clone();
  h2PhyLabCorr = (TH2F*)hphytmp->Clone();
  h2PhyLab->SetTitle( Form("Physics Simulation, E_{beam} vs. #vartheta_{lab}") );
  h2PhyLabCorr->SetTitle( Form("Physics Simulation, Acceptance corrected, E_{beam} vs. #vartheta_{lab}") );
  

  ////sprintf(tmp1, "simLightThetaCM:genBeamEnergy*132.0:simLightThetaLab>>hmat(%f, %f, %f, %f, %f, %f, %f, %f, %f)", thBins, thMin, thMax, eBins, eMin, eMax, thBins, thMin, thMax);
  ////sprintf(tmp2, "anaMissingMass>%f && anaMissingMass<%f", mmcMin, mmcMax);
  ////analysis2_1->Draw(tmp1, tmp2,"");
  ////h3LabE_CM_Mat = (TH3F*)hmat->Clone();
  ////h2LabE_CM_Mat = (TH2F*)h3LabE_CM_Mat->Project3D("yx");
  //sprintf(tmp1, "simLightThetaCM:simLightThetaLab>>hmat(%f, %f, %f, %f, %f, %f)", thBins, thMin, thMax, thBins, thMin, thMax);
  //sprintf(tmp2, "anaMissingMass>%f && anaMissingMass<%f", mmcMin, mmcMax);
  //analysis2_1->Draw(tmp1, tmp2,"");
  //h2LabE_CM_Mat = (TH2F*)hmat->Clone();

  
  // acceptence events
  fileNor->cd();
  sprintf(tmp1, "beamEnergy*%f:lightTheta/TMath::Pi()*180.0>>hnortmp(%f, %f, %f, %f, %f, %f)", massNr, thBins, thMin, thMax, eBins, eMin, eMax);
  //sprintf(tmp2, "missingMass>%f && missingMass<%f", mmcMin, mmcMax);
  sprintf(tmp2, "state==%d", state);
  events->Draw(tmp1, tmp2,"");
  h2NorLab = (TH2F*)hnortmp->Clone();
  h2NorLab->SetTitle( Form("Normalization, E_{beam} vs. #vartheta_{lab}") );
  

  // acceptacne simulation
  fileAcc->cd();
  sprintf(tmp1, "genBeamEnergy*%f:simLightThetaLab/TMath::Pi()*180.0>>hacctmp(%f, %f, %f, %f, %f, %f)", massNr, thBins, thMin, thMax, eBins, eMin, eMax);
  //sprintf(tmp2, "anaMissingMass>%f && anaMissingMass<%f", mmcMin, mmcMax);
  sprintf(tmp2, "genState==%d", state);
  analysis2_1->Draw(tmp1, tmp2,"");
  h2AccLab = (TH2F*)hacctmp->Clone();
  h2AccLabNor = (TH2F*)hacctmp->Clone();
  h2AccLab->SetTitle( Form("Acceptance (NOT normalized), E_{beam} vs. #vartheta_{lab}") );
  h2AccLabNor->SetTitle( Form("Acceptance (normalized), E_{beam} vs. #vartheta_{lab}") );

  h2AccLabNor->Divide(h2NorLab);
  
  
  // acceptance correction of physics data 
  h2PhyLabCorr->Divide(h2AccLabNor);

  //h2PhyCMCorr = (TH2F*)h2PhyLabCorr->Clone(0); // to get the same binning

  //// convert to CM with the obtained 3d matrix
  //for(Int_t t=0; t<thBins; t++){
  //for(Int_t e=0; e<eBins; e++){
  //  Float_t cs = h2PhyLabCorr->GetBinContent(t,e);
  //  Float_t thcm = h2LabE_CM_Mat->GetBinContent(t,e);   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  //  Int_t bcm = h2PhyCMCorr->GetXaxis()->FindBin(thcm);
  //  //h2PhyCMCorr->SetBinContent(bcm, e, cs);
  //  h2PhyCMCorr->Fill(bcm, e, cs); // increment by cs
  //}
  //}




  // *****************************************************************************
  // ************************* CM System *****************************************
  // *****************************************************************************
  
  // fresco events
  fileFre->cd();
  
  sprintf(tmp1, "beamEnergy*%f:lightThetaCM/TMath::Pi()*180.0>>hfretmp(%f, %f, %f, %f, %f, %f)", massNr, thBins, thMin, thMax, eBins, eMin, eMax);
  //sprintf(tmp2, "missingMass>%f && missingMass<%f", mmcMin, mmcMax);
  sprintf(tmp2, "state==%d", state);
  events->Draw(tmp1, tmp2,"");
  h2GenCM = (TH2F*)hfretmp->Clone();
  h2GenCM->SetTitle( Form("Generated events, #vartheta_{CM}, E_{beam} vs. #vartheta_{CM}") );

  // physics/fresco simulation
  filePhy->cd();
  sprintf(tmp1, "genBeamEnergy*%f:simLightThetaCM/TMath::Pi()*180.0>>hphytmp(%f, %f, %f, %f, %f, %f)", massNr, thBins, thMin, thMax, eBins, eMin, eMax);
  //sprintf(tmp2, "anaMissingMass>%f && anaMissingMass<%f", mmcMin, mmcMax);
  sprintf(tmp2, "genState==%d", state);
  analysis2_1->Draw(tmp1, tmp2,"");
  h2PhyCM = (TH2F*)hphytmp->Clone();
  h2PhyCMCorr = (TH2F*)hphytmp->Clone();
  h2PhyCM->SetTitle( Form("Physics Simulation, E_{beam} vs. #vartheta_{CM}") );
  h2PhyCMCorr->SetTitle( Form("Phys Sim, Acc corr, E_{beam} vs. #vartheta_{CM}") );
  

  // acceptence events
  fileNor->cd();
  sprintf(tmp1, "beamEnergy*%f:lightThetaCM/TMath::Pi()*180.0>>hnortmp(%f, %f, %f, %f, %f, %f)", massNr, thBins, thMin, thMax, eBins, eMin, eMax);
  //sprintf(tmp2, "missingMass>%f && missingMass<%f", mmcMin, mmcMax);
  sprintf(tmp2, "state==%d", state);
  events->Draw(tmp1, tmp2,"");
  h2NorCM = (TH2F*)hnortmp->Clone();
  h2NorCM->SetTitle( Form("Normalization, E_{beam} vs. #vartheta_{CM}") );
  

  // acceptacne simulation
  fileAcc->cd();
  sprintf(tmp1, "genBeamEnergy*%f:simLightThetaCM/TMath::Pi()*180.0>>hacctmp(%f, %f, %f, %f, %f, %f)", massNr, thBins, thMin, thMax, eBins, eMin, eMax);
  //sprintf(tmp2, "anaMissingMass>%f && anaMissingMass<%f", mmcMin, mmcMax);
  sprintf(tmp2, "genState==%d", state);
  analysis2_1->Draw(tmp1, tmp2,"");
  h2AccCM = (TH2F*)hacctmp->Clone();
  h2AccCMNor = (TH2F*)hacctmp->Clone();
  h2AccCM->SetTitle( Form("Acceptance (NOT normalized), E_{beam} vs. #vartheta_{CM}") );
  h2AccCMNor->SetTitle( Form("Acceptance (normalized), E_{beam} vs. #vartheta_{CM}") );

  h2AccCMNor->Divide(h2NorCM);
  
  
  // acceptance correction of physics data 
  h2PhyCMCorr->Divide(h2AccCMNor);























  // *****************************************************************************
  // ****************************** Plot *****************************************
  // *****************************************************************************

  sprintf(tmp1,"colz");

  TCanvas* can3 = new TCanvas("can3", "Analysis, CM");
  can3->Divide(3,1);
  can3->cd(1);
  h2FreCM->Draw(tmp1);
  can3->cd(2);
  h2GenCM->Draw(tmp1);
  //h2LabE_CM_Mat->Draw(tmp1);
  can3->cd(3);
  h2GenLab->Draw(tmp1);
  //h2PhyCMCorr->Draw(tmp1);

  TCanvas* can4 = new TCanvas("can4", "Sim, Lab");
  can4->Divide(3,2);
  can4->cd(1);
  h2GenLab->Draw(tmp1);
  can4->cd(2);
  h2PhyLab->Draw(tmp1);
  can4->cd(3);
  h2PhyLabCorr->Draw(tmp1);
  
  can4->cd(4);
  h2NorLab->Draw(tmp1);
  can4->cd(5);
  h2AccLab->Draw(tmp1);
  can4->cd(6);
  h2AccLabNor->Draw(tmp1);

  TCanvas* can5 = new TCanvas("can5", "Sim, CM");
  can5->Divide(3,2);
  can5->cd(1);
  h2GenCM->Draw(tmp1);
  can5->cd(2);
  h2PhyCM->Draw(tmp1);
  can5->cd(3);
  h2PhyCMCorr->Draw(tmp1);
  
  can5->cd(4);
  h2NorCM->Draw(tmp1);
  can5->cd(5);
  h2AccCM->Draw(tmp1);
  can5->cd(6);
  h2AccCMNor->Draw(tmp1);


//  TCanvas* can6 = new TCanvas("can6", "CS, CM");
//  can6->cd();
//  TH1D* h1cscm_pr = h2PhyCMCorr->ProjectionX("h1cscm_pr", 0, -1);
//  h1cscm_pr->Sumw2();
//  h1cscm_pr->Draw();

  TCanvas* can7 = new TCanvas("can7", "CS, Lab");
  can7->cd();
  TH1D* h1cslab_pr = h2PhyLabCorr->ProjectionX("h1cslab_pr", 0, -1);
  for(Int_t b=0; b<h1cslab_pr->GetXaxis()->GetNbins(); b++){
    Double_t cont=h1cslab_pr->GetBinContent(b);
    Double_t angl=h1cslab_pr->GetXaxis()->GetBinCenter(b)/180.0*TMath::Pi();
    //cont/=TMath::Sin(angl);

    h1cslab_pr->SetBinContent(b, cont);
  }
  h1cslab_pr->Sumw2();
  h1cslab_pr->Draw();
  h1cslab_pr->Rebin(rebinnr);
  h1cslab_pr->GetXaxis()->SetRangeUser(80,180);
  h1cslab_pr->GetXaxis()->SetTitle("#vartheta_lab");
  h1cslab_pr->GetYaxis()->SetTitle("d#sigma/d#vartheta_lab in mb/deg");

  h1cslab_pr->Scale(scalefactor*specFact);













  Double_t p[8];
p[0]                        =      180.761;
p[1]                        =     -3.25234;
p[2]                        =    0.0405916;
p[3]                        = -0.000238102;
p[4]                        =  -1.0551e-06;
p[5]                        =  2.29259e-08;
p[6]                        = -1.18655e-10;
p[7]                        =  2.10794e-13;

p[0]                        =      3.16169;
p[1]                        =     -3.33532;
p[2]                        =      2.43578;
p[3]                        =    -0.740411;
p[4]                        =     -0.36691;
p[5]                        =     0.368405;
p[6]                        =    -0.109112;
p[7]                        =    0.0113104;

p[0]                        =      3.16639; 
p[1]                        =      -3.3419; 
p[2]                        =      2.52711; 
p[3]                        =    -0.963284; 
p[4]                        =    -0.139066; 
p[5]                        =     0.252848; 
p[6]                        =    -0.080406; 
p[7]                        =   0.00852544; 

  //TF1* cm2lab = new TF1*("cm2lab", "3.16571 - x*3.34246 + x*x*2.52951 - x*x*x*0.955381 - x*x*x*x*0.158001 + x*x*x*x*x*0.267008 - x*x*x*x*x*x*0.0848927 + x*x*x*x*x*x*x*0.00904326", 0, TMath::Pi());
  //TF1* cm2lab = new TF1*("cm2lab", "par[0] - x*par[1] + x*x*par[2] - x*x*x*par[3] - x*x*x*x*par[4] + x*x*x*x*x*par[5] - x*x*x*x*x*x*par[6] + x*x*x*x*x*x*x*par[7]", 0, TMath::Pi());
  //TF1* cm2lab = new TF1*("cm2lab", "pol7", 0, TMath::Pi());

  TFile* file[3];
  //file[0]=TFile::Open("/mnt/raid/OEDO/philipp/events/56Cr_10.0.root", "read");
  //file[1]=TFile::Open("/mnt/raid/OEDO/philipp/events/56Cr_09.2.root", "read");
  //file[2]=TFile::Open("/mnt/raid/OEDO/philipp/events/56Cr_11.0.root", "read");
  file[0]=TFile::Open("/mnt/raid/OEDO/philipp/events/56Cr_exp_dL_10AMeV.root", "read");

  TH1F* hist[6];

  // Cr, dL = 1
  //hist[0]=(TH1F*)file[0]->Get("state02/histCScmFresco_00");
  //hist[1]=(TH1F*)file[1]->Get("state02/histCScmFresco_00");
  //hist[2]=(TH1F*)file[2]->Get("state02/histCScmFresco_00");
  //hist[1]->SetLineStyle(2);
  //hist[2]->SetLineStyle(3);
  hist[0]=(TH1F*)file[0]->Get("state03/histCScmFresco_00");

  // Cr, dL = 3
  //hist[3]=(TH1F*)file[0]->Get("state04/histCScmFresco_00");
  //hist[4]=(TH1F*)file[1]->Get("state04/histCScmFresco_00");
  //hist[5]=(TH1F*)file[2]->Get("state04/histCScmFresco_00");
  hist[3]=(TH1F*)file[0]->Get("state04/histCScmFresco_00");

  //TH1F* histlab[6];
  TGraph* graflab[6];
  Int_t nbins = hist[0]->GetXaxis()->GetNbins();
  Double_t binmin=0;
  Double_t binmax=180.0;

  for(Int_t h=0; h<6; h++){
    if(!(h==0 || h==3)){
      continue;
    }
    //histlab[h]=new TH1F(Form("histlab%d", h), Form("histlab%d", h), 180,0,180);
    graflab[h]=new TGraph();

    for(Int_t b=0; b<nbins; b++){
      Double_t tcm = hist[h]->GetXaxis()->GetBinCenter(b);
      Double_t cs = hist[h]->GetBinContent(b);
      Double_t tlab;// = cm2lab->Eval(tcm)/TMath::Pi()*180.0;

      Double_t x = tcm;
      //tlab=3.16571 - x*3.34246 + x*x*2.52951 - x*x*x*0.955381 - x*x*x*x*0.158001 + x*x*x*x*x*0.267008 - x*x*x*x*x*x*0.0848927 + x*x*x*x*x*x*x*0.00904326;
      tlab=p[0] + x*p[1] + x*x*p[2] + x*x*x*p[3] + x*x*x*x*p[4] + x*x*x*x*x*p[5] + x*x*x*x*x*x*p[6] + x*x*x*x*x*x*x*p[7];
      tlab*= 180.0/TMath::Pi();

      //Int_t binlab = histlab[h]->GetXaxis()->FindBin(tlab);
      //histlab[h]->SetBinContent(binlab, cs);
      
      if(h==0){
        graflab[h]->SetPoint(b, tlab, cs*specFact*4.0);
      }
      if(h==3){
        graflab[h]->SetPoint(b, tlab, cs*specFact);
      }

      //graflab[h]->SetPoint(b, tlab, cs*TMath::Sin(tcm*TMath::Pi()/180.0));
      //printf("tcm %f, tlab %f, cs %f\n", tcm, tlab, cs);
    }

  }

  //graflab[1]->SetLineStyle(2);
  //graflab[2]->SetLineStyle(3);
  //graflab[4]->SetLineStyle(2);
  //graflab[5]->SetLineStyle(3);

  graflab[0]->SetLineWidth(2);
  graflab[0]->SetLineColor(kBlack);
  //graflab[1]->SetLineColor(kBlack);
  //graflab[2]->SetLineColor(kBlack);

  graflab[3]->SetLineWidth(2);
  graflab[3]->SetLineColor(kRed);
  //graflab[4]->SetLineColor(kRed);
  //graflab[5]->SetLineColor(kRed);

  graflab[0]->GetXaxis()->SetRangeUser(80, 175);
  graflab[0]->GetXaxis()->SetTitle("#vartheta_{lab}");
  graflab[0]->GetYaxis()->SetRangeUser(0.11, 20);
  graflab[0]->GetYaxis()->SetTitle("d#sigma/d#vartheta in mb/sr");

  //gPad->SetLogy();

  can7->cd();
//  graflab[0]->Draw("Csame");
//  //graflab[1]->Draw("Csame");
//  //graflab[2]->Draw("Csame");
//  graflab[3]->Draw("Csame");
//  //graflab[4]->Draw("Csame");
//  //graflab[5]->Draw("Csame");

































}
