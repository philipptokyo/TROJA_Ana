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
#include "Analysis.hh"

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
  cout << endl;

  DetectorInfo* detInfo = new DetectorInfo();
  detInfo->Parse(info->fInFileNameGeometry);
  cout << endl;
  
  TRandom3* randomizer = new TRandom3();
  randomizer->SetSeed(0); 
  


  TStopwatch* watch = new TStopwatch();
  watch->Start();

  Analysis* ana = new Analysis(info, detInfo);
  cout << endl;
  

 
  ana->Analysis1();
  
  Bool_t haveCuts = ana->GetCuts();
  
  if(haveCuts){
    cout << "Found dE-E cuts! Proceeding with analysis2 (missing mass)" << endl;
    ana->Analysis2();
  }else{
    cout << "\033[1;34mGraphical cuts on dE-E plot not yet set! Please create the cuts from the dE-E output of Analysis1!\033[0m" << endl;
    //cout << "\033[1;31mbold red text\033[0m\n" << endl;
    cout << endl;
  }
  
  watch->Stop();
  cout << "Took: real time " << watch->RealTime() << "sec., CPU time " << watch->CpuTime() << " sec." << endl;
  cout << endl;
  
//  Bool_t plot = false;
//  //  plot = true;
//
//  // plot histograms
//  if(plot){  
//
//    TCanvas* candEE = new TCanvas();
//    candEE->cd();
//    hdEE->Draw("colz");
//  
//  
//    TCanvas* canTheta = new TCanvas();
//    canTheta->Divide(1,2);
//    canTheta->cd(1);
//    hThetaLab->Draw();
//    canTheta->cd(2);
//    hThetaCM->Draw();
//  
//    
//  
//    TCanvas* canMissMass = new TCanvas();
//    canMissMass->cd();
//    hMiss->GetXaxis()->SetTitle("E_{miss} / MeV");
//    hMiss->GetYaxis()->SetTitle("cts");
//    gStyle->SetOptStat(0);
//    hMiss->Draw();
//    //hMissTheta->Draw();
//  
//    //FitTwoGaus(hMiss);  
//     
//    TF1* fit1 = new TF1("gaus1", "gaus(0)", -5, 1);
//    TF1* fit2 = new TF1("gaus2", "gaus(0)", -5, 1);
//  
//    fit1->SetParameters(100, -2, 0.1);
//    fit2->SetParameters(100, 0, 0.1);
//  
//    fit1->SetParLimits(0,0,1000000);
//    fit1->SetParLimits(1,-4.5,-1.5);
//    fit1->SetParLimits(2,0.0,1.0  );
//    fit2->SetParLimits(0,0,1000000);
//    fit2->SetParLimits(1,-0.5,0.5);
//    fit2->SetParLimits(2,0.0,1.0  );
//  
//    hMiss->Fit(fit1 , "","", -2.3, -1.8);
//    hMiss->Fit(fit2 , "","", -0.4, 0.4);
//  
//    fit1->Draw("same"); 
//    
//  
//    TCanvas* canEth = new TCanvas();
//    canEth->cd();
//    hEth->Draw("colz");
//
//  }

  //printf("Integrals: Fit1 = %f, Fit2 = %f\n", fit1->Integral(), fit2->Integral()); 
  

  
  
  
  
  
  // if histograms shall be plotted, run theApp
  // otherwise program closes
//  if(plot){
//    theApp->Run();
//  }else{
//    fileAnalysis->Close();
    delete theApp;  
//  }
  return 0;
}


