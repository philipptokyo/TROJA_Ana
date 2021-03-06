#include "EnLoss.hh"


EnLoss::EnLoss(){

}

EnLoss::EnLoss(InputInfo *i, DetectorInfo* d){

  fInfo = i;
  fDetInfo = d;

}


EnLoss::~EnLoss(){
}



void EnLoss::CollectData(string filename, Int_t index){
  
  printf("\n\nParsing energy loss data file\n");


  //check if in/out file exists
  std::ifstream fin;
  std::ifstream FileTest(filename.c_str());
  if(!FileTest){
    printf("Can not open file '%s'! Does it exist?\n", filename.c_str());
    abort();
  }else{
    printf("Opening file '%s'\n", filename.c_str());
    fin.open(filename.c_str());
  }


  //parse input file line-wise
  std::string line;
  Int_t counter=-1;
  const Int_t stopper=2*maxEnLossEntriesE;

  char temp[maxEnLossEntriesX+2][100];

  fEnLoss[index].ebins=0;
  fEnLoss[index].xbins=0;

  while(fin.good())  {

    counter++;

    //reset values
    for(Int_t i=0; i< maxEnLossEntriesX+2; i++){temp[i][0]='\0';}

    getline(fin,line);
    //printf("Parsing line '%s'\n",line.c_str());

    //skip comments
    line.copy(temp[0],2,0); //first two characters
    temp[0][2]='\0';
    //printf("temp1=%s\n",temp1);
    if(strcmp(temp[0],"//")==0){continue;}
    
    std::istringstream iss(line);
    for(Int_t i=0; i<maxEnLossEntriesX+2; i++){iss >> temp[i];}

    //skip empty lines
    if(strcmp(temp[0],"")==0){continue;}

    if(strcmp(temp[0],"E_X")==0){ // first line with x/distance values
      printf("Got dX values: ");
      for(Int_t i=0; i<maxEnLossEntriesX+2; i++){
        if(strcmp(temp[i+1],"micron")==0){
          break;
        }else{
          fEnLoss[index].x[i]=atof(temp[i+1]);
          printf("%f (%d), ", fEnLoss[index].x[i], i);
          fEnLoss[index].xbins++;
          if(fEnLoss[index].xbins>=maxEnLossEntriesX){
            printf("Seg fault warning: xbins are at the limit of %d! Increase the 'maxEnLossEntriesX' in EnLoss.hh!\n", maxEnLossEntriesX);
            abort();
          }
        }
      }
      printf("\nFound %d entries for the distance (x) through the material\n", fEnLoss[index].xbins);

    }else if(strcmp(temp[0],"MeV")!=0){
      fEnLoss[index].e[fEnLoss[index].ebins]=atof(temp[0]);
      printf("Got particle energy %f, energy losses: ", fEnLoss[index].e[fEnLoss[index].ebins]);

      for(Int_t i=0; i<fEnLoss[index].xbins; i++){
        fEnLoss[index].enLoss[fEnLoss[index].ebins][i]=atof(temp[i+1]);
        printf("%f ", fEnLoss[index].enLoss[fEnLoss[index].ebins][i]);
      }
      printf("(entry %d)\n", fEnLoss[index].ebins);
      fEnLoss[index].ebins++;
      if(fEnLoss[index].ebins>=maxEnLossEntriesE){
        printf("Seg fault warning: ebins are at the limit of %d! Increase the 'maxEnLossEntriesE' in EnLoss.hh!\n", maxEnLossEntriesX);
        abort();
      }
    }else{break;}

    //count lines
    if(counter>stopper){
      printf("Reached %d lines in file '%s'! Stopping!\n", counter, filename.c_str());
      abort();
    }



  } // end of while loop

  
  printf("EnLoss::CollectData finished! Inverting matrix...\n");
  InvertData(index, 0.001);
  printf("... finished inverting matrix!\n");

}








void EnLoss::InvertData(Int_t index, Double_t prec=0.001){
  
  // copy energy and thickness binnings:
  fEnLossInv[index].ebins=fEnLoss[index].ebins;
  fEnLossInv[index].xbins=fEnLoss[index].xbins;

  for(Int_t eb=0; eb<fEnLoss[index].ebins; eb++){
    fEnLossInv[index].e[eb] = fEnLoss[index].e[eb];
    
    printf("Det energy %f MeV (thickness[micron] / particle energy[MeV]): \n", fEnLossInv[index].e[eb]);
    
    for(Int_t xb=0; xb<fEnLoss[index].xbins; xb++){
      fEnLossInv[index].x[xb] = fEnLoss[index].x[xb];

      Double_t eninv=fEnLossInv[index].e[eb];

      Int_t counter=0;
      
      //printf("  Start iterating: eninv %f, in thickness %f, EnLoss %f\n", eninv, fEnLossInv[index].x[xb], CalcEnLoss(eninv, fEnLossInv[index].x[xb], index, false));
      while( (eninv-CalcEnLoss(eninv, fEnLossInv[index].x[xb], index, false)) < fEnLossInv[index].e[eb] ){
        
        //eninv += CalcEnLoss(eninv, fEnLossInv[index].x[xb], index);
        eninv += prec; 
      
        //printf("    iterating: eninv %f, EnLoss %f\n", eninv, CalcEnLoss(eninv, fEnLossInv[index].x[xb], index));
        
        counter++;
        if(counter>(Int_t)(10000.0/prec)){
          printf("Error while creating InvertData (EnLoss.cc): more than %d iterations for energy %f and thickness %f! Aborting...\n", (Int_t)(10000.0/prec), fEnLossInv[index].e[eb], fEnLossInv[index].x[xb]);
          abort();
        }

      }


      fEnLossInv[index].enLoss[eb][xb] = eninv;

      printf("%3.1f / %2.3f,  ", fEnLossInv[index].x[xb], fEnLossInv[index].enLoss[eb][xb]);

    }

    printf("\n");

  }

    

}





Double_t EnLoss::CalcEnLoss(Double_t en, Double_t dist, Int_t index){
  return CalcEnLoss(en, dist, index, true);
}


Double_t EnLoss::CalcEnLoss(Double_t en, Double_t dist, Int_t index, Bool_t warnings=true){
  
  en-=0.000001;
  dist-=0.000001;

  //printf("Calc En Loss for en %f, dist %f\n", en, dist);
  //printf("x[xbins-1] %f, e[ebins-1] %f\n", fEnLoss[index].x[fEnLoss[index].xbins-1], fEnLoss[index].e[fEnLoss[index].ebins-1]);

  //determine energy bin
  Int_t ebin=0;
  if(en>fEnLoss[index].e[fEnLoss[index].ebins-1]){
    if(warnings){
      printf("paricle energy %f MeV above maximum given value in list (%f) - extrapolating from bins %d (at %f mu) and %d (at %f mu) \n", en, fEnLoss[index].e[fEnLoss[index].ebins-1], fEnLoss[index].ebins-2, fEnLoss[index].e[fEnLoss[index].ebins-2], fEnLoss[index].ebins-1, fEnLoss[index].e[fEnLoss[index].ebins-1]);
    }
    ebin=fEnLoss[index].ebins-2;
  }else{
    for(Int_t eb=0; eb<fEnLoss[index].ebins-1; eb++){
      if((en>fEnLoss[index].e[eb]) && (en<fEnLoss[index].e[eb+1] )){
        ebin = eb;
        break;
      }
    }
  }
  
  //determine x bin
  Int_t xbin=0;
  if(dist>fEnLoss[index].x[fEnLoss[index].xbins-1]){
    if(warnings){
      printf("dX %f micron above maximum value given in list (%f) - extrapolating from bins %d (at %f mu) and %d (at %f mu) \n", dist, fEnLoss[index].x[fEnLoss[index].xbins-1], fEnLoss[index].xbins-2, fEnLoss[index].x[fEnLoss[index].xbins-2], fEnLoss[index].xbins-1, fEnLoss[index].x[fEnLoss[index].xbins-1]);
    }
    xbin=fEnLoss[index].xbins-2;
  }else{
    for(Int_t xb=0; xb<fEnLoss[index].xbins-1; xb++){
      if((dist>fEnLoss[index].x[xb]) && (dist<fEnLoss[index].x[xb+1] )){
        xbin = xb;
        break;
      }
    }
  }
  
  //printf("using xbin %d, ebin %d\n", xbin, ebin);

  //interpolating at xbin and xbin+1
  TGraph *gx1 = new TGraph();
  gx1->SetPoint(0, fEnLoss[index].e[ebin], fEnLoss[index].enLoss[ebin][xbin]);
  gx1->SetPoint(1, fEnLoss[index].e[ebin+1], fEnLoss[index].enLoss[ebin+1][xbin]);
  
  //printf("Graph 1: e %f %f, dE %f %f\n", fEnLoss[index].e[ebin], fEnLoss[index].e[ebin+1], fEnLoss[index].enLoss[ebin][xbin], fEnLoss[index].enLoss[ebin+1][xbin]);

  TGraph *gx2 = new TGraph();
  gx2->SetPoint(0, fEnLoss[index].e[ebin], fEnLoss[index].enLoss[ebin][xbin+1]);
  gx2->SetPoint(1, fEnLoss[index].e[ebin+1], fEnLoss[index].enLoss[ebin+1][xbin+1]);

  //printf("Graph 2: e %f %f, dE %f %f\n", fEnLoss[index].e[ebin], fEnLoss[index].e[ebin+1], fEnLoss[index].enLoss[ebin][xbin+1], fEnLoss[index].enLoss[ebin+1][xbin+1]);

  Double_t dE1=gx1->Eval(en);
  Double_t dE2=gx2->Eval(en);

  TGraph *gem = new TGraph();
  gem->SetPoint(0, fEnLoss[index].x[xbin], dE1);
  gem->SetPoint(1, fEnLoss[index].x[xbin+1], dE2);
  
  //printf("Graph M: x %f %f, dE %f %f\n", fEnLoss[index].x[xbin], fEnLoss[index].x[xbin+1], dE1, dE2);

  Double_t enloca = gem->Eval(dist);
 
  //printf("Resulting En Loss is %f\n", enloca);

  return enloca;
}






Double_t EnLoss::CalcParticleEnergy(Double_t en, Double_t dist, Int_t index){
  
  en-=0.000001;
  dist-=0.000001;

  //printf("Calc En Loss for en %f, dist %f\n", en, dist);
  //printf("x[xbins-1] %f, e[ebins-1] %f\n", fEnLoss[index].x[fEnLoss[index].xbins-1], fEnLoss[index].e[fEnLoss[index].ebins-1]);

  //determine energy bin
  Int_t ebin=0;
  if(en>fEnLossInv[index].e[fEnLossInv[index].ebins-1]){
    printf("paricle energy %f MeV above maximum given value in list (%f) - extrapolating from bins %d (at %f mu) and %d (at %f mu) \n", en, fEnLossInv[index].e[fEnLossInv[index].ebins-1], fEnLossInv[index].ebins-2, fEnLossInv[index].e[fEnLossInv[index].ebins-2], fEnLossInv[index].ebins-1, fEnLossInv[index].e[fEnLossInv[index].ebins-1]);
    ebin=fEnLossInv[index].ebins-2;
  }else{
    for(Int_t eb=0; eb<fEnLossInv[index].ebins-1; eb++){
      if((en>fEnLossInv[index].e[eb]) && (en<fEnLossInv[index].e[eb+1] )){
        ebin = eb;
        break;
      }
    }
  }
  
  //determine x bin
  Int_t xbin=0;
  if(dist>fEnLossInv[index].x[fEnLossInv[index].xbins-1]){
    printf("dX %f micron above maximum value given in list (%f) - extrapolating from bins %d (at %f mu) and %d (at %f mu) \n", dist, fEnLossInv[index].x[fEnLossInv[index].xbins-1], fEnLossInv[index].xbins-2, fEnLossInv[index].x[fEnLossInv[index].xbins-2], fEnLossInv[index].xbins-1, fEnLossInv[index].x[fEnLossInv[index].xbins-1]);
    xbin=fEnLossInv[index].xbins-2;
  }else{
    for(Int_t xb=0; xb<fEnLossInv[index].xbins-1; xb++){
      if((dist>fEnLossInv[index].x[xb]) && (dist<fEnLossInv[index].x[xb+1] )){
        xbin = xb;
        break;
      }
    }
  }
  
  //printf("using xbin %d, ebin %d\n", xbin, ebin);

  //interpolating at xbin and xbin+1
  TGraph *gx1 = new TGraph();
  gx1->SetPoint(0, fEnLossInv[index].e[ebin], fEnLossInv[index].enLoss[ebin][xbin]);
  gx1->SetPoint(1, fEnLossInv[index].e[ebin+1], fEnLossInv[index].enLoss[ebin+1][xbin]);
  
  //printf("Graph 1: e %f %f, dE %f %f\n", fEnLossInv[index].e[ebin], fEnLossInv[index].e[ebin+1], fEnLossInv[index].enLoss[ebin][xbin], fEnLossInv[index].enLoss[ebin+1][xbin]);

  TGraph *gx2 = new TGraph();
  gx2->SetPoint(0, fEnLossInv[index].e[ebin], fEnLossInv[index].enLoss[ebin][xbin+1]);
  gx2->SetPoint(1, fEnLossInv[index].e[ebin+1], fEnLossInv[index].enLoss[ebin+1][xbin+1]);

  //printf("Graph 2: e %f %f, dE %f %f\n", fEnLossInv[index].e[ebin], fEnLossInv[index].e[ebin+1], fEnLossInv[index].enLoss[ebin][xbin+1], fEnLossInv[index].enLoss[ebin+1][xbin+1]);

  Double_t dE1=gx1->Eval(en);
  Double_t dE2=gx2->Eval(en);

  TGraph *gem = new TGraph();
  gem->SetPoint(0, fEnLossInv[index].x[xbin], dE1);
  gem->SetPoint(1, fEnLossInv[index].x[xbin+1], dE2);
  
  //printf("Graph M: x %f %f, dE %f %f\n", fEnLossInv[index].x[xbin], fEnLossInv[index].x[xbin+1], dE1, dE2);

  Double_t enca = gem->Eval(dist);
 
  //printf("Resulting En Loss is %f\n", enloca);

  return enca;
}

