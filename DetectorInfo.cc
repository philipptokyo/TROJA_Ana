#include "DetectorInfo.hh"
#include "DetectorGlobals.hh"

//#include "G4PhysicalConstants.hh"
//#include "G4SystemOfUnits.hh"
//
//#include "G4Material.hh"
//#include "TVector3.hh"

//#include "TMath.h"
//
//#include <iostream>
//#include <sstream>
//#include <fstream>
//#include <string>

#include "LibPerso.h"

const Double_t mm = 1.0;
const Double_t deg = TMath::Pi()/180.0;

DetectorInfo::DetectorInfo()
{ 
  
  ClearGeometry();
  ResetData();
  
}

DetectorInfo::~DetectorInfo()
{ }


void DetectorInfo::ClearGeometry()
{
  
  for(Int_t d=0; d<GetMaxNoDetectors(); d++){
      SetCenterX(d, NAN);
      SetCenterY(d, NAN);
      SetCenterZ(d, NAN);

    SetRotationX(d, NAN);
    SetRotationY(d, NAN);
    SetRotationZ(d, NAN);

        SetSize0(d, NAN);
        SetSize1(d, NAN);
        SetSize2(d, NAN);
        SetSize3(d, NAN);
        SetSize4(d, NAN);

    SetNoStripsX(d, 0);
    SetNoStripsY(d, 0);

         SetName(d, "");
         SetType(d, "");
  }

}

void DetectorInfo::ResetData()
{

  // keep first interaction point information
  // for debugging...
  detData.fIDetID = -1; // first interaction point detector ID

  detData.fIX = NAN; // first interaction point, x coordinate (origin: center of target)
  detData.fIY = NAN;
  detData.fIZ = NAN;
  
  for(Int_t d=0; d<maxDetectors; d++){
    detData.haveHit[d]= 0;
    detData.haveHitID[d]= -1;
    detData.stripX[d] = -1;
    detData.stripY[d] = -1;
    detData.energy[d] = NAN;

    detData.hitPositionX[d]=NAN;
    detData.hitPositionY[d]=NAN;
    detData.hitPositionZ[d]=NAN;
  }

  
  
}





Bool_t DetectorInfo::IsInFront(Int_t id1, Int_t id2)
{
  // this routine checks if detector 1 is closer to the center of target than detector 2
  // if both detectors have the same distance to target, than false is returned

  TVector3 detCenter1, detCenter2;

  detCenter1.SetXYZ(GetCenterX(id1), GetCenterY(id1), GetCenterZ(id1));
  detCenter2.SetXYZ(GetCenterX(id2), GetCenterY(id2), GetCenterZ(id2));

  if(detCenter2.Mag() > detCenter1.Mag()){
    return true;
  }else{
    return false;
  }

}











void DetectorInfo::CalcStripNumbers(Int_t detID, Double_t hx, Double_t hy, Double_t hz, Int_t &stripx, Int_t &stripy)
{
  //printf("inside: calculating strip no x\n");
  TVector3 vhit(hx, hy, hz); // vector from origin to FI
  TVector3 vcen(GetCenterX(detID), GetCenterY(detID), GetCenterZ(detID)); // vector from origin to center detector

  TVector3 vdet=vhit-vcen; // hit with respect to center of detector
  
  vdet.RotateX(GetRotationX(detID));
  vdet.RotateY(GetRotationY(detID));
  vdet.RotateZ(GetRotationZ(detID));

  if(strcmp(det[detID].type.c_str(), "DSSDbox")==0){
    
    TVector3 vori(GetSize0(detID)/2.0, GetSize1(detID)/2.0, 0.0); // vector to move origin to center of lower left corner of detector
  
    TVector3 vpix=vdet+vori;

    stripx = (Int_t)(vpix.X()/(GetSize0(detID)/GetNoStripsX(detID)));
    stripy = (Int_t)(vpix.Y()/(GetSize1(detID)/GetNoStripsY(detID)));

  } else if(strcmp(det[detID].type.c_str(), "DSSDtube")==0){
    
    Double_t phips = GetSize4(detID)/((Double_t)GetNoStripsX(detID));  // angle coverage phi per x strip
    Double_t radps = (GetSize1(detID)-GetSize0(detID))/((Double_t)GetNoStripsY(detID)); // radius coverage per y strip
    
    Double_t hitPhi = vdet.Phi(); // vector: phi = -pi ... pi
    if(hitPhi<0.0){ // detector: phi = 0 ... 2pi
      hitPhi+=2.0*TMath::Pi();
    }
    Double_t hitRad = vdet.Mag();
    
    stripx = (Int_t)((hitPhi-GetSize3(detID))/phips);
    stripy = (Int_t)((hitRad-GetSize0(detID))/radps);
    
    //printf("FIglob: %f %f %f\n", hx, hy, hz);
    //printf("FIdet:  %f %f %f\n", vdet.x(), vdet.y(), vdet.z());
    //printf("hitPhi %f, phi_s %f, phips %f;      hitRad %f, radMin %f, radps %f\n", hitPhi, GetSize3(detID), phips, hitRad, GetSize0(detID), radps);
    //printf("CalcStrip: got strip x = %d, strip y = %d\n\n", stripx, stripy);

  } else {
    stripx = -1;
    stripy = -1;
  }
  
}






void DetectorInfo::CalcHitPosition(Int_t detID, Int_t stripx, Int_t stripy, Double_t &hx, Double_t &hy, Double_t &hz)
{

  TVector3 vcen(GetCenterX(detID), GetCenterY(detID), GetCenterZ(detID)); // vector from origin to center detector
  TVector3 vdet(0.0, 1.0, 0.0);

  if(strcmp(det[detID].type.c_str(), "DSSDbox")==0){

    Double_t pixx = (Double_t)stripx * (GetSize0(detID)/GetNoStripsX(detID)) + (GetSize0(detID)/GetNoStripsX(detID))/2.0;
    Double_t pixy = (Double_t)stripy * (GetSize1(detID)/GetNoStripsY(detID)) + (GetSize1(detID)/GetNoStripsY(detID))/2.0;

    TVector3 vpix(pixx, pixy, 0.0);

    TVector3 vori(GetSize0(detID)/2.0, GetSize1(detID)/2.0, 0.0);

    vdet = vpix-vori;


  } else if(strcmp(det[detID].type.c_str(), "DSSDtube")==0){

    Double_t phips = GetSize4(detID)/((Double_t)GetNoStripsX(detID));  // angle coverage phi per x strip
    Double_t radps = (GetSize1(detID)-GetSize0(detID))/((Double_t)GetNoStripsY(detID)); // radius coverage per y strip

    Double_t hitPhi = (Double_t)(stripx) * phips + GetSize3(detID) + phips/2.0 ;
    Double_t hitRad = (Double_t)(stripy) * radps + GetSize0(detID) + radps/2.0;

    if(hitPhi>TMath::Pi()){
      hitPhi-=2.0*TMath::Pi();
    }

    // vector from origin of detector to center of strip
    vdet.SetMag(hitRad);
    vdet.SetPhi(hitPhi);


  } else {
    printf("<DetectorInfo::CalcHitPosition> invalid geometry type '%s'! Cannot determine position!\n", det[detID].type.c_str());
    abort();
  }

  vdet.RotateZ(-GetRotationZ(detID));
  vdet.RotateY(-GetRotationY(detID));
  vdet.RotateX(-GetRotationX(detID));

  TVector3 vhit=vdet+vcen;

  hx = vhit.X();
  hy = vhit.Y();
  hz = vhit.Z();
}





void DetectorInfo::CalcHitPosition(Int_t detID, Int_t stripx, Int_t stripy)
{
  CalcHitPosition(detID, stripx, stripy, detData.hitPositionX[detID], detData.hitPositionY[detID], detData.hitPositionZ[detID]);
}









void DetectorInfo::Parse(string filename)
{
  
  printf("\n\nParsing of detector geometry file\n");


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
  Int_t counter=0;
  const Int_t stopper=10000;

  const Int_t maxArg=7;
  char temp[maxArg][500];



  while(fin.good())  {

    //reset values
    for(Int_t i=0; i< maxArg; i++){temp[i][0]='\0';}


    getline(fin,line);
    //printf("Parsing line '%s'\n",line.c_str());



    //skip comments
    line.copy(temp[0],2,0); //first two characters
    temp[0][2]='\0';
    //printf("temp1=%s\n",temp1);
    if(strcmp(temp[0],"//")==0){continue;}

    //parse the line
    //write each word to a char
    //delimiter is whitespace
    std::istringstream iss(line);
    for(Int_t i=0; i<maxArg; i++){iss >> temp[i];}

    //skip empty lines
    if(strcmp(temp[0],"")==0){continue;}



    //get input parameter:

    // file names   
    if(strcmp(temp[0],"number_of_detectors")==0)  {
      fNoOfDet=atoi(temp[1]);
      printf("Number of detectors %d\n", fNoOfDet);
    }
    else if(strcmp(temp[0],"position")==0){
      Int_t index = atoi(temp[1]);
      SetCenterX(index, atof(temp[2])*mm);
      SetCenterY(index, atof(temp[3])*mm);
      SetCenterZ(index, atof(temp[4])*mm);
      printf("Got position of detector %d: x = %f, y = %f, z = %f\n", index, GetCenterX(index), GetCenterY(index), GetCenterZ(index));
    }
    else if(strcmp(temp[0],"rotation")==0){
      Int_t index = atoi(temp[1]);
      SetRotationX(index, atof(temp[2])*deg);
      SetRotationY(index, atof(temp[3])*deg);
      SetRotationZ(index, atof(temp[4])*deg);
      printf("Got rotation of detector %d: x = %f, y = %f, z = %f\n", index, GetRotationX(index), GetRotationY(index), GetRotationZ(index));
    }
    else if(strcmp(temp[0],"size")==0){
      Int_t index = atoi(temp[1]);
      SetSize0(index, atof(temp[2])*mm);
      SetSize1(index, atof(temp[3])*mm);
      SetSize2(index, atof(temp[4])*mm);
      SetSize3(index, atof(temp[5])*deg);
      SetSize4(index, atof(temp[6])*deg);
      printf("Got full size of detector %d: 1st (x or r_min) = %f, 2nd (y or r_max) = %f, 3rd (z) = %f, 4th (void or phi_start) = %f, 5th (void or phi_D) = %f\n", index, GetSize0(index), GetSize1(index), GetSize2(index), GetSize3(index), GetSize4(index));
    }
    else if(strcmp(temp[0],"number_of_stripes")==0){
      Int_t index = atoi(temp[1]);
      SetNoStripsX(index, atoi(temp[2]));
      SetNoStripsY(index, atoi(temp[3]));
      printf("Got number of strips of detector %d: x = %d, y = %d\n", index, GetNoStripsX(index), GetNoStripsY(index));
    }
    else if(strcmp(temp[0],"name")==0){
      printf("Info: 'name' not yet implemented (is 'logical').\n");
    }
    else if(strcmp(temp[0],"type")==0){
      Int_t index = atoi(temp[1]);
      SetType(index, temp[2]);
      printf("Got type '%s'.\n", GetType(index).c_str());
    }
    else {
      printf("Cannot read your input keyword '%s'. Aborting program.\n", temp[0]);
      abort();
    }


    //count lines
    counter++;
    if(counter>stopper){
      printf("Reached %d lines in file '%s'! Stopping!\n", counter, filename.c_str());
      abort();
    }
  } // end of reading input file

  
  printf("DetectorInfo::Parse done! Checking input...\n");
  CheckInput();


}












void DetectorInfo::CheckInput()
{
  
  if(fNoOfDet>maxDetectors){
    printf("Error: more detectors in input file (%d) than allowed in DetectorGlobals.hh (%d)\n", fNoOfDet, maxDetectors);
    abort();
  }

  if(fNoOfDet==0){
    printf("Info: number of detectors is 0!\n");
  }else if(fNoOfDet<0){
    printf("Error: invalid (negative) number of detectors: %d\n", fNoOfDet);
    abort();
  }

  for(Int_t d=0; d<fNoOfDet; d++){
    Double_t checkSum=0.0;
    
    if( !((strcmp(det[d].type.c_str(),"DSSDbox")==0) || (strcmp(det[d].type.c_str(),"DSSDtube")==0)) ){
      printf("Error: invalid detector type '%s'!\n", GetType(d).c_str());
      abort();
    }
    
    checkSum = det[d].center[0] * det[d].center[1] * det[d].center[2];
    if(TMath::IsNaN(checkSum)){
      printf("Error: invalid position for detector %d: x = %f, y = %f, z = %f\n", d, det[d].center[0], det[d].center[1], det[d].center[2]);
      abort();
    }

    checkSum = det[d].rotation[0] * det[d].rotation[1] * det[d].rotation[2];
    if(TMath::IsNaN(checkSum)){
      printf("Error: invalid rotation for detector %d: x = %f, y = %f, z = %f\n", d, det[d].rotation[0], det[d].rotation[1], det[d].rotation[2]);
      abort();
    } 

    checkSum = det[d].size[0] * det[d].size[1] * det[d].size[2];
    if(TMath::IsNaN(checkSum) || checkSum==0){
      printf("Error: invalid size for detector %d: 1st = %f, 2nd = %f, 3rd = %f\n", d, det[d].size[0], det[d].size[1], det[d].size[2]);
      abort();
    }
    checkSum = det[d].size[3] * det[d].size[4];
    if( (strcmp(det[d].type.c_str(),"DSSDtube")==0) && TMath::IsNaN(checkSum) ){
      printf("Error: invalid angle for type 'DSSDtube': phi_start = %f, phi_D = %f\n", det[d].size[3], det[d].size[4]);
      abort();
    }
    
    if( (det[d].noStrips[0] < 0) || (det[d].noStrips[1] < 0) ){
      printf("Error: invalid number of strips in detector %d: %d x, %d y (have to be nonnegative; 0 and 1 is same geometry)\n", d, det[d].noStrips[0], det[d].noStrips[1]);
      abort();
    }

  }



  //for testing:
  printf("CheckInput completed sucessfully!\n");
  //abort();
}



