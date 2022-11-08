#include <TF1.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TTimer.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TMath.h>
#include <TH1D.h>

#include <TROOT.h>
#include <TApplication.h>


#include <stdio.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

Double_t sigma[] = { 660, 920, 1400, 1900, 2500, 3700, 4800, 6200, 8700, 11000, 14000, 19000, 24000, 30000, 39000, 48000, 59000, 75000 }; //pb cc cross section cooper-sarkar 2011
Double_t sigmaNC[] = { 240, 350, 530, 730, 980, 1400, 1900, 2400, 3400, 4400, 5600, 7600, 9600, 12000, 16000, 20000, 24000, 31000 }; //pb
Double_t Esig[] = { 1e6, 2e6, 5e6, 1e7, 2e7, 5e7, 1e8, 2e8, 5e8, 1e9, 2e9, 5e9, 1e10, 2e10, 5e10, 1e11, 2e11, 5e11 }; //GeV

//make plot of length distribution of showers
//make plot of azimuth and elevation distribution as function of distance and wi
//need good estimate of length of shower


//Ghandi 1996
//Double_t sigma[] = { 634, 960, 1412, 1749, 2554, 3630, 4436, 6283, 8700, 10490, 14660, 20100, 23790, 32890, 44270, 53570, 73200, 99270, 117900 }; //pb
//Double_t sigmaNC[] = { 260, 402, 600, 748, 1104, 1581, 1939, 2763, 3837, 4641, 6490, 8931, 10660, 14650, 19950, 23770, 32470, 43770, 51960 }; //pb
//Double_t Esig[] = { 1e6, 2.5e6, 6e6, 1e7, 2.5e7, 6e7, 1e8, 2.5e8, 6e8, 1e9, 2.5e9, 6e9, 1e10, 2.5e10, 6e10, 1e11, 2.5e11, 6e11, 1e12 }; //GeV


int iColors[] = {kBlue-3,kCyan-3,kGreen-3,kYellow-3,kRed-3,kRed+3,kMagenta-3};
int marker[] = { 23, 22, 29, 21, 20, 28, 25};
double markerSize[] = { 0.9, 0.9, 1.2, 0.7, 0.9, 1.3, 0.9};

TF1 *fPE = 0;

Double_t c = 299792; //km/s
Double_t pi = 3.14159265359;
Double_t DecayTime = 0.290e-12; //s
Double_t Mtau = 1.7768; //GeV
Double_t REarth = 6371; //km

//All the coefficients to get the right PE intensity
///////////////////////////////

Int_t iConfig = 2;
double DetectorAltitude[] = { 0, 1, 2, 3};

// Timing information
Int_t iTrigWin = 2;
double TriggerWindow[] = {15.00, 30.00, 45.00, 60.00, 75.00}; // ns
 
//obtained from 3e4 GeV gamma rays
double lincorr[] = { 0.0509251, 0.0522854, 0.0595455, 0.0642221};
double scalefirst[] = { 1.00001, 1.03346, 1.53535, 2.36961};
double eleScaling[] = { 0.00163848, 0.00191408, 0.0071185, 0.0513182};
double absorptionlength[] = { 16.7049, 16.6106, 17.9808, 19.0274};
//Parameterization of PE distribution at 50km, 0ele, and 0 Altitude
//par[0]*exp(-xx/par[1])+par[2]*exp(-xx/(par[3]+xx*par[4]))
double parPEF[] = { 0.000332267, 0.580756, 4.25751e-05, 1.9491, 0.0427249};

/*
//obtained from 1e6 GeV gamma rays
double lincorr[] = { 0.0497124, 0.0493805, 0.0552341, 0.0604423};
double scalefirst[] = { 1.00001, 1.02504, 1.42546, 1.73371};
double eleScaling[] = { 0.00143344, 0.00164536, 0.00419904, 0.00479111};
double absorptionlength[] = { 16.7647, 16.8572, 18.3202, 19.4597};
//Parameterization of PE distribution at 50km, 0ele, and 0 Altitude
//par[0]*exp(-xx/par[1])+par[2]*exp(-xx/(par[3]+xx*par[4]))
double parPEF[] = { 0.00038827, 0.555588, 4.66631e-05, 1.90266, 0.0426453};
*/

// Parameterization of PE distribution with timing at 40km,0 ele,and 2km altitude
//obtained from 1e6 GeV gamma rays with various ns trigger window (15,30,45,60,75)
double DetectorAltitude2[] = { 2, 2, 2, 2, 2};
double lincorr2[] = { 0.04, 0.04, 0.04, 0.04, 0.0400001};
double scalefirst2[] = { 1, 1, 1, 1, 1};
double eleScaling2[] = { 0.00499707, 0.00501046, 0.0049938, 0.00499608, 0.00500326};
double absorptionlength2[] = { 20, 20, 20, 20, 20};
double parPEF2[5][5] = {{ 6.77183e-05, 0.599578, 0.000109705, 1.04491, 0.0346329}, { 2.48094e-09, 0.971087, 0.000211414, 1.04562, 0.0346362}, { 4.93947e-09, 0.384049, 0.000275539, 1.07095, 0.0344304}, { 2.52288e-10, 0.379361, 0.000306272, 1.10666, 0.0341286}, { 2.03024e-09, 0.346843, 0.000314862, 1.14685, 0.0337836}};

// Parameterization of PE distribution with timing at 80km,0 ele,and 2km altitude
//obtained from 1e6 GeV gamma rays with various trigger windows (15,30,45,60,75)
double DetectorAltitude3[] = { 2, 2, 2, 2, 2};
double lincorr3[] = { 0.0792027, 0.0786967, 0.0782509, 0.0779326, 0.0777705};
double scalefirst3[] = { 1, 1, 1, 1, 1};
double eleScaling3[] = { 0.00169174, 0.00170477, 0.00172341, 0.00173431, 0.0017268};
double absorptionlength3[] = { 20, 20, 20, 20, 20};
double parPEF3[5][5] = {{ 4.8878e-05, 0.765578, 7.53198e-06, 1.22237, 0.0324224},{ 3.41281e-05, 0.911229, 1.30399e-05, 1.24591, 0.032254}, { 1.84094e-05, 1, 2.07172e-05, 1.23037, 0.0323958}, { 1.0245e-09, 0.321118, 3.32187e-05, 1.19095, 0.0326761}, { 2.21387e-05, 0.138689, 3.20128e-05, 1.24981, 0.0321155}};

///////////////////////

Double_t dMinEnu = 8.5;
Double_t dMaxEnu = 9.5;
Double_t dST = 10; //km max height of shower tip above ground;
Double_t yMin = 1;
Double_t yMax = 60;
Double_t yDelta = 1;
Double_t MaxElevation = 10; //elevation angle (determines path through Earth;
Double_t DeltaAngle = 0.1; //steps in azimuth and elevation 
Double_t DeltaAngleAz = 0.3; //steps in azimuth  
Double_t nuIndex = 2; //power law index of the neutrino spectrum the minus sign is added later
Double_t dMaxCherenkovAzimuthAngle = 40.0; //maximum azimuth angle for cherenkov 
Double_t dMaxFluorescenceDistance = 100;

//next three parameters are key to the instrument
//plus the telescop height above ground which can be selected above
Double_t tanFoV = tan(5/180.*pi); //Field of view of telescope above the horizon
Double_t dFoVBelow = 5/180.*pi; //Field of view of telescope  below horizon 
Double_t dMinLength = 0.3; //mimnimum length a shower has to have in the camera, in degrees. This is a conservative estimate because it assumes that the shower starts at a distance l from the detector, which is not necessarily tru for showers with shallow elevation angles.
Double_t dMinimumNumberPhotoelectrons = 20; 

Int_t iMirrorSize = 1;
Double_t dMirrorA[] = {1.0, 5.0, 10.0, 100.0}; //m^2 
//Double_t dThreshold[] = {8*3, 19*3, 22*3, 120*3}; //pe //three fold coincidence. 
//Double_t dThreshold[] = {10*2, 22*2, 24*2, 155*2}; //pe //two fold coincidence
Double_t dThreshold[] = {10*2, 22*2, 24*2, 155*2}; //pe //two fold coincidence

Bool_t bFluorescence = kFALSE;
Bool_t bCombined = kFALSE;
Bool_t bMonoNu = kFALSE; //simulate monoenergetic neutrinos, only good for acceptance calculation. For all other simulations set it to kFALSE
TGraph *grsCC;
TGraph *grsNC;

TH1D *hTriggeredAzimuthAngles;

string Hold()
{
  string input;

  //hold the code;
  TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
  timer.TurnOn();
  cout<<"Press Enter to continue:"<<endl;
  getline(cin,input);
  timer.TurnOff();
  return input;
  
}

Double_t myPEfunction(Double_t *x, Double_t *par)
{
   //par[0] Distance to where tau comes out in km
   //par[1] Elevation in rad
   //azimuth angle is our x
   Float_t xx =x[0]; //angle is rad here 
   //if angle is larger 40 degrees return 0
   if(xx>0.69813170)
     return 0;

   //Calculate azimuth angle in frame of master pe distribution (50km, 0ele,
   //0altitude
   Double_t dTelAngle = atan(DetectorAltitude[iConfig]*1e-3/par[0]);
   Double_t dAngle = sqrt(xx*xx + (par[1]-dTelAngle)*(par[1]-dTelAngle))*57.295780; //in deg
   //calculate how many PEs / per m2 per GeV
   Double_t f = 0;
   if(dAngle<1.3)
     f = parPEF[0]*exp(-1.1/parPEF[1])+parPEF[2]*exp(-1.1/(parPEF[3]+1.1*parPEF[4]));
   else
     f = parPEF[0]*exp(-dAngle/parPEF[1])+parPEF[2]*exp(-dAngle/(parPEF[3]+dAngle*parPEF[4]));

//other parameters to get PE intensity for different distance, azimuth and
//elevation
   //scale PE distribution to first PE distribution at 50km distance
   f*=   scalefirst[iConfig];
   //Get elevation dependence
   f*= (2-exp(-par[1]/eleScaling[iConfig]));
   //Get Distance dependence
   f*=  exp(-(par[0]-55)/
           (absorptionlength[iConfig]+(par[0]-55)*lincorr[iConfig])); //55km is the distance for which the normalized PE distribution is extracted

   return f;
}

/*
 Modeled after Nepomuk's PE function. This parameterization includes trigger timing using a 15,30,45 ns trigger window around the median of the timing distribution. This function is only valid for 2km and uses two master distributions. One taken at 45 km and 85km. 
*/ 
Double_t myPEfunction_trigwindow(Double_t *x, Double_t *par)
{
  // par[0] is distance 
  Double_t l = par[0];
  // par[1] is elevation 
  Double_t elv = par[1];
  // par[2] is trig time window index
  Int_t ti = par[2];
  // x is azimuth angle in rad
  Double_t azi = x[0];
  
  //cout << par[0]<< " " << par[1]<< " " << par[2] << endl;
  //cout << "ti" << ti << "TW: " << TriggerWindow[ti] << " Fit Parameters: " << parPEF2[ti][0] << endl;
  
  if (azi > 0.69813170) // if > 40 degrees then return 0
    return 0; 
  
  
  Double_t dTelAngle = atan(DetectorAltitude[iConfig]*1e-3/l);
  Double_t dAngle = sqrt(azi*azi + (elv-dTelAngle)*(elv-dTelAngle))*57.295780; //in deg
  //calculate how many PEs / per m2 per GeV
  Double_t f = 0;

  //Calculate azimuth angle in frame of master pe distribution (40km, 0ele,2km altitude)
  // if l < 50 km
  if (l < 50){
    if(dAngle<1.3)
      f = parPEF2[ti][0]*exp(-1.1/parPEF2[ti][1])+parPEF2[ti][2]*exp(-1.1/(parPEF2[ti][3]+1.1*parPEF2[ti][4]));
    else
      f = parPEF2[ti][0]*exp(-(dAngle+0.19*azi)/parPEF2[ti][1])+parPEF2[ti][2]*exp(-(dAngle+0.19*azi)/(parPEF2[ti][3]+(dAngle+0.19*azi)*parPEF2[ti][4]));
    
    //scale PE distribution to first PE distribution at 50km distance
    f*=  3.3*scalefirst2[ti]; // Scale is 1 at master distribution height
    //Get elevation dependence
    f*= (2-exp(-elv));
    //Get Distance dependence
    f*=  exp(-(l-45)/
	     (absorptionlength2[ti]+(l-45)*lincorr2[ti])); //40km is the distance for which the normalized PE distribution is extracted
  } else {
    //Calculate azimuth angle in frame of master pe distribution (80km, 0ele,2km altitude)
    // if >= 50 km
    if(dAngle<1.3)
      f = parPEF3[ti][0]*exp(-1.1/parPEF3[ti][1])+parPEF3[ti][2]*exp(-1.1/(parPEF3[ti][3]+1.1*parPEF3[ti][4]));
    else
      f = parPEF3[ti][0]*exp(-dAngle/parPEF3[ti][1])+parPEF3[ti][2]*exp(-dAngle/(parPEF3[ti][3]+dAngle*parPEF3[ti][4]));
    
    //scale PE distribution to first PE distribution at 50km distance
    f*=   scalefirst3[ti]; // Scale is 1 at master distribution height
    //Get elevation dependence
    f*= (2-exp(-elv/eleScaling3[ti]));
    //Get Distance dependence
    f*=  exp(-(l-85)/
	     (absorptionlength3[ti]+(l-85)*lincorr3[ti])); //80km is the distance for which the normalized PE distribution is extracted
    }
  
  
  return f; // return number of photoelectrons arriving at camera in trigger window
}  

Double_t DistanceThroughEarth(Double_t y, Double_t elevation, Double_t azimuth)
{

  elevation = elevation/180*pi; //elevation angle (determines path through Earth;
  azimuth = azimuth/180.*pi;  //azimuth angle

  Double_t l = y; //Distance from detector to where the tau comes out detector is always at z=0

  Double_t v = sqrt((REarth+DetectorAltitude[iConfig])*(REarth+DetectorAltitude[iConfig])-REarth*REarth);

  //shortest distance d between tau trajectory and detector
  Double_t nproj = y*sqrt( 1 + tan(azimuth)*tan(azimuth) ); //projection of trajectory to x-y plane
  Double_t denomsquared= y*tan(azimuth)*y*tan(azimuth) + y*y + nproj*nproj*tan(elevation)*tan(elevation)  ;
  
  //normalized trajectory vector of tau
  Double_t dNormalize = y/sqrt(denomsquared);
  //Double_t dNx = dNormalize * tan(azimuth);
  Double_t dNy = -dNormalize;
  Double_t dNz = dNormalize * sqrt( 1 + tan(azimuth)*tan(azimuth) ) * tan(elevation);


  Double_t p = 2 * ( REarth*dNz - (v-l)*dNy );
  Double_t q = (v-l)*(v-l);

  if(q-p*p/4>=0) //trajectory does not intersect with Earth
      return 0;
        
  Double_t i1 = p/2. - sqrt(p*p/4.-q);
  Double_t i2 = p/2. + sqrt(p*p/4.-q);

return fabs(i2-i1);

}


string star ;
double number;
int distanceNumber;
vector<double> enerNu,enerTau,prob,dist,EtauNorm;
unsigned rem = 0 ;
void removeDuplicates()
{
	while(rem < enerNu.size()-1){
		if(enerNu[rem]==enerNu[rem+1]){
			enerNu.erase(enerNu.begin()+rem+1);
			removeDuplicates();
		}else if(enerNu[rem+1]==enerNu[rem+2]){
			rem ++;
			removeDuplicates();
		}
	}
}
void readFromTable(){
	ifstream ifs("table_with_e_05_a_1.txt") ;
	if(ifs.is_open()){
	ifs>>star;
		while(ifs.good()){
			ifs>>number;
			enerNu.push_back( pow(10,number-9.0) );
			ifs>>number;
			//angle.push_back(number);
			dist.push_back(cos((180-number)*pi/180)*2*REarth);
			for(int i=0;i<100;i++){
				ifs>>number;
				enerTau.push_back( pow(10,number-9.0) );
                                EtauNorm.push_back( pow(10,4+i*0.07) );
				ifs>>number;
				prob.push_back(number);
			}
			ifs>>number;
			enerTau.push_back( pow(10,number-9.0) );
                        EtauNorm.push_back( pow(10,4+100*0.07) );
			ifs>>star;
		}
		cout << "size: " << enerNu.size() << endl;
		removeDuplicates() ;
	}
}
void findAngleNumber(){
	for(unsigned i=1;i<enerNu.size();i++){
		if(dist[0]==dist[i]){
			distanceNumber = i;
			break;
		}else{
			continue ;
		}
	}
}
double biLinearInterpolation(double a1,double n1,double q11,double a2,double n2,double q22,double x,double y){
	double q12 = (q11 + q22)/2 ;
	double q21 = q12 ;
	double p = (n2-y)/(n2-n1)*( (a2-x)/(a2-a1)*q11 + (x-a1)/(a2-a1)*q21 ) + (y-n1)/(n2-n1)*( (a2-x)/(a2-a1)*q12 + (x-a1)/(a2-a1)*q22 ) ;

return p ;
}

//Find index in lookuptable
int FindLion(double dValue, vector<double> &vData,int iSize)
{
          int iWidth = iSize/2;
          int index = iWidth;
          //cout<<dValue<<": ";
          while(iWidth>1 && vData[index]!=dValue && index > 0  && index <iSize)
            {
              iWidth = iWidth/2+ iWidth%2;
              index = vData[index]<dValue ? index + iWidth :  index - iWidth;
              //cout<<index<<" "<<iWidth<<"  "<<vData[index]<<"; ";
            }
          while((dValue>vData[index] || index<0) && index < iSize)
             index++;
          index--;

          if(dValue==vData[index] || index>=iSize)
             index--;

          if(index<0)
            index=0; 
 
         //cout<<index<<endl;

return index;
}
//Calculates the probability of tau emergence using NuTauSim LUT
Double_t PEtauNTauSim(Double_t D,Double_t Etau, Double_t Enu, TH1D *hTau) 
{
	//Enu = log10(Enu) + 9.0;
	//Etau = log10(Etau) + 9.0 ;
	//double zenithAngle = 180 - acos(D/2/REarth)/M_PI*180 ;
//	if(zenithAngle >= angle[0] && zenithAngle <= angle[angleNumber-1] && Enu>=enerNu[0] && Enu <=enerNu[enerNu.size()-1]){
	if(D >= dist[0] && D <= dist[distanceNumber-1] && Enu>=enerNu[0] && Enu <=enerNu[enerNu.size()-1] && Etau>=enerTau[0] && Etau <=enerTau[100]){

         int indexEnu = FindLion(Enu,enerNu,enerNu.size());

	 int indexDistance = FindLion(D,dist,distanceNumber);

	 int indexEtau = FindLion(Etau,enerTau,100);
		
                int indexProb1 = indexEnu*distanceNumber*100+indexDistance*100+indexEtau;
		double p1 = prob[indexProb1] ;
		int indexProb2 = (indexEnu+1)*distanceNumber*100 + (indexDistance+1)*100 + indexEtau ;
		double p2 = prob[indexProb2] ;

		double Prob = biLinearInterpolation(dist[indexDistance],enerNu[indexEnu],p1,dist[indexDistance+1],enerNu[indexEnu+1],p2,D,Enu)/
						(EtauNorm[indexEtau+1]-EtauNorm[indexEtau]);
						//(pow(10,4+(indexEtau+1)*0.07)-pow(10,4+indexEtau*0.07));
		//cout<<Prob<<endl;
		return Prob ;
	}else{
		return 0 ;
	}
}



//Probability that Tau with Energy Etau emerges for initial nu energy Enu
//Thickness of matter d 
//Does not use energy loss of tau in Earth.
//Assumes the energy of the tau is 0.8*Enu
//Double_t PEtauNoTauEnergyLoss(Double_t D,Double_t Etau, Double_t Enu, TH1D *hTau) 
//Double_t PEtauNoTauEnergyLoss(Double_t D,Double_t Etau, Double_t Enu, TH1D *hTau) 
Double_t PEtau(Double_t D,Double_t Etau, Double_t Enu, TH1D *hTau) 
{

int n = hTau->FindBin(Etau);
if(hTau->GetBinLowEdge(n)>0.8*Enu  || hTau->GetBinLowEdge(n+1)<0.8*Enu )
  return 0;

Double_t sCC = grsCC->Eval(Enu); //crossection in pB
Double_t sNC = grsNC->Eval(Enu); //crossection in pB
Double_t rho = 2.65; //density in g/cm3
Double_t NA = 6.022142e23;

Double_t dInvConvCC = sCC*rho*NA*1e-31; // 1/km 1e-12*1e-28*1e4*1e5
Double_t db = (sNC+sCC)*rho*NA*1e-31; // 1/km
Double_t da = Mtau/(DecayTime*c*0.8*Enu); //1/km

//cout<<dInvConvCC<<endl;
//cout<<"neutrino interaction: "<<db<<" "<<exp(-1.0*db*D)<<endl;
//cout<<"tau survival: "<<da<<"  "<<exp(-1.0*da*D)<<endl;

Double_t Prob = dInvConvCC / (da-db) * ( exp(-1.0*db*D)-exp(-1.0*da*D) );  

if(Prob<0)
  return 0;

Prob /= (hTau->GetBinLowEdge(n+1)-hTau->GetBinLowEdge(n));

return Prob;  
}

//Probability that Tau with Energy Etau emerges for initial nu energy Enu
//Thickness of matter d 
//follows description in Dutta 2005 in particular equation 28 with
//parameterization of beta in equation 13 case II
//energies in GeV distances in km at inptut
//Fails <1e8 GeV because energy loss (Beta) becoms <0
Double_t PEtauDutta(Double_t D,Double_t Etau, Double_t Enu,TH1D *hTau) 
{

if(Etau>0.8*Enu)
  return 0;

Double_t sCC = grsCC->Eval(Enu); //crossection in pB
Double_t sNC = grsNC->Eval(Enu); //crossection in pB
Double_t rho = 2.65; //density in g/cm3
Double_t NA = 6.022142e23;

Double_t dInvConvCC = sCC*rho*NA*1e-31; // 1/km 1e-12*1e-28*1e4*1e5
Double_t dInvConvtotal = (sNC+sCC)*rho*NA*1e-31; // 1/km
Double_t beta = 1.2e-6 + 0.16e-6 * log(Etau/1e10); //cm2/g Equation 13 case II in Dutta 
//cout<<"beta "<<beta<<endl;
Double_t prefactor = Mtau/(DecayTime*c*1e5*beta*rho*Etau); //dimensionless
//cout<<"Prefactor: "<<prefactor<<endl;

Double_t xDelta = log(Etau/(0.8*Enu))/(beta*rho)*1e-5+D; //where delta function is non zero; 1e-5 convert from cm to km 

if(xDelta<0)
  return 0;

//cout<<"beta: "<<beta<<" rho: "<<rho<<" Etau: "<<Etau<<endl;
Double_t Prob = 1.e-5/(beta*rho*Etau); //km/GeV
//cout<<"Prob: "<<Prob<<endl;
Double_t Pnu = dInvConvCC*exp(-xDelta*dInvConvtotal);  // 1/km
//cout<<"Pnu: "<<Pnu<<endl;
Prob *= Pnu; 
//cout<<"Prob2: "<<Prob<<endl;
if(Prob<0)
  return 0;

if(prefactor<0)
  return 0;
  
Double_t Ptau = exp(-prefactor*(1.0-exp(-beta*rho*(D-xDelta)*1e5))); //1e5 to convert from km to cm 
//cout<<"Ptau: "<<Ptau<<endl;
//cout<<1.0-exp(-beta*rho*(D-xDelta)*1e5)<<endl;
Prob *= Ptau;

return Prob;  
}

Double_t PDecayFluorescence(Double_t Etau, Double_t y, Double_t elevation, Double_t azimuth)
{
  //cout<<endl<<"elevation: "<<elevation<<" azimuth: "<<azimuth<<" distance: "<<y<<endl;
  elevation = elevation/180*pi; //elevation angle (determines path through Earth;
  azimuth = azimuth/180.*pi;  //azimuth angle

  Double_t l = y; //Distance from detector to where the tau comes out detector is always at z=0

  //shortest distance d between tau trajectory and detector
  Double_t nproj = y*sqrt( 1 + tan(azimuth)*tan(azimuth) ); //projection of trajectory to x-y plane
  Double_t denomsquared= y*tan(azimuth)*y*tan(azimuth) + y*y + nproj*nproj*tan(elevation)*tan(elevation)  ;
  
  //normalized trajectory vector of tau
  Double_t dNormalize = y/sqrt(denomsquared);
  Double_t dNx = dNormalize * tan(azimuth);
  Double_t dNy = -dNormalize;
  Double_t dNz = dNormalize * sqrt( 1 + tan(azimuth)*tan(azimuth) ) * tan(elevation);

  if(azimuth>=pi/2.)
    {
      dNx*=-1;
      dNy*=-1;
    }

  //cout<<"trajectory vector normalized x: "<<dNx<<" y: "<<dNy<<" z: "<<dNz<<" normalization: "<<dNormalize<<endl;

  //crossproduct of trajectory vector and vector of where tau emerges. Gives the
  //distance between the to perpendicular to the trajecotory vector
  //Double_t dx = y*dNz; 
  //Double_t dy = 0;
  //Double_t dz = y*dNx;

  //Double_t d = sqrt(dx*dx+dy*dy+dz*dz); //shortest distance d between tau trajectory and detector

  //Double_t dem = sqrt(l*l-d*d);

  //maximum length of trajectory above horizon befor track leaves atmosphere (dST above ground)
  //calculation is not entirely correct but we do not max out on this distance
  //anyway
  Double_t v = sqrt((REarth+DetectorAltitude[iConfig])*(REarth+DetectorAltitude[iConfig])-REarth*REarth);

  Double_t phi = elevation + asin( REarth/sqrt(REarth*REarth+(l-v)*(l-v)) ); //if azimuth > 90

  Double_t alpha = asin( sin(phi) * sqrt(REarth*REarth+(l-v)*(l-v)) / (REarth+dST)  );

  Double_t gamma = pi - alpha - phi;

  Double_t dMaxDist = (REarth+dST) *  sin(gamma)/sin(phi);


 //trim the path length to be fully inside the atmosphere, should always be the
 //case
  Double_t dED = 0; //extra distance
  Double_t dInFoV=0; //Distance in between lower edge of FoV and line of sight to horizon

  //If the camera FoV below the horizon is larger than the maximum angle needed to
  //cover all the solid angle below the horizon, set the FoV below the horizon
  //to the maximum angle needed.
  Double_t dMaxFoVBelow = asin(REarth/(REarth+DetectorAltitude[iConfig]));
  if(dFoVBelow>dMaxFoVBelow)
      dFoVBelow=dMaxFoVBelow; 


  //length of the visible trajectory between the plane to horizon and the lower edge of the camera
  //FoV below the horizon
  dInFoV = l * sin(dFoVBelow) / ( dNy * sin(dFoVBelow) +  dNz * cos(dFoVBelow)  );

  //length of tau trajectory below the horizon and earth surface
  Double_t p = 2 * ( REarth*dNz - (v-l)*dNy );
  Double_t q = (v-l)*(v-l);

  if(q-p*p/4>=0) //trajectory does not intersect with Earth
      return 0;
        
  Double_t i1 = p/2. - sqrt(p*p/4.-q);
  Double_t i2 = p/2. + sqrt(p*p/4.-q);
            
  if(i1<0 || i2<0)
    cout<<"PDecay: i1 or i2 less than 0: "<<i1<<"  "<<i2<<endl;

  dED= i1<i2 ? i1 : i2;

  if( (dInFoV>dED && dInFoV>0 ) || dInFoV<0)
          dInFoV=dED;
 
  if(l>v) //if shower emerges beyond horizon
     dInFoV = 0;


  dMaxDist = dST/sin(elevation) + dInFoV; //takes into account distance visible below plane to horizon
 // cout<<"dMaxDist: "<<dMaxDist<<endl;

  Double_t dTermInSquareRoot = cos(elevation)*cos(azimuth)*cos(elevation)*cos(azimuth)
                               + sin(elevation)*sin(elevation)/tanFoV/tanFoV
                               - cos(elevation)*cos(elevation);

  if(dTermInSquareRoot>0) //if it is negative the shower is contained in the FoV anywhere along the track
    {
        //maximum trajectory length before the tip of the shower is not contained in the camera anymore
       Double_t dMaxDisttoSatisfyFovReq = l / ( cos(elevation)*cos(azimuth) + sqrt( dTermInSquareRoot ));
       //if the value is negative, the elevation angle is less than the Max FoV would have to point below the
       //horizon
       if(dMaxDist>dMaxDisttoSatisfyFovReq+dInFoV && dMaxDisttoSatisfyFovReq>0) //correct max. trajectory length
           dMaxDist = dMaxDisttoSatisfyFovReq+dInFoV;
   //   cout<<" dMaxDisttoSatisfyFovReq: "<<dMaxDisttoSatisfyFovReq<<endl;
    }

  
  //length of shower in camera plane
  Double_t dShwrLgth = 0.304 * log(Etau*0.5/0.088)/log(2); //0.304km radiation length at see level, 0.088GeV critical energy of electrons in air,only 0.5 of the energy goes into the electromagnetic shower
  //cout<<"Etau "<<Etau<<" size of shower in km: "<<dShwrLgth<<endl;
  //cout<<" dMaxDist: "<<dMaxDist<<endl;
  
  //ok have taken all requirements into account. Lets see if the trajectory
  //length allows for a full shower development. If not we quit.
  if(dMaxDist<dShwrLgth)
   return 0;


  //make sure the shower does not develop past the point where more than 90% of the
  //taus have decayed
  Double_t DecayLength = Etau * c * DecayTime / Mtau;
  Double_t d90PctDecayLength = -log(0.1)*DecayLength;
  //cout<<"90% of taus decayed after: "<<d90PctDecayLength<<endl;

  if(d90PctDecayLength+dShwrLgth<dMaxDist)
    dMaxDist = d90PctDecayLength+dShwrLgth;

  //cout<<" dMaxDist: "<<dMaxDist<<endl;
  //calculate how far away from the detector the shower can be to be still
  //detected
  //emitted light intensity
  Double_t dLight = 5.95844e3*Etau; //in photons. 5.958 comes from the macro FluorescenceDetectionYield.C and includes PDE of S14520-6050CN, it is the integral from 300 to 430nm 
  dLight /= 4 * pi; //so we do not have to do it in every loop below
  //and absorption 0.9 at 337nm is used in condition below
  Double_t dFluorescenceMaximumDistance = 10; //km

  while(1)
   {
      //the 1e-6 is for a 1m^2 mirror
      if(dLight*1e-6/dFluorescenceMaximumDistance/dFluorescenceMaximumDistance*exp(dFluorescenceMaximumDistance*log(0.9))>dMinimumNumberPhotoelectrons)
      dFluorescenceMaximumDistance++; 
     else
      break;
     //cout<<dFluorescenceMaximumDistance<<"   "<<dLight*1e-6/dFluorescenceMaximumDistance/dFluorescenceMaximumDistance*exp(dFluorescenceMaximumDistance*log(0.9))<<endl;

   }  
  dFluorescenceMaximumDistance--;
  //cout<<"Maximum Distance between Detector and Shower: "<<dFluorescenceMaximumDistance<<endl;
  //don't know if this is actually good
 
  //check if we need to increase limit 
  if(dMaxDist>dFluorescenceMaximumDistance) //too speed up calculations
    dMaxDist=dFluorescenceMaximumDistance;
  

  //now lets check if the size of the shower is fullfilling the minimum length
  //requirement

  //calculating length of shower in the camera assuming the shower happens late
  //and develops up to the maximumg possible point along the trajectory

  Double_t dLength = 0.0;
  Double_t m = dMaxDist;
  Double_t B = sqrt( m* dNx * m* dNx + (m*dNy+y) * (m*dNy+y) + m*dNz * m*dNz  );
  Double_t A = 0.0;


  //This is if the shower passed and develops behind
  while(dMaxDist>dShwrLgth && (dLength<dMinLength || B>dFluorescenceMaximumDistance ))
   {
     Double_t n = dMaxDist - dShwrLgth -dInFoV;
              m = dMaxDist;

     Double_t A = sqrt( n* dNx * n* dNx + (n*dNy+y) * (n*dNy+y) + n*dNz * n*dNz  );
              B = sqrt( m* dNx * m* dNx + (m*dNy+y) * (m*dNy+y) + m*dNz * m*dNz  );

     Double_t costheta = n*dNx * m*dNx + (n*dNy+y)*(m*dNy+y) + m*dNz * n*dNz;
     costheta = costheta / (A*B);
     dLength = acos(costheta)*180/pi;
 //    cout<<"l"<<l<<"az"<<azimuth*180/pi<<"size of shower in degrees: "<<dLength<<" cos of angle:  "<<costheta<<" dMaxDist "<<dMaxDist<<"dShwrLgth "<<dShwrLgth<<" minimum shower length in deg  "<<dMinLength<<" Distance of tip of shower to telescope  "<<B<<endl;
     if((dLength<dMinLength && B>l) || B>dFluorescenceMaximumDistance )
       dMaxDist-=1.0;
     else
       break;
   }
  dMaxDist+=1.0;

  if(dLength<dMinLength  || B>dFluorescenceMaximumDistance || A>dFluorescenceMaximumDistance) //shower is not long enough in the camera or either end of the shower does not produce sufficient intensity in the telescope
      return 0;

//  need to check if the shower is pointing away from the telescope and if we find a distance in which the tau can decay and develop a shower which appears larger. Need to adjust the distance so the image has the minimal required size.

  //Ok finally we are there. Lets decay the tau in the remaining distance we
  //have
 
  Double_t ProbTauDecay = exp(-(dED-dInFoV)/DecayLength); //Tau has to to survive before it becomes visible to the detector

  ProbTauDecay *= 1-exp(-(dMaxDist-dShwrLgth)/DecayLength); // then it has to decay before it is out of the FoV

  ProbTauDecay*=0.8;//only 80% of taus make a shower

//if(azimuth*180/pi>90 && azimuth*180/pi<8.1)
//cout<<"Etau: "<<Etau<<" el: "<<elevation*180/pi<<" az: "<<azimuth*180/pi<<" alpha "<<alpha<<" beta "<<180/pi*(pi - asin(sinbeta))<<" Prob: "<<ProbTauDecay<<" l "<<l<<" v "<<v<<" below horizon: "<<dBH<<" MaxDistOfTrack: "<<dMaxDist<<" minimal distance between trajectory and telescope "<<d<<endl;
return ProbTauDecay;
}



Double_t PDecay(Double_t Etau, Double_t y, Double_t elevation, Double_t azimuth)
{
  elevation = elevation/180*pi; //elevation angle (determines path through Earth;
  azimuth = azimuth/180.*pi;  //azimuth angle

  Double_t l = y; //Distance between the detector and the point where the tau emerges from the ground. The detector is always at z=0
  //Distance between telescope and horizon
  Double_t v = sqrt((REarth+DetectorAltitude[iConfig])*(REarth+DetectorAltitude[iConfig])-REarth*REarth);

  //Below: the shortest distance d between tau trajectory and detector d
  Double_t nproj = y*sqrt( 1 + tan(azimuth)*tan(azimuth) ); //projection of trajectory to x-y plane
  Double_t denomsquared= y*tan(azimuth)*y*tan(azimuth) + y*y + nproj*nproj*tan(elevation)*tan(elevation)  ;
  
  //normalized trajectory vector of tau
  Double_t dNormalize = y/sqrt(denomsquared);
  Double_t dNx = dNormalize * tan(azimuth);
  Double_t dNy = -dNormalize;
  Double_t dNz = dNormalize * sqrt( 1 + tan(azimuth)*tan(azimuth) ) * tan(elevation);

  //cout<<"trajectory vector normalized x: "<<dNx<<" y: "<<dNy<<" z: "<<dNz<<" normalization: "<<dNormalize<<endl;

  //crossproduct of the trajectory vector with the vector pointing to where the tau emerged from the ground. 
  //The magnitude of the cross product gives the distance of closest approach of the tau to the detector as it travels along 
  //its trajectory
  Double_t dx = y*dNz; 
  Double_t dy = 0;
  Double_t dz = y*dNx;

  Double_t d = sqrt(dx*dx+dy*dy+dz*dz); //shortest distance d between tau trajectory and detector

  //add some extra distance if the tau passes plane in between telescope and
  //horizon
  //using dED term assumes we see shower in camera but lower edge of FoV is aligned with line of sight to the horizon
  Double_t dED = 0; //extra distance
  Double_t dInFoV=0; //Distance in between lower edge of FoV and line of sight to horizon

  //Reset FoV below to maximum possible if it is larger
  Double_t dMaxFoVBelow = asin(REarth/(REarth+DetectorAltitude[iConfig]));
  if(dFoVBelow>dMaxFoVBelow)
      dFoVBelow=dMaxFoVBelow; 


  //length of trajectory between plane to horizon and lower edge of camera
  //FoV
  dInFoV = l * sin(dFoVBelow) / ( dNy * sin(dFoVBelow) +  dNz * cos(dFoVBelow)  );

  //length of trajectory between plane to horizon and earth surface
  Double_t p = 2 * ( REarth*dNz - (v-l)*dNy );
  Double_t q = (v-l)*(v-l);

  if(q-p*p/4>=0) //trajectory does not intersect with Earth
      return 0;
        
  Double_t i1 = p/2. - sqrt(p*p/4.-q);
  Double_t i2 = p/2. + sqrt(p*p/4.-q);
            
  if(i1<0 || i2<0)
    cout<<"PDecay: i1 or i2 less than 0: "<<i1<<"  "<<i2<<endl;

  dED= i1<i2 ? i1 : i2;

  if( (dInFoV>dED && dInFoV>0 ) || dInFoV<0)
          dInFoV=dED;
 
  if(l>v) //if shower emerges beyond horizon
     dInFoV = 0;

  //maximum distance from where tau emerges to that plane
  //modify this to take into account extra length if shower emerges before the
  //horizon
  //add extra length if trajectory crosses plane between telescope and horizon
  Double_t dem = sqrt(l*l-d*d) + dInFoV;


  //fix this to use elevation measured when shower emerges from ground
  fPE->FixParameter(1,elevation ); //Shower elevation in rad

  //length of shower in camera plane
  Double_t dShwrLgth = 0.304 * log(Etau*0.5/0.088)/log(2); //0.304km radiation length at see level, 0.088GeV critical energy of electrons in air,only 0.5 of the energy goes into the electromagnetic shower
  //cout<<"Etau "<<Etau<<" size of shower in km: "<<dShwrLgth<<endl;
  //cout<<" dMaxDist: "<<dMaxDist<<endl;

  //minimal distance from tip of shower to the plane that is normal to trajectory and goes through the origin (where the detector is located), constrained by maximum angle a sufficient Cherenkov light reaches the detector. values below are 90-a and a in the sines. Need to have functions of tau energy and distance for a 
  //dd is first used as the distance to the start of the shower not hte tip of
  //the shower
  Double_t dd = dShwrLgth+5; //we need to at least have the shower develop. note we neglect the decay length here, which does not really matter. That is taken care of later.
 
  if(dd>dem)
     dd = dShwrLgth;
 
  while(dd<dem)
   {
  //move along trajectory and find spot where MaxCherenkovAngle condition is
  //fullfilled
  Double_t dDistanceToWhereTauStarts = sqrt(d*d+dd*dd);
  fPE->FixParameter(0,dDistanceToWhereTauStarts); //Distance to where the tau starts shower
   //get new azimuth
   Double_t g = (dem-dd) * cos(elevation);
   Double_t az = 0;
   if(g>0 && azimuth>0)
     {
        Double_t xi = ( l - g * cos(azimuth)) / (g * sin(azimuth));
        az =  pi*0.5 + azimuth - atan(xi);
        
        if (dem-dd-dInFoV < 0 ) //if we are below the plane
          az = azimuth;
     }
    else
      cout<<"in PDecay, azimuth is zero or g is smaller zero"<<endl;
   //get PE for new azimuth
     //double az = asin(d/dDistanceToWhereTauStarts); //azimuth for that distance
   if(fPE->Eval(az)*Etau*0.5<dMinimumNumberPhotoelectrons)
       dd++;
     else
      break;

   }

  //cout<<"maximum available length for decay and shower to happen (dem:) "<<dem<<" distance between trajectory and telescope d: "<<d<<" dd: "<<dd<<endl;
  //tip of the shower has to be inside the atmosphere. Check if that is the case. if not adjust dd
  //cout<<"need to takeaway dd so we can see all the Cherenkov light: "<<dd<<endl; 
  //if dem is less then dd, which means Cherenkov light will not hit the
  //telescope. return 0
  if(dd>dem) // the shower cannot be seen by the telescope because the cherenkov cone does not illuminate the telescope anywhere along the track
   return 0;

  //dd below is used as the distence between the plane perp. to the trajectory
  //and the tip of the shower so lets subtract the length of the shower and the
  //5 km again
  if(dd<dShwrLgth+5)
    dd -= dShwrLgth;
  else
    dd -= (dShwrLgth+5);
  
  
  //maximum length of trajectory above horizon befor track leaves atmosphere (dST above ground)
  Double_t phi = elevation + asin( REarth/sqrt(REarth*REarth+(l-v)*(l-v)) );
  
  Double_t alpha = asin( sin(phi) * sqrt(REarth*REarth+(l-v)*(l-v)) / (REarth+dST)  );

  Double_t gamma = pi - alpha - phi;

  Double_t dMaxDist = (REarth+dST) *  sin(gamma)/sin(phi);

 //trim the path length to be fully inside the atmosphere, should always be the
 //case
 if(l>v)//if the shower is not seen over the entire track in the atmosphere and l>v. Reset the maximum possible track length to the portion that can be seen
  {
    dMaxDist = dMaxDist>dem-dd ? dem-dd : dMaxDist;
  }
 else//do the same if the shower emerges l<v from telescope v is where the tangent touches earth.
  {
    dMaxDist = dST/sin(elevation) + dInFoV;
    dMaxDist = dMaxDist>dem-dd ? dem-dd : dMaxDist;
  }

 // Double_t delta = acos(-dNy);
 // cout<<"Delta: "<<delta*180/pi<<" angle from tip of shower to telescope: "<<asin(y*sin(delta)/sqrt(y*y+dMaxDist*dMaxDist-2*dMaxDist*y*cos(delta)))*180/pi  <<" shortest distance so Cherenkov light goes to camera: "<<dd<<endl;


  Double_t dTermInSquareRoot = cos(elevation)*cos(azimuth)*cos(elevation)*cos(azimuth)
                               + sin(elevation)*sin(elevation)/tanFoV/tanFoV
                               - cos(elevation)*cos(elevation);

  if(dTermInSquareRoot>0) //if it is negative the shower is contained in the FoV anywhere along the track
    {
        //maximum trajectory length before the tip of the shower is not contained in the camera anymore
       Double_t dMaxDisttoSatisfyFovReq = l / ( cos(elevation)*cos(azimuth) + sqrt( dTermInSquareRoot )) + dInFoV;

       if(dMaxDist>dMaxDisttoSatisfyFovReq+dInFoV && dMaxDisttoSatisfyFovReq>0) //correct max. trajectory length
           dMaxDist = dMaxDisttoSatisfyFovReq+dInFoV;
   //    cout<<" dMaxDisttoSatisfyFovReq: "<<dMaxDisttoSatisfyFovReq<<endl;
    }

 
 
  
  //ok have taken all requirements into account. Lets see if the trajectory
  //length allows for a full shower development.
  //If not we quit.
  if(dMaxDist<dShwrLgth)
   return 0;
  //make sure the shower does not develop past the point where more than 90% of the
  //taus have decayed
  Double_t DecayLength = Etau * c * DecayTime / Mtau;
  Double_t d90PctDecayLength = -log(0.1)*DecayLength;
  //cout<<"90% of taus decayed after: "<<d90PctDecayLength<<endl;

  if(d90PctDecayLength+dShwrLgth<dMaxDist)
    dMaxDist = d90PctDecayLength+dShwrLgth;

  //now lets check if the size of the shower is fullfilling the minimum length
  //requirement

  //calculating length of shower in the camera assuming the shower happens late
  //and develops up to the maximumg possible point along the trajectory
  //add that distance is within maximum distance like we do for the Fluorescence
  //part
  Double_t n = dMaxDist - dShwrLgth -dInFoV;
  Double_t m = dMaxDist;

  Double_t A = sqrt( n* dNx * n* dNx + (n*dNy+y) * (n*dNy+y) + n*dNz * n*dNz  );
  Double_t B = sqrt( m* dNx * m* dNx + (m*dNy+y) * (m*dNy+y) + m*dNz * m*dNz  );

  Double_t costheta = n*dNx * m*dNx + (n*dNy+y)*(m*dNy+y) + m*dNz * n*dNz;
  costheta = costheta / (A*B);
  Double_t dLength = acos(costheta)*180/pi;
  //cout<<"size of shower in degrees: "<<dLength<<" cos of angle:  "<<costheta<<endl;


  if(dLength<dMinLength) //shower image is too short
   return 0;


  //Ok finally we are there. Lets decay the tau in the remaining distance we
  //have
 
  //Probabilty that tau survives if it is not in the field of view
  //below works for tau emerging beyond horizon and before horizon 
  Double_t ProbTauDecay = exp(-(dED-dInFoV)/DecayLength); //Tau has to to survive before it becomes visible to the detector

  
  ProbTauDecay *= 1-exp(-(dMaxDist-dShwrLgth)/DecayLength); // then it has to decay before it is out of the FoV

  ProbTauDecay*=0.8;//only 80% of taus make a shower


//if(azimuth*180/pi>90 && azimuth*180/pi<8.1)
//cout<<"Etau: "<<Etau<<" el: "<<elevation*180/pi<<" az: "<<azimuth*180/pi<<" alpha "<<alpha<<" beta "<<180/pi*(pi - asin(sinbeta))<<" Prob: "<<ProbTauDecay<<" l "<<l<<" v "<<v<<" below horizon: "<<dBH<<" MaxDistOfTrack: "<<dMaxDist<<" minimal distance between trajectory and telescope "<<d<<endl;
return ProbTauDecay;
}


void PlotEmergenceProbability()
{

  TH1D *hTau = new TH1D("hTauS","",70,4,11);
  //hTau->SetMaximum(1);
  hTau->GetXaxis()->SetTitle("energy [GeV]");
  hTau->GetYaxis()->SetTitle("F_tau/F_nu");
  TAxis *axis = hTau->GetXaxis();
  int bins = axis->GetNbins();
  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];
  for (int i = 0; i <= bins; i++) {
     new_bins[i] = TMath::Power(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  TMultiGraph *mg = new TMultiGraph();
  TLegend *leg = new TLegend(0.7,0.4,0.89,0.88,"neutrino energy");


  double Enulog = 11;
  double Enuminlog = 5.9;
  double Enusteplog = 0.5;
  int s=0;
  while(Enulog>Enuminlog)
    {
      //cout<<Enulog<<endl;
      double Enu = pow(10,Enulog);
      TGraph *grProb = new TGraph(); 
      grProb->SetMarkerStyle(20+s);
      TString title;
      title.Form("%0.0e GeV",Enu);
      leg->AddEntry(grProb,title.Data(),"p"); 
      mg->Add(grProb,"lp");
      s++;

      // loop over target thickness
      double d = 0; //in 10^dmin km
      double dmax = 4;
      double dstep = 0.2;
      int p=0;
      while(d<dmax)
       {          
        double targetthickness = pow(10,d); 
        double dSumProb = 0;
        hTau->Reset();
        for(int i=1;i<hTau->GetNbinsX();i++)
           {

               Double_t Etau = hTau->GetBinCenter(i+1); 
               if(hTau->GetBinLowEdge(i+2)<=Enu)
                 {
                    Double_t P = PEtau(targetthickness,Etau,Enu,hTau);
                    //cout<<i+1<<" "<<targetthickness<<" . "<<Etau<<"  "<<Enu<<" P "<<P<<endl;
                  P *= (hTau->GetBinLowEdge(i+1)-hTau->GetBinLowEdge(i));

                    hTau->Fill(Etau,P); 
                    hTau->SetBinError(i,0);
                    dSumProb+=P;
                 }

            }//got the energy spectrum of the taus for this azimuth and elevation

         grProb->SetPoint(p,targetthickness,dSumProb);
         p++;
         d+=dstep;            
        }    
      Enulog-=Enusteplog;
     }
  TCanvas *cProbOfEmergence = new TCanvas("cProbOfEmergence","Probability of emergence",750,500);
  cProbOfEmergence->Draw();
  cProbOfEmergence->SetLogx();
  cProbOfEmergence->SetLogy();
  mg->Draw("a");
  mg->GetXaxis()->SetTitle("target thickness [km]");
  mg->GetXaxis()->SetTitleSize(0.045);
  mg->GetXaxis()->SetTitleOffset(1.1);
  mg->GetXaxis()->SetLabelSize(0.045);
  mg->GetYaxis()->SetTitle("probability of #tau emergence");
  mg->GetYaxis()->SetTitleOffset(1.0);
  mg->GetYaxis()->SetTitleSize(0.045);
  mg->GetYaxis()->SetLabelSize(0.045);
  mg->GetYaxis()->SetRangeUser(1e-6,1);
 // mg->GetXaxis()->SetRangeUser(1,5e3);
  leg->Draw();
  //TF1 *fdeg=new TF1("fdeg","90-180/3.1415*TMath::ASin(0.5*x/6371)",0,1e4);
  cout<<mg->GetXaxis()->GetXmin()<<endl;
  cout<<180/3.1415*TMath::ASin(0.5*mg->GetXaxis()->GetXmin()/6371)<<"  "<<180/3.1415*TMath::ASin(0.5*mg->GetXaxis()->GetXmax()/6371)<<endl;
  TF1 *fdeg=new TF1("fdeg","x",180/3.1415*TMath::ASin(0.5*mg->GetXaxis()->GetXmin()/6371),180/3.1415*TMath::ASin(0.5*mg->GetXaxis()->GetXmax()/6371));
  fdeg->Eval(0);
  TGaxis *degaxis = new TGaxis(mg->GetXaxis()->GetXmin(),1,mg->GetXaxis()->GetXmax(),1,"fdeg",510,"-G");
  degaxis->SetTitle("elevation angle [degrees]");
  degaxis->SetTitleFont(42);
  degaxis->SetLabelFont(42);
  degaxis->Draw();
}


void GetTauDistribution(TH1D *hTauSpec, Double_t d, Double_t Enumin = 1e9, Double_t Enumax = 3.16e9)
{
   hTauSpec->Reset(); //get the energy spectrum of taus coming out of the Earth, starting with nus in the range expmin expMax
   //int nEnuSteps = 20;
   //Double_t DeltaEnu = (Enumax - Enumin)/nEnuSteps;
   Double_t DeltaEnu = 0.1; //logscale
   int nEnuSteps = (0.0001+log10(Enumax) - log10(Enumin))/DeltaEnu; //the 0.0001 is due to small uncertainties making sure we get the right number of steps
   Double_t Normalization = (nuIndex-1)/(pow(Enumin,1-nuIndex)-pow(Enumax,1-nuIndex)); 
   vector<double> Enu;
   vector<double> EnuWeight;
   for(int i=0;i<nEnuSteps;i++)
      {
        //Enu.push_back(Enumin+i*DeltaEnu+0.5*DeltaEnu);
        //NO Enu.push_back(pow(10,(log10(Enumin+i*DeltaEnu)+log10(Enumin+(i+1)*DeltaEnu))*0.5));
        Enu.push_back(pow(10,log10(Enumin)+i*DeltaEnu+DeltaEnu*0.5));
        //EnuWeight.push_back(pow(Enumin+i*DeltaEnu+0.5*DeltaEnu,-nuIndex)*Normalization*DeltaEnu);
        //NO EnuWeight.push_back(pow(Enu[Enu.size()-1],-nuIndex)*Normalization*DeltaEnu);
        EnuWeight.push_back(pow(Enu[Enu.size()-1],-nuIndex)*Normalization*
                       (pow(10,log10(Enumin)+(i+1)*DeltaEnu)-pow(10,log10(Enumin)+i*DeltaEnu)));
      }

   for(int i=1;i<=hTauSpec->GetNbinsX();i++)
       {

          Double_t Etau = hTauSpec->GetBinCenter(i); 

         //loop over all Enu in this energy bin to calculate the sensitivity
          if(bMonoNu) //if we want to simulate only monoenergetic neutrinos (only good for acceptance calculations
            {
               //NO removed because that is transferred into PEtau 
               //if(Etau<=0.8*Enu) //0.8 because I need to account for part of the energy being transferred to the also produced neutrinos
                 {
                    Double_t Enu = pow(10,(log10(Enumin)+log10(Enumax))*0.5);
                    Double_t P = PEtau(d,Etau,Enu,hTauSpec);
                    hTauSpec->Fill(Etau,P);
                 } 
            }
          else
            {
              //while(Enu<Enumax)//loop over nu energy bin
              for(int n = 0; n<nEnuSteps;n++)
                {
                  //NO removed because that is transferred into PEtau 
                  //if(Etau<=0.8*Enu) //0.8 because I need to account for part of the energy being transferred to the also produced neutrinos
                    {
                      Double_t P = 0;
                      //if(Enu==Enumin)  //mod
                      P = PEtau(d,Etau,Enu[n],hTauSpec);
                      //if(Enu==Enumin) //mod
                      //cout<<P<<"  "<<d<<"  "<<Etau<<" "<<Enu[n]<<"  "<<EnuWeight[n]<<endl; //mod
                      //P *= DeltaEnu/(Enumax-Enumin); 
                      P *= EnuWeight[n];

                      //do not know how the below is calculated
                      //P *= ( pow(Enu,1-nuIndex) - pow(Enu+DeltaEnu,1-nuIndex) ) 
                      // / ( pow(Enumin,1-nuIndex) - pow(Enumax,1-nuIndex) );//multiplying in the with
                      hTauSpec->SetBinContent(i,hTauSpec->GetBinContent(i)+P); 
                      //hTauSpec->Fill(Etau,P); 
                    }
                //Enu += DeltaEnu;  
               }
           }    
      }//got the energy spectrum of the taus for this azimuth and elevation

 /*
               for(int i=0;i<hTauSpec->GetNbinsX();i++)
                  {
cout<<i+1<<"  "<<hTauSpec->GetBinCenter(i+1)  <<" taus cont: "<<hTauSpec->GetBinContent(i+1)<<endl;
}
*/
      //finding lowest and highest bin that is non-zero
      //bin with maximum content
      Int_t iBinMax = hTauSpec->GetMaximumBin();
      Double_t dContInMax = hTauSpec->GetBinContent(iBinMax);
      //   cout<<"iBinMax: "<<iBinMax<<" maxcontent: "<<dContInMax<<endl;       
   
      //Double_t dHighE = hTauSpec->GetBinLowEdge(iBinMax+1);
      //finding lowest bin
      Int_t iBinMin = iBinMax;
      while(hTauSpec->GetBinContent(iBinMin)>1e-4*dContInMax)
        {
           iBinMin--;   
        }
      // cout<<"iBinMin: "<<iBinMin<<endl;

     //finding highest bin
      while(hTauSpec->GetBinContent(iBinMax+1)>0)
         iBinMax++;
 
      //cout<<"iBinMax: "<<iBinMax<<endl;       


       //Double_t dLowE = hTauSpec->GetBinLowEdge(iBinMin);
       //Double_t dDeltaE = dHighE-dLowE;

       for(int i=1;i<hTauSpec->GetNbinsX();i++) //weight each bin with energy scale need to do that because of logarithmic energy scale
          {
            Double_t dContent = hTauSpec->GetBinContent(i);
            dContent *= (hTauSpec->GetBinLowEdge(i+1)-hTauSpec->GetBinLowEdge(i));
            if(i<iBinMin || i>iBinMax)
                 dContent = 0;
            hTauSpec->SetBinContent(i,dContent);
            hTauSpec->SetBinError(i,0);
          }

}


//Calculates the acceptance by looping over distance, azimuth and elevation for
//a given energy bin
//returns the integral acceptance
//and returns a graph that holds the acceptance for each step in distance in the
//loop
Double_t CalculateAcceptance(Double_t dMinEnu, Double_t dMaxEnu,TGraph *grDiffAcceptance,TH1D *hTau, Int_t tIndex)
{
  if (tIndex==5){
    fPE = new TF1("fPE",myPEfunction,0,40,2);
  }else{
    cout<<"/////" << tIndex <<endl;
    fPE = new TF1("fPE",myPEfunction_trigwindow,0,40,3);
    fPE->FixParameter(2,tIndex);
  }
  
  dMinEnu = pow(10,dMinEnu);
  dMaxEnu = pow(10,dMaxEnu);
  
  //the hTau passed to the function is ignored. Create a new one below
  //hTau = new TH1D("hTauNew","",dMaxEnu/(dMinEnu/10.),dMinEnu/10.,dMaxEnu);
  /*
  //make logarithmic axis
  TAxis *axis = hTauNew->GetXaxis();
  int bins = axis->GetNbins();
  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];
  for (int i = 0; i <= bins; i++) {
  new_bins[i] = TMath::Power(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  axis->Delete();
  */
  //area of cell
  Double_t dConversion=yDelta*2*pi; //multiply area of cell taking into account that we have a 360 degree FoV
  dConversion*=1e10; //from km2 to cm2
  //solid angle
  dConversion*=DeltaAngleAz/180.*pi*DeltaAngle/180.*pi; //multiply area of solidangle cell
  
  //time Do that in the sensitivity calculation. Acceptance is calculated
  //without the observing time
  //dConversion*=3*365*24*3600*0.20; //exposure time 3 years in seconds with 20% duty cycle
  
  
  dConversion*=2; //because we only calculate for azimuth angles 0 to azimuth max. There are also negative azimuth values due to symmetry of the problem 
  
  
  Double_t dIntegratedAcceptance=0;
  Int_t p = 0;
  Double_t y = yMin; //y distance from telescope where tau comes out of the ground;
  while(y<yMax) //loop over distance to telescope
    {
      //cout<<"Distance from Detector: "<<y<<endl;
      //calculate for given elevation the length of the trajectory through earth.
      Double_t dAcceptance=0.0;
      
      Double_t MaxAzimuth = dMaxCherenkovAzimuthAngle;
      if( bFluorescence || (bCombined && y<dMaxFluorescenceDistance) ) // so we can make full use of fluoresence events
	MaxAzimuth = 180;
      
      Double_t elevation=DeltaAngle*0.5;    
      while(elevation<MaxElevation) //loop over elevation
	{
	  Double_t dWeightForTriggeredAzimuth = sin(elevation/180.*pi)*y;
	  Double_t azimuth = DeltaAngleAz*0.5;
	  while(azimuth<MaxAzimuth) // loop over azimuth
	    {
	      
	      if(azimuth>MaxAzimuth && y>dMaxFluorescenceDistance)
		cout<<" Azimuth:" <<azimuth<<" should not be here "<<endl;
	      
	      Double_t dEarth = DistanceThroughEarth(y,elevation,azimuth);
	      //cout<<"   length of trajectory in Earth: "<<dEarth<<" km"<<endl;
	      
	      //cout<<dEarth<<"  "<<dMinEnu<<"  "<<dMaxEnu<<endl; 
	      GetTauDistribution(hTau,dEarth,dMinEnu,dMaxEnu);                
	      //cout<<"done GetTau"<<endl;
	      
	      Double_t dDeltaAcceptance=0;
	      
	      //Calculate probability that taus with E convert before they are 150 km away from detector when they make it out of the earth. 
	      //150km is for an angle of 5 degrees assuming 10 degree opening angle. If the angle is free the distance of closest approach depends on
	      // where the tau comes out of the earth and under what angle (vertical and horizontal. 
	      //That is probably dependend on the initial nu energy
	      Double_t dP = 0;
	      //     if(elevation>0.5 && elevation<1 && azimuth>9 && azimuth<9.2)
	      for(int i=0;i<hTau->GetNbinsX();i++)
		{
		  //cout<<dMinEnu<<"  "<<dMaxEnu<<endl;
		  //cout<<"taus cont: "<<hTau->GetBinContent(i+1)<<endl;
		  if(hTau->GetBinContent(i+1)>0)
		    {
		      Double_t dPFluorescence = 0.0;
		      Double_t dPCherenkov = 0.0;
		      //cout<<bFluorescence<<" "<<bCombined<<"  "<<y<<"<"<<dMaxFluorescenceDistance<<endl;
		      if( bFluorescence || (bCombined && y<dMaxFluorescenceDistance) )
			dPFluorescence = PDecayFluorescence(hTau->GetBinCenter(i+1),y,elevation,azimuth);
		      if( (!bFluorescence || bCombined) && azimuth<dMaxCherenkovAzimuthAngle  )
			dPCherenkov = PDecay(hTau->GetBinCenter(i+1),y,elevation,azimuth);
		      
		      if(bCombined)
			dP = dPFluorescence > dPCherenkov ? dPFluorescence : dPCherenkov;
		      else if(bFluorescence)
			dP = dPFluorescence;
		      else
			dP = dPCherenkov;
		      dDeltaAcceptance+=hTau->GetBinContent(i+1)*dP;
		      //dDeltaAcceptance+=hTau->GetBinContent(i+1); //use above
		      if(!bMonoNu)
			hTriggeredAzimuthAngles->Fill(azimuth,hTau->GetBinContent(i+1)*dP*dWeightForTriggeredAzimuth);
		    }
		  //if(hTau->GetBinContent(i+1)*dP>0 )
		  //cout<<hTau->GetBinCenter(i+1)<<"  "<<hTau->GetBinContent(i+1)<<" y:  "<<y<<"  el: "<<elevation<<" az: "<<azimuth<<" dp: "<<dP<<" dDeltaAccept: "<<dDeltaAcceptance<<" prod: "<<hTau->GetBinContent(i+1)*dP<<endl;
		}
	      if(dDeltaAcceptance<1e-10 && y>50) //won't get any more acceptance. The >60 is to make sure we do not miss fluorescence events whic can be seen from the back
		break;
	      
	      dAcceptance+=dDeltaAcceptance*sin(elevation/180.*pi); //projection of area cell to trajectory
	      //   cout<<"distance "<<y<<" prob"<<dDeltaAcceptance<<" elevation  "<<elevation<<endl;
	      azimuth+=DeltaAngleAz;
	      
	      //Add absorption in the atmosphere between shower and observer
	      //Go over target area and calculate acceptance angle for each dA. Integrate over energy spectrum of taus coming out of the earth at that point. multiplied with detection efficiency(absorption).
	    }//finished looping over all azimuth angles
	  elevation+=DeltaAngle;    
	  //cout<<azimuth<<"  "<<dAcceptance<<endl;
	  //multiply with dOmega  DeltaAngle*DeltaAngle
        }//finished looping over all elevation angles
      //dAcceptance*=yDelta*DeltaAngle/180*pi*y; //multiply area of cell
      //dAcceptance*=DeltaAngle/180*pi*DeltaAngle/180*pi; //multiply area of solidangle cell
      dAcceptance*=y; //multiply with area of cell (note that the yDelta*DeltaAngle/180*pi is included in dConversion
      dIntegratedAcceptance+=dAcceptance;
      
      grDiffAcceptance->SetPoint(p,y,dAcceptance*dConversion); 
      p++;
      cout<<"distance: "<<y<<" differ. acceptance "<<(dAcceptance*dConversion)<<" integr. acceptance.: "<<(dIntegratedAcceptance*dConversion)<<endl; 
      cout<<dAcceptance<<"  "<<dConversion<<endl;
      y+=yDelta;
      if(dAcceptance<1e-10 && y>50) //no sense to increase in distance if we can't see any showers now
	break;
    }//end looping over distances
  return (dIntegratedAcceptance*dConversion);
}

void CalculateAcceptanceVsLowerFoV(TH1D *hTau)
{

  bCombined = kFALSE;

  MaxElevation = 10; //elevation angle (determines path through Earth; 
  DeltaAngle = 0.05; //steps in azimuth and elevation 
  yMin = 0;
  yDelta = 5; 


  iMirrorSize = 2;
  dMinimumNumberPhotoelectrons = dThreshold[iMirrorSize]/dMirrorA[iMirrorSize]; 

  dMinLength = 0.3; //mimnimum length a shower has to have in the camera, in degrees. This is a conservative estimate because it assumes that the shower starts at a distance l from the detector, which is not necessarily tru for showers with shallow elevation angles.

  tanFoV = tan(10/180.*pi); //Field of view of telescope above the horizon
  

  Double_t dPreFactor = pow(10,(dMaxEnu+dMinEnu)) / (pow(10,dMaxEnu)-pow(10,dMinEnu));
  dPreFactor *= 3;//three neutrino flavours

  iConfig = 2; //telescope altitude


  TMultiGraph *mgDiffAcceptance = new TMultiGraph();
  TMultiGraph *mgAcceptance = new TMultiGraph();
  TLegend *legAcceptance = new TLegend(0.6,0.54,0.91,0.91,"FoV below horizon");
  legAcceptance->SetTextSize(0.05);

  for(int b = 0;b<2;b++)
  {
    bFluorescence=b;
    TGraph *grAcceptvsFoV = new TGraph();
    grAcceptvsFoV->SetLineWidth(3);
    grAcceptvsFoV->SetLineColor(kBlue+3);
    mgAcceptance->Add(grAcceptvsFoV);

     if(bFluorescence)
      {
        
        yMax = 100;
        grAcceptvsFoV->SetLineStyle(9);
      }
    else
      {
         yMax = 400;
       }

    for(int f = 0;f<6;f++)
     {
        dFoVBelow = (f+0.01)/180.*pi;

       TGraph *grDiffAcceptance = new TGraph();

       if(bFluorescence)
        {
           grDiffAcceptance->SetLineStyle(2);
           TString title;
           title.Form("%i#circ",(f));
           legAcceptance->AddEntry(grDiffAcceptance,title.Data(),"p"); 
        }


       grDiffAcceptance->SetMarkerStyle(marker[f]);
       grDiffAcceptance->SetMarkerSize(markerSize[f]);
       grDiffAcceptance->SetMarkerColor(iColors[f]);
       grDiffAcceptance->SetLineColor(iColors[f]);
       mgDiffAcceptance->Add(grDiffAcceptance,"lp");
    
       Double_t dAcceptance = CalculateAcceptance(dMinEnu,dMaxEnu,grDiffAcceptance,hTau,iTrigWin);
     //  cout<<"Sensitivity "<<dPreFactor/dAcceptance<<endl;
       grAcceptvsFoV->SetPoint(f,f*2+3,dAcceptance);
     }//looping over different FoV
   } //looping over Fluo and Cherenkov
  //Done with calculating the sensitivity

  TCanvas *cAcceptvsFoV = new TCanvas("cAcceptvsFoV","Acceptance vs. lower FoV",750,500);
  cAcceptvsFoV->Draw();
  mgAcceptance->Draw("alp");
  mgAcceptance->GetXaxis()->SetTitle("vertical field of view [degrees]");
  mgAcceptance->GetYaxis()->SetTitle("acceptance [cm^{2} sr]");
  mgAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgAcceptance->GetXaxis()->SetLabelSize(0.04);
  

  TCanvas *cDiffAcceptDistance = new TCanvas("cDiffAcceptDistance","Acceptance vs. Distance",750,500);
  cDiffAcceptDistance->Draw();
  mgDiffAcceptance->Draw("a");
  mgDiffAcceptance->GetXaxis()->SetTitle("distance of tau emergence from telescope [km]");
  mgDiffAcceptance->GetYaxis()->SetTitle("radial acceptance [cm^{2} sr]");
  mgDiffAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetLabelSize(0.04);
  legAcceptance->Draw();

  
}


void CalculateAcceptanceVsUpperFoV(TH1D *hTau)
{

  bCombined = kFALSE;

  MaxElevation = 10; //elevation angle (determines path through Earth; 
  DeltaAngle = 0.05; //steps in azimuth and elevation 
  yMin = 0;
  yDelta = 5; 


  iMirrorSize = 2;
  dMinimumNumberPhotoelectrons = dThreshold[iMirrorSize]/dMirrorA[iMirrorSize]; 

  dMinLength = 0.3; //mimnimum length a shower has to have in the camera, in degrees. This is a conservative estimate because it assumes that the shower starts at a distance l from the detector, which is not necessarily tru for showers with shallow elevation angles.

  Double_t dPreFactor = pow(10,(dMaxEnu+dMinEnu)) / (pow(10,dMaxEnu)-pow(10,dMinEnu));
  dPreFactor *= 3;//three neutrino flavours

  iConfig = 2; //telescope altitude

  dFoVBelow =   asin(REarth/(REarth+DetectorAltitude[iConfig]));
  //dFoVBelow =  3.0/180.*pi; 

  TMultiGraph *mgDiffAcceptance = new TMultiGraph();
  TMultiGraph *mgAcceptance = new TMultiGraph();
  TLegend *legAcceptance = new TLegend(0.6,0.54,0.91,0.91,"FoV above horizon");
  legAcceptance->SetTextSize(0.05);

  for(int b = 0;b<2;b++)
  {
    bFluorescence=b;
    TGraph *grAcceptvsFoV = new TGraph();
    grAcceptvsFoV->SetLineWidth(3);
    grAcceptvsFoV->SetLineColor(kBlue+3);
    mgAcceptance->Add(grAcceptvsFoV);

     if(bFluorescence)
      {
        
        yMax = 100;
        grAcceptvsFoV->SetLineStyle(9);
      }
    else
      {
         yMax = 400;
       }

    for(int f = 0;f<6;f++)
     {
       tanFoV = tan((f*0.5+0.01)/180.*pi);

       TGraph *grDiffAcceptance = new TGraph();

       if(bFluorescence)
        {
           grDiffAcceptance->SetLineStyle(2);
           TString title;
           title.Form("%0.1f#circ",(f*0.5));
           legAcceptance->AddEntry(grDiffAcceptance,title.Data(),"p"); 
        }


       grDiffAcceptance->SetMarkerStyle(marker[f]);
       grDiffAcceptance->SetMarkerSize(markerSize[f]);
       grDiffAcceptance->SetMarkerColor(iColors[f]);
       grDiffAcceptance->SetLineColor(iColors[f]);
       mgDiffAcceptance->Add(grDiffAcceptance,"lp");
    
       Double_t dAcceptance = CalculateAcceptance(dMinEnu,dMaxEnu,grDiffAcceptance,hTau,iTrigWin);
       cout<<"Sensitivity "<<dPreFactor/dAcceptance<<endl;
       grAcceptvsFoV->SetPoint(f,f*2+3,dAcceptance);
     }//looping over different FoV
   } //looping over Fluo and Cherenkov
  //Done with calculating the sensitivity

  TCanvas *cAcceptvsFoV = new TCanvas("cAcceptvsFoV","Acceptance vs. FoV",750,500);
  cAcceptvsFoV->Draw();
  mgAcceptance->Draw("alp");
  mgAcceptance->GetXaxis()->SetTitle("vertical field of view [degrees]");
  mgAcceptance->GetYaxis()->SetTitle("acceptance [cm^{2} sr]");
  mgAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgAcceptance->GetXaxis()->SetLabelSize(0.04);
  

  TCanvas *cDiffAcceptDistance = new TCanvas("cDiffAcceptDistance","Acceptance vs. Distance",750,500);
  cDiffAcceptDistance->Draw();
  mgDiffAcceptance->Draw("a");
  mgDiffAcceptance->GetXaxis()->SetTitle("distance of tau emergence from telescope [km]");
  mgDiffAcceptance->GetYaxis()->SetTitle("radial acceptance [cm^{2} sr]");
  mgDiffAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetLabelSize(0.04);
  legAcceptance->Draw();

  
}

void CalculateAcceptanceVsEnergy(TH1D *hTau)
{

  bCombined = kFALSE;

  MaxElevation = 10; //elevation angle (determines path through Earth; 
  DeltaAngle = 0.05; //steps in azimuth and elevation 
  yMin = 0;
  yDelta = 5; 

  iConfig = 2; //telescope altitude
  dFoVBelow = asin(REarth/(REarth+DetectorAltitude[iConfig]));
  tanFoV = tan(5./180.*pi);
  dMinLength = 0.3; //mimnimum length a shower has to have in the camera, in degrees. This is a conservative estimate because it assumes that the shower starts at a distance l from the detector, which is not necessarily tru for showers with shallow elevation angles.

  iMirrorSize = 2;
  dMinimumNumberPhotoelectrons = dThreshold[iMirrorSize]/dMirrorA[iMirrorSize]; 

  Double_t dLogEnergyStep = 0.2;
  Double_t dHalfEnergyBinWidth =1/2.; //in log
  Double_t logEmin = 10.0; //was 7




  TMultiGraph *mgDiffAcceptance = new TMultiGraph();
  TMultiGraph *mgAcceptance = new TMultiGraph();
  TLegend *legAcceptance = new TLegend(0.64,0.64,0.95,0.95,"minimum image length");
  legAcceptance->SetTextSize(0.05);

 for(int b = 0;b<2;b++)
  {
    bFluorescence=b;
    TGraph *grAcceptvsImageSize = new TGraph();
    grAcceptvsImageSize->SetLineWidth(3);
    grAcceptvsImageSize->SetLineColor(kBlue+3);
    mgAcceptance->Add(grAcceptvsImageSize);

     if(bFluorescence)
      {
        
        yMax = 100;
        grAcceptvsImageSize->SetLineStyle(9);
      }
    else
      {
         yMax = 400;
       }

    for(int f = 0;f<1;f++) 
      {
       dMinEnu=logEmin+f*dLogEnergyStep-dHalfEnergyBinWidth;
       dMaxEnu=logEmin+f*dLogEnergyStep+dHalfEnergyBinWidth;
       Double_t dPreFactor = pow(10,dMaxEnu+dMinEnu) / (pow(10,dMaxEnu)-pow(10,dMinEnu));
       dPreFactor *= 3;//three neutrino flavours
       

       TGraph *grDiffAcceptance = new TGraph();
       if(bFluorescence)
        {
           grDiffAcceptance->SetLineStyle(2);
           TString title;
           title.Form("%0.1f GeV",logEmin+f*dLogEnergyStep);
           legAcceptance->AddEntry(grDiffAcceptance,title.Data(),"p"); 
        }
       
       grDiffAcceptance->SetMarkerStyle(marker[f]);
       grDiffAcceptance->SetMarkerSize(markerSize[f]);
       grDiffAcceptance->SetMarkerColor(iColors[f]);
       grDiffAcceptance->SetLineColor(iColors[f]);
       grDiffAcceptance->SetLineWidth(2);
       mgDiffAcceptance->Add(grDiffAcceptance,"lp");
    
       Double_t dAcceptance = CalculateAcceptance(dMinEnu,dMaxEnu,grDiffAcceptance,hTau,iTrigWin);
       cout<<"Sensitivity "<<dPreFactor/dAcceptance<<endl;
       grAcceptvsImageSize->SetPoint(f,dMinLength,dAcceptance);
      } //looping over different energies
   } //looping over Fluo and Cherenkov

  TCanvas *cAcceptvsImageSize = new TCanvas("cAcceptvsImageSize","Acceptance vs. Image Size",750,500);
  cAcceptvsImageSize->Draw();
  mgAcceptance->Draw("alp");
  mgAcceptance->GetXaxis()->SetTitle("Energy [log10(GeV)]");
  mgAcceptance->GetYaxis()->SetTitle("acceptance [cm^{2} sr]");
  mgAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgAcceptance->GetXaxis()->SetLabelSize(0.04);
  

  TCanvas *cDiffAcceptDistance = new TCanvas("cDiffAcceptDistance","Acceptance vs. Distance",750,500);
  cDiffAcceptDistance->Draw();
  mgDiffAcceptance->Draw("alp");
  mgDiffAcceptance->GetXaxis()->SetTitle("distance of tau emergence from telescope [km]");
  mgDiffAcceptance->GetYaxis()->SetTitle("radial acceptance [cm^{2} sr]");
  mgDiffAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetLabelSize(0.04);
  legAcceptance->Draw();


}
void CalculateAcceptanceVsImageLength(TH1D *hTau)
{

  bCombined = kFALSE;

  MaxElevation = 10; //elevation angle (determines path through Earth; 
  DeltaAngle = 0.05; //steps in azimuth and elevation 
  yMin = 0;
  yDelta = 5; 

  iConfig = 2; //telescope altitude
  dFoVBelow = 3./180*pi; //asin(REarth/(REarth+DetectorAltitude[iConfig]));
  tanFoV = tan(2./180.*pi);

  iMirrorSize = 2;
  dMinimumNumberPhotoelectrons = dThreshold[iMirrorSize]/dMirrorA[iMirrorSize]; 



  Double_t dPreFactor = pow(10,dMaxEnu+dMinEnu) / (pow(10,dMaxEnu)-pow(10,dMinEnu));
  dPreFactor *= 3;//three neutrino flavours


  TMultiGraph *mgDiffAcceptance = new TMultiGraph();
  TMultiGraph *mgAcceptance = new TMultiGraph();
  TLegend *legAcceptance = new TLegend(0.64,0.64,0.95,0.95,"minimum image length");
  legAcceptance->SetTextSize(0.05);

  for(int b = 0;b<2;b++)
  {
    bFluorescence=b;
    TGraph *grAcceptvsImageSize = new TGraph();
    grAcceptvsImageSize->SetLineWidth(3);
    grAcceptvsImageSize->SetLineColor(kBlue+3);
    mgAcceptance->Add(grAcceptvsImageSize);

     if(bFluorescence)
      {
        
        yMax = 100;
        grAcceptvsImageSize->SetLineStyle(9);
      }
    else
      {
         yMax = 400;
       }

    for(int f = 0;f<5;f++) 
      {
        dMinLength = 0.2*f+0.1; //mimnimum length a shower has to have in the camera, in degrees. This is a conservative estimate because it assumes that the shower starts at a distance l from the detector, which is not necessarily tru for showers with shallow elevation angles.
       TGraph *grDiffAcceptance = new TGraph();
       if(bFluorescence)
        {
           grDiffAcceptance->SetLineStyle(2);
           TString title;
           title.Form("%0.1f#circ",dMinLength);
           legAcceptance->AddEntry(grDiffAcceptance,title.Data(),"p"); 
        }
       
       grDiffAcceptance->SetMarkerStyle(marker[f]);
       grDiffAcceptance->SetMarkerSize(markerSize[f]);
       grDiffAcceptance->SetMarkerColor(iColors[f]);
       grDiffAcceptance->SetLineColor(iColors[f]);
       grDiffAcceptance->SetLineWidth(2);
       mgDiffAcceptance->Add(grDiffAcceptance,"lp");
    
       Double_t dAcceptance = CalculateAcceptance(dMinEnu,dMaxEnu,grDiffAcceptance,hTau,iTrigWin);
       cout<<"Sensitivity "<<dPreFactor/dAcceptance<<endl;
       grAcceptvsImageSize->SetPoint(f,dMinLength,dAcceptance);
      } //looping over different minimum image size
   } //looping over Fluo and Cherenkov

  TCanvas *cAcceptvsImageSize = new TCanvas("cAcceptvsImageSize","Acceptance vs. Image Size",750,500);
  cAcceptvsImageSize->Draw();
  mgAcceptance->Draw("alp");
  mgAcceptance->GetXaxis()->SetTitle("minimum image length [degrees]");
  mgAcceptance->GetYaxis()->SetTitle("acceptance [cm^{2} sr]");
  mgAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgAcceptance->GetXaxis()->SetLabelSize(0.04);
  

  TCanvas *cDiffAcceptDistance = new TCanvas("cDiffAcceptDistance","Acceptance vs. Distance",750,500);
  cDiffAcceptDistance->Draw();
  mgDiffAcceptance->Draw("alp");
  mgDiffAcceptance->GetXaxis()->SetTitle("distance of tau emergence from telescope [km]");
  mgDiffAcceptance->GetYaxis()->SetTitle("radial acceptance [cm^{2} sr]");
  mgDiffAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetLabelSize(0.04);
  legAcceptance->Draw();


}

void CalculateAcceptanceVsTelescopeHeight(TH1D *hTau)
{

  bCombined = kFALSE;

  MaxElevation = 10; //elevation angle (determines path through Earth; 
  DeltaAngle = 0.05; //steps in azimuth and elevation 
  yMin = 0;
  yDelta = 5; 

  tanFoV = tan(10./180.*pi);
  dFoVBelow = asin(REarth/(REarth+DetectorAltitude[3]));

  iMirrorSize = 2;
  dMinimumNumberPhotoelectrons = dThreshold[iMirrorSize]/dMirrorA[iMirrorSize]; 

  dMinLength = 0.3; //mimnimum length a shower has to have in the camera, in degrees. This is a conservative estimate because it assumes that the shower starts at a distance l from the detector, which is not necessarily tru for showers with shallow elevation angles.

  Double_t dPreFactor = pow(10,(dMaxEnu+dMinEnu)) / (pow(10,dMaxEnu)-pow(10,dMinEnu));
  dPreFactor *= 3;//three neutrino flavours


  TMultiGraph *mgDiffAcceptance = new TMultiGraph();
  TMultiGraph *mgAcceptance = new TMultiGraph();
  TLegend *legAcceptance = new TLegend(0.6,0.54,0.91,0.91,"telescope height");
  legAcceptance->SetTextSize(0.05);

  for(int b = 0;b<2;b++)
  {
    bFluorescence=b;
    TGraph *grAcceptvsHeight = new TGraph();
    grAcceptvsHeight->SetLineWidth(3);
    grAcceptvsHeight->SetLineColor(kBlue+3);
    mgAcceptance->Add(grAcceptvsHeight);

     if(bFluorescence)
      {
        
        yMax = 100;
        grAcceptvsHeight->SetLineStyle(9);
      }
    else
      {
         yMax = 400;
       }

    for(unsigned int f = 0;f<sizeof(DetectorAltitude)/sizeof(*DetectorAltitude);f++)
     {
      iConfig = f; //telescope altitude
      TGraph *grDiffAcceptance = new TGraph();
         if(bFluorescence)
          {
            grDiffAcceptance->SetLineStyle(2);
            TString title;
            title.Form("%i km",f);
            legAcceptance->AddEntry(grDiffAcceptance,title.Data(),"p"); 
          }

      grDiffAcceptance->SetMarkerStyle(marker[f]);
      grDiffAcceptance->SetMarkerSize(markerSize[f]);
      grDiffAcceptance->SetMarkerColor(iColors[f]);
      grDiffAcceptance->SetLineColor(iColors[f]);
      mgDiffAcceptance->Add(grDiffAcceptance,"lp");
    
      Double_t dAcceptance = CalculateAcceptance(dMinEnu,dMaxEnu,grDiffAcceptance,hTau,iTrigWin);
      cout<<"Sensitivity "<<dPreFactor/dAcceptance<<endl;
      grAcceptvsHeight->SetPoint(f,DetectorAltitude[f],dAcceptance);
     }//looping over different FoV
   } //looping over Fluo and Cherenkov
  //Done with calculating the sensitivity

  TCanvas *cAcceptvsHeight = new TCanvas("cAcceptvsHeight","Acceptance vs. Telescope Height",750,500);
  cAcceptvsHeight->Draw();
  mgAcceptance->Draw("alp");
  mgAcceptance->GetXaxis()->SetTitle("telescope height [km]");
  mgAcceptance->GetYaxis()->SetTitle("acceptance [cm^{2} sr]");
  mgAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgAcceptance->GetXaxis()->SetLabelSize(0.04);
  

  TCanvas *cDiffAcceptDistance = new TCanvas("cDiffAcceptDistance","Acceptance vs. Distance",750,500);
  cDiffAcceptDistance->Draw();
  mgDiffAcceptance->Draw("a");
  mgDiffAcceptance->GetXaxis()->SetTitle("distance of tau emergence from telescope [km]");
  mgDiffAcceptance->GetYaxis()->SetTitle("radial acceptance [cm^{2} sr]");
  mgDiffAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetLabelSize(0.04);
  legAcceptance->Draw();


}

void CalculateAcceptanceVsTriggerWindow(TH1D *hTau)
{

  bCombined = kFALSE;

  MaxElevation = 10; //elevation angle (determines path through Earth; 
  DeltaAngle = 0.05; //steps in azimuth and elevation 
  yMin = 0;
  yDelta = 5; 

  tanFoV = tan(10./180.*pi);
  dFoVBelow = asin(REarth/(REarth+DetectorAltitude[3]));

  iMirrorSize = 2;
  dMinimumNumberPhotoelectrons = dThreshold[iMirrorSize]/dMirrorA[iMirrorSize]; 

  dMinLength = 0.3; //mimnimum length a shower has to have in the camera, in degrees. This is a conservative estimate because it assumes that the shower starts at a distance l from the detector, which is not necessarily tru for showers with shallow elevation angles.

  Double_t dPreFactor = pow(10,(dMaxEnu+dMinEnu)) / (pow(10,dMaxEnu)-pow(10,dMinEnu));
  dPreFactor *= 3;//three neutrino flavours


  TMultiGraph *mgDiffAcceptance = new TMultiGraph();
  TMultiGraph *mgAcceptance = new TMultiGraph();
  TLegend *legAcceptance = new TLegend(0.6,0.54,0.91,0.91,"trigger window");
  legAcceptance->SetTextSize(0.05);

  for(int b = 0;b<2;b++)
  {
    bFluorescence=b;
    TGraph *grAcceptvsTrigWindow = new TGraph();
    grAcceptvsTrigWindow->SetLineWidth(3);
    grAcceptvsTrigWindow->SetLineColor(kBlue+3);
    mgAcceptance->Add(grAcceptvsTrigWindow);

     if(bFluorescence)
      {
        
        yMax = 100;
        grAcceptvsTrigWindow->SetLineStyle(9);
      }
    else
      {
         yMax = 400;
       }

    for(unsigned int f = 0;f<6;f++)
     {
      iTrigWin = f; //telescope trigger window
      TGraph *grDiffAcceptance = new TGraph();
         if(bFluorescence)
          {
            grDiffAcceptance->SetLineStyle(2);
            TString title;
	    if (f == 5 ){
	      title.Form("all time");
	    }else{
	      title.Form("%f",TriggerWindow[f]);
	    }
            legAcceptance->AddEntry(grDiffAcceptance,title.Data(),"p"); 
          }

      grDiffAcceptance->SetMarkerStyle(marker[f]);
      grDiffAcceptance->SetMarkerSize(markerSize[f]);
      grDiffAcceptance->SetMarkerColor(iColors[f]);
      grDiffAcceptance->SetLineColor(iColors[f]);
      mgDiffAcceptance->Add(grDiffAcceptance,"lp");
    
      Double_t dAcceptance = CalculateAcceptance(dMinEnu,dMaxEnu,grDiffAcceptance,hTau,f);
      cout<<"Sensitivity "<<dPreFactor/dAcceptance<<endl;
      if (f < 3)
	grAcceptvsTrigWindow->SetPoint(f,TriggerWindow[f],dAcceptance);
     }//looping over different trigger window widths
   } //looping over Fluo and Cherenkov
  //Done with calculating the sensitivity

  TCanvas *cAcceptvsTrigWindow = new TCanvas("cAcceptvsTrigWindow","Acceptance vs. Trigger Window",750,500);
  cAcceptvsTrigWindow->Draw();
  mgAcceptance->Draw("alp");
  mgAcceptance->GetXaxis()->SetTitle("trigger window [ns]");
  mgAcceptance->GetYaxis()->SetTitle("acceptance [cm^{2} sr]");
  mgAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgAcceptance->GetXaxis()->SetLabelSize(0.04);
  

  TCanvas *cDiffAcceptDistance = new TCanvas("cDiffAcceptDistance","Acceptance vs. Distance",750,500);
  cDiffAcceptDistance->Draw();
  mgDiffAcceptance->Draw("a");
  mgDiffAcceptance->GetXaxis()->SetTitle("distance of tau emergence from telescope [km]");
  mgDiffAcceptance->GetYaxis()->SetTitle("radial acceptance [cm^{2} sr]");
  mgDiffAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetLabelSize(0.04);
  legAcceptance->Draw();
}

void CalculateAcceptanceVsThreshold(TH1D *hTau)
{

  bCombined = kFALSE;

  MaxElevation = 10; //elevation angle (determines path through Earth; 
  DeltaAngle = 0.05; //steps in azimuth and elevation 
  yMin = 0;
  yDelta = 5; 

  tanFoV = tan(2./180.*pi);

  iConfig = 2; //telescope altitude
  dFoVBelow = 3./180*pi; 
 // dFoVBelow = asin(REarth/(REarth+DetectorAltitude[iConfig]));

  dMinLength = 0.3; //mimnimum length a shower has to have in the camera, in degrees. This is a conservative estimate because it assumes that the shower starts at a distance l from the detector, which is not necessarily tru for showers with shallow elevation angles.

  Double_t dPreFactor = pow(10,(dMaxEnu+dMinEnu)) / (pow(10,dMaxEnu)-pow(10,dMinEnu));
  dPreFactor *= 3;//three neutrino flavours


  TMultiGraph *mgDiffAcceptance = new TMultiGraph();
  TMultiGraph *mgAcceptance = new TMultiGraph();
  TLegend *legAcceptance = new TLegend(0.6,0.54,0.91,0.91,"mirror size");
  legAcceptance->SetTextSize(0.05);

  for(int b = 0;b<2;b++)
  {
    bFluorescence=b;
    TGraph *grAcceptvsMirror = new TGraph();
    grAcceptvsMirror->SetLineWidth(3);
    grAcceptvsMirror->SetLineColor(kBlue+3);
    mgAcceptance->Add(grAcceptvsMirror);

     if(bFluorescence)
      {
        
        yMax = 100;
        grAcceptvsMirror->SetLineStyle(9);
      }
    else
      {
         yMax = 400;
       }

  for(unsigned int f = 0;f<4;f++)
   {
     iMirrorSize = f;
     dMinimumNumberPhotoelectrons = dThreshold[iMirrorSize]/dMirrorA[iMirrorSize]; 

     TGraph *grDiffAcceptance = new TGraph();

     if(bFluorescence)
        {
           grDiffAcceptance->SetLineStyle(2);
           TString title;
           title.Form("%0.0f m^{2}",dMirrorA[f]);
           legAcceptance->AddEntry(grDiffAcceptance,title.Data(),"p"); 
        }

     
    grDiffAcceptance->SetMarkerStyle(marker[f]);
    grDiffAcceptance->SetMarkerSize(markerSize[f]);
    grDiffAcceptance->SetMarkerColor(iColors[f]);
    grDiffAcceptance->SetLineColor(iColors[f]);
    mgDiffAcceptance->Add(grDiffAcceptance,"lp");
    
    Double_t dAcceptance = CalculateAcceptance(dMinEnu,dMaxEnu,grDiffAcceptance,hTau,iTrigWin);
    cout<<"Sensitivity "<<dPreFactor/dAcceptance<<endl;
    grAcceptvsMirror->SetPoint(f,dMirrorA[iMirrorSize],dAcceptance);
   }//looping over different FoV
  } //looping over Fluo and Cherenkov
  //Done with calculating the sensitivity

  TCanvas *cAcceptvsMirrorArea = new TCanvas("cAcceptvsImageSize","Acceptance vs. Mirror Area",750,500);
  cAcceptvsMirrorArea->Draw();
  mgAcceptance->Draw("alp");
  mgAcceptance->GetXaxis()->SetTitle("mirror area [m^{2}]");
  mgAcceptance->GetYaxis()->SetTitle("acceptance [cm^{2} sr]");
  mgAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgAcceptance->GetXaxis()->SetLabelSize(0.04);
  

  TCanvas *cDiffAcceptDistance = new TCanvas("cDiffAcceptDistance","Acceptance vs. Distance",750,500);
  cDiffAcceptDistance->Draw();
  mgDiffAcceptance->Draw("a");
  mgDiffAcceptance->GetXaxis()->SetTitle("distance of tau emergence from telescope [km]");
  mgDiffAcceptance->GetYaxis()->SetTitle("radial acceptance [cm^{2} sr]");
  mgDiffAcceptance->GetYaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetYaxis()->SetLabelSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetTitleSize(0.04);
  mgDiffAcceptance->GetXaxis()->SetLabelSize(0.04);
  legAcceptance->Draw();


}
////////////////////////////////////////////////////////////////////
//
//
// Calculate Integral Sensitiivty
//
//
void CalculateIntegralSensitivity(TH1D *hTau)
{

    Double_t dLogEnergyStep = 0.2;
    Double_t logEmin = 6;
    Double_t logEmax = 10;

    if(bFluorescence)
       {
          yMin = 1;
          yMax = 100;
          yDelta = 5;
          MaxElevation = 10; //elevation angle (determines path through Earth;
          DeltaAngle = 0.3; //steps in azimuth and elevation 
       }
    else
       {
          yMin = 10;
          yMax = 400;
          yDelta = 10;
          MaxElevation = 10; //elevation angle (determines path through Earth;
          DeltaAngle = 0.3; //steps in azimuth and elevation 
       }


    dMinimumNumberPhotoelectrons = 20;
    tanFoV = tan(5./180.*pi);
    dMinLength = 0.3; //mimnimum length a shower has to have in the camera, in degrees. This is a conservative estimate because it assumes that the shower starts at a distance l from the detector, which is not necessarily tru for showers with shallow elevation angles.

    TGraph *grSensitivity = new TGraph();

    TGraph *grDiffAcceptance = new TGraph();

  TCanvas *cIntSensitivity = new TCanvas("cIntSensitivity","Integral Sensitivity",750,500);
  cIntSensitivity->Draw();
  cIntSensitivity->SetLogy();
  cIntSensitivity->SetLogx();
  

    //Move in steps from highest energy to lowest adding acceptance
    Double_t dLogE = logEmax;
    Int_t n=0;
    Double_t dIntegratedAcceptance=0;
    while(dLogE>logEmin)
       {
	 Double_t dAcceptance = CalculateAcceptance(dLogE,dLogE+dLogEnergyStep,grDiffAcceptance,hTau,iTrigWin);
          Double_t dEIndexed = pow(10,-dLogE*nuIndex);
          dIntegratedAcceptance+=dAcceptance
             *dEIndexed*(pow(10,dLogE+dLogEnergyStep)-pow(10,dLogE));

          Double_t dnuFnu = 3 * pow(10,dLogE*2) * dEIndexed / dIntegratedAcceptance;
          cout<<"Energy "<<dLogE<<" integrated acceptance:  "<<dIntegratedAcceptance<<" converted to nuFnu: "<<dnuFnu<<" for power law with index -"<<nuIndex<<endl;
          grSensitivity->SetPoint(n,pow(10,dLogE),dnuFnu);

          n++;
          dLogE-=dLogEnergyStep;
          cIntSensitivity->cd();
          grSensitivity->Draw("alp");
          cIntSensitivity->Modified();
          cIntSensitivity->Update();
       }

  grSensitivity->SetLineWidth(3);
  grSensitivity->SetLineColor(kBlue+3);
  grSensitivity->GetXaxis()->SetTitle("energy [GeV]");
  grSensitivity->GetYaxis()->SetTitle("E^{2} dN/dE [ GeV cm^{-2} s^{-1} sr^{-1} ]");
  grSensitivity->GetYaxis()->SetTitleSize(0.04);
  grSensitivity->GetYaxis()->SetLabelSize(0.04);
  grSensitivity->GetXaxis()->SetTitleSize(0.04);
  grSensitivity->GetXaxis()->SetLabelSize(0.04);
}

////////////////////////////////////////////////////////////////////
//
//
// Calculate Differential Sensitiivty
//
//
void CalculateDifferentialSensitivity(TH1D *hTau)
{

    Double_t dLogEnergyStep = 0.2; //0.2
    Double_t dHalfEnergyBinWidth =1/2.; //in log was 1/2
    Double_t logEmin = 5; //7
    Double_t logEmax = 10.5; //11

    bCombined = kTRUE;
    yMin = 5; //5
    yMax = 500; //500
    yDelta = 5; //5
    MaxElevation = 40; //elevation angle (determines path through Earth;
    DeltaAngle = 0.05; //steps in azimuth and elevation 

    iConfig = 2; //telescope altitude
  
    //exposure
    Double_t dExposure=10*365*24*3600*0.20; //exposure time 10 years in seconds with 20% duty cycle

    Double_t dFoV = 89;  //test 0, 1, 2, 10
    tanFoV = tan(dFoV/180.*pi);
    //dFoVBelow = asin(REarth/(REarth+DetectorAltitude[iConfig]));
    dFoVBelow =  89/180.*pi; 
    iMirrorSize = 3;
    dMinimumNumberPhotoelectrons = dThreshold[iMirrorSize]/dMirrorA[iMirrorSize]; 

    dMinLength = 0.3; //mimnimum length a shower has to have in the camera, in degrees. This is a conservative estimate because it assumes that the shower starts at a distance l from the detector, which is not necessarily tru for showers with shallow elevation angles.
    TGraph *grDiffAcceptance = new TGraph();

    TGraph *grSensitivity = new TGraph();
    grSensitivity->SetLineWidth(3);
    grSensitivity->SetLineColor(kBlue+3);

    TGraph *grSensitivityNuCommunity = new TGraph();
    grSensitivityNuCommunity->SetLineWidth(3);
    grSensitivityNuCommunity->SetLineColor(kRed+3);



    TCanvas *cDiffSensitivity = new TCanvas("cDiffSensitivity","Differential Sensitivity",750,500);
    cDiffSensitivity->Draw();
    cDiffSensitivity->SetLogy();
    cDiffSensitivity->SetLogx();

    TGraph *grAcceptance = new TGraph();
    grAcceptance->SetLineWidth(3);
    grAcceptance->SetLineColor(kBlue+3);

    TGraph *grAcceptanceMonoEnergy = new TGraph();
    grAcceptanceMonoEnergy->SetLineWidth(3);
    grAcceptanceMonoEnergy->SetLineColor(kRed+3);

    TCanvas *cAcceptance = new TCanvas("cAcceptance","Acceptance",750,500);
    cAcceptance->Draw();
    cAcceptance->SetLogy();
    cAcceptance->SetLogx();
  
    TCanvas *cTriggeredAzimuthAngles = new TCanvas("cTriggeredAzimuthAngles","Triggered Azimuth Angles",750,500);
    hTriggeredAzimuthAngles->Draw("HIST");


    //Move in steps from lowest to highes energy
    Double_t dLogE = logEmin;
    Int_t n=0;
    while(dLogE<=logEmax)
       {
          bMonoNu = kFALSE; //Calculate the acceptance averaged over Enu bin
          Double_t dAcceptance = CalculateAcceptance(dLogE-dHalfEnergyBinWidth,dLogE+dHalfEnergyBinWidth,grDiffAcceptance,hTau,iTrigWin);

          if(dAcceptance>1)
           {
              grAcceptance->SetPoint(n,pow(10,dLogE),dAcceptance);

              Double_t dnuFnu = 3 * 2.44 * pow(10,dLogE*2) / dAcceptance / dExposure / (pow(10,dLogE+dHalfEnergyBinWidth)- pow(10,dLogE-dHalfEnergyBinWidth));
              cout<<"Energy "<<dLogE-dHalfEnergyBinWidth<<" to "<<dLogE+dHalfEnergyBinWidth<<" acceptance:  "<<dAcceptance<<" nuFnu: "<<dnuFnu<<" for power law with index -"<<nuIndex<<endl;
              grSensitivity->SetPoint(n,pow(10,dLogE),dnuFnu);
              cDiffSensitivity->cd();
              grSensitivity->Draw("alp");

              //Calculate the Sensitivity as done by the Nu Community
              bMonoNu = kTRUE; //calculate the acceptance at Enu bin center
              dAcceptance = CalculateAcceptance(dLogE-dHalfEnergyBinWidth,dLogE+dHalfEnergyBinWidth,grDiffAcceptance,hTau,iTrigWin);
              grAcceptanceMonoEnergy->SetPoint(n,pow(10,dLogE),dAcceptance);
              bMonoNu = kFALSE;
              dnuFnu = 3 * 2.44 / dAcceptance / dExposure / log(10) / (2*dHalfEnergyBinWidth) * pow(10,dLogE); //2.44 is from Feldman Cousin 90% confidence upper limit
              cout<<"NuCommunity Definition: "<<dnuFnu<<" acceptance at center energy: "<<dAcceptance<<endl;
              grSensitivityNuCommunity->SetPoint(n,pow(10,dLogE),dnuFnu);
              grSensitivityNuCommunity->Draw("lp");

              cDiffSensitivity->Modified();
              cDiffSensitivity->Update();

              cAcceptance->cd();
              grAcceptance->Draw("alp");
              grAcceptanceMonoEnergy->Draw("lp");
              cAcceptance->Modified();
              cAcceptance->Update();

              cTriggeredAzimuthAngles->cd();
              cTriggeredAzimuthAngles->Modified();
              cTriggeredAzimuthAngles->Update();

              n++;
            }
          dLogE+=dLogEnergyStep;
       }


  cDiffSensitivity->cd();
  TLegend *legend = new TLegend(0.53,0.7,0.75,0.88);
  TString legstr;
  legend->AddEntry(grSensitivity,"Integration","l");
  legend->AddEntry(grSensitivityNuCommunity,"Approximation","l");
  legend->Draw();

  grSensitivity->GetXaxis()->SetTitle("energy [GeV]");
  grSensitivity->GetYaxis()->SetTitle("E^{2} dN/dE [ GeV cm^{-2} s^{-1} sr^{-1} ]");
  grSensitivity->GetYaxis()->SetTitleSize(0.04);
  grSensitivity->GetYaxis()->SetLabelSize(0.04);
  grSensitivity->GetXaxis()->SetTitleSize(0.04);
  grSensitivity->GetXaxis()->SetLabelSize(0.04);

  cAcceptance->cd();
  legend = new TLegend(0.755,0.49,0.84,0.63);
  legend->AddEntry(grAcceptance,"Averaged","l");
  legend->AddEntry(grAcceptanceMonoEnergy,"Mono","l");
  legend->Draw();

  grAcceptance->GetXaxis()->SetTitle("energy [GeV]");
  grAcceptance->GetYaxis()->SetTitle("acceptance [ cm^{2} sr ]");
  grAcceptance->GetYaxis()->SetTitleSize(0.04);
  grAcceptance->GetYaxis()->SetLabelSize(0.04);
  grAcceptance->GetXaxis()->SetTitleSize(0.04);
  grAcceptance->GetXaxis()->SetLabelSize(0.04);

  hTriggeredAzimuthAngles->SetLineWidth(3);
  hTriggeredAzimuthAngles->SetLineColor(kBlue+3);
  hTriggeredAzimuthAngles->GetYaxis()->SetTitleSize(0.04);
  hTriggeredAzimuthAngles->GetYaxis()->SetLabelSize(0.04);
  hTriggeredAzimuthAngles->GetXaxis()->SetTitleSize(0.04);
  hTriggeredAzimuthAngles->GetXaxis()->SetLabelSize(0.04);
  hTriggeredAzimuthAngles->Scale(1.0/hTriggeredAzimuthAngles->Integral(),"nosw2");

  TString Filename;
  Filename.Form("SensitivityResults/StudyHighestEnergies/DifferentialSensitivityTrinity_NuTauSim_%ikmAboveGround_%0.0fsqrmMirror_%0.1fdegUpperFoV_%0.1fdegLowerFoV_%0.1fdegMinShowerLength.root",iConfig,dMirrorA[iMirrorSize],dFoV,dFoVBelow/pi*180.,dMinLength);
  cDiffSensitivity->SaveAs(Filename.Data());
  Filename.Form("SensitivityResults/StudyHighestEnergies/AcceptanceTrinity_NuTauSim_%ikmAboveGround_%0.0fsqrmMirror_%0.1fdegUpperFoV_%0.1fdegLowerFoV_%0.1fdegMinShowerLength.root",iConfig,dMirrorA[iMirrorSize],dFoV,dFoVBelow/pi*180.,dMinLength);
  cAcceptance->SaveAs(Filename.Data());
  Filename.Form("SensitivityResults/StudyHighestEnergies/TriggeredAzimuthTrinity_NuTauSim_%ikmAboveGround_%0.0fsqrmMirror_%0.1fdegUpperFoV_%0.1fdegLowerFoV_%0.1fdegMinShowerLength.root",iConfig,dMirrorA[iMirrorSize],dFoV,dFoVBelow/pi*180.,dMinLength);
  cTriggeredAzimuthAngles->SaveAs(Filename.Data());

  Double_t *dE = grSensitivity->GetX(); 
  Double_t *dF = grSensitivity->GetY(); 
  for(int i=0;i<grSensitivity->GetN();i++)
     {
       cout<<dE[i]<<"  "<<dF[i]<<endl;;
     }

}

////////////////////////////////////////////////////////////////////
// Main Program
//
//////////////////////////////////////////////////////////////////
int main (int argc, char **argv) {


  //initiate root
  TROOT root("DisplayEvts","Display Results");
  TApplication *theApp = new TApplication("App",&argc,argv);
  gROOT->ProcessLine("#include <vector>"); //need this otherwise we cannot save vectors in the root file

  //read tables to use NuTauSimResults	
  readFromTable() ;
  findAngleNumber() ;
  
  //cout<<PDecayFluorescence(1e9,10,3,170)<<endl;;
  
  //fPE = new TF1("fPE",myPEfunction_trigwindow,0,40,2);
  //fix the trigger time window time index
  //fPE->FixParameter(2,0); // 0 = 15 ns, 1 = 30 ns, 2 = 45ns
  //fPE = new TF1("fPE",myPEfunction,0,40,2);
  
  //TLegend *leg = new TLegend(0.14,0.6,0.31,0.87);
  grsCC = new TGraph(19, Esig, sigma);
  grsNC = new TGraph(19, Esig, sigmaNC);
  
  hTriggeredAzimuthAngles = new TH1D("hTriggeredAzimuthAngles","",int(180/DeltaAngleAz),0,180);
  hTriggeredAzimuthAngles->GetXaxis()->SetTitle("angle #alpha [degrees]");
  hTriggeredAzimuthAngles->GetYaxis()->SetTitle("probability");
  hTriggeredAzimuthAngles->SetLineWidth(2);
  hTriggeredAzimuthAngles->SetLineColor(kBlue+3);
  
  //plot the probability that a tau emerges when a neutrino hits the Earth
  PlotEmergenceProbability();
  cout<<"Done with Emergence Probability Plot"<<endl;
  //Double_t dmin = 1;
  //Double_t dmax = 6001.0; //thickness in km
  //Double_t xStepCoarse = 0.1;
  
  gStyle->SetOptStat(0);
  
  TH1D *hTau = new TH1D("hTau","",70,4,11);
  hTau->GetXaxis()->SetTitle("energy [GeV]");
  hTau->GetYaxis()->SetTitle("F_tau/F_nu");
  hTau->SetLineWidth(3);
  hTau->GetXaxis()->SetLabelSize(0.045);
  hTau->GetXaxis()->SetTitleOffset(1.5);
  hTau->GetYaxis()->SetTitleOffset(1.3);
  hTau->GetYaxis()->SetTitleSize(0.04);
  hTau->GetXaxis()->SetTitleSize(0.04);
  hTau->GetYaxis()->SetLabelSize(0.045);
  
  TAxis *axis = hTau->GetXaxis();
  int bins = axis->GetNbins();
  Axis_t from = axis->GetXmin();
  Axis_t to = axis->GetXmax();
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];
  for (int i = 0; i <= bins; i++) {
    new_bins[i] = TMath::Power(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  
  //making Fta/Fnu plot from dutta paper
  for(int i=0;i<hTau->GetNbinsX();i++)
    {
      Double_t El = hTau->GetBinLowEdge(i+1);
      Double_t Eh = hTau->GetBinLowEdge(i+2);
      Double_t Enu=hTau->GetBinCenter(i+1);
      
      Double_t weight = log(Eh)-log(El) ;
      
      Double_t Etau = hTau->GetBinCenter(1); 
      int n = 0;
      while(Etau<=Eh)
	{
	  //  cout<<"Energy "<<E<<endl;
	  Double_t P = PEtau(100,Etau,Enu,hTau) ;
	  hTau->Fill(Etau,weight*P);
	  n++;
	  Etau=  hTau->GetBinCenter(n+1);
	}
    }
  
  //convert to different units
  for(int i=0;i<hTau->GetNbinsX();i++)
    {
      //Double_t Enu = hTau->GetBinCenter(i+1);
      Double_t El = hTau->GetBinLowEdge(i+1);
      Double_t Eh = hTau->GetBinLowEdge(i+2);
      hTau->SetBinContent(i+1,hTau->GetBinContent(i+1)/(log(Eh)-log(El)));
      hTau->SetBinError(i+1,0);
    }
  
  //plot clone of hTau
  
  
  //produce plot with tau energy distribution for four different target
  //thicknesses
  
  TCanvas *cTauSpeMonoEn = new TCanvas("cTauSpeMonoEn","Energy Distribution of Taus produced by 10^9 GeV Nus",750,500);
  cTauSpeMonoEn->Draw();
  cTauSpeMonoEn->SetLogx();
  cTauSpeMonoEn->SetLogy();
  
  hTau->GetYaxis()->SetTitle("probability of #tau emergence");
  
  TLegend *legend = new TLegend(0.75,0.6,0.97,0.95);
  TString legstr;
  for(float i=0;i<4;i++)
    {
      GetTauDistribution(hTau,pow(10,i),pow(10,9.0),pow(10,9.1001));//9.0...9.1
      TH1D *h = (TH1D*)hTau->Clone("h");
      h->SetLineStyle(i+1);
      h->Draw("same");
      legstr.Form("%0.0f km",pow(10,i));
      legend->AddEntry(h,legstr.Data(),"l");
    }
  legend->Draw();
 
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //Sensitivity calculation starts here
  bFluorescence = kFALSE;
  
  //CalculateAcceptanceVsImageLength(hTau);
  
  //CalculateAcceptanceVsUpperFoV(hTau);
  //
  //CalculateAcceptanceVsLowerFoV(hTau);
  
  //CalculateAcceptanceVsTelescopeHeight(hTau);
  //
  CalculateAcceptanceVsTriggerWindow(hTau);
  //
  //CalculateAcceptanceVsThreshold(hTau);
  //
  //CalculateAcceptanceVsEnergy(hTau);
  //
  //CalculateIntegralSensitivity(hTau);
  //CalculateDifferentialSensitivity(hTau);
  //
 
/*
  cout<<"DEBUGGING"<<endl;
  Double_t dEarth = DistanceThroughEarth(5,0.75);
  //GetTauDistribution(hTau,dEarth,9.3,10.3);                
  GetTauDistribution(hTau,dEarth,9.5,10.5);                
  for(int i=0;i<hTau->GetNbinsX();i++)
  {
  cout<<hTau->GetBinCenter(i+1)  <<" taus cont: "<<hTau->GetBinContent(i+1)<<endl;
  }
  
*/
/*
  yMin = 5;
  yMax = 5.1;
  yDelta = 5;
  MaxElevation = 30;
  DeltaAngle = 0.05;
  dMaxCherenkovAzimuthAngle = 30;
  dMinEnu = 5.5;
  dMaxEnu = 6.5;
  
  Double_t dFoV = 3;  //test 0, 1, 2, 10
  tanFoV = tan(dFoV/180.*pi);
  //dFoVBelow = asin(REarth/(REarth+DetectorAltitude[iConfig]));
  dFoVBelow =  2/180.*pi; 
  iMirrorSize = 2;
  dMinimumNumberPhotoelectrons = dThreshold[iMirrorSize]/dMirrorA[iMirrorSize]; 
  
  dMinLength = 0.3; //mimnimum length a shower has to have in the camera, in degrees. This is a conservative estimate because it assumes that the shower starts at a distance l from the detector, which is not necessarily tru for showers with shallow elevation angles.
  TGraph *grTMP = new TGraph();
  CalculateAcceptance(dMinEnu,dMaxEnu,grTMP,hTau);
*/
 
//loop over threshold energy//
//
 
//loop over nu energy bins weighting with power law with index g and fill spectrum with tau energies
//Generate Tau spectrum for neutrinos in this bin emerging under a given angle
//
/*
  TCanvas *cSens = new TCanvas("cSens","Sensitivity",750,500);
  cSens->Draw();
  TH1D *hSens = new TH1D("hSens","",1000,0.9*dmin,dmax*1.1);
  hSens->Draw();
  hSens->SetMaximum(1e-3);
  hSens->SetMinimum(1e-10);
  hSens->GetXaxis()->SetTitle("Target thickness [km]");
  hSens->GetYaxis()->SetTitle("Sensitivity [GeV cm^{-2} s^{-1} sr^{-1}]");
  cSens->SetLogx();
  cSens->SetLogy();
  TMultiGraph *mGrSens = new TMultiGraph("mGrSens","MGrSens");
  
  
  cSens->cd();
  mGrSens->Draw("PL");
  leg->Draw();
  leg->Draw();
*/
 
 
/*
  TTimer timer("gSystem->ProcessEvents();", 50, kFALSE);
  timer.TurnOn();
  Getline("Type <return> to go on: ");
  timer.TurnOff();
*/
 cout<<"done"<<endl;
 theApp->Run();
 
 return 0;
}
