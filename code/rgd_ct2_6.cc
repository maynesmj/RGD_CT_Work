//******************************************************************
//*  ██╗  ██╗██╗██████╗  ██████╗     ██╗  ██╗    ██████╗
//*  ██║  ██║██║██╔══██╗██╔═══██╗    ██║  ██║   ██╔═████╗
//*  ███████║██║██████╔╝██║   ██║    ███████║   ██║██╔██║
//*  ██╔══██║██║██╔═══╝ ██║   ██║    ╚════██║   ████╔╝██║
//*  ██║  ██║██║██║     ╚██████╔╝         ██║██╗╚██████╔╝
//*  ╚═╝  ╚═╝╚═╝╚═╝      ╚═════╝          ╚═╝╚═╝ ╚═════╝
//************************ Jefferson National Lab (2017) ***********
//******************************************************************

/*

The latest version of the main analysis code
This code takes a list file of locations to hipo files. An example of such list has been added to this branch.
The code produces a root file containing several useful graphs and histograms

input file is located at line: 863
name of root file is located at line: 861

Lines (868,  887 - 890 ) allow you to set which target specific cuts and sector of the drift chamber you wish to look at

the loop per hipo file starts at line: 991 

the loop per event starts at line: 1143

lines ( 1448 - 1537) take the information for e, pip, & pim per event

lines ( 1588 - 1961) link the detector banks and make calulations related to fiducal cuts

lines ( 2516 - 2601 ) make all calculations related to rho0 and related invarient mass plots

*/

#include <cstdlib>
#include <iostream>
#include <string>
#include "reader.h" 
#include "TFile.h"
//#include "TTree.h"
//#include "TChain.h"
#include "TLorentzVector.h"
#include <filesystem>
#include <sys/stat.h>
#include <vector>

#include "TTree.h"
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include <THStack.h>
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include <Math/Vector3D.h>


//namespace fs = std::filesystem;
int rho_total = 0;
int d_rho_counter = 0;
int nun_rho = 0;
int event_counter = 0;
int pT = 0;
int nT = 0;


int t_e_num = 0;
   int t_pip_num = 0;
   int t_pim_num = 0;
   int t_rho_num = 0;
   float t_bC = 0.0;

   float e_yield = 0.0;
   float pip_yield = 0.0;
   float pim_yield = 0.0;
   float rho_yield = 0.0;

   float e_pf = 0.0;
   float pip_pf = 0.0;
   float pim_pf = 0.0;
   float rho_pf = 0.0;

Double_t fitf(Double_t *x, Double_t *par){
   Double_t arg = 0;
   if (par[2]!= 0) arg = (x[0] - par[1])/par[2];
   Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg) + par[3]*(x[0])*(x[0]) + par[4]*(x[0]) + par[5];
   //Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg) - 1000*(x[0])*(x[0]) + 500*(x[0]) + 700;  for vx 19125
   return fitval;
}

Double_t fitfy(Double_t *x, Double_t *par){
   Double_t arg = 0;
   if (par[2]!= 0) arg = (x[0] - par[1])/par[2];
   Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg) + par[3]*(x[0]-0.5)*(x[0]-0.5) + par[4]*(x[0]-0.5) + par[5];
   //Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg) - 1000*(x[0])*(x[0]) + 500*(x[0]) + 700;  for vx 19125
   return fitval;
}

//Breit-Wigner function
Double_t fitfbw(Double_t* x, Double_t* par)
{
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = par[1]*par[1]*par[2]*par[2]; //Gamma=par[1]  M=par[2]
  Double_t arg3 = ((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = x[0]*x[0]*x[0]*x[0]*((par[1]*par[1])/(par[2]*par[2]));
  return par[0]*arg1*arg2/(arg3 + arg4) + par[3]*(x[0])*(x[0])*(x[0]) + par[4]*(x[0])*(x[0]) + par[5]*(x[0]) + par[6];
}



inline auto angle_plane(ROOT::Math::XYZVector vA, ROOT::Math::XYZVector vB, ROOT::Math::XYZVector vC, ROOT::Math::XYZVector vD) -> double {



ROOT::Math::XYZVector crossAB = vA.Cross(vB); // AxB

ROOT::Math::XYZVector crossCD = vC.Cross(vD); // CxD



double sgn = crossAB.Dot(vD); // (AxB).D

if(std::abs(sgn)<0.00001) return -10000;

sgn /= std::abs(sgn);



// calculate numerator and denominator 

double numer = crossAB.Dot(crossCD); // (AxB).(CxD)

double denom = crossAB.R() * crossCD.R(); // |AxB|*|CxD|

if(std::abs(denom) < 0.00001) return -10000;



// return angle

return sgn * std::acos(numer/denom);

}

//This code is build to get the kinematics and invarient mass distributions. seperate sectors and well as yields comparison per run number graphs.
//This code also has the samplpe fraction code.

int main(){

   auto *c1 = new TCanvas("c1", "c1");
   TStyle *gstyle = new TStyle();
   //gstyle->SetPalette(kRainBow); 

   auto *c2 = new TCanvas("c2", "c2");
   //TStyle *gstyle = new TStyle();
   //gstyle->SetPalette(kRainBow); 

   auto *c3 = new TCanvas("c3", "c3");
   //gStyle->SetOptStat("");
	gStyle->SetCanvasPreferGL();
   TH1F *z_hist = new TH1F("Zh", " Zh ", 200, 0, 1.2);
   //bC_hist->SetStats(0);

   TH2F *qz_hist = new TH2F("lc vs Q2", "lc vs Q2", 200, 0, 10, 200, 0, 7);
   
   TH2F *xy_hist = new TH2F("px vs py", "px vs py", 200, -5, 5, 200, -5, 5);

   TH1F *q_hist = new TH1F("Q2", " Q2", 200, 0, 10);
   
   TH1F *lc_hist = new TH1F("lc", " lc", 200, 0, 3);
   
   TH1F *e_hist = new TH1F("electrons", " electrons ", 200, -300, 300);
   
   TH1F *pip_hist = new TH1F("pip", "pip", 200, -300, 300);
   
   TH1F *pim_hist = new TH1F("pim", "pim", 200, -300, 300);

   //TH2F *qlc_hist = new TH2F("Q2 vs lc", "Q2 vs lc", 200, 0, 10, 200, 0, 20);

   TH1F *w_hist = new TH1F("w (Invariant Mass)", " ", 200, 0, 5);

   //TH2F *tz_hist = new TH2F("-t vs rho_pz", "-t vs rho_pz", 200, 0, 10, 200, 0, 20);

   TH1F *t_hist = new TH1F("-t", "-t", 200, 0, 1);

   TH1F *wonphe_hist = new TH1F("wonphe", " ", 500, 0, 0.5);

   TH1F *wnphe_hist = new TH1F("wnphe", " ", 500, 0, 0.5);

   TH1F *wopcal_hist = new TH1F("Total Energy", " ", 500, 0, 1);

   TH1F *wpcal_hist = new TH1F("PCAL", " ", 500, 0, 1);

   TH2F *wocut_hist = new TH2F("wocut", "wocut", 500, 0, 20, 500, 0, 1);

   TH2F *wcut_hist = new TH2F("Energy", "", 500, 0.6, 11, 500, 0, 0.3);
   
   TH1F *e_pht_hist = new TH1F("e phi vs theta", "e phi vs theta", 200, 0, 40);
   
   TH2F *pip_pht_hist = new TH2F("pip phi vs theta", "pip phi vs theta", 200, -180.0, 180.0, 200, 0, 40);
   
   TH2F *pim_pht_hist = new TH2F("pim phi vs theta", "pim phi vs theta", 200, -180.0, 180.0, 200, 0, 40);
   
   TH1F *sam_hist = new TH1F("sm", "Sample Fractions", 200, 0, 0.4);

   TH2F *samp_hist = new TH2F("smvsp", "Sample Fraction vs momentum", 200, 0, 11, 200, 0.01, 0.4);
   samp_hist->SetStats(0);
   
   TH2F *samp1_hist = new TH2F("samp1_hist", "Sample Fraction vs momentum", 200, 0, 11, 200, 0.01, 0.4);
   //samp1_hist->SetStats(0);
   
   TH2F *samp1_cut_hist = new TH2F("samp1_cut_hist", "Sample Fraction vs momentum", 200, 0, 11, 200, 0.01, 0.4);
   //samp1_hist->SetStats(0);
   
   TH2F *samp2_hist = new TH2F("smvsp2", "Sample Fraction vs momentum", 200, 0, 11, 200, 0.01, 0.4);
   samp2_hist->SetStats(0);
   
   TH2F *samp3_hist = new TH2F("smvsp3", "Sample Fraction vs momentum", 200, 0, 11, 200, 0.01, 0.4);
   samp3_hist->SetStats(0);
   
   TH2F *samp4_hist = new TH2F("smvsp4", "Sample Fraction vs momentum", 200, 0, 11, 200, 0.01, 0.4);
   samp4_hist->SetStats(0);
   
   TH2F *samp5_hist = new TH2F("smvsp5", "Sample Fraction vs momentum", 200, 0, 11, 200, 0.01, 0.4);
   samp5_hist->SetStats(0);
   
   TH2F *samp6_hist = new TH2F("smvsp6", "Sample Fraction vs momentum", 200, 0, 11, 200, 0.01, 0.4);
   samp6_hist->SetStats(0);
   
   
   TH2F *samv_hist = new TH2F("smvsv", "Sample Fraction vs v", 200, 0, 70, 200, 0, 0.4);
   samv_hist->SetStats(0);
   TH2F *samw_hist = new TH2F("smvsw", "Sample Fraction vs w", 200, 0, 70, 200, 0, 0.4);
   samw_hist->SetStats(0);
   
   
    TH2F *thve_hist = new TH2F("thve", "#theta_{DC} vs edge", 200, 0, 60, 200, 5, 35);
   thve_hist->SetStats(0);
   
   TH2F *thve_hist1 = new TH2F("thve1", "#theta_{DC} vs edge", 200, 0, 60, 200, 5, 35);
   thve_hist1->SetStats(0);
   TH2F *thve_hist2 = new TH2F("thve2", "#theta_{DC} vs edge", 200, 0, 60, 200, 5, 35);
   thve_hist2->SetStats(0);
   TH2F *thve_hist3 = new TH2F("thve3", "#theta_{DC} vs edge", 200, 0, 60, 200, 5, 35);
   thve_hist3->SetStats(0);
   TH2F *thve_hist4 = new TH2F("thve4", "sector 4 #chi^{2} / NDF vs edge", 200, 0, 30, 200, 0, 100);
   thve_hist4->SetStats(0);
   TH2F *thve_hist5 = new TH2F("thve5", "sector 5 #chi^{2} / NDF vs edge", 200, 0, 30, 200, 0, 100);
   thve_hist5->SetStats(0);
   TH2F *thve_hist6 = new TH2F("thve6", "sector 6 #chi^{2} / NDF vs edge", 200, 0, 30, 200, 0, 100);
   thve_hist6->SetStats(0);
   
   
   TProfile *chin_hist = new TProfile("chin", "#chi^{2} / NDF vs edge", 200, 0, 30, 0, 100);
   chin_hist->SetStats(0);
   
   TProfile *chin_hist1 = new TProfile("chin1", "sector 1 #chi^{2} / NDF vs edge", 200, 0, 30, 0, 100);
   chin_hist1->SetStats(0);
   TProfile *chin_hist2 = new TProfile("chin2", "sector 2 #chi^{2} / NDF vs edge", 200, 0, 30, 0, 100);
   chin_hist2->SetStats(0);
   TProfile *chin_hist3 = new TProfile("chin3", "sector 3 #chi^{2} / NDF vs edge", 200, 0, 30, 0, 100);
   chin_hist3->SetStats(0);
   TProfile *chin_hist4 = new TProfile("chin4", "sector 4 #chi^{2} / NDF vs edge", 200, 0, 30, 0, 100);
   chin_hist4->SetStats(0);
   TProfile *chin_hist5 = new TProfile("chin5", "sector 5 #chi^{2} / NDF vs edge", 200, 0, 30, 0, 100);
   chin_hist5->SetStats(0);
   TProfile *chin_hist6 = new TProfile("chin6", "sector 6 #chi^{2} / NDF vs edge", 200, 0, 30, 0, 100);
   chin_hist6->SetStats(0);
   
   
   

   auto *vxs = new THStack("vxs", "vxs");
   auto *vys = new THStack("vys", "vys");
   auto *etots = new THStack("etots", "etots");
   auto *ps = new THStack("ps" , "ps");

int e_count = 0;
int pip_count = 0;
int pim_count = 0;
int zerosp = 0;
int total_r = 0;
int r_count = 0;
int rw_count = 0;
int rwt_count = 0;
int rwtz_count = 0;
int r1_count = 0;
int r2_count = 0;
int r3_count = 0;
int r4_count = 0;
int r5_count = 0;
int r6_count = 0;
int bin1_count = 0;
int bin2_count = 0;
int bin3_count = 0;
int bin4_count = 0;
int bin5_count = 0;
int bin6_count = 0;


TH1F *ep_hist = new TH1F("Electron Momentum", "Electron Momentum", 200, 0, 11);
       ep_hist->SetLineColor(kBlack);
TH1F *ewp_hist = new TH1F("Electron Momentum With W cut", "Electron Momentum With W cut", 200, 0, 11);
       ewp_hist->SetLineColor(kBlack);
TH1F *ewtp_hist = new TH1F("Electron Momentum With W & -t cut", "Electron Momentum With W & -t cut", 200, 0, 11);
       ewtp_hist->SetLineColor(kBlack);
TH1F *ewtzp_hist = new TH1F("Electron Momentum With W, -t & Zh cut", "Electron Momentum With W, -t & Zh cut", 200, 0, 11);
       ewtzp_hist->SetLineColor(kBlack);
       
       
TH1F *pipp_hist = new TH1F("#Pi^{+} Momentum", "#Pi^{+} Momentum", 200, 0, 10);
       pipp_hist->SetLineColor(kBlack);
TH1F *pipwp_hist = new TH1F("#Pi^{+} Momentum With W cut", "#Pi^{+} Momentum With W cut", 200, 0, 10);
       pipwp_hist->SetLineColor(kBlack);
TH1F *pipwtp_hist = new TH1F("#Pi^{+} Momentum With W & -t cut", "#Pi^{+} Momentum With W & -t cut", 200, 0, 10);
       pipwtp_hist->SetLineColor(kBlack);
TH1F *pipwtzp_hist = new TH1F("#Pi^{+} Momentum With W, -t & Zh cut", "#Pi^{+} Momentum With W, -t & Zh cut", 200, 0, 10);
       pipwtzp_hist->SetLineColor(kBlack);
      
TH1F *pimp_hist = new TH1F("#Pi^{-} Momentum", "#Pi^{-} Momentum", 200, 0, 10);
       pimp_hist->SetLineColor(kBlack);
TH1F *pimwp_hist = new TH1F("#Pi^{-} Momentum With W cut", "#Pi^{-} Momentum With W cut", 200, 0, 10);
       pimwp_hist->SetLineColor(kBlack);
TH1F *pimwtp_hist = new TH1F("#Pi^{-} Momentum With W & -t cut", "#Pi^{-} Momentum With W & -t cut", 200, 0, 10);
       pimwtp_hist->SetLineColor(kBlack);
TH1F *pimwtzp_hist = new TH1F("#Pi^{-} Momentum With W, -t & Zh cut", "#Pi^{-} Momentum With W, -t & Zh cut", 200, 0, 10);
       pimwtzp_hist->SetLineColor(kBlack);   





   TH1F *vx_hist = new TH1F(" vx", " vx", 200, -5, 5);
       vx_hist->SetLineColor(kBlack);
   TH1F *vy_hist = new TH1F(" vy", " vy", 200, -5, 5);
       vy_hist->SetLineColor(kBlack);
   //TH2F *xy_hist = new TH2F(" px .vs py", " px .vs py", 200, -1, 1, 200, -1, 1);
   TH1F *vz_hist = new TH1F(" vz", " vz", 200, -30, 20);
          vz_hist->SetLineColor(kBlack);
   TH1F *fmt_vz_hist = new TH1F("fmt vz", "fmt vz", 200, -30, 20);
         fmt_vz_hist->SetLineColor(kBlack);    
     
     TH1F *fmt_vz1_hist = new TH1F("fmt vz1", "fmt vz1", 200, -30, 20);
         fmt_vz1_hist->SetLineColor(kRed); 
         
     TH1F *fmt_vz2_hist = new TH1F("fmt vz2", "fmt vz2", 200, -30, 20);
         fmt_vz2_hist->SetLineColor(kOrange);
         
         TH1F *fmt_vz3_hist = new TH1F("fmt vz3", "fmt vz3", 200, -30, 20);
         fmt_vz3_hist->SetLineColor(kYellow); 
         
         TH1F *fmt_vz4_hist = new TH1F("fmt vz4", "fmt vz4", 200, -30, 20);
         fmt_vz4_hist->SetLineColor(kGreen); 
         
         TH1F *fmt_vz5_hist = new TH1F("fmt vz5", "fmt vz5", 200, -30, 20);
         fmt_vz5_hist->SetLineColor(kBlue); 
         
         TH1F *fmt_vz6_hist = new TH1F("fmt vz6", "fmt vz6", 200, -30, 20);
         fmt_vz6_hist->SetLineColor(kMagenta); 
          
          TH1F *vxx_hist = new TH1F(" vxx", " vxx", 200, -5, 5);
       vxx_hist->SetLineColor(kBlack);
   TH1F *vyy_hist = new TH1F(" vyy", " vyy", 200, -5, 5);
       vyy_hist->SetLineColor(kBlack);
   //TH2F *xy_hist = new TH2F(" px .vs py", " px .vs py", 200, -1, 1, 200, -1, 1);
   TH1F *vzz_hist = new TH1F(" vzz", " vzz", 200, -30, 20);
          vzz_hist->SetLineColor(kBlack);
   
   TH1F *pip_vx_hist = new TH1F("pip vx", "pip vx", 200, -5, 5);
       pip_vx_hist->SetLineColor(kBlack);
   TH1F *pip_vy_hist = new TH1F("pip vy", "pip vy", 200, -5, 5);
       pip_vy_hist->SetLineColor(kBlack);
   TH2F *pip_xy_hist = new TH2F("sector 1 px .vs py", "sector 1 px .vs py", 200, -1, 1, 200, -1, 1);
   TH1F *pip_vz_hist = new TH1F("pip vz", "pip vz", 200, -30, 20);
          pip_vz_hist->SetLineColor(kBlack);
   
   TH1F *pim_vx_hist = new TH1F("pim vx", "pim vx", 200, -5, 5);
       pim_vx_hist->SetLineColor(kBlack);
   TH1F *pim_vy_hist = new TH1F("pim vy", "pim vy", 200, -5, 5);
       pim_vy_hist->SetLineColor(kBlack);
   TH2F *pim_xy_hist = new TH2F("sector 1 px .vs py", "sector 1 px .vs py", 200, -1, 1, 200, -1, 1);
   TH1F *pim_vz_hist = new TH1F("pim vz", "pim vz", 200, -30, 20);
          pim_vz_hist->SetLineColor(kBlack);
   
   
   

   TH1F *ivx_hist = new TH1F("sector 1 vx", "Sector 1 vx", 200, -5, 5);
       ivx_hist->SetLineColor(kCyan+1);
   TH1F *ivy_hist = new TH1F("sector 1 vy", "Sector 1 vy", 200, -5, 5);
       ivy_hist->SetLineColor(kCyan+1);
   TH2F *ixy_hist = new TH2F("sector 1 px .vs py", "sector 1 px .vs py", 200, -1, 1, 200, -1, 1);
   TH1F *ivz_hist = new TH1F("sector 1 vz", "Sector 1 vz", 200, -30, 20);
   ivz_hist->SetLineColor(kRed+1);

   TH1F *iivx_hist = new TH1F("sector 2 vx", "Sector 2 vx", 200, -5, 5);
	iivx_hist->SetLineColor(kOrange+1);
   TH1F *iivy_hist = new TH1F("sector 2 vy", "Sector 2 vy", 200, -5, 5);
	iivy_hist->SetLineColor(kOrange+1);
   TH2F *iixy_hist = new TH2F("sector 2 vx .vs vy", "" , 500, 0.6, 11, 500, 0, 0.3);
   TH1F *iivz_hist = new TH1F("sector 2 vz", "Sector 2 vz", 200, -30, 20);
   	iivz_hist->SetLineColor(kOrange+1);

   TH1F *iiivx_hist = new TH1F("sector 3 vx", "Sector 3 vx", 200, -5, 5);
	iiivx_hist->SetLineColor(kYellow+1);
   TH1F *iiivy_hist = new TH1F("sector 3 vy", "Sector 3 vy", 200, -5, 5);
	iiivy_hist->SetLineColor(kYellow+1);
   TH2F *iiixy_hist = new TH2F("sector 3 vx .vs vy", "sector 3 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *iiivz_hist = new TH1F("sector 3 vz", "Sector 3 vz", 200, -30, 20);
   	iiivz_hist->SetLineColor(kYellow+1);

   TH1F *ivvx_hist = new TH1F("sector 4 vx", "Sector 4 vx", 200, -5, 5);
	ivvx_hist->SetLineColor(kOrange + 2);
   TH1F *ivvy_hist = new TH1F("sector 4 vy", "Sector 4 vy", 200, -5, 5);
	ivvy_hist->SetLineColor(kOrange + 2);
   TH2F *ivxy_hist = new TH2F("sector 4 vx .vs vy", "sector 4 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *ivvz_hist = new TH1F("sector 4 vz", "Sector 4 vz", 200, -30, 20);
   	ivvz_hist->SetLineColor(kGreen + 2);

   TH1F *vvx_hist = new TH1F("sector 5 vx", "Sector 5 vx", 200, -5, 5);
	vvx_hist->SetLineColor(kOrange + 2);
   TH1F *vvy_hist = new TH1F("sector 5 vy", "Sector 5 vy", 200, -5, 5);
	vvy_hist->SetLineColor(kOrange + 2);
   TH2F *vxy_hist = new TH2F("sector 5 vx .vs vy", "sector 5 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *vvz_hist = new TH1F("sector 5 vz", "Sector 5 vz", 200, -30, 20);
   	vvz_hist->SetLineColor(kBlue + 2);

   TH1F *vivx_hist = new TH1F("sector 6 vx", "Sector 6 vx", 200, -5, 5);
	vivx_hist->SetLineColor(kMagenta+1);
   TH1F *vivy_hist = new TH1F("sector 6 vy", "Sector 6 vy", 200, -5, 5);
        vivy_hist->SetLineColor(kMagenta+1);
   TH2F *vixy_hist = new TH2F("sector 6 vx .vs vy", "sector 6 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *vivz_hist = new TH1F("sector 6 vz", "Sector 6 vz", 200, -30, 20);
   	vivz_hist->SetLineColor(kMagenta+1);
   
   
   
   
   TH1F *pip_ivx_hist = new TH1F("pip_sector 1 vx", "pip_Sector 1 vx", 200, -5, 5);
       pip_ivx_hist->SetLineColor(kCyan+1);
   TH1F *pip_ivy_hist = new TH1F("pip_sector 1 vy", "pip_Sector 1 vy", 200, -5, 5);
       pip_ivy_hist->SetLineColor(kCyan+1);
   TH2F *pip_ixy_hist = new TH2F("pip_sector 1 px .vs py", "pip_sector 1 px .vs py", 200, -1, 1, 200, -1, 1);
   TH1F *pip_ivz_hist = new TH1F("pip_sector 1 vz", "pip_Sector 1 vz", 200, -30, 20);
   	pip_ivz_hist->SetLineColor(kCyan+1);

   TH1F *pip_iivx_hist = new TH1F("pip_sector 2 vx", "pip_Sector 2 vx", 200, -5, 5);
	pip_iivx_hist->SetLineColor(kOrange+1);
   TH1F *pip_iivy_hist = new TH1F("pip_sector 2 vy", "pip_Sector 2 vy", 200, -5, 5);
	pip_iivy_hist->SetLineColor(kOrange+1);
   TH2F *pip_iixy_hist = new TH2F("pip_sector 2 vx .vs vy", "pip_vxvy" , 500, 0.6, 11, 500, 0, 0.3);
   TH1F *pip_iivz_hist = new TH1F("pip_sector 2 vz", "pip_Sector 2 vz", 200, -30, 20);
   	pip_iivz_hist->SetLineColor(kOrange+1);

   TH1F *pip_iiivx_hist = new TH1F("pip_sector 3 vx", "pip_Sector 3 vx", 200, -5, 5);
	pip_iiivx_hist->SetLineColor(kYellow+1);
   TH1F *pip_iiivy_hist = new TH1F("pip_sector 3 vy", "pip_Sector 3 vy", 200, -5, 5);
	pip_iiivy_hist->SetLineColor(kYellow+1);
   TH2F *pip_iiixy_hist = new TH2F("pip_sector 3 vx .vs vy", "pip_sector 3 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *pip_iiivz_hist = new TH1F("pip_sector 3 vz", "pip_Sector 3 vz", 200, -30, 20);
   	pip_iiivz_hist->SetLineColor(kYellow+1);

   TH1F *pip_ivvx_hist = new TH1F("pip_sector 4 vx", "pip_Sector 4 vx", 200, -5, 5);
	pip_ivvx_hist->SetLineColor(kOrange + 2);
   TH1F *pip_ivvy_hist = new TH1F("pip_sector 4 vy", "pip_Sector 4 vy", 200, -5, 5);
	pip_ivvy_hist->SetLineColor(kOrange + 2);
   TH2F *pip_ivxy_hist = new TH2F("pip_sector 4 vx .vs vy", "pip_sector 4 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *pip_ivvz_hist = new TH1F("pip_sector 4 vz", "pip_Sector 4 vz", 200, -30, 20);
   	pip_ivvz_hist->SetLineColor(kOrange + 2);

   TH1F *pip_vvx_hist = new TH1F("pip_sector 5 vx", "pip_Sector 5 vx", 200, -5, 5);
	pip_vvx_hist->SetLineColor(kOrange + 2);
   TH1F *pip_vvy_hist = new TH1F("pip_sector 5 vy", "pip_Sector 5 vy", 200, -5, 5);
	pip_vvy_hist->SetLineColor(kOrange + 2);
   TH2F *pip_vxy_hist = new TH2F("pip_sector 5 vx .vs vy", "pip_sector 5 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *pip_vvz_hist = new TH1F("pip_sector 5 vz", "pip_Sector 5 vz", 200, -30, 20);
   	pip_vvz_hist->SetLineColor(kOrange + 2);

   TH1F *pip_vivx_hist = new TH1F("pip_sector 6 vx", "pip_Sector 6 vx", 200, -5, 5);
	pip_vivx_hist->SetLineColor(kMagenta+1);
   TH1F *pip_vivy_hist = new TH1F("pip_sector 6 vy", "pip_Sector 6 vy", 200, -5, 5);
        pip_vivy_hist->SetLineColor(kMagenta+1);
   TH2F *pip_vixy_hist = new TH2F("pip_sector 6 vx .vs vy", "pip_sector 6 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *pip_vivz_hist = new TH1F("pip_sector 6 vz", "pip_Sector 6 vz", 200, -30, 20);
   	pip_vivz_hist->SetLineColor(kMagenta+1);
   
   
   
   
   
   TH1F *pim_ivx_hist = new TH1F("pim_sector 1 vx", "pim_Sector 1 vx", 200, -5, 5);
         pim_ivx_hist->SetLineColor(kCyan+1);
   TH1F *pim_ivy_hist = new TH1F("pim_sector 1 vy", "pim_Sector 1 vy", 200, -5, 5);
         pim_ivy_hist->SetLineColor(kCyan+1);
   TH2F *pim_ixy_hist = new TH2F("pim_sector 1 px .vs py", "pim_sector 1 px .vs py", 200, -1, 1, 200, -1, 1);
   TH1F *pim_ivz_hist = new TH1F("pim_sector 1 vz", "pim_Sector 1 vz", 200, -30, 20);
   	pim_ivz_hist->SetLineColor(kCyan+1);

   TH1F *pim_iivx_hist = new TH1F("pim_sector 2 vx", "pim_Sector 2 vx", 200, -5, 5);
	 pim_iivx_hist->SetLineColor(kOrange+1);
   TH1F *pim_iivy_hist = new TH1F("pim_sector 2 vy", "pim_Sector 2 vy", 200, -5, 5);
	 pim_iivy_hist->SetLineColor(kOrange+1);
   TH2F *pim_iixy_hist = new TH2F("pim_sector 2 vx .vs vy", "pim_vxvy" , 500, 0.6, 11, 500, 0, 0.3);
   TH1F *pim_iivz_hist = new TH1F("pim_sector 2 vz", "pim_Sector 2 vz", 200, -30, 20);
   	pim_iivz_hist->SetLineColor(kOrange+1);

   TH1F *pim_iiivx_hist = new TH1F("pim_sector 3 vx", "pim_Sector 3 vx", 200, -5, 5);
	 pim_iiivx_hist->SetLineColor(kYellow+1);
   TH1F *pim_iiivy_hist = new TH1F("pim_sector 3 vy", "pim_Sector 3 vy", 200, -5, 5);
	 pim_iiivy_hist->SetLineColor(kYellow+1);
   TH2F *pim_iiixy_hist = new TH2F("pim_sector 3 vx .vs vy", "pim_sector 3 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *pim_iiivz_hist = new TH1F("pim_sector 3 vz", "pim_Sector 3 vz", 200, -30, 20);
   	pim_iiivz_hist->SetLineColor(kYellow+1);

   TH1F *pim_ivvx_hist = new TH1F("pim_sector 4 vx", "pim_Sector 4 vx", 200, -5, 5);
	 pim_ivvx_hist->SetLineColor(kOrange + 2);
   TH1F *pim_ivvy_hist = new TH1F("pim_sector 4 vy", "pim_Sector 4 vy", 200, -5, 5);
	 pim_ivvy_hist->SetLineColor(kOrange + 2);
   TH2F *pim_ivxy_hist = new TH2F("pim_sector 4 vx .vs vy", "pim_sector 4 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *pim_ivvz_hist = new TH1F("pim_sector 4 vz", "pim_Sector 4 vz", 200, -30, 20);
   	pim_ivvz_hist->SetLineColor(kOrange + 2);

   TH1F *pim_vvx_hist = new TH1F("pim_sector 5 vx", "pim_Sector 5 vx", 200, -5, 5);
	 pim_vvx_hist->SetLineColor(kOrange + 2);
   TH1F *pim_vvy_hist = new TH1F("pim_sector 5 vy", "pim_Sector 5 vy", 200, -5, 5);
	 pim_vvy_hist->SetLineColor(kOrange + 2);
   TH2F *pim_vxy_hist = new TH2F("pim_sector 5 vx .vs vy", "pim_sector 5 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *pim_vvz_hist = new TH1F("pim_sector 5 vz", "pim_Sector 5 vz", 200, -30, 20);
   	pim_vvz_hist->SetLineColor(kOrange + 2);

   TH1F *pim_vivx_hist = new TH1F("pim_sector 6 vx", "pim_Sector 6 vx", 200, -5, 5);
	 pim_vivx_hist->SetLineColor(kMagenta+1);
   TH1F *pim_vivy_hist = new TH1F("pim_sector 6 vy", "pim_Sector 6 vy", 200, -5, 5);
         pim_vivy_hist->SetLineColor(kMagenta+1);
   TH2F *pim_vixy_hist = new TH2F("pim_sector 6 vx .vs vy", "pim_sector 6 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *pim_vivz_hist = new TH1F("pim_sector 6 vz", "pim_Sector 6 vz", 200, -30, 20);
   	pim_vivz_hist->SetLineColor(kMagenta+1);
   
   TH1F *echi_hist = new TH1F("echi", "electron #Chi^{2}", 200, -30, 30);
      echi_hist->SetStats(0);
   TH1F *pipchi_hist = new TH1F("pipchi", "#pi^{+} #Chi^{2}", 200, -30, 30);
      pipchi_hist->SetStats(0);
   TH1F *pimchi_hist = new TH1F("pimchi", "#pi^{-} #Chi^{2}", 200, -30, 30);
      pimchi_hist->SetStats(0);

   TH2F *tz_hist = new TH2F("vz vs phi", "vz vs phi", 200, -180, 180, 200, -30, 20);

   

   TH1F *rm_no = new TH1F("Without Cuts", "Without Cuts", 200, 0.2, 1.4);
      rm_no->SetStats(0);
   TH1F *rm_w = new TH1F("With W Cut", "With W Cut ", 200, 0.2, 1.4);
      rm_w->SetStats(0);
   TH1F *rm_wt = new TH1F("With W & -t Cut", "With W & -t Cut ", 200, 0.2, 1.4);
      rm_wt->SetStats(0);
   TH1F *rm_wtz = new TH1F("With W, -t & z Cut", " With W, -t & z Cut", 200, 0.2, 1.4);
      rm_wtz->SetStats(0);

   //auto sdaq_ev = new TGraph();
   //   sdaq_ev->SetLineColor(kRed);
   
   auto edaq_ev = new TGraph();
      edaq_ev->SetMarkerColor(kOrange + 2);
      edaq_ev->SetMarkerStyle(kFullCircle);
      edaq_ev->SetLineWidth(0);


   auto ietot_p = new TGraph();
      ietot_p->SetMarkerColor(kRed);
      ietot_p->SetMarkerStyle(kFullCircle);
      ietot_p->SetLineWidth(0);

   auto iietot_p = new TGraph();
      iietot_p->SetMarkerColor(kOrange);
      iietot_p->SetMarkerStyle(kFullCircle);
      iietot_p->SetLineWidth(0);

   auto iiietot_p = new TGraph();
      iiietot_p->SetMarkerColor(kYellow);
      iiietot_p->SetMarkerStyle(kFullCircle);
      iiietot_p->SetLineWidth(0);

   auto ivetot_p = new TGraph();
      ivetot_p->SetMarkerColor(kOrange + 2);
      ivetot_p->SetMarkerStyle(kFullCircle);
      ivetot_p->SetLineWidth(0);

   auto vetot_p = new TGraph();
      vetot_p->SetMarkerColor(kOrange + 2);
      vetot_p->SetMarkerStyle(kFullCircle);
      vetot_p->SetLineWidth(0);

   auto vietot_p = new TGraph();
      vietot_p->SetMarkerColor(kBlack);
      vietot_p->SetMarkerStyle(kFullCircle);
      vietot_p->SetLineWidth(0);









   
   auto ebeam_ev = new TGraph();
      ebeam_ev->SetMarkerColor(kRed);
      ebeam_ev->SetMarkerStyle(kFullCircle);
      ebeam_ev->SetLineWidth(0);

   auto fcup_ev = new TGraph();
      fcup_ev->SetMarkerColor(kRed);
      fcup_ev->SetMarkerStyle(kFullCircle);
      fcup_ev->SetLineWidth(0);

   auto fcupg_ev = new TGraph();
      fcupg_ev->SetMarkerColor(kRed);
      fcupg_ev->SetMarkerStyle(kFullCircle);
      fcupg_ev->SetLineWidth(0);

   auto hel_fcup = new TGraph();
      hel_fcup->SetMarkerColor(kOrange + 2);
      hel_fcup->SetMarkerStyle(kFullCircle);
      hel_fcup->SetLineWidth(0);

   auto hel_fcupgated = new TGraph();
      hel_fcupgated->SetMarkerColor(kBlack);
      hel_fcupgated->SetMarkerStyle(kFullCircle);
      hel_fcupgated->SetLineWidth(0);

   auto hel_clock = new TGraph();
      hel_clock->SetMarkerColor(kOrange + 2);
      hel_clock->SetMarkerStyle(kFullCircle);
      hel_clock->SetLineWidth(0);

   auto hel_clockgated = new TGraph();
      hel_clockgated->SetMarkerColor(kBlack);
      hel_clockgated->SetMarkerStyle(kFullCircle);
      hel_clockgated->SetLineWidth(0);
      
      auto tvz1_g = new TGraph();
      tvz1_g->SetName("tvz1_g");
      tvz1_g->SetTitle("Sector 1");
      tvz1_g->SetMarkerColor(kMagenta);
      tvz1_g->SetMarkerStyle(kFullCircle);
      tvz1_g->SetLineWidth(0);
      
      auto tvz2_g = new TGraph();
      tvz2_g->SetName("tvz2_g");
      tvz2_g->SetTitle("Sector 2");
      tvz2_g->SetMarkerColor(kBlue);
      tvz2_g->SetMarkerStyle(kFullCircle);
      tvz2_g->SetLineWidth(0);
      
      auto tvz3_g = new TGraph();
      tvz3_g->SetName("tvz3_g");
      tvz3_g->SetTitle("Sector 3");
      tvz3_g->SetMarkerColor(kOrange + 2);
      tvz3_g->SetMarkerStyle(kFullCircle);
      tvz3_g->SetLineWidth(0);
      
      
      auto tvz4_g = new TGraph();
      tvz4_g->SetName("tvz4_g");
      tvz4_g->SetTitle("Sector 4");
      tvz4_g->SetMarkerColor(kYellow);
      tvz4_g->SetMarkerStyle(kFullCircle);
      tvz4_g->SetLineWidth(0);
      
      auto tvz5_g = new TGraph();
      tvz5_g->SetName("tvz5_g");
      tvz5_g->SetTitle("Sector 5");
      tvz5_g->SetMarkerColor(kOrange);
      tvz5_g->SetMarkerStyle(kFullCircle);
      tvz5_g->SetLineWidth(0);
      
      auto tvz6_g = new TGraph();
      tvz6_g->SetName("tvz6_g");
      tvz6_g->SetTitle("Sector 6");
      tvz6_g->SetMarkerColor(kRed);
      tvz6_g->SetMarkerStyle(kFullCircle);
      tvz6_g->SetLineWidth(0);
      
      
      auto te_g = new TGraph();
      te_g->SetName("te_g");
      te_g->SetTitle("Electron Yield");
      te_g->SetMarkerColor(kYellow);
      te_g->SetMarkerStyle(kFullCircle);
      te_g->SetLineWidth(0);
      
      auto tpip_g = new TGraph();
      tpip_g->SetName("tpip_g");
      tpip_g->SetTitle("#pi_{+} Yield");
      tpip_g->SetMarkerColor(kOrange);
      tpip_g->SetMarkerStyle(kFullCircle);
      tpip_g->SetLineWidth(0);
      
      auto tpim_g = new TGraph();
      tpim_g->SetName("tpim_g");
      tpim_g->SetTitle("#pi_{-} Yield");
      tpim_g->SetMarkerColor(kRed);
      tpim_g->SetMarkerStyle(kFullCircle);
      tpim_g->SetLineWidth(0);
      
      
      auto tr1_g = new TGraph();
      tr1_g->SetName("tr1_g");
      tr1_g->SetTitle("Q^{2} bin 1");
      tr1_g->SetMarkerColor(kMagenta);
      tr1_g->SetMarkerStyle(kFullCircle);
      tr1_g->SetLineWidth(0);
      
      auto tr2_g = new TGraph();
      tr2_g->SetName("tr2_g");
      tr2_g->SetTitle("Q^{2} bin 2");
      tr2_g->SetMarkerColor(kBlue);
      tr2_g->SetMarkerStyle(kFullCircle);
      tr2_g->SetLineWidth(0);
      
      auto tr3_g = new TGraph();
      tr3_g->SetName("tr3_g");
      tr3_g->SetTitle("Q^{2} bin 3");
      tr3_g->SetMarkerColor(kOrange + 2);
      tr3_g->SetMarkerStyle(kFullCircle);
      tr3_g->SetLineWidth(0);
      
      
      auto tr4_g = new TGraph();
      tr4_g->SetName("tr4_g");
      tr4_g->SetTitle("Q^{2} bin 4");
      tr4_g->SetMarkerColor(kYellow);
      tr4_g->SetMarkerStyle(kFullCircle);
      tr4_g->SetLineWidth(0);
      
      auto tr5_g = new TGraph();
      tr5_g->SetName("tr5_g");
      tr5_g->SetTitle("Q^{2} bin 5");
      tr5_g->SetMarkerColor(kOrange);
      tr5_g->SetMarkerStyle(kFullCircle);
      tr5_g->SetLineWidth(0);
      
      auto tr6_g = new TGraph();
      tr6_g->SetName("tr6_g");
      tr6_g->SetTitle("Q^{2} bin 6");
      tr6_g->SetMarkerColor(kRed);
      tr6_g->SetMarkerStyle(kFullCircle);
      tr6_g->SetLineWidth(0);
      
      
      
      
      
      TH1D *phih_hist = new TH1D("phih_hist","#phi_{h}",200,-180,180);
      TH1D *xbj_hist = new TH1D("xbj_hist","X",200,-4,4);
      TH1D *ybj_hist = new TH1D("ybj_hist","Y",200,-4,4);

   TH1F *ietot_hist = new TH1F("sector 1 E_{tot}", "Sector 1 E_{tot}", 200, 0, 0.4);
       ietot_hist->SetFillColor(kCyan+1);
   TH1F *ip_hist = new TH1F("Q2 (2->2.5)", "Q2 (2->2.5)", 200, 0.2, 2);
       ip_hist->SetFillColor(kOrange + 2);
   TH1F *p_hist = new TH1F("Q2 (1->2)", "Q2 (1->2)", 200, 0.2, 2);
       p_hist->SetFillColor(kOrange + 2);

   TH1F *iietot_hist = new TH1F("sector 2 E_{tot}", "Sector 2 E_{tot}", 200, 0, 0.4);
	iietot_hist->SetFillColor(kOrange+1);
   TH1F *iip_hist = new TH1F("Q2 (2.5->3)", "Q2 (2.5->3)", 200, 0.2, 2);
	iip_hist->SetFillColor(kOrange + 2);
   

   TH1F *iiietot_hist = new TH1F("sector 3 E_{tot}", "Sector 3 E_{tot}", 200, 0, 0.4);
	iiietot_hist->SetFillColor(kYellow+1);
   TH1F *iiip_hist = new TH1F("Q2 (3->3.5)", "Q2 (3->3.5)", 200, 0.2, 2);
	iiip_hist->SetFillColor(kOrange + 2);
   

   TH1F *ivetot_hist = new TH1F("sector 4 E_{tot}", "Sector 4 E_{tot}", 200, 0, 0.4);
	ivetot_hist->SetFillColor(kMagenta);
   TH1F *ivp_hist = new TH1F("Q2 (3.5->4.5)", "Q2 (3.5->4.5)", 200, 0.2, 2);
	ivp_hist->SetFillColor(kOrange + 2);
   

   TH1F *vetot_hist = new TH1F("sector 5 E_{tot}", "Sector 5 E_{tot}", 200, 0, 0.4);
	vetot_hist->SetFillColor(kOrange + 2);
   TH1F *vp_hist = new TH1F("Q2 (4.5->6)", "Q2 (4.5->6)", 200, 0.2, 2);
	vp_hist->SetFillColor(kOrange + 2);
	
   TH1F *e_phi = new TH1F("e phi", "electron #phi", 200, -180, 180);
   
   TH1F *e_theta = new TH1F("e theta", "electron #theta", 200, 0, 90);
   
   TH1F *pip_phi = new TH1F("pip phi", "#pi^{+} #phi", 200, -180, 180);
   
   TH1F *pip_theta = new TH1F("pip theta", "#pi^{+} #theta", 200, 0, 90);
   
   TH1F *pim_phi = new TH1F("pim phi", "#pi^{-} #phi", 200, -180, 180);
   
   TH1F *pim_theta = new TH1F("pim theta", "#pi^{-} #theta", 200, 0, 90);
   
   //TH1F *pim_theta = new TH1F("pim theta", "#pi^{-} #theta", 200, 0, 90);
   
   
   TH2F *xy_uncut1 = new TH2F("xy_uncut1", "OB LD_{2} DC Region 1 Electron Uncut", 200, -150, 150, 200, -150, 150);
   TH2F *xy_cut1 = new TH2F("xy_cut1", "OB LD_{2} DC Region 1 Electron Cut", 200, -150, 150, 200, -150, 150);
   
   TH2F *xy_uncut2 = new TH2F("xy_uncut2", "OB LD_{2} DC Region 1 Electron Uncut", 200, -150, 150, 200, -150, 150);
   TH2F *xy_cut2 = new TH2F("xy_cut2", "OB LD_{2} DC Region 1 Electron Cut", 200, -150, 150, 200, -150, 150);
   
   TH2F *xy_uncut3 = new TH2F("xy_uncut3", "OB LD_{2} DC Region 1 Electron Uncut", 200, -150, 150, 200, -150, 150);
   TH2F *xy_cut3 = new TH2F("xy_cut3", "OB LD_{2} DC Region 1 Electron Cut", 200, -150, 150, 200, -150, 150);
   
   TH2F *xy_uncut4 = new TH2F("xy_uncut4", "OB LD_{2} DC Region 1 Electron Uncut", 200, -150, 150, 200, -150, 150);
   TH2F *xy_cut4 = new TH2F("xy_cut4", "OB LD_{2} DC Region 1 Electron Cut", 200, -150, 150, 200, -150, 150);
   
   TH2F *xy_uncut5 = new TH2F("xy_uncut5", "OB LD_{2} DC Region 1 Electron Uncut", 200, -150, 150, 200, -150, 150);
   TH2F *xy_cut5 = new TH2F("xy_cut5", "OB LD_{2} DC Region 1 Electron Cut", 200, -150, 150, 200, -150, 150);
   
   TH2F *xy_uncut6 = new TH2F("xy_uncut6", "OB LD_{2} DC Region 1 Electron Uncut", 200, -150, 150, 200, -150, 150);
   TH2F *xy_cut6 = new TH2F("xy_cut6", "OB LD_{2} DC Region 1 Electron Cut", 200, -150, 150, 200, -150, 150);
   
   
   

   TH1F *vietot_hist = new TH1F("sector 6 E_{tot}", "Sector 6 E_{tot}", 200, 0, 0.4);
	vietot_hist->SetFillColor(kMagenta+1);
   TH1F *vip_hist = new TH1F("sector 6 P", "Q2 (4.5->6)", 200, 0, 2);
        vip_hist->SetFillColor(kOrange + 2);

     float zl_b = -300000000000000.0 ; //-10.0 
     float zu_b = 1000000000000.0; // -5.0
     float pip_zl_b = -20000000000000.0 ; //-10.0 
     float pip_zu_b = 100000000000000.0; // -5.0
     float pim_zl_b = -20000000000000.0 ; //-10.0 
     float pim_zu_b = 100000000000000.0; // -5.0
   //TH2F *z_hist = new TH2F("px vs py", "px vs py", 200, -2, 2, 200, -2, 2);

   //auto cgr = new TGraph();
  // std::cout << "1 th" << std::endl;

   std::vector<std::string> v3_a;
   
   std::vector<int> v4_a;
   std::vector<double> v5_a;
   std::vector<int> v6_a;
   int num = 0;
   
   int ecount_ful = 0;
   int ecount_q = 0;
   
//47, 48, 50, 52, 53, 54, 55, 56, 58, 59, 60
//670, 671, 672, 673, 675, 676, 677, 678, 679
//9110, 9111, 9112, 9113, 9114, 9115, 9117, 9118, 9119

   TFile f("cctest_pass0v11_ldob_48_r1.root", "update"); //cusn_pass0v7_63_test.root

    std::ifstream v3("ld_pass0v11_cskim.list");  // all_pass0v11_ib_cskim.list  rge_align.list    cusn_calv9_skim.list   cxc_pass0v10_cskim.list   ld_pass0v10_cskim.list    cusn_pass0v10_cskim.list
    std::ifstream v4("ld_pass1_runnum.list");
    std::ifstream v5("ld_pass1_charge.list");
    std::ifstream v6("ld_pass1_fn.list");
    
    int dc_region = 6;
        
        std::vector<float> sam_v;
        std::vector<float> p_v;
        std::vector<float> lv_v;
        std::vector<float> lw_v;
        
        
         std::string line;
   while(std::getline(v3, line)){
	std::stringstream ss(line);
        std::string hip;
        std::getline(ss, hip, ',');
        v3_a.push_back(hip);
        num++;
}



 int targ = 3;   //  0 = Cu ob;   1 = Sn ob;  2 = CxC  ob;   3 = LD2 ob&ib;   4 = Cu ib;   5 = Sn ib;   6 = CxC ib;
    
    int start_p = 13;
    int end_p = v3_a.size();
    
    float ed1_cut = 2.1;
    float ed2_cut = 2.4;
    float ed3_cut = 1.35;
    float ed4_cut = 1.5;
    float ed5_cut = 2.1;
    float ed6_cut = 2.1;
    
    float sc1l[4] = {0.126274, 0.0335092, -0.00458362, 0.000199021};
    float sc1u[4] = {0.295144, -0.00215462, -0.0000984186, -0.00000817584};







while(std::getline(v4, line)){
	std::stringstream ss(line);
        std::string hip;
        std::getline(ss, hip, ',');
        int nhip = std::stoi(hip);
        v4_a.push_back(nhip);
        num++;
}
while(std::getline(v5, line)){
	std::stringstream ss(line);
        std::string hip;
        std::getline(ss, hip, ',');
        //std::cout << hip << std::endl;
        double nhip = std::stod(hip);
      //  std::cout << nhip << std::endl;
        v5_a.push_back(nhip);
        num++;
}
while(std::getline(v6, line)){
	std::stringstream ss(line);
        std::string hip;
        std::getline(ss, hip, ',');
        int nhip = std::stoi(hip);
        v6_a.push_back(nhip);
        num++;
}


int e_tot_c = 0;
   //std::ofstream fw("actual_2_18318_inbend_yield.txt", std::ofstream::out);   //011270 fall & 011361 spring

   std::string dir_path[4] = { "/volatile/clas12/rg-d/production/lumi/v0_LD2_aicv/dcalign/recon/", "/volatile/clas12/rg-d/production/prod/v0_ob_CxC/dst/recon/", "/volatile/clas12/rg-d/production/pass0.2/mon/recon/", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/"};
//std::cout << "1 th" << std::endl;
   
   std::string run_num[4] = {"018321", "018451", "018490", "019128"}; //"011283"};
   int file_num[4] = {56, 200, 10, 260};//need to add the files for 18309       361 for 18451   200 for v0

   int first_num = 0;
   int second_num = 4;

   int tt = 0;
   
   int bad_e_run [50];
   int bad_pip_run [50];
   int bad_pim_run [50];
   int bad_r_run [50];
   int ber = 0;
   int bpipr = 0;
   int bpimr = 0;
   int brr = 0;
   
   
 int e_num = 0;
	   int e_pip_num = 0;
	   int e_pim_num = 0;
	   int e_pip_pim_num = 0; 
	   
	   int count_1 = 0;
      int count_2 = 0;
      int count_3 = 0;

   
      t_e_num = 0;
      t_pip_num = 0;
      t_pim_num = 0;
      t_rho_num = 0;
      t_bC = 0.0;

      e_yield = 0.0;
      pip_yield = 0.0;
      pim_yield = 0.0;
      rho_yield = 0.0;

      e_pf = 0.0;
      pip_pf = 0.0;
      pim_pf = 0.0;
      rho_pf = 0.0;
  // for(int redo = 0; redo < 4; redo++){
  
  double t_charge = 0.0;
  double nt_charge = 0.0;
  
                       
      for(int p = start_p; p < end_p; p++){ // for(int p = 0; p < v3_a.size(); p++){  cusn(34 37)33  cxc(18, 20)13  LD2(18,)
        std::string inFile;
            inFile =  v3_a[p];
         std::cout << v3_a[p] << std::endl;
       //  std::cout << v5_a[p] << std::endl;
       
   //int pid = 200;
   //float px = 200.0;
   //float py = 200.0;
   //float pz = 200.0;
   //float vx = 200.0;
   //float vy = 200.0;
   //float vz = 200.0;
   //float energy = 200.0;

   //int pip_pid = 200;
   //float pip_px = 200.0;
   //float pip_py = 200.0;                          //These values will be fill the root tree
   //float pip_pz = 200.0;
   //float pip_vx = 200.0;
   //float pip_vy = 200.0;                          
   //float pip_vz = 200.0;
   //float pip_energy = 200.0;
   //TH1F *tvz1_hist = new TH1F(" vz1", " vz", 200, -30, 20);
   //TH1F *tvz2_hist = new TH1F(" vz2", " vz", 200, -30, 20);
   //TH1F *tvz3_hist = new TH1F(" vz3", " vz", 200, -30, 20);
   //TH1F *tvz4_hist = new TH1F(" vz4", " vz", 200, -30, 20);
   //TH1F *tvz5_hist = new TH1F(" vz5", " vz", 200, -30, 20);
   //TH1F *tvz6_hist = new TH1F(" vz6", " vz", 200, -30, 20);
   
   int rr1_count = 0;
   int rr2_count = 0;
   int rr3_count = 0;
   int rr4_count = 0;
   int rr5_count = 0;
   int rr6_count = 0;
   
   // in GeV/c^2
   float m_pro = 0.938; //in Gev/c^2
   //float m_ele = 0.000511; // in GeV/c^2
   float m_pi = 0.1396;
   float m_rho = 0.770;
   
   int fe_count = 0;
   int fpip_count = 0;
   int fpim_count = 0;

   //int pim_pid = 200;
   //float pim_px = 200.0;
   //float pim_py = 200.0;
   //float pim_pz = 200.0;
   //float pim_vx = 200.0;
   //float pim_vy = 200.0;
   //float pim_vz = 200.0;
   //float pim_energy = 200.0;

   float rho_px = 200.0;
   float rho_py = 200.0;
   float rho_pz = 200.0;
   float rho_energy = 200.0;
   float Q2 = 200.0;
   float dQ2 = 200.0;
   float w = 200.0;
   float t = 200.0;
   float lc = 200.0;
   float Zh = 200.0;
   float rho_mass = 200.0;
   float beam_energy = 200.0;

	int r = 0;
	int i = 0;
	int c = 0;

   
      beam_energy = 10.56;
   

	
   std::vector<float> c_bC;
   std::vector<int> c_ev;
   std::vector<float> bC;
   std::vector<int> ev;

   

   std::vector<std::string> eve_bC;

   
      char* inputFile_char = new char[inFile.length() + 1];   // give inFile.c_str() in reader.open(
      strcpy(inputFile_char, inFile.c_str());
      const char* inputFile = inputFile_char;
   

      hipo::reader  reader;  								//hipo is a custom class defined in hipo4, it is used to define functions
      reader.open(inputFile);
      hipo::dictionary  factory; 								// This is stating that it is declaring a dictionary of class hipo named factory. 
      reader.readDictionary(factory);							// We are then using the hipo function reader to read this dictionary

      //factory.show();

      hipo::bank  particle(factory.getSchema("REC::Particle"));     //in class hipo particles is a type of bank, as well as detectors. the events correlate to the rows of the data  
      hipo::bank config(factory.getSchema("RUN::config"));
      hipo::bank events(factory.getSchema("REC::Event"));
      hipo::bank scaler(factory.getSchema("RUN::scaler"));
      hipo::bank hel(factory.getSchema("HEL::scaler"));
      hipo::bank hel_a(factory.getSchema("HEL::online"));
      hipo::bank track(factory.getSchema("REC::Track"));
      hipo::bank cher(factory.getSchema("REC::Cherenkov"));
      hipo::bank cal(factory.getSchema("REC::Calorimeter"));
      hipo::bank traj(factory.getSchema("REC::Traj"));

      hipo::event      event; 								//This declares a variable of class hipo named event

      //int counter_p = 
      int counter = 0;
      

      TLorentzVector e_lv;
      e_lv.SetPxPyPzE(0.0, 0.0, beam_energy, beam_energy);

      TLorentzVector p_lv;
      p_lv.SetPxPyPzE(0.0, 0.0, 0.0, m_pro);

      TLorentzVector te_lv;
      TLorentzVector pip_te_lv;
      TLorentzVector pim_te_lv;
      //te_lv.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

	  

	   int pid_e_count = 0;
      int pid_pip_count = 0;
      int pid_pim_count = 0;

      float bigger = 0.0;
      float diff_bC = 0.0;
      float total_bC = 0.0;
      
      float sec_dec = 200;
      float pip_sec_dec = 200;
      float pim_sec_dec = 200;
      
      float e_theta = 200;
      float pip_theta = 200;
      float pim_theta = 200;
      
      
      
      
      
      

      while(reader.next() == true){
         reader.read(event);
      
         event.getStructure(particle);
	 event.getStructure(events);
         event.getStructure(config);
         event.getStructure(scaler);
         event.getStructure(hel);
         event.getStructure(hel_a);
         event.getStructure(cher);
         event.getStructure(cal);
         event.getStructure(track);
         event.getStructure(traj);
         
       // particle.show();
        // events.show();
        //hel.show();
        //cher.show();
        //cal.show();
        //fmt.show();
        
        //track.show();
        //traj.show();
        
        int pid = 200;
   float px = 200.0;
   float py = 200.0;
   float pz = 200.0;
   float vx = 200.0;
   float vy = 200.0;
   float vz = 200.0;
   float energy = 200.0;

   int pip_pid = 200;
   float pip_px = 200.0;
   float pip_py = 200.0;                          //These values will be fill the root tree
   float pip_pz = 200.0;
   float pip_vx = 200.0;
   float pip_vy = 200.0;                          
   float pip_vz = 200.0;
   float pip_energy = 200.0;
   
   // in GeV/c^2
   //float m_pro = 0.938; //in Gev/c^2
   //float m_ele = 0.000511; // in GeV/c^2
   //float m_pi = 0.1396;
   //float m_rho = 0.770;

   int pim_pid = 200;
   float pim_px = 200.0;
   float pim_py = 200.0;
   float pim_pz = 200.0;
   float pim_vx = 200.0;
   float pim_vy = 200.0;
   float pim_vz = 200.0;
   float pim_energy = 200.0;

   //float rho_px = 200.0;
   //float rho_py = 200.0;
   //float rho_pz = 200.0;
   //float rho_energy = 200.0;
   //float Q2 = 200.0;
   //float w = 200.0;
   //float t = 200.0;
   //float lc = 200.0;
   //float Zh = 200.0;
   //float rho_mass = 200.0;
   //float beam_energy = 200.0;

         int nrows = 0; 
         int a_pid[200] = {200};
         float a_px[200] = {200}; 
         float a_py[200] = {200};
         float a_pz[200] = {200};
         float a_vx[200] = {200};
         float a_vy[200] = {200};
         float a_vz[200] = {200};
         float a_energy[200] = {200};
         float a_pin[30] = {30};
         float a_sect[30] = {30};
         float ac_px[30] = {200};
         float ac_py[30] = {200};
         float ac_pz[30] = {200};
         float ac_energy[30] = {200};

         float a_rho_mass[200] = {200};
         float a_rho_px[200] = {200};
         float a_rho_py[200] = {200};
         float a_rho_pz[200] = {200};
         float a_rho_energy[200] = {200};

         int a_pip_pid[200] = {200};
         float a_pip_px[200] = {200};
         float a_pip_py[200] = {200};
         float a_pip_pz[200] = {200};
         float a_pip_vx[200] = {200};
         float a_pip_vy[200] = {200};
         float a_pip_vz[200] = {200};
         float a_pip_energy[200] = {200};
         float a_pip_pin[30] = {30};
         float a_pip_sect[30] = {30};
         float ac_pip_px[30] = {200};
         float ac_pip_py[30] = {200};
         float ac_pip_pz[30] = {200};
         float ac_pip_energy[30] = {200};

         int a_pim_pid[200] = {200};                                                     //This section declares the arrays that saves the data of one event. Each will be rewritten for each event.
         float a_pim_px[200] = {200};
         float a_pim_py[200] = {200};
         float a_pim_pz[200] = {200};
         float a_pim_vx[200] = {200};
         float a_pim_vy[200] = {200};
         float a_pim_vz[200] = {200};
         float a_pim_energy[200] = {200};
         float a_pim_pin[30] = {30};
         float a_pim_sect[30] = {30};
         
         int b_pip_pid[30] = {300};
         float b_pip_px[30] = {300};
         float b_pip_py[30] = {300};
         float b_pip_pz[30] = {300};
         float b_pip_vx[30] = {300};
         float b_pip_vy[30] = {300};
         float b_pip_vz[30] = {300};
         float b_pip_energy[30] = {300};
         float b_pip_pin[30] = {300};
         float b_pip_sect[30] = {300};
         float bc_pip_px[30] = {300};
         float bc_pip_py[30] = {300};
         float bc_pip_pz[30] = {300};
         float bc_pip_energy[30] = {300};

         int b_pim_pid[30] = {30};                                                     //This section declares the arrays that saves the data of one event. Each will be rewritten for each event.
         float b_pim_px[30] = {30};
         float b_pim_py[30] = {30};
         float b_pim_pz[30] = {30};
         float b_pim_vx[30] = {30};
         float b_pim_vy[30] = {30};
         float b_pim_vz[30] = {30};
         float b_pim_energy[30] = {30};
         float b_pim_pin[30] = {30};
         float b_pim_sect[30] = {30};
         
         float ac_pim_px[30] = {200};
         float ac_pim_py[30] = {200};
         float ac_pim_pz[30] = {200};
         float ac_pim_energy[30] = {200};

   	 float ivx_sigma = 0.054;
   	 float ivy_sigma = 0.163737;
   	 float iivx_sigma = 0.0577;
   	 float iivy_sigma = 0.162782;
   	 float iiivx_sigma = 0.05309;
   	 float iiivy_sigma = 0.166;
   	 float ivvx_sigma = 0.055688;
   	 float ivvy_sigma = 0.1614;
   	 float vvx_sigma = 0.0544952;
   	 float vvy_sigma = 0.1658;
   	 float vivx_sigma = 0.062778;
   	 float vivy_sigma = 0.16245;
   	 
      
        std::vector<int> a_row;

         nrows = particle.getRows(); 
          pid_e_count = 0;
          pid_pip_count = 0;
          pid_pim_count = 0;
         int address = -1;	                   //the rows start at 0 and increase by into the positive. Initiation at -1 was choosen to show that no row was asigned

         float e_final = 0.0;

         TLorentzVector pip_lv;// = new TLorentzVector(0., 0., 0., 0.);                //TLorentzVector(px,py,pz,e)
         TLorentzVector pim_lv;// = new TLorentzVector(0., 0., 0., 0.);            //This declares the lorentze four vector for pions and negative pions
        
         TLorentzVector rho_lv;
         TLorentzVector vpho_lv;
         TLorentzVector t_lv;
         TLorentzVector el_lv;
         TLorentzVector w_lv;
         TLorentzVector del_lv;
         TLorentzVector dvpho_lv;
         
         float e_total = 0.0;
         float e1_total = 0.0;
         float e2_total = 0.0;
         float e3_total = 0.0;
         float e4_total = 0.0;
         float e5_total = 0.0;
         float e6_total = 0.0;
         float ppcal = 0.0;
         float pcal = 0.0;
         float nphe = 0.0;
    

         
         
         //outbending
     /*
         double ez_lcut = -10.0;  // outbending
         double ez_ucut = 5;
         double pipz_lcut = -10.0;  
         double pipz_ucut = 5;
         double pimz_lcut = -10.0;  
         double pimz_ucut = 5;
         */
      /*
         double ez_lcut = -9.6;  // Inbending
         double ez_ucut = 5;
         double pipz_lcut = -9.6;  
         double pipz_ucut = 5;
         double pimz_lcut = -9.6;  
         double pimz_ucut = 5;
        */ 
         
       //  double ez_lcut = -15.0;  // LD2
       //  double ez_ucut = 5;
       //  double pipz_lcut = -15.0;  
       //  double pipz_ucut = 5;
       //  double pimz_lcut = -15.0;  
       //  double pimz_ucut = 5;
     
     
     
         //double ez_lcut = -10.49509;  //for Copper
         //double ez_ucut = -6.50215;
         //double ez_lcut = -6.01642;  //for Tin
         //double ez_ucut = 5;
         //double ez_lcut = -15;  //for Liquid Deuterium
         //double ez_ucut = 5;
         //double ez_lcut = -10.52183; //or Carbon
         //double ez_ucut = 5; 
         //double ez_lcut = -10000000000;  //for no cuts
         //double ez_ucut =  10000000000;
         
         //double pipz_lcut = -11.03832;  //for Copper
         //double pipz_ucut = -5.85492;
         //double pipz_lcut = -4.74283;  //for Tin
         //double pipz_ucut = 5;
         //double pipz_lcut = -15;  //for Liquid Deuterium
         //double pipz_ucut = 5;
         //double pipz_lcut = -9.67505;  //for Carbon
         //double pipz_ucut = 5;
         //double pipz_lcut = -10000000000;  //for no cuts
         //double pipz_ucut =  10000000000;
         
         //double pimz_lcut = -11.20497;  //for Copper
         //double pimz_ucut = -6.29121;
         //double pimz_lcut = -5.74827;  //for Tin
         //double pimz_ucut = 5;
         //double pimz_lcut = -15;  //for Liquid Deuterium
         //double pimz_ucut = 5;
         //double pimz_lcut = -10.55386;  //for Carbon
         //double pimz_ucut = 5;
         //double pimz_lcut = -10000000000;  //for no cuts
         //double pimz_ucut =  10000000000;

//inbending

         //double ez_lcut = -9.20876;  //for Copper
         //double ez_ucut = -5.75756;
         //double ez_lcut = -4.08899;  //for Tin
         //double ez_ucut = 5;
         //double ez_lcut = -15;  //for Liquid Deuterium
         //double ez_ucut = 5;
         //double ez_lcut = -9.10818;  //for Carbon
         //double ez_ucut = 5;
         //double ez_lcut = -10000000000;  //for no cuts
         //double ez_ucut =  10000000000;
         
         //double pipz_lcut = -9.95123;  //for Copper
         //double pipz_ucut = -6.15124;
         //double pipz_lcut = -5.26301;  //for Tin
         //double pipz_ucut = 5;
         //double pipz_lcut = -15;  //for Liquid Deuterium
         //double pipz_ucut = 5;
         //double pipz_lcut = -10.42521;  //for Carbon
         //double pipz_ucut = 5;
         //double pipz_lcut = -10000000000;  //for no cuts
         //double pipz_ucut =  10000000000;
         
         //double pimz_lcut = -9.95123;  //for Copper
         //double pimz_ucut = -5.26301;
         //double pimz_lcut = -4.86721;  //for Tin
         //double pimz_ucut = 5;
         //double pimz_lcut = -15;  //for Liquid Deuterium
         //double pimz_ucut = 5;
         //double pimz_lcut = -9.67362;  //for Carbon
         //double pimz_ucut = 5;
         //double pimz_lcut = -10000000000;  //for no cuts
         //double pimz_ucut =  10000000000;
         
         double ez_lcut[7] = {-10.49509, -6.01642, -10.52183, -15, -9.20876, -4.08899, -9.10818};
         double ez_ucut[7] = {-6.50215, 5, 5, 5, -5.75756, 5, 5};
         
         double pipz_lcut[7] = {-11.03832, -4.74283, -9.67505, -15, -9.95123, -5.26301, -10.42521};
         double pipz_ucut[7] = {-5.85492, 5, 5, 5, -6.15124, 5, 5};
         
         double pimz_lcut[7] = {-11.20497, -5.74827, -10.55386, -15, -9.95123, -4.86721, -9.67362};
         double pimz_ucut[7] = {-6.29121, 5, 5, 5, -5.26301, 5, 5};
         
         
       
       
         //std::cout << nrows << std::endl;
         for(int row = 0; row <= nrows; row++){
     //    std::cout << row << " | " << nrows << " | " << ((p + 1.0) / v3_a.size()) * 100 << "%" <<std::endl;
            float chi2 = 0.0;
            int temp_s = 0;
            float temp_vz = 0.0;
            temp_vz = particle.getFloat("vz", row);
            temp_s = particle.getInt("status", row);
           // if(abs(temp_s)/2000 == 1){
            if(particle.getInt("pid",row) == 11){
               
		        // float temp_vz = 0.0;
    		      //int temp_s = 0;
    		      
        

                  //   temp_vz = particle.getFloat("vz", row);
		     //temp_s = particle.getInt("status", row);
		     chi2 = particle.getFloat("chi2pid",row);
                         count_1++;
                         
                        // if(abs(temp_s)/2000 == 1 && temp_s < 0){
                          echi_hist->Fill(chi2);
                 //    if(temp_vz >= -11.5 && temp_vz <= 5){
                        //std::cout << nrows << std::endl;
			               
                        if (abs(temp_s)/2000 == 1 && temp_s < 0 && temp_vz >= ez_lcut[targ] && temp_vz <= ez_ucut[targ]){
			if(chi2 > -5 && chi2 < 5){	   
				a_row.push_back(row);
				count_2++;
           	      		   a_pid[pid_e_count] = 11;
           	      		   a_px[pid_e_count] = particle.getFloat("px",row);
            	   		   a_py[pid_e_count] = particle.getFloat("py",row);
              	   		   a_pz[pid_e_count] = particle.getFloat("pz",row);
                  		   a_vx[pid_e_count] = particle.getFloat("vx",row);
           	      		   a_vy[pid_e_count] = particle.getFloat("vy",row);
           	      		   a_vz[pid_e_count] = temp_vz;
                  		   a_energy[pid_e_count] = sqrt(pow(a_px[pid_e_count],2) + pow(a_py[pid_e_count],2) + pow(a_pz[pid_e_count],2));
                  		   a_pin[pid_e_count] = row;
                  		   ep_hist->Fill(sqrt(pow(a_px[pid_e_count],2) + pow(a_py[pid_e_count],2) + pow(a_pz[pid_e_count],2)));
                  		//   std::cout << fmt.getInt("layer", row);
                       		   pid_e_count++;
					//te_lv.SetPxPyPzE(a_px[pid_e_count],a_py[pid_e_count],a_pz[pid_e_count],a_energy[pid_e_count]);
					ixy_hist->Fill(px, py);
					e_tot_c++;
                                   e_count++;
         // }
	               }
               }
            } else if(particle.getInt("pid",row) == 211){
            chi2 = particle.getFloat("chi2pid",row);
            pipchi_hist->Fill(chi2);
            if(chi2 > -10 && chi2 < 10 && temp_vz >= pipz_lcut[targ] && temp_vz <= pipz_ucut[targ] && abs(temp_s)/2000 == 1){	
              	a_pip_pid[pid_pip_count] = 211;
              	a_pip_px[pid_pip_count] = particle.getFloat("px",row);
           	   a_pip_py[pid_pip_count] = particle.getFloat("py",row);
              	a_pip_pz[pid_pip_count] = particle.getFloat("pz",row);
               a_pip_vx[pid_pip_count] = particle.getFloat("vx",row);
              	a_pip_vy[pid_pip_count] = particle.getFloat("vy",row);
           	   a_pip_vz[pid_pip_count] = particle.getFloat("vz",row);
               a_pip_energy[pid_pip_count] = sqrt(pow(a_pip_px[pid_pip_count],2) + pow(a_pip_py[pid_pip_count],2) + pow(a_pip_pz[pid_pip_count],2) + pow(m_pi,2));
              // pipp_hist->Fill(sqrt(pow(a_pip_px[pid_pip_count],2) + pow(a_pip_py[pid_pip_count],2) + pow(a_pip_pz[pid_pip_count],2)));
               a_pip_pin[pid_pip_count] = row;
               pid_pip_count++;
               pip_count++;
               //std::cout << a_pip_vx[150] << " | " << a_pip_vx[151] << " | " << a_pip_vx[152] << " | " << a_pip_vx[153] << " | " << a_pip_vx[154] << " | " << std::endl;
               
               if (particle.getFloat("vx",row) == 0){
                   zerosp = zerosp + 1;
               }
           }
            } else if(particle.getInt("pid",row) == -211){
            chi2 = particle.getFloat("chi2pid",row);
            pimchi_hist->Fill(chi2);
            if(chi2 > -10 && chi2 < 10 && temp_vz >= pimz_lcut[targ] && temp_vz <= pimz_ucut[targ] && abs(temp_s)/2000 == 1){	
              	a_pim_pid[pid_pim_count] = -211;
              	a_pim_px[pid_pim_count] = particle.getFloat("px",row);
              	a_pim_py[pid_pim_count] = particle.getFloat("py",row);
           	   a_pim_pz[pid_pim_count] = particle.getFloat("pz",row);
               a_pim_vx[pid_pim_count] = particle.getFloat("vx",row);
              	a_pim_vy[pid_pim_count] = particle.getFloat("vy",row);
              	a_pim_vz[pid_pim_count] = particle.getFloat("vz",row);
               a_pim_energy[pid_pim_count] = sqrt(pow(a_pim_px[pid_pim_count],2) + pow(a_pim_py[pid_pim_count],2) + pow(a_pim_pz[pid_pim_count],2) + pow(m_pi,2));
            //   pimp_hist->Fill(sqrt(pow(a_pim_px[pid_pim_count],2) + pow(a_pim_py[pid_pim_count],2) + pow(a_pim_pz[pid_pim_count],2)));
               a_pim_pin[pid_pim_count] = row;
               pid_pim_count++;
               pim_count++;
               }
            }
           // }
         }

	int psize = a_row.size();
//ietot_p->AddPoint(x,y);
        int sect = 0;
      
	      int w_cut = 0;
         int z_cut = 0;
         int t_cut = 0;
         
         if (pid_e_count > 0){
         e_num++;
         }
         
          if (pid_e_count > 0 && pid_pip_count > 0){
         e_pip_num++;
         }
         
          if (pid_e_count > 0 && pid_pim_count > 0){
         e_pim_num++;
         }
         
          if (pid_e_count > 0 && pid_pip_count > 0 && pid_pim_count > 0){
         e_pip_pim_num++;
         }
         
         
      // std::cout << pid_pip_count << std::endl;
      
      
    //     float tvz = a_vz[0];
        //    float tuez = (-8.477 + 3*0.675);      //Cu for 34-37
        //    float tlez = (-8.477 - 3*0.675);
//float tuez = (-3.620 + 3*0.795);        //Sn for 34-37
  //          float tlez = (-3.620 - 3*0.795);
        //    float tuez =  500000000; //for no vz cut
        //    float tlez = -100000000;
          //  float tuez =  5; //for cxc vz cut
          //  float tlez = (-13.688906 + 3*0.738610);
       //     if(vz >= lez && vz <= uez){
      
      
      
         if( pid_e_count > 0 ){               //This ensures that only the events with all three particle will be looked at

            int ne = pid_e_count;
            int npip = pid_pip_count;
            int npim = pid_pim_count;                                                 //Since the arrays are filled from the first position until all of that particle is in
            int trows = -1;
            int frows = -1;
            
            trows = track.getRows(); 
           int ctrows = cal.getRows(); 
            for(int ec = 0; ec < ne; ec++){
                for(int row = 0; row < ctrows; row++){
                    if(cal.getInt("pindex",row) == a_pin[ec]){
                        a_sect[ec] = cal.getInt("sector",row);
                        
         	        
                       // std::cout << a_sect[ec] << std::endl;
                     //  ac_px[ec] = cher.getFloat("px", row);
                   /*
                     for(int frow = 0; frow < frows; frow++){
                     if(fmt.getInt("index",frow) == row){
                     
                     float fmtt = 100000;
                     fmtt = fmt.getFloat("Vtx0_z",frow);
                     
                     if(fmtt != 0){
                     fmt_vz_hist->Fill(fmtt);
                     if (a_sect[ec] == 1){
                         fmt_vz1_hist->Fill(fmtt);
                     }else if(a_sect[ec] ==2){
                         fmt_vz2_hist->Fill(fmtt);
                     }else if(a_sect[ec] ==3){
                         fmt_vz3_hist->Fill(fmtt);
                     }else if(a_sect[ec] ==4){
                         fmt_vz4_hist->Fill(fmtt);
                     }else if(a_sect[ec] ==5){
                         fmt_vz5_hist->Fill(fmtt);
                     }else if(a_sect[ec] ==6){
                         fmt_vz6_hist->Fill(fmtt);
                     
                     }
                     }
                     }
                     }
                     */
                    }
                }
            }
            
            for(int pipc = 0; pipc < npip; pipc++){
                for(int row = 0; row < nrows; row++){
                    if(cal.getInt("pindex",row) == a_pip_pin[pipc]){
                        a_pip_sect[pipc] = cal.getInt("sector",row);
                        //std::cout << a_pip_sect[pipc] << std::endl;
                    }
                }
            }
            
            for(int pimc = 0; pimc < npim; pimc++){
                for(int row = 0; row < nrows; row++){
                    if(cal.getInt("pindex",row) == a_pim_pin[pimc]){
                        a_pim_sect[pimc] = cal.getInt("sector",row);
                    }
                }
            }
            
           // 	int psize = a_row.size();
//ietot_p->AddPoint(x,y);
       // int sect = 0;

int carows = cal.getRows(); 
int trkrows = track.getRows();
int trjrows = traj.getRows();
float lv_cal = 0.0;
float lw_cal = 0.0;
float edge = 0.0;


double xtr = 0.0;
double ytr = 0.0;
double ztr = 0.0;
			
			
double angle = 0.0;
double s = 0.0;
double c = 0.0;
double x0 = 0.0;
double y0 = 0.0;
			
int track_sect = 0;	
			
double the_n = 0.0;
//int trkc = 0;
/*
for(int ec = 0; ec < ne; ec++){
                for(int row = 0; row < trkrows; row++){
                    if(track.getInt("pindex",row) == a_pin[ec]){
                    
			
			track_sect = track.getInt("sector", row);
		//trkc++;	
			 
			}
		}
           } 
 
*/

edge = 0.0;
float chindf = 0.0;
            
            for(int ec = 0; ec < ne; ec++){
                for(int row = 0; row < trjrows; row++){
                    if(traj.getInt("pindex",row) == a_pin[ec] && traj.getInt("detector", row) == 6 && traj.getInt("layer", row) == dc_region){
                    
                    edge = 0.0;
                 edge = traj.getFloat("edge", row);
                 
                    for(int tkrow = 0; tkrow < trkrows; tkrow++){
                    	if(track.getInt("pindex",tkrow) == a_pin[ec] ){
                    
			 
				track_sect = track.getInt("sector", tkrow);   //use Tprofile function
				
				
				//if(track.getInt("NDF", tkrow) != 0.0){
				chindf = (track.getFloat("chi2", tkrow) * 1.0) / (track.getInt("NDF", tkrow) * 1.0) ;
				
				//std::cout << edge << " | " << track_sect << " | " << chindf << std::endl;
				
		         if(track_sect == 1){
			 	chin_hist1->Fill(edge, chindf);
			 } else if(track_sect == 2){
			 	chin_hist2->Fill(edge, chindf);
			 } else if(track_sect == 3){
			 	chin_hist3->Fill(edge, chindf);
			 } else if(track_sect == 4){
			 	chin_hist4->Fill(edge, chindf);
			 } else if(track_sect == 5){
			 	chin_hist5->Fill(edge, chindf);
			 } else if(track_sect == 6){
			 	chin_hist6->Fill(edge, chindf);
			 }
		         //}
			
                    
                    
                    xtr = 1110.0;
                    ytr = 1110.0;
                    ztr = 1110.0;
                    x0 = 1110.0;
                    y0 = 1110.0;
                    s = 0.0;
                    c = 0.0;
                    angle = 0.0;
                    the_n = 0.0;
                    
                    
                 xtr = traj.getFloat("x", row);
			 ytr = traj.getFloat("y", row);
			 ztr = traj.getFloat("z", row);
			
			
			 angle = (track_sect - 1.0) * M_PI / 3.0;
			 s = std::sin( angle);
			 c = std::cos( angle);
			 x0 = xtr;
			 y0 = ytr;
			
			xtr = c * x0 - s * y0;
			ytr = s * x0 + c * y0;
			
			 the_n = std::atan2( std::sqrt(xtr*xtr + ytr*ytr), ztr ) * 180.0 / M_PI;
			 
			 thve_hist->Fill(edge, the_n);
			 
			 
			 
			 if(track_sect == 1){
			 	xy_uncut1->Fill(x0, y0);
			 	if( edge >= ed1_cut){
			 
			     		xy_cut1->Fill(x0, y0);
			 
			 	}
			 } else if(track_sect == 2){
			 	xy_uncut2->Fill(x0, y0);
			 	if( edge >= ed2_cut){
			 
			     		xy_cut2->Fill(x0, y0);
			 
			 	}
			 } else if(track_sect == 3){
			 	xy_uncut3->Fill(x0, y0);
			 	if( edge >= ed3_cut){
			 
			     		xy_cut3->Fill(x0, y0);
			 
			 	}
			 } else if(track_sect == 4){
			 	xy_uncut4->Fill(x0, y0);
			 	if( edge >= ed4_cut){
			 
			     		xy_cut4->Fill(x0, y0);
			 
			 	}
			 } else if(track_sect == 5){
			 	xy_uncut5->Fill(x0, y0);
			 	if( edge >= ed5_cut){
			 
			     		xy_cut5->Fill(x0, y0);
			 
			 	}
			 } else if(track_sect == 6){
			 	xy_uncut6->Fill(x0, y0);
			 	if( edge >= ed6_cut){
			 
			     		xy_cut6->Fill(x0, y0);
			 
			 	}
			 }
			 
			 }
		    }
                    
                    
			 
			 
			 //}
			}
		   }
               } 
               
               
            for(int ec = 0; ec < npim; ec++){
                for(int row = 0; row < trjrows; row++){
                    if(traj.getInt("pindex",row) == a_pim_pin[ec] && traj.getInt("detector", row) == 6 && traj.getInt("layer", row) == dc_region){
                    //std::cout << "tracks e = " << traj.getInt("layer", row) << std::endl;
                 //   if(a_pip_sect[i] == 1 || a_pip_sect[i] == 2 || a_pip_sect[i] == 3 || a_pip_sect[i] == 4 || a_pip_sect[i] == 5 || a_pip_sect[i] == 6){
                 
                 
                    for(int tkrow = 0; tkrow < trkrows; tkrow++){
                    	if(track.getInt("pindex",tkrow) == traj.getInt("index",row) ){
                    
			
			}
		    }
                    
                    edge = 0.0;
                   
                    edge = traj.getFloat("edge", row);
                    
                    xtr = 1110.0;
                    ytr = 1110.0;
                    ztr = 1110.0;
                    x0 = 1110.0;
                    y0 = 1110.0;
                    s = 0.0;
                    c = 0.0;
                    angle = 0.0;
                    the_n = 0.0;
                    
                    
                 xtr = traj.getFloat("x", row);
			 ytr = traj.getFloat("y", row);
			 ztr = traj.getFloat("z", row);
			
			
			 angle = (track_sect - 1.0) * M_PI / 3.0;
			 s = std::sin( angle);
			 c = std::cos( angle);
			 x0 = xtr;
			 y0 = ytr;
			
			xtr = c * x0 - s * y0;
			ytr = s * x0 + c * y0;
			
			 the_n = std::atan2( std::sqrt(xtr*xtr + ytr*ytr), ztr ) * 180.0 / M_PI;
			 
			 
			 thve_hist->Fill(edge, the_n);
			 
			 
			// xy_uncut->Fill(x0, y0);
			 if( edge >= 3.71){
			 
			     //xy_cut->Fill(x0, y0);
			 
			 }
			 
			 
			 //}
			}
		   }
               } 
                
                float e_inner = 0.0;
                float e_outer = 0.0;
                for(int ec = 0; ec < ne; ec++){
                for(int row = 0; row < carows; row++){
                    if(cal.getInt("pindex",row) == a_pin[ec]){
                    
        		float pcal = cal.getFloat("energy", row);
        		lv_cal = lv_cal + cal.getFloat("lv", row);   //only layer one
        		lw_cal = lw_cal + cal.getFloat("lw", row);
         	        e_total = e_total + pcal;
         	        
         	        

                        sect = cal.getInt("sector",row);
				
			if( cal.getInt("layer", row) == 1){
				ppcal = pcal;
			}
			if( cal.getInt("layer", row) == 4){
				e_inner = pcal;
			}
			if( cal.getInt("layer", row) == 7){
				e_outer = pcal;
			}
			
			if (sect == 1){
			    e1_total = e1_total + pcal;
			}else if (sect == 2){
			    e2_total = e2_total + pcal;
			}else if (sect == 3){
			    e3_total = e3_total + pcal;
			}else if (sect == 4){
			    e4_total = e4_total + pcal;
			}else if (sect == 5){
			    e5_total = e5_total + pcal;
			}else if (sect == 6){
			    e6_total = e6_total + pcal;
			}
			
			

                    }
               }
           }
           
           TLorentzVector en_lv;
           en_lv.SetPxPyPzE(a_px[0],a_py[0],a_pz[0],a_energy[0]);
           float en_theta = en_lv.Theta() * 180.0/M_PI;
           
           for(int ec = 0; ec < ne; ec++){
           	//if(ppcal != 0){
           	//if((e_total/a_energy[ec]) <= (0.0137 * 3) + 0.24  && (e_total/a_energy[ec]) >= (-0.0137 * 3) + 0.24){
		sam_hist->Fill(e_total/a_energy[ec]);
		
        		 samp_hist->Fill(a_energy[ec], e_total/a_energy[ec]);
        		 samp1_hist->Fill(a_energy[ec], e1_total/a_energy[ec]);
        		 samp2_hist->Fill(a_energy[ec], e2_total/a_energy[ec]);
        		 samp3_hist->Fill(a_energy[ec], e3_total/a_energy[ec]);
        		 samp4_hist->Fill(a_energy[ec], e4_total/a_energy[ec]);
        		 samp5_hist->Fill(a_energy[ec], e5_total/a_energy[ec]);
        		 samp6_hist->Fill(a_energy[ec], e6_total/a_energy[ec]);
        		 
        		 if ( (e1_total/a_energy[ec]) >= (sc1l[0] + sc1l[1]*a_energy[ec] + sc1l[2]*pow(a_energy[ec],2) + sc1l[3]*pow(a_energy[ec], 3) )  && (e1_total/a_energy[ec]) <= (sc1u[0] + sc1u[1]*a_energy[ec] + sc1u[2]*pow(a_energy[ec],2) + sc1u[3]*pow(a_energy[ec], 3) )){
        		       samp1_cut_hist->Fill(a_energy[ec], e1_total/a_energy[ec]);
        		 }
        		 
        		 samv_hist->Fill(lv_cal, e_total/a_energy[ec]);
        		 samw_hist->Fill(lw_cal, e_total/a_energy[ec]);
        		 
        		 
        		 
        		 sam_v.push_back(e_total/a_energy[ec]);
        		 p_v.push_back(a_energy[ec]);
        		 lv_v.push_back(lv_cal);
        		 lw_v.push_back(lw_cal);
        	//}	 
        		 
			// wnphe_hist->Fill(nphe);
       		        //   wpcal_hist->Fill(ppcal);
                       //if((e_total/a_energy[ec]) <= (0.0137 * 3) + 0.24  && (e_total/a_energy[ec]) >= (-0.0137 * 3) + 0.24){
			//	iixy_hist->Fill(a_energy[ec], e_total/a_energy[ec]);
			//}

                
         //}
         }
         
            
            for(int i = 0; i < ne; i++){
               if(e_final < a_energy[i]){                                 //This saves the highest energy electron
                  e_final = a_energy[i];
                  address = i;
               }
            } 
            
            
         
            if(address > -1 ){
               pid = a_pid[address];
               px = a_px[address];
               py = a_py[address];
               pz = a_pz[address];
               vx = a_vx[address];
               vy = a_vy[address];
               vz = a_vz[address];
               energy = e_final;
               
             //  ep_hist->Fill(sqrt(pow(a_px[pid_e_count],2) + pow(a_py[pid_e_count],2) + pow(a_pz[pid_e_count],2)));
               e_hist->Fill(pid);
               xy_hist->Fill(px,py);
               te_lv.SetPxPyPzE(px,py,pz,energy);
			count_3++;		
                                   sec_dec =  te_lv.Phi() * 180/M_PI;
                                   e_theta = te_lv.Theta() * 180/M_PI;
                                   e_pht_hist->Fill(e_theta);
                                   
               del_lv.SetPxPyPzE(a_px[0],a_py[0],a_pz[0],a_energy[0]);
                                 
                                 e_phi->Fill(sec_dec);
                                 
                                // e_theta->Fill(e_theta);
                                 
               dvpho_lv = e_lv - del_lv;

               dQ2 = - dvpho_lv.Mag2(); 
             //  q_hist->Fill(dQ2);
             
             fe_count++;
                             if(dQ2 >= 1 && dQ2 <2){
         			bin1_count++;
                            }else if(dQ2 >= 2 && dQ2 <2.5){
                                bin2_count++;
                            }else if(dQ2 >= 2.5 && dQ2 <3){
    			        bin3_count++;
                           }else if(dQ2 >= 3 && dQ2 <3.5){
                                bin4_count++;
                           }else if(dQ2 >= 3.5 && dQ2 <4){
                                bin5_count++;
                           }else if(dQ2 >= 4 && dQ2 <6){
                                bin6_count++;
                           } 
                          
                                
                 float mul = 110000.0; // 3.0 4.0 5.0
                 float lxlim = 0.00656862 - mul*(0.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float uxlim = 0.00656862 + mul*(0.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float lylim = -0.1793 - mul*(0.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float uylim = -0.1793 + mul*(0.46594);  //1.5056; // 0.83372 1.16966 1.5056       
                 //      if( vx >= -6 && vx <= 6 && vy >= -6 && vy <= 6)  {   
                 
                               if( vx >= lxlim && vx <= uxlim )  {  
				vxx_hist->Fill(vx);
				//vx1_c++;
				}
				if(  vy >= lylim && vy <= uylim)  {  
   				vyy_hist->Fill(vy);
   			//	vy1_c++;
   				}
   				if( vx >= lxlim && vx <= uxlim && vy >= lylim && vy <= uylim && vz >= zl_b && vz <= zu_b)  {  
   				vzz_hist->Fill(vz);
   			//	vz1_c++;
				}
				
				
				             
                          tz_hist->Fill(sec_dec, vz);
				if(a_sect[0] == 1) {
				if( vx >= lxlim && vx <= uxlim )  {  
				ivx_hist->Fill(vx);
				//vx1_c++;
				}
				if(  vy >= lylim && vy <= uylim)  {  
   				ivy_hist->Fill(vy);
   			//	vy1_c++;
   				}
   				if( vx >= lxlim && vx <= uxlim && vy >= lylim && vy <= uylim && vz >= zl_b && vz <= zu_b)  {  
   				ivz_hist->Fill(vz);
   				//tvz1_hist->Fill(vz);
   			//	vz1_c++;
				}
				
                        } else if(a_sect[0] == 2) {
                                if( vx >= lxlim && vx <= uxlim )  {  
				iivx_hist->Fill(vx);
			//	vx2_c++;
				}
				if(  vy >= lylim && vy <= uylim)  {  
   				iivy_hist->Fill(vy);
   			//	vy2_c++;
   				}
				iixy_hist->Fill(px, py);
				if( vx >= lxlim && vx <= uxlim && vy >= lylim && vy <= uylim && vz >= zl_b && vz <= zu_b)  {  
				iivz_hist->Fill(vz);
				   				//tvz2_hist->Fill(vz);
			//	vz2_c++;
				}
				
                        } else if(a_sect[0] == 3) {
                                if( vx >= lxlim && vx <= uxlim )  {  
				iiivx_hist->Fill(vx);
			//	vx3_c++;
				}
				if(  vy >= lylim && vy <= uylim)  {  
   				iiivy_hist->Fill(vy);
   			//	vy3_c++;
   				}
				iiixy_hist->Fill(px, py);
				if( vx >= lxlim && vx <= uxlim && vy >= lylim && vy <= uylim && vz >= zl_b && vz <= zu_b)  {  
				iiivz_hist->Fill(vz);
				   				//tvz3_hist->Fill(vz);
			//	vz3_c++;
				}
				
				
                        } else if(a_sect[0] == 4) {
                                if( vx >= lxlim && vx <= uxlim )  {  
				ivvx_hist->Fill(vx);
			//	vx4_c++;
				//vx5_c++;
				}
				if(  vy >= lylim && vy <= uylim)  {  
   				ivvy_hist->Fill(vy);
   			//	vy4_c++;
   				}
				ivxy_hist->Fill(px, py);
				if( vx >= lxlim && vx <= uxlim && vy >= lylim && vy <= uylim && vz >= zl_b && vz <= zu_b)  {  
				ivvz_hist->Fill(vz);
				   				//tvz4_hist->Fill(vz);
			//	vz4_c++;
				}
				
                        } else if(a_sect[0] == 5) {
                                if( vx >= lxlim && vx <= uxlim )  {  
				vvx_hist->Fill(vx);
			//	vx5_c++;
				}
				if( vy >= lylim && vy <= uylim)  {  
   				vvy_hist->Fill(vy);
   			//	vy5_c++;
   				}
				vxy_hist->Fill(px, py);
				if( vx >= lxlim && vx <= uxlim && vy >= lylim && vy <= uylim && vz >= zl_b && vz <= zu_b)  {  
				vvz_hist->Fill(vz);
				   				//tvz5_hist->Fill(vz);
			//	vz5_c++;
				}
				
				
                        } else if(a_sect[0] == 6) {
                                if( vx >= lxlim && vx <= uxlim )  {  
				vivx_hist->Fill(vx);
			//	vx6_c++;
				}
				if( vy >= lylim && vy <= uylim)  {  
   				vivy_hist->Fill(vy);
   			//	vy6_c++;
   				}
				vixy_hist->Fill(px, py);
				if( vx >= lxlim && vx <= uxlim && vy >= lylim && vy <= uylim && vz >= zl_b && vz <= zu_b)  {  
				vivz_hist->Fill(vz);
				   				//tvz6_hist->Fill(vz);
			///	vz6_c++;
				}
				
				}
                       }
            
            
         int bp = 0;
            
            for(int i = 0; i < npip; i++){
            //if(a_pip_vz[i] >= pipz_lcut && a_pip_vz[0] <= pipz_ucut){
               pip_pid = a_pip_pid[i];
               pip_px = a_pip_px[i];
               pip_py = a_pip_py[i];
               pip_pz = a_pip_pz[i];
               pip_vx = a_pip_vx[i];
               pip_vy = a_pip_vy[i];
               pip_vz = a_pip_vz[i];
               pip_energy = a_pip_energy[i];
	       pip_hist->Fill(a_pip_sect[i]); 
	       
	       pipp_hist->Fill(sqrt(pow(a_pip_px[i],2) + pow(a_pip_py[i],2) + pow(a_pip_pz[i],2)));
	       
                              pip_te_lv.SetPxPyPzE(pip_px,pip_py,pip_pz,pip_energy);
			//count_3++;		
                                   pip_sec_dec =  pip_te_lv.Phi() * 180/M_PI;
                                   pip_theta = pip_te_lv.Theta() * 180/M_PI;
                                   pip_pht_hist->Fill(pip_sec_dec, pip_theta);
                                   
                                   pip_phi->Fill(pip_sec_dec);
                                  // pip_theta->Fill(pip_theta);
                                
                 float pip_mul = 110000.0; // 3.0 4.0 5.0
                 float pip_lxlim = 0.00656862 - pip_mul*(0.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float pip_uxlim = 0.00656862 + pip_mul*(0.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float pip_lylim = -0.1793 - pip_mul*(0.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float pip_uylim = -0.1793 + pip_mul*(0.46594);  //1.5056; // 0.83372 1.16966 1.5056       
                   //    if( pip_vx >= -6 && pip_vx <= 6 && pip_vy >= -6 && pip_vy <= 6 )  {                
                          //tz_hist->Fill(sec_dec, vz);
                          if(a_pip_sect[i] == 1 || a_pip_sect[i] == 2 || a_pip_sect[i] == 3 || a_pip_sect[i] == 4 || a_pip_sect[i] == 5 || a_pip_sect[i] == 6){
                          b_pip_pid[bp] = a_pip_pid[i];
                          b_pip_px[bp]  = a_pip_px[i];
            		  b_pip_py[bp]  = a_pip_py[i];
            		  b_pip_pz[bp]  = a_pip_pz[i];
           		  b_pip_vx[bp]  = a_pip_vx[i];
          		  b_pip_vy[bp]  = a_pip_vy[i];
          		  b_pip_vz[bp]  = a_pip_vz[i];
                          fpip_count++;
                          bp++;
                          }
		        if(pip_sec_dec > -30 && pip_sec_dec <= 30) {
				if( pip_vx >= pip_lxlim && pip_vx <= pip_uxlim )  {  
				pip_ivx_hist->Fill(pip_vx);
				//vx1_c++;
				}
				if(  pip_vy >= pip_lylim && pip_vy <= pip_uylim)  {  
   				pip_ivy_hist->Fill(pip_vy);
   				//vy1_c++;
   				}
   				if( pip_vx >= pip_lxlim && pip_vx <= pip_uxlim && pip_vy >= pip_lylim && pip_vy <= pip_uylim && pip_vz >= pip_zl_b && pip_vz <= pip_zu_b)  {  
   				pip_ivz_hist->Fill(pip_vz);
   				//vz1_c++;
				}
				
		         } else if(pip_sec_dec > 30 && pip_sec_dec <= 90) {
                                if( pip_vx >= pip_lxlim && pip_vx <= pip_uxlim )  {  
				pip_iivx_hist->Fill(pip_vx);
				//vx2_c++;
				}
				if(  pip_vy >= pip_lylim && pip_vy <= pip_uylim)  {  
   				pip_iivy_hist->Fill(pip_vy);
   				//vy2_c++;
   				}
				//iixy_hist->Fill(px, py);
				if( pip_vx >= pip_lxlim && pip_vx <= pip_uxlim && pip_vy >= pip_lylim && pip_vy <= pip_uylim && pip_vz >= pip_zl_b && pip_vz <= pip_zu_b)  {  
				pip_iivz_hist->Fill(pip_vz);
				//vz2_c++;
				}
				
			} else if(pip_sec_dec > 90 && pip_sec_dec <= 150) {
                                if( pip_vx >= pip_lxlim && pip_vx <= pip_uxlim )  {  
				pip_iiivx_hist->Fill(pip_vx);
				//vx3_c++;
				}
				if(  pip_vy >= pip_lylim && pip_vy <= pip_uylim)  {  
   				pip_iiivy_hist->Fill(pip_vy);
   				//vy3_c++;
   				}
				//iiixy_hist->Fill(px, py);
				if( pip_vx >= pip_lxlim && pip_vx <= pip_uxlim && pip_vy >= pip_lylim && pip_vy <= pip_uylim && pip_vz >= pip_zl_b && pip_vz <= pip_zu_b)  {  
				pip_iiivz_hist->Fill(pip_vz);
				//vz3_c++;
				}
				
			 } else if(pip_sec_dec > 150 || pip_sec_dec <= -150) {
                                if( pip_vx >= pip_lxlim && pip_vx <= pip_uxlim )  {  
				pip_ivvx_hist->Fill(pip_vx);
				//vx4_c++;
				//vx5_c++;
				}
				if(  pip_vy >= pip_lylim && pip_vy <= pip_uylim)  {  
   				pip_ivvy_hist->Fill(pip_vy);
   				//vy4_c++;
   				}
				//ivxy_hist->Fill(px, py);
				if( pip_vx >= pip_lxlim && pip_vx <= pip_uxlim && pip_vy >= pip_lylim && pip_vy <= pip_uylim && pip_vz >= pip_zl_b && pip_vz <= pip_zu_b)  {  
				pip_ivvz_hist->Fill(pip_vz);
				//vz4_c++;
				}
				
			 } else if(pip_sec_dec > -150 && pip_sec_dec <= -90) {
                                if( pip_vx >= pip_lxlim && pip_vx <= pip_uxlim )  {  
				pip_vvx_hist->Fill(pip_vx);
				//vx5_c++;
				}
				if( pip_vy >= pip_lylim && pip_vy <= pip_uylim)  {  
   				pip_vvy_hist->Fill(pip_vy);
   				//vy5_c++;
   				}
				//vxy_hist->Fill(px, py);
				if( pip_vx >= pip_lxlim && pip_vx <= pip_uxlim && pip_vy >= pip_lylim && pip_vy <= pip_uylim && pip_vz >= pip_zl_b && pip_vz <= pip_zu_b)  {  
				pip_vvz_hist->Fill(pip_vz);
				//vz5_c++;
				}
				
			 } else if(pip_sec_dec > -90 && pip_sec_dec <= -30) {
                                if( pip_vx >= pip_lxlim && pip_vx <= pip_uxlim )  {  
				pip_vivx_hist->Fill(pip_vx);
				//vx6_c++;
				}
				if( pip_vy >= pip_lylim && pip_vy <= pip_uylim)  {  
   				pip_vivy_hist->Fill(pip_vy);
   				//vy6_c++;
   				}
				//pim_vixy_hist->Fill(px, py);
				if( pip_vx >= pip_lxlim && pip_vx <= pip_uxlim && pip_vy >= pip_lylim && pip_vy <= pip_uylim && pip_vz >= pip_zl_b && pip_vz <= pip_zu_b)  {  
				pip_vivz_hist->Fill(pip_vz);
				//vz6_c++;
				}
				
            
                        }
                 // }
            }
       int bm = 0;
            for(int j = 0; j < npim; j++){ 
               pim_pid = a_pim_pid[j];
               pim_px = a_pim_px[j];
               pim_py = a_pim_py[j];
               pim_pz = a_pim_pz[j];
               pim_vx = a_pim_vx[j];
               pim_vy = a_pim_vy[j];
               pim_vz = a_pim_vz[j];
               pim_energy = a_pim_energy[j];
               pim_hist->Fill(pim_pid);
               
               //std::cout << pim_vx << std::endl;
               //std::cout << npim << std::endl;
                              pim_te_lv.SetPxPyPzE(pim_px,pim_py,pim_pz,energy);
			//count_3++;		
                                   pim_sec_dec =  pim_te_lv.Phi() * 180/M_PI;
                                   pim_theta = pim_te_lv.Theta() * 180/M_PI;
                                   pim_pht_hist->Fill(pim_sec_dec, pim_theta);
                                   
                                   pim_phi->Fill(pim_sec_dec);
                              //    std::cout << pim_sec_dec << std::endl; 
                                //   pim_theta->Fill(pim_theta);
                                   
                              pimp_hist->Fill(sqrt(pow(a_pim_px[j],2) + pow(a_pim_py[j],2) + pow(a_pim_pz[j],2)));     
                                
                 float pim_mul = 1100000.0; // 3.0 4.0 5.0
                 float pim_lxlim = 0.00656862 - pim_mul*(0.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float pim_uxlim = 0.00656862 + pim_mul*(0.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float pim_lylim = -0.1793 - pim_mul*(0.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float pim_uylim = -0.1793 + pim_mul*(0.46594);  //1.5056; // 0.83372 1.16966 1.5056       
                     //  if( pim_vx >= -6 && pim_vx <= 6 && pim_vy >= -6 && pim_vy <= 6)  {                
                          //tz_hist->Fill(sec_dec, vz);
                          if(a_pim_sect[j] == 1 || a_pim_sect[j] == 2 || a_pim_sect[j] == 3 || a_pim_sect[j] == 4 || a_pim_sect[j] == 5 || a_pim_sect[j] == 6){
                          b_pim_pid[bm] = a_pim_pid[j];
                          b_pim_px[bm]  = a_pim_px[j];
            		   b_pim_py[bm]  = a_pim_py[j];
            		   b_pim_pz[bm]  = a_pim_pz[j];
           		    b_pim_vx[bm]  = a_pim_vx[j];
          		     b_pim_vy[bm]  = a_pim_vy[j];
          		     b_pim_vz[bm]  = a_pim_vz[j];
          		     b_pim_energy[bm]  = a_pim_energy[j];
                          fpim_count++;
                          bm++;
                          }
                          
		        if(pim_sec_dec > -30 && pim_sec_dec <= 30) {
				if( pim_vx >= pim_lxlim && pim_vx <= pim_uxlim )  {  
				pim_ivx_hist->Fill(pim_vx);
				//vx1_c++;
				}
				if(  pim_vy >= pim_lylim && pim_vy <= pim_uylim)  {  
   				pim_ivy_hist->Fill(pim_vy);
   				//vy1_c++;
   				}
   				if( pim_vx >= pim_lxlim && pim_vx <= pim_uxlim && pim_vy >= pim_lylim && pim_vy <= pim_uylim && pim_vz >= pim_zl_b && pim_vz <= pim_zu_b)  {  
   				pim_ivz_hist->Fill(pim_vz);
   				//vz1_c++;
				}
				
		         } else if(pim_sec_dec > 30 && pim_sec_dec <= 90) {
                                if( pim_vx >= pim_lxlim && pim_vx <= pim_uxlim )  {  
				pim_iivx_hist->Fill(pim_vx);
				//vx2_c++;
				}
				if(  pim_vy >= pim_lylim && pim_vy <= pim_uylim)  {  
   				pim_iivy_hist->Fill(pim_vy);
   				//vy2_c++;
   				}
				//iixy_hist->Fill(px, py);
				if( pim_vx >= pim_lxlim && pim_vx <= pim_uxlim && pim_vy >= pim_lylim && pim_vy <= pim_uylim && pim_vz >= pim_zl_b && pim_vz <= pim_zu_b)  {  
				pim_iivz_hist->Fill(pim_vz);
				//vz2_c++;
				}
			
			} else if(pim_sec_dec > 90 && pim_sec_dec <= 150) {
                                if( pim_vx >= pim_lxlim && pim_vx <= pim_uxlim )  {  
				pim_iiivx_hist->Fill(pim_vx);
				//vx3_c++;
				}
				if(  pim_vy >= pim_lylim && pim_vy <= pim_uylim)  {  
   				pim_iiivy_hist->Fill(pim_vy);
   				//vy3_c++;
   				}
				//iiixy_hist->Fill(px, py);
				if( pim_vx >= pim_lxlim && pim_vx <= pim_uxlim && pim_vy >= pim_lylim && pim_vy <= pim_uylim && pim_vz >= pim_zl_b && pim_vz <= pim_zu_b)  {  
				pim_iiivz_hist->Fill(pim_vz);
				//vz3_c++;
				}
				
			 } else if(pim_sec_dec > 150 || pim_sec_dec <= -150) {
                                if( pim_vx >= pim_lxlim && pim_vx <= pim_uxlim )  {  
				pim_ivvx_hist->Fill(pim_vx);
				//vx4_c++;
				//vx5_c++;
				}
				if(  pim_vy >= pim_lylim && pim_vy <= pim_uylim)  {  
   				pim_ivvy_hist->Fill(pim_vy);
   				//vy4_c++;
   				}
				//ivxy_hist->Fill(px, py);
				if( pim_vx >= pim_lxlim && pim_vx <= pim_uxlim && pim_vy >= pim_lylim && pim_vy <= pim_uylim && pim_vz >= pim_zl_b && pim_vz <= pim_zu_b)  {  
				pim_ivvz_hist->Fill(pim_vz);
				//vz4_c++;
				}
				
			 } else if(pim_sec_dec > -150 && pim_sec_dec <= -90) {
                                if( pim_vx >= pim_lxlim && pim_vx <= pim_uxlim )  {  
				pim_vvx_hist->Fill(pim_vx);
				//vx5_c++;
				}
				if( pim_vy >= pim_lylim && pim_vy <= pim_uylim)  {  
   				pim_vvy_hist->Fill(pim_vy);
   				//vy5_c++;
   				}
				//vxy_hist->Fill(px, py);
				if( pim_vx >= pim_lxlim && pim_vx <= pim_uxlim && pim_vy >= pim_lylim && pim_vy <= pim_uylim && pim_vz >= pim_zl_b && pim_vz <= pim_zu_b)  {  
				pim_vvz_hist->Fill(pim_vz);
				//vz5_c++;
				}
				
			 } else if(pim_sec_dec > -90 && pim_sec_dec <= -30) {
                                if( pim_vx >= pim_lxlim && pim_vx <= pim_uxlim )  {  
				pim_vivx_hist->Fill(pim_vx);
				//vx6_c++;
				}
				if( pim_vy >= pim_lylim && pim_vy <= pim_uylim)  {  
   				pim_vivy_hist->Fill(pim_vy);
   				//vy6_c++;
   				}
				//pim_vixy_hist->Fill(px, py);
				if( pim_vx >= pim_lxlim && pim_vx <= pim_uxlim && pim_vy >= pim_lylim && pim_vy <= pim_uylim && pim_vz >= pim_zl_b && pim_vz <= pim_zu_b)  {  
				pim_vivz_hist->Fill(pim_vz);
				//vz6_c++;
				}
				
                       }
                 // }
            }
       //  }
            
            

            int rho_counter = 0;
          
            ROOT::Math::XYZVector ep_v;
            ep_v.SetX(a_px[0]);
            ep_v.SetY(a_py[0]);
            ep_v.SetZ(a_pz[0]);
            
            

           el_lv.SetPxPyPzE(a_px[0],a_py[0],a_pz[0],a_energy[0]);

               vpho_lv = e_lv - el_lv;

               Q2 = - vpho_lv.Mag2(); //((e_lv - el_lv) * (e_lv - el_lv));
               w_lv = vpho_lv + p_lv; 
               w = w_lv.Mag();
      // if(npim > 0 && npip > 0){
               w_hist->Fill(w);
        //   }
               ROOT::Math::XYZVector vpp_v;
            vpp_v.SetX(vpho_lv.Px());
            vpp_v.SetY(vpho_lv.Py());
            vpp_v.SetZ(vpho_lv.Pz());
               double x_bj = Q2 / (2 * m_pro * vpho_lv.E());
               double y_bj = vpho_lv.E() / beam_energy;
               
               xbj_hist->Fill(x_bj);
               ybj_hist->Fill(y_bj);

            for(int i = 0; i < npip; i++){                               //This double for loop will compare each pip to each pim in an event
                     ROOT::Math::XYZVector pip_v;
            pip_v.SetX(a_pip_px[i]);
            pip_v.SetY(a_pip_py[i]);
            pip_v.SetZ(a_pip_pz[i]);
            
 //           ROOT::Math::XYZVector pim_v;
  //          pim_v.SetX(a_pim_px[j]);
  //          pim_v.SetY(a_pim_py[j]);
   //         pim_v.SetZ(a_pim_pz[j]);
            
            double phi_h = angle_plane(vpp_v, ep_v, vpp_v, pip_v) * 180/M_PI;
        //    std::cout << "phi_h " << phi_h << "Q2 " << Q2 <<"nu " << vpho_lv.E() <<std::endl;
            phih_hist->Fill(phi_h);
            
               for(int j = 0; j < npim; j++){                            //This is primarly done using TLorentzsVector (the commented parts are my original attempt to calculate)
               
       //        ROOT::Math::XYZVector pip_v;
       //     pip_v.SetX(a_pip_px[i]);
       //     pip_v.SetY(a_pip_py[i]);
        //    pip_v.SetZ(a_pip_pz[i]);
            
            ROOT::Math::XYZVector pim_v;
            pim_v.SetX(a_pim_px[j]);
            pim_v.SetY(a_pim_py[j]);
            pim_v.SetZ(a_pim_pz[j]);
            
     //       double phi_h = angle_plane(ep_v, vpp_v, pip_v, pim_v);
     //       std::cout << "phi_h " << phi_h << "Q2 " << Q2 <<"nu " << vpho_lv.E() <<std::endl;
     //       phih_hist->Fill(phi_h);
            
            
                    if(a_pip_sect[i] == 1 || a_pip_sect[i] == 2 || a_pip_sect[i] == 3 || a_pip_sect[i] == 4 || a_pip_sect[i] == 5 || a_pip_sect[i] == 6 ){
                    if (a_pim_sect[j] == 1 || a_pim_sect[j] == 2 || a_pim_sect[j] == 3 || a_pim_sect[j] == 4 || a_pim_sect[j] == 5 || a_pim_sect[j] == 6){
                  pip_lv.SetPxPyPzE(a_pip_px[i], a_pip_py[i], a_pip_pz[i], a_pip_energy[i]);
                  pim_lv.SetPxPyPzE(a_pim_px[j], a_pim_py[j], a_pim_pz[j], a_pim_energy[j]);
                    
                  rho_lv = pip_lv + pim_lv;
                   // }
                    
                    
                  t_lv = (vpho_lv - rho_lv);
                 
                  
                  


                  
                 // if(0.6 <= rho_lv.Mag() && rho_lv.Mag() <= 1.0){   //should be 0.77 with +- 1.50[ first standard deviation from hist] error
                     
                     a_rho_mass[rho_counter] = rho_lv.Mag();
                     a_rho_px[rho_counter] = rho_lv.Px();
                     a_rho_py[rho_counter] = rho_lv.Py();
                     a_rho_pz[rho_counter] = rho_lv.Pz();
                     a_rho_energy[rho_counter] = rho_lv.E();
                     rho_counter++;
                  }
                  }
               } 
            }

            event_counter++;
         
            if(rho_counter > 0){
               for(int i = 0; i < rho_counter; i++){
               
               //sqrt(pow(a_pip_px[i],2) + pow(a_pip_py[i],2) + pow(a_pip_pz[i],2))

                  rho_lv.SetPxPyPzE(a_rho_px[i], a_rho_py[i], a_rho_pz[i], a_rho_energy[i]);

                  t_lv = (vpho_lv - rho_lv);
                  t = t_lv.Mag2();

                  //double tl = rho_lv.Mag2() - Q2 - 2*rho_lv.E()*vpho_lv.E() + 2*(rho_lv.Vect()*vpho_lv.Vect()*cos(vpho_lv.Angle(rho_lv.Vect())) );
                  lc = 0.1973*(2 * vpho_lv.E()) / (Q2 + pow(m_rho,2));           //rho_lv.Mag2());  (0.1973)

                  Zh = rho_lv.E()/vpho_lv.E();

                  rho_mass = a_rho_mass[i];
                  rho_px = a_rho_px[i];
                  rho_py = a_rho_py[i];
                  rho_pz = a_rho_pz[i];
                  rho_energy = a_rho_energy[i];

                  t_hist->Fill(-t);
                  z_hist->Fill(Zh);
                  q_hist->Fill(Q2);
                  qz_hist->Fill(Q2,lc);
                  lc_hist->Fill(lc);
                 //   w_hist->Fill(w);
                  
                  rm_no->Fill(rho_mass);
                  //std::cout << rho_lv.Mag() << std::endl;
                  //std::cout << w << " | " << Q2 << std::endl;
                  
                  r_count++;
                  
                  if (w > 2 && Q2 >= 1){
                  ecount_q++;
                  }

                  if(w > 2 && lc <= 1.0 ){
               rw_count++;
              rm_w->Fill(rho_mass); 
                 ewp_hist->Fill(sqrt(pow(a_px[0],2) + pow(a_py[0],2) + pow(a_pz[0],2)));
                 pipwp_hist->Fill(sqrt(pow(a_pip_px[i],2) + pow(a_pip_py[i],2) + pow(a_pip_pz[i],2)));
                 pimwp_hist->Fill(sqrt(pow(a_pim_px[i],2) + pow(a_pim_py[i],2) + pow(a_pim_pz[i],2)));
                 
                   
                     if(-t > 0.1 && -t < 0.5){
                       rm_wt->Fill(rho_mass); 
                       rwt_count++;
                       
                       ewtp_hist->Fill(sqrt(pow(a_px[0],2) + pow(a_py[0],2) + pow(a_pz[0],2)));
                 pipwtp_hist->Fill(sqrt(pow(a_pip_px[i],2) + pow(a_pip_py[i],2) + pow(a_pip_pz[i],2)));
                 pimwtp_hist->Fill(sqrt(pow(a_pim_px[i],2) + pow(a_pim_py[i],2) + pow(a_pim_pz[i],2)));
                       
                        if(Zh > 0.8 && Zh < 1.0 && Q2 >= 1){
                          rm_wtz->Fill(rho_mass);
                         // q_hist->Fill(Q2);
                         ecount_ful++;
                      
                      ewtzp_hist->Fill(sqrt(pow(a_px[0],2) + pow(a_py[0],2) + pow(a_pz[0],2)));
                 pipwtzp_hist->Fill(sqrt(pow(a_pip_px[i],2) + pow(a_pip_py[i],2) + pow(a_pip_pz[i],2)));
                 pimwtzp_hist->Fill(sqrt(pow(a_pim_px[i],2) + pow(a_pim_py[i],2) + pow(a_pim_pz[i],2)));   
                 rr1_count++;
                           if(Q2 >= 1 && Q2 <2){
         			r1_count++;
         			//rr1_count++;
                              		p_hist->Fill(rho_mass);
				
                            }else if(Q2 >= 2 && Q2 <2.5){
                              ip_hist->Fill(rho_mass);
                              r2_count++;
                              rr2_count++;
                            }else if(Q2 >= 2.5 && Q2 <3){
    			r3_count++;
    			rr3_count++;
                              		iip_hist->Fill(rho_mass);
				
                           }else if(Q2 >= 3 && Q2 <3.5){
                              iiip_hist->Fill(rho_mass);
                              r4_count++;
                              rr4_count++;
                           }else if(Q2 >= 3.5 && Q2 < 4.5){
                              ivp_hist->Fill(rho_mass);
                              r5_count++;
                              rr5_count++;
                           }else if(Q2 >= 4.5 && Q2 <6){
                              vp_hist->Fill(rho_mass);
                              r6_count++;
                              rr6_count++;
                        
                           } 
                          
                        
                     }
                  }
                 //w_cut = 1;
                 //z_cut = 1;
                 //t_cut = 1;
                     
               }
            }
            //if(w_cut > 0 && z_cut > 0 && t_cut > 0){
              // if(rho_counter > 0){
                //  e_neg_pi_num++;
              // }
            //}
         
         }
      
         counter++;
      
      
   }
    
   }
   
   /*
   if (p < 147){
   
   
   
   if (v4_a[p] < 18337){
       
       
       tvz1_hist->Scale(1.0 / tvz1_hist->Integral());
      //TF1 *tvz1_fit = new TF1("tvz1_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
      //    tvz1_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
      //    tvz1_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz12_fit = new TF1("tvz12_fit", "gaus(0)", -10.5, -2);
          tvz12_fit->SetParameters(0.01, -5, 1);
          tvz12_fit->SetParNames("Copper const", "Copper mean", "Copper sigma");
   //tvz1_hist->Fit(tvz1_fit, "R");
   tvz1_hist->Fit(tvz12_fit, "R");
   
   tvz2_hist->Scale(1.0 / tvz2_hist->Integral());
   //   TF1 *tvz2_fit = new TF1("tvz2_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
   //       tvz2_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
   //       tvz2_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz22_fit = new TF1("tvz22_fit", "gaus(0)", -10.5, -2);
          tvz22_fit->SetParameters(0.01, -5, 1);
          tvz22_fit->SetParNames("Copper const", "Copper mean", "Copper sigma");
   //tvz2_hist->Fit(tvz2_fit, "R");
   tvz2_hist->Fit(tvz22_fit, "R");
     
     
     tvz3_hist->Scale(1.0 / tvz3_hist->Integral());
     // TF1 *tvz3_fit = new TF1("tvz3_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
     //     tvz3_fit->SetParameters(0.01, -18 0.5, 0.01, -14, 0.5);
     //     tvz3_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz32_fit = new TF1("tvz32_fit", "gaus(0)", -10.5, -2);
          tvz32_fit->SetParameters(0.01, -5, 1);
          tvz32_fit->SetParNames("Copper const", "Copper mean", "Copper sigma");
   //tvz3_hist->Fit(tvz3_fit, "R");
   tvz3_hist->Fit(tvz32_fit, "R");
     
     
     tvz4_hist->Scale(1.0 / tvz4_hist->Integral());
     // TF1 *tvz4_fit = new TF1("tvz4_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
     //     tvz4_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
     //     tvz4_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz42_fit = new TF1("tvz42_fit", "gaus(0)", -10.5, -2);
          tvz42_fit->SetParameters(0.01, -5, 1);
          tvz42_fit->SetParNames("Copper const", "Copper mean", "Copper sigma");
   //tvz4_hist->Fit(tvz4_fit, "R");
   tvz4_hist->Fit(tvz42_fit, "R");
     
     
     tvz5_hist->Scale(1.0 / tvz5_hist->Integral());
     // TF1 *tvz5_fit = new TF1("tvz5_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
     //     tvz5_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
     //     tvz5_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz52_fit = new TF1("tvz52_fit", "gaus(0)", -10.5, -2);
          tvz52_fit->SetParameters(0.01, -5, 1);
          tvz52_fit->SetParNames("Copper const", "Copper mean", "Copper sigma");
   //tvz5_hist->Fit(tvz5_fit, "R");
   tvz5_hist->Fit(tvz52_fit, "R");
     
     
     tvz6_hist->Scale(1.0 / tvz6_hist->Integral());
     // TF1 *tvz6_fit = new TF1("tvz6_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
     //     tvz6_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
     //     tvz6_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz62_fit = new TF1("tvz62_fit", "gaus(0)", -10.5, -2);
          tvz62_fit->SetParameters(0.01, -5, 1);
          tvz62_fit->SetParNames("Copper const", "Copper mean", "Copper sigma");
   //tvz6_hist->Fit(tvz6_fit, "R");
   tvz6_hist->Fit(tvz62_fit, "R");
   
   tvz1_g->AddPoint(v4_a[p], tvz12_fit->GetParameter(1));
   //tvz1_g->AddPoint(v4_a[p], tvz12_fit->GetParameter(4));
   
   tvz2_g->AddPoint(v4_a[p], tvz22_fit->GetParameter(1));
   //tvz2_g->AddPoint(v4_a[p], tvz22_fit->GetParameter(4));
   
   tvz3_g->AddPoint(v4_a[p], tvz32_fit->GetParameter(1));
   //tvz3_g->AddPoint(v4_a[p], tvz32_fit->GetParameter(4));
   
   tvz4_g->AddPoint(v4_a[p], tvz42_fit->GetParameter(1));
   //tvz4_g->AddPoint(v4_a[p], tvz42_fit->GetParameter(4));
   
   tvz5_g->AddPoint(v4_a[p], tvz52_fit->GetParameter(1));
   //tvz5_g->AddPoint(v4_a[p], tvz52_fit->GetParameter(4));
   
   tvz6_g->AddPoint(v4_a[p], tvz62_fit->GetParameter(1));
   //tvz6_g->AddPoint(v4_a[p], tvz62_fit->GetParameter(4));
   /////////
   tvz1_hist->GetXaxis()->SetTitle(std::to_string(v4_a[p]).c_str());
       tvz1_hist->Draw();
       tvz1_hist->Write();
    
  tvz2_hist->GetXaxis()->SetTitle(std::to_string(v4_a[p]).c_str());     
       tvz2_hist->Draw();
       tvz2_hist->Write();
       
  tvz3_hist->GetXaxis()->SetTitle(std::to_string(v4_a[p]).c_str());     
       tvz3_hist->Draw();
       tvz3_hist->Write();
       
  tvz4_hist->GetXaxis()->SetTitle(std::to_string(v4_a[p]).c_str());     
       tvz4_hist->Draw();
       tvz4_hist->Write();
       
  tvz5_hist->GetXaxis()->SetTitle(std::to_string(v4_a[p]).c_str());    
       tvz5_hist->Draw();
       tvz5_hist->Write();
       
  tvz6_hist->GetXaxis()->SetTitle(std::to_string(v4_a[p]).c_str());    
       tvz6_hist->Draw();
       tvz6_hist->Write();
   ////////
   } else {
    tvz1_hist->Scale(1.0 / tvz1_hist->Integral());
      //TF1 *tvz1_fit = new TF1("tvz1_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
      //    tvz1_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
      //    tvz1_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz12_fit = new TF1("tvz12_fit", "gaus(0)", -10.5, -2);
          tvz12_fit->SetParameters(0.01, -7, 1);
          tvz12_fit->SetParNames("Copper const", "Copper mean", "Copper sigma");
   //tvz1_hist->Fit(tvz1_fit, "R");
   tvz1_hist->Fit(tvz12_fit, "R");
   
   tvz2_hist->Scale(1.0 / tvz2_hist->Integral());
   //   TF1 *tvz2_fit = new TF1("tvz2_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
   //       tvz2_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
   //       tvz2_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz22_fit = new TF1("tvz22_fit", "gaus(0)", -10.5, -2);
          tvz22_fit->SetParameters(0.01, -7, 1);
          tvz22_fit->SetParNames("Copper const", "Copper mean", "Copper sigma");
   //tvz2_hist->Fit(tvz2_fit, "R");
   tvz2_hist->Fit(tvz22_fit, "R");
     
     
     tvz3_hist->Scale(1.0 / tvz3_hist->Integral());
     // TF1 *tvz3_fit = new TF1("tvz3_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
     //     tvz3_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
     //     tvz3_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz32_fit = new TF1("tvz32_fit", "gaus(0)", -10.5, -2);
          tvz32_fit->SetParameters(0.01, -7, 1);
          tvz32_fit->SetParNames("Copper const", "Copper mean", "Copper sigma");
   //tvz3_hist->Fit(tvz3_fit, "R");
   tvz3_hist->Fit(tvz32_fit, "R");
     
     
     tvz4_hist->Scale(1.0 / tvz4_hist->Integral());
     // TF1 *tvz4_fit = new TF1("tvz4_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
     //     tvz4_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
     //     tvz4_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz42_fit = new TF1("tvz42_fit", "gaus(0)", -10.5, -2);
          tvz42_fit->SetParameters(0.01, -7, 1);
          tvz42_fit->SetParNames("Copper const", "Copper mean", "Copper sigma");
   //tvz4_hist->Fit(tvz4_fit, "R");
   tvz4_hist->Fit(tvz42_fit, "R");
     
     
     tvz5_hist->Scale(1.0 / tvz5_hist->Integral());
     // TF1 *tvz5_fit = new TF1("tvz5_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
     //     tvz5_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
     //     tvz5_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz52_fit = new TF1("tvz52_fit", "gaus(0)", -10.5, -2);
          tvz52_fit->SetParameters(0.01, -7, 1);
          tvz52_fit->SetParNames("Copper const", "Copper mean", "Copper sigma");
   //tvz5_hist->Fit(tvz5_fit, "R");
   tvz5_hist->Fit(tvz52_fit, "R");
     
     
     tvz6_hist->Scale(1.0 / tvz6_hist->Integral());
     // TF1 *tvz6_fit = new TF1("tvz6_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
     //     tvz6_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
     //     tvz6_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz62_fit = new TF1("tvz62_fit", "gaus(0)", -10.5, -2);
          tvz62_fit->SetParameters(0.01, -7, 1);
          tvz62_fit->SetParNames("Copper const", "Copper mean", "Copper sigma");
   //tvz6_hist->Fit(tvz6_fit, "R");
   tvz6_hist->Fit(tvz62_fit, "R");
   
   tvz1_g->AddPoint(v4_a[p], tvz12_fit->GetParameter(1));
   //tvz1_g->AddPoint(v4_a[p], tvz12_fit->GetParameter(4));
   
   tvz2_g->AddPoint(v4_a[p], tvz22_fit->GetParameter(1));
   //tvz2_g->AddPoint(v4_a[p], tvz22_fit->GetParameter(4));
   
   tvz3_g->AddPoint(v4_a[p], tvz32_fit->GetParameter(1));
   //tvz3_g->AddPoint(v4_a[p], tvz32_fit->GetParameter(4));
   
   tvz4_g->AddPoint(v4_a[p], tvz42_fit->GetParameter(1));
   //tvz4_g->AddPoint(v4_a[p], tvz42_fit->GetParameter(4));
   
   tvz5_g->AddPoint(v4_a[p], tvz52_fit->GetParameter(1));
   //tvz5_g->AddPoint(v4_a[p], tvz52_fit->GetParameter(4));
   
   tvz6_g->AddPoint(v4_a[p], tvz62_fit->GetParameter(1));
   //tvz6_g->AddPoint(v4_a[p], tvz62_fit->GetParameter(4));
   
   
   }
   
   
   }else if (v4_a[p] == 18448  ||  v4_a[p] == 19090 ){
   tvz1_hist->Scale(1.0 / tvz1_hist->Integral());
     // TF1 *tvz1_fit = new TF1("tvz1_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
     //     tvz1_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
     //     tvz1_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz12_fit = new TF1("tvz12_fit", "gaus(0)+gaus(3)+pol1(6)", -12.5, -2);
          tvz12_fit->SetParameters(0.01, -9, 1, 0.01, -3, 1, 0.01, 0);
          tvz12_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
  // tvz1_hist->Fit(tvz1_fit, "R");
   tvz1_hist->Fit(tvz12_fit, "R");
   
   tvz2_hist->Scale(1.0 / tvz2_hist->Integral());
   //   TF1 *tvz2_fit = new TF1("tvz2_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
   //       tvz2_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
   //       tvz2_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz22_fit = new TF1("tvz22_fit", "gaus(0)+gaus(3)+pol1(6)", -12.5, -2);
          tvz22_fit->SetParameters(0.01, -9, 1, 0.01, -3, 1, 0.01, 0);
          tvz22_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
   //tvz2_hist->Fit(tvz2_fit, "R");
   tvz2_hist->Fit(tvz22_fit, "R");
     
     
     tvz3_hist->Scale(1.0 / tvz3_hist->Integral());
   //   TF1 *tvz3_fit = new TF1("tvz3_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
   //       tvz3_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
   //       tvz3_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz32_fit = new TF1("tvz32_fit", "gaus(0)+gaus(3)+pol1(6)", -12.5, -2);
          tvz32_fit->SetParameters(0.01, -9, 1, 0.01, -3, 1, 0.01, 0);
          tvz32_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
  // tvz3_hist->Fit(tvz3_fit, "R");
   tvz3_hist->Fit(tvz32_fit, "R");
     
     
     tvz4_hist->Scale(1.0 / tvz4_hist->Integral());
  //    TF1 *tvz4_fit = new TF1("tvz4_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
  //        tvz4_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
  //        tvz4_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz42_fit = new TF1("tvz42_fit", "gaus(0)+gaus(3)+pol1(6)", -12.5, -2);
          tvz42_fit->SetParameters(0.01, -9, 1, 0.01, -3, 1, 0.01, 0);
          tvz42_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
  // tvz4_hist->Fit(tvz4_fit, "R");
   tvz4_hist->Fit(tvz42_fit, "R");
     
     
     tvz5_hist->Scale(1.0 / tvz5_hist->Integral());
  //    TF1 *tvz5_fit = new TF1("tvz5_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
  //        tvz5_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
  //        tvz5_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz52_fit = new TF1("tvz52_fit", "gaus(0)+gaus(3)+pol1(6)", -12.5, -2);
          tvz52_fit->SetParameters(0.01, -9, 1, 0.01, -3, 1, 0.01, 0);
          tvz52_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
 //  tvz5_hist->Fit(tvz5_fit, "R");
   tvz5_hist->Fit(tvz52_fit, "R");
     
     
     tvz6_hist->Scale(1.0 / tvz6_hist->Integral());
 //     TF1 *tvz6_fit = new TF1("tvz6_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
 //         tvz6_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
 //         tvz6_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz62_fit = new TF1("tvz62_fit", "gaus(0)+gaus(3)+pol1(6)", -12.5, -2);
          tvz62_fit->SetParameters(0.01, -9, 1, 0.01, -3, 1, 0.01, 0);
          tvz62_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
  // tvz6_hist->Fit(tvz6_fit, "R");
   tvz6_hist->Fit(tvz62_fit, "R");
   
   /////////////
   //if(tvz22_fit->GetParameter(4) < 20000 ){
   tvz1_g->AddPoint(v4_a[p], tvz12_fit->GetParameter(1));
   tvz1_g->AddPoint(v4_a[p], tvz12_fit->GetParameter(4));
   
   tvz2_g->AddPoint(v4_a[p], tvz22_fit->GetParameter(1));
   tvz2_g->AddPoint(v4_a[p], tvz22_fit->GetParameter(4));
   
   tvz3_g->AddPoint(v4_a[p], tvz32_fit->GetParameter(1));
   tvz3_g->AddPoint(v4_a[p], tvz32_fit->GetParameter(4));
   
   tvz4_g->AddPoint(v4_a[p], tvz42_fit->GetParameter(1));
   tvz4_g->AddPoint(v4_a[p], tvz42_fit->GetParameter(4));
   
   tvz5_g->AddPoint(v4_a[p], tvz52_fit->GetParameter(1));
   tvz5_g->AddPoint(v4_a[p], tvz52_fit->GetParameter(4));
   
   tvz6_g->AddPoint(v4_a[p], tvz62_fit->GetParameter(1));
   tvz6_g->AddPoint(v4_a[p], tvz62_fit->GetParameter(4));
   /////////////
   //} else {
   tvz1_hist->GetXaxis()->SetTitle(std::to_string(v4_a[p]).c_str());
       tvz1_hist->Draw();
       tvz1_hist->Write();
    
  tvz2_hist->GetXaxis()->SetTitle(std::to_string(v4_a[p]).c_str());     
       tvz2_hist->Draw();
       tvz2_hist->Write();
       
  tvz3_hist->GetXaxis()->SetTitle(std::to_string(v4_a[p]).c_str());     
       tvz3_hist->Draw();
       tvz3_hist->Write();
       
  tvz4_hist->GetXaxis()->SetTitle(std::to_string(v4_a[p]).c_str());     
       tvz4_hist->Draw();
       tvz4_hist->Write();
       
  tvz5_hist->GetXaxis()->SetTitle(std::to_string(v4_a[p]).c_str());    
  
       tvz5_hist->Draw();
       tvz5_hist->Write();
       
  tvz6_hist->GetXaxis()->SetTitle(std::to_string(v4_a[p]).c_str());    
       tvz6_hist->Draw();
       tvz6_hist->Write();
   
   //}
   
     } else{
   tvz1_hist->Scale(1.0 / tvz1_hist->Integral());
      TF1 *tvz1_fit = new TF1("tvz1_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
          tvz1_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
          tvz1_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz12_fit = new TF1("tvz12_fit", "gaus(0)+gaus(3)+pol1(6)", -10.5, -2);
          tvz12_fit->SetParameters(0.01, -8, 1, 0.01, -2, 1, 0.01, 0);
          tvz12_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
   tvz1_hist->Fit(tvz1_fit, "R");
   tvz1_hist->Fit(tvz12_fit, "R+");
   
   tvz2_hist->Scale(1.0 / tvz2_hist->Integral());
      TF1 *tvz2_fit = new TF1("tvz2_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
          tvz2_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
          tvz2_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz22_fit = new TF1("tvz22_fit", "gaus(0)+gaus(3)+pol1(6)", -10.5, -2);
          tvz22_fit->SetParameters(0.01, -8, 1, 0.01, -2, 1, 0.01, 0);
          tvz22_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
   tvz2_hist->Fit(tvz2_fit, "R");
   tvz2_hist->Fit(tvz22_fit, "R+");
     
     
     tvz3_hist->Scale(1.0 / tvz3_hist->Integral());
      TF1 *tvz3_fit = new TF1("tvz3_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
          tvz3_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
          tvz3_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz32_fit = new TF1("tvz32_fit", "gaus(0)+gaus(3)+pol1(6)", -10.5, -2);
          tvz32_fit->SetParameters(0.01, -8, 1, 0.01, -2, 1, 0.01, 0);
          tvz32_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
   tvz3_hist->Fit(tvz3_fit, "R");
   tvz3_hist->Fit(tvz32_fit, "R+");
     
     
     tvz4_hist->Scale(1.0 / tvz4_hist->Integral());
      TF1 *tvz4_fit = new TF1("tvz4_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
          tvz4_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
          tvz4_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz42_fit = new TF1("tvz42_fit", "gaus(0)+gaus(3)+pol1(6)", -10.5, -2);
          tvz42_fit->SetParameters(0.01, -8, 1, 0.01, -2, 1, 0.01, 0);
          tvz42_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
   tvz4_hist->Fit(tvz4_fit, "R");
   tvz4_hist->Fit(tvz42_fit, "R+");
     
     
     tvz5_hist->Scale(1.0 / tvz5_hist->Integral());
      TF1 *tvz5_fit = new TF1("tvz5_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
          tvz5_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
          tvz5_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz52_fit = new TF1("tvz52_fit", "gaus(0)+gaus(3)+pol1(6)", -10.5, -2);
          tvz52_fit->SetParameters(0.01, -8, 1, 0.01, -2, 1, 0.01, 0);
          tvz52_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
   tvz5_hist->Fit(tvz5_fit, "R");
   tvz5_hist->Fit(tvz52_fit, "R+");
     
     
     tvz6_hist->Scale(1.0 / tvz6_hist->Integral());
      TF1 *tvz6_fit = new TF1("tvz6_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
          tvz6_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
          tvz6_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *tvz62_fit = new TF1("tvz62_fit", "gaus(0)+gaus(3)+pol1(6)", -10.5, -2);
          tvz62_fit->SetParameters(0.01, -8, 1, 0.01, -2, 1, 0.01, 0);
          tvz62_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
   tvz6_hist->Fit(tvz6_fit, "R");
   tvz6_hist->Fit(tvz62_fit, "R+");
   
   if(tvz22_fit->GetParameter(4) < -2.2 ){
   tvz1_g->AddPoint(v4_a[p], tvz12_fit->GetParameter(1));
   tvz1_g->AddPoint(v4_a[p], tvz12_fit->GetParameter(4));
   
   tvz2_g->AddPoint(v4_a[p], tvz22_fit->GetParameter(1));
   tvz2_g->AddPoint(v4_a[p], tvz22_fit->GetParameter(4));
   
   tvz3_g->AddPoint(v4_a[p], tvz32_fit->GetParameter(1));
   tvz3_g->AddPoint(v4_a[p], tvz32_fit->GetParameter(4));
   
   tvz4_g->AddPoint(v4_a[p], tvz42_fit->GetParameter(1));
   tvz4_g->AddPoint(v4_a[p], tvz42_fit->GetParameter(4));
   
   tvz5_g->AddPoint(v4_a[p], tvz52_fit->GetParameter(1));
   tvz5_g->AddPoint(v4_a[p], tvz52_fit->GetParameter(4));
   
   tvz6_g->AddPoint(v4_a[p], tvz62_fit->GetParameter(1));
   tvz6_g->AddPoint(v4_a[p], tvz62_fit->GetParameter(4));
   } 
   
     }
     */
      if( v6_a[p] < 40){ 
          double cn = (v5_a[p] * 1000000.0);
          double ne_count = fe_count  / cn;
          double npip_count = fpip_count  / cn;
          double npim_count = fpim_count  / cn;
          
          t_charge = t_charge + (v5_a[p]* 1000000.0);
          nt_charge = nt_charge + cn;
          
        //  std::cout << "electron count = " << v5_a[p] << std::endl;
        //  std::cout << "pip count = " << v6_a[p] << std::endl;
        //  std::cout << "pim count = " << npim_count << std::endl;
          
          if( ne_count < 3000000000000000000000){
          te_g->AddPoint(v4_a[p], ne_count);
          } else {
          bad_e_run[ber] = v4_a[p];
          ber++;
          }
          if ( npip_count < 6000000000000000000000000 ){
          tpip_g->AddPoint(v4_a[p], npip_count);
          } else {
          bad_pip_run[bpipr] = v4_a[p];
          bpipr++;
          }
          if (npim_count < 60000000000000000000){
          tpim_g->AddPoint(v4_a[p], npim_count);
          } else {
          bad_pim_run[bpimr] = v4_a[p];
          bpimr++;
          }
          double arr1_count = rr1_count / cn;
          
          if( arr1_count < 2000000000000000){
          tr1_g->AddPoint(v4_a[p], arr1_count);
          } else {
          bad_r_run[brr] = v4_a[p];
          brr++;
          }
          tr2_g->AddPoint(v4_a[p], rr2_count);
          tr3_g->AddPoint(v4_a[p], rr3_count);
          tr4_g->AddPoint(v4_a[p], rr4_count);
          tr5_g->AddPoint(v4_a[p], rr5_count);
          tr6_g->AddPoint(v4_a[p], rr6_count);
      } else {
          double cnn = ((v5_a[p] * 1000000.0) * (40.0 / v6_a[p]));
          double ne_count = fe_count / cnn;
          double npip_count = fpip_count / cnn;
          double npim_count = fpim_count / cnn;
          
          t_charge = t_charge + (v5_a[p]* 1000000.0);
          nt_charge = nt_charge + cnn;
          
         // std::cout << v6_a[p] << std::endl;
         // std::cout << cnn << std::endl;
         // std::cout << v5_a[p] << std::endl;
         // std::cout << "electron count = " << fe_count << std::endl;
         // std::cout << "pip count = " << fpip_count << std::endl;
         // std::cout << "pim count = " << fpim_count << std::endl;
          
          double arr1_count = rr1_count / cnn;
          
          if( ne_count < 30000000000000){
          te_g->AddPoint(v4_a[p], ne_count);
          } else {
          bad_e_run[ber] = v4_a[p];
          ber++;
          }
          if( npip_count < 60000000000000){
          tpip_g->AddPoint(v4_a[p], npip_count);
          } else {
          bad_pip_run[bpipr] = v4_a[p];
          bpipr++;
          }
          if( npim_count < 600000000000000000){
          tpim_g->AddPoint(v4_a[p], npim_count);
          } else {
          bad_pim_run[bpimr] = v4_a[p];
          bpimr++;
          }
          if( arr1_count < 2000000000000000){
          tr1_g->AddPoint(v4_a[p], arr1_count);
          } else {
          bad_r_run[brr] = v4_a[p];
          brr++;
          }
          tr2_g->AddPoint(v4_a[p], rr2_count);
          tr3_g->AddPoint(v4_a[p], rr3_count);
          tr4_g->AddPoint(v4_a[p], rr4_count);
          tr5_g->AddPoint(v4_a[p], rr5_count);
          tr6_g->AddPoint(v4_a[p], rr6_count);
          
          
      }
     
     
     
     
     
      std::cout << "p = " << ((p + 1.0 - start_p) / (end_p - start_p)) * 100 << "%" <<std::endl;
      
      
     
      
   }
//}
   std::cout << "e with pid = " << count_1 << std::endl;
   std::cout << "e with status, vz, chi2pid = " << count_2 << std::endl;
   std::cout << "e wit all = " << count_3 << std::endl;
   std::cout << "e pip pim events = " << e_pip_pim_num << std::endl;
   
   
   
   //std::cout << "e count " << e_count << std::endl;
   
   
   
   
   e_hist->GetXaxis()->SetTitle("e");
   e_hist->Draw();
   e_hist->Write();

   pip_hist->GetXaxis()->SetTitle("pip");
   pip_hist->Draw();
   pip_hist->Write();

   pim_hist->GetXaxis()->SetTitle("pim");
   pim_hist->Draw();
   pim_hist->Write();
   
   e_pht_hist->GetXaxis()->SetTitle("#phi_{e} [deg]");
   e_pht_hist->GetYaxis()->SetTitle("#theta_{e} [deg]");
   e_pht_hist->Draw();
   e_pht_hist->Write();
   
   pip_pht_hist->GetXaxis()->SetTitle("#phi_{#pi^{+}} [deg]");
   pip_pht_hist->GetYaxis()->SetTitle("#theta_{#pi^{+}} [deg]");
   pip_pht_hist->Draw();
   pip_pht_hist->Write();
   
   pim_pht_hist->GetXaxis()->SetTitle("#phi_{#pi^{-}} [deg]");
   pim_pht_hist->GetYaxis()->SetTitle("#theta_{#pi^{-}} [deg]");
   pim_pht_hist->Draw();
   pim_pht_hist->Write();
   
   phih_hist->GetXaxis()->SetTitle("#phi_{h} [deg]");
   //phih_hist->GetYaxis()->SetTitle("#theta_{#pi^{-}} [deg]");
   phih_hist->Draw();
   phih_hist->Write();
   
   xbj_hist->GetXaxis()->SetTitle("X");
   //phih_hist->GetYaxis()->SetTitle("#theta_{#pi^{-}} [deg]");
   xbj_hist->Draw();
   xbj_hist->Write();
   
   ybj_hist->GetXaxis()->SetTitle("Y");
   //phih_hist->GetYaxis()->SetTitle("#theta_{#pi^{-}} [deg]");
   ybj_hist->Draw();
   ybj_hist->Write();
   
   
   tz_hist->GetXaxis()->SetTitle("#phi_{e} [deg]");
   tz_hist->GetYaxis()->SetTitle("vz [cm]");
   tz_hist->Draw();
   tz_hist->Write();
   
   echi_hist->GetXaxis()->SetTitle("e #Chi^{2}");
   echi_hist->Draw();
   echi_hist->Write();
   
   pipchi_hist->GetXaxis()->SetTitle("#pie^{+} #Chi^{2}");
   pipchi_hist->Draw();
   pipchi_hist->Write();
   
   pimchi_hist->GetXaxis()->SetTitle("#pie^{-} #Chi^{2}");
   pimchi_hist->Draw();
   pimchi_hist->Write();

   q_hist->GetXaxis()->SetTitle("Q2");
   q_hist->Draw();
   q_hist->Write();
   
   qz_hist->GetXaxis()->SetTitle("Q2 [GeV^{2}]");
   qz_hist->GetYaxis()->SetTitle("lc [fm]");
   qz_hist->Draw();
   qz_hist->Write();
   
   sam_hist->GetXaxis()->SetTitle("E_{dep,tot}/p [GeV]");
   sam_hist->Draw();
   sam_hist->Write();
   
   samp_hist->GetXaxis()->SetTitle("p [GeV/c]");
   samp_hist->GetYaxis()->SetTitle("E_{dep,tot}/p [GeV]");
   samp_hist->Draw();
   samp_hist->Write();
   
   samv_hist->GetXaxis()->SetTitle("lv [cm]");
   samv_hist->GetYaxis()->SetTitle("E_{dep,tot}/p [GeV]");
   samv_hist->Draw();
   samv_hist->Write();
   
   samw_hist->GetXaxis()->SetTitle("lw [cm]");
   samw_hist->GetYaxis()->SetTitle("E_{dep,tot}/p [GeV]");
   samw_hist->Draw();
   samw_hist->Write();
   
   
   auto samp_g = new TGraph();
      samp_g->SetName("samp_g");
      samp_g->SetTitle("samp mean");
      samp_g->SetMarkerColor(kRed);
      samp_g->SetMarkerStyle(kFullCircle);
      samp_g->SetLineWidth(0);
      
  auto sampu_g = new TGraph();
      sampu_g->SetName("sampu_g");
      sampu_g->SetTitle("samp up mean");
      sampu_g->SetMarkerColor(kYellow );
      sampu_g->SetMarkerStyle(kFullCircle);
      sampu_g->SetLineWidth(0);
    
   auto sampd_g = new TGraph();
      sampd_g->SetName("sampd_g");
      sampd_g->SetTitle("samp down mean");
      sampd_g->SetMarkerColor(kYellow);
      sampd_g->SetMarkerStyle(kFullCircle);
      sampd_g->SetLineWidth(0);
      
      TH1F *samp01_hist = new TH1F();
      
      
      auto samp1_g = new TGraph();
      samp1_g->SetName("samp1_g");
      samp1_g->SetTitle("samp1 mean");
      samp1_g->SetMarkerColor(kRed);
      samp1_g->SetMarkerStyle(kFullCircle);
      samp1_g->SetLineWidth(0);
      
  auto sampu1_g = new TGraph();
      sampu1_g->SetName("sampu1_g");
      sampu1_g->SetTitle("samp1 up mean");
      sampu1_g->SetMarkerColor(kYellow );
      sampu1_g->SetMarkerStyle(kFullCircle);
      sampu1_g->SetLineWidth(0);
    
   auto sampd1_g = new TGraph();
      sampd1_g->SetName("sampd1_g");
      sampd1_g->SetTitle("samp1 down mean");
      sampd1_g->SetMarkerColor(kYellow);
      sampd1_g->SetMarkerStyle(kFullCircle);
      sampd1_g->SetLineWidth(0);
      
      TH1F *samp11_hist = new TH1F();
      
      
      auto samp2_g = new TGraph();
      samp2_g->SetName("samp2_g");
      samp2_g->SetTitle("samp2 mean");
      samp2_g->SetMarkerColor(kRed);
      samp2_g->SetMarkerStyle(kFullCircle);
      samp2_g->SetLineWidth(0);
      
  auto sampu2_g = new TGraph();
      sampu2_g->SetName("sampu2_g");
      sampu2_g->SetTitle("samp2 up mean");
      sampu2_g->SetMarkerColor(kYellow );
      sampu2_g->SetMarkerStyle(kFullCircle);
      sampu2_g->SetLineWidth(0);
    
   auto sampd2_g = new TGraph();
      sampd2_g->SetName("sampd2_g");
      sampd2_g->SetTitle("samp2 down mean");
      sampd2_g->SetMarkerColor(kYellow);
      sampd2_g->SetMarkerStyle(kFullCircle);
      sampd2_g->SetLineWidth(0);
      
      TH1F *samp12_hist = new TH1F();
      
      
      auto samp3_g = new TGraph();
      samp3_g->SetName("samp3_g");
      samp3_g->SetTitle("samp3 mean");
      samp3_g->SetMarkerColor(kRed);
      samp3_g->SetMarkerStyle(kFullCircle);
      samp3_g->SetLineWidth(0);
      
  auto sampu3_g = new TGraph();
      sampu3_g->SetName("sampu3_g");
      sampu3_g->SetTitle("samp3 up mean");
      sampu3_g->SetMarkerColor(kYellow );
      sampu3_g->SetMarkerStyle(kFullCircle);
      sampu3_g->SetLineWidth(0);
    
   auto sampd3_g = new TGraph();
      sampd3_g->SetName("sampd3_g");
      sampd3_g->SetTitle("samp3 down mean");
      sampd3_g->SetMarkerColor(kYellow);
      sampd3_g->SetMarkerStyle(kFullCircle);
      sampd3_g->SetLineWidth(0);
      
      TH1F *samp13_hist = new TH1F();
      
      
      auto samp4_g = new TGraph();
      samp4_g->SetName("samp4_g");
      samp4_g->SetTitle("samp4 mean");
      samp4_g->SetMarkerColor(kRed);
      samp4_g->SetMarkerStyle(kFullCircle);
      samp4_g->SetLineWidth(0);
      
  auto sampu4_g = new TGraph();
      sampu4_g->SetName("sampu4_g");
      sampu4_g->SetTitle("samp4 up mean");
      sampu4_g->SetMarkerColor(kYellow );
      sampu4_g->SetMarkerStyle(kFullCircle);
      sampu4_g->SetLineWidth(0);
    
   auto sampd4_g = new TGraph();
      sampd4_g->SetName("sampd4_g");
      sampd4_g->SetTitle("samp4 down mean");
      sampd4_g->SetMarkerColor(kYellow);
      sampd4_g->SetMarkerStyle(kFullCircle);
      sampd4_g->SetLineWidth(0);
      
      TH1F *samp14_hist = new TH1F();
      
      
      auto samp5_g = new TGraph();
      samp5_g->SetName("samp5_g");
      samp5_g->SetTitle("samp5 mean");
      samp5_g->SetMarkerColor(kRed);
      samp5_g->SetMarkerStyle(kFullCircle);
      samp5_g->SetLineWidth(0);
      
  auto sampu5_g = new TGraph();
      sampu5_g->SetName("sampu5_g");
      sampu5_g->SetTitle("samp5 up mean");
      sampu5_g->SetMarkerColor(kYellow );
      sampu5_g->SetMarkerStyle(kFullCircle);
      sampu5_g->SetLineWidth(0);
    
   auto sampd5_g = new TGraph();
      sampd5_g->SetName("sampd5_g");
      sampd5_g->SetTitle("samp5 down mean");
      sampd5_g->SetMarkerColor(kYellow);
      sampd5_g->SetMarkerStyle(kFullCircle);
      sampd5_g->SetLineWidth(0);
      
      TH1F *samp15_hist = new TH1F();
      
      
      auto samp6_g = new TGraph();
      samp6_g->SetName("samp6_g");
      samp6_g->SetTitle("samp6 mean");
      samp6_g->SetMarkerColor(kRed);
      samp6_g->SetMarkerStyle(kFullCircle);
      samp6_g->SetLineWidth(0);
      
  auto sampu6_g = new TGraph();
      sampu6_g->SetName("sampu6_g");
      sampu6_g->SetTitle("samp6 up mean");
      sampu6_g->SetMarkerColor(kYellow );
      sampu6_g->SetMarkerStyle(kFullCircle);
      sampu6_g->SetLineWidth(0);
    
   auto sampd6_g = new TGraph();
      sampd6_g->SetName("sampd6_g");
      sampd6_g->SetTitle("samp6 down mean");
      sampd6_g->SetMarkerColor(kYellow);
      sampd6_g->SetMarkerStyle(kFullCircle);
      sampd6_g->SetLineWidth(0);
      
      TH1F *samp16_hist = new TH1F();
      
  
      
      
   float lwpp = 0.0;
   float upp = 0.0; 
   
   //spline
   
   int pvsize = p_v.size();
   /*
   TH1F *sam1_vhist = new TH1F("sam1", "E_{pcal}/P", 200, 0, 0.4);
   TH1F *sam2_vhist = new TH1F("sam2", "E_{pcal}/P", 200, 0, 0.4);
   TH1F *sam3_vhist = new TH1F("sam3", "E_{pcal}/P", 200, 0, 0.4);
   TH1F *sam4_vhist = new TH1F("sam4", "E_{pcal}/P", 200, 0, 0.4);
   TH1F *sam5_vhist = new TH1F("sam5", "E_{pcal}/P", 200, 0, 0.4);
   TH1F *sam6_vhist = new TH1F("sam6", "E_{pcal}/P", 200, 0, 0.4);
   TH1F *sam7_vhist = new TH1F("sam7", "E_{pcal}/P", 200, 0, 0.4);
   TH1F *sam8_vhist = new TH1F("sam8", "E_{pcal}/P", 200, 0, 0.4);
   TH1F *sam9_vhist = new TH1F("sam9", "E_{pcal}/P", 200, 0, 0.4);
   TH1F *sam10_vhist = new TH1F("sam10", "E_{pcal}/P", 200, 0, 0.4);
   */
   //TF1 *samp_fit = new TF1("samp_fit", "gaus(0)", 0.05, 0.35);
           // samp_fit->SetParameters(0.1, 0.24, 0.03);
            //samp_fit->SetParNames("const", "mean", "sigma");
   /*
   for(int i = 0; i < 11; i++){
   
       upp = lwpp + 1.0;
    
        TH1F *sam_vhist = new TH1F();
        std::cout << "sam_v = " << sam_v[i] << std::endl;
       for(int j = 0; j < pvsize; j++){
           if(p_v[j] >= 0.0 && p_v[j] < 0.55){
                 sam1_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= 0.55 && p_v[j] < 1.1){
                 sam2_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= 1.1 && p_v[j] < 1.65){
                 sam3_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= 1.65 && p_v[j] < 2.2){
                 sam4_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= 2.2 && p_v[j] < 2.75){
                 sam5_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= 2.75 && p_v[j] < 3.3){
                 sam6_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= 3.3 && p_v[j] < 3.85){
                 sam7_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= 3.85 && p_v[j] < 4){
                 sam8_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= lwpp && p_v[j] < upp){
                 sam9_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= lwpp && p_v[j] < upp){
                 sam10_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= lwpp && p_v[j] < upp){
                 sam_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= lwpp && p_v[j] < upp){
                 sam_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= lwpp && p_v[j] < upp){
                 sam_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= lwpp && p_v[j] < upp){
                 sam_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= lwpp && p_v[j] < upp){
                 sam_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= lwpp && p_v[j] < upp){
                 sam_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= lwpp && p_v[j] < upp){
                 sam_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= lwpp && p_v[j] < upp){
                 sam_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= lwpp && p_v[j] < upp){
                 sam_vhist->Fill(sam_v[j]);
           } else if(p_v[j] >= lwpp && p_v[j] < upp){
                 sam_vhist->Fill(sam_v[j]);
           }
       }
       
       
       sam_vhist->Fit(samp_fit, "R");
       sam_vhist->Draw();
       sam_vhist->Write();
       
       
       samp_g->AddPoint( ((lwpp + upp)/2.0) , samp_fit->GetParameter(1));
       sampd_g->AddPoint( ((lwpp + upp)/2.0) , samp_fit->GetParameter(1) - (3 * samp_fit->GetParameter(2)));
       sampu_g->AddPoint( ((lwpp + upp)/2.0) , samp_fit->GetParameter(1) + (3 * samp_fit->GetParameter(2)));
       
       
       lwpp = lwpp + 1.0;
   
   }
   */
   /*
       sam1_vhist->Fit(samp_fit, "R");
       sam1_vhist->Draw();
       sam1_vhist->Write();
       
       sam2_vhist->Fit(samp_fit, "R");
       sam2_vhist->Draw();
       sam2_vhist->Write();
       
       sam3_vhist->Fit(samp_fit, "R");
       sam3_vhist->Draw();
       sam3_vhist->Write();
       
       sam4_vhist->Fit(samp_fit, "R");
       sam4_vhist->Draw();
       sam4_vhist->Write();
       
       sam5_vhist->Fit(samp_fit, "R");
       sam5_vhist->Draw();
       sam5_vhist->Write();
       
       sam6_vhist->Fit(samp_fit, "R");
       sam6_vhist->Draw();
       sam6_vhist->Write();
       
       sam7_vhist->Fit(samp_fit, "R");
       sam7_vhist->Draw();
       sam7_vhist->Write();
       
       sam8_vhist->Fit(samp_fit, "R");
       sam8_vhist->Draw();
       sam8_vhist->Write();
       
       sam9_vhist->Fit(samp_fit, "R");
       sam9_vhist->Draw();
       sam9_vhist->Write();
       
       sam10_vhist->Fit(samp_fit, "R");
       sam10_vhist->Draw();
       sam10_vhist->Write();
   */
   
   
   for(int i = 0; i < 50; i++){
       upp = lwpp + 4.0;
       
       samp01_hist = (TH1F*)samp_hist->ProjectionY("", lwpp, upp);
       samp11_hist = (TH1F*)samp1_hist->ProjectionY("", lwpp, upp);
       samp12_hist = (TH1F*)samp2_hist->ProjectionY("", lwpp, upp);
       samp13_hist = (TH1F*)samp3_hist->ProjectionY("", lwpp, upp);
       samp14_hist = (TH1F*)samp4_hist->ProjectionY("", lwpp, upp);
       samp15_hist = (TH1F*)samp5_hist->ProjectionY("", lwpp, upp);
       samp16_hist = (TH1F*)samp6_hist->ProjectionY("", lwpp, upp);
       TF1 *samp_fit = new TF1("samp_fit", "gaus(0)", 0.05, 0.35);
           // samp_fit->SetParameters(0.1, 0.24, 0.03);
            //samp_fit->SetParNames("const", "mean", "sigma");
       samp1_hist->Fit(samp_fit, "R");
      // samp1_hist->Draw();
      // samp1_hist->Write();
      
      TF1 *samp1_fit = new TF1("samp1_fit", "gaus(0)", 0.05, 0.35);
       samp11_hist->Fit(samp1_fit, "R");
      TF1 *samp2_fit = new TF1("samp2_fit", "gaus(0)", 0.05, 0.35);
       samp12_hist->Fit(samp2_fit, "R");
      TF1 *samp3_fit = new TF1("samp3_fit", "gaus(0)", 0.05, 0.35);
       samp13_hist->Fit(samp3_fit, "R");
      TF1 *samp4_fit = new TF1("samp4_fit", "gaus(0)", 0.05, 0.35);
       samp14_hist->Fit(samp4_fit, "R");
      TF1 *samp5_fit = new TF1("samp5_fit", "gaus(0)", 0.05, 0.35);
       samp15_hist->Fit(samp5_fit, "R");
      TF1 *samp6_fit = new TF1("samp6_fit", "gaus(0)", 0.05, 0.35);
       samp16_hist->Fit(samp6_fit, "R");
       
       samp11_hist->Draw();
       samp11_hist->Write();
       
       if(lwpp >= 28 && lwpp <= 180){
       samp_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0) , samp_fit->GetParameter(1));
       sampd_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp_fit->GetParameter(1) - (3 * samp_fit->GetParameter(2)));
       sampu_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp_fit->GetParameter(1) + (3 * samp_fit->GetParameter(2)));
       
       samp1_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0) , samp1_fit->GetParameter(1));
       sampd1_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp1_fit->GetParameter(1) - (3 * samp1_fit->GetParameter(2)));
       sampu1_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp1_fit->GetParameter(1) + (3 * samp1_fit->GetParameter(2)));
       
       samp2_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0) , samp2_fit->GetParameter(1));
       sampd2_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp2_fit->GetParameter(1) - (3 * samp2_fit->GetParameter(2)));
       sampu2_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp2_fit->GetParameter(1) + (3 * samp2_fit->GetParameter(2)));
       
       samp3_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0) , samp3_fit->GetParameter(1));
       sampd3_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp3_fit->GetParameter(1) - (3 * samp3_fit->GetParameter(2)));
       sampu3_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp3_fit->GetParameter(1) + (3 * samp3_fit->GetParameter(2)));
       
       samp4_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0) , samp4_fit->GetParameter(1));
       sampd4_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp4_fit->GetParameter(1) - (3 * samp4_fit->GetParameter(2)));
       sampu4_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp4_fit->GetParameter(1) + (3 * samp4_fit->GetParameter(2)));
       
       samp5_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0) , samp5_fit->GetParameter(1));
       sampd5_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp5_fit->GetParameter(1) - (3 * samp5_fit->GetParameter(2)));
       sampu5_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp5_fit->GetParameter(1) + (3 * samp5_fit->GetParameter(2)));
       
       samp6_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0) , samp6_fit->GetParameter(1));
       sampd6_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp6_fit->GetParameter(1) - (3 * samp6_fit->GetParameter(2)));
       sampu6_g->AddPoint( (11.0/200.0)  * ((lwpp + upp)/2.0), samp6_fit->GetParameter(1) + (3 * samp6_fit->GetParameter(2)));
       }
       
       lwpp = lwpp + 4.0;
   
   }
   
   
   TF1 *samp_gfit = new TF1("", "pol3", 2, 10);
   samp_g->Fit(samp_gfit);
   samp_g->Draw();
   samp_g->Write();
   TF1 *sampd_gfit = new TF1("", "pol3", 2, 10);
   sampd_g->Fit(sampd_gfit);
   sampd_g->Draw();
   sampd_g->Write();
   TF1 *sampu_gfit = new TF1("", "pol3", 2, 10);
   sampu_g->Fit(sampu_gfit);
   sampu_g->Draw();
   sampu_g->Write();
   
   
   
   TF1 *samp1_gfit = new TF1("", "pol3", 2, 10);
   samp1_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   samp1_g->Fit(samp1_gfit);
   TF1 *sampd1_gfit = new TF1("", "pol3", 2, 10);
   sampd1_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   sampd1_g->Fit(sampd1_gfit);
   TF1 *sampu1_gfit = new TF1("", "pol3", 2, 10);
   sampu1_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   sampu1_g->Fit(sampu1_gfit);
   
   TF1 *samp2_gfit = new TF1("", "pol3", 2, 10);
   samp2_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   samp2_g->Fit(samp2_gfit);
   TF1 *sampd2_gfit = new TF1("", "pol3", 2, 10);
   sampd2_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   sampd2_g->Fit(sampd2_gfit);
   TF1 *sampu2_gfit = new TF1("", "pol3", 2, 10);
   sampu2_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   sampu2_g->Fit(sampu2_gfit);
   
   TF1 *samp3_gfit = new TF1("", "pol3", 2, 10);
   samp3_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   samp3_g->Fit(samp3_gfit);
   TF1 *sampd3_gfit = new TF1("", "pol3", 2, 10);
   sampd3_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   sampd3_g->Fit(sampd3_gfit);
   TF1 *sampu3_gfit = new TF1("", "pol3", 2, 10);
   sampu3_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   sampu3_g->Fit(sampu3_gfit);
   
   TF1 *samp4_gfit = new TF1("", "pol3", 2, 10);
   samp4_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   samp4_g->Fit(samp4_gfit);
   TF1 *sampd4_gfit = new TF1("", "pol3", 2, 10);
   sampd4_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   sampd4_g->Fit(sampd4_gfit);
   TF1 *sampu4_gfit = new TF1("", "pol3", 2, 10);
   sampu4_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   sampu4_g->Fit(sampu4_gfit);
   
   TF1 *samp5_gfit = new TF1("", "pol3", 2, 10);
   samp5_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   samp5_g->Fit(samp5_gfit);
   TF1 *sampd5_gfit = new TF1("", "pol3", 2, 10);
   sampd5_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   sampd5_g->Fit(sampd5_gfit);
   TF1 *sampu5_gfit = new TF1("", "pol3", 2, 10);
   sampu5_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   sampu5_g->Fit(sampu5_gfit);
   
   TF1 *samp6_gfit = new TF1("", "pol3", 2, 10);
   samp6_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   samp6_g->Fit(samp6_gfit);
   TF1 *sampd6_gfit = new TF1("", "pol3", 2, 10);
   sampd6_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   sampd6_g->Fit(sampd6_gfit);
   TF1 *sampu6_gfit = new TF1("", "pol3", 2, 10);
   sampu6_gfit->SetParameters(0.2, 0.004, -0.0007, 0.000017);
   sampu6_g->Fit(sampu6_gfit);
   
   
   samp1_g->Draw();
   samp1_g->Write();
   sampd1_g->Draw();
   sampd1_g->Write();
   sampu1_g->Draw();
   sampu1_g->Write();
   
   samp2_g->Draw();
   samp2_g->Write();
   sampd2_g->Draw();
   sampd2_g->Write();
   sampu2_g->Draw();
   sampu2_g->Write();
   
   samp3_g->Draw();
   samp3_g->Write();
   sampd3_g->Draw();
   sampd3_g->Write();
   sampu3_g->Draw();
   sampu3_g->Write();
   
   samp4_g->Draw();
   samp4_g->Write();
   sampd4_g->Draw();
   sampd4_g->Write();
   sampu4_g->Draw();
   sampu4_g->Write();
   
   samp5_g->Draw();
   samp5_g->Write();
   sampd5_g->Draw();
   sampd5_g->Write();
   sampu5_g->Draw();
   sampu5_g->Write();
   
   samp6_g->Draw();
   samp6_g->Write();
   sampd6_g->Draw();
   sampd6_g->Write();
   sampu6_g->Draw();
   sampu6_g->Write();
   
   
   
     auto samv_g = new TGraph();
      samv_g->SetName("samv_g");
      samv_g->SetTitle("samv mean");
      samv_g->SetMarkerColor(kRed);
      samv_g->SetMarkerStyle(kFullCircle);
      samv_g->SetLineWidth(0);
      
  auto samvu_g = new TGraph();
      samvu_g->SetName("samvu_g");
      samvu_g->SetTitle("samv up mean");
      samvu_g->SetMarkerColor(kYellow );
      samvu_g->SetMarkerStyle(kFullCircle);
      samvu_g->SetLineWidth(0);
    
   auto samvd_g = new TGraph();
      samvd_g->SetName("samvd_g");
      samvd_g->SetTitle("samv down mean");
      samvd_g->SetMarkerColor(kYellow );
      samvd_g->SetMarkerStyle(kFullCircle);
      samvd_g->SetLineWidth(0);
      
      
   float lwpv = 0.0;
   float upv = 0.0; 
   
   
   /*
   for(int i = 0; i < 50; i++){
       upv = lwpv + 4.0;
       TH1F *samv1_hist = (TH1F*)samv_hist->ProjectionY("", lwpv, upv);
       TF1 *samv_fit = new TF1("samv_fit", "gaus(0)", 0.05, 0.35);
          //  samv_fit->SetParameters(0.01, 0.24, 0.03);
            samv_fit->SetParNames("const", "mean", "sigma");
       samv1_hist->Fit(samv_fit, "R");
       
       
       
       samv_g->AddPoint( (70.0/200.0)  * ((lwpv + upv)/2.0), samv_fit->GetParameter(1));                       //samp6_gfit->GetParameter(0)
       samvd_g->AddPoint( (70.0/200.0)  * ((lwpv + upv)/2.0), samv_fit->GetParameter(1) - (3 * samv_fit->GetParameter(2)));
       samvu_g->AddPoint( (70.0/200.0)  * ((lwpv + upv)/2.0), samv_fit->GetParameter(1) + (3 * samv_fit->GetParameter(2)));
       
       lwpv = lwpv + 4.0;
   
   }
   */
   
 //  TF1 *samv_gfit = new TF1("", "pol3", 0, 70);
 //  samv_g->Fit(samv_gfit);
   samv_g->Draw();
   samv_g->Write();
 //  TF1 *samvd_gfit = new TF1("", "pol3", 0, 70);
 //  samvd_g->Fit(samvd_gfit);
   samvd_g->Draw();
   samvd_g->Write();
 //  TF1 *samvu_gfit = new TF1("", "pol3", 0, 70);
 //  samvu_g->Fit(samvu_gfit);
   samvu_g->Draw();
   samvu_g->Write();
   
   
   
     auto samw_g = new TGraph();
      samw_g->SetName("samw_g");
      samw_g->SetTitle("samw mean");
      samw_g->SetMarkerColor(kRed);
      samw_g->SetMarkerStyle(kFullCircle);
      samw_g->SetLineWidth(0);
      
  auto samwu_g = new TGraph();
      samwu_g->SetName("samwu_g");
      samwu_g->SetTitle("samw up mean");
      samwu_g->SetMarkerColor(kYellow );
      samwu_g->SetMarkerStyle(kFullCircle);
      samwu_g->SetLineWidth(0);
    
   auto samwd_g = new TGraph();
      samwd_g->SetName("samwd_g");
      samwd_g->SetTitle("samw down mean");
      samwd_g->SetMarkerColor(kYellow );
      samwd_g->SetMarkerStyle(kFullCircle);
      samwd_g->SetLineWidth(0);
      
      
   float lwpw = 0.0;
   float upw = 0.0; 
   /*
   for(int i = 0; i < 50; i++){
       upw = lwpw + 4.0;
       TH1F *samw1_hist = (TH1F*)samw_hist->ProjectionY("", lwpw, upw);
       TF1 *samw_fit = new TF1("samw_fit", "gaus(0)", 0.05, 0.35);
         //   samw_fit->SetParameters(0.01, 0.24, 0.03);
            samw_fit->SetParNames("const", "mean", "sigma");
       samw1_hist->Fit(samw_fit, "R");
       
       samw_g->AddPoint( (70.0/200.0)  * ((lwpw + upw)/2.0), samw_fit->GetParameter(1));
       samwd_g->AddPoint( (70.0/200.0)  * ((lwpw + upw)/2.0), samw_fit->GetParameter(1) - (3 * samw_fit->GetParameter(2)));
       samwu_g->AddPoint( (70.0/200.0)  * ((lwpw + upw)/2.0), samw_fit->GetParameter(1) + (3 * samw_fit->GetParameter(2)));
       
       lwpw = lwpw + 4.0;
   
   }
   */
   
 //  TF1 *samw_gfit = new TF1("", "pol3", 0, 70);
 //  samw_g->Fit(samw_gfit);
   samw_g->Draw();
   samw_g->Write();
//   TF1 *samwd_gfit = new TF1("", "pol3", 0, 70);
//   samwd_g->Fit(samwd_gfit);
   samwd_g->Draw();
   samwd_g->Write();
 //  TF1 *samwu_gfit = new TF1("", "pol3", 0, 70);
 //  samwu_g->Fit(samwu_gfit);
   samwu_g->Draw();
   samwu_g->Write();
   
   thve_hist->GetXaxis()->SetTitle("edge [cm]");
   thve_hist->GetYaxis()->SetTitle("#theta [deg]");
   thve_hist->Draw();
   thve_hist->Write();
   
   
   
   
   TF1 *chin1_fit = new TF1("chin1_fit", "pol0", 12, 30);
       chin_hist1->Fit(chin1_fit, "R");
       
   TF1 *chin2_fit = new TF1("chin2_fit", "pol0", 12, 30);
       chin_hist2->Fit(chin2_fit, "R");
       
   TF1 *chin3_fit = new TF1("chin3_fit", "pol0", 12, 30);
       chin_hist3->Fit(chin3_fit, "R");
       
   TF1 *chin4_fit = new TF1("chin4_fit", "pol0", 12, 30);
       chin_hist4->Fit(chin4_fit, "R");
       
   TF1 *chin5_fit = new TF1("chin5_fit", "pol0", 12, 30);
       chin_hist5->Fit(chin5_fit, "R");
       
   TF1 *chin6_fit = new TF1("chin6_fit", "pol0", 12, 30);
       chin_hist6->Fit(chin6_fit, "R");
       
       float chin1_std = chin1_fit->GetParError(0) * 50;
       float chin2_std = chin2_fit->GetParError(0) * 50;
       float chin3_std = chin3_fit->GetParError(0) * 50;
       float chin4_std = chin4_fit->GetParError(0) * 50;
       float chin5_std = chin5_fit->GetParError(0) * 50;
       float chin6_std = chin6_fit->GetParError(0) * 50;
       
       std::cout << chin1_std << std::endl;
       /*
       float chih1 = chin1_std + chin1_fit->GetParameter(0);
       float chih2 = chin2_std + chin2_fit->GetParameter(0);
       float chih3 = chin3_std + chin3_fit->GetParameter(0);
       float chih4 = chin4_std + chin4_fit->GetParameter(0);
       float chih5 = chin5_std + chin5_fit->GetParameter(0);
       float chih6 = chin6_std + chin6_fit->GetParameter(0);
       */
       int chin1c = 0;
       int chin2c = 0;
       int chin3c = 0;
       int chin4c = 0;
       int chin5c = 0;
       int chin6c = 0;
       
       float chin1_cut = 0.0;
       float chin2_cut = 0.0;
       float chin3_cut = 0.0;
       float chin4_cut = 0.0;
       float chin5_cut = 0.0;
       float chin6_cut = 0.0;
       
       float chin1t = 0.0;
       float chin2t = 0.0;
       float chin3t = 0.0;
       float chin4t = 0.0;
       float chin5t = 0.0;
       float chin6t = 0.0;
       
       for(int i = 100; i < 200; i++){
       
       float chin1i = std::abs(chin1_fit->GetParameter(0) - chin_hist1->GetBinContent(i));
       float chin2i = std::abs(chin2_fit->GetParameter(0) - chin_hist2->GetBinContent(i));
       float chin3i = std::abs(chin3_fit->GetParameter(0) - chin_hist3->GetBinContent(i));
       float chin4i = std::abs(chin4_fit->GetParameter(0) - chin_hist4->GetBinContent(i));
       float chin5i = std::abs(chin5_fit->GetParameter(0) - chin_hist5->GetBinContent(i));
       float chin6i = std::abs(chin6_fit->GetParameter(0) - chin_hist6->GetBinContent(i));
       
       if(chin1t < chin1i){
          chin1t = chin1i;
       }
       
       if(chin2t < chin2i){
          chin2t = chin2i;
       }
       
       if(chin3t < chin3i){
          chin3t = chin3i;
       }
       
       if(chin4t < chin4i){
          chin4t = chin4i;
       }
       
       if(chin5t < chin5i){
          chin5t = chin5i;
       }
       
       if(chin6t < chin6i){
          chin6t = chin6i;
       }
       
       }
       
       float chih1 = chin1t + chin1_fit->GetParameter(0);
       float chih2 = chin2t + chin2_fit->GetParameter(0);
       float chih3 = chin3t + chin3_fit->GetParameter(0);
       float chih4 = chin4t + chin4_fit->GetParameter(0);
       float chih5 = chin5t + chin5_fit->GetParameter(0);
       float chih6 = chin6t + chin6_fit->GetParameter(0);
       
       for(int i = 0; i < 95; i++){
       
       if( chin_hist1->GetBinContent(i) < chih1 && chin_hist1->GetBinContent(i+1) < chih1 && chin_hist1->GetBinContent(i+2) < chih1 && chin_hist1->GetBinContent(i+3) < chih1 && chin_hist1->GetBinContent(i+4) < chih1 && chin1c == 0  ){
       
       chin1_cut = (i) * (30.0/200.0);
       
       chin1c++;
       }
       
       
       if( chin_hist2->GetBinContent(i) < chih2 && chin_hist2->GetBinContent(i+1) < chih2 && chin_hist2->GetBinContent(i+2) < chih2 && chin_hist2->GetBinContent(i+3) < chih2 && chin_hist2->GetBinContent(i+4) < chih2 && chin2c == 0  ){
       
       chin2_cut = (i)* (30.0/200.0);
       
       chin2c++;
       }
       
       
       if( chin_hist3->GetBinContent(i) < chih3 && chin_hist3->GetBinContent(i+1) < chih3 && chin_hist3->GetBinContent(i+2) < chih3 && chin_hist3->GetBinContent(i+3) < chih3 && chin_hist3->GetBinContent(i+4) < chih3 && chin3c == 0  ){
       
       chin3_cut = (i)* (30.0/200.0);
       
       chin3c++;
       }
       
       
       if( chin_hist4->GetBinContent(i) < chih4 && chin_hist4->GetBinContent(i+1) < chih4 && chin_hist4->GetBinContent(i+2) < chih4 && chin_hist4->GetBinContent(i+3) < chih4 && chin_hist4->GetBinContent(i+4) < chih4 && chin4c == 0  ){
       
       chin4_cut = (i)* (30.0/200.0);
       
       chin4c++;
       }
       
       
       if( chin_hist5->GetBinContent(i) < chih5 && chin_hist5->GetBinContent(i+1) < chih5 && chin_hist5->GetBinContent(i+2) < chih5 && chin_hist5->GetBinContent(i+3) < chih5 && chin_hist5->GetBinContent(i+4) < chih5 && chin5c == 0  ){
       
       chin5_cut = (i)* (30.0/200.0);
       
       chin5c++;
       }
       
       
       if( chin_hist6->GetBinContent(i) < chih6 && chin_hist6->GetBinContent(i+1) < chih6 && chin_hist6->GetBinContent(i+2) < chih6 && chin_hist6->GetBinContent(i+3) < chih6 && chin_hist6->GetBinContent(i+4) < chih6 && chin6c == 0  ){
       
       chin6_cut = (i)* (30.0/200.0);
       
       chin6c++;
       }
       
           
       
       
       }
       
       TLine *chin1_x = new TLine(0.1, chin1_fit->GetParameter(0), 30, chin1_fit->GetParameter(0));
       TLine *chin2_x = new TLine(0.1, chin2_fit->GetParameter(0), 30, chin2_fit->GetParameter(0));
       TLine *chin3_x = new TLine(0.1, chin3_fit->GetParameter(0), 30, chin3_fit->GetParameter(0));
       TLine *chin4_x = new TLine(0.1, chin4_fit->GetParameter(0), 30, chin4_fit->GetParameter(0));
       TLine *chin5_x = new TLine(0.1, chin5_fit->GetParameter(0), 30, chin5_fit->GetParameter(0));
       TLine *chin6_x = new TLine(0.1, chin6_fit->GetParameter(0), 30, chin6_fit->GetParameter(0));
       
       TLine *chin1_x2 = new TLine(0.1, chih1, 30, chih1);
       TLine *chin2_x2 = new TLine(0.1, chih2, 30, chih2);
       TLine *chin3_x2 = new TLine(0.1, chih3, 30, chih3);
       TLine *chin4_x2 = new TLine(0.1, chih4, 30, chih4);
       TLine *chin5_x2 = new TLine(0.1, chih5, 30, chih5);
       TLine *chin6_x2 = new TLine(0.1, chih6, 30, chih6);
       
       TLine *chin1_y = new TLine(chin1_cut, chin_hist1->GetYmin(), chin1_cut, chin_hist1->GetYmax());
       TLine *chin2_y = new TLine(chin2_cut, chin_hist2->GetYmin(), chin2_cut, chin_hist2->GetYmax());
       TLine *chin3_y = new TLine(chin3_cut, chin_hist3->GetYmin(), chin3_cut, chin_hist3->GetYmax());
       TLine *chin4_y = new TLine(chin4_cut, chin_hist4->GetYmin(), chin4_cut, chin_hist4->GetYmax());
       TLine *chin5_y = new TLine(chin5_cut, chin_hist5->GetYmin(), chin5_cut, chin_hist5->GetYmax());
       TLine *chin6_y = new TLine(chin6_cut, chin_hist6->GetYmin(), chin6_cut, chin_hist6->GetYmax());
       
       chin1_x->SetLineWidth(3);
       chin2_x->SetLineWidth(3);
       chin3_x->SetLineWidth(3);
       chin4_x->SetLineWidth(3);
       chin5_x->SetLineWidth(3);
       chin6_x->SetLineWidth(3);
       
       chin1_x2->SetLineWidth(3);
       chin2_x2->SetLineWidth(3);
       chin3_x2->SetLineWidth(3);
       chin4_x2->SetLineWidth(3);
       chin5_x2->SetLineWidth(3);
       chin6_x2->SetLineWidth(3);
       
       chin1_y->SetLineWidth(3);
       chin2_y->SetLineWidth(3);
       chin3_y->SetLineWidth(3);
       chin4_y->SetLineWidth(3);
       chin5_y->SetLineWidth(3);
       chin6_y->SetLineWidth(3);
       
       chin1_x->SetLineColor(kBlue);
       chin2_x->SetLineColor(kBlue);
       chin3_x->SetLineColor(kBlue);
       chin4_x->SetLineColor(kBlue);
       chin5_x->SetLineColor(kBlue);
       chin6_x->SetLineColor(kBlue);
       
       chin1_x2->SetLineColor(kRed);
       chin2_x2->SetLineColor(kRed);
       chin3_x2->SetLineColor(kRed);
       chin4_x2->SetLineColor(kRed);
       chin5_x2->SetLineColor(kRed);
       chin6_x2->SetLineColor(kRed);
       
       chin1_y->SetLineColor(kGreen+2);
       chin2_y->SetLineColor(kGreen+2);
       chin3_y->SetLineColor(kGreen+2);
       chin4_y->SetLineColor(kGreen+2);
       chin5_y->SetLineColor(kGreen+2);
       chin6_y->SetLineColor(kGreen+2);
       
   
   chin_hist1->GetXaxis()->SetTitle("edge [cm]");
   chin_hist1->GetYaxis()->SetTitle("<#chi^{2} / NDF>");
   chin_hist1->Draw("APL");
   chin_hist1->Write();
   
   chin_hist2->GetXaxis()->SetTitle("edge [cm]");
   chin_hist2->GetYaxis()->SetTitle("<#chi^{2} / NDF>");
   chin_hist2->Draw("APL");
   chin_hist2->Write();
   
   chin_hist3->GetXaxis()->SetTitle("edge [cm]");
   chin_hist3->GetYaxis()->SetTitle("<#chi^{2} / NDF>");
   chin_hist3->Draw("APL");
   chin_hist3->Write();
   
   chin_hist4->GetXaxis()->SetTitle("edge [cm]");
   chin_hist4->GetYaxis()->SetTitle("<#chi^{2} / NDF>");
   chin_hist4->Draw("APL");
   chin_hist4->Write();
   
   chin_hist5->GetXaxis()->SetTitle("edge [cm]");
   chin_hist5->GetYaxis()->SetTitle("<#chi^{2} / NDF>");
   chin_hist5->Draw("APL");
   chin_hist5->Write();
   
   chin_hist6->GetXaxis()->SetTitle("edge [cm]");
   chin_hist6->GetYaxis()->SetTitle("<#chi^{2} / NDF>");
   chin_hist6->Draw("APL");
   chin_hist6->Write();
   
   
   
    auto thve_g = new TGraph();
      thve_g->SetName("Sector Independent");
      thve_g->SetTitle("samw mean");
      thve_g->SetMarkerColor(kRed);
      thve_g->SetMarkerStyle(kFullCircle);
      thve_g->SetLineWidth(0);
      
  auto thveu_g = new TGraph();
      thveu_g->SetName("sector Independent");
      thveu_g->SetTitle("#mu + 3#sigma");
      thveu_g->SetMarkerColor(kBlack );
      thveu_g->SetMarkerStyle(kFullCircle);
      thveu_g->SetLineWidth(0);
      
      
      auto thveu_g1 = new TGraph();
      thveu_g1->SetName("Sector 1");
      thveu_g1->SetTitle("samw up mean");
      thveu_g1->SetMarkerColor(kRed );
      thveu_g1->SetMarkerStyle(kFullCircle);
      thveu_g1->SetLineWidth(0);
      
      
      auto thveu_g2 = new TGraph();
      thveu_g2->SetName("Sector 2");
      thveu_g2->SetTitle("samw up mean");
      thveu_g2->SetMarkerColor(kOrange );
      thveu_g2->SetMarkerStyle(kFullCircle);
      thveu_g2->SetLineWidth(0);
      
      
      auto thveu_g3 = new TGraph();
      thveu_g3->SetName("Sector 3");
      thveu_g3->SetTitle("samw up mean");
      thveu_g3->SetMarkerColor(kYellow );
      thveu_g3->SetMarkerStyle(kFullCircle);
      thveu_g3->SetLineWidth(0);
      
      
      auto thveu_g4 = new TGraph();
      thveu_g4->SetName("Sector 4");
      thveu_g4->SetTitle("samw up mean");
      thveu_g4->SetMarkerColor(kGreen );
      thveu_g4->SetMarkerStyle(kFullCircle);
      thveu_g4->SetLineWidth(0);
      
      
      auto thveu_g5 = new TGraph();
      thveu_g5->SetName("Sector 5");
      thveu_g5->SetTitle("samw up mean");
      thveu_g5->SetMarkerColor(kBlue );
      thveu_g5->SetMarkerStyle(kFullCircle);
      thveu_g5->SetLineWidth(0);
      
      
      auto thveu_g6 = new TGraph();
      thveu_g6->SetName("Sector 6");
      thveu_g6->SetTitle("samw up mean");
      thveu_g6->SetMarkerColor(kMagenta );
      thveu_g6->SetMarkerStyle(kFullCircle);
      thveu_g6->SetLineWidth(0);
    
   auto thved_g = new TGraph();
      thved_g->SetName("samwd_g");
      thved_g->SetTitle("samw down mean");
      thved_g->SetMarkerColor(kBlack );
      thved_g->SetMarkerStyle(kFullCircle);
      thved_g->SetLineWidth(0);
      
      
      
      
   float lwpt = 0.0;
   float upt = 0.0; 
   for(int i = 0; i < 30; i++){
       upt = lwpt + 6;
       TH1F *thve1_hist = (TH1F*)thve_hist->ProjectionX("edge", lwpt, upt);
       /*
       TH1F *thve1_hist1 = (TH1F*)thve_hist1->ProjectionX("edge1", lwpt, upt);
       TH1F *thve1_hist2 = (TH1F*)thve_hist2->ProjectionX("edge2", lwpt, upt);
       TH1F *thve1_hist3 = (TH1F*)thve_hist3->ProjectionX("edge3", lwpt, upt);
       TH1F *thve1_hist4 = (TH1F*)thve_hist4->ProjectionX("edge4", lwpt, upt);
       TH1F *thve1_hist5 = (TH1F*)thve_hist5->ProjectionX("edge5", lwpt, upt);
       TH1F *thve1_hist6 = (TH1F*)thve_hist6->ProjectionX("edge6", lwpt, upt);
       */
       
       TH1F *thve_grad_hist = new TH1F("thve_grad", "gradient", 200, 0, 60);
             	    thve_grad_hist->SetLineColor(kBlack);
             	    /*
       TH1F *thve_grad_hist1 = new TH1F("thve_grad1", "sector 1 gradient", 200, 0, 20);
                    thve_grad_hist1->SetLineColor(kRed);
       TH1F *thve_grad_hist2 = new TH1F("thve_grad2", "sector 2 gradient", 200, 0, 20);
                    thve_grad_hist2->SetLineColor(kOrange);
       TH1F *thve_grad_hist3 = new TH1F("thve_grad3", "sector 3 gradient", 200, 0, 20);
                    thve_grad_hist3->SetLineColor(kYellow);
       TH1F *thve_grad_hist4 = new TH1F("thve_grad4", "sector 4 gradient", 200, 0, 20);
                    thve_grad_hist4->SetLineColor(kGreen);
       TH1F *thve_grad_hist5 = new TH1F("thve_grad5", "sector 5 gradient", 200, 0, 20);
                    thve_grad_hist5->SetLineColor(kBlue);
       TH1F *thve_grad_hist6 = new TH1F("thve_grad6", "sector 6 gradient", 200, 0, 20);
                    thve_grad_hist6->SetLineColor(kMagenta);
       */ 
       for(int k = 1; k < 199; k++){
       
       //if(thve1_hist->GetBinContent(k+1) - thve1_hist->GetBinContent(k)){
       double slope = (thve1_hist->GetBinContent(k+1) - thve1_hist->GetBinContent(k)) / ((60.0/200.0)*(k+1) - (60.0/200.0)*k);
       //}
       thve_grad_hist->SetBinContent(k, slope);
       
       }
       
      // TH1F *thve_grad_vhist = new TH1F("thve_grad", "gradient", 200, 0, 20.0);
       //TH1F *thvel_grad_hist = (TH1F*)thve1_hist->GetGradient();
       
       //thve_grad_hist->Scale(1.0 / thve_grad_hist->Integral());
       /*
       thve_grad_hist1->Scale(1.0 / thve_grad_hist1->Integral());
       thve_grad_hist2->Scale(1.0 / thve_grad_hist2->Integral());
       thve_grad_hist3->Scale(1.0 / thve_grad_hist3->Integral());
       thve_grad_hist4->Scale(1.0 / thve_grad_hist4->Integral());
       thve_grad_hist5->Scale(1.0 / thve_grad_hist5->Integral());
       thve_grad_hist6->Scale(1.0 / thve_grad_hist6->Integral());
       */ 
       
       TF1 *thve_fit = new TF1("thve_fit", "gaus(0)", 0, 4);
            //thve_fit->SetParameters(1, 3);
            thve_fit->SetParNames("const", "mean", "sigma");
       thve_grad_hist->Fit(thve_fit, "R");
      /* 
       TF1 *thve_fit1 = new TF1("thve_fit1", "gaus(0)", 0, 20);
            thve_fit1->SetParameters(1, 3);
            thve_fit1->SetParNames("const", "mean", "sigma");
       thve_grad_hist1->Fit(thve_fit1, "R");
       
       TF1 *thve_fit2 = new TF1("thve_fit2", "gaus(0)", 0, 20);
            thve_fit2->SetParameters(1, 3);
            thve_fit2->SetParNames("const", "mean", "sigma");
       thve_grad_hist2->Fit(thve_fit2, "R");
       
       TF1 *thve_fit3 = new TF1("thve_fit3", "gaus(0)", 0, 20);
            thve_fit3->SetParameters(1, 3);
            thve_fit3->SetParNames("const", "mean", "sigma");
       thve_grad_hist3->Fit(thve_fit3, "R");
       
       TF1 *thve_fit4 = new TF1("thve_fit4", "gaus(0)", 0, 20);
            thve_fit4->SetParameters(1, 3);
            thve_fit4->SetParNames("const", "mean", "sigma");
       thve_grad_hist4->Fit(thve_fit4, "R");
       
       TF1 *thve_fit5 = new TF1("thve_fit5", "gaus(0)", 0, 20);
            thve_fit5->SetParameters(1, 3);
            thve_fit5->SetParNames("const", "mean", "sigma");
       thve_grad_hist5->Fit(thve_fit5, "R");
       
       TF1 *thve_fit6 = new TF1("thve_fit6", "gaus(0)", 0, 20);
            thve_fit6->SetParameters(1, 3);
            thve_fit6->SetParNames("const", "mean", "sigma");
       thve_grad_hist6->Fit(thve_fit6, "R");
       */
      // thve_grad_hist->Draw();
      // thve_grad_hist->Write();
      
      //if ((30.0/200.0)  * ((lwpt + upt)/2.0) + 5.0 <= 21 && thve_fit->GetParameter(1) + (3 * thve_fit->GetParameter(2)) > 0 && thve_fit->GetParameter(1) + (3 * thve_fit->GetParameter(2)) < 10 ){ //&& thve_fit2->GetParameter(1) + (3 * thve_fit2->GetParameter(2)) > 0 && thve_fit2->GetParameter(1) + (3 * thve_fit2->GetParameter(2)) < 10 && thve_fit3->GetParameter(1) + (3 * thve_fit3->GetParameter(2)) > 0 && thve_fit3->GetParameter(1) + (3 * thve_fit3->GetParameter(2)) < 10){
      if((thve_fit->GetParameter(1) + (3 * thve_fit->GetParameter(2))) < 4 && (thve_fit->GetParameter(1) + (3 * thve_fit->GetParameter(2))) > 2){
       thve_g->AddPoint( thve_fit->GetParameter(1),  (30.0/200.0)  * ((lwpt + upt)/2.0) + 5.0);
       thved_g->AddPoint( thve_fit->GetParameter(1) - (3 * thve_fit->GetParameter(2)) , (30.0/200.0)  * ((lwpt + upt)/2.0) + 5.0 );
       thveu_g->AddPoint( thve_fit->GetParameter(1) + (3 * thve_fit->GetParameter(2)) , (30.0/200.0)  * ((lwpt + upt)/2.0) + 5.0);
       }
       /*
       thveu_g1->AddPoint( (30.0/200.0)  * ((lwpt + upt)/2.0) + 5.0 ,  thve_fit1->GetParameter(1) + (3 * thve_fit1->GetParameter(2)));
       thveu_g2->AddPoint( (30.0/200.0)  * ((lwpt + upt)/2.0) + 5.0 ,  thve_fit2->GetParameter(1) + (3 * thve_fit2->GetParameter(2)));
       thveu_g3->AddPoint( (30.0/200.0)  * ((lwpt + upt)/2.0) + 5.0 ,  thve_fit3->GetParameter(1) + (3 * thve_fit3->GetParameter(2)));
       thveu_g4->AddPoint( (30.0/200.0)  * ((lwpt + upt)/2.0) + 5.0 ,  thve_fit4->GetParameter(1) + (3 * thve_fit4->GetParameter(2)));
       thveu_g5->AddPoint( (30.0/200.0)  * ((lwpt + upt)/2.0) + 5.0 ,  thve_fit5->GetParameter(1) + (3 * thve_fit5->GetParameter(2)));
       thveu_g6->AddPoint( (30.0/200.0)  * ((lwpt + upt)/2.0) + 5.0 ,  thve_fit6->GetParameter(1) + (3 * thve_fit6->GetParameter(2)));
       */
       //std::cout << thve_fit1->GetParameter(1) << " | " << thve_fit2->GetParameter(1) << " | " << thve_fit3->GetParameter(1) << " | " << thve_fit4->GetParameter(1) << std::endl;
       //}
       
       lwpt = lwpt + 6.0;
      // if(lwpt == 80){
      thve_fit->Draw();
      thve_fit->Write();
      
      //thve1_hist->Scale(100.0 / thve1_hist->Integral());
      
      //thve_grad_hist->Scale(1.0 / thve_grad_hist->Integral());
      
       //TCanvas *thve_c = new TCanvas("thve_c","thve_c",0,0,640,480);
       thve1_hist->Draw();
       thve1_hist->Write();
       thve_grad_hist->Draw();
       thve_grad_hist->Write();
       //thve_grad_hist1->Draw("HIST Same");
       //thve_grad_hist2->Draw("HIST Same");
       //thve_grad_hist3->Draw("HIST Same");
       //thve_grad_hist4->Draw("HIST Same");
       //thve_grad_hist5->Draw("HIST Same");
       //thve_grad_hist6->Draw("HIST Same");
       //thve_fit->Draw();
       //thve_c->BuildLegend();
       //thve_c->Write();
   //}
   
   
   }
   
   TF1 *thveu_fit = new TF1("thve_fit", "pol1", 2, 4);
            thveu_fit->SetParNames("const", "mean");
       thveu_g->Fit(thveu_fit);
   
   //TF1 *samw_gfit = new TF1("", "pol3", 0, 70);
   //samw_g->Fit(samw_gfit);
   thve_g->Draw();
   thve_g->Write();
   //TF1 *samwd_gfit = new TF1("", "pol3", 0, 70);
   //samwd_g->Fit(samwd_gfit);
   thved_g->Draw();
   thved_g->Write();
   //TF1 *samwu_gfit = new TF1("", "pol3", 0, 70);
   //samwu_g->Fit(samwu_gfit);
   thveu_g->Draw();
   thveu_g->Write();
   
   xy_uncut1->GetXaxis()->SetTitle("x [cm]");
   xy_uncut1->GetYaxis()->SetTitle("y [cm]");
   xy_uncut1->Draw();
   xy_uncut1->Write();
   
   
   xy_cut1->GetXaxis()->SetTitle("x [cm]");
   xy_cut1->GetYaxis()->SetTitle("y [cm]");
   xy_cut1->Draw();
   xy_cut1->Write();
   
   xy_cut2->GetXaxis()->SetTitle("x [cm]");
   xy_cut2->GetYaxis()->SetTitle("y [cm]");
   xy_cut2->Draw();
   xy_cut2->Write();
   
   
   
   
   lc_hist->GetXaxis()->SetTitle("lc");
   lc_hist->Draw();
   lc_hist->Write();
   
   t_hist->GetXaxis()->SetTitle("-t");
   t_hist->Draw();
   t_hist->Write();
   
   w_hist->GetXaxis()->SetTitle("W");
   w_hist->Draw();
   w_hist->Write();
   
   z_hist->GetXaxis()->SetTitle("Zh");
   z_hist->Draw();
   z_hist->Write();

  
   rm_no->GetXaxis()->SetTitle("M_{#rho} [GeV]");
   rm_no->Draw();
   rm_no->Write();

   rm_w->GetXaxis()->SetTitle("M_{#rho} [GeV]");
   rm_w->Draw();
   rm_w->Write();

   rm_wt->GetXaxis()->SetTitle("M_{#rho} [GeV]");
   rm_wt->Draw();
   rm_wt->Write();

   rm_wtz->GetXaxis()->SetTitle("M_{#rho} [GeV]");
   rm_wtz->Draw();
   rm_wtz->Write();
   
   
   e_phi->GetXaxis()->SetTitle("e #phi [degree]");
   e_phi->Draw();
   e_phi->Write();
   
   pip_phi->GetXaxis()->SetTitle("#pi^{+} #phi [degree]");
   pip_phi->Draw();
   pip_phi->Write();
   
   pim_phi->GetXaxis()->SetTitle("#pi^{-} #phi [degree]");
   pim_phi->Draw();
   pim_phi->Write();
   
   e_theta->GetXaxis()->SetTitle("e #theta [degree]");
   e_theta->Draw();
   e_theta->Write();
   
   pip_theta->GetXaxis()->SetTitle("#pi^{+} #theta [degree]");
   pip_theta->Draw();
   pip_theta->Write();
   
   pim_theta->GetXaxis()->SetTitle("#pi^{-} #theta [degree]");
   pim_theta->Draw();
   pim_theta->Write();
   
   
   ep_hist->GetXaxis()->SetTitle("P_{e} [GeV]");
   ep_hist->Draw();
   ep_hist->Write();
   
   pipp_hist->GetXaxis()->SetTitle("P_{#Pi^{+}} [GeV]");
   pipp_hist->Draw();
   pipp_hist->Write();
   
   pimp_hist->GetXaxis()->SetTitle("P_{#Pi^{-}} [GeV]");
   pimp_hist->Draw();
   pimp_hist->Write();
   
   ewp_hist->GetXaxis()->SetTitle("P_{e} [GeV]");
   ewp_hist->Draw();
   ewp_hist->Write();
   
   pipwp_hist->GetXaxis()->SetTitle("P_{#Pi^{+}} [GeV]");
   pipwp_hist->Draw();
   pipwp_hist->Write();
   
   pimwp_hist->GetXaxis()->SetTitle("P_{#Pi^{-}} [GeV]");
   pimwp_hist->Draw();
   pimwp_hist->Write();
   
   ewtp_hist->GetXaxis()->SetTitle("P_{e} [GeV]");
   ewtp_hist->Draw();
   ewtp_hist->Write();
   
   pipwtp_hist->GetXaxis()->SetTitle("P_{#Pi^{+}} [GeV]");
   pipwtp_hist->Draw();
   pipwtp_hist->Write();
   
   pimwtp_hist->GetXaxis()->SetTitle("P_{#Pi^{-}} [GeV]");
   pimwtp_hist->Draw();
   pimwtp_hist->Write();
   
   ewtzp_hist->GetXaxis()->SetTitle("P_{e} [GeV]");
   ewtzp_hist->Draw();
   ewtzp_hist->Write();
   
   pipwtzp_hist->GetXaxis()->SetTitle("P_{#Pi^{+}} [GeV]");
   pipwtzp_hist->Draw();
   pipwtzp_hist->Write();
   
   pimwtzp_hist->GetXaxis()->SetTitle("P_{#Pi^{-}} [GeV]");
   pimwtzp_hist->Draw();
   pimwtzp_hist->Write();
   /*
   vzz_hist->SetLineColor(kBlack);
   fmt_vz_hist->SetLineColor(kGreen+2);
   
   vzz_hist->SetLineWidth(3);
   fmt_vz_hist->SetLineWidth(3);
   
   vzz_hist->Scale(1.0 / vzz_hist->Integral());
   fmt_vz_hist->Scale(1.0 / fmt_vz_hist->Integral());
   
   TCanvas *fmt_vzt = new TCanvas("fmt_vzt","fmt vz",0,0,640,480);
   vzz_hist->Draw("HIST");
   fmt_vz_hist->Draw("HIST Same");
   fmt_vzt->BuildLegend();
   fmt_vzt->Write();
   
   fmt_vz_hist->Draw();
   fmt_vz_hist->Write();
   */
   
   
      ivx_hist->Scale(1.0 / ivx_hist->Integral());
   //   ivx_hist->SetMinimum(0.0);
   //   ivx_hist->SetMaximum(0.19);
   TF1 *ivx_fit = new TF1("ivx_fit", "gaus(0)+pol4(3)",  - 1.5,   1.5);
   ivx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   ivx_fit->SetParNames("vx1 const", "vx1 mean", "vx1 sigma");
   ivx_hist->Fit(ivx_fit, "R");
   ivx_hist->Draw();
   ivx_hist->Write();
   iivx_hist->Scale(1.0 / iivx_hist->Integral());
      //iivx_hist->SetMinimum(0.0);
      //iivx_hist->SetMaximum(0.);
   TF1 *iivx_fit = new TF1("iivx_fit", "gaus(0)+pol4(3)", - 1.5, 1.5);
   iivx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   iivx_fit->SetParNames("vx2 const", "vx2 mean", "vx2 sigma");
   iivx_hist->Fit(iivx_fit, "R");
   iivx_hist->Draw();
   iivx_hist->Write();
   iiivx_hist->Scale(1.0 / iiivx_hist->Integral());
   //   iiivx_hist->SetMinimum(0.0);
   //   iiivx_hist->SetMaximum(0.19);
   TF1 *iiivx_fit = new TF1("iiivx_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   iiivx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   iiivx_fit->SetParNames("vx3 const", "vx3 mean", "vx3 sigma");
   iiivx_hist->Fit(iiivx_fit, "R");
   iiivx_hist->Draw();
   iiivx_hist->Write();
   ivvx_hist->Scale(1.0 / ivvx_hist->Integral());
    //  ivvx_hist->SetMinimum(0.0);
    //  ivvx_hist->SetMaximum(0.19);
   TF1 *ivvx_fit = new TF1("ivvx_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   ivvx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   ivvx_fit->SetParNames("vx4 const", "vx4 mean", "vx4 sigma");
   ivvx_hist->Fit(ivvx_fit, "R");
   ivvx_hist->Draw();
   ivvx_hist->Write();
   vvx_hist->Scale(1.0 / vvx_hist->Integral());
   //   vvx_hist->SetMinimum(0.0);
   //   vvx_hist->SetMaximum(0.19);
   TF1 *vvx_fit = new TF1("vvx_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   vvx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   vvx_fit->SetParNames("vx5 const", "vx5 mean", "vx5 sigma");
   vvx_hist->Fit(vvx_fit, "R");
   vvx_hist->Draw();
   vvx_hist->Write();
   vivx_hist->Scale(1.0 / vivx_hist->Integral());
   //   vivx_hist->SetMinimum(0.0);
   //   vivx_hist->SetMaximum(0.19);
   TF1 *vivx_fit = new TF1("vivx_fit", "gaus(0)+pol4(3)",  - 1.5, 1.5);
   vivx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   vivx_fit->SetParNames("vx6 const", "vx6 mean", "vx6 sigma");
   vivx_hist->Fit(vivx_fit, "R");
   vivx_hist->Draw();
   vivx_hist->Write();
   
   
   ivy_hist->Scale(1.0 / ivy_hist->Integral());
    //  ivy_hist->SetMinimum(0.0);
    //  ivy_hist->SetMaximum(0.12);
   TF1 *ivy_fit = new TF1("ivy_fit", "gaus(0)+pol4(3)",  - 1.5, 1.5);
   ivy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   ivy_fit->SetParNames("vy1 const", "vy1 mean", "vy1 sigma");
   ivy_hist->Fit(ivy_fit, "R");
   ivy_hist->Draw();
   ivy_hist->Write();
   iivy_hist->Scale(1.0 / iivy_hist->Integral());
    //  iivy_hist->SetMinimum(0.0);
    //  iivy_hist->SetMaximum(0.12);
   TF1 *iivy_fit = new TF1("iivy_fit", "gaus(0)+pol4(3)",  - 1.5,   1.5);
   iivy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   iivy_fit->SetParNames("vy2 const", "vy2 mean", "vy2 sigma");
   iivy_hist->Fit(iivy_fit, "R");
   iivy_hist->Draw();
   iivy_hist->Write();
   iiivy_hist->Scale(1.0 / iiivy_hist->Integral());
    //  iiivy_hist->SetMinimum(0.0);
    //  iiivy_hist->SetMaximum(0.12);
   TF1 *iiivy_fit = new TF1("iiivy_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   iiivy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   iiivy_fit->SetParNames("vy3 const", "vy3 mean", "vy3 sigma");
   iiivy_hist->Fit(iiivy_fit, "R");
   iiivy_hist->Draw();
   iiivy_hist->Write();
   ivvy_hist->Scale(1.0 / ivvy_hist->Integral());
    //  ivvy_hist->SetMinimum(0.0);
    //  ivvy_hist->SetMaximum(0.55);
   TF1 *ivvy_fit = new TF1("ivvy_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   ivvy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   ivvy_fit->SetParNames("vy4 const", "vy4 mean", "vy4 sigma");
   ivvy_hist->Fit(ivvy_fit, "R");
   ivvy_hist->Draw();
   ivvy_hist->Write();
   vvy_hist->Scale(1.0 / vvy_hist->Integral());
   //   vvy_hist->SetMinimum(0.0);
   //   vvy_hist->SetMaximum(0.12);
   TF1 *vvy_fit = new TF1("vvy_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   vvy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   vvy_fit->SetParNames("vy5 const", "vy5 mean", "vy5 sigma");
   vvy_hist->Fit(vvy_fit, "R");
   vvy_hist->Draw();
   vvy_hist->Write();
   vivy_hist->Scale(1.0 / vivy_hist->Integral());
    //  vivy_hist->SetMinimum(0.0);
    //  vivy_hist->SetMaximum(0.12);
   TF1 *vivy_fit = new TF1("vivy_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   vivy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   vivy_fit->SetParNames("vy6 const", "vy6 mean", "vy6 sigma");
   vivy_hist->Fit(vivy_fit, "R");
   vivy_hist->Draw();
   vivy_hist->Write();
   
   ivz_hist->Scale(1.0 / ivz_hist->Integral());
    //  ivz_hist->SetMinimum(0.0);
    //  ivz_hist->SetMaximum(0.04);
   ivz_hist->Draw("L");
   ivz_hist->Write();
   iivz_hist->Scale(1.0 / iivz_hist->Integral());
    //  iivz_hist->SetMinimum(0.0);
    //  iivz_hist->SetMaximum(0.04);
   iivz_hist->Draw();
   iivz_hist->Write();
   iiivz_hist->Scale(1.0 / iiivz_hist->Integral());
    //  iiivz_hist->SetMinimum(0.0);
    //  iiivz_hist->SetMaximum(0.04);
   iiivz_hist->Draw();
   iiivz_hist->Write();
   ivvz_hist->Scale(1.0 / ivvz_hist->Integral());
    //  ivvz_hist->SetMinimum(0.0);
    //  ivvz_hist->SetMaximum(0.04);
   ivvz_hist->Draw();
   ivvz_hist->Write();
   vvz_hist->Scale(1.0 / vvz_hist->Integral());
    //  vvz_hist->SetMinimum(0.0);
    //  vvz_hist->SetMaximum(0.04);
   vvz_hist->Draw();
   vvz_hist->Write();
   vivz_hist->Scale(1.0 / vivz_hist->Integral());
    //  vivz_hist->SetMinimum(0.0);
    //  vivz_hist->SetMaximum(0.04);
   vivz_hist->Draw();
   vivz_hist->Write();
   
   
   pip_ivx_hist->Scale(1.0 / pip_ivx_hist->Integral());
    //  pip_ivx_hist->SetMinimum(0.0);
    //  pip_ivx_hist->SetMaximum(0.14);
   TF1 *pip_ivx_fit = new TF1("pip_ivx_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   pip_ivx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pip_ivx_hist->Fit(pip_ivx_fit, "R");
   pip_ivx_hist->Draw();
   pip_ivx_hist->Write();
   pip_iivx_hist->Scale(1.0 / pip_iivx_hist->Integral());
    //  pip_iivx_hist->SetMinimum(0.0);
    //  pip_iivx_hist->SetMaximum(0.14);
   TF1 *pip_iivx_fit = new TF1("pip_iivx_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   pip_iivx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pip_iivx_hist->Fit(pip_iivx_fit, "R");
   pip_iivx_hist->Draw();
   pip_iivx_hist->Write();
   pip_iiivx_hist->Scale(1.0 / pip_iiivx_hist->Integral());
    //  pip_iiivx_hist->SetMinimum(0.0);
    //  pip_iiivx_hist->SetMaximum(0.14);
   TF1 *pip_iiivx_fit = new TF1("pip_iiivx_fit", "gaus(0)+pol4(3)",  - 2.5,  2.5);
   pip_iiivx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pip_iiivx_hist->Fit(pip_iiivx_fit, "R");
   pip_iiivx_hist->Draw();
   pip_iiivx_hist->Write();
   pip_ivvx_hist->Scale(1.0 / pip_ivvx_hist->Integral());
    //  pip_ivvx_hist->SetMinimum(0.0);
    //  pip_ivvx_hist->SetMaximum(0.14);
   TF1 *pip_ivvx_fit = new TF1("pip_ivvx_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   pip_ivvx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pip_ivvx_hist->Fit(pip_ivvx_fit, "R");
   pip_ivvx_hist->Draw();
   pip_ivvx_hist->Write();
   pip_vvx_hist->Scale(1.0 / pip_vvx_hist->Integral());
    //  pip_vvx_hist->SetMinimum(0.0);
    //  pip_vvx_hist->SetMaximum(0.14);
   TF1 *pip_vvx_fit = new TF1("pip_vvx_fit", "gaus(0)+pol4(3)",  - 2.0,  2.0);
   pip_vvx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pip_vvx_hist->Fit(pip_vvx_fit, "R");
   pip_vvx_hist->Draw();
   pip_vvx_hist->Write();
   pip_vivx_hist->Scale(1.0 / pip_vivx_hist->Integral());
    //  pip_vivx_hist->SetMinimum(0.0);
    //  pip_vivx_hist->SetMaximum(0.14);
   TF1 *pip_vivx_fit = new TF1("pip_vivx_fit", "gaus(0)+pol4(3)",  - 2.0,  2.0);
   pip_vivx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pip_vivx_hist->Fit(pip_vivx_fit, "R");
   pip_vivx_hist->Draw();
   pip_vivx_hist->Write();
   
   
   pip_ivy_hist->Scale(1.0 / pip_ivy_hist->Integral());
    //  pip_ivy_hist->SetMinimum(0.0);
    //  pip_ivy_hist->SetMaximum(0.1);
   TF1 *pip_ivy_fit = new TF1("pip_ivy_fit", "gaus(0)+pol3(3)",  - 1.5,  1.5);
   pip_ivy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1, 0.001);
   pip_ivy_hist->Fit(pip_ivy_fit, "R");
   pip_ivy_hist->Draw();
   pip_ivy_hist->Write();
   pip_iivy_hist->Scale(1.0 / pip_iivy_hist->Integral());
    //  pip_iivy_hist->SetMinimum(0.0);
    //  pip_iivy_hist->SetMaximum(0.1);
   TF1 *pip_iivy_fit = new TF1("pip_iivy_fit", "gaus(0)+pol4(3)",  - 1.5, 1.5);
   pip_iivy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pip_iivy_hist->Fit(pip_iivy_fit, "R");
   pip_iivy_hist->Draw();
   pip_iivy_hist->Write();
   pip_iiivy_hist->Scale(1.0 / pip_iiivy_hist->Integral());
    //  pip_iiivy_hist->SetMinimum(0.0);
    //  pip_iiivy_hist->SetMaximum(0.1);
   TF1 *pip_iiivy_fit = new TF1("pip_iiivy_fit", "gaus(0)+pol4(3)", - 1.5, 1.5);
   pip_iiivy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pip_iiivy_hist->Fit(pip_iiivy_fit, "R");
   pip_iiivy_hist->Draw();
   pip_iiivy_hist->Write();
   pip_ivvy_hist->Scale(1.0 / pip_ivvy_hist->Integral());
     // pip_ivvy_hist->SetMinimum(0.0);
     // pip_ivvy_hist->SetMaximum(0.1);
   TF1 *pip_ivvy_fit = new TF1("pip_ivvy_fit", "gaus(0)+pol4(3)",  - 2.5,  2.5);
   pip_ivvy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1, 0.001);
   pip_ivvy_hist->Fit(pip_ivvy_fit, "R");
   pip_ivvy_hist->Draw();
   pip_ivvy_hist->Write();
   pip_vvy_hist->Scale(1.0 / pip_vvy_hist->Integral());
     // pip_vvy_hist->SetMinimum(0.0);
     // pip_vvy_hist->SetMaximum(0.1);
   TF1 *pip_vvy_fit = new TF1("pip_vvy_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   pip_vvy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pip_vvy_hist->Fit(pip_vvy_fit, "R");
   pip_vvy_hist->Draw();
   pip_vvy_hist->Write();
   pip_vivy_hist->Scale(1.0 / pip_vivy_hist->Integral());
     // pip_vivy_hist->SetMinimum(0.0);
     // pip_vivy_hist->SetMaximum(0.1);
   TF1 *pip_vivy_fit = new TF1("pip_vivy_fit", "gaus(0)+pol4(3)",  - 1.5, 1.5);
   pip_vivy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pip_vivy_hist->Fit(pip_vivy_fit, "R");
   pip_vivy_hist->Draw();
   pip_vivy_hist->Write();
   
   pip_ivz_hist->Scale(1.0 / pip_ivz_hist->Integral());
   //   pip_ivz_hist->SetMinimum(0.0);
   //   pip_ivz_hist->SetMaximum(0.03);
   TF1 *pip_ivz_fit = new TF1("pip_ivz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
   pip_ivz_fit->SetParameters(0.01, -6, 1, 0.01, -3, 0.5);
   pip_ivz_hist->Fit(pip_ivz_fit, "R");
   pip_ivz_hist->Draw();
   pip_ivz_hist->Write();
   pip_iivz_hist->Scale(1.0 / pip_iivz_hist->Integral());
   //   pip_iivz_hist->SetMinimum(0.0);
   //   pip_iivz_hist->SetMaximum(0.03);
   TF1 *pip_iivz_fit = new TF1("pip_iivz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
   pip_iivz_fit->SetParameters(0.01, -6, 1, 0.01, -3, 1);
   pip_iivz_hist->Fit(pip_iivz_fit, "R");
   pip_iivz_hist->Draw();
   pip_iivz_hist->Write();
   pip_iiivz_hist->Scale(1.0 / pip_iiivz_hist->Integral());
   //   pip_iiivz_hist->SetMinimum(0.0);
   //   pip_iiivz_hist->SetMaximum(0.03);
   TF1 *pip_iiivz_fit = new TF1("pip_iiivz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
   pip_iiivz_fit->SetParameters(0.01, -6, 0.5, 0.01, -3, 1);
   pip_iiivz_hist->Fit(pip_iiivz_fit, "R");
   pip_iiivz_hist->Draw();
   pip_iiivz_hist->Write();
   pip_ivvz_hist->Scale(1.0 / pip_ivvz_hist->Integral());
   //   pip_ivvz_hist->SetMinimum(0.0);
   //   pip_ivvz_hist->SetMaximum(0.03);
   TF1 *pip_ivvz_fit = new TF1("pip_ivvz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
   pip_ivvz_fit->SetParameters(0.01, -6, 1, 0.01, -3, 1);
   pip_ivvz_hist->Fit(pip_ivvz_fit, "R");
   pip_ivvz_hist->Draw();
   pip_ivvz_hist->Write();
   pip_vvz_hist->Scale(1.0 / pip_vvz_hist->Integral());
   //   pip_vvz_hist->SetMinimum(0.0);
   //   pip_vvz_hist->SetMaximum(0.03);
   TF1 *pip_vvz_fit = new TF1("pip_vvz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
   pip_vvz_fit->SetParameters(0.01, -6, 1, 0.01, -3, 1);
   pip_vvz_hist->Fit(pip_vvz_fit, "R");
   pip_vvz_hist->Draw();
   pip_vvz_hist->Write();
   pip_vivz_hist->Scale(1.0 / pip_vivz_hist->Integral());
   //   pip_vivz_hist->SetMinimum(0.0);
   //   pip_vivz_hist->SetMaximum(0.03);
   TF1 *pip_vivz_fit = new TF1("pip_vivz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
   pip_vivz_fit->SetParameters(0.01, -6, 1, 0.01, -3, 1);
   pip_vivz_hist->Fit(pip_vivz_fit, "R");
   pip_vivz_hist->Draw();
   pip_vivz_hist->Write();
   
   
   pim_ivx_hist->Scale(1.0 / pim_ivx_hist->Integral());
   //   pim_ivx_hist->SetMinimum(0.0);
   //   pim_ivx_hist->SetMaximum(0.14);
   TF1 *pim_ivx_fit = new TF1("pim_ivx_fit", "gaus(0)+pol4(3)",  - 2.0,  2.0);
   pim_ivx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pim_ivx_hist->Fit(pim_ivx_fit, "R");
   pim_ivx_hist->Draw();
   pim_ivx_hist->Write();
   pim_iivx_hist->Scale(1.0 / pim_iivx_hist->Integral());
    //  pim_iivx_hist->SetMinimum(0.0);
    //  pim_iivx_hist->SetMaximum(0.14);
   TF1 *pim_iivx_fit = new TF1("pim_iivx_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   pim_iivx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pim_iivx_hist->Fit(pim_iivx_fit, "R");
   pim_iivx_hist->Draw();
   pim_iivx_hist->Write();
   pim_iiivx_hist->Scale(1.0 / pim_iiivx_hist->Integral());
    //  pim_iiivx_hist->SetMinimum(0.0);
    //  pim_iiivx_hist->SetMaximum(0.14);
   TF1 *pim_iiivx_fit = new TF1("pim_iiivx_fit", "gaus(0)+pol4(3)",  - 2.0,  2.0);
   pim_iiivx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pim_iiivx_hist->Fit(pim_iiivx_fit, "R");
   pim_iiivx_hist->Draw();
   pim_iiivx_hist->Write();
   pim_ivvx_hist->Scale(1.0 / pim_ivvx_hist->Integral());
    //  pim_ivvx_hist->SetMinimum(0.0);
    //  pim_ivvx_hist->SetMaximum(0.14);
   TF1 *pim_ivvx_fit = new TF1("pim_ivvx_fit", "gaus(0)+pol4(3)", - 2.0,  2.0);
   pim_ivvx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pim_ivvx_hist->Fit(pim_ivvx_fit, "R");
   pim_ivvx_hist->Draw();
   pim_ivvx_hist->Write();
   pim_vvx_hist->Scale(1.0 / pim_vvx_hist->Integral());
    //  pim_vvx_hist->SetMinimum(0.0);
    //  pim_vvx_hist->SetMaximum(0.14);
   TF1 *pim_vvx_fit = new TF1("pim_vvx_fit", "gaus(0)+pol4(3)",  - 1.5, 1.5);
   pim_vvx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pim_vvx_hist->Fit(pim_vvx_fit, "R");
   pim_vvx_hist->Draw();
   pim_vvx_hist->Write();
   pim_vivx_hist->Scale(1.0 / pim_vivx_hist->Integral());
    //  pim_vivx_hist->SetMinimum(0.0);
    //  pim_vivx_hist->SetMaximum(0.14);
   TF1 *pim_vivx_fit = new TF1("pim_vivx_fit", "gaus(0)+pol4(3)",  - 2.0,  2.0);
   pim_vivx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pim_vivx_hist->Fit(pim_vivx_fit, "R");
   pim_vivx_hist->Draw();
   pim_vivx_hist->Write();
   
   
   pim_ivy_hist->Scale(1.0 / pim_ivy_hist->Integral());
      //pim_ivy_hist->SetMinimum(0.0);
      //pim_ivy_hist->SetMaximum(0.1);
   TF1 *pim_ivy_fit = new TF1("pim_ivy_fit", "gaus(0)+pol4(3)",  - 2.0,  2.0);
   pim_ivy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1, 0.001, 0.001);
   pim_ivy_hist->Fit(pim_ivy_fit, "R");
   pim_ivy_hist->Draw();
   pim_ivy_hist->Write();
   pim_iivy_hist->Scale(1.0 / pim_iivy_hist->Integral());
   //   pim_iivy_hist->SetMinimum(0.0);
   //   pim_iivy_hist->SetMaximum(0.1);
   TF1 *pim_iivy_fit = new TF1("pim_iivy_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   pim_iivy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pim_iivy_hist->Fit(pim_iivy_fit, "R");
   pim_iivy_hist->Draw();
   pim_iivy_hist->Write();
   pim_iiivy_hist->Scale(1.0 / pim_iiivy_hist->Integral());
    //  pim_iiivy_hist->SetMinimum(0.0);
   //   pim_iiivy_hist->SetMaximum(0.1);
   TF1 *pim_iiivy_fit = new TF1("pim_iiivy_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   pim_iiivy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pim_iiivy_hist->Fit(pim_iiivy_fit, "R");
   pim_iiivy_hist->Draw();
   pim_iiivy_hist->Write();
   pim_ivvy_hist->Scale(1.0 / pim_ivvy_hist->Integral());
    //  pim_ivvy_hist->SetMinimum(0.0);
    //  pim_ivvy_hist->SetMaximum(0.1);
   TF1 *pim_ivvy_fit = new TF1("pim_ivvy_fit", "gaus(0)+pol4(3)",  - 2.0,  2.5);
   pim_ivvy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1, 0.01, 0.01);
   pim_ivvy_hist->Fit(pim_ivvy_fit, "R");
   pim_ivvy_hist->Draw();
   pim_ivvy_hist->Write();
   pim_vvy_hist->Scale(1.0 / pim_vvy_hist->Integral());
   //   pim_vvy_hist->SetMinimum(0.0);
   //   pim_vvy_hist->SetMaximum(0.1);
   TF1 *pim_vvy_fit = new TF1("pim_vvy_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   pim_vvy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pim_vvy_hist->Fit(pim_vvy_fit, "R");
   pim_vvy_hist->Draw();
   pim_vvy_hist->Write();
   pim_vivy_hist->Scale(1.0 / pim_vivy_hist->Integral());
   //   pim_vivy_hist->SetMinimum(0.0);
   //   pim_vivy_hist->SetMaximum(0.1);
   TF1 *pim_vivy_fit = new TF1("pim_vivy_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   pim_vivy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pim_vivy_hist->Fit(pim_vivy_fit, "R");
   pim_vivy_hist->Draw();
   pim_vivy_hist->Write();
   
   
   pim_ivz_hist->Scale(1.0 / pim_ivz_hist->Integral());
    //  pim_ivz_hist->SetMinimum(0.0);
    //  pim_ivz_hist->SetMaximum(0.03);
    TF1 *pim_ivz_fit = new TF1("pim_ivz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
   pim_ivz_fit->SetParameters(0.01, -8, 1, 0.01, -3, 1);
   pim_ivz_hist->Fit(pim_ivz_fit, "R");
   pim_ivz_hist->Draw();
   pim_ivz_hist->Write();
   pim_iivz_hist->Scale(1.0 / pim_iivz_hist->Integral());
    //  pim_iivz_hist->SetMinimum(0.0);
    //  pim_iivz_hist->SetMaximum(0.03);
     TF1 *pim_iivz_fit = new TF1("pim_iivz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
   pim_iivz_fit->SetParameters(0.01, -8, 1, 0.01, -3, 1);
   pim_iivz_hist->Fit(pim_iivz_fit, "R");
   pim_iivz_hist->Draw();
   pim_iivz_hist->Write();
   pim_iiivz_hist->Scale(1.0 / pim_iiivz_hist->Integral());
    //  pim_iiivz_hist->SetMinimum(0.0);
    //  pim_iiivz_hist->SetMaximum(0.03);
     TF1 *pim_iiivz_fit = new TF1("pim_iiivz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
   pim_iiivz_fit->SetParameters(0.01, -8, 1, 0.01, -3, 1);
   pim_iiivz_hist->Fit(pim_iiivz_fit, "R");
   pim_iiivz_hist->Draw();
   pim_iiivz_hist->Write();
   pim_ivvz_hist->Scale(1.0 / pim_ivvz_hist->Integral());
    //  pim_ivvz_hist->SetMinimum(0.0);
    //  pim_ivvz_hist->SetMaximum(0.03);
     TF1 *pim_ivvz_fit = new TF1("pim_ivvz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
   pim_ivvz_fit->SetParameters(0.01, -8, 1, 0.01, -3, 1);
   pim_ivvz_hist->Fit(pim_ivvz_fit, "R");
   pim_ivvz_hist->Draw();
   pim_ivvz_hist->Write();
   pim_vvz_hist->Scale(1.0 / pim_vvz_hist->Integral());
    //  pim_vvz_hist->SetMinimum(0.0);
    //  pim_vvz_hist->SetMaximum(0.03);
     TF1 *pim_vvz_fit = new TF1("pim_vvz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
   pim_vvz_fit->SetParameters(0.01, -8, 1, 0.01, -3, 1);
   pim_vvz_hist->Fit(pim_vvz_fit, "R");
   pim_vvz_hist->Draw();
   pim_vvz_hist->Write();
   pim_vivz_hist->Scale(1.0 / pim_vivz_hist->Integral());
     // pim_vivz_hist->SetMinimum(0.0);
     // pim_vivz_hist->SetMaximum(0.03);
      TF1 *pim_vivz_fit = new TF1("pim_vivz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
   pim_vivz_fit->SetParameters(0.01, -8, 1, 0.01, -3, 1);
   pim_vivz_hist->Fit(pim_vivz_fit, "R");
   pim_vivz_hist->Draw();
   pim_vivz_hist->Write();
   
   
      vx_hist->Add(ivx_hist);
      vx_hist->Add(iivx_hist);
      vx_hist->Add(iiivx_hist);
      vx_hist->Add(ivvx_hist);
      vx_hist->Add(vvx_hist);
      vx_hist->Add(vivx_hist);
      vx_hist->Scale(1.0 / vx_hist->Integral());
    //  vx_hist->SetMinimum(0.0);
    //  vx_hist->SetMaximum(0.19);
      TF1 *vx_fit = new TF1("vx_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   vx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1, 0.1, 0.1);
   vx_fit->SetParNames("vx const", "vx mean", "vx sigma");
   vx_hist->Fit(vx_fit, "R");
   vx_hist->Draw();
   vx_hist->Write();
      vy_hist->Add(ivy_hist);
      vy_hist->Add(iivy_hist);
      vy_hist->Add(iiivy_hist);
      vy_hist->Add(ivvy_hist);
      vy_hist->Add(vvy_hist);
      vy_hist->Add(vivy_hist);
      vy_hist->Scale(1.0 / vy_hist->Integral());
    //  vy_hist->SetMinimum(0.0);
    //  vy_hist->SetMaximum(0.12);
      TF1 *vy_fit = new TF1("vy_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   vy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   vy_fit->SetParNames("vy const", "vy mean", "vy sigma");
   vy_hist->Fit(vy_fit, "R");
   vy_hist->Draw();
   vy_hist->Write();
      vz_hist->Add(ivz_hist);
      vz_hist->Add(iivz_hist);
      vz_hist->Add(iiivz_hist);
      vz_hist->Add(ivvz_hist);
      vz_hist->Add(vvz_hist);
      vz_hist->Add(vivz_hist);
      vz_hist->Scale(1.0 / vz_hist->Integral());
   //   vz_hist->SetMinimum(0.0);
   //   vz_hist->SetMaximum(0.04);
      //TF1 *vz_fit = new TF1("vz_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -12.0);
        //  vz_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
        //  vz_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
           TF1 *vz_fit = new TF1("vz_fit", "[0] -[1]/([2]*x + [3]) + gaus(4) + gaus(7)", -25.0, -10.5);
          vz_fit->SetParameters(0, 0.1, 20, 142.687, 0.01, -19, 0.5, 0.01, -14, 0.5);
          vz_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *vz2_fit = new TF1("vz2_fit", "gaus(0)+gaus(3)+pol1(6)", -10.5, -2);
          vz2_fit->SetParameters(0.01, -8, 1, 0.01, -2, 1, 0.01, 0);
          vz2_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
      //   TF1 *vz2_fit = new TF1("vz2_fit", "gaus(0)+pol1(3)", -11, 0);
      //    vz2_fit->SetParameters(0.01, -7, 1);
      //    vz2_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma"); 
   vz_hist->Fit(vz_fit, "R");
   vz_hist->Fit(vz2_fit, "R+");
   vz_hist->Draw("L");
   vz_hist->Write();
   
   
      vxx_hist->Scale(1.0 / vxx_hist->Integral());
  //    vxx_hist->SetMinimum(0.0);
  //    vxx_hist->SetMaximum(0.4);
      TF1 *vxx_fit = new TF1("vxx_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   vxx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   vxx_fit->SetParNames("vxx const", "vxx mean", "vxx sigma");
   vxx_hist->Fit(vxx_fit, "R");
   vxx_hist->Draw();
   vxx_hist->Write();
      
      vyy_hist->Scale(1.0 / vyy_hist->Integral());
  //    vyy_hist->SetMinimum(0.0);
  //    vyy_hist->SetMaximum(0.12);
      TF1 *vyy_fit = new TF1("vyy_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   vyy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   vyy_fit->SetParNames("vyy const", "vyy mean", "vyy sigma");
   vyy_hist->Fit(vyy_fit, "R");
   vyy_hist->Draw();
   vyy_hist->Write();
   
   vzz_hist->Scale(1.0 / vzz_hist->Integral());
      TF1 *vzz_fit = new TF1("vzz_fit", "gaus(0)+gaus(3)+pol1(6)", -23.0, -11.0);
          vzz_fit->SetParameters(0.01, -18, 0.5, 0.01, -14, 0.5);
          vzz_fit->SetParNames("enter const", "enter mean", "enter sigma", "exit const", "exit mean", "exit sigma");
       TF1 *vzz2_fit = new TF1("vzz2_fit", "gaus(0)+gaus(3)+pol1(6)", -10.5, -2);
          vzz2_fit->SetParameters(0.01, -8, 1, 0.01, -2, 1, 0.01, 0);
          vzz2_fit->SetParNames("Copper const", "Copper mean", "Copper sigma", "Tin const", "Tin mean", "Tin sigma");
   //vzz_hist->Fit(vzz_fit, "R");
   //vzz_hist->Fit(vzz2_fit, "R+");
   vzz_hist->SetLineWidth(4);
   vzz_hist->Draw("L");
   vzz_hist->Write();
      
      //std::cout << "this is where vz is -------------------------------------" << std::endl;
      pip_vx_hist->Add(pip_ivx_hist);
      pip_vx_hist->Add(pip_iivx_hist);
      pip_vx_hist->Add(pip_iiivx_hist);
      pip_vx_hist->Add(pip_ivvx_hist);
      pip_vx_hist->Add(pip_vvx_hist);
      pip_vx_hist->Add(pip_vivx_hist);
      pip_vx_hist->Scale(1.0 / pip_vx_hist->Integral());
      //pip_vx_hist->SetMinimum(0.0);
      //pip_vx_hist->SetMaximum(0.14);
      TF1 *pip_vx_fit = new TF1("pip_vx_fit", "gaus(0)+pol4(3)",  - 1.5, 1.5);
   pip_vx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pip_vx_hist->Fit(pip_vx_fit, "R");
   pip_vx_hist->Draw();
   pip_vx_hist->Write();
      pip_vy_hist->Add(pip_ivy_hist);
      pip_vy_hist->Add(pip_iivy_hist);
      pip_vy_hist->Add(pip_iiivy_hist);
      pip_vy_hist->Add(pip_ivvy_hist);
      pip_vy_hist->Add(pip_vvy_hist);
      pip_vy_hist->Add(pip_vivy_hist);
      pip_vy_hist->Scale(1.0 / pip_vy_hist->Integral());
      //pip_vy_hist->SetMinimum(0.0);
      //pip_vy_hist->SetMaximum(0.1);
      TF1 *pip_vy_fit = new TF1("pip_vy_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   pip_vy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pip_vy_hist->Fit(pip_vy_fit, "R");
   pip_vy_hist->Draw();
   pip_vy_hist->Write();
      pip_vz_hist->Add(pip_ivz_hist);
      pip_vz_hist->Add(pip_iivz_hist);
      pip_vz_hist->Add(pip_iiivz_hist);
      pip_vz_hist->Add(pip_ivvz_hist);
      pip_vz_hist->Add(pip_vvz_hist);
      pip_vz_hist->Add(pip_vivz_hist);
      pip_vz_hist->Scale(1.0 / pip_vz_hist->Integral());
      //pip_vz_hist->SetMinimum(0.0);
      //pip_vz_hist->SetMaximum(0.03);
      TF1 *pip_vz_fit = new TF1("pip_vz_fit", "gaus(0) + gaus(3) + pol1(6)",  - 11, -1);
   pip_vz_fit->SetParameters(0.01, -7, 1, 0.01, -2, 1);
   pip_vz_hist->Fit(pip_vz_fit, "R");
   pip_vz_hist->Draw();
   pip_vz_hist->Write();
      
      
      pim_vx_hist->Add(pim_ivx_hist);
      pim_vx_hist->Add(pim_iivx_hist);
      pim_vx_hist->Add(pim_iiivx_hist);
      pim_vx_hist->Add(pim_ivvx_hist);
      pim_vx_hist->Add(pim_vvx_hist);
      pim_vx_hist->Add(pim_vivx_hist);
      pim_vx_hist->Scale(1.0 / pim_vx_hist->Integral());
      //pim_vx_hist->SetMinimum(0.0);
      //pim_vx_hist->SetMaximum(0.14);
      TF1 *pim_vx_fit = new TF1("pim_vx_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   pim_vx_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pim_vx_hist->Fit(pim_vx_fit, "R");
      pim_vx_hist->Draw();
   pim_vx_hist->Write();
      pim_vy_hist->Add(pim_ivy_hist);
      pim_vy_hist->Add(pim_iivy_hist);
      pim_vy_hist->Add(pim_iiivy_hist);
      pim_vy_hist->Add(pim_ivvy_hist);
      pim_vy_hist->Add(pim_vvy_hist);
      pim_vy_hist->Add(pim_vivy_hist);
      pim_vy_hist->Scale(1.0 / pim_vy_hist->Integral());
      //pim_vy_hist->SetMinimum(0.0);
      //pim_vy_hist->SetMaximum(0.1);
      TF1 *pim_vy_fit = new TF1("pim_vy_fit", "gaus(0)+pol4(3)",  - 1.5,  1.5);
   pim_vy_fit->SetParameters(0.1, 0, 0.1, 0.001, 0.1, 0.1);
   pim_vy_hist->Fit(pim_vy_fit, "R");
      pim_vy_hist->Draw();
   pim_vy_hist->Write();
      pim_vz_hist->Add(pim_ivz_hist);
      pim_vz_hist->Add(pim_iivz_hist);
      pim_vz_hist->Add(pim_iiivz_hist);
      pim_vz_hist->Add(pim_ivvz_hist);
      pim_vz_hist->Add(pim_vvz_hist);
      pim_vz_hist->Add(pim_vivz_hist);
      pim_vz_hist->Scale(1.0 / pim_vz_hist->Integral());
      //pim_vz_hist->SetMinimum(0.0);
      //pim_vz_hist->SetMaximum(0.03);
      TF1 *pim_vz_fit = new TF1("pim_vz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, 0);
      pim_vz_fit->SetParameters(0.01, -8, 1, 0.01, -3, 1);
      pim_vz_hist->Fit(pim_vz_fit, "R");
      pim_vz_hist->Draw();
      pim_vz_hist->Write();
      
      
   TCanvas *vxt = new TCanvas("vxt","vx",0,0,640,480);
     // TLine *llxlim5=new TLine(lxlim5,vxt->GetUymin(),lxlim5,vxt->GetUymax());
    //  llxlim5->SetLineColor(kCyan-1);
    //  TLine *luxlim5=new TLine(uxlim5,vxt->GetUymin(),uxlim5,vxt->GetUymax());
    //  luxlim5->SetLineColor(kCyan-1);
   ivx_hist->Draw("HIST");
   vx_hist->Draw("HIST Same");
   iivx_hist->Draw("HIST Same");
   iiivx_hist->Draw("HIST Same");
   ivvx_hist->Draw("HIST Same");
   vvx_hist->Draw("HIST Same");
   vivx_hist->Draw("HIST Same");
  // llxlim5->Draw("HIST Same");
  // luxlim5->Draw("HIST Same");
   vxt->Write();
   
   
   TCanvas *vyt = new TCanvas("vyt","vy",0,0,640,480);
    //  TLine *llylim5=new TLine(lylim5,vyt->GetUymin(),lylim5,vyt->GetUymax());
    //  llylim5->SetLineColor(kCyan-1);
    //  TLine *luylim5=new TLine(uylim5,vyt->GetUymin(),uylim5,vyt->GetUymax());
    //  luylim5->SetLineColor(kCyan-1);
   vivy_hist->Draw("HIST");
   ivy_hist->Draw("HIST Same");
   iivy_hist->Draw("HIST Same");
   iiivy_hist->Draw("HIST Same");
   ivvy_hist->Draw("HIST Same");
   vy_hist->Draw("HIST Same");
   vvy_hist->Draw("HIST Same");
  // llylim5->Draw("HIST Same");
  // luylim5->Draw("HIST Same");
   vyt->Write();
   
   /*
   
   //fmt_vz_hist->Scale(1.0 / fmt_vz_hist->Integral());
   fmt_vz1_hist->Scale(1.0 / fmt_vz1_hist->Integral());
   fmt_vz2_hist->Scale(1.0 / fmt_vz2_hist->Integral());
   fmt_vz3_hist->Scale(1.0 / fmt_vz3_hist->Integral());
   fmt_vz4_hist->Scale(1.0 / fmt_vz4_hist->Integral());
   fmt_vz5_hist->Scale(1.0 / fmt_vz5_hist->Integral());
   fmt_vz6_hist->Scale(1.0 / fmt_vz6_hist->Integral());
   
   TCanvas *vzt = new TCanvas("vzt","vz",0,0,640,480);
   //   TLine *llclim5=new TLine(lvy1c,vzt->GetUymin(),lvy1c,vzt->GetUymax());
   //   llclim5->SetLineColor(kCyan-1);
   //   llclim5->SetLineWidth(3);
    //  TLine *luclim5=new TLine(uvy1c,vzt->GetUymin(),uvy1c,vzt->GetUymax());
    //  luclim5->SetLineColor(kCyan-1);
    //  luclim5->SetLineWidth(3);
    //  TLine *llslim5=new TLine(ls5,vzt->GetUymin(),ls5,vzt->GetUymax());
    //  llslim5->SetLineColor(kMagenta-7);
    //  llslim5->SetLineWidth(3);
    //  llslim5->SetLineStyle(10);
    //  TLine *luslim5=new TLine(us5,vzt->GetUymin(),us5,vzt->GetUymax());
    //  luslim5->SetLineColor(kMagenta-7);
    //  luslim5->SetLineWidth(3);
    //  luslim5->SetLineStyle(10);
    fmt_vz_hist->SetLineWidth(3);
    fmt_vz1_hist->SetLineWidth(3);
    fmt_vz2_hist->SetLineWidth(3);
    fmt_vz3_hist->SetLineWidth(3);
    fmt_vz4_hist->SetLineWidth(3);
    fmt_vz5_hist->SetLineWidth(3);
    fmt_vz6_hist->SetLineWidth(3);
    vz_hist->SetLineWidth(3);
    ivz_hist->SetLineWidth(3);
    iivz_hist->SetLineWidth(3);
    iiivz_hist->SetLineWidth(3);
    ivvz_hist->SetLineWidth(3);
    vvz_hist->SetLineWidth(3);
    vivz_hist->SetLineWidth(3);
   //vz_hist->Draw("HIST");
   //ivz_hist->Draw("HIST Same");
   //iivz_hist->Draw("HIST Same");
   //iiivz_hist->Draw("HIST Same");
   //ivvz_hist->Draw("HIST Same");
   //vvz_hist->Draw("HIST Same");
   //vivz_hist->Draw("HIST Same");
   fmt_vz_hist->Draw("HIST");
   fmt_vz1_hist->Draw("HIST Same");
   fmt_vz2_hist->Draw("HIST Same");
   fmt_vz3_hist->Draw("HIST Same");
   fmt_vz4_hist->Draw("HIST Same");
   fmt_vz5_hist->Draw("HIST Same");
   fmt_vz6_hist->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vzt->BuildLegend();
   vzt->Write();
   
   TCanvas *vz1t = new TCanvas("vz1t","vz",0,0,640,480);
    fmt_vz1_hist->SetLineWidth(3);
    ivz_hist->SetLineWidth(3);
    ivz_hist->Draw("HIST Same");
   fmt_vz1_hist->Draw("HIST Same");
   vz1t->BuildLegend();
   vz1t->Write();
   
   TCanvas *vz2t = new TCanvas("vz2t","vz",0,0,640,480);
    fmt_vz2_hist->SetLineWidth(3);
    iivz_hist->SetLineWidth(3);
   iivz_hist->Draw("HIST Same");
   fmt_vz2_hist->Draw("HIST Same");
   vz2t->BuildLegend();
   vz2t->Write();
   
   TCanvas *vz3t = new TCanvas("vz3t","vz",0,0,640,480);
    fmt_vz3_hist->SetLineWidth(3);
    iiivz_hist->SetLineWidth(3);
   iiivz_hist->Draw("HIST Same");
   fmt_vz3_hist->Draw("HIST Same");
   vz3t->BuildLegend();
   vz3t->Write();
   
   TCanvas *vz4t = new TCanvas("vz4t","vz",0,0,640,480);
    fmt_vz4_hist->SetLineWidth(3);
    ivvz_hist->SetLineWidth(3);
   ivvz_hist->Draw("HIST Same");
   fmt_vz4_hist->Draw("HIST Same");
   vz4t->BuildLegend();
   vz4t->Write();
   
   TCanvas *vz5t = new TCanvas("vz5t","vz",0,0,640,480);
    fmt_vz5_hist->SetLineWidth(3);
    vvz_hist->SetLineWidth(3);
   vvz_hist->Draw("HIST Same");
   fmt_vz5_hist->Draw("HIST Same");
   vz5t->BuildLegend();
   vz5t->Write();
   
   TCanvas *vz6t = new TCanvas("vz6t","vz",0,0,640,480);
    fmt_vz6_hist->SetLineWidth(3);
    vivz_hist->SetLineWidth(3);
   vivz_hist->Draw("HIST Same");
   fmt_vz6_hist->Draw("HIST Same");
   vz6t->BuildLegend();
   vz6t->Write();
   */
   
   te_g->Draw();
   te_g->Write();
   
   tpip_g->Draw();
   tpip_g->Write();
   
   tpim_g->Draw();
   tpim_g->Write();
   
   TCanvas *pip_vxt = new TCanvas("pip_vxt","pip_vx",0,0,640,480);
   //   TLine *pip_llxlim5=new TLine(pip_lxlim5,pip_vxt->GetUymin(),pip_lxlim5,pip_vxt->GetUymax());
   //   pip_llxlim5->SetLineColor(kCyan-1);
   //   TLine *pip_luxlim5=new TLine(pip_uxlim5,pip_vxt->GetUymin(),pip_uxlim5,pip_vxt->GetUymax());
   //   pip_luxlim5->SetLineColor(kCyan-1);
   pip_ivx_hist->Draw("HIST");
   pip_vx_hist->Draw("HIST Same");
   pip_iivx_hist->Draw("HIST Same");
   pip_iiivx_hist->Draw("HIST Same");
   pip_ivvx_hist->Draw("HIST Same");
   pip_vvx_hist->Draw("HIST Same");
   pip_vivx_hist->Draw("HIST Same");
   //pip_llxlim5->Draw("HIST Same");
   //pip_luxlim5->Draw("HIST Same");
   pip_vxt->Write();
   
   
   TCanvas *pip_vyt = new TCanvas("pip_vyt","pip_vy",0,0,640,480);
   //   TLine *pip_llylim5=new TLine(pip_lylim5,pip_vyt->GetUymin(),pip_lylim5,pip_vyt->GetUymax());
   //   pip_llylim5->SetLineColor(kCyan-1);
   //   TLine *pip_luylim5=new TLine(pip_uylim5,pip_vyt->GetUymin(),pip_uylim5,pip_vyt->GetUymax());
   //   pip_luylim5->SetLineColor(kCyan-1);
   pip_vivy_hist->Draw("HIST");
   pip_ivy_hist->Draw("HIST Same");
   pip_iivy_hist->Draw("HIST Same");
   pip_iiivy_hist->Draw("HIST Same");
   pip_ivvy_hist->Draw("HIST Same");
   pip_vy_hist->Draw("HIST Same");
   pip_vvy_hist->Draw("HIST Same");
   //pip_llylim5->Draw("HIST Same");
   //pip_luylim5->Draw("HIST Same");
   pip_vyt->Write();
   
   
   TCanvas *pip_vzt = new TCanvas("pip_vzt","pip_vz",0,0,640,480);
      /*
      TLine *llclim5=new TLine(lvy1c,vz1c->GetUymin(),lvy1c,ivz_hist->GetMaximum());
      llclim5->SetLineColor(kCyan-1);
      llclim5->SetLineWidth(3);
      TLine *luclim5=new TLine(uvy1c,vz1c->GetUymin(),uvy1c,ivz_hist->GetMaximum());
      luclim5->SetLineColor(kCyan-1);
      luclim5->SetLineWidth(3);
      TLine *llslim5=new TLine(ls5,vz1c->GetUymin(),ls5,ivz_hist->GetMaximum());
      llslim5->SetLineColor(kMagenta-7);
      llslim5->SetLineWidth(3);
      llslim5->SetLineStyle(10);
      TLine *luslim5=new TLine(us5,vz1c->GetUymin(),us5,ivz_hist->GetMaximum());
      luslim5->SetLineColor(kMagenta-7);
      luslim5->SetLineWidth(3);
      luslim5->SetLineStyle(10);
      */
   pip_vz_hist->Draw("HIST");
   pip_ivz_hist->Draw("HIST Same");
   pip_iivz_hist->Draw("HIST Same");
   pip_iiivz_hist->Draw("HIST Same");
   pip_ivvz_hist->Draw("HIST Same");
   pip_vvz_hist->Draw("HIST Same");
   pip_vivz_hist->Draw("HIST Same");
   //pip_llzlim5->Draw("HIST Same");
   //pip_luzlim5->Draw("HIST Same");
   pip_vzt->Write();
   
   
   
   TCanvas *pim_vxt = new TCanvas("pim_vxt","pim_vx",0,0,640,480);
   //   TLine *pim_llxlim5=new TLine(pim_lxlim5,pim_vxt->GetUymin(),pim_lxlim5,pim_vxt->GetUymax());
   //   pim_llxlim5->SetLineColor(kCyan-1);
   //   TLine *pim_luxlim5=new TLine(pim_uxlim5,pim_vxt->GetUymin(),pim_uxlim5,pim_vxt->GetUymax());
   //   pim_luxlim5->SetLineColor(kCyan-1);
   pim_ivx_hist->Draw("HIST");
   pim_vx_hist->Draw("HIST Same");
   pim_iivx_hist->Draw("HIST Same");
   pim_iiivx_hist->Draw("HIST Same");
   pim_ivvx_hist->Draw("HIST Same");
   pim_vvx_hist->Draw("HIST Same");
   pim_vivx_hist->Draw("HIST Same");
   //pim_llxlim5->Draw("HIST Same");
   //pim_luxlim5->Draw("HIST Same");
   pim_vxt->Write();
   
   
   TCanvas *pim_vyt = new TCanvas("pim_vyt","pim_vy",0,0,640,480);
   //   TLine *pim_llylim5=new TLine(pim_lylim5,pim_vyt->GetUymin(),pim_lylim5,pim_vyt->GetUymax());
   //   pim_llylim5->SetLineColor(kCyan-1);
   //   TLine *pim_luylim5=new TLine(pim_uylim5,pim_vyt->GetUymin(),pim_uylim5,pim_vyt->GetUymax());
   //   pim_luylim5->SetLineColor(kCyan-1);
   pim_vivy_hist->Draw("HIST");
   pim_ivy_hist->Draw("HIST Same");
   pim_iivy_hist->Draw("HIST Same");
   pim_iiivy_hist->Draw("HIST Same");
   pim_ivvy_hist->Draw("HIST Same");
   pim_vy_hist->Draw("HIST Same");
   pim_vvy_hist->Draw("HIST Same");
   //pim_llylim5->Draw("HIST Same");
   //pim_luylim5->Draw("HIST Same");
   pim_vyt->Write();
   
   
   TCanvas *pim_vzt = new TCanvas("pim_vzt","pim_vz",0,0,640,480);
   /*
      TLine *pim_llzlim5=new TLine(pim_lzlim5,pim_vzt->GetUymin(),pim_lzlim5,pim_vzt->GetUymax());
      pim_llzlim5->SetLineColor(kCyan-1);
      TLine *pim_luzlim5=new TLine(pim_uzlim5,pim_vzt->GetUymin(),pim_uzlim5,pim_vzt->GetUymax());
      pim_luzlim5->SetLineColor(kCyan-1);
      */
   pim_vz_hist->Draw("HIST");
   pim_ivz_hist->Draw("HIST Same");
   pim_iivz_hist->Draw("HIST Same");
   pim_iiivz_hist->Draw("HIST Same");
   pim_ivvz_hist->Draw("HIST Same");
   pim_vvz_hist->Draw("HIST Same");
   pim_vivz_hist->Draw("HIST Same");
   //pim_llzlim5->Draw("HIST Same");
   //pim_luzlim5->Draw("HIST Same");
   pim_vzt->Write();
   
   
   TCanvas *sampc = new TCanvas("sampc","sampc",0,0,640,480);
   auto mgsampc = new TMultiGraph;
   mgsampc->Add(samp_g);
   mgsampc->Add(sampu_g);
   mgsampc->Add(sampd_g);
   mgsampc->SetMaximum(0.4);
   mgsampc->SetMinimum(0.0);
   samp_hist->Draw("colz");
   mgsampc->Draw("P");
   sampc->Modified();
   sampc->Update();
   sampc->Write();
   mgsampc->Write();
   
   
   TCanvas *sampc1 = new TCanvas("sampc1","sampc1",0,0,640,480);
   auto mgsampc1 = new TMultiGraph;
   mgsampc1->Add(samp1_g);
   mgsampc1->Add(sampu1_g);
   mgsampc1->Add(sampd1_g);
   mgsampc1->SetMaximum(0.4);
   mgsampc1->SetMinimum(0.0);
   samp1_hist->Draw("colz");
   mgsampc1->Draw("P");
   sampc1->Modified();
   sampc1->Update();
   sampc1->Write();
   mgsampc1->Write();
   
   
   TCanvas *sampc1cut = new TCanvas("sampc1cut","sampc1cut",0,0,640,480);
   auto mgsampc1cut = new TMultiGraph;
   mgsampc1cut->Add(samp1_g);
   mgsampc1cut->Add(sampu1_g);
   mgsampc1cut->Add(sampd1_g);
   mgsampc1cut->SetMaximum(0.4);
   mgsampc1cut->SetMinimum(0.0);
   samp1_cut_hist->Draw("colz");
   mgsampc1cut->Draw("P");
   sampc1cut->Modified();
   sampc1cut->Update();
   sampc1cut->Write();
   mgsampc1cut->Write();
   
   
   TCanvas *sampc2 = new TCanvas("sampc2","sampc2",0,0,640,480);
   auto mgsampc2 = new TMultiGraph;
   mgsampc2->Add(samp2_g);
   mgsampc2->Add(sampu2_g);
   mgsampc2->Add(sampd2_g);
   mgsampc2->SetMaximum(0.4);
   mgsampc2->SetMinimum(0.0);
   samp2_hist->Draw("colz");
   mgsampc2->Draw("P");
   sampc2->Modified();
   sampc2->Update();
   sampc2->Write();
   mgsampc2->Write();
   
   
   TCanvas *sampc3 = new TCanvas("sampc3","sampc3",0,0,640,480);
   auto mgsampc3 = new TMultiGraph;
   mgsampc3->Add(samp3_g);
   mgsampc3->Add(sampu3_g);
   mgsampc3->Add(sampd3_g);
   mgsampc3->SetMaximum(0.4);
   mgsampc3->SetMinimum(0.0);
   samp3_hist->Draw("colz");
   mgsampc3->Draw("P");
   sampc3->Modified();
   sampc3->Update();
   sampc3->Write();
   mgsampc3->Write();
   
   
   TCanvas *sampc4 = new TCanvas("sampc4","sampc4",0,0,640,480);
   auto mgsampc4 = new TMultiGraph;
   mgsampc4->Add(samp4_g);
   mgsampc4->Add(sampu4_g);
   mgsampc4->Add(sampd4_g);
   mgsampc4->SetMaximum(0.4);
   mgsampc4->SetMinimum(0.0);
   samp4_hist->Draw("colz");
   mgsampc4->Draw("P");
   sampc4->Modified();
   sampc4->Update();
   sampc4->Write();
   mgsampc4->Write();
   
   
   
   TCanvas *sampc5 = new TCanvas("sampc5","sampc5",0,0,640,480);
   auto mgsampc5 = new TMultiGraph;
   mgsampc5->Add(samp5_g);
   mgsampc5->Add(sampu5_g);
   mgsampc5->Add(sampd5_g);
   mgsampc5->SetMaximum(0.4);
   mgsampc5->SetMinimum(0.0);
   samp5_hist->Draw("colz");
   mgsampc5->Draw("P");
   sampc5->Modified();
   sampc5->Update();
   sampc5->Write();
   mgsampc5->Write();
   
   
   
   TCanvas *sampc6 = new TCanvas("sampc6","sampc6",0,0,640,480);
   auto mgsampc6 = new TMultiGraph;
   mgsampc6->Add(samp6_g);
   mgsampc6->Add(sampu6_g);
   mgsampc6->Add(sampd6_g);
   mgsampc6->SetMaximum(0.4);
   mgsampc6->SetMinimum(0.0);
   samp6_hist->Draw("colz");
   mgsampc6->Draw("P");
   sampc6->Modified();
   sampc6->Update();
   sampc6->Write();
   mgsampc6->Write();
   
   TCanvas *samvc = new TCanvas("samvc","samvc",0,0,640,480);
   auto mgsamvc = new TMultiGraph;
   mgsamvc->Add(samv_g);
   mgsamvc->Add(samvu_g);
   mgsamvc->Add(samvd_g);
   mgsamvc->SetMaximum(0.4);
   mgsamvc->SetMinimum(0.0);
   samv_hist->Draw("colz");
   mgsamvc->Draw("P");
   samvc->Modified();
   samvc->Update();
   samvc->Write();
   mgsamvc->Write();
   
   TCanvas *samwc = new TCanvas("samwc","samwc",0,0,640,480);
   auto mgsamwc = new TMultiGraph;
   mgsamwc->Add(samw_g);
   mgsamwc->Add(samwu_g);
   mgsamwc->Add(samwd_g);
   mgsamwc->SetMaximum(0.4);
   mgsamwc->SetMinimum(0.0);
   samw_hist->Draw("colz");
   mgsamwc->Draw("P");
   samwc->Modified();
   samwc->Update();
   samwc->Write();
   mgsamwc->Write();
   
   thve_g->Draw();
   thve_g->Write();
   thveu_g->Draw();
   thveu_g->Write();
   thved_g->Draw();
   thved_g->Write();
   
   auto mgthvec = new TMultiGraph;
   //mgthvec->Add(thve_g);
   mgthvec->Add(thveu_g);
   //mgthvec->Add(thved_g);
   //mgthvec->SetMaximum(0.4);
   //mgthvec->SetMinimum(0.0);
   
   
   TCanvas *thvec = new TCanvas("thvec","thvec",0,0,640,480);
   thve_hist->Draw("colz");
   mgthvec->Draw("P");
   thvec->Modified();
   thvec->Update();
   thvec->Write();
   mgthvec->Write();
   
   TCanvas *chin1cav = new TCanvas("chin1cav","chin1cav",0,0,640,480);
   chin1cav->Range(0.0, 0.0, 30.0, 11.0);
   chin_hist1->Draw();
   chin1_x->Draw();
   chin1_x2->Draw();
   chin1_y->Draw();
   chin1cav->Write();
   
   TCanvas *chin2cav = new TCanvas("chin2cav","chin2cav",0,0,640,480);
   chin2cav->Range(0.0, 0.0, 30.0, 11.0);
   chin_hist2->Draw();
   chin2_x->Draw();
   chin2_x2->Draw();
   chin2_y->Draw();
   chin2cav->Write();
   
   TCanvas *chin3cav = new TCanvas("chin3cav","chin3cav",0,0,640,480);
   chin3cav->Range(0.0, 0.0, 30.0, 11.0);
   chin_hist3->Draw();
   chin3_x->Draw();
   chin3_x2->Draw();
   chin3_y->Draw();
   chin3cav->Write();
   
   TCanvas *chin4cav = new TCanvas("chin4cav","chin4cav",0,0,640,480);
   chin4cav->Range(0.0, 0.0, 30.0, 11.0);
   chin_hist4->Draw();
   chin4_x->Draw();
   chin4_x2->Draw();
   chin4_y->Draw();
   chin4cav->Write();
   
   TCanvas *chin5cav = new TCanvas("chin5cav","chin5cav",0,0,640,480);
   chin5cav->Range(0.0, 0.0, 30.0, 11.0);
   chin_hist5->Draw();
   chin5_x->Draw();
   chin5_x2->Draw();
   chin5_y->Draw();
   chin5cav->Write();
   
   TCanvas *chin6cav = new TCanvas("chin6cav","chin6cav",0,0,640,480);
   chin6cav->Range(0.0, 0.0, 30.0, 11.0);
   chin_hist6->Draw();
   chin6_x->Draw();
   chin6_x2->Draw();
   chin6_y->Draw();
   chin6cav->Write();
   
   thveu_g1->Draw();
   thveu_g1->Write();
   
   thveu_g2->Draw();
   thveu_g2->Write();
   
   thveu_g3->Draw();
   thveu_g3->Write();
   
   thveu_g4->Draw();
   thveu_g4->Write();
   
   thveu_g5->Draw();
   thveu_g5->Write();
   
   thveu_g6->Draw();
   thveu_g6->Write();
   
   TCanvas *means_sect = new TCanvas("means_sect", "mean_sect", 0, 0, 640, 480);
   thveu_g->Draw();
   thveu_g1->Draw("SAME");
   thveu_g2->Draw("SAME");
   thveu_g3->Draw("SAME");
   thveu_g4->Draw("SAME");
   thveu_g5->Draw("SAME");
   thveu_g6->Draw("SAME");
   means_sect->BuildLegend();
   means_sect->Draw();
   means_sect->Write();
   
   
   
   
  auto mgthvec22 = new TMultiGraph;
   mgthvec22->Add(thveu_g);
   mgthvec22->Add(thveu_g1);
   mgthvec22->Add(thveu_g2);
   mgthvec22->Add(thveu_g3);
   mgthvec22->Add(thveu_g4);
   mgthvec22->Add(thveu_g5);
   mgthvec22->Add(thveu_g6);
   mgthvec22->Draw("AP");
   mgthvec22->Write();




 
   ietot_hist->GetXaxis()->SetTitle("E_{tot}/P");
   TF1 *ietot_func = new TF1("ietot_fit",fitf,-1,1,3);
   ietot_func->SetParameters(500,ietot_hist->GetMean(),(ietot_hist->GetRMS()));
   ietot_func->SetParNames("Constant","Mean_Value","Sigma");
   //ietot_hist->Fit("ietot_fit","","",0.21,0.26);
   //ietot_hist->Draw();
   //ietot_hist->Write();

   p_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   TF1 *p_func = new TF1("p_fit",fitfbw,-1,1,7);
   p_func->SetParameters(0,0,0, 1,-142857, 1, 70);  //500,p_hist->GetMean(),p_hist->GetRMS(), -142857, 0, 70);
   //p_func->SetParLimits(0,500, 1500);
   //p_func->SetParLimits(1,-0.2,0.2);
   //p_func->SetParLimits(2, );
   p_func->SetParLimits(3,-20000, 20000);
   p_func->SetParLimits(4,-2000000, 0);
   p_func->SetParLimits(5,0, 2000000);
   p_func->SetParLimits(6,-20000, 20000);
   p_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //p_hist->Fit("p_fit","","",0.3,1.7);
   p_hist->Draw();
   p_hist->Write();

   ip_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   TF1 *ip_func = new TF1("ip_fit",fitfbw,-1,1,7);
   ip_func->SetParameters(500,0.77,(150.0/2), 1,-142857, 1, 70);
   ip_func->SetParLimits(0,-20000, 20000);
   //ip_func->SetParLimits(1,-20000, 20000);
   //ip_func->SetParLimits(2,-20000, 20000);
   ip_func->SetParLimits(3,-20000, 20000);
   ip_func->SetParLimits(4,-2000000,0);
   ip_func->SetParLimits(5,0,2000000);
   ip_func->SetParLimits(6,-20000,20000);
   ip_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //ip_hist->Fit("ip_fit","","",0.3,1.7);
   ip_hist->Draw();
   ip_hist->Write();



   iietot_hist->GetXaxis()->SetTitle("E_{tot}/P");
   TF1 *iietot_func = new TF1("iietot_fit",fitf,-1,1,3);
   iietot_func->SetParameters(500,iietot_hist->GetMean(),iietot_hist->GetRMS());
   iietot_func->SetParNames("Constant","Mean_Value","Sigma");
   //iietot_hist->Fit("iietot_fit","","",0.21,0.26);
   //iietot_hist->Draw();
   //iietot_hist->Write();

   iip_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   TF1 *iip_func = new TF1("iip_fit",fitfbw,-1,1,7);
   iip_func->SetParameters(500,0.77,(150.0/2), 0,0, 0, 0);
   iip_func->SetParLimits(0,-20000, 20000);
   //iip_func->SetParLimits(1,-20000, 20000);
   //iip_func->SetParLimits(2,-20000, 20000);
   //iip_func->SetParLimits(3,-20000, 20000);
   //iip_func->SetParLimits(4,-200000,20000);
   //iip_func->SetParLimits(5,-20000,20000);
   //iip_func->SetParLimits(6,-20000,20000);
   iip_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //iip_hist->Fit("iip_fit","","",0.6,1);
   iip_hist->Draw();
   iip_hist->Write();




   iiietot_hist->GetXaxis()->SetTitle("E_{tot}/P");
   TF1 *iiietot_func = new TF1("iiietot_fit",fitf,-1,1,3);
   iiietot_func->SetParameters(500,iiietot_hist->GetMean(),iiietot_hist->GetRMS());
   iiietot_func->SetParNames("Constant","Mean_Value","Sigma");
   //iiietot_hist->Fit("iiietot_fit","","",0.21,0.26);
   //iiietot_hist->Draw();
   //iiietot_hist->Write();

   iiip_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   TF1 *iiip_func = new TF1("iiip_fit",fitfbw,-1,1,7);
   iiip_func->SetParameters(500,0.77,(150.0/2), 0,-18000, 0, 45);     //iiip_hist->GetMean(),iiip_hist->GetRMS(
   iiip_func->SetParLimits(0,-20000, 20000);
   //iiip_func->SetParLimits(1,-20000, 20000);
   //iiip_func->SetParLimits(2,-20000, 20000);
   //iiip_func->SetParLimits(3,-20000, 20000);
   iiip_func->SetParLimits(4,-2000000,0);
   iiip_func->SetParLimits(5,0,2000000);
   iiip_func->SetParLimits(6,-20000,20000);
   iiip_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
  // iiip_hist->Fit("iiip_fit","","",0.3,1.7);
   iiip_hist->Draw();
   iiip_hist->Write();




   ivetot_hist->GetXaxis()->SetTitle("E_{tot}/P");
   TF1 *ivetot_func = new TF1("ivetot_fit",fitf,-1,1,3);
   ivetot_func->SetParameters(500,ivetot_hist->GetMean(),ivetot_hist->GetRMS());
   ivetot_func->SetParNames("Constant","Mean_Value","Sigma");
   //ivetot_hist->Fit("ivetot_fit","","",0.21,0.26);
   //ivetot_hist->Draw();
   //ivetot_hist->Write();

   ivp_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   TF1 *ivp_func = new TF1("ivp_fit",fitfbw,-1,1,7);
   ivp_func->SetParameters(500,0.77,(150.0/2), 0,-18000, 0, 45);
   ivp_func->SetParLimits(0,-20000, 20000);
   //ivp_func->SetParLimits(1,-20000, 20000);
   //ivp_func->SetParLimits(2,-20000, 20000);
   ivp_func->SetParLimits(3,-20000, 20000);
   ivp_func->SetParLimits(4,-200000,20000);
   ivp_func->SetParLimits(5,-20000,20000);
   ivp_func->SetParLimits(6,-20000,20000);
   ivp_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //ivp_hist->Fit("ivp_fit","","",0.6,1.0);
   ivp_hist->Draw();
   ivp_hist->Write();




   vetot_hist->GetXaxis()->SetTitle("E_{tot}/P");
   TF1 *vetot_func = new TF1("vetot_fit",fitf,-1,1,3);
   vetot_func->SetParameters(500,vetot_hist->GetMean(),vetot_hist->GetRMS());
   vetot_func->SetParNames("Constant","Mean_Value","Sigma");
   //vetot_hist->Fit("vetot_fit","","",0.21,0.26);
   //vetot_hist->Draw();
   //vetot_hist->Write();

   vp_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   TF1 *vp_func = new TF1("vp_fit",fitfbw,-1,1,7);
   vp_func->SetParameters(500,0.77,(150.0/2), 0,-18000, 0, 45);
   vp_func->SetParLimits(0,-20000, 20000);
   //vp_func->SetParLimits(1,-20000, 20000);
   //vp_func->SetParLimits(2,-20000, 20000);
   vp_func->SetParLimits(3,-20000, 20000);
   vp_func->SetParLimits(4,-200000,20000);
   vp_func->SetParLimits(5,-20000,20000);
   vp_func->SetParLimits(6,-20000,20000);
   vp_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //vp_hist->Fit("vp_fit","","",0.6,1.0);
   vp_hist->Draw();
   vp_hist->Write();



   vietot_hist->GetXaxis()->SetTitle("E_{tot}/P");
   TF1 *vietot_func = new TF1("vietot_fit",fitf,-1,1,3);
   vietot_func->SetParameters(500,vietot_hist->GetMean(),vietot_hist->GetRMS());
   vietot_func->SetParNames("Constant","Mean_Value","Sigma");
   //vietot_hist->Fit("vietot_fit","","",0.21,0.26);
   //vietot_hist->Draw();
   //vietot_hist->Write();

   vip_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   TF1 *vip_func = new TF1("vip_fit",fitfbw,-1,1,7);
   vip_func->SetParameters(500,0.77,(150.0/2), 0,-18000, 0, 45);
   vip_func->SetParLimits(0,-20000, 20000);
   //vip_func->SetParLimits(1,-20000, 20000);
   //vip_func->SetParLimits(2,-20000, 20000);
   vip_func->SetParLimits(3,-20000, 20000);
   vip_func->SetParLimits(4,-200000,20000);
   vip_func->SetParLimits(5,-20000,20000);
   vip_func->SetParLimits(6,-20000,20000);
   vip_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //vip_hist->Fit("vip_fit","","",0.6,1.0);
   vip_hist->Draw();
   vip_hist->Write();

/*
   //etots->Add(ietot_hist);
   //etots->Add(iietot_hist);
   //etots->Add(iiietot_hist);
   etots->Add(iixy_hist);
   //etots->Add(vetot_hist);
   //etots->Add(vietot_hist);
   //etots->GetXaxis()->SetTitle("E_{tot} / P");
   etots->Draw("colz");
   c3->SetLogz();
   c3->Modified(); c3->Update();
   c3->SaveAs("zz_19125_etots_2.pdf");
   //etots->Write("nostack 0lego1 PFC");

   //ps->Add(ip_hist);
   //ps->Add(iip_hist);
   //ps->Add(iiip_hist);
   //ps->Add(ivp_hist);
   //ps->Add(vp_hist);
   //ps->Add(vip_hist);
   //ps->GetXaxis()->SetTitle("P [GeV]");
   //ps->Draw("nostack 0lego1 PFC");
   //ps->Write("nostack 0lego1 PFC");
/*
   ietot_p->SetTitle("Sector 1 P .vs E_{tot};P [GeV];E_{tot} [GeV]");
   ietot_p->Draw("AP");
   ietot_p->Write();

   iietot_p->SetTitle("Sector 2 P .vs E_{tot};P [GeV];E_{tot} [GeV]");
   iietot_p->Draw("AP");
   iietot_p->Write();

   iiietot_p->SetTitle("Sector 3 P .vs E_{tot};P [GeV];E_{tot} [GeV]");
   iiietot_p->Draw("AP");
   iiietot_p->Write();

   ivetot_p->SetTitle("Sector 4 P .vs E_{tot};P [GeV];E_{tot} [GeV]");
   ivetot_p->Draw("AP");
   ivetot_p->Write();

   vetot_p->SetTitle("Sector 5 P .vs E_{tot};P [GeV];E_{tot} [GeV]");
   vetot_p->Draw("AP");
   vetot_p->Write();

   vietot_p->SetTitle("Sector 6 P .vs E_{tot};P [GeV];E_{tot} [GeV]");
   vietot_p->Draw("AP");
   vietot_p->Write();

   TMultiGraph *both_hfcup = new TMultiGraph();
   both_hfcup->SetTitle("P .vs E_{tot}/P;P [GeV];E_{tot}/P");
   both_hfcup->Add(ietot_p);
   both_hfcup->Add(iietot_p);
   both_hfcup->Add(iiietot_p);
   both_hfcup->Add(ivetot_p);
   both_hfcup->Add(vetot_p);
   both_hfcup->Add(vietot_p);
   both_hfcup->Draw(); 
   //both_hfcup->Write();
*/
   

   //ebeam_ev->SetTitle("REC::event beamCharge;Timestamp [ns];REC::event beamCharge[nC]");
   //ebeam_ev->Draw("AP");
   //ebeam_ev->Write();

   //fcup_ev->SetTitle("Timestamp v RUN::scaler fcup;Timestamp;RUN::scaler fcup[nC]");
   //fcup_ev->Draw();
   //fcup_ev->Write();

   //fcupg_ev->SetTitle("RUN::scaler fcupgated[nC];Timestamp [ns];RUN::scaler fcupgated[nC]");
   //fcupg_ev->Draw("AP");
   //fcupg_ev->Write();

   //hel_fcup->SetTitle("HEL::scaler fcup[nC];Timestamp [ns];Hel::scaler fcup[nC]");
   //hel_fcup->Draw("AP");
   //hel_fcup->Write();

   //hel_fcupgated->SetTitle("HEL::scaler fcupgated[nC];Timestamp [ns];Hel::scaler fcupgated[nC]");
   //hel_fcupgated->Draw("AP");
   //hel_fcupgated->Write();

   //hel_clock->SetTitle("HEL::scaler clock;Timestamp [ns];Hel::scaler clock");
   //hel_clock->Draw("AP");
   //hel_clock->Write();

   //hel_clockgated->SetTitle("HEL::scaler clockgated;Timestamp [ns];Hel::scaler clockgated");
   //hel_clockgated->Draw("AP");
   //hel_clockgated->Write();

   //TMultiGraph *both = new TMultiGraph();
   //both->SetTitle("Scaler lifetime and event lifetime; Event number;DAQ Lifetime [s]");
   //both->Add(sdaq_ev);
   //both->Add(edaq_ev);
   //both->Draw(); 
   //both->Write();

   //TMultiGraph *both_hfcup = new TMultiGraph();
   //both_hfcup->SetTitle("Hel::scaler; Timestamp [s];Faraday Cup [nC]");
   //both_hfcup->Add(hel_fcup);
   //both_hfcup->Add(hel_fcupgated);
   //both_hfcup->Draw(); 
   //both_hfcup->Write();

   //TMultiGraph *both_hclock = new TMultiGraph();
   //both_hclock->SetTitle("Hel::scaler; Timestamp [s];Clock");
   //both_hclock->Add(hel_clock);
   //both_hclock->Add(hel_clockgated);
   //both_hclock->Draw(); 
   //both_hclock->Write();
   std::cout << e_tot_c << std::endl;
   std::cout << "---------------------" << std::endl;
   std::cout << "vx   mean  :   sigma" << std::endl;
   std::cout << "ind: " << vx_fit->GetParameter(1) << " : " << vx_fit->GetParameter(2) << std::endl;
   std::cout << "1: " << ivx_fit->GetParameter(1) << " : " << ivx_fit->GetParameter(2) << std::endl;
   std::cout << "2: " << iivx_fit->GetParameter(1) << " : " << iivx_fit->GetParameter(2) << std::endl;
   std::cout << "3: " << iiivx_fit->GetParameter(1) << " : " << iiivx_fit->GetParameter(2) << std::endl;
   std::cout << "4: " << ivvx_fit->GetParameter(1) << " : " << ivvx_fit->GetParameter(2) << std::endl;
   std::cout << "5: " << vvx_fit->GetParameter(1) << " : " << vvx_fit->GetParameter(2) << std::endl;
   std::cout << "6: " << vivx_fit->GetParameter(1) << " : " << vivx_fit->GetParameter(2) << std::endl;
   std::cout << "---------------------" << std::endl;
   std::cout << "vy   mean  :   sigma" << std::endl;
   std::cout << "ind: " << vy_fit->GetParameter(1) << " : " << vy_fit->GetParameter(2) << std::endl;
   std::cout << "1: " << ivy_fit->GetParameter(1) << " : " << ivy_fit->GetParameter(2) << std::endl;
   std::cout << "2: " << iivy_fit->GetParameter(1) << " : " << iivy_fit->GetParameter(2) << std::endl;
   std::cout << "3: " << iiivy_fit->GetParameter(1) << " : " << iiivy_fit->GetParameter(2) << std::endl;
   std::cout << "4: " << ivvy_fit->GetParameter(1) << " : " << ivvy_fit->GetParameter(2) << std::endl;
   std::cout << "5: " << vvy_fit->GetParameter(1) << " : " << vvy_fit->GetParameter(2) << std::endl;
   std::cout << "6: " << vivy_fit->GetParameter(1) << " : " << vivy_fit->GetParameter(2) << std::endl;
   std::cout << "---------------------" << std::endl;
   std::cout << "vz   mean  :   sigma" << std::endl;
//   std::cout << "enter: " << vz_fit->GetParameter(5) << " : " << vz_fit->GetParameter(6) <<   std::endl;
//   std::cout << "exit: "  << vz_fit->GetParameter(8) << " : " << vz_fit->GetParameter(9) <<  std::endl;
   std::cout << "peak: " << vz2_fit->GetParameter(1) << " : " << vz2_fit->GetParameter(2) <<  std::endl;
  // std::cout << "Sn: " << vz2_fit->GetParameter(4) << " : " << vz2_fit->GetParameter(5) <<  std::endl;
   std::cout << "---------------------" << std::endl;
   std::cout << "pip vx   mean  :   sigma" << std::endl;
   std::cout << "ind: " << pip_vx_fit->GetParameter(1) << " : " << pip_vx_fit->GetParameter(2) << std::endl;
   std::cout << "1: " << pip_ivx_fit->GetParameter(1) << " : " << pip_ivx_fit->GetParameter(2) << std::endl;
   std::cout << "2: " << pip_iivx_fit->GetParameter(1) << " : " << pip_iivx_fit->GetParameter(2) << std::endl;
   std::cout << "3: " << pip_iiivx_fit->GetParameter(1) << " : " << pip_iiivx_fit->GetParameter(2) << std::endl;
   std::cout << "4: " << pip_ivvx_fit->GetParameter(1) << " : " << pip_ivvx_fit->GetParameter(2) << std::endl;
   std::cout << "5: " << pip_vvx_fit->GetParameter(1) << " : " << pip_vvx_fit->GetParameter(2) << std::endl;
   std::cout << "6: " << pip_vivx_fit->GetParameter(1) << " : " << pip_vivx_fit->GetParameter(2) << std::endl;
   std::cout << "---------------------" << std::endl;
   std::cout << "pip vy   mean  :   sigma" << std::endl;
   std::cout << "ind: " << pip_vy_fit->GetParameter(1) << " : " << pip_vy_fit->GetParameter(2) << std::endl;
   std::cout << "1: " << pip_ivy_fit->GetParameter(1) << " : " << pip_ivy_fit->GetParameter(2) << std::endl;
   std::cout << "2: " << pip_iivy_fit->GetParameter(1) << " : " << pip_iivy_fit->GetParameter(2) << std::endl;
   std::cout << "3: " << pip_iiivy_fit->GetParameter(1) << " : " << pip_iiivy_fit->GetParameter(2) << std::endl;
   std::cout << "4: " << pip_ivvy_fit->GetParameter(1) << " : " << pip_ivvy_fit->GetParameter(2) << std::endl;
   std::cout << "5: " << pip_vvy_fit->GetParameter(1) << " : " << pip_vvy_fit->GetParameter(2) << std::endl;
   std::cout << "6: " << pip_vivy_fit->GetParameter(1) << " : " << pip_vivy_fit->GetParameter(2) << std::endl;
   std::cout << "---------------------" << std::endl;
   std::cout << "pip vz   mean  :   sigma" << std::endl;
   std::cout << "ind: " << pip_vz_fit->GetParameter(1) << " : " << pip_vz_fit->GetParameter(2) << " : " << pip_vz_fit->GetParameter(4) << " : " << pip_vz_fit->GetParameter(5) <<  std::endl;
   //std::cout << "1: " << pip_ivz_fit->GetParameter(1) << " : " << pip_ivz_fit->GetParameter(2) << " : " << pip_ivz_fit->GetParameter(4) << " : " << pip_ivz_fit->GetParameter(5) << std::endl;
   //std::cout << "2: " << pip_iivz_fit->GetParameter(1) << " : " << pip_iivz_fit->GetParameter(2) << " : " << pip_iivz_fit->GetParameter(4) << " : " << pip_iivz_fit->GetParameter(5) <<  std::endl;
   //std::cout << "3: " << pip_iiivz_fit->GetParameter(1) << " : " << pip_iiivz_fit->GetParameter(2) << " : " << pip_iiivz_fit->GetParameter(4) << " : " << pip_iiivz_fit->GetParameter(5) <<  std::endl;
   //std::cout << "4: " << pip_ivvz_fit->GetParameter(1) << " : " << pip_ivvz_fit->GetParameter(2) << " : " << pip_ivvz_fit->GetParameter(4) << " : " << pip_ivvz_fit->GetParameter(5) <<  std::endl;
   //std::cout << "5: " << pip_vvz_fit->GetParameter(1) << " : " << pip_vvz_fit->GetParameter(2) << " : " << pip_vvz_fit->GetParameter(4) << " : " << pip_vvz_fit->GetParameter(5) <<  std::endl;
   //std::cout << "6: " << pip_vivz_fit->GetParameter(1) << " : " << pip_vivz_fit->GetParameter(2) << " : " << pip_vivz_fit->GetParameter(4) << " : " << pip_vivz_fit->GetParameter(5) <<  std::endl;
   std::cout << "---------------------" << std::endl;
   std::cout << "pim vx   mean  :   sigma" << std::endl;
   std::cout << "ind: " << pim_vx_fit->GetParameter(1) << " : " << pim_vx_fit->GetParameter(2) << std::endl;
   std::cout << "1: " << pim_ivx_fit->GetParameter(1) << " : " << pim_ivx_fit->GetParameter(2) << std::endl;
   std::cout << "2: " << pim_iivx_fit->GetParameter(1) << " : " << pim_iivx_fit->GetParameter(2) << std::endl;
   std::cout << "3: " << pim_iiivx_fit->GetParameter(1) << " : " << pim_iiivx_fit->GetParameter(2) << std::endl;
   std::cout << "4: " << pim_ivvx_fit->GetParameter(1) << " : " << pim_ivvx_fit->GetParameter(2) << std::endl;
   std::cout << "5: " << pim_vvx_fit->GetParameter(1) << " : " << pim_vvx_fit->GetParameter(2) << std::endl;
   std::cout << "6: " << pim_vivx_fit->GetParameter(1) << " : " << pim_vivx_fit->GetParameter(2) << std::endl;
   std::cout << "---------------------" << std::endl;
   std::cout << "pim vy   mean  :   sigma" << std::endl;
   std::cout << "ind: " << pim_vy_fit->GetParameter(1) << " : " << pim_vy_fit->GetParameter(2) << std::endl;
   std::cout << "1: " << pim_ivy_fit->GetParameter(1) << " : " << pim_ivy_fit->GetParameter(2) << std::endl;
   std::cout << "2: " << pim_iivy_fit->GetParameter(1) << " : " << pim_iivy_fit->GetParameter(2) << std::endl;
   std::cout << "3: " << pim_iiivy_fit->GetParameter(1) << " : " << pim_iiivy_fit->GetParameter(2) << std::endl;
   std::cout << "4: " << pim_ivvy_fit->GetParameter(1) << " : " << pim_ivvy_fit->GetParameter(2) << std::endl;
   std::cout << "5: " << pim_vvy_fit->GetParameter(1) << " : " << pim_vvy_fit->GetParameter(2) << std::endl;
   std::cout << "6: " << pim_vivy_fit->GetParameter(1) << " : " << pim_vivy_fit->GetParameter(2) << std::endl;
   std::cout << "---------------------" << std::endl;
   std::cout << "pim vz   mean  :   sigma" << std::endl;
   std::cout << "ind: " << pim_vz_fit->GetParameter(1) << " : " << pim_vz_fit->GetParameter(2) << " : " << pim_vz_fit->GetParameter(4) << " : " << pim_vz_fit->GetParameter(5) <<  std::endl;
   //std::cout << "1: " << pim_ivz_fit->GetParameter(1) << " : " << pim_ivz_fit->GetParameter(2) << " : " << pim_ivz_fit->GetParameter(4) << " : " << pim_ivz_fit->GetParameter(5) << std::endl;
   //std::cout << "2: " << pim_iivz_fit->GetParameter(1) << " : " << pim_iivz_fit->GetParameter(2) << " : " << pim_iivz_fit->GetParameter(4) << " : " << pim_iivz_fit->GetParameter(5) <<  std::endl;
   //std::cout << "3: " << pim_iiivz_fit->GetParameter(1) << " : " << pim_iiivz_fit->GetParameter(2) << " : " << pim_iiivz_fit->GetParameter(4) << " : " << pim_iiivz_fit->GetParameter(5) <<  std::endl;
   //std::cout << "4: " << pim_ivvz_fit->GetParameter(1) << " : " << pim_ivvz_fit->GetParameter(2) << " : " << pim_ivvz_fit->GetParameter(4) << " : " << pim_ivvz_fit->GetParameter(5) <<  std::endl;
   //std::cout << "5: " << pim_vvz_fit->GetParameter(1) << " : " << pim_vvz_fit->GetParameter(2) << " : " << pim_vvz_fit->GetParameter(4) << " : " << pim_vvz_fit->GetParameter(5) <<  std::endl;
   //std::cout << "6: " << pim_vivz_fit->GetParameter(1) << " : " << pim_vivz_fit->GetParameter(2) << " : " << pim_vivz_fit->GetParameter(4) << " : " << pim_vivz_fit->GetParameter(5) <<  std::endl;
   std::cout << "---------------------" << std::endl;
   
   tr1_g->Draw();
   tr1_g->Write();
   
   tr2_g->Draw();
   tr2_g->Write();
   
   tr3_g->Draw();
   tr3_g->Write();
   
   tr4_g->Draw();
   tr4_g->Write();
   
   tr5_g->Draw();
   tr5_g->Write();
   
   tr6_g->Draw();
   tr6_g->Write();
   
   
   TCanvas *mgzc = new TCanvas("mgzc","mgzc",0,0,640,480);
   auto mgz = new TMultiGraph;
   	mgz->SetTitle("Vz vs Run Number; Run Number; Vz Peak Mean [cm]");
   	mgz->Add(tvz1_g);
   	mgz->Add(tvz2_g);
   	mgz->Add(tvz3_g);
   	mgz->Add(tvz4_g);
   	mgz->Add(tvz5_g);
   	mgz->Add(tvz6_g);
   	mgz->SetMinimum(-10.0);
   	mgz->SetMaximum(0);
   	//TAxis *axis = mgz->GetXaxis();
   	//axis->SetLimits(1,6);
   	mgz->Draw("AP");
   	mgzc->BuildLegend();
   	mgzc->Write();
   	
   	
   	TCanvas *mgrc = new TCanvas("mgrc","mgrc",0,0,640,480);
   auto mgr = new TMultiGraph;
   	mgr->SetTitle("#rho_{0} Yield vs Run Number; Run Number; #rho_{0} Yield");
   	mgr->Add(tr1_g);
   	mgr->Add(tr2_g);
   	mgr->Add(tr3_g);
   	mgr->Add(tr4_g);
   	mgr->Add(tr5_g);
   	mgr->Add(tr6_g);
   	//mgr->SetMinimum(-10.0);
   	//mgr->SetMaximum(0);
   	//TAxis *axis = mgz->GetXaxis();
   	//axis->SetLimits(1,6);
   	mgr->Draw("AP");
   	mgrc->BuildLegend();
   	mgrc->Write();
   
   std::cout << "e count " << e_count << std::endl;
   std::cout << "pip count " << pip_count << std::endl;
   std::cout << "pim count " << pim_count << std::endl;
   std::cout << "rho w/o cuts " << r_count << std::endl;
   std::cout << "rho w,lc cuts " << rw_count << std::endl;
   std::cout << "rho w,lc,t cuts " << rwt_count << std::endl;
   std::cout << "rho w,lc,t,z cuts " << rwtz_count << std::endl;
   std::cout << "bin1 " << r1_count << std::endl;
   std::cout << "bin2 " << r2_count << std::endl;
   std::cout << "bin3 " << r3_count << std::endl;
   std::cout << "bin4 " << r4_count << std::endl;
   std::cout << "bin5 " << r5_count << std::endl;
   std::cout << "bin6 " << r6_count << std::endl;
   std::cout << "ebin1 " << bin1_count << std::endl;
   std::cout << "ebin2 " << bin2_count << std::endl;
   std::cout << "ebin3 " << bin3_count << std::endl;
   std::cout << "ebin4 " << bin4_count << std::endl;
   std::cout << "ebin5 " << bin5_count << std::endl;
   std::cout << "ebin6 " << bin6_count << std::endl;
   
   std::cout << "e count q " << ecount_q << std::endl;
   std::cout << "e count full " << ecount_ful << std::endl;
   
   std::cout << " Total Charge = " << t_charge << std::endl;
   std::cout << "normalized charge = " << nt_charge << std::endl;
   
   std::cout << "Electron Outlier runs" << std::endl;
   for(int i = 0; i < ber; i++){
   std::cout << bad_e_run[i] << std::endl;
   }
   std::cout << "______________________" << std::endl;
   std::cout << "Pip Outlier runs" << std::endl;
   for(int i = 0; i < bpipr; i++){
   std::cout << bad_pip_run[i] << std::endl;
   }
   std::cout << "______________________" << std::endl;
   std::cout << "Pim Outlier runs" << std::endl;
   for(int i = 0; i < bpimr; i++){
   std::cout << bad_pim_run[i] << std::endl;
   }
   std::cout << "______________________" << std::endl;
   std::cout << "Rho Outlier runs" << std::endl;
   for(int i = 0; i < brr; i++){
   std::cout << bad_r_run[i] << std::endl;
   }
   
   std::cout << "___________________________________________________________________" << std::endl;
   std::cout << " Sector      |        a       |     bx        |       cx^2    |     dx^3   " << std::endl; 
   std::cout << "   1         | " << samp1_gfit->GetParameter(0) << " | " << samp1_gfit->GetParameter(1) << " | " << samp1_gfit->GetParameter(2) << " | " << samp1_gfit->GetParameter(3) << std::endl;
   std::cout << "   1u         | " << sampu1_gfit->GetParameter(0) << " | " << sampu1_gfit->GetParameter(1) << " | " << sampu1_gfit->GetParameter(2) << " | " << sampu1_gfit->GetParameter(3) << std::endl;
   std::cout << "   1d         | " << sampd1_gfit->GetParameter(0) << " | " << sampd1_gfit->GetParameter(1) << " | " << sampd1_gfit->GetParameter(2) << " | " << sampd1_gfit->GetParameter(3) << std::endl;
   std::cout << "___________________________________________________________________" << std::endl;
   std::cout << " Sector      |        a       |     bx        |       cx^2    |     dx^3   " << std::endl; 
   std::cout << "   2         | " << samp2_gfit->GetParameter(0) << " | " << samp2_gfit->GetParameter(1) << " | " << samp2_gfit->GetParameter(2) << " | " << samp2_gfit->GetParameter(3) << std::endl;
   std::cout << "   2u         | " << sampu2_gfit->GetParameter(0) << " | " << sampu2_gfit->GetParameter(1) << " | " << sampu2_gfit->GetParameter(2) << " | " << sampu2_gfit->GetParameter(3) << std::endl;
   std::cout << "   2d         | " << sampd2_gfit->GetParameter(0) << " | " << sampd2_gfit->GetParameter(1) << " | " << sampd2_gfit->GetParameter(2) << " | " << sampd2_gfit->GetParameter(3) << std::endl;
   std::cout << "___________________________________________________________________" << std::endl;
   std::cout << " Sector      |        a       |     bx        |       cx^2    |     dx^3   " << std::endl; 
   std::cout << "   3         | " << samp3_gfit->GetParameter(0) << " | " << samp3_gfit->GetParameter(1) << " | " << samp3_gfit->GetParameter(2) << " | " << samp3_gfit->GetParameter(3) << std::endl;
   std::cout << "   3u         | " << sampu3_gfit->GetParameter(0) << " | " << sampu3_gfit->GetParameter(1) << " | " << sampu3_gfit->GetParameter(2) << " | " << sampu3_gfit->GetParameter(3) << std::endl;
   std::cout << "   3d         | " << sampd3_gfit->GetParameter(0) << " | " << sampd3_gfit->GetParameter(1) << " | " << sampd3_gfit->GetParameter(2) << " | " << sampd3_gfit->GetParameter(3) << std::endl;
   std::cout << "___________________________________________________________________" << std::endl;
   std::cout << " Sector      |        a       |     bx        |       cx^2    |     dx^3   " << std::endl; 
   std::cout << "   4         | " << samp4_gfit->GetParameter(0) << " | " << samp4_gfit->GetParameter(1) << " | " << samp4_gfit->GetParameter(2) << " | " << samp4_gfit->GetParameter(3) << std::endl;
   std::cout << "   4u         | " << sampu4_gfit->GetParameter(0) << " | " << sampu4_gfit->GetParameter(1) << " | " << sampu4_gfit->GetParameter(2) << " | " << sampu4_gfit->GetParameter(3) << std::endl;
   std::cout << "   4d         | " << sampd4_gfit->GetParameter(0) << " | " << sampd4_gfit->GetParameter(1) << " | " << sampd4_gfit->GetParameter(2) << " | " << sampd4_gfit->GetParameter(3) << std::endl;
   std::cout << "___________________________________________________________________" << std::endl;
   std::cout << " Sector      |        a       |     bx        |       cx^2    |     dx^3   " << std::endl; 
   std::cout << "   5         | " << samp5_gfit->GetParameter(0) << " | " << samp5_gfit->GetParameter(1) << " | " << samp5_gfit->GetParameter(2) << " | " << samp5_gfit->GetParameter(3) << std::endl;
   std::cout << "   5u         | " << sampu5_gfit->GetParameter(0) << " | " << sampu5_gfit->GetParameter(1) << " | " << sampu5_gfit->GetParameter(2) << " | " << sampu5_gfit->GetParameter(3) << std::endl;
   std::cout << "   5d         | " << sampd5_gfit->GetParameter(0) << " | " << sampd5_gfit->GetParameter(1) << " | " << sampd5_gfit->GetParameter(2) << " | " << sampd5_gfit->GetParameter(3) << std::endl;
   std::cout << "___________________________________________________________________" << std::endl;
   std::cout << " Sector      |        a       |     bx        |       cx^2    |     dx^3   " << std::endl; 
   std::cout << "   6         | " << samp6_gfit->GetParameter(0) << " | " << samp6_gfit->GetParameter(1) << " | " << samp6_gfit->GetParameter(2) << " | " << samp6_gfit->GetParameter(3) << std::endl;
   std::cout << "   6u         | " << sampu6_gfit->GetParameter(0) << " | " << sampu6_gfit->GetParameter(1) << " | " << sampu6_gfit->GetParameter(2) << " | " << sampu6_gfit->GetParameter(3) << std::endl;
   std::cout << "   6d         | " << sampd6_gfit->GetParameter(0) << " | " << sampd6_gfit->GetParameter(1) << " | " << sampd6_gfit->GetParameter(2) << " | " << sampd6_gfit->GetParameter(3) << std::endl;
   std::cout << "___________________________________________________________________" << std::endl;
   std::cout << "   Thve         | " << thveu_fit->GetParameter(0) << " | " << thveu_fit->GetParameter(1) << std::endl;
   std::cout << "   sampd1         | " << sampd1_gfit->GetParameter(0) << " | " << sampd1_gfit->GetParameter(1) << " | " << sampd1_gfit->GetParameter(2) << " | " << sampd1_gfit->GetParameter(3) << std::endl;
   std::cout << "   sampu1         | " << sampu1_gfit->GetParameter(0) << " | " << sampu1_gfit->GetParameter(1) << " | " << sampu1_gfit->GetParameter(2) << " | " << sampu1_gfit->GetParameter(3) << std::endl;
   std::cout << "   chin 1-6      | " << chin1_cut << " | " << chin2_cut << " | " << chin3_cut << " | " << chin4_cut << " | " << chin5_cut << " | " << chin6_cut << " | " << std::endl;
   
   return 0;
}

/*
if(sect == 1) {
				ivx_hist->Fill(a_vx[row]);
   				ivy_hist->Fill(a_vy[row]);
				ixy_hist->Fill(a_px[row], a_py[row]);
                                if( (0.07204 * -3) <= a_vx[row] && (0.07204 * 3) >= a_vx[row] && (0.2734 * -3) <= a_vy[row] && (0.2734 * 3) >= a_vy[row]){
				ivz_hist->Fill(a_vz[row]);
				}
                        } else if(sect == 2) {
				iivx_hist->Fill(a_vx[row]);
   				iivy_hist->Fill(a_vy[row]);
				iixy_hist->Fill(a_px[row], a_py[row]);
				if( (0.111 * -3) <= a_vx[row] && (0.111 * 3) >= a_vx[row] && (0.2025 * -3) <= a_vy[row] && (0.2025 * 3) >= a_vy[row]){
				iivz_hist->Fill(a_vz[row]);
				}
                        } else if(sect == 3) {
				iiivx_hist->Fill(a_vx[row]);
   				iiivy_hist->Fill(a_vy[row]);
				iiixy_hist->Fill(a_px[row], a_py[row]);
				if( (0.1082 * -3) <= a_vx[row] && (0.1082 * 3) >= a_vx[row] && (0.2431 * -3) <= a_vy[row] && (0.2431 * 3) >= a_vy[row]){
				iiivz_hist->Fill(a_vz[row]);
				}
                        } else if(sect == 4) {
				ivvx_hist->Fill(a_vx[row]);
   				ivvy_hist->Fill(a_vy[row]);
				ivxy_hist->Fill(a_px[row], a_py[row]);
				if( (0.08755 * -3) <= a_vx[row] && (0.08755 * 3) >= a_vx[row] && (0.2828 * -3) <= a_vy[row] && (0.2828 * 3) >= a_vy[row]){
				ivvz_hist->Fill(a_vz[row]);
				}
                        } else if(sect == 5) {
				vvx_hist->Fill(a_vx[row]);
   				vvy_hist->Fill(a_vy[row]);
				vxy_hist->Fill(a_px[row], a_py[row]);
				if( (0.111 * -3) <= a_vx[row] && (0.111 * 3) >= a_vx[row] && (0.2089 * -3) <= a_vy[row] && (0.2089 * 3) >= a_vy[row]){
				vvz_hist->Fill(a_vz[row]);
				}
                        } else if(sect == 6) {
				vivx_hist->Fill(a_vx[row]);
   				vivy_hist->Fill(a_vy[row]);
				vixy_hist->Fill(a_px[row], a_py[row]);
				if( (0.1097 * -3) <= a_vx[row] && (0.1097 * 3) >= a_vx[row] && (0.242 * -3) <= a_vy[row] && (0.242 * 3) >= a_vy[row]){
				vivz_hist->Fill(a_vz[row]);
				}
*/

