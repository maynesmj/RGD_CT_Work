//******************************************************************
//*  ██╗  ██╗██╗██████╗  ██████╗     ██╗  ██╗    ██████╗
//*  ██║  ██║██║██╔══██╗██╔═══██╗    ██║  ██║   ██╔═████╗
//*  ███████║██║██████╔╝██║   ██║    ███████║   ██║██╔██║
//*  ██╔══██║██║██╔═══╝ ██║   ██║    ╚════██║   ████╔╝██║
//*  ██║  ██║██║██║     ╚██████╔╝         ██║██╗╚██████╔╝
//*  ╚═╝  ╚═╝╚═╝╚═╝      ╚═════╝          ╚═╝╚═╝ ╚═════╝
//************************ Jefferson National Lab (2017) ***********
//******************************************************************
//* Example program for reading HIPO-4 Files..
//* Reads the file created by writeFile program
//*--
//* Author: G.Gavalian
//*

// This version of hipo file reader will take 25,000 events and save highest energy electron, +/- pi, and mass of rho in root trees


#include <cstdlib>
#include <iostream>
#include <string>
//#include "reader.h" 
//#include "TFile.h"
//#include "TTree.h"
//#include "TChain.h"
//#include "TLorentzVector.h"
//#include <filesystem>
//#include <sys/stat.h>
//#include <vector>

//#include "TTree.h"
//#include <TH1D.h>
//#include <TH1F.h>
//#include <TH2F.h>
//#include "TGraph.h"
//#include "TMultiGraph.h"
//#include "TCanvas.h"
//#include "TStyle.h"
//#include "TLine.h"
//#include "TLatex.h"
//#include <THStack.h>
//#include "TF1.h"
//#include "TMath.h"


//namespace fs = std::filesystem;

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
  //Double_t arg1 = 1.0/4.0; // 2 over pi
  Double_t arg2 = 0.150; //Gamma=par[1]  M=par[2]
  Double_t arg3 = (par[2] - x[0])*(par[2] - x[0]); //((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = (0.150*0.150);
  return par[0]*(arg2/(4*arg3 + arg4)) + par[3] + par[4]*pow(x[0],1) + par[5]*pow(x[0],2) + par[6]*pow(x[0],3); // + par[7]*pow(x[0],4) + par[8]*pow(x[0],5);
}

//"BreitWinger"

void fitslc(){


   TFile input("zzzzz_rgd_02_ld.root", "read");

   
    TFile f("zzzzzzfits_ld_lc_02.root", "update");

       

        TH1D *ilc_hist = (TH1D*)input.Get("lc (0,0.5) [fm]");
      TH1D *iilc_hist = (TH1D*)input.Get("lc (0.5,0.75) [fm]");
      TH1D *iiilc_hist = (TH1D*)input.Get("lc (0.75,1) [fm]");
      TH1D *ivlc_hist = (TH1D*)input.Get("lc (1,1.25) [fm]");
      TH1D *vlc_hist = (TH1D*)input.Get("lc (1.25,1.5) [fm]");
      TH1D *vilc_hist = (TH1D*)input.Get("lc (1.5,1.75) [fm]");
      TH1D *viilc_hist = (TH1D*)input.Get("lc (1.75,2.0) [fm]");
      TH1D *viiilc_hist = (TH1D*)input.Get("lc (2,2.5) [fm]");
      TH1D *ixlc_hist = (TH1D*)input.Get("lc (2.5,3) [fm]");
      TH1D *xlc_hist = (TH1D*)input.Get("lc (3,4) [fm]");
      TH1D *xilc_hist = (TH1D*)input.Get("lc (4,6) [fm]");
      


   //TH2F *z_hist = new TH2F("px vs py", "px vs py", 200, -2, 2, 200, -2, 2);

   //auto cgr = new TGraph();
   //std::cout << "1 th" << std::endl;

   std::vector<std::string> v3_a;
   int num = 0;





   
  

   ilc_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *ilc_func = new TF1("ilc_fit",fitfbw,0.2,1.4,9);
   ilc_func->SetParameter(1,0.150);
   ilc_func->SetParameter(2,0.77);
   ilc_func->SetParLimits(2,0.65,0.89);
   ilc_func->SetParNames("double ilcConstant","double ilcgamma","double ilcmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   ilc_hist->Fit("ilc_fit","","",0.2,1.4);
   ilc_hist->Draw();
   ilc_hist->Write();

   iilc_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *iilc_func = new TF1("iilc_fit",fitfbw,0.2,1.4,9);
   iilc_func->SetParameter(1,0.150);
   iilc_func->SetParameter(2,0.77);
   iilc_func->SetParLimits(2,0.65,0.89);
   iilc_func->SetParNames("double iilcConstant","double iilcgamma","double iilcmean","double ipa1","double ipa2","double ipa3", "double ipa4", "double ipa5", "double ipa6");
   iilc_hist->Fit("iilc_fit","","",0.2,1.4);
   iilc_hist->Draw();
   iilc_hist->Write();


   iiilc_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *iiilc_func = new TF1("iiilc_fit",fitfbw,0.2,1.4,9);
   iiilc_func->SetParameter(1,0.150);
   iiilc_func->SetParameter(2,0.77);
   iiilc_func->SetParLimits(2,0.65,0.89);
   iiilc_func->SetParNames("double iiilcConstant","double iiilcgamma","double iiilcmean","double iipa1","double iipa2","double iipa3", "double iipa4", "double iipa5", "double iipa6");
   iiilc_hist->Fit("iiilc_fit","","",0.2,1.4);
   iiilc_hist->Draw();
   iiilc_hist->Write();





   ivlc_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *ivlc_func = new TF1("ivlc_fit",fitfbw,0.2,1.4,9);
   ivlc_func->SetParameter(1,0.150);
   ivlc_func->SetParameter(2,0.77);
   ivlc_func->SetParLimits(2,0.65,0.89);
   ivlc_func->SetParNames("double ivlcConstant","double ivlcgamma","double ivlcmean","double ivpa1","double ivpa2","double ivpa3", "double ivpa4", "double ivpa5", "double ivpa6");
   ivlc_hist->Fit("ivlc_fit","","",0.2,1.4);
   ivlc_hist->Draw();
   ivlc_hist->Write();


   vlc_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *vlc_func = new TF1("vlc_fit",fitfbw,0.2,1.4,9);
   vlc_func->SetParameter(1,0.150);
   vlc_func->SetParameter(2,0.77);
   vlc_func->SetParLimits(2,0.65,0.89);
   vlc_func->SetParNames("double vlcConstant","double vlcgamma","double vlcmean","double vpa1","double vpa2","double vpa3", "double vpa4", "double vpa5", "double vpa6");
   vlc_hist->Fit("vlc_fit","","",0.2,1.4);
   vlc_hist->Draw();
   vlc_hist->Write();


   vilc_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *vilc_func = new TF1("vilc_fit",fitfbw,0.2,1.4,9);
   vilc_func->SetParameter(1,0.150);
   vilc_func->SetParameter(2,0.77);
   vilc_func->SetParLimits(2,0.65,0.89);
   vilc_func->SetParNames("double vilcConstant","double vilcgamma","double vilcmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   vilc_hist->Fit("vilc_fit","","",0.2,1.4);
   vilc_hist->Draw();
   vilc_hist->Write();

   viilc_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *viilc_func = new TF1("viilc_fit",fitfbw,0.2,1.4,9);
   viilc_func->SetParameter(1,0.150);
   viilc_func->SetParameter(2,0.77);
   viilc_func->SetParLimits(2,0.65,0.89);
   viilc_func->SetParNames("double viilcConstant","double viilcgamma","double viilcmean","double ipa1","double ipa2","double ipa3", "double ipa4", "double ipa5", "double ipa6");
   viilc_hist->Fit("viilc_fit","","",0.2,1.4);
   viilc_hist->Draw();
   viilc_hist->Write();


   viiilc_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *viiilc_func = new TF1("viiilc_fit",fitfbw,0.2,1.4,9);
   viiilc_func->SetParameter(1,0.150);
   viiilc_func->SetParameter(2,0.77);
   viiilc_func->SetParLimits(2,0.65,0.89);
   viiilc_func->SetParNames("double viiilcConstant","double viiilcgamma","double viiilcmean","double iipa1","double iipa2","double iipa3", "double iipa4", "double iipa5", "double iipa6");
   viiilc_hist->Fit("viiilc_fit","","",0.2,1.4);
   viiilc_hist->Draw();
   viiilc_hist->Write();





   ixlc_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *ixlc_func = new TF1("ixlc_fit",fitfbw,0.2,1.4,9);
   ixlc_func->SetParameter(1,0.150);
   ixlc_func->SetParameter(2,0.77);
   ixlc_func->SetParLimits(2,0.65,0.89);
   ixlc_func->SetParNames("double ixlcConstant","double ixlcgamma","double ixlcmean","double ivpa1","double ivpa2","double ivpa3", "double ivpa4", "double ivpa5", "double ivpa6");
   ixlc_hist->Fit("ixlc_fit","","",0.2,1.4);
   ixlc_hist->Draw();
   ixlc_hist->Write();


   xlc_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *xlc_func = new TF1("xlc_fit",fitfbw,0.2,1.4,9);
   xlc_func->SetParameter(1,0.150);
   xlc_func->SetParameter(2,0.77);
   xlc_func->SetParLimits(2,0.65,0.89);
   xlc_func->SetParNames("double xlcConstant","double xlcgamma","double xlcmean","double vpa1","double vpa2","double vpa3", "double vpa4", "double vpa5", "double vpa6");
   xlc_hist->Fit("xlc_fit","","",0.2,1.4);
   xlc_hist->Draw();
   xlc_hist->Write();
   
   xilc_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *xilc_func = new TF1("xilc_fit",fitfbw,0.2,1.4,9);
   xilc_func->SetParameter(1,0.150);
   xilc_func->SetParameter(2,0.77);
   xilc_func->SetParLimits(2,0.65,0.89);
   xilc_func->SetParNames("double xilcConstant","double xilcgamma","double xilcmean","double vpa1","double vpa2","double vpa3", "double vpa4", "double vpa5", "double vpa6");
   xilc_hist->Fit("xilc_fit","","",0.2,1.4);
   xilc_hist->Draw();
   xilc_hist->Write();


}





// lc fits ld -------------------------------------------------------------------------------

/*
****************************************
 
double ilcConstant        =      116.244   +/-   0           
double ilcgamma           =         0.15   +/-   0           
double ilcmean            =     0.751406   +/-   0            	 (limited)
double pa1                =     -710.328   +/-   0           
double pa2                =      3659.76   +/-   0           
double pa3                =     -4014.87   +/-   0           
double pa4                =      1303.09   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
       
double iilcConstant       =      838.656   +/-   0           
double iilcgamma          =         0.15   +/-   0           
double iilcmean           =     0.763414   +/-   0            	 (limited)
double ipa1               =     -1423.66   +/-   0           
double ipa2               =      5947.95   +/-   0           
double ipa3               =      -2950.1   +/-   0           
double ipa4               =     -266.445   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
       
double iiilcConstant      =      460.093   +/-   0           
double iiilcgamma         =         0.15   +/-   0           
double iiilcmean          =     0.770876   +/-   0            	 (limited)
double iipa1              =     -11.5613   +/-   0           
double iipa2              =     -932.826   +/-   0           
double iipa3              =      3533.57   +/-   0           
double iipa4              =     -1867.03   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         
double ivlcConstant       =      172.829   +/-   0           
double ivlcgamma          =         0.15   +/-   0           
double ivlcmean           =     0.772275   +/-   0            	 (limited)
double ivpa1              =      27.9105   +/-   0           
double ivpa2              =     -544.065   +/-   0           
double ivpa3              =      1506.67   +/-   0           
double ivpa4              =     -696.404   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         
double vlcConstant        =      59.3194   +/-   0           
double vlcgamma           =         0.15   +/-   0           
double vlcmean            =      0.77407   +/-   0            	 (limited)
double vpa1               =      38.6272   +/-   0           
double vpa2               =     -330.635   +/-   0           
double vpa3               =      673.783   +/-   0           
double vpa4               =     -271.197   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         
double vilcConstant       =      17.8224   +/-   0           
double vilcgamma          =         0.15   +/-   0           
double vilcmean           =     0.780049   +/-   0            	 (limited)
double pa1                =      39.3967   +/-   0           
double pa2                =     -226.468   +/-   0           
double pa3                =      354.394   +/-   0           
double pa4                =     -128.448   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
        
double viilcConstant      =      3.72275   +/-   0           
double viilcgamma         =         0.15   +/-   0           
double viilcmean          =     0.792714   +/-   0            	 (limited)
double ipa1               =      23.0316   +/-   0           
double ipa2               =      -107.73   +/-   0           
double ipa3               =      144.845   +/-   0           
double ipa4               =     -48.2236   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         
double viiilcConstant     =     0.426149   +/-   0           
double viiilcgamma        =         0.15   +/-   0           
double viiilcmean         =      0.79274   +/-   0            	 (limited)
double iipa1              =      3.08569   +/-   0           
double iipa2              =     -5.09386   +/-   0           
double iipa3              =     -2.31915   +/-   0           
double iipa4              =       6.2969   +/-   0           
Warning in <Fit>: Fit data is empty 
Warning in <Fit>: Fit data is empty 
Warning in <Fit>: Fit data is empty 
*/


