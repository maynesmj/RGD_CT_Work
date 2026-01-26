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
  return par[0]*(1/(2*3.141592653))*(arg2/(4*arg3 + (arg4))) + par[3] + par[4]*pow(x[0],1) + par[5]*pow(x[0],2) + par[6]*pow(x[0],3); // + par[7]*pow(x[0],4) + par[8]*pow(x[0],5);
}

//"BreitWinger"

void fits(){


   TFile input("zzzzz_rgd_02_ld.root", "read");

   
    TFile f("zzzzzfits_ld_lc(0,0_5)_05.root", "update");

       

        TH1D *p_hist = (TH1D*)input.Get("Q^{2} (1,2) [GeV^{2}]");
      TH1D *ip_hist = (TH1D*)input.Get("Q^{2} (2,2.5) [GeV^{2}]");
      TH1D *iip_hist = (TH1D*)input.Get("Q^{2} (2.5,3) [GeV^{2}]");
      TH1D *iiip_hist = (TH1D*)input.Get("Q^{2} (3,3.5) [GeV^{2}]");
      TH1D *ivp_hist = (TH1D*)input.Get("Q^{2} (3.5,4.5) [GeV^{2}]");
      TH1D *vp_hist = (TH1D*)input.Get("Q^{2} (4.5,6) [GeV^{2}]");



   //TH2F *z_hist = new TH2F("px vs py", "px vs py", 200, -2, 2, 200, -2, 2);

   //auto cgr = new TGraph();
   //std::cout << "1 th" << std::endl;

   std::vector<std::string> v3_a;
   int num = 0;





   
  

   p_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *p_func = new TF1("p_fit",fitfbw,0.2,1.4,8);
   p_func->SetParameter(1,0.150);
   p_func->SetParameter(2,0.77);
   p_func->SetParLimits(2,0.65,0.89);
   p_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   p_hist->Fit("p_fit","","",0.2,1.4);
   p_hist->Draw();
   p_hist->Write();

   ip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *ip_func = new TF1("ip_fit",fitfbw,0.2,1.4,8);
   ip_func->SetParameter(1,0.150);
   ip_func->SetParameter(2,0.77);
   ip_func->SetParLimits(2,0.65,0.89);
   ip_func->SetParNames("double ipConstant","double ipgamma","double ipmean","double ipa1","double ipa2","double ipa3", "double ipa4", "double ipa5", "double ipa6");
   ip_hist->Fit("ip_fit","","",0.2,1.4);
   ip_hist->Draw();
   ip_hist->Write();


   iip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *iip_func = new TF1("iip_fit",fitfbw,0.2,1.4,8);
   iip_func->SetParameter(1,0.150);
   iip_func->SetParameter(2,0.77);
   iip_func->SetParLimits(2,0.65,0.89);
   iip_func->SetParNames("double iipConstant","double iipgamma","double iipmean","double iipa1","double iipa2","double iipa3", "double iipa4", "double iipa5", "double iipa6");
   iip_hist->Fit("iip_fit","","",0.2,1.4);
   iip_hist->Draw();
   iip_hist->Write();




 
   iiip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *iiip_func = new TF1("iiip_fit",fitfbw,0.2,1.4,8);
   iiip_func->SetParameter(1,0.150);
   iiip_func->SetParameter(2,0.77);
   iiip_func->SetParLimits(2,0.65,0.89);
   iiip_func->SetParNames("double iiipConstant","double iiipgamma","double iiipmean","double iiipa1","double iiipa2","double iiipa3", "double iiipa4", "double iiipa5", "double iiipa6");
   iiip_hist->Fit("iiip_fit","","",0.2,1.4);
   iiip_hist->Draw();
   iiip_hist->Write();



   ivp_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *ivp_func = new TF1("ivp_fit",fitfbw,0.2,1.4,8);
   ivp_func->SetParameter(1,0.150);
   ivp_func->SetParameter(2,0.77);
   ivp_func->SetParLimits(2,0.65,0.89);
   ivp_func->SetParNames("double ivpConstant","double ivpgamma","double ivpmean","double ivpa1","double ivpa2","double ivpa3", "double ivpa4", "double ivpa5", "double ivpa6");
   ivp_hist->Fit("ivp_fit","","",0.2,1.4);
   ivp_hist->Draw();
   ivp_hist->Write();


   vp_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}[GeV]");
   TF1 *vp_func = new TF1("vp_fit",fitfbw,0.2,1.4,8);
   vp_func->SetParameter(1,0.150);
   vp_func->SetParameter(2,0.77);
   vp_func->SetParLimits(2,0.65,0.89);
   vp_func->SetParNames("double vpConstant","double vpgamma","double vpmean","double vpa1","double vpa2","double vpa3", "double vpa4", "double vpa5", "double vpa6");
   vp_hist->Fit("vp_fit","","",0.2,1.4);
   vp_hist->Draw();
   vp_hist->Write();



}





/*
****************************************
         Invalid FitResult  (status = 2 )
****************************************

double pConstant          =      2168.77   +/-   0           
double pgamma             =         0.15   +/-   0           
double pmean              =     0.764356   +/-   0            	 (limited)
double pa1                =       221.98   +/-   0                  2608.33
double pa2                =     -4250.32   +/-   0           
double pa3                =      20014.4   +/-   0           
double pa4                =     -12098.3   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         Invalid FitResult  (status = 2 )
****************************************

double ipConstant         =      264.688   +/-   0           
double ipgamma            =         0.15   +/-   0           
double ipmean             =     0.765416   +/-   0            	 (limited)
double ipa1               =     -248.165   +/-   0           
double ipa2               =      1034.45   +/-   0                     318.479
double ipa3               =     -82.9882   +/-   0           
double ipa4               =     -333.572   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         Invalid FitResult  (status = 2 )
****************************************

double iipConstant        =      113.183   +/-   0           
double iipgamma           =         0.15   +/-   0           
double iipmean            =     0.765586   +/-   0            	 (limited)
double iipa1              =     -199.568   +/-   0           
double iipa2              =      904.387   +/-   0                       136.191
double iipa3              =     -675.768   +/-   0           
double iipa4              =       125.29   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         Invalid FitResult  (status = 2 )
****************************************

double iiipConstant       =      51.6768   +/-   0           
double iiipgamma          =         0.15   +/-   0           
double iiipmean           =     0.765298   +/-   0            	 (limited)
double iiipa1             =     -53.4725   +/-   0           
double iiipa2             =      237.568   +/-   0                     62.1762
double iiipa3             =      -106.97   +/-   0           
double iiipa4             =     -9.93169   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         Invalid FitResult  (status = 2 )
****************************************

double bivpConstant       =      34.6234   +/-   0           
double bivpgamma          =         0.15   +/-   0           
double bivpmean           =     0.766156   +/-   0            	 (limited)
double bivpa1             =     -59.6859   +/-   0           
double bivpa2             =      274.967   +/-   0                     41.674
double bivpa3             =     -234.613   +/-   0           
double bivpa4             =      65.1402   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         Invalid FitResult  (status = 2 )
****************************************

double vipConstant        =      6.58297   +/-   0           
double vipgamma           =         0.15   +/-   0           
double vipmean            =     0.759759   +/-   0            	 (limited)
double vipa1              =      -13.577   +/-   0           
double vipa2              =      41.0654   +/-   0                     7.89982
double vipa3              =      20.5305   +/-   0           
double vipa4              =     -32.2752   +/-   0  



cxc 






Warning in <Fit>: Abnormal termination of minimization.
Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
****************************************
         Invalid FitResult  (status = 2 )
****************************************

double pConstant          =       1644.4   +/-   0           
double pgamma             =         0.15   +/-   0           
double pmean              =     0.766751   +/-   0            	 (limited)
double pa1                =       565.07   +/-   0           
double pa2                =        -5248   +/-   0                     1979.67
double pa3                =      17053.4   +/-   0           
double pa4                =      -9567.3   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         Invalid FitResult  (status = 2 )
****************************************

double ipConstant         =       220.74   +/-   0           
double ipgamma            =         0.15   +/-   0           
double ipmean             =     0.766556   +/-   0            	 (limited)
double ipa1               =     -162.492   +/-   0           
double ipa2               =      587.376   +/-   0                     265.73
double ipa3               =      267.765   +/-   0           
double ipa4               =     -413.079   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         Invalid FitResult  (status = 2 )
****************************************

double iipConstant        =      95.3823   +/-   0           
double iipgamma           =         0.15   +/-   0           
double iipmean            =     0.764952   +/-   0            	 (limited)
double iipa1              =     -137.634   +/-   0           
double iipa2              =      572.158   +/-   0           
double iipa3              =     -312.751   +/-   0                      114.744
double iipa4              =      1.45923   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         Invalid FitResult  (status = 2 )
****************************************

double iiipConstant       =      44.6379   +/-   0           
double iiipgamma          =         0.15   +/-   0           
double iiipmean           =     0.765613   +/-   0            	 (limited)
double iiipa1             =      -89.693   +/-   0           
double iiipa2             =      413.154   +/-   0                     53.7143
double iiipa3             =     -399.421   +/-   0           
double iiipa4             =      129.067   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         Invalid FitResult  (status = 2 )
****************************************

double bivpConstant       =      31.8808   +/-   0           
double bivpgamma          =         0.15   +/-   0           
double bivpmean           =     0.768581   +/-   0            	 (limited)
double bivpa1             =     -44.9654   +/-   0           
double bivpa2             =      190.866   +/-   0                     38.4098
double bivpa3             =     -138.588   +/-   0           
double bivpa4             =      27.5065   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
         Invalid FitResult  (status = 2 )
****************************************

double vipConstant        =       7.2101   +/-   0           
double vipgamma           =         0.15   +/-   0           
double vipmean            =     0.762585   +/-   0            	 (limited)
double vipa1              =     -12.6212   +/-   0           
double vipa2              =      34.4824   +/-   0                     8.66438
double vipa3              =      19.1262   +/-   0           
double vipa4              =     -27.3664   +/-   0 
*/



//    fits from zzzzz_cxc_01_02    ---------------------------------------------------------------------------------------------


/*
                                                                        32.232
double pConstant          =      26.9952   +/-   0           
double pgamma             =         0.15   +/-   0           
double pmean              =     0.750966   +/-   0            	 (limited)
double pa1                =     -174.071   +/-   0           
double pa2                =      888.866   +/-   0           
double pa3                =     -921.516   +/-   0           
double pa4                =      278.481   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
                                                                       23.2851
double ipConstant         =      19.4452   +/-   0           
double ipgamma            =         0.15   +/-   0           
double ipmean             =     0.755786   +/-   0            	 (limited)
double ipa1               =     -85.6629   +/-   0           
double ipa2               =      466.017   +/-   0           
double ipa3               =     -512.568   +/-   0           
double ipa4               =        166.3   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
                                                                      14.9409
double iipConstant        =      12.5008   +/-   0           
double iipgamma           =         0.15   +/-   0           
double iipmean            =     0.752569   +/-   0            	 (limited)
double iipa1              =     -60.1673   +/-   0           
double iipa2              =      302.643   +/-   0           
double iipa3              =     -331.836   +/-   0           
double iipa4              =      108.301   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
                                                                      8.40333
double iiipConstant       =      7.00644   +/-   0           
double iiipgamma          =         0.15   +/-   0           
double iiipmean           =     0.758687   +/-   0            	 (limited)
double iiipa1             =     -23.6006   +/-   0           
double iiipa2             =      123.858   +/-   0           
double iiipa3             =     -125.047   +/-   0           
double iiipa4             =      35.9879   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
                                                                     9.222
double ivpConstant        =       7.6562   +/-   0           
double ivpgamma           =         0.15   +/-   0           
double ivpmean            =     0.767991   +/-   0            	 (limited)
double ivpa1              =     -23.9404   +/-   0           
double ivpa2              =      112.399   +/-   0           
double ivpa3              =     -116.806   +/-   0           
double ivpa4              =      37.3145   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
                                                                      3.49978
double vpConstant         =      2.91422   +/-   0           
double vpgamma            =         0.15   +/-   0           
double vpmean             =     0.761251   +/-   0            	 (limited)
double vpa1               =     -6.88773   +/-   0           
double vpa2               =      29.6491   +/-   0           
double vpa3               =     -21.2073   +/-   0           
double vpa4               =      3.42055   +/-   0  
*/

//    fits from zzzzz_ld_01_02    ---------------------------------------------------------------------------------------------

                                             
/*
double pConstant          =        47.97   +/-   0                    57.1249       
double pgamma             =         0.15   +/-   0           
double pmean              =     0.747041   +/-   0            	 (limited)
double pa1                =     -352.474   +/-   0           
double pa2                =      1808.31   +/-   0           
double pa3                =     -1974.76   +/-   0           
double pa4                =       633.27   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
    
double ipConstant         =      29.0397   +/-   0                     34.6447     
double ipgamma            =         0.15   +/-   0           
double ipmean             =     0.749702   +/-   0            	 (limited)
double ipa1               =     -166.354   +/-   0           
double ipa2               =      864.303   +/-   0           
double ipa3               =      -954.72   +/-   0           
double ipa4               =      311.225   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
                                                                      20.9128
double iipConstant        =      17.4752   +/-   0           
double iipgamma           =         0.15   +/-   0           
double iipmean            =     0.754685   +/-   0            	 (limited)
double iipa1              =     -98.0208   +/-   0           
double iipa2              =      508.792   +/-   0           
double iipa3              =     -579.644   +/-   0           
double iipa4              =      196.833   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
                                                                      11.862
double iiipConstant       =      9.89416   +/-   0           
double iiipgamma          =         0.15   +/-   0           
double iiipmean           =     0.757929   +/-   0            	 (limited)
double iiipa1             =     -40.3156   +/-   0           
double iiipa2             =      206.736   +/-   0           
double iiipa3             =     -226.197   +/-   0           
double iiipa4             =      73.9983   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
                                                                     10.7233
double ivpConstant        =      8.93473   +/-   0           
double ivpgamma           =         0.15   +/-   0           
double ivpmean            =     0.760009   +/-   0            	 (limited)
double ivpa1              =       -35.68   +/-   0           
double ivpa2              =      171.031   +/-   0           
double ivpa3              =     -177.432   +/-   0           
double ivpa4              =      55.8486   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
                                                                         3.80744
double vpConstant         =      3.17238   +/-   0           
double vpgamma            =         0.15   +/-   0           
double vpmean             =     0.760066   +/-   0            	 (limited)
double vpa1               =     -9.53226   +/-   0           
double vpa2               =      41.0644   +/-   0           
double vpa3               =     -34.8885   +/-   0           
double vpa4               =      8.89049   +/-   0  
*/


