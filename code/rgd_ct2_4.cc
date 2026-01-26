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
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"
#include <THStack.h>
#include "TF1.h"
#include "TMath.h"


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
   TH1F *z_hist = new TH1F("Zh", " Zh ", 200, 0, 10);
   //bC_hist->SetStats(0);

   TH2F *qz_hist = new TH2F("lc vs Q2", "lc vs Q2", 200, 0, 10, 200, 0, 7);
   
   TH2F *xy_hist = new TH2F("px vs py", "px vs py", 200, -5, 5, 200, -5, 5);

   TH1F *q_hist = new TH1F("Q2", " Q2", 200, 0, 10);
   
   TH1F *lc_hist = new TH1F("lc", " lc", 200, 0, 10);
   
   TH1F *e_hist = new TH1F("electrons", " electrons ", 200, -300, 300);
   
   TH1F *pip_hist = new TH1F("pip", "pip", 200, -300, 300);
   
   TH1F *pim_hist = new TH1F("pim", "pim", 200, -300, 300);

   //TH2F *qlc_hist = new TH2F("Q2 vs lc", "Q2 vs lc", 200, 0, 10, 200, 0, 20);

   TH1F *w_hist = new TH1F("w (Invariant Mass)", " ", 200, 0, 5);

   //TH2F *tz_hist = new TH2F("-t vs rho_pz", "-t vs rho_pz", 200, 0, 10, 200, 0, 20);

   TH1F *t_hist = new TH1F("-t", "-t", 200, 0, 20);

   TH1F *wonphe_hist = new TH1F("wonphe", " ", 500, 0, 0.5);

   TH1F *wnphe_hist = new TH1F("wnphe", " ", 500, 0, 0.5);

   TH1F *wopcal_hist = new TH1F("Total Energy", " ", 500, 0, 1);

   TH1F *wpcal_hist = new TH1F("PCAL", " ", 500, 0, 1);

   TH2F *wocut_hist = new TH2F("wocut", "wocut", 500, 0, 20, 500, 0, 1);

   TH2F *wcut_hist = new TH2F("Energy", "", 500, 0.6, 11, 500, 0, 0.3);
   
   TH2F *e_pht_hist = new TH2F("e phi vs theta", "e phi vs theta", 200, -180.0, 180.0, 200, 0, 40);
   
   TH2F *pip_pht_hist = new TH2F("pip phi vs theta", "pip phi vs theta", 200, -180.0, 180.0, 200, 0, 40);
   
   TH2F *pim_pht_hist = new TH2F("pim phi vs theta", "pim phi vs theta", 200, -180.0, 180.0, 200, 0, 40);

   auto *vxs = new THStack("vxs", "vxs");
   auto *vys = new THStack("vys", "vys");
   auto *etots = new THStack("etots", "etots");
   auto *ps = new THStack("ps" , "ps");

float lx13 = 0.00617297 - 3*(0.0695994);   
float ux13 = 0.00617297 + 3*(0.0695994); 
float lx14 = 0.00617297 - 4*(0.0695994);   
float ux14 = 0.00617297 + 4*(0.0695994); 
float lx15 = 0.00617297 - 5*(0.0695994);   
float ux15 = 0.00617297 + 5*(0.0695994); 


float ly13 = -0.319247 - 3*(0.287735);   
float uy13 = -0.319247 + 3*(0.287735); 
float ly14 = -0.319247 - 4*(0.287735);   
float uy14 = -0.319247 + 4*(0.287735); 
float ly15 = -0.319247 - 5*(0.287735);   
float uy15 = -0.319247 + 5*(0.287735); 


float lx23 = 0.301558 - 3*(0.324802);   
float ux23 = 0.301558 + 3*(0.324802); 
float lx24 = 0.301558 - 4*(0.324802);   
float ux24 = 0.301558 + 4*(0.324802); 
float lx25 = 0.301558 - 5*(0.324802);   
float ux25 = 0.301558 + 5*(0.324802); 
   

float ly23 = -0.309686 - 3*(0.55718);   
float uy23 = -0.309686 + 3*(0.55718); 
float ly24 = -0.309686 - 4*(0.55718);   
float uy24 = -0.309686 + 4*(0.55718); 
float ly25 = -0.309686 - 5*(0.55718);   
float uy25 = -0.309686 + 5*(0.55718); 


float lx33 = 0.0687002 - 3*(0.28029);   
float ux33 = 0.0687002 + 3*(0.28029); 
float lx34 = 0.0687002 - 4*(0.28029);   
float ux34 = 0.0687002 + 4*(0.28029); 
float lx35 = 0.0687002 - 5*(0.28029);   
float ux35 = 0.0687002 + 5*(0.28029); 

 
float ly33 = -0.5418 - 3*(0.195704);   
float uy33 = -0.5418 + 3*(0.195704); 
float ly34 = -0.5418 - 4*(0.195704);   
float uy34 = -0.5418 + 4*(0.195704);
float ly35 = -0.5418 - 5*(0.195704);   
float uy35 = -0.5418 + 5*(0.195704);


float lx43 = -0.0365462 - 3*(0.117509);   
float ux43 = -0.0365462 + 3*(0.117509);
float lx44 = -0.0365462 - 4*(0.117509);   
float ux44 = -0.0365462 + 4*(0.117509);
float lx45 = -0.0365462 - 5*(0.117509);   
float ux45 = -0.0365462 + 5*(0.117509);


float ly43 = 0.267435 - 3*(0.402091);   
float uy43 = 0.267435 + 3*(0.402091);
float ly44 = 0.267435 - 4*(0.402091);   
float uy44 = 0.267435 + 4*(0.402091);
float ly45 = 0.267435 - 5*(0.402091);   
float uy45 = 0.267435 + 5*(0.402091);


float lx53 = -0.277040 - 3*(0.289325);   
float ux53 = -0.277040 + 3*(0.289325);
float lx54 = -0.277040 - 4*(0.289325);   
float ux54 = -0.277040 + 4*(0.289325);
float lx55 = -0.277040 - 5*(0.289325);   
float ux55 = -0.277040 + 5*(0.289325);


float ly53 = -0.0392163 - 3*(0.161128);   
float uy53 = -0.0392163 + 3*(0.161128);
float ly54 = -0.0392163 - 4*(0.161128);   
float uy54 = -0.0392163 + 4*(0.161128);
float ly55 = -0.0392163 - 5*(0.161128);   
float uy55 = -0.0392163 + 5*(0.161128);


float lx63 = -0.0359976 - 3*(0.266055);   
float ux63 = -0.0359976 + 3*(0.266055);
float lx64 = -0.0359976 - 4*(0.266055);   
float ux64 = -0.0359976 + 4*(0.266055);
float lx65 = -0.0359976 - 5*(0.266055);   
float ux65 = -0.0359976 + 5*(0.266055);


float ly63 = -0.184030  - 3*(0.191211);   
float uy63 = -0.184030  + 3*(0.191211);
float ly64 = -0.184030  - 4*(0.191211);   
float uy64 = -0.184030  + 4*(0.191211);
float ly65 = -0.184030  - 5*(0.191211);   
float uy65 = -0.184030  + 5*(0.191211);
                 float lxlim3 = 0.00656862 - 3*(0.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float uxlim3 = 0.00656862 + 3*(0.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float lylim3 = -0.1793 - 3*(0.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float uylim3 = -0.1793 + 3*(0.46594);
                 float lxlim4 = 0.00656862 - 4*(0.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float uxlim4 = 0.00656862 + 4*(0.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float lylim4 = -0.1793 - 4*(0.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float uylim4 = -0.1793 + 4*(0.46594);
                 float lxlim5 = 0.00656862 - 5*(0.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float uxlim5 = 0.00656862 + 5*(0.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float lylim5 = -0.1793 - 5*(0.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float uylim5 = -0.1793 + 5*(0.46594);

float vx13_c = 0.0;
float vx14_c = 0.0;
float vx15_c = 0.0;
float vxlim13_c = 0.0;
float vxlim14_c = 0.0;
float vxlim15_c = 0.0;
float vy13_c = 0.0;
float vy14_c = 0.0;
float vy15_c = 0.0;
float vylim13_c = 0.0;
float vylim14_c = 0.0;
float vylim15_c = 0.0;
float vz13_c = 0.0;
float vz14_c = 0.0;
float vz15_c = 0.0;
float vzlim13_c = 0.0;
float vzlim14_c = 0.0;
float vzlim15_c = 0.0;

float vx23_c = 0.0;
float vx24_c = 0.0;
float vx25_c = 0.0;
float vxlim23_c = 0.0;
float vxlim24_c = 0.0;
float vxlim25_c = 0.0;
float vy23_c = 0.0;
float vy24_c = 0.0;
float vy25_c = 0.0;
float vylim23_c = 0.0;
float vylim24_c = 0.0;
float vylim25_c = 0.0;
float vz23_c = 0.0;
float vz24_c = 0.0;
float vz25_c = 0.0;
float vzlim23_c = 0.0;
float vzlim24_c = 0.0;
float vzlim25_c = 0.0;

float vx33_c = 0.0;
float vx34_c = 0.0;
float vx35_c = 0.0;
float vxlim33_c = 0.0;
float vxlim34_c = 0.0;
float vxlim35_c = 0.0;
float vy33_c = 0.0;
float vy34_c = 0.0;
float vy35_c = 0.0;
float vylim33_c = 0.0;
float vylim34_c = 0.0;
float vylim35_c = 0.0;
float vz33_c = 0.0;
float vz34_c = 0.0;
float vz35_c = 0.0;
float vzlim33_c = 0.0;
float vzlim34_c = 0.0;
float vzlim35_c = 0.0;

float vx43_c = 0.0;
float vx44_c = 0.0;
float vx45_c = 0.0;
float vxlim43_c = 0.0;
float vxlim44_c = 0.0;
float vxlim45_c = 0.0;
float vy43_c = 0.0;
float vy44_c = 0.0;
float vy45_c = 0.0;
float vylim43_c = 0.0;
float vylim44_c = 0.0;
float vylim45_c = 0.0;
float vz43_c = 0.0;
float vz44_c = 0.0;
float vz45_c = 0.0;
float vzlim43_c = 0.0;
float vzlim44_c = 0.0;
float vzlim45_c = 0.0;

float vx53_c = 0.0;
float vx54_c = 0.0;
float vx55_c = 0.0;
float vxlim53_c = 0.0;
float vxlim54_c = 0.0;
float vxlim55_c = 0.0;
float vy53_c = 0.0;
float vy54_c = 0.0;
float vy55_c = 0.0;
float vylim53_c = 0.0;
float vylim54_c = 0.0;
float vylim55_c = 0.0;
float vz53_c = 0.0;
float vz54_c = 0.0;
float vz55_c = 0.0;
float vzlim53_c = 0.0;
float vzlim54_c = 0.0;
float vzlim55_c = 0.0;

float vx63_c = 0.0;
float vx64_c = 0.0;
float vx65_c = 0.0;
float vxlim63_c = 0.0;
float vxlim64_c = 0.0;
float vxlim65_c = 0.0;
float vy63_c = 0.0;
float vy64_c = 0.0;
float vy65_c = 0.0;
float vylim63_c = 0.0;
float vylim64_c = 0.0;
float vylim65_c = 0.0;
float vz63_c = 0.0;
float vz64_c = 0.0;
float vz65_c = 0.0;
float vzlim63_c = 0.0;
float vzlim64_c = 0.0;
float vzlim65_c = 0.0;

float vx1_c = 0.0;
float vy1_c = 0.0;
float vz1_c = 0.0;
float vx2_c = 0.0;
float vy2_c = 0.0;
float vz2_c = 0.0;
float vx3_c = 0.0;
float vy3_c = 0.0;
float vz3_c = 0.0;
float vx4_c = 0.0;
float vy4_c = 0.0;
float vz4_c = 0.0;
float vx5_c = 0.0;
float vy5_c = 0.0;
float vz5_c = 0.0;
float vx6_c = 0.0;
float vy6_c = 0.0;
float vz6_c = 0.0;

float lvy5c = -7.58642 - 1.3 - 3*(0.527298);   
float uvy5c = -7.58642 - 1.3 + 3*(0.527298); 
float lvx6c = -7.58642 - 1.3 - 4*(0.527298);   
float uvx6c = -7.58642 - 1.3 + 4*(0.527298); 
float lvy6c = -7.58642 - 1.3 - 5*(0.527298);   
float uvy6c = -7.58642 - 1.3 + 5*(0.527298); 


float ls13 = -2.65302 - 1.3 - 3*(0.544745);   
float us13 = -2.65302 - 1.3 + 3*(0.544745); 
float ls14 = -2.65302 - 1.3 - 4*(0.544745);   
float us14 = -2.65302 - 1.3 + 4*(0.544745); 
float ls15 = -2.65302 - 1.3 - 5*(0.544745);   
float us15 = -2.65302 - 1.3 + 5*(0.544745); 


float lc23 = -7.42106 - 1.3 - 3*(0.622596);   
float uc23 = -7.42106 - 1.3 + 3*(0.622596); 
float lc24 = -7.42106 - 1.3 - 4*(0.622596);   
float uc24 = -7.42106 - 1.3 + 4*(0.622596); 
float lc25 = -7.42106 - 1.3 - 5*(0.622596);   
float uc25 = -7.42106 - 1.3 + 5*(0.622596); 
   

float ls23 = -2.51326 - 1.3 - 3*(0.59944);   
float us23 = -2.51326 - 1.3 + 3*(0.59944); 
float ls24 = -2.51326 - 1.3 - 4*(0.59944);   
float us24 = -2.51326 - 1.3 + 4*(0.59944); 
float ls25 = -2.51326 - 1.3 - 5*(0.59944);   
float us25 = -2.51326 - 1.3 + 5*(0.59944); 


float lc33 = -7.41344 - 1.3 - 3*(0.572976);   
float uc33 = -7.41344 - 1.3 + 3*(0.572976); 
float lc34 = -7.41344 - 1.3 - 4*(0.572976);   
float uc34 = -7.41344 - 1.3 + 4*(0.572976); 
float lc35 = -7.41344 - 1.3 - 5*(0.572976);   
float uc35 = -7.41344 - 1.3 + 5*(0.572976); 

 
float ls33 = -2.48824 - 1.3 - 3*(0.54446);   
float us33 = -2.48824 - 1.3 + 3*(0.54446); 
float ls34 = -2.48824 - 1.3 - 4*(0.54446);   
float us34 = -2.48824 - 1.3 + 4*(0.54446);
float ls35 = -2.48824 - 1.3 - 5*(0.54446);   
float us35 = -2.48824 - 1.3 + 5*(0.54446);


float lvx1c3 = -7.60018 - 1.3 - 3*(0.571731);   
float uvx1c3 = -7.60018 - 1.3 + 3*(0.571731);
float lvx1c4 = -7.60018 - 1.3 - 4*(0.571731);   
float uvx1c4 = -7.60018 - 1.3 + 4*(0.571731);
float lvx1vy1c = -7.60018 - 1.3 - 5*(0.571731);   
float uvx1vy1c = -7.60018 - 1.3 + 5*(0.571731);


float ls43 = -2.66478 - 1.3 - 3*(0.559822);   
float us43 = -2.66478 - 1.3 + 3*(0.559822);
float ls44 = -2.66478 - 1.3 - 4*(0.559822);   
float us44 = -2.66478 - 1.3 + 4*(0.559822);
float ls45 = -2.66478 - 1.3 - 5*(0.559822);   
float us45 = -2.66478 - 1.3 + 5*(0.559822);


float lvy1c3 = -7.44553 - 1.3 - 3*(0.569488);   
float uvy1c3 = -7.44553 - 1.3 + 3*(0.569488);
float lvy1c4 = -7.44553 - 1.3 - 4*(0.569488);   
float uvy1c4 = -7.44553 - 1.3 + 4*(0.569488);
float lvy1c5 = -7.44553 - 1.3 - 5*(0.569488);   
float uvy1c5 = -7.44553 - 1.3 + 5*(0.569488);


float ls53 = -2.52733 - 1.3 - 3*(0.563511);   
float us53 = -2.52733 - 1.3 + 3*(0.563511);
float ls54 = -2.52733 - 1.3 - 4*(0.563511);   
float us54 = -2.52733 - 1.3 + 4*(0.563511);
float ls55 = -2.52733 - 1.3 - 5*(0.563511);   
float us55 = -2.52733 - 1.3 + 5*(0.563511);


float lvx2c3 = -7.57334 - 1.3 - 3*(0.515616);   
float uvx2c3 = -7.57334 - 1.3 + 3*(0.515616);
float lvx2c4 = -7.57334 - 1.3 - 4*(0.515616);   
float uvx2c4 = -7.57334 - 1.3 + 4*(0.515616);
float lvx2c5 = -7.57334 - 1.3 - 5*(0.515616);   
float uvx2c5 = -7.57334 - 1.3 + 5*(0.515616);

float ls63 = -2.64149  - 1.3 - 3*(0.508674);   
float us63 = -2.64149  - 1.3 + 3*(0.508674);
float ls64 = -2.64149  - 1.3 - 4*(0.508674);   
float us64 = -2.64149  - 1.3 + 4*(0.508674);
float ls65 = -2.64149  - 1.3 - 5*(0.508674);   
float us65 = -2.64149  - 1.3 + 5*(0.508674);

float lc3 = -7.50864 - 1.3 - 3*(0.5685);   
float uc3 = -7.50864 - 1.3 + 3*(0.5685);
float lvx1c = -7.50864 - 1.3 - 4*(0.5685);   
float uvx1c = -7.50864 - 1.3 + 4*(0.5685);
float lvy1c = -7.50864 - 1.3 - 5*(0.5685);   
float uvy1c = -7.50864 - 1.3 + 5*(0.5685);

float ls3 = -2.58252 - 1.3  - 3*(0.559604);   
float us3 = -2.58252 - 1.3  + 3*(0.559604);
float ls4 = -2.58252 - 1.3  - 4*(0.559604);   
float us4 = -2.58252 - 1.3  + 4*(0.559604);
float ls5 = -2.58252 - 1.3  - 5*(0.559604);   
float us5 = -2.58252 - 1.3  + 5*(0.559604);


float pip_lx13 = 10.00617297 - 3*(10.0695994);   
float pip_ux13 = 10.00617297 + 3*(10.0695994); 
float pip_lx14 = 10.00617297 - 4*(10.0695994);   
float pip_ux14 = 10.00617297 + 4*(10.0695994); 
float pip_lx15 = 10.00617297 - 5*(10.0695994);   
float pip_ux15 = 10.00617297 + 5*(10.0695994); 


float pip_ly13 = -10.319247 - 3*(10.287735);   
float pip_uy13 = -10.319247 + 3*(10.287735); 
float pip_ly14 = -10.319247 - 4*(10.287735);   
float pip_uy14 = -10.319247 + 4*(10.287735); 
float pip_ly15 = -10.319247 - 5*(10.287735);   
float pip_uy15 = -10.319247 + 5*(10.287735); 


float pip_lx23 = 10.1301558 - 3*(10.324802);   
float pip_ux23 = 10.1301558 + 3*(10.324802); 
float pip_lx24 = 10.1301558 - 4*(10.324802);   
float pip_ux24 = 10.1301558 + 4*(10.324802); 
float pip_lx25 = 10.1301558 - 5*(10.324802);   
float pip_ux25 = 10.1301558 + 5*(10.324802); 
   

float pip_ly23 = -10.309686 - 3*(10.55718);   
float pip_uy23 = -10.309686 + 3*(10.55718); 
float pip_ly24 = -10.309686 - 4*(10.55718);   
float pip_uy24 = -10.309686 + 4*(10.55718); 
float pip_ly25 = -10.309686 - 5*(10.55718);   
float pip_uy25 = -10.309686 + 5*(10.55718); 


float pip_lx33 = 10.0687002 - 3*(10.28029);   
float pip_ux33 = 10.0687002 + 3*(10.28029); 
float pip_lx34 = 10.0687002 - 4*(10.28029);   
float pip_ux34 = 10.0687002 + 4*(10.28029); 
float pip_lx35 = 10.0687002 - 5*(10.28029);   
float pip_ux35 = 10.0687002 + 5*(10.28029); 

 
float pip_ly33 = -10.5418 - 3*(10.195704);   
float pip_uy33 = -10.5418 + 3*(10.195704); 
float pip_ly34 = -10.5418 - 4*(10.195704);   
float pip_uy34 = -10.5418 + 4*(10.195704);
float pip_ly35 = -10.5418 - 5*(10.195704);   
float pip_uy35 = -10.5418 + 5*(10.195704);


float pip_lx43 = -10.0365462 - 3*(10.117509);   
float pip_ux43 = -10.0365462 + 3*(10.117509);
float pip_lx44 = -10.0365462 - 4*(10.117509);   
float pip_ux44 = -10.0365462 + 4*(10.117509);
float pip_lx45 = -10.0365462 - 5*(10.117509);   
float pip_ux45 = -10.0365462 + 5*(10.117509);


float pip_ly43 = 10.267435 - 3*(10.402091);   
float pip_uy43 = 10.267435 + 3*(10.402091);
float pip_ly44 = 10.267435 - 4*(10.402091);   
float pip_uy44 = 10.267435 + 4*(10.402091);
float pip_ly45 = 10.267435 - 5*(10.402091);   
float pip_uy45 = 10.267435 + 5*(10.402091);


float pip_lx53 = -10.277040 - 3*(10.289325);   
float pip_ux53 = -10.277040 + 3*(10.289325);
float pip_lx54 = -10.277040 - 4*(10.289325);   
float pip_ux54 = -10.277040 + 4*(10.289325);
float pip_lx55 = -10.277040 - 5*(10.289325);   
float pip_ux55 = -10.277040 + 5*(10.289325);


float pip_ly53 = -10.0392163 - 3*(10.161128);   
float pip_uy53 = -10.0392163 + 3*(10.161128);
float pip_ly54 = -10.0392163 - 4*(10.161128);   
float pip_uy54 = -10.0392163 + 4*(10.161128);
float pip_ly55 = -10.0392163 - 5*(10.161128);   
float pip_uy55 = -10.0392163 + 5*(10.161128);


float pip_lx63 = -10.0359976 - 3*(10.266055);   
float pip_ux63 = -10.0359976 + 3*(10.266055);
float pip_lx64 = -10.0359976 - 4*(10.266055);   
float pip_ux64 = -10.0359976 + 4*(10.266055);
float pip_lx65 = -10.0359976 - 5*(10.266055);   
float pip_ux65 = -10.0359976 + 5*(10.266055);


float pip_ly63 = -10.184030  - 3*(10.191211);   
float pip_uy63 = -10.184030  + 3*(10.191211);
float pip_ly64 = -10.184030  - 4*(10.191211);   
float pip_uy64 = -10.184030  + 4*(10.191211);
float pip_ly65 = -10.184030  - 5*(10.191211);   
float pip_uy65 = -10.184030  + 5*(10.191211);
                 float pip_lxlim3 = 10.00656862 - 3*(10.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float pip_uxlim3 = 10.00656862 + 3*(10.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float pip_lylim3 = -10.1793 - 3*(10.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float pip_uylim3 = -10.1793 + 3*(10.46594);
                 float pip_lxlim4 = 10.00656862 - 4*(10.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float pip_uxlim4 = 10.00656862 + 4*(10.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float pip_lylim4 = -10.1793 - 4*(10.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float pip_uylim4 = -10.1793 + 4*(10.46594);
                 float pip_lxlim5 = 10.00656862 - 5*(10.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float pip_uxlim5 = 10.00656862 + 5*(10.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float pip_lylim5 = -10.1793 - 5*(10.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float pip_uylim5 = -10.1793 + 5*(10.46594);

float pip_vx13_c = 0.0;
float pip_vx14_c = 0.0;
float pip_vx15_c = 0.0;
float pip_vxlim13_c = 0.0;
float pip_vxlim14_c = 0.0;
float pip_vxlim15_c = 0.0;
float pip_vy13_c = 0.0;
float pip_vy14_c = 0.0;
float pip_vy15_c = 0.0;
float pip_vylim13_c = 0.0;
float pip_vylim14_c = 0.0;
float pip_vylim15_c = 0.0;
float pip_vz13_c = 0.0;
float pip_vz14_c = 0.0;
float pip_vz15_c = 0.0;
float pip_vzlim13_c = 0.0;
float pip_vzlim14_c = 0.0;
float pip_vzlim15_c = 0.0;

float pip_vx23_c = 0.0;
float pip_vx24_c = 0.0;
float pip_vx25_c = 0.0;
float pip_vxlim23_c = 0.0;
float pip_vxlim24_c = 0.0;
float pip_vxlim25_c = 0.0;
float pip_vy23_c = 0.0;
float pip_vy24_c = 0.0;
float pip_vy25_c = 0.0;
float pip_vylim23_c = 0.0;
float pip_vylim24_c = 0.0;
float pip_vylim25_c = 0.0;
float pip_vz23_c = 0.0;
float pip_vz24_c = 0.0;
float pip_vz25_c = 0.0;
float pip_vzlim23_c = 0.0;
float pip_vzlim24_c = 0.0;
float pip_vzlim25_c = 0.0;

float pip_vx33_c = 0.0;
float pip_vx34_c = 0.0;
float pip_vx35_c = 0.0;
float pip_vxlim33_c = 0.0;
float pip_vxlim34_c = 0.0;
float pip_vxlim35_c = 0.0;
float pip_vy33_c = 0.0;
float pip_vy34_c = 0.0;
float pip_vy35_c = 0.0;
float pip_vylim33_c = 0.0;
float pip_vylim34_c = 0.0;
float pip_vylim35_c = 0.0;
float pip_vz33_c = 0.0;
float pip_vz34_c = 0.0;
float pip_vz35_c = 0.0;
float pip_vzlim33_c = 0.0;
float pip_vzlim34_c = 0.0;
float pip_vzlim35_c = 0.0;

float pip_vx43_c = 0.0;
float pip_vx44_c = 0.0;
float pip_vx45_c = 0.0;
float pip_vxlim43_c = 0.0;
float pip_vxlim44_c = 0.0;
float pip_vxlim45_c = 0.0;
float pip_vy43_c = 0.0;
float pip_vy44_c = 0.0;
float pip_vy45_c = 0.0;
float pip_vylim43_c = 0.0;
float pip_vylim44_c = 0.0;
float pip_vylim45_c = 0.0;
float pip_vz43_c = 0.0;
float pip_vz44_c = 0.0;
float pip_vz45_c = 0.0;
float pip_vzlim43_c = 0.0;
float pip_vzlim44_c = 0.0;
float pip_vzlim45_c = 0.0;

float pip_vx53_c = 0.0;
float pip_vx54_c = 0.0;
float pip_vx55_c = 0.0;
float pip_vxlim53_c = 0.0;
float pip_vxlim54_c = 0.0;
float pip_vxlim55_c = 0.0;
float pip_vy53_c = 0.0;
float pip_vy54_c = 0.0;
float pip_vy55_c = 0.0;
float pip_vylim53_c = 0.0;
float pip_vylim54_c = 0.0;
float pip_vylim55_c = 0.0;
float pip_vz53_c = 0.0;
float pip_vz54_c = 0.0;
float pip_vz55_c = 0.0;
float pip_vzlim53_c = 0.0;
float pip_vzlim54_c = 0.0;
float pip_vzlim55_c = 0.0;

float pip_vx63_c = 0.0;
float pip_vx64_c = 0.0;
float pip_vx65_c = 0.0;
float pip_vxlim63_c = 0.0;
float pip_vxlim64_c = 0.0;
float pip_vxlim65_c = 0.0;
float pip_vy63_c = 0.0;
float pip_vy64_c = 0.0;
float pip_vy65_c = 0.0;
float pip_vylim63_c = 0.0;
float pip_vylim64_c = 0.0;
float pip_vylim65_c = 0.0;
float pip_vz63_c = 0.0;
float pip_vz64_c = 0.0;
float pip_vz65_c = 0.0;
float pip_vzlim63_c = 0.0;
float pip_vzlim64_c = 0.0;
float pip_vzlim65_c = 0.0;

float pip_vx1_c = 0.0;
float pip_vy1_c = 0.0;
float pip_vz1_c = 0.0;
float pip_vx2_c = 0.0;
float pip_vy2_c = 0.0;
float pip_vz2_c = 0.0;
float pip_vx3_c = 0.0;
float pip_vy3_c = 0.0;
float pip_vz3_c = 0.0;
float pip_vx4_c = 0.0;
float pip_vy4_c = 0.0;
float pip_vz4_c = 0.0;
float pip_vx5_c = 0.0;
float pip_vy5_c = 0.0;
float pip_vz5_c = 0.0;
float pip_vx6_c = 0.0;
float pip_vy6_c = 0.0;
float pip_vz6_c = 0.0;

float pip_lvy5c = -17.58642 - 1.3 - 3*(10.527298);   
float pip_uvy5c = -17.58642 - 1.3 + 3*(10.527298); 
float pip_lvx6c = -17.58642 - 1.3 - 4*(10.527298);   
float pip_uvx6c = -17.58642 - 1.3 + 4*(10.527298); 
float pip_lvy6c = -17.58642 - 1.3 - 5*(10.527298);   
float pip_uvy6c = -17.58642 - 1.3 + 5*(10.527298); 


float pip_ls13 = -12.65302 - 1.3 - 3*(10.544745);   
float pip_us13 = -12.65302 - 1.3 + 3*(10.544745); 
float pip_ls14 = -12.65302 - 1.3 - 4*(10.544745);   
float pip_us14 = -12.65302 - 1.3 + 4*(10.544745); 
float pip_ls15 = -12.65302 - 1.3 - 5*(10.544745);   
float pip_us15 = -12.65302 - 1.3 + 5*(10.544745); 


float pip_lc23 = -17.42106 - 1.3 - 3*(10.622596);   
float pip_uc23 = -17.42106 - 1.3 + 3*(10.622596); 
float pip_lc24 = -17.42106 - 1.3 - 4*(10.622596);   
float pip_uc24 = -17.42106 - 1.3 + 4*(10.622596); 
float pip_lc25 = -17.42106 - 1.3 - 5*(10.622596);   
float pip_uc25 = -17.42106 - 1.3 + 5*(10.622596); 
   

float pip_ls23 = -12.51326 - 1.3 - 3*(10.59944);   
float pip_us23 = -12.51326 - 1.3 + 3*(10.59944); 
float pip_ls24 = -12.51326 - 1.3 - 4*(10.59944);   
float pip_us24 = -12.51326 - 1.3 + 4*(10.59944); 
float pip_ls25 = -12.51326 - 1.3 - 5*(10.59944);   
float pip_us25 = -12.51326 - 1.3 + 5*(10.59944); 


float pip_lc33 = -17.41344 - 1.3 - 3*(10.572976);   
float pip_uc33 = -17.41344 - 1.3 + 3*(10.572976); 
float pip_lc34 = -17.41344 - 1.3 - 4*(10.572976);   
float pip_uc34 = -17.41344 - 1.3 + 4*(10.572976); 
float pip_lc35 = -17.41344 - 1.3 - 5*(10.572976);   
float pip_uc35 = -17.41344 - 1.3 + 5*(10.572976); 

 
float pip_ls33 = -12.48824 - 1.3 - 3*(10.54446);   
float pip_us33 = -12.48824 - 1.3 + 3*(10.54446); 
float pip_ls34 = -12.48824 - 1.3 - 4*(10.54446);   
float pip_us34 = -12.48824 - 1.3 + 4*(10.54446);
float pip_ls35 = -12.48824 - 1.3 - 5*(10.54446);   
float pip_us35 = -12.48824 - 1.3 + 5*(10.54446);


float pip_lvx1c3 = -17.60018 - 1.3 - 3*(10.571731);   
float pip_uvx1c3 = -17.60018 - 1.3 + 3*(10.571731);
float pip_lvx1c4 = -17.60018 - 1.3 - 4*(10.571731);   
float pip_uvx1c4 = -17.60018 - 1.3 + 4*(10.571731);
float pip_lvx1c5 = -17.60018 - 1.3 - 5*(10.571731);   
float pip_uvx1c5 = -17.60018 - 1.3 + 5*(10.571731);


float pip_ls43 = -12.66478 - 1.3 - 3*(10.559822);   
float pip_us43 = -12.66478 - 1.3 + 3*(10.559822);
float pip_ls44 = -12.66478 - 1.3 - 4*(10.559822);   
float pip_us44 = -12.66478 - 1.3 + 4*(10.559822);
float pip_ls45 = -12.66478 - 1.3 - 5*(10.559822);   
float pip_us45 = -12.66478 - 1.3 + 5*(10.559822);


float pip_lvy1c3 = -17.44553 - 1.3 - 3*(10.569488);   
float pip_uvy1c3 = -17.44553 - 1.3 + 3*(10.569488);
float pip_lvy1c4 = -17.44553 - 1.3 - 4*(10.569488);   
float pip_uvy1c4 = -17.44553 - 1.3 + 4*(10.569488);
float pip_lvy1c5 = -17.44553 - 1.3 - 5*(10.569488);   
float pip_uvy1c5 = -17.44553 - 1.3 + 5*(10.569488);


float pip_ls53 = -12.52733 - 1.3 - 3*(10.563511);   
float pip_us53 = -12.52733 - 1.3 + 3*(10.563511);
float pip_ls54 = -12.52733 - 1.3 - 4*(10.563511);   
float pip_us54 = -12.52733 - 1.3 + 4*(10.563511);
float pip_ls55 = -12.52733 - 1.3 - 5*(10.563511);   
float pip_us55 = -12.52733 - 1.3 + 5*(10.563511);


float pip_lvx2c3 = -17.57334 - 1.3 - 3*(10.515616);   
float pip_uvx2c3 = -17.57334 - 1.3 + 3*(10.515616);
float pip_lvx2c4 = -17.57334 - 1.3 - 4*(10.515616);   
float pip_uvx2c4 = -17.57334 - 1.3 + 4*(10.515616);
float pip_lvx2c5 = -17.57334 - 1.3 - 5*(10.515616);   
float pip_uvx2c5 = -17.57334 - 1.3 + 5*(10.515616);

float pip_ls63 = -12.64149  - 1.3 - 3*(10.508674);   
float pip_us63 = -12.64149  - 1.3 + 3*(10.508674);
float pip_ls64 = -12.64149  - 1.3 - 4*(10.508674);   
float pip_us64 = -12.64149  - 1.3 + 4*(10.508674);
float pip_ls65 = -12.64149  - 1.3 - 5*(10.508674);   
float pip_us65 = -12.64149  - 1.3 + 5*(10.508674);

float pip_lc3 = -17.50864 - 1.3 - 3*(10.5685);   
float pip_uc3 = -17.50864 - 1.3 + 3*(10.5685);
float pip_lvx1c = -17.50864 - 1.3 - 4*(10.5685);   
float pip_uvx1c = -17.50864 - 1.3 + 4*(10.5685);
float pip_lvy1c = -17.50864 - 1.3 - 5*(10.5685);   
float pip_uvy1c = -17.50864 - 1.3 + 5*(10.5685);

float pip_ls3 = -12.58252 - 1.3  - 3*(10.559604);   
float pip_us3 = -12.58252 - 1.3  + 3*(10.559604);
float pip_ls4 = -12.58252 - 1.3  - 4*(10.559604);   
float pip_us4 = -12.58252 - 1.3  + 4*(10.559604);
float pip_ls5 = -12.58252 - 1.3  - 5*(10.559604);   
float pip_us5 = -12.58252 - 1.3  + 5*(10.559604);




float pim_lx13 = 20.00617297 - 3*(20.0695994);   
float pim_ux13 = 20.00617297 + 3*(20.0695994); 
float pim_lx14 = 20.00617297 - 4*(20.0695994);   
float pim_ux14 = 20.00617297 + 4*(20.0695994); 
float pim_lx15 = 20.00617297 - 5*(20.0695994);   
float pim_ux15 = 20.00617297 + 5*(20.0695994); 


float pim_ly13 = -20.319247 - 3*(20.287735);   
float pim_uy13 = -20.319247 + 3*(20.287735); 
float pim_ly14 = -20.319247 - 4*(20.287735);   
float pim_uy14 = -20.319247 + 4*(20.287735); 
float pim_ly15 = -20.319247 - 5*(20.287735);   
float pim_uy15 = -20.319247 + 5*(20.287735); 


float pim_lx23 = 20.301558 - 3*(20.324802);   
float pim_ux23 = 20.301558 + 3*(20.324802); 
float pim_lx24 = 20.301558 - 4*(20.324802);   
float pim_ux24 = 20.301558 + 4*(20.324802); 
float pim_lx25 = 20.301558 - 5*(20.324802);   
float pim_ux25 = 20.301558 + 5*(20.324802); 
   

float pim_ly23 = -20.309686 - 3*(20.55718);   
float pim_uy23 = -20.309686 + 3*(20.55718); 
float pim_ly24 = -20.309686 - 4*(20.55718);   
float pim_uy24 = -20.309686 + 4*(20.55718); 
float pim_ly25 = -20.309686 - 5*(20.55718);   
float pim_uy25 = -20.309686 + 5*(20.55718); 


float pim_lx33 = 20.0687002 - 3*(20.28029);   
float pim_ux33 = 20.0687002 + 3*(20.28029); 
float pim_lx34 = 20.0687002 - 4*(20.28029);   
float pim_ux34 = 20.0687002 + 4*(20.28029); 
float pim_lx35 = 20.0687002 - 5*(20.28029);   
float pim_ux35 = 20.0687002 + 5*(20.28029); 

 
float pim_ly33 = -20.5418 - 3*(20.195704);   
float pim_uy33 = -20.5418 + 3*(20.195704); 
float pim_ly34 = -20.5418 - 4*(20.195704);   
float pim_uy34 = -20.5418 + 4*(20.195704);
float pim_ly35 = -20.5418 - 5*(20.195704);   
float pim_uy35 = -20.5418 + 5*(20.195704);


float pim_lx43 = -20.0365462 - 3*(20.117509);   
float pim_ux43 = -20.0365462 + 3*(20.117509);
float pim_lx44 = -20.0365462 - 4*(20.117509);   
float pim_ux44 = -20.0365462 + 4*(20.117509);
float pim_lx45 = -20.0365462 - 5*(20.117509);   
float pim_ux45 = -20.0365462 + 5*(20.117509);


float pim_ly43 = 20.267435 - 3*(20.402091);   
float pim_uy43 = 20.267435 + 3*(20.402091);
float pim_ly44 = 20.267435 - 4*(20.402091);   
float pim_uy44 = 20.267435 + 4*(20.402091);
float pim_ly45 = 20.267435 - 5*(20.402091);   
float pim_uy45 = 20.267435 + 5*(20.402091);


float pim_lx53 = -20.277040 - 3*(20.289325);   
float pim_ux53 = -20.277040 + 3*(20.289325);
float pim_lx54 = -20.277040 - 4*(20.289325);   
float pim_ux54 = -20.277040 + 4*(20.289325);
float pim_lx55 = -20.277040 - 5*(20.289325);   
float pim_ux55 = -20.277040 + 5*(20.289325);


float pim_ly53 = -20.0392163 - 3*(20.161128);   
float pim_uy53 = -20.0392163 + 3*(20.161128);
float pim_ly54 = -20.0392163 - 4*(20.161128);   
float pim_uy54 = -20.0392163 + 4*(20.161128);
float pim_ly55 = -20.0392163 - 5*(20.161128);   
float pim_uy55 = -20.0392163 + 5*(20.161128);


float pim_lx63 = -20.0359976 - 3*(20.266055);   
float pim_ux63 = -20.0359976 + 3*(20.266055);
float pim_lx64 = -20.0359976 - 4*(20.266055);   
float pim_ux64 = -20.0359976 + 4*(20.266055);
float pim_lx65 = -20.0359976 - 5*(20.266055);   
float pim_ux65 = -20.0359976 + 5*(20.266055);


float pim_ly63 = -20.184030  - 3*(20.191211);   
float pim_uy63 = -20.184030  + 3*(20.191211);
float pim_ly64 = -20.184030  - 4*(20.191211);   
float pim_uy64 = -20.184030  + 4*(20.191211);
float pim_ly65 = -20.184030  - 5*(20.191211);   
float pim_uy65 = -20.184030  + 5*(20.191211);
                 float pim_lxlim3 = 20.00656862 - 3*(20.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float pim_uxlim3 = 20.00656862 + 3*(20.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float pim_lylim3 = -20.1793 - 3*(20.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float pim_uylim3 = -20.1793 + 3*(20.46594);
                 float pim_lxlim4 = 20.00656862 - 4*(20.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float pim_uxlim4 = 20.00656862 + 4*(20.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float pim_lylim4 = -20.1793 - 4*(20.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float pim_uylim4 = -20.1793 + 4*(20.46594);
                 float pim_lxlim5 = 20.00656862 - 5*(20.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float pim_uxlim5 = 20.00656862 + 5*(20.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float pim_lylim5 = -20.1793 - 5*(20.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float pim_uylim5 = -20.1793 + 5*(20.46594);

float pim_vx13_c = 0.0;
float pim_vx14_c = 0.0;
float pim_vx15_c = 0.0;
float pim_vxlim13_c = 0.0;
float pim_vxlim14_c = 0.0;
float pim_vxlim15_c = 0.0;
float pim_vy13_c = 0.0;
float pim_vy14_c = 0.0;
float pim_vy15_c = 0.0;
float pim_vylim13_c = 0.0;
float pim_vylim14_c = 0.0;
float pim_vylim15_c = 0.0;
float pim_vz13_c = 0.0;
float pim_vz14_c = 0.0;
float pim_vz15_c = 0.0;
float pim_vzlim13_c = 0.0;
float pim_vzlim14_c = 0.0;
float pim_vzlim15_c = 0.0;

float pim_vx23_c = 0.0;
float pim_vx24_c = 0.0;
float pim_vx25_c = 0.0;
float pim_vxlim23_c = 0.0;
float pim_vxlim24_c = 0.0;
float pim_vxlim25_c = 0.0;
float pim_vy23_c = 0.0;
float pim_vy24_c = 0.0;
float pim_vy25_c = 0.0;
float pim_vylim23_c = 0.0;
float pim_vylim24_c = 0.0;
float pim_vylim25_c = 0.0;
float pim_vz23_c = 0.0;
float pim_vz24_c = 0.0;
float pim_vz25_c = 0.0;
float pim_vzlim23_c = 0.0;
float pim_vzlim24_c = 0.0;
float pim_vzlim25_c = 0.0;

float pim_vx33_c = 0.0;
float pim_vx34_c = 0.0;
float pim_vx35_c = 0.0;
float pim_vxlim33_c = 0.0;
float pim_vxlim34_c = 0.0;
float pim_vxlim35_c = 0.0;
float pim_vy33_c = 0.0;
float pim_vy34_c = 0.0;
float pim_vy35_c = 0.0;
float pim_vylim33_c = 0.0;
float pim_vylim34_c = 0.0;
float pim_vylim35_c = 0.0;
float pim_vz33_c = 0.0;
float pim_vz34_c = 0.0;
float pim_vz35_c = 0.0;
float pim_vzlim33_c = 0.0;
float pim_vzlim34_c = 0.0;
float pim_vzlim35_c = 0.0;

float pim_vx43_c = 0.0;
float pim_vx44_c = 0.0;
float pim_vx45_c = 0.0;
float pim_vxlim43_c = 0.0;
float pim_vxlim44_c = 0.0;
float pim_vxlim45_c = 0.0;
float pim_vy43_c = 0.0;
float pim_vy44_c = 0.0;
float pim_vy45_c = 0.0;
float pim_vylim43_c = 0.0;
float pim_vylim44_c = 0.0;
float pim_vylim45_c = 0.0;
float pim_vz43_c = 0.0;
float pim_vz44_c = 0.0;
float pim_vz45_c = 0.0;
float pim_vzlim43_c = 0.0;
float pim_vzlim44_c = 0.0;
float pim_vzlim45_c = 0.0;

float pim_vx53_c = 0.0;
float pim_vx54_c = 0.0;
float pim_vx55_c = 0.0;
float pim_vxlim53_c = 0.0;
float pim_vxlim54_c = 0.0;
float pim_vxlim55_c = 0.0;
float pim_vy53_c = 0.0;
float pim_vy54_c = 0.0;
float pim_vy55_c = 0.0;
float pim_vylim53_c = 0.0;
float pim_vylim54_c = 0.0;
float pim_vylim55_c = 0.0;
float pim_vz53_c = 0.0;
float pim_vz54_c = 0.0;
float pim_vz55_c = 0.0;
float pim_vzlim53_c = 0.0;
float pim_vzlim54_c = 0.0;
float pim_vzlim55_c = 0.0;

float pim_vx63_c = 0.0;
float pim_vx64_c = 0.0;
float pim_vx65_c = 0.0;
float pim_vxlim63_c = 0.0;
float pim_vxlim64_c = 0.0;
float pim_vxlim65_c = 0.0;
float pim_vy63_c = 0.0;
float pim_vy64_c = 0.0;
float pim_vy65_c = 0.0;
float pim_vylim63_c = 0.0;
float pim_vylim64_c = 0.0;
float pim_vylim65_c = 0.0;
float pim_vz63_c = 0.0;
float pim_vz64_c = 0.0;
float pim_vz65_c = 0.0;
float pim_vzlim63_c = 0.0;
float pim_vzlim64_c = 0.0;
float pim_vzlim65_c = 0.0;

float pim_vx1_c = 0.0;
float pim_vy1_c = 0.0;
float pim_vz1_c = 0.0;
float pim_vx2_c = 0.0;
float pim_vy2_c = 0.0;
float pim_vz2_c = 0.0;
float pim_vx3_c = 0.0;
float pim_vy3_c = 0.0;
float pim_vz3_c = 0.0;
float pim_vx4_c = 0.0;
float pim_vy4_c = 0.0;
float pim_vz4_c = 0.0;
float pim_vx5_c = 0.0;
float pim_vy5_c = 0.0;
float pim_vz5_c = 0.0;
float pim_vx6_c = 0.0;
float pim_vy6_c = 0.0;
float pim_vz6_c = 0.0;

float pim_lvy5c = -27.58642 - 1.3 - 3*(20.527298);   
float pim_uvy5c = -27.58642 - 1.3 + 3*(20.527298); 
float pim_lvx6c = -27.58642 - 1.3 - 4*(20.527298);   
float pim_uvx6c = -27.58642 - 1.3 + 4*(20.527298); 
float pim_lvy6c = -27.58642 - 1.3 - 5*(20.527298);   
float pim_uvy6c = -27.58642 - 1.3 + 5*(20.527298); 


float pim_ls13 = -22.65302 - 1.3 - 3*(20.544745);   
float pim_us13 = -22.65302 - 1.3 + 3*(20.544745); 
float pim_ls14 = -22.65302 - 1.3 - 4*(20.544745);   
float pim_us14 = -22.65302 - 1.3 + 4*(20.544745); 
float pim_ls15 = -22.65302 - 1.3 - 5*(20.544745);   
float pim_us15 = -22.65302 - 1.3 + 5*(20.544745); 


float pim_lc23 = -27.42106 - 1.3 - 3*(20.622596);   
float pim_uc23 = -27.42106 - 1.3 + 3*(20.622596); 
float pim_lc24 = -27.42106 - 1.3 - 4*(20.622596);   
float pim_uc24 = -27.42106 - 1.3 + 4*(20.622596); 
float pim_lc25 = -27.42106 - 1.3 - 5*(20.622596);   
float pim_uc25 = -27.42106 - 1.3 + 5*(20.622596); 
   

float pim_ls23 = -22.51326 - 1.3 - 3*(20.59944);   
float pim_us23 = -22.51326 - 1.3 + 3*(20.59944); 
float pim_ls24 = -22.51326 - 1.3 - 4*(20.59944);   
float pim_us24 = -22.51326 - 1.3 + 4*(20.59944); 
float pim_ls25 = -22.51326 - 1.3 - 5*(20.59944);   
float pim_us25 = -22.51326 - 1.3 + 5*(20.59944); 


float pim_lc33 = -27.41344 - 1.3 - 3*(20.572976);   
float pim_uc33 = -27.41344 - 1.3 + 3*(20.572976); 
float pim_lc34 = -27.41344 - 1.3 - 4*(20.572976);   
float pim_uc34 = -27.41344 - 1.3 + 4*(20.572976); 
float pim_lc35 = -27.41344 - 1.3 - 5*(20.572976);   
float pim_uc35 = -27.41344 - 1.3 + 5*(20.572976); 

 
float pim_ls33 = -22.48824 - 1.3 - 3*(20.54446);   
float pim_us33 = -22.48824 - 1.3 + 3*(20.54446); 
float pim_ls34 = -22.48824 - 1.3 - 4*(20.54446);   
float pim_us34 = -22.48824 - 1.3 + 4*(20.54446);
float pim_ls35 = -22.48824 - 1.3 - 5*(20.54446);   
float pim_us35 = -22.48824 - 1.3 + 5*(20.54446);


float pim_lvx1c3 = -27.60018 - 1.3 - 3*(20.571731);   
float pim_uvx1c3 = -27.60018 - 1.3 + 3*(20.571731);
float pim_lvx1c4 = -27.60018 - 1.3 - 4*(20.571731);   
float pim_uvx1c4 = -27.60018 - 1.3 + 4*(20.571731);
float pim_lvx1vy1c = -27.60018 - 1.3 - 5*(20.571731);   
float pim_uvx1vy1c = -27.60018 - 1.3 + 5*(20.571731);


float pim_ls43 = -22.66478 - 1.3 - 3*(20.559822);   
float pim_us43 = -22.66478 - 1.3 + 3*(20.559822);
float pim_ls44 = -22.66478 - 1.3 - 4*(20.559822);   
float pim_us44 = -22.66478 - 1.3 + 4*(20.559822);
float pim_ls45 = -22.66478 - 1.3 - 5*(20.559822);   
float pim_us45 = -22.66478 - 1.3 + 5*(20.559822);


float pim_lvy1c3 = -27.44553 - 1.3 - 3*(20.569488);   
float pim_uvy1c3 = -27.44553 - 1.3 + 3*(20.569488);
float pim_lvy1c4 = -27.44553 - 1.3 - 4*(20.569488);   
float pim_uvy1c4 = -27.44553 - 1.3 + 4*(20.569488);
float pim_lvy1c5 = -27.44553 - 1.3 - 5*(20.569488);   
float pim_uvy1c5 = -27.44553 - 1.3 + 5*(20.569488);


float pim_ls53 = -22.52733 - 1.3 - 3*(20.563511);   
float pim_us53 = -22.52733 - 1.3 + 3*(20.563511);
float pim_ls54 = -22.52733 - 1.3 - 4*(20.563511);   
float pim_us54 = -22.52733 - 1.3 + 4*(20.563511);
float pim_ls55 = -22.52733 - 1.3 - 5*(20.563511);   
float pim_us55 = -22.52733 - 1.3 + 5*(20.563511);


float pim_lvx2c3 = -27.57334 - 1.3 - 3*(20.515616);   
float pim_uvx2c3 = -27.57334 - 1.3 + 3*(20.515616);
float pim_lvx2c4 = -27.57334 - 1.3 - 4*(20.515616);   
float pim_uvx2c4 = -27.57334 - 1.3 + 4*(20.515616);
float pim_lvx2c5 = -27.57334 - 1.3 - 5*(20.515616);   
float pim_uvx2c5 = -27.57334 - 1.3 + 5*(20.515616);

float pim_ls63 = -22.64149  - 1.3 - 3*(20.508674);   
float pim_us63 = -22.64149  - 1.3 + 3*(20.508674);
float pim_ls64 = -22.64149  - 1.3 - 4*(20.508674);   
float pim_us64 = -22.64149  - 1.3 + 4*(20.508674);
float pim_ls65 = -22.64149  - 1.3 - 5*(20.508674);   
float pim_us65 = -22.64149  - 1.3 + 5*(20.508674);

float pim_lc3 = -27.50864 - 1.3 - 3*(20.5685);   
float pim_uc3 = -27.50864 - 1.3 + 3*(20.5685);
float pim_lvx1c = -27.50864 - 1.3 - 4*(20.5685);   
float pim_uvx1c = -27.50864 - 1.3 + 4*(20.5685);
float pim_lvy1c = -27.50864 - 1.3 - 5*(20.5685);   
float pim_uvy1c = -27.50864 - 1.3 + 5*(20.5685);

float pim_ls3 = -22.58252 - 1.3  - 3*(20.559604);   
float pim_us3 = -22.58252 - 1.3  + 3*(20.559604);
float pim_ls4 = -22.58252 - 1.3  - 4*(20.559604);   
float pim_us4 = -22.58252 - 1.3  + 4*(20.559604);
float pim_ls5 = -22.58252 - 1.3  - 5*(20.559604);   
float pim_us5 = -22.58252 - 1.3  + 5*(20.559604);

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
   TH1F *vx_hist = new TH1F(" vx", " vx", 200, -5, 5);
       vx_hist->SetLineColor(kBlack);
   TH1F *vy_hist = new TH1F(" vy", " vy", 200, -5, 5);
       vy_hist->SetLineColor(kBlack);
   //TH2F *xy_hist = new TH2F(" px .vs py", " px .vs py", 200, -1, 1, 200, -1, 1);
   TH1F *vz_hist = new TH1F(" vz", " vz", 200, -30, 20);
          vz_hist->SetLineColor(kBlack);
          
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
   ivz_hist->SetLineColor(kCyan+1);

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
	ivvx_hist->SetLineColor(kGreen);
   TH1F *ivvy_hist = new TH1F("sector 4 vy", "Sector 4 vy", 200, -5, 5);
	ivvy_hist->SetLineColor(kGreen);
   TH2F *ivxy_hist = new TH2F("sector 4 vx .vs vy", "sector 4 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *ivvz_hist = new TH1F("sector 4 vz", "Sector 4 vz", 200, -30, 20);
   	ivvz_hist->SetLineColor(kGreen);

   TH1F *vvx_hist = new TH1F("sector 5 vx", "Sector 5 vx", 200, -5, 5);
	vvx_hist->SetLineColor(kGreen);
   TH1F *vvy_hist = new TH1F("sector 5 vy", "Sector 5 vy", 200, -5, 5);
	vvy_hist->SetLineColor(kGreen);
   TH2F *vxy_hist = new TH2F("sector 5 vx .vs vy", "sector 5 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *vvz_hist = new TH1F("sector 5 vz", "Sector 5 vz", 200, -30, 20);
   	vvz_hist->SetLineColor(kGreen);

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
	pip_ivvx_hist->SetLineColor(kGreen);
   TH1F *pip_ivvy_hist = new TH1F("pip_sector 4 vy", "pip_Sector 4 vy", 200, -5, 5);
	pip_ivvy_hist->SetLineColor(kGreen);
   TH2F *pip_ivxy_hist = new TH2F("pip_sector 4 vx .vs vy", "pip_sector 4 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *pip_ivvz_hist = new TH1F("pip_sector 4 vz", "pip_Sector 4 vz", 200, -30, 20);
   	pip_ivvz_hist->SetLineColor(kGreen);

   TH1F *pip_vvx_hist = new TH1F("pip_sector 5 vx", "pip_Sector 5 vx", 200, -5, 5);
	pip_vvx_hist->SetLineColor(kGreen);
   TH1F *pip_vvy_hist = new TH1F("pip_sector 5 vy", "pip_Sector 5 vy", 200, -5, 5);
	pip_vvy_hist->SetLineColor(kGreen);
   TH2F *pip_vxy_hist = new TH2F("pip_sector 5 vx .vs vy", "pip_sector 5 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *pip_vvz_hist = new TH1F("pip_sector 5 vz", "pip_Sector 5 vz", 200, -30, 20);
   	pip_vvz_hist->SetLineColor(kGreen);

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
	 pim_ivvx_hist->SetLineColor(kGreen);
   TH1F *pim_ivvy_hist = new TH1F("pim_sector 4 vy", "pim_Sector 4 vy", 200, -5, 5);
	 pim_ivvy_hist->SetLineColor(kGreen);
   TH2F *pim_ivxy_hist = new TH2F("pim_sector 4 vx .vs vy", "pim_sector 4 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *pim_ivvz_hist = new TH1F("pim_sector 4 vz", "pim_Sector 4 vz", 200, -30, 20);
   	pim_ivvz_hist->SetLineColor(kGreen);

   TH1F *pim_vvx_hist = new TH1F("pim_sector 5 vx", "pim_Sector 5 vx", 200, -5, 5);
	 pim_vvx_hist->SetLineColor(kGreen);
   TH1F *pim_vvy_hist = new TH1F("pim_sector 5 vy", "pim_Sector 5 vy", 200, -5, 5);
	 pim_vvy_hist->SetLineColor(kGreen);
   TH2F *pim_vxy_hist = new TH2F("pim_sector 5 vx .vs vy", "pim_sector 5 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *pim_vvz_hist = new TH1F("pim_sector 5 vz", "pim_Sector 5 vz", 200, -30, 20);
   	pim_vvz_hist->SetLineColor(kGreen);

   TH1F *pim_vivx_hist = new TH1F("pim_sector 6 vx", "pim_Sector 6 vx", 200, -5, 5);
	 pim_vivx_hist->SetLineColor(kMagenta+1);
   TH1F *pim_vivy_hist = new TH1F("pim_sector 6 vy", "pim_Sector 6 vy", 200, -5, 5);
         pim_vivy_hist->SetLineColor(kMagenta+1);
   TH2F *pim_vixy_hist = new TH2F("pim_sector 6 vx .vs vy", "pim_sector 6 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1F *pim_vivz_hist = new TH1F("pim_sector 6 vz", "pim_Sector 6 vz", 200, -30, 20);
   	pim_vivz_hist->SetLineColor(kMagenta+1);
   
   TH1F *echi_hist = new TH1F("echi", "electron chi2 pid", 200, -40, 40);
   TH1F *pipchi_hist = new TH1F("pipchi", "pip chi2 pid", 200, -40, 40);
   TH1F *pimchi_hist = new TH1F("pimchi", "pim chi2 pid", 200, -40, 40);

   TH2F *tz_hist = new TH2F("vz vs phi", "vz vs phi", 200, -180, 180, 200, -30, 20);

   

   TH1F *rm_no = new TH1F("Without Cuts", "Without Cuts", 200, 0.2, 1.4);
      //rm_no->SetStats(0);
   TH1F *rm_w = new TH1F("With W Cut", "With W Cut ", 200, 0.2, 1.4);
      //rm_w->SetStats(0);
   TH1F *rm_wt = new TH1F("With W Cut & -t Cut", "With W Cut & -t Cut ", 200, 0.2, 1.4);
      //rm_wt->SetStats(0);
   TH1F *rm_wtz = new TH1F("With W Cut & -t Cut & z Cut", " With W Cut & -t Cut & z Cut", 200, 0.2, 1.4);
      //rm_wtz->SetStats(0);

   //auto sdaq_ev = new TGraph();
   //   sdaq_ev->SetLineColor(kRed);
   
   auto edaq_ev = new TGraph();
      edaq_ev->SetMarkerColor(kGreen);
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
      ivetot_p->SetMarkerColor(kGreen);
      ivetot_p->SetMarkerStyle(kFullCircle);
      ivetot_p->SetLineWidth(0);

   auto vetot_p = new TGraph();
      vetot_p->SetMarkerColor(kGreen);
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
      hel_fcup->SetMarkerColor(kGreen);
      hel_fcup->SetMarkerStyle(kFullCircle);
      hel_fcup->SetLineWidth(0);

   auto hel_fcupgated = new TGraph();
      hel_fcupgated->SetMarkerColor(kBlack);
      hel_fcupgated->SetMarkerStyle(kFullCircle);
      hel_fcupgated->SetLineWidth(0);

   auto hel_clock = new TGraph();
      hel_clock->SetMarkerColor(kGreen);
      hel_clock->SetMarkerStyle(kFullCircle);
      hel_clock->SetLineWidth(0);

   auto hel_clockgated = new TGraph();
      hel_clockgated->SetMarkerColor(kBlack);
      hel_clockgated->SetMarkerStyle(kFullCircle);
      hel_clockgated->SetLineWidth(0);

   TH1F *ietot_hist = new TH1F("sector 1 E_{tot}", "Sector 1 E_{tot}", 200, 0, 0.4);
       ietot_hist->SetFillColor(kCyan+1);
   TH1F *ip_hist = new TH1F("Q2 (2->2.5)", "Q2 (2->2.5)", 200, 0, 2);
       ip_hist->SetFillColor(kGreen);
   TH1F *p_hist = new TH1F("Q2 (1->2)", "Q2 (1->2)", 200, 0, 2);
       p_hist->SetFillColor(kGreen);

   TH1F *iietot_hist = new TH1F("sector 2 E_{tot}", "Sector 2 E_{tot}", 200, 0, 0.4);
	iietot_hist->SetFillColor(kOrange+1);
   TH1F *iip_hist = new TH1F("Q2 (2.5->3)", "Q2 (2.5->3)", 200, 0, 2);
	iip_hist->SetFillColor(kGreen);
   

   TH1F *iiietot_hist = new TH1F("sector 3 E_{tot}", "Sector 3 E_{tot}", 200, 0, 0.4);
	iiietot_hist->SetFillColor(kYellow+1);
   TH1F *iiip_hist = new TH1F("Q2 (3->3.5)", "Q2 (3->3.5)", 200, 0, 2);
	iiip_hist->SetFillColor(kGreen);
   

   TH1F *ivetot_hist = new TH1F("sector 4 E_{tot}", "Sector 4 E_{tot}", 200, 0, 0.4);
	ivetot_hist->SetFillColor(kMagenta);
   TH1F *ivp_hist = new TH1F("Q2 (3.5->4)", "Q2 (3.5->4)", 200, 0, 2);
	ivp_hist->SetFillColor(kGreen);
   

   TH1F *vetot_hist = new TH1F("sector 5 E_{tot}", "Sector 5 E_{tot}", 200, 0, 0.4);
	vetot_hist->SetFillColor(kGreen);
   TH1F *vp_hist = new TH1F("Q2 (4->6)", "Q2 (4->6)", 200, 0, 2);
	vp_hist->SetFillColor(kGreen);
   

   TH1F *vietot_hist = new TH1F("sector 6 E_{tot}", "Sector 6 E_{tot}", 200, 0, 0.4);
	vietot_hist->SetFillColor(kMagenta+1);
   TH1F *vip_hist = new TH1F("sector 6 P", "Q2 (4.5->6)", 200, 0, 2);
        vip_hist->SetFillColor(kGreen);

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
   int num = 0;
//47, 48, 50, 52, 53, 54, 55, 56, 58, 59, 60
//670, 671, 672, 673, 675, 676, 677, 678, 679
//9110, 9111, 9112, 9113, 9114, 9115, 9117, 9118, 9119

   TFile f("cu_pass0v7_nt_21.root", "update"); //cusn_pass0v7_63_test.root

    std::ifstream v3("cusn_calv9_skim.list");  //   rge_align.list    cusn_calv9_skim.list   cxc_pass0v7_skim.list   ld_pass0v7_skim.list
         std::string line;
   while(std::getline(v3, line)){
	std::stringstream ss(line);
        std::string hip;
        std::getline(ss, hip, ',');
        v3_a.push_back(hip);
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
      for(int p = 33; p < v3_a.size(); p++){ // for(int p = 0; p < v3_a.size(); p++){  cusn(34 37)33  cxc(18, 20)13  LD2(18,)
        std::string inFile;
            inFile =  v3_a[p];
         std::cout << v3_a[p] << std::endl;
       
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
   
   // in GeV/c^2
   float m_pro = 0.938; //in Gev/c^2
   //float m_ele = 0.000511; // in GeV/c^2
   float m_pi = 0.1396;
   float m_rho = 0.770;

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

   
      char* inputFile_char = new char[inFile.length() + 1];
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
      //hipo::bank res(factory.getSchema("REC::Response"));

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
         
       // particle.show();
        // events.show();
        //hel.show();
        //cher.show();
        //cal.show();
        //track.show();
        
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
         
         float e_total = 0.0;
         float ppcal = 0.0;
         float pcal = 0.0;
         float nphe = 0.0;
         
         double ez_lcut = -10.502317;  //for Copper
         double ez_ucut = -6.452905;
         //double ez_lcut = -6.008624;  //for Tin
         //double ez_ucut = 5;
         //double ez_lcut = -20;  //for Liquid Deuterium
         //double ez_ucut = 5;
         //double ez_lcut = -10.531918;  //for Carbon
         //double ez_ucut = 5;
         //double ez_lcut = -10000000000;  //for no cuts
         //double ez_ucut =  10000000000;
         
         double pipz_lcut = -9.760247;  //for Copper
         double pipz_ucut = -5.167685;
         //double pipz_lcut = -4.641313;  //for Tin
         //double pipz_ucut = 5;
         //double pipz_lcut = -20;  //for Liquid Deuterium
         //double pipz_ucut = 5;
         //double pipz_lcut = -9.658612;  //for Carbon
         //double pipz_ucut = 5;
         //double pipz_lcut = -10000000000;  //for no cuts
         //double pipz_ucut =  10000000000;
         
         double pimz_lcut = -10.955292;  //for Copper
         double pimz_ucut = -5.8338;
         //double pimz_lcut = -5.462858;  //for Tin
         //double pimz_ucut = 5;
         //double pimz_lcut = -20;  //for Liquid Deuterium
         //double pimz_ucut = 5;
         //double pimz_lcut = -10.350112;  //for Carbon
         //double pimz_ucut = 5;
         //double pimz_lcut = -10000000000;  //for no cuts
         //double pimz_ucut =  10000000000;

        
       
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
			               
                        if (abs(temp_s)/2000 == 1 && temp_s < 0 && temp_vz >= ez_lcut && temp_vz <= ez_ucut){
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
            if(chi2 > -10 && chi2 < 10 && temp_vz >= pipz_lcut && temp_vz <= pipz_ucut){	
              	a_pip_pid[pid_pip_count] = 211;
              	a_pip_px[pid_pip_count] = particle.getFloat("px",row);
           	   a_pip_py[pid_pip_count] = particle.getFloat("py",row);
              	a_pip_pz[pid_pip_count] = particle.getFloat("pz",row);
               a_pip_vx[pid_pip_count] = particle.getFloat("vx",row);
              	a_pip_vy[pid_pip_count] = particle.getFloat("vy",row);
           	   a_pip_vz[pid_pip_count] = particle.getFloat("vz",row);
               a_pip_energy[pid_pip_count] = sqrt(pow(a_pip_px[pid_pip_count],2) + pow(a_pip_py[pid_pip_count],2) + pow(a_pip_pz[pid_pip_count],2) + pow(m_pi,2));
               
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
            if(chi2 > -10 && chi2 < 10 && temp_vz >= pimz_lcut && temp_vz <= pimz_ucut){	
              	a_pim_pid[pid_pim_count] = -211;
              	a_pim_px[pid_pim_count] = particle.getFloat("px",row);
              	a_pim_py[pid_pim_count] = particle.getFloat("py",row);
           	   a_pim_pz[pid_pim_count] = particle.getFloat("pz",row);
               a_pim_vx[pid_pim_count] = particle.getFloat("vx",row);
              	a_pim_vy[pid_pim_count] = particle.getFloat("vy",row);
              	a_pim_vz[pid_pim_count] = particle.getFloat("vz",row);
               a_pim_energy[pid_pim_count] = sqrt(pow(a_pim_px[pid_pim_count],2) + pow(a_pim_py[pid_pim_count],2) + pow(a_pim_pz[pid_pim_count],2) + pow(m_pi,2));
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
            
            for(int ec = 0; ec < ne; ec++){
                for(int row = 0; row < nrows; row++){
                    if(track.getInt("pindex",row) == a_pin[ec]){
                        a_sect[ec] = track.getInt("sector",row);
                       // std::cout << a_sect[ec] << std::endl;
                     //  ac_px[ec] = cher.getFloat("px", row);
                    }
                }
            }
            
            for(int pipc = 0; pipc < npip; pipc++){
                for(int row = 0; row < nrows; row++){
                    if(track.getInt("pindex",row) == a_pip_pin[pipc]){
                        a_pip_sect[pipc] = track.getInt("sector",row);
                        //std::cout << a_pip_sect[pipc] << std::endl;
                    }
                }
            }
            
            for(int pimc = 0; pimc < npim; pimc++){
                for(int row = 0; row < nrows; row++){
                    if(track.getInt("pindex",row) == a_pim_pin[pimc]){
                        a_pim_sect[pimc] = track.getInt("sector",row);
                    }
                }
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
               
               
               e_hist->Fill(pid);
               xy_hist->Fill(px,py);
               te_lv.SetPxPyPzE(px,py,pz,energy);
			count_3++;		
                                   sec_dec =  te_lv.Phi() * 180/M_PI;
                                   e_theta = te_lv.Theta() * 180/M_PI;
                                   e_pht_hist->Fill(sec_dec, e_theta);
                                
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
				if(sec_dec > -30 && sec_dec <= 30) {
				if( vx >= lxlim && vx <= uxlim )  {  
				ivx_hist->Fill(vx);
				vx1_c++;
				}
				if(  vy >= lylim && vy <= uylim)  {  
   				ivy_hist->Fill(vy);
   				vy1_c++;
   				}
   				if( vx >= lxlim && vx <= uxlim && vy >= lylim && vy <= uylim && vz >= zl_b && vz <= zu_b)  {  
   				ivz_hist->Fill(vz);
   				vz1_c++;
				}
				
                        } else if(sec_dec > 30 && sec_dec <= 90) {
                                if( vx >= lxlim && vx <= uxlim )  {  
				iivx_hist->Fill(vx);
				vx2_c++;
				}
				if(  vy >= lylim && vy <= uylim)  {  
   				iivy_hist->Fill(vy);
   				vy2_c++;
   				}
				iixy_hist->Fill(px, py);
				if( vx >= lxlim && vx <= uxlim && vy >= lylim && vy <= uylim && vz >= zl_b && vz <= zu_b)  {  
				iivz_hist->Fill(vz);
				vz2_c++;
				}
				
                        } else if(sec_dec > 90 && sec_dec <= 150) {
                                if( vx >= lxlim && vx <= uxlim )  {  
				iiivx_hist->Fill(vx);
				vx3_c++;
				}
				if(  vy >= lylim && vy <= uylim)  {  
   				iiivy_hist->Fill(vy);
   				vy3_c++;
   				}
				iiixy_hist->Fill(px, py);
				if( vx >= lxlim && vx <= uxlim && vy >= lylim && vy <= uylim && vz >= zl_b && vz <= zu_b)  {  
				iiivz_hist->Fill(vz);
				vz3_c++;
				}
				
				
                        } else if(sec_dec > 150 || sec_dec <= -150) {
                                if( vx >= lxlim && vx <= uxlim )  {  
				ivvx_hist->Fill(vx);
				vx4_c++;
				//vx5_c++;
				}
				if(  vy >= lylim && vy <= uylim)  {  
   				ivvy_hist->Fill(vy);
   				vy4_c++;
   				}
				ivxy_hist->Fill(px, py);
				if( vx >= lxlim && vx <= uxlim && vy >= lylim && vy <= uylim && vz >= zl_b && vz <= zu_b)  {  
				ivvz_hist->Fill(vz);
				vz4_c++;
				}
				
                        } else if(sec_dec > -150 && sec_dec <= -90) {
                                if( vx >= lxlim && vx <= uxlim )  {  
				vvx_hist->Fill(vx);
				vx5_c++;
				}
				if( vy >= lylim && vy <= uylim)  {  
   				vvy_hist->Fill(vy);
   				vy5_c++;
   				}
				vxy_hist->Fill(px, py);
				if( vx >= lxlim && vx <= uxlim && vy >= lylim && vy <= uylim && vz >= zl_b && vz <= zu_b)  {  
				vvz_hist->Fill(vz);
				vz5_c++;
				}
				
				
                        } else if(sec_dec > -90 && sec_dec <= -30) {
                                if( vx >= lxlim && vx <= uxlim )  {  
				vivx_hist->Fill(vx);
				vx6_c++;
				}
				if( vy >= lylim && vy <= uylim)  {  
   				vivy_hist->Fill(vy);
   				vy6_c++;
   				}
				vixy_hist->Fill(px, py);
				if( vx >= lxlim && vx <= uxlim && vy >= lylim && vy <= uylim && vz >= zl_b && vz <= zu_b)  {  
				vivz_hist->Fill(vz);
				vz6_c++;
				}
				
				}
                       }
            
            
         
            
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
	       
	       
	       
                              pip_te_lv.SetPxPyPzE(pip_px,pip_py,pip_pz,pip_energy);
			//count_3++;		
                                   pip_sec_dec =  pip_te_lv.Phi() * 180/M_PI;
                                   pip_theta = pip_te_lv.Theta() * 180/M_PI;
                                   pip_pht_hist->Fill(pip_sec_dec, pip_theta);
                                
                 float pip_mul = 110000.0; // 3.0 4.0 5.0
                 float pip_lxlim = 0.00656862 - pip_mul*(0.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float pip_uxlim = 0.00656862 + pip_mul*(0.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float pip_lylim = -0.1793 - pip_mul*(0.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float pip_uylim = -0.1793 + pip_mul*(0.46594);  //1.5056; // 0.83372 1.16966 1.5056       
                   //    if( pip_vx >= -6 && pip_vx <= 6 && pip_vy >= -6 && pip_vy <= 6 )  {                
                          //tz_hist->Fill(sec_dec, vz);
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
                                
                 float pim_mul = 1100000.0; // 3.0 4.0 5.0
                 float pim_lxlim = 0.00656862 - pim_mul*(0.4567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float pim_uxlim = 0.00656862 + pim_mul*(0.4567);  //1.5741; // 0.9525 1.2633 1.5741
                 float pim_lylim = -0.1793 - pim_mul*(0.46594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float pim_uylim = -0.1793 + pim_mul*(0.46594);  //1.5056; // 0.83372 1.16966 1.5056       
                     //  if( pim_vx >= -6 && pim_vx <= 6 && pim_vy >= -6 && pim_vy <= 6)  {                
                          //tz_hist->Fill(sec_dec, vz);
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
            
            


            for(int i = 0; i < npip; i++){                               //This double for loop will compare each pip to each pim in an event
               for(int j = 0; j < npim; j++){                            //This is primarly done using TLorentzsVector (the commented parts are my original attempt to calculate)
               
                    if(a_vz[0] >= -100000000000000 && a_vz[0] <= 100000000000000000 ){
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

            event_counter++;
         
            if(rho_counter > 0){
               for(int i = 0; i < rho_counter; i++){
               
               el_lv.SetPxPyPzE(a_px[0],a_py[0],a_pz[0],a_energy[0]);

               vpho_lv = e_lv - el_lv;

               Q2 = - vpho_lv.Mag2(); //((e_lv - el_lv) * (e_lv - el_lv));
               w_lv = vpho_lv + p_lv; 
               w = w_lv.Mag();

               w_hist->Fill(w);

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
                  
                  rm_no->Fill(rho_mass);
                  //std::cout << rho_lv.Mag() << std::endl;
                  //std::cout << w << " | " << Q2 << std::endl;
                  
                  r_count++;

                  if(w > 2 && lc <= 1.0 ){
               rw_count++;
              rm_w->Fill(rho_mass); 
                   
                     if(-t > 0.1 && -t < 0.5){
                       rm_wt->Fill(rho_mass); 
                       rwt_count++;
                        if(Zh > 0.9 && Zh < 1.0){
                          rm_wtz->Fill(rho_mass);
                         rwtz_count++; 
                           if(Q2 >= 1 && Q2 <2){
         			r1_count++;
                              		p_hist->Fill(rho_mass);
				
                            }else if(Q2 >= 2 && Q2 <2.5){
                              ip_hist->Fill(rho_mass);
                              r2_count++;
                            }else if(Q2 >= 2.5 && Q2 <3){
    			r3_count++;
                              		iip_hist->Fill(rho_mass);
				
                           }else if(Q2 >= 3 && Q2 <3.5){
                              iiip_hist->Fill(rho_mass);
                              r4_count++;
                           }else if(Q2 >= 3.5 && Q2 <4){
                              ivp_hist->Fill(rho_mass);
                              r5_count++;
                           }else if(Q2 >= 4 && Q2 <6){
                              vp_hist->Fill(rho_mass);
                              r6_count++;
                        
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
     
      std::cout << "p = " << ((p + 1.0) / v3_a.size()) * 100 << "%" <<std::endl;
   }
//}
   std::cout << "e with pid = " << count_1 << std::endl;
   std::cout << "e with status, vz, chi2pid = " << count_2 << std::endl;
   std::cout << "e wit all = " << count_3 << std::endl;
   std::cout << "e pip pim events = " << e_pip_pim_num << std::endl;
   
   //std::cout << "e count " << e_count << std::endl;
   
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vx1 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vx13_c/vx1_c) * 100 << " | " << (vxlim13_c/vx1_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vx14_c/vx1_c) * 100 << " | " << (vxlim14_c/vx1_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vx15_c/vx1_c) * 100 << " | " << (vxlim15_c/vx1_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vx2 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vx23_c/vx2_c) * 100 << " | " << (vxlim23_c/vx2_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vx24_c/vx2_c) * 100 << " | " << (vxlim24_c/vx2_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vx25_c/vx2_c) * 100 << " | " << (vxlim25_c/vx2_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vx3 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vx33_c/vx3_c) * 100 << " | " << (vxlim33_c/vx3_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vx34_c/vx3_c) * 100 << " | " << (vxlim34_c/vx3_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vx35_c/vx3_c) * 100 << " | " << (vxlim35_c/vx3_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vx4 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vx43_c/vx4_c) * 100 << " | " << (vxlim43_c/vx4_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vx44_c/vx4_c) * 100 << " | " << (vxlim44_c/vx4_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vx45_c/vx4_c) * 100 << " | " << (vxlim45_c/vx4_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vx5 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vx53_c/vx5_c) * 100 << " | " << (vxlim53_c/vx5_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vx54_c/vx5_c) * 100 << " | " << (vxlim54_c/vx5_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vx55_c/vx5_c) * 100 << " | " << (vxlim55_c/vx5_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vx6 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vx63_c/vx6_c) * 100 << " | " << (vxlim63_c/vx6_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vx64_c/vx6_c) * 100 << " | " << (vxlim64_c/vx6_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vx65_c/vx6_c) * 100 << " | " << (vxlim65_c/vx6_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vy1 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vy13_c/vy1_c) * 100 << " | " << (vylim13_c/vy1_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vy14_c/vy1_c) * 100 << " | " << (vylim14_c/vy1_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vy15_c/vy1_c) * 100 << " | " << (vylim15_c/vy1_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vy2 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vy23_c/vy2_c) * 100 << " | " << (vylim23_c/vy2_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vy24_c/vy2_c) * 100 << " | " << (vylim24_c/vy2_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vy25_c/vy2_c) * 100 << " | " << (vylim25_c/vy2_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vy3 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vy33_c/vy3_c) * 100 << " | " << (vylim33_c/vy3_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vy34_c/vy3_c) * 100 << " | " << (vylim34_c/vy3_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vy35_c/vy3_c) * 100 << " | " << (vylim35_c/vy3_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vy4 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vy43_c/vy4_c) * 100 << " | " << (vylim43_c/vy4_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vy44_c/vy4_c) * 100 << " | " << (vylim44_c/vy4_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vy45_c/vy4_c) * 100 << " | " << (vylim45_c/vy4_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vy5 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vy53_c/vy5_c) * 100 << " | " << (vylim53_c/vy5_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vy54_c/vy5_c) * 100 << " | " << (vylim54_c/vy5_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vy55_c/vy5_c) * 100 << " | " << (vylim55_c/vy5_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vy6 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vy63_c/vy6_c) * 100 << " | " << (vylim63_c/vy6_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vy64_c/vy6_c) * 100 << " | " << (vylim64_c/vy6_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vy65_c/vy6_c) * 100 << " | " << (vylim65_c/vy6_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vz1 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vz13_c/vz1_c) * 100 << " | " << (vzlim13_c/vz1_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vz14_c/vz1_c) * 100 << " | " << (vzlim14_c/vz1_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vz15_c/vz1_c) * 100 << " | " << (vzlim15_c/vz1_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vz2 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vz23_c/vz2_c) * 100 << " | " << (vzlim23_c/vz2_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vz24_c/vz2_c) * 100 << " | " << (vzlim24_c/vz2_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vz25_c/vz2_c) * 100 << " | " << (vzlim25_c/vz2_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vz3 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vz33_c/vz3_c) * 100 << " | " << (vzlim33_c/vz3_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vz34_c/vz3_c) * 100 << " | " << (vzlim34_c/vz3_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vz35_c/vz3_c) * 100 << " | " << (vzlim35_c/vz3_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vz4 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vz43_c/vz4_c) * 100 << " | " << (vzlim43_c/vz4_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vz44_c/vz4_c) * 100 << " | " << (vzlim44_c/vz4_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vz45_c/vz4_c) * 100 << " | " << (vzlim45_c/vz4_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vz5 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vz53_c/vz5_c) * 100 << " | " << (vzlim53_c/vz5_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vz54_c/vz5_c) * 100 << " | " << (vzlim54_c/vz5_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vz55_c/vz5_c) * 100 << " | " << (vzlim55_c/vz5_c) * 100 << std::endl;
   std::cout << "---------------------------------" << std::endl;
   std::cout << "vz6 : sector fit % | total fit %" << std::endl;
   std::cout << "3 sigma : These don't mean shit yet " << (vz63_c/vz6_c) * 100 << " | " << (vzlim63_c/vz6_c) * 100 << std::endl;
   std::cout << "4 sigma : These don't mean shit yet " << (vz64_c/vz6_c) * 100 << " | " << (vzlim64_c/vz6_c) * 100 << std::endl;
   std::cout << "5 sigma : These don't mean shit yet " << (vz65_c/vz6_c) * 100 << " | " << (vzlim65_c/vz6_c) * 100 << std::endl;
   
   
   
   
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
   
   
   tz_hist->GetXaxis()->SetTitle("#phi_{e} [deg]");
   tz_hist->GetYaxis()->SetTitle("vz [cm]");
   tz_hist->Draw();
   tz_hist->Write();
   
   echi_hist->GetXaxis()->SetTitle("e chi2");
   echi_hist->Draw();
   echi_hist->Write();
   
   pipchi_hist->GetXaxis()->SetTitle("pip chi2");
   pipchi_hist->Draw();
   pipchi_hist->Write();
   
   pimchi_hist->GetXaxis()->SetTitle("pim chi2");
   pimchi_hist->Draw();
   pimchi_hist->Write();

   q_hist->GetXaxis()->SetTitle("Q2");
   q_hist->Draw();
   q_hist->Write();
   
   qz_hist->GetXaxis()->SetTitle("Q2 [GeV^{2}]");
   qz_hist->GetYaxis()->SetTitle("lc [fm]");
   qz_hist->Draw();
   qz_hist->Write();
   
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

   //wonphe_hist->GetXaxis()->SetTitle("# of Photoelectrons");
   //wonphe_hist->Draw();
   //wonphe_hist->Write();

   //wnphe_hist->GetXaxis()->SetTitle("# of Photoelectrons");
   //wnphe_hist->Draw();
   //wnphe_hist->Write();

   //wopcal_hist->GetXaxis()->SetTitle("Total Energy [GeV]/Momentum");
   //wopcal_hist->Draw();
   //wopcal_hist->Write();

   //wpcal_hist->GetXaxis()->SetTitle("PCAL Energy [GeV]");
   //wpcal_hist->Draw();
   //wpcal_hist->Write();

   //wocut_hist->GetXaxis()->SetTitle("# of PhotoElectron");
   //wocut_hist->GetYaxis()->SetTitle("PCAL [GeV]");
   //wocut_hist->Draw();
   //wocut_hist->Write();
/*
   wcut_hist->GetXaxis()->SetTitle("P [GeV]");
   wcut_hist->GetYaxis()->SetTitle("E_{tot} / P");
   wcut_hist->Draw();
   wcut_hist->Write();
*/
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

   //sdaq_ev->SetTitle("Event number v scaler lifetime;Event number;DAQ lifetime[s]");
   //sdaq_ev->Draw();
   //sdaq_ev->Write();

   //edaq_ev->SetTitle("Nuclear Transparency;Q2 [GeV];# of rho");
   //edaq_ev->Draw("AP");
   //edaq_ev->Write();
/*
   ivx_hist->Draw();
   c1->Modified(); c1->Update();
   iivx_hist->Draw("SAME");
   iiivx_hist->Draw("SAME");
   ivvx_hist->Draw("SAME");
   vvx_hist->Draw("SAME");
   vivx_hist->Draw("SAME");
   c1->SaveAs("zz_19125_vx.pdf");

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
   ivz_hist->Draw();
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
   vz_hist->Fit(vz_fit, "R");
   vz_hist->Fit(vz2_fit, "R+");
   vz_hist->Draw();
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
   vzz_hist->Fit(vzz_fit, "R");
   vzz_hist->Fit(vzz2_fit, "R+");
   vzz_hist->Draw();
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
      TF1 *pip_vz_fit = new TF1("pip_vz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
   pip_vz_fit->SetParameters(0.01, -8, 1, 0.01, -3, 1);
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
      TF1 *pim_vz_fit = new TF1("pim_vz_fit", "gaus(0)+ gaus(3) + pol1(6)",  - 10, -0.5);
      pim_vz_fit->SetParameters(0.01, -8, 1, 0.01, -3, 1);
      pim_vz_hist->Fit(pim_vz_fit, "R");
      pim_vz_hist->Draw();
      pim_vz_hist->Write();
      
      
   TCanvas *vxt = new TCanvas("vxt","vx",0,0,640,480);
      TLine *llxlim5=new TLine(lxlim5,vxt->GetUymin(),lxlim5,vxt->GetUymax());
      llxlim5->SetLineColor(kCyan-1);
      TLine *luxlim5=new TLine(uxlim5,vxt->GetUymin(),uxlim5,vxt->GetUymax());
      luxlim5->SetLineColor(kCyan-1);
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
      TLine *llylim5=new TLine(lylim5,vyt->GetUymin(),lylim5,vyt->GetUymax());
      llylim5->SetLineColor(kCyan-1);
      TLine *luylim5=new TLine(uylim5,vyt->GetUymin(),uylim5,vyt->GetUymax());
      luylim5->SetLineColor(kCyan-1);
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
   
   
   TCanvas *vzt = new TCanvas("vzt","vz",0,0,640,480);
      TLine *llclim5=new TLine(lvy1c,vzt->GetUymin(),lvy1c,vzt->GetUymax());
      llclim5->SetLineColor(kCyan-1);
      llclim5->SetLineWidth(3);
      TLine *luclim5=new TLine(uvy1c,vzt->GetUymin(),uvy1c,vzt->GetUymax());
      luclim5->SetLineColor(kCyan-1);
      luclim5->SetLineWidth(3);
      TLine *llslim5=new TLine(ls5,vzt->GetUymin(),ls5,vzt->GetUymax());
      llslim5->SetLineColor(kMagenta-7);
      llslim5->SetLineWidth(3);
      llslim5->SetLineStyle(10);
      TLine *luslim5=new TLine(us5,vzt->GetUymin(),us5,vzt->GetUymax());
      luslim5->SetLineColor(kMagenta-7);
      luslim5->SetLineWidth(3);
      luslim5->SetLineStyle(10);
   vz_hist->Draw("HIST");
   ivz_hist->Draw("HIST Same");
   iivz_hist->Draw("HIST Same");
   iiivz_hist->Draw("HIST Same");
   ivvz_hist->Draw("HIST Same");
   vvz_hist->Draw("HIST Same");
   vivz_hist->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vzt->Write();
   
   
   TCanvas *pip_vxt = new TCanvas("pip_vxt","pip_vx",0,0,640,480);
      TLine *pip_llxlim5=new TLine(pip_lxlim5,pip_vxt->GetUymin(),pip_lxlim5,pip_vxt->GetUymax());
      pip_llxlim5->SetLineColor(kCyan-1);
      TLine *pip_luxlim5=new TLine(pip_uxlim5,pip_vxt->GetUymin(),pip_uxlim5,pip_vxt->GetUymax());
      pip_luxlim5->SetLineColor(kCyan-1);
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
      TLine *pip_llylim5=new TLine(pip_lylim5,pip_vyt->GetUymin(),pip_lylim5,pip_vyt->GetUymax());
      pip_llylim5->SetLineColor(kCyan-1);
      TLine *pip_luylim5=new TLine(pip_uylim5,pip_vyt->GetUymin(),pip_uylim5,pip_vyt->GetUymax());
      pip_luylim5->SetLineColor(kCyan-1);
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
      TLine *pim_llxlim5=new TLine(pim_lxlim5,pim_vxt->GetUymin(),pim_lxlim5,pim_vxt->GetUymax());
      pim_llxlim5->SetLineColor(kCyan-1);
      TLine *pim_luxlim5=new TLine(pim_uxlim5,pim_vxt->GetUymin(),pim_uxlim5,pim_vxt->GetUymax());
      pim_luxlim5->SetLineColor(kCyan-1);
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
      TLine *pim_llylim5=new TLine(pim_lylim5,pim_vyt->GetUymin(),pim_lylim5,pim_vyt->GetUymax());
      pim_llylim5->SetLineColor(kCyan-1);
      TLine *pim_luylim5=new TLine(pim_uylim5,pim_vyt->GetUymin(),pim_uylim5,pim_vyt->GetUymax());
      pim_luylim5->SetLineColor(kCyan-1);
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
   
   
   
/*
TCanvas *vx1c = new TCanvas("vx1c","vx1",0,0,640,480);
      TLine *llx13=new TLine(lx13,vx1c->GetUymin(),lx13,ivx_hist->GetMaximum());
      llx13->SetLineColor(kRed);
      TLine *lux13=new TLine(ux13,vx1c->GetUymin(),ux13,ivx_hist->GetMaximum());
      lux13->SetLineColor(kRed);
      TLine *llx14=new TLine(lx14,vx1c->GetUymin(),lx14,ivx_hist->GetMaximum());
      llx14->SetLineColor(kBlack);
      TLine *lux14=new TLine(ux14,vx1c->GetUymin(),ux14,ivx_hist->GetMaximum());
      lux14->SetLineColor(kBlack);
      TLine *llx15=new TLine(lx15,vx1c->GetUymin(),lx15,ivx_hist->GetMaximum());
      llx15->SetLineColor(kRed+2);
      TLine *lux15=new TLine(ux15,vx1c->GetUymin(),ux15,ivx_hist->GetMaximum());
      lux15->SetLineColor(kRed+2);
      TLine *llxlim3=new TLine(lxlim3,vx1c->GetUymin(),lxlim3,ivx_hist->GetMaximum());
      llxlim3->SetLineColor(kRed-7);
      TLine *luxlim3=new TLine(uxlim3,vx1c->GetUymin(),uxlim3,ivx_hist->GetMaximum());
      luxlim3->SetLineColor(kRed-7);
      TLine *llxlim4=new TLine(lxlim4,vx1c->GetUymin(),lxlim4,ivx_hist->GetMaximum());
      llxlim4->SetLineColor(kGray+1);
      TLine *luxlim4=new TLine(uxlim4,vx1c->GetUymin(),uxlim4,ivx_hist->GetMaximum());
      luxlim4->SetLineColor(kGray+1);
      TLine *llxlim5=new TLine(lxlim5,vx1c->GetUymin(),lxlim5,ivx_hist->GetMaximum());
      llxlim5->SetLineColor(kCyan-1);
      TLine *luxlim5=new TLine(uxlim5,vx1c->GetUymin(),uxlim5,ivx_hist->GetMaximum());
      luxlim5->SetLineColor(kCyan-1);
   ivx_hist->Draw("HIST");
   llx13->Draw("HIST Same");
   lux13->Draw("HIST Same");
   llx14->Draw("HIST Same");
   lux14->Draw("HIST Same");
   llx15->Draw("HIST Same");
   lux15->Draw("HIST Same");
   llxlim3->Draw("HIST Same");
   luxlim3->Draw("HIST Same");
   llxlim4->Draw("HIST Same");
   luxlim4->Draw("HIST Same");
   llxlim5->Draw("HIST Same");
   luxlim5->Draw("HIST Same");
   vx1c->Write();


   TCanvas *vy1c = new TCanvas("vy1c","vy1",0,0,640,480);
      TLine *lly13=new TLine(ly13,vy1c->GetUymin(),ly13,ivy_hist->GetMaximum());
      lly13->SetLineColor(kRed);
      TLine *luy13=new TLine(uy13,vy1c->GetUymin(),uy13,ivy_hist->GetMaximum());
      luy13->SetLineColor(kRed);
      TLine *lly14=new TLine(ly14,vy1c->GetUymin(),ly14,ivy_hist->GetMaximum());
      lly14->SetLineColor(kBlack);
      TLine *luy14=new TLine(uy14,vy1c->GetUymin(),uy14,ivy_hist->GetMaximum());
      luy14->SetLineColor(kBlack);
      TLine *lly15=new TLine(ly15,vy1c->GetUymin(),ly15,ivy_hist->GetMaximum());
      lly15->SetLineColor(kRed+2);
      TLine *luy15=new TLine(uy15,vy1c->GetUymin(),uy15,ivy_hist->GetMaximum());
      luy15->SetLineColor(kRed+2);
      TLine *llylim3=new TLine(lylim3,vy1c->GetUymin(),lylim3,ivy_hist->GetMaximum());
      llylim3->SetLineColor(kRed-7);
      TLine *luylim3=new TLine(uylim3,vy1c->GetUymin(),uylim3,ivy_hist->GetMaximum());
      luylim3->SetLineColor(kRed-7);
      TLine *llylim4=new TLine(lylim4,vy1c->GetUymin(),lylim4,ivy_hist->GetMaximum());
      llylim4->SetLineColor(kGray+1);
      TLine *luylim4=new TLine(uylim4,vy1c->GetUymin(),uylim4,ivy_hist->GetMaximum());
      luylim4->SetLineColor(kGray+1);
      TLine *llylim5=new TLine(lylim5,vy1c->GetUymin(),lylim5,ivy_hist->GetMaximum());
      llylim5->SetLineColor(kCyan-1);
      TLine *luylim5=new TLine(uylim5,vy1c->GetUymin(),uylim5,ivy_hist->GetMaximum());
      luylim5->SetLineColor(kCyan-1);
   ivy_hist->Draw("HIST");
   lly13->Draw("HIST Same");
   luy13->Draw("HIST Same");
   lly14->Draw("HIST Same");
   luy14->Draw("HIST Same");
   lly15->Draw("HIST Same");
   luy15->Draw("HIST Same");
   llylim3->Draw("HIST Same");
   luylim3->Draw("HIST Same");
   llylim4->Draw("HIST Same");
   luylim4->Draw("HIST Same");
   llylim5->Draw("HIST Same");
   luylim5->Draw("HIST Same");
   vy1c->Write();


   TCanvas *vz1c = new TCanvas("vz1c","vz1",0,0,640,480);
      TLine *llvy5c=new TLine(lvy5c,vz1c->GetUymin(),lvy5c,ivz_hist->GetMaximum());
      llvy5c->SetLineColor(kRed);
      llvy5c->SetLineWidth(3);
      TLine *luvy5c=new TLine(uvy5c,vz1c->GetUymin(),uvy5c,ivz_hist->GetMaximum());
      luvy5c->SetLineColor(kRed);
      luvy5c->SetLineWidth(3);
      TLine *llvx6c=new TLine(lvx6c,vz1c->GetUymin(),lvx6c,ivz_hist->GetMaximum());
      llvx6c->SetLineColor(kBlack);
      llvx6c->SetLineWidth(3);
      TLine *luvx6c=new TLine(uvx6c,vz1c->GetUymin(),uvx6c,ivz_hist->GetMaximum());
      luvx6c->SetLineColor(kBlack);
      luvx6c->SetLineWidth(3);
      TLine *llvy6c=new TLine(lvy6c,vz1c->GetUymin(),lvy6c,ivz_hist->GetMaximum());
      llvy6c->SetLineColor(kRed+2);
      llvy6c->SetLineWidth(3);
      TLine *luvy6c=new TLine(uvy6c,vz1c->GetUymin(),uvy6c,ivz_hist->GetMaximum());
      luvy6c->SetLineColor(kRed+2);
      luvy6c->SetLineWidth(3);
      TLine *lls13=new TLine(ls13,vz1c->GetUymin(),ls13,ivz_hist->GetMaximum());
      lls13->SetLineColor(kGreen);
      lls13->SetLineWidth(3);
      lls13->SetLineStyle(10);
      TLine *lus13=new TLine(us13,vz1c->GetUymin(),us13,ivz_hist->GetMaximum());
      lus13->SetLineColor(kGreen);
      lus13->SetLineWidth(3);
      lus13->SetLineStyle(10);
      TLine *lls14=new TLine(ls14,vz1c->GetUymin(),ls14,ivz_hist->GetMaximum());
      lls14->SetLineColor(kMagenta);
      lls14->SetLineWidth(3);
      lls14->SetLineStyle(10);
      TLine *lus14=new TLine(us14,vz1c->GetUymin(),us14,ivz_hist->GetMaximum());
      lus14->SetLineColor(kMagenta);
      lus14->SetLineWidth(3);
      lus14->SetLineStyle(10);
      TLine *lls15=new TLine(ls15,vz1c->GetUymin(),ls15,ivz_hist->GetMaximum());
      lls15->SetLineColor(kGreen+2);
      lls15->SetLineWidth(3);
      lls15->SetLineStyle(10);
      TLine *lus15=new TLine(us15,vz1c->GetUymin(),us15,ivz_hist->GetMaximum());
      lus15->SetLineColor(kGreen+2);
      lus15->SetLineWidth(3);
      lus15->SetLineStyle(10);
      TLine *llslim3=new TLine(ls3,vz1c->GetUymin(),ls3,ivz_hist->GetMaximum());
      llslim3->SetLineColor(kGreen+2);
      llslim3->SetLineWidth(3);
      llslim3->SetLineStyle(10);
      TLine *luslim3=new TLine(us3,vz1c->GetUymin(),us3,ivz_hist->GetMaximum());
      luslim3->SetLineColor(kGreen+2);
      luslim3->SetLineWidth(3);
      luslim3->SetLineStyle(10);
      TLine *llslim4=new TLine(ls4,vz1c->GetUymin(),ls4,ivz_hist->GetMaximum());
      llslim4->SetLineColor(kMagenta-7);
      llslim4->SetLineWidth(3);
      llslim4->SetLineStyle(10);
      TLine *luslim4=new TLine(us4,vz1c->GetUymin(),us4,ivz_hist->GetMaximum());
      luslim4->SetLineColor(kMagenta-7);
      luslim4->SetLineWidth(3);
      luslim4->SetLineStyle(10);
      TLine *llslim5=new TLine(ls5,vz1c->GetUymin(),ls5,ivz_hist->GetMaximum());
      llslim5->SetLineColor(kMagenta-7);
      llslim5->SetLineWidth(3);
      llslim5->SetLineStyle(10);
      TLine *luslim5=new TLine(us5,vz1c->GetUymin(),us5,ivz_hist->GetMaximum());
      luslim5->SetLineColor(kMagenta-7);
      luslim5->SetLineWidth(3);
      luslim5->SetLineStyle(10);
      TLine *llclim3=new TLine(lc3,vz1c->GetUymin(),lc3,ivz_hist->GetMaximum());
      llclim3->SetLineColor(kRed-7);
      llclim3->SetLineWidth(3);
      TLine *luclim3=new TLine(uc3,vz1c->GetUymin(),uc3,ivz_hist->GetMaximum());
      luclim3->SetLineColor(kRed-7);
      luclim3->SetLineWidth(3);
      TLine *llclim4=new TLine(lvx1c,vz1c->GetUymin(),lvx1c,ivz_hist->GetMaximum());
      llclim4->SetLineColor(kGray+1);
      llclim4->SetLineWidth(3);
      TLine *luclim4=new TLine(uvx1c,vz1c->GetUymin(),uvx1c,ivz_hist->GetMaximum());
      luclim4->SetLineColor(kGray+1);
      luclim4->SetLineWidth(3);
      TLine *llclim5=new TLine(lvy1c,vz1c->GetUymin(),lvy1c,ivz_hist->GetMaximum());
      llclim5->SetLineColor(kCyan-1);
      llclim5->SetLineWidth(3);
      TLine *luclim5=new TLine(uvy1c,vz1c->GetUymin(),uvy1c,ivz_hist->GetMaximum());
      luclim5->SetLineColor(kCyan-1);
      luclim5->SetLineWidth(3);
   ivz_hist->Draw("HIST");
   lls13->Draw("HIST Same");
   lus13->Draw("HIST Same");
   lls14->Draw("HIST Same");
   lus14->Draw("HIST Same");
   //lls15->Draw("HIST Same");
   //lus15->Draw("HIST Same");
   llvy5c->Draw("HIST Same");
   luvy5c->Draw("HIST Same");
   llvx6c->Draw("HIST Same");
   luvx6c->Draw("HIST Same");
   //llvy6c->Draw("HIST Same");
   //luvy6c->Draw("HIST Same");
   llclim3->Draw("HIST Same");
   luclim3->Draw("HIST Same");
   llclim4->Draw("HIST Same");
   luclim4->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   llslim3->Draw("HIST Same");
   luslim3->Draw("HIST Same");
   llslim4->Draw("HIST Same");
   luslim4->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vz1c->Write();
   



//std::cout << iivx_hist->GetRMS() << std::endl;

TCanvas *vx2c = new TCanvas("vx2c","vx2",0,0,640,480);
      TLine *llx23=new TLine(lx23,vx2c->GetUymin(),lx23,iivx_hist->GetMaximum());
      llx23->SetLineColor(kRed);
      TLine *lux23=new TLine(ux23,vx2c->GetUymin(),ux23,iivx_hist->GetMaximum());
      lux23->SetLineColor(kRed);
      TLine *llx24=new TLine(lx24,vx2c->GetUymin(),lx24,iivx_hist->GetMaximum());
      llx24->SetLineColor(kBlack);
      TLine *lux24=new TLine(ux24,vx2c->GetUymin(),ux24,iivx_hist->GetMaximum());
      lux24->SetLineColor(kBlack);
      TLine *llx25=new TLine(lx25,vx2c->GetUymin(),lx25,iivx_hist->GetMaximum());
      llx25->SetLineColor(kRed+2);
      TLine *lux25=new TLine(ux25,vx2c->GetUymin(),ux25,iivx_hist->GetMaximum());
      lux25->SetLineColor(kRed+2);
   iivx_hist->GetXaxis()->SetTitle("vx [cm]");
   iivx_hist->Draw("HIST");
   llx23->Draw("HIST Same");
   lux23->Draw("HIST Same");
   llx24->Draw("HIST Same");
   lux24->Draw("HIST Same");
   llx25->Draw("HIST Same");
   lux25->Draw("HIST Same");
   llxlim3->Draw("HIST Same");
   luxlim3->Draw("HIST Same");
   llxlim4->Draw("HIST Same");
   luxlim4->Draw("HIST Same");
   llxlim5->Draw("HIST Same");
   luxlim5->Draw("HIST Same");
   vx2c->Write();

TCanvas *vy2c = new TCanvas("vy2c","vy2",0,0,640,480);
      TLine *lly23=new TLine(ly23,vy2c->GetUymin(),ly23,iivy_hist->GetMaximum());
      lly23->SetLineColor(kRed);
      TLine *luy23=new TLine(uy23,vy2c->GetUymin(),uy23,iivy_hist->GetMaximum());
      luy23->SetLineColor(kRed);
      TLine *lly24=new TLine(ly24,vy2c->GetUymin(),ly24,iivy_hist->GetMaximum());
      lly24->SetLineColor(kBlack);
      TLine *luy24=new TLine(uy24,vy2c->GetUymin(),uy24,iivy_hist->GetMaximum());
      luy24->SetLineColor(kBlack);
      TLine *lly25=new TLine(ly25,vy2c->GetUymin(),ly25,iivy_hist->GetMaximum());
      lly25->SetLineColor(kRed+2);
      TLine *luy25=new TLine(uy25,vy2c->GetUymin(),uy25,iivy_hist->GetMaximum());
      luy25->SetLineColor(kRed+2);
   iivy_hist->GetXaxis()->SetTitle("vy [cm]");
   iivy_hist->Draw("HIST");
   lly23->Draw("HIST Same");
   luy23->Draw("HIST Same");
   lly24->Draw("HIST Same");
   luy24->Draw("HIST Same");
   lly25->Draw("HIST Same");
   luy25->Draw("HIST Same");
   llylim3->Draw("HIST Same");
   luylim3->Draw("HIST Same");
   llylim4->Draw("HIST Same");
   luylim4->Draw("HIST Same");
   llylim5->Draw("HIST Same");
   luylim5->Draw("HIST Same");
   vy2c->Write();


TCanvas *vz2c = new TCanvas("vz2c","vz2",0,0,640,480);
      TLine *llc23=new TLine(lc23,vz2c->GetUymin(),lc23,iivz_hist->GetMaximum());
      llc23->SetLineColor(kRed);
      TLine *luc23=new TLine(uc23,vz2c->GetUymin(),uc23,iivz_hist->GetMaximum());
      luc23->SetLineColor(kRed);
      TLine *llc24=new TLine(lc24,vz2c->GetUymin(),lc24,iivz_hist->GetMaximum());
      llc24->SetLineColor(kBlack);
      TLine *luc24=new TLine(uc24,vz2c->GetUymin(),uc24,iivz_hist->GetMaximum());
      luc24->SetLineColor(kBlack);
      TLine *llc25=new TLine(lc25,vz2c->GetUymin(),lc25,iivz_hist->GetMaximum());
      llc25->SetLineColor(kRed+2);
      TLine *luc25=new TLine(uc25,vz2c->GetUymin(),uc25,iivz_hist->GetMaximum());
      luc25->SetLineColor(kRed+2);
      TLine *lls23=new TLine(ls23,vz2c->GetUymin(),ls23,iivz_hist->GetMaximum());
      lls23->SetLineColor(kGreen);
      lls23->SetLineWidth(3);
      lls23->SetLineStyle(10);
      TLine *lus23=new TLine(us23,vz2c->GetUymin(),us23,iivz_hist->GetMaximum());
      lus23->SetLineColor(kGreen);
      lus23->SetLineWidth(3);
      lus23->SetLineStyle(10);
      TLine *lls24=new TLine(ls24,vz2c->GetUymin(),ls24,iivz_hist->GetMaximum());
      lls24->SetLineColor(kMagenta);
      lls24->SetLineWidth(3);
      lls24->SetLineStyle(10);
      TLine *lus24=new TLine(us24,vz2c->GetUymin(),us24,iivz_hist->GetMaximum());
      lus24->SetLineColor(kMagenta);
      lus24->SetLineWidth(3);
      lus24->SetLineStyle(10);
      TLine *lls25=new TLine(ls25,vz2c->GetUymin(),ls25,iivz_hist->GetMaximum());
      lls25->SetLineColor(kGreen+2);
      lls25->SetLineWidth(3);
      lls25->SetLineStyle(10);
      TLine *lus25=new TLine(us25,vz2c->GetUymin(),us25,iivz_hist->GetMaximum());
      lus25->SetLineColor(kGreen+2);
      lus25->SetLineWidth(3);
      lus25->SetLineStyle(10);
      iivz_hist->Draw("HIST");
   lls23->Draw("HIST Same");
   lus23->Draw("HIST Same");
   lls24->Draw("HIST Same");
   lus24->Draw("HIST Same");
   //lls25->Draw("HIST Same");
   //lus25->Draw("HIST Same");
   llc23->Draw("HIST Same");
   luc23->Draw("HIST Same");
   llc24->Draw("HIST Same");
   luc24->Draw("HIST Same");
   //llc25->Draw("HIST Same");
   //luc25->Draw("HIST Same");
   llclim3->Draw("HIST Same");
   luclim3->Draw("HIST Same");
   llclim4->Draw("HIST Same");
   luclim4->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   llslim3->Draw("HIST Same");
   luslim3->Draw("HIST Same");
   llslim4->Draw("HIST Same");
   luslim4->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vz2c->Write();
   




    TCanvas *vx3c = new TCanvas("vx3c","vx3",0,0,640,480);
      TLine *llx33=new TLine(lx33,vx3c->GetUymin(),lx33,iiivx_hist->GetMaximum());
      llx33->SetLineColor(kRed);
      TLine *lux33=new TLine(ux33,vx3c->GetUymin(),ux33,iiivx_hist->GetMaximum());
      lux33->SetLineColor(kRed);
      TLine *llx34=new TLine(lx34,vx3c->GetUymin(),lx34,iiivx_hist->GetMaximum());
      llx34->SetLineColor(kBlack);
      TLine *lux34=new TLine(ux34,vx3c->GetUymin(),ux34,iiivx_hist->GetMaximum());
      lux34->SetLineColor(kBlack);
      TLine *llx35=new TLine(lx35,vx3c->GetUymin(),lx35,iiivx_hist->GetMaximum());
      llx35->SetLineColor(kRed+2);
      TLine *lux35=new TLine(ux35,vx3c->GetUymin(),ux35,iiivx_hist->GetMaximum());
      lux35->SetLineColor(kRed+2);
   iiivx_hist->GetXaxis()->SetTitle("vx [cm]");
   iiivx_hist->Draw("HIST");
   llx33->Draw("HIST Same");
   lux33->Draw("HIST Same");
   llx34->Draw("HIST Same");
   lux34->Draw("HIST Same");
   llx35->Draw("HIST Same");
   lux35->Draw("HIST Same");
   llxlim3->Draw("HIST Same");
   luxlim3->Draw("HIST Same");
   llxlim4->Draw("HIST Same");
   luxlim4->Draw("HIST Same");
   llxlim5->Draw("HIST Same");
   luxlim5->Draw("HIST Same");
   vx3c->Write();

TCanvas *vy3c = new TCanvas("vy3c","vy3",0,0,640,480);
      TLine *lly33=new TLine(ly33,vy3c->GetUymin(),ly33,iiivy_hist->GetMaximum());
      lly33->SetLineColor(kRed);
      TLine *luy33=new TLine(uy33,vy3c->GetUymin(),uy33,iiivy_hist->GetMaximum());
      luy33->SetLineColor(kRed);
      TLine *lly34=new TLine(ly34,vy3c->GetUymin(),ly34,iiivy_hist->GetMaximum());
      lly34->SetLineColor(kBlack);
      TLine *luy34=new TLine(uy34,vy3c->GetUymin(),uy34,iiivy_hist->GetMaximum());
      luy34->SetLineColor(kBlack);
      TLine *lly35=new TLine(ly35,vy3c->GetUymin(),ly35,iiivy_hist->GetMaximum());
      lly35->SetLineColor(kRed+2);
      TLine *luy35=new TLine(uy35,vy3c->GetUymin(),uy35,iiivy_hist->GetMaximum());
      luy35->SetLineColor(kRed+2);
   iiivy_hist->Draw("HIST");
   lly33->Draw("HIST Same");
   luy33->Draw("HIST Same");
   lly34->Draw("HIST Same");
   luy34->Draw("HIST Same");
   lly35->Draw("HIST Same");
   luy35->Draw("HIST Same");
   llylim3->Draw("HIST Same");
   luylim3->Draw("HIST Same");
   llylim4->Draw("HIST Same");
   luylim4->Draw("HIST Same");
   llylim5->Draw("HIST Same");
   luylim5->Draw("HIST Same");
   vy3c->Write();

   

TCanvas *vz3c = new TCanvas("vz3c","vz3",0,0,640,480);
      TLine *llc33=new TLine(lc33,vz3c->GetUymin(),lc33,iiivz_hist->GetMaximum());
      llc33->SetLineColor(kRed);
      TLine *luc33=new TLine(uc33,vz3c->GetUymin(),uc33,iiivz_hist->GetMaximum());
      luc33->SetLineColor(kRed);
      TLine *llc34=new TLine(lc34,vz3c->GetUymin(),lc34,iiivz_hist->GetMaximum());
      llc34->SetLineColor(kBlack);
      TLine *luc34=new TLine(uc34,vz3c->GetUymin(),uc34,iiivz_hist->GetMaximum());
      luc34->SetLineColor(kBlack);
      TLine *llc35=new TLine(lc35,vz3c->GetUymin(),lc35,iiivz_hist->GetMaximum());
      llc35->SetLineColor(kRed+2);
      TLine *luc35=new TLine(uc35,vz3c->GetUymin(),uc35,iiivz_hist->GetMaximum());
      luc35->SetLineColor(kRed+2);
      TLine *lls33=new TLine(ls33,vz3c->GetUymin(),ls33,iiivz_hist->GetMaximum());
      lls33->SetLineColor(kGreen);
      lls33->SetLineWidth(3);
      lls33->SetLineStyle(10);
      TLine *lus33=new TLine(us33,vz3c->GetUymin(),us33,iiivz_hist->GetMaximum());
      lus33->SetLineColor(kGreen);
      lus33->SetLineWidth(3);
      lus33->SetLineStyle(10);
      TLine *lls34=new TLine(ls34,vz3c->GetUymin(),ls34,iiivz_hist->GetMaximum());
      lls34->SetLineColor(kMagenta);
      lls34->SetLineWidth(3);
      lls34->SetLineStyle(10);
      TLine *lus34=new TLine(us34,vz3c->GetUymin(),us34,iiivz_hist->GetMaximum());
      lus34->SetLineColor(kMagenta);
      lus34->SetLineWidth(3);
      lus34->SetLineStyle(10);
      TLine *lls35=new TLine(ls35,vz3c->GetUymin(),ls35,iiivz_hist->GetMaximum());
      lls35->SetLineColor(kGreen+2);
      lls35->SetLineWidth(3);
      lls35->SetLineStyle(10);
      TLine *lus35=new TLine(us35,vz3c->GetUymin(),us35,iiivz_hist->GetMaximum());
      lus35->SetLineColor(kGreen+2);
      lus35->SetLineWidth(3);
      lus35->SetLineStyle(10);
      iiivz_hist->Draw("HIST");
   lls33->Draw("HIST Same");
   lus33->Draw("HIST Same");
   lls34->Draw("HIST Same");
   lus34->Draw("HIST Same");
   //lls35->Draw("HIST Same");
   //lus35->Draw("HIST Same");
   llc33->Draw("HIST Same");
   luc33->Draw("HIST Same");
   llc34->Draw("HIST Same");
   luc34->Draw("HIST Same");
   //llc35->Draw("HIST Same");
   //luc35->Draw("HIST Same");
   llclim3->Draw("HIST Same");
   luclim3->Draw("HIST Same");
   llclim4->Draw("HIST Same");
   luclim4->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   llslim3->Draw("HIST Same");
   luslim3->Draw("HIST Same");
   llslim4->Draw("HIST Same");
   luslim4->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vz3c->Write();





TCanvas *vx4c = new TCanvas("vx4c","vx4",0,0,640,480);
      TLine *llx43=new TLine(lx43,vx4c->GetUymin(),lx43,ivvx_hist->GetMaximum());
      llx43->SetLineColor(kRed);
      TLine *lux43=new TLine(ux43,vx4c->GetUymin(),ux43,ivvx_hist->GetMaximum());
      lux43->SetLineColor(kRed);
      TLine *llx44=new TLine(lx44,vx4c->GetUymin(),lx44,ivvx_hist->GetMaximum());
      llx44->SetLineColor(kBlack);
      TLine *lux44=new TLine(ux44,vx4c->GetUymin(),ux44,ivvx_hist->GetMaximum());
      lux44->SetLineColor(kBlack);
      TLine *llx45=new TLine(lx45,vx4c->GetUymin(),lx45,ivvx_hist->GetMaximum());
      llx45->SetLineColor(kRed+2);
      TLine *lux45=new TLine(ux45,vx4c->GetUymin(),ux45,ivvx_hist->GetMaximum());
      lux45->SetLineColor(kRed+2);
   ivvx_hist->Draw("HIST");
   llx43->Draw("HIST Same");
   lux43->Draw("HIST Same");
   llx44->Draw("HIST Same");
   lux44->Draw("HIST Same");
   llx45->Draw("HIST Same");
   lux45->Draw("HIST Same");
   llxlim3->Draw("HIST Same");
   luxlim3->Draw("HIST Same");
   llxlim4->Draw("HIST Same");
   luxlim4->Draw("HIST Same");
   llxlim5->Draw("HIST Same");
   luxlim5->Draw("HIST Same");
   vx4c->Write();

TCanvas *vy4c = new TCanvas("vy4c","vy4",0,0,640,480);
      TLine *lly43=new TLine(ly43,vy4c->GetUymin(),ly43,ivvy_hist->GetMaximum());
      lly43->SetLineColor(kRed);
      TLine *luy43=new TLine(uy43,vy4c->GetUymin(),uy43,ivvy_hist->GetMaximum());
      luy43->SetLineColor(kRed);
      TLine *lly44=new TLine(ly44,vy4c->GetUymin(),ly44,ivvy_hist->GetMaximum());
      lly44->SetLineColor(kBlack);
      TLine *luy44=new TLine(uy44,vy4c->GetUymin(),uy44,ivvy_hist->GetMaximum());
      luy44->SetLineColor(kBlack);
      TLine *lly45=new TLine(ly45,vy4c->GetUymin(),ly45,ivvy_hist->GetMaximum());
      lly45->SetLineColor(kRed+2);
      TLine *luy45=new TLine(uy45,vy4c->GetUymin(),uy45,ivvy_hist->GetMaximum());
      luy45->SetLineColor(kRed+2);
   ivvy_hist->Draw("HIST");
   lly43->Draw("HIST Same");
   luy43->Draw("HIST Same");
   lly44->Draw("HIST Same");
   luy44->Draw("HIST Same");
   lly45->Draw("HIST Same");
   luy45->Draw("HIST Same");
   llylim3->Draw("HIST Same");
   luylim3->Draw("HIST Same");
   llylim4->Draw("HIST Same");
   luylim4->Draw("HIST Same");
   llylim5->Draw("HIST Same");
   luylim5->Draw("HIST Same");
   vy4c->Write();
  

   

TCanvas *vz4c = new TCanvas("vz4c","vz4",0,0,640,480);
      TLine *llvx1c3=new TLine(lvx1c3,vz4c->GetUymin(),lvx1c3,ivvz_hist->GetMaximum());
      llvx1c3->SetLineColor(kRed);
      TLine *luvx1c3=new TLine(uvx1c3,vz4c->GetUymin(),uvx1c3,ivvz_hist->GetMaximum());
      luvx1c3->SetLineColor(kRed);
      TLine *llvx1c4=new TLine(lvx1c4,vz4c->GetUymin(),lvx1c4,ivvz_hist->GetMaximum());
      llvx1c4->SetLineColor(kBlack);
      TLine *luvx1c4=new TLine(uvx1c4,vz4c->GetUymin(),uvx1c4,ivvz_hist->GetMaximum());
      luvx1c4->SetLineColor(kBlack);
      TLine *llvx1vy1c=new TLine(lvx1vy1c,vz4c->GetUymin(),lvx1vy1c,ivvz_hist->GetMaximum());
      llvx1vy1c->SetLineColor(kRed+2);
      TLine *luvx1vy1c=new TLine(uvx1vy1c,vz4c->GetUymin(),uvx1vy1c,ivvz_hist->GetMaximum());
      luvx1vy1c->SetLineColor(kRed+2);
      TLine *lls43=new TLine(ls43,vz4c->GetUymin(),ls43,ivvz_hist->GetMaximum());
      lls43->SetLineColor(kGreen);
      lls43->SetLineWidth(3);
      lls43->SetLineStyle(10);
      TLine *lus43=new TLine(us43,vz4c->GetUymin(),us43,ivvz_hist->GetMaximum());
      lus43->SetLineColor(kGreen);
      lus43->SetLineWidth(3);
      lus43->SetLineStyle(10);
      TLine *lls44=new TLine(ls44,vz4c->GetUymin(),ls44,ivvz_hist->GetMaximum());
      lls44->SetLineColor(kMagenta);
      lls44->SetLineWidth(3);
      lls44->SetLineStyle(10);
      TLine *lus44=new TLine(us44,vz4c->GetUymin(),us44,ivvz_hist->GetMaximum());
      lus44->SetLineColor(kMagenta);
      lus44->SetLineWidth(3);
      lus44->SetLineStyle(10);
      TLine *lls45=new TLine(ls45,vz4c->GetUymin(),ls45,ivvz_hist->GetMaximum());
      lls45->SetLineColor(kGreen+2);
      lls45->SetLineWidth(3);
      lls45->SetLineStyle(10);
      TLine *lus45=new TLine(us45,vz4c->GetUymin(),us45,ivvz_hist->GetMaximum());
      lus45->SetLineColor(kGreen+2);
      lus45->SetLineWidth(3);
      lus45->SetLineStyle(10);
      ivvz_hist->Draw("HIST");
   lls43->Draw("HIST Same");
   lus43->Draw("HIST Same");
   lls44->Draw("HIST Same");
   lus44->Draw("HIST Same");
   //lls45->Draw("HIST Same");
   //lus45->Draw("HIST Same");
   llvx1c3->Draw("HIST Same");
   luvx1c3->Draw("HIST Same");
   llvx1c4->Draw("HIST Same");
   luvx1c4->Draw("HIST Same");
   //llvx1vy1c->Draw("HIST Same");
   //luvx1vy1c->Draw("HIST Same");
   llclim3->Draw("HIST Same");
   luclim3->Draw("HIST Same");
   llclim4->Draw("HIST Same");
   luclim4->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   llslim3->Draw("HIST Same");
   luslim3->Draw("HIST Same");
   llslim4->Draw("HIST Same");
   luslim4->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vz4c->Write();
   




TCanvas *vx5c = new TCanvas("vx5c","vx5",0,0,640,480);
      TLine *llx53=new TLine(lx53,vx5c->GetUymin(),lx53,vvx_hist->GetMaximum());
      llx53->SetLineColor(kRed);
      TLine *lux53=new TLine(ux53,vx5c->GetUymin(),ux53,vvx_hist->GetMaximum());
      lux53->SetLineColor(kRed);
      TLine *llx54=new TLine(lx54,vx5c->GetUymin(),lx54,vvx_hist->GetMaximum());
      llx54->SetLineColor(kBlack);
      TLine *lux54=new TLine(ux54,vx5c->GetUymin(),ux54,vvx_hist->GetMaximum());
      lux54->SetLineColor(kBlack);
      TLine *llx55=new TLine(lx55,vx5c->GetUymin(),lx55,vvx_hist->GetMaximum());
      llx55->SetLineColor(kRed+2);
      TLine *lux55=new TLine(ux55,vx5c->GetUymin(),ux55,vvx_hist->GetMaximum());
      lux55->SetLineColor(kRed+2);
   vvx_hist->GetXaxis()->SetTitle("vx [cm]");
   vvx_hist->Draw("HIST");
   llx53->Draw("HIST Same");
   lux53->Draw("HIST Same");
   llx54->Draw("HIST Same");
   lux54->Draw("HIST Same");
   llx55->Draw("HIST Same");
   lux55->Draw("HIST Same");
   llxlim3->Draw("HIST Same");
   luxlim3->Draw("HIST Same");
   llxlim4->Draw("HIST Same");
   luxlim4->Draw("HIST Same");
   llxlim5->Draw("HIST Same");
   luxlim5->Draw("HIST Same");
   vx5c->Write();
  

TCanvas *vy5c = new TCanvas("vy5c","vy5",0,0,640,480);
      TLine *lly53=new TLine(ly53,vy5c->GetUymin(),ly53,vvy_hist->GetMaximum());
      lly53->SetLineColor(kRed);
      TLine *luy53=new TLine(uy53,vy5c->GetUymin(),uy53,vvy_hist->GetMaximum());
      luy53->SetLineColor(kRed);
      TLine *lly54=new TLine(ly54,vy5c->GetUymin(),ly54,vvy_hist->GetMaximum());
      lly54->SetLineColor(kBlack);
      TLine *luy54=new TLine(uy54,vy5c->GetUymin(),uy54,vvy_hist->GetMaximum());
      luy54->SetLineColor(kBlack);
      TLine *lly55=new TLine(ly55,vy5c->GetUymin(),ly55,vvy_hist->GetMaximum());
      lly55->SetLineColor(kRed+2);
      TLine *luy55=new TLine(uy55,vy5c->GetUymin(),uy55,vvy_hist->GetMaximum());
      luy55->SetLineColor(kRed+2);
   vvy_hist->GetXaxis()->SetTitle("vy [cm]");
   vvy_hist->Draw("HIST");
   lly53->Draw("HIST Same");
   luy53->Draw("HIST Same");
   lly54->Draw("HIST Same");
   luy54->Draw("HIST Same");
   lly55->Draw("HIST Same");
   luy55->Draw("HIST Same");
   llylim3->Draw("HIST Same");
   luylim3->Draw("HIST Same");
   llylim4->Draw("HIST Same");
   luylim4->Draw("HIST Same");
   llylim5->Draw("HIST Same");
   luylim5->Draw("HIST Same");
   vy5c->Write();
   


TCanvas *vz5c = new TCanvas("vz5c","vz5",0,0,640,480);
      TLine *llvy1c3=new TLine(lvy1c3,vz5c->GetUymin(),lvy1c3,vvz_hist->GetMaximum());
      llvy1c3->SetLineColor(kRed);
      TLine *luvy1c3=new TLine(uvy1c3,vz5c->GetUymin(),uvy1c3,vvz_hist->GetMaximum());
      luvy1c3->SetLineColor(kRed);
      TLine *llvy1c4=new TLine(lvy1c4,vz5c->GetUymin(),lvy1c4,vvz_hist->GetMaximum());
      llvy1c4->SetLineColor(kBlack);
      TLine *luvy1c4=new TLine(uvy1c4,vz5c->GetUymin(),uvy1c4,vvz_hist->GetMaximum());
      luvy1c4->SetLineColor(kBlack);
      TLine *llvy1c5=new TLine(lvy1c5,vz5c->GetUymin(),lvy1c5,vvz_hist->GetMaximum());
      llvy1c5->SetLineColor(kRed+2);
      TLine *luvy1c5=new TLine(uvy1c5,vz5c->GetUymin(),uvy1c5,vvz_hist->GetMaximum());
      luvy1c5->SetLineColor(kRed+2);
      TLine *lls53=new TLine(ls53,vz5c->GetUymin(),ls53,vvz_hist->GetMaximum());
      lls53->SetLineColor(kGreen);
      lls53->SetLineWidth(3);
      lls53->SetLineStyle(10);
      TLine *lus53=new TLine(us53,vz5c->GetUymin(),us53,vvz_hist->GetMaximum());
      lus53->SetLineColor(kGreen);
      lus53->SetLineWidth(3);
      lus53->SetLineStyle(10);
      TLine *lls54=new TLine(ls54,vz5c->GetUymin(),ls54,vvz_hist->GetMaximum());
      lls54->SetLineColor(kMagenta);
      lls54->SetLineWidth(3);
      lls54->SetLineStyle(10);
      TLine *lus54=new TLine(us54,vz5c->GetUymin(),us54,vvz_hist->GetMaximum());
      lus54->SetLineColor(kMagenta);
      lus54->SetLineWidth(3);
      lus54->SetLineStyle(10);
      TLine *lls55=new TLine(ls55,vz5c->GetUymin(),ls55,vvz_hist->GetMaximum());
      lls55->SetLineColor(kGreen+2);
      lls55->SetLineWidth(3);
      lls55->SetLineStyle(10);
      TLine *lus55=new TLine(us55,vz5c->GetUymin(),us55,vvz_hist->GetMaximum());
      lus55->SetLineColor(kGreen+2);
      lus55->SetLineWidth(3);
      lus55->SetLineStyle(10);
      vvz_hist->Draw("HIST");
   lls53->Draw("HIST Same");
   lus53->Draw("HIST Same");
   lls54->Draw("HIST Same");
   lus54->Draw("HIST Same");
   //lls55->Draw("HIST Same");
   //lus55->Draw("HIST Same");
   llvy1c3->Draw("HIST Same");
   luvy1c3->Draw("HIST Same");
   llvy1c4->Draw("HIST Same");
   luvy1c4->Draw("HIST Same");
   //llvy1c5->Draw("HIST Same");
   //luvy1c5->Draw("HIST Same");
   llclim3->Draw("HIST Same");
   luclim3->Draw("HIST Same");
   llclim4->Draw("HIST Same");
   luclim4->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   llslim3->Draw("HIST Same");
   luslim3->Draw("HIST Same");
   llslim4->Draw("HIST Same");
   luslim4->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vz5c->Write();
  

 


TCanvas *vx6c = new TCanvas("vx6c","vx6",0,0,640,480);
      TLine *llx63=new TLine(lx63,vx6c->GetUymin(),lx63,vivx_hist->GetMaximum());
      llx63->SetLineColor(kRed);
      TLine *lux63=new TLine(ux63,vx6c->GetUymin(),ux63,vivx_hist->GetMaximum());
      lux63->SetLineColor(kRed);
      TLine *llx64=new TLine(lx64,vx6c->GetUymin(),lx64,vivx_hist->GetMaximum());
      llx64->SetLineColor(kBlack);
      TLine *lux64=new TLine(ux64,vx6c->GetUymin(),ux64,vivx_hist->GetMaximum());
      lux64->SetLineColor(kBlack);
      TLine *llx65=new TLine(lx65,vx6c->GetUymin(),lx65,vivx_hist->GetMaximum());
      llx65->SetLineColor(kRed+2);
      TLine *lux65=new TLine(ux65,vx6c->GetUymin(),ux65,vivx_hist->GetMaximum());
      lux65->SetLineColor(kRed+2);
   vivx_hist->GetXaxis()->SetTitle("vx [cm]");
   vivx_hist->Draw("HIST");
   llx63->Draw("HIST Same");
   lux63->Draw("HIST Same");
   llx64->Draw("HIST Same");
   lux64->Draw("HIST Same");
   llx65->Draw("HIST Same");
   lux65->Draw("HIST Same");
   llxlim3->Draw("HIST Same");
   luxlim3->Draw("HIST Same");
   llxlim4->Draw("HIST Same");
   luxlim4->Draw("HIST Same");
   llxlim5->Draw("HIST Same");
   luxlim5->Draw("HIST Same");
   vx6c->Write();
   

TCanvas *vy6c = new TCanvas("vy6c","vy6",0,0,640,480);
      TLine *lly63=new TLine(ly63,vy6c->GetUymin(),ly63,vivy_hist->GetMaximum());
      lly63->SetLineColor(kRed);
      TLine *luy63=new TLine(uy63,vy6c->GetUymin(),uy63,vivy_hist->GetMaximum());
      luy63->SetLineColor(kRed);
      TLine *lly64=new TLine(ly64,vy6c->GetUymin(),ly64,vivy_hist->GetMaximum());
      lly64->SetLineColor(kBlack);
      TLine *luy64=new TLine(uy64,vy6c->GetUymin(),uy64,vivy_hist->GetMaximum());
      luy64->SetLineColor(kBlack);
      TLine *lly65=new TLine(ly65,vy6c->GetUymin(),ly65,vivy_hist->GetMaximum());
      lly65->SetLineColor(kRed+2);
      TLine *luy65=new TLine(uy65,vy6c->GetUymin(),uy65,vivy_hist->GetMaximum());
      luy65->SetLineColor(kRed+2);
   vivy_hist->GetXaxis()->SetTitle("vy [cm]");
   vivy_hist->Draw("HIST");
   lly63->Draw("HIST Same");
   luy63->Draw("HIST Same");
   lly64->Draw("HIST Same");
   luy64->Draw("HIST Same");
   lly65->Draw("HIST Same");
   luy65->Draw("HIST Same");
   llylim3->Draw("HIST Same");
   luylim3->Draw("HIST Same");
   llylim4->Draw("HIST Same");
   luylim4->Draw("HIST Same");
   llylim5->Draw("HIST Same");
   luylim5->Draw("HIST Same");
   vy6c->Write();
   
   
   
   TCanvas *vz6c = new TCanvas("vz6c","vz6",0,0,640,480);
      TLine *llvx2c3=new TLine(lvx2c3,vz6c->GetUymin(),lvx2c3,vivz_hist->GetMaximum());
      llvx2c3->SetLineColor(kRed);
      TLine *luvx2c3=new TLine(uvx2c3,vz6c->GetUymin(),uvx2c3,vivz_hist->GetMaximum());
      luvx2c3->SetLineColor(kRed);
      TLine *llvx2c4=new TLine(lvx2c4,vz6c->GetUymin(),lvx2c4,vivz_hist->GetMaximum());
      llvx2c4->SetLineColor(kBlack);
      TLine *luvx2c4=new TLine(uvx2c4,vz6c->GetUymin(),uvx2c4,vivz_hist->GetMaximum());
      luvx2c4->SetLineColor(kBlack);
      TLine *llvx2c5=new TLine(lvx2c5,vz6c->GetUymin(),lvx2c5,vivz_hist->GetMaximum());
      llvx2c5->SetLineColor(kRed+2);
      TLine *luvx2c5=new TLine(uvx2c5,vz6c->GetUymin(),uvx2c5,vivz_hist->GetMaximum());
      luvx2c5->SetLineColor(kRed+2);
      TLine *lls63=new TLine(ls63,vz6c->GetUymin(),ls63,vivz_hist->GetMaximum());
      lls63->SetLineColor(kGreen);
      lls63->SetLineWidth(3);
      lls63->SetLineStyle(10);
      TLine *lus63=new TLine(us63,vz6c->GetUymin(),us63,vivz_hist->GetMaximum());
      lus63->SetLineColor(kGreen);
      lus63->SetLineWidth(3);
      lus63->SetLineStyle(10);
      TLine *lls64=new TLine(ls64,vz6c->GetUymin(),ls64,vivz_hist->GetMaximum());
      lls64->SetLineColor(kMagenta);
      lls64->SetLineWidth(3);
      lls64->SetLineStyle(10);
      TLine *lus64=new TLine(us64,vz6c->GetUymin(),us64,vivz_hist->GetMaximum());
      lus64->SetLineColor(kMagenta);
      lus64->SetLineWidth(3);
      lus64->SetLineStyle(10);
      TLine *lls65=new TLine(ls65,vz6c->GetUymin(),ls65,vivz_hist->GetMaximum());
      lls65->SetLineColor(kGreen+2);
      lls65->SetLineWidth(3);
      lls65->SetLineStyle(10);
      TLine *lus65=new TLine(us65,vz6c->GetUymin(),us65,vivz_hist->GetMaximum());
      lus65->SetLineColor(kGreen+2);
      lus65->SetLineWidth(3);
      lus65->SetLineStyle(10);
      vivz_hist->Draw("HIST");
   lls63->Draw("HIST Same");
   lus63->Draw("HIST Same");
   lls64->Draw("HIST Same");
   lus64->Draw("HIST Same");
   //lls65->Draw("HIST Same");
   //lus65->Draw("HIST Same");
   llvx2c3->Draw("HIST Same");
   luvx2c3->Draw("HIST Same");
   llvx2c4->Draw("HIST Same");
   luvx2c4->Draw("HIST Same");
   //llvx2c5->Draw("HIST Same");
   //luvx2c5->Draw("HIST Same");
   llclim3->Draw("HIST Same");
   luclim3->Draw("HIST Same");
   llclim4->Draw("HIST Same");
   luclim4->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   llslim3->Draw("HIST Same");
   luslim3->Draw("HIST Same");
   llslim4->Draw("HIST Same");
   luslim4->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vz6c->Write();


TCanvas *pip_vx1c = new TCanvas("pip_vx1c","pip_vx1",0,0,640,480);
      TLine *pip_llx13=new TLine(pip_lx13,pip_vx1c->GetUymin(),pip_lx13,pip_ivx_hist->GetMaximum());
      pip_llx13->SetLineColor(kRed);
      TLine *pip_lux13=new TLine(pip_ux13,pip_vx1c->GetUymin(),pip_ux13,pip_ivx_hist->GetMaximum());
      pip_lux13->SetLineColor(kRed);
      TLine *pip_llx14=new TLine(pip_lx14,pip_vx1c->GetUymin(),pip_lx14,pip_ivx_hist->GetMaximum());
      pip_llx14->SetLineColor(kBlack);
      TLine *pip_lux14=new TLine(pip_ux14,pip_vx1c->GetUymin(),pip_ux14,pip_ivx_hist->GetMaximum());
      pip_lux14->SetLineColor(kBlack);
      TLine *pip_llx15=new TLine(pip_lx15,pip_vx1c->GetUymin(),pip_lx15,pip_ivx_hist->GetMaximum());
      pip_llx15->SetLineColor(kRed+2);
      TLine *pip_lux15=new TLine(pip_ux15,pip_vx1c->GetUymin(),pip_ux15,pip_ivx_hist->GetMaximum());
      pip_lux15->SetLineColor(kRed+2);
      TLine *pip_llxlim3=new TLine(pip_lxlim3,pip_vx1c->GetUymin(),pip_lxlim3,pip_ivx_hist->GetMaximum());
      pip_llxlim3->SetLineColor(kRed-7);
      TLine *pip_luxlim3=new TLine(pip_uxlim3,pip_vx1c->GetUymin(),pip_uxlim3,pip_ivx_hist->GetMaximum());
      pip_luxlim3->SetLineColor(kRed-7);
      TLine *pip_llxlim4=new TLine(pip_lxlim4,pip_vx1c->GetUymin(),pip_lxlim4,pip_ivx_hist->GetMaximum());
      pip_llxlim4->SetLineColor(kGray+1);
      TLine *pip_luxlim4=new TLine(pip_uxlim4,pip_vx1c->GetUymin(),pip_uxlim4,pip_ivx_hist->GetMaximum());
      pip_luxlim4->SetLineColor(kGray+1);
      TLine *pip_llxlim5=new TLine(pip_lxlim5,pip_vx1c->GetUymin(),pip_lxlim5,pip_ivx_hist->GetMaximum());
      pip_llxlim5->SetLineColor(kCyan-1);
      TLine *pip_luxlim5=new TLine(pip_uxlim5,pip_vx1c->GetUymin(),pip_uxlim5,pip_ivx_hist->GetMaximum());
      pip_luxlim5->SetLineColor(kCyan-1);
   pip_ivx_hist->Draw("HIST");
   pip_llx13->Draw("HIST Same");
   pip_lux13->Draw("HIST Same");
   pip_llx14->Draw("HIST Same");
   pip_lux14->Draw("HIST Same");
   pip_llx15->Draw("HIST Same");
   pip_lux15->Draw("HIST Same");
   pip_llxlim3->Draw("HIST Same");
   pip_luxlim3->Draw("HIST Same");
   pip_llxlim4->Draw("HIST Same");
   pip_luxlim4->Draw("HIST Same");
   pip_llxlim5->Draw("HIST Same");
   pip_luxlim5->Draw("HIST Same");
   pip_vx1c->Write();


   TCanvas *pip_vy1c = new TCanvas("pip_vy1c","pip_vy1",0,0,640,480);
      TLine *pip_lly13=new TLine(pip_ly13,pip_vy1c->GetUymin(),pip_ly13,pip_ivy_hist->GetMaximum());
      pip_lly13->SetLineColor(kRed);
      TLine *pip_luy13=new TLine(pip_uy13,pip_vy1c->GetUymin(),pip_uy13,pip_ivy_hist->GetMaximum());
      pip_luy13->SetLineColor(kRed);
      TLine *pip_lly14=new TLine(pip_ly14,pip_vy1c->GetUymin(),pip_ly14,pip_ivy_hist->GetMaximum());
      pip_lly14->SetLineColor(kBlack);
      TLine *pip_luy14=new TLine(pip_uy14,pip_vy1c->GetUymin(),pip_uy14,pip_ivy_hist->GetMaximum());
      pip_luy14->SetLineColor(kBlack);
      TLine *pip_lly15=new TLine(pip_ly15,pip_vy1c->GetUymin(),pip_ly15,pip_ivy_hist->GetMaximum());
      pip_lly15->SetLineColor(kRed+2);
      TLine *pip_luy15=new TLine(pip_uy15,pip_vy1c->GetUymin(),pip_uy15,pip_ivy_hist->GetMaximum());
      pip_luy15->SetLineColor(kRed+2);
      TLine *pip_llylim3=new TLine(pip_lylim3,pip_vy1c->GetUymin(),pip_lylim3,pip_ivy_hist->GetMaximum());
      pip_llylim3->SetLineColor(kRed-7);
      TLine *pip_luylim3=new TLine(pip_uylim3,pip_vy1c->GetUymin(),pip_uylim3,pip_ivy_hist->GetMaximum());
      pip_luylim3->SetLineColor(kRed-7);
      TLine *pip_llylim4=new TLine(pip_lylim4,pip_vy1c->GetUymin(),pip_lylim4,pip_ivy_hist->GetMaximum());
      pip_llylim4->SetLineColor(kGray+1);
      TLine *pip_luylim4=new TLine(pip_uylim4,pip_vy1c->GetUymin(),pip_uylim4,pip_ivy_hist->GetMaximum());
      pip_luylim4->SetLineColor(kGray+1);
      TLine *pip_llylim5=new TLine(pip_lylim5,pip_vy1c->GetUymin(),pip_lylim5,pip_ivy_hist->GetMaximum());
      pip_llylim5->SetLineColor(kCyan-1);
      TLine *pip_luylim5=new TLine(pip_uylim5,pip_vy1c->GetUymin(),pip_uylim5,pip_ivy_hist->GetMaximum());
      pip_luylim5->SetLineColor(kCyan-1);
   pip_ivy_hist->Draw("HIST");
   pip_lly13->Draw("HIST Same");
   pip_luy13->Draw("HIST Same");
   pip_lly14->Draw("HIST Same");
   pip_luy14->Draw("HIST Same");
   pip_lly15->Draw("HIST Same");
   pip_luy15->Draw("HIST Same");
   pip_llylim3->Draw("HIST Same");
   pip_luylim3->Draw("HIST Same");
   pip_llylim4->Draw("HIST Same");
   pip_luylim4->Draw("HIST Same");
   pip_llylim5->Draw("HIST Same");
   pip_luylim5->Draw("HIST Same");
   pip_vy1c->Write();

std::cout << "zeros = " << zerosp << std::endl;
*/


/*
   TCanvas *vz1c = new TCanvas("vz1c","vz1",0,0,640,480);
      TLine *llvy5c=new TLine(lvy5c,vz1c->GetUymin(),lvy5c,ivz_hist->GetMaximum());
      llvy5c->SetLineColor(kRed);
      llvy5c->SetLineWidth(3);
      TLine *luvy5c=new TLine(uvy5c,vz1c->GetUymin(),uvy5c,ivz_hist->GetMaximum());
      luvy5c->SetLineColor(kRed);
      luvy5c->SetLineWidth(3);
      TLine *llvx6c=new TLine(lvx6c,vz1c->GetUymin(),lvx6c,ivz_hist->GetMaximum());
      llvx6c->SetLineColor(kBlack);
      llvx6c->SetLineWidth(3);
      TLine *luvx6c=new TLine(uvx6c,vz1c->GetUymin(),uvx6c,ivz_hist->GetMaximum());
      luvx6c->SetLineColor(kBlack);
      luvx6c->SetLineWidth(3);
      TLine *llvy6c=new TLine(lvy6c,vz1c->GetUymin(),lvy6c,ivz_hist->GetMaximum());
      llvy6c->SetLineColor(kRed+2);
      llvy6c->SetLineWidth(3);
      TLine *luvy6c=new TLine(uvy6c,vz1c->GetUymin(),uvy6c,ivz_hist->GetMaximum());
      luvy6c->SetLineColor(kRed+2);
      luvy6c->SetLineWidth(3);
      TLine *lls13=new TLine(ls13,vz1c->GetUymin(),ls13,ivz_hist->GetMaximum());
      lls13->SetLineColor(kGreen);
      lls13->SetLineWidth(3);
      lls13->SetLineStyle(10);
      TLine *lus13=new TLine(us13,vz1c->GetUymin(),us13,ivz_hist->GetMaximum());
      lus13->SetLineColor(kGreen);
      lus13->SetLineWidth(3);
      lus13->SetLineStyle(10);
      TLine *lls14=new TLine(ls14,vz1c->GetUymin(),ls14,ivz_hist->GetMaximum());
      lls14->SetLineColor(kMagenta);
      lls14->SetLineWidth(3);
      lls14->SetLineStyle(10);
      TLine *lus14=new TLine(us14,vz1c->GetUymin(),us14,ivz_hist->GetMaximum());
      lus14->SetLineColor(kMagenta);
      lus14->SetLineWidth(3);
      lus14->SetLineStyle(10);
      TLine *lls15=new TLine(ls15,vz1c->GetUymin(),ls15,ivz_hist->GetMaximum());
      lls15->SetLineColor(kGreen+2);
      lls15->SetLineWidth(3);
      lls15->SetLineStyle(10);
      TLine *lus15=new TLine(us15,vz1c->GetUymin(),us15,ivz_hist->GetMaximum());
      lus15->SetLineColor(kGreen+2);
      lus15->SetLineWidth(3);
      lus15->SetLineStyle(10);
      TLine *llslim3=new TLine(ls3,vz1c->GetUymin(),ls3,ivz_hist->GetMaximum());
      llslim3->SetLineColor(kGreen+2);
      llslim3->SetLineWidth(3);
      llslim3->SetLineStyle(10);
      TLine *luslim3=new TLine(us3,vz1c->GetUymin(),us3,ivz_hist->GetMaximum());
      luslim3->SetLineColor(kGreen+2);
      luslim3->SetLineWidth(3);
      luslim3->SetLineStyle(10);
      TLine *llslim4=new TLine(ls4,vz1c->GetUymin(),ls4,ivz_hist->GetMaximum());
      llslim4->SetLineColor(kMagenta-7);
      llslim4->SetLineWidth(3);
      llslim4->SetLineStyle(10);
      TLine *luslim4=new TLine(us4,vz1c->GetUymin(),us4,ivz_hist->GetMaximum());
      luslim4->SetLineColor(kMagenta-7);
      luslim4->SetLineWidth(3);
      luslim4->SetLineStyle(10);
      TLine *llslim5=new TLine(ls5,vz1c->GetUymin(),ls5,ivz_hist->GetMaximum());
      llslim5->SetLineColor(kMagenta-7);
      llslim5->SetLineWidth(3);
      llslim5->SetLineStyle(10);
      TLine *luslim5=new TLine(us5,vz1c->GetUymin(),us5,ivz_hist->GetMaximum());
      luslim5->SetLineColor(kMagenta-7);
      luslim5->SetLineWidth(3);
      luslim5->SetLineStyle(10);
      TLine *llclim3=new TLine(lc3,vz1c->GetUymin(),lc3,ivz_hist->GetMaximum());
      llclim3->SetLineColor(kRed-7);
      llclim3->SetLineWidth(3);
      TLine *luclim3=new TLine(uc3,vz1c->GetUymin(),uc3,ivz_hist->GetMaximum());
      luclim3->SetLineColor(kRed-7);
      luclim3->SetLineWidth(3);
      TLine *llclim4=new TLine(lvx1c,vz1c->GetUymin(),lvx1c,ivz_hist->GetMaximum());
      llclim4->SetLineColor(kGray+1);
      llclim4->SetLineWidth(3);
      TLine *luclim4=new TLine(uvx1c,vz1c->GetUymin(),uvx1c,ivz_hist->GetMaximum());
      luclim4->SetLineColor(kGray+1);
      luclim4->SetLineWidth(3);
      TLine *llclim5=new TLine(lvy1c,vz1c->GetUymin(),lvy1c,ivz_hist->GetMaximum());
      llclim5->SetLineColor(kCyan-1);
      llclim5->SetLineWidth(3);
      TLine *luclim5=new TLine(uvy1c,vz1c->GetUymin(),uvy1c,ivz_hist->GetMaximum());
      luclim5->SetLineColor(kCyan-1);
      luclim5->SetLineWidth(3);
   ivz_hist->Draw("HIST");
   lls13->Draw("HIST Same");
   lus13->Draw("HIST Same");
   lls14->Draw("HIST Same");
   lus14->Draw("HIST Same");
   //lls15->Draw("HIST Same");
   //lus15->Draw("HIST Same");
   llvy5c->Draw("HIST Same");
   luvy5c->Draw("HIST Same");
   llvx6c->Draw("HIST Same");
   luvx6c->Draw("HIST Same");
   //llvy6c->Draw("HIST Same");
   //luvy6c->Draw("HIST Same");
   llclim3->Draw("HIST Same");
   luclim3->Draw("HIST Same");
   llclim4->Draw("HIST Same");
   luclim4->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   llslim3->Draw("HIST Same");
   luslim3->Draw("HIST Same");
   llslim4->Draw("HIST Same");
   luslim4->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vz1c->Write();
   



//std::cout << iivx_hist->GetRMS() << std::endl;

TCanvas *vx2c = new TCanvas("vx2c","vx2",0,0,640,480);
      TLine *llx23=new TLine(lx23,vx2c->GetUymin(),lx23,iivx_hist->GetMaximum());
      llx23->SetLineColor(kRed);
      TLine *lux23=new TLine(ux23,vx2c->GetUymin(),ux23,iivx_hist->GetMaximum());
      lux23->SetLineColor(kRed);
      TLine *llx24=new TLine(lx24,vx2c->GetUymin(),lx24,iivx_hist->GetMaximum());
      llx24->SetLineColor(kBlack);
      TLine *lux24=new TLine(ux24,vx2c->GetUymin(),ux24,iivx_hist->GetMaximum());
      lux24->SetLineColor(kBlack);
      TLine *llx25=new TLine(lx25,vx2c->GetUymin(),lx25,iivx_hist->GetMaximum());
      llx25->SetLineColor(kRed+2);
      TLine *lux25=new TLine(ux25,vx2c->GetUymin(),ux25,iivx_hist->GetMaximum());
      lux25->SetLineColor(kRed+2);
   iivx_hist->GetXaxis()->SetTitle("vx [cm]");
   iivx_hist->Draw("HIST");
   llx23->Draw("HIST Same");
   lux23->Draw("HIST Same");
   llx24->Draw("HIST Same");
   lux24->Draw("HIST Same");
   llx25->Draw("HIST Same");
   lux25->Draw("HIST Same");
   llxlim3->Draw("HIST Same");
   luxlim3->Draw("HIST Same");
   llxlim4->Draw("HIST Same");
   luxlim4->Draw("HIST Same");
   llxlim5->Draw("HIST Same");
   luxlim5->Draw("HIST Same");
   vx2c->Write();

TCanvas *vy2c = new TCanvas("vy2c","vy2",0,0,640,480);
      TLine *lly23=new TLine(ly23,vy2c->GetUymin(),ly23,iivy_hist->GetMaximum());
      lly23->SetLineColor(kRed);
      TLine *luy23=new TLine(uy23,vy2c->GetUymin(),uy23,iivy_hist->GetMaximum());
      luy23->SetLineColor(kRed);
      TLine *lly24=new TLine(ly24,vy2c->GetUymin(),ly24,iivy_hist->GetMaximum());
      lly24->SetLineColor(kBlack);
      TLine *luy24=new TLine(uy24,vy2c->GetUymin(),uy24,iivy_hist->GetMaximum());
      luy24->SetLineColor(kBlack);
      TLine *lly25=new TLine(ly25,vy2c->GetUymin(),ly25,iivy_hist->GetMaximum());
      lly25->SetLineColor(kRed+2);
      TLine *luy25=new TLine(uy25,vy2c->GetUymin(),uy25,iivy_hist->GetMaximum());
      luy25->SetLineColor(kRed+2);
   iivy_hist->GetXaxis()->SetTitle("vy [cm]");
   iivy_hist->Draw("HIST");
   lly23->Draw("HIST Same");
   luy23->Draw("HIST Same");
   lly24->Draw("HIST Same");
   luy24->Draw("HIST Same");
   lly25->Draw("HIST Same");
   luy25->Draw("HIST Same");
   llylim3->Draw("HIST Same");
   luylim3->Draw("HIST Same");
   llylim4->Draw("HIST Same");
   luylim4->Draw("HIST Same");
   llylim5->Draw("HIST Same");
   luylim5->Draw("HIST Same");
   vy2c->Write();


TCanvas *vz2c = new TCanvas("vz2c","vz2",0,0,640,480);
      TLine *llc23=new TLine(lc23,vz2c->GetUymin(),lc23,iivz_hist->GetMaximum());
      llc23->SetLineColor(kRed);
      TLine *luc23=new TLine(uc23,vz2c->GetUymin(),uc23,iivz_hist->GetMaximum());
      luc23->SetLineColor(kRed);
      TLine *llc24=new TLine(lc24,vz2c->GetUymin(),lc24,iivz_hist->GetMaximum());
      llc24->SetLineColor(kBlack);
      TLine *luc24=new TLine(uc24,vz2c->GetUymin(),uc24,iivz_hist->GetMaximum());
      luc24->SetLineColor(kBlack);
      TLine *llc25=new TLine(lc25,vz2c->GetUymin(),lc25,iivz_hist->GetMaximum());
      llc25->SetLineColor(kRed+2);
      TLine *luc25=new TLine(uc25,vz2c->GetUymin(),uc25,iivz_hist->GetMaximum());
      luc25->SetLineColor(kRed+2);
      TLine *lls23=new TLine(ls23,vz2c->GetUymin(),ls23,iivz_hist->GetMaximum());
      lls23->SetLineColor(kGreen);
      lls23->SetLineWidth(3);
      lls23->SetLineStyle(10);
      TLine *lus23=new TLine(us23,vz2c->GetUymin(),us23,iivz_hist->GetMaximum());
      lus23->SetLineColor(kGreen);
      lus23->SetLineWidth(3);
      lus23->SetLineStyle(10);
      TLine *lls24=new TLine(ls24,vz2c->GetUymin(),ls24,iivz_hist->GetMaximum());
      lls24->SetLineColor(kMagenta);
      lls24->SetLineWidth(3);
      lls24->SetLineStyle(10);
      TLine *lus24=new TLine(us24,vz2c->GetUymin(),us24,iivz_hist->GetMaximum());
      lus24->SetLineColor(kMagenta);
      lus24->SetLineWidth(3);
      lus24->SetLineStyle(10);
      TLine *lls25=new TLine(ls25,vz2c->GetUymin(),ls25,iivz_hist->GetMaximum());
      lls25->SetLineColor(kGreen+2);
      lls25->SetLineWidth(3);
      lls25->SetLineStyle(10);
      TLine *lus25=new TLine(us25,vz2c->GetUymin(),us25,iivz_hist->GetMaximum());
      lus25->SetLineColor(kGreen+2);
      lus25->SetLineWidth(3);
      lus25->SetLineStyle(10);
      iivz_hist->Draw("HIST");
   lls23->Draw("HIST Same");
   lus23->Draw("HIST Same");
   lls24->Draw("HIST Same");
   lus24->Draw("HIST Same");
   //lls25->Draw("HIST Same");
   //lus25->Draw("HIST Same");
   llc23->Draw("HIST Same");
   luc23->Draw("HIST Same");
   llc24->Draw("HIST Same");
   luc24->Draw("HIST Same");
   //llc25->Draw("HIST Same");
   //luc25->Draw("HIST Same");
   llclim3->Draw("HIST Same");
   luclim3->Draw("HIST Same");
   llclim4->Draw("HIST Same");
   luclim4->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   llslim3->Draw("HIST Same");
   luslim3->Draw("HIST Same");
   llslim4->Draw("HIST Same");
   luslim4->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vz2c->Write();
   




    TCanvas *vx3c = new TCanvas("vx3c","vx3",0,0,640,480);
      TLine *llx33=new TLine(lx33,vx3c->GetUymin(),lx33,iiivx_hist->GetMaximum());
      llx33->SetLineColor(kRed);
      TLine *lux33=new TLine(ux33,vx3c->GetUymin(),ux33,iiivx_hist->GetMaximum());
      lux33->SetLineColor(kRed);
      TLine *llx34=new TLine(lx34,vx3c->GetUymin(),lx34,iiivx_hist->GetMaximum());
      llx34->SetLineColor(kBlack);
      TLine *lux34=new TLine(ux34,vx3c->GetUymin(),ux34,iiivx_hist->GetMaximum());
      lux34->SetLineColor(kBlack);
      TLine *llx35=new TLine(lx35,vx3c->GetUymin(),lx35,iiivx_hist->GetMaximum());
      llx35->SetLineColor(kRed+2);
      TLine *lux35=new TLine(ux35,vx3c->GetUymin(),ux35,iiivx_hist->GetMaximum());
      lux35->SetLineColor(kRed+2);
   iiivx_hist->GetXaxis()->SetTitle("vx [cm]");
   iiivx_hist->Draw("HIST");
   llx33->Draw("HIST Same");
   lux33->Draw("HIST Same");
   llx34->Draw("HIST Same");
   lux34->Draw("HIST Same");
   llx35->Draw("HIST Same");
   lux35->Draw("HIST Same");
   llxlim3->Draw("HIST Same");
   luxlim3->Draw("HIST Same");
   llxlim4->Draw("HIST Same");
   luxlim4->Draw("HIST Same");
   llxlim5->Draw("HIST Same");
   luxlim5->Draw("HIST Same");
   vx3c->Write();

TCanvas *vy3c = new TCanvas("vy3c","vy3",0,0,640,480);
      TLine *lly33=new TLine(ly33,vy3c->GetUymin(),ly33,iiivy_hist->GetMaximum());
      lly33->SetLineColor(kRed);
      TLine *luy33=new TLine(uy33,vy3c->GetUymin(),uy33,iiivy_hist->GetMaximum());
      luy33->SetLineColor(kRed);
      TLine *lly34=new TLine(ly34,vy3c->GetUymin(),ly34,iiivy_hist->GetMaximum());
      lly34->SetLineColor(kBlack);
      TLine *luy34=new TLine(uy34,vy3c->GetUymin(),uy34,iiivy_hist->GetMaximum());
      luy34->SetLineColor(kBlack);
      TLine *lly35=new TLine(ly35,vy3c->GetUymin(),ly35,iiivy_hist->GetMaximum());
      lly35->SetLineColor(kRed+2);
      TLine *luy35=new TLine(uy35,vy3c->GetUymin(),uy35,iiivy_hist->GetMaximum());
      luy35->SetLineColor(kRed+2);
   iiivy_hist->Draw("HIST");
   lly33->Draw("HIST Same");
   luy33->Draw("HIST Same");
   lly34->Draw("HIST Same");
   luy34->Draw("HIST Same");
   lly35->Draw("HIST Same");
   luy35->Draw("HIST Same");
   llylim3->Draw("HIST Same");
   luylim3->Draw("HIST Same");
   llylim4->Draw("HIST Same");
   luylim4->Draw("HIST Same");
   llylim5->Draw("HIST Same");
   luylim5->Draw("HIST Same");
   vy3c->Write();

   

TCanvas *vz3c = new TCanvas("vz3c","vz3",0,0,640,480);
      TLine *llc33=new TLine(lc33,vz3c->GetUymin(),lc33,iiivz_hist->GetMaximum());
      llc33->SetLineColor(kRed);
      TLine *luc33=new TLine(uc33,vz3c->GetUymin(),uc33,iiivz_hist->GetMaximum());
      luc33->SetLineColor(kRed);
      TLine *llc34=new TLine(lc34,vz3c->GetUymin(),lc34,iiivz_hist->GetMaximum());
      llc34->SetLineColor(kBlack);
      TLine *luc34=new TLine(uc34,vz3c->GetUymin(),uc34,iiivz_hist->GetMaximum());
      luc34->SetLineColor(kBlack);
      TLine *llc35=new TLine(lc35,vz3c->GetUymin(),lc35,iiivz_hist->GetMaximum());
      llc35->SetLineColor(kRed+2);
      TLine *luc35=new TLine(uc35,vz3c->GetUymin(),uc35,iiivz_hist->GetMaximum());
      luc35->SetLineColor(kRed+2);
      TLine *lls33=new TLine(ls33,vz3c->GetUymin(),ls33,iiivz_hist->GetMaximum());
      lls33->SetLineColor(kGreen);
      lls33->SetLineWidth(3);
      lls33->SetLineStyle(10);
      TLine *lus33=new TLine(us33,vz3c->GetUymin(),us33,iiivz_hist->GetMaximum());
      lus33->SetLineColor(kGreen);
      lus33->SetLineWidth(3);
      lus33->SetLineStyle(10);
      TLine *lls34=new TLine(ls34,vz3c->GetUymin(),ls34,iiivz_hist->GetMaximum());
      lls34->SetLineColor(kMagenta);
      lls34->SetLineWidth(3);
      lls34->SetLineStyle(10);
      TLine *lus34=new TLine(us34,vz3c->GetUymin(),us34,iiivz_hist->GetMaximum());
      lus34->SetLineColor(kMagenta);
      lus34->SetLineWidth(3);
      lus34->SetLineStyle(10);
      TLine *lls35=new TLine(ls35,vz3c->GetUymin(),ls35,iiivz_hist->GetMaximum());
      lls35->SetLineColor(kGreen+2);
      lls35->SetLineWidth(3);
      lls35->SetLineStyle(10);
      TLine *lus35=new TLine(us35,vz3c->GetUymin(),us35,iiivz_hist->GetMaximum());
      lus35->SetLineColor(kGreen+2);
      lus35->SetLineWidth(3);
      lus35->SetLineStyle(10);
      iiivz_hist->Draw("HIST");
   lls33->Draw("HIST Same");
   lus33->Draw("HIST Same");
   lls34->Draw("HIST Same");
   lus34->Draw("HIST Same");
   //lls35->Draw("HIST Same");
   //lus35->Draw("HIST Same");
   llc33->Draw("HIST Same");
   luc33->Draw("HIST Same");
   llc34->Draw("HIST Same");
   luc34->Draw("HIST Same");
   //llc35->Draw("HIST Same");
   //luc35->Draw("HIST Same");
   llclim3->Draw("HIST Same");
   luclim3->Draw("HIST Same");
   llclim4->Draw("HIST Same");
   luclim4->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   llslim3->Draw("HIST Same");
   luslim3->Draw("HIST Same");
   llslim4->Draw("HIST Same");
   luslim4->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vz3c->Write();





TCanvas *vx4c = new TCanvas("vx4c","vx4",0,0,640,480);
      TLine *llx43=new TLine(lx43,vx4c->GetUymin(),lx43,ivvx_hist->GetMaximum());
      llx43->SetLineColor(kRed);
      TLine *lux43=new TLine(ux43,vx4c->GetUymin(),ux43,ivvx_hist->GetMaximum());
      lux43->SetLineColor(kRed);
      TLine *llx44=new TLine(lx44,vx4c->GetUymin(),lx44,ivvx_hist->GetMaximum());
      llx44->SetLineColor(kBlack);
      TLine *lux44=new TLine(ux44,vx4c->GetUymin(),ux44,ivvx_hist->GetMaximum());
      lux44->SetLineColor(kBlack);
      TLine *llx45=new TLine(lx45,vx4c->GetUymin(),lx45,ivvx_hist->GetMaximum());
      llx45->SetLineColor(kRed+2);
      TLine *lux45=new TLine(ux45,vx4c->GetUymin(),ux45,ivvx_hist->GetMaximum());
      lux45->SetLineColor(kRed+2);
   ivvx_hist->Draw("HIST");
   llx43->Draw("HIST Same");
   lux43->Draw("HIST Same");
   llx44->Draw("HIST Same");
   lux44->Draw("HIST Same");
   llx45->Draw("HIST Same");
   lux45->Draw("HIST Same");
   llxlim3->Draw("HIST Same");
   luxlim3->Draw("HIST Same");
   llxlim4->Draw("HIST Same");
   luxlim4->Draw("HIST Same");
   llxlim5->Draw("HIST Same");
   luxlim5->Draw("HIST Same");
   vx4c->Write();

TCanvas *vy4c = new TCanvas("vy4c","vy4",0,0,640,480);
      TLine *lly43=new TLine(ly43,vy4c->GetUymin(),ly43,ivvy_hist->GetMaximum());
      lly43->SetLineColor(kRed);
      TLine *luy43=new TLine(uy43,vy4c->GetUymin(),uy43,ivvy_hist->GetMaximum());
      luy43->SetLineColor(kRed);
      TLine *lly44=new TLine(ly44,vy4c->GetUymin(),ly44,ivvy_hist->GetMaximum());
      lly44->SetLineColor(kBlack);
      TLine *luy44=new TLine(uy44,vy4c->GetUymin(),uy44,ivvy_hist->GetMaximum());
      luy44->SetLineColor(kBlack);
      TLine *lly45=new TLine(ly45,vy4c->GetUymin(),ly45,ivvy_hist->GetMaximum());
      lly45->SetLineColor(kRed+2);
      TLine *luy45=new TLine(uy45,vy4c->GetUymin(),uy45,ivvy_hist->GetMaximum());
      luy45->SetLineColor(kRed+2);
   ivvy_hist->Draw("HIST");
   lly43->Draw("HIST Same");
   luy43->Draw("HIST Same");
   lly44->Draw("HIST Same");
   luy44->Draw("HIST Same");
   lly45->Draw("HIST Same");
   luy45->Draw("HIST Same");
   llylim3->Draw("HIST Same");
   luylim3->Draw("HIST Same");
   llylim4->Draw("HIST Same");
   luylim4->Draw("HIST Same");
   llylim5->Draw("HIST Same");
   luylim5->Draw("HIST Same");
   vy4c->Write();
  

   

TCanvas *vz4c = new TCanvas("vz4c","vz4",0,0,640,480);
      TLine *llvx1c3=new TLine(lvx1c3,vz4c->GetUymin(),lvx1c3,ivvz_hist->GetMaximum());
      llvx1c3->SetLineColor(kRed);
      TLine *luvx1c3=new TLine(uvx1c3,vz4c->GetUymin(),uvx1c3,ivvz_hist->GetMaximum());
      luvx1c3->SetLineColor(kRed);
      TLine *llvx1c4=new TLine(lvx1c4,vz4c->GetUymin(),lvx1c4,ivvz_hist->GetMaximum());
      llvx1c4->SetLineColor(kBlack);
      TLine *luvx1c4=new TLine(uvx1c4,vz4c->GetUymin(),uvx1c4,ivvz_hist->GetMaximum());
      luvx1c4->SetLineColor(kBlack);
      TLine *llvx1vy1c=new TLine(lvx1vy1c,vz4c->GetUymin(),lvx1vy1c,ivvz_hist->GetMaximum());
      llvx1vy1c->SetLineColor(kRed+2);
      TLine *luvx1vy1c=new TLine(uvx1vy1c,vz4c->GetUymin(),uvx1vy1c,ivvz_hist->GetMaximum());
      luvx1vy1c->SetLineColor(kRed+2);
      TLine *lls43=new TLine(ls43,vz4c->GetUymin(),ls43,ivvz_hist->GetMaximum());
      lls43->SetLineColor(kGreen);
      lls43->SetLineWidth(3);
      lls43->SetLineStyle(10);
      TLine *lus43=new TLine(us43,vz4c->GetUymin(),us43,ivvz_hist->GetMaximum());
      lus43->SetLineColor(kGreen);
      lus43->SetLineWidth(3);
      lus43->SetLineStyle(10);
      TLine *lls44=new TLine(ls44,vz4c->GetUymin(),ls44,ivvz_hist->GetMaximum());
      lls44->SetLineColor(kMagenta);
      lls44->SetLineWidth(3);
      lls44->SetLineStyle(10);
      TLine *lus44=new TLine(us44,vz4c->GetUymin(),us44,ivvz_hist->GetMaximum());
      lus44->SetLineColor(kMagenta);
      lus44->SetLineWidth(3);
      lus44->SetLineStyle(10);
      TLine *lls45=new TLine(ls45,vz4c->GetUymin(),ls45,ivvz_hist->GetMaximum());
      lls45->SetLineColor(kGreen+2);
      lls45->SetLineWidth(3);
      lls45->SetLineStyle(10);
      TLine *lus45=new TLine(us45,vz4c->GetUymin(),us45,ivvz_hist->GetMaximum());
      lus45->SetLineColor(kGreen+2);
      lus45->SetLineWidth(3);
      lus45->SetLineStyle(10);
      ivvz_hist->Draw("HIST");
   lls43->Draw("HIST Same");
   lus43->Draw("HIST Same");
   lls44->Draw("HIST Same");
   lus44->Draw("HIST Same");
   //lls45->Draw("HIST Same");
   //lus45->Draw("HIST Same");
   llvx1c3->Draw("HIST Same");
   luvx1c3->Draw("HIST Same");
   llvx1c4->Draw("HIST Same");
   luvx1c4->Draw("HIST Same");
   //llvx1vy1c->Draw("HIST Same");
   //luvx1vy1c->Draw("HIST Same");
   llclim3->Draw("HIST Same");
   luclim3->Draw("HIST Same");
   llclim4->Draw("HIST Same");
   luclim4->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   llslim3->Draw("HIST Same");
   luslim3->Draw("HIST Same");
   llslim4->Draw("HIST Same");
   luslim4->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vz4c->Write();
   




TCanvas *vx5c = new TCanvas("vx5c","vx5",0,0,640,480);
      TLine *llx53=new TLine(lx53,vx5c->GetUymin(),lx53,vvx_hist->GetMaximum());
      llx53->SetLineColor(kRed);
      TLine *lux53=new TLine(ux53,vx5c->GetUymin(),ux53,vvx_hist->GetMaximum());
      lux53->SetLineColor(kRed);
      TLine *llx54=new TLine(lx54,vx5c->GetUymin(),lx54,vvx_hist->GetMaximum());
      llx54->SetLineColor(kBlack);
      TLine *lux54=new TLine(ux54,vx5c->GetUymin(),ux54,vvx_hist->GetMaximum());
      lux54->SetLineColor(kBlack);
      TLine *llx55=new TLine(lx55,vx5c->GetUymin(),lx55,vvx_hist->GetMaximum());
      llx55->SetLineColor(kRed+2);
      TLine *lux55=new TLine(ux55,vx5c->GetUymin(),ux55,vvx_hist->GetMaximum());
      lux55->SetLineColor(kRed+2);
   vvx_hist->GetXaxis()->SetTitle("vx [cm]");
   vvx_hist->Draw("HIST");
   llx53->Draw("HIST Same");
   lux53->Draw("HIST Same");
   llx54->Draw("HIST Same");
   lux54->Draw("HIST Same");
   llx55->Draw("HIST Same");
   lux55->Draw("HIST Same");
   llxlim3->Draw("HIST Same");
   luxlim3->Draw("HIST Same");
   llxlim4->Draw("HIST Same");
   luxlim4->Draw("HIST Same");
   llxlim5->Draw("HIST Same");
   luxlim5->Draw("HIST Same");
   vx5c->Write();
  

TCanvas *vy5c = new TCanvas("vy5c","vy5",0,0,640,480);
      TLine *lly53=new TLine(ly53,vy5c->GetUymin(),ly53,vvy_hist->GetMaximum());
      lly53->SetLineColor(kRed);
      TLine *luy53=new TLine(uy53,vy5c->GetUymin(),uy53,vvy_hist->GetMaximum());
      luy53->SetLineColor(kRed);
      TLine *lly54=new TLine(ly54,vy5c->GetUymin(),ly54,vvy_hist->GetMaximum());
      lly54->SetLineColor(kBlack);
      TLine *luy54=new TLine(uy54,vy5c->GetUymin(),uy54,vvy_hist->GetMaximum());
      luy54->SetLineColor(kBlack);
      TLine *lly55=new TLine(ly55,vy5c->GetUymin(),ly55,vvy_hist->GetMaximum());
      lly55->SetLineColor(kRed+2);
      TLine *luy55=new TLine(uy55,vy5c->GetUymin(),uy55,vvy_hist->GetMaximum());
      luy55->SetLineColor(kRed+2);
   vvy_hist->GetXaxis()->SetTitle("vy [cm]");
   vvy_hist->Draw("HIST");
   lly53->Draw("HIST Same");
   luy53->Draw("HIST Same");
   lly54->Draw("HIST Same");
   luy54->Draw("HIST Same");
   lly55->Draw("HIST Same");
   luy55->Draw("HIST Same");
   llylim3->Draw("HIST Same");
   luylim3->Draw("HIST Same");
   llylim4->Draw("HIST Same");
   luylim4->Draw("HIST Same");
   llylim5->Draw("HIST Same");
   luylim5->Draw("HIST Same");
   vy5c->Write();
   


TCanvas *vz5c = new TCanvas("vz5c","vz5",0,0,640,480);
      TLine *llvy1c3=new TLine(lvy1c3,vz5c->GetUymin(),lvy1c3,vvz_hist->GetMaximum());
      llvy1c3->SetLineColor(kRed);
      TLine *luvy1c3=new TLine(uvy1c3,vz5c->GetUymin(),uvy1c3,vvz_hist->GetMaximum());
      luvy1c3->SetLineColor(kRed);
      TLine *llvy1c4=new TLine(lvy1c4,vz5c->GetUymin(),lvy1c4,vvz_hist->GetMaximum());
      llvy1c4->SetLineColor(kBlack);
      TLine *luvy1c4=new TLine(uvy1c4,vz5c->GetUymin(),uvy1c4,vvz_hist->GetMaximum());
      luvy1c4->SetLineColor(kBlack);
      TLine *llvy1c5=new TLine(lvy1c5,vz5c->GetUymin(),lvy1c5,vvz_hist->GetMaximum());
      llvy1c5->SetLineColor(kRed+2);
      TLine *luvy1c5=new TLine(uvy1c5,vz5c->GetUymin(),uvy1c5,vvz_hist->GetMaximum());
      luvy1c5->SetLineColor(kRed+2);
      TLine *lls53=new TLine(ls53,vz5c->GetUymin(),ls53,vvz_hist->GetMaximum());
      lls53->SetLineColor(kGreen);
      lls53->SetLineWidth(3);
      lls53->SetLineStyle(10);
      TLine *lus53=new TLine(us53,vz5c->GetUymin(),us53,vvz_hist->GetMaximum());
      lus53->SetLineColor(kGreen);
      lus53->SetLineWidth(3);
      lus53->SetLineStyle(10);
      TLine *lls54=new TLine(ls54,vz5c->GetUymin(),ls54,vvz_hist->GetMaximum());
      lls54->SetLineColor(kMagenta);
      lls54->SetLineWidth(3);
      lls54->SetLineStyle(10);
      TLine *lus54=new TLine(us54,vz5c->GetUymin(),us54,vvz_hist->GetMaximum());
      lus54->SetLineColor(kMagenta);
      lus54->SetLineWidth(3);
      lus54->SetLineStyle(10);
      TLine *lls55=new TLine(ls55,vz5c->GetUymin(),ls55,vvz_hist->GetMaximum());
      lls55->SetLineColor(kGreen+2);
      lls55->SetLineWidth(3);
      lls55->SetLineStyle(10);
      TLine *lus55=new TLine(us55,vz5c->GetUymin(),us55,vvz_hist->GetMaximum());
      lus55->SetLineColor(kGreen+2);
      lus55->SetLineWidth(3);
      lus55->SetLineStyle(10);
      vvz_hist->Draw("HIST");
   lls53->Draw("HIST Same");
   lus53->Draw("HIST Same");
   lls54->Draw("HIST Same");
   lus54->Draw("HIST Same");
   //lls55->Draw("HIST Same");
   //lus55->Draw("HIST Same");
   llvy1c3->Draw("HIST Same");
   luvy1c3->Draw("HIST Same");
   llvy1c4->Draw("HIST Same");
   luvy1c4->Draw("HIST Same");
   //llvy1c5->Draw("HIST Same");
   //luvy1c5->Draw("HIST Same");
   llclim3->Draw("HIST Same");
   luclim3->Draw("HIST Same");
   llclim4->Draw("HIST Same");
   luclim4->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   llslim3->Draw("HIST Same");
   luslim3->Draw("HIST Same");
   llslim4->Draw("HIST Same");
   luslim4->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vz5c->Write();
  

 


TCanvas *vx6c = new TCanvas("vx6c","vx6",0,0,640,480);
      TLine *llx63=new TLine(lx63,vx6c->GetUymin(),lx63,vivx_hist->GetMaximum());
      llx63->SetLineColor(kRed);
      TLine *lux63=new TLine(ux63,vx6c->GetUymin(),ux63,vivx_hist->GetMaximum());
      lux63->SetLineColor(kRed);
      TLine *llx64=new TLine(lx64,vx6c->GetUymin(),lx64,vivx_hist->GetMaximum());
      llx64->SetLineColor(kBlack);
      TLine *lux64=new TLine(ux64,vx6c->GetUymin(),ux64,vivx_hist->GetMaximum());
      lux64->SetLineColor(kBlack);
      TLine *llx65=new TLine(lx65,vx6c->GetUymin(),lx65,vivx_hist->GetMaximum());
      llx65->SetLineColor(kRed+2);
      TLine *lux65=new TLine(ux65,vx6c->GetUymin(),ux65,vivx_hist->GetMaximum());
      lux65->SetLineColor(kRed+2);
   vivx_hist->GetXaxis()->SetTitle("vx [cm]");
   vivx_hist->Draw("HIST");
   llx63->Draw("HIST Same");
   lux63->Draw("HIST Same");
   llx64->Draw("HIST Same");
   lux64->Draw("HIST Same");
   llx65->Draw("HIST Same");
   lux65->Draw("HIST Same");
   llxlim3->Draw("HIST Same");
   luxlim3->Draw("HIST Same");
   llxlim4->Draw("HIST Same");
   luxlim4->Draw("HIST Same");
   llxlim5->Draw("HIST Same");
   luxlim5->Draw("HIST Same");
   vx6c->Write();
   

TCanvas *vy6c = new TCanvas("vy6c","vy6",0,0,640,480);
      TLine *lly63=new TLine(ly63,vy6c->GetUymin(),ly63,vivy_hist->GetMaximum());
      lly63->SetLineColor(kRed);
      TLine *luy63=new TLine(uy63,vy6c->GetUymin(),uy63,vivy_hist->GetMaximum());
      luy63->SetLineColor(kRed);
      TLine *lly64=new TLine(ly64,vy6c->GetUymin(),ly64,vivy_hist->GetMaximum());
      lly64->SetLineColor(kBlack);
      TLine *luy64=new TLine(uy64,vy6c->GetUymin(),uy64,vivy_hist->GetMaximum());
      luy64->SetLineColor(kBlack);
      TLine *lly65=new TLine(ly65,vy6c->GetUymin(),ly65,vivy_hist->GetMaximum());
      lly65->SetLineColor(kRed+2);
      TLine *luy65=new TLine(uy65,vy6c->GetUymin(),uy65,vivy_hist->GetMaximum());
      luy65->SetLineColor(kRed+2);
   vivy_hist->GetXaxis()->SetTitle("vy [cm]");
   vivy_hist->Draw("HIST");
   lly63->Draw("HIST Same");
   luy63->Draw("HIST Same");
   lly64->Draw("HIST Same");
   luy64->Draw("HIST Same");
   lly65->Draw("HIST Same");
   luy65->Draw("HIST Same");
   llylim3->Draw("HIST Same");
   luylim3->Draw("HIST Same");
   llylim4->Draw("HIST Same");
   luylim4->Draw("HIST Same");
   llylim5->Draw("HIST Same");
   luylim5->Draw("HIST Same");
   vy6c->Write();
   
   
   
   TCanvas *vz6c = new TCanvas("vz6c","vz6",0,0,640,480);
      TLine *llvx2c3=new TLine(lvx2c3,vz6c->GetUymin(),lvx2c3,vivz_hist->GetMaximum());
      llvx2c3->SetLineColor(kRed);
      TLine *luvx2c3=new TLine(uvx2c3,vz6c->GetUymin(),uvx2c3,vivz_hist->GetMaximum());
      luvx2c3->SetLineColor(kRed);
      TLine *llvx2c4=new TLine(lvx2c4,vz6c->GetUymin(),lvx2c4,vivz_hist->GetMaximum());
      llvx2c4->SetLineColor(kBlack);
      TLine *luvx2c4=new TLine(uvx2c4,vz6c->GetUymin(),uvx2c4,vivz_hist->GetMaximum());
      luvx2c4->SetLineColor(kBlack);
      TLine *llvx2c5=new TLine(lvx2c5,vz6c->GetUymin(),lvx2c5,vivz_hist->GetMaximum());
      llvx2c5->SetLineColor(kRed+2);
      TLine *luvx2c5=new TLine(uvx2c5,vz6c->GetUymin(),uvx2c5,vivz_hist->GetMaximum());
      luvx2c5->SetLineColor(kRed+2);
      TLine *lls63=new TLine(ls63,vz6c->GetUymin(),ls63,vivz_hist->GetMaximum());
      lls63->SetLineColor(kGreen);
      lls63->SetLineWidth(3);
      lls63->SetLineStyle(10);
      TLine *lus63=new TLine(us63,vz6c->GetUymin(),us63,vivz_hist->GetMaximum());
      lus63->SetLineColor(kGreen);
      lus63->SetLineWidth(3);
      lus63->SetLineStyle(10);
      TLine *lls64=new TLine(ls64,vz6c->GetUymin(),ls64,vivz_hist->GetMaximum());
      lls64->SetLineColor(kMagenta);
      lls64->SetLineWidth(3);
      lls64->SetLineStyle(10);
      TLine *lus64=new TLine(us64,vz6c->GetUymin(),us64,vivz_hist->GetMaximum());
      lus64->SetLineColor(kMagenta);
      lus64->SetLineWidth(3);
      lus64->SetLineStyle(10);
      TLine *lls65=new TLine(ls65,vz6c->GetUymin(),ls65,vivz_hist->GetMaximum());
      lls65->SetLineColor(kGreen+2);
      lls65->SetLineWidth(3);
      lls65->SetLineStyle(10);
      TLine *lus65=new TLine(us65,vz6c->GetUymin(),us65,vivz_hist->GetMaximum());
      lus65->SetLineColor(kGreen+2);
      lus65->SetLineWidth(3);
      lus65->SetLineStyle(10);
      vivz_hist->Draw("HIST");
   lls63->Draw("HIST Same");
   lus63->Draw("HIST Same");
   lls64->Draw("HIST Same");
   lus64->Draw("HIST Same");
   //lls65->Draw("HIST Same");
   //lus65->Draw("HIST Same");
   llvx2c3->Draw("HIST Same");
   luvx2c3->Draw("HIST Same");
   llvx2c4->Draw("HIST Same");
   luvx2c4->Draw("HIST Same");
   //llvx2c5->Draw("HIST Same");
   //luvx2c5->Draw("HIST Same");
   llclim3->Draw("HIST Same");
   luclim3->Draw("HIST Same");
   llclim4->Draw("HIST Same");
   luclim4->Draw("HIST Same");
   //llclim5->Draw("HIST Same");
   //luclim5->Draw("HIST Same");
   llslim3->Draw("HIST Same");
   luslim3->Draw("HIST Same");
   llslim4->Draw("HIST Same");
   luslim4->Draw("HIST Same");
   //llslim5->Draw("HIST Same");
   //luslim5->Draw("HIST Same");
   vz6c->Write();

  
*/
   
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
   std::cout << "enter: " << vz_fit->GetParameter(5) << " : " << vz_fit->GetParameter(6) <<   std::endl;
   std::cout << "exit: "  << vz_fit->GetParameter(8) << " : " << vz_fit->GetParameter(9) <<  std::endl;
   std::cout << "Cu: " << vz2_fit->GetParameter(1) << " : " << vz2_fit->GetParameter(2) <<  std::endl;
   std::cout << "Sn: " << vz2_fit->GetParameter(4) << " : " << vz2_fit->GetParameter(5) <<  std::endl;
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
   return 0;
}



