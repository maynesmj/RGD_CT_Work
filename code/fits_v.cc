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

void fits_v(){


   TFile input_1("cusn_18624_hist_2-2.root", "read");
   TFile input_2("cusn_18625_hist_2-2.root", "read");
   TFile input_3("cusn_18626_hist_2-2.root", "read");
   TFile input_4("cusn_18627_hist_2-2.root", "read");
   TFile input_5("cusn_18628_hist_2-2.root", "read");
   TFile input_6("cusn_18629_hist_2-2.root", "read");
   TFile input_7("cusn_18630_hist_2-2.root", "read");
   TFile input_8("cusn_18631_hist_2-2.root", "read");
   TFile input_9("cusn_18632_hist_2-2.root", "read");
   TFile input_10("cusn_18636_hist_2-2.root", "read");

   
    TFile f("cusnfits_cache_sectors_02-2.root", "update");

       TH1D *tvx1_hist = new TH1D("sector 1 vx", "sector 1 vx", 300, -1, 1);
	tvx1_hist->SetStats(0);
	
	TH1D *tvy1_hist = new TH1D("sector 1 vy", "sector 1 vy", 300, -1, 1);
	tvy1_hist->SetStats(0);
	
	TH1D *tvz1_hist = new TH1D("sector 1 vz", "sector 1 vz", 300, -20, 20);
	tvz1_hist->SetStats(0);

      TH1D *vx1_1_hist = (TH1D*)input_1.Get("sector 1 vx");
      TH1D *vy1_1_hist = (TH1D*)input_1.Get("sector 1 vy");
      TH1D *vz1_1_hist = (TH1D*)input_1.Get("sector 1 vz");
      TH1D *vx1_2_hist = (TH1D*)input_2.Get("sector 1 vx");
      TH1D *vy1_2_hist = (TH1D*)input_2.Get("sector 1 vy");
      TH1D *vz1_2_hist = (TH1D*)input_2.Get("sector 1 vz");
      TH1D *vx1_3_hist = (TH1D*)input_3.Get("sector 1 vx");
      TH1D *vy1_3_hist = (TH1D*)input_3.Get("sector 1 vy");
      TH1D *vz1_3_hist = (TH1D*)input_3.Get("sector 1 vz");
      TH1D *vx1_4_hist = (TH1D*)input_4.Get("sector 1 vx");
      TH1D *vy1_4_hist = (TH1D*)input_4.Get("sector 1 vy");
      TH1D *vz1_4_hist = (TH1D*)input_4.Get("sector 1 vz");
      TH1D *vx1_5_hist = (TH1D*)input_5.Get("sector 1 vx");
      TH1D *vy1_5_hist = (TH1D*)input_5.Get("sector 1 vy");
      TH1D *vz1_5_hist = (TH1D*)input_5.Get("sector 1 vz");
      TH1D *vx1_6_hist = (TH1D*)input_6.Get("sector 1 vx");
      TH1D *vy1_6_hist = (TH1D*)input_6.Get("sector 1 vy");
      TH1D *vz1_6_hist = (TH1D*)input_6.Get("sector 1 vz");
      TH1D *vx1_7_hist = (TH1D*)input_7.Get("sector 1 vx");
      TH1D *vy1_7_hist = (TH1D*)input_7.Get("sector 1 vy");
      TH1D *vz1_7_hist = (TH1D*)input_7.Get("sector 1 vz");
      TH1D *vx1_8_hist = (TH1D*)input_8.Get("sector 1 vx");
      TH1D *vy1_8_hist = (TH1D*)input_8.Get("sector 1 vy");
      TH1D *vz1_8_hist = (TH1D*)input_8.Get("sector 1 vz");
      TH1D *vx1_9_hist = (TH1D*)input_9.Get("sector 1 vx");
      TH1D *vy1_9_hist = (TH1D*)input_9.Get("sector 1 vy");
      TH1D *vz1_9_hist = (TH1D*)input_9.Get("sector 1 vz");
      TH1D *vx1_10_hist = (TH1D*)input_10.Get("sector 1 vx");
      TH1D *vy1_10_hist = (TH1D*)input_10.Get("sector 1 vy");
      TH1D *vz1_10_hist = (TH1D*)input_10.Get("sector 1 vz");
      
      
      tvx1_hist->Add(vx1_1_hist);
      tvx1_hist->Add(vx1_2_hist);
      tvx1_hist->Add(vx1_3_hist);
      tvx1_hist->Add(vx1_4_hist);
      tvx1_hist->Add(vx1_5_hist);
      tvx1_hist->Add(vx1_6_hist);
      tvx1_hist->Add(vx1_7_hist);
      tvx1_hist->Add(vx1_8_hist);
      tvx1_hist->Add(vx1_9_hist);
      tvx1_hist->Add(vx1_10_hist);
      tvy1_hist->Add(vy1_1_hist);
      tvy1_hist->Add(vy1_2_hist);
      tvy1_hist->Add(vy1_3_hist);
      tvy1_hist->Add(vy1_4_hist);
      tvy1_hist->Add(vy1_5_hist);
      tvy1_hist->Add(vy1_6_hist);
      tvy1_hist->Add(vy1_7_hist);
      tvy1_hist->Add(vy1_8_hist);
      tvy1_hist->Add(vy1_9_hist);
      tvy1_hist->Add(vy1_10_hist);
      tvz1_hist->Add(vz1_1_hist);
      tvz1_hist->Add(vz1_2_hist);
      tvz1_hist->Add(vz1_3_hist);
      tvz1_hist->Add(vz1_4_hist);
      tvz1_hist->Add(vz1_5_hist);
      tvz1_hist->Add(vz1_6_hist);
      tvz1_hist->Add(vz1_7_hist);
      tvz1_hist->Add(vz1_8_hist);
      tvz1_hist->Add(vz1_9_hist);
      tvz1_hist->Add(vz1_10_hist);
      
        TH1D *tvx2_hist = new TH1D("sector 2 vx", "sector 2 vx", 300, -1, 1);
	tvx2_hist->SetStats(0);
	
	TH1D *tvy2_hist = new TH1D("sector 2 vy", "sector 2 vy", 300, -1, 1);
	tvy2_hist->SetStats(0);
	
	TH1D *tvz2_hist = new TH1D("sector 2 vz", "sector 2 vz", 300, -20, 20);
	tvz2_hist->SetStats(0);

      TH1D *vx2_1_hist = (TH1D*)input_1.Get("sector 2 vx");
      TH1D *vy2_1_hist = (TH1D*)input_1.Get("sector 2 vy");
      TH1D *vz2_1_hist = (TH1D*)input_1.Get("sector 2 vz");
      TH1D *vx2_2_hist = (TH1D*)input_2.Get("sector 2 vx");
      TH1D *vy2_2_hist = (TH1D*)input_2.Get("sector 2 vy");
      TH1D *vz2_2_hist = (TH1D*)input_2.Get("sector 2 vz");
      TH1D *vx2_3_hist = (TH1D*)input_3.Get("sector 2 vx");
      TH1D *vy2_3_hist = (TH1D*)input_3.Get("sector 2 vy");
      TH1D *vz2_3_hist = (TH1D*)input_3.Get("sector 2 vz");
      TH1D *vx2_4_hist = (TH1D*)input_4.Get("sector 2 vx");
      TH1D *vy2_4_hist = (TH1D*)input_4.Get("sector 2 vy");
      TH1D *vz2_4_hist = (TH1D*)input_4.Get("sector 2 vz");
      TH1D *vx2_5_hist = (TH1D*)input_5.Get("sector 2 vx");
      TH1D *vy2_5_hist = (TH1D*)input_5.Get("sector 2 vy");
      TH1D *vz2_5_hist = (TH1D*)input_5.Get("sector 2 vz");
      TH1D *vx2_6_hist = (TH1D*)input_6.Get("sector 2 vx");
      TH1D *vy2_6_hist = (TH1D*)input_6.Get("sector 2 vy");
      TH1D *vz2_6_hist = (TH1D*)input_6.Get("sector 2 vz");
      TH1D *vx2_7_hist = (TH1D*)input_7.Get("sector 2 vx");
      TH1D *vy2_7_hist = (TH1D*)input_7.Get("sector 2 vy");
      TH1D *vz2_7_hist = (TH1D*)input_7.Get("sector 2 vz");
      TH1D *vx2_8_hist = (TH1D*)input_8.Get("sector 2 vx");
      TH1D *vy2_8_hist = (TH1D*)input_8.Get("sector 2 vy");
      TH1D *vz2_8_hist = (TH1D*)input_8.Get("sector 2 vz");
      TH1D *vx2_9_hist = (TH1D*)input_9.Get("sector 2 vx");
      TH1D *vy2_9_hist = (TH1D*)input_9.Get("sector 2 vy");
      TH1D *vz2_9_hist = (TH1D*)input_9.Get("sector 2 vz");
      TH1D *vx2_10_hist = (TH1D*)input_10.Get("sector 2 vx");
      TH1D *vy2_10_hist = (TH1D*)input_10.Get("sector 2 vy");
      TH1D *vz2_10_hist = (TH1D*)input_10.Get("sector 2 vz");

      tvx2_hist->Add(vx2_1_hist);
      tvx2_hist->Add(vx2_2_hist);
      tvx2_hist->Add(vx2_3_hist);
      tvx2_hist->Add(vx2_4_hist);
      tvx2_hist->Add(vx2_5_hist);
      tvx2_hist->Add(vx2_6_hist);
      tvx2_hist->Add(vx2_7_hist);
      tvx2_hist->Add(vx2_8_hist);
      tvx2_hist->Add(vx2_9_hist);
      tvx2_hist->Add(vx2_10_hist);
      tvy2_hist->Add(vy2_1_hist);
      tvy2_hist->Add(vy2_2_hist);
      tvy2_hist->Add(vy2_3_hist);
      tvy2_hist->Add(vy2_4_hist);
      tvy2_hist->Add(vy2_5_hist);
      tvy2_hist->Add(vy2_6_hist);
      tvy2_hist->Add(vy2_7_hist);
      tvy2_hist->Add(vy2_8_hist);
      tvy2_hist->Add(vy2_9_hist);
      tvy2_hist->Add(vy2_10_hist);
      tvz2_hist->Add(vz2_1_hist);
      tvz2_hist->Add(vz2_2_hist);
      tvz2_hist->Add(vz2_3_hist);
      tvz2_hist->Add(vz2_4_hist);
      tvz2_hist->Add(vz2_5_hist);
      tvz2_hist->Add(vz2_6_hist);
      tvz2_hist->Add(vz2_7_hist);
      tvz2_hist->Add(vz2_8_hist);
      tvz2_hist->Add(vz2_9_hist);
      tvz2_hist->Add(vz2_10_hist);
      
      
      
      
      
      TH1D *tvx3_hist = new TH1D("sector 3 vx", "sector 3 vx", 300, -1, 1);
	tvx3_hist->SetStats(0);
	
	TH1D *tvy3_hist = new TH1D("sector 3 vy", "sector 3 vy", 300, -1, 1);
	tvy3_hist->SetStats(0);
	
	TH1D *tvz3_hist = new TH1D("sector 3 vz", "sector 3 vz", 300, -20, 20);
	tvz3_hist->SetStats(0);

      TH1D *vx3_1_hist = (TH1D*)input_1.Get("sector 3 vx");
      TH1D *vy3_1_hist = (TH1D*)input_1.Get("sector 3 vy");
      TH1D *vz3_1_hist = (TH1D*)input_1.Get("sector 3 vz");
      TH1D *vx3_2_hist = (TH1D*)input_2.Get("sector 3 vx");
      TH1D *vy3_2_hist = (TH1D*)input_2.Get("sector 3 vy");
      TH1D *vz3_2_hist = (TH1D*)input_2.Get("sector 3 vz");
      TH1D *vx3_3_hist = (TH1D*)input_3.Get("sector 3 vx");
      TH1D *vy3_3_hist = (TH1D*)input_3.Get("sector 3 vy");
      TH1D *vz3_3_hist = (TH1D*)input_3.Get("sector 3 vz");
      TH1D *vx3_4_hist = (TH1D*)input_4.Get("sector 3 vx");
      TH1D *vy3_4_hist = (TH1D*)input_4.Get("sector 3 vy");
      TH1D *vz3_4_hist = (TH1D*)input_4.Get("sector 3 vz");
      TH1D *vx3_5_hist = (TH1D*)input_5.Get("sector 3 vx");
      TH1D *vy3_5_hist = (TH1D*)input_5.Get("sector 3 vy");
      TH1D *vz3_5_hist = (TH1D*)input_5.Get("sector 3 vz");
      TH1D *vx3_6_hist = (TH1D*)input_6.Get("sector 3 vx");
      TH1D *vy3_6_hist = (TH1D*)input_6.Get("sector 3 vy");
      TH1D *vz3_6_hist = (TH1D*)input_6.Get("sector 3 vz");
      TH1D *vx3_7_hist = (TH1D*)input_7.Get("sector 3 vx");
      TH1D *vy3_7_hist = (TH1D*)input_7.Get("sector 3 vy");
      TH1D *vz3_7_hist = (TH1D*)input_7.Get("sector 3 vz");
      TH1D *vx3_8_hist = (TH1D*)input_8.Get("sector 3 vx");
      TH1D *vy3_8_hist = (TH1D*)input_8.Get("sector 3 vy");
      TH1D *vz3_8_hist = (TH1D*)input_8.Get("sector 3 vz");
      TH1D *vx3_9_hist = (TH1D*)input_9.Get("sector 3 vx");
      TH1D *vy3_9_hist = (TH1D*)input_9.Get("sector 3 vy");
      TH1D *vz3_9_hist = (TH1D*)input_9.Get("sector 3 vz");
      TH1D *vx3_10_hist = (TH1D*)input_10.Get("sector 3 vx");
      TH1D *vy3_10_hist = (TH1D*)input_10.Get("sector 3 vy");
      TH1D *vz3_10_hist = (TH1D*)input_10.Get("sector 3 vz");

      tvx3_hist->Add(vx3_1_hist);
      tvx3_hist->Add(vx3_2_hist);
      tvx3_hist->Add(vx3_3_hist);
      tvx3_hist->Add(vx3_4_hist);
      tvx3_hist->Add(vx3_5_hist);
      tvx3_hist->Add(vx3_6_hist);
      tvx3_hist->Add(vx3_7_hist);
      tvx3_hist->Add(vx3_8_hist);
      tvx3_hist->Add(vx3_9_hist);
      tvx3_hist->Add(vx3_10_hist);
      tvy3_hist->Add(vy3_1_hist);
      tvy3_hist->Add(vy3_2_hist);
      tvy3_hist->Add(vy3_3_hist);
      tvy3_hist->Add(vy3_4_hist);
      tvy3_hist->Add(vy3_5_hist);
      tvy3_hist->Add(vy3_6_hist);
      tvy3_hist->Add(vy3_7_hist);
      tvy3_hist->Add(vy3_8_hist);
      tvy3_hist->Add(vy3_9_hist);
      tvy3_hist->Add(vy3_10_hist);
      tvz3_hist->Add(vz3_1_hist);
      tvz3_hist->Add(vz3_2_hist);
      tvz3_hist->Add(vz3_3_hist);
      tvz3_hist->Add(vz3_4_hist);
      tvz3_hist->Add(vz3_5_hist);
      tvz3_hist->Add(vz3_6_hist);
      tvz3_hist->Add(vz3_7_hist);
      tvz3_hist->Add(vz3_8_hist);
      tvz3_hist->Add(vz3_9_hist);
      tvz3_hist->Add(vz3_10_hist);
      
      TH1D *tvx4_hist = new TH1D("sector 4 vx", "sector 4 vx", 300, -1, 1);
	tvx4_hist->SetStats(0);
	
	TH1D *tvy4_hist = new TH1D("sector 4 vy", "sector 4 vy", 300, -1, 1);
	tvy4_hist->SetStats(0);
	
	TH1D *tvz4_hist = new TH1D("sector 4 vz", "sector 4 vz", 300, -20, 20);
	tvz4_hist->SetStats(0);

      TH1D *vx4_1_hist = (TH1D*)input_1.Get("sector 4 vx");
      TH1D *vy4_1_hist = (TH1D*)input_1.Get("sector 4 vy");
      TH1D *vz4_1_hist = (TH1D*)input_1.Get("sector 4 vz");
      TH1D *vx4_2_hist = (TH1D*)input_2.Get("sector 4 vx");
      TH1D *vy4_2_hist = (TH1D*)input_2.Get("sector 4 vy");
      TH1D *vz4_2_hist = (TH1D*)input_2.Get("sector 4 vz");
      TH1D *vx4_3_hist = (TH1D*)input_3.Get("sector 4 vx");
      TH1D *vy4_3_hist = (TH1D*)input_3.Get("sector 4 vy");
      TH1D *vz4_3_hist = (TH1D*)input_3.Get("sector 4 vz");
      TH1D *vx4_4_hist = (TH1D*)input_4.Get("sector 4 vx");
      TH1D *vy4_4_hist = (TH1D*)input_4.Get("sector 4 vy");
      TH1D *vz4_4_hist = (TH1D*)input_4.Get("sector 4 vz");
      TH1D *vx4_5_hist = (TH1D*)input_5.Get("sector 4 vx");
      TH1D *vy4_5_hist = (TH1D*)input_5.Get("sector 4 vy");
      TH1D *vz4_5_hist = (TH1D*)input_5.Get("sector 4 vz");
      TH1D *vx4_6_hist = (TH1D*)input_6.Get("sector 4 vx");
      TH1D *vy4_6_hist = (TH1D*)input_6.Get("sector 4 vy");
      TH1D *vz4_6_hist = (TH1D*)input_6.Get("sector 4 vz");
      TH1D *vx4_7_hist = (TH1D*)input_7.Get("sector 4 vx");
      TH1D *vy4_7_hist = (TH1D*)input_7.Get("sector 4 vy");
      TH1D *vz4_7_hist = (TH1D*)input_7.Get("sector 4 vz");
      TH1D *vx4_8_hist = (TH1D*)input_8.Get("sector 4 vx");
      TH1D *vy4_8_hist = (TH1D*)input_8.Get("sector 4 vy");
      TH1D *vz4_8_hist = (TH1D*)input_8.Get("sector 4 vz");
      TH1D *vx4_9_hist = (TH1D*)input_9.Get("sector 4 vx");
      TH1D *vy4_9_hist = (TH1D*)input_9.Get("sector 4 vy");
      TH1D *vz4_9_hist = (TH1D*)input_9.Get("sector 4 vz");
      TH1D *vx4_10_hist = (TH1D*)input_10.Get("sector 4 vx");
      TH1D *vy4_10_hist = (TH1D*)input_10.Get("sector 4 vy");
      TH1D *vz4_10_hist = (TH1D*)input_10.Get("sector 4 vz");

      tvx4_hist->Add(vx4_1_hist);
      tvx4_hist->Add(vx4_2_hist);
      tvx4_hist->Add(vx4_3_hist);
      tvx4_hist->Add(vx4_4_hist);
      tvx4_hist->Add(vx4_5_hist);
      tvx4_hist->Add(vx4_6_hist);
      tvx4_hist->Add(vx4_7_hist);
      tvx4_hist->Add(vx4_8_hist);
      tvx4_hist->Add(vx4_9_hist);
      tvx4_hist->Add(vx4_10_hist);
      tvy4_hist->Add(vy4_1_hist);
      tvy4_hist->Add(vy4_2_hist);
      tvy4_hist->Add(vy4_3_hist);
      tvy4_hist->Add(vy4_4_hist);
      tvy4_hist->Add(vy4_5_hist);
      tvy4_hist->Add(vy4_6_hist);
      tvy4_hist->Add(vy4_7_hist);
      tvy4_hist->Add(vy4_8_hist);
      tvy4_hist->Add(vy4_9_hist);
      tvy4_hist->Add(vy4_10_hist);
      tvz4_hist->Add(vz4_1_hist);
      tvz4_hist->Add(vz4_2_hist);
      tvz4_hist->Add(vz4_3_hist);
      tvz4_hist->Add(vz4_4_hist);
      tvz4_hist->Add(vz4_5_hist);
      tvz4_hist->Add(vz4_6_hist);
      tvz4_hist->Add(vz4_7_hist);
      tvz4_hist->Add(vz4_8_hist);
      tvz4_hist->Add(vz4_9_hist);
      tvz4_hist->Add(vz4_10_hist);
      
      
      
      TH1D *tvx5_hist = new TH1D("sector 5 vx", "sector 5 vx", 300, -1, 1);
	tvx5_hist->SetStats(0);
	
	TH1D *tvy5_hist = new TH1D("sector 5 vy", "sector 5 vy", 300, -1, 1);
	tvy5_hist->SetStats(0);
	
	TH1D *tvz5_hist = new TH1D("sector 5 vz", "sector 5 vz", 300, -20, 20);
	tvz5_hist->SetStats(0);

      TH1D *vx5_1_hist = (TH1D*)input_1.Get("sector 5 vx");
      TH1D *vy5_1_hist = (TH1D*)input_1.Get("sector 5 vy");
      TH1D *vz5_1_hist = (TH1D*)input_1.Get("sector 5 vz");
      TH1D *vx5_2_hist = (TH1D*)input_2.Get("sector 5 vx");
      TH1D *vy5_2_hist = (TH1D*)input_2.Get("sector 5 vy");
      TH1D *vz5_2_hist = (TH1D*)input_2.Get("sector 5 vz");
      TH1D *vx5_3_hist = (TH1D*)input_3.Get("sector 5 vx");
      TH1D *vy5_3_hist = (TH1D*)input_3.Get("sector 5 vy");
      TH1D *vz5_3_hist = (TH1D*)input_3.Get("sector 5 vz");
      TH1D *vx5_4_hist = (TH1D*)input_4.Get("sector 5 vx");
      TH1D *vy5_4_hist = (TH1D*)input_4.Get("sector 5 vy");
      TH1D *vz5_4_hist = (TH1D*)input_4.Get("sector 5 vz");
      TH1D *vx5_5_hist = (TH1D*)input_5.Get("sector 5 vx");
      TH1D *vy5_5_hist = (TH1D*)input_5.Get("sector 5 vy");
      TH1D *vz5_5_hist = (TH1D*)input_5.Get("sector 5 vz");
      TH1D *vx5_6_hist = (TH1D*)input_6.Get("sector 5 vx");
      TH1D *vy5_6_hist = (TH1D*)input_6.Get("sector 5 vy");
      TH1D *vz5_6_hist = (TH1D*)input_6.Get("sector 5 vz");
      TH1D *vx5_7_hist = (TH1D*)input_7.Get("sector 5 vx");
      TH1D *vy5_7_hist = (TH1D*)input_7.Get("sector 5 vy");
      TH1D *vz5_7_hist = (TH1D*)input_7.Get("sector 5 vz");
      TH1D *vx5_8_hist = (TH1D*)input_8.Get("sector 5 vx");
      TH1D *vy5_8_hist = (TH1D*)input_8.Get("sector 5 vy");
      TH1D *vz5_8_hist = (TH1D*)input_8.Get("sector 5 vz");
      TH1D *vx5_9_hist = (TH1D*)input_9.Get("sector 5 vx");
      TH1D *vy5_9_hist = (TH1D*)input_9.Get("sector 5 vy");
      TH1D *vz5_9_hist = (TH1D*)input_9.Get("sector 5 vz");
      TH1D *vx5_10_hist = (TH1D*)input_10.Get("sector 5 vx");
      TH1D *vy5_10_hist = (TH1D*)input_10.Get("sector 5 vy");
      TH1D *vz5_10_hist = (TH1D*)input_10.Get("sector 5 vz");

      tvx5_hist->Add(vx5_1_hist);
      tvx5_hist->Add(vx5_2_hist);
      tvx5_hist->Add(vx5_3_hist);
      tvx5_hist->Add(vx5_4_hist);
      tvx5_hist->Add(vx5_5_hist);
      tvx5_hist->Add(vx5_6_hist);
      tvx5_hist->Add(vx5_7_hist);
      tvx5_hist->Add(vx5_8_hist);
      tvx5_hist->Add(vx5_9_hist);
      tvx5_hist->Add(vx5_10_hist);
      tvy5_hist->Add(vy5_1_hist);
      tvy5_hist->Add(vy5_2_hist);
      tvy5_hist->Add(vy5_3_hist);
      tvy5_hist->Add(vy5_4_hist);
      tvy5_hist->Add(vy5_5_hist);
      tvy5_hist->Add(vy5_6_hist);
      tvy5_hist->Add(vy5_7_hist);
      tvy5_hist->Add(vy5_8_hist);
      tvy5_hist->Add(vy5_9_hist);
      tvy5_hist->Add(vy5_10_hist);
      tvz5_hist->Add(vz5_1_hist);
      tvz5_hist->Add(vz5_2_hist);
      tvz5_hist->Add(vz5_3_hist);
      tvz5_hist->Add(vz5_4_hist);
      tvz5_hist->Add(vz5_5_hist);
      tvz5_hist->Add(vz5_6_hist);
      tvz5_hist->Add(vz5_7_hist);
      tvz5_hist->Add(vz5_8_hist);
      tvz5_hist->Add(vz5_9_hist);
      tvz5_hist->Add(vz5_10_hist);
      
      
      
      
      TH1D *tvx6_hist = new TH1D("sector 6 vx", "sector 6 vx", 300, -1, 1);
	tvx6_hist->SetStats(0);
	
	TH1D *tvy6_hist = new TH1D("sector 6 vy", "sector 6 vy", 300, -1, 1);
	tvy6_hist->SetStats(0);
	
	TH1D *tvz6_hist = new TH1D("sector 6 vz", "sector 6 vz", 300, -20, 20);
	tvz6_hist->SetStats(0);

      TH1D *vx6_1_hist = (TH1D*)input_1.Get("sector 6 vx");
      TH1D *vy6_1_hist = (TH1D*)input_1.Get("sector 6 vy");
      TH1D *vz6_1_hist = (TH1D*)input_1.Get("sector 6 vz");
      TH1D *vx6_2_hist = (TH1D*)input_2.Get("sector 6 vx");
      TH1D *vy6_2_hist = (TH1D*)input_2.Get("sector 6 vy");
      TH1D *vz6_2_hist = (TH1D*)input_2.Get("sector 6 vz");
      TH1D *vx6_3_hist = (TH1D*)input_3.Get("sector 6 vx");
      TH1D *vy6_3_hist = (TH1D*)input_3.Get("sector 6 vy");
      TH1D *vz6_3_hist = (TH1D*)input_3.Get("sector 6 vz");
      TH1D *vx6_4_hist = (TH1D*)input_4.Get("sector 6 vx");
      TH1D *vy6_4_hist = (TH1D*)input_4.Get("sector 6 vy");
      TH1D *vz6_4_hist = (TH1D*)input_4.Get("sector 6 vz");
      TH1D *vx6_5_hist = (TH1D*)input_5.Get("sector 6 vx");
      TH1D *vy6_5_hist = (TH1D*)input_5.Get("sector 6 vy");
      TH1D *vz6_5_hist = (TH1D*)input_5.Get("sector 6 vz");
      TH1D *vx6_6_hist = (TH1D*)input_6.Get("sector 6 vx");
      TH1D *vy6_6_hist = (TH1D*)input_6.Get("sector 6 vy");
      TH1D *vz6_6_hist = (TH1D*)input_6.Get("sector 6 vz");
      TH1D *vx6_7_hist = (TH1D*)input_7.Get("sector 6 vx");
      TH1D *vy6_7_hist = (TH1D*)input_7.Get("sector 6 vy");
      TH1D *vz6_7_hist = (TH1D*)input_7.Get("sector 6 vz");
      TH1D *vx6_8_hist = (TH1D*)input_8.Get("sector 6 vx");
      TH1D *vy6_8_hist = (TH1D*)input_8.Get("sector 6 vy");
      TH1D *vz6_8_hist = (TH1D*)input_8.Get("sector 6 vz");
      TH1D *vx6_9_hist = (TH1D*)input_9.Get("sector 6 vx");
      TH1D *vy6_9_hist = (TH1D*)input_9.Get("sector 6 vy");
      TH1D *vz6_9_hist = (TH1D*)input_9.Get("sector 6 vz");
      TH1D *vx6_10_hist = (TH1D*)input_10.Get("sector 6 vx");
      TH1D *vy6_10_hist = (TH1D*)input_10.Get("sector 6 vy");
      TH1D *vz6_10_hist = (TH1D*)input_10.Get("sector 6 vz");

      tvx6_hist->Add(vx6_1_hist);
      tvx6_hist->Add(vx6_2_hist);
      tvx6_hist->Add(vx6_3_hist);
      tvx6_hist->Add(vx6_4_hist);
      tvx6_hist->Add(vx6_5_hist);
      tvx6_hist->Add(vx6_6_hist);
      tvx6_hist->Add(vx6_7_hist);
      tvx6_hist->Add(vx6_8_hist);
      tvx6_hist->Add(vx6_9_hist);
      tvx6_hist->Add(vx6_10_hist);
      tvy6_hist->Add(vy6_1_hist);
      tvy6_hist->Add(vy6_2_hist);
      tvy6_hist->Add(vy6_3_hist);
      tvy6_hist->Add(vy6_4_hist);
      tvy6_hist->Add(vy6_5_hist);
      tvy6_hist->Add(vy6_6_hist);
      tvy6_hist->Add(vy6_7_hist);
      tvy6_hist->Add(vy6_8_hist);
      tvy6_hist->Add(vy6_9_hist);
      tvy6_hist->Add(vy6_10_hist);
      tvz6_hist->Add(vz6_1_hist);
      tvz6_hist->Add(vz6_2_hist);
      tvz6_hist->Add(vz6_3_hist);
      tvz6_hist->Add(vz6_4_hist);
      tvz6_hist->Add(vz6_5_hist);
      tvz6_hist->Add(vz6_6_hist);
      tvz6_hist->Add(vz6_7_hist);
      tvz6_hist->Add(vz6_8_hist);
      tvz6_hist->Add(vz6_9_hist);
      tvz6_hist->Add(vz6_10_hist);
      
      
      
      
      
      
      TH1D *tvx_hist = new TH1D("vx", "vx", 300, -1, 1);
	tvx_hist->SetStats(0);
	
	TH1D *tvy_hist = new TH1D("vy", "vy", 300, -1, 1);
	tvy_hist->SetStats(0);
	
	TH1D *tvz_hist = new TH1D("vz", "vz", 300, -20, 20);
	tvz_hist->SetStats(0);


      tvx_hist->Add(tvx1_hist);
      tvx_hist->Add(tvx2_hist);
      tvx_hist->Add(tvx3_hist);
      tvx_hist->Add(tvx4_hist);
      tvx_hist->Add(tvx5_hist);
      tvx_hist->Add(tvx6_hist);
      
      
      tvy_hist->Add(tvy1_hist);
      tvy_hist->Add(tvy2_hist);
      tvy_hist->Add(tvy3_hist);
      tvy_hist->Add(tvy4_hist);
      tvy_hist->Add(tvy5_hist);
      tvy_hist->Add(tvy6_hist);
      
      
      
      tvz_hist->Add(tvz1_hist);
      tvz_hist->Add(tvz2_hist);
      tvz_hist->Add(tvz3_hist);
      tvz_hist->Add(tvz4_hist);
      tvz_hist->Add(tvz5_hist);
      tvz_hist->Add(tvz6_hist);






   //TH2F *z_hist = new TH2F("px vs py", "px vs py", 200, -2, 2, 200, -2, 2);

   //auto cgr = new TGraph();
   //std::cout << "1 th" << std::endl;

   std::vector<std::string> v3_a;
   int num = 0;


   tvx1_hist->GetXaxis()->SetTitle("vx [cm]");
   TF1 *tvx1_func = new TF1("tvx1_fit",fitf,-0.5,0.5,6);
   //ivp_func->SetParameter(1,0.150);
   //ivp_func->SetParameter(2,0.77);
   //ivp_func->SetParLimits(2,0.65,0.89);
   //ivp_func->SetParNames("double ivpConstant","double ivpgamma","double ivpmean","double ivpa1","double ivpa2","double ivpa3", "double ivpa4", "double ivpa5", "double ivpa6");
   tvx1_hist->Fit("tvx1_fit","","",-0.5,0.5);
   tvx1_hist->Draw();
   tvx1_hist->Write();
   
   tvy1_hist->GetXaxis()->SetTitle("vy [cm]");
   tvy1_hist->Draw();
   tvy1_hist->Write();
   
   tvz1_hist->GetXaxis()->SetTitle("vz [cm]");
   tvz1_hist->Draw();
   tvz1_hist->Write();



   tvx2_hist->GetXaxis()->SetTitle("vx [cm]");
   tvx2_hist->Draw();
   tvx2_hist->Write();
   
   tvy2_hist->GetXaxis()->SetTitle("vy [cm]");
   tvy2_hist->Draw();
   tvy2_hist->Write();
   
   tvz2_hist->GetXaxis()->SetTitle("vz [cm]");
   tvz2_hist->Draw();
   tvz2_hist->Write();
   
   
   
   tvx3_hist->GetXaxis()->SetTitle("vx [cm]");
   tvx3_hist->Draw();
   tvx3_hist->Write();
   
   tvy3_hist->GetXaxis()->SetTitle("vy [cm]");
   tvy3_hist->Draw();
   tvy3_hist->Write();
   
   tvz3_hist->GetXaxis()->SetTitle("vz [cm]");
   tvz3_hist->Draw();
   tvz3_hist->Write();
   
   
   
   tvx4_hist->GetXaxis()->SetTitle("vx [cm]");
   tvx4_hist->Draw();
   tvx4_hist->Write();
   
   tvy4_hist->GetXaxis()->SetTitle("vy [cm]");
   tvy4_hist->Draw();
   tvy4_hist->Write();
   
   tvz4_hist->GetXaxis()->SetTitle("vz [cm]");
   tvz4_hist->Draw();
   tvz4_hist->Write();
   
   
   
   tvx5_hist->GetXaxis()->SetTitle("vx [cm]");
   tvx5_hist->Draw();
   tvx5_hist->Write();
   
   tvy5_hist->GetXaxis()->SetTitle("vy [cm]");
   tvy5_hist->Draw();
   tvy5_hist->Write();
   
   tvz5_hist->GetXaxis()->SetTitle("vz [cm]");
   tvz5_hist->Draw();
   tvz5_hist->Write();
   
   
   
   
   tvx6_hist->GetXaxis()->SetTitle("vx [cm]");
   tvx6_hist->Draw();
   tvx6_hist->Write();
  
   tvy6_hist->GetXaxis()->SetTitle("vy [cm]");
   tvy6_hist->Draw();
   tvy6_hist->Write();
   
   tvz6_hist->GetXaxis()->SetTitle("vz [cm]");
   tvz6_hist->Draw();
   tvz6_hist->Write();
   
   
   
   tvx_hist->GetXaxis()->SetTitle("vx [cm]");
   tvx_hist->Draw();
   tvx_hist->Write();
  
   tvy_hist->GetXaxis()->SetTitle("vy [cm]");
   tvy_hist->Draw();
   tvy_hist->Write();
   
   tvz_hist->GetXaxis()->SetTitle("vz [cm]");
   //TF1 *tvz1_func = new TF1("m1","gaus",-7,-5);
   //TF1 *tvz2_func = new TF1("m2","gaus",-5,-3);
   //TF1 *total = new TF1("totalm","gaus(0)+gaus(3)",-7,-3);
   //tvz_hist->Fit("total");
   tvz_hist->Draw();
   tvz_hist->Write();

   
 /* 

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
*/


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


