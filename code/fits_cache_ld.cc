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
  return par[0]*(1/(2*3.141592653))*(arg2/(4*arg3 + (arg4))) + par[3] + par[4]*pow(x[0],1) + par[5]*pow(x[0],2) + par[6]*pow(x[0],3) + par[7]*pow(x[0],4) + par[8]*pow(x[0],5);
}

//"BreitWinger"

void fits_cache_ld(){


   TFile input_1("ld_ib_kin_18312.root", "read");
   TFile input_2("ld_ib_kin_18324.root", "read");
   TFile input_3("ld_ib_kin_18325.root", "read");
   TFile input_4("ld_ib_kin_18326.root", "read");
   TFile input_5("ld_ib_kin_18330.root", "read");
   TFile input_6("ld_ib_kin_18331.root", "read");
   TFile input_7("ld_ib_kin_18335.root", "read");
   TFile input_8("ld_ib_kin_18336.root", "read");
   //TFile input_9("cu_su_18632_sep_hist_3.root", "read");
   //TFile input_10("cu_su_18636_sep_hist_3.root", "read");
   //TFile input_11("cu_su_19100_sep_hist_3.root", "read");
   //TFile input_12("cu_su_19101_sep_hist_3.root", "read");
   //TFile input_13("cu_su_19102_sep_hist_3.root", "read");
   //TFile input_14("cu_su_19103_sep_hist_3.root", "read");
   //TFile input_15("cu_su_19107_sep_hist_3.root", "read");

   
    TFile f("aaaaaaaaaaaaaaaald_fits_sep_8_3.root", "update");

   TH1D *cu1_hist = new TH1D("sector 0 P", "Q^{2} (1-2) [GeV^{2}]", 200, 0, 2);
	cu1_hist->SetStats(0);
	cu1_hist->SetLineColor(kOrange-3);
	
	TH1D *cu2_hist = new TH1D("sector 1 P", "Q^{2} (2-2.5) [GeV^{2}]", 200, 0, 2);
	cu2_hist->SetStats(0);
	cu2_hist->SetLineColor(kOrange-3);
	
	TH1D *cu3_hist = new TH1D("sector 2 P", "Q^{2} (2.5-3) [GeV^{2}]", 200, 0, 2);
	cu3_hist->SetStats(0);
	cu3_hist->SetLineColor(kOrange-3);
	
	TH1D *cu4_hist = new TH1D("sector 3 P", "Q^{2} (3-3.5) [GeV^{2}]", 200, 0, 2);
	cu4_hist->SetStats(0);
	cu4_hist->SetLineColor(kOrange-3);
	
	TH1D *cu5_hist = new TH1D("sector 4 P", "Q^{2} (3.5-6) [GeV^{2}]", 200, 0, 2);
	cu5_hist->SetStats(0);
	cu5_hist->SetLineColor(kOrange-3);
	
	TH1D *cu6_hist = new TH1D("sector 5 P", "sector 5 P", 200, 0, 2);
	cu6_hist->SetStats(0);
        cu6_hist->SetLineColor(kOrange-3);
        
      TH1D *cu1_1_hist = (TH1D*)input_1.Get("sector 0 P");
      TH1D *cu2_1_hist = (TH1D*)input_1.Get("sector 1 P");
      TH1D *cu3_1_hist = (TH1D*)input_1.Get("sector 2 P");
      TH1D *cu4_1_hist = (TH1D*)input_1.Get("sector 3 P");
      TH1D *cu5_1_hist = (TH1D*)input_1.Get("sector 4 P");
      TH1D *cu6_1_hist = (TH1D*)input_1.Get("sector 5 P");
      TH1D *cu1_2_hist = (TH1D*)input_2.Get("sector 0 P");
      TH1D *cu2_2_hist = (TH1D*)input_2.Get("sector 1 P");
      TH1D *cu3_2_hist = (TH1D*)input_2.Get("sector 2 P");
      TH1D *cu4_2_hist = (TH1D*)input_2.Get("sector 3 P");
      TH1D *cu5_2_hist = (TH1D*)input_2.Get("sector 4 P");
      TH1D *cu6_2_hist = (TH1D*)input_2.Get("sector 5 P");
      TH1D *cu1_3_hist = (TH1D*)input_3.Get("sector 0 P");
      TH1D *cu2_3_hist = (TH1D*)input_3.Get("sector 1 P");
      TH1D *cu3_3_hist = (TH1D*)input_3.Get("sector 2 P");
      TH1D *cu4_3_hist = (TH1D*)input_3.Get("sector 3 P");
      TH1D *cu5_3_hist = (TH1D*)input_3.Get("sector 4 P");
      TH1D *cu6_3_hist = (TH1D*)input_3.Get("sector 5 P");
      TH1D *cu1_4_hist = (TH1D*)input_4.Get("sector 0 P");
      TH1D *cu2_4_hist = (TH1D*)input_4.Get("sector 1 P");
      TH1D *cu3_4_hist = (TH1D*)input_4.Get("sector 2 P");
      TH1D *cu4_4_hist = (TH1D*)input_4.Get("sector 3 P");
      TH1D *cu5_4_hist = (TH1D*)input_4.Get("sector 4 P");
      TH1D *cu6_4_hist = (TH1D*)input_4.Get("sector 5 P");
      TH1D *cu1_5_hist = (TH1D*)input_5.Get("sector 0 P");
      TH1D *cu2_5_hist = (TH1D*)input_5.Get("sector 1 P");
      TH1D *cu3_5_hist = (TH1D*)input_5.Get("sector 2 P");
      TH1D *cu4_5_hist = (TH1D*)input_5.Get("sector 3 P");
      TH1D *cu5_5_hist = (TH1D*)input_5.Get("sector 4 P");
      TH1D *cu6_5_hist = (TH1D*)input_5.Get("sector 5 P");
      TH1D *cu1_6_hist = (TH1D*)input_6.Get("sector 0 P");
      TH1D *cu2_6_hist = (TH1D*)input_6.Get("sector 1 P");
      TH1D *cu3_6_hist = (TH1D*)input_6.Get("sector 2 P");
      TH1D *cu4_6_hist = (TH1D*)input_6.Get("sector 3 P");
      TH1D *cu5_6_hist = (TH1D*)input_6.Get("sector 4 P");
      TH1D *cu6_6_hist = (TH1D*)input_6.Get("sector 5 P");
      TH1D *cu1_7_hist = (TH1D*)input_7.Get("sector 0 P");
      TH1D *cu2_7_hist = (TH1D*)input_7.Get("sector 1 P");
      TH1D *cu3_7_hist = (TH1D*)input_7.Get("sector 2 P");
      TH1D *cu4_7_hist = (TH1D*)input_7.Get("sector 3 P");
      TH1D *cu5_7_hist = (TH1D*)input_7.Get("sector 4 P");
      TH1D *cu6_7_hist = (TH1D*)input_7.Get("sector 5 P");
      TH1D *cu1_8_hist = (TH1D*)input_8.Get("sector 0 P");
      TH1D *cu2_8_hist = (TH1D*)input_8.Get("sector 1 P");
      TH1D *cu3_8_hist = (TH1D*)input_8.Get("sector 2 P");
      TH1D *cu4_8_hist = (TH1D*)input_8.Get("sector 3 P");
      TH1D *cu5_8_hist = (TH1D*)input_8.Get("sector 4 P");
      TH1D *cu6_8_hist = (TH1D*)input_8.Get("sector 5 P");
      //TH1D *cu1_9_hist = (TH1D*)input_9.Get("sector 0 P");
      //TH1D *cu2_9_hist = (TH1D*)input_9.Get("sector 1 P");
      //TH1D *cu3_9_hist = (TH1D*)input_9.Get("sector 2 P");
      //TH1D *cu4_9_hist = (TH1D*)input_9.Get("sector 3 P");
      //TH1D *cu5_9_hist = (TH1D*)input_9.Get("sector 4 P");
      //TH1D *cu6_9_hist = (TH1D*)input_9.Get("sector 5 P");
      //TH1D *cu1_10_hist = (TH1D*)input_10.Get("sector 0 P");
      //TH1D *cu2_10_hist = (TH1D*)input_10.Get("sector 1 P");
      //TH1D *cu3_10_hist = (TH1D*)input_10.Get("sector 2 P");
      //TH1D *cu4_10_hist = (TH1D*)input_10.Get("sector 3 P");
      //TH1D *cu5_10_hist = (TH1D*)input_10.Get("sector 4 P");
      //TH1D *cu6_10_hist = (TH1D*)input_10.Get("sector 5 P");
      //TH1D *cu1_11_hist = (TH1D*)input_11.Get("sector 0 P");
      //TH1D *cu2_11_hist = (TH1D*)input_11.Get("sector 1 P");
      //TH1D *cu3_11_hist = (TH1D*)input_11.Get("sector 2 P");
      //TH1D *cu4_11_hist = (TH1D*)input_11.Get("sector 3 P");
      //TH1D *cu5_11_hist = (TH1D*)input_11.Get("sector 4 P");
      //TH1D *cu6_11_hist = (TH1D*)input_11.Get("sector 5 P");
      //TH1D *cu1_12_hist = (TH1D*)input_12.Get("sector 0 P");
      //TH1D *cu2_12_hist = (TH1D*)input_12.Get("sector 1 P");
      //TH1D *cu3_12_hist = (TH1D*)input_12.Get("sector 2 P");
      //TH1D *cu4_12_hist = (TH1D*)input_12.Get("sector 3 P");
      //TH1D *cu5_12_hist = (TH1D*)input_12.Get("sector 4 P");
      //TH1D *cu6_12_hist = (TH1D*)input_12.Get("sector 5 P");
      //TH1D *cu1_13_hist = (TH1D*)input_13.Get("sector 0 P");
      //TH1D *cu2_13_hist = (TH1D*)input_13.Get("sector 1 P");
      //TH1D *cu3_13_hist = (TH1D*)input_13.Get("sector 2 P");
      //TH1D *cu4_13_hist = (TH1D*)input_13.Get("sector 3 P");
      //TH1D *cu5_13_hist = (TH1D*)input_13.Get("sector 4 P");
      //TH1D *cu6_13_hist = (TH1D*)input_13.Get("sector 5 P");
      //TH1D *cu1_14_hist = (TH1D*)input_14.Get("sector 0 P");
      //TH1D *cu2_14_hist = (TH1D*)input_14.Get("sector 1 P");
      //TH1D *cu3_14_hist = (TH1D*)input_14.Get("sector 2 P");
      //TH1D *cu4_14_hist = (TH1D*)input_14.Get("sector 3 P");
      //TH1D *cu5_14_hist = (TH1D*)input_14.Get("sector 4 P");
      //TH1D *cu6_14_hist = (TH1D*)input_14.Get("sector 5 P");
      //TH1D *cu1_15_hist = (TH1D*)input_15.Get("sector 0 P");
      //TH1D *cu2_15_hist = (TH1D*)input_15.Get("sector 1 P");
      //TH1D *cu3_15_hist = (TH1D*)input_15.Get("sector 2 P");
      //TH1D *cu4_15_hist = (TH1D*)input_15.Get("sector 3 P");
      //TH1D *cu5_15_hist = (TH1D*)input_15.Get("sector 4 P");
      //TH1D *cu6_15_hist = (TH1D*)input_15.Get("sector 5 P");
      

      cu1_hist->Add(cu1_1_hist);
      cu1_hist->Add(cu1_2_hist);
      cu1_hist->Add(cu1_3_hist);
      cu1_hist->Add(cu1_4_hist);
      cu1_hist->Add(cu1_5_hist);
      cu1_hist->Add(cu1_6_hist);
      cu1_hist->Add(cu1_7_hist);
      cu1_hist->Add(cu1_8_hist);
      //cu1_hist->Add(cu1_9_hist);
      //cu1_hist->Add(cu1_10_hist);
      //cu1_hist->Add(cu1_11_hist);
      //cu1_hist->Add(cu1_12_hist);
      //cu1_hist->Add(cu1_13_hist);
      //cu1_hist->Add(cu1_14_hist);
      //cu1_hist->Add(cu1_15_hist);

      cu2_hist->Add(cu2_1_hist);
      cu2_hist->Add(cu2_2_hist);
      cu2_hist->Add(cu2_3_hist);
      cu2_hist->Add(cu2_4_hist);
      cu2_hist->Add(cu2_5_hist);
      cu2_hist->Add(cu2_6_hist);
      cu2_hist->Add(cu2_7_hist);
      cu2_hist->Add(cu2_8_hist);
      //cu2_hist->Add(cu2_9_hist);
      //cu2_hist->Add(cu2_10_hist);
      //cu2_hist->Add(cu2_11_hist);
      //cu2_hist->Add(cu2_12_hist);
      //cu2_hist->Add(cu2_13_hist);
      //cu2_hist->Add(cu2_14_hist);
      //cu2_hist->Add(cu2_15_hist);

      cu3_hist->Add(cu3_1_hist);
      cu3_hist->Add(cu3_2_hist);
      cu3_hist->Add(cu3_3_hist);
      cu3_hist->Add(cu3_4_hist);
      cu3_hist->Add(cu3_5_hist);
      cu3_hist->Add(cu3_6_hist);
      cu3_hist->Add(cu3_7_hist);
      cu3_hist->Add(cu3_8_hist);
      //cu3_hist->Add(cu3_9_hist);
      //cu3_hist->Add(cu3_10_hist);
      //cu3_hist->Add(cu3_11_hist);
      //cu3_hist->Add(cu3_12_hist);
      //cu3_hist->Add(cu3_13_hist);
      //cu3_hist->Add(cu3_14_hist);
      //cu3_hist->Add(cu3_15_hist);

      cu4_hist->Add(cu4_1_hist);
      cu4_hist->Add(cu4_2_hist);
      cu4_hist->Add(cu4_3_hist);
      cu4_hist->Add(cu4_4_hist);
      cu4_hist->Add(cu4_5_hist);
      cu4_hist->Add(cu4_6_hist);
      cu4_hist->Add(cu4_7_hist);
      cu4_hist->Add(cu4_8_hist);
      //cu4_hist->Add(cu4_9_hist);
      //cu4_hist->Add(cu4_10_hist);
      //cu4_hist->Add(cu4_11_hist);
      //cu4_hist->Add(cu4_12_hist);
      //cu4_hist->Add(cu4_13_hist);
      //cu4_hist->Add(cu4_14_hist);
      //cu4_hist->Add(cu4_15_hist);

      cu5_hist->Add(cu5_1_hist);
      cu5_hist->Add(cu5_2_hist);
      cu5_hist->Add(cu5_3_hist);
      cu5_hist->Add(cu5_4_hist);
      cu5_hist->Add(cu5_5_hist);
      cu5_hist->Add(cu5_6_hist);
      cu5_hist->Add(cu5_7_hist);
      cu5_hist->Add(cu5_8_hist);
      //cu5_hist->Add(cu5_9_hist);
      //cu5_hist->Add(cu5_10_hist);
      //cu5_hist->Add(cu5_11_hist);
      //cu5_hist->Add(cu5_12_hist);
      //cu5_hist->Add(cu5_13_hist);
      //cu5_hist->Add(cu5_14_hist);
      //cu5_hist->Add(cu5_15_hist);

      cu5_hist->Add(cu6_1_hist);
      cu5_hist->Add(cu6_2_hist);
      cu5_hist->Add(cu6_3_hist);
      cu5_hist->Add(cu6_4_hist);
      cu5_hist->Add(cu6_5_hist);
      cu5_hist->Add(cu6_6_hist);
      cu5_hist->Add(cu6_7_hist);
      cu5_hist->Add(cu6_8_hist);
      //cu5_hist->Add(cu6_9_hist);
      //cu5_hist->Add(cu6_10_hist);
      //cu5_hist->Add(cu6_11_hist);
      //cu5_hist->Add(cu6_12_hist);
      //cu5_hist->Add(cu6_13_hist);
      //cu5_hist->Add(cu6_14_hist);
      //cu5_hist->Add(cu6_15_hist);
      

      
  


   std::vector<std::string> v3_a;
   int num = 0;
   
   
   cu1_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //cu1_hist->Scale(1.0 / cu1_hist->Integral());
   TF1 *cu1_func = new TF1("cu1_fit",fitfbw,0.4,1.4,8);
   cu1_func->SetParameter(1,0.150);
   cu1_func->SetParameter(2,0.77);
   //cu1_func->SetParLimits(2,0.65,0.89);
   cu1_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   cu1_hist->Fit("cu1_fit","","",0.4,1.4);
   cu1_hist->Draw();
   cu1_hist->Write();

   cu2_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //cu2_hist->Scale(1.0 / cu2_hist->Integral());
   TF1 *cu2_func = new TF1("cu2_fit",fitfbw,0.4,1.4,8);
   cu2_func->SetParameter(1,0.150);
   cu2_func->SetParameter(2,0.77);
   //cu2_func->SetParLimits(2,0.65,0.89);
   cu2_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   cu2_hist->Fit("cu2_fit","","",0.4,1.4);
   cu2_hist->Draw();
   cu2_hist->Write();

   cu3_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //cu3_hist->Scale(1.0 / cu3_hist->Integral());
   TF1 *cu3_func = new TF1("cu3_fit",fitfbw,0.4,1.4,8);
   cu3_func->SetParameter(1,0.150);
   cu3_func->SetParameter(2,0.77);
   //cu3_func->SetParLimits(2,0.65,0.89);
   cu3_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   cu3_hist->Fit("cu3_fit","","",0.4,1.4);
   cu3_hist->Draw();
   cu3_hist->Write();

   cu4_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //cu4_hist->Scale(1.0 / cu4_hist->Integral());
   TF1 *cu4_func = new TF1("cu4_fit",fitfbw,0.4,1.4,8);
   cu4_func->SetParameter(1,0.150);
   cu4_func->SetParameter(2,0.77);
   //cu4_func->SetParLimits(2,0.65,0.89);
   cu4_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   cu4_hist->Fit("cu4_fit","","",0.4,1.4);
   cu4_hist->Draw();
   cu4_hist->Write();

   cu5_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //cu5_hist->Scale(1.0 / cu5_hist->Integral());
   TF1 *cu5_func = new TF1("cu5_fit",fitfbw,0.4,1.4,8);
   cu5_func->SetParameter(1,0.150);
   cu5_func->SetParameter(2,0.77);
   //cu5_func->SetParLimits(2,0.65,0.89);
   cu5_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   cu5_hist->Fit("cu5_fit","","",0.4,1.4);
   cu5_hist->Draw();
   cu5_hist->Write();

   cu6_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //cu6_hist->Scale(1.0 / cu6_hist->Integral());
   TF1 *cu6_func = new TF1("cu6_fit",fitfbw,0.4,1.4,8);
   cu6_func->SetParameter(1,0.150);
   cu6_func->SetParameter(2,0.77);
   //cu6_func->SetParLimits(2,0.65,0.89);
   cu6_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   cu6_hist->Fit("cu6_fit","","",0.4,1.4);
   cu6_hist->Draw();
   cu6_hist->Write();
   
   



   
   
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


