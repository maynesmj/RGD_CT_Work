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
  return par[0]*(arg2/(4*arg3 + (arg4))) + par[3] + par[4]*pow(x[0],1) + par[5]*pow(x[0],2) + par[6]*pow(x[0],3); // + par[7]*pow(x[0],4) + par[8]*pow(x[0],5);
}

//"BreitWinger"

void fits_snv2(){


   TFile input_1("cu_su_18624_sep_hist_5.root", "read");
   TFile input_2("cu_su_18625_sep_hist_5.root", "read");
   TFile input_3("cu_su_18626_sep_hist_5.root", "read");
   TFile input_4("cu_su_18627_sep_hist_5.root", "read");
   TFile input_5("cu_su_18628_sep_hist_5.root", "read");
   TFile input_6("cu_su_18629_sep_hist_5.root", "read");
   TFile input_7("cu_su_18630_sep_hist_5.root", "read");
   TFile input_8("cu_su_18631_sep_hist_5.root", "read");
   TFile input_9("cu_su_18632_sep_hist_5.root", "read");
   TFile input_10("cu_su_18636_sep_hist_5.root", "read");
   TFile input_11("cu_su_19100_sep_hist_5.root", "read");
   TFile input_12("cu_su_19101_sep_hist_5.root", "read");
   TFile input_13("cu_su_19102_sep_hist_5.root", "read");
   TFile input_14("cu_su_19103_sep_hist_5.root", "read");
   TFile input_15("cu_su_19107_sep_hist_5.root", "read");

   
    TFile f("cu_sn_fits_sep_7.root", "update");

   TH1D *sn1_hist = new TH1D("sn bin1", "Sn bin1", 300, 0, 2);
	sn1_hist->SetStats(0);
	
	TH1D *sn2_hist = new TH1D("sn bin2", "Sn bin2", 300, 0, 2);
	sn2_hist->SetStats(0);
	
	TH1D *sn3_hist = new TH1D("sn bin3", "Sn bin3", 300, 0, 2);
	sn3_hist->SetStats(0);
	
	TH1D *sn4_hist = new TH1D("sn bin4", "Sn bin4", 300, 0, 2);
	sn4_hist->SetStats(0);
	
	TH1D *sn5_hist = new TH1D("sn bin5", "Sn bin5", 300, 0, 2);
	sn5_hist->SetStats(0);
	
	TH1D *sn6_hist = new TH1D("sn bin6", "Sn bin6", 300, 0, 2);
	sn6_hist->SetStats(0);

      TH1D *sn1_1_hist = (TH1D*)input_1.Get("Sn bin1");
      TH1D *sn2_1_hist = (TH1D*)input_1.Get("Sn bin2");
      TH1D *sn3_1_hist = (TH1D*)input_1.Get("Sn bin3");
      TH1D *sn4_1_hist = (TH1D*)input_1.Get("Sn bin4");
      TH1D *sn5_1_hist = (TH1D*)input_1.Get("Sn bin5");
      TH1D *sn6_1_hist = (TH1D*)input_1.Get("Sn bin6");
      TH1D *sn1_2_hist = (TH1D*)input_2.Get("Sn bin1");
      TH1D *sn2_2_hist = (TH1D*)input_2.Get("Sn bin2");
      TH1D *sn3_2_hist = (TH1D*)input_2.Get("Sn bin3");
      TH1D *sn4_2_hist = (TH1D*)input_2.Get("Sn bin4");
      TH1D *sn5_2_hist = (TH1D*)input_2.Get("Sn bin5");
      TH1D *sn6_2_hist = (TH1D*)input_2.Get("Sn bin6");
      TH1D *sn1_3_hist = (TH1D*)input_3.Get("Sn bin1");
      TH1D *sn2_3_hist = (TH1D*)input_3.Get("Sn bin2");
      TH1D *sn3_3_hist = (TH1D*)input_3.Get("Sn bin3");
      TH1D *sn4_3_hist = (TH1D*)input_3.Get("Sn bin4");
      TH1D *sn5_3_hist = (TH1D*)input_3.Get("Sn bin5");
      TH1D *sn6_3_hist = (TH1D*)input_3.Get("Sn bin6");
      TH1D *sn1_4_hist = (TH1D*)input_4.Get("Sn bin1");
      TH1D *sn2_4_hist = (TH1D*)input_4.Get("Sn bin2");
      TH1D *sn3_4_hist = (TH1D*)input_4.Get("Sn bin3");
      TH1D *sn4_4_hist = (TH1D*)input_4.Get("Sn bin4");
      TH1D *sn5_4_hist = (TH1D*)input_4.Get("Sn bin5");
      TH1D *sn6_4_hist = (TH1D*)input_4.Get("Sn bin6");
      TH1D *sn1_5_hist = (TH1D*)input_5.Get("Sn bin1");
      TH1D *sn2_5_hist = (TH1D*)input_5.Get("Sn bin2");
      TH1D *sn3_5_hist = (TH1D*)input_5.Get("Sn bin3");
      TH1D *sn4_5_hist = (TH1D*)input_5.Get("Sn bin4");
      TH1D *sn5_5_hist = (TH1D*)input_5.Get("Sn bin5");
      TH1D *sn6_5_hist = (TH1D*)input_5.Get("Sn bin6");
      TH1D *sn1_6_hist = (TH1D*)input_6.Get("Sn bin1");
      TH1D *sn2_6_hist = (TH1D*)input_6.Get("Sn bin2");
      TH1D *sn3_6_hist = (TH1D*)input_6.Get("Sn bin3");
      TH1D *sn4_6_hist = (TH1D*)input_6.Get("Sn bin4");
      TH1D *sn5_6_hist = (TH1D*)input_6.Get("Sn bin5");
      TH1D *sn6_6_hist = (TH1D*)input_6.Get("Sn bin6");
      TH1D *sn1_7_hist = (TH1D*)input_7.Get("Sn bin1");
      TH1D *sn2_7_hist = (TH1D*)input_7.Get("Sn bin2");
      TH1D *sn3_7_hist = (TH1D*)input_7.Get("Sn bin3");
      TH1D *sn4_7_hist = (TH1D*)input_7.Get("Sn bin4");
      TH1D *sn5_7_hist = (TH1D*)input_7.Get("Sn bin5");
      TH1D *sn6_7_hist = (TH1D*)input_7.Get("Sn bin6");
      TH1D *sn1_8_hist = (TH1D*)input_8.Get("Sn bin1");
      TH1D *sn2_8_hist = (TH1D*)input_8.Get("Sn bin2");
      TH1D *sn3_8_hist = (TH1D*)input_8.Get("Sn bin3");
      TH1D *sn4_8_hist = (TH1D*)input_8.Get("Sn bin4");
      TH1D *sn5_8_hist = (TH1D*)input_8.Get("Sn bin5");
      TH1D *sn6_8_hist = (TH1D*)input_8.Get("Sn bin6");
      TH1D *sn1_9_hist = (TH1D*)input_1.Get("Sn bin1");
      TH1D *sn2_9_hist = (TH1D*)input_1.Get("Sn bin2");
      TH1D *sn3_9_hist = (TH1D*)input_1.Get("Sn bin3");
      TH1D *sn4_9_hist = (TH1D*)input_1.Get("Sn bin4");
      TH1D *sn5_9_hist = (TH1D*)input_1.Get("Sn bin5");
      TH1D *sn6_9_hist = (TH1D*)input_1.Get("Sn bin6");
      TH1D *sn1_10_hist = (TH1D*)input_10.Get("Sn bin1");
      TH1D *sn2_10_hist = (TH1D*)input_10.Get("Sn bin2");
      TH1D *sn3_10_hist = (TH1D*)input_10.Get("Sn bin3");
      TH1D *sn4_10_hist = (TH1D*)input_10.Get("Sn bin4");
      TH1D *sn5_10_hist = (TH1D*)input_10.Get("Sn bin5");
      TH1D *sn6_10_hist = (TH1D*)input_10.Get("Sn bin6");
      TH1D *sn1_11_hist = (TH1D*)input_11.Get("Sn bin1");
      TH1D *sn2_11_hist = (TH1D*)input_11.Get("Sn bin2");
      TH1D *sn3_11_hist = (TH1D*)input_11.Get("Sn bin3");
      TH1D *sn4_11_hist = (TH1D*)input_11.Get("Sn bin4");
      TH1D *sn5_11_hist = (TH1D*)input_11.Get("Sn bin5");
      TH1D *sn6_11_hist = (TH1D*)input_11.Get("Sn bin6");
      TH1D *sn1_12_hist = (TH1D*)input_12.Get("Sn bin1");
      TH1D *sn2_12_hist = (TH1D*)input_12.Get("Sn bin2");
      TH1D *sn3_12_hist = (TH1D*)input_12.Get("Sn bin3");
      TH1D *sn4_12_hist = (TH1D*)input_12.Get("Sn bin4");
      TH1D *sn5_12_hist = (TH1D*)input_12.Get("Sn bin5");
      TH1D *sn6_12_hist = (TH1D*)input_12.Get("Sn bin6");
      TH1D *sn1_13_hist = (TH1D*)input_13.Get("Sn bin1");
      TH1D *sn2_13_hist = (TH1D*)input_13.Get("Sn bin2");
      TH1D *sn3_13_hist = (TH1D*)input_13.Get("Sn bin3");
      TH1D *sn4_13_hist = (TH1D*)input_13.Get("Sn bin4");
      TH1D *sn5_13_hist = (TH1D*)input_13.Get("Sn bin5");
      TH1D *sn6_13_hist = (TH1D*)input_13.Get("Sn bin6");
      TH1D *sn1_14_hist = (TH1D*)input_14.Get("Sn bin1");
      TH1D *sn2_14_hist = (TH1D*)input_14.Get("Sn bin2");
      TH1D *sn3_14_hist = (TH1D*)input_14.Get("Sn bin3");
      TH1D *sn4_14_hist = (TH1D*)input_14.Get("Sn bin4");
      TH1D *sn5_14_hist = (TH1D*)input_14.Get("Sn bin5");
      TH1D *sn6_14_hist = (TH1D*)input_14.Get("Sn bin6");
      TH1D *sn1_15_hist = (TH1D*)input_15.Get("Sn bin1");
      TH1D *sn2_15_hist = (TH1D*)input_15.Get("Sn bin2");
      TH1D *sn3_15_hist = (TH1D*)input_15.Get("Sn bin3");
      TH1D *sn4_15_hist = (TH1D*)input_15.Get("Sn bin4");
      TH1D *sn5_15_hist = (TH1D*)input_15.Get("Sn bin5");
      TH1D *sn6_15_hist = (TH1D*)input_15.Get("Sn bin6");
      

      sn1_hist->Add(sn1_1_hist);
      sn1_hist->Add(sn1_2_hist);
      sn1_hist->Add(sn1_3_hist);
      sn1_hist->Add(sn1_4_hist);
      sn1_hist->Add(sn1_5_hist);
      sn1_hist->Add(sn1_6_hist);
      sn1_hist->Add(sn1_7_hist);
      sn1_hist->Add(sn1_8_hist);
      sn1_hist->Add(sn1_9_hist);
      sn1_hist->Add(sn1_10_hist);
      sn1_hist->Add(sn1_11_hist);
      sn1_hist->Add(sn1_12_hist);
      sn1_hist->Add(sn1_13_hist);
      sn1_hist->Add(sn1_14_hist);
      sn1_hist->Add(sn1_15_hist);

      sn2_hist->Add(sn2_1_hist);
      sn2_hist->Add(sn2_2_hist);
      sn2_hist->Add(sn2_3_hist);
      sn2_hist->Add(sn2_4_hist);
      sn2_hist->Add(sn2_5_hist);
      sn2_hist->Add(sn2_6_hist);
      sn2_hist->Add(sn2_7_hist);
      sn2_hist->Add(sn2_8_hist);
      sn2_hist->Add(sn2_9_hist);
      sn2_hist->Add(sn2_10_hist);
      sn2_hist->Add(sn2_11_hist);
      sn2_hist->Add(sn2_12_hist);
      sn2_hist->Add(sn2_13_hist);
      sn2_hist->Add(sn2_14_hist);
      sn2_hist->Add(sn2_15_hist);

      sn3_hist->Add(sn3_1_hist);
      sn3_hist->Add(sn3_2_hist);
      sn3_hist->Add(sn3_3_hist);
      sn3_hist->Add(sn3_4_hist);
      sn3_hist->Add(sn3_5_hist);
      sn3_hist->Add(sn3_6_hist);
      sn3_hist->Add(sn3_7_hist);
      sn3_hist->Add(sn3_8_hist);
      sn3_hist->Add(sn3_9_hist);
      sn3_hist->Add(sn3_10_hist);
      sn3_hist->Add(sn3_11_hist);
      sn3_hist->Add(sn3_12_hist);
      sn3_hist->Add(sn3_13_hist);
      sn3_hist->Add(sn3_14_hist);
      sn3_hist->Add(sn3_15_hist);

      sn4_hist->Add(sn4_1_hist);
      sn4_hist->Add(sn4_2_hist);
      sn4_hist->Add(sn4_3_hist);
      sn4_hist->Add(sn4_4_hist);
      sn4_hist->Add(sn4_5_hist);
      sn4_hist->Add(sn4_6_hist);
      sn4_hist->Add(sn4_7_hist);
      sn4_hist->Add(sn4_8_hist);
      sn4_hist->Add(sn4_9_hist);
      sn4_hist->Add(sn4_10_hist);
      sn4_hist->Add(sn4_11_hist);
      sn4_hist->Add(sn4_12_hist);
      sn4_hist->Add(sn4_13_hist);
      sn4_hist->Add(sn4_14_hist);
      sn4_hist->Add(sn4_15_hist);

      sn5_hist->Add(sn5_1_hist);
      sn5_hist->Add(sn5_2_hist);
      sn5_hist->Add(sn5_3_hist);
      sn5_hist->Add(sn5_4_hist);
      sn5_hist->Add(sn5_5_hist);
      sn5_hist->Add(sn5_6_hist);
      sn5_hist->Add(sn5_7_hist);
      sn5_hist->Add(sn5_8_hist);
      sn5_hist->Add(sn5_9_hist);
      sn5_hist->Add(sn5_10_hist);
      sn5_hist->Add(sn5_11_hist);
      sn5_hist->Add(sn5_12_hist);
      sn5_hist->Add(sn5_13_hist);
      sn5_hist->Add(sn5_14_hist);
      sn5_hist->Add(sn5_15_hist);

      sn5_hist->Add(sn6_1_hist);
      sn5_hist->Add(sn6_2_hist);
      sn5_hist->Add(sn6_3_hist);
      sn5_hist->Add(sn6_4_hist);
      sn5_hist->Add(sn6_5_hist);
      sn5_hist->Add(sn6_6_hist);
      sn5_hist->Add(sn6_7_hist);
      sn5_hist->Add(sn6_8_hist);
      sn5_hist->Add(sn6_9_hist);
      sn5_hist->Add(sn6_10_hist);
      sn5_hist->Add(sn6_11_hist);
      sn5_hist->Add(sn6_12_hist);
      sn5_hist->Add(sn6_13_hist);
      sn5_hist->Add(sn6_14_hist);
      sn5_hist->Add(sn6_15_hist);
      

      
  


   std::vector<std::string> v3_a;
   int num = 0;
   
   
   sn1_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   TF1 *sn1_func = new TF1("sn1_fit",fitfbw,0.2,1.4,7);
   sn1_func->SetParameter(1,0.150);
   sn1_func->SetParameter(2,0.77);
   sn1_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   sn1_hist->Fit("sn1_fit","","",0.2,1.4);
   sn1_hist->Draw();
   sn1_hist->Write();

   sn2_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   TF1 *sn2_func = new TF1("sn2_fit",fitfbw,0.2,1.4,7);
   sn2_func->SetParameter(1,0.150);
   sn2_func->SetParameter(2,0.77);
   sn2_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   sn2_hist->Fit("sn2_fit","","",0.2,1.4);
   sn2_hist->Draw();
   sn2_hist->Write();

   sn3_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   TF1 *sn3_func = new TF1("sn3_fit",fitfbw,0.2,1.4,7);
   sn3_func->SetParameter(1,0.150);
   sn3_func->SetParameter(2,0.77);
   sn3_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   sn3_hist->Fit("sn3_fit","","",0.2,1.4);
   sn3_hist->Draw();
   sn3_hist->Write();

   sn4_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   TF1 *sn4_func = new TF1("sn4_fit",fitfbw,0.2,1.4,7);
   sn4_func->SetParameter(1,0.150);
   sn4_func->SetParameter(2,0.77);
   sn4_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   sn4_hist->Fit("sn4_fit","","",0.2,1.4);
   sn4_hist->Draw();
   sn4_hist->Write();

   sn5_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   TF1 *sn5_func = new TF1("sn5_fit",fitfbw,0.2,1.4,7);
   sn5_func->SetParameter(1,0.150);
   sn5_func->SetParameter(2,0.77);
   sn5_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   sn5_hist->Fit("sn5_fit","","",0.2,1.4);
   sn5_hist->Draw();
   sn5_hist->Write();

   sn6_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   TF1 *sn6_func = new TF1("sn6_fit",fitfbw,0.2,1.4,7);
   sn6_func->SetParameter(1,0.150);
   sn6_func->SetParameter(2,0.77);
   sn6_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   sn6_hist->Fit("sn6_fit","","",0.2,1.4);
   sn6_hist->Draw();
   sn6_hist->Write();
   
   



   
   
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


