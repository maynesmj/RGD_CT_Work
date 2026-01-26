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
   Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg) + par[3]*(x[0])*(x[0]) + par[4]*(x[0]) + par[5];
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
  return par[0]*(arg2/(4*arg3 + arg4)) + par[3] + par[4]*pow(x[0],1) + par[5]*pow(x[0],2) + par[6]*pow(x[0],3); //+ par[7]*pow(x[0],4) + par[8]*pow(x[0],5);
}

//"BreitWinger"

void fitsaps_4(){


   TFile input("cusn_pass0v4_kin5_300_no_sigma.root", "read");    ///fits_ld2.root
   TFile input2("cusn_ob_bspot_kin4_5_sigma.root", "read"); 
   
    TFile f("aaa_cusn_vol_stacked_5sig_vertex_pass0v4bspot_4.root", "update");
   gStyle->SetOptStat("");
	gStyle->SetCanvasPreferGL();

   TH1F *vx1_hist = new TH1F("Sector 1 vx", "Sector 1 vx", 300, -5, 5);
       vx1_hist->SetStats(0);
       vx1_hist->SetLineColor(kRed+1);
   TH1D *vy1_hist = new TH1D("Sector 1 vy", "Sector 1 vy", 300, -5, 5);
       vy1_hist->SetStats(0);
       vy1_hist->SetLineColor(kRed+1);
   TH1D *vz1_hist = new TH1D("Sector 1 vz", "Sector 1 vz", 300, -20, 20);
       vz1_hist->SetStats(0);
       vz1_hist->SetLineColor(kRed+1);
       
  TH1F *vx2_hist = new TH1F("Sector 2 vx", "Sector 2 vx", 300, -5, 5);
       vx2_hist->SetStats(0);
       vx2_hist->SetLineColor(kRed+1);
   TH1D *vy2_hist = new TH1D("Sector 2 vy", "Sector 2 vy", 300, -5, 5);
       vy2_hist->SetStats(0);
       vy2_hist->SetLineColor(kRed+1);
   TH1D *vz2_hist = new TH1D("Sector 2 vz", "Sector 2 vz", 300, -20, 20);
       vz2_hist->SetStats(0);
       vz2_hist->SetLineColor(kRed+1);
             
   TH1F *vx3_hist = new TH1F("Sector 3 vx", "Sector 3 vx", 300, -5, 5);
       vx3_hist->SetStats(0);
       vx3_hist->SetLineColor(kRed+1);
   TH1D *vy3_hist = new TH1D("Sector 3 vy", "Sector 3 vy", 300, -5, 5);
       vy3_hist->SetStats(0);
       vy3_hist->SetLineColor(kRed+1);
   TH1D *vz3_hist = new TH1D("Sector 3 vz", "Sector 3 vz", 300, -20, 20);
       vz3_hist->SetStats(0);
       vz3_hist->SetLineColor(kRed+1);
       
   TH1F *vx4_hist = new TH1F("Sector 4 vx", "Sector 4 vx", 300, -5, 5);
       vx4_hist->SetStats(0);
       vx4_hist->SetLineColor(kRed+1);
   TH1D *vy4_hist = new TH1D("Sector 4 vy", "Sector 4 vy", 300, -5, 5);
       vy4_hist->SetStats(0);
       vy4_hist->SetLineColor(kRed+1);
   TH1D *vz4_hist = new TH1D("Sector 4 vz", "Sector 4 vz", 300, -20, 20);
       vz4_hist->SetStats(0);
       vz4_hist->SetLineColor(kRed+1);
       
   TH1F *vx5_hist = new TH1F("Sector 5 vx", "Sector 5 vx", 300, -5, 5);
       vx5_hist->SetStats(0);
       vx5_hist->SetLineColor(kRed+1);
   TH1D *vy5_hist = new TH1D("Sector 5 vy", "Sector 5 vy", 300, -5, 5);
       vy5_hist->SetStats(0);
       vy5_hist->SetLineColor(kRed+1);
   TH1D *vz5_hist = new TH1D("Sector 5 vz", "Sector 5 vz", 300, -20, 20);
       vz5_hist->SetStats(0);
       vz5_hist->SetLineColor(kRed+1);
       
   TH1F *vx6_hist = new TH1F("Sector 6 vx", "Sector 6 vx", 300, -5, 5);
       vx6_hist->SetStats(0);
       vx6_hist->SetLineColor(kRed+1);
   TH1D *vy6_hist = new TH1D("Sector 6 vy", "Sector 6 vy", 300, -5, 5);
       vy6_hist->SetStats(0);
       vy6_hist->SetLineColor(kRed+1);
   TH1D *vz6_hist = new TH1D("Sector 6 vz", "Sector 6 vz", 300, -20, 20);
       vz6_hist->SetStats(0);
       vz6_hist->SetLineColor(kRed+1);
       
       TH1F *yvx1_hist = new TH1F("Sector 1 vx", "Sector 1 vx", 300, -5, 5);
       yvx1_hist->SetStats(0);
       yvx1_hist->SetLineColor(kCyan+1);
   TH1D *yvy1_hist = new TH1D("Sector 1 vy", "Sector 1 vy", 300, -5, 5);
       yvy1_hist->SetStats(0);
       yvy1_hist->SetLineColor(kCyan-1);
   TH1D *yvz1_hist = new TH1D("Sector 1 vz", "Sector 1 vz", 300, -20, 20);
       yvz1_hist->SetStats(0);
       yvz1_hist->SetLineColor(kCyan+1);
       
  TH1F *yvx2_hist = new TH1F("Sector 2 vx", "Sector 2 vx", 300, -5, 5);
       yvx2_hist->SetStats(0);
       yvx2_hist->SetLineColor(kOrange+1);
   TH1D *yvy2_hist = new TH1D("Sector 2 vy", "Sector 2 vy", 300, -5, 5);
       yvy2_hist->SetStats(0);
       yvy2_hist->SetLineColor(kOrange-1);
   TH1D *yvz2_hist = new TH1D("Sector 2 vz", "Sector 2 vz", 300, -20, 20);
       yvz2_hist->SetStats(0);
       yvz2_hist->SetLineColor(kOrange+1);
             
   TH1F *yvx3_hist = new TH1F("Sector 3 vx", "Sector 3 vx", 300, -5, 5);
       yvx3_hist->SetStats(0);
       yvx3_hist->SetLineColor(kYellow+1);
   TH1D *yvy3_hist = new TH1D("Sector 3 vy", "Sector 3 vy", 300, -5, 5);
       yvy3_hist->SetStats(0);
       yvy3_hist->SetLineColor(kYellow-1);
   TH1D *yvz3_hist = new TH1D("Sector 3 vz", "Sector 3 vz", 300, -20, 20);
       yvz3_hist->SetStats(0);
       yvz3_hist->SetLineColor(kYellow+1);
       
   TH1F *yvx4_hist = new TH1F("Sector 4 vx", "Sector 4 vx", 300, -5, 5);
       yvx4_hist->SetStats(0);
       yvx4_hist->SetLineColor(kGreen+1);
   TH1D *yvy4_hist = new TH1D("Sector 4 vy", "Sector 4 vy", 300, -5, 5);
       yvy4_hist->SetStats(0);
       yvy4_hist->SetLineColor(kGreen-1);
   TH1D *yvz4_hist = new TH1D("Sector 4 vz", "Sector 4 vz", 300, -20, 20);
       yvz4_hist->SetStats(0);
       yvz4_hist->SetLineColor(kGreen+1);
       
   TH1F *yvx5_hist = new TH1F("Sector 5 vx", "Sector 5 vx", 300, -5, 5);
       yvx5_hist->SetStats(0);
       yvx5_hist->SetLineColor(kBlue+1);
   TH1D *yvy5_hist = new TH1D("Sector 5 vy", "Sector 5 vy", 300, -5, 5);
       yvy5_hist->SetStats(0);
       yvy5_hist->SetLineColor(kBlue-1);
   TH1D *yvz5_hist = new TH1D("Sector 5 vz", "Sector 5 vz", 300, -20, 20);
       yvz5_hist->SetStats(0);
       yvz5_hist->SetLineColor(kBlue+1);
       
   TH1F *yvx6_hist = new TH1F("Sector 6 vx", "Sector 6 vx", 300, -5, 5);
       yvx6_hist->SetStats(0);
       yvx6_hist->SetLineColor(kMagenta+1);
   TH1D *yvy6_hist = new TH1D("Sector 6 vy", "Sector 6 vy", 300, -5, 5);
       yvy6_hist->SetStats(0);
       yvy6_hist->SetLineColor(kMagenta-1);
   TH1D *yvz6_hist = new TH1D("Sector 6 vz", "Sector 6 vz", 300, -20, 20);
       yvz6_hist->SetStats(0);
       yvz6_hist->SetLineColor(kMagenta+1);

  

      TH1F *rvx1_hist = (TH1F*)input.Get("sector 1 vx"); 
      TH1D *rvy1_hist = (TH1D*)input.Get("sector 1 vy");
      TH1D *rvz1_hist = (TH1D*)input.Get("sector 1 vz");
      
      TH1F *rvx2_hist = (TH1F*)input.Get("sector 2 vx"); 
      TH1D *rvy2_hist = (TH1D*)input.Get("sector 2 vy");
      TH1D *rvz2_hist = (TH1D*)input.Get("sector 2 vz");
      
      TH1F *rvx3_hist = (TH1F*)input.Get("sector 3 vx"); 
      TH1D *rvy3_hist = (TH1D*)input.Get("sector 3 vy");
      TH1D *rvz3_hist = (TH1D*)input.Get("sector 3 vz");
      
      TH1F *rvx4_hist = (TH1F*)input.Get("sector 4 vx"); 
      TH1D *rvy4_hist = (TH1D*)input.Get("sector 4 vy");
      TH1D *rvz4_hist = (TH1D*)input.Get("sector 4 vz");
      
      TH1F *rvx5_hist = (TH1F*)input.Get("sector 5 vx"); 
      TH1D *rvy5_hist = (TH1D*)input.Get("sector 5 vy");
      TH1D *rvz5_hist = (TH1D*)input.Get("sector 5 vz");
      
      TH1F *rvx6_hist = (TH1F*)input.Get("sector 6 vx"); 
      TH1D *rvy6_hist = (TH1D*)input.Get("sector 6 vy");
      TH1D *rvz6_hist = (TH1D*)input.Get("sector 6 vz");
      
      
      vx1_hist->Add(rvx1_hist);
      vx1_hist->Scale(1.0 / vx1_hist->Integral());
      vx1_hist->SetMinimum(0.0);
      vx1_hist->SetMaximum(0.35);
      vy1_hist->Add(rvy1_hist);
      vy1_hist->Scale(1.0 / vy1_hist->Integral());
      vy1_hist->SetMinimum(0.0);
      vy1_hist->SetMaximum(0.15);
      vz1_hist->Add(rvz1_hist);
      vz1_hist->Scale(1.0 / vz1_hist->Integral());
      vz1_hist->SetMinimum(0.0);
      vz1_hist->SetMaximum(0.05);
      
      vx2_hist->Add(rvx2_hist);
      vx2_hist->Scale(1.0 / vx2_hist->Integral());
      vx2_hist->SetMinimum(0.0);
      vx2_hist->SetMaximum(0.35);
      vy2_hist->Add(rvy2_hist);
      vy2_hist->Scale(1.0 / vy2_hist->Integral());
      vy2_hist->SetMinimum(0.0);
      vy2_hist->SetMaximum(0.15);
      vz2_hist->Add(rvz2_hist);
      vz2_hist->Scale(1.0 / vz2_hist->Integral());
      vz2_hist->SetMinimum(0.0);
      vz2_hist->SetMaximum(0.05);
      
      vx3_hist->Add(rvx3_hist);
      vx3_hist->Scale(1.0 / vx3_hist->Integral());
      vx3_hist->SetMinimum(0.0);
      vx3_hist->SetMaximum(0.35);
      vy3_hist->Add(rvy3_hist);
      vy3_hist->Scale(1.0 / vy3_hist->Integral());
      vy3_hist->SetMinimum(0.0);
      vy3_hist->SetMaximum(0.15);
      vz3_hist->Add(rvz3_hist);
      vz3_hist->Scale(1.0 / vz3_hist->Integral());
      vz3_hist->SetMinimum(0.0);
      vz3_hist->SetMaximum(0.05);
      
      vx4_hist->Add(rvx4_hist);
      vx4_hist->Scale(1.0 / vx4_hist->Integral());
      vx4_hist->SetMinimum(0.0);
      vx4_hist->SetMaximum(0.35);
      vy4_hist->Add(rvy4_hist);
      vy4_hist->Scale(1.0 / vy4_hist->Integral());
      vy4_hist->SetMinimum(0.0);
      vy4_hist->SetMaximum(0.15);
      vz4_hist->Add(rvz4_hist);
      vz4_hist->Scale(1.0 / vz4_hist->Integral());
      vz4_hist->SetMinimum(0.0);
      vz4_hist->SetMaximum(0.05);
      
      vx5_hist->Add(rvx5_hist);
      vx5_hist->Scale(1.0 / vx5_hist->Integral());
      vx5_hist->SetMinimum(0.0);
      vx5_hist->SetMaximum(0.35);
      vy5_hist->Add(rvy5_hist);
      vy5_hist->Scale(1.0 / vy5_hist->Integral());
      vy5_hist->SetMinimum(0.0);
      vy5_hist->SetMaximum(0.15);
      vz5_hist->Add(rvz5_hist);
      vz5_hist->Scale(1.0 / vz5_hist->Integral());
      vz5_hist->SetMinimum(0.0);
      vz5_hist->SetMaximum(0.05);
      
      vx6_hist->Add(rvx6_hist);
      vx6_hist->Scale(1.0 / vx6_hist->Integral());
      vx6_hist->SetMinimum(0.0);
      vx6_hist->SetMaximum(0.35);
      vy6_hist->Add(rvy6_hist);
      vy6_hist->Scale(1.0 / vy6_hist->Integral());
      vy6_hist->SetMinimum(0.0);
      vy6_hist->SetMaximum(0.15);
      vz6_hist->Add(rvz6_hist);
      vz6_hist->Scale(1.0 / vz6_hist->Integral());
      vz6_hist->SetMinimum(0.0);
      vz6_hist->SetMaximum(0.05);
      
      yvx1_hist->Add(rvx1_hist);
 //     yvx1_hist->SetMinimum(0.0);
 //     yvx1_hist->SetMaximum(0.35);
      yvy1_hist->Add(rvy1_hist);
 //     yvy1_hist->SetMinimum(0.0);
 //     yvy1_hist->SetMaximum(0.15);
      yvz1_hist->Add(rvz1_hist);
 //     yvz1_hist->SetMinimum(0.0);
 //     yvz1_hist->SetMaximum(0.05);
      
      yvx2_hist->Add(rvx2_hist);
 //     yvx2_hist->SetMinimum(0.0);
 //     yvx2_hist->SetMaximum(0.35);
      yvy2_hist->Add(rvy2_hist);
 //     yvy2_hist->SetMinimum(0.0);
 //     yvy2_hist->SetMaximum(0.15);
      yvz2_hist->Add(rvz2_hist);
 //     yvz2_hist->SetMinimum(0.0);
 //     yvz2_hist->SetMaximum(0.05);
      
      yvx3_hist->Add(rvx3_hist);
 //     yvx3_hist->SetMinimum(0.0);
 //     yvx3_hist->SetMaximum(0.35);
      yvy3_hist->Add(rvy3_hist);
 //     yvy3_hist->SetMinimum(0.0);
 //     yvy3_hist->SetMaximum(0.15);
      yvz3_hist->Add(rvz3_hist);
 //     yvz3_hist->SetMinimum(0.0);
 //     yvz3_hist->SetMaximum(0.05);
      
      yvx4_hist->Add(rvx4_hist);
 //     yvx4_hist->SetMinimum(0.0);
  //    yvx4_hist->SetMaximum(0.35);
      yvy4_hist->Add(rvy4_hist);
  //    yvy4_hist->SetMinimum(0.0);
  //    yvy4_hist->SetMaximum(0.15);
      yvz4_hist->Add(rvz4_hist);
  //    yvz4_hist->SetMinimum(0.0);
  //    yvz4_hist->SetMaximum(0.05);
      
      yvx5_hist->Add(rvx5_hist);
//      yvx5_hist->SetMinimum(0.0);
//      yvx5_hist->SetMaximum(0.35);
      yvy5_hist->Add(rvy5_hist);
//      yvy5_hist->SetMinimum(0.0);
//      yvy5_hist->SetMaximum(0.15);
      yvz5_hist->Add(rvz5_hist);
//      yvz5_hist->SetMinimum(0.0);
//      yvz5_hist->SetMaximum(0.05);
      
      yvx6_hist->Add(rvx6_hist);
//      yvx6_hist->SetMinimum(0.0);
//      yvx6_hist->SetMaximum(0.35);
      yvy6_hist->Add(rvy6_hist);
//      yvy6_hist->SetMinimum(0.0);
//      yvy6_hist->SetMaximum(0.15);
      yvz6_hist->Add(rvz6_hist);
     // yvz6_hist->SetMinimum(0.0);
    //  yvz6_hist->SetMaximum(0.05);
      
      
      TH1F *tvx1_hist = new TH1F("Sector 1 vx", "Sector 1 vx", 300, -5, 5);
       tvx1_hist->SetStats(0);
       tvx1_hist->SetLineColor(kBlue+1);
   TH1D *tvy1_hist = new TH1D("Sector 1 vy", "Sector 1 vy", 300, -5, 5);
       tvy1_hist->SetStats(0);
       tvy1_hist->SetLineColor(kBlue+1);
   TH1D *tvz1_hist = new TH1D("Sector 1 vz", "Sector 1 vz", 300, -20, 20);
       tvz1_hist->SetStats(0);
       tvz1_hist->SetLineColor(kBlue+1);
       
  TH1F *tvx2_hist = new TH1F("Sector 2 vx", "Sector 2 vx", 300, -5, 5);
       tvx2_hist->SetStats(0);
       tvx2_hist->SetLineColor(kBlue+1);
   TH1D *tvy2_hist = new TH1D("Sector 2 vy", "Sector 2 vy", 300, -5, 5);
       tvy2_hist->SetStats(0);
       tvy2_hist->SetLineColor(kBlue+1);
   TH1D *tvz2_hist = new TH1D("Sector 2 vz", "Sector 2 vz", 300, -20, 20);
       tvz2_hist->SetStats(0);
       tvz2_hist->SetLineColor(kBlue+1);
             
   TH1F *tvx3_hist = new TH1F("Sector 3 vx", "Sector 3 vx", 300, -5, 5);
       tvx3_hist->SetStats(0);
       tvx3_hist->SetLineColor(kBlue+1);
   TH1D *tvy3_hist = new TH1D("Sector 3 vy", "Sector 3 vy", 300, -5, 5);
       tvy3_hist->SetStats(0);
       tvy3_hist->SetLineColor(kBlue+1);
   TH1D *tvz3_hist = new TH1D("Sector 3 vz", "Sector 3 vz", 300, -20, 20);
       tvz3_hist->SetStats(0);
       tvz3_hist->SetLineColor(kBlue+1);
       
   TH1F *tvx4_hist = new TH1F("Sector 4 vx", "Sector 4 vx", 300, -5, 5);
       tvx4_hist->SetStats(0);
       tvx4_hist->SetLineColor(kBlue+1);
   TH1D *tvy4_hist = new TH1D("Sector 4 vy", "Sector 4 vy", 300, -5, 5);
       tvy4_hist->SetStats(0);
       tvy4_hist->SetLineColor(kBlue+1);
   TH1D *tvz4_hist = new TH1D("Sector 4 vz", "Sector 4 vz", 300, -20, 20);
       tvz4_hist->SetStats(0);
       tvz4_hist->SetLineColor(kBlue+1);
       
   TH1F *tvx5_hist = new TH1F("Sector 5 vx", "Sector 5 vx", 300, -5, 5);
       tvx5_hist->SetStats(0);
       tvx5_hist->SetLineColor(kBlue+1);
   TH1D *tvy5_hist = new TH1D("Sector 5 vy", "Sector 5 vy", 300, -5, 5);
       tvy5_hist->SetStats(0);
       tvy5_hist->SetLineColor(kBlue+1);
   TH1D *tvz5_hist = new TH1D("Sector 5 vz", "Sector 5 vz", 300, -20, 20);
       tvz5_hist->SetStats(0);
       tvz5_hist->SetLineColor(kBlue+1);
       
   TH1F *tvx6_hist = new TH1F("Sector 6 vx", "Sector 6 vx", 300, -5, 5);
       tvx6_hist->SetStats(0);
       tvx6_hist->SetLineColor(kBlue+1);
   TH1D *tvy6_hist = new TH1D("Sector 6 vy", "Sector 6 vy", 300, -5, 5);
       tvy6_hist->SetStats(0);
       tvy6_hist->SetLineColor(kBlue+1);
   TH1D *tvz6_hist = new TH1D("Sector 6 vz", "Sector 6 vz", 300, -20, 20);
       tvz6_hist->SetStats(0);
       tvz6_hist->SetLineColor(kBlue+1);

  

      TH1F *ttvx1_hist = (TH1F*)input2.Get("sector 1 vx"); 
      TH1D *ttvy1_hist = (TH1D*)input2.Get("sector 1 vy");
      TH1D *ttvz1_hist = (TH1D*)input2.Get("sector 1 vz");
      
      TH1F *ttvx2_hist = (TH1F*)input2.Get("sector 2 vx"); 
      TH1D *ttvy2_hist = (TH1D*)input2.Get("sector 2 vy");
      TH1D *ttvz2_hist = (TH1D*)input2.Get("sector 2 vz");
      
      TH1F *ttvx3_hist = (TH1F*)input2.Get("sector 3 vx"); 
      TH1D *ttvy3_hist = (TH1D*)input2.Get("sector 3 vy");
      TH1D *ttvz3_hist = (TH1D*)input2.Get("sector 3 vz");
      
      TH1F *ttvx4_hist = (TH1F*)input2.Get("sector 4 vx"); 
      TH1D *ttvy4_hist = (TH1D*)input2.Get("sector 4 vy");
      TH1D *ttvz4_hist = (TH1D*)input2.Get("sector 4 vz");
      
      TH1F *ttvx5_hist = (TH1F*)input2.Get("sector 5 vx"); 
      TH1D *ttvy5_hist = (TH1D*)input2.Get("sector 5 vy");
      TH1D *ttvz5_hist = (TH1D*)input2.Get("sector 5 vz");
      
      TH1F *ttvx6_hist = (TH1F*)input2.Get("sector 6 vx"); 
      TH1D *ttvy6_hist = (TH1D*)input2.Get("sector 6 vy");
      TH1D *ttvz6_hist = (TH1D*)input2.Get("sector 6 vz");
      
      
      tvx1_hist->Add(ttvx1_hist);
      tvx1_hist->Scale(1.0 / tvx1_hist->Integral());
      tvx1_hist->SetMinimum(0.0);
      tvx1_hist->SetMaximum(0.35);
      tvy1_hist->Add(ttvy1_hist);
      tvy1_hist->Scale(1.0 / tvy1_hist->Integral());
      tvy1_hist->SetMinimum(0.0);
      tvy1_hist->SetMaximum(0.15);
      tvz1_hist->Add(ttvz1_hist);
      tvz1_hist->Scale(1.0 / tvz1_hist->Integral());
      tvz1_hist->SetMinimum(0.0);
      tvz1_hist->SetMaximum(0.05);
      
      tvx2_hist->Add(ttvx2_hist);
      tvx2_hist->Scale(1.0 / tvx2_hist->Integral());
      tvx2_hist->SetMinimum(0.0);
      tvx2_hist->SetMaximum(0.35);
      tvy2_hist->Add(ttvy2_hist);
      tvy2_hist->Scale(1.0 / tvy2_hist->Integral());
      tvy2_hist->SetMinimum(0.0);
      tvy2_hist->SetMaximum(0.15);
      tvz2_hist->Add(ttvz2_hist);
      tvz2_hist->Scale(1.0 / tvz2_hist->Integral());
      tvz2_hist->SetMinimum(0.0);
      tvz2_hist->SetMaximum(0.05);
      
      tvx3_hist->Add(ttvx3_hist);
      tvx3_hist->Scale(1.0 / tvx3_hist->Integral());
      tvx3_hist->SetMinimum(0.0);
      tvx3_hist->SetMaximum(0.35);
      tvy3_hist->Add(ttvy3_hist);
      tvy3_hist->Scale(1.0 / tvy3_hist->Integral());
      tvy3_hist->SetMinimum(0.0);
      tvy3_hist->SetMaximum(0.15);
      tvz3_hist->Add(ttvz3_hist);
      tvz3_hist->Scale(1.0 / tvz3_hist->Integral());
      tvz3_hist->SetMinimum(0.0);
      tvz3_hist->SetMaximum(0.05);
      
      tvx4_hist->Add(ttvx4_hist);
      tvx4_hist->Scale(1.0 / tvx4_hist->Integral());
      tvx4_hist->SetMinimum(0.0);
      tvx4_hist->SetMaximum(0.35);
      tvy4_hist->Add(ttvy4_hist);
      tvy4_hist->Scale(1.0 / tvy4_hist->Integral());
      tvy4_hist->SetMinimum(0.0);
      tvy4_hist->SetMaximum(0.15);
      tvz4_hist->Add(ttvz4_hist);
      tvz4_hist->Scale(1.0 / tvz4_hist->Integral());
      tvz4_hist->SetMinimum(0.0);
      tvz4_hist->SetMaximum(0.05);
      
      tvx5_hist->Add(ttvx5_hist);
      tvx5_hist->Scale(1.0 / tvx5_hist->Integral());
      tvx5_hist->SetMinimum(0.0);
      tvx5_hist->SetMaximum(0.35);
      tvy5_hist->Add(ttvy5_hist);
      tvy5_hist->Scale(1.0 / tvy5_hist->Integral());
      tvy5_hist->SetMinimum(0.0);
      tvy5_hist->SetMaximum(0.15);
      tvz5_hist->Add(ttvz5_hist);
      tvz5_hist->Scale(1.0 / tvz5_hist->Integral());
      tvz5_hist->SetMinimum(0.0);
      tvz5_hist->SetMaximum(0.05);
      
      tvx6_hist->Add(ttvx6_hist);
      tvx6_hist->Scale(1.0 / tvx6_hist->Integral());
      tvx6_hist->SetMinimum(0.0);
      tvx6_hist->SetMaximum(0.35);
      tvy6_hist->Add(ttvy6_hist);
      tvy6_hist->Scale(1.0 / tvy6_hist->Integral());
      tvy6_hist->SetMinimum(0.0);
      tvy6_hist->SetMaximum(0.15);
      tvz6_hist->Add(ttvz6_hist);
      tvz6_hist->Scale(1.0 / tvz6_hist->Integral());
      tvz6_hist->SetMinimum(0.0);
      tvz6_hist->SetMaximum(0.05);
      
      
      
      //gStyle->SetCanvasPreferGL();
      TCanvas *c1 = new TCanvas("c1","vx",0,0,640,480);
      vx1_hist->Draw("HIST");
      vx2_hist->Draw("HIST Same");
      vx3_hist->Draw("HIST Same");
      vx4_hist->Draw("HIST Same");
      vx5_hist->Draw("HIST Same");
      vx6_hist->Draw("HIST Same");
      tvx1_hist->Draw("HIST Same");
      tvx2_hist->Draw("HIST Same");
      tvx3_hist->Draw("HIST Same");
      tvx4_hist->Draw("HIST Same");
      tvx5_hist->Draw("HIST Same");
      tvx6_hist->Draw("HIST Same");
      c1->Write();
      
      TCanvas *c2 = new TCanvas("c2","vy",0,0,640,480);
      vy1_hist->Draw("HIST");
      vy2_hist->Draw("HIST Same");
      vy3_hist->Draw("HIST Same");
      vy4_hist->Draw("HIST Same");
      vy5_hist->Draw("HIST Same");
      vy6_hist->Draw("HIST Same");
      tvy1_hist->Draw("HIST Same");
      tvy2_hist->Draw("HIST Same");
      tvy3_hist->Draw("HIST Same");
      tvy4_hist->Draw("HIST Same");
      tvy5_hist->Draw("HIST Same");
      tvy6_hist->Draw("HIST Same");
      c2->Write();
      
      
      TCanvas *c3 = new TCanvas("c3","vz",0,0,640,480);
      vz1_hist->Draw("HIST");
      vz2_hist->Draw("HIST Same");
      vz3_hist->Draw("HIST Same");
      vz4_hist->Draw("HIST Same");
      vz5_hist->Draw("HIST Same");
      vz6_hist->Draw("HIST Same");
      tvz1_hist->Draw("HIST Same");
      tvz2_hist->Draw("HIST Same");
      tvz3_hist->Draw("HIST Same");
      tvz4_hist->Draw("HIST Same");
      tvz5_hist->Draw("HIST Same");
      tvz6_hist->Draw("HIST Same");
      c3->Write();
      
      TCanvas *c4 = new TCanvas("c4","Sector 1",0,0,640,480);
      yvx1_hist->Draw("HIST");
      yvy1_hist->Draw("HIST Same");
      c4->Write();
      
      TCanvas *c5 = new TCanvas("c5","Sector 2",0,0,640,480);
      yvy2_hist->Draw("HIST");
      yvx2_hist->Draw("HIST Same");
      c5->Write();
      
      TCanvas *c6 = new TCanvas("c6","Sector 3",0,0,640,480);
      yvy3_hist->Draw("HIST");
      yvx3_hist->Draw("HIST Same");
      c6->Write();
      
      TCanvas *c7 = new TCanvas("c7","Sector 4",0,0,640,480);
      yvx4_hist->Draw("HIST");
      yvy4_hist->Draw("HIST Same");
      c7->Write();
      
      TCanvas *c8 = new TCanvas("c8","Sector 5",0,0,640,480);
      yvy5_hist->Draw("HIST");
      yvx5_hist->Draw("HIST Same");
      c8->Write();
      
      TCanvas *c9 = new TCanvas("c9","Sector 6",0,0,640,480);
      yvy6_hist->Draw("HIST");
      yvx6_hist->Draw("HIST Same");
      c9->Write();
      
     
  


   
  
 


}


