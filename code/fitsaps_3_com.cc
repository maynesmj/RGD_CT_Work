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

void fitsaps_3_com(){


   TFile input("cu_pass0v7_ghp_022.root", "read");    ///fits_ld2.root    cusn_kin5_dnp_4_sigma_tzS.root  cusn_kin5_pass0v5_dnp_200_no_sigma_tzS.root
   TFile input22("cusn_pass0v11_gcom_good_ibnt_03.root", "read");
   TFile input2("cusn_pass0v11_gcom_bad_ibnt_03.root", "read");    ///fits_ld2.root
   TFile input3("ld_pass0v11_gcom_bad_ibnt_03.root", "read");    ///fits_ld2.root
   TFile input4("cxc_pass0v11_ccom_obnt_02.root", "read");    ///fits_ld2.root
    TFile f("sol_pass0v10_obnt_24_fits.root", "update");
   
   gStyle->SetOptStat("");
	gStyle->SetCanvasPreferGL();
/*	
   TH1F *vx_hist = new TH1F("vx", "All Sectors vx", 200, -5, 5);
       //vx_hist->SetStats(0);
       vx_hist->SetLineColor(kBlack);
       vx_hist->GetXaxis()->SetTitle("vx [cm]");
   TH1D *vy_hist = new TH1D("vy", "All Sectors vy", 200, -5, 5);
       //vy_hist->SetStats(0);
       vy_hist->SetLineColor(kBlack);
       vy_hist->GetXaxis()->SetTitle("vy [cm]");
   TH1D *vz_hist = new TH1D("vz", "All Sectors vz", 200, -30, 20);
       //vz_hist->SetStats(0);
       vz_hist->SetLineColor(kBlack);
       vz_hist->GetXaxis()->SetTitle("vz [cm]");

   TH1F *vx1_hist = new TH1F("Sector 1 vx", "Sector 1 vx", 200, -5, 5);
      // vx1_hist->SetStats(0);
       vx1_hist->SetLineColor(kCyan+1);
   TH1D *vy1_hist = new TH1D("Sector 1 vy", "Sector 1 vy", 200, -5, 5);
     //  vy1_hist->SetStats(0);
       vy1_hist->SetLineColor(kCyan+1);
   TH1D *vz1_hist = new TH1D("Sector 1 vz", "Sector 1 vz", 200, -30, 20);
    //   vz1_hist->SetStats(0);
       vz1_hist->SetLineColor(kWhite);
       
  TH1F *vx2_hist = new TH1F("Sector 2 vx", "Sector 2 vx", 200, -5, 5);
   //    vx2_hist->SetStats(0);
       vx2_hist->SetLineColor(kOrange + 2+1);
   TH1D *vy2_hist = new TH1D("Sector 2 vy", "Sector 2 vy", 200, -5, 5);
    //   vy2_hist->SetStats(0);
       vy2_hist->SetLineColor(kOrange + 2+1);
   TH1D *vz2_hist = new TH1D("Sector 2 vz", "Sector 2 vz", 200, -30, 20);
    //   vz2_hist->SetStats(0);
       vz2_hist->SetLineColor(kOrange + 2+1);
             
   TH1F *vx3_hist = new TH1F("Sector 3 vx", "Sector 3 vx", 200, -5, 5);
    //   vx3_hist->SetStats(0);
       vx3_hist->SetLineColor(kYellow+1);
   TH1D *vy3_hist = new TH1D("Sector 3 vy", "Sector 3 vy", 200, -5, 5);
    //   vy3_hist->SetStats(0);
       vy3_hist->SetLineColor(kYellow+1);
   TH1D *vz3_hist = new TH1D("Sector 3 vz", "Sector 3 vz", 200, -30, 20);
   //    vz3_hist->SetStats(0);
       vz3_hist->SetLineColor(kYellow+1);
       
   TH1F *vx4_hist = new TH1F("Sector 4 vx", "Sector 4 vx", 200, -5, 5);
   //    vx4_hist->SetStats(0);
       vx4_hist->SetLineColor(kGreen + 2+1);
   TH1D *vy4_hist = new TH1D("Sector 4 vy", "Sector 4 vy", 200, -5, 5);
   //    vy4_hist->SetStats(0);
       vy4_hist->SetLineColor(kGreen + 2+1);
   TH1D *vz4_hist = new TH1D("Sector 4 vz", "Sector 4 vz", 200, -30, 20);
  //     vz4_hist->SetStats(0);
       vz4_hist->SetLineColor(kGreen + 2+1);
       
   TH1F *vx5_hist = new TH1F("Sector 5 vx", "Sector 5 vx", 200, -5, 5);
   //    vx5_hist->SetStats(0);
       vx5_hist->SetLineColor(kBlue + 2+1);
   TH1D *vy5_hist = new TH1D("Sector 5 vy", "Sector 5 vy", 200, -5, 5);
    //   vy5_hist->SetStats(0);
       vy5_hist->SetLineColor(kBlue + 2+1);
   TH1D *vz5_hist = new TH1D("Sector 5 vz", "Sector 5 vz", 200, -30, 20);
    //   vz5_hist->SetStats(0);
       vz5_hist->SetLineColor(kBlue + 2+1);
       
   TH1F *vx6_hist = new TH1F("Sector 6 vx", "Sector 6 vx", 200, -5, 5);
    //   vx6_hist->SetStats(0);
       vx6_hist->SetLineColor(kMagenta + 2+1);
   TH1D *vy6_hist = new TH1D("Sector 6 vy", "Sector 6 vy", 200, -5, 5);
   //    vy6_hist->SetStats(0);
       vy6_hist->SetLineColor(kMagenta + 2+1);
   TH1D *vz6_hist = new TH1D("Sector 6 vz", "Sector 6 vz", 200, -30, 20);
  //     vz6_hist->SetStats(0);
       vz6_hist->SetLineColor(kMagenta + 2+1);
       
       TH1F *vxtest_hist = new TH1F("Sector 6 vx", "Sector 6 vx", 200, -5, 5);
    //   vx6_hist->SetStats(0);
       vxtest_hist->SetLineColor(kMagenta + 2+1);
   TH1D *vytest_hist = new TH1D("Sector 6 vy", "Sector 6 vy", 200, -5, 5);
   //    vy6_hist->SetStats(0);
       vytest_hist->SetLineColor(kMagenta + 2+1);
   TH1D *vztest_hist = new TH1D("Sector 6 vz", "Sector 6 vz", 200, -30, 20);
  //     vz6_hist->SetStats(0);
       vztest_hist->SetLineColor(kMagenta + 2+1);
  */     
    TH1F *vcu_hist = new TH1F("Sector 6 vx", "LD_{2} Inbending", 200, -30, 20);
       vcu_hist->SetStats(0);
       vcu_hist->SetLineColor(kOrange + 2);
       vcu_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   TH1F *vsn_hist = new TH1F("Sector 6 vy", "LD_{2} Outbending", 200, -30, 20);
       vsn_hist->SetStats(0);
       vsn_hist->SetLineColor(kMagenta + 2);
   TH1F *vcxc_hist = new TH1F("Sector 6 vx", "Good DCProblematic DC Outbending", 200, -30, 20);
       vcxc_hist->SetStats(0);
       vcxc_hist->SetLineColor(kGreen + 2);
   TH1F *vld_hist = new TH1F("Sector 6 vy", "CxC Outbending", 200, -30, 20);
       vld_hist->SetStats(0);
       vld_hist->SetLineColor(kBlue + 2);
       
       
   TH1F *cu_q2_hist = new TH1F("Good DC Q2", "Good DC", 200, 0, 10);
       cu_q2_hist->SetStats(0);
       cu_q2_hist->SetLineColor(kOrange + 2);
       cu_q2_hist->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
   TH1F *sn_q2_hist = new TH1F("Good DC Q2", "Problematic DC", 200, 0, 10);
       sn_q2_hist->SetStats(0);
       sn_q2_hist->SetLineColor(kMagenta + 2);
   TH1F *cu_w_hist = new TH1F("Good DC w", "Good DC", 200, 0, 5);
       cu_w_hist->SetStats(0);
       cu_w_hist->SetLineColor(kOrange + 2);
       cu_w_hist->GetXaxis()->SetTitle("W [GeV]");
   TH1F *sn_w_hist = new TH1F("Good DC w", "Problematic DC", 200, 0, 5);
       sn_w_hist->SetStats(0);
       sn_w_hist->SetLineColor(kMagenta + 2);
       sn_w_hist->GetXaxis()->SetTitle("W [GeV]");
   TH1F *ld_q2_hist = new TH1F("Good DC Q2", "LD_{2}", 200, 0, 10);
       ld_q2_hist->SetStats(0);
       ld_q2_hist->SetLineColor(kBlue + 2);
       ld_q2_hist->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
       //cu_w_hist->GetXaxis()->SetTitle("W [GeV]");
   TH1F *cxc_q2_hist = new TH1F("Good DC Q2", "CxC", 200, 0, 10);
       cxc_q2_hist->SetStats(0);
       cxc_q2_hist->SetLineColor(kGreen + 2);
       //cu_w_hist->GetXaxis()->SetTitle("W [GeV]");
   TH1F *ld_w_hist = new TH1F("Good DC w", "LD_{2}", 200, 0, 5);
       ld_w_hist->SetStats(0);
       ld_w_hist->SetLineColor(kBlue + 2);
       ld_w_hist->GetXaxis()->SetTitle("W [GeV]");
   TH1F *cxc_w_hist = new TH1F("Good DC w", "CxC", 200, 0, 5);
       cxc_w_hist->SetStats(0);
       cxc_w_hist->SetLineColor(kGreen + 2);
       cxc_w_hist->GetXaxis()->SetTitle("W [GeV]");
       
   TH1F *cu_lc_hist = new TH1F("Good DC bad", "Good DC", 200, 0, 3);
       cu_lc_hist->SetStats(0);
       cu_lc_hist->SetLineColor(kOrange + 2);
       cu_lc_hist->GetXaxis()->SetTitle("lc [fm]");
   TH1F *sn_lc_hist = new TH1F("Problematic DC lc", "Problematic DC", 200, 0, 3);
       sn_lc_hist->SetStats(0);
       sn_lc_hist->SetLineColor(kMagenta + 2);
       sn_lc_hist->GetXaxis()->SetTitle("lc [fm]");
   TH1F *cu_t_hist = new TH1F("Good DC t", "Good DC", 200, 0, 1);
       cu_t_hist->SetStats(0);
       cu_t_hist->SetLineColor(kOrange + 2);
       cu_t_hist->GetXaxis()->SetTitle("-t [GeV]");
   TH1F *sn_t_hist = new TH1F("Problematic DC t", "Problematic DC", 200, 0, 1);
       sn_t_hist->SetStats(0);
       sn_t_hist->SetLineColor(kMagenta + 2);
       sn_t_hist->GetXaxis()->SetTitle("-t [GeV]");
   TH1F *ld_lc_hist = new TH1F("Good DC bad", "LD_{2}", 200, 0, 3);
       ld_lc_hist->SetStats(0);
       ld_lc_hist->SetLineColor(kBlue + 2);
       ld_lc_hist->GetXaxis()->SetTitle("lc [fm]");
       ld_lc_hist->GetXaxis()->SetTitle("lc [fm]");
   TH1F *cxc_lc_hist = new TH1F("Problematic DC lc", "CxC", 200, 0, 3);
       cxc_lc_hist->SetStats(0);
       cxc_lc_hist->SetLineColor(kGreen + 2);
       cxc_lc_hist->GetXaxis()->SetTitle("lc [fm]");
   TH1F *ld_t_hist = new TH1F("Good DC t", "LD_{2}", 200, 0, 1);
       ld_t_hist->SetStats(0);
       ld_t_hist->SetLineColor(kBlue + 2);
       ld_t_hist->GetXaxis()->SetTitle("-t [GeV]");
   TH1F *cxc_t_hist = new TH1F("Problematic DC t", "CxC", 200, 0, 1);
       cxc_t_hist->SetStats(0);
       cxc_t_hist->SetLineColor(kGreen + 2);
       cxc_t_hist->GetXaxis()->SetTitle("-t [GeV]");
       
   TH1F *cu_z_hist = new TH1F("Good DC z", "Good DC", 200, 0, 1.2);
       cu_z_hist->SetStats(0);
       cu_z_hist->SetLineColor(kOrange + 2);
       cu_z_hist->GetXaxis()->SetTitle("Z_{h}");
   TH1F *sn_z_hist = new TH1F("Problematic DC z", "Problematic DC", 200, 0, 1.2);
       sn_z_hist->SetStats(0);
       sn_z_hist->SetLineColor(kMagenta + 2);
       sn_z_hist->GetXaxis()->SetTitle("Z_{h}");
   TH1F *cu_ep_hist = new TH1F("Good DC electron momentum", "Good DC", 200, 0, 11);
       cu_ep_hist->SetStats(0);
       cu_ep_hist->SetLineColor(kOrange + 2);
       cu_ep_hist->GetXaxis()->SetTitle("P_{e} [GeV]");
   TH1F *sn_ep_hist = new TH1F("Problematic DC electron momentum", "Problematic DC", 200, 0, 11);
       sn_ep_hist->SetStats(0);
       sn_ep_hist->SetLineColor(kMagenta + 2);
       sn_ep_hist->GetXaxis()->SetTitle("P_{e} [GeV]");
   TH1F *ld_z_hist = new TH1F("Good DC z", "LD_{2}", 200, 0, 1.2);
       ld_z_hist->SetStats(0);
       ld_z_hist->SetLineColor(kBlue + 2);
       ld_z_hist->GetXaxis()->SetTitle("Z_{h}");
   TH1F *cxc_z_hist = new TH1F("Problematic DC z", "CxC", 200, 0, 1.2);
       cxc_z_hist->SetStats(0);
       cxc_z_hist->SetLineColor(kGreen + 2);
       cxc_z_hist->GetXaxis()->SetTitle("Z_{h}");
   TH1F *ld_ep_hist = new TH1F("Good DC electron momentum", "LD_{2}", 200, 0, 11);
       ld_ep_hist->SetStats(0);
       ld_ep_hist->SetLineColor(kBlue + 2);
       ld_ep_hist->GetXaxis()->SetTitle("P_{e} [GeV]");
   TH1F *cxc_ep_hist = new TH1F("Problematic DC electron momentum", "CxC", 200, 0, 11);
       cxc_ep_hist->SetStats(0);
       cxc_ep_hist->SetLineColor(kGreen + 2);
       cxc_ep_hist->GetXaxis()->SetTitle("P_{e} [GeV]");
       
   TH1F *cu_pipp_hist = new TH1F("Good DC pip momentum", "Good DC", 200, 0, 10);
       cu_pipp_hist->SetStats(0);
       cu_pipp_hist->SetLineColor(kOrange + 2);
       cu_pipp_hist->GetXaxis()->SetTitle("P_{#pi^{+}} [GeV]");
   TH1F *sn_pipp_hist = new TH1F("Problematic DC pip momentum", "Problematic DC", 200, 0, 10);
       sn_pipp_hist->SetStats(0);
       sn_pipp_hist->SetLineColor(kMagenta + 2);
       sn_pipp_hist->GetXaxis()->SetTitle("P_{#pi^{+}} [GeV]");
   TH1F *cu_pimp_hist = new TH1F("Good DC pim momentum", "Good DC", 200, 0, 10);
       cu_pimp_hist->SetStats(0);
       cu_pimp_hist->SetLineColor(kOrange + 2);
       cu_pimp_hist->GetXaxis()->SetTitle("P_{#pi^{-}} [GeV]");
   TH1F *sn_pimp_hist = new TH1F("Problematic DC pim momentum", "Problematic DC", 200, 0, 10);
       sn_pimp_hist->SetStats(0);
       sn_pimp_hist->SetLineColor(kMagenta + 2);
       sn_pimp_hist->GetXaxis()->SetTitle("P_{#pi^{-}} [GeV]");
   TH1F *ld_pipp_hist = new TH1F("Good DC pip momentum", "LD_{2}", 200, 0, 10);
       ld_pipp_hist->SetStats(0);
       ld_pipp_hist->SetLineColor(kBlue + 2);
       ld_pipp_hist->GetXaxis()->SetTitle("P_{#pi^{+}} [GeV]");
   TH1F *cxc_pipp_hist = new TH1F("Problematic DC pip momentum", "CxC", 200, 0, 10);
       cxc_pipp_hist->SetStats(0);
       cxc_pipp_hist->SetLineColor(kGreen + 2);
       cxc_pipp_hist->GetXaxis()->SetTitle("P_{#pi^{+}} [GeV]");
   TH1F *ld_pimp_hist = new TH1F("Good DC pim momentum", "LD_{2}", 200, 0, 10);
       ld_pimp_hist->SetStats(0);
       ld_pimp_hist->SetLineColor(kBlue + 2);
       ld_pimp_hist->GetXaxis()->SetTitle("P_{#pi^{-}} [GeV]");
   TH1F *cxc_pimp_hist = new TH1F("Problematic DC pim momentum", "CxC", 200, 0, 10);
       cxc_pimp_hist->SetStats(0);
       cxc_pimp_hist->SetLineColor(kGreen + 2);
       cxc_pimp_hist->GetXaxis()->SetTitle("P_{#pi^{-}} [GeV]");
       
       
   TH1F *cu_eph_hist = new TH1F("Good DC electron phi", "Good DC", 200, -180, 180);
       cu_eph_hist->SetStats(0);
       cu_eph_hist->SetLineColor(kOrange + 2);
       cu_eph_hist->GetXaxis()->SetTitle("#phi_{e} [Degree]");
   TH1F *sn_eph_hist = new TH1F("Problematic DC electron phi", "Problematic DC", 200, -180, 180);
       sn_eph_hist->SetStats(0);
       sn_eph_hist->SetLineColor(kMagenta + 2);
       sn_eph_hist->GetXaxis()->SetTitle("#phi_{e} [Degree]");
       
   TH1F *cu_pipph_hist = new TH1F("Good DC pip phi", "Good DC", 200, 0, 40);
       cu_pipph_hist->SetStats(0);
       cu_pipph_hist->SetLineColor(kOrange + 2);
       cu_pipph_hist->GetXaxis()->SetTitle("#theta_{e} [Degree]");
       //cu_pipph_hist->GetXaxis()->SetTitle("#pi^{+} #phi [Degree]");
   TH1F *sn_pipph_hist = new TH1F("Problematic DC pip phi", "Problematic DC", 200, 0, 40);
       sn_pipph_hist->SetStats(0);
       sn_pipph_hist->SetLineColor(kMagenta + 2);
       sn_pipph_hist->GetXaxis()->SetTitle("#theta_{e} [Degree]");
   TH1F *cu_pimph_hist = new TH1F("Good DC pim phi", "Good DC", 200, -180, 180);
       cu_pimph_hist->SetStats(0);
       cu_pimph_hist->SetLineColor(kOrange + 2);
       cu_pimph_hist->GetXaxis()->SetTitle("#pi^{-} #phi [Degree]");
   TH1F *sn_pimph_hist = new TH1F("Problematic DC pim phi", "Problematic DC", 200, -180, 180);
       sn_pimph_hist->SetStats(0);
       sn_pimph_hist->SetLineColor(kMagenta + 2);
       
       TH1F *cu_mass_hist = new TH1F("Good DC mass", "Good DC", 200, 0.2, 1.4);
       cu_mass_hist->SetStats(0);
       cu_mass_hist->SetLineColor(kOrange + 2);
       cu_mass_hist->GetXaxis()->SetTitle("M_{#rho} [GeV]");
   TH1F *sn_mass_hist = new TH1F("Problematic DC mass", "Problematic DC", 200, 0.2, 1.4);
       sn_mass_hist->SetStats(0);
       sn_mass_hist->SetLineColor(kMagenta + 2);
       
       
   TH1F *ld_eph_hist = new TH1F("Good DC electron phi", "LD_{2}", 200, -180, 180);
       ld_eph_hist->SetStats(0);
       ld_eph_hist->SetLineColor(kBlue + 2);
       ld_eph_hist->GetXaxis()->SetTitle("#phi_{e} [Degree]");
   TH1F *cxc_eph_hist = new TH1F("Problematic DC electron phi", "CxC", 200, -180, 180);
       cxc_eph_hist->SetStats(0);
       cxc_eph_hist->SetLineColor(kGreen + 2);
       cxc_eph_hist->GetXaxis()->SetTitle("#phi_{e} [Degree]");
       
   TH1F *ld_pipph_hist = new TH1F("Good DC pip phi", "LD_{2}", 200, 0, 40);
       ld_pipph_hist->SetStats(0);
       ld_pipph_hist->SetLineColor(kBlue + 2);
       ld_pipph_hist->GetXaxis()->SetTitle("#theta_{e} [Degree]");
   TH1F *cxc_pipph_hist = new TH1F("Problematic DC pip phi", "CxC", 200, 0, 40);
       cxc_pipph_hist->SetStats(0);
       cxc_pipph_hist->SetLineColor(kGreen + 2);
       cxc_pipph_hist->GetXaxis()->SetTitle("#theta_{e} [Degree]");
   TH1F *ld_pimph_hist = new TH1F("Good DC pim phi", "LD_{2}", 200, -180, 180);
       ld_pimph_hist->SetStats(0);
       ld_pimph_hist->SetLineColor(kBlue + 2);
       ld_pimph_hist->GetXaxis()->SetTitle("#pi^{-} #phi [Degree]");
   TH1F *cxc_pimph_hist = new TH1F("Problematic DC pim phi", "CxC", 200, -180, 180);
       cxc_pimph_hist->SetStats(0);
       cxc_pimph_hist->SetLineColor(kGreen + 2);
       
       TH1F *ld_mass_hist = new TH1F("Good DC mass", "LD_{2}", 200, 0.2, 1.4);
       ld_mass_hist->SetStats(0);
       ld_mass_hist->SetLineColor(kBlue + 2);
       ld_mass_hist->GetXaxis()->SetTitle("M_{#rho} [GeV]");
   TH1F *cxc_mass_hist = new TH1F("Problematic DC mass", "CxC", 200, 0.2, 1.4);
       cxc_mass_hist->SetStats(0);
       cxc_mass_hist->SetLineColor(kGreen + 2);




  /*

      TH1F *tvx1_hist = (TH1F*)input.Get("sector 1 vx"); 
      TH1D *tvy1_hist = (TH1D*)input.Get("sector 1 vy");
      TH1D *tvz1_hist = (TH1D*)input.Get("sector 1 vz");
      
      TH1F *tvx2_hist = (TH1F*)input.Get("sector 2 vx"); 
      TH1D *tvy2_hist = (TH1D*)input.Get("sector 2 vy");
      TH1D *tvz2_hist = (TH1D*)input.Get("sector 2 vz");
      
      TH1F *tvx3_hist = (TH1F*)input.Get("sector 3 vx"); 
      TH1D *tvy3_hist = (TH1D*)input.Get("sector 3 vy");
      TH1D *tvz3_hist = (TH1D*)input.Get("sector 3 vz");
      
      TH1F *tvx4_hist = (TH1F*)input.Get("sector 4 vx"); 
      TH1D *tvy4_hist = (TH1D*)input.Get("sector 4 vy");
      TH1D *tvz4_hist = (TH1D*)input.Get("sector 4 vz");
      
      TH1F *tvx5_hist = (TH1F*)input.Get("sector 5 vx"); 
      TH1D *tvy5_hist = (TH1D*)input.Get("sector 5 vy");
      TH1D *tvz5_hist = (TH1D*)input.Get("sector 5 vz");
      
      TH1F *tvx6_hist = (TH1F*)input.Get("sector 6 vx"); 
      TH1D *tvy6_hist = (TH1D*)input.Get("sector 6 vy");
      TH1D *tvz6_hist = (TH1D*)input.Get("sector 6 vz");
  */    
      TH1F *cu_hist = (TH1F*)input22.Get(" vz"); 
      TH1F *sn_hist = (TH1F*)input2.Get(" vz");
      TH1F *cxc_hist = (TH1F*)input3.Get(" vz");
      TH1F *ld_hist = (TH1F*)input4.Get(" vz"); 
      
      TH1F *vcu_q2_hist = (TH1F*)input22.Get("Q2"); 
      TH1F *vsn_q2_hist = (TH1F*)input2.Get("Q2");
      TH1F *vld_q2_hist = (TH1F*)input3.Get("Q2"); 
      TH1F *vcxc_q2_hist = (TH1F*)input4.Get("Q2");
      
      TH1F *vcu_w_hist = (TH1F*)input22.Get("w (Invariant Mass)"); 
      TH1F *vsn_w_hist = (TH1F*)input2.Get("w (Invariant Mass)");
      TH1F *vld_w_hist = (TH1F*)input3.Get("w (Invariant Mass)"); 
      TH1F *vcxc_w_hist = (TH1F*)input4.Get("w (Invariant Mass)");
      
      TH1F *vcu_lc_hist = (TH1F*)input22.Get("lc"); 
      TH1F *vsn_lc_hist = (TH1F*)input2.Get("lc");
      TH1F *vld_lc_hist = (TH1F*)input3.Get("lc"); 
      TH1F *vcxc_lc_hist = (TH1F*)input4.Get("lc");
      
      TH1F *vcu_t_hist = (TH1F*)input22.Get("-t"); 
      TH1F *vsn_t_hist = (TH1F*)input2.Get("-t");
      TH1F *vld_t_hist = (TH1F*)input3.Get("-t"); 
      TH1F *vcxc_t_hist = (TH1F*)input4.Get("-t");
      
      TH1F *vcu_z_hist = (TH1F*)input22.Get("Zh"); 
      TH1F *vsn_z_hist = (TH1F*)input2.Get("Zh");
      TH1F *vld_z_hist = (TH1F*)input3.Get("Zh"); 
      TH1F *vcxc_z_hist = (TH1F*)input4.Get("Zh");
      
      TH1F *vcu_ep_hist = (TH1F*)input22.Get("Electron Momentum"); 
      TH1F *vsn_ep_hist = (TH1F*)input2.Get("Electron Momentum");
      TH1F *vld_ep_hist = (TH1F*)input3.Get("Electron Momentum"); 
      TH1F *vcxc_ep_hist = (TH1F*)input4.Get("Electron Momentum");
      
      TH1F *vcu_pipp_hist = (TH1F*)input22.Get("#Pi^{+} Momentum"); 
      TH1F *vsn_pipp_hist = (TH1F*)input2.Get("#Pi^{+} Momentum");
      TH1F *vld_pipp_hist = (TH1F*)input3.Get("#Pi^{+} Momentum"); 
      TH1F *vcxc_pipp_hist = (TH1F*)input4.Get("#Pi^{+} Momentum");
      
      TH1F *vcu_pimp_hist = (TH1F*)input22.Get("#Pi^{-} Momentum"); 
      TH1F *vsn_pimp_hist = (TH1F*)input2.Get("#Pi^{-} Momentum");
      TH1F *vld_pimp_hist = (TH1F*)input3.Get("#Pi^{-} Momentum"); 
      TH1F *vcxc_pimp_hist = (TH1F*)input4.Get("#Pi^{-} Momentum");
      
      TH1F *vcu_eph_hist = (TH1F*)input22.Get("e phi"); 
      TH1F *vsn_eph_hist = (TH1F*)input2.Get("e phi");
      TH1F *vld_eph_hist = (TH1F*)input3.Get("e phi"); 
      TH1F *vcxc_eph_hist = (TH1F*)input4.Get("e phi");
      
      TH1F *vcu_pipph_hist = (TH1F*)input22.Get("e phi vs theta"); 
      TH1F *vsn_pipph_hist = (TH1F*)input2.Get("e phi vs theta");
      TH1F *vld_pipph_hist = (TH1F*)input3.Get("e phi vs theta"); 
      TH1F *vcxc_pipph_hist = (TH1F*)input4.Get("e phi vs theta");
      
      TH1F *vcu_pimph_hist = (TH1F*)input22.Get("pim phi"); 
      TH1F *vsn_pimph_hist = (TH1F*)input2.Get("pim phi");
      TH1F *vld_pimph_hist = (TH1F*)input3.Get("pim phi"); 
      TH1F *vcxc_pimph_hist = (TH1F*)input4.Get("pim phi");
      
      TH1F *vcu_mass_hist = (TH1F*)input22.Get("With W, -t & z Good DCt"); 
      TH1F *vsn_mass_hist = (TH1F*)input2.Get("With W, -t & z Good DCt");
      TH1F *vld_mass_hist = (TH1F*)input3.Get("With W, -t & z Good DCt"); 
      TH1F *vcxc_mass_hist = (TH1F*)input4.Get("With W, -t & z Good DCt");
      
      cu_q2_hist->Add(vcu_q2_hist);
      sn_q2_hist->Add(vsn_q2_hist);
      
      cu_w_hist->Add(vcu_w_hist);
      sn_w_hist->Add(vsn_w_hist);
      
      cu_lc_hist->Add(vcu_lc_hist);
      sn_lc_hist->Add(vsn_lc_hist);
      
      cu_t_hist->Add(vcu_t_hist);
      sn_t_hist->Add(vsn_t_hist);
      
      cu_z_hist->Add(vcu_z_hist);
      sn_z_hist->Add(vsn_z_hist);
      
      cu_ep_hist->Add(vcu_ep_hist);
      sn_ep_hist->Add(vsn_ep_hist);
      
      cu_pipp_hist->Add(vcu_pipp_hist);
      sn_pipp_hist->Add(vsn_pipp_hist);
      
      cu_pimp_hist->Add(vcu_pimp_hist);
      sn_pimp_hist->Add(vsn_pimp_hist);
      
      cu_eph_hist->Add(vcu_eph_hist);
      sn_eph_hist->Add(vsn_eph_hist);
      
      cu_pipph_hist->Add(vcu_pipph_hist);
      sn_pipph_hist->Add(vsn_pipph_hist);
      
      cu_pimph_hist->Add(vcu_pimph_hist);
      sn_pimph_hist->Add(vsn_pimph_hist);
      
      cu_mass_hist->Add(vcu_mass_hist);
      sn_mass_hist->Add(vsn_mass_hist);
      
      ld_q2_hist->Add(vld_q2_hist);
      cxc_q2_hist->Add(vcxc_q2_hist);
      
      ld_w_hist->Add(vld_w_hist);
      cxc_w_hist->Add(vcxc_w_hist);
      
      ld_lc_hist->Add(vld_lc_hist);
      cxc_lc_hist->Add(vcxc_lc_hist);
      
      ld_t_hist->Add(vld_t_hist);
      cxc_t_hist->Add(vcxc_t_hist);
      
      ld_z_hist->Add(vld_z_hist);
      cxc_z_hist->Add(vcxc_z_hist);
      
      ld_ep_hist->Add(vld_ep_hist);
      cxc_ep_hist->Add(vcxc_ep_hist);
      
      ld_pipp_hist->Add(vld_pipp_hist);
      cxc_pipp_hist->Add(vcxc_pipp_hist);
      
      ld_pimp_hist->Add(vld_pimp_hist);
      cxc_pimp_hist->Add(vcxc_pimp_hist);
      
      ld_eph_hist->Add(vld_eph_hist);
      cxc_eph_hist->Add(vcxc_eph_hist);
      
      ld_pipph_hist->Add(vld_pipph_hist);
      cxc_pipph_hist->Add(vcxc_pipph_hist);
      
      ld_pimph_hist->Add(vld_pimph_hist);
      cxc_pimph_hist->Add(vcxc_pimph_hist);
      
      ld_mass_hist->Add(vld_mass_hist);
      cxc_mass_hist->Add(vcxc_mass_hist);
      
      cu_q2_hist->Scale(1.0/ cu_q2_hist->Integral("width"));
      sn_q2_hist->Scale(1.0/ sn_q2_hist->Integral("width"));
      
      cu_w_hist->Scale(1.0/ cu_w_hist->Integral("width"));
      sn_w_hist->Scale(1.0/ sn_w_hist->Integral("width"));
      
      cu_lc_hist->Scale(1.0/ cu_lc_hist->Integral("width"));
      sn_lc_hist->Scale(1.0/ sn_lc_hist->Integral("width"));
      
      cu_t_hist->Scale(1.0/ cu_t_hist->Integral("width"));
      sn_t_hist->Scale(1.0/ sn_t_hist->Integral("width"));
      
      cu_z_hist->Scale(1.0/ cu_z_hist->Integral("width"));
      sn_z_hist->Scale(1.0/ sn_z_hist->Integral("width"));
      
      cu_ep_hist->Scale(1.0/ cu_ep_hist->Integral("width"));
      sn_ep_hist->Scale(1.0/ sn_ep_hist->Integral("width"));
      
      cu_pipp_hist->Scale(1.0/ cu_pipp_hist->Integral("width"));
      sn_pipp_hist->Scale(1.0/ sn_pipp_hist->Integral("width"));
      
      cu_pimp_hist->Scale(1.0/ cu_pimp_hist->Integral("width"));
      sn_pimp_hist->Scale(1.0/ sn_pimp_hist->Integral("width"));
      
      cu_eph_hist->Scale(1.0/ cu_eph_hist->Integral("width"));
      sn_eph_hist->Scale(1.0/ sn_eph_hist->Integral("width"));
      
      cu_pipph_hist->Scale(1.0/ cu_pipph_hist->Integral("width"));
      sn_pipph_hist->Scale(1.0/ sn_pipph_hist->Integral("width"));
      
      cu_pimph_hist->Scale(1.0/ cu_pimph_hist->Integral("width"));
      sn_pimph_hist->Scale(1.0/ sn_pimph_hist->Integral("width"));
      
      cu_mass_hist->Scale(1.0/ cu_mass_hist->Integral("width"));
      sn_mass_hist->Scale(1.0/ sn_mass_hist->Integral("width"));
      
      
      ld_q2_hist->Scale(1.0/ ld_q2_hist->Integral("width"));
      cxc_q2_hist->Scale(1.0/ cxc_q2_hist->Integral("width"));
      
      ld_w_hist->Scale(1.0/ ld_w_hist->Integral("width"));
      cxc_w_hist->Scale(1.0/ cxc_w_hist->Integral("width"));
      
      ld_lc_hist->Scale(1.0/ ld_lc_hist->Integral("width"));
      cxc_lc_hist->Scale(1.0/ cxc_lc_hist->Integral("width"));
      
      ld_t_hist->Scale(1.0/ ld_t_hist->Integral("width"));
      cxc_t_hist->Scale(1.0/ cxc_t_hist->Integral("width"));
      
      ld_z_hist->Scale(1.0/ ld_z_hist->Integral("width"));
      cxc_z_hist->Scale(1.0/ cxc_z_hist->Integral("width"));
      
      ld_ep_hist->Scale(1.0/ ld_ep_hist->Integral("width"));
      cxc_ep_hist->Scale(1.0/ cxc_ep_hist->Integral("width"));
      
      ld_pipp_hist->Scale(1.0/ ld_pipp_hist->Integral("width"));
      cxc_pipp_hist->Scale(1.0/ cxc_pipp_hist->Integral("width"));
      
      ld_pimp_hist->Scale(1.0/ ld_pimp_hist->Integral("width"));
      cxc_pimp_hist->Scale(1.0/ cxc_pimp_hist->Integral("width"));
      
      ld_eph_hist->Scale(1.0/ ld_eph_hist->Integral("width"));
      cxc_eph_hist->Scale(1.0/ cxc_eph_hist->Integral("width"));
      
      ld_pipph_hist->Scale(1.0/ ld_pipph_hist->Integral("width"));
      cxc_pipph_hist->Scale(1.0/ cxc_pipph_hist->Integral("width"));
      
      ld_pimph_hist->Scale(1.0/ ld_pimph_hist->Integral("width"));
      cxc_pimph_hist->Scale(1.0/ cxc_pimph_hist->Integral("width"));
      
      ld_mass_hist->Scale(1.0/ ld_mass_hist->Integral("width"));
      cxc_mass_hist->Scale(1.0/ cxc_mass_hist->Integral("width"));
      
      
      
      vcu_hist->Add(cu_hist);
      vcu_hist->Scale(1.0);                        //1.0 / vcu_hist->Integral("width"));
      vcu_hist->SetMinimum(0.0);
     // vcu_hist->SetMaximum(50);
      vcu_hist->SetLineColor(kGreen + 2);
      vcu_hist->SetLineWidth(3);
      vsn_hist->Add(sn_hist);
      vsn_hist->Scale(1.0);/// vsn_hist->Integral("width"));
      vsn_hist->SetMinimum(0.0);
     // vsn_hist->SetMaximum(50);
      vsn_hist->SetLineColor(kBlue + 2);
      vsn_hist->SetLineWidth(3);
      vcxc_hist->Add(cxc_hist);
      vcxc_hist->Scale(1.0/vcxc_hist->Integral("width")); // / vcxc_hist->Integral("width"));
      vcxc_hist->SetMinimum(0.0);
     // vcxc_hist->SetMaximum(50);
      vcxc_hist->SetLineColor(kGreen + 2);
      vcxc_hist->SetLineWidth(3);
      vld_hist->Add(ld_hist);
      vld_hist->Scale(1.0/vld_hist->Integral("width"));   /// vld_hist->Integral("width"));
      vld_hist->SetMinimum(0.0);
    //  vld_hist->SetMaximum(50);
      vld_hist->SetLineColor(kBlue + 2);
      vld_hist->SetLineWidth(3);
      
    /* 
      vx_hist->Add(tvx1_hist);
      vx_hist->Add(tvx2_hist);
      vx_hist->Add(tvx3_hist);
      vx_hist->Add(tvx4_hist);
      vx_hist->Add(tvx5_hist);
      vx_hist->Add(tvx6_hist);
      vx_hist->Scale(1.0 / vx_hist->Integral());
      vx_hist->SetMinimum(0.0);
      vx_hist->SetMaximum(0.19);
      vy_hist->Add(tvy1_hist);
      vy_hist->Add(tvy2_hist);
      vy_hist->Add(tvy3_hist);
      vy_hist->Add(tvy4_hist);
      vy_hist->Add(tvy5_hist);
      vy_hist->Add(tvy6_hist);
      vy_hist->Scale(1.0 / vy_hist->Integral());
      vy_hist->SetMinimum(0.0);
      vy_hist->SetMaximum(0.07);
      vz_hist->Add(tvz1_hist);
      vz_hist->Add(tvz2_hist);
      vz_hist->Add(tvz3_hist);
      vz_hist->Add(tvz4_hist);
      vz_hist->Add(tvz5_hist);
      vz_hist->Add(tvz6_hist);
      vz_hist->Scale(1.0 / vz_hist->Integral());
      vz_hist->SetMinimum(0.0);
      vz_hist->SetMaximum(0.05);
      
      
      
      vx1_hist->Add(tvx1_hist);
      vx1_hist->Scale(1.0 / vx1_hist->Integral());
      vx1_hist->SetMinimum(0.0);
      vx1_hist->SetMaximum(0.19);
      vy1_hist->Add(tvy1_hist);
      vy1_hist->Scale(1.0 / vy1_hist->Integral());
      vy1_hist->SetMinimum(0.0);
      vy1_hist->SetMaximum(0.07);
      vz1_hist->Add(tvz1_hist);
      vz1_hist->Scale(1.0 / vz1_hist->Integral());
      vz1_hist->SetMinimum(0.0);
      vz1_hist->SetMaximum(0.05);
      
      vx2_hist->Add(tvx2_hist);
      vx2_hist->Scale(1.0 / vx2_hist->Integral());
      vx2_hist->SetMinimum(0.0);
      vx2_hist->SetMaximum(0.19);
      vy2_hist->Add(tvy2_hist);
      vy2_hist->Scale(1.0 / vy2_hist->Integral());
      vy2_hist->SetMinimum(0.0);
      vy2_hist->SetMaximum(0.07);
      vz2_hist->Add(tvz2_hist);
      vz2_hist->Scale(1.0 / vz2_hist->Integral());
      vz2_hist->SetMinimum(0.0);
      vz2_hist->SetMaximum(0.05);
      
      vx3_hist->Add(tvx3_hist);
      vx3_hist->Scale(1.0 / vx3_hist->Integral());
      vx3_hist->SetMinimum(0.0);
      vx3_hist->SetMaximum(0.19);
      vy3_hist->Add(tvy3_hist);
      vy3_hist->Scale(1.0 / vy3_hist->Integral());
      vy3_hist->SetMinimum(0.0);
      vy3_hist->SetMaximum(0.07);
      vz3_hist->Add(tvz3_hist);
      vz3_hist->Scale(1.0 / vz3_hist->Integral());
      vz3_hist->SetMinimum(0.0);
      vz3_hist->SetMaximum(0.05);
      
      vx4_hist->Add(tvx4_hist);
      vx4_hist->Scale(1.0 / vx4_hist->Integral());
      vx4_hist->SetMinimum(0.0);
      vx4_hist->SetMaximum(0.19);
      vy4_hist->Add(tvy4_hist);
      vy4_hist->Scale(1.0 / vy4_hist->Integral());
      vy4_hist->SetMinimum(0.0);
      vy4_hist->SetMaximum(0.07);
      vz4_hist->Add(tvz4_hist);
      vz4_hist->Scale(1.0 / vz4_hist->Integral());
      vz4_hist->SetMinimum(0.0);
      vz4_hist->SetMaximum(0.07);
      
      vx5_hist->Add(tvx5_hist);
      vx5_hist->Scale(1.0 / vx5_hist->Integral());
      vx5_hist->SetMinimum(0.0);
      vx5_hist->SetMaximum(0.19);
      vy5_hist->Add(tvy5_hist);
      vy5_hist->Scale(1.0 / vy5_hist->Integral());
      vy5_hist->SetMinimum(0.0);
      vy5_hist->SetMaximum(0.07);
      vz5_hist->Add(tvz5_hist);
      vz5_hist->Scale(1.0 / vz5_hist->Integral());
      vz5_hist->SetMinimum(0.0);
      vz5_hist->SetMaximum(0.05);
      
      vx6_hist->Add(tvx6_hist);
      vx6_hist->Scale(1.0 / vx6_hist->Integral());
      vx6_hist->SetMinimum(0.0);
      vx6_hist->SetMaximum(0.19);
      vy6_hist->Add(tvy6_hist);
      vy6_hist->Scale(1.0 / vy6_hist->Integral());
      vy6_hist->SetMinimum(0.0);
      vy6_hist->SetMaximum(0.07);
      vz6_hist->Add(tvz6_hist);
      vz6_hist->Scale(1.0 / vz6_hist->Integral());
      vz6_hist->SetMinimum(0.0);
      vz6_hist->SetMaximum(0.05);
      */
      /*
      This section was already committed out
      vx_hist->Add(vx1_hist);
      vx_hist->Add(vx2_hist);
      vx_hist->Add(vx3_hist);
      vx_hist->Add(vx4_hist);
      vx_hist->Add(vx5_hist);
      vx_hist->Add(vx6_hist);
      //vx_hist->Scale(1.0 / vx_hist->Integral());
      vx_hist->SetMinimum(0.0);
      vx_hist->SetMaximum(0.19);
      vy_hist->Add(vy1_hist);
      vy_hist->Add(vy2_hist);
      vy_hist->Add(vy3_hist);
      vy_hist->Add(vy4_hist);
      vy_hist->Add(vy5_hist);
      vy_hist->Add(vy6_hist);
      //vy_hist->Scale(1.0 / vy_hist->Integral());
      vy_hist->SetMinimum(0.0);
      vy_hist->SetMaximum(0.19);
      vz_hist->Add(vz1_hist);
      vz_hist->Add(vz2_hist);
      vz_hist->Add(vz3_hist);
      vz_hist->Add(vz4_hist);
      vz_hist->Add(vz5_hist);
      vz_hist->Add(vz6_hist);
      //vz_hist->Scale(1.0 / vz_hist->Integral());
      vz_hist->SetMinimum(0.0);
      vz_hist->SetMaximum(0.05);
      */
      
      
      /*
      //float mul = 200000000000000000000000000000000000000000.0; // 3.0 4.0 5.0
                 float ilxlim = 0.00656862 - 3*(0.3567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float iuxlim = 0.00656862 + 3*(0.3567);  //1.5741; // 0.9525 1.2633 1.5741
                 float ilylim = -0.1793 - 3*(0.33594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float iuylim = -0.1793 + 3*(0.33594); 
                 float iilxlim = 0.00656862 - 4*(0.3567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float iiuxlim = 0.00656862 + 4*(0.3567);  //1.5741; // 0.9525 1.2633 1.5741
                 float iilylim = -0.1793 - 4*(0.33594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float iiuylim = -0.1793 + 4*(0.33594); 
                 float iiilxlim = 0.00656862 - 5*(0.3567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float iiiuxlim = 0.00656862 + 5*(0.3567);  //1.5741; // 0.9525 1.2633 1.5741
                 float iiilylim = -0.1793 - 5*(0.33594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float iiiuylim = -0.1793 + 5*(0.33594); 
      //gStyle->SetCanvasPreferGL();
      TCanvas *c1 = new TCanvas("c1","vx",0,0,640,480);
      std::cout << "max height " << c1->GetUymax() << std::endl;
      TLine *ilx=new TLine(ilxlim,c1->GetUymin(),ilxlim,c1->GetUymax()-0.2);
      ilx->SetLineColor(kBlue + 2);
      TLine *iux=new TLine(iuxlim,c1->GetUymin(),iuxlim,c1->GetUymax()-0.2);
      iux->SetLineColor(kBlue + 2);
      TLine *iilx=new TLine(iilxlim,c1->GetUymin(),iilxlim,c1->GetUymax()-0.2);
      iilx->SetLineColor(kGreen + 2);
      TLine *iiux=new TLine(iiuxlim,c1->GetUymin(),iiuxlim,c1->GetUymax()-0.2);
      iiux->SetLineColor(kGreen + 2);
      TLine *iiilx=new TLine(iiilxlim,c1->GetUymin(),iiilxlim,c1->GetUymax()-0.2);
      iiilx->SetLineColor(kMagenta + 2);
      TLine *iiiux=new TLine(iiiuxlim,c1->GetUymin(),iiiuxlim,c1->GetUymax()-0.2);
      iiiux->SetLineColor(kMagenta + 2);
      
      vx_hist->Draw("HIST");
      vx1_hist->Draw("HIST Same");
      vx2_hist->Draw("HIST Same");
      vx3_hist->Draw("HIST Same");
      vx4_hist->Draw("HIST Same");
      vx5_hist->Draw("HIST Same");
      vx6_hist->Draw("HIST Same");
      c1->BuildLegend();
      ilx->Draw("HIST Same");
      iux->Draw("HIST Same");
      iilx->Draw("HIST Same");
      iiux->Draw("HIST Same");
      iiilx->Draw("HIST Same");
      iiiux->Draw("HIST Same");
      c1->Write();
      
      TCanvas *c2 = new TCanvas("c2","vy",0,0,640,480);
      //std::cout << "max height " << c2->GetUymax() << std::endl;
      TLine *ily=new TLine(ilylim,c2->GetUymin(),ilylim,c2->GetUymax()-0.2);
      ily->SetLineColor(kBlue + 2);
      TLine *iuy=new TLine(iuylim,c2->GetUymin(),iuylim,c2->GetUymax()-0.2);
      iuy->SetLineColor(kBlue + 2);
      TLine *iily=new TLine(iilylim,c2->GetUymin(),iilylim,c2->GetUymax()-0.3);
      iily->SetLineColor(kGreen + 2);
      TLine *iiuy=new TLine(iiuylim,c2->GetUymin(),iiuylim,c2->GetUymax()-0.3);
      iiuy->SetLineColor(kGreen + 2);
      TLine *iiily=new TLine(iiilylim,c2->GetUymin(),iiilylim,c2->GetUymax()-0.4);
      iiily->SetLineColor(kMagenta + 2);
      TLine *iiiuy=new TLine(iiiuylim,c2->GetUymin(),iiiuylim,c2->GetUymax()-0.4);
      iiiuy->SetLineColor(kMagenta + 2);
      vy_hist->Draw("HIST");
      vy1_hist->Draw("HIST Same");
      vy2_hist->Draw("HIST Same");
      vy3_hist->Draw("HIST Same");
      vy4_hist->Draw("HIST Same");
      vy5_hist->Draw("HIST Same");
      vy6_hist->Draw("HIST Same");
      c2->BuildLegend();
      ily->Draw("HIST Same");
      iuy->Draw("HIST Same");
      iily->Draw("HIST Same");
      iiuy->Draw("HIST Same");
      iiily->Draw("HIST Same");
      iiiuy->Draw("HIST Same");
      c2->Write();
      
      
      TCanvas *c3 = new TCanvas("c3","vz",0,0,640,480);
      vz_hist->Draw("HIST");
      vz1_hist->Draw("HIST Same");
      vz2_hist->Draw("HIST Same");
      vz3_hist->Draw("HIST Same");
      vz4_hist->Draw("HIST Same");
      vz5_hist->Draw("HIST Same");
      vz6_hist->Draw("HIST Same");
      c3->BuildLegend();
      c3->Write();
      */
      TCanvas *c4 = new TCanvas("Invariant Mass","Invariant Mass",0,0,640,480);
      //vcu_hist->Draw("HIST");
      //vsn_hist->Draw("HIST Same");
      vcxc_hist->Draw("HIST");
      vld_hist->Draw("HIST Same");
      c4->BuildLegend();
      c4->Write();
      
       TCanvas *cq2 = new TCanvas("Q2c","Q2",0,0,640,480);
      cu_q2_hist->Draw("HIST");
      sn_q2_hist->Draw("HIST Same");
      //ld_q2_hist->Draw("HIST Same");
      //cxc_q2_hist->Draw("HIST Same");
      cq2->BuildLegend();
      cq2->Write();
      
             TCanvas *cw = new TCanvas("Wc","W",0,0,640,480);
      //cxc_w_hist->Draw("HIST");
      sn_w_hist->Draw("HIST Same");
      //ld_w_hist->Draw("HIST Same");
      cu_w_hist->Draw("HIST Same");
      cw->BuildLegend();
      cw->Write();
      
             TCanvas *clc = new TCanvas("lcc","lc",0,0,640,480);
      cu_lc_hist->Draw("HIST");
      sn_lc_hist->Draw("HIST Same");
      //ld_lc_hist->Draw("HIST Same");
      //cxc_lc_hist->Draw("HIST Same");
      clc->BuildLegend();
      clc->Write();
      
             TCanvas *ct = new TCanvas("-tc","-t",0,0,640,480);
      cu_t_hist->Draw("HIST");
      sn_t_hist->Draw("HIST Same");
      //ld_t_hist->Draw("HIST Same");
      //cxc_t_hist->Draw("HIST Same");
      ct->BuildLegend();
      ct->Write();
      
             TCanvas *cz = new TCanvas("Zh","Zh",0,0,640,480);
      //cxc_z_hist->Draw("HIST");
      sn_z_hist->Draw("HIST Same");
      //ld_z_hist->Draw("HIST Same");
      cu_z_hist->Draw("HIST Same");
      cz->BuildLegend();
      cz->Write();
      
             TCanvas *cep = new TCanvas("cep","ep",0,0,640,480);
      //cxc_ep_hist->Draw("HIST");
      sn_ep_hist->Draw("HIST Same");
      //ld_ep_hist->Draw("HIST Same");
      cu_ep_hist->Draw("HIST Same");
      cep->BuildLegend();
      cep->Write();
      
             TCanvas *cpipp = new TCanvas("cpipp","pipp",0,0,640,480);
      sn_pipp_hist->Draw("HIST");
      //cxc_pipp_hist->Draw("HIST Same");
      //ld_pipp_hist->Draw("HIST Same");
      cu_pipp_hist->Draw("HIST Same");
      cpipp->BuildLegend();
      cpipp->Write();
      
             TCanvas *cpimp = new TCanvas("cpimp","pimp",0,0,640,480);
      sn_pimp_hist->Draw("HIST");
      cu_pimp_hist->Draw("HIST Same");
      //ld_pimp_hist->Draw("HIST Same");
      //cxc_pimp_hist->Draw("HIST Same");
      cpimp->BuildLegend();
      cpimp->Write();
      
      
             TCanvas *ceph = new TCanvas("ceph","eph",0,0,640,480);
      cu_eph_hist->Draw("HIST");
      sn_eph_hist->Draw("HIST Same");
      //ld_eph_hist->Draw("HIST Same");
      //cxc_eph_hist->Draw("HIST Same");
      ceph->BuildLegend();
      ceph->Write();
      
             TCanvas *cpipph = new TCanvas("cpipph","pipph",0,0,640,480);
      cu_pipph_hist->Draw("HIST");
      sn_pipph_hist->Draw("HIST Same");
      //ld_pipph_hist->Draw("HIST Same");
      //cxc_pipph_hist->Draw("HIST Same");
      cpipph->BuildLegend();
      cpipph->Write();
      
             TCanvas *cpimph = new TCanvas("cpimph","pimph",0,0,640,480);
      cu_pimph_hist->Draw("HIST");
      sn_pimph_hist->Draw("HIST Same");
      //ld_pimph_hist->Draw("HIST Same");
      //cxc_pimph_hist->Draw("HIST Same");
      cpimph->BuildLegend();
      cpimph->Write();
      
             TCanvas *cmass = new TCanvas("cmass","mass",0,0,640,480);
      cu_mass_hist->Draw("HIST");
      sn_mass_hist->Draw("HIST Same");
      //ld_mass_hist->Draw("HIST Same");
      //cxc_mass_hist->Draw("HIST Same");
      cmass->BuildLegend();
      cmass->Write();
      /*
      TCanvas *c5 = new TCanvas("c5","vx",0,0,640,480);
      vx_hist->Draw("HIST");
      vxtest_hist->Draw("HIST Same");
      //vx_hist->Write();
      c5->Write();
      
      TCanvas *c6 = new TCanvas("c6","vy",0,0,640,480);
      vy_hist->Draw("HIST");
      vytest_hist->Draw("HIST Same");
      //vy_hist->Write();
      c6->Write();
      
      TCanvas *c7 = new TCanvas("c7","vz",0,0,640,480);
      vz_hist->Draw("HIST");
      vztest_hist->Draw("HIST Same");
      //vz_hist->Write();
      c7->Write();
      
      
      vx_hist->Write();
      vy_hist->Write();
      
       TF1 *g1 = new TF1 ("m1", "gaus", -10.0, -5);
  g1->SetLineColor(kRed);
  TF1 *g2 = new TF1 ("m2", "gaus", -5, -0.0);
  g2->SetLineColor(kGreen + 2);
  TF1 *f1 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f1->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f1->SetLineColor(kBlue + 2);
  
  gStyle->SetOptFit(1);
  
  vz_hist->Fit(g1, "R");
  vz_hist->Fit(g2, "R+");
  
  Double_t par[8];
  g1->GetParameters(&par[0]);
  g2->GetParameters(&par[3]);
  f1->SetParameters(par);
  
  vz_hist->Fit(f1, "R+");
  // h->Fit(f1, "+", "e1", 4, 16);
  vz_hist->Draw("e1");
      vz_hist->Write();
     
  TF1 *g11 = new TF1 ("m11", "gaus", -10.0, -5);
  g11->SetLineColor(kRed);
  TF1 *g21 = new TF1 ("m21", "gaus", -5, -0.0);
  g21->SetLineColor(kGreen + 2);
  TF1 *f11 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f11->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f11->SetLineColor(kBlue + 2);
  
  gStyle->SetOptFit(1);
  
  tvz1_hist->Fit(g11, "R");
  tvz1_hist->Fit(g21, "R+");
  
  Double_t par1[8];
  g11->GetParameters(&par1[0]);
  g21->GetParameters(&par1[3]);
  f11->SetParameters(par1);
  
  tvz1_hist->Fit(f11, "R+");
  // h->Fit(f1, "+", "e1", 4, 16);
  tvz1_hist->Draw("e1");
      tvz1_hist->Write();
      
      
   TF1 *g12 = new TF1 ("m12", "gaus", -10.0, -5);
  g12->SetLineColor(kRed);
  TF1 *g22 = new TF1 ("m22", "gaus", -5, -0.0);
  g22->SetLineColor(kGreen + 2);
  TF1 *f12 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f12->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f12->SetLineColor(kBlue + 2);
  
  gStyle->SetOptFit(1);
  
  tvz2_hist->Fit(g12, "R");
  tvz2_hist->Fit(g22, "R+");
  
  Double_t par2[8];
  g12->GetParameters(&par2[0]);
  g22->GetParameters(&par2[3]);
  f12->SetParameters(par2);
  
  tvz2_hist->Fit(f12, "R+");
  // h->Fit(f1, "+", "e1", 4, 16);
  tvz2_hist->Draw("e1");
      tvz2_hist->Write();


  
  
  TF1 *g13 = new TF1 ("m13", "gaus", -10.0, -5);
  g13->SetLineColor(kRed);
  TF1 *g23 = new TF1 ("m23", "gaus", -5, -0.0);
  g23->SetLineColor(kGreen + 2);
  TF1 *f13 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f13->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f13->SetLineColor(kBlue + 2);
  
  gStyle->SetOptFit(1);
  
  tvz3_hist->Fit(g13, "R");
  tvz3_hist->Fit(g23, "R+");
  
  Double_t par3[8];
  g13->GetParameters(&par3[0]);
  g23->GetParameters(&par3[3]);
  f13->SetParameters(par3);
  
  tvz3_hist->Fit(f13, "R+");
  tvz3_hist->Draw("e1");
  tvz3_hist->Write(); 
  
  
  
  TF1 *g14 = new TF1 ("m14", "gaus", -10.0, -5);
  g14->SetLineColor(kRed);
  TF1 *g24 = new TF1 ("m24", "gaus", -5, -0.0);
  g24->SetLineColor(kGreen + 2);
  TF1 *f14 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f14->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f14->SetLineColor(kBlue + 2);
  
  gStyle->SetOptFit(1);
  
  tvz4_hist->Fit(g14, "R");
  tvz4_hist->Fit(g24, "R+");
  
  Double_t par4[8];
  g14->GetParameters(&par4[0]);
  g24->GetParameters(&par4[3]);
  f14->SetParameters(par4);
  
  tvz4_hist->Fit(f14, "R+");
  tvz4_hist->Draw("e1");
  tvz4_hist->Write(); 
  
  
  
  
  TF1 *g15 = new TF1 ("m15", "gaus", -10.0, -5);
  g15->SetLineColor(kRed);
  TF1 *g25 = new TF1 ("m25", "gaus", -5, -0.0);
  g25->SetLineColor(kGreen + 2);
  TF1 *f15 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f15->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f15->SetLineColor(kBlue + 2);
  
  gStyle->SetOptFit(1);
  
  tvz5_hist->Fit(g15, "R");
  tvz5_hist->Fit(g25, "R+");
  
  Double_t par5[8];
  g15->GetParameters(&par5[0]);
  g25->GetParameters(&par5[3]);
  f15->SetParameters(par5);
  
  tvz5_hist->Fit(f15, "R+");
  tvz5_hist->Draw("e1");
  tvz5_hist->Write(); 
 


  TF1 *g16 = new TF1 ("m16", "gaus", -10.0, -5);
  g16->SetLineColor(kRed);
  TF1 *g26 = new TF1 ("m23", "gaus", -5, -0.0);
  g26->SetLineColor(kGreen + 2);
  TF1 *f16 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f16->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f16->SetLineColor(kBlue + 2);
  
  gStyle->SetOptFit(1);
  
  tvz6_hist->Fit(g16, "R");
  tvz6_hist->Fit(g26, "R+");
  
  Double_t par6[8];
  g16->GetParameters(&par6[0]);
  g26->GetParameters(&par6[3]);
  f16->SetParameters(par6);
  
  tvz6_hist->Fit(f16, "R+");
  tvz6_hist->Draw("e1");
  tvz6_hist->Write(); 
  */

}


