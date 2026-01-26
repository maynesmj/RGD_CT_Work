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

void fitsaps_3_comm(){


   TFile input("cusn_pass0v11_zzgg_obnt_27.root", "read"); 
   TFile input2("sn_pass0v11_zzgg_obnt_09.root", "read");    ///fits_ld2.root
   TFile input3("ld_pass0v11_zzgg_obnt_27.root", "read");    ///fits_ld2.root
   TFile input4("cxc_pass0v11_zzgg_obnt_27.root", "read");    ///fits_ld2.root
   TFile input5("cusn_pass0v11_zzgg_ibnt_26.root", "read"); 
   TFile input6("sn_pass0v11_zzgg_ibnt_09.root", "read");    ///fits_ld2.root
   TFile input7("ld_pass0v11_zzgg_ibnt_11.root", "read");    ///fits_ld2.root
   TFile input8("cxc_pass0v11_zzgg_ibnt_26.root", "read");    ///fits_ld2.root
   TFile input9("graph_normalized_electron_yield_all.root", "read"); 
   TFile input10("graph_normalized_pip_yield_all.root", "read");    ///fits_ld2.root
   TFile input11("graph_normalized_pim_yield_all.root", "read");    ///fits_ld2.root
   TFile input12("graph_normalized_rho0_yield_all.root", "read");    ///fits_ld2.root
    TFile f("zzzzzsol_pass0v11_obnt_63_fits.root", "update");
   
   gStyle->SetOptStat("");
	gStyle->SetCanvasPreferGL();

       
       
   
      TGraph *cu_ob_e_g = (TGraph*)input.Get("te_g");
              cu_ob_e_g->SetMarkerColor(kOrange); 
      TGraph *sn_ob_e_g = (TGraph*)input2.Get("te_g");
              sn_ob_e_g->SetMarkerColor(kMagenta);
      TGraph *ld_ob_e_g = (TGraph*)input3.Get("te_g");
              ld_ob_e_g->SetMarkerColor(kBlue);
      TGraph *cxc_ob_e_g = (TGraph*)input4.Get("te_g");
              cxc_ob_e_g->SetMarkerColor(kGreen);
      TGraph *cu_ib_e_g = (TGraph*)input5.Get("te_g");
              cu_ib_e_g->SetMarkerColor(kOrange); 
      TGraph *sn_ib_e_g = (TGraph*)input6.Get("te_g");
              sn_ib_e_g->SetMarkerColor(kMagenta);
      TGraph *ld_ib_e_g = (TGraph*)input7.Get("te_g");
              ld_ib_e_g->SetMarkerColor(kBlue);
      TGraph *cxc_ib_e_g = (TGraph*)input8.Get("te_g");
              cxc_ib_e_g->SetMarkerColor(kGreen);
      
      TGraph *cu_ob_pip_g = (TGraph*)input.Get("tpip_g"); 
              cu_ob_pip_g->SetMarkerColor(kOrange);
      TGraph *sn_ob_pip_g = (TGraph*)input2.Get("tpip_g");
              sn_ob_pip_g->SetMarkerColor(kMagenta);
      TGraph *ld_ob_pip_g = (TGraph*)input3.Get("tpip_g");
              ld_ob_pip_g->SetMarkerColor(kBlue);
      TGraph *cxc_ob_pip_g = (TGraph*)input4.Get("tpip_g");
              cxc_ob_pip_g->SetMarkerColor(kGreen);
      TGraph *cu_ib_pip_g = (TGraph*)input5.Get("tpip_g");
              cu_ib_pip_g->SetMarkerColor(kOrange); 
      TGraph *sn_ib_pip_g = (TGraph*)input6.Get("tpip_g");
              sn_ib_pip_g->SetMarkerColor(kMagenta);
      TGraph *ld_ib_pip_g = (TGraph*)input7.Get("tpip_g");
              ld_ib_pip_g->SetMarkerColor(kBlue);
      TGraph *cxc_ib_pip_g = (TGraph*)input8.Get("tpip_g");
              cxc_ib_pip_g->SetMarkerColor(kGreen); 
      
      TGraph *cu_ob_pim_g = (TGraph*)input.Get("tpim_g"); 
              cu_ob_pim_g->SetMarkerColor(kOrange);
      TGraph *sn_ob_pim_g = (TGraph*)input2.Get("tpim_g");
              sn_ob_pim_g->SetMarkerColor(kMagenta);
      TGraph *ld_ob_pim_g = (TGraph*)input3.Get("tpim_g");
              ld_ob_pim_g->SetMarkerColor(kBlue);
      TGraph *cxc_ob_pim_g = (TGraph*)input4.Get("tpim_g");
              cxc_ob_pim_g->SetMarkerColor(kGreen);
      TGraph *cu_ib_pim_g = (TGraph*)input5.Get("tpim_g"); 
              cu_ib_pim_g->SetMarkerColor(kOrange);
      TGraph *sn_ib_pim_g = (TGraph*)input6.Get("tpim_g");
              sn_ib_pim_g->SetMarkerColor(kMagenta);
      TGraph *ld_ib_pim_g = (TGraph*)input7.Get("tpim_g");
              ld_ib_pim_g->SetMarkerColor(kBlue);
      TGraph *cxc_ib_pim_g = (TGraph*)input8.Get("tpim_g");  
              cxc_ib_pim_g->SetMarkerColor(kGreen);
      
      TGraph *cu_ob_r1_g = (TGraph*)input.Get("tr1_g");
              cu_ob_r1_g->SetMarkerColor(kOrange); 
      TGraph *sn_ob_r1_g = (TGraph*)input2.Get("tr1_g");
              sn_ob_r1_g->SetMarkerColor(kMagenta);
      TGraph *ld_ob_r1_g = (TGraph*)input3.Get("tr1_g");
              ld_ob_r1_g->SetMarkerColor(kBlue);
      TGraph *cxc_ob_r1_g = (TGraph*)input4.Get("tr1_g");
              cxc_ob_r1_g->SetMarkerColor(kGreen);
      TGraph *cu_ib_r1_g = (TGraph*)input5.Get("tr1_g"); 
              cu_ib_r1_g->SetMarkerColor(kOrange);
      TGraph *sn_ib_r1_g = (TGraph*)input6.Get("tr1_g");
              sn_ib_r1_g->SetMarkerColor(kMagenta);
      TGraph *ld_ib_r1_g = (TGraph*)input7.Get("tr1_g");
              ld_ib_r1_g->SetMarkerColor(kBlue);
      TGraph *cxc_ib_r1_g = (TGraph*)input8.Get("tr1_g");
              cxc_ib_r1_g->SetMarkerColor(kGreen);
      
      TGraph *cu_ob_r2_g = (TGraph*)input.Get("tr2_g");
              cu_ob_r2_g->SetMarkerColor(kOrange); 
      TGraph *sn_ob_r2_g = (TGraph*)input2.Get("tr2_g");
              sn_ob_r2_g->SetMarkerColor(kMagenta);
      TGraph *ld_ob_r2_g = (TGraph*)input3.Get("tr2_g");
              ld_ob_r2_g->SetMarkerColor(kBlue);
      TGraph *cxc_ob_r2_g = (TGraph*)input4.Get("tr2_g");
              cxc_ob_r2_g->SetMarkerColor(kGreen);
      TGraph *cu_ib_r2_g = (TGraph*)input5.Get("tr2_g"); 
              cu_ib_r2_g->SetMarkerColor(kOrange);
      TGraph *sn_ib_r2_g = (TGraph*)input6.Get("tr2_g");
              sn_ib_r2_g->SetMarkerColor(kMagenta);
      TGraph *ld_ib_r2_g = (TGraph*)input7.Get("tr2_g");
              ld_ib_r2_g->SetMarkerColor(kBlue);
      TGraph *cxc_ib_r2_g = (TGraph*)input8.Get("tr2_g");
              cxc_ib_r2_g->SetMarkerColor(kGreen); 
      
      TGraph *cu_ob_r3_g = (TGraph*)input.Get("tr3_g"); 
              cu_ob_r3_g->SetMarkerColor(kOrange);
      TGraph *sn_ob_r3_g = (TGraph*)input2.Get("tr3_g");
              sn_ob_r3_g->SetMarkerColor(kMagenta);
      TGraph *ld_ob_r3_g = (TGraph*)input3.Get("tr3_g");
              ld_ob_r3_g->SetMarkerColor(kBlue);
      TGraph *cxc_ob_r3_g = (TGraph*)input4.Get("tr3_g");
              cxc_ob_r3_g->SetMarkerColor(kGreen);
      TGraph *cu_ib_r3_g = (TGraph*)input5.Get("tr3_g"); 
              cu_ib_r3_g->SetMarkerColor(kOrange);
      TGraph *sn_ib_r3_g = (TGraph*)input6.Get("tr3_g");
              sn_ib_r3_g->SetMarkerColor(kMagenta);
      TGraph *ld_ib_r3_g = (TGraph*)input7.Get("tr3_g");
              ld_ib_r3_g->SetMarkerColor(kBlue);
      TGraph *cxc_ib_r3_g = (TGraph*)input8.Get("tr3_g");
              cxc_ib_r3_g->SetMarkerColor(kGreen);
      
      TGraph *cu_ob_r4_g = (TGraph*)input.Get("tr4_g"); 
              cu_ob_r4_g->SetMarkerColor(kOrange);
      TGraph *sn_ob_r4_g = (TGraph*)input2.Get("tr4_g");
              sn_ob_r4_g->SetMarkerColor(kMagenta);
      TGraph *ld_ob_r4_g = (TGraph*)input3.Get("tr4_g");
              ld_ob_r4_g->SetMarkerColor(kBlue);
      TGraph *cxc_ob_r4_g = (TGraph*)input4.Get("tr4_g");
              cxc_ob_r4_g->SetMarkerColor(kGreen);
      TGraph *cu_ib_r4_g = (TGraph*)input5.Get("tr4_g"); 
              cu_ib_r4_g->SetMarkerColor(kOrange);
      TGraph *sn_ib_r4_g = (TGraph*)input6.Get("tr4_g");
              sn_ib_r4_g->SetMarkerColor(kMagenta);
      TGraph *ld_ib_r4_g = (TGraph*)input7.Get("tr4_g");
              ld_ib_r4_g->SetMarkerColor(kBlue);
      TGraph *cxc_ib_r4_g = (TGraph*)input8.Get("tr4_g");
              cxc_ib_r4_g->SetMarkerColor(kGreen);
      
      TGraph *cu_ob_r5_g = (TGraph*)input.Get("tr5_g"); 
              cu_ob_r5_g->SetMarkerColor(kOrange);
      TGraph *sn_ob_r5_g = (TGraph*)input2.Get("tr5_g");
              sn_ob_r5_g->SetMarkerColor(kMagenta);
      TGraph *ld_ob_r5_g = (TGraph*)input3.Get("tr5_g");
              ld_ob_r5_g->SetMarkerColor(kBlue);
      TGraph *cxc_ob_r5_g = (TGraph*)input4.Get("tr5_g");
              cxc_ob_r5_g->SetMarkerColor(kGreen);
      TGraph *cu_ib_r5_g = (TGraph*)input5.Get("tr5_g"); 
              cu_ib_r5_g->SetMarkerColor(kOrange);
      TGraph *sn_ib_r5_g = (TGraph*)input6.Get("tr5_g");
              sn_ib_r5_g->SetMarkerColor(kMagenta);
      TGraph *ld_ib_r5_g = (TGraph*)input7.Get("tr5_g");
              ld_ib_r5_g->SetMarkerColor(kBlue);
      TGraph *cxc_ib_r5_g = (TGraph*)input8.Get("tr5_g");
              cxc_ib_r5_g->SetMarkerColor(kGreen);
      
      TGraph *cu_ob_r6_g = (TGraph*)input.Get("tr6_g");
              cu_ob_r6_g->SetMarkerColor(kOrange); 
      TGraph *sn_ob_r6_g = (TGraph*)input2.Get("tr6_g");
              sn_ob_r6_g->SetMarkerColor(kMagenta);
      TGraph *ld_ob_r6_g = (TGraph*)input3.Get("tr6_g");
              ld_ob_r6_g->SetMarkerColor(kBlue);
      TGraph *cxc_ob_r6_g = (TGraph*)input4.Get("tr6_g");
              cxc_ob_r6_g->SetMarkerColor(kGreen);
      TGraph *cu_ib_r6_g = (TGraph*)input5.Get("tr6_g"); 
              cu_ib_r6_g->SetMarkerColor(kOrange);
      TGraph *sn_ib_r6_g = (TGraph*)input6.Get("tr6_g");
              sn_ib_r6_g->SetMarkerColor(kMagenta);
      TGraph *ld_ib_r6_g = (TGraph*)input7.Get("tr6_g");
              ld_ib_r6_g->SetMarkerColor(kBlue);
      TGraph *cxc_ib_r6_g = (TGraph*)input8.Get("tr6_g");
              cxc_ib_r6_g->SetMarkerColor(kGreen);
              
              
      TGraph *mat_e_g = (TGraph*)input9.Get("");
              mat_e_g->SetMarkerStyle(8);
              mat_e_g->SetMarkerColor(kBlack);
              mat_e_g->SetMarkerSize(0.7);
              mat_e_g->SetLineWidth(0.0); 
      TGraph *mat_pip_g = (TGraph*)input10.Get("");
              mat_pip_g->SetMarkerStyle(8);
              mat_pip_g->SetMarkerColor(kBlack);
              mat_pip_g->SetMarkerSize(0.7);
              mat_pip_g->SetLineWidth(0.0);
      TGraph *mat_pim_g = (TGraph*)input11.Get("");
              mat_pim_g->SetMarkerStyle(8);
              mat_pim_g->SetMarkerColor(kBlack);
              mat_pim_g->SetMarkerSize(0.7);
              mat_pim_g->SetLineWidth(0.0);
      TGraph *mat_rho_g = (TGraph*)input12.Get("");
              mat_rho_g->SetMarkerStyle(8);
              mat_rho_g->SetMarkerColor(kBlack);
              mat_rho_g->SetMarkerSize(0.7);
              mat_rho_g->SetLineWidth(0.0);
              
     int mates = mat_e_g->GetN();
     int matpips = mat_pip_g->GetN();
     int matpims = mat_pim_g->GetN();
     int matrhos = mat_rho_g->GetN();
     
     int cuoes = cu_ob_e_g->GetN();
     int cuopips = cu_ob_pip_g->GetN();
     int cuopims = cu_ob_pim_g->GetN();
     int cuorhos = cu_ob_r1_g->GetN();
     int cuies = cu_ib_e_g->GetN();
     int cuipips = cu_ib_pip_g->GetN();
     int cuipims = cu_ib_pim_g->GetN();
     int cuirhos = cu_ib_r1_g->GetN();
     
     int cxcoes = cxc_ob_e_g->GetN();
     int cxcopips = cxc_ob_pip_g->GetN();
     int cxcopims = cxc_ob_pim_g->GetN();
     int cxcorhos = cxc_ob_r1_g->GetN();
     int cxcies = cxc_ib_e_g->GetN();
     int cxcipips = cxc_ib_pip_g->GetN();
     int cxcipims = cxc_ib_pim_g->GetN();
     int cxcirhos = cxc_ib_r1_g->GetN();
     
     int ldoes = ld_ob_e_g->GetN();
     int ldopips = ld_ob_pip_g->GetN();
     int ldopims = ld_ob_pim_g->GetN();
     int ldorhos = ld_ob_r1_g->GetN();
     
     TGraph *diff_e_g = new TGraph();
              diff_e_g->SetMarkerStyle(8);
              diff_e_g->SetMarkerColor(kBlack);
              diff_e_g->SetMarkerSize(0.7);
              diff_e_g->SetLineWidth(0.0);
     TGraph *diff_pip_g = new TGraph();
              diff_pip_g->SetMarkerStyle(8);
              diff_pip_g->SetMarkerColor(kBlack);
              diff_pip_g->SetMarkerSize(0.7);
              diff_pip_g->SetLineWidth(0.0);
     TGraph *diff_pim_g = new TGraph();
              diff_pim_g->SetMarkerStyle(8);
              diff_pim_g->SetMarkerColor(kBlack);
              diff_pim_g->SetMarkerSize(0.7);
              diff_pim_g->SetLineWidth(0.0);
     TGraph *diff_rho_g = new TGraph();
              diff_rho_g->SetMarkerStyle(8);
              diff_rho_g->SetMarkerColor(kBlack);
              diff_rho_g->SetMarkerSize(0.7);
              diff_rho_g->SetLineWidth(0.0);
     
      
      for(int i = 0; i < mates ; i++){
          for(int j = 0; j < cuoes; j++){
              if( mat_e_g->GetPointX(i) == cu_ob_e_g->GetPointX(j) && (mat_e_g->GetPointY(i) - cu_ob_e_g->GetPointY(j)) < 20){
              		diff_e_g->AddPoint(mat_e_g->GetPointX(i), (mat_e_g->GetPointY(i) - cu_ob_e_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < cuies; j++){
              if( mat_e_g->GetPointX(i) == cu_ib_e_g->GetPointX(j) && (mat_e_g->GetPointY(i) - cu_ib_e_g->GetPointY(j)) < 20){
              		diff_e_g->AddPoint(mat_e_g->GetPointX(i), (mat_e_g->GetPointY(i) - cu_ib_e_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < cxcoes; j++){
              if( mat_e_g->GetPointX(i) == cxc_ob_e_g->GetPointX(j) && (mat_e_g->GetPointY(i) - cxc_ob_e_g->GetPointY(j)) < 20){
              		diff_e_g->AddPoint(mat_e_g->GetPointX(i), (mat_e_g->GetPointY(i) - cxc_ob_e_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < cxcies; j++){
              if( mat_e_g->GetPointX(i) == cxc_ib_e_g->GetPointX(j) && (mat_e_g->GetPointY(i) - cxc_ib_e_g->GetPointY(j)) < 20){
              		diff_e_g->AddPoint(mat_e_g->GetPointX(i), (mat_e_g->GetPointY(i) - cxc_ib_e_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < ldoes; j++){
              if( mat_e_g->GetPointX(i) == ld_ob_e_g->GetPointX(j) && (mat_e_g->GetPointY(i) - ld_ob_e_g->GetPointY(j)) < 20){
              		diff_e_g->AddPoint(mat_e_g->GetPointX(i), (mat_e_g->GetPointY(i) - ld_ob_e_g->GetPointY(j)));
              }
          }
      
      
      }
      
      for(int i = 0; i < matpips ; i++){
          for(int j = 0; j < cuopips; j++){
              if( mat_pip_g->GetPointX(i) == cu_ob_pip_g->GetPointX(j) && (mat_pip_g->GetPointY(i) - cu_ob_pip_g->GetPointY(j)) < 3){
              		diff_pip_g->AddPoint(mat_pip_g->GetPointX(i), (mat_pip_g->GetPointY(i) - cu_ob_pip_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < cuipips; j++){
              if( mat_pip_g->GetPointX(i) == cu_ib_pip_g->GetPointX(j) && (mat_pip_g->GetPointY(i) - cu_ib_pip_g->GetPointY(j)) < 3){
              		diff_pip_g->AddPoint(mat_pip_g->GetPointX(i), (mat_pip_g->GetPointY(i) - cu_ib_pip_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < cxcopips; j++){
              if( mat_pip_g->GetPointX(i) == cxc_ob_pip_g->GetPointX(j) && (mat_pip_g->GetPointY(i) - cxc_ob_pip_g->GetPointY(j)) < 3){
              		diff_pip_g->AddPoint(mat_pip_g->GetPointX(i), (mat_pip_g->GetPointY(i) - cxc_ob_pip_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < cxcipips; j++){
              if( mat_pip_g->GetPointX(i) == cxc_ib_pip_g->GetPointX(j) && (mat_pip_g->GetPointY(i) - cxc_ib_pip_g->GetPointY(j)) < 3){
              		diff_pip_g->AddPoint(mat_pip_g->GetPointX(i), (mat_pip_g->GetPointY(i) - cxc_ib_pip_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < ldopips; j++){
              if( mat_pip_g->GetPointX(i) == ld_ob_pip_g->GetPointX(j) && (mat_pip_g->GetPointY(i) - ld_ob_pip_g->GetPointY(j)) < 3){
              		diff_pip_g->AddPoint(mat_pip_g->GetPointX(i), (mat_pip_g->GetPointY(i) - ld_ob_pip_g->GetPointY(j)));
              }
          }
      
      
      }
      
      for(int i = 0; i < matpims ; i++){
          for(int j = 0; j < cuopims; j++){
              if( mat_pim_g->GetPointX(i) == cu_ob_pim_g->GetPointX(j) && (mat_pim_g->GetPointY(i) - cu_ob_pim_g->GetPointY(j)) < 3){
              		diff_pim_g->AddPoint(mat_pim_g->GetPointX(i), (mat_pim_g->GetPointY(i) - cu_ob_pim_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < cuipims; j++){
              if( mat_pim_g->GetPointX(i) == cu_ib_pim_g->GetPointX(j) && (mat_pim_g->GetPointY(i) - cu_ib_pim_g->GetPointY(j)) < 3){
              		diff_pim_g->AddPoint(mat_pim_g->GetPointX(i), (mat_pim_g->GetPointY(i) - cu_ib_pim_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < cxcopims; j++){
              if( mat_pim_g->GetPointX(i) == cxc_ob_pim_g->GetPointX(j) && (mat_pim_g->GetPointY(i) - cxc_ob_pim_g->GetPointY(j)) < 3){
              		diff_pim_g->AddPoint(mat_pim_g->GetPointX(i), (mat_pim_g->GetPointY(i) - cxc_ob_pim_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < cxcipims; j++){
              if( mat_pim_g->GetPointX(i) == cxc_ib_pim_g->GetPointX(j) && (mat_pim_g->GetPointY(i) - cxc_ib_pim_g->GetPointY(j)) < 3){
              		diff_pim_g->AddPoint(mat_pim_g->GetPointX(i), (mat_pim_g->GetPointY(i) - cxc_ib_pim_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < ldopims; j++){
              if( mat_pim_g->GetPointX(i) == ld_ob_pim_g->GetPointX(j) && (mat_pim_g->GetPointY(i) - ld_ob_pim_g->GetPointY(j)) < 3){
              		diff_pim_g->AddPoint(mat_pim_g->GetPointX(i), (mat_pim_g->GetPointY(i) - ld_ob_pim_g->GetPointY(j)));
              }
          }
      
      
      }
      
      
      for(int i = 0; i < matrhos ; i++){
          for(int j = 0; j < cuorhos; j++){
              if( mat_rho_g->GetPointX(i) == cu_ob_r1_g->GetPointX(j) ){
              		diff_rho_g->AddPoint(mat_rho_g->GetPointX(i), (mat_rho_g->GetPointY(i) - cu_ob_r1_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < cuirhos; j++){
              if( mat_rho_g->GetPointX(i) == cu_ib_r1_g->GetPointX(j) ){
              		diff_rho_g->AddPoint(mat_rho_g->GetPointX(i), (mat_rho_g->GetPointY(i) - cu_ib_r1_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < cxcorhos; j++){
              if( mat_rho_g->GetPointX(i) == cxc_ob_r1_g->GetPointX(j) ){
              		diff_rho_g->AddPoint(mat_rho_g->GetPointX(i), (mat_rho_g->GetPointY(i) - cxc_ob_r1_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < cxcirhos; j++){
              if( mat_rho_g->GetPointX(i) == cxc_ib_r1_g->GetPointX(j) ){
              		diff_rho_g->AddPoint(mat_rho_g->GetPointX(i), (mat_rho_g->GetPointY(i) - cxc_ib_r1_g->GetPointY(j)));
              }
          }
          
          for(int j = 0; j < ldorhos; j++){
              if( mat_rho_g->GetPointX(i) == ld_ob_r1_g->GetPointX(j) ){
              		diff_rho_g->AddPoint(mat_rho_g->GetPointX(i), (mat_rho_g->GetPointY(i) - ld_ob_r1_g->GetPointY(j)));
              }
          }
      
      
      }
      
      
     TCanvas *c2 = new TCanvas("c2","",0,0,640,480);
     auto legend = new TLegend(0.1,0.7,0.48,0.9);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   legend->AddEntry(ld_ob_e_g,"LD_{2}","p");
   legend->AddEntry(cxc_ob_e_g,"CxC","p");
   legend->AddEntry(cu_ob_e_g,"CuSn","p");
   legend->Draw();
   	//cu_ob_e_g->Draw();
   	//ld_ob_e_g->Draw("SAME");
   	//cxc_ob_e_g->Draw("SAME");
   	//c1->BuildLedgend();
   	//c2->BuildLegend();
   	c2->Draw();
   	c2->Write();
      
      
     
      auto mge = new TMultiGraph;
   	mge->SetTitle("Electron Yield;Run Number; Y_{e} / Q [nC^{-1}]");
   	mge->Add(cu_ob_e_g);
   	//mge->Add(sn_ob_e_g);
   	mge->Add(ld_ob_e_g);
   	mge->Add(cxc_ob_e_g);
   	mge->Add(cu_ib_e_g);
   	//mge->Add(sn_ib_e_g);
   	//mge->Add(ld_ib_e_g);
   	mge->Add(cxc_ib_e_g);
   	//mge->Add(mat_e_g);
   	//mge->SetMinimum(0.0);
   	mge->SetMaximum(600);
   	//TAxis *axis = mg->GetXaxis();
   	//axis->SetLimits(1,6);
   	mge->Draw("AP");
   	mge->Write();
   	
   	auto mgpp = new TMultiGraph;
   	mgpp->SetTitle("#pi^{+} Yield;Run Number; Y_{#pi^{+}} / Q [nC^{-1}]");
   	mgpp->Add(cu_ob_pip_g);
   	//mgpp->Add(sn_ob_pip_g);
   	mgpp->Add(ld_ob_pip_g);
   	mgpp->Add(cxc_ob_pip_g);
   	mgpp->Add(cu_ib_pip_g);
   	//mgpp->Add(sn_ib_pip_g);
   	//mgpp->Add(ld_ib_pip_g);
   	mgpp->Add(cxc_ib_pip_g);
   	//mgpp->Add(mat_pip_g);
   	mgpp->SetMaximum(60);
   	mgpp->Draw("AP");
   	mgpp->Write();
   	
   	auto mgpm = new TMultiGraph;
   	mgpm->SetTitle("#pi^{-} Yield;Run Number; Y_{#pi^{-}} / Q [nC^{-1}]");
   	mgpm->Add(cu_ob_pim_g);
   	//mgpm->Add(sn_ob_pim_g);
   	mgpm->Add(ld_ob_pim_g);
   	mgpm->Add(cxc_ob_pim_g);
   	mgpm->Add(cu_ib_pim_g);
   	//mgpm->Add(sn_ib_pim_g);
   	//mgpm->Add(ld_ib_pim_g);
   	mgpm->Add(cxc_ib_pim_g);
   	//mgpm->Add(mat_pim_g);
   	mgpm->SetMaximum(60);
   	mgpm->Draw("AP");
   	mgpm->Write();
   	
   	auto mgr1 = new TMultiGraph;
   	mgr1->SetTitle("#rho^{0} Yield;Run Number; Y_{#rho^{0}} / Q [nC^{-1}]");
   	mgr1->Add(cu_ob_r1_g);
   	//mgr1->Add(sn_ob_r1_g);
   	mgr1->Add(ld_ob_r1_g);
   	mgr1->Add(cxc_ob_r1_g);
   	mgr1->Add(cu_ib_r1_g);
   	//mgr1->Add(sn_ib_r1_g);
   	//mgr1->Add(ld_ib_r1_g);
   	mgr1->Add(cxc_ib_r1_g);
   	//mgr1->Add(mat_rho_g);
   	mgr1->SetMaximum(0.4);
   	mgr1->Draw("AP");
   	mgr1->Write();
   	//std::cout << cu_ob_r1_g->GetPointX(40) << std::endl;
   	auto mgr2 = new TMultiGraph;
   	mgr2->SetTitle("#rho^{0} Yield for Q^{2} bin 2;Run Number; Y_{#rho^{0}}");
   	mgr2->Add(cu_ob_r2_g);
   	mgr2->Add(sn_ob_r2_g);
   	mgr2->Add(ld_ob_r2_g);
   	mgr2->Add(cxc_ob_r2_g);
   	mgr2->Add(cu_ib_r2_g);
   	mgr2->Add(sn_ib_r2_g);
   	mgr2->Add(ld_ib_r2_g);
   	mgr2->Add(cxc_ib_r2_g);
   	mgr2->Draw("AP");
   	mgr2->Write();
   	
   	auto mgr3 = new TMultiGraph;
   	mgr3->SetTitle("#rho^{0} Yield for Q^{2} bin 3;Run Number; Y_{#rho^{0}}");
   	mgr3->Add(cu_ob_r3_g);
   	mgr3->Add(sn_ob_r3_g);
   	mgr3->Add(ld_ob_r3_g);
   	mgr3->Add(cxc_ob_r3_g);
   	mgr3->Add(cu_ib_r3_g);
   	mgr3->Add(sn_ib_r3_g);
   	mgr3->Add(ld_ib_r3_g);
   	mgr3->Add(cxc_ib_r3_g);
   	mgr3->Draw("AP");
   	mgr3->Write();
   	
   	auto mgr4 = new TMultiGraph;
   	mgr4->SetTitle("#rho^{0} Yield for Q^{2} bin 4;Run Number; Y_{#rho^{0}}");
   	mgr4->Add(cu_ob_r4_g);
   	mgr4->Add(sn_ob_r4_g);
   	mgr4->Add(ld_ob_r4_g);
   	mgr4->Add(cxc_ob_r4_g);
   	mgr4->Add(cu_ib_r4_g);
   	mgr4->Add(sn_ib_r4_g);
   	mgr4->Add(ld_ib_r4_g);
   	mgr4->Add(cxc_ib_r4_g);
   	mgr4->Draw("AP");
   	mgr4->Write();
   	
   	auto mgr5 = new TMultiGraph;
   	mgr5->SetTitle("#rho^{0} Yield for Q^{2} bin 5;Run Number; Y_{#rho^{0}}");
   	mgr5->Add(cu_ob_r5_g);
   	mgr5->Add(sn_ob_r5_g);
   	mgr5->Add(ld_ob_r5_g);
   	mgr5->Add(cxc_ob_r5_g);
   	mgr5->Add(cu_ib_r5_g);
   	mgr5->Add(sn_ib_r5_g);
   	mgr5->Add(ld_ib_r5_g);
   	mgr5->Add(cxc_ib_r5_g);
   	mgr5->Draw("AP");
   	mgr5->Write();
   	
   	auto mgr6 = new TMultiGraph;
   	mgr6->SetTitle("#rho^{0} Yield for Q^{2} bin 6;Run Number; Y_{#rho^{0}}");
   	mgr6->Add(cu_ob_r6_g);
   	//mgr6->Add(sn_ob_r6_g);
   	mgr6->Add(ld_ob_r6_g);
   	mgr6->Add(cxc_ob_r6_g);
   	mgr6->Add(cu_ib_r6_g);
   	//mgr6->Add(sn_ib_r6_g);
   	mgr6->Add(ld_ib_r6_g);
   	mgr6->Add(cxc_ib_r6_g);
   	mgr6->Draw("AP");
   	mgr6->Write();
   	
   	diff_e_g->SetTitle("Electron Yield Difference between Mathieu and Matt; Run Number; Y_{Mathieu} - Y_{Matt} [nC^{-1}]");
   	diff_e_g->SetMaximum(100);
   	diff_e_g->SetMinimum(-100);
   	diff_e_g->Draw();
   	diff_e_g->Write();
   	
   	diff_pip_g->SetTitle("#pi^{+} Yield Difference between Mathieu and Matt; Run Number; Y_{Mathieu} - Y_{Matt} [nC^{-1}]");
   	diff_pip_g->SetMaximum(100);
   	diff_pip_g->SetMinimum(-100);
   	diff_pip_g->Draw();
   	diff_pip_g->Write();
   	
   	diff_pim_g->SetTitle("#pi^{-} Yield Difference between Mathieu and Matt; Run Number; Y_{Mathieu} - Y_{Matt} [nC^{-1}]");
   	diff_pim_g->SetMaximum(100);
   	diff_pim_g->SetMinimum(-100);
   	diff_pim_g->Draw();
   	diff_pim_g->Write();
   	
   	diff_rho_g->SetTitle("#rho^{0} Yield Difference between Mathieu and Matt; Run Number; Y_{Mathieu} - Y_{Matt} [nC^{-1}]");
   	diff_rho_g->SetMaximum(3);
   	diff_rho_g->SetMinimum(-3);
   	diff_rho_g->Draw();
   	diff_rho_g->Write();
     
}


