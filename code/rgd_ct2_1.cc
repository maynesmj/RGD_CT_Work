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
#include <TH1D.h>
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
   gstyle->SetPalette(kRainBow); 

   auto *c2 = new TCanvas("c2", "c2");
   //TStyle *gstyle = new TStyle();
   //gstyle->SetPalette(kRainBow); 

   auto *c3 = new TCanvas("c3", "c3");
   
   TH1D *z_hist = new TH1D("Zh", " Zh ", 200, 0, 10);
   //bC_hist->SetStats(0);

   TH2F *qz_hist = new TH2F("Q2 vs lc", "Q2 vs lc", 200, 0, 10, 200, 0, 10);
   
   TH2F *xy_hist = new TH2F("px vs py", "px vs py", 200, -3, 3, 200, -3, 3);

   TH1D *q_hist = new TH1D("Q2", " Q2", 200, 0, 10);
   
   TH1D *lc_hist = new TH1D("lc", " lc", 200, 0, 10);
   
   TH1D *e_hist = new TH1D("electrons", " electrons ", 200, -300, 300);
   
   TH1D *pip_hist = new TH1D("pip", "pip", 200, -300, 300);
   
   TH1D *pim_hist = new TH1D("pim", "pim", 200, -300, 300);

   //TH2F *qlc_hist = new TH2F("Q2 vs lc", "Q2 vs lc", 200, 0, 10, 200, 0, 20);

   TH1D *w_hist = new TH1D("w (Invariant Mass)", " ", 200, 0, 5);

   //TH2F *tz_hist = new TH2F("-t vs rho_pz", "-t vs rho_pz", 200, 0, 10, 200, 0, 20);

   TH1D *t_hist = new TH1D("-t", "-t", 200, 0, 20);

   TH1D *wonphe_hist = new TH1D("wonphe", " ", 500, 0, 0.5);

   TH1D *wnphe_hist = new TH1D("wnphe", " ", 500, 0, 0.5);

   TH1D *wopcal_hist = new TH1D("Total Energy", " ", 500, 0, 1);

   TH1D *wpcal_hist = new TH1D("PCAL", " ", 500, 0, 1);

   TH2F *wocut_hist = new TH2F("wocut", "wocut", 500, 0, 20, 500, 0, 1);

   TH2F *wcut_hist = new TH2F("Energy", "", 500, 0.6, 11, 500, 0, 0.3);

   auto *vxs = new THStack("vxs", "vxs");
   auto *vys = new THStack("vys", "vys");
   auto *etots = new THStack("etots", "etots");
   auto *ps = new THStack("ps" , "ps");

   


   TH1F *ivx_hist = new TH1F("sector 1 vx", "Sector 1 vx", 300, -5, 5);
       ivx_hist->SetFillColor(kCyan+1);
   TH1D *ivy_hist = new TH1D("sector 1 vy", "Sector 1 vy", 300, -5, 5);
       ivy_hist->SetFillColor(kCyan+1);
   TH2F *ixy_hist = new TH2F("sector 1 px .vs py", "sector 1 px .vs py", 200, -1, 1, 200, -1, 1);
   TH1D *ivz_hist = new TH1D("sector 1 vz", "Sector 1 vz", 300, -20, 20);

   TH1F *iivx_hist = new TH1F("sector 2 vx", "Sector 2 vx", 300, -5, 5);
	iivx_hist->SetFillColor(kOrange+1);
   TH1D *iivy_hist = new TH1D("sector 2 vy", "Sector 2 vy", 300, -5, 5);
	iivy_hist->SetFillColor(kOrange+1);
   TH2F *iixy_hist = new TH2F("sector 2 vx .vs vy", "" , 500, 0.6, 11, 500, 0, 0.3);
   TH1D *iivz_hist = new TH1D("sector 2 vz", "Sector 2 vz", 300, -20, 20);

   TH1F *iiivx_hist = new TH1F("sector 3 vx", "Sector 3 vx", 300, -5, 5);
	iiivx_hist->SetFillColor(kYellow+1);
   TH1D *iiivy_hist = new TH1D("sector 3 vy", "Sector 3 vy", 300, -5, 5);
	iiivy_hist->SetFillColor(kYellow+1);
   TH2F *iiixy_hist = new TH2F("sector 3 vx .vs vy", "sector 3 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1D *iiivz_hist = new TH1D("sector 3 vz", "Sector 3 vz", 300, -20, 20);

   TH1F *ivvx_hist = new TH1F("sector 4 vx", "Sector 4 vx", 300, -5, 5);
	ivvx_hist->SetFillColor(kGreen+1);
   TH1D *ivvy_hist = new TH1D("sector 4 vy", "Sector 4 vy", 300, -5, 5);
	ivvy_hist->SetFillColor(kGreen+1);
   TH2F *ivxy_hist = new TH2F("sector 4 vx .vs vy", "sector 4 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1D *ivvz_hist = new TH1D("sector 4 vz", "Sector 4 vz", 300, -20, 20);

   TH1F *vvx_hist = new TH1F("sector 5 vx", "Sector 5 vx", 300, -5, 5);
	vvx_hist->SetFillColor(kBlue);
   TH1D *vvy_hist = new TH1D("sector 5 vy", "Sector 5 vy", 300, -5, 5);
	vvy_hist->SetFillColor(kBlue);
   TH2F *vxy_hist = new TH2F("sector 5 vx .vs vy", "sector 5 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1D *vvz_hist = new TH1D("sector 5 vz", "Sector 5 vz", 300, -20, 20);

   TH1F *vivx_hist = new TH1F("sector 6 vx", "Sector 6 vx", 300, -5, 5);
	vivx_hist->SetFillColor(kMagenta+1);
   TH1D *vivy_hist = new TH1D("sector 6 vy", "Sector 6 vy", 300, -5, 5);
        vivy_hist->SetFillColor(kMagenta+1);
   TH2F *vixy_hist = new TH2F("sector 6 vx .vs vy", "sector 6 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1D *vivz_hist = new TH1D("sector 6 vz", "Sector 6 vz", 300, -20, 20);


   

   TH1D *rm_no = new TH1D("Without Cuts", "Without Cuts", 200, 0.2, 1.4);
      //rm_no->SetStats(0);
   TH1D *rm_w = new TH1D("With W Cut", "With W Cut ", 200, 0.2, 1.4);
      //rm_w->SetStats(0);
   TH1D *rm_wt = new TH1D("With W Cut & -t Cut", "With W Cut & -t Cut ", 200, 0.2, 1.4);
      //rm_wt->SetStats(0);
   TH1D *rm_wtz = new TH1D("With W Cut & -t Cut & z Cut", " With W Cut & -t Cut & z Cut", 200, 0.2, 1.4);
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
      vetot_p->SetMarkerColor(kBlue);
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
      hel_fcup->SetMarkerColor(kBlue);
      hel_fcup->SetMarkerStyle(kFullCircle);
      hel_fcup->SetLineWidth(0);

   auto hel_fcupgated = new TGraph();
      hel_fcupgated->SetMarkerColor(kBlack);
      hel_fcupgated->SetMarkerStyle(kFullCircle);
      hel_fcupgated->SetLineWidth(0);

   auto hel_clock = new TGraph();
      hel_clock->SetMarkerColor(kBlue);
      hel_clock->SetMarkerStyle(kFullCircle);
      hel_clock->SetLineWidth(0);

   auto hel_clockgated = new TGraph();
      hel_clockgated->SetMarkerColor(kBlack);
      hel_clockgated->SetMarkerStyle(kFullCircle);
      hel_clockgated->SetLineWidth(0);

   TH1D *ietot_hist = new TH1D("sector 1 E_{tot}", "Sector 1 E_{tot}", 200, 0, 0.4);
       ietot_hist->SetFillColor(kCyan+1);
   TH1D *ip_hist = new TH1D("sector 1 P", "Q2 (2->2.5)", 200, 0, 2);
       ip_hist->SetFillColor(kBlue);
   TH1D *p_hist = new TH1D("sector 0 P", "Q2 (1->2)", 200, 0, 2);
       p_hist->SetFillColor(kBlue);

   TH1D *iietot_hist = new TH1D("sector 2 E_{tot}", "Sector 2 E_{tot}", 200, 0, 0.4);
	iietot_hist->SetFillColor(kOrange+1);
   TH1D *iip_hist = new TH1D("sector 2 P", "Q2 (2.5->3)", 200, 0, 2);
	iip_hist->SetFillColor(kBlue);
   

   TH1D *iiietot_hist = new TH1D("sector 3 E_{tot}", "Sector 3 E_{tot}", 200, 0, 0.4);
	iiietot_hist->SetFillColor(kYellow+1);
   TH1D *iiip_hist = new TH1D("sector 3 P", "Q2 (3->3.5)", 200, 0, 2);
	iiip_hist->SetFillColor(kBlue);
   

   TH1D *ivetot_hist = new TH1D("sector 4 E_{tot}", "Sector 4 E_{tot}", 200, 0, 0.4);
	ivetot_hist->SetFillColor(kGreen+1);
   TH1D *ivp_hist = new TH1D("sector 4 P", "Q2 (3.5->6)", 200, 0, 2);
	ivp_hist->SetFillColor(kBlue);
   

   TH1D *vetot_hist = new TH1D("sector 5 E_{tot}", "Sector 5 E_{tot}", 200, 0, 0.4);
	vetot_hist->SetFillColor(kBlue);
   TH1D *vp_hist = new TH1D("sector 5 P", "Q2 (4->4.5)", 200, 0, 2);
	vp_hist->SetFillColor(kBlue);
   

   TH1D *vietot_hist = new TH1D("sector 6 E_{tot}", "Sector 6 E_{tot}", 200, 0, 0.4);
	vietot_hist->SetFillColor(kMagenta+1);
   TH1D *vip_hist = new TH1D("sector 6 P", "Q2 (4.5->6)", 200, 0, 2);
        vip_hist->SetFillColor(kBlue);



   //TH2F *z_hist = new TH2F("px vs py", "px vs py", 200, -2, 2, 200, -2, 2);

   //auto cgr = new TGraph();
  // std::cout << "1 th" << std::endl;

   std::vector<std::string> v3_a;
   int num = 0;


   TFile f("ld_ob_test_kin4.root", "update");

    std::ifstream v3("ld_ob.list");
         std::string line;
   while(std::getline(v3, line)){
	std::stringstream ss(line);
        std::string hip;
        std::getline(ss, hip, ',');
        v3_a.push_back(hip);
        num++;
}


   //std::ofstream fw("actual_2_18318_inbend_yield.txt", std::ofstream::out);   //011270 fall & 011361 spring

   std::string dir_path[4] = { "/volatile/clas12/rg-d/production/lumi/v0_LD2_aicv/dcalign/recon/", "/volatile/clas12/rg-d/production/prod/v0_ob_CxC/dst/recon/", "/volatile/clas12/rg-d/production/pass0.2/mon/recon/", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/"};
std::cout << "1 th" << std::endl;
   
   std::string run_num[4] = {"018321", "018451", "018490", "019128"}; //"011283"};
   int file_num[4] = {56, 200, 10, 260};//need to add the files for 18309       361 for 18451   200 for v0

   int first_num = 0;
   int second_num = 4;

   int tt = 0;


   for(int i = 3; i < 4; i++){
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
      for(int p = 0; p < v3_a.size(); p++){
         int new_first_num = 0;
         int new_second_num = 0;
         if (p <= 81 && i == 1){
            new_first_num = first_num + (5*p);
            new_second_num = second_num + (5*p);
         } else {
            new_first_num = first_num + (5*p) + 10;
            new_second_num = second_num + (5*p) + 10;
         }

         std::string str_first_num;
         std::string str_second_num;
        // std::cout << "first num = " << new_first_num << " | second num = " << new_second_num << " " << i << std::endl;

       

         if(new_first_num < 10){
            str_first_num = std::to_string(0) + std::to_string(0) + std::to_string(0)  + std::to_string(0) + std::to_string(new_first_num);
         } else if(new_first_num >= 10 && new_first_num < 100){
            str_first_num = std::to_string(0) + std::to_string(0) + std::to_string(0) + std::to_string(new_first_num);
         }else if(new_first_num >= 100 && new_first_num < 1000){
               str_first_num = std::to_string(0) + std::to_string(0) + std::to_string(new_first_num);
         }else if(new_first_num >= 1000){
            str_first_num = std::to_string(0) + std::to_string(new_first_num);
         } 

         if(new_second_num < 10){
            str_second_num = std::to_string(0) + std::to_string(0) + std::to_string(0) + std::to_string(0) + std::to_string(new_second_num);
         } else if(new_second_num >= 10 && new_second_num < 100){
               if(i == 0 && new_first_num == 40){      // 90 & 90 for 18319
                   str_second_num = "00041";
               } else {
                   str_second_num = std::to_string(0) + std::to_string(0) + std::to_string(0) + std::to_string(new_second_num);
               }
         }else if(new_second_num >= 100 && new_second_num < 1000){
            //if(i == 1 && new_first_num == 280){      // 90 & 90 for 18319
            //      str_second_num = "00283";
            //   } else {
            str_second_num =  std::to_string(0) + std::to_string(0) + std::to_string(new_second_num);
            //   }  
         }else if(new_second_num >= 1000){
            if(i == 1 && new_first_num == 1815){      // 90 & 90 for 18319
                  str_second_num = "01816";
               } else {
                  str_second_num = std::to_string(0) + std::to_string(new_second_num);
               }  
         } 
         std::string inFile;
         if (i < 3){
            inFile = dir_path[i] + run_num[i] + "/rec_clas_" + run_num[i] + ".evio." + str_first_num + "-" + str_second_num + ".hipo";
         } else {
            //if( p != 3 && p!= 7){
/*
	        std::string direr_path[8] =    {"/volatile/clas12/rg-d/production/pass0.3/mon/recon/019128/rec_clas_019128.evio.00040.hipo","/volatile/clas12/rg-d/production/pass0.3/mon/recon/019128/rec_clas_019128.evio.00041.hipo", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/019128/rec_clas_019128.evio.00042.hipo", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/019128/rec_clas_019128.evio.00044.hipo", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/019128/rec_clas_019128.evio.00045.hipo", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/019128/rec_clas_019128.evio.00046.hipo", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/019128/rec_clas_019128.evio.00048.hipo", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/019128/rec_clas_019128.evio.00049.hipo"};
*/
            //int bnb = p;
            //std::string bnbs = std::to_string(bnb);


	     std::string direr_path[10] =    {"/volatile/clas12/rg-d/production/pass0.3/mon/recon/019125/rec_clas_019125.evio.00040.hipo","/volatile/clas12/rg-d/production/pass0.3/mon/recon/019125/rec_clas_019125.evio.00041.hipo", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/019125/rec_clas_019125.evio.00042.hipo", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/019125/rec_clas_019125.evio.00043.hipo","/volatile/clas12/rg-d/production/pass0.3/mon/recon/019125/rec_clas_019125.evio.00044.hipo", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/019125/rec_clas_019125.evio.00045.hipo", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/019125/rec_clas_019125.evio.00046.hipo", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/019125/rec_clas_019125.evio.00047.hipo", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/019125/rec_clas_019125.evio.00048.hipo", "/volatile/clas12/rg-d/production/pass0.3/mon/recon/019125/rec_clas_019125.evio.00049.hipo"};


//  /volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/
        std::string direrer_path[20] = {"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00050-00054.hipo",   
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00055-00059.hipo",  
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00090-00094.hipo",   
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00095-00099.hipo",   
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00270-00274.hipo",   
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00275-00279.hipo",  
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00280-00284.hipo",   
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00285-00289.hipo",  
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00310-00314.hipo",  
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00315-00319.hipo",   
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00320-00324.hipo",   
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00325-00329.hipo",  
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00330-00334.hipo",   
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00335-00339.hipo",   
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00340-00344.hipo",   
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00345-00349.hipo",  
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00350-00354.hipo",  
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00355-00359.hipo",  
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00360-00364.hipo",  
					"/volatile/clas12/rg-d/production/prod/v3_ob_CxC/dst/recon/018463/rec_clas_018463.evio.00365-00369.hipo"};



          //std::string dirererer_path[] = {};
        


            inFile =  v3_a[p];//std::getline(v3, line);     //  dir_path[i] + run_num[i] + "/rec_clas_" + run_num[i] + ".evio.0004" + bnbs + ".hipo";
            //}
         }
       
int pid = 0;
   float px = 0.0;
   float py = 0.0;
   float pz = 0.0;
   float vx = 0.0;
   float vy = 0.0;
   float vz = 0.0;
   float energy = 0.0;

   int pip_pid = 0;
   float pip_px = 0.0;
   float pip_py = 0.0;                          //These values will be fill the root tree
   float pip_pz = 0.0;
   float pip_vx = 0.0;
   float pip_vy = 0.0;                          
   float pip_vz = 0.0;
   float pip_energy = 0.0;
   
   // in GeV/c^2
   float m_pro = 0.938; //in Gev/c^2
   //float m_ele = 0.000511; // in GeV/c^2
   float m_pi = 0.1396;
   float m_rho = 0.770;

   int pim_pid = 0;
   float pim_px = 0.0;
   float pim_py = 0.0;
   float pim_pz = 0.0;
   float pim_vx = 0.0;
   float pim_vy = 0.0;
   float pim_vz = 0.0;
   float pim_energy = 0.0;

   float rho_px = 0.0;
   float rho_py = 0.0;
   float rho_pz = 0.0;
   float rho_energy = 0.0;
   float Q2 = 0.0;
   float w = 0.0;
   float t = 0.0;
   float lc = 0.0;
   float Zh = 0.0;
   float rho_mass = 0.0;
   float beam_energy = 0.0;

	int r = 0;
	int i = 0;
	int c = 0;

   if(tt < 1){
      beam_energy = 10.5;
   } else {
      beam_energy = 10.5;
   }

	
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

      hipo::bank cher(factory.getSchema("REC::Cherenkov"));
      hipo::bank cal(factory.getSchema("REC::Calorimeter"));

      hipo::event      event; 								//This declares a variable of class hipo named event

      //int counter_p = 
      int counter = 0;
      

      TLorentzVector e_lv;
      e_lv.SetPxPyPzE(0.0, 0.0, beam_energy, beam_energy);

      TLorentzVector p_lv;
      p_lv.SetPxPyPzE(0.0, 0.0, 0.0, m_pro);

      TLorentzVector te_lv;
      te_lv.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

	   int e_num = 0;
	   int e_pi_num = 0;
	   int e_neg_num = 0;
	   int e_neg_pi_num = 0; 

	   int pid_e_count = 0;
      int pid_pip_count = 0;
      int pid_pim_count = 0;

      float bigger = 0.0;
      float diff_bC = 0.0;
      float total_bC = 0.0;


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
         
       // particle.show();
        // events.show();
        //hel.show();
        //cher.show();
        //cal.show();

         int nrows = 0; 
         int a_pid[200];
         float a_px[200]; 
         float a_py[200];
         float a_pz[200];
         float a_vx[200];
         float a_vy[200];
         float a_vz[200];
         float a_energy[200];
      

         float a_rho_mass[200];
         float a_rho_px[200];
         float a_rho_py[200];
         float a_rho_pz[200];
         float a_rho_energy[200];

         int a_pip_pid[200];
         float a_pip_px[200];
         float a_pip_py[200];
         float a_pip_pz[200];
         float a_pip_vx[200];
         float a_pip_vy[200];
         float a_pip_vz[200];
         float a_pip_energy[200];

         int a_pim_pid[200];                                                     //This section declares the arrays that saves the data of one event. Each will be rewritten for each event.
         float a_pim_px[200];
         float a_pim_py[200];
         float a_pim_pz[200];
         float a_pim_vx[200];
         float a_pim_vy[200];
         float a_pim_vz[200];
         float a_pim_energy[200];

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
        
       
         //std::cout << nrows << std::endl;
         for(int row = 0; row < nrows; row++){
            
            if(particle.getInt("pid",row) == 11){
               
		         float temp_vz = 0.0;
    		      int temp_s = 0;
        

                     temp_vz = particle.getFloat("vz", row);
		     temp_s = particle.getInt("status", row);

                     if(temp_vz >= -20 && temp_vz <= 20){
                        //std::cout << nrows << std::endl;
			               
                        if(abs(temp_s)/2000 == 1){
                          
			if(temp_s < 0){	   
				a_row.push_back(row);
				
           	      		   a_pid[pid_e_count] = 11;
           	      		   a_px[pid_e_count] = particle.getFloat("px",row);
            	   		   a_py[pid_e_count] = particle.getFloat("py",row);
              	   		   a_pz[pid_e_count] = particle.getFloat("pz",row);
                  		   a_vx[pid_e_count] = particle.getFloat("vx",row);
           	      		   a_vy[pid_e_count] = particle.getFloat("vy",row);
           	      		   a_vz[pid_e_count] = temp_vz;
                  		   a_energy[pid_e_count] = sqrt(pow(a_px[pid_e_count],2) + pow(a_py[pid_e_count],2) + pow(a_pz[pid_e_count],2));
                       		   pid_e_count++;
					te_lv.SetPxPyPzE(a_px[pid_e_count],a_py[pid_e_count],a_pz[pid_e_count],a_energy[pid_e_count]);
					ixy_hist->Fill(a_px[row], a_py[row]);
                                   float sec_dec =  te_lv.Phi() * 180/M_PI;
                                   std::cout << sec_dec << std::endl;
                                    //el_lv.SetPxPyPzE(px,py,pz,energy);
                                    //phi = el_lv.Phi() * 180/M_PI;
                                   //  ivx_hist->Fill(a_vx[row]);
                          
				if(sec_dec >= -30 && sec_dec < 30) {
				ivx_hist->Fill(a_vx[pid_e_count]);
   				ivy_hist->Fill(a_vy[pid_e_count]);
				
                                if( (ivx_sigma * -3) <= a_vx[pid_e_count] && (ivx_sigma * 3) >= a_vx[pid_e_count] && (ivy_sigma * -3) <= a_vy[pid_e_count] && (ivy_sigma * 3) >= a_vy[pid_e_count]){
				ivz_hist->Fill(a_vz[pid_e_count]);
				}
                        } else if(sec_dec >= 30 && sec_dec < 90) {
				iivx_hist->Fill(a_vx[row]);
   				iivy_hist->Fill(a_vy[row]);
				iixy_hist->Fill(a_px[row], a_py[row]);
				if( (iivx_sigma * -3) <= a_vx[row] && (iivx_sigma * 3) >= a_vx[row] && (iivy_sigma * -3) <= a_vy[row] && (iivy_sigma * 3) >= a_vy[row]){
				iivz_hist->Fill(a_vz[row]);
				}
                        } else if(sec_dec >= 90 && sec_dec < 150) {
				iiivx_hist->Fill(a_vx[row]);
   				iiivy_hist->Fill(a_vy[row]);
				iiixy_hist->Fill(a_px[row], a_py[row]);
				if( (iiivx_sigma * -3) <= a_vx[row] && (iiivx_sigma * 3) >= a_vx[row] && (iiivy_sigma * -3) <= a_vy[row] && (iiivy_sigma * 3) >= a_vy[row]){
				iiivz_hist->Fill(a_vz[row]);
				}
                        } else if(sec_dec >= 150 && sec_dec < -150) {
				ivvx_hist->Fill(a_vx[row]);
   				ivvy_hist->Fill(a_vy[row]);
				ivxy_hist->Fill(a_px[row], a_py[row]);
				if( (ivvx_sigma * -3) <= a_vx[row] && (ivvx_sigma * 3) >= a_vx[row] && (ivvy_sigma * -3) <= a_vy[row] && (ivvy_sigma * 3) >= a_vy[row]){
				ivvz_hist->Fill(a_vz[row]);
				}
                        } else if(sec_dec >= -150 && sec_dec < -90) {
				vvx_hist->Fill(a_vx[row]);
   				vvy_hist->Fill(a_vy[row]);
				vxy_hist->Fill(a_px[row], a_py[row]);
				if( (vvx_sigma * -3) <= a_vx[row] && (vvx_sigma * 3) >= a_vx[row] && (vvy_sigma * -3) <= a_vy[row] && (vvy_sigma * 3) >= a_vy[row]){
				vvz_hist->Fill(a_vz[row]);
				}
                        } else if(sec_dec >= -90 && sec_dec < -30) {
				vivx_hist->Fill(a_vx[row]);
   				vivy_hist->Fill(a_vy[row]);
				vixy_hist->Fill(a_px[row], a_py[row]);
				if( (vivx_sigma * -3) <= a_vx[row] && (vivx_sigma * 3) >= a_vx[row] && (vivy_sigma * -3) <= a_vy[row] && (vivy_sigma * 3) >= a_vy[row]){
				vivz_hist->Fill(a_vz[row]);
				}
                         }
                       		
          //            } 
          //         }  
          }
	               }
               }
            } else if(particle.getInt("pid",row) == 211){
              	a_pip_pid[pid_pip_count] = 211;
              	a_pip_px[pid_pip_count] = particle.getFloat("px",row);
           	   a_pip_py[pid_pip_count] = particle.getFloat("py",row);
              	a_pip_pz[pid_pip_count] = particle.getFloat("pz",row);
               a_pip_vx[pid_pip_count] = particle.getFloat("vx",row);
              	a_pip_vy[pid_pip_count] = particle.getFloat("vy",row);
           	   a_pip_vz[pid_pip_count] = particle.getFloat("vz",row);
               a_pip_energy[pid_pip_count] = sqrt(pow(a_pip_px[pid_pip_count],2) + pow(a_pip_py[pid_pip_count],2) + pow(a_pip_pz[pid_pip_count],2) + pow(m_pi,2));
               pid_pip_count++;
           
            } else if(particle.getInt("pid",row) == -211){
              	a_pim_pid[pid_pim_count] = -211;
              	a_pim_px[pid_pim_count] = particle.getFloat("px",row);
              	a_pim_py[pid_pim_count] = particle.getFloat("py",row);
           	   a_pim_pz[pid_pim_count] = particle.getFloat("pz",row);
               a_pim_vx[pid_pim_count] = particle.getFloat("vx",row);
              	a_pim_vy[pid_pim_count] = particle.getFloat("vy",row);
              	a_pim_vz[pid_pim_count] = particle.getFloat("vz",row);
               a_pim_energy[pid_pim_count] = sqrt(pow(a_pim_px[pid_pim_count],2) + pow(a_pim_py[pid_pim_count],2) + pow(a_pim_pz[pid_pim_count],2) + pow(m_pi,2));
               pid_pim_count++;
            }
         }

	int psize = a_row.size();
//ietot_p->AddPoint(x,y);
        int sect = 0;
/*
         for (int row = 0; row < psize; row++){
           for(int pindex = 0; pindex < 20; pindex++){
		int pin = cal.getInt("pindex",pindex);
                if(pin == a_row[row]){
        		float pcal = cal.getFloat("energy", pindex);
         	        e_total = e_total + pcal;

                        sect = cal.getInt("sector",pindex);
				
			if( cal.getInt("layer", pindex) == 1){
				ppcal = pcal;
			}
 

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

                        //}
 


                }
           }
           	//if(ppcal != 0){
		//wopcal_hist->Fill(e_total/a_energy[row]);
        	//	 wcut_hist->Fill(a_energy[row], e_total/a_energy[row]);
		//	 wnphe_hist->Fill(nphe);
       		//           wpcal_hist->Fill(ppcal);
                //       if((e_total/a_energy[row]) <= (0.0137 * 3) + 0.24  && (e_total/a_energy[row]) >= (-0.0137 * 3) + 0.24){
		//		iixy_hist->Fill(a_energy[row], e_total/a_energy[row]);
		//	}

                  for(int pindex = 0; pindex < 20; pindex++){
		
                        sect = cal.getInt("sector",pindex);

		    if(e_total/a_energy[row] <= 1){

		        if(sect == 1) {
				ietot_p->AddPoint(a_energy[row], e_total/a_energy[row]);
				ietot_hist->Fill(e_total/a_energy[row]);
				//ip_hist->Fill(a_energy[row]);
                        } else if(sect == 2) {
				iietot_p->AddPoint(a_energy[row], e_total/a_energy[row]);
				iietot_hist->Fill(e_total/a_energy[row]);
				//iip_hist->Fill(a_energy[row]);
                        } else if(sect == 3) {
				iiietot_p->AddPoint(a_energy[row], e_total/a_energy[row]);
				iiietot_hist->Fill(e_total/a_energy[row]);
				//iiip_hist->Fill(a_energy[row]);
                        } else if(sect == 4) {
				ivetot_p->AddPoint(a_energy[row], e_total/a_energy[row]);
				ivetot_hist->Fill(e_total/a_energy[row]);
				//ivp_hist->Fill(a_energy[row]);
                        } else if(sect == 5) {
				vetot_p->AddPoint(a_energy[row], e_total/a_energy[row]);
				vetot_hist->Fill(e_total/a_energy[row]);
				//vp_hist->Fill(a_energy[row]);
                        } else if(sect == 6) {
				vietot_p->AddPoint(a_energy[row], e_total/a_energy[row]);
				vietot_hist->Fill(e_total/a_energy[row]);
				//vip_hist->Fill(a_energy[row]);
                        }
                     }
                 }
                
         }
  */       
	      int w_cut = 0;
         int z_cut = 0;
         int t_cut = 0;
       
         if( pid_e_count > 0  && pid_pip_count > 0  && pid_pim_count > 0){               //This ensures that only the events with all three particle will be looked at

            int ne = pid_e_count;
            int npip = pid_pip_count;
            int npim = pid_pim_count;                                                 //Since the arrays are filled from the first position until all of that particle is in
         
         
            for(int i = 0; i < npip; i++){
               pip_pid = a_pip_pid[i];
               pip_px = a_pip_px[i];
               pip_py = a_pip_py[i];
               pip_pz = a_pip_pz[i];
               pip_vx = a_pip_vx[i];
               pip_vy = a_pip_vy[i];
               pip_vz = a_pip_vz[i];
               pip_energy = a_pip_energy[i];
	       pip_hist->Fill(pip_pid);    
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
            }
         
            for(int i = 0; i < ne; i++){
               if(e_final < a_energy[i]){                                 //This saves the highest energy electron
                  e_final = a_energy[i];
                  address = i;
               }
            } 
         
            if(address > -1){
               pid = a_pid[address];
               px = a_px[address];
               py = a_py[address];
               pz = a_pz[address];
               vx = a_vx[address];
               vy = a_vy[address];
               vz = a_vz[address];
               energy = e_final;
               
               el_lv.SetPxPyPzE(px,py,pz,energy);

               vpho_lv = e_lv - el_lv;

               Q2 = - vpho_lv.Mag2(); //((e_lv - el_lv) * (e_lv - el_lv));
               w_lv = vpho_lv + p_lv; 
               w = w_lv.Mag();

              // q_hist->Fill(Q2);
               w_hist->Fill(w);
               e_hist->Fill(pid);
               xy_hist->Fill(px,py);
            }

            int rho_counter = 0;

            for(int i = 0; i < npip; i++){                               //This double for loop will compare each pip to each pim in an event
               for(int j = 0; j < npim; j++){                            //This is primarly done using TLorentzsVector (the commented parts are my original attempt to calculate)

                  pip_lv.SetPxPyPzE(a_pip_px[i], a_pip_py[i], a_pip_pz[i], a_pip_energy[i]);
                  pim_lv.SetPxPyPzE(a_pim_px[j], a_pim_py[j], a_pim_pz[j], a_pim_energy[j]);

                  rho_lv = pip_lv + pim_lv;

                  t_lv = (vpho_lv - rho_lv);
                 
                  
                  rm_no->Fill(rho_lv.Mag());
                  //std::cout << rho_lv.Mag() << std::endl;
                  //std::cout << w << " | " << Q2 << std::endl;
                  if(w > 2 && Q2 > 1){
                     rm_w->Fill(rho_lv.Mag());
                     float tb = t_lv.Mag2();

                     if(-tb > 0.1 && -tb < 0.5){
                        float Zhb = rho_lv.E()/vpho_lv.E();
                        rm_wt->Fill(rho_lv.Mag());

                        if(Zhb > 0.9){
                           rm_wtz->Fill(rho_lv.Mag());
                        }
                     }
                  }


                  
                 // if(0.6 <= rho_lv.Mag() && rho_lv.Mag() <= 1.0){   //should be 0.77 with +- 0.20[ first standard deviation from hist] error
                     
                     a_rho_mass[rho_counter] = rho_lv.Mag();
                     a_rho_px[rho_counter] = rho_lv.Px();
                     a_rho_py[rho_counter] = rho_lv.Py();
                     a_rho_pz[rho_counter] = rho_lv.Pz();
                     a_rho_energy[rho_counter] = rho_lv.E();
                     rho_counter++;
                  //}
               } 
            }

            event_counter++;
         
            if(rho_counter > 0){
               for(int i = 0; i < rho_counter; i++){

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

                  if(w > 2 ){//&& lc <= 0.5 ){
                     w_cut++;
                     if(-t > 0.1 && -t < 0.5){
                        //z_hist->Fill(Zh);
                        t_cut++;
                        if(Zh > 0.9){
                           z_cut++;
		           //edaq_ev->AddPoint(Q2,rho_counter);
                           if(Q2 >= 1 && Q2 <2){
         			//if(rho_mass <= 0.6 || rho_mass >= 1.0){
                              		p_hist->Fill(rho_mass);
				//}
                            }else if(Q2 >= 2 && Q2 <2.5){
                              ip_hist->Fill(rho_mass);
                            }else if(Q2 >= 2.5 && Q2 <3){
    				//if(rho_mass >= 0.6 && rho_mass <= 1.0){
                              		iip_hist->Fill(rho_mass);
				//}
                           }else if(Q2 >= 3 && Q2 <3.5){
                              iiip_hist->Fill(rho_mass);
                           }else if(Q2 >= 3.5 && Q2 <4.5){
                              ivp_hist->Fill(rho_mass);
                           }else if(Q2 >= 4.5 && Q2 <6){
                              vp_hist->Fill(rho_mass);
                           //} else if(Q2 >= 4.5 && Q2 <6){
                           //   vip_hist->Fill(rho_mass);
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
      /*
      float fcup = 0.0;
      float daq = 0.0;
      float daq_e = 0.0;
      float beamCharge = 0.0;
      int event_num = 0;

      float helfcup = hel.getFloat("fcup", 0);
      float helfcupgated = hel.getFloat("fcupgated", 0);
      float helclock = hel.getFloat("clock", 0);
      float helclockgated = hel.getFloat("clockgated", 0);
      
      
      fcup = scaler.getFloat("fcup",0);
      daq = scaler.getFloat("livetime", 0);
      daq_e = events.getDouble("liveTime", 0);
      int daq_ev = config.getInt("event",0);
      float g_beam = scaler.getFloat("fcupgated",0);

      //if (g_beam != 0){
      //	hel_a.show();
      //}

      long int timestamp = config.getLong("timestamp",0); // / pow(10,9);

      float e_beamch = events.getFloat("beamCharge",0);
      ebeam_ev->AddPoint(timestamp,e_beamch);

      hel_fcup->AddPoint(timestamp, helfcup);
      hel_fcupgated->AddPoint(timestamp, helfcupgated);
      hel_clock->AddPoint(timestamp, helclock);
      hel_clockgated->AddPoint(timestamp, helclockgated);

      //std::cout << "fcupgated = " << g_beam << " | beamcharge = " << e_beamch << " | fcup = " << fcup << std::endl;
      //std::cout << daq_ev << std::endl;
      //if (daq >= 0.09 && daq <= 1.0){
        // beamCharge = fcup * daq;
        // event_num = config.getInt("event",0);
         //sdaq_ev->AddPoint(timestamp,daq);

//         edaq_ev->AddPoint(timestamp,daq_e);    currently being used for Nuclear Transparency graph

         //fcup_ev->AddPoint(event_num, beamCharge);
         //fcupg_ev->AddPoint(event_num, g_beam);
      //} else if(fcup == daq) {
         beamCharge = fcup; 
         event_num = config.getInt("event",0); 
         //fcup_ev->AddPoint(event_num, beamCharge);
         //fcupg_ev->AddPoint(event_num, g_beam); 
      //}
      
	   
        // if(pid_e_count == 1){
        //    e_num++;
        // }
        // if(pid_e_count == 1 && pid_pip_count > 0){
         //   e_pi_num++;
        // }
        // if(pid_e_count == 1 && pid_pim_count > 0){
        //    e_neg_num++;
        // }

      r++;
      i++;
      beamCharge = fcup;
      if(beamCharge != 0){
         std::string s_ev = std::to_string(event_num);
         std::string s_bC = std::to_string(beamCharge);
         std::string evbc = s_ev + "a" + s_bC;
         eve_bC.push_back(evbc);
         //fcup_ev->AddPoint(event_num, beamCharge);
         fcupg_ev->AddPoint(timestamp, g_beam);
      }
     */ 
      
   }
    /*
   int n = eve_bC.size();
   for(int i = 0; i < n; i++){
      std::string temp_i_name = eve_bC[i];
      int a_pos = temp_i_name.find("a");
      int leg_i = temp_i_name.length();
      //std::string st_e = temp_i_name.erase(a_pos, leg_i);
      int i_eve = std::stoi(temp_i_name.erase(a_pos, leg_i));

      //std::cout << i_eve << std::endl;

      for(int j = i+1; j < 200; j++){
         std::string temp_j_name = eve_bC[j];
         int aj_pos = temp_j_name.find("a");
         int leg_j = temp_j_name.length();
         //std::string st_ej = temp_j_name.erase(aj_pos, leg_j);
         int j_eve = std::stoi(temp_j_name.erase(aj_pos, leg_j));
         
         if(i_eve > j_eve){
            std::string temp = eve_bC[i];
            eve_bC[i] = eve_bC[j];
            eve_bC[j] = temp;
         }
      }
   }
   

   bigger = -10000000000.0;
   total_bC = 0.0;
   for(int d = 0; d < eve_bC.size(); d++){
      float p_d = 0.0;

      std::string temp_bc_name = eve_bC[d];
      int ab_pos = temp_bc_name.find("a");
      //std::string st_bC = temp_bc_name.erase(0,ab_pos+1);
      float d_bC = std::stof(temp_bc_name.erase(0,ab_pos+1));

      //std::cout << d_bC << std::endl;

      std::string temp_ev_name = eve_bC[d];
      //int ae_pos = temp_ev_name.find("a");
      int leg_ev = temp_ev_name.length();
      //std::string st_ev = temp_ev_name.erase(ab_pos, leg_ev);
      int d_eve = std::stoi(temp_ev_name.erase(ab_pos, leg_ev));

      if(d_bC > bigger && d_bC != 0){
         if(c == 0){
            bigger = d_bC;
            c_bC.push_back(d_bC);
          
            c++;
         } else {
            bigger = d_bC;
            c_bC.push_back(d_bC);
            
            p_d = c_bC[c] - c_bC[c-1];
	         
            total_bC = total_bC + p_d;
               
            c++;
         }
      }
   */
   }
     
      std::cout << "p = " << p + 1 << std::endl;
   }
//}
   
  
   e_hist->GetXaxis()->SetTitle("e");
   e_hist->Draw();
   e_hist->Write();

   pip_hist->GetXaxis()->SetTitle("pip");
   pip_hist->Draw();
   pip_hist->Write();

   pim_hist->GetXaxis()->SetTitle("pim");
   pim_hist->Draw();
   pim_hist->Write();

   q_hist->GetXaxis()->SetTitle("Q2");
   q_hist->Draw();
   q_hist->Write();
   
   qz_hist->GetXaxis()->SetTitle("Q2");
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
   xy_hist->GetXaxis()->SetTitle("px [GeV]");
   xy_hist->GetYaxis()->SetTitle("py [GeV]");
   xy_hist->Draw();
   xy_hist->Write();

   ivx_hist->GetXaxis()->SetTitle("vx [cm]");
   TF1 *ivx_func = new TF1("ivx_fit",fitf,-1,1,6);
   ivx_func->SetParameters(1200,ivx_hist->GetMean(),ivx_hist->GetRMS(), -40000000, 4000000, 0);
   ivx_func->SetParLimits(0,500, 1500);
   ivx_func->SetParLimits(1,-0.2,0.2);
   //ivx_func->SetParLimits(2, );
   ivx_func->SetParLimits(3,-50000000, 50000000);
   ivx_func->SetParLimits(4,-50000000, 50000000);
   ivx_func->SetParLimits(5,-50000000, 50000000);
   ivx_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //ivx_func->GetYaxis()->SetLogy();
   //ivx_hist->Fit("ivx_fit","","",-0.8,0.8);
   ivx_hist->Draw();
   ivx_hist->Write();


   
   ivy_hist->GetXaxis()->SetTitle("vy [cm]");
   TF1 *ivy_func = new TF1("ivy_fit",fitfy,-1,1,6);
   ivy_func->SetParameters(1200,ivy_hist->GetMean(),ivy_hist->GetRMS(), -13400000, 0, 75000);
   ivy_func->SetParLimits(0,500, 1500);
   //ivy_func->SetParLimits(1,-0.2,0.2);
   //ivy_func->SetParLimits(2, );
   ivy_func->SetParLimits(3,-2000000, -15000);
   ivy_func->SetParLimits(4,-150,200);
   ivy_func->SetParLimits(5,0,150000);
   ivy_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //ivy_hist->Fit("ivy_fit","","",0.4,0.6);
   ivy_hist->Draw();
   ivy_hist->Write();

   ixy_hist->GetXaxis()->SetTitle("px [GeV]");
   ixy_hist->GetYaxis()->SetTitle("py [GeV]");
   ixy_hist->Draw();
   ixy_hist->Write();

   ivz_hist->GetXaxis()->SetTitle("vz [cm]");
   ivz_hist->Draw();
   ivz_hist->Write();



//std::cout << iivx_hist->GetRMS() << std::endl;

   iivx_hist->GetXaxis()->SetTitle("vx [cm]");
   TF1 *iivx_func = new TF1("iivx_fit",fitf,-1,1,6);
   iivx_func->SetParameters(1200,iivx_hist->GetMean(),iivx_hist->GetRMS(), -100000, 50000, 70000);
   iivx_func->SetParLimits(0,500, 1500);
   iivx_func->SetParLimits(1,-0.2,0.2);
   //iivx_func->SetParLimits(2, );
   iivx_func->SetParLimits(3,-150000, -500);
   iivx_func->SetParLimits(4,250,75000);
   iivx_func->SetParLimits(5,500,150000);
   iivx_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //iivx_hist->Fit("iivx_fit","","",-0.1,0.1);
   iivx_hist->Draw();
   iivx_hist->Write();

   iivy_hist->GetXaxis()->SetTitle("vy [cm]");
   TF1 *iivy_func = new TF1("iivy_fit",fitfy,-1,1,6);
   iivy_func->SetParameters(1200,iivy_hist->GetMean(),iivy_hist->GetRMS(), -134000000, 0, 75000);
   iivy_func->SetParLimits(0,500, 1500);
   iivy_func->SetParLimits(1,-0.2,0.2);
   //iivy_func->SetParLimits(2, );
   iivy_func->SetParLimits(3,-150000000, -500);
   iivy_func->SetParLimits(4,0,750);
   iivy_func->SetParLimits(5,100,100000);
   iivy_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //iivy_func->GetYaxis()->SetLogy();
   //iivy_hist->Fit("iivy_fit","","",0.3,0.7);
   iivy_hist->Draw();
   iivy_hist->Write();

   iixy_hist->GetXaxis()->SetTitle("P [GeV]");
   iixy_hist->GetYaxis()->SetTitle("E_{tot} / P");
   iixy_hist->Draw();
   iixy_hist->Write();

   iivz_hist->GetXaxis()->SetTitle("vz [cm]");
   iivz_hist->Draw();
   iivz_hist->Write();





   iiivx_hist->GetXaxis()->SetTitle("vx [cm]");
   TF1 *iiivx_func = new TF1("iiivx_fit",fitf,-1,1,6);
   iiivx_func->SetParameters(1200,iiivx_hist->GetMean(),iiivx_hist->GetRMS(), -100000, 50000, 70000);
   iiivx_func->SetParLimits(0,500, 1500);
   iiivx_func->SetParLimits(1,-0.2,0.2);
   //iiivx_func->SetParLimits(2, );
   iiivx_func->SetParLimits(3,-150000, -500);
   iiivx_func->SetParLimits(4,250,75000);
   iiivx_func->SetParLimits(5,500,150000);
   iiivx_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //iiivx_func->GetYaxis()->SetLogy();
   //iiivx_hist->Fit("iiivx_fit","","",-0.1,0.1);
   iiivx_hist->Draw();
   iiivx_hist->Write();

   iiivy_hist->GetXaxis()->SetTitle("vy [cm]");
   TF1 *iiivy_func = new TF1("iiivy_fit",fitfy,-1,1,6);
   iiivy_func->SetParameters(1200,iiivy_hist->GetMean(),iiivy_hist->GetRMS(), -13400000, 0, 75000);
   iiivy_func->SetParLimits(0,500, 1500);
   iiivy_func->SetParLimits(1,-0.2,0.2);
   //iiivy_func->SetParLimits(2, );
   iiivy_func->SetParLimits(3,-15000000, -500);
   iiivy_func->SetParLimits(4,0,750);
   iiivy_func->SetParLimits(5,100,100000);
   iiivy_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //iiivy_func->GetYaxis()->SetLogy();
   //iiivy_hist->Fit("iiivy_fit","","",0.3,0.7);
   iiivy_hist->Draw();
   iiivy_hist->Write();

   iiixy_hist->GetXaxis()->SetTitle("vx [cm]");
   iiixy_hist->GetYaxis()->SetTitle("vy [cm]");
   iiixy_hist->Draw();
   iiixy_hist->Write();

   iiivz_hist->GetXaxis()->SetTitle("vz [cm]");
   iiivz_hist->Draw();
   iiivz_hist->Write();





   ivvx_hist->GetXaxis()->SetTitle("vx [cm]");
   TF1 *ivvx_func = new TF1("ivvx_fit",fitf,-1,1,6);
   ivvx_func->SetParameters(1200,ivvx_hist->GetMean(),ivvx_hist->GetRMS(), -100000, 50000, 70000);
   ivvx_func->SetParLimits(0,500, 1500);
   ivvx_func->SetParLimits(1,-0.2,0.2);
   //ivvx_func->SetParLimits(2, );
   ivvx_func->SetParLimits(3,-150000, -500);
   ivvx_func->SetParLimits(4,250,75000);
   ivvx_func->SetParLimits(5,500,150000);
   ivvx_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //ivvx_func->GetYaxis()->SetLogy();
   //ivvx_hist->Fit("ivvx_fit","","",-0.1,0.1);
   ivvx_hist->Draw();
   ivvx_hist->Write();

   ivvy_hist->GetXaxis()->SetTitle("vy [cm]");
   TF1 *ivvy_func = new TF1("ivvy_fit",fitfy,-1,1,6);
   ivvy_func->SetParameters(1200,ivvy_hist->GetMean(),ivvy_hist->GetRMS(), -13400000, 0, 75000);
   ivvy_func->SetParLimits(0,500, 1500);
   ivvy_func->SetParLimits(1,-0.2,0.2);
   //ivvy_func->SetParLimits(2, );
   ivvy_func->SetParLimits(3,-15000000, -500);
   ivvy_func->SetParLimits(4,0,750);
   ivvy_func->SetParLimits(5,100,100000);
   ivvy_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //ivvy_func->GetYaxis()->SetLogy();
   //ivvy_hist->Fit("ivvy_fit","","",0.3,0.7);
   ivvy_hist->Draw();
   ivvy_hist->Write();

   ivxy_hist->GetXaxis()->SetTitle("vx [cm]");
   ivxy_hist->GetYaxis()->SetTitle("vy [cm]");
   ivxy_hist->Draw();
   ivxy_hist->Write();

   ivvz_hist->GetXaxis()->SetTitle("vz [cm]");
   ivvz_hist->Draw();
   ivvz_hist->Write();





   vvx_hist->GetXaxis()->SetTitle("vx [cm]");
   TF1 *vvx_func = new TF1("vvx_fit",fitf,-1,1,6);
   vvx_func->SetParameters(1200,vvx_hist->GetMean(),vvx_hist->GetRMS(), -100000, 50000, 70000);
   vvx_func->SetParLimits(0,500, 1500);
   vvx_func->SetParLimits(1,-0.2,0.2);
   //vvx_func->SetParLimits(2, );
   vvx_func->SetParLimits(3,-150000, -500);
   vvx_func->SetParLimits(4,250,75000);
   vvx_func->SetParLimits(5,500,150000);
   vvx_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //vvx_func->GetYaxis()->SetLogy();
   //vvx_hist->Fit("vvx_fit","","",-0.1,0.1);
   vvx_hist->Draw();
   vvx_hist->Write();

   vvy_hist->GetXaxis()->SetTitle("vy [cm]");
   TF1 *vvy_func = new TF1("vvy_fit",fitfy,-1,1,6);
   vvy_func->SetParameters(1200,vvy_hist->GetMean(),vvy_hist->GetRMS(), -13400000, 0, 75000);
   vvy_func->SetParLimits(0,500, 1500);
   vvy_func->SetParLimits(1,-0.2,0.2);
   //vvy_func->SetParLimits(2, );
   vvy_func->SetParLimits(3,-15000000, -500);
   vvy_func->SetParLimits(4,0,750);
   vvy_func->SetParLimits(5,100,100000);
   vvy_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //vvy_func->GetYaxis()->SetLogy();
   //vvy_hist->Fit("vvy_fit","","",0.3,0.7);
   vvy_hist->Draw();
   vvy_hist->Write();

   //vxy_hist->GetXaxis()->SetTitle("vx [cm]");
   //vxy_hist->GetYaxis()->SetTitle("vy [cm]");
   //vxy_hist->Draw();
   //vxy_hist->Write();

   vvz_hist->GetXaxis()->SetTitle("vz [cm]");
   vvz_hist->Draw();
   vvz_hist->Write();

 



   vivx_hist->GetXaxis()->SetTitle("vx [cm]");
   TF1 *vivx_func = new TF1("vivx_fit",fitf,-1,1,6);
   vivx_func->SetParameters(1200,vivx_hist->GetMean(),vivx_hist->GetRMS(), -100000, 50000, 70000);
   vivx_func->SetParLimits(0,500, 1500);
   vivx_func->SetParLimits(1,-0.2,0.2);
   //vivx_func->SetParLimits(2, );
   vivx_func->SetParLimits(3,-150000, -500);
   vivx_func->SetParLimits(4,250,75000);
   vivx_func->SetParLimits(5,500,150000);
   vivx_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //vivx_func->GetYaxis()->SetLogy();
   //vivx_hist->Fit("vivx_fit","","",-0.1,0.1);
   vivx_hist->Draw();
   vivx_hist->Write();

   vivy_hist->GetXaxis()->SetTitle("vy [cm]");
   TF1 *vivy_func = new TF1("vivy_fit",fitfy,-1,1,6);
   vivy_func->SetParameters(1200,vivy_hist->GetMean(),vivy_hist->GetRMS(), -13400000, 0, 75000);
   vivy_func->SetParLimits(0,500, 1500);
   vivy_func->SetParLimits(1,-0.2,0.2);
   //vivy_func->SetParLimits(2, );
   vivy_func->SetParLimits(3,-15000000, -500);
   vivy_func->SetParLimits(4,0,750);
   vivy_func->SetParLimits(5,100,100000);
   vivy_func->SetParNames("Constant","Mean_Value","Sigma","c1","c2","c3");
   //vivy_func->GetYaxis()->SetLogy();
   //vivy_hist->Fit("ivx_fit","","",0.3,0.7);
   vivy_hist->Draw();
   vivy_hist->Write();

   vixy_hist->GetXaxis()->SetTitle("vx [cm]");
   vixy_hist->GetYaxis()->SetTitle("vy [cm]");
   vixy_hist->Draw();
   vixy_hist->Write();

   vivz_hist->GetXaxis()->SetTitle("vz [cm]");
   vivz_hist->Draw();
   vivz_hist->Write();

   //ivx_hist->Draw();
   //iivx_hist->Draw("SAME");
   //iiivx_hist->Draw("SAME");
   //ivvx_hist->Draw("SAME");
   //vvx_hist->Draw("SAME");
   //vivx_hist->Draw("SAME");
   //c1->Modified(); c1->Update();
   //c1->SaveAs("zz_19125_vx.pdf");
   //vxs->Write("nostack 0lego1 PFC");


/*
   vxs->Add(ivx_hist);
   vxs->Add(iivx_hist);
   vxs->Add(iiivx_hist);
   vxs->Add(ivvx_hist);
   vxs->Add(vvx_hist);
   vxs->Add(vivx_hist);
   vxs->Draw("nostack");
   c1->Modified(); c1->Update();
   c1->SaveAs("zz_19125_vx.pdf");




   ivy_hist->Draw();
   iivy_hist->Draw("SAME");
   iiivy_hist->Draw("SAME");
   ivvy_hist->Draw("SAME");
   vvy_hist->Draw("SAME");
   vivy_hist->Draw("SAME");
   c2->Modified(); c2->Update();
   c2->SaveAs("zz_19125_vy.pdf");
   //vys->Write("nostack 0lego1 PFC");

   
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
   
   }
   return 0;
}



