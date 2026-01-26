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

// This version of hipo file reader will list the first 10 events and save in root tree


#include <cstdlib>
#include <iostream>
#include "TLorentzVector.h"
//#include "reader.h" 
#include "TFile.h"
#include "TTree.h"
#include <TH1D.h>
#include <TH2F.h>
#include "TCanvas.h"
#include "TStyle.h"

									//This header file is found in /group/clas12/package/hipo/1.8/hipo4/



											//argc stands for argument count and argv stands for argument vector
											//argc is the number of arguments pointed to by argv
int main() {                   //int argc char** argv


   TCanvas *c1 = new TCanvas(); 
   TStyle *gstyle = new TStyle();
   gstyle->SetPalette(kRainBow); 
   
   TH2F *w_hist = new TH2F("Q2 vs w", "Q2 vs w", 100, 0, 10, 100, 1, 5);
   w_hist->SetStats(0);

   TH2F *t_hist = new TH2F("Q2 vs -t", "Q2 vs -t", 100, 0, 10, 100, -5, 5);
   t_hist->SetStats(0);

   TH2F *lc_hist = new TH2F("Q2 vs lc", "Q2 vs lc", 100, 0, 10, 100, 0, 5);
   lc_hist->SetStats(0);

   TH2F *Zh_hist = new TH2F("Q2 vs Zh", "Q2 vs Zh", 100, 0, 10, 100, 0, 2);
   Zh_hist->SetStats(0);

   TH2F *xy_hist = new TH2F("px vs py", "px vs py", 100, -3, 3, 100, -3, 3);
   xy_hist->SetStats(0);
   
   TH1F *cuq_hist = new TH1F("Cu Q2", "Cu Q2", 300, 0, 10);
   TH1F *cuw_hist = new TH1F("Cu w", "Cu w", 300, 1, 5);
   TH1F *cut_hist = new TH1F("Cu t", "Cu t", 300, -10, 10);
   TH1F *culc_hist = new TH1F("Cu lc", "Cu lc", 300, 0, 5);
   TH1F *cuzh_hist = new TH1F("Cu zh", "Cu zh", 300, 0, 5);
   
   std::cout << "1" << std::endl;
   TH1F *ivx_hist = new TH1F("sector 1 vx", "Sector 1 vx", 300, -1, 1);
       ivx_hist->SetFillColor(kCyan+1);
       ivx_hist->SetStats(0);
   TH1D *ivy_hist = new TH1D("sector 1 vy", "Sector 1 vy", 300, -1, 1);
       ivy_hist->SetFillColor(kCyan+1);
       ivy_hist->SetStats(0);
   TH2F *ixy_hist = new TH2F("sector 1 px .vs py", "sector 1 px .vs py", 200, -1, 1, 200, -1, 1);
   TH1D *ivz_hist = new TH1D("sector 1 vz", "Sector 1 vz", 300, -20, 20);

   TH1F *iivx_hist = new TH1F("sector 2 vx", "Sector 2 vx", 300, -1, 1);
	iivx_hist->SetFillColor(kOrange+1);
	iivx_hist->SetStats(0);
   TH1D *iivy_hist = new TH1D("sector 2 vy", "Sector 2 vy", 300, -1, 1);
	iivy_hist->SetFillColor(kOrange+1);
	iivy_hist->SetStats(0);
   TH2F *iixy_hist = new TH2F("sector 2 vx .vs vy", "" , 500, 0.6, 11, 500, 0, 0.3);
   TH1D *iivz_hist = new TH1D("sector 2 vz", "Sector 2 vz", 300, -20, 20);

   TH1F *iiivx_hist = new TH1F("sector 3 vx", "Sector 3 vx", 300, -1, 1);
	iiivx_hist->SetFillColor(kYellow+1);
	iiivx_hist->SetStats(0);
   TH1D *iiivy_hist = new TH1D("sector 3 vy", "Sector 3 vy", 300, -1, 1);
	iiivy_hist->SetFillColor(kYellow+1);
	iiivy_hist->SetStats(0);
   TH2F *iiixy_hist = new TH2F("sector 3 vx .vs vy", "sector 3 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1D *iiivz_hist = new TH1D("sector 3 vz", "Sector 3 vz", 300, -20, 20);

   TH1F *ivvx_hist = new TH1F("sector 4 vx", "Sector 4 vx", 300, -1, 1);
	ivvx_hist->SetFillColor(kGreen+1);
	ivvx_hist->SetStats(0);
   TH1D *ivvy_hist = new TH1D("sector 4 vy", "Sector 4 vy", 300, -1, 1);
	ivvy_hist->SetFillColor(kGreen+1);
	ivvy_hist->SetStats(0);
   TH2F *ivxy_hist = new TH2F("sector 4 vx .vs vy", "sector 4 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1D *ivvz_hist = new TH1D("sector 4 vz", "Sector 4 vz", 300, -20, 20);

   TH1F *vvx_hist = new TH1F("sector 5 vx", "Sector 5 vx", 300, -1, 1);
	vvx_hist->SetFillColor(kBlue);
	vvx_hist->SetStats(0);
   TH1D *vvy_hist = new TH1D("sector 5 vy", "Sector 5 vy", 300, -1, 1);
	vvy_hist->SetFillColor(kBlue);
	vvy_hist->SetStats(0);
   TH2F *vxy_hist = new TH2F("sector 5 vx .vs vy", "sector 5 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1D *vvz_hist = new TH1D("sector 5 vz", "Sector 5 vz", 300, -20, 20);

   TH1F *vivx_hist = new TH1F("sector 6 vx", "Sector 6 vx", 300, -1, 1);
	vivx_hist->SetFillColor(kMagenta+1);
	vivx_hist->SetStats(0);
   TH1D *vivy_hist = new TH1D("sector 6 vy", "Sector 6 vy", 300, -1, 1);
        vivy_hist->SetFillColor(kMagenta+1);
        vivy_hist->SetStats(0);
   TH2F *vixy_hist = new TH2F("sector 6 vx .vs vy", "sector 6 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1D *vivz_hist = new TH1D("sector 6 vz", "Sector 6 vz", 300, -20, 20);

   //int row; 
   //int pid;
   //float q2;
   float t;
   float Q2;
   float w;
   float lc;
   float Zh;

   float rho_px;
   float rho_py;
   
    // in GeV/c^2
   float m_pro = 0.938; //in Gev/c^2
   //float m_ele = 0.000511; // in GeV/c^2
   float m_pi = 0.1396;
   float m_rho = 0.770;
   
   float a_pip_vx[10];
   float a_pip_vy[10];
   float a_pip_vz[10];
   float a_pip_px[10];
   float a_pip_py[10];
   float a_pip_pz[10];
   float a_pip_energy[10];
   
   float a_pim_vx[10];
   float a_pim_vy[10];
   float a_pim_vz[10];
   float a_pim_px[10];
   float a_pim_py[10];
   float a_pim_pz[10];
   float a_pim_energy[10];

   float phi;
   float vx;
   float vy;
   float vz;
   float px;
   float py;
   float pz;
   float energy;

std::cout << "1.1" << std::endl;

   TFile *input = new TFile("cusn_skimp_cache_18624.root", "read");
std::cout << "1.2" << std::endl;
   TTree *tree = (TTree*)input->Get("electron_tree");
   TTree *piptree = (TTree*)input->Get("pip_tree");
   TTree *pimtree = (TTree*)input->Get("pim_tree");
   //TFile f("cusn_18624_hist.root", "update");
std::cout << "1.3" << std::endl;
TFile *output = new TFile("cusn_18624_hist_2-2.root", "recreate");
std::cout << "1.4" << std::endl;
   tree->SetBranchAddress("phi" ,&phi);
   tree->SetBranchAddress("vx", &vx);
   tree->SetBranchAddress("vy", &vy);
   tree->SetBranchAddress("vz", &vz);
   tree->SetBranchAddress("px", &px);
   tree->SetBranchAddress("py", &py);
   tree->SetBranchAddress("pz", &pz);
   tree->SetBranchAddress("energy", &energy);
   
   piptree->SetBranchAddress("a_pip_px", &a_pip_px);
   piptree->SetBranchAddress("a_pip_py", &a_pip_py);
   piptree->SetBranchAddress("a_pip_pz", &a_pip_pz);
   piptree->SetBranchAddress("a_pip_energy", &a_pip_energy);
   
   pimtree->SetBranchAddress("a_pim_px", &a_pim_px);
   pimtree->SetBranchAddress("a_pim_py", &a_pim_py);
   pimtree->SetBranchAddress("a_pim_pz", &a_pim_pz);
   pimtree->SetBranchAddress("a_pim_energy", &a_pim_energy);
   
   std::cout << "2" << std::endl;
   
    int nentries = tree->GetEntries();


   std::vector<float> sn_vx;
   std::vector<float> sn_vy;
   std::vector<float> sn_vz;
   std::vector<float> sn_px;
   std::vector<float> sn_py;
   std::vector<float> sn_pz;
   std::vector<float> sn_energy;
   
   float sn_pip_vx[nentries][10];
   float sn_pip_vy[nentries][10];
   float sn_pip_vz[nentries][10];
   float sn_pip_px[nentries][10];
   float sn_pip_py[nentries][10];
   float sn_pip_pz[nentries][10];
   float sn_pip_energy[nentries][10];
   
   float sn_pim_vx[nentries][10];
   float sn_pim_vy[nentries][10];
   float sn_pim_vz[nentries][10];
   float sn_pim_px[nentries][10];
   float sn_pim_py[nentries][10];
   float sn_pim_pz[nentries][10];
   float sn_pim_energy[nentries][10];
   
   
   std::vector<float> cu_vx;
   std::vector<float> cu_vy;
   std::vector<float> cu_vz;
   std::vector<float> cu_px;
   std::vector<float> cu_py;
   std::vector<float> cu_pz;
   std::vector<float> cu_energy;
   
   float cu_pip_vx[nentries][10];
   float cu_pip_vy[nentries][10];
   float cu_pip_vz[nentries][10];
   float cu_pip_px[nentries][10];
   float cu_pip_py[nentries][10];
   float cu_pip_pz[nentries][10];
   float cu_pip_energy[nentries][10];
   
   float cu_pim_vx[nentries][10];
   float cu_pim_vy[nentries][10];
   float cu_pim_vz[nentries][10];
   float cu_pim_px[nentries][10];
   float cu_pim_py[nentries][10];
   float cu_pim_pz[nentries][10];
   float cu_pim_energy[nentries][10];
   
      TLorentzVector e_lv;
      e_lv.SetPxPyPzE(0.0, 0.0, 10.5, 10.5);
      
      TLorentzVector el_lv;
      TLorentzVector vpho_lv;
      TLorentzVector w_lv;
      TLorentzVector pip_lv;
      TLorentzVector pim_lv;
      TLorentzVector rho_lv;
      TLorentzVector t_lv;

      TLorentzVector p_lv;
      p_lv.SetPxPyPzE(0.0, 0.0, 0.0, 0.938);

      TLorentzVector te_lv;
      te_lv.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
      
      float ivx_sigma = 0.0558;
      float ivy_sigma = 0.2362;
      float iivx_sigma = 0.2508;
      float iivy_sigma = 0.0739;
      float iiivx_sigma = 0.236;
      float iiivy_sigma = 0.144;
      float ivvx_sigma = 0.0811;
      float ivvy_sigma = 0.2803;
      float vvx_sigma = 0.272;
      float vvy_sigma = 0.1098;
      float vivx_sigma = 0.204;
      float vivy_sigma = 0.171;

      float ivx_mean = 0.011;
      float ivy_mean = -0.168;
      float iivx_mean = 0.138;
      float iivy_mean = -0.294;
      float iiivx_mean = 0.1394;
      float iiivy_mean = -0.172;
      float ivvx_mean = 0.003;
      float ivvy_mean = 0.235;
      float vvx_mean = -0.1323;
      float vvy_mean = -0.1589;
      float vivx_mean = -0.0378;
      float vivy_mean = -0.290;
      
      int cu_count = 0;
      int sn_count = 0;
      
      
   /*   
      
     // TFile *output = new TFile( outFile_char, "recreate"); 

   TTree *cu_tree = new TTree("cu_electron_tree", "cu_electron_tree");
   TTree *cu_pip_tree = new TTree("cu_pip_tree",  "cu_pip_tree");
   TTree *cu_pim_tree = new TTree("cu_pim_tree",  "cu_pim_tree");               //This is the names of each tree
   


   
   cu_tree->Branch("cu_px", &cu_px, "cu_px/F");
   cu_tree->Branch("cu_py", &cu_py, "cu_py/F");
   cu_tree->Branch("cu_pz", &cu_pz, "cu_pz/F");
   cu_tree->Branch("cu_vx", &cu_vx, "cu_vx/F");
   cu_tree->Branch("cu_vy", &cu_vy, "cu_vy/F");
   cu_tree->Branch("cu_vz", &cu_vz, "cu_vz/F");
   cu_tree->Branch("cu_energy", &cu_energy, "cu_energy/F");
   
   
   cu_pip_tree->Branch("cu_pip_px", &cu_pip_px, "cu_pip_px/F");
   cu_pip_tree->Branch("cu_pip_py", &cu_pip_py, "cu_pip_py/F");
   cu_pip_tree->Branch("cu_pip_pz", &cu_pip_pz, "cu_pip_pz/F");
   cu_pip_tree->Branch("cu_pip_vx", &cu_pip_vx, "cu_pip_vx/F");
   cu_pip_tree->Branch("cu_pip_vy", &cu_pip_vy, "cu_pip_vy/F");
   cu_pip_tree->Branch("cu_pip_vz", &cu_pip_vz, "cu_pip_vz/F");
   cu_pip_tree->Branch("cu_pip_energy", &cu_pip_energy, "cu_pip_energy/F");

   
   cu_pim_tree->Branch("cu_pim_px", &cu_pim_px, "cu_pim_px/F");
   cu_pim_tree->Branch("cu_pim_py", &cu_pim_py, "cu_pim_py/F");
   cu_pim_tree->Branch("cu_pim_pz", &cu_pim_pz, "cu_pim_pz/F");
   cu_pim_tree->Branch("cu_pim_vx", &cu_pim_vx, "cu_pim_vx/F");
   cu_pim_tree->Branch("cu_pim_vy", &cu_pim_vy, "cu_pim_vy/F");
   cu_pim_tree->Branch("cu_pim_vz", &cu_pim_vz, "cu_pim_vz/F");
   cu_pim_tree->Branch("cu_pim_energy", &cu_pim_energy, "cu_pim_energy/F");
   
   
   
   
   TTree *sn_tree = new TTree("sn_electron_tree", "sn_electron_tree");
   TTree *sn_pip_tree = new TTree("sn_pip_tree",  "sn_pip_tree");
   TTree *sn_pim_tree = new TTree("sn_pim_tree",  "sn_pim_tree");               //This is the names of each tree
   


   sn_tree->Branch("sn_px", &sn_px, "sn_px/F");
   sn_tree->Branch("sn_py", &sn_py, "sn_py/F");
   sn_tree->Branch("sn_pz", &sn_pz, "sn_pz/F");
   sn_tree->Branch("sn_vx", &sn_vx, "sn_vx/F");
   sn_tree->Branch("sn_vy", &sn_vy, "sn_vy/F");
   sn_tree->Branch("sn_vz", &sn_vz, "sn_vz/F");
   sn_tree->Branch("sn_energy", &sn_energy, "sn_energy/F");
   
   
   sn_pip_tree->Branch("sn_pip_px", &sn_pip_px, "sn_pip_px/F");
   sn_pip_tree->Branch("sn_pip_py", &sn_pip_py, "sn_pip_py/F");
   sn_pip_tree->Branch("sn_pip_pz", &sn_pip_pz, "sn_pip_pz/F");
   sn_pip_tree->Branch("sn_pip_vx", &sn_pip_vx, "sn_pip_vx/F");
   sn_pip_tree->Branch("sn_pip_vy", &sn_pip_vy, "sn_pip_vy/F");
   sn_pip_tree->Branch("sn_pip_vz", &sn_pip_vz, "sn_pip_vz/F");
   sn_pip_tree->Branch("sn_pip_energy", &sn_pip_energy, "sn_pip_energy/F");

   
   sn_pim_tree->Branch("sn_pim_px", &sn_pim_px, "sn_pim_px/F");
   sn_pim_tree->Branch("sn_pim_py", &sn_pim_py, "sn_pim_py/F");
   sn_pim_tree->Branch("sn_pim_pz", &sn_pim_pz, "sn_pim_pz/F");
   sn_pim_tree->Branch("sn_pim_vx", &sn_pim_vx, "sn_pim_vx/F");
   sn_pim_tree->Branch("sn_pim_vy", &sn_pim_vy, "sn_pim_vy/F");
   sn_pim_tree->Branch("sn_pim_vz", &sn_pim_vz, "sn_pim_vz/F");
   sn_pim_tree->Branch("sn_pim_energy", &sn_pim_energy, "sn_pim_energy/F");
   */


  
      for(int i = 0; i < nentries; i++){
            vx = 30;
            vy = 30;
            vz = 30;
            px = 30;
            py = 30;
            pz = 30;
            energy = 30;
      
        std::cout << i << "/" << nentries << std::endl;
            tree->GetEntry(i);
            piptree->GetEntry(i);
            pimtree->GetEntry(i);
            
            el_lv.SetPxPyPzE(px,py,pz,energy);
            
            float sec_dec =  phi; //el_lv.Phi() * 180/M_PI;  //atan(sqrt(pow(px,2) + pow(py,2))/pz);// * 180/M_PI;
                                   //  ivx_hist->Fill(vx);
                 //  std::cout << sec_dec << std::endl; 
                 if(vx >= -5 && vx < 5 && vy >= -5 && vy <= 5) {      
				if(sec_dec >= -20 && sec_dec < 40) {
				ivx_hist->Fill(vx);
   				ivy_hist->Fill(vy);
				
                                if( ((ivx_sigma * -3) + ivx_mean) <= vx && ((ivx_sigma * 3)+ ivx_mean) >= vx && ((ivy_sigma * -3)+ ivy_mean) <= vy && ((ivy_sigma * 3)+ ivy_mean) >= vy){
				ivz_hist->Fill(vz);
				    /* if((-9.5) <= vz && (-6.5) >= vz){			      
    						cu_vx.push_back(vx);
   						cu_vy.push_back(vy);
   						cu_vz.push_back(vz);
   						cu_px.push_back(px);
   						cu_py.push_back(py);
   						cu_pz.push_back(pz);
   						cu_energy.push_back(energy);
                                             for(int m = 0; i < 10; i++){
   						cu_pip_vx[cu_count][m] = a_pip_vx[m];
   						cu_pip_vy[cu_count][m] = a_pip_vy[m];
   						cu_pip_vz[cu_count][m] = a_pip_vz[m];
   						cu_pip_px[cu_count][m] = a_pip_px[m];
   						cu_pip_py[cu_count][m] = a_pip_py[m];
   						cu_pip_pz[cu_count][m] = a_pip_pz[m];
   						cu_pip_energy[cu_count][m] = a_pip_energy[m];
   
   						cu_pim_vx[cu_count][m] = a_pim_vx[m];
   						cu_pim_vy[cu_count][m] = a_pim_vy[m];
   						cu_pim_vz[cu_count][m] = a_pim_vz[m];
   						cu_pim_px[cu_count][m] = a_pim_px[m];
   						cu_pim_py[cu_count][m] = a_pim_py[m];
   						cu_pim_pz[cu_count][m] = a_pim_pz[m];
   						cu_pim_energy[cu_count][m] = a_pim_energy[m];
                                              }
                                              cu_count++;
				     } else if ((-4.5) <= vz && (-1.5) >= vz){
    						sn_vx.push_back(vx);
   						sn_vy.push_back(vy);
   						sn_vz.push_back(vz);
   						sn_px.push_back(px);
   						sn_py.push_back(py);
   						sn_pz.push_back(pz);
   						sn_energy.push_back(energy);
   
   						  for(int m = 0; i < 10; i++){
   						sn_pip_vx[sn_count][m] = a_pip_vx[m];
   						sn_pip_vy[sn_count][m] = a_pip_vy[m];
   						sn_pip_vz[sn_count][m] = a_pip_vz[m];
   						sn_pip_px[sn_count][m] = a_pip_px[m];
   						sn_pip_py[sn_count][m] = a_pip_py[m];
   						sn_pip_pz[sn_count][m] = a_pip_pz[m];
   						sn_pip_energy[sn_count][m] = a_pip_energy[m];
   
   						sn_pim_vx[sn_count][m] = a_pim_vx[m];
   						sn_pim_vy[sn_count][m] = a_pim_vy[m];
   						sn_pim_vz[sn_count][m] = a_pim_vz[m];
   						sn_pim_px[sn_count][m] = a_pim_px[m];
   						sn_pim_py[sn_count][m] = a_pim_py[m];
   						sn_pim_pz[sn_count][m] = a_pim_pz[m];
   						sn_pim_energy[sn_count][m] = a_pim_energy[m];
                                              }
                                              sn_count++;
				     } */
				}
				
                        } else if(sec_dec >= 40 && sec_dec < 100) {
				iivx_hist->Fill(vx);
   				iivy_hist->Fill(vy);
				//iixy_hist->Fill(px, py);
				if( ((iivx_sigma * -3) + iivx_mean) <= vx && ((iivx_sigma * 3)+ iivx_mean) >= vx && ((iivy_sigma * -3)+ iivy_mean) <= vy && ((iivy_sigma * 3)+ iivy_mean) >= vy){
				iivz_hist->Fill(vz);
				/*if((-9.5) <= vz && (-6.5) >= vz){			      
    						cu_vx.push_back(vx);
   						cu_vy.push_back(vy);
   						cu_vz.push_back(vz);
   						cu_px.push_back(px);
   						cu_py.push_back(py);
   						cu_pz.push_back(pz);
   						cu_energy.push_back(energy);
   
   						  for(int m = 0; i < 10; i++){
   						cu_pip_vx[cu_count][m] = a_pip_vx[m];
   						cu_pip_vy[cu_count][m] = a_pip_vy[m];
   						cu_pip_vz[cu_count][m] = a_pip_vz[m];
   						cu_pip_px[cu_count][m] = a_pip_px[m];
   						cu_pip_py[cu_count][m] = a_pip_py[m];
   						cu_pip_pz[cu_count][m] = a_pip_pz[m];
   						cu_pip_energy[cu_count][m] = a_pip_energy[m];
   
   						cu_pim_vx[cu_count][m] = a_pim_vx[m];
   						cu_pim_vy[cu_count][m] = a_pim_vy[m];
   						cu_pim_vz[cu_count][m] = a_pim_vz[m];
   						cu_pim_px[cu_count][m] = a_pim_px[m];
   						cu_pim_py[cu_count][m] = a_pim_py[m];
   						cu_pim_pz[cu_count][m] = a_pim_pz[m];
   						cu_pim_energy[cu_count][m] = a_pim_energy[m];
                                              }
                                              cu_count++;
				     } else if ((-4.5) <= vz && (-1.5) >= vz){
    						sn_vx.push_back(vx);
   						sn_vy.push_back(vy);
   						sn_vz.push_back(vz);
   						sn_px.push_back(px);
   						sn_py.push_back(py);
   						sn_pz.push_back(pz);
   						sn_energy.push_back(energy);
   
   						  for(int m = 0; i < 10; i++){
   						sn_pip_vx[sn_count][m] = a_pip_vx[m];
   						sn_pip_vy[sn_count][m] = a_pip_vy[m];
   						sn_pip_vz[sn_count][m] = a_pip_vz[m];
   						sn_pip_px[sn_count][m] = a_pip_px[m];
   						sn_pip_py[sn_count][m] = a_pip_py[m];
   						sn_pip_pz[sn_count][m] = a_pip_pz[m];
   						sn_pip_energy[sn_count][m] = a_pip_energy[m];
   
   						sn_pim_vx[sn_count][m] = a_pim_vx[m];
   						sn_pim_vy[sn_count][m] = a_pim_vy[m];
   						sn_pim_vz[sn_count][m] = a_pim_vz[m];
   						sn_pim_px[sn_count][m] = a_pim_px[m];
   						sn_pim_py[sn_count][m] = a_pim_py[m];
   						sn_pim_pz[sn_count][m] = a_pim_pz[m];
   						sn_pim_energy[sn_count][m] = a_pim_energy[m];
                                              }
                                              sn_count++;
				     }  */
				}
                        } else if(sec_dec >= 100 && sec_dec < 160) {
				iiivx_hist->Fill(vx);
   				iiivy_hist->Fill(vy);
				iiixy_hist->Fill(px, py);
				if( ((iiivx_sigma * -3) + iiivx_mean) <= vx && ((iiivx_sigma * 3)+ iiivx_mean) >= vx && ((iiivy_sigma * -3)+ iiivy_mean) <= vy && ((iiivy_sigma * 3)+ iiivy_mean) >= vy){
				iiivz_hist->Fill(vz);
			/*	if((-9.5) <= vz && (-6.5) >= vz){			      
    						cu_vx.push_back(vx);
   						cu_vy.push_back(vy);
   						cu_vz.push_back(vz);
   						cu_px.push_back(px);
   						cu_py.push_back(py);
   						cu_pz.push_back(pz);
   						cu_energy.push_back(energy);
   
   						  for(int m = 0; i < 10; i++){
   						cu_pip_vx[cu_count][m] = a_pip_vx[m];
   						cu_pip_vy[cu_count][m] = a_pip_vy[m];
   						cu_pip_vz[cu_count][m] = a_pip_vz[m];
   						cu_pip_px[cu_count][m] = a_pip_px[m];
   						cu_pip_py[cu_count][m] = a_pip_py[m];
   						cu_pip_pz[cu_count][m] = a_pip_pz[m];
   						cu_pip_energy[cu_count][m] = a_pip_energy[m];
   
   						cu_pim_vx[cu_count][m] = a_pim_vx[m];
   						cu_pim_vy[cu_count][m] = a_pim_vy[m];
   						cu_pim_vz[cu_count][m] = a_pim_vz[m];
   						cu_pim_px[cu_count][m] = a_pim_px[m];
   						cu_pim_py[cu_count][m] = a_pim_py[m];
   						cu_pim_pz[cu_count][m] = a_pim_pz[m];
   						cu_pim_energy[cu_count][m] = a_pim_energy[m];
                                              }
                                           cu_count++;
				     } else if ((-4.5) <= vz && (-1.5) >= vz){
    						sn_vx.push_back(vx);
   						sn_vy.push_back(vy);
   						sn_vz.push_back(vz);
   						sn_px.push_back(px);
   						sn_py.push_back(py);
   						sn_pz.push_back(pz);
   						sn_energy.push_back(energy);
   
   						  for(int m = 0; i < 10; i++){
   						sn_pip_vx[sn_count][m] = a_pip_vx[m];
   						sn_pip_vy[sn_count][m] = a_pip_vy[m];
   						sn_pip_vz[sn_count][m] = a_pip_vz[m];
   						sn_pip_px[sn_count][m] = a_pip_px[m];
   						sn_pip_py[sn_count][m] = a_pip_py[m];
   						sn_pip_pz[sn_count][m] = a_pip_pz[m];
   						sn_pip_energy[sn_count][m] = a_pip_energy[m];
   
   						sn_pim_vx[sn_count][m] = a_pim_vx[m];
   						sn_pim_vy[sn_count][m] = a_pim_vy[m];
   						sn_pim_vz[sn_count][m] = a_pim_vz[m];
   						sn_pim_px[sn_count][m] = a_pim_px[m];
   						sn_pim_py[sn_count][m] = a_pim_py[m];
   						sn_pim_pz[sn_count][m] = a_pim_pz[m];
   						sn_pim_energy[sn_count][m] = a_pim_energy[m];
                                              }
                                              sn_count++;
				     }   */
				}
                        } else if(sec_dec >= 160 || sec_dec < -140) {
				ivvx_hist->Fill(vx);
   				ivvy_hist->Fill(vy);
				ivxy_hist->Fill(px, py);
				if( ((ivvx_sigma * -3) + ivvx_mean) <= vx && ((ivvx_sigma * 3)+ ivvx_mean) >= vx && ((ivvy_sigma * -3)+ ivvy_mean) <= vy && ((ivvy_sigma * 3)+ ivvy_mean) >= vy){
				ivvz_hist->Fill(vz);
			/*	if((-9.5) <= vz && (-6.5) >= vz){			      
    						cu_vx.push_back(vx);
   						cu_vy.push_back(vy);
   						cu_vz.push_back(vz);
   						cu_px.push_back(px);
   						cu_py.push_back(py);
   						cu_pz.push_back(pz);
   						cu_energy.push_back(energy);
   
   						  for(int m = 0; i < 10; i++){
   						cu_pip_vx[cu_count][m] = a_pip_vx[m];
   						cu_pip_vy[cu_count][m] = a_pip_vy[m];
   						cu_pip_vz[cu_count][m] = a_pip_vz[m];
   						cu_pip_px[cu_count][m] = a_pip_px[m];
   						cu_pip_py[cu_count][m] = a_pip_py[m];
   						cu_pip_pz[cu_count][m] = a_pip_pz[m];
   						cu_pip_energy[cu_count][m] = a_pip_energy[m];
   
   						cu_pim_vx[cu_count][m] = a_pim_vx[m];
   						cu_pim_vy[cu_count][m] = a_pim_vy[m];
   						cu_pim_vz[cu_count][m] = a_pim_vz[m];
   						cu_pim_px[cu_count][m] = a_pim_px[m];
   						cu_pim_py[cu_count][m] = a_pim_py[m];
   						cu_pim_pz[cu_count][m] = a_pim_pz[m];
   						cu_pim_energy[cu_count][m] = a_pim_energy[m];
                                              }
                                               cu_count++;
				     } else if ((-4.5) <= vz && (-1.5) >= vz){
    						sn_vx.push_back(vx);
   						sn_vy.push_back(vy);
   						sn_vz.push_back(vz);
   						sn_px.push_back(px);
   						sn_py.push_back(py);
   						sn_pz.push_back(pz);
   						sn_energy.push_back(energy);
   
   						  for(int m = 0; i < 10; i++){
   						sn_pip_vx[sn_count][m] = a_pip_vx[m];
   						sn_pip_vy[sn_count][m] = a_pip_vy[m];
   						sn_pip_vz[sn_count][m] = a_pip_vz[m];
   						sn_pip_px[sn_count][m] = a_pip_px[m];
   						sn_pip_py[sn_count][m] = a_pip_py[m];
   						sn_pip_pz[sn_count][m] = a_pip_pz[m];
   						sn_pip_energy[sn_count][m] = a_pip_energy[m];
   
   						sn_pim_vx[sn_count][m] = a_pim_vx[m];
   						sn_pim_vy[sn_count][m] = a_pim_vy[m];
   						sn_pim_vz[sn_count][m] = a_pim_vz[m];
   						sn_pim_px[sn_count][m] = a_pim_px[m];
   						sn_pim_py[sn_count][m] = a_pim_py[m];
   						sn_pim_pz[sn_count][m] = a_pim_pz[m];
   						sn_pim_energy[sn_count][m] = a_pim_energy[m];
                                              }
                                              sn_count++;
				     } */
				}
                        } else if(sec_dec >= -140 && sec_dec < -80) {
				vvx_hist->Fill(vx);
   				vvy_hist->Fill(vy);
				vxy_hist->Fill(px, py);
				if( ((vvx_sigma * -3) + vvx_mean) <= vx && ((vvx_sigma * 3)+ vvx_mean) >= vx && ((vvy_sigma * -3)+ vvy_mean) <= vy && ((vvy_sigma * 3)+ vvy_mean) >= vy){
				vvz_hist->Fill(vz);
			/*	if((-9.5) <= vz && (-6.5) >= vz){			      
    						cu_vx.push_back(vx);
   						cu_vy.push_back(vy);
   						cu_vz.push_back(vz);
   						cu_px.push_back(px);
   						cu_py.push_back(py);
   						cu_pz.push_back(pz);
   						cu_energy.push_back(energy);
   
   						  for(int m = 0; i < 10; i++){
   						cu_pip_vx[cu_count][m] = a_pip_vx[m];
   						cu_pip_vy[cu_count][m] = a_pip_vy[m];
   						cu_pip_vz[cu_count][m] = a_pip_vz[m];
   						cu_pip_px[cu_count][m] = a_pip_px[m];
   						cu_pip_py[cu_count][m] = a_pip_py[m];
   						cu_pip_pz[cu_count][m] = a_pip_pz[m];
   						cu_pip_energy[cu_count][m] = a_pip_energy[m];
   
   						cu_pim_vx[cu_count][m] = a_pim_vx[m];
   						cu_pim_vy[cu_count][m] = a_pim_vy[m];
   						cu_pim_vz[cu_count][m] = a_pim_vz[m];
   						cu_pim_px[cu_count][m] = a_pim_px[m];
   						cu_pim_py[cu_count][m] = a_pim_py[m];
   						cu_pim_pz[cu_count][m] = a_pim_pz[m];
   						cu_pim_energy[cu_count][m] = a_pim_energy[m];
                                              }
                                               cu_count++;
				     } else if ((-4.5) <= vz && (-1.5) >= vz){
    						sn_vx.push_back(vx);
   						sn_vy.push_back(vy);
   						sn_vz.push_back(vz);
   						sn_px.push_back(px);
   						sn_py.push_back(py);
   						sn_pz.push_back(pz);
   						sn_energy.push_back(energy);
   
   						  for(int m = 0; i < 10; i++){
   						sn_pip_vx[sn_count][m] = a_pip_vx[m];
   						sn_pip_vy[sn_count][m] = a_pip_vy[m];
   						sn_pip_vz[sn_count][m] = a_pip_vz[m];
   						sn_pip_px[sn_count][m] = a_pip_px[m];
   						sn_pip_py[sn_count][m] = a_pip_py[m];
   						sn_pip_pz[sn_count][m] = a_pip_pz[m];
   						sn_pip_energy[sn_count][m] = a_pip_energy[m];
   
   						sn_pim_vx[sn_count][m] = a_pim_vx[m];
   						sn_pim_vy[sn_count][m] = a_pim_vy[m];
   						sn_pim_vz[sn_count][m] = a_pim_vz[m];
   						sn_pim_px[sn_count][m] = a_pim_px[m];
   						sn_pim_py[sn_count][m] = a_pim_py[m];
   						sn_pim_pz[sn_count][m] = a_pim_pz[m];
   						sn_pim_energy[sn_count][m] = a_pim_energy[m];
                                              }
                                              sn_count++;
				     } */
				}
                        } else if(sec_dec >= -80 && sec_dec < -20) {
				vivx_hist->Fill(vx);
   				vivy_hist->Fill(vy);
				vixy_hist->Fill(px, py);
				if( ((vivx_sigma * -3) + vivx_mean) <= vx && ((vivx_sigma * 3)+ vivx_mean) >= vx && ((vivy_sigma * -3)+ vivy_mean) <= vy && ((vivy_sigma * 3)+ vivy_mean) >= vy){
				vivz_hist->Fill(vz);
			/*	if((-9.5) <= vz && (-6.5) >= vz){			      
    						cu_vx.push_back(vx);
   						cu_vy.push_back(vy);
   						cu_vz.push_back(vz);
   						cu_px.push_back(px);
   						cu_py.push_back(py);
   						cu_pz.push_back(pz);
   						cu_energy.push_back(energy);
   
   						  for(int m = 0; i < 10; i++){
   						cu_pip_vx[cu_count][m] = a_pip_vx[m];
   						cu_pip_vy[cu_count][m] = a_pip_vy[m];
   						cu_pip_vz[cu_count][m] = a_pip_vz[m];
   						cu_pip_px[cu_count][m] = a_pip_px[m];
   						cu_pip_py[cu_count][m] = a_pip_py[m];
   						cu_pip_pz[cu_count][m] = a_pip_pz[m];
   						cu_pip_energy[cu_count][m] = a_pip_energy[m];
   
   						cu_pim_vx[cu_count][m] = a_pim_vx[m];
   						cu_pim_vy[cu_count][m] = a_pim_vy[m];
   						cu_pim_vz[cu_count][m] = a_pim_vz[m];
   						cu_pim_px[cu_count][m] = a_pim_px[m];
   						cu_pim_py[cu_count][m] = a_pim_py[m];
   						cu_pim_pz[cu_count][m] = a_pim_pz[m];
   						cu_pim_energy[cu_count][m] = a_pim_energy[m];
                                              }
                                              cu_count++;
   
				     } else if ((-4.5) <= vz && (-1.5) >= vz){
    						sn_vx.push_back(vx);
   						sn_vy.push_back(vy);
   						sn_vz.push_back(vz);
   						sn_px.push_back(px);
   						sn_py.push_back(py);
   						sn_pz.push_back(pz);
   						sn_energy.push_back(energy);
   
   						  for(int m = 0; i < 10; i++){
   						sn_pip_vx[sn_count][m] = a_pip_vx[m];
   						sn_pip_vy[sn_count][m] = a_pip_vy[m];
   						sn_pip_vz[sn_count][m] = a_pip_vz[m];
   						sn_pip_px[sn_count][m] = a_pip_px[m];
   						sn_pip_py[sn_count][m] = a_pip_py[m];
   						sn_pip_pz[sn_count][m] = a_pip_pz[m];
   						sn_pip_energy[sn_count][m] = a_pip_energy[m];
   
   						sn_pim_vx[sn_count][m] = a_pim_vx[m];
   						sn_pim_vy[sn_count][m] = a_pim_vy[m];
   						sn_pim_vz[sn_count][m] = a_pim_vz[m];
   						sn_pim_px[sn_count][m] = a_pim_px[m];
   						sn_pim_py[sn_count][m] = a_pim_py[m];
   						sn_pim_pz[sn_count][m] = a_pim_pz[m];
   						sn_pim_energy[sn_count][m] = a_pim_energy[m];
                                              }
                                              sn_count++;
				     }   */
				}
                         }
                         }
            
              
               
	}
	
	//TLorentzVector e_lv;
      //e_lv.SetPxPyPzE(0.0, 0.0, 10.5, 10.5);
      std::cout << "2" << std::endl;
      TLorentzVector cu_el_lv;
      TLorentzVector cu_vpho_lv;
      TLorentzVector cu_w_lv;
      TLorentzVector cu_pip_lv;
      TLorentzVector cu_pim_lv;
      TLorentzVector cu_rho_lv;
      TLorentzVector cu_t_lv;
      
      TLorentzVector sn_el_lv;
      TLorentzVector sn_vpho_lv;
      TLorentzVector sn_w_lv;
      TLorentzVector sn_pip_lv;
      TLorentzVector sn_pim_lv;
      TLorentzVector sn_rho_lv;
      TLorentzVector sn_t_lv;

      //TLorentzVector p_lv;
      //p_lv.SetPxPyPzE(0.0, 0.0, 0.0, 0.938);

      //TLorentzVector te_lv;
      //te_lv.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
	/*
	for(int i = 0; i < cu_count; i++){
	
	       cu_el_lv.SetPxPyPzE( cu_px[i], cu_py[i], cu_pz[i], cu_energy[i] );
	
	       cu_vpho_lv = e_lv - cu_el_lv;

               float Q2 = - cu_vpho_lv.Mag2(); //((e_lv - el_lv) * (e_lv - el_lv));
               cu_w_lv = cu_vpho_lv + p_lv; 
               float w = cu_w_lv.Mag();
               
               int npip = 0;
               int npim = 0;

               std::cout << "3" << std::endl;
               for(int n = 0; n < 10; n++){
               	if( cu_pip_px[i][n] != 0){
               		npip++;
               	}
               }
               
               for(int n = 0; n < 10; n++){
               	if( cu_pim_px[i][n] != 0){
               		npim++;
               	}
               }

              // for(int l = 0; l < npip; l++){                               //This double for loop will compare each pip to each pim in an event
              // 	for(int j = 0; j < npim; j++){                            //This is primarly done using TLorentzsVector (the commented parts are my original attempt to calculate)

                  //cu_pip_lv.SetPxPyPzE(cu_pip_px[i][l], cu_pip_py[i][l], cu_pip_pz[i][l], cu_pip_energy[i][l]);
                  //cu_pim_lv.SetPxPyPzE(cu_pim_px[i][j], cu_pim_py[i][j], cu_pim_pz[i][j], cu_pim_energy[i][j]);

                  //cu_rho_lv = cu_pip_lv + cu_pim_lv;

                 
                  //cu_t_lv = (cu_vpho_lv - cu_rho_lv);
                  //t = cu_t_lv.Mag2();

                  //lc = 0.1973*(2 * cu_vpho_lv.E()) / (Q2 + pow(m_rho,2));           //rho_lv.Mag2());  (0.1973)

                  //Zh = cu_rho_lv.E()/cu_vpho_lv.E();

            w_hist->Fill(Q2, w);    
            t_hist->Fill(Q2,-t);
            lc_hist->Fill(lc,Q2);
            Zh_hist->Fill(Q2,Zh);
            
            cuq_hist->Fill(Q2);
            cuw_hist->Fill(w);
            cut_hist->Fill(-t);
            culc_hist->Fill(lc);
            cuzh_hist->Fill(Zh);
               //
//}
//}
               }
               
       for(int i = 0; i < sn_count; i++){
       
       
       }

       
    */
          

   cuq_hist->Draw();
   cuq_hist->Write();
   
   cuw_hist->Draw();
   cuw_hist->Write();
   
   cut_hist->Draw();
   cut_hist->Write();
   
   culc_hist->Draw();
   culc_hist->Write();
   
   cuzh_hist->Draw();
   cuzh_hist->Write();
   
      
   w_hist->GetXaxis()->SetTitle("Q2 [GeV^2]");
   w_hist->GetYaxis()->SetTitle("w [GeV]");
   w_hist->GetZaxis()->SetTitle("Gradiant");
   w_hist->Draw("colz");
   w_hist->Write();

   t_hist->GetXaxis()->SetTitle("Q2 [GeV^2]");
   t_hist->GetYaxis()->SetTitle("-t [GeV^2]");
   t_hist->GetZaxis()->SetTitle("Gradiant");
   t_hist->Draw("colz");
   t_hist->Write();

   lc_hist->GetXaxis()->SetTitle("Q2 [GeV^2]");
   lc_hist->GetYaxis()->SetTitle("lc [fm]");
   lc_hist->GetZaxis()->SetTitle("Gradiant");
   lc_hist->Draw("colz");
   lc_hist->Write();

   Zh_hist->GetXaxis()->SetTitle("Q2 [GeV^2]");
   Zh_hist->GetYaxis()->SetTitle("Zh [GeV]");
   Zh_hist->GetZaxis()->SetTitle("Gradiant");
   Zh_hist->Draw("colz");
   Zh_hist->Write();

   xy_hist->GetXaxis()->SetTitle("Rho_px");
   xy_hist->GetYaxis()->SetTitle("Rho_py");
   xy_hist->GetZaxis()->SetTitle("Gradiant");
   
   
   ivx_hist->GetXaxis()->SetTitle("vx [cm]");
   
   ivx_hist->Draw();
   ivx_hist->Write();


   
   ivy_hist->GetXaxis()->SetTitle("vy [cm]");
   
   ivy_hist->Draw();
   ivy_hist->Write();


   ivz_hist->GetXaxis()->SetTitle("vz [cm]");
   ivz_hist->Draw();
   ivz_hist->Write();



//std::cout << iivx_hist->GetRMS() << std::endl;

   iivx_hist->GetXaxis()->SetTitle("vx [cm]");
   
   iivx_hist->Draw();
   iivx_hist->Write();

   iivy_hist->GetXaxis()->SetTitle("vy [cm]");
   
   iivy_hist->Draw();
   iivy_hist->Write();


   iivz_hist->GetXaxis()->SetTitle("vz [cm]");
   iivz_hist->Draw();
   iivz_hist->Write();





   iiivx_hist->GetXaxis()->SetTitle("vx [cm]");
   
   iiivx_hist->Draw();
   iiivx_hist->Write();

   iiivy_hist->GetXaxis()->SetTitle("vy [cm]");
   
   iiivy_hist->Draw();
   iiivy_hist->Write();

   

   iiivz_hist->GetXaxis()->SetTitle("vz [cm]");
   iiivz_hist->Draw();
   iiivz_hist->Write();





   ivvx_hist->GetXaxis()->SetTitle("vx [cm]");
   
   ivvx_hist->Draw();
   ivvx_hist->Write();

   ivvy_hist->GetXaxis()->SetTitle("vy [cm]");
 
   ivvy_hist->Draw();
   ivvy_hist->Write();

  

   ivvz_hist->GetXaxis()->SetTitle("vz [cm]");
   ivvz_hist->Draw();
   ivvz_hist->Write();





   vvx_hist->GetXaxis()->SetTitle("vx [cm]");
 
   vvx_hist->Draw();
   vvx_hist->Write();

   vvy_hist->GetXaxis()->SetTitle("vy [cm]");

   vvy_hist->Draw();
   vvy_hist->Write();

  
   vvz_hist->GetXaxis()->SetTitle("vz [cm]");
   vvz_hist->Draw();
   vvz_hist->Write();

 



   vivx_hist->GetXaxis()->SetTitle("vx [cm]");
   
   vivx_hist->Draw();
   vivx_hist->Write();

   vivy_hist->GetXaxis()->SetTitle("vy [cm]");
  
   vivy_hist->Draw();
   vivy_hist->Write();

  
   vivz_hist->GetXaxis()->SetTitle("vz [cm]");
   vivz_hist->Draw();
   vivz_hist->Write();
   //xy_hist->Draw("colz");
   
   output->Write();
   output->Close();
   return 0;
   //xy_hist->Write();
}
//### END OF GENERATED CODE
