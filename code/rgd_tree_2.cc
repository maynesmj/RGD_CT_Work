#include <cstdlib>
#include <iostream>
#include <string>
#include "reader.h" 
#include "TFile.h"
#include "TTree.h"
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

int main(){


std::vector<std::string> v3_a;
   int num = 0;


   //TFile f("zzzzz_rgd_03b_cxc.root", "update");

    std::ifstream v3("cusn_calv9_skim.list");
         std::string line;
   while(std::getline(v3, line)){
	std::stringstream ss(line);
        std::string hip;
        std::getline(ss, hip, ',');
        v3_a.push_back(hip);
        num++;
}

 TFile *output = new TFile( "cusn_skim_pass0v7_tree_test2.root", "recreate"); 


int pid = 0;
   float phi = 200.0;
   float px = 333.0;
   float py = 333.0;
   float pz = 333.0;
   float vx = 333.0;
   float vy = 333.0;
   float vz = 333.0;
   float vt = 333.0;
   float beta = 333.0;
   float energy = 333.0;
   
   float px1 = 333.0;
   float py1 = 333.0;
   float pz1 = 333.0;
   float vx1 = 333.0;
   float vy1 = 333.0;
   float vz1 = 333.0;
   
   float px2 = 333.0;
   float py2 = 333.0;
   float pz2 = 333.0;
   float vx2 = 333.0;
   float vy2 = 333.0;
   float vz2 = 333.0;
   
   float px3 = 333.0;
   float py3 = 333.0;
   float pz3 = 333.0;
   float vx3 = 333.0;
   float vy3 = 333.0;
   float vz3 = 333.0;
   
   float px4 = 333.0;
   float py4 = 333.0;
   float pz4 = 333.0;
   float vx4 = 333.0;
   float vy4 = 333.0;
   float vz4 = 333.0;
   
   float px5 = 333.0;
   float py5 = 333.0;
   float pz5 = 333.0;
   float vx5 = 333.0;
   float vy5 = 333.0;
   float vz5 = 333.0;
   
   float px6 = 333.0;
   float py6 = 333.0;
   float pz6 = 333.0;
   float vx6 = 333.0;
   float vy6 = 333.0;
   float vz6 = 333.0;

   float pip_phi = 200.0;
   int pip_pid = 0;
   float pip_px = 333.0;
   float pip_py = 333.0;                          //These values will be fill the root tree
   float pip_pz = 333.0;
   float pip_vx = 333.0;
   float pip_vy = 333.0;                          
   float pip_vz = 333.0;
   float pip_vt = 333.0;
   float pip_beta = 333.0;
   float pip_energy = 333.0;
   
   // in GeV/c^2
   float m_pro = 0.938; //in Gev/c^2
   //float m_ele = 0.000511; // in GeV/c^2
   float m_pi = 0.1396;
   float m_rho = 0.770;

   float pim_phi = 200.0;
   int pim_pid = 0;
   float pim_px = 333.0;
   float pim_py = 333.0;
   float pim_pz = 333.0;
   float pim_vx = 333.0;
   float pim_vy = 333.0;
   float pim_vz = 333.0;
   float pim_vt = 333.0;
   float pim_beta = 333.0;
   float pim_energy = 333.0;

float rho_mass = 333.0;
   float rho_px = 333.0;
   float rho_py = 333.0;
   float rho_pz = 333.0;
   float rho_vx = 333.0;
   float rho_vy = 333.0;
   float rho_vz = 333.0;
   float rho_energy = 333.0;

   float Q2 = 0.0;
   float t = 0.0;
   float lc = 0.0;
   float Zh = 0.0;

   int npip = 0;
   int npim = 0;

  int a_pid[10];
         float a_px[10]; 
         float a_py[10];
         float a_pz[10];
         float a_vx[10];
         float a_vy[10];
         float a_vz[10];
         float a_vt[10];
         float a_beta[10];
         float a_energy[10];
         
         

         float a_rho_mass[10];
         float a_rho_px[10];
         float a_rho_py[10];
         float a_rho_pz[10];
         float a_rho_energy[10];

         int a_pip_pid[10];
         float a_pip_px[10];
         float a_pip_py[10];
         float a_pip_pz[10];
         float a_pip_vx[10];
         float a_pip_vy[10];
         float a_pip_vz[10];
         float a_pip_vt[10];
         float a_pip_beta[10];
         float a_pip_energy[10];
         
         float a_pip_px1[10];
         float a_pip_py1[10];
         float a_pip_pz1[10];
         float a_pip_vx1[10];
         float a_pip_vy1[10];
         float a_pip_vz1[10];
         
         float a_pip_px2[10];
         float a_pip_py2[10];
         float a_pip_pz2[10];
         float a_pip_vx2[10];
         float a_pip_vy2[10];
         float a_pip_vz2[10];
         
         float a_pip_px3[10];
         float a_pip_py3[10];
         float a_pip_pz3[10];
         float a_pip_vx3[10];
         float a_pip_vy3[10];
         float a_pip_vz3[10];
         
         float a_pip_px4[10];
         float a_pip_py4[10];
         float a_pip_pz4[10];
         float a_pip_vx4[10];
         float a_pip_vy4[10];
         float a_pip_vz4[10];
         
         float a_pip_px5[10];
         float a_pip_py5[10];
         float a_pip_pz5[10];
         float a_pip_vx5[10];
         float a_pip_vy5[10];
         float a_pip_vz5[10];
         
         float a_pip_px6[10];
         float a_pip_py6[10];
         float a_pip_pz6[10];
         float a_pip_vx6[10];
         float a_pip_vy6[10];
         float a_pip_vz6[10];

         int a_pim_pid[10];                                                     //This section declares the arrays that saves the data of one event. Each will be rewritten for each event.
         float a_pim_px[10];
         float a_pim_py[10];
         float a_pim_pz[10];
         float a_pim_vx[10];
         float a_pim_vy[10];
         float a_pim_vz[10];
         float a_pim_vt[10];
         float a_pim_beta[10];
         float a_pim_energy[10];
         
         float a_pim_px1[10];
         float a_pim_py1[10];
         float a_pim_pz1[10];
         float a_pim_vx1[10];
         float a_pim_vy1[10];
         float a_pim_vz1[10];
         
         float a_pim_px2[10];
         float a_pim_py2[10];
         float a_pim_pz2[10];
         float a_pim_vx2[10];
         float a_pim_vy2[10];
         float a_pim_vz2[10];
         
         float a_pim_px3[10];
         float a_pim_py3[10];
         float a_pim_pz3[10];
         float a_pim_vx3[10];
         float a_pim_vy3[10];
         float a_pim_vz3[10];
         
         float a_pim_px4[10];
         float a_pim_py4[10];
         float a_pim_pz4[10];
         float a_pim_vx4[10];
         float a_pim_vy4[10];
         float a_pim_vz4[10];
         
         float a_pim_px5[10];
         float a_pim_py5[10];
         float a_pim_pz5[10];
         float a_pim_vx5[10];
         float a_pim_vy5[10];
         float a_pim_vz5[10];
         
         float a_pim_px6[10];
         float a_pim_py6[10];
         float a_pim_pz6[10];
         float a_pim_vx6[10];
         float a_pim_vy6[10];
         float a_pim_vz6[10];


      TLorentzVector e_lv;
      e_lv.SetPxPyPzE(0.0, 0.0, 10.5, 10.5);

      TLorentzVector p_lv;
      p_lv.SetPxPyPzE(0.0, 0.0, 0.0, m_pro);

      TLorentzVector te_lv;
      te_lv.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

TLorentzVector pip_lv;// = new TLorentzVector(0., 0., 0., 0.);                //TLorentzVector(px,py,pz,e)
         TLorentzVector pim_lv;// = new TLorentzVector(0., 0., 0., 0.);            //This declares the lorentze four vector for pions and negative pions
        
         TLorentzVector rho_lv;
         TLorentzVector vpho_lv;
         TLorentzVector t_lv;
         TLorentzVector el_lv;
         TLorentzVector w_lv;

   //int file_num[7] = {0, 1000, 2000, 3000, 4000, 5000, 6724};
   //const char* file_name[7] = {"rgd_cxc_tree_01.root","rgd_cxc_tree_02.root","rgd_cxc_tree_03.root","rgd_cxc_tree_04.root","rgd_cxc_tree_05.root","rgd_cxc_tree_06.root","rgd_cxc_tree_07.root"};

  // int file_num[17] = {0, 500, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500};                                                                        
  // const char* file_name[17] = {"rgd_ld_tree_01.root","rgd_ld_tree_02.root","rgd_ld_tree_03.root","rgd_ld_tree_04.root","rgd_ld_\
tree_05.root","rgd_ld_tree_06.root","rgd_ld_tree_07.root", "rgd_ld_tree_08.root","rgd_ld_tree_09.root","rgd_ld_tree_10.root","rgd_ld_tree_11.root","rgd_ld_\
tree_12.root","rgd_ld_tree_13.root","rgd_ld_tree_14.root","rgd_ld_tree_15.root","rgd_ld_tree_16.root","rgd_ld_tree_17.root"};  

//for(int w = 0; w < 15; w++){
   
      //char* outFile_char = new char[file_name[w].length() + 1];
      //strcpy(outFile_char, file_name[w].c_str());
      //const char* outFile = outFile_char;

 //TH1D *wonphe_hist = new TH1D("Without Cuts", "Without Cuts", 200, 0.2, 1.4);

  // TH1D *wnphe_hist = new TH1D("With W cut", "With W cut", 200, 0.2, 1.4);

 //  TH1D *wopcal_hist = new TH1D("With W & -t Cuts", "With W & -t Cuts", 200, 0.2, 1.4);

 //  TH1D *wpcal_hist = new TH1D("With W, -t, & z Cuts", "With W, -t, & z Cuts", 200, 0.2, 1.4);

 //  TFile f("zzzrm_rg_b_11_cxcsm.root", "update");


  

   TTree *tree = new TTree("electron_tree", "electron_tree");
   TTree *pip_tree = new TTree("pip_tree", "pip_tree");
   TTree *pim_tree = new TTree("pim_tree", "pim_tree");
   //TTree *rho_tree = new TTree("rho_tree",  "rho_tree");
   //TTree *kin_tree = new TTree("kin_tree",  "kin_tree"); 


  

   
   tree->Branch("phi", &phi, "phi/F");
   tree->Branch("px1", &px1, "px1/F");
   tree->Branch("py1", &py1, "py1/F");
   tree->Branch("pz1", &pz1, "pz1/F");
   tree->Branch("vx1", &vx1, "vx1/F");
   tree->Branch("vy1", &vy1, "vy1/F");
   tree->Branch("vz1", &vz1, "vz1/F");
   //tree->Branch("vt1", &vt1, "vt/F");
   //tree->Branch("energy1", &energy1, "energy/F");
   tree->Branch("px2", &px2, "px2/F");
   tree->Branch("py2", &py2, "py2/F");
   tree->Branch("pz2", &pz2, "pz2/F");
   tree->Branch("vx2", &vx2, "vx2/F");
   tree->Branch("vy2", &vy2, "vy2/F");
   tree->Branch("vz2", &vz2, "vz2/F");
   
   tree->Branch("px3", &px3, "px3/F");
   tree->Branch("py3", &py3, "py3/F");
   tree->Branch("pz3", &pz3, "pz3/F");
   tree->Branch("vx3", &vx3, "vx3/F");
   tree->Branch("vy3", &vy3, "vy3/F");
   tree->Branch("vz3", &vz3, "vz3/F");
   
   tree->Branch("px4", &px4, "px4/F");
   tree->Branch("py4", &py4, "py4/F");
   tree->Branch("pz4", &pz4, "pz4/F");
   tree->Branch("vx4", &vx4, "vx4/F");
   tree->Branch("vy4", &vy4, "vy4/F");
   tree->Branch("vz4", &vz4, "vz4/F");
   
   tree->Branch("px5", &px5, "px5/F");
   tree->Branch("py5", &py5, "py5/F");
   tree->Branch("pz5", &pz5, "pz5/F");
   tree->Branch("vx5", &vx5, "vx5/F");
   tree->Branch("vy5", &vy5, "vy5/F");
   tree->Branch("vz5", &vz5, "vz5/F");
   
   tree->Branch("px6", &px6, "px6/F");
   tree->Branch("py6", &py6, "py6/F");
   tree->Branch("pz6", &pz6, "pz6/F");
   tree->Branch("vx6", &vx6, "vx6/F");
   tree->Branch("vy6", &vy6, "vy6/F");
   tree->Branch("vz6", &vz6, "vz6/F");
   

   pip_tree->Branch("npip", &npip, "npip/I");
   pip_tree->Branch("a_pip_px1", a_pip_px1, TString::Format("a_pip_px1[%d]/F", 10));
   pip_tree->Branch("a_pip_py1", a_pip_py1, TString::Format("a_pip_py1[%d]/F", 10));
   pip_tree->Branch("a_pip_pz1", a_pip_pz1, TString::Format("a_pip_pz1[%d]/F", 10));
   pip_tree->Branch("a_pip_vx1", a_pip_vx1, TString::Format("a_pip_vx1[%d]/F", 10));
   pip_tree->Branch("a_pip_vy1", a_pip_vy1, TString::Format("a_pip_vy1[%d]/F", 10));
   pip_tree->Branch("a_pip_vz1", a_pip_vz1, TString::Format("a_pip_vz1[%d]/F", 10));
   pip_tree->Branch("a_pip_vt", a_pip_vt, TString::Format("a_pip_vt[%d]/F", 10));
   pip_tree->Branch("a_pip_beta", a_pip_beta, TString::Format("a_pip_beta[%d]/F", 10));
   pip_tree->Branch("a_pip_energy", a_pip_energy, TString::Format("a_pip_energy[%d]/F", 10));
   
   pip_tree->Branch("a_pip_px2", a_pip_px2, TString::Format("a_pip_px2[%d]/F", 10));
   pip_tree->Branch("a_pip_py2", a_pip_py2, TString::Format("a_pip_py2[%d]/F", 10));
   pip_tree->Branch("a_pip_pz2", a_pip_pz2, TString::Format("a_pip_pz2[%d]/F", 10));
   pip_tree->Branch("a_pip_vx2", a_pip_vx2, TString::Format("a_pip_vx2[%d]/F", 10));
   pip_tree->Branch("a_pip_vy2", a_pip_vy2, TString::Format("a_pip_vy2[%d]/F", 10));
   pip_tree->Branch("a_pip_vz2", a_pip_vz2, TString::Format("a_pip_vz2[%d]/F", 10));
   
   pip_tree->Branch("a_pip_px3", a_pip_px3, TString::Format("a_pip_px3[%d]/F", 10));
   pip_tree->Branch("a_pip_py3", a_pip_py3, TString::Format("a_pip_py3[%d]/F", 10));
   pip_tree->Branch("a_pip_pz3", a_pip_pz3, TString::Format("a_pip_pz3[%d]/F", 10));
   pip_tree->Branch("a_pip_vx3", a_pip_vx3, TString::Format("a_pip_vx3[%d]/F", 10));
   pip_tree->Branch("a_pip_vy3", a_pip_vy3, TString::Format("a_pip_vy3[%d]/F", 10));
   pip_tree->Branch("a_pip_vz3", a_pip_vz3, TString::Format("a_pip_vz3[%d]/F", 10));
   
   pip_tree->Branch("a_pip_px4", a_pip_px4, TString::Format("a_pip_px4[%d]/F", 10));
   pip_tree->Branch("a_pip_py4", a_pip_py4, TString::Format("a_pip_py4[%d]/F", 10));
   pip_tree->Branch("a_pip_pz4", a_pip_pz4, TString::Format("a_pip_pz4[%d]/F", 10));
   pip_tree->Branch("a_pip_vx4", a_pip_vx4, TString::Format("a_pip_vx4[%d]/F", 10));
   pip_tree->Branch("a_pip_vy4", a_pip_vy4, TString::Format("a_pip_vy4[%d]/F", 10));
   pip_tree->Branch("a_pip_vz4", a_pip_vz4, TString::Format("a_pip_vz4[%d]/F", 10));
   
   pip_tree->Branch("a_pip_px5", a_pip_px5, TString::Format("a_pip_px5[%d]/F", 10));
   pip_tree->Branch("a_pip_py5", a_pip_py5, TString::Format("a_pip_py5[%d]/F", 10));
   pip_tree->Branch("a_pip_pz5", a_pip_pz5, TString::Format("a_pip_pz5[%d]/F", 10));
   pip_tree->Branch("a_pip_vx5", a_pip_vx5, TString::Format("a_pip_vx5[%d]/F", 10));
   pip_tree->Branch("a_pip_vy5", a_pip_vy5, TString::Format("a_pip_vy5[%d]/F", 10));
   pip_tree->Branch("a_pip_vz5", a_pip_vz5, TString::Format("a_pip_vz5[%d]/F", 10));
   
   pip_tree->Branch("a_pip_px6", a_pip_px6, TString::Format("a_pip_px6[%d]/F", 10));
   pip_tree->Branch("a_pip_py6", a_pip_py6, TString::Format("a_pip_py6[%d]/F", 10));
   pip_tree->Branch("a_pip_pz6", a_pip_pz6, TString::Format("a_pip_pz6[%d]/F", 10));
   pip_tree->Branch("a_pip_vx6", a_pip_vx6, TString::Format("a_pip_vx6[%d]/F", 10));
   pip_tree->Branch("a_pip_vy6", a_pip_vy6, TString::Format("a_pip_vy6[%d]/F", 10));
   pip_tree->Branch("a_pip_vz6", a_pip_vz6, TString::Format("a_pip_vz6[%d]/F", 10));

   pim_tree->Branch("npim", &npim, "npim/I");
   pim_tree->Branch("a_pim_px1", a_pim_px1, TString::Format("a_pim_px1[%d]/F", 10));
   pim_tree->Branch("a_pim_py1", a_pim_py1, TString::Format("a_pim_py1[%d]/F", 10));
   pim_tree->Branch("a_pim_pz1", a_pim_pz1, TString::Format("a_pim_pz1[%d]/F", 10));
   pim_tree->Branch("a_pim_vx1", a_pim_vx1, TString::Format("a_pim_vx1[%d]/F", 10));
   pim_tree->Branch("a_pim_vy1", a_pim_vy1, TString::Format("a_pim_vy1[%d]/F", 10));
   pim_tree->Branch("a_pim_vz1", a_pim_vz1, TString::Format("a_pim_vz1[%d]/F", 10));
   pim_tree->Branch("a_pim_vt", a_pim_vt, TString::Format("a_pim_vt[%d]/F", 10));
   pim_tree->Branch("a_pim_beta", a_pim_beta, TString::Format("a_pim_beta[%d]/F", 10));
   pim_tree->Branch("a_pim_energy", a_pim_energy, TString::Format("a_pim_energy[%d]/F", 10));
   
   pim_tree->Branch("a_pim_px2", a_pim_px2, TString::Format("a_pim_px2[%d]/F", 10));
   pim_tree->Branch("a_pim_py2", a_pim_py2, TString::Format("a_pim_py2[%d]/F", 10));
   pim_tree->Branch("a_pim_pz2", a_pim_pz2, TString::Format("a_pim_pz2[%d]/F", 10));
   pim_tree->Branch("a_pim_vx2", a_pim_vx2, TString::Format("a_pim_vx2[%d]/F", 10));
   pim_tree->Branch("a_pim_vy2", a_pim_vy2, TString::Format("a_pim_vy2[%d]/F", 10));
   pim_tree->Branch("a_pim_vz2", a_pim_vz2, TString::Format("a_pim_vz2[%d]/F", 10));
   
   pim_tree->Branch("a_pim_px3", a_pim_px3, TString::Format("a_pim_px3[%d]/F", 10));
   pim_tree->Branch("a_pim_py3", a_pim_py3, TString::Format("a_pim_py3[%d]/F", 10));
   pim_tree->Branch("a_pim_pz3", a_pim_pz3, TString::Format("a_pim_pz3[%d]/F", 10));
   pim_tree->Branch("a_pim_vx3", a_pim_vx3, TString::Format("a_pim_vx3[%d]/F", 10));
   pim_tree->Branch("a_pim_vy3", a_pim_vy3, TString::Format("a_pim_vy3[%d]/F", 10));
   pim_tree->Branch("a_pim_vz3", a_pim_vz3, TString::Format("a_pim_vz3[%d]/F", 10));
   
   pim_tree->Branch("a_pim_px4", a_pim_px4, TString::Format("a_pim_px4[%d]/F", 10));
   pim_tree->Branch("a_pim_py4", a_pim_py4, TString::Format("a_pim_py4[%d]/F", 10));
   pim_tree->Branch("a_pim_pz4", a_pim_pz4, TString::Format("a_pim_pz4[%d]/F", 10));
   pim_tree->Branch("a_pim_vx4", a_pim_vx4, TString::Format("a_pim_vx4[%d]/F", 10));
   pim_tree->Branch("a_pim_vy4", a_pim_vy4, TString::Format("a_pim_vy4[%d]/F", 10));
   pim_tree->Branch("a_pim_vz4", a_pim_vz4, TString::Format("a_pim_vz4[%d]/F", 10));
   
   pim_tree->Branch("a_pim_px5", a_pim_px5, TString::Format("a_pim_px5[%d]/F", 10));
   pim_tree->Branch("a_pim_py5", a_pim_py5, TString::Format("a_pim_py5[%d]/F", 10));
   pim_tree->Branch("a_pim_pz5", a_pim_pz5, TString::Format("a_pim_pz5[%d]/F", 10));
   pim_tree->Branch("a_pim_vx5", a_pim_vx5, TString::Format("a_pim_vx5[%d]/F", 10));
   pim_tree->Branch("a_pim_vy5", a_pim_vy5, TString::Format("a_pim_vy5[%d]/F", 10));
   pim_tree->Branch("a_pim_vz5", a_pim_vz5, TString::Format("a_pim_vz5[%d]/F", 10));
   
   pim_tree->Branch("a_pim_px6", a_pim_px6, TString::Format("a_pim_px6[%d]/F", 10));
   pim_tree->Branch("a_pim_py6", a_pim_py6, TString::Format("a_pim_py6[%d]/F", 10));
   pim_tree->Branch("a_pim_pz6", a_pim_pz6, TString::Format("a_pim_pz6[%d]/F", 10));
   pim_tree->Branch("a_pim_vx6", a_pim_vx6, TString::Format("a_pim_vx6[%d]/F", 10));
   pim_tree->Branch("a_pim_vy6", a_pim_vy6, TString::Format("a_pim_vy6[%d]/F", 10));
   pim_tree->Branch("a_pim_vz6", a_pim_vz6, TString::Format("a_pim_vz6[%d]/F", 10));

   //rho_tree->Branch("rho_mass", &rho_mass, "rho_mass/F");
   //rho_tree->Branch("rho_px", &rho_px, "rho_px/F");
   //rho_tree->Branch("rho_py", &rho_py, "rho_py/F");
   //rho_tree->Branch("rho_pz", &rho_pz, "rho_pz/F");
   //rho_tree->Branch("rho_energy", &rho_energy, "rho_energy/F");
   //rho_tree->Branch("Q2", &Q2, "Q2/F");
   //rho_tree->Branch("w", &w, "w/F");
   //kin_tree->Branch("t", &t, "t/F");
   //rho_tree->Branch("lc", &lc, "lc/F");
   //kin_tree->Branch("Zh", &Zh, "Zh/F");




//int file_num[4] = {56, 1, 10, 6724}; //3930}; //4164};//need to add the files for 18309       361 for 18451   200 for v0


for(int p = 34; p < 37; p++){

std::string inFile =  v3_a[p];
std::cout << p << " | " << num-1 << std::endl;

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

     // hipo::bank hel(factory.getSchema("HEL::scaler"));
     // hipo::bank hel_a(factory.getSchema("HEL::online"));

      hipo::bank cher(factory.getSchema("REC::Cherenkov"));
      hipo::bank cal(factory.getSchema("REC::Calorimeter"));

      hipo::event      event; 			



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
         
         int nrows = 0; 
         std::fill_n(a_pid, 10, -330);
         std::fill_n(a_px, 10, -330); 
         std::fill_n(a_py, 10, -330);
         std::fill_n(a_pz, 10, -330);
         std::fill_n(a_vx, 10, -330);
         std::fill_n(a_vy, 10, -330);
         std::fill_n(a_vz, 10, -330);
         std::fill_n(a_vt, 10, -330);
         std::fill_n(a_beta, 10, -330);
         std::fill_n(a_energy, 10, -330);
         
         
      

         std::fill_n(a_pip_pid, 10, -330);
         std::fill_n(a_pip_px, 10, -330);
         std::fill_n(a_pip_py, 10, -330);
         std::fill_n(a_pip_pz, 10, -330);
         std::fill_n(a_pip_vx, 10, -330);
         std::fill_n(a_pip_vy, 10, -330);
         std::fill_n(a_pip_vz, 10, -330);
         std::fill_n(a_pip_vt, 10, -330);
         std::fill_n(a_pip_beta, 10, -330);
         std::fill_n(a_pip_energy, 10, -330);
         
         std::fill_n(a_pip_px1, 10, -330);
         std::fill_n(a_pip_py1, 10, -330);
         std::fill_n(a_pip_pz1, 10, -330);
         std::fill_n(a_pip_vx1, 10, -330);
         std::fill_n(a_pip_vy1, 10, -330);
         std::fill_n(a_pip_vz1, 10, -330);
         
         std::fill_n(a_pip_px2, 10, -330);
         std::fill_n(a_pip_py2, 10, -330);
         std::fill_n(a_pip_pz2, 10, -330);
         std::fill_n(a_pip_vx2, 10, -330);
         std::fill_n(a_pip_vy2, 10, -330);
         std::fill_n(a_pip_vz2, 10, -330);
         
         std::fill_n(a_pip_px3, 10, -330);
         std::fill_n(a_pip_py3, 10, -330);
         std::fill_n(a_pip_pz3, 10, -330);
         std::fill_n(a_pip_vx3, 10, -330);
         std::fill_n(a_pip_vy3, 10, -330);
         std::fill_n(a_pip_vz3, 10, -330);
         
         std::fill_n(a_pip_px4, 10, -330);
         std::fill_n(a_pip_py4, 10, -330);
         std::fill_n(a_pip_pz4, 10, -330);
         std::fill_n(a_pip_vx4, 10, -330);
         std::fill_n(a_pip_vy4, 10, -330);
         std::fill_n(a_pip_vz4, 10, -330);
         
         std::fill_n(a_pip_px5, 10, -330);
         std::fill_n(a_pip_py5, 10, -330);
         std::fill_n(a_pip_pz5, 10, -330);
         std::fill_n(a_pip_vx5, 10, -330);
         std::fill_n(a_pip_vy5, 10, -330);
         std::fill_n(a_pip_vz5, 10, -330);
         
         std::fill_n(a_pip_px6, 10, -330);
         std::fill_n(a_pip_py6, 10, -330);
         std::fill_n(a_pip_pz6, 10, -330);
         std::fill_n(a_pip_vx6, 10, -330);
         std::fill_n(a_pip_vy6, 10, -330);
         std::fill_n(a_pip_vz6, 10, -330);

         std::fill_n(a_pim_pid, 10, -330);                
         std::fill_n(a_pim_px, 10, -330);
         std::fill_n(a_pim_py, 10, -330);
         std::fill_n(a_pim_pz, 10, -330);
         std::fill_n(a_pim_vx, 10, -330);
         std::fill_n(a_pim_vy, 10, -330);
         std::fill_n(a_pim_vz, 10, -330);
         std::fill_n(a_pim_vt, 10, -330);
         std::fill_n(a_pim_beta, 10, -330);
         std::fill_n(a_pim_energy, 10, -330);
         
         std::fill_n(a_pim_px1, 10, -330);
         std::fill_n(a_pim_py1, 10, -330);
         std::fill_n(a_pim_pz1, 10, -330);
         std::fill_n(a_pim_vx1, 10, -330);
         std::fill_n(a_pim_vy1, 10, -330);
         std::fill_n(a_pim_vz1, 10, -330);
         
         std::fill_n(a_pim_px2, 10, -330);
         std::fill_n(a_pim_py2, 10, -330);
         std::fill_n(a_pim_pz2, 10, -330);
         std::fill_n(a_pim_vx2, 10, -330);
         std::fill_n(a_pim_vy2, 10, -330);
         std::fill_n(a_pim_vz2, 10, -330);
         
         std::fill_n(a_pim_px3, 10, -330);
         std::fill_n(a_pim_py3, 10, -330);
         std::fill_n(a_pim_pz3, 10, -330);
         std::fill_n(a_pim_vx3, 10, -330);
         std::fill_n(a_pim_vy3, 10, -330);
         std::fill_n(a_pim_vz3, 10, -330);
         
         std::fill_n(a_pim_px4, 10, -330);
         std::fill_n(a_pim_py4, 10, -330);
         std::fill_n(a_pim_pz4, 10, -330);
         std::fill_n(a_pim_vx4, 10, -330);
         std::fill_n(a_pim_vy4, 10, -330);
         std::fill_n(a_pim_vz4, 10, -330);
         
         std::fill_n(a_pim_px5, 10, -330);
         std::fill_n(a_pim_py5, 10, -330);
         std::fill_n(a_pim_pz5, 10, -330);
         std::fill_n(a_pim_vx5, 10, -330);
         std::fill_n(a_pim_vy5, 10, -330);
         std::fill_n(a_pim_vz5, 10, -330);
         
         std::fill_n(a_pim_px6, 10, -330);
         std::fill_n(a_pim_py6, 10, -330);
         std::fill_n(a_pim_pz6, 10, -330);
         std::fill_n(a_pim_vx6, 10, -330);
         std::fill_n(a_pim_vy6, 10, -330);
         std::fill_n(a_pim_vz6, 10, -330);

        

         nrows = particle.getRows(); 
          pid_e_count = 0;
          pid_pip_count = 0;
          pid_pim_count = 0;
         int address = -1;	                   //the rows start at 0 and increase by into the positive. Initiation at -1 was choosen to show that no row was asigned

         float e_final = 0.0;

         float e_total = 0.0;
         float ppcal = 0.0;
         float pcal = 0.0;
         float nphe = 0.0;
        
       
         //std::cout << nrows << std::endl;
         for(int row = 0; row < nrows; row++){
           // float temp_chi = particle.getFloat("chi2pid",row);
           int temp_s = 0;
           float chi2 = 0.0;
           chi2 = particle.getFloat("chi2pid", row);
           temp_s = particle.getInt("status", row);
           if(abs(temp_s)/2000 == 1){
            if(particle.getInt("pid",row) == 11){
               
		         float temp_vz = 0.0;
    		      
        
		

                     temp_vz = particle.getFloat("vz", row);
		     

                     if(temp_vz >= -11.5 && temp_vz <= 5){
                        if(chi2 > -5 && chi2 < 5){
			               
                        
                          
				
				
           	      		   a_pid[pid_e_count] = 11;
           	      		   a_px[pid_e_count] = particle.getFloat("px",row);
            	   		   a_py[pid_e_count] = particle.getFloat("py",row);
              	   		   a_pz[pid_e_count] = particle.getFloat("pz",row);
                  		   a_vx[pid_e_count] = particle.getFloat("vx",row);
           	      		   a_vy[pid_e_count] = particle.getFloat("vy",row);
           	      		   a_vz[pid_e_count] = temp_vz;
				   a_vt[pid_e_count] = particle.getFloat("vt",row);
           	      		   a_beta[pid_e_count] = particle.getFloat("beta",row);
                  		   a_energy[pid_e_count] = sqrt(pow(a_px[pid_e_count],2) + pow(a_py[pid_e_count],2) + pow(a_pz[pid_e_count],2));
                       		   pid_e_count++;

	               }
               }
            } else if(particle.getInt("pid",row) == 211){
            if(chi2 > -10 && chi2 < 10){
              	a_pip_pid[pid_pip_count] = 211;
              	a_pip_px[pid_pip_count] = particle.getFloat("px",row);
           	   a_pip_py[pid_pip_count] = particle.getFloat("py",row);
              	a_pip_pz[pid_pip_count] = particle.getFloat("pz",row);
               a_pip_vx[pid_pip_count] = particle.getFloat("vx",row);
              	a_pip_vy[pid_pip_count] = particle.getFloat("vy",row);
           	   a_pip_vz[pid_pip_count] = particle.getFloat("vz",row);
		a_pip_vt[pid_pip_count] = particle.getFloat("vt",row);
           	   a_pip_beta[pid_pip_count] = particle.getFloat("beta",row);
               a_pip_energy[pid_pip_count] = sqrt(pow(a_pip_px[pid_pip_count],2) + pow(a_pip_py[pid_pip_count],2) + pow(a_pip_pz[pid_pip_count],2) + pow(m_pi,2));
               pid_pip_count++;
               }
           
            } else if(particle.getInt("pid",row) == -211){
            if(chi2 > -10 && chi2 < 10){
              	a_pim_pid[pid_pim_count] = -211;
              	a_pim_px[pid_pim_count] = particle.getFloat("px",row);
              	a_pim_py[pid_pim_count] = particle.getFloat("py",row);
           	   a_pim_pz[pid_pim_count] = particle.getFloat("pz",row);
               a_pim_vx[pid_pim_count] = particle.getFloat("vx",row);
              	a_pim_vy[pid_pim_count] = particle.getFloat("vy",row);
              	a_pim_vz[pid_pim_count] = particle.getFloat("vz",row);
		a_pim_vt[pid_pim_count] = particle.getFloat("vt",row);
              	a_pim_beta[pid_pim_count] = particle.getFloat("beta",row);
               a_pim_energy[pid_pim_count] = sqrt(pow(a_pim_px[pid_pim_count],2) + pow(a_pim_py[pid_pim_count],2) + pow(a_pim_pz[pid_pim_count],2) + pow(m_pi,2));
               pid_pim_count++;
               }
            }
            }
         }

if( (pid_e_count > 0  && pid_pip_count > 0 && pid_pim_count > 0)){               //This ensures that only the events with all three particle will be looked at

            int ne = pid_e_count;
             npip = pid_pip_count;
             npim = pid_pim_count; 
             
            for(int i = 0; i < npip; i++){
                pip_lv.SetPxPyPzE(a_pip_px[i], a_pip_py[i], a_pip_pz[i], a_pip_energy[i]);
                pip_phi = pip_lv.Phi() * 180/M_PI;
                if(pip_phi > -30 && pip_phi <= 30){
                    a_pip_px1[i] = a_pip_px[i];
                    a_pip_py1[i] = a_pip_py[i];
                    a_pip_pz1[i] = a_pip_pz[i];
                    a_pip_vx1[i] = a_pip_vx[i];
                    a_pip_vy1[i] = a_pip_vy[i];
                    a_pip_vz1[i] = a_pip_vz[i];
                } else if (pip_phi > 30 && pip_phi <= 90){
                    a_pip_px2[i] = a_pip_px[i];
                    a_pip_py2[i] = a_pip_py[i];
                    a_pip_pz2[i] = a_pip_pz[i];
                    a_pip_vx2[i] = a_pip_vx[i];
                    a_pip_vy2[i] = a_pip_vy[i];
                    a_pip_vz2[i] = a_pip_vz[i];
                } else if (pip_phi > 90 && pip_phi <= 150){
                    a_pip_px3[i] = a_pip_px[i];
                    a_pip_py3[i] = a_pip_py[i];
                    a_pip_pz3[i] = a_pip_pz[i];
                    a_pip_vx3[i] = a_pip_vx[i];
                    a_pip_vy3[i] = a_pip_vy[i];
                    a_pip_vz3[i] = a_pip_vz[i];
                } else if (pip_phi > 150 || pip_phi <= -150){
                    a_pip_px4[i] = a_pip_px[i];
                    a_pip_py4[i] = a_pip_py[i];
                    a_pip_pz4[i] = a_pip_pz[i];
                    a_pip_vx4[i] = a_pip_vx[i];
                    a_pip_vy4[i] = a_pip_vy[i];
                    a_pip_vz4[i] = a_pip_vz[i];
                } else if (pip_phi > -150 && pip_phi <= -90){
                    a_pip_px5[i] = a_pip_px[i];
                    a_pip_py5[i] = a_pip_py[i];
                    a_pip_pz5[i] = a_pip_pz[i];
                    a_pip_vx5[i] = a_pip_vx[i];
                    a_pip_vy5[i] = a_pip_vy[i];
                    a_pip_vz5[i] = a_pip_vz[i];
                } else if (pip_phi > -90 && pip_phi <= -30){
                    a_pip_px6[i] = a_pip_px[i];
                    a_pip_py6[i] = a_pip_py[i];
                    a_pip_pz6[i] = a_pip_pz[i];
                    a_pip_vx6[i] = a_pip_vx[i];
                    a_pip_vy6[i] = a_pip_vy[i];
                    a_pip_vz6[i] = a_pip_vz[i];
                }
     
           }  
           
           for(int i = 0; i < npim; i++){
                pim_lv.SetPxPyPzE(a_pim_px[i], a_pim_py[i], a_pim_pz[i], a_pim_energy[i]);
                pim_phi = pim_lv.Phi() * 180/M_PI;
                if(pim_phi > -30 && pim_phi <= 30){
                    a_pim_px1[i] = a_pim_px[i];
                    a_pim_py1[i] = a_pim_py[i];
                    a_pim_pz1[i] = a_pim_pz[i];
                    a_pim_vx1[i] = a_pim_vx[i];
                    a_pim_vy1[i] = a_pim_vy[i];
                    a_pim_vz1[i] = a_pim_vz[i];
                } else if (pim_phi > 30 && pim_phi <= 90){
                    a_pim_px2[i] = a_pim_px[i];
                    a_pim_py2[i] = a_pim_py[i];
                    a_pim_pz2[i] = a_pim_pz[i];
                    a_pim_vx2[i] = a_pim_vx[i];
                    a_pim_vy2[i] = a_pim_vy[i];
                    a_pim_vz2[i] = a_pim_vz[i];
                } else if (pim_phi > 90 && pim_phi <= 150){
                    a_pim_px3[i] = a_pim_px[i];
                    a_pim_py3[i] = a_pim_py[i];
                    a_pim_pz3[i] = a_pim_pz[i];
                    a_pim_vx3[i] = a_pim_vx[i];
                    a_pim_vy3[i] = a_pim_vy[i];
                    a_pim_vz3[i] = a_pim_vz[i];
                } else if (pim_phi > 150 || pim_phi <= -150){
                    a_pim_px4[i] = a_pim_px[i];
                    a_pim_py4[i] = a_pim_py[i];
                    a_pim_pz4[i] = a_pim_pz[i];
                    a_pim_vx4[i] = a_pim_vx[i];
                    a_pim_vy4[i] = a_pim_vy[i];
                    a_pim_vz4[i] = a_pim_vz[i];
                } else if (pim_phi > -150 && pim_phi <= -90){
                    a_pim_px5[i] = a_pim_px[i];
                    a_pim_py5[i] = a_pim_py[i];
                    a_pim_pz5[i] = a_pim_pz[i];
                    a_pim_vx5[i] = a_pim_vx[i];
                    a_pim_vy5[i] = a_pim_vy[i];
                    a_pim_vz5[i] = a_pim_vz[i];
                } else if (pim_phi > -90 && pim_phi <= -30){
                    a_pim_px6[i] = a_pim_px[i];
                    a_pim_py6[i] = a_pim_py[i];
                    a_pim_pz6[i] = a_pim_pz[i];
                    a_pim_vx6[i] = a_pim_vx[i];
                    a_pim_vy6[i] = a_pim_vy[i];
                    a_pim_vz6[i] = a_pim_vz[i];
                }
     
           }                                           
         
         pip_tree->Fill();
         
         pim_tree->Fill();
            
       
         
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
               vt = a_vt[address];
               beta = a_beta[address];
               energy = e_final;
               el_lv.SetPxPyPzE(px,py,pz,energy);
               phi = el_lv.Phi() * 180/M_PI;
               
               if (phi > -30 && phi <= 30){
                   px1 = px;
                   py1 = py;
                   pz1 = pz;
                   vx1 = vx;
                   vy1 = vy;
                   vz1 = vz;
               } else if (phi > 30 && phi <= 90){
                   px2 = px;
                   py2 = py;
                   pz2 = pz;
                   vx2 = vx;
                   vy2 = vy;
                   vz2 = vz;
               } else if (phi > 90 && phi <= 150){
                   px3 = px;
                   py3 = py;
                   pz3 = pz;
                   vx3 = vx;
                   vy3 = vy;
                   vz3 = vz;
               } else if (phi > 150 || phi <= -150){
                   px4 = px;
                   py4 = py;
                   pz4 = pz;
                   vx4 = vx;
                   vy4 = vy;
                   vz4 = vz;
               } else if (phi > -150 && phi <= -90){
                   px5 = px;
                   py5 = py;
                   pz5 = pz;
                   vx5 = vx;
                   vy5 = vy;
                   vz5 = vz;
               } else if (phi > -90 && phi <= -30){
                   px6 = px;
                   py6 = py;
                   pz6 = pz;
                   vx6 = vx;
                   vy6 = vy;
                   vz6 = vz;
               }
               
               tree->Fill();
            }

               
}
}
}
   
   output->Write();
   output->Close();


std::fclose (stdout);
/*
wonphe_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   wonphe_hist->Draw();
   wonphe_hist->Write();

wnphe_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   wnphe_hist->Draw();
   wnphe_hist->Write();

wopcal_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   wopcal_hist->Draw();
   wopcal_hist->Write();

wpcal_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   wpcal_hist->Draw();
   wpcal_hist->Write();
*/
//}

   return 0;
}

