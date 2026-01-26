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
   
   
   TH1F *ivx_hist = new TH1F("sector 1 vx", "Sector 1 vx", 300, -5, 5);
       ivx_hist->SetFillColor(kCyan+1);
       ivx_hist->SetStats(0);
   TH1D *ivy_hist = new TH1D("sector 1 vy", "Sector 1 vy", 300, -5, 5);
       ivy_hist->SetFillColor(kCyan+1);
       ivy_hist->SetStats(0);
   TH2F *ixy_hist = new TH2F("sector 1 px .vs py", "sector 1 px .vs py", 200, -1, 1, 200, -1, 1);
   TH1D *ivz_hist = new TH1D("sector 1 vz", "Sector 1 vz", 300, -20, 20);

   TH1F *iivx_hist = new TH1F("sector 2 vx", "Sector 2 vx", 300, -5, 5);
	iivx_hist->SetFillColor(kOrange+1);
	iivx_hist->SetStats(0);
   TH1D *iivy_hist = new TH1D("sector 2 vy", "Sector 2 vy", 300, -5, 5);
	iivy_hist->SetFillColor(kOrange+1);
	iivy_hist->SetStats(0);
   TH2F *iixy_hist = new TH2F("sector 2 vx .vs vy", "" , 500, 0.6, 11, 500, 0, 0.3);
   TH1D *iivz_hist = new TH1D("sector 2 vz", "Sector 2 vz", 300, -20, 20);

   TH1F *iiivx_hist = new TH1F("sector 3 vx", "Sector 3 vx", 300, -5, 5);
	iiivx_hist->SetFillColor(kYellow+1);
	iiivx_hist->SetStats(0);
   TH1D *iiivy_hist = new TH1D("sector 3 vy", "Sector 3 vy", 300, -5, 5);
	iiivy_hist->SetFillColor(kYellow+1);
	iiivy_hist->SetStats(0);
   TH2F *iiixy_hist = new TH2F("sector 3 vx .vs vy", "sector 3 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1D *iiivz_hist = new TH1D("sector 3 vz", "Sector 3 vz", 300, -20, 20);

   TH1F *ivvx_hist = new TH1F("sector 4 vx", "Sector 4 vx", 300, -5, 5);
	ivvx_hist->SetFillColor(kGreen+1);
	ivvx_hist->SetStats(0);
   TH1D *ivvy_hist = new TH1D("sector 4 vy", "Sector 4 vy", 300, -5, 5);
	ivvy_hist->SetFillColor(kGreen+1);
	ivvy_hist->SetStats(0);
   TH2F *ivxy_hist = new TH2F("sector 4 vx .vs vy", "sector 4 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1D *ivvz_hist = new TH1D("sector 4 vz", "Sector 4 vz", 300, -20, 20);

   TH1F *vvx_hist = new TH1F("sector 5 vx", "Sector 5 vx", 300, -5, 5);
	vvx_hist->SetFillColor(kBlue);
	vvx_hist->SetStats(0);
   TH1D *vvy_hist = new TH1D("sector 5 vy", "Sector 5 vy", 300, -5, 5);
	vvy_hist->SetFillColor(kBlue);
	vvy_hist->SetStats(0);
   TH2F *vxy_hist = new TH2F("sector 5 vx .vs vy", "sector 5 vx .vs vy", 200, -1, 1, 200, -1, 1);
   TH1D *vvz_hist = new TH1D("sector 5 vz", "Sector 5 vz", 300, -20, 20);

   TH1F *vivx_hist = new TH1F("sector 6 vx", "Sector 6 vx", 300, -5, 5);
	vivx_hist->SetFillColor(kMagenta+1);
	vivx_hist->SetStats(0);
   TH1D *vivy_hist = new TH1D("sector 6 vy", "Sector 6 vy", 300, -5, 5);
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
   
   float a_pip_px[10];
   float a_pip_py[10];
   float a_pip_pz[10];
   float a_pip_energy[10];
   
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



   TFile *input = new TFile("cusn_skimp_18564.root", "read");

   TTree *tree = (TTree*)input->Get("electron_tree");
   TTree *piptree = (TTree*)input->Get("pip_tree");
   TTree *pimtree = (TTree*)input->Get("pim_tree");
   TFile f("hcusn_18564_hist.root", "update");

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
   
   
    int nentries = tree->GetEntries();

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

      for(int i = 0; i < nentries; i++){
            vx = 0;
            vy = 0;
            vz = 0;
            px = 0;
            py = 0;
            pz = 0;
            energy = 0;
      
         
            tree->GetEntry(i);
            piptree->GetEntry(i);
            pimtree->GetEntry(i);
            
            el_lv.SetPxPyPzE(px,py,pz,energy);
            
            float sec_dec =  phi; //el_lv.Phi() * 180/M_PI;  //atan(sqrt(pow(px,2) + pow(py,2))/pz);// * 180/M_PI;
                                   //  ivx_hist->Fill(vx);
                 //  std::cout << sec_dec << std::endl;       
				if(vx >= -5 && vx < 5 && vy >= -5 && vy <= 5) {
				ivx_hist->Fill(vx);
   				ivy_hist->Fill(vy);
				
                                //if( (ivx_sigma * -3) <= vx && (ivx_sigma * 3) >= vx && (ivy_sigma * -3) <= vy && (ivy_sigma * 3) >= vy){
				ivz_hist->Fill(vz);
				//}
                        } else if(sec_dec >= 40 && sec_dec < 100) {
				iivx_hist->Fill(vx);
   				iivy_hist->Fill(vy);
				//iixy_hist->Fill(px, py);
				//if( (iivx_sigma * -3) <= vx && (iivx_sigma * 3) >= vx && (iivy_sigma * -3) <= vy && (iivy_sigma * 3) >= vy){
				iivz_hist->Fill(vz);
				//}
                        } else if(sec_dec >= 100 && sec_dec < 160) {
				iiivx_hist->Fill(vx);
   				iiivy_hist->Fill(vy);
				iiixy_hist->Fill(px, py);
				//if( (iiivx_sigma * -3) <= vx && (iiivx_sigma * 3) >= vx && (iiivy_sigma * -3) <= vy && (iiivy_sigma * 3) >= vy){
				iiivz_hist->Fill(vz);
				//}
                        } else if(sec_dec >= 160 && sec_dec < -140) {
				ivvx_hist->Fill(vx);
   				ivvy_hist->Fill(vy);
				ivxy_hist->Fill(px, py);
				//if( (ivvx_sigma * -3) <= vx && (ivvx_sigma * 3) >= vx && (ivvy_sigma * -3) <= vy && (ivvy_sigma * 3) >= vy){
				ivvz_hist->Fill(vz);
				//}
                        } else if(sec_dec >= -140 && sec_dec < -80) {
				vvx_hist->Fill(vx);
   				vvy_hist->Fill(vy);
				vxy_hist->Fill(px, py);
				//if( (vvx_sigma * -3) <= vx && (vvx_sigma * 3) >= vx && (vvy_sigma * -3) <= vy && (vvy_sigma * 3) >= vy){
				vvz_hist->Fill(vz);
				//}
                        } else if(sec_dec >= -80 && sec_dec < -20) {
				vivx_hist->Fill(vx);
   				vivy_hist->Fill(vy);
				vixy_hist->Fill(px, py);
				//if( (vivx_sigma * -3) <= vx && (vivx_sigma * 3) >= vx && (vivy_sigma * -3) <= vy && (vivy_sigma * 3) >= vy){
				vivz_hist->Fill(vz);
				//}
                         }
            
            
               

               vpho_lv = e_lv - el_lv;

               float Q2 = - vpho_lv.Mag2(); //((e_lv - el_lv) * (e_lv - el_lv));
               w_lv = vpho_lv + p_lv; 
               float w = w_lv.Mag();

               //q_hist->Fill(Q2);
              // w_hist->Fill(w);
            
                 // pip_lv.SetPxPyPzE(pip_px, pip_py, pip_pz, pip_energy);
                  //pim_lv.SetPxPyPzE(pim_px, pim_py, pim_pz, pim_energy);

                  //rho_lv = pip_lv + pim_lv;

                 
                  //t_lv = (vpho_lv - rho_lv);
                  //t = t_lv.Mag2();

                  //lc = 0.1973*(2 * vpho_lv.E()) / (Q2 + pow(m_rho,2));           //rho_lv.Mag2());  (0.1973)

                  //eZh = rho_lv.E()/vpho_lv.E();

            //xy_hist->Fill(rho_px, rho_py);
            w_hist->Fill(Q2,w);
            t_hist->Fill(Q2,-t);
            lc_hist->Fill(lc,Q2);
            Zh_hist->Fill(Q2,Zh);
	}
      
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
   
   return 0;
   //xy_hist->Write();
}
