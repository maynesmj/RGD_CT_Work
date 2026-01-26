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
   
   TH2F *cutz_hist = new TH2F("Theta vs Cu vz", "Theta vs Cu vz", 100, -14, -6, 100, 0, 60);
   cutz_hist->SetStats(0);
   
   TH2F *sntz_hist = new TH2F("Theta vs Sn vz", "Theta vs Sn vz", 100, -6, 0, 100, 0, 60);
   sntz_hist->SetStats(0);
   
   TH1F *cuq_hist = new TH1F("Cu Q2", "", 300, 0, 10);
	cuq_hist->SetStats(0);
   TH1F *cuw_hist = new TH1F("Cu w", "", 300, 1, 5);
	cuw_hist->SetStats(0);
   TH1F *cut_hist = new TH1F("Cu t", "", 300, 0, 1);
	cut_hist->SetStats(0);
   TH1F *culc_hist = new TH1F("Cu lc", "", 300, 0, 5);
	culc_hist->SetStats(0);
   TH1F *cuzh_hist = new TH1F("Cu zh", "", 300, 0, 5);
	cuzh_hist->SetStats(0);
   
   TH1F *bin1_cu = new TH1F("Cu bin1", "Q^{2} (1-2) [GeV^{2}]", 300, 0, 2);
	bin1_cu->SetStats(0);
   TH1F *bin2_cu = new TH1F("Cu bin2", "Q^{2} (2-2.5) [GeV^{2}]", 300, 0, 2);
	bin2_cu->SetStats(0);
   TH1F *bin3_cu = new TH1F("Cu bin3", "Q^{2} (2.5-3) [GeV^{2}]", 300, 0, 2);
	bin3_cu->SetStats(0);
   TH1F *bin4_cu = new TH1F("Cu bin4", "Q^{2} (3-3.5) [GeV^{2}]", 300, 0, 2);
	bin4_cu->SetStats(0);
   TH1F *bin5_cu = new TH1F("Cu bin5", "Q^{2} (3.5-4.5) [GeV^{2}]", 300, 0, 2);
	bin5_cu->SetStats(0);
   TH1F *bin6_cu = new TH1F("Cu bin6", "Q^{2} (4.5-6) [GeV^{2}]", 300, 0, 2);
	bin6_cu->SetStats(0);

   TH1F *snq_hist = new TH1F("Sn Q2", "", 300, 0, 10);
	snq_hist->SetStats(0);
   TH1F *snw_hist = new TH1F("Sn w", "", 300, 1, 5);
	snw_hist->SetStats(0);
   TH1F *snt_hist = new TH1F("Sn t", "", 300, 0, 1);
	snt_hist->SetStats(0);
   TH1F *snlc_hist = new TH1F("Sn lc", "", 300, 0, 5);
	snlc_hist->SetStats(0);
   TH1F *snzh_hist = new TH1F("Sn zh", "", 300, 0, 5);
	snzh_hist->SetStats(0);   

   TH1F *bin1_sn = new TH1F("Sn bin1", "Q^{2} (1-2) [GeV^{2}]", 300, 0, 2);
	bin1_sn->SetStats(0);
   TH1F *bin2_sn = new TH1F("Sn bin2", "Q^{2} (2-2.5) [GeV^{2}]", 300, 0, 2);
	bin2_sn->SetStats(0);
   TH1F *bin3_sn = new TH1F("Sn bin3", "Q^{2} (2.5-3) [GeV^{2}]", 300, 0, 2);
	bin3_sn->SetStats(0);
   TH1F *bin4_sn = new TH1F("Sn bin4", "Q^{2} (3-3.5) [GeV^{2}]", 300, 0, 2);
	bin4_sn->SetStats(0);
   TH1F *bin5_sn = new TH1F("Sn bin5", "Q^{2} (3.5-4.5) [GeV^{2}]", 300, 0, 2);
	bin5_sn->SetStats(0);
   TH1F *bin6_sn = new TH1F("Sn bin6", "Q^{2} (4.5-6) [GeV^{2}]", 300, 0, 2);
	bin6_sn->SetStats(0);
	
    TH1F *sn_full = new TH1F("sn_full", "sn", 200, 0, 2);
	sn_full->SetStats(0);
    TH1F *cu_full = new TH1F("cu_full", "cu", 200, 0, 2);
	cu_full->SetStats(0);

   TH1F *vz_hist = new TH1F("vz_hist", "vz", 300, -20, 0);
        vz_hist->SetStats(0);


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



   TFile *input = new TFile("Cu_Sn_bspot_tree_2.root", "read"); //Cu_Sn_pass0v4_300_tree_2.root


   TFile f("cu_sn_bspot_hist_6.root", "update");

   

   TTree *cu_electron_tree = (TTree*)input->Get("cu_electron_tree");
   TTree *cu_pip_tree = (TTree*)input->Get("cu_pip_tree");
   TTree *cu_pim_tree = (TTree*)input->Get("cu_pim_tree");               //This is the names of each tree
   
   cu_electron_tree->SetBranchAddress("phi" ,&phi);
   //tree->SetBranchAddress("vx", &vx);
   //tree->SetBranchAddress("vy", &vy);
   cu_electron_tree->SetBranchAddress("vz", &vz);
   cu_electron_tree->SetBranchAddress("px", &px);
   cu_electron_tree->SetBranchAddress("py", &py);
   cu_electron_tree->SetBranchAddress("pz", &pz);
   cu_electron_tree->SetBranchAddress("energy", &energy);
   
   //piptree->SetBranchAddress("npip", &npip, "npip/I");
   cu_pip_tree->SetBranchAddress("a_pip_px", &a_pip_px);
   cu_pip_tree->SetBranchAddress("a_pip_py", &a_pip_py);
   cu_pip_tree->SetBranchAddress("a_pip_pz", &a_pip_pz);
   cu_pip_tree->SetBranchAddress("a_pip_energy", &a_pip_energy);
   
   //pimtree->SetBranchAddress("npim", &npim, "npim/I");
   cu_pim_tree->SetBranchAddress("a_pim_px", &a_pim_px);
   cu_pim_tree->SetBranchAddress("a_pim_py", &a_pim_py);
   cu_pim_tree->SetBranchAddress("a_pim_pz", &a_pim_pz);
   cu_pim_tree->SetBranchAddress("a_pim_energy", &a_pim_energy);

   TTree *sn_electron_tree = (TTree*)input->Get("sn_electron_tree");
   TTree *sn_pip_tree = (TTree*)input->Get("sn_pip_tree");
   TTree *sn_pim_tree = (TTree*)input->Get("sn_pim_tree");               //This is the names of each tree
   
   sn_electron_tree->SetBranchAddress("phi" ,&phi);
   //tree->SetBranchAddress("vx", &vx);
   //tree->SetBranchAddress("vy", &vy);
   sn_electron_tree->SetBranchAddress("vz", &vz);
   sn_electron_tree->SetBranchAddress("px", &px);
   sn_electron_tree->SetBranchAddress("py", &py);
   sn_electron_tree->SetBranchAddress("pz", &pz);
   sn_electron_tree->SetBranchAddress("energy", &energy);
   
   //piptree->SetBranchAddress("npip", &npip, "npip/I");
   sn_pip_tree->SetBranchAddress("a_pip_px", &a_pip_px);
   sn_pip_tree->SetBranchAddress("a_pip_py", &a_pip_py);
   sn_pip_tree->SetBranchAddress("a_pip_pz", &a_pip_pz);
   sn_pip_tree->SetBranchAddress("a_pip_energy", &a_pip_energy);
   
   //pimtree->SetBranchAddress("npim", &npim, "npim/I");
   sn_pim_tree->SetBranchAddress("a_pim_px", &a_pim_px);
   sn_pim_tree->SetBranchAddress("a_pim_py", &a_pim_py);
   sn_pim_tree->SetBranchAddress("a_pim_pz", &a_pim_pz);
   sn_pim_tree->SetBranchAddress("a_pim_energy", &a_pim_energy);
   

   //TTree *sn_tree = (TTree*)input->Get("sn_electron_tree");
   //TTree *sn_pip_tree = (TTree*)input->Get("sn_pip_tree");
   //TTree *sn_pim_tree = (TTree*)input->Get("sn_pim_tree");               //This is the names of each tree
   
   
    int nentries = cu_electron_tree->GetEntries();


   
 
   
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
      
      float ivx_sigma = 0.142252;
      float ivy_sigma = 0.1312978;
      float iivx_sigma = 0.1329112;
      float iivy_sigma = 0.1142967;
      float iiivx_sigma = 0.1382390;
      float iiivy_sigma = 0.1418667;
      float ivvx_sigma = 0.1367910;
      float ivvy_sigma = 0.1181533;
      float vvx_sigma = 0.1154537;
      float vvy_sigma = 0.1336475;
      float vivx_sigma = 0.1279556;
      float vivy_sigma = 0.1220071;
      
      int cu_count = 0;
      int sn_count = 0;
      
      
      
      
     


   
  
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
      int tik = 0;
      for(int i = 0; i < nentries; i++){
            //vx = 30;
            //vy = 30;
            //vz = 30;
            px = 30;
            py = 30;
            pz = 30;
            energy = 30;

         std::cout << i << "/" << nentries << std::endl;
         
            cu_electron_tree->GetEntry(i);
            cu_pip_tree->GetEntry(i);
            cu_pim_tree->GetEntry(i);
            tik++;
            el_lv.SetPxPyPzE(px,py,pz,energy);
            vz_hist->Fill(vz);
            float cut = el_lv.Theta() * 180/M_PI;
            std::cout << cut << std::endl;
            cutz_hist->Fill(vz, cut);
            //std::cout << tik << " |py " << py << " |pz " << pz << " |energy " << energy << std::endl;
            //std::cout << a_pip_px[1] << " |pip_py " << a_pip_py[1] << " |pip_pz " << a_pip_pz[1] << " |pip_energy " << a_pip_energy[1] << std::endl;
            //std::cout << a_pim_px[1] << " |pim_py " << a_pim_py[1] << " |pim_pz " << a_pim_pz[1] << " |pim_energy " << a_pim_energy[1] << std::endl;

   //pip_lv.SetPxPyPzE(a_pip_px[l], a_pip_py[l], a_pip_pz[l], a_pip_energy[l]);
   //               pim_lv.SetPxPyPzE(a_pim_px[j], a_pim_py[j], a_pim_pz[j], a_pim_energy[j]);
	       vpho_lv = e_lv - el_lv;

               float Q2 = - vpho_lv.Mag2(); //((e_lv - el_lv) * (e_lv - el_lv));
               w_lv = vpho_lv + p_lv; 
               float w = w_lv.Mag();
               //std::cout << Q2 << " | " << w << std::endl;
               int npip = 0;
               int npim = 0;

               //std::cout << "3" << std::endl;
               for(int n = 0; n < 10; n++){
               	if( a_pip_px[n] != 0){
               		npip++;
               	}
               }
               
               for(int n = 0; n < 10; n++){
               	if( a_pim_px[n] != 0){
               		npim++;
               	}
               }

               for(int l = 0; l < npip; l++){                               //This double for loop will compare each pip to each pim in an event
               	for(int j = 0; j < npim; j++){                            //This is primarly done using TLorentzsVector (the commented parts are my original attempt to calculate)

                  pip_lv.SetPxPyPzE(a_pip_px[l], a_pip_py[l], a_pip_pz[l], a_pip_energy[l]);
                  pim_lv.SetPxPyPzE(a_pim_px[j], a_pim_py[j], a_pim_pz[j], a_pim_energy[j]);

                  rho_lv = pip_lv + pim_lv;

                  float rho_mass = rho_lv.Mag();
                 
                  t_lv = (vpho_lv - rho_lv);
                  t = t_lv.Mag2();

                  lc = 0.1973*(2 * vpho_lv.E()) / (Q2 + pow(m_rho,2));           //rho_lv.Mag2());  (0.1973)

                  Zh = rho_lv.E()/vpho_lv.E();

            w_hist->Fill(Q2, w);    
            t_hist->Fill(Q2,-t);
            lc_hist->Fill(lc,Q2);
            Zh_hist->Fill(Q2,Zh);
            
            cuq_hist->Fill(Q2);
            cuw_hist->Fill(w);
            cut_hist->Fill(-t);
            culc_hist->Fill(lc);
            cuzh_hist->Fill(Zh);

            

            if(w > 2 && lc <= 1.0){
                     if(-t > 0.1 && -t < 0.5){
                        if(Zh > 0.9){
			   cu_full->Fill(rho_mass);
                           if(Q2 >= 1 && Q2 < 2){
                              		bin1_cu->Fill(rho_mass);
                            }else if(Q2 >= 2 && Q2 < 2.5){
                              bin2_cu->Fill(rho_mass);
                            }else if(Q2 >= 2.5 && Q2 < 3){
                              		bin3_cu->Fill(rho_mass);
                           }else if(Q2 >= 3 && Q2 < 3.5){
                              bin4_cu->Fill(rho_mass);
                           }else if(Q2 >= 3.5 && Q2 < 4.5){
                              bin5_cu->Fill(rho_mass);
                     
                           } else if(Q2 >= 4.5 && Q2 < 6){
                              bin6_cu->Fill(rho_mass);
                           } 
                          
                        
                     }
                  }
               }
            //   std::cout << Q2 << " |w " << w << " |t " << t << " |lc " << lc << " |Zh " << Zh << std::endl;
}
}
               
            
              
               
	}
	nentries = sn_electron_tree->GetEntries();
for(int i = 0; i < nentries; i++){
            //vx = 30;
            //vy = 30;
            //vz = 30;
            px = 30;
            py = 30;
            pz = 30;
            energy = 30;

         std::cout << i << "/" << nentries << std::endl;
         
            sn_electron_tree->GetEntry(i);
            sn_pip_tree->GetEntry(i);
            sn_pim_tree->GetEntry(i);
            tik++;
            sn_el_lv.SetPxPyPzE(px,py,pz,energy);
            vz_hist->Fill(vz);
            float snt = sn_el_lv.Theta() * 180/M_PI;
            sntz_hist->Fill(vz, snt);
            //std::cout << tik << " |py " << py << " |pz " << pz << " |energy " << energy << std::endl;
            //std::cout << a_pip_px[1] << " |pip_py " << a_pip_py[1] << " |pip_pz " << a_pip_pz[1] << " |pip_energy " << a_pip_energy[1] << std::endl;
            //std::cout << a_pim_px[1] << " |pim_py " << a_pim_py[1] << " |pim_pz " << a_pim_pz[1] << " |pim_energy " << a_pim_energy[1] << std::endl;

   //pip_lv.SetPxPyPzE(a_pip_px[l], a_pip_py[l], a_pip_pz[l], a_pip_energy[l]);
   //               pim_lv.SetPxPyPzE(a_pim_px[j], a_pim_py[j], a_pim_pz[j], a_pim_energy[j]);
	       sn_vpho_lv = e_lv - sn_el_lv;

               float Q2 = - sn_vpho_lv.Mag2(); //((e_lv - el_lv) * (e_lv - el_lv));
               sn_w_lv = sn_vpho_lv + p_lv; 
               float w = sn_w_lv.Mag();
               //std::cout << Q2 << " | " << w << std::endl;
               int npip = 0;
               int npim = 0;

               //std::cout << "3" << std::endl;
               for(int n = 0; n < 10; n++){
               	if( a_pip_px[n] != 0){
               		npip++;
               	}
               }
               
               for(int n = 0; n < 10; n++){
               	if( a_pim_px[n] != 0){
               		npim++;
               	}
               }

               for(int l = 0; l < npip; l++){                               //This double for loop will compare each pip to each pim in an event
               	for(int j = 0; j < npim; j++){                            //This is primarly done using TLorentzsVector (the commented parts are my original attempt to calculate)

                  pip_lv.SetPxPyPzE(a_pip_px[l], a_pip_py[l], a_pip_pz[l], a_pip_energy[l]);
                  pim_lv.SetPxPyPzE(a_pim_px[j], a_pim_py[j], a_pim_pz[j], a_pim_energy[j]);

                  rho_lv = pip_lv + pim_lv;

                  float rho_mass = rho_lv.Mag();
                 
                  t_lv = (sn_vpho_lv - rho_lv);
                  t = t_lv.Mag2();

                  lc = 0.1973*(2 * sn_vpho_lv.E()) / (Q2 + pow(m_rho,2));           //rho_lv.Mag2());  (0.1973)

                  Zh = rho_lv.E()/sn_vpho_lv.E();

            w_hist->Fill(Q2, w);    
            t_hist->Fill(Q2,-t);
            lc_hist->Fill(lc,Q2);
            Zh_hist->Fill(Q2,Zh);
            
            snq_hist->Fill(Q2);
            snw_hist->Fill(w);
            snt_hist->Fill(-t);
            snlc_hist->Fill(lc);
            snzh_hist->Fill(Zh);
            if(w > 2 && lc <= 1.0){
                     if(-t > 0.1 && -t < 0.5){
                        if(Zh > 0.9){
			   sn_full->Fill(rho_mass);
                           if(Q2 >= 1 && Q2 < 2){
                              		bin1_sn->Fill(rho_mass);
                            }else if(Q2 >= 2 && Q2 < 2.5){
                              bin2_sn->Fill(rho_mass);
                            }else if(Q2 >= 2.5 && Q2 < 3){
                              	bin3_sn->Fill(rho_mass);
                           }else if(Q2 >= 3 && Q2 < 3.5){
                              bin4_sn->Fill(rho_mass);
                           }else if(Q2 >= 3.5 && Q2 < 4.5){
                              bin5_sn->Fill(rho_mass);
                     
                           } else if(Q2 >= 4.5 && Q2 < 6){
                              bin6_sn->Fill(rho_mass);
                           } 
                          
                        
                     }
                  }
               }
            //   std::cout << Q2 << " |w " << w << " |t " << t << " |lc " << lc << " |Zh " << Zh << std::endl;
}
}
               
            
              
               
	}	



   vz_hist->GetXaxis()->SetTitle("vz [cm]");
   vz_hist->Draw();
   vz_hist->Write();
    
          
   cuq_hist->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
   cuq_hist->Draw();
   cuq_hist->Write();
   
   cuw_hist->GetXaxis()->SetTitle("W [GeV]");
   cuw_hist->Draw();
   cuw_hist->Write();
   
   cut_hist->GetXaxis()->SetTitle("-t [GeV^{2}]");
   cut_hist->Draw();
   cut_hist->Write();
   
   culc_hist->GetXaxis()->SetTitle("lc [fm]");
   culc_hist->Draw();
   culc_hist->Write();
   
   cuzh_hist->Draw();
   cuzh_hist->Write();

   bin1_cu->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   bin1_cu->Draw();
   bin1_cu->Write();

   bin2_cu->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   bin2_cu->Draw();
   bin2_cu->Write();

   bin3_cu->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   bin3_cu->Draw();
   bin3_cu->Write();

   bin4_cu->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   bin4_cu->Draw();
   bin4_cu->Write();

   bin5_cu->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   bin5_cu->Draw();
   bin5_cu->Write();

   bin6_cu->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   bin6_cu->Draw();
   bin6_cu->Write();

   snq_hist->Draw();
   snq_hist->Write();
   
   snw_hist->Draw();
   snw_hist->Write();
   
   snt_hist->Draw();
   snt_hist->Write();
   
   snlc_hist->Draw();
   snlc_hist->Write();
   
   snzh_hist->Draw();
   snzh_hist->Write();

   bin1_sn->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   bin1_sn->Draw();
   bin1_sn->Write();

   bin2_sn->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   bin2_sn->Draw();
   bin2_sn->Write();

   bin3_sn->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   bin3_sn->Draw();
   bin3_sn->Write();

   bin4_sn->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   bin4_sn->Draw();
   bin4_sn->Write();

   bin5_sn->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   bin5_sn->Draw();
   bin5_sn->Write();

   bin6_sn->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   bin6_sn->Draw();
   bin6_sn->Write();
   
      
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
   
   cutz_hist->GetXaxis()->SetTitle("vz [cm]");
   cutz_hist->GetYaxis()->SetTitle("Theta [deg]");
   cutz_hist->GetZaxis()->SetTitle("Gradiant");
   cutz_hist->Draw("colz");
   cutz_hist->Write();
   
   sntz_hist->GetXaxis()->SetTitle("vz [cm]");
   sntz_hist->GetYaxis()->SetTitle("Theta [deg]");
   sntz_hist->GetZaxis()->SetTitle("Gradiant");
   sntz_hist->Draw("colz");
   sntz_hist->Write();
   
   cu_full->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   cu_full->Draw();
   cu_full->Write();
   
   sn_full->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   sn_full->Draw();
   sn_full->Write();
   
   
   
   
   return 0;
   //xy_hist->Write();
}
//### END OF GENERATED CODE
