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

int main(){

   TH1D *z_hist = new TH1D("Zh", " z ", 200, 0, 10);
   //bC_hist->SetStats(0);

   int p_ent = 0;

   //TH2F *qz_hist = new TH2F("Q2 vs pz", "Q2 vs pz", 200, 0, 10, 200, 0, 10);

   TH1D *q_hist = new TH1D("Q^{2}", " Q^{2} ", 200, 0, 10);

   //TH2F *qlc_hist = new TH2F("Q2 vs lc", "Q2 vs lc", 200, 0, 10, 200, 0, 20);

   TH1D *w_hist = new TH1D("W (Invariant Mass)", "W ", 200, 0, 5);

   //TH2F *tz_hist = new TH2F("-t vs rho_pz", "-t vs rho_pz", 200, 0, 10, 200, 0, 20);

   TH1D *t_hist = new TH1D("-t", "-t", 200, 0, 20);

   TH1D *wonphe_hist = new TH1D("wonphe", " ", 500, 0, 0.5);

   TH1D *wnphe_hist = new TH1D("wnphe", " ", 500, 0, 0.5);

   TH1D *wopcal_hist = new TH1D("Total Energy", " ", 500, 0, 1);

   TH1D *wpcal_hist = new TH1D("PCAL", " ", 500, 0, 1);

   TH2F *wocut_hist = new TH2F("wocut", "wocut", 500, 0, 20, 500, 0, 1);

   TH2F *wcut_hist = new TH2F("Energy", "", 500, 0.6, 11, 500, 0, 0.3);
   
    
  



   
   TH1D *ip_hist = new TH1D("Q^{2} (2->2.5) GeV^{2}", "Q^{2} (2->2.5) GeV^{2}", 200, 0, 2);
       ip_hist->SetFillColor(kCyan+1);
   TH1D *p_hist = new TH1D("Q^{2} (1->2) GeV^{2}", "Q^{2} (1->2) GeV^{2}", 200, 0, 2);
       p_hist->SetFillColor(kGreen-1);

   
   TH1D *iip_hist = new TH1D("Q^{2} (2.5->3) GeV^{2}", "Q^{2} (2.5->3) GeV^{2}", 200, 0, 2);
	iip_hist->SetFillColor(kOrange+1);
   

   
   TH1D *iiip_hist = new TH1D("Q^{2} (3->3.5) GeV^{2}", "Q^{2} (3->3.5) GeV^{2}", 200, 0, 2);
	iiip_hist->SetFillColor(kYellow+1);
   

  
   TH1D *ivp_hist = new TH1D("Q^{2} (3.5->4.5) GeV^{2}", "Q^{2} (3.5->4.5) GeV^{2}", 200, 0, 2);
	ivp_hist->SetFillColor(kGreen+1);
   

   
   TH1D *vip_hist = new TH1D("Q^{2} (4.5->6) GeV^{2}", "Q^{2} (4.5->6) GeV^{2}", 200, 0, 2);
        vip_hist->SetFillColor(kMagenta+1);



   std::vector<std::string> v3_a;
   int num = 0;


   TFile f("zzzzrm_rg_b_11_cxc.root", "update");

    std::ifstream v3("v3.list");
         std::string line;
   while(std::getline(v3, line)){
	std::stringstream ss(line);
        std::string hip;
        std::getline(ss, hip, ',');
        v3_a.push_back(hip);
        num++;
}


   //std::ofstream fw("actual_2_18318_inbend_yield.txt", std::ofstream::out);   //011270 fall & 011361 spring

  
   int file_num[4] = {56, 1, 10, 6700};//need to add the files for 18309       361 for 18451   200 for v0

   
//std::string inFile;

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

      for(int p = 0; p < 6700; p++){
         
std::string inFile;
            inFile =  v3_a[p];
         
       
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


	int r = 0;
	int i = 0;
	int c = 0;

      beam_energy = 10.5;
   

	
   std::vector<float> c_bC;
   std::vector<int> c_ev;
   std::vector<float> bC;
   std::vector<int> ev;

   

   //std::vector<std::string> eve_bC;

   
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
         int nrows = 0;
        


      while(reader.next() == true){
         reader.read(event);
      
         event.getStructure(particle);
	      event.getStructure(events);
         event.getStructure(config);
         event.getStructure(scaler);
         //event.getStructure(hel);
        // event.getStructure(hel_a);
         event.getStructure(cher);
         event.getStructure(cal);
         
       // particle.show();
        // events.show();
        //hel.show();
        //cher.show();
        //cal.show();

         nrows = 0; 
        
   	 
       nrows = particle.getRows();
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
               //std::cout << address << std::endl;
               el_lv.SetPxPyPzE(px,py,pz,energy);

               vpho_lv = e_lv - el_lv;

               Q2 = - vpho_lv.Mag2(); //((e_lv - el_lv) * (e_lv - el_lv));
               w_lv = vpho_lv + p_lv; 
               w = w_lv.Mag();

               //q_hist->Fill(Q2);
               w_hist->Fill(w);
               
            }

            int rho_counter = 0;

            for(int i = 0; i < npip; i++){                               //This double for loop will compare each pip to each pim in an event
               for(int j = 0; j < npim; j++){                            //This is primarly done using TLorentzsVector (the commented parts are my original attempt to calculate)

                  pip_lv.SetPxPyPzE(a_pip_px[i], a_pip_py[i], a_pip_pz[i], a_pip_energy[i]);
                  pim_lv.SetPxPyPzE(a_pim_px[j], a_pim_py[j], a_pim_pz[j], a_pim_energy[j]);

                  rho_lv = pip_lv + pim_lv;

                 
                  t_lv = (vpho_lv - rho_lv);
                  t = t_lv.Mag2();

                  lc = 0.1973*(2 * vpho_lv.E()) / (Q2 + pow(m_rho,2));           //rho_lv.Mag2());  (0.1973)

                  Zh = rho_lv.E()/vpho_lv.E();

                  rho_mass = a_rho_mass[i];
                  rho_px = a_rho_px[i];
                  rho_py = a_rho_py[i];
                  rho_pz = a_rho_pz[i];
                  rho_energy = a_rho_energy[i];

                  t_hist->Fill(-t);
                  //z_hist->Fill(Zh);
                  q_hist->Fill(Q2);


                  



                  if(w > 2 && lc <= 1.0){
                     if(-t > 0.1 && -t < 0.5){
                        if(Zh > 0.9){
			   
                           if(Q2 >= 1 && Q2 < 2){
                              		p_hist->Fill(rho_mass);
                            }else if(Q2 >= 2 && Q2 < 2.5){
                              ip_hist->Fill(rho_mass);
                            }else if(Q2 >= 2.5 && Q2 < 3){
                              		iip_hist->Fill(rho_mass);
                           }else if(Q2 >= 3 && Q2 < 3.5){
                              iiip_hist->Fill(rho_mass);
                           }else if(Q2 >= 3.5 && Q2 < 4){
                              ivp_hist->Fill(rho_mass);
                     
                           } else if(Q2 >= 4.5 && Q2 < 6){
                              vip_hist->Fill(rho_mass);
                           } 
                          
                        
                     }
                  }
               }
 

            }
 	       }

            
         } 
         counter++;
      //std::cout << "1" << std::endl; 
   }
   
   std::cout << "2" << std::endl; 
     
     // std::cout << "p = " << p << std::endl;
   }
 //}
   
  
   t_hist->GetXaxis()->SetTitle("-t");
   t_hist->Draw();
   t_hist->Write();

   //z_hist->GetXaxis()->SetTitle("Zh");
   //z_hist->Draw();
   //z_hist->Write();

   w_hist->GetXaxis()->SetTitle("W [GeV]");
   w_hist->Draw();
   w_hist->Write();

   q_hist->GetXaxis()->SetTitle("Q^{2} [GeV^{2}/c]");
   q_hist->Draw();
   q_hist->Write();

  /*
   wcut_hist->GetXaxis()->SetTitle("P [GeV]");
   wcut_hist->GetYaxis()->SetTitle("E_{tot} / P");
   wcut_hist->Draw();
   wcut_hist->Write();

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

   


   ivx_hist->GetXaxis()->SetTitle("vx [cm]");

   ivx_hist->Draw();
   ivx_hist->Write();


   
   ivy_hist->GetXaxis()->SetTitle("vy [cm]");
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
   iivx_hist->Draw();
   iivx_hist->Write();

   iivy_hist->GetXaxis()->SetTitle("vy [cm]");
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

   */
   
   p_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}");
   p_hist->Draw();
   p_hist->Write();

   ip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}");
   ip_hist->Draw();
   ip_hist->Write();


   iip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}");
   iip_hist->Draw();
   iip_hist->Write();


   iiip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}");
   iiip_hist->Draw();
   iiip_hist->Write();

   ivp_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}");
   ivp_hist->Draw();
   ivp_hist->Write();


   vip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}");
   vip_hist->Draw();
   vip_hist->Write();


   
   return 0;
}


