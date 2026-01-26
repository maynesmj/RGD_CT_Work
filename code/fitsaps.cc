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

void fitsaps(){


   TFile input("ld_ob_bspot_kin1.root", "read");    ///fits_ld2.root

   
    TFile f("aaaafits_ld_ob_bspot_for_CLAS_5.root", "update");
/*
        std::cout << (0.82 * 1966.09) / (0.88 * 2604.77) << " " << (0.82 * 1966.09) / (0.88 * 2497.7) << std::endl;
	std::cout << (0.82 * 266.725) / (0.88 * 325.119) << " " << (0.82 * 266.725) / (0.88 * 317.708) << std::endl;
	std::cout << (0.82 * 1966.09) / (0.88 * 2604.77) << " " << (0.82 * 1966.09) / (0.88 * 2497.7) << std::endl;
	std::cout << "" << std::endl;
	std::cout << (0.88 * 1966.09) / (0.82 * 2604.77) << " " << (0.88 * 1966.09) / (0.82 * 2497.7) << std::endl;
	std::cout << (0.88 * 266.725) / (0.82 * 325.119) << " " << (0.88 * 266.725) / (0.82 * 317.708) << std::endl;
	std::cout << (0.88 * 1966.09) / (0.82 * 2604.77) << " " << (0.88 * 1966.09) / (0.82 * 2497.7) << std::endl;
*/    

   TH1D *ip_hist = new TH1D("Q^{2} (2-2.5) [GeV^{2}]", "Q^{2} (2-2.5) [GeV^{2}]", 200, 0, 2);
       ip_hist->SetFillColor(kBlue);
       ip_hist->SetStats(0);
   TH1D *p_hist = new TH1D("Q^{2} (1-2) [GeV^{2}]", "Q^{2} (1-2) [GeV^{2}]", 200, 0, 2);
       p_hist->SetFillColor(kBlue);
       p_hist->SetStats(0);

   TH1D *iip_hist = new TH1D("Q^{2} (2.5-3) [GeV^{2}]", "Q^{2} (2.5-3) [GeV^{2}]", 200, 0, 2);
	iip_hist->SetFillColor(kBlue);
	iip_hist->SetStats(0);
   

   TH1D *iiip_hist = new TH1D("Q^{2} (3-3.5) [GeV^{2}]", "Q^{2} (3-3.5) [GeV^{2}]", 200, 0, 2);
	iiip_hist->SetFillColor(kBlue);
	iiip_hist->SetStats(0);
   

   TH1D *ivp_hist = new TH1D("Q^{2} (3.5-4.5) [GeV^{2}]", "Q^{2} (3.5-4.5) [GeV^{2}]", 200, 0,2);
	ivp_hist->SetFillColor(kBlue);
	ivp_hist->SetStats(0);
	//ivp_hist->SetFillStyle(4000);
   
   TH1D *vp_hist = new TH1D("Q^{2} (4.5-6) [GeV^{2}]", "Q^{2} (4.5-6) [GeV^{2}]", 200, 0, 2);
	vp_hist->SetFillColor(kBlue);
	vp_hist->SetStats(0);
   
   TH1D *bivp_hist = new TH1D("bbbQ2 (3.5 4.5)", "Q2 (3.5 4.5)", 200, 0.2, 1.4);
	bivp_hist->SetFillColor(kBlue-1);
   
   TH1D *vip_hist = new TH1D("bbbQ2 (4.5 6)", "Q2 (4.5 6)", 200, 0.2, 1.4);
        vip_hist->SetFillColor(kBlue-1);
        
       // TH1D *vizero_hist = new TH1D("","",200,,);
       // TFile b("zrm_rg_b_18451_4_15_cxc.root", "update");

        TH1D *tp_hist = (TH1D*)input.Get("sector 0 P");   //Q^{2} (1 2) [GeV^{2}]
      TH1D *tip_hist = (TH1D*)input.Get("sector 1 P");
      TH1D *tiip_hist = (TH1D*)input.Get("sector 2 P");
      TH1D *tiiip_hist = (TH1D*)input.Get("sector 3 P");
      TH1D *tivp_hist = (TH1D*)input.Get("sector 4 P");
      TH1D *tvp_hist = (TH1D*)input.Get("sector 5 P");
      TH1D *tvip_hist = (TH1D*)input.Get("sector 6 P");
      p_hist->Add(tp_hist);
      ip_hist->Add(tip_hist);
      iip_hist->Add(tiip_hist);
      iiip_hist->Add(tiiip_hist);
      ivp_hist->Add(tivp_hist);
      vp_hist->Add(tvp_hist);
      //vp_hist->Add(tvip_hist);
     //bivp_hist->Add(vp_hist);
     //vip_hist->Add(ivp_hist);
     //vip_hist->Add(vizero_hist);
     //TFile input("zzzzrm_rg_b__4_16_cxc.root", "read");
     
     TH1D *q_hist = new TH1D("Q^{2} [GeV^{2}]", "", 200, 0, 10);
       q_hist->SetStats(0);
   TH1D *w_hist = new TH1D("W [GeV^{2}]", "", 200, 0, 5);
       w_hist->SetStats(0);

   TH1D *t_hist = new TH1D("t [GeV^{2}]", "", 200, 0, 20);
	t_hist->SetStats(0);
   

   TH1D *z_hist = new TH1D("z", "", 200, 0, 10);
	z_hist->SetStats(0);
   

   TH1D *wo_hist = new TH1D("Without Cuts", "Without Cuts", 200, 0.2,1.4);
	wo_hist->SetStats(0);
   
   TH1D *ww_hist = new TH1D("With W Cut", "With W Cut", 200, 0.2, 1.4);
	ww_hist->SetStats(0);
   
   TH1D *wwt_hist = new TH1D("With W & -t Cuts", "With W & -t Cuts", 200, 0.2, 1.4);
	wwt_hist->SetStats(0);
   
   TH1D *wwtz_hist = new TH1D("With W, -t, & z Cuts", "With W, -t, & z Cuts", 200, 0.2, 1.4);
        wwtz_hist->SetStats(0);
        

   TH1D *ivx_hist = (TH1D*)input.Get("lc");
   TH1D *ivy_hist = (TH1D*)input.Get("w (Invariant Mass)");
   
   TH1D *iivx_hist = (TH1D*)input.Get("-t");
   TH1D *iivy_hist = (TH1D*)input.Get("Zh");
   
   TH1D *iiivx_hist = (TH1D*)input.Get("Without Cuts");
   TH1D *iiivy_hist = (TH1D*)input.Get("With W Cut");
   
   TH1D *ivvx_hist = (TH1D*)input.Get("With W Cut & -t Cut");
   TH1D *ivvy_hist = (TH1D*)input.Get("With W Cut & -t Cut & z Cut");
    q_hist->Add(ivx_hist);
      w_hist->Add(ivy_hist);
      t_hist->Add(iivx_hist);
      z_hist->Add(iivy_hist);
      wo_hist->Add(iiivx_hist);
      ww_hist->Add(iiivy_hist);
      wwt_hist->Add(ivvx_hist);
      wwtz_hist->Add(ivvy_hist);
   
   

   //TH2F *z_hist = new TH2F("px vs py", "px vs py", 200, -2, 2, 200, -2, 2);

   //auto cgr = new TGraph();
   //std::cout << "1 th" << std::endl;

   std::vector<std::string> v3_a;
   int num = 0;





   
  

   p_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //p_hist->Scale(1.0 / p_hist->Integral());
   TF1 *p_func = new TF1("p_fit",fitfbw,0.4,1.4,7);
   //p_func->SetParameters(0,0,0, 1,-142857, 1, 70);  //500,p_hist->GetMean(),p_hist->GetRMS(), -142857, 0, 70);
   //p_func->SetParameter(0,7000);
   p_func->SetParameter(1,0.150);
   p_func->SetParameter(2,0.77);
   //p_func->SetParLimits(2,0.65,0.89);
   //p_func->SetParameter(8,0);
   //p_func->SetParLimits(3,-20000, 20000);
   //p_func->SetParLimits(4,-2000000, 0);
   //p_func->SetParLimits(5,0, 2000000);
   //p_func->SetParLimits(6,-20000, 20000);
   p_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   p_hist->Fit("p_fit","","",0.4,1.4);
   p_hist->Draw();
   p_hist->Write();

   ip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //ip_hist->Scale(1.0 / ip_hist->Integral());
   TF1 *ip_func = new TF1("ip_fit",fitfbw,0.4,1.4,7);
   //ip_func->SetParameters(500,0.77,(150.0/2), 1,-142857, 1, 70);
   //ip_func->SetParLimits(0,-20000, 20000);
   //ip_func->SetParLimits(0,480, 480);
   ip_func->SetParameter(1,0.150);
   ip_func->SetParameter(2,0.77);
   //ip_func->SetParLimits(2,0.65,0.89);
   //ip_func->SetParLimits(1,-20000, 20000);
   //ip_func->SetParLimits(2,-20000, 20000);
   //ip_func->SetParLimits(3,-20000, 20000);
   //ip_func->SetParLimits(4,-2000000,0);
   //ip_func->SetParLimits(5,0,2000000);
   //ip_func->SetParLimits(6,-20000,20000);
   ip_func->SetParNames("double ipConstant","double ipgamma","double ipmean","double ipa1","double ipa2","double ipa3", "double ipa4", "double ipa5", "double ipa6");
   ip_hist->Fit("ip_fit","","",0.4,1.4);
   ip_hist->Draw();
   ip_hist->Write();


   iip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //iip_hist->Scale(1.0 / iip_hist->Integral());
   TF1 *iip_func = new TF1("iip_fit",fitfbw,0.4,1.4,7);
   //iip_func->SetParameters(500,0.77,(150.0/2), 0,0, 0, 0);
   //iip_func->SetParLimits(0,-20000, 20000);
   iip_func->SetParameter(1,0.150);
   iip_func->SetParameter(2,0.77);
   //iip_func->SetParLimits(2,0.65,0.89);
   //iip_func->SetParLimits(1,-20000, 20000);
   //iip_func->SetParLimits(2,-20000, 20000);
   //iip_func->SetParLimits(3,-20000, 20000);
   //iip_func->SetParLimits(4,-200000,20000);
   //iip_func->SetParLimits(5,-20000,20000);
   //iip_func->SetParLimits(6,-20000,20000);
   iip_func->SetParNames("double iipConstant","double iipgamma","double iipmean","double iipa1","double iipa2","double iipa3", "double iipa4", "double iipa5", "double iipa6");
   iip_hist->Fit("iip_fit","","",0.4,1.4);
   iip_hist->Draw();
   iip_hist->Write();




 
   iiip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //iiip_hist->Scale(1.0 / iiip_hist->Integral());
   TF1 *iiip_func = new TF1("iiip_fit",fitfbw,0.4,1.4,7);
   //iiip_func->SetParameters(500,0.77,(150.0/2), 0,-18000, 0, 45);     //iiip_hist->GetMean(),iiip_hist->GetRMS(
   //iiip_func->SetParLimits(0,-20000, 20000);
   //iiip_func->SetParLimits(1,-20000, 20000);
   //iiip_func->SetParLimits(2,-20000, 20000);
   //iiip_func->SetParLimits(3,-20000, 20000);
   //iiip_func->SetParLimits(4,-2000000,0);
   //iiip_func->SetParLimits(5,0,2000000);
   //iiip_func->SetParLimits(6,-20000,20000);
   iiip_func->SetParameter(1,0.150);
   iiip_func->SetParameter(2,0.77);
   //iiip_func->SetParLimits(2,0.65,0.89);
   iiip_func->SetParNames("double iiipConstant","double iiipgamma","double iiipmean","double iiipa1","double iiipa2","double iiipa3", "double iiipa4", "double iiipa5", "double iiipa6");
   iiip_hist->Fit("iiip_fit","","",0.4,1.4);
   iiip_hist->Draw();
   iiip_hist->Write();



   ivp_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //ivp_hist->Scale(1.0 / ivp_hist->Integral());
   TF1 *ivp_func = new TF1("ivp_fit",fitfbw,0.4,1.4,7);
   //ivp_func->SetParameters(500,0.77,(150.0/2), 0,-18000, 0, 45);
   //ivp_func->SetParLimits(0,-20000, 20000);
   //ivp_func->SetParLimits(1,-20000, 20000);
   //ivp_func->SetParLimits(2,-20000, 20000);
   //ivp_func->SetParLimits(3,-20000, 20000);
   //ivp_func->SetParLimits(4,-200000,20000);
   //ivp_func->SetParLimits(5,-20000,20000);
   //ivp_func->SetParLimits(6,-20000,20000);
   ivp_func->SetParameter(1,0.150);
   ivp_func->SetParameter(2,0.77);
   //ivp_func->SetParLimits(2,0.65,0.89);
   ivp_func->SetParNames("double ivpConstant","double ivpgamma","double ivpmean","double ivpa1","double ivpa2","double ivpa3", "double ivpa4", "double ivpa5", "double ivpa6");
   ivp_hist->Fit("ivp_fit","","",0.4,1.4);
   ivp_hist->Draw();
   ivp_hist->Write();



   vp_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //vp_hist->Scale(1.0 / vp_hist->Integral());
   TF1 *vp_func = new TF1("vp_fit",fitfbw,0.4,1.4,7);
   //vp_func->SetParameters(500,0.77,(150.0/2), 0,-18000, 0, 45);
   //vp_func->SetParLimits(0,-20000, 20000);
   //vp_func->SetParLimits(1,-20000, 20000);
   //vp_func->SetParLimits(2,-20000, 20000);
   //vp_func->SetParLimits(3,-20000, 20000);
   //vp_func->SetParLimits(4,-200000,20000);
   //vp_func->SetParLimits(5,-20000,20000);
   //vp_func->SetParLimits(6,-20000,20000);
   vp_func->SetParameter(1,0.150);
   vp_func->SetParameter(2,0.77);
   //vp_func->SetParLimits(2,0.65,0.89);
   vp_func->SetParNames("double vpConstant","double vpgamma","double vpmean","double vpa1","double vpa2","double vpa3", "double vpa4", "double vpa5", "double vpa6");
   vp_hist->Fit("vp_fit","","",0.4,1.4);
   vp_hist->Draw();
   vp_hist->Write();
   
   
   bivp_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   //bivp_hist->Scale(1.0 / bivp_hist->Integral());
   TF1 *bivp_func = new TF1("bivp_fit",fitfbw,0,2,9);
   //ivp_func->SetParameters(500,0.77,(150.0/2), 0,-18000, 0, 45);
   //ivp_func->SetParLimits(0,-20000, 20000);
   //ivp_func->SetParLimits(1,-20000, 20000);
   //ivp_func->SetParLimits(2,-20000, 20000);
   //ivp_func->SetParLimits(3,-20000, 20000);
   //ivp_func->SetParLimits(4,-200000,20000);
   //ivp_func->SetParLimits(5,-20000,20000);
   //ivp_func->SetParLimits(6,-20000,20000);
   bivp_func->SetParameter(1,0.150);
   bivp_func->SetParameter(2,0.77);
   //bivp_func->SetParLimits(2,0.65,0.89);
   bivp_func->SetParNames("double bivpConstant","double bivpgamma","double bivpmean","double bivpa1","double bivpa2","double bivpa3", "double bivpa4", "double bivpa5", "double bivpa6");
   bivp_hist->Fit("bivp_fit","","",0.3,1.2);
   bivp_hist->Draw();
   bivp_hist->Write();



   vip_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   //vip_hist->Scale(1.0 / vip_hist->Integral());
   TF1 *vip_func = new TF1("vip_fit",fitfbw,0,2,9);
   //vip_func->SetParameters(500,0.77,(150.0/2), 0,-18000, 0, 45);
   //vip_func->SetParLimits(0,-20000, 20000);
   //vip_func->SetParLimits(1,-20000, 20000);
   //vip_func->SetParLimits(2,-20000, 20000);
   //vip_func->SetParLimits(3,-20000, 20000);
   //vip_func->SetParLimits(4,-200000,20000);
   //vip_func->SetParLimits(5,-20000,20000);
   //vip_func->SetParLimits(6,-20000,20000);
   vip_func->SetParameter(1,0.150);
   vip_func->SetParameter(2,0.77);
   //vip_func->SetParLimits(2,0.65,0.89);
   vip_func->SetParNames("double vipConstant","double vipgamma","double vipmean","double vipa1","double vipa2","double vipa3", "double vipa4", "double vipa5", "double vipa6");
   vip_hist->Fit("vip_fit","","",0.3,1.2);
   vip_hist->Draw();
   vip_hist->Write();
   
   q_hist->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
   q_hist->Draw();
   q_hist->Write();
   
   w_hist->GetXaxis()->SetTitle("W [GeV]");
   w_hist->Draw();
   w_hist->Write();
 
 
 t_hist->GetXaxis()->SetTitle("-t [GeV^{2}]");
   t_hist->Draw();
   t_hist->Write();
 
 
 //z_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   z_hist->Draw();
   z_hist->Write();
 
 
 wo_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   wo_hist->Draw();
   wo_hist->Write();
 
 
 ww_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
  ww_hist->Draw();
   ww_hist->Write();
 
 
 wwt_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   wwt_hist->Draw();
   wwt_hist->Write();
 
 
 wwtz_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   wwtz_hist->Draw();
   wwtz_hist->Write();
 
 
/*  
    ivx_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   TF1 *ivx_func = new TF1("ivx_fit",fitf,-1,1,6);
   //p_func->SetParameters(0,0,0, 1,-142857, 1, 70);  //500,p_hist->GetMean(),p_hist->GetRMS(), -142857, 0, 70);
   //p_func->SetParameter(2,0.77);
   //p_func->SetParLimits(2,0.65,0.89);
   ivx_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3");
   ivx_hist->Fit("ivx_fit","","",-0.4,0.4);
   ivx_hist->Draw();
   ivx_hist->Write();
   
   ivy_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   TF1 *ivy_func = new TF1("ivy_fit",fitf,-1,1,6);
   //p_func->SetParameters(0,0,0, 1,-142857, 1, 70);  //500,p_hist->GetMean(),p_hist->GetRMS(), -142857, 0, 70);
   //p_func->SetParameter(2,0.77);
   //p_func->SetParLimits(2,0.65,0.89);
   ivy_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3");
   ivy_hist->Fit("ivy_fit","","",0.2,0.8);
   ivy_hist->Draw();
   ivy_hist->Write();

iivx_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   TF1 *iivx_func = new TF1("iivx_fit",fitfy,-1,1,6);
   //p_func->SetParameters(00,0, 1,-142857, 1, 70);  //500,p_hist->GetMean(),p_hist->GetRMS(), -142857, 0, 70);
   //p_func->SetParameter(2,0.77);
   //p_func->SetParLimits(2,0.65,0.89);
   iivx_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3");
   iivx_hist->Fit("iivx_fit","","",-0.4,0.4);
   iivx_hist->Draw();
   iivx_hist->Write();

*/

}


