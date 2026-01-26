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
  Double_t arg3 = (par[1] - x[0])*(par[1] - x[0]); //((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = (0.150*0.150);
  return par[0]*(arg2/(4*arg3 + arg4)) + par[2] + par[3]*pow(x[0],1) + par[4]*pow(x[0],2) + par[5]*pow(x[0],3); //+ par[7]*pow(x[0],4) + par[8]*pow(x[0],5);
}

//"BreitWinger"

void fitsaps_2(){


   TFile input("ld_bspot_kin5_dnp_5_sigma_tzS.root", "read");    ///fits_ld2.root

   
    TFile f("aaaaafits_ld_ob_bspoot_matt_kin3_6.root", "update");
   

   TH1D *ip_hist = new TH1D("Q^{2} (2-2.5) [GeV^{2}]", "Q^{2} (2-2.5) [GeV^{2}]", 200, 0, 2);
       ip_hist->SetLineColor(kBlue);
       //ip_hist->SetStats(0);
   TH1D *p_hist = new TH1D("Q^{2} (1-2) [GeV^{2}]", "Q^{2} (1-2) [GeV^{2}]", 200, 0, 2);
       p_hist->SetLineColor(kBlue);
       //p_hist->SetStats(0);

   TH1D *iip_hist = new TH1D("Q^{2} (2.5-3) [GeV^{2}]", "Q^{2} (2.5-3) [GeV^{2}]", 200, 0, 2);
	iip_hist->SetLineColor(kBlue);
	//iip_hist->SetStats(0);
   

   TH1D *iiip_hist = new TH1D("Q^{2} (3-3.5) [GeV^{2}]", "Q^{2} (3-3.5) [GeV^{2}]", 200, 0, 2);
	iiip_hist->SetLineColor(kBlue);
	//iiip_hist->SetStats(0);
   

   TH1D *ivp_hist = new TH1D("Q^{2} (3.5-4.5) [GeV^{2}]", "Q^{2} (3.5-6) [GeV^{2}]", 200, 0,2);
	ivp_hist->SetLineColor(kBlue);
	//ivp_hist->SetStats(0);
	//ivp_hist->SetFillStyle(4000);
   
   TH1D *vp_hist = new TH1D("Q^{2} (4.5-6) [GeV^{2}]", "Q^{2} (4.5-6) [GeV^{2}]", 200, 0, 2);
	vp_hist->SetLineColor(kBlue);
	//vp_hist->SetStats(0);
   /*
   TH1D *bivp_hist = new TH1D("bbbQ2 (3.5 4.5)", "Q2 (3.5 4.5)", 200, 0.2, 1.4);
	bivp_hist->SetLineColor(kGreen-1);
   
   TH1D *vip_hist = new TH1D("bbbQ2 (4.5 6)", "Q2 (4.5 6)", 200, 0.2, 1.4);
        vip_hist->SetLineColor(kGreen-1);
     */   
       // TH1D *vizero_hist = new TH1D("","",200,,);
       // TFile b("zrm_rg_b_18451_4_15_cxc.root", "update");

        TH1D *tp_hist = (TH1D*)input.Get("sector 0 P");  //"Q^{2} (1 2) [GeV^{2}]
      TH1D *tip_hist = (TH1D*)input.Get("sector 1 P");
      TH1D *tiip_hist = (TH1D*)input.Get("sector 2 P");
      TH1D *tiiip_hist = (TH1D*)input.Get("sector 3 P");
      TH1D *tivp_hist = (TH1D*)input.Get("sector 4 P");
      TH1D *tvp_hist = (TH1D*)input.Get("sector 5 P");
      //TH1D *tvip_hist = (TH1D*)input.Get("sector 0 P");
      p_hist->Add(tp_hist);
      ip_hist->Add(tip_hist);
      iip_hist->Add(tiip_hist);
      iiip_hist->Add(tiiip_hist);
      ivp_hist->Add(tivp_hist);
      vp_hist->Add(tvp_hist);
      //ivp_hist->Add(tvip_hist);
     
  


   
  

   p_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //p_hist->Scale(1.0 / p_hist->Integral());
   TF1 *p_func = new TF1("p_fit",fitfbw,0.3,1.4,6);
   //p_func->SetParameter(1,0.150);
   p_func->SetParameter(1,0.77);
   p_func->SetParNames("double pConstant","double pgamma","double pmean","double pa1","double pa2","double pa3", "double pa4", "double pa5", "double pa6");
   p_hist->Fit("p_fit","","",0.3,1.4);
   p_hist->Draw();
   p_hist->Write();

   ip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //ip_hist->Scale(1.0 / ip_hist->Integral());
   TF1 *ip_func = new TF1("ip_fit",fitfbw,0.3,1.4,6);
   //ip_func->SetParameter(1,0.150);
   ip_func->SetParameter(1,0.77);
   ip_func->SetParNames("double ipConstant","double ipgamma","double ipmean","double ipa1","double ipa2","double ipa3", "double ipa4", "double ipa5", "double ipa6");
   ip_hist->Fit("ip_fit","","",0.3,1.4);
   ip_hist->Draw();
   ip_hist->Write();


   iip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //iip_hist->Scale(1.0 / iip_hist->Integral());
   TF1 *iip_func = new TF1("iip_fit",fitfbw,0.3,1.4,6);
   //iip_func->SetParameter(1,0.150);
   iip_func->SetParameter(1,0.77);
   iip_func->SetParNames("double iipConstant","double iipgamma","double iipmean","double iipa1","double iipa2","double iipa3", "double iipa4", "double iipa5", "double iipa6");
   iip_hist->Fit("iip_fit","","",0.3,1.4);
   iip_hist->Draw();
   iip_hist->Write();




 
   iiip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //iiip_hist->Scale(1.0 / iiip_hist->Integral());
   TF1 *iiip_func = new TF1("iiip_fit",fitfbw,0.3,1.4,6);
   //iiip_func->SetParameter(1,0.150);
   iiip_func->SetParameter(1,0.77);
   iiip_func->SetParNames("double iiipConstant","double iiipgamma","double iiipmean","double iiipa1","double iiipa2","double iiipa3", "double iiipa4", "double iiipa5", "double iiipa6");
   iiip_hist->Fit("iiip_fit","","",0.3,1.4);
   iiip_hist->Draw();
   iiip_hist->Write();



   ivp_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //ivp_hist->Scale(1.0 / ivp_hist->Integral());
   TF1 *ivp_func = new TF1("ivp_fit",fitfbw,0.3,1.4,6);
   //ivp_func->SetParameter(1,0.150);
   ivp_func->SetParameter(1,0.77);
   ivp_func->SetParNames("double ivpConstant","double ivpgamma","double ivpmean","double ivpa1","double ivpa2","double ivpa3", "double ivpa4", "double ivpa5", "double ivpa6");
   ivp_hist->Fit("ivp_fit","","",0.3,1.4);
   ivp_hist->Draw();
   ivp_hist->Write();



   vp_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //vp_hist->Scale(1.0 / vp_hist->Integral());
   TF1 *vp_func = new TF1("vp_fit",fitfbw,0.3,1.4,6);
   //vp_func->SetParameter(1,0.150);
   vp_func->SetParameter(1,0.77);
   vp_func->SetParNames("double vpConstant","double vpgamma","double vpmean","double vpa1","double vpa2","double vpa3", "double vpa4", "double vpa5", "double vpa6");
   vp_hist->Fit("vp_fit","","",0.3,1.4);
   vp_hist->Draw();
   vp_hist->Write();
   
  /* 
   bivp_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   //bivp_hist->Scale(1.0 / bivp_hist->Integral());
   TF1 *bivp_func = new TF1("bivp_fit",fitfbw,0.4,1.4,7);
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
   TF1 *vip_func = new TF1("vip_fit",fitfbw,0.4,1.4,7);
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
   
  */ 
 


}


