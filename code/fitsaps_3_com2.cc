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

void fitsaps_3_com2(){


   TFile input("cusn_pass0v11_zzgg_obnt_11.root", "read"); 
   TFile input2("sn_pass0v11_zzgg_obnt_09.root", "read");    ///fits_ld2.root
   TFile input3("ld_pass0v11_zzgg_obnt_11.root", "read");    ///fits_ld2.root
   TFile input4("cxc_pass0v11_zzgg_obnt_11.root", "read");    ///fits_ld2.root
   TFile input5("cusn_pass0v11_zzgg_ibnt_11.root", "read"); 
   TFile input6("sn_pass0v11_zzgg_ibnt_09.root", "read");    ///fits_ld2.root
   TFile input7("ld_pass0v11_zzgg_ibnt_11.root", "read");    ///fits_ld2.root
   TFile input8("cxc_pass0v11_zzgg_ibnt_11.root", "read");    ///fits_ld2.root
    TFile f("sol_pass0v11_obnt_28_fits.root", "update");
   
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
      
      
      
     
      auto mge = new TMultiGraph;
   	mge->SetTitle("Electron Yield;Run Number; Y_{e} / Q [nC^{-1}]");
   	mge->Add(cu_ob_e_g);
   	mge->Add(sn_ob_e_g);
   	mge->Add(ld_ob_e_g);
   	mge->Add(cxc_ob_e_g);
   	mge->Add(cu_ib_e_g);
   	mge->Add(sn_ib_e_g);
   	mge->Add(ld_ib_e_g);
   	mge->Add(cxc_ib_e_g);
   	//mge->SetMinimum(0.0);
   	//mge->SetMaximum(1.2);
   	//TAxis *axis = mg->GetXaxis();
   	//axis->SetLimits(1,6);
   	mge->Draw("AP");
   	mge->Write();
   	
   	auto mgpp = new TMultiGraph;
   	mgpp->SetTitle("#pi^{+} Yield;Run Number; Y_{#pi^{+}} / Q [nC^{-1}]");
   	mgpp->Add(cu_ob_pip_g);
   	mgpp->Add(sn_ob_pip_g);
   	mgpp->Add(ld_ob_pip_g);
   	mgpp->Add(cxc_ob_pip_g);
   	mgpp->Add(cu_ib_pip_g);
   	mgpp->Add(sn_ib_pip_g);
   	mgpp->Add(ld_ib_pip_g);
   	mgpp->Add(cxc_ib_pip_g);
   	mgpp->Draw("AP");
   	mgpp->Write();
   	
   	auto mgpm = new TMultiGraph;
   	mgpm->SetTitle("#pi^{-} Yield;Run Number; Y_{#pi^{-}} / Q [nC^{-1}]");
   	mgpm->Add(cu_ob_pim_g);
   	mgpm->Add(sn_ob_pim_g);
   	mgpm->Add(ld_ob_pim_g);
   	mgpm->Add(cxc_ob_pim_g);
   	mgpm->Add(cu_ib_pim_g);
   	mgpm->Add(sn_ib_pim_g);
   	mgpm->Add(ld_ib_pim_g);
   	mgpm->Add(cxc_ib_pim_g);
   	mgpm->Draw("AP");
   	mgpm->Write();
   	
   	auto mgr1 = new TMultiGraph;
   	mgr1->SetTitle("#rho^{0} Yield for Q^{2} bin 1;Run Number; Y_{#rho^{0}}");
   	mgr1->Add(cu_ob_r1_g);
   	mgr1->Add(sn_ob_r1_g);
   	mgr1->Add(ld_ob_r1_g);
   	mgr1->Add(cxc_ob_r1_g);
   	mgr1->Add(cu_ib_r1_g);
   	mgr1->Add(sn_ib_r1_g);
   	mgr1->Add(ld_ib_r1_g);
   	mgr1->Add(cxc_ib_r1_g);
   	mgr1->Draw("AP");
   	mgr1->Write();
   	
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
     
}


