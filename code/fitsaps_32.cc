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

void fitsaps_32(){


   TFile input("cusn_kin5_pass0v5_dnp_300_no_sigma_tzS1_fits.root", "read");    ///fits_ld2.root    cusn_kin5_dnp_4_sigma_tzS.root  cusn_kin5_pass0v5_dnp_300_no_sigma_tzS.root
   TFile input2("cusn_kin5_pass0v4_dnp_300_no_sigma_tzS1_fits.root", "read");
    TFile f("cusn_kin5_pass0v4&5_dnp_300_no_sigma_tzS1_fits.root", "update");
   gStyle->SetOptStat("");
	gStyle->SetCanvasPreferGL();
	
   TH1F *vx_hist = new TH1F("vx", "total vx", 300, -3, 3);
       //vx_hist->SetStats(0);
       vx_hist->SetLineColor(kBlack);
       vx_hist->GetXaxis()->SetTitle("vx [cm]");
   TH1D *vy_hist = new TH1D("vy", "total vy", 300, -3, 3);
       //vy_hist->SetStats(0);
       vy_hist->SetLineColor(kBlack);
       vy_hist->GetXaxis()->SetTitle("vy [cm]");
   TH1D *vz_hist = new TH1D("vz", "total vz", 300, -20, 10);
       //vz_hist->SetStats(0);
       vz_hist->SetLineColor(kBlack);
       vz_hist->GetXaxis()->SetTitle("vz [cm]");

   TH1F *vx1_hist = new TH1F("Sector 1 vx", "Pass0v5 vx", 300, -3, 3);
      // vx1_hist->SetStats(0);
       vx1_hist->SetLineColor(kCyan+1);
   TH1D *vy1_hist = new TH1D("Sector 1 vy", "Pass0v5 vy", 300, -3, 3);
     //  vy1_hist->SetStats(0);
       vy1_hist->SetLineColor(kCyan+1);
   TH1D *vz1_hist = new TH1D("Sector 1 vz", "Pass0v5 vz", 300, -20, 10);
    //   vz1_hist->SetStats(0);
       vz1_hist->SetLineColor(kCyan+1);
       
  TH1F *vx2_hist = new TH1F("Sector 2 vx", "Pass0v4 vx", 300, -3, 3);
   //    vx2_hist->SetStats(0);
       vx2_hist->SetLineColor(kOrange+1);
   TH1D *vy2_hist = new TH1D("Sector 2 vy", "Pass0v4 vy", 300, -3, 3);
    //   vy2_hist->SetStats(0);
       vy2_hist->SetLineColor(kOrange+1);
   TH1D *vz2_hist = new TH1D("Sector 2 vz", "Pass0v4 vz", 300, -20, 10);
    //   vz2_hist->SetStats(0);
       vz2_hist->SetLineColor(kOrange+1);
             
  
  

      TH1F *tvx1_hist = (TH1F*)input.Get("vx"); 
      TH1D *tvy1_hist = (TH1D*)input.Get("vy");
      TH1D *tvz1_hist = (TH1D*)input.Get("vz");
      
      TH1F *tvx2_hist = (TH1F*)input2.Get("vx"); 
      TH1D *tvy2_hist = (TH1D*)input2.Get("vy");
      TH1D *tvz2_hist = (TH1D*)input2.Get("vz");
     
      
      vx_hist->Add(tvx1_hist);
      vx_hist->Add(tvx2_hist);
      vx_hist->Scale(1.0 / vx_hist->Integral());
      vx_hist->SetMinimum(0.0);
      vx_hist->SetMaximum(0.08);
      vy_hist->Add(tvy1_hist);
      vy_hist->Add(tvy2_hist);
      vy_hist->Scale(1.0 / vy_hist->Integral());
      vy_hist->SetMinimum(0.0);
      vy_hist->SetMaximum(0.04);
      vz_hist->Add(tvz1_hist);
      vz_hist->Add(tvz2_hist);
      vz_hist->Scale(1.0 / vz_hist->Integral());
      vz_hist->SetMinimum(0.0);
      vz_hist->SetMaximum(0.04);
      
      
      
      vx1_hist->Add(tvx1_hist);
      vx1_hist->Scale(1.0 / vx1_hist->Integral());
      vx1_hist->SetMinimum(0.0);
      vx1_hist->SetMaximum(0.08);
      vy1_hist->Add(tvy1_hist);
      vy1_hist->Scale(1.0 / vy1_hist->Integral());
      vy1_hist->SetMinimum(0.0);
      vy1_hist->SetMaximum(0.04);
      vz1_hist->Add(tvz1_hist);
      vz1_hist->Scale(1.0 / vz1_hist->Integral());
      vz1_hist->SetMinimum(0.0);
      vz1_hist->SetMaximum(0.04);
      
      vx2_hist->Add(tvx2_hist);
      vx2_hist->Scale(1.0 / vx2_hist->Integral());
      vx2_hist->SetMinimum(0.0);
      vx2_hist->SetMaximum(0.08);
      vy2_hist->Add(tvy2_hist);
      vy2_hist->Scale(1.0 / vy2_hist->Integral());
      vy2_hist->SetMinimum(0.0);
      vy2_hist->SetMaximum(0.04);
      vz2_hist->Add(tvz2_hist);
      vz2_hist->Scale(1.0 / vz2_hist->Integral());
      vz2_hist->SetMinimum(0.0);
      vz2_hist->SetMaximum(0.04);
     
      //float mul = 300000000000000000000000000000000000000000.0; // 3.0 4.0 5.0
                 float ilxlim = 0.00656862 - 3*(0.3567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float iuxlim = 0.00656862 + 3*(0.3567);  //1.5741; // 0.9525 1.2633 1.5741
                 float ilylim = -0.1793 - 3*(0.33594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float iuylim = -0.1793 + 3*(0.33594); 
                 float iilxlim = 0.00656862 - 4*(0.3567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float iiuxlim = 0.00656862 + 4*(0.3567);  //1.5741; // 0.9525 1.2633 1.5741
                 float iilylim = -0.1793 - 4*(0.33594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float iiuylim = -0.1793 + 4*(0.33594); 
                 float iiilxlim = 0.00656862 - 5*(0.3567);   //-1.5339; // -0.9123 -1.2231 -1.5339 
                 float iiiuxlim = 0.00656862 + 5*(0.3567);  //1.5741; // 0.9525 1.2633 1.5741
                 float iiilylim = -0.1793 - 5*(0.33594);  //-1.8538; // -1.18192 -1.51786 -1.8538
                 float iiiuylim = -0.1793 + 5*(0.33594); 
      //gStyle->SetCanvasPreferGL();
      TCanvas *c1 = new TCanvas("c1","vx",0,0,640,480);
      std::cout << "max height " << c1->GetUymax() << std::endl;
      TLine *ilx=new TLine(ilxlim,c1->GetUymin(),ilxlim,c1->GetUymax()-0.2);
      ilx->SetLineColor(kBlue);
      TLine *iux=new TLine(iuxlim,c1->GetUymin(),iuxlim,c1->GetUymax()-0.2);
      iux->SetLineColor(kBlue);
      TLine *iilx=new TLine(iilxlim,c1->GetUymin(),iilxlim,c1->GetUymax()-0.2);
      iilx->SetLineColor(kGreen);
      TLine *iiux=new TLine(iiuxlim,c1->GetUymin(),iiuxlim,c1->GetUymax()-0.2);
      iiux->SetLineColor(kGreen);
      TLine *iiilx=new TLine(iiilxlim,c1->GetUymin(),iiilxlim,c1->GetUymax()-0.2);
      iiilx->SetLineColor(kMagenta);
      TLine *iiiux=new TLine(iiiuxlim,c1->GetUymin(),iiiuxlim,c1->GetUymax()-0.2);
      iiiux->SetLineColor(kMagenta);
      
      vx_hist->Draw("HIST");
      vx1_hist->Draw("HIST Same");
      vx2_hist->Draw("HIST Same");
      c1->BuildLegend();
      ilx->Draw("HIST Same");
      iux->Draw("HIST Same");
      iilx->Draw("HIST Same");
      iiux->Draw("HIST Same");
      iiilx->Draw("HIST Same");
      iiiux->Draw("HIST Same");
      c1->Write();
      
      TCanvas *c2 = new TCanvas("c2","vy",0,0,640,480);
      std::cout << "max height " << c2->GetUymax() << std::endl;
      TLine *ily=new TLine(ilylim,c2->GetUymin(),ilylim,c2->GetUymax()-0.2);
      ily->SetLineColor(kBlue);
      TLine *iuy=new TLine(iuylim,c2->GetUymin(),iuylim,c2->GetUymax()-0.2);
      iuy->SetLineColor(kBlue);
      TLine *iily=new TLine(iilylim,c2->GetUymin(),iilylim,c2->GetUymax()-0.3);
      iily->SetLineColor(kGreen);
      TLine *iiuy=new TLine(iiuylim,c2->GetUymin(),iiuylim,c2->GetUymax()-0.3);
      iiuy->SetLineColor(kGreen);
      TLine *iiily=new TLine(iiilylim,c2->GetUymin(),iiilylim,c2->GetUymax()-0.4);
      iiily->SetLineColor(kMagenta);
      TLine *iiiuy=new TLine(iiiuylim,c2->GetUymin(),iiiuylim,c2->GetUymax()-0.4);
      iiiuy->SetLineColor(kMagenta);
      vy_hist->Draw("HIST");
      vy1_hist->Draw("HIST Same");
      vy2_hist->Draw("HIST Same");
      c2->BuildLegend();
      ily->Draw("HIST Same");
      iuy->Draw("HIST Same");
      iily->Draw("HIST Same");
      iiuy->Draw("HIST Same");
      iiily->Draw("HIST Same");
      iiiuy->Draw("HIST Same");
      c2->Write();
      
      
      TCanvas *c3 = new TCanvas("c3","vz",0,0,640,480);
      vz_hist->Draw("HIST");
      vz1_hist->Draw("HIST Same");
      vz2_hist->Draw("HIST Same");
      c3->BuildLegend();
      c3->Write();
      
      /*
      TCanvas *c5 = new TCanvas("c5","vx",0,0,640,480);
      vx_hist->Draw("HIST");
      //vx_hist->Write();
      c5->Write();
      
      TCanvas *c6 = new TCanvas("c6","vy",0,0,640,480);
      vy_hist->Draw("HIST");
      //vy_hist->Write();
      c6->Write();
      
      TCanvas *c7 = new TCanvas("c7","vz",0,0,640,480);
      vz_hist->Draw("HIST");
      //vz_hist->Write();
      c7->Write();
      */
      vx_hist->SetMinimum(0.0001);
      vx_hist->Write();
      vy_hist->SetMinimum(0.0001);
      vy_hist->Write();
      vz_hist->SetMinimum(0.0001);
      vz_hist->Write();
     
  


   
  
 


}


