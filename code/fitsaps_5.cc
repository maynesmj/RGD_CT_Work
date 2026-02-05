//******************************************************************
//*  ██╗  ██╗██╗██████╗  ██████╗     ██╗  ██╗    ██████╗
//*  ██║  ██║██║██╔══██╗██╔═══██╗    ██║  ██║   ██╔═████╗
//*  ███████║██║██████╔╝██║   ██║    ███████║   ██║██╔██║
//*  ██╔══██║██║██╔═══╝ ██║   ██║    ╚════██║   ████╔╝██║
//*  ██║  ██║██║██║     ╚██████╔╝         ██║██╗╚██████╔╝
//*  ╚═╝  ╚═╝╚═╝╚═╝      ╚═════╝          ╚═╝╚═╝ ╚═════╝
//************************ Jefferson National Lab (2017) ***********
//******************************************************************


/*

This code applies a Briet-Winger + 3D poly fit to the invarent mass plots of the root file. The code will produce a separet root file of the invarient mass plots with the 
fits. This is what I used for presentations. It will also print out the values of the yield based on each fit.

*/

#include <cstdlib>
#include <iostream>
#include <string>
#include <TMath.h>



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
  //Double_t arg1 = 2.0/4.0; // 2 over pi
  Double_t arg2 = 0.150/2; //Gamma=par[1]  M=par[2]
  Double_t arg3 = (par[1] - x[0])*(par[1] - x[0]); //((x[0]*x[0]) - (par[2]*par[2]))*((x[0]*x[0]) - (par[2]*par[2]));
  Double_t arg4 = (arg2*arg2);
  return par[0]*(1.0/(2* TMath::Pi()))*(arg2/(arg3 + arg4)) + par[2] + par[3]*pow(x[0],1) + par[4]*pow(x[0],2) + par[5]*pow(x[0],3);// + par[6]*pow(x[0],4) + par[7]*pow(x[0],5);
}

//"BreitWinger"

double funcy(double *xf, double *par){
         double bwy = par[0]*(1.0/(2* TMath::Pi()))*((0.15/2)/(((par[1] - xf[0])*(par[1] - xf[0])) + ((0.15/2)*(0.15/2))));     
     return bwy;
     }
     
double funcp(double *xp, double *par){
         double ypol = par[0] + par[1]*pow(xp[0],1) + par[2]*pow(xp[0],2) + par[3]*pow(xp[0],3);//  + par[4]*pow(xp[0],4)  + par[5]*pow(xp[0],5);    
     return ypol;
     }


void fitsaps_5(){


   TFile input("bcolloquial_pass1_ldob_02_r1.root", "read");    ///ld_kin5_dnp_5_sigma_tzS.root  ld_calvp_no_sigma_tzS02.root

   
    TFile f("Invarient_mass_plots_for_Mathieu.root", "update");
   

   TH1D *ip_hist = new TH1D("Q^{2} (2-2.5) [GeV^{2}]", "Q^{2} (2-2.5) [GeV^{2}]", 200, 0.2, 2);
       ip_hist->SetFillColor(kBlue+2);
       ip_hist->SetStats(0);
   TH1D *p_hist = new TH1D("Q^{2} (1-2) [GeV^{2}]", "Q^{2} (1-2) [GeV^{2}]", 200, 0.2, 2);
       p_hist->SetFillColor(kBlue+2);
       p_hist->SetStats(0);

   TH1D *iip_hist = new TH1D("Q^{2} (2.5-3) [GeV^{2}]", "Q^{2} (2.5-3) [GeV^{2}]", 200, 0.2, 2);
	iip_hist->SetFillColor(kBlue+2);
	iip_hist->SetStats(0);
   

   TH1D *iiip_hist = new TH1D("Q^{2} (3-3.5) [GeV^{2}]", "Q^{2} (3-3.5) [GeV^{2}]", 200, 0.2, 2);
	iiip_hist->SetFillColor(kBlue+2);
	iiip_hist->SetStats(0);
   

   TH1D *ivp_hist = new TH1D("Q^{2} (3.5-4.5) [GeV^{2}]", "Q^{2} (3.5-4.5) [GeV^{2}]", 200, 0.2, 2);
	ivp_hist->SetFillColor(kBlue+2);
	ivp_hist->SetStats(0);
	//ivp_hist->SetFillStyle(4000);
   
   TH1D *vp_hist = new TH1D("Q^{2} (4.5-6) [GeV^{2}]", "Q^{2} (4.5-6) [GeV^{2}]", 200, 0.2, 2);
	vp_hist->SetFillColor(kBlue+2);
	vp_hist->SetStats(0);
   /*
   TH1D *bivp_hist = new TH1D("bbbQ2 (3.5 4.5)", "Q2 (3.5 4.5)", 200, 0.3, 2.0);
	bivp_hist->SetFillColor(kBlue-1);
   
   TH1D *vip_hist = new TH1D("bbbQ2 (4.56)", "Q2 (4.56)", 200, 0.3, 2.0);
        vip_hist->SetFillColor(kBlue-1);
     */   
       // TH1D *vizero_hist = new TH1D("","",200,,);
       // TFile b("zrm_rg_b_18451_4_15_cxc.root", "update");

        TH1D *tp_hist = (TH1D*)input.Get("Q2 (1->2)");  //"Q^{2} (1 2) [GeV^{2}]
      TH1D *tip_hist = (TH1D*)input.Get("Q2 (2->2.5)");
      TH1D *tiip_hist = (TH1D*)input.Get("Q2 (2.5->3)");
      TH1D *tiiip_hist = (TH1D*)input.Get("Q2 (3->3.5)");
      TH1D *tivp_hist = (TH1D*)input.Get("Q2 (3.5->4.5)");
      TH1D *tvp_hist = (TH1D*)input.Get("Q2 (4.5->6)");
     // TH1D *tvip_hist = (TH1D*)input.Get("Q2 (1->2)");
      p_hist->Add(tp_hist);
      //p_hist->Scale(2.0 / p_hist->Integral("width"));
      ip_hist->Add(tip_hist);
      //ip_hist->Scale(2.0 / ip_hist->Integral("width"));
      iip_hist->Add(tiip_hist);
      //iip_hist->Scale(2.0 / iip_hist->Integral("width"));
      iiip_hist->Add(tiiip_hist);
      //iiip_hist->Scale(2.0 / iiip_hist->Integral("width"));
      ivp_hist->Add(tivp_hist);
      //ivp_hist->Scale(2.0 / ivp_hist->Integral("width"));
      vp_hist->Add(tvp_hist);
      //vp_hist->Scale(2.0 / vp_hist->Integral("width"));
    //  vp_hist->Add(tvip_hist);
     
  


   double x1 = 0.6;
   double x2 = 1.0;
  

   p_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //p_hist->Scale(2.0 / p_hist->Integral());
   TF1 *p_func = new TF1("Signal + Background",fitfbw,0.6,1.0,6);
   p_func->SetParameter(0,100);
   p_func->SetParameter(1,0.77);
   p_func->SetParNames("double pConstant","double pmean","double pa1","double pa2","double pa3", "double pa4");
   p_func->SetLineColor(kRed);
   p_hist->Fit("Signal + Background","","",0.6,1.0);
   //p_hist->Draw();
   //p_hist->Write();
   TF1 *fbw1 = new TF1("Breit-Wigner",funcy, 0.6,1.0, 2);
   fbw1->SetParameter(0, p_func->GetParameter(0));
   fbw1->SetParameter(1, p_func->GetParameter(1));
   fbw1->SetLineColor(kCyan-9);
   TF1 *fp1 = new TF1("fp1",funcp, 0.6,1.0,4);
   fp1->SetParameter(0, p_func->GetParameter(2));
   fp1->SetParameter(1, p_func->GetParameter(3));
   fp1->SetParameter(2, p_func->GetParameter(4));
   fp1->SetParameter(3, p_func->GetParameter(5));
   fp1->SetLineColor(kBlue);
   TFitResultPtr r1 = p_hist->Fit(p_func, "S", "",0.6,1.0); 
   TMatrixDSym c1_exp = r1->GetCovarianceMatrix(); 
   //signalFcn->IntegralError(xmin,xmax, signalFcn->GetParameters(),c1_exp.GetMatrixArray())))/h->GetBinWidth(1);
   
   double ybwy1 = fbw1->Integral(fbw1->GetParameter(1) - 3*(0.15) , fbw1->GetParameter(1) + 3*(0.15)) / 0.009;
   
   
   double eybwy1 = fbw1->IntegralError(0.6,1.0, fbw1->GetParameters(), c1_exp.GetMatrixArray()) / 0.009;
   //std::cout << "LD2 bin1 = " << ybwy1 << " +/- " << eybwy1 << std::endl;
   
   auto c1le = new TLegend(0.6,0.75,0.9,0.9);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   c1le->AddEntry(p_func,"Signal + Background","l");
   c1le->AddEntry(fbw1,"Breit-Wigner","l");
   
   TCanvas *c1 = new TCanvas("LD2 bin 1","",0,0,640,480);
      TLine *c1l=new TLine(fbw1->GetParameter(1) - 3*(0.15), c1->GetUymin(), fbw1->GetParameter(1) - 3*(0.15), c1->GetUymax()-0.2);
      c1l->SetLineColor(kBlue - 10);
      TLine *c1u=new TLine(fbw1->GetParameter(1) + 3*(0.15), c1->GetUymin(), fbw1->GetParameter(1) + 3*(0.15), c1->GetUymax()-0.2);
      c1u->SetLineColor(kBlue - 10);
      p_hist->Draw();
      fbw1->Draw("Same");
     // c1->BuildLegend();
     // c1l->Draw("Same");
     // c1u->Draw("Same");
      c1le->Draw();
      c1->Write();
      
     

   ip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //ip_hist->Scale(2.0 / ip_hist->Integral());
   TF1 *ip_func = new TF1("Signal + Background",fitfbw,0.6,1.0,6);
   ip_func->SetParameter(0,100);
   ip_func->SetParameter(1,0.77);
   ip_func->SetParNames("double ipConstant","double ipmean","double ipa1","double ipa2","double ipa3", "double ipa4");
   ip_hist->Fit("Signal + Background","","",0.6,1.0);
   //ip_hist->Draw();
   //ip_hist->Write();
   TF1 *fbw2 = new TF1("Breit-Wigner",funcy, x1, x2, 2);
   fbw2->SetParameter(0, ip_func->GetParameter(0));
   fbw2->SetParameter(1, ip_func->GetParameter(1));
   fbw2->SetLineColor(kCyan-9);
   TF1 *fp2 = new TF1("fp2",funcp, 0.6,1.0,4);
   fp2->SetParameter(0, ip_func->GetParameter(2));
   fp2->SetParameter(1, ip_func->GetParameter(3));
   fp2->SetParameter(2, ip_func->GetParameter(4));
   fp2->SetParameter(3, ip_func->GetParameter(5));
   TFitResultPtr r2 = ip_hist->Fit(ip_func, "S", "",0.6,1.0); 
   TMatrixDSym c2_exp = r2->GetCovarianceMatrix(); 
   //signalFcn->IntegralError(xmin,xmax, signalFcn->GetParameters(),c1_exp.GetMatrixArray())))/h->GetBinWidth(1);
   double ybwy2 = fbw2->Integral(fbw2->GetParameter(1) - 3*(0.15) , fbw2->GetParameter(1) + 3*(0.15)) / 0.009;
   double eybwy2 = fbw2->IntegralError(0.6,1.0, fbw2->GetParameters(), c2_exp.GetMatrixArray()) / 0.009;
   
     auto c2le = new TLegend(0.6,0.75,0.9,0.9);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   c2le->AddEntry(ip_func,"Signal + Background","l");
   c2le->AddEntry(fbw2,"Breit-Wigner","l");
   
   TCanvas *c2 = new TCanvas("LD2 bin 2","Q^{2} (2-2.5) [GeV^{2}]",0,0,640,480);
   TLine *c2l=new TLine(fbw2->GetParameter(1) - 3*(0.15), c2->GetUymin(), fbw2->GetParameter(1) - 3*(0.15), c2->GetUymax()-0.2);
      c2l->SetLineColor(kBlue - 10);
      TLine *c2u=new TLine(fbw2->GetParameter(1) + 3*(0.15), c2->GetUymin(), fbw2->GetParameter(1) + 3*(0.15), c2->GetUymax()-0.2);
      c2u->SetLineColor(kBlue - 10);
      ip_hist->Draw();
      fbw2->Draw("Same");
      //c2->BuildLegend();
     // c2l->Draw("Same");
     // c2u->Draw("Same");
      c2le->Draw();
      c2->Write();


   iip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //iip_hist->Scale(2.0 / iip_hist->Integral());
   TF1 *iip_func = new TF1("Signal + Background",fitfbw,0.6,1.0,6);
   iip_func->SetParameter(0,100);
   iip_func->SetParameter(1,0.77);
   iip_func->SetParNames("double iipConstant","double iipmean","double iipa1","double iipa2","double iipa3", "double iipa4");
   iip_hist->Fit("Signal + Background","","",0.6,1.0);
   //iip_hist->Draw();
   //iip_hist->Write();
   TF1 *fbw3 = new TF1("Breit-Wigner",funcy, x1, x2, 2);
   fbw3->SetParameter(0, iip_func->GetParameter(0));
   fbw3->SetParameter(1, iip_func->GetParameter(1));
   fbw3->SetLineColor(kCyan-9);
   TF1 *fp3 = new TF1("fp3",funcp, 0.6,1.0,4);
   fp3->SetParameter(0, iip_func->GetParameter(2));
   fp3->SetParameter(1, iip_func->GetParameter(3));
   fp3->SetParameter(2, iip_func->GetParameter(4));
   fp3->SetParameter(3, iip_func->GetParameter(5));
   TFitResultPtr r3 = iip_hist->Fit(iip_func, "S", "",0.6,1.0); 
   TMatrixDSym c3_exp = r3->GetCovarianceMatrix(); 
   //signalFcn->IntegralError(xmin,xmax, signalFcn->GetParameters(),c1_exp.GetMatrixArray())))/h->GetBinWidth(1);
   double ybwy3 = fbw3->Integral(fbw3->GetParameter(1) - 3*(0.15) , fbw3->GetParameter(1) + 3*(0.15)) / 0.009;
   double eybwy3 = fbw3->IntegralError(0.6,1.0, fbw3->GetParameters(), c3_exp.GetMatrixArray()) / 0.009;


  auto c3le = new TLegend(0.6,0.75,0.9,0.9);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   c3le->AddEntry(iip_func,"Signal + Background","l");
   c3le->AddEntry(fbw3,"Breit-Wigner","l");

   TCanvas *c3 = new TCanvas("LD2 bin 3","Q^{2} (2.5-3) [GeV^{2}]",0,0,640,480);
   TLine *c3l=new TLine(fbw3->GetParameter(1) - 3*(0.15), c3->GetUymin(), fbw3->GetParameter(1) - 3*(0.15), c3->GetUymax()-0.2);
      c3l->SetLineColor(kBlue - 10);
      TLine *c3u=new TLine(fbw3->GetParameter(1) + 3*(0.15), c3->GetUymin(), fbw3->GetParameter(1) + 3*(0.15), c3->GetUymax()-0.2);
      c3u->SetLineColor(kBlue - 10);
      iip_hist->Draw();
      fbw3->Draw("Same");
      //c3->BuildLegend();
    //  c3l->Draw("Same");
    //  c3u->Draw("Same");
      c3le->Draw();
      c3->Write();



 
   iiip_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //iiip_hist->Scale(2.0 / iiip_hist->Integral());
   TF1 *iiip_func = new TF1("Signal + Background",fitfbw,0.6,1.0,6);
   iiip_func->SetParameter(0,100);
   iiip_func->SetParameter(1,0.77);
   iiip_func->SetParNames("double iiipConstant","double iiipmean","double iiipa1","double iiipa2","double iiipa3", "double iiipa4");
   iiip_hist->Fit("Signal + Background","","",0.6,1.0);
   //iiip_hist->Draw();
   //iiip_hist->Write();
   TF1 *fbw4 = new TF1("Breit-Wigner",funcy, x1, x2, 2);
   fbw4->SetParameter(0, iiip_func->GetParameter(0));
   fbw4->SetParameter(1, iiip_func->GetParameter(1));
   fbw4->SetLineColor(kCyan-9);
   TF1 *fp4 = new TF1("fp4",funcp, 0.6,1.0,4);
   fp4->SetParameter(0, iiip_func->GetParameter(2));
   fp4->SetParameter(1, iiip_func->GetParameter(3));
   fp4->SetParameter(2, iiip_func->GetParameter(4));
   fp4->SetParameter(3, iiip_func->GetParameter(5));
   TFitResultPtr r4 = iiip_hist->Fit(iiip_func, "S", "",0.6,1.0); 
   TMatrixDSym c4_exp = r4->GetCovarianceMatrix(); 
   //signalFcn->IntegralError(xmin,xmax, signalFcn->GetParameters(),c1_exp.GetMatrixArray())))/h->GetBinWidth(1);
   double ybwy4 = fbw4->Integral(fbw4->GetParameter(1) - 3*(0.15) , fbw4->GetParameter(1) + 3*(0.15)) / 0.009;
   double eybwy4 = fbw4->IntegralError(0.6,1.0, fbw4->GetParameters(), c4_exp.GetMatrixArray()) / 0.009;
   
     auto c4le = new TLegend(0.6,0.75,0.9,0.9);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   c4le->AddEntry(iiip_func,"Signal + Background","l");
   c4le->AddEntry(fbw4,"Breit-Wigner","l");
   
   TCanvas *c4 = new TCanvas("LD2 bin 4","Q^{2} (3-3.5) [GeV^{2}]",0,0,640,480);
   TLine *c4l=new TLine(fbw4->GetParameter(1) - 3*(0.15), c4->GetUymin(), fbw4->GetParameter(1) - 3*(0.15), c4->GetUymax()-0.2);
      c4l->SetLineColor(kBlue - 10);
      TLine *c4u=new TLine(fbw4->GetParameter(1) + 3*(0.15), c4->GetUymin(), fbw4->GetParameter(1) + 3*(0.15), c4->GetUymax()-0.2);
      c4u->SetLineColor(kBlue - 10);
      iiip_hist->Draw();
      fbw4->Draw("Same");
      //c4->BuildLegend();
     // c4l->Draw("Same");
     // c4u->Draw("Same");
      c4le->Draw();
      c4->Write();



   ivp_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //ivp_hist->Scale(2.0 / ivp_hist->Integral());
   TF1 *ivp_func = new TF1("Signal + Background",fitfbw,0.6,1.0,6);
   ivp_func->SetParameter(0,100);
   ivp_func->SetParameter(1,0.77);
   ivp_func->SetParNames("double ivpConstant","double ivpmean","double ivpa1","double ivpa2","double ivpa3", "double ivpa4", "double ivpa5", "double ivpa6");
   ivp_hist->Fit("Signal + Background","","",0.6,1.0);
   //ivp_hist->Draw();
   //ivp_hist->Write();
   TF1 *fbw5 = new TF1("Breit-Wigner",funcy, x1, x2, 2);
   fbw5->SetParameter(0, ivp_func->GetParameter(0));
   fbw5->SetParameter(1, ivp_func->GetParameter(1));
   fbw5->SetLineColor(kCyan-9);
   TF1 *fp5 = new TF1("fp5",funcp, 0.6,1.0,4);
   fp5->SetParameter(0, ivp_func->GetParameter(2));
   fp5->SetParameter(1, ivp_func->GetParameter(3));
   fp5->SetParameter(2, ivp_func->GetParameter(4));
   fp5->SetParameter(3, ivp_func->GetParameter(5));
   TFitResultPtr r5 = ivp_hist->Fit(ivp_func, "S", "",0.6,1.0); 
   TMatrixDSym c5_exp = r5->GetCovarianceMatrix(); 
   //signalFcn->IntegralError(xmin,xmax, signalFcn->GetParameters(),c1_exp.GetMatrixArray())))/h->GetBinWidth(1);
   double ybwy5 = fbw5->Integral(fbw5->GetParameter(1) - 3*(0.15) , fbw5->GetParameter(1) + 3*(0.15)) / 0.009;
   double eybwy5 = fbw5->IntegralError(0.6,1.0, fbw5->GetParameters(), c5_exp.GetMatrixArray()) / 0.009;
   
   
     auto c5le = new TLegend(0.6,0.75,0.9,0.9);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   c5le->AddEntry(ivp_func,"Signal + Background","l");
   c5le->AddEntry(fbw5,"Breit-Wigner","l");
   
   TCanvas *c5 = new TCanvas("LD2 bin 5","Q^{2} (3.5-4.5) [GeV^{2}]",0,0,640,480);
   //TLine *c5l=new TLine(fbw5->GetParameter(1) - 3*(0.15), c5->GetUymin(), fbw5->GetParameter(1) - 3*(0.15), c5->GetUymax()-0.2);
     // c5l->SetLineColor(kBlue - 10);
     // TLine *c5u=new TLine(fbw5->GetParameter(1) + 3*(0.15), c5->GetUymin(), fbw5->GetParameter(1) + 3*(0.15), c5->GetUymax()-0.2);
     // c5u->SetLineColor(kBlue - 10);
      ivp_hist->Draw();
      fbw5->Draw("Same");
      //c5->BuildLegend();
      //c5l->Draw("Same");
      //c5u->Draw("Same");
      c5le->Draw();
      c5->Write();
      
      
      
      vp_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   //ivp_hist->Scale(2.0 / ivp_hist->Integral());
   TF1 *vp_func = new TF1("Signal + Background",fitfbw,0.6,1.0,6);
   vp_func->SetParameter(0,100);
   vp_func->SetParameter(1,0.77);
   vp_func->SetParNames("double vpConstant","double vpmean","double vpa1","double vpa2","double vpa3", "double vpa4", "double vpa5", "double vpa6");
   vp_hist->Fit("Signal + Background","","",0.6,1.0);
   //vp_hist->Draw();
   //vp_hist->Write();
   TF1 *fbw6 = new TF1("Breit-Wigner",funcy, x1, x2, 2);
   fbw6->SetParameter(0, vp_func->GetParameter(0));
   fbw6->SetParameter(1, vp_func->GetParameter(1));
   fbw6->SetLineColor(kCyan-9);
   TF1 *fp6 = new TF1("fp6",funcp, 0.6,1.0,4);
   fp6->SetParameter(0, vp_func->GetParameter(2));
   fp6->SetParameter(1, vp_func->GetParameter(3));
   fp6->SetParameter(2, vp_func->GetParameter(4));
   fp6->SetParameter(3, vp_func->GetParameter(5));
   TFitResultPtr r6 = vp_hist->Fit(vp_func, "S", "",0.6,1.0); 
   TMatrixDSym c6_exp = r6->GetCovarianceMatrix(); 
   //signalFcn->IntegralError(xmin,xmax, signalFcn->GetParameters(),c1_exp.GetMatrixArray())))/h->GetBinWidth(1);
   double ybwy6 = fbw6->Integral(fbw6->GetParameter(1) - 3*(0.15) , fbw6->GetParameter(1) + 3*(0.15)) / 0.009;
   double eybwy6 = fbw6->IntegralError(0.6,1.0, fbw6->GetParameters(), c6_exp.GetMatrixArray()) / 0.009;
   
   
     auto c6le = new TLegend(0.6,0.75,0.9,0.9);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   c6le->AddEntry(vp_func,"Signal + Background","l");
   c6le->AddEntry(fbw6,"Breit-Wigner","l");
   
   TCanvas *c6 = new TCanvas("LD2 bin 6","Q^{2} (4.5-6) [GeV^{2}]",0,0,640,480);
   //TLine *c5l=new TLine(fbw5->GetParameter(1) - 3*(0.15), c5->GetUymin(), fbw5->GetParameter(1) - 3*(0.15), c5->GetUymax()-0.2);
     // c5l->SetLineColor(kBlue - 10);
     // TLine *c5u=new TLine(fbw5->GetParameter(1) + 3*(0.15), c5->GetUymin(), fbw5->GetParameter(1) + 3*(0.15), c5->GetUymax()-0.2);
     // c5u->SetLineColor(kBlue - 10);
      vp_hist->Draw();
      fbw6->Draw("Same");
      //c5->BuildLegend();
      //c5l->Draw("Same");
      //c5u->Draw("Same");
      c6le->Draw();
      c6->Write();





   
   
   std::cout << "Cu" << std::endl;
   std::cout << "bin1 = " << ybwy1 << " +/- " << eybwy1 << std::endl;
   std::cout << "bin2 = " << ybwy2 << " +/- " << eybwy2 << std::endl;
   std::cout << "bin3 = " << ybwy3 << " +/- " << eybwy3 << std::endl;
   std::cout << "bin4 = " << ybwy4 << " +/- " << eybwy4 << std::endl;
   std::cout << "bin5 = " << ybwy5 << " +/- " << eybwy5 << std::endl;
   std::cout << "bin6 = " << ybwy6 << " +/- " << eybwy6 << std::endl;
  // std::cout << p_hist->GetBinContent(50) << std::endl;
   
  /* 
   bivp_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   //bivp_hist->Scale(2.0 / bivp_hist->Integral());
   TF1 *bivp_func = new TF1("bivbin1_fit",fitfbw,0.6,2.0,7);
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
   //bivp_func->SetParLimits(2,0.35,0.89);
   bivp_func->SetParNames("double bivpConstant","double bivpgamma","double bivpmean","double bivpa1","double bivpa2","double bivpa3", "double bivpa4", "double bivpa5", "double bivpa6");
   bivp_hist->Fit("bivbin1_fit","","",0.3,1.0);
   bivp_hist->Draw();
   bivp_hist->Write();



   vip_hist->GetXaxis()->SetTitle("#pi^{+} + #pi^{-} mass");
   //vip_hist->Scale(2.0 / vip_hist->Integral());
   TF1 *vip_func = new TF1("vibin1_fit",fitfbw,0.6,2.0,7);
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
   //vip_func->SetParLimits(2,0.35,0.89);
   vip_func->SetParNames("double vipConstant","double vipgamma","double vipmean","double vipa1","double vipa2","double vipa3", "double vipa4", "double vipa5", "double vipa6");
   vip_hist->Fit("vibin1_fit","","",0.3,1.0);
   vip_hist->Draw();
   vip_hist->Write();
   
  */ 
 


}


