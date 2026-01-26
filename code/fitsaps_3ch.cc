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

void fitsaps_3ch(){


   TFile input("cusn_calvp_no_sigma_tzS022.root", "read");    ///fits_ld2.root    cusn_kin5_dnp_4_sigma_tzS.root  cusn_kin5_pass0v5_dnp_300_no_sigma_tzS.root
   TFile input22("cu_pass0v7_cc_02.root", "read");
   TFile input2("sn_pass0v7_cc_02.root", "read");    ///fits_ld2.root
   TFile input3("cxc_pass0v7_cc_02.root", "read");    ///fits_ld2.root
   TFile input4("ld_pass0v7_cc_02.root", "read");    ///fits_ld2.root
    TFile f("sol_pass0v7_ghp_03_fits.root", "update");
   
   gStyle->SetOptStat("");
	gStyle->SetCanvasPreferGL();
	
   TH1F *vx_hist = new TH1F("vx", "All Sectors vx", 300, -5.5, 5.5);
       //vx_hist->SetStats(0);
       vx_hist->SetLineColor(kBlack);
       vx_hist->GetXaxis()->SetTitle("vx [cm]");
   TH1D *vy_hist = new TH1D("vy", "All Sectors vy", 300, -5.5, 5.5);
       //vy_hist->SetStats(0);
       vy_hist->SetLineColor(kBlack);
       vy_hist->GetXaxis()->SetTitle("vy [cm]");
   TH1D *vz_hist = new TH1D("vz", "All Sectors vz", 300, -20, 10);
       //vz_hist->SetStats(0);
       vz_hist->SetLineColor(kBlack);
       vz_hist->GetXaxis()->SetTitle("vz [cm]");

   TH1F *vx1_hist = new TH1F("Sector 1 vx", "Sector 1 vx", 300, -5.5, 5.5);
      // vx1_hist->SetStats(0);
       vx1_hist->SetLineColor(kCyan+1);
   TH1D *vy1_hist = new TH1D("Sector 1 vy", "Sector 1 vy", 300, -5.5, 5.5);
     //  vy1_hist->SetStats(0);
       vy1_hist->SetLineColor(kCyan+1);
   TH1D *vz1_hist = new TH1D("Sector 1 vz", "Sector 1 vz", 300, -20, 10);
    //   vz1_hist->SetStats(0);
       vz1_hist->SetLineColor(kWhite);
       
  TH1F *vx2_hist = new TH1F("Sector 2 vx", "Sector 2 vx", 300, -5.5, 5.5);
   //    vx2_hist->SetStats(0);
       vx2_hist->SetLineColor(kOrange+1);
   TH1D *vy2_hist = new TH1D("Sector 2 vy", "Sector 2 vy", 300, -5.5, 5.5);
    //   vy2_hist->SetStats(0);
       vy2_hist->SetLineColor(kOrange+1);
   TH1D *vz2_hist = new TH1D("Sector 2 vz", "Sector 2 vz", 300, -20, 10);
    //   vz2_hist->SetStats(0);
       vz2_hist->SetLineColor(kOrange+1);
             
   TH1F *vx3_hist = new TH1F("Sector 3 vx", "Sector 3 vx", 300, -5.5, 5.5);
    //   vx3_hist->SetStats(0);
       vx3_hist->SetLineColor(kYellow+1);
   TH1D *vy3_hist = new TH1D("Sector 3 vy", "Sector 3 vy", 300, -5.5, 5.5);
    //   vy3_hist->SetStats(0);
       vy3_hist->SetLineColor(kYellow+1);
   TH1D *vz3_hist = new TH1D("Sector 3 vz", "Sector 3 vz", 300, -20, 10);
   //    vz3_hist->SetStats(0);
       vz3_hist->SetLineColor(kYellow+1);
       
   TH1F *vx4_hist = new TH1F("Sector 4 vx", "Sector 4 vx", 300, -5.5, 5.5);
   //    vx4_hist->SetStats(0);
       vx4_hist->SetLineColor(kGreen+1);
   TH1D *vy4_hist = new TH1D("Sector 4 vy", "Sector 4 vy", 300, -5.5, 5.5);
   //    vy4_hist->SetStats(0);
       vy4_hist->SetLineColor(kGreen+1);
   TH1D *vz4_hist = new TH1D("Sector 4 vz", "Sector 4 vz", 300, -20, 10);
  //     vz4_hist->SetStats(0);
       vz4_hist->SetLineColor(kGreen+1);
       
   TH1F *vx5_hist = new TH1F("Sector 5 vx", "Sector 5 vx", 300, -5.5, 5.5);
   //    vx5_hist->SetStats(0);
       vx5_hist->SetLineColor(kBlue+1);
   TH1D *vy5_hist = new TH1D("Sector 5 vy", "Sector 5 vy", 300, -5.5, 5.5);
    //   vy5_hist->SetStats(0);
       vy5_hist->SetLineColor(kBlue+1);
   TH1D *vz5_hist = new TH1D("Sector 5 vz", "Sector 5 vz", 300, -20, 10);
    //   vz5_hist->SetStats(0);
       vz5_hist->SetLineColor(kBlue+1);
       
   TH1F *vx6_hist = new TH1F("Sector 6 vx", "Sector 6 vx", 300, -5.5, 5.5);
    //   vx6_hist->SetStats(0);
       vx6_hist->SetLineColor(kMagenta+1);
   TH1D *vy6_hist = new TH1D("Sector 6 vy", "Sector 6 vy", 300, -5.5, 5.5);
   //    vy6_hist->SetStats(0);
       vy6_hist->SetLineColor(kMagenta+1);
   TH1D *vz6_hist = new TH1D("Sector 6 vz", "Sector 6 vz", 300, -20, 10);
  //     vz6_hist->SetStats(0);
       vz6_hist->SetLineColor(kMagenta+1);
       
       TH1F *vxtest_hist = new TH1F("Sector 6 vx", "Sector 6 vx", 300, -5.5, 5.5);
    //   vx6_hist->SetStats(0);
       vxtest_hist->SetLineColor(kMagenta+1);
   TH1D *vytest_hist = new TH1D("Sector 6 vy", "Sector 6 vy", 300, -5.5, 5.5);
   //    vy6_hist->SetStats(0);
       vytest_hist->SetLineColor(kMagenta+1);
   TH1D *vztest_hist = new TH1D("Sector 6 vz", "Sector 6 vz", 300, -20, 10);
  //     vz6_hist->SetStats(0);
       vztest_hist->SetLineColor(kMagenta+1);
       
       TH1F *vcu_hist = new TH1F("Sector 6 vx", "Cu", 200, 0, 2);
       vcu_hist->SetStats(0);
       vcu_hist->SetLineColor(kOrange);
       vcu_hist->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}} [GeV]");
   TH1F *vsn_hist = new TH1F("Sector 6 vy", "Sn", 200, 0, 2);
       vsn_hist->SetStats(0);
       vsn_hist->SetLineColor(kMagenta+1);
       TH1F *vcxc_hist = new TH1F("Sector 6 vx", "CxC", 200, 0, 2);
       vcxc_hist->SetStats(0);
       vcxc_hist->SetLineColor(kBlack);
   TH1F *vld_hist = new TH1F("Sector 6 vy", "LD_{2}", 200, 0, 2);
       vld_hist->SetStats(0);
       vz6_hist->SetLineColor(kBlue-1);
       
       double ld_charge = -0.10371671
-0.026121232
-0.165757921
-0.175941958
-0.13575266
-0.160131587
-0.159366704
-0.145902379
-0.022946483
-0.051999923
-0.070096738
-0.052477074
-0.205105837
-0.159515121
-0.239047584
-0.236212043
-0.246259843
-0.218683524
-0.242181835
-0.148569817
-0.096008422
-0.28293769
-0.29412128
-0.294545202
-0.283528661
-0.305232182
-0.087577924
-0.080814581
-0.291759627
-0.284328937
-0.126206992
-0.302094818
-0.305426796
-0.300111809
-0.120576336
-0.292568426
-0.192820571
-0.056433115
-0.289657247
-0.2961642
-0.059440749
-0.286692186
-0.297766804
-0.065970377
-0.296426577
-0.297741217
-0.045058183
-0.03096454
-0.12938844
-0.285818958
-0.268611001
-0.289796754
-0.285739976
-0.286428448
-0.279755714
-0.288796953
-0.286882865
-0.286637881
-0.285423941
-0.063304549
-0.218542515
-0.117558183
-0.282269415
-0.055336435
-0.145719576
-0.282985363
-0.284213858
-0.287869904
-0.284637246
-0.09565723
-0.28395258
-0.280417403
-0.285385957
-0.279160035
-0.280456838
-0.28212898
-0.264316782
-0.281476517
-0.275941838
-0.159076612
-0.29320495
-0.282956072
-0.192941946
-0.283967467
-0.280668282
-0.303445092
-0.121589385
-0.260221945
-0.278845495
-0.281114361
-0.054206235
-0.129841942
-0.01906473
-0.27630575
-0.277764724
-0.279372002
-0.286091061
-0.13072159
-0.278777322
-0.275540005
-0.06957881
-0.283827709
-0.279896475
-0.141193201
-0.097773592
-0.128341349
-0.109611244
-0.123720383
-0.279053331
-0.134378774
-0.275809045
-0.276550879
-0.275094819
-0.090895018
-0.264984804
-0.286145886
-0.27477323
-0.278662732
-0.017872181
-0.117857679
-0.17354754
-0.280087273
-0.279844503
-0.280479798
-0.285327566
-0.273136111
-0.27880999
-0.278239916
-0.144163883
-0.234343511;

double cxc_charge = -0.021409655
-0.058376027
-0.059119247
-0.039595818
-0.04523127
-0.269821781
-0.253627025
-0.072084335
-0.141775711
-0.068509891
-0.16050247
-0.251116421
-0.25488341
-0.264614433
-0.254030006
-0.153630405
-0.082118515
-0.252363164
-0.253391037
-0.209178119
-0.266346409
-0.264891987
-0.256974357
-0.252957803
-0.262653302
-0.18305665
-0.135885772
-0.29355055
-0.299880335
-0.289616591
-0.297309913
-0.292706486
-0.296692881
-0.09993449
-0.293670146
-0.287469465
-0.319354027
-0.276430211
-0.291709284
-0.310829913
-0.29781652
-0.188613166
-0.30865846
-0.303587395
-0.135310337
-0.313866528
-0.305721715
-0.077346001
-0.10282063
-0.054454093
-0.214878129
-0.314215376
-0.31359421
-0.31604368
-0.064592577
-0.309300082
-0.30789066
-0.313746799
-0.296974843
-0.309476135
-0.308858337
-0.312464919
-0.096273723
-0.324944378
-0.340697249
-0.026176505
-0.327986534
-0.305314544
-0.308238744
-0.307641907
-0.141273018
-0.262791658
-0.302953285
-0.139992699
-0.028507818
-0.307454789
-0.308425591
-0.309768761
-0.204611167
-0.308332821
-0.121746468
-0.309039574
-0.309433726
-0.309431106
-0.311541426
-0.308599214
-0.309131203
-0.056640892
-0.312225748
-0.307727552
-0.326574027
-0.307579657
-0.326233115
-0.30528374
-0.307009488
-0.048047967
-0.14034691
-0.306676354
-0.307545579
-0.306179133
-0.023077337
-0.066913182
-0.309520164
-0.027657723
-0.307674878
-0.307371019
-0.300433437
-0.306470044
-0.305567195
-0.053722379
-0.114178446
-0.304888779;



double cu_charge = -0.772779703
-0.249950477
-0.845492767
-0.268954742
-0.469632849
-0.103529074
-0.07810441
-0.205324211
-0.100862341
-0.846908404
-0.89980439
-0.887859192
-0.898181153
-0.517005422
-0.211752181
-0.55201839
-0.75267132
-0.392714073
-0.264253646
-0.897019754
-0.122577029
-0.894324224
-0.894126547
-0.948261026
-0.577190948
-0.451316229
-0.176581788
-0.248013517
-0.555830285
-0.906314773
-0.087538428
-0.891223264
-0.889506114
-0.626579809
-0.080886648
-0.36193337
-0.444310359
-0.114206273
-0.191995906
-0.581194407
-0.175837083
-0.884277844
-0.211850836
-0.166422671
-0.112903694
-0.800874636
-0.894654573
-0.885933843
-0.89821577
-0.899391129
-0.89538211
-0.900416686
-0.896066703
-0.884408786
-0.980575204
-0.896307108
-0.753334719
-0.885681611
-0.023755564
-0.227654033
-0.869113911
-0.063542174
-0.044110816
-0.524458642
-0.117674537
-0.87538487
-0.87525709
-0.872092982
-0.872227782
-0.873642362
-0.32216525
-0.875888933
-0.869908082
-0.695738262
-0.746645864
-0.697946739
-0.875340883
-0.765746714
-0.512288772
-0.03786725
-0.900965214
-0.075069676
-0.123296565
-0.847240424
-0.877391596
-0.59994788
-0.440605854
-0.38053088
-0.374799672
-0.183335804
-0.674366322
-0.322992972
-0.264128331
-0.833975982
-0.840634099
-0.875483865
-0.411608518
-0.871769185
-0.801875982
-0.338178358
-0.760534103
-0.871565657
-0.873773857
-0.507866486
-0.88220889
-0.875497753
-0.049807006
-0.868283039
-0.052466931
-0.092699522
-0.519269643
-0.577278799
-0.203436031
-0.592757596
-0.882562102
-0.112113036
-0.063006983
-0.866117965
-0.867142237
-0.531255538
-0.034012117
-0.878006035
-0.594205339
-0.871105369
-0.868976504
-0.607893305
-0.673859952
-0.512973932
-0.400380502
-0.872330451
-0.866955182
-0.864326771
-0.458816089
-0.287305109
-0.412841234
-0.629889437
-0.849793838
-0.861943979
-0.852433875
-0.135693393
-0.404606722
-0.100744135
-0.848334326
-0.851198395
-0.097955765
-0.850228333
-0.066726527
-0.578626908
-0.369794684
-0.289287654
-0.233011581
-0.047072217
-0.855675581
-0.853163361
-0.806725441
-0.851977414
-0.854858001
-0.436479072
-0.852381568
-0.855308329
-0.85939577
-0.853335147
-0.854832498
-0.858654441
-0.8680825
-0.844213107
-0.85111702
-0.6426627
-0.786236694
-0.055094584
-0.850576963
-0.110120617
-0.846933328
-0.860790947
-0.413028081
-0.848320956
-0.848270198
-0.476841995
-0.072204425
-0.278563301
-0.847406713
-0.686723343
-0.850438464
-0.247774063
-0.878832592
-0.848596646
-0.00902137
-0.4410995
-0.60330036
-0.683728063
-0.267373425
-0.050192388
-0.867195291
-0.866697806
-0.070060311
-0.078879342
-0.85304265
-0.536719399
-0.692351993
-0.622786126
-0.10925462
-0.729003758
-0.620293817
-0.034714503
-0.066678223
-0.386663961
-0.852700702
-0.15915104
-0.598025811
-0.406242211
-0.850206514
-0.755313806
-0.848144078
-0.094972386
-0.071231317
-0.847684947
-0.830035548
-0.544179577
-0.398426273
-0.584913882
-0.215312189
-0.856295944
-0.849033989
-0.264906969
-0.220213047
-0.317187304
-0.847307432
-0.842323074
-0.833560918
-0.662795445
-0.845398333
-0.803479566
-0.532264289
-0.620020806
-0.21407956
-0.843652907
-0.848707909
-0.127971485
-0.846991104
-0.847227982
-0.846835781
-0.433436291
-0.452869003
-0.849080402
-0.850315202
-0.848416355
-0.84846497
-0.793466952
-0.755031339
-0.845448343
-0.846484882
-0.534493512
-0.748317093
-0.851595684
-0.189127736
-0.599037275
-0.656142944
-0.787058418
-0.847393339
-0.848700137
-0.753504037
-0.849742484
-0.086242404
-0.823693197
-0.856907792
-0.278421074
-0.441536449
-0.847410573
-0.591330742;





std::cout << "ld charge = " << ld_charge << std::endl;
std::cout << "cxc charge = " << cxc_charge << std::endl;
std::cout << "cusn charge = " << cu_charge << std::endl;



  

      TH1F *tvx1_hist = (TH1F*)input.Get("sector 1 vx"); 
      TH1D *tvy1_hist = (TH1D*)input.Get("sector 1 vy");
      TH1D *tvz1_hist = (TH1D*)input.Get("sector 1 vz");
      
      TH1F *tvx2_hist = (TH1F*)input.Get("sector 2 vx"); 
      TH1D *tvy2_hist = (TH1D*)input.Get("sector 2 vy");
      TH1D *tvz2_hist = (TH1D*)input.Get("sector 2 vz");
      
      TH1F *tvx3_hist = (TH1F*)input.Get("sector 3 vx"); 
      TH1D *tvy3_hist = (TH1D*)input.Get("sector 3 vy");
      TH1D *tvz3_hist = (TH1D*)input.Get("sector 3 vz");
      
      TH1F *tvx4_hist = (TH1F*)input.Get("sector 4 vx"); 
      TH1D *tvy4_hist = (TH1D*)input.Get("sector 4 vy");
      TH1D *tvz4_hist = (TH1D*)input.Get("sector 4 vz");
      
      TH1F *tvx5_hist = (TH1F*)input.Get("sector 5 vx"); 
      TH1D *tvy5_hist = (TH1D*)input.Get("sector 5 vy");
      TH1D *tvz5_hist = (TH1D*)input.Get("sector 5 vz");
      
      TH1F *tvx6_hist = (TH1F*)input.Get("sector 6 vx"); 
      TH1D *tvy6_hist = (TH1D*)input.Get("sector 6 vy");
      TH1D *tvz6_hist = (TH1D*)input.Get("sector 6 vz");
      
      TH1F *cu_hist = (TH1F*)input22.Get("Q2 (3.5->4)"); 
      TH1F *sn_hist = (TH1F*)input2.Get("Q2 (3.5->4)");
      TH1F *cxc_hist = (TH1F*)input3.Get("Q2 (3.5->4)");
      TH1F *ld_hist = (TH1F*)input4.Get("Q2 (3.5->4)"); 
      
      vcu_hist->Add(cu_hist);
      vcu_hist->Scale(1.0 / vcu_hist->Integral("width"));
      vcu_hist->SetMinimum(0.0);
      vcu_hist->SetMaximum(5);
      vcu_hist->SetLineColor(kOrange + 1);
      vcu_hist->SetLineWidth(3);
      vsn_hist->Add(sn_hist);
      vsn_hist->Scale(1.0 / vsn_hist->Integral("width"));
      vsn_hist->SetMinimum(0.0);
      vsn_hist->SetMaximum(5);
      vsn_hist->SetLineColor(kMagenta + 1);
      vsn_hist->SetLineWidth(3);
      vcxc_hist->Add(cxc_hist);
      vcxc_hist->Scale(1.0 / vcxc_hist->Integral("width"));
      vcxc_hist->SetMinimum(0.0);
      vcxc_hist->SetMaximum(5);
      vcxc_hist->SetLineColor(kGreen + 1);
      vcxc_hist->SetLineWidth(3);
      vld_hist->Add(ld_hist);
      vld_hist->Scale(1.0 / vld_hist->Integral("width"));
      vld_hist->SetMinimum(0.0);
      vld_hist->SetMaximum(5);
      vld_hist->SetLineColor(kBlue + 1);
      vld_hist->SetLineWidth(3);
      
      
      vx_hist->Add(tvx1_hist);
      vx_hist->Add(tvx2_hist);
      vx_hist->Add(tvx3_hist);
      vx_hist->Add(tvx4_hist);
      vx_hist->Add(tvx5_hist);
      vx_hist->Add(tvx6_hist);
      vx_hist->Scale(1.0 / vx_hist->Integral());
      vx_hist->SetMinimum(0.0);
      vx_hist->SetMaximum(0.19);
      vy_hist->Add(tvy1_hist);
      vy_hist->Add(tvy2_hist);
      vy_hist->Add(tvy3_hist);
      vy_hist->Add(tvy4_hist);
      vy_hist->Add(tvy5_hist);
      vy_hist->Add(tvy6_hist);
      vy_hist->Scale(1.0 / vy_hist->Integral());
      vy_hist->SetMinimum(0.0);
      vy_hist->SetMaximum(0.07);
      vz_hist->Add(tvz1_hist);
      vz_hist->Add(tvz2_hist);
      vz_hist->Add(tvz3_hist);
      vz_hist->Add(tvz4_hist);
      vz_hist->Add(tvz5_hist);
      vz_hist->Add(tvz6_hist);
      vz_hist->Scale(1.0 / vz_hist->Integral());
      vz_hist->SetMinimum(0.0);
      vz_hist->SetMaximum(0.05);
      
      
      
      vx1_hist->Add(tvx1_hist);
      vx1_hist->Scale(1.0 / vx1_hist->Integral());
      vx1_hist->SetMinimum(0.0);
      vx1_hist->SetMaximum(0.19);
      vy1_hist->Add(tvy1_hist);
      vy1_hist->Scale(1.0 / vy1_hist->Integral());
      vy1_hist->SetMinimum(0.0);
      vy1_hist->SetMaximum(0.07);
      vz1_hist->Add(tvz1_hist);
      vz1_hist->Scale(1.0 / vz1_hist->Integral());
      vz1_hist->SetMinimum(0.0);
      vz1_hist->SetMaximum(0.05);
      
      vx2_hist->Add(tvx2_hist);
      vx2_hist->Scale(1.0 / vx2_hist->Integral());
      vx2_hist->SetMinimum(0.0);
      vx2_hist->SetMaximum(0.19);
      vy2_hist->Add(tvy2_hist);
      vy2_hist->Scale(1.0 / vy2_hist->Integral());
      vy2_hist->SetMinimum(0.0);
      vy2_hist->SetMaximum(0.07);
      vz2_hist->Add(tvz2_hist);
      vz2_hist->Scale(1.0 / vz2_hist->Integral());
      vz2_hist->SetMinimum(0.0);
      vz2_hist->SetMaximum(0.05);
      
      vx3_hist->Add(tvx3_hist);
      vx3_hist->Scale(1.0 / vx3_hist->Integral());
      vx3_hist->SetMinimum(0.0);
      vx3_hist->SetMaximum(0.19);
      vy3_hist->Add(tvy3_hist);
      vy3_hist->Scale(1.0 / vy3_hist->Integral());
      vy3_hist->SetMinimum(0.0);
      vy3_hist->SetMaximum(0.07);
      vz3_hist->Add(tvz3_hist);
      vz3_hist->Scale(1.0 / vz3_hist->Integral());
      vz3_hist->SetMinimum(0.0);
      vz3_hist->SetMaximum(0.05);
      
      vx4_hist->Add(tvx4_hist);
      vx4_hist->Scale(1.0 / vx4_hist->Integral());
      vx4_hist->SetMinimum(0.0);
      vx4_hist->SetMaximum(0.19);
      vy4_hist->Add(tvy4_hist);
      vy4_hist->Scale(1.0 / vy4_hist->Integral());
      vy4_hist->SetMinimum(0.0);
      vy4_hist->SetMaximum(0.07);
      vz4_hist->Add(tvz4_hist);
      vz4_hist->Scale(1.0 / vz4_hist->Integral());
      vz4_hist->SetMinimum(0.0);
      vz4_hist->SetMaximum(0.07);
      
      vx5_hist->Add(tvx5_hist);
      vx5_hist->Scale(1.0 / vx5_hist->Integral());
      vx5_hist->SetMinimum(0.0);
      vx5_hist->SetMaximum(0.19);
      vy5_hist->Add(tvy5_hist);
      vy5_hist->Scale(1.0 / vy5_hist->Integral());
      vy5_hist->SetMinimum(0.0);
      vy5_hist->SetMaximum(0.07);
      vz5_hist->Add(tvz5_hist);
      vz5_hist->Scale(1.0 / vz5_hist->Integral());
      vz5_hist->SetMinimum(0.0);
      vz5_hist->SetMaximum(0.05);
      
      vx6_hist->Add(tvx6_hist);
      vx6_hist->Scale(1.0 / vx6_hist->Integral());
      vx6_hist->SetMinimum(0.0);
      vx6_hist->SetMaximum(0.19);
      vy6_hist->Add(tvy6_hist);
      vy6_hist->Scale(1.0 / vy6_hist->Integral());
      vy6_hist->SetMinimum(0.0);
      vy6_hist->SetMaximum(0.07);
      vz6_hist->Add(tvz6_hist);
      vz6_hist->Scale(1.0 / vz6_hist->Integral());
      vz6_hist->SetMinimum(0.0);
      vz6_hist->SetMaximum(0.05);
      
      /*
      vx_hist->Add(vx1_hist);
      vx_hist->Add(vx2_hist);
      vx_hist->Add(vx3_hist);
      vx_hist->Add(vx4_hist);
      vx_hist->Add(vx5_hist);
      vx_hist->Add(vx6_hist);
      //vx_hist->Scale(1.0 / vx_hist->Integral());
      vx_hist->SetMinimum(0.0);
      vx_hist->SetMaximum(0.19);
      vy_hist->Add(vy1_hist);
      vy_hist->Add(vy2_hist);
      vy_hist->Add(vy3_hist);
      vy_hist->Add(vy4_hist);
      vy_hist->Add(vy5_hist);
      vy_hist->Add(vy6_hist);
      //vy_hist->Scale(1.0 / vy_hist->Integral());
      vy_hist->SetMinimum(0.0);
      vy_hist->SetMaximum(0.19);
      vz_hist->Add(vz1_hist);
      vz_hist->Add(vz2_hist);
      vz_hist->Add(vz3_hist);
      vz_hist->Add(vz4_hist);
      vz_hist->Add(vz5_hist);
      vz_hist->Add(vz6_hist);
      //vz_hist->Scale(1.0 / vz_hist->Integral());
      vz_hist->SetMinimum(0.0);
      vz_hist->SetMaximum(0.05);
      */
      
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
      vx3_hist->Draw("HIST Same");
      vx4_hist->Draw("HIST Same");
      vx5_hist->Draw("HIST Same");
      vx6_hist->Draw("HIST Same");
      c1->BuildLegend();
      ilx->Draw("HIST Same");
      iux->Draw("HIST Same");
      iilx->Draw("HIST Same");
      iiux->Draw("HIST Same");
      iiilx->Draw("HIST Same");
      iiiux->Draw("HIST Same");
      c1->Write();
      
      TCanvas *c2 = new TCanvas("c2","vy",0,0,640,480);
      //std::cout << "max height " << c2->GetUymax() << std::endl;
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
      vy3_hist->Draw("HIST Same");
      vy4_hist->Draw("HIST Same");
      vy5_hist->Draw("HIST Same");
      vy6_hist->Draw("HIST Same");
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
      vz3_hist->Draw("HIST Same");
      vz4_hist->Draw("HIST Same");
      vz5_hist->Draw("HIST Same");
      vz6_hist->Draw("HIST Same");
      c3->BuildLegend();
      c3->Write();
      
      TCanvas *c4 = new TCanvas("Invariant Mass","Invariant Mass",0,0,640,480);
      vcu_hist->Draw("HIST");
      vsn_hist->Draw("HIST Same");
      vcxc_hist->Draw("HIST Same");
      vld_hist->Draw("HIST Same");
      c4->BuildLegend();
      c4->Write();
      
      TCanvas *c5 = new TCanvas("c5","vx",0,0,640,480);
      vx_hist->Draw("HIST");
      vxtest_hist->Draw("HIST Same");
      //vx_hist->Write();
      c5->Write();
      
      TCanvas *c6 = new TCanvas("c6","vy",0,0,640,480);
      vy_hist->Draw("HIST");
      vytest_hist->Draw("HIST Same");
      //vy_hist->Write();
      c6->Write();
      
      TCanvas *c7 = new TCanvas("c7","vz",0,0,640,480);
      vz_hist->Draw("HIST");
      vztest_hist->Draw("HIST Same");
      //vz_hist->Write();
      c7->Write();
      
      
      vx_hist->Write();
      vy_hist->Write();
      
       TF1 *g1 = new TF1 ("m1", "gaus", -10.0, -5);
  g1->SetLineColor(kRed);
  TF1 *g2 = new TF1 ("m2", "gaus", -5, -0.0);
  g2->SetLineColor(kGreen);
  TF1 *f1 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f1->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f1->SetLineColor(kBlue);
  
  gStyle->SetOptFit(1);
  
  vz_hist->Fit(g1, "R");
  vz_hist->Fit(g2, "R+");
  
  Double_t par[8];
  g1->GetParameters(&par[0]);
  g2->GetParameters(&par[3]);
  f1->SetParameters(par);
  
  vz_hist->Fit(f1, "R+");
  // h->Fit(f1, "+", "e1", 4, 16);
  vz_hist->Draw("e1");
      vz_hist->Write();
     
  TF1 *g11 = new TF1 ("m11", "gaus", -10.0, -5);
  g11->SetLineColor(kRed);
  TF1 *g21 = new TF1 ("m21", "gaus", -5, -0.0);
  g21->SetLineColor(kGreen);
  TF1 *f11 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f11->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f11->SetLineColor(kBlue);
  
  gStyle->SetOptFit(1);
  
  tvz1_hist->Fit(g11, "R");
  tvz1_hist->Fit(g21, "R+");
  
  Double_t par1[8];
  g11->GetParameters(&par1[0]);
  g21->GetParameters(&par1[3]);
  f11->SetParameters(par1);
  
  tvz1_hist->Fit(f11, "R+");
  // h->Fit(f1, "+", "e1", 4, 16);
  tvz1_hist->Draw("e1");
      tvz1_hist->Write();
      
      
   TF1 *g12 = new TF1 ("m12", "gaus", -10.0, -5);
  g12->SetLineColor(kRed);
  TF1 *g22 = new TF1 ("m22", "gaus", -5, -0.0);
  g22->SetLineColor(kGreen);
  TF1 *f12 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f12->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f12->SetLineColor(kBlue);
  
  gStyle->SetOptFit(1);
  
  tvz2_hist->Fit(g12, "R");
  tvz2_hist->Fit(g22, "R+");
  
  Double_t par2[8];
  g12->GetParameters(&par2[0]);
  g22->GetParameters(&par2[3]);
  f12->SetParameters(par2);
  
  tvz2_hist->Fit(f12, "R+");
  // h->Fit(f1, "+", "e1", 4, 16);
  tvz2_hist->Draw("e1");
      tvz2_hist->Write();


  
  
  TF1 *g13 = new TF1 ("m13", "gaus", -10.0, -5);
  g13->SetLineColor(kRed);
  TF1 *g23 = new TF1 ("m23", "gaus", -5, -0.0);
  g23->SetLineColor(kGreen);
  TF1 *f13 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f13->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f13->SetLineColor(kBlue);
  
  gStyle->SetOptFit(1);
  
  tvz3_hist->Fit(g13, "R");
  tvz3_hist->Fit(g23, "R+");
  
  Double_t par3[8];
  g13->GetParameters(&par3[0]);
  g23->GetParameters(&par3[3]);
  f13->SetParameters(par3);
  
  tvz3_hist->Fit(f13, "R+");
  tvz3_hist->Draw("e1");
  tvz3_hist->Write(); 
  
  
  
  TF1 *g14 = new TF1 ("m14", "gaus", -10.0, -5);
  g14->SetLineColor(kRed);
  TF1 *g24 = new TF1 ("m24", "gaus", -5, -0.0);
  g24->SetLineColor(kGreen);
  TF1 *f14 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f14->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f14->SetLineColor(kBlue);
  
  gStyle->SetOptFit(1);
  
  tvz4_hist->Fit(g14, "R");
  tvz4_hist->Fit(g24, "R+");
  
  Double_t par4[8];
  g14->GetParameters(&par4[0]);
  g24->GetParameters(&par4[3]);
  f14->SetParameters(par4);
  
  tvz4_hist->Fit(f14, "R+");
  tvz4_hist->Draw("e1");
  tvz4_hist->Write(); 
  
  
  
  
  TF1 *g15 = new TF1 ("m15", "gaus", -10.0, -5);
  g15->SetLineColor(kRed);
  TF1 *g25 = new TF1 ("m25", "gaus", -5, -0.0);
  g25->SetLineColor(kGreen);
  TF1 *f15 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f15->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f15->SetLineColor(kBlue);
  
  gStyle->SetOptFit(1);
  
  tvz5_hist->Fit(g15, "R");
  tvz5_hist->Fit(g25, "R+");
  
  Double_t par5[8];
  g15->GetParameters(&par5[0]);
  g25->GetParameters(&par5[3]);
  f15->SetParameters(par5);
  
  tvz5_hist->Fit(f15, "R+");
  tvz5_hist->Draw("e1");
  tvz5_hist->Write(); 
 


  TF1 *g16 = new TF1 ("m16", "gaus", -10.0, -5);
  g16->SetLineColor(kRed);
  TF1 *g26 = new TF1 ("m23", "gaus", -5, -0.0);
  g26->SetLineColor(kGreen);
  TF1 *f16 = new TF1("double_gaus", "gaus(0) + gaus(3) + [6] + [7]*x", -10, 0.0);
  f16->SetParNames("Constant 1", "Mean 1", "Sigma 1",
                  "Constant 2", "Mean 2", "Sigma 2");
  f16->SetLineColor(kBlue);
  
  gStyle->SetOptFit(1);
  
  tvz6_hist->Fit(g16, "R");
  tvz6_hist->Fit(g26, "R+");
  
  Double_t par6[8];
  g16->GetParameters(&par6[0]);
  g26->GetParameters(&par6[3]);
  f16->SetParameters(par6);
  
  tvz6_hist->Fit(f16, "R+");
  tvz6_hist->Draw("e1");
  tvz6_hist->Write(); 

}


