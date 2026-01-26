#include <cstdlib>
#include <iostream>
#include <string>


void graph_cusn(){
        //double con = 1651.83;
	//double gamma = 0.155213;
	//double mean = 0.76275;
	//double a1 = -30445.3;
	//double a2 = 256975.0;
	//double a3 = -817571.0;
	//double a4 = 1241440.0;
	//double a5 = -887689.0;
	//double a6 = 240180.0;
	double x = 0.2;
	
	double pbw = 0.0;
	double ppol = 0.0;
	
	double ipbw = 0.0;
	double ippol = 0.0;
	
	double iipbw = 0.0;
	double iippol = 0.0;
	
	double iiipbw = 0.0;
	double iiippol = 0.0;
	
	double ivpbw = 0.0;
	double ivppol = 0.0;
	
	double vpbw = 0.0;
	double vppol = 0.0;
	
	double vipbw = 0.0;
	double vippol = 0.0;
	
	
double pConstant                 =      1651.83; //  +/-   28.354      
double pgamma                    =     0.155213; //  +/-   0.00195911    int = 1966.09   rat = 0.0039
double pmean                     =      0.76275; //  +/-   0.000340017  	 (limited)      integral = 
double pa1                       =     -30445.3; //  +/-   1062.1      
double pa2                       =       256975; //  +/-   8961.76     
double pa3                       =      -817571; //  +/-   28911       
double pa4                       =  1.24144e+06; //  +/-   44811.6     
double pa5                       =      -887689; //  +/-   33276.5     
double pa6                       =       240180; //  +/-   9443.89     
//****************************************
//Minimizer is Minuit2 / Migrad
//Chi2                      =      118.362
//NDf                       =           81
//Edm                       =  5.88706e-07
//NCalls                    =         1002
double ipConstant                =      223.097; //  +/-   4.01896          int = 266.725   rat = 0.00429191
double ipgamma                   =     0.152923; //  +/-   0.00290336  
double ipmean                    =     0.763404; //  +/-   0.000724172  	 (limited)
double ipa1                      =     -3481.96; //  +/-   8.44367     
double ipa2                      =      28327.9; //  +/-   27.6155     
double ipa3                      =     -86848.1; //  +/-   30.4475     
double ipa4                      =       128325; //  +/-   30.1977     
double ipa5                      =     -90110.3; //  +/-   27.3555     
double ipa6                      =      24076.4; //  +/-   19.859      
//***************************************;*
//Minimizer is Minuit2 / Migra;d
//Chi2                      = ;     82.5062
//NDf                       = ;          81
//Edm                       = ; 1.68512e-07
//NCalls                    =          842
double iipConstant               =      81.9288; //  +/-   4.92204     
double iipgamma                  =     0.133513; //  +/-   0.00655357      int = 101.439        rat = 0.00366974
double iipmean                   =     0.761094;  // +/-   0.0012187    	 (limited)
double iipa1                     =      -1591.4;  // +/-   240.088     
double iipa2                     =      13176.4;  // +/-   2017.03     
double iipa3                     =     -41558.8;  // +/-   6443.83     
double iipa4                     =      63867.2;  // +/-   9840.07     
double iipa5                     =     -46701.8; ///  +/-   7195.42     
double iipa6                     =      12954.3; //  +/-   2015.9     
/*
150

Minimizer is Minuit2 / Migrad
Chi2                      =      129.131
NDf                       =           81
Edm                       =  6.59787e-06
NCalls                    =          371
double iipConstant        =      94.7681   +/-   0           
double iipgamma           =         0.15   +/-   0                     int = 114.002     rat  0.00412423
double iipmean            =      0.76488   +/-   0            	 (limited)
double iipa1              =     -105.086   +/-   0           
double iipa2              =       375.01   +/-   0           
double iipa3              =      54.6388   +/-   0           
double iipa4              =     -170.135   +/-   0           
double iipa5              =     -114.987   +/-   0           
double iipa6              =      84.6158   +/-   0           
Warning in <Fit>: Abnormal termination of minimization.
****************************************
 */
//***************************************;*
//Minimizer is Minuit2 / Migrad
//Chi2                      =      107.218
//NDf                       =           81
//Edm                       =  1.83593e-06
//NCalls                    =          494
double iiipConstant              =       41.902; //  +/-   2.03196            int = 51.2197   rat 0.00396039
double iiipgamma                 =     0.141403; //  +/-   0.00708898  
double iiipmean                  =     0.765213; //  +/-   0.00167836   	 (limited)
double iiipa1                    =     -60.4562; //  +/-   9.70804     
double iiipa2                    =      220.031; //  +/-   32.3952     
double iiipa3                    =      12.9697; //  +/-   8.57654     
double iiipa4                    =     -119.769; //  +/-   24.965      
double iiipa5                    =     -68.8728; //  +/-   10.773      
double iiipa6                    =      73.1755; //  +/-   19.5442  


/*
double iiipConstant       =      44.0643   +/-   0                  int = 53.0227       rat = 0.00410019
double iiipgamma          =         0.15   +/-   0           
double iiipmean           =     0.765549   +/-   0            	 (limited)
double iiipa1             =     -56.3565   +/-   0           
double iiipa2             =      205.865   +/-   0           
double iiipa3             =      10.0496   +/-   0           
double iiipa4             =     -114.084   +/-   0           
double iiipa5             =     -64.0078   +/-   0           
double iiipa6             =      72.1904   +/-   0     
//***************************************;*
*/
double ivpConstant               =      34.2567;  // +/-   1.21327             int = 26.3094    rat = 0.00427934
double ivpgamma                  =     0.15;  // +/-   0.00868573  
double ivpmean                   =     0.766416;  // +/-   0.00239354   	 (limited)
double ivpa1                     =     -339.148;  // +/-   2.52223     
double ivpa2                     =      2793.88;  // +/-   8.25338     
double ivpa3                     =     -8727.89;  // +/-   9.14226     
double ivpa4                     =      13136.7;  // +/-   9.08193     
double ivpa5                     =     -9410.21;  // +/-   8.23323     
double ivpa6                     =      2569.54;  // +/-   6.02937     

/*
double ivpConstant        =      21.7925   +/-   0                          int = 26.2605    rat = 0.00427139
double ivpgamma           =         0.15   +/-   0           
double ivpmean            =      0.76908   +/-   0            	 (limited)
double ivpa1              =     -18.2665   +/-   0           
double ivpa2              =      65.5027   +/-   0           
double ivpa3              =      10.2114   +/-   0           
double ivpa4              =     -31.4208   +/-   0           
double ivpa5              =     -21.6524   +/-   0           
double ivpa6              =      18.0993   +/-   0  
*/
//****************************************
//Minimizer is Minuit2 / Migrad
//Chi2                      =      76.0156
//NDf                       =           81
//Edm                       =   4.4825e-06
//NCalls                    =          875
double vpConstant                =      7.29832;   //+/-   0.595321            int = 9.37705  rat = 0.00314772
double vpgamma                   =     0.114038;  // +/-   0.0104113   
double vpmean                    =     0.763698;   //+/-   0.00321902   	 (limited)
double vpa1                      =     -110.789;  // +/-   1.55833     
double vpa2                      =      949.323;  // +/-   5.4053      
double vpa3                      =     -3190.33;  // +/-   6.37889     
double vpa4                      =      5292.02;  // +/-   6.47822     
double vpa5                      =     -4129.01;  // +/-   5.88213     
double vpa6                      =      1204.55;  // +/-   4.21301  
/*
double vpConstant         =      9.98792                                   int =  12.0251    rat = 0.00403662
double vpgamma            =         0.15   +/-   0           
double vpmean             =      0.76686   +/-   0            	 (limited)
double vpa1               =     -11.9996   +/-   0           
double vpa2               =      36.3603   +/-   0           
double vpa3               =      6.32251   +/-   0           
double vpa4               =     -14.5877   +/-   0           
double vpa5               =     -12.3923   +/-   0           
double vpa6               =      8.23894   +/-   0  
*/      
//****************************************
//Minimizer is Minuit2 / Migrad
//Chi2                      =       69.526
//NDf                       =           80
//Edm                       =  6.59736e-06
//NCalls                    =          735
double vipConstant               =      7.77908;  // +/-   0.794348           int = 9.21326
double vipgamma                  =     0.157837;  // +/-   0.0161658   
double vipmean                   =     0.762542;  // +/-   0.00434969   	 (limited)
double vipa1                     =     -21.9119;  // +/-   1.57196     
double vipa2                     =      94.0523;  // +/-   5.38083     
double vipa3                     =     -97.0014;  // +/-   6.21599     
double vipa4                     =      13.5484;  // +/-   6.23649     
double vipa5                     =      65.7141;  // +/-   5.66268     
double vipa6                     =     -41.2746;  // +/-   4.10677 

/*
double vipConstant        =      7.20368   +/-   0           
double vipgamma           =         0.15   +/-   0                     int = 8.65494
double vipmean            =     0.762163   +/-   0            	 (limited)
double vipa1              =     -11.9756   +/-   0           
double vipa2              =      33.3373   +/-   0           
double vipa3              =      12.7825   +/-   0           
double vipa4              =     -10.9316   +/-   0           
double vipa5              =     -11.3066   +/-   0           
double vipa6              =      1.87752   +/-   0  
*/
/*
double pConstant          =      2195.85 ;//  +/-   13.7669           intergal = 2604.77 
double pgamma             =     0.156555 ;//  +/-   0.00100132  
double pmean              =     0.760923 ;//  +/-   0.00024804   	 (limited)
double pa1                =     -35668.3 ;//  +/-   29.2293     
double pa2                =       300512 ;//  +/-   96.0103     
double pa3                =      -954487 ;//  +/-   108.442     
double pa4                =  1.45498e+06 ;//  +/-   108.719     
double pa5                = -1.04542e+06 ;//  +/-   99.2377     
double pa6                =       284013 ;//  +/-   71.8779 


double pConstant          =      2080.41   +/-   0                  inte = 2497.7
double pgamma             =         0.15   +/-   0           
double pmean              =     0.760658   +/-   0            	 (limited)
double pa1                =     -35547.9   +/-   0           
double pa2                =       302200   +/-   0           
double pa3                =      -970251   +/-   0           
double pa4                =  1.49772e+06   +/-   0           
double pa5                = -1.08893e+06   +/-   0           
double pa6                =       298895   +/-   0      
//****************************************
//Minimizer is Minuit2 / Migrad
//Chi2                      =      121.683
//NDf                       =           81
//Edm                       =   1.5548e-06
//NCalls                    =         1087
double ipConstant         =       273.51 ;//  +/-   4.58083          intergal = 325.119
double ipgamma            =     0.155766 ;//  +/-   0.00270433  
double ipmean             =     0.762155 ;//  +/-   0.000674074  	 (limited)
double ipa1               =     -4490.02 ;//  +/-   9.73157     
double ipa2               =      36377.8 ;//  +/-   31.6159     
double ipa3               =      -110641 ;//  +/-   34.7326     
double ipa4               =       162202 ;//  +/-   34.5117     
double ipa5               =      -113169 ;//  +/-   31.321      
double ipa6               =      30095.6 ;//  +/-   22.7192  


double ipConstant         =      264.056   +/-   0              int3 = 317.708
double ipgamma            =         0.15   +/-   0           
double ipmean             =     0.765322   +/-   0            	 (limited)
double ipa1               =     -212.794   +/-   0           
double ipa2               =      844.947   +/-   0           
double ipa3               =      189.811   +/-   0           
double ipa4               =     -313.095   +/-   0           
double ipa5               =     -253.349   +/-   0           
double ipa6               =      115.906   +/-   0     
//****************************************
//Minimizer is Minuit2 / Migrad
//Chi2                      =      117.137
//NDf                       =           81
//Edm                       =  6.20721e-07
//NCalls                    =         1042
double iipConstant        =      117.349 ;//  +/-   2.91732            intergal = 140.402
double iipgamma           =      0.15234 ;//  +/-   0.00394914  
double iipmean            =     0.762491 ;//  +/-   0.0010042    	 (limited)
double iipa1              =     -1982.07 ;//  +/-   6.34314     
double iipa2              =      15432.9 ;//  +/-   20.6716     
double iipa3              =     -44995.4 ;//  +/-   22.6164     
double iipa4              =      63562.1 ;//  +/-   22.4724     
double iipa5              =     -43005.3 ;//  +/-   20.3675     
double iipa6              =      11146.1 ;//  +/-   14.7415  


double iipConstant        =      117.616   +/-   0               inte = 141.514
double iipgamma           =         0.15   +/-   0           
double iipmean            =     0.765327   +/-   0            	 (limited)
double iipa1              =     -397.021   +/-   0           
double iipa2              =      2081.02   +/-   0           
double iipa3              =     -2664.24   +/-   0           
double iipa4              =      330.482   +/-   0           
double iipa5              =      1728.46   +/-   0           
double iipa6              =     -925.592   +/-   0    
//****************************************
//Minimizer is Minuit2 / Migrad
//Chi2                      =      96.3548
//NDf                       =           81
//Edm                       =  3.55089e-06
//NCalls                    =          781
double iiipConstant       =      44.2737  ;// +/-   3.4819          intergal = 54.8568  
double iiipgamma          =     0.133124  ;// +/-   0.00837219  
double iiipmean           =     0.761119  ;// +/-   0.00172218   	 (limited)
double iiipa1             =     -885.526  ;// +/-   186.35      
double iiipa2             =      7389.84  ;// +/-   1545.48     
double iiipa3             =     -23354.9  ;// +/-   4881.76     
double iiipa4             =      35802.9  ;// +/-   7383.36     
double iiipa5             =       -26089  ;// +/-   5357.32     
double iiipa6             =      7216.52  ;// +/-   1492.08     


double iiipConstant       =      51.4535   +/-   0                   inte = 61.9055
double iiipgamma          =         0.15   +/-   0           
double iiipmean           =     0.765224   +/-   0            	 (limited)
double iiipa1             =     -40.5626   +/-   0           
double iiipa2             =      161.567   +/-   0           
double iiipa3             =      29.4999   +/-   0           
double iiipa4             =     -67.0568   +/-   0           
double iiipa5             =      -48.571   +/-   0           
double iiipa6             =      32.6829   +/-   0 
//****************************************
//Minimizer is Minuit2 / Migrad
//Chi2                      =      73.0677
//NDf                       =           81
//Edm                       =  1.55236e-05
//NCalls                    =          745
double ivpConstant        =      23.1862  ;// +/-   1.33819            intergal = 27.8058
double ivpgamma           =      0.15171  ;// +/-   0.00924073  
double ivpmean            =     0.765026  ;// +/-   0.00226144   	 (limited)
double ivpa1              =     -12.3947  ;// +/-   2.89305     
double ivpa2              =      27.1792  ;// +/-   9.36471     
double ivpa3              =      87.9707  ;// +/-   10.3272     
double ivpa4              =      79.8874  ;// +/-   10.2777     
double ivpa5              =     -299.582  ;// +/-   9.35313     
double ivpa6              =        147.5  ;// +/-   6.80119  


 double ivpConstant        =      23.3662   +/-   0               inte = 28.1099
double ivpgamma           =         0.15   +/-   0           
double ivpmean            =     0.764994   +/-   0            	 (limited)
double ivpa1              =     -27.4377   +/-   0           
double ivpa2              =      105.055   +/-   0           
double ivpa3              =      9.40464   +/-   0           
double ivpa4              =     -56.7378   +/-   0           
double ivpa5              =     -34.8899   +/-   0           
double ivpa6              =      35.2898   +/-   0     
//****************************************
//Minimizer is Minuit2 / Migrad
//Chi2                      =      100.442
//NDf                       =           81
//Edm                       =  3.82867e-06
//NCalls                    =          905
double vpConstant         =      10.2765  ;// +/-   0.916189           intergal = 12.2504
double vpgamma            =     0.155252  ;// +/-   0.014151    
double vpmean             =     0.766629  ;// +/-   0.00361535   	 (limited)
double vpa1               =     -154.129  ;// +/-   1.8358      
double vpa2               =      1297.13  ;// +/-   6.20827     
double vpa3               =     -4236.24  ;// +/-   7.15024     
double vpa4               =      6773.61  ;// +/-   7.20096     
double vpa5               =     -5151.76  ;// +/-   6.60256     
double vpa6               =      1485.78  ;// +/-   4.83569  


double vpConstant         =      10.2406   +/-   0                   inte  = 12.3444
double vpgamma            =         0.15   +/-   0           
double vpmean             =     0.769987   +/-   0            	 (limited)
double vpa1               =      14.7936   +/-   0           
double vpa2               =     -126.021   +/-   0           
double vpa3               =      286.845   +/-   0           
double vpa4               =     -6.39728   +/-   0           
double vpa5               =     -335.874   +/-   0           
double vpa6               =      180.614   +/-   0    
//****************************************
//Minimizer is Minuit2 / Migrad
//Chi2                      =      98.0827
//NDf                       =           81
//Edm                       =  3.00513e-07
//NCalls                    =          952
double vipConstant        =      9.83009  ;// +/-   1.08424            intergal =  10.9478
double vipgamma           =      0.19281  ;// +/-   0.0197536   
double vipmean            =     0.768585  ;// +/-   0.00442418   	 (limited)
double vipa1              =      201.453  ;// +/-   1.70241     
double vipa2              =     -1847.88   ;//+/-   5.81726     
double vipa3              =      6301.85   ;//+/-   6.49115     
double vipa4              =     -9937.05  ;// +/-   6.39344     
double vipa5              =      7360.22  ;// +/-   5.79301     
double vipa6              =      -2067.7  ;// +/-   4.25135 

double vipConstant        =      6.57723   +/-   0                    inte  =  7.89291
double vipgamma           =         0.15   +/-   0           
double vipmean            =     0.759755   +/-   0            	 (limited)
double vipa1              =     -13.7721   +/-   0           
double vipa2              =      43.3038   +/-   0           
double vipa3              =      11.2283   +/-   0           
double vipa4              =      -14.741   +/-   0           
double vipa5              =     -14.9668   +/-   0           
double vipa6              =      4.68957   +/-   0  
*/	
	TFile f("zzzzzzzzzzzzzzzz_graph_lc_3D.root", "update");
	
	auto pbw_g = new TGraph();
      pbw_g->SetLineColor(kGreen);
     
      auto ppol_g = new TGraph();
      ppol_g->SetLineColor(kBlue);
      
      auto ipbw_g = new TGraph();
      ipbw_g->SetLineColor(kGreen);
     
      auto ippol_g = new TGraph();
      ippol_g->SetLineColor(kBlue);
      
      auto iipbw_g = new TGraph();
      iipbw_g->SetLineColor(kGreen);
     
      auto iippol_g = new TGraph();
      iippol_g->SetLineColor(kBlue);
      
      auto iiipbw_g = new TGraph();
      iiipbw_g->SetLineColor(kGreen);
     
      auto iiippol_g = new TGraph();
      iiippol_g->SetLineColor(kBlue);
      
      auto ivpbw_g = new TGraph();
      ivpbw_g->SetLineColor(kGreen);
      
     
      auto ivppol_g = new TGraph();
      ivppol_g->SetLineColor(kBlue);
      
      auto vpbw_g = new TGraph();
      vpbw_g->SetLineColor(kGreen);
     
      auto vppol_g = new TGraph();
      vppol_g->SetLineColor(kBlue);
      
      auto vipbw_g = new TGraph();
      vipbw_g->SetLineColor(kGreen);
     
     
     
      auto vippol_g = new TGraph();
      vippol_g->SetLineColor(kBlue);

	  double cucon [5] = { 0.0166, 0.0189, 0.0174, 0.771877, 0.0178265};
      double cumean [5] = { 0.767057, 0.767482, 0.76673, 0.771877,0.763332};
      
      double ccon [5] = { 0.00343964, 0.00373877, 0.00362233, 0.0034664,0.00347056};
      double cmean [5] = { 0.765765, 0.766401, 0.764707, 0.766536,0.768765};
      
      double dcon [5] = { 0.00347077, 0.003624, 0.00346489, 0.00338155,0.003222};
      double dmean [5] = { 0.76475, 0.76502, 0.765078, 0.766053,0.766426};
      
      //double yc [6] = {32.232, 23.2851, 14.9409, 8.4033, 9.222, 3.49978};
      //double yd [6] = {57.1249, 34.6447, 20.9128, 11.862, 10.7233, 3.80744};
      
      //double yc [6] = {1979.67, 265.73, 114.744, 53.7143, 38.4098, 8.66438};
      //double yd [6] = {2608.33, 318.479, 136.191, 62.1762, 41.674, 7.89982};
      
      //double yc [6] = {1901.48, 270.564, 114.002, 53.0227, 38.1232, 8.65};
      //double yd [6] = {2497.7, 317.708, 141.514, 61.9, 41.24, 7.89291};
      
      //double yc [6] = {1966.09, 266.725, 101.439, 51.22, 38.1232, 9.21};
      //double yd [6] = {2604.77, 325.119, 140.402, 54.8, 41.24, 10.9478};
      double q2[5] = { 1.125, 1.325, 1.625, 1.825, 2.125};
      double eq[5] = {0,0,0,0,0};
      double ta[6];
      double eta[6];
      
      for(int i = 0; i < 5; i++){
         double yc = 0.0;
         double yd = 0.0;
         double xs = 0.6;
         for(int j = 0; j <= 100000000; j++){
            yc = yc + cucon[i]*(1/(2*3.14159))*(0.15/(4*((cumean[i] - xs)*(cumean[i] - xs)) + (0.15*0.15)));
            yd = yd + dcon[i]*(1/(2*3.14159))*(0.15/(4*((dmean[i] - xs)*(dmean[i] - xs)) + (0.15*0.15)));
            xs = xs + 0.000000004;
          }
      
         ta[i] = (yc * 0.82) / (yd * 0.08064);  
         std::cout << ta[i] << std::endl;
         
         eta[i] = (ta[i] * sqrt((1/yc)+ (1/yd))) / 2;
         //nt->AddPoint(q2[i],ta);
      }
      
      auto nt = new TGraphErrors(6,q2,ta,eq,eta);
      nt->SetLineColor(kGreen-1);
      //nt->SetLineWidth(0);
      nt->SetMarkerColor(kBlue-1);
      nt->SetMarkerStyle(21);
      /*
	for( int i = 0; i <= 100; i++){
		
		pbw = pConstant * (pgamma / (4*((pmean - x)*(pmean - x)) + pow(pgamma,2))); 
		ppol = pa1 + pa2*x + pa3*pow(x,2) + pa4*pow(x,3) + pa5*pow(x,4) + pa6*pow(x,5);
		
		pbw_g->AddPoint(x,pbw);
		ppol_g->AddPoint(x,ppol);
		
		
		
		ipbw = ipConstant * (ipgamma / (4*((ipmean - x)*(ipmean - x)) + pow(ipgamma,2))); 
		ippol = ipa1 + ipa2*x + ipa3*pow(x,2) + ipa4*pow(x,3) + ipa5*pow(x,4) + ipa6*pow(x,5);
		
		ipbw_g->AddPoint(x,ipbw);
		ippol_g->AddPoint(x,ippol);
		
		
		
		iipbw = iipConstant * (iipgamma / (4*((iipmean - x)*(iipmean - x)) + pow(iipgamma,2))); 
		iippol = iipa1 + iipa2*x + iipa3*pow(x,2) + iipa4*pow(x,3) + iipa5*pow(x,4) + iipa6*pow(x,5);
		
		iipbw_g->AddPoint(x,iipbw);
		iippol_g->AddPoint(x,iippol);
		
		
		iiipbw = iiipConstant * (iiipgamma / (4*((iiipmean - x)*(iiipmean - x)) + pow(iiipgamma,2))); 
		iiippol = iiipa1 + iiipa2*x + iiipa3*pow(x,2) + iiipa4*pow(x,3) + iiipa5*pow(x,4) + iiipa6*pow(x,5);
		
		iiipbw_g->AddPoint(x,iiipbw);
		iiippol_g->AddPoint(x,iiippol);
		
		
		
		ivpbw = ivpConstant * (ivpgamma / (4*((ivpmean - x)*(ivpmean - x)) + pow(ivpgamma,2))); 
		ivppol = ivpa1 + ivpa2*x + ivpa3*pow(x,2) + ivpa4*pow(x,3) + ivpa5*pow(x,4) + ivpa6*pow(x,5);
		
		ivpbw_g->AddPoint(x,ivpbw);
		ivppol_g->AddPoint(x,ivppol);
		
		
		
		
		vpbw = vpConstant * (vpgamma / (4*((vpmean - x)*(vpmean - x)) + pow(vpgamma,2))); 
		vppol = vpa1 + vpa2*x + vpa3*pow(x,2) + vpa4*pow(x,3) + vpa5*pow(x,4) + vpa6*pow(x,5);
		
		vpbw_g->AddPoint(x,vpbw);
		vppol_g->AddPoint(x,vppol);
		
		
		
		vipbw = vipConstant * (vipgamma / (4*((vipmean - x)*(vipmean - x)) + pow(vipgamma,2))); 
		vippol = vipa1 + vipa2*x + vipa3*pow(x,2) + vipa4*pow(x,3) + vipa5*pow(x,4) + vipa6*pow(x,5);
		
		vipbw_g->AddPoint(x,vipbw);
		vippol_g->AddPoint(x,vippol);
		
		x = x + 0.01;
	} 
	*/
	
	ppol_g->SetTitle("Q^{2} (1->2) Polynomial;#pi^{+} + #pi^{-} mass;");
   	ppol_g->Draw("AP");
   	ppol_g->Write();

   	pbw_g->SetTitle("Q^{2} (1->2) Breit-Wigner;#pi^{+} + #pi^{-} mass;");
   	pbw_g->Draw("AP");
   	pbw_g->Write();

   	TMultiGraph *both_p = new TMultiGraph();
   	both_p->SetTitle("Q^{2} (1->2) fits;#pi^{+} + #pi^{-} mass;");
   	both_p->Add(ppol_g);
   	both_p->Add(pbw_g);
   	both_p->Draw(); 
        both_p->Write();
        
        ippol_g->SetTitle("Q^{2} (2->2.5) Polynomial;#pi^{+} + #pi^{-} mass;");
   	ippol_g->Draw("AP");
   	ippol_g->Write();

   	ipbw_g->SetTitle("Q^{2} (2->2.5) Breit-Wigner;#pi^{+} + #pi^{-} mass;");
   	ipbw_g->Draw("AP");
   	ipbw_g->Write();

   	TMultiGraph *both_ip = new TMultiGraph();
   	both_ip->SetTitle("Q^{2} (2->2.5) fits;#pi^{+} + #pi^{-} mass;");
   	both_ip->Add(ippol_g);
   	both_ip->Add(ipbw_g);
   	both_ip->Draw(); 
        both_ip->Write();
        
        
        iippol_g->SetTitle("Q^{2} (2.5->3) Polynomial;#pi^{+} + #pi^{-} mass;");
   	iippol_g->Draw("AP");
   	iippol_g->Write();

   	iipbw_g->SetTitle("Q^{2} (2.5->3) Breit-Wigner;#pi^{+} + #pi^{-} mass;");
   	iipbw_g->Draw("AP");
   	iipbw_g->Write();

   	TMultiGraph *both_iip = new TMultiGraph();
   	both_iip->SetTitle("Q^{2} (2.5->3) fits;#pi^{+} + #pi^{-} mass;");
   	both_iip->Add(iippol_g);
   	both_iip->Add(iipbw_g);
   	both_iip->Draw(); 
        both_iip->Write();
        
        
        
        iiippol_g->SetTitle("Q^{2} (3->3.5) Polynomial;#pi^{+} + #pi^{-} mass;");
   	iiippol_g->Draw("AP");
   	iiippol_g->Write();

   	iiipbw_g->SetTitle("Q^{2} (3->3.5) Breit-Wigner;#pi^{+} + #pi^{-} mass;");
   	iiipbw_g->Draw("AP");
   	iiipbw_g->Write();

   	TMultiGraph *both_iiip = new TMultiGraph();
   	both_iiip->SetTitle("Q^{2} (3->3.5) fits;#pi^{+} + #pi^{-} mass;");
   	both_iiip->Add(iiippol_g);
   	both_iiip->Add(iiipbw_g);
   	both_iiip->Draw(); 
        both_iiip->Write();
        
        
        
        ivppol_g->SetTitle("Q^{2} (3.5->4) Polynomial;#pi^{+} + #pi^{-} mass;");
   	ivppol_g->Draw("AP");
   	ivppol_g->Write();

   	ivpbw_g->SetTitle("Q^{2} (3.5->4.5) Breit-Wigner;#pi^{+} + #pi^{-} mass;");
   	ivpbw_g->Draw("AP");
   	ivpbw_g->Write();

   	TMultiGraph *both_ivp = new TMultiGraph();
   	both_ivp->SetTitle("Q^{2} (3.5->4) fits;#pi^{+} + #pi^{-} mass;");
   	both_ivp->Add(ivppol_g);
   	both_ivp->Add(ivpbw_g);
   	both_ivp->Draw(); 
        both_ivp->Write();
        
        
        
        
        vppol_g->SetTitle("Q^{2} (4->4.5) Polynomial;#pi^{+} + #pi^{-} mass;");
   	vppol_g->Draw("AP");
   	vppol_g->Write();

   	vpbw_g->SetTitle("Q^{2} (4->4.5) Breit-Wigner;#pi^{+} + #pi^{-} mass;");
   	vpbw_g->Draw("AP");
   	vpbw_g->Write();

   	TMultiGraph *both_vp = new TMultiGraph();
   	both_vp->SetTitle("Q^{2} (4->4.5) fits;#pi^{+} + #pi^{-} mass;");
   	both_vp->Add(vppol_g);
   	both_vp->Add(vpbw_g);
   	both_vp->Draw(); 
        both_vp->Write();
        
        
        
        vippol_g->SetTitle("Q^{2} (4.5->6) Polynomial;#pi^{+} + #pi^{-} mass;");
   	vippol_g->Draw("AP");
   	vippol_g->Write();

   	vipbw_g->SetTitle("Q^{2} (4.5->6) Breit-Wigner;#pi^{+} + #pi^{-} mass;");
   	vipbw_g->Draw("AP");
   	vipbw_g->Write();

   	TMultiGraph *both_vip = new TMultiGraph();
   	both_vip->SetTitle("Q^{2} (4.5->6) fits;#pi^{+} + #pi^{-} mass;");
   	both_vip->Add(vippol_g);
   	both_vip->Add(vipbw_g);
   	both_vip->Draw(); 
        both_vip->Write();
        
        nt->SetTitle("Nuclear Transparency;Q^{2} [GeV^{2}];T_{A}");
   	nt->Draw("AP");
   	nt->Write();
        
        /*
        std::cout << (0.82 * 1966.09) / (0.88 * 2604.77) << " " << (0.82 * 1966.09) / (0.88 * 2497.7) << std::endl;
	std::cout << (0.82 * 266.725) / (0.88 * 325.119) << " " << (0.82 * 266.725) / (0.88 * 317.708) << std::endl;
	std::cout << (0.82 * 1966.09) / (0.88 * 2604.77) << " " << (0.82 * 1966.09) / (0.88 * 2497.7) << << std::endl;
	std::cout << "" << std::endl;
	std::cout << (0.88 * 1966.09) / (0.82 * 2604.77) << " " << (0.88 * 1966.09) / (0.82 * 2497.7) << std::endl;
	std::cout << (0.88 * 266.725) / (0.82 * 325.119) << " " << (0.88 * 266.725) / (0.82 * 317.708) << std::endl;
	std::cout << (0.88 * 1966.09) / (0.82 * 2604.77) << " " << (0.88 * 1966.09) / (0.82 * 2497.7) << << std::endl;
	*/
}


/*

cxcsm_11

Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
****************************************
Minimizer is Minuit2 / Migrad
Chi2                      =        328.4
NDf                       =           81
Edm                       =  1.62272e-08
NCalls                    =          784
Constant                  =      1651.83   +/-   28.354      
gamma                     =     0.155213   +/-   0.00195911  
mean                      =      0.76275   +/-   0.000340017  	 (limited)
c1                        =     -30445.3   +/-   1062.1      
c2                        =       256975   +/-   8961.76     
c3                        =      -817571   +/-   28911       
p6                        =  1.24144e+06   +/-   44811.6     
p7                        =      -887689   +/-   33276.5     
p8                        =       240180   +/-   9443.89     
****************************************
Minimizer is Minuit2 / Migrad
Chi2                      =      118.362
NDf                       =           81
Edm                       =  5.88706e-07
NCalls                    =         1002
Constant                  =      223.097   +/-   4.01896     
Mean_Value                =     0.152923   +/-   0.00290336  
Sigma                     =     0.763404   +/-   0.000724172  	 (limited)
c1                        =     -3481.96   +/-   8.44367     
c2                        =      28327.9   +/-   27.6155     
c3                        =     -86848.1   +/-   30.4475     
p6                        =       128325   +/-   30.1977     
p7                        =     -90110.3   +/-   27.3555     
p8                        =      24076.4   +/-   19.859      
****************************************
Minimizer is Minuit2 / Migrad
Chi2                      =      82.5062
NDf                       =           81
Edm                       =  1.11688e-06
NCalls                    =          821
Constant                  =      81.9284   +/-   4.8523      
Mean_Value                =     0.133514   +/-   0.00652244  
Sigma                     =     0.761094   +/-   0.00123559   	 (limited)
c1                        =     -1591.51   +/-   241.331     
c2                        =      13177.4   +/-   2026.68     
c3                        =     -41561.8   +/-   6469.43     
p6                        =      63871.8   +/-   9867.1      
p7                        =       -46705   +/-   7204.77     
p8                        =      12955.2   +/-   2015.69     
****************************************
Minimizer is Minuit2 / Migrad
Chi2                      =      107.217
NDf                       =           81
Edm                       =  1.41234e-06
NCalls                    =          468
Constant                  =      41.9022   +/-   1.9635      
Mean_Value                =     0.141402   +/-   0.0067432   
Sigma                     =     0.765213   +/-   0.00165915   	 (limited)
c1                        =     -60.4709   +/-   9.53705     
c2                        =      220.121   +/-   31.6024     
c3                        =      12.8184   +/-   8.61823     
p6                        =     -119.759   +/-   24.586      
p7                        =     -68.7354   +/-   10.6112     
p8                        =      73.1039   +/-   19.3291     
****************************************
Minimizer is Minuit2 / Migrad
Chi2                      =      68.6856
NDf                       =           81
Edm                       =  3.29553e-06
NCalls                    =          856
Constant                  =      22.1872   +/-   1.21333     
Mean_Value                =     0.157922   +/-   0.00868613  
Sigma                     =     0.765588   +/-   0.00239276   	 (limited)
c1                        =     -339.075   +/-   2.52223     
c2                        =      2793.14   +/-   8.25339     
c3                        =     -8725.08   +/-   9.14226     
p6                        =      13131.9   +/-   9.0819      
p7                        =     -9406.35   +/-   8.23324     
p8                        =      2568.39   +/-   6.0294      
****************************************
Minimizer is Minuit2 / Migrad
Chi2                      =      76.0156
NDf                       =           81
Edm                       =  2.57121e-07
NCalls                    =          943
Constant                  =      7.29928   +/-   0.595416    
Mean_Value                =     0.114044   +/-   0.0104119   
Sigma                     =     0.763702   +/-   0.00321698   	 (limited)
c1                        =     -110.777   +/-   1.55834     
c2                        =      949.182   +/-   5.40532     
c3                        =     -3189.78   +/-   6.37893     
p6                        =      5291.03   +/-   6.47821     
p7                        =     -4128.21   +/-   5.88216     
p8                        =      1204.32   +/-   4.21307     
****************************************
Minimizer is Minuit2 / Migrad
Chi2                      =      69.5259
NDf                       =           80
Edm                       =  3.45208e-06
NCalls                    =          755
Constant                  =      7.77472   +/-   0.793977    
Mean_Value                =     0.157769   +/-   0.0161612   
Sigma                     =     0.762545   +/-   0.00434566   	 (limited)
c1                        =     -21.8893   +/-   1.57197     
c2                        =      93.8516   +/-   5.38063     
c3                        =     -96.5137   +/-   6.21595     
p6                        =      13.3987   +/-   6.23655     
p7                        =      65.2951   +/-   5.66271     
p8                        =     -41.0131   +/-   4.10672   

*/
