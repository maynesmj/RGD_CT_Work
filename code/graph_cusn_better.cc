#include <cstdlib>
#include <iostream>
#include <string>


     double cxc4mean [6] = { 0.764032, 0.765215, 0.763842, 0.764708,0.762921, 0.761806};
     double cxc4con [6] =  {331.234 , 63.5916, 30.509, 12.8173, 9.40861, 2.29666};                 //1 bspot  vz -12 - 5
     
     double ld4mean [6] = { 0.762435, 0.763823, 0.763654, 0.765803,0.758345, 0.769388};             //1 bspot
     double ld4con [6] = { 851.676, 122.189, 50.2311, 22.5951, 14.0952, 2.98078};


     double cxc3mean [6] = { 0.76385, 0.765267, 0.763896, 0.764276,0.763048, 0.764755};
     double cxc3con [6] = {332.749 , 65.6906, 31.5598, 13.543, 9.663595, 2.35268};                 //1 bspot vz -20 5
     
     double ld3mean [6] = { 0.762462, 0.76382, 0.763731, 0.765909,0.758326, 0.769388};             //1 bspot
     double ld3con [6] = { 854.012, 122.392, 50.3537, 22.6104, 14.0904, 2.98078};

     double cxc2mean [6] = { 0.746603, 0.753011, 0.752014, 0.750565,0.760022, 0.763643};
     double cxc2con [6] = { 15.9309, 11.887, 8.61356, 4.98051, 3.97888, 1.50378};                 //0.5 bspot
     
     double ld2mean [6] = { 0.737535, 0.747182, 0.74725, 0.75425,0.751862, 0.761226};             //0.5 bspot
     double ld2con [6] = { 43.3353, 27.7255, 13.5687, 7.54917, 6.43461, 2.14625};
     
     double cxc1con [6] = { 637.701, 100.633, 48.7874, 21.4715, 14.0962, 4.24545};                 //cache c
     double cxc1mean [6] = { 0.766028, 0.766471, 0.766974, 0.767237,0.7696, 0.76571};
     
     double ld1con [6] = { 1163.01, 169.421, 72.1004, 32.3734, 20.6719, 4.94969};                 //only ob cache files
     double ld1mean [6] = { 0.76441, 0.763454, 0.764234, 0.764745,0.767267, 0.76603};
     
     double cchcon [6] = { cxc3con[0], cxc3con[1], cxc3con[2], cxc3con[3], cxc3con[4], cxc3con[5]};         
     double cchmean [6] = { cxc3mean[0], cxc3mean[1], cxc3mean[2], cxc3mean[3],cxc3mean[4], cxc3mean[5]};
     
     double dchocon [6] = { ld3con[0], ld3con[1], ld3con[2], ld3con[3], ld3con[4], ld3con[5]};            
     double dchomean [6] = {ld3mean[0], ld3mean[1], ld3mean[2], ld3mean[3],ld3mean[4], ld3mean[5]};
     
double funcy(double *xf, double *par){
         double bwy = par[0]*(0.15/(4*((par[1] - xf[0])*(par[1] - xf[0])) + (0.15*0.15)));     
     return bwy;
     }
     
//TFitResultPtr r = h->Fit(fitfunc, "S", "",xmin,xmax); 
// TMatrixDSym c_exp = r->GetCovarianceMatrix(); 
//signalFcn->IntegralError(xmin,xmax, signalFcn->GetParameters(),c_exp.GetMatrixArray())))/h->GetBinWidth(1);
 
void graph_cusn_better(){

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
	double x1;
	double x2;
	double xf;
	
       TF1 *fbw1 = new TF1("fbw1",funcy, x1, x2, 2);
    	fbw1->SetParameter(0,cchcon[0]);
    	fbw1->SetParameter(1,cchmean[0]);
    	
    	double ybwy1 = fbw1->Integral(0.3, 1.4) / 0.01;
    	//double eybwy1 = fbw1->IntegralError(0.3, 1.4, fbw1->GetParameters()) / 0.01;
    //    std::cout << "CxC bin1 = " << ybwy1 << std::endl;
        
        TF1 *fbw2 = new TF1("fbw2",funcy, x1, x2, 2);
    	fbw2->SetParameter(0,cchcon[1]);
    	fbw2->SetParameter(1,cchmean[1]);
    	double ybwy2 = fbw2->Integral(0.3, 1.4);
     //   std::cout << "CxC bin2 = " << ybwy2/0.01 << std::endl;
        
        TF1 *fbw3 = new TF1("fbw3",funcy, x1, x2, 2);
    	fbw3->SetParameter(0,cchcon[2]);
    	fbw3->SetParameter(1,cchmean[2]);
    	double ybwy3 = fbw3->Integral(0.3, 1.4);
        std::cout << "CxC bin3 = " << ybwy3/0.01 << std::endl;
        
        TF1 *fbw4 = new TF1("fbw4",funcy, x1, x2, 2);
    	fbw4->SetParameter(0,cchcon[3]);
    	fbw4->SetParameter(1,cchmean[3]);
    	double ybwy4 = fbw4->Integral(0.3, 1.4);
     //   std::cout << "CxC bin4 = " << ybwy4/0.01 << std::endl;
        
        TF1 *fbw5 = new TF1("fbw5",funcy, x1, x2, 2);
    	fbw5->SetParameter(0,cchcon[4]);
    	fbw5->SetParameter(1,cchmean[4]);
    	double ybwy5 = fbw5->Integral(0.3, 1.4);
        std::cout << "CxC bin5 = " << ybwy5/0.01 << std::endl;
        
        TF1 *fbw6 = new TF1("fbw6",funcy, x1, x2, 2);
    	fbw6->SetParameter(0,cchcon[5]);
    	fbw6->SetParameter(1,cchmean[5]);
    	double ybwy6 = fbw6->Integral(0.3, 1.4);
     //   std::cout << "CxC bin6 = " << ybwy6/0.01 << std::endl;
        
        TF1 *lfbw1 = new TF1("lfbw1",funcy, x1, x2, 2);
    	lfbw1->SetParameter(0,dchocon[0]);
    	lfbw1->SetParameter(1,dchomean[0]);
    	double lybwy1 = lfbw1->Integral(0.3, 1.4);
        std::cout << "LD2 bin1 = " << lybwy1/0.01 << std::endl;
        
        TF1 *lfbw2 = new TF1("lfbw2",funcy, x1, x2, 2);
    	lfbw2->SetParameter(0,dchocon[1]);
    	lfbw2->SetParameter(1,dchomean[1]);
    	double lybwy2 = lfbw2->Integral(0.3, 1.4);
      //  std::cout << "LD2 bin2 = " << lybwy2/0.01 << std::endl;
        
        TF1 *lfbw3 = new TF1("lfbw3",funcy, x1, x2, 2);
    	lfbw3->SetParameter(0,dchocon[2]);
    	lfbw3->SetParameter(1,dchomean[2]);
    	double lybwy3 = lfbw3->Integral(0.3, 1.4);
      //  std::cout << "LD2 bin3 = " << lybwy3/0.01 << std::endl;
        
        TF1 *lfbw4 = new TF1("lfbw4",funcy, x1, x2, 2);
    	lfbw4->SetParameter(0,dchocon[3]);
    	lfbw4->SetParameter(1,dchomean[3]);
    	double lybwy4 = lfbw4->Integral(0.3, 1.4);
       // std::cout << "LD2 bin4 = " << lybwy4/0.01 << std::endl;
        
        TF1 *lfbw5 = new TF1("lfbw5",funcy, x1, x2, 2);
    	lfbw5->SetParameter(0,dchocon[4]);
    	lfbw5->SetParameter(1,dchomean[4]);
    	double lybwy5 = lfbw5->Integral(0.3, 1.4);
       // std::cout << "LD2 bin5 = " << lybwy5/0.01 << std::endl;
        
        TF1 *lfbw6 = new TF1("lfbw6",funcy, x1, x2, 2);
    	lfbw6->SetParameter(0,dchocon[5]);
    	lfbw6->SetParameter(1,dchomean[5]);
    	double lybwy6 = lfbw6->Integral(0.3, 1.4);
       // std::cout << "LD2 bin6 = " << lybwy6/0.01 << std::endl;
    	
    	double ydcon [6] = { lybwy1, lybwy2, lybwy3, lybwy4, lybwy5, lybwy6};                 
        double yccon [6] = { ybwy1, ybwy2, ybwy3, ybwy4, ybwy5, ybwy6};   
 
	TFile f("Nuclear_transparency_for_Matheiu.root", "update");
	
	auto pbw_g = new TGraph();
      pbw_g->SetLineColor(kYellow);
     pbw_g->SetLineWidth(5);
      auto ppol_g = new TGraph();
      ppol_g->SetLineColor(kBlue);
      
      auto ipbw_g = new TGraph();
      ipbw_g->SetLineColor(kYellow);
       ipbw_g->SetLineWidth(5);
     
      auto ippol_g = new TGraph();
      ippol_g->SetLineColor(kBlue);
      
      auto iipbw_g = new TGraph();
      iipbw_g->SetLineColor(kYellow);
       iipbw_g->SetLineWidth(5);
     
      auto iippol_g = new TGraph();
      iippol_g->SetLineColor(kBlue);
      
      auto iiipbw_g = new TGraph();
      iiipbw_g->SetLineColor(kYellow);
       iiipbw_g->SetLineWidth(5);
     
      auto iiippol_g = new TGraph();
      iiippol_g->SetLineColor(kBlue);
      
      auto ivpbw_g = new TGraph();
      ivpbw_g->SetLineColor(kYellow);
       ivpbw_g->SetLineWidth(5);
      
     
      auto ivppol_g = new TGraph();
      ivppol_g->SetLineColor(kBlue);
      
      auto vpbw_g = new TGraph();
      vpbw_g->SetLineColor(kYellow);
       vpbw_g->SetLineWidth(5);
     
      auto vppol_g = new TGraph();
      vppol_g->SetLineColor(kBlue);
      
      auto vipbw_g = new TGraph();
      vipbw_g->SetLineColor(kYellow);
     
     
     
      auto vippol_g = new TGraph();
      vippol_g->SetLineColor(kBlue);
      
      //double ccon [6] = {169.616, 122.177, 78.5446, 44.0228, 48.1053, 18.3106};
      //double cmean [6] = {0.750966, 0.755786, 0.752569, 0.758687, 0.767991, 0.761251};
      
      //double dcon [6] = {301.404, 182.462, 109.8, 62.1669, 56.1386, 19.9327};
      //double dmean [6] = {0.747041, 0.749702, 0.754685, 0.757929, 0.760009, 0.760066};
      
      ////double yc [6] = {2330.27, 377.003, 169.162, 76.5722, 56.3637, 13.8727};
      ////double yd [6] = {4190.71, 604.823, 264.057, 117.085, 78.1214, 17.7632};
      
      //double yc [6] = {1979.67, 265.73, 114.744, 53.7143, 38.4098, 8.66438};
      //double yd [6] = {2608.33, 318.479, 136.191, 62.1762, 41.674, 7.89982};
      
      //double yc [6] = {1901.48, 270.564, 114.002, 53.0227, 38.1232, 8.65};
      //double yd [6] = {2497.7, 317.708, 141.514, 61.9, 41.24, 7.89291};
      
      //double yc [6] = {1966.09, 266.725, 101.439, 51.22, 38.1232, 9.21};
      //double yd [6] = {2604.77, 325.119, 140.402, 54.8, 41.24, 10.9478};
      double q2[5] = {1.5, 2.25, 2.75, 3.25, 4.75};
      double eq[5] = {0,0,0,0,0};
      double q2_2[6] = {1.5, 2.25, 2.75, 3.25, 4, 5.25};
      double eq_2[6] = {0,0,0,0,0,0};
      double ta[6];
      double eta[6];
      double ta2[6];
      double eta2[6];
      double ta3[6];
      double eta3[6];
      double ta4[5];
      double eta4[5];
      double ta5[6];
      double eta5[6];
      double eta6[6];
      
      double ta_ib[4];
      double eta_ib[4];
      double q2_ib[4] = {1.5, 2.25, 2.75, 4.5};
      double eq_ib[4] = {0,0,0,0};
//	  double cucon [5] = { 966.446/6, 170.258/6, 72.0879/6, 31.8191/6, 26.42/6};
//      double cumean [5] = { 0.766316, 0.7666626, 0.765916, 0.771877,0.76253};
      
      //double ccon [5] = { 0.00343964, 0.00373877, 0.00362233, 0.0034664,0.00347056};
      //double cmean [5] = { 0.765765, 0.766401, 0.764707, 0.766536,0.768765};
      
      //double dcon [5] = { 0.00347077, 0.003624, 0.00346489, 0.00338155,0.003222};
      //double dmean [5] = { 0.76475, 0.76502, 0.765078, 0.766053,0.766426};
      
      //double ccon [5] = { 4208.22, 673.941, 326.329, 145.753,117.703};
      //double cmean [5] = { 0.764721, 0.765424, 0.766017, 0.766229,0.76793};
      
      //double ccon [5] = { 669.758, 107.261, 51.9368, 23.1973, 18.7331};
      //double cmean [5] = { 0.764721, 0.765424, 0.766017, 0.766229,0.76793};
      
      //double dcon [5] = { 1815.87, 874.913, 475.313, 231.402, 149.866};
      //double dmean [5] = { 0.760601, 0.760077, 0.761425, 0.761183,0.761669};
      
    //  double dcon [5] = { 289.004, 139.247, 75.6484, 36.8287, 23.8519};
    //  double dmean [5] = { 0.760601, 0.760077, 0.761425, 0.761183,0.761669};
      
      
 //     double dcon [5] = { 1420.66, 289.965, 136.054, 63.518, 50.1955};
 //     double dmean [5] = { 0.763922, 0.762454, 0.763573, 0.76372,0.765136};
 
  //    double dcon [5] = { 1319.32, 175.398, 72.357, 32.3734, 25.2654};
  //    double dmean [5] = { 0.766345, 0.764046, 0.764443, 0.764745, 0.767623};
      
    //  double dcon [5] = { 56.7176, 32.9886, 17.6794, 9.68959, 11.7712};
    //  double dmean [5] = { 0.739508, 0.741838, 0.748383, 0.75439, 0.762814};
      int ce = 43464260;
      int de = 11608090;
      double sncon [6] = { 148.209, 24.1306, 11.2566, 5.48538, 3.9852,0};    //based on 3 with lc cut of 1.0
      double snmean [6] = { 0.765354, 0.760853, 0.762455, 0.762381,0.771106,0};
      
      double sn1con [6] = { 148.209, 24.1306, 11.2566, 5.48538, 3.9852,0};    //based on 3 with lc cut of 1.0
      double sn1mean [6] = { 0.765354, 0.760853, 0.762455, 0.762381,0.771106,0};
      
      double cucon [6] = { 145.718, 24.7849, 10.328, 5.04566, 3.70732,0};     // based on 4 with lc cut of 0.5
      double cumean [6] = { 0.766989, 0.767591, 0.766961, 0.769011,0.764005,0};
      
      double cu10con [6] = { 73.5858, 14.9832, 6.46788, 3.03546, 2.232, 0};
      double cu5con [6] = { 41.23358, 9.28734, 3.56906, 1.432, 1.12194, 0};
      
      double ccon [6] = { 1891.33, 313.193, 140.633, 63.6085, 32.2802, 14.6037};                    //aps c
      double cmean [6] = { 0.766488, 0.766401, 0.764707, 0.766536,0.7696, 0.76803};
      
     double dcon [6] = { 3405.05, 493.659, 219.487, 97.2822, 44.2606, 20.7415};                    //aps ld2 with 6 bins
     double dmean [6] = { 0.765431, 0.765536, 0.765078, 0.766051, 0.766135, 0.767852};
     
     double ddcon [6] = { 3405.05, 493.659, 219.487, 97.2822, 64.8886,0};                         //aps ld2 with 5 bins
     double ddmean [6] = { 0.765431, 0.765536, 0.765078, 0.766051, 0.766796,0};
     
  //   double cchcon [6] = { 637.701, 100.633, 48.7874, 21.4715, 14.0962, 4.24545};                 //cache c
  //    double cchmean [6] = { 0.766028, 0.766471, 0.766974, 0.767237,0.7696, 0.76571};
     
 //    double dchcon [6] = { 1420.66, 289.965, 136.054, 63.518, 41.3791, 9.06455};                   //ib and ob files
 //    double dchmean [6] = { 0.763922, 0.762454, 0.763573, 0.76372,0.765072, 0.764828};
     
     //double dchocon [6] = { 1163.01, 169.421, 72.1004, 32.3734, 20.6719, 4.94969};                 //only ob files
    // double dchomean [6] = { 0.76441, 0.763454, 0.764234, 0.764745,0.767267, 0.76603};
     
      double dicon [6] = { 264.511/de, 123.318/de, 65.4501/de, 31.5651/de, 25.0419/de, 0};                        //only ib files
     double dimean [6] = { 0.761009, 0.760497, 0.762143, 0.761927,0.762542, 0};
     
     double dbcon [6] = { 49.9254, 34.5609, 18.2289, 10.3275, 11.933, 0};                        //only bspot files
     double dbmean [6] = { 0.761009, 0.760497, 0.762143, 0.761927,0.762542, 0};
     double dilcon [6] = { 49.9254, 34.5609, 18.2289, 10.3275, 11.933, 0};
     
   //  double cxcmean [6] = { 0.761009, 0.760497, 0.752014, 0.750565,0.760022, 0.763643};
   //  double cxccon [6] = { 49.9254, 34.5609, 8.61356, 4.98051, 3.97888, 1.50378};
     
 //    double func(double xf, int i){
 //        double bwy = cxccon[i]*(0.15/(4*((cxcmean[i] - xf)*(cxcmean[i] - xf)) + (0.15*0.15)));     
 //    return bwy;
 //    }
 
               //only ob files
   double cxc_yields [6] = {46000.2, 9106.79, 4368.79, 1835.48, 1347.21, 328.836};
   double cxc_uyields [6] = {351.549, 159.63, 112.09, 74.1935, 62.0573, 32.2385};
   
   double ld_yields [6] = {121948, 17497, 7192.85, 3235.88, 2017.76, 426.963};
   double ld_uyields [6] = {580.718, 222.984, 147.447, 95.627, 76.2998, 34.5363};
   
   double ld2_yields [6] = {116140, 16604, 6863.53, 3086.41, 1930.6, 396.535};
   double ld2_uyields [6] = {566.409, 217.318, 143.953, 93.5444, 74.8588, 33.8865};
 
 //pass0v456
   double cxcp_yields [6] = {32208, 4893.4 , 2142.68 , 1240.45 , 865.131 , 104.004 };
   double cee = 21516150;
   double ldp_yields [6] = {64319.8 , 8549.12 , 4645.1 , 1911.34 , 1415.64 , 361.293 };
   double lee = 13249680;  
   
   double cxcpn_yields [6] = {60.8703, 56.1697 , 57.1813 , 56.7851 , 57.4128 , 62.7105 };
   double ldpn_yields [6] = {58.5845, 57.9553 , 58.4485 , 56.8284 , 54.5783 , 51.9748 };
   double cxcpnn_yields [6] = {73207.9, 10889.7 , 5133 , 2316.96 , 1564.32 , 425.176 };
   double ldpnn_yields [6] = {60507.2, 8810.56 , 3992.05 , 1665.95 , 1003.53 , 210.533 };
   //double cxcp_yields [6] = {32208, 4893.4 , 2142.68 , 1240.45 , 865.131 , 104.004 };
   double ld0v7_yields [6] = {41110.3, 6639.52, 3287.51, 1428.5, 934.702, 187.972};
   double cu0v7_yields [6] = {13268.6, 2511.84, 1125.64, 534.488, 360.462, 38.5781};
   double sn0v7_yields [6] = {15668.2, 2985.32, 1312.92, 679.581, 470.503, 66.6071};
   double ld0v72_yields [6] = {30240.5, 4714.07, 2556.2, 1377.07, 864.756, 173.899};
   double cu0v72_yields [6] = {9366.29, 1889.73, 937.684, 432.933, 247.585, 34.5971};
   double sn0v72_yields [6] = {11831, 2274.38, 1039.49, 648.933, 306.313, 18.8931};
   double ld0v73_yields [6] = {37181.7, 5946.84, 3007.59, 1344.07, 908.613, 173.938};
   double cu0v73_yields [6] = {11961.3, 2325.25, 1060.78, 497.718, 305.199, 39.6263};
   double sn0v73_yields [6] = {13870.4, 2717.26, 1233.39, 660.689, 428.192, 51.9666};
   double cxc0v73_yields [6] = {18189.8, 3092.56, 1444.3, 706.914, 459.012, 108.623};
   double ld0v74_yields [6] = {39337.9, 6140.38, 3041.53, 1356.88, 862.118, 168.687};
   double cu0v74_yields [6] = {11961.3, 2325.25, 1060.78, 497.718, 305.199, 39.6263};
   double sn0v74_yields [6] = {8991.61, 1812.28, 884.527, 482.567, 302.699, 40.449};
   double cxc0v74_yields [6] = {18189.8, 3092.56, 1444.3, 706.914, 459.012, 108.623};
   
   double ld0v75_yields [6] = {44639.6, 5826.2, 2655.02, 1064.06, 631.896, 110.587};
   double cu0v75_yields [6] = {9513.67, 1699.3, 719.257, 336.357, 200.708, 17.6941};
   double sn0v75_yields [6] = {11731.6, 1937.99, 849.92, 443.909, 289.909, 41.457};
   double cxc0v75_yields [6] = {19007.2, 2888.13, 1209.33, 533.43, 354.043, 63.2046};
   
   double cu07_e [6] = {30791600.0, 4621160.0, 2451230.0, 1355330.0, 777351.0, 981338.0};
   double sn07_e [6] = {43294894.0, 6585955.0, 3501706.0, 1956465.0, 1131479.0, 1460798.0};
   double ld07_e [6] = {45161979.0, 6838216.0, 3608882.0, 1991569.0, 1138923.0, 1437183.0};
   double cxc07_e [6] = {32361815.0, 4937021.0, 2624547.0, 1464742.0, 844143.0, 1088323.0};
   
   double ldc = 29.57;
   double ccc = 26.07;
   double snc = 158.52;
   
   /*
   if (target == Core::Target::CxC){
        thinkness_A = 0.4;
        density_A = 2.2;
        charge_A = 26.07;
    }
    if (target == Core::Target::Cu){
        thinkness_A = 0.0093;
        density_A = 8.96;
        charge_A = 158.52;
    }    
    if (target == Core::Target::Sn){
        thinkness_A = 0.0171;
        density_A = 7.31;
        charge_A = 158.52;
    }

    const double thinkness_LD2 = 5;
    const double density_LD2 = 0.164;
    const double charge_LD2 = 29.57;
   */
   
   double cxc0v76_yields [5] = {17140.6/ccc, 2625.63/ccc, 1111.99/ccc, 457.388/ccc, 412.923/ccc};
   double ld0v76_yields [5] = {40192.5/ldc, 5241.59/ldc, 2388.26/ldc, 934.584/ldc, 689.907/ldc};
   double cu0v76_yields [5] = {8566.15/snc, 1577.19/snc, 645.729/snc, 304.614/snc, 197.094/snc};
   
   
   double cxc0v77_yields [5] = {21411.9/ccc, 3293.23/ccc, 1353.86/ccc, 564.885/ccc, 517.77/ccc};
   double ld0v77_yields [5] = {51217/ldc, 6690.41/ldc, 3008.08/ldc, 1168.81/ldc, 828.395/ldc};
   double cu0v77_yields [5] = {10480.8/snc, 1894.98/snc, 788.38/snc, 370.17/snc, 246.2/snc};
   double sn0v77_yields [5] = {12981.9/snc, 2153.55/snc, 972.767/snc, 435.312/snc, 402.977/snc};
   
   double l7e = 133315400;
   double cu7e = 85025880;
   double sn7e = 110102700;
   
   //double cu0v11_ob_yields [5] = {19718.6, 3214.13, 1522.7, 719.155, 541.62};
   //double sn0v11_ob_yields [5] = {30879.9, 4651.89, 2108.62, 1099.52, 908.345};
   //double ld0v11_ob_yields [5] = {97207.7, 12724.8, 5593.14, 2399.47, 1807.75};
   //double cxc0v11_ob_yields [5] = {44515.1, 6402.48, 2759.63, 1262.86, 1131.46};
   
   //double cu0v11_ib_yields [5] = {328.582, 157.214, 65.3127, 39.2983, 17.9676};
   //double sn0v11_ib_yields [5] = {556.698, 312.982, 165.762, 31.9786, 51.5133};
   //double ld0v11_ib_yields [5] = {2470.04, 1233.09, 612.818, 269.875, 191.861};
   //double cxc0v11_ib_yields [5] = {1102.83, 666.887, 371.888, 146.63, 113.569};
   
   double cu0v11_ib_yields[4] = {187.582, 72.3252, 73.9344, 41.0114};    //with Calorimeter cut
   double sn0v11_ib_yields[4] = {277.802, 163.472, 59.1719, 53.536};    //with Calorimeter cut
   double ld0v11_ib_yields[4] = {1199.31, 644.467, 297.125, 258.592};    //with Calorimeter cut
   double cxc0v11_ib_yields[4] = {586.044, 357.761, 194.761, 183.845};    //with Calorimeter cut
   
   double cu0v11_ob_yields[4] = {10135.3, 1962.96, 953.221, 779.065};   //with Calorimeter cut
   double sn0v11_ob_yields[4] = {12532.4, 2243.06, 1104.24, 1055.32};    //with Calorimeter cut
   double ld0v11_ob_yields[4] = {42209.2, 6549.08, 2967.73, 2276.7};    //with Calorimeter cut
   double cxc0v11_ob_yields[4] = {20521.3, 3424.07, 1636.56, 1415.24};    //with Calorimeter cut
   
   //double cu0v11_ob_yields[5] = {10135.3, 1962.96, 953.221, 427.509, 290.244};   //with Calorimeter cut
   //double sn0v11_ob_yields[5] = {12532.4, 2243.06, 1104.24, 569.271, 496.545};    //with Calorimeter cut
   //double ld0v11_ob_yields[5] = {42209.2, 6549.08, 2967.73, 1336.75, 948.481};    //with Calorimeter cut
   //double cxc0v11_ob_yields[5] = {20521.3, 3424.07, 1636.56, 732.548, 695.664};    //with Calorimeter cut
   
   double ldp1_ob_yields[5] = {78526.4, 12742.4, 5464.44, 2283.9, 2005.01};    //with Calorimeter cut
   double cxcp1_ob_yields[5] = {25411.3, 4399.93, 2030.53, 820.623, 839.247};    //with Calorimeter cut
   
   double cu1_ob_yields[6] = {36423.3, 6466.49, 3097.23, 1322.53, 1127.71, 250.55};   //with Calorimeter cut
   double sn1_ob_yields[6] = {44144.7, 7812.47, 3506.69, 1782.51, 1219.8, 251.392};    //with Calorimeter cut
   double ld1_ob_yields[6] = {164149, 26376.6, 11399.2, 5085.1, 3356.22, 757.455};    //with Calorimeter cut
   double cxc1_ob_yields[6] = {79359.2, 13675.5, 6318.63, 2765.32, 1934.3, 492.539};    //with Calorimeter cut
   
   
   /*
   double cuto = 157.905 * pow(10,6);
   double cunto = 4.70315 * pow(10,9);
   
   double snto = 157.905 * pow(10,6);
   double snnto = 4.70315 * pow(10,9);
   
   double ldto = 27.574 * pow(10,6);
   double ldnto = 8.74571 * pow(10,8);
   
   double cxcto = 26.0234 * pow(10,6);
   double cxcnto = 9.13238 * pow(10,8);
   
   double cuti = 11.6958 * pow(10,6);
   double cunti = 3.24905 * pow(10,8);
   
   double snti = 11.6958 * pow(10,6);
   double snnti = 3.24905 * pow(10,8);
   
   double ldti = 1.76664 * pow(10,6);
   double ldnti = 3.80188 * pow(10,7);
   
   double cxcti = 1.65141 * pow(10,6);
   double cxcnti = 4.00005 * pow(10,7);
   */
   
   double cunti  = 831845.0;
   double snnti  = 831845.0;
   double cxcnti = 135249.0;
   double ldnti  = 132491.0;
   
   double cunto  = 6852860.0;
   double snnto  = 6852860.0;
   double cxcnto =  901523.0;
   double ldnto  = 1036990.0;
   
   double ldob1 = 389439;
   double cxcob1 = 209752;
   
   double ldp1 = 4331578;
   double cxcp1 = 3304196;
   double cup1 = 24934751;
   // cusn inbending Total Charge = 1.16958e+07
   //           normalized charge = 831845


   //cxc inbending  Total Charge = 1.65141e+06
   //          normalized charge = 135249
   
   
   
   //LD2 inbending  Total Charge = 1.76664e+06
   //          normalized charge = 132491



   //CuSn outbending  Total Charge = 1.57905e+08
   //            normalized charge = 6.85286e+06
   
   //LD2 outbending  Total Charge = 2.7574e+07
   //           normalized charge = 1.03699e+06
   
   //CxC outbending  Total Charge = 2.60234e+07
   //           normalized charge = 901523
   
   
	  for(int i = 0; i < 6; i++){
         double yc = 0.0;
         double yc2 = 0.0;
         double ycu = 0.0;
         double ysn = 0.0;
         double yd = 0.0;
         double ydd = 0.0;
         double yd2 = 0.0;
         double xs = 0.3;
         /*
         for(int j = 0; j < 10000; j++){
            //yc = yc + cxccon[i]*(0.15/(4*((cxcmean[i] - xs)*(cxcmean[i] - xs)) + (0.15*0.15)));
            yc2 = yc2 + cchcon[i]*(0.15/(4*((cchmean[i] - xs)*(cchmean[i] - xs)) + (0.15*0.15)));
            if(i < 5){
            ycu = ycu + cucon[i]*(0.15/(4*((cumean[i] - xs)*(cumean[i] - xs)) + (0.15*0.15)));
            ysn = ysn + sncon[i]*(0.15/(4*((snmean[i] - xs)*(snmean[i] - xs)) + (0.15*0.15)));
            }
            yd = yd + dcon[i]*(0.15/(4*((dmean[i] - xs)*(dmean[i] - xs)) + (0.15*0.15)));                        //dmean
            ydd = ydd + ddcon[i]*(0.15/(4*((ddmean[i] - xs)*(ddmean[i] - xs)) + (0.15*0.15)));                   //ddmean
            yd2 = yd2 + dchocon[i]*(0.15/(4*((dchomean[i] - xs)*(dchomean[i] - xs)) + (0.15*0.15)));             //dchomean
            xs = xs + 0.00011;
            //std::cout << xs << std::endl;
          }
          
         */
     
    // double ybwy = fbw->Integral(0.3, 1.4);
    //      std::cout << "ybwy = " << ybwy/0.01 << std::endl;
          //std::cout << "cxc bin " << yc << std::endl;
         int ycc = 1.9;
         int ydc = 4.12;
         ta[i] = ((cxc1_ob_yields[i] * 0.82 ) / (ld1_ob_yields[i] * 0.80)) * ((ldp1)/(cxcp1));   //* ((27.5931)/(26.0234));  //(ld07_e[i]/cxc07_e[i]);       //0.0833 for cu   ccon[i] dcon[i]
         ta2[i] = ((cu1_ob_yields[i] * 0.82) / (ld1_ob_yields[i] * 0.0813 )) * ((ldp1)/(cup1));  //* ((27.5931)/(157.914));//(ld07_e[i]/(cu07_e[i] + sn07_e[i]));    // * (cu7e/l7e);       //0.0833 for cu   * (l7e/cu7e)
         ta3[i] = ((sn1_ob_yields[i] * 0.82 ) / (ld1_ob_yields[i] * 0.117 )) * ((ldp1)/(cup1));  //* ((27.5931)/(157.914));//(ld07_e[i]/(cu07_e[i] + sn07_e[i]));      //* (sn7e/l7e);          //ddcon
         ta4[i] = ((yccon[i]) * 0.82) / ((ydcon[i]) * 0.88);        //dchocon
         //std::cout << yc << std::endl;
         
         ta5[i] = ((cxcpnn_yields[i]) * 0.82) / ((ldpnn_yields[i]) * 0.88) ;//* (lee/cee); 
         double cxc_per = (cxc_uyields[i] / cxc_yields[i]) * 100.0;
         double ld2_per = (ld2_uyields[i] / ld2_yields[i]) * 100.0;
         eta5[i] = ta5[i] * (sqrt(pow(cxc_per,2) + pow(ld2_per,2)) / 100);
         
         
         eta[i] = (ta[i] * sqrt((1.0/cxc1_ob_yields[i]) + (1.0/ld1_ob_yields[i]))) * 1.5 ;
         eta2[i] = (ta2[i] * sqrt((1.0/cu1_ob_yields[i]) + (1.0/ld1_ob_yields[i]))) * 1.5 ;
         eta3[i] = (ta3[i] * sqrt((1.0/sn1_ob_yields[i]) + (1.0/ld1_ob_yields[i]))) * 1.5 ;
         eta4[i] = (ta4[i] * sqrt((1/(yccon[i]))+ (1/(ydcon[i])))) / 2;
         //nt->AddPoint(q2[i],ta);
         eta6[i] = ta5[i] * (sqrt((1/cxcpnn_yields[i]) + (1/ldpnn_yields[i])) );
         std::cout << "---------------------------------------------------------------------------------------------" << std::endl;
         std::cout << "bin" << i+1 << ": cxc = " << cxc_yields[i] << " | sqrt(cxc) = " << sqrt(cxc_yields[i]) << " | uncertanty = " << cxc_uyields[i] << std::endl;
         std::cout << "bin" << i+1 << ": LD2 = " << ld2_yields[i] << " | sqrt(LD2) = " << sqrt(ld2_yields[i]) << " | uncertanty = " << ld2_uyields[i] << std::endl;
         std::cout << "NT = " << ta5[i] << " | un"  << eta5[i] << " | LE_un = " << ta5[i] * (sqrt((1/cxc_yields[i]) + (1/ld2_yields[i])) ) << std::endl;
         
         
         std::cout << ta[i] << std::endl;
         //std::cout << (1/ccon[i])+ (1/dcon[i]) << std::endl;
         //std::cout << sqrt((1/ccon[i])+ (1/dcon[i])) << std::endl;
      }
      double ta_aps[6] = {0.52, 0.5801, 0.598, 0.6095, 0.673, 0.73};
      double eta_aps[6] = {0.006, 0.012, 0.024, 0.045, 0.056, 0.13};
      auto nt = new TGraphErrors(6,q2_2,ta,eq_2,eta);
      nt->SetLineColor(kGreen-1);
      //nt->SetLineWidth(0)
      nt->SetMarkerColor(kGreen-1);
      nt->SetMarkerStyle(21);
      
      auto nt2 = new TGraphErrors(6,q2_2,ta2,eq_2,eta2);
      nt2->SetLineColor(kGreen-1);
      //nt2->SetLineWidth(0);
      nt2->SetMarkerColor(kOrange+2);
      nt2->SetMarkerStyle(21);
      
      auto nt3 = new TGraphErrors(6,q2_2,ta3,eq_2,eta3);
      nt3->SetLineColor(kGreen-1);
      //nt3->SetLineWidth(0);
      nt3->SetMarkerColor(kMagenta+2);
      nt3->SetMarkerStyle(21);
      
      auto nt4 = new TGraphErrors(5,q2,ta5,eq,eta5);
      nt4->SetLineColor(kGreen-1);
      //nt4->SetLineWidth(0);
      nt4->SetMarkerColor(kGreen-1);
      nt4->SetMarkerStyle(21);
      
	
		
        
        nt->SetTitle("Carbon;Q^{2} [GeV^{2}];T_{A}");
        nt->SetMinimum(0.2);
   	nt->SetMaximum(1.1);
   	
   	nt->Draw("AP");
   	//nt->Write();
   	
   	nt2->SetTitle("Copper;Q^{2} [GeV^{2}];T_{A}");
   	nt2->SetMinimum(0.2);
   	nt2->SetMaximum(1.1);
   	nt2->Draw("AP");
   	//nt2->Write();
   	
   	nt3->SetTitle("Tin;Q^{2} [GeV^{2}];T_{A}");
   	nt3->SetMinimum(0.2);
   	nt3->SetMaximum(1.1);
   	nt3->Draw("AP");
   	//nt3->Write();
   	
   	nt4->SetTitle("Cache Carbon Nuclear Transparency;Q^{2} [GeV^{2}];T_{A}");
   	nt4->SetMinimum(0.0);
   	nt4->SetMaximum(1.1);
   	nt4->Draw("AP");
   	//nt4->Write();
   	
   	
   	auto mg = new TMultiGraph;
   	mg->SetTitle("Nuclear Transparency;Q^{2} [GeV^{2}]; T_{A}");
   	mg->Add(nt);
   	mg->Add(nt2);
   	mg->Add(nt3);
   	//mg->Add(nt4);
   	mg->SetMinimum(0.0);
   	mg->SetMaximum(1.2);
   	TAxis *axis = mg->GetXaxis();
   	axis->SetLimits(1,6);
   	mg->Draw("AP");
   	mg->Write();
   	
   	auto c2le = new TLegend(0.1,0.75,0.25,0.9);
  // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   c2le->AddEntry(nt,"Carbon","p");
   c2le->AddEntry(nt2,"Copper","p");
   c2le->AddEntry(nt3,"Tin","p");
   	
   	TCanvas *c2 = new TCanvas("c2","",0,0,640,480);
   	//nt->Draw();
   	//nt2->Draw("SAME");
   	//nt3->Draw("SAME");
   	//c1->BuildLedgend();
   	mg->Draw("AP");
   	c2le->Draw();
   	c2->Draw();
   	//c2->Write();
        
        
}


