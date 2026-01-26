#include <cstdlib>
#include <iostream>
#include <string>


void graph_runs(){


    	
	TFile f("graph_cache_runs_comparison_3.root", "update");
	
	
      
      double ldob_run[5] = {18528, 18529, 18530, 18532, 18533};
      double ldib_run[4] = { 18324, 18331, 18335, 18336}; //{18312, 18324, 18325, 18326, 18330, 18331, 18335, 18336}; 324 , 331, 335, 336
      double cxcob_run[5] = {18451, 18452, 18453, 18454, 18459};
      
      double ldob_file[5] = {255, 343, 344, 258, 359};
      double ldib_file[4] = { 65, 242, 257, 171}; //{224, 65, 164, 181, 158, 242, 257, 171};
      double cxcob_file[5] = {364, 367, 381, 358, 360};
      
      double ldob_e[5] = {68302395/255, 68283413/343, 68222517/344, 68168033/258, 70473349/359};
      double ldib_e[4] =  {3676025/65,  14690120/242, 14711142/257, 7230476/171}; //{11844766/224, 3676025/65, 6873782/164, 5465134/181, 9741124/158, 14690120/242, 14711142/257, 7230476/171};
      double cxcob_e[5] = {62569699/364, 63432912/367, 65799302/381, 63163949/358, 62754364/360};
      
      double ldob_e_pip[5] = {20851172/255, 21090036/343, 21044711/344, 20915433/258, 21554806/359};
      double ldib_e_pip[4] = {1789154/65, 6464115/242, 7162946/257, 3524934/171}; //{5669219/224, 1789154/65, 3480564/164, 2858217/181, 4679576/158, 6464115/242, 7162946/257, 3524934/171};
      double cxcob_e_pip[5] = {18551257/364, 18810430/367, 19475729/381, 18558476/358, 18521082/360};
      
      double ldob_e_pim[5] = {24649835/255, 25045870/343, 24981986/344, 24838776/258, 25636357/359};
      double ldib_e_pim[4] = { 959957/65, 3138927/242, 3841196/257, 1893108/171}; //{3053802/224, 959957/65, 1883295/164, 1584606/181, 2502066/158, 3138927/242, 3841196/257, 1893108/171};
      double cxcob_e_pim[5] = {19774858/364, 20038695/367, 20693999/381, 19871823/358, 19661376/360};
      
      double ldob_e_pip_pim[5] = {8856480/255, 9107285/343, 9075802/344, 8975841/258, 9238326/395};
      double ldib_e_pip_pim[4] = { 536061/65, 1582132/242, 2144575/257, 1057504/171}; //{1666889/224, 536061/65, 1107113/164, 971595/181, 1377180/158, 1582132/242, 2144575/257, 1057504/171};
      double cxcob_e_pip_pim[5] = {7548047/364, 7642735/367, 7870472/381, 7497963/358, 7450893/360};
      
      
      double eq[5] = {0,0,0,0,0};
      double eta[5] = {0,0,0,0,0};
      double eq2[4] = {0,0,0,0};
      double eta2[4] = {0,0,0,0};

      
      
      auto ldobe = new TGraphErrors(5,ldob_run, ldob_e,eq,eta);
      ldobe->SetLineColor(kGreen-1);
      //nt->SetLineWidth(0);
      ldobe->SetMarkerColor(kBlue);
      ldobe->SetMarkerStyle(21);
      
      auto ldibe = new TGraphErrors(4,ldib_run, ldib_e,eq2,eta2);
      ldibe->SetLineColor(kGreen-1);
      //nt2->SetLineWidth(0);
      ldibe->SetMarkerColor(kCyan+1);
      ldibe->SetMarkerStyle(21);
      
      auto cxcobe = new TGraphErrors(5,cxcob_run, cxcob_e,eq,eta);
      cxcobe->SetLineColor(kGreen-1);
      //nt3->SetLineWidth(0);
      cxcobe->SetMarkerColor(kBlack);
      cxcobe->SetMarkerStyle(21);
      
      auto ldobep = new TGraphErrors(5,ldob_run, ldob_e_pip,eq,eta);
      ldobep->SetLineColor(kGreen-1);
      //nt->SetLineWidth(0);
      ldobep->SetMarkerColor(kBlue);
      ldobep->SetMarkerStyle(21);
      
      auto ldibep = new TGraphErrors(4,ldib_run, ldib_e_pip,eq2,eta2);
      ldibep->SetLineColor(kGreen-1);
      //nt2->SetLineWidth(0);
      ldibep->SetMarkerColor(kCyan+1);
      ldibep->SetMarkerStyle(21);
      
      auto cxcobep = new TGraphErrors(5,cxcob_run, cxcob_e_pip,eq,eta);
      cxcobep->SetLineColor(kGreen-1);
      //nt3->SetLineWidth(0);
      cxcobep->SetMarkerColor(kBlack);
      cxcobep->SetMarkerStyle(21);
      
      auto ldobem = new TGraphErrors(5,ldob_run, ldob_e_pim,eq,eta);
      ldobem->SetLineColor(kGreen-1);
      //nt->SetLineWidth(0);
      ldobem->SetMarkerColor(kBlue);
      ldobem->SetMarkerStyle(21);
      
      auto ldibem = new TGraphErrors(4,ldib_run, ldib_e_pim,eq2,eta2);
      ldibem->SetLineColor(kGreen-1);
      //nt2->SetLineWidth(0);
      ldibem->SetMarkerColor(kCyan+1);
      ldibem->SetMarkerStyle(21);
     
      auto cxcobem = new TGraphErrors(5,cxcob_run, cxcob_e_pim,eq,eta);
      cxcobem->SetLineColor(kGreen-1);
      //nt3->SetLineWidth(0);
      cxcobem->SetMarkerColor(kBlack);
      cxcobem->SetMarkerStyle(21);
      
      auto ldobepm = new TGraphErrors(5,ldob_run, ldob_e_pip_pim,eq,eta);
      ldobepm->SetLineColor(kGreen-1);
      //nt->SetLineWidth(0);
      ldobepm->SetMarkerColor(kBlue);
      ldobepm->SetMarkerStyle(21);
      
      auto ldibepm = new TGraphErrors(4,ldib_run, ldib_e_pip_pim,eq2,eta2);
      ldibepm->SetLineColor(kGreen-1);
      //nt2->SetLineWidth(0);
      ldibepm->SetMarkerColor(kCyan+1);
      ldibepm->SetMarkerStyle(21);
      
      auto cxcobepm = new TGraphErrors(5,cxcob_run, cxcob_e_pip_pim,eq,eta);
      cxcobepm->SetLineColor(kGreen-1);
      //nt3->SetLineWidth(0);
      cxcobepm->SetMarkerColor(kBlack);
      cxcobepm->SetMarkerStyle(21);
     
      
     
	
		
        
        ldobe->SetTitle(" ;run #;e events");
        ldobe->SetMinimum(000000);
        ldobe->SetMaximum(200000);
   	ldobe->Draw("AP");
   	ldobe->Write();
   	
   	ldibe->SetTitle(" ;run #;e events");
   	ldibe->SetMinimum(000000);
   	ldibe->SetMaximum(200000);
   	ldibe->Draw("AP");
   	ldibe->Write();
   	
   	cxcobe->SetTitle(" ;run #;e events");
   	cxcobe->SetMinimum(000000);
   	cxcobe->SetMaximum(200000);
   	cxcobe->Draw("AP");
   	cxcobe->Write();
   	
   	
   	ldobep->SetTitle(" ;run #;e events");
        ldobep->SetMinimum(0000);
        ldobep->SetMaximum(100000);
   	ldobep->Draw("AP");
   	ldobep->Write();
   	
   	ldibep->SetTitle(" ;run #;e events");
   	ldibep->SetMinimum(0000);
   	ldibep->SetMaximum(100000);
   	ldibep->Draw("AP");
   	ldibep->Write();
   	
   	cxcobep->SetTitle(" ;run #;e events");
   	cxcobep->SetMinimum(0000);
   	cxcobep->SetMaximum(100000);
   	cxcobep->Draw("AP");
   	cxcobep->Write();
   	
   	
   	ldobem->SetTitle(" ;run #;e events");
        ldobem->SetMinimum(0000);
        ldobem->SetMaximum(100000);
   	ldobem->Draw("AP");
   	ldobem->Write();
   	
   	ldibem->SetTitle(" ;run #;e events");
   	ldibem->SetMinimum(0000);
   	ldibem->SetMaximum(100000);
   	ldibem->Draw("AP");
   	ldibem->Write();
   	
   	cxcobem->SetTitle(" ;run #;e events");
   	cxcobem->SetMinimum(0000);
   	cxcobem->SetMaximum(100000);
   	cxcobem->Draw("AP");
   	cxcobem->Write();
   	
   	
   	ldobepm->SetTitle(" ;run #;e events");
        ldobepm->SetMinimum(0000);
        ldobepm->SetMaximum(100000);
   	ldobepm->Draw("AP");
   	ldobepm->Write();
   	
   	ldibepm->SetTitle(" ;run #;e events");
   	ldibepm->SetMinimum(0000);
   	ldibepm->SetMaximum(100000);
   	ldibepm->Draw("AP");
   	ldibepm->Write();
   	
   	cxcobepm->SetTitle(" ;run #;e events");
   	cxcobepm->SetMinimum(0);
   	cxcobepm->SetMaximum(100000);
   	cxcobepm->Draw("AP");
   	cxcobepm->Write();
   	
   	
   	
   	auto mg = new TMultiGraph;
   	mg->SetTitle(" ;run num; e events / per file");
   	mg->Add(ldobe);
   	mg->Add(ldibe);
   //	mg->Add(cxcobe);
   	mg->SetMinimum(000000);
   	mg->SetMaximum(200000);
   	mg->Draw("AP");
   	mg->Write();
   	
   	auto mg2 = new TMultiGraph;
   	mg2->SetTitle(" ;run num; e pip events / per file");
   	mg2->Add(ldobep);
   	mg2->Add(ldibep);
   //	mg2->Add(cxcobep);
   	mg2->SetMinimum(0000);
   	mg2->SetMaximum(100000);
   	mg2->Draw("AP");
   	mg2->Write();
   	
   	auto mg3 = new TMultiGraph;
   	mg3->SetTitle(" ;run num; e pim events / per file");
   	mg3->Add(ldobem);
   	mg3->Add(ldibem);
   //	mg3->Add(cxcobem);
   	mg3->SetMinimum(0000);
   	mg3->SetMaximum(100000);
   	mg3->Draw("AP");
   	mg3->Write();
   	
   	auto mg4 = new TMultiGraph;
   	mg4->SetTitle(" ;run num; e pip pim events / per file");
   	mg4->Add(ldobepm);
   	mg4->Add(ldibepm);
  // 	mg4->Add(cxcobepm);
   	mg4->SetMinimum(0000);
   	mg4->SetMaximum(100000);
   	mg4->Draw("AP");
   	mg4->Write();
        
        
}


