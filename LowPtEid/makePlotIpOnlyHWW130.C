{

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  TFile *_file0 = TFile::Open("testIP_baby_EGMon2010B.root");
  //TFile *_file0 = TFile::Open("baby_EGMonitor.root");
  TFile *_file1 = TFile::Open("testIP_baby_GluGluToHToWWTo2L2Nu_M-130_7TeV_LOPU.root");

  TTree* tree0 = (TTree*) _file0->Get("tree");
  TTree* tree1 = (TTree*) _file1->Get("tree");

  TCanvas *c1 = new TCanvas();
  c1->cd();

  TCut bar = "abs(el_etaSC)<1.479";
  TCut end = "abs(el_etaSC)>1.479";
  TCut cutNoPt = "met<20&&el_Mt<20";
  TCut base0 = "el_MZ<0&&met<20&&el_Mt<20&&el10_sw&&el_pt<20";
  TCut base1 = "el_pt<20";
  
  TCut v0id("el_VBTF80 && (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95&&abs(el_dPhiIn*el_q)<0.006))") ;
  TCut eta("abs(el_etaSC)<2.2");
  TCut mva("el_mva>0.4");
  TCut newconv("abs(el_newconv_dist)>0.05||abs(el_newconv_dcot)>0.02||abs(el_newconv_dmh)>1");
  TCut iso("20*el_relIso/el_pt<0.1");
  //iso = "(el_CICt&2)>0";
  TCut ip("abs(el_d0corr)<0.02");
  TCut conv("(abs(el_conv_dist)>0.02||abs(el_conv_dcot)>0.02) && el_innerlayer39X==0");
 
  TCut vbtf90_id  ("el_VBTF90");
  TCut vbtf85_id  ("el_VBTF90 && (abs(el_etaSC)<1.479&&abs(el_dPhiIn)<0.06&&abs(el_dEtaIn)<0.006&&el_hOverE<0.04 || abs(el_etaSC)>1.479&&abs(el_dPhiIn)<0.04&&abs(el_dEtaIn)<0.007&&el_hOverE<0.025)");
  TCut vbtf80_id ("el_VBTF80");
  TCut vbtf70_id ("el_VBTF70");
  TCut vbtf60_id ("el_VBTF70 && (abs(el_etaSC)<1.479&&abs(el_dPhiIn)<0.025 || abs(el_etaSC)>1.479)");
  TCut vbtf80p_id("el_VBTF80 && (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95))") ;
  TCut vbtf70p_id("el_VBTF70 && (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95))") ;

  TCut vbtf80_nohoeend_id ("el_VBTF80 & abs(el_etaSC)<1.479 || abs(el_etaSC)>1.479&&el_sigmaIEtaIEta<0.03&&abs(el_dPhiIn)<0.03&&abs(el_dEtaIn)<0.007");
  TCut vbtf80p_nohoeend_id ("(el_VBTF80 & abs(el_etaSC)<1.479 || abs(el_etaSC)>1.479&&el_sigmaIEtaIEta<0.03&&abs(el_dPhiIn)<0.03&&abs(el_dEtaIn)<0.007) && (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95))");

  TCut vbtf70_nohoeend_id ("el_VBTF70 & abs(el_etaSC)<1.479 || abs(el_etaSC)>1.479&&el_sigmaIEtaIEta<0.03&&abs(el_dPhiIn)<0.02&&abs(el_dEtaIn)<0.005");
  TCut vbtf70p_nohoeend_id ("(el_VBTF70 & abs(el_etaSC)<1.479 || abs(el_etaSC)>1.479&&el_sigmaIEtaIEta<0.03&&abs(el_dPhiIn)<0.02&&abs(el_dEtaIn)<0.005) && (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95))");

  TCut cic_ht1_id  ("(el_CICht1&1)>0");
  TCut cic_ht1_iso ("(el_CICht1&2)>0");
  TCut cic_ht1_conv("(el_CICht1&4)>0");
  TCut cic_ht1_ip  ("(el_CICht1&8)>0");
  TCut cic_ht2_id  ("(el_CICht2&1)>0");
  TCut cic_ht2_iso ("(el_CICht2&2)>0");
  TCut cic_ht2_conv("(el_CICht2&4)>0");
  TCut cic_ht2_ip  ("(el_CICht2&8)>0");
  TCut cic_ht3_id  ("(el_CICht3&1)>0");
  TCut cic_ht3_iso ("(el_CICht3&2)>0");
  TCut cic_ht3_conv("(el_CICht3&4)>0");
  TCut cic_ht3_ip  ("(el_CICht3&8)>0");
  TCut cic_ht4_id  ("(el_CICht4&1)>0");
  TCut cic_ht4_iso ("(el_CICht4&2)>0");
  TCut cic_ht4_conv("(el_CICht4&4)>0");
  TCut cic_ht4_ip  ("(el_CICht4&8)>0");
  TCut cic_st_id  ("(el_CICst&1)>0");
  TCut cic_st_iso ("(el_CICst&2)>0");
  TCut cic_st_conv("(el_CICst&4)>0");
  TCut cic_st_ip  ("(el_CICst&8)>0");
  TCut cic_t_id  ("(el_CICt&1)>0");
  TCut cic_t_iso ("(el_CICt&2)>0");
  TCut cic_t_conv("(el_CICt&4)>0");
  TCut cic_t_ip  ("(el_CICt&8)>0");

  TCut lik_t_id ("abs(el_etaSC)<1.479&&el_fbrem<0.15&&el_LL>0.820 || abs(el_etaSC)<1.479&&el_fbrem>0.15&&el_LL>0.925 || abs(el_etaSC)>1.479&&el_fbrem<0.15&&el_LL>0.930 || abs(el_etaSC)>1.479&&el_fbrem>0.15&&el_LL>0.979");
  TCut pfmva3_id ("el_mva>0.3");
  TCut pfmva5_id ("el_mva>0.5");
  TCut pfmva6_id ("el_mva>0.6");
  TCut pfmva7_id ("el_mva>0.7");
  TCut pfmva75_id ("el_mva>0.75");
  TCut pfmva8_id ("el_mva>0.8");

  base0=base0+vbtf80_id+conv+iso;
  base1=base1+vbtf80_id+conv+iso;

  float ptbins[] = {15.,20.};
  float etabins[] = {0,2.5};
  TH2F* pt_vs_eta = new TH2F("pt_vs_eta","pt_vs_eta",1,etabins,1,ptbins);

  TGraph gr_iso(2);
  gr_iso.SetTitle("");
  gr_iso.SetPoint(0, 0, 0);
  gr_iso.SetPoint(1, 1.5, 3.0);
  gr_iso.SetMarkerStyle(1);

  TCut ip_d02d_0008("ip_d02d_0008","abs(el_d0pv2d)<0.008");
  TCut ip_d02d_001("ip_d02d_001","abs(el_d0pv2d)<0.01");
  TCut ip_d02d_0015("ip_d02d_0015","abs(el_d0pv2d)<0.015");
  TCut ip_d02d_002("ip_d02d_002","abs(el_d0pv2d)<0.02");
  TCut ip_d02d_005("ip_d02d_005","abs(el_d0pv2d)<0.05");

  TCut ip_d03d_002("ip_d03d_002","abs(el_d0pv3d)<0.02");
  TCut ip_d03d_003("ip_d03d_003","abs(el_d0pv3d)<0.03");
  TCut ip_d03d_005("ip_d03d_005","abs(el_d0pv3d)<0.05");
  TCut ip_d03d_010("ip_d03d_010","abs(el_d0pv3d)<0.10");
  TCut ip_d03d_020("ip_d03d_020","abs(el_d0pv3d)<0.20");

  TCut ip_d02ds_12("ip_d02ds_12","abs(el_d0pv2d/el_d0pv2dErr)<1.2");
  TCut ip_d02ds_14("ip_d02ds_14","abs(el_d0pv2d/el_d0pv2dErr)<1.4");
  TCut ip_d02ds_17("ip_d02ds_17","abs(el_d0pv2d/el_d0pv2dErr)<1.7");
  TCut ip_d02ds_2("ip_d02ds_2","abs(el_d0pv2d/el_d0pv2dErr)<2");
  TCut ip_d02ds_5("ip_d02ds_5","abs(el_d0pv2d/el_d0pv2dErr)<5");

  TCut ip_d03ds_16("ip_d03ds_16","abs(el_d0pv3d/el_d0pv3dErr)<1.6");
  TCut ip_d03ds_18("ip_d03ds_18","abs(el_d0pv3d/el_d0pv3dErr)<1.8");
  TCut ip_d03ds_2("ip_d03ds_2","abs(el_d0pv3d/el_d0pv3dErr)<2");
  TCut ip_d03ds_3("ip_d03ds_3","abs(el_d0pv3d/el_d0pv3dErr)<3");
  TCut ip_d03ds_9("ip_d03ds_9","abs(el_d0pv3d/el_d0pv3dErr)<9");

  TCut cuts[] = {ip_d02d_0008,ip_d02d_001,ip_d02d_0015,ip_d02d_002,ip_d02d_005,
		 ip_d03d_002,ip_d03d_003,ip_d03d_005,ip_d03d_010,ip_d03d_020,
		 ip_d02ds_12,ip_d02ds_14,ip_d02ds_17,ip_d02ds_2,ip_d02ds_5,
  		 ip_d03ds_16,ip_d03ds_18,ip_d03ds_2,ip_d03ds_3,ip_d03ds_9};
  int ncuts = sizeof(cuts)/sizeof(TCut);


  tree0->Draw("el_pt:abs(el_etaSC)>>pt_vs_eta",base0+ip_d02d_002);
  float all0 = pt_vs_eta->GetBinContent(1,1);
  tree1->Draw("el_pt:abs(el_etaSC)>>pt_vs_eta",base1+ip_d02d_002);
  float all1 = pt_vs_eta->GetBinContent(1,1);

  TLegend *l1 = new TLegend(0.2, 0.7, 0.6, 0.9);
  l1->SetNColumns(2);
  l1->SetFillColor(kWhite);

  gr_iso.Draw("AP");
  gr_iso.GetXaxis()->SetRangeUser(0.85,1.05);
  gr_iso.GetYaxis()->SetRangeUser(0.6,1.3);
  gr_iso.GetXaxis()->SetTitle("Electron IP rel eff (HWW130)");
  gr_iso.GetYaxis()->SetTitle("Electron IP rel eff (EGMon2010B)");

  TGraph gr_d02d(ncuts/4);
  gr_d02d.SetMarkerStyle(20); 
  gr_d02d.SetMarkerColor(kBlue);
  gr_d02d.SetLineColor(kBlue);

  TGraph gr_d03d(ncuts/4);
  gr_d03d.SetMarkerStyle(20); 
  gr_d03d.SetMarkerColor(kRed);
  gr_d03d.SetLineColor(kRed);

  TGraph gr_d02ds(ncuts/4);
  gr_d02ds.SetMarkerStyle(21); 
  gr_d02ds.SetMarkerColor(kBlue);
  gr_d02ds.SetLineColor(kBlue);

  TGraph gr_d03ds(ncuts/4);
  gr_d03ds.SetMarkerStyle(21); 
  gr_d03ds.SetMarkerColor(kRed);
  gr_d03ds.SetLineColor(kRed);

  int n_d02d=0,n_d03d=0,n_d02ds=0,n_d03ds=0;
  for (unsigned int ic=0;ic<ncuts;++ic){
    tree0->Draw("el_pt:abs(el_etaSC)>>pt_vs_eta",base0+cuts[ic],"goff");
    float cut0 = pt_vs_eta->GetBinContent(1,1);
    tree1->Draw("el_pt:abs(el_etaSC)>>pt_vs_eta",base1+cuts[ic],"goff");
    float cut1 = pt_vs_eta->GetBinContent(1,1);
    if (TString(cuts[ic].GetName()).Contains("ip_d02d_")) gr_d02d.SetPoint(n_d02d++, cut1/all1, cut0/all0);
    if (TString(cuts[ic].GetName()).Contains("ip_d03d_")) gr_d03d.SetPoint(n_d03d++, cut1/all1, cut0/all0);
    if (TString(cuts[ic].GetName()).Contains("ip_d02ds_")) gr_d02ds.SetPoint(n_d02ds++, cut1/all1, cut0/all0);
    if (TString(cuts[ic].GetName()).Contains("ip_d03ds_")) gr_d03ds.SetPoint(n_d03ds++, cut1/all1, cut0/all0);
    cout << cuts[ic].GetName() << " hww: " << Form("%4.2f",cut1/all1) << " - eg2010a: " << Form("%4.2f",cut0/all0) << endl;

  }

  l1->AddEntry(&gr_d02d, "d0 2D", "p");
  gr_d02d.Draw("PL");
  l1->AddEntry(&gr_d03d, "d0 3D", "p");
  gr_d03d.Draw("PL");
  l1->AddEntry(&gr_d02ds, "d0 2Ds", "p");
  gr_d02ds.Draw("PL");
  l1->AddEntry(&gr_d03ds, "d0 3Ds", "p");
  gr_d03ds.Draw("PL");
  l1->Draw();
  c1->SaveAs("ipEffic.png");  

}
