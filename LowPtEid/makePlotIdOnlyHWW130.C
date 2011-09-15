{

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  TFile *_file0 = TFile::Open("baby_JunkRun2011A.root");
  TFile *_file1 = TFile::Open("baby_GluGluToHToWWTo2L2Nu_M-130_7TeV.root");

  TTree* tree0 = (TTree*) _file0->Get("tree");
  TTree* tree1 = (TTree*) _file1->Get("tree");

  TCanvas *c1 = new TCanvas();
  c1->cd();

  TCut bar = "abs(el_etaSC)<1.479";
  TCut end = "abs(el_etaSC)>1.479";
  TCut cutNoPt = "met<20&&el_Mt<20";
  TCut base0 = "el_MZ<0&&met<20&&el_Mt<20&&(el8_v1||el8_v2||el8_v3)&&el_pt<200";
  TCut base1 = "el_pt<200";
  
  TCut v0id("el_VBTF80 && (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95&&abs(el_dPhiIn*el_q)<0.006))") ;
  TCut eta("abs(el_etaSC)<2.2");
  TCut mva("el_mva>0.4");
  TCut iso("20*el_relIso/el_pt<0.1");
  //iso = "(el_CICt&2)>0";
  TCut ip("abs(el_d0corr)<0.02");
  //TCut newconv("abs(el_newconv_dist)>0.05||abs(el_newconv_dcot)>0.02||abs(el_newconv_dmh)>1");
  //TCut conv("(abs(el_conv_dist)>0.02||abs(el_conv_dcot)>0.02) && el_innerlayer==0");
  TCut conv("el_mitconv==0 && el_innerlayer==0");
  TCut newconv = conv;
 
  TCut vbtf90_id  ("VBTF90"  ,"el_VBTF90");
  TCut vbtf85_id  ("VBTF85"  ,"el_VBTF90 && (abs(el_etaSC)<1.479&&abs(el_dPhiIn)<0.06&&abs(el_dEtaIn)<0.006&&el_hOverE<0.04 || abs(el_etaSC)>1.479&&abs(el_dPhiIn)<0.04&&abs(el_dEtaIn)<0.007&&el_hOverE<0.025)");
  TCut vbtf80_id  ("VBTF80"  ,"el_VBTF80");
  TCut vbtf70_id  ("VBTF70"  ,"el_VBTF70");
  TCut vbtf60_id  ("VBTF60"  ,"el_VBTF70 && (abs(el_etaSC)<1.479&&abs(el_dPhiIn)<0.025 || abs(el_etaSC)>1.479)");
  TCut vbtf70p_id ("VBTF70+" ,"el_VBTF70 && (el_pt>20 || (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95)))");
  TCut vbtf70pp_id("VBTF70++","el_VBTF70 && (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95))");
  TCut vbtf80p_id ("VBTF80+" ,"el_VBTF80 && (el_pt>20 || (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95)))");
  TCut vbtf80pp_id("VBTF80++","el_VBTF80 && (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95))");

  TCut vbtf80_nohoeend_id ("VBTF80-","el_VBTF80 & abs(el_etaSC)<1.479 || abs(el_etaSC)>1.479&&el_sigmaIEtaIEta<0.03&&abs(el_dPhiIn)<0.03&&abs(el_dEtaIn)<0.007");
  TCut vbtf80p_nohoeend_id ("VBTF80-+","(el_VBTF80 & abs(el_etaSC)<1.479 || abs(el_etaSC)>1.479&&el_sigmaIEtaIEta<0.03&&abs(el_dPhiIn)<0.03&&abs(el_dEtaIn)<0.007) && (el_pt>20 || (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95)))");
  TCut vbtf80pp_nohoeend_id ("VBTF80-++","(el_VBTF80 & abs(el_etaSC)<1.479 || abs(el_etaSC)>1.479&&el_sigmaIEtaIEta<0.03&&abs(el_dPhiIn)<0.03&&abs(el_dEtaIn)<0.007) && (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95))");
  TCut vbtf70_nohoeend_id ("VBTF70-","el_VBTF70 & abs(el_etaSC)<1.479 || abs(el_etaSC)>1.479&&el_sigmaIEtaIEta<0.03&&abs(el_dPhiIn)<0.02&&abs(el_dEtaIn)<0.005");
  TCut vbtf70p_nohoeend_id ("VBTF70-+","(el_VBTF70 & abs(el_etaSC)<1.479 || abs(el_etaSC)>1.479&&el_sigmaIEtaIEta<0.03&&abs(el_dPhiIn)<0.02&&abs(el_dEtaIn)<0.005) && (el_pt>20 || (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95)))");
  TCut vbtf70pp_nohoeend_id ("VBTF70-++","(el_VBTF70 & abs(el_etaSC)<1.479 || abs(el_etaSC)>1.479&&el_sigmaIEtaIEta<0.03&&abs(el_dPhiIn)<0.02&&abs(el_dEtaIn)<0.005) && (el_fbrem>0.15 || (abs(el_etaSC)<1.&&el_eOverPIn>0.95))");

  TCut cic_ht1_id  ("cic_ht1","(el_CICht1&1)>0");
  TCut cic_ht1_iso ("cic_ht1_iso","(el_CICht1&2)>0");
  TCut cic_ht1_conv("cic_ht1_conv","(el_CICht1&4)>0");
  TCut cic_ht1_ip  ("cic_ht1_ip","(el_CICht1&8)>0");
  TCut cic_ht2_id  ("cic_ht2","(el_CICht2&1)>0");
  TCut cic_ht2_iso ("cic_ht2_iso","(el_CICht2&2)>0");
  TCut cic_ht2_conv("cic_ht2_conv","(el_CICht2&4)>0");
  TCut cic_ht2_ip  ("cic_ht2_ip","(el_CICht2&8)>0");
  TCut cic_ht3_id  ("cic_ht3","(el_CICht3&1)>0");
  TCut cic_ht3_iso ("cic_ht3_iso","(el_CICht3&2)>0");
  TCut cic_ht3_conv("cic_ht3_conv","(el_CICht3&4)>0");
  TCut cic_ht3_ip  ("cic_ht3_ip","(el_CICht3&8)>0");
  TCut cic_ht4_id  ("cic_ht4","(el_CICht4&1)>0");
  TCut cic_ht4_iso ("cic_ht4_iso","(el_CICht4&2)>0");
  TCut cic_ht4_conv("cic_ht4_conv","(el_CICht4&4)>0");
  TCut cic_ht4_ip  ("cic_ht4_ip","(el_CICht4&8)>0");
  TCut cic_st_id   ("cic_st","(el_CICst&1)>0");
  TCut cic_st_iso  ("cic_st_iso","(el_CICst&2)>0");
  TCut cic_st_conv ("cic_st_conv","(el_CICst&4)>0");
  TCut cic_st_ip   ("cic_st_ip","(el_CICst&8)>0");
  TCut cic_t_id    ("cic_t","(el_CICt&1)>0");
  TCut cic_t_iso   ("cic_t_iso","(el_CICt&2)>0");
  TCut cic_t_conv  ("cic_t_conv","(el_CICt&4)>0");
  TCut cic_t_ip    ("cic_t_ip","(el_CICt&8)>0");

  TCut lik_t_id ("lht","abs(el_etaSC)<1.479&&el_fbrem<0.15&&el_LL>0.820 || abs(el_etaSC)<1.479&&el_fbrem>0.15&&el_LL>0.925 || abs(el_etaSC)>1.479&&el_fbrem<0.15&&el_LL>0.930 || abs(el_etaSC)>1.479&&el_fbrem>0.15&&el_LL>0.979");

  TCut lik95_id ("lh95","abs(el_etaSC)<1.479&&el_nbrem==0&&el_lh>-4.274 || abs(el_etaSC)<1.479&&el_nbrem>=1&&el_lh>-3.773 || abs(el_etaSC)>1.479&&el_nbrem==0&&el_lh>-5.092 || abs(el_etaSC)>1.479&&el_nbrem>=1&&el_lh>-2.796");
  TCut lik90_id ("lh90","abs(el_etaSC)<1.479&&el_nbrem==0&&el_lh>-1.497 || abs(el_etaSC)<1.479&&el_nbrem>=1&&el_lh>-1.521 || abs(el_etaSC)>1.479&&el_nbrem==0&&el_lh>-2.571 || abs(el_etaSC)>1.479&&el_nbrem>=1&&el_lh>-0.657");
  TCut lik85_id ("lh85","abs(el_etaSC)<1.479&&el_nbrem==0&&el_lh>+0.163 || abs(el_etaSC)<1.479&&el_nbrem>=1&&el_lh>+0.065 || abs(el_etaSC)>1.479&&el_nbrem==0&&el_lh>-0.683 || abs(el_etaSC)>1.479&&el_nbrem>=1&&el_lh>+1.564");
  TCut lik80_id ("lh80","abs(el_etaSC)<1.479&&el_nbrem==0&&el_lh>+1.193 || abs(el_etaSC)<1.479&&el_nbrem>=1&&el_lh>+1.345 || abs(el_etaSC)>1.479&&el_nbrem==0&&el_lh>+0.810 || abs(el_etaSC)>1.479&&el_nbrem>=1&&el_lh>+3.021");
  TCut lik70_id ("lh70","abs(el_etaSC)<1.479&&el_nbrem==0&&el_lh>+1.781 || abs(el_etaSC)<1.479&&el_nbrem>=1&&el_lh>+2.397 || abs(el_etaSC)>1.479&&el_nbrem==0&&el_lh>+2.361 || abs(el_etaSC)>1.479&&el_nbrem>=1&&el_lh>+4.052");


  TCut pfmva3_id ("mva3","el_mva>0.3");
  TCut pfmva5_id ("mva5","el_mva>0.5");
  TCut pfmva6_id ("mva6","el_mva>0.6");
  TCut pfmva7_id ("mva7","el_mva>0.7");
  TCut pfmva8_id ("mva8","el_mva>0.8");

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


  base0=base0+ip+conv+iso;
  base1=base1+ip+conv+iso;

  float ptbins[] = {20.,999.};
  float etabins[] = {0,2.5};
  TH2F* pt_vs_eta = new TH2F("pt_vs_eta","pt_vs_eta",1,etabins,1,ptbins);

  TGraph gr(2);
  gr.SetTitle("");
  gr.SetPoint(0, 0, 0);
  gr.SetPoint(1, 1.3, 2.0);
  gr.SetMarkerStyle(1);

  TCut cuts[] = {vbtf90_id,vbtf85_id,vbtf80_id,vbtf70_id,vbtf80p_id,vbtf70p_id, 
		 vbtf80_nohoeend_id,vbtf70_nohoeend_id,vbtf80p_nohoeend_id,vbtf70p_nohoeend_id,
		 cic_t_id,cic_st_id,cic_ht1_id,cic_ht2_id,cic_ht3_id,cic_ht4_id,
		 lik95_id,lik90_id,lik85_id,lik80_id,lik70_id,
		 pfmva3_id,pfmva5_id,pfmva6_id,pfmva7_id,pfmva8_id};
  int ncuts = sizeof(cuts)/sizeof(TCut);

  tree0->Draw("el_pt:abs(el_etaSC)>>pt_vs_eta",base0+vbtf80_id,"goff");
  float all0 = pt_vs_eta->GetBinContent(1,1);
  tree1->Draw("el_pt:abs(el_etaSC)>>pt_vs_eta",base1+vbtf80_id,"goff");
  float all1 = pt_vs_eta->GetBinContent(1,1);

  TLegend *l1 = new TLegend(0.2, 0.7, 0.6, 0.9);
  l1->SetNColumns(2);
  l1->SetFillColor(kWhite);

  gr.Draw("AP");
  gr.GetXaxis()->SetRangeUser(0.4,1.3);
  gr.GetYaxis()->SetRangeUser(0.,2.);
  gr.GetXaxis()->SetTitle("Electron ID/VBTF80 eff (HWW130 MC)");
  gr.GetYaxis()->SetTitle("Electron ID/VBTF80 eff (2011 data)");

  TGraph gr_vbtf(6);
  gr_vbtf.SetMarkerStyle(20); 
  gr_vbtf.SetMarkerColor(kBlue);
  gr_vbtf.SetLineColor(kBlue);

  TGraph gr_vbtf_nohoeend(4);
  gr_vbtf_nohoeend.SetMarkerStyle(20); 
  gr_vbtf_nohoeend.SetMarkerColor(kMagenta);
  gr_vbtf_nohoeend.SetLineColor(kMagenta);

  TGraph gr_cic(6);
  gr_cic.SetMarkerStyle(21); 
  gr_cic.SetMarkerColor(kRed);
  gr_cic.SetLineColor(kRed);

  TGraph gr_mva(5);
  gr_mva.SetMarkerStyle(22); 
  gr_mva.SetMarkerColor(kGreen);
  gr_mva.SetLineColor(kGreen);

  TGraph gr_lik(3);
  gr_lik.SetMarkerStyle(23); 
  gr_lik.SetMarkerColor(kOrange);
  gr_lik.SetLineColor(kOrange);

  int n_vbtf=0,n_vbtf_nohoeend=0,n_cic=0,n_mva=0,n_lik=0;
  for (unsigned int ic=0;ic<ncuts;++ic){
    tree0->Draw("el_pt:abs(el_etaSC)>>pt_vs_eta",base0+cuts[ic],"goff");
    float cut0 = pt_vs_eta->GetBinContent(1,1);
    tree1->Draw("el_pt:abs(el_etaSC)>>pt_vs_eta",base1+cuts[ic],"goff");
    float cut1 = pt_vs_eta->GetBinContent(1,1);
    if (TString(cuts[ic].GetName()).Contains("VBTF")&&!TString(cuts[ic].GetName()).Contains("-")) gr_vbtf.SetPoint(n_vbtf++, cut1/all1, cut0/all0);
    if (TString(cuts[ic].GetName()).Contains("VBTF")&&TString(cuts[ic].GetName()).Contains("-")) gr_vbtf_nohoeend.SetPoint(n_vbtf_nohoeend++, cut1/all1, cut0/all0);
    if (TString(cuts[ic].GetName()).Contains("cic")) gr_cic.SetPoint(n_cic++, cut1/all1, cut0/all0);
    if (TString(cuts[ic].GetName()).Contains("mva")) gr_mva.SetPoint(n_mva++, cut1/all1, cut0/all0);
    if (TString(cuts[ic].GetName()).Contains("lh")) gr_lik.SetPoint(n_lik++, cut1/all1, cut0/all0);
    cout << cuts[ic].GetName() << " hww: " << Form("%4.2f",cut1/all1) << " - eg2010a: " << Form("%4.2f",cut0/all0) << endl;

  }

  l1->AddEntry(&gr_vbtf, "VBTF", "p");
  gr_vbtf.Draw("PL");
  l1->AddEntry(&gr_vbtf_nohoeend, "VBTF-", "p");
  gr_vbtf_nohoeend.Draw("PL");
  l1->AddEntry(&gr_cic, "CIC", "p");
  gr_cic.Draw("PL");
  l1->AddEntry(&gr_mva, "MVA", "p");
  gr_mva.Draw("PL");
  l1->AddEntry(&gr_lik, "LH", "p");
  gr_lik.Draw("PL");
  l1->Draw();
  c1->SaveAs("idEffic_20ptup.png");  

}
