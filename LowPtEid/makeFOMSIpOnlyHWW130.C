{

  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  TFile *_fileB = TFile::Open("testIP_baby_EGMon2010B.root");
  TFile *_fileS = TFile::Open("testIP_baby_GluGluToHToWWTo2L2Nu_M-130_7TeV_LOPU.root");

  TCanvas c1;

  TTree* treeB = (TTree*) _fileB->Get("tree");
  TTree* treeS = (TTree*) _fileS->Get("tree");

  TCut bar = "abs(el_etaSC)<1.479";
  TCut end = "abs(el_etaSC)>1.479";
  TCut ptcut = "";
  TCut baseB = "el_MZ<0&&met<20&&el_Mt<20&&el10_sw"+ptcut;
  TCut baseS = ""+ptcut;

  TCut ip("abs(el_d0corr)<0.02");
  TCut conv("(abs(el_conv_dist)>0.02||abs(el_conv_dcot)>0.02) && el_innerlayer39X==0");

  TCut iso("20*el_relIso/el_pt<0.1");
  
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

  //iso = cic_ht2_iso;

  TCut lik_t_id ("lht","abs(el_etaSC)<1.479&&el_fbrem<0.15&&el_LL>0.820 || abs(el_etaSC)<1.479&&el_fbrem>0.15&&el_LL>0.925 || abs(el_etaSC)>1.479&&el_fbrem<0.15&&el_LL>0.930 || abs(el_etaSC)>1.479&&el_fbrem>0.15&&el_LL>0.979");
  TCut pfmva3_id ("mva3","el_mva>0.3");
  TCut pfmva5_id ("mva5","el_mva>0.5");
  TCut pfmva6_id ("mva6","el_mva>0.6");
  TCut pfmva7_id ("mva7","el_mva>0.7");
  TCut pfmva8_id ("mva8","el_mva>0.8");

  TCut ip_d02d_0008("d02d_0008","abs(el_d0pv2d)<0.008");
  TCut ip_d02d_001("d02d_001","abs(el_d0pv2d)<0.01");
  TCut ip_d02d_0015("d02d_0015","abs(el_d0pv2d)<0.015");
  TCut ip_d02d_002("d02d_002","abs(el_d0pv2d)<0.02");
  TCut ip_d02d_005("d02d_005","abs(el_d0pv2d)<0.05");

  TCut ip_d03d_002("d03d_002","abs(el_d0pv3d)<0.02");
  TCut ip_d03d_003("d03d_003","abs(el_d0pv3d)<0.03");
  TCut ip_d03d_005("d03d_005","abs(el_d0pv3d)<0.05");
  TCut ip_d03d_010("d03d_010","abs(el_d0pv3d)<0.10");
  TCut ip_d03d_020("d03d_020","abs(el_d0pv3d)<0.20");

  TCut ip_d02ds_12("d02ds_12","abs(el_d0pv2d/el_d0pv2dErr)<1.2");
  TCut ip_d02ds_14("d02ds_14","abs(el_d0pv2d/el_d0pv2dErr)<1.4");
  TCut ip_d02ds_17("d02ds_17","abs(el_d0pv2d/el_d0pv2dErr)<1.7");
  TCut ip_d02ds_2("d02ds_2","abs(el_d0pv2d/el_d0pv2dErr)<2");
  TCut ip_d02ds_5("d02ds_5","abs(el_d0pv2d/el_d0pv2dErr)<5");

  TCut ip_d03ds_16("d03ds_16","abs(el_d0pv3d/el_d0pv3dErr)<1.6");
  TCut ip_d03ds_18("d03ds_18","abs(el_d0pv3d/el_d0pv3dErr)<1.8");
  TCut ip_d03ds_2("d03ds_2","abs(el_d0pv3d/el_d0pv3dErr)<2");
  TCut ip_d03ds_3("d03ds_3","abs(el_d0pv3d/el_d0pv3dErr)<3");
  TCut ip_d03ds_9("d03ds_9","abs(el_d0pv3d/el_d0pv3dErr)<9");

  TCut ip_d02d_sli("d02d_sli","el_pt<20&&abs(el_d0pv2d)<0.010 || el_pt>20&&abs(el_d0pv2d)<0.020");


  baseB=baseB+conv+iso+vbtf80_id;
  baseS=baseS+conv+iso+vbtf80_id;

  float do1015 = 0.;
  float do1520 = 1.;
  float do20up = 1.;

  TString ptrange = "";
  if (do1015&&do1520&&do20up) ptrange = "pt10";
  else if (do1520&&do20up) ptrange = "pt15";
  else if (do20up) ptrange = "pt20";

//   TCut cuts[] = {ip_d02d_sli,ip_d03d_010,ip_d02d_002}

  TCut cuts[] = {ip_d02d_0008,ip_d02d_001,ip_d02d_0015,ip_d02d_002,ip_d02d_005,
		 ip_d03d_002,ip_d03d_003,ip_d03d_005,ip_d03d_010,ip_d03d_020,
		 ip_d02ds_12,ip_d02ds_14,ip_d02ds_17,ip_d02ds_2,ip_d02ds_5,
  		 ip_d03ds_16,ip_d03ds_18,ip_d03ds_2,ip_d03ds_3,ip_d03ds_9};

  float ptbins[] =  {10., 15., 20., 9999.};
  float etabins[] = {0.,  1.,  1.5, 2.5};
  TH2F* pt_vs_eta = new TH2F("pt_vs_eta","pt_vs_eta",3,etabins,3,ptbins);

  treeB->Draw("el_pt:abs(el_etaSC)>>pt_vs_eta",baseB+ip_d02d_002,"goff");
  float allB_11 = pt_vs_eta->GetBinContent(1,1);
  float allB_21 = pt_vs_eta->GetBinContent(2,1);
  float allB_31 = pt_vs_eta->GetBinContent(3,1);
  float allB_12 = pt_vs_eta->GetBinContent(1,2);
  float allB_22 = pt_vs_eta->GetBinContent(2,2);
  float allB_32 = pt_vs_eta->GetBinContent(3,2);
  float allB_13 = pt_vs_eta->GetBinContent(1,3);
  float allB_23 = pt_vs_eta->GetBinContent(2,3);
  float allB_33 = pt_vs_eta->GetBinContent(3,3);
  treeS->Draw("el_pt:abs(el_etaSC)>>pt_vs_eta",baseS+ip_d02d_002,"goff");
  float allS_11 = pt_vs_eta->GetBinContent(1,1);
  float allS_21 = pt_vs_eta->GetBinContent(2,1);
  float allS_31 = pt_vs_eta->GetBinContent(3,1);
  float allS_12 = pt_vs_eta->GetBinContent(1,2);
  float allS_22 = pt_vs_eta->GetBinContent(2,2);
  float allS_32 = pt_vs_eta->GetBinContent(3,2);
  float allS_13 = pt_vs_eta->GetBinContent(1,3);
  float allS_23 = pt_vs_eta->GetBinContent(2,3);
  float allS_33 = pt_vs_eta->GetBinContent(3,3);

  //FAKES
  // pT 10 - 15 : 13.70: 9.11, 3.42, 1.18
  // pT 15-20   :  4.15: 3.22, 0.0 , 0.93
  // pT 20-Inf  :  5.14: 1.24, 2.85, 1.05
  float totB_11 = 9.11;
  float totB_21 = 3.42;
  float totB_31 = 1.18;
  float totB_12 = 3.22;
  float totB_22 = 0.00;
  float totB_32 = 0.93;
  float totB_13 = 1.24;
  float totB_23 = 2.85;
  float totB_33 = 1.05;
  float totB = totB_11+totB_21+totB_31+totB_12+totB_22+totB_32+totB_13+totB_23+totB_33;

  //ALL BACKGROUNDS EXCEPT WJETS
  // pT 10 - 15 : 4.00 : 2.15,  0.81, 1.03
  // pT 15-20   : 6.87 : 4.18,  1.42, 1.26
  // pT 20-Inf  : 18.71: 10.26, 4.89, 3.57
  float totOB_11 = 2.15;
  float totOB_21 = 0.81;
  float totOB_31 = 1.03;
  float totOB_12 = 4.18;
  float totOB_22 = 1.42;
  float totOB_32 = 1.26;
  float totOB_13 = 10.26;
  float totOB_23 = 4.89;
  float totOB_33 = 3.57;
  float totOB = totOB_11+totOB_21+totOB_31+totOB_12+totOB_22+totOB_32+totOB_13+totOB_23+totOB_33;

  //HWW130
  // pT 10 - 15 : 1.19: 0.79, 0.22, 0.18
  // pT 15-20   : 1.77: 1.09, 0.37, 0.32
  // pT 20-Inf  : 3.41: 2.18, 0.69, 0.54
  float totS_11 = 0.79;
  float totS_21 = 0.22;
  float totS_31 = 0.18;
  float totS_12 = 1.09;
  float totS_22 = 0.37;
  float totS_32 = 0.32;
  float totS_13 = 2.18;
  float totS_23 = 0.69;
  float totS_33 = 0.54;
  float totS = totS_11+totS_21+totS_31+totS_12+totS_22+totS_32+totS_13+totS_23+totS_33;

  vector<TString> names;
  vector<float> fom1s;
  vector<float> fom2s;
  vector<float> fom3s;

  int ncuts = sizeof(cuts)/sizeof(TCut);
  for (int id_it=0;id_it<ncuts;++id_it) {

    TCut idTest = cuts[id_it];
	      
    treeB->Draw("el_pt:abs(el_etaSC)>>pt_vs_eta",baseB+idTest,"goff");
    float id_cutB_11 = pt_vs_eta->GetBinContent(1,1);
    float id_cutB_21 = pt_vs_eta->GetBinContent(2,1);
    float id_cutB_31 = pt_vs_eta->GetBinContent(3,1);
    float id_cutB_12 = pt_vs_eta->GetBinContent(1,2);
    float id_cutB_22 = pt_vs_eta->GetBinContent(2,2);
    float id_cutB_32 = pt_vs_eta->GetBinContent(3,2);
    float id_cutB_13 = pt_vs_eta->GetBinContent(1,3);
    float id_cutB_23 = pt_vs_eta->GetBinContent(2,3);
    float id_cutB_33 = pt_vs_eta->GetBinContent(3,3);
    treeS->Draw("el_pt:abs(el_etaSC)>>pt_vs_eta",baseS+idTest,"goff");
    float id_cutS_11 = pt_vs_eta->GetBinContent(1,1);
    float id_cutS_21 = pt_vs_eta->GetBinContent(2,1);
    float id_cutS_31 = pt_vs_eta->GetBinContent(3,1);
    float id_cutS_12 = pt_vs_eta->GetBinContent(1,2);
    float id_cutS_22 = pt_vs_eta->GetBinContent(2,2);
    float id_cutS_32 = pt_vs_eta->GetBinContent(3,2);
    float id_cutS_13 = pt_vs_eta->GetBinContent(1,3);
    float id_cutS_23 = pt_vs_eta->GetBinContent(2,3);
    float id_cutS_33 = pt_vs_eta->GetBinContent(3,3);
    
    //compute projected yields for fake bkg
    float B_11 = do1015*totB_11*id_cutB_11/allB_11;
    float B_21 = do1015*totB_21*id_cutB_21/allB_21;
    float B_31 = do1015*totB_31*id_cutB_31/allB_31;
    float B_12 = do1520*totB_12*id_cutB_12/allB_12;
    float B_22 = do1520*totB_22*id_cutB_22/allB_22;
    float B_32 = do1520*totB_32*id_cutB_32/allB_32;
    float B_13 = do20up*totB_13*id_cutB_13/allB_13;
    float B_23 = do20up*totB_23*id_cutB_23/allB_23;
    float B_33 = do20up*totB_33*id_cutB_33/allB_33;
    
    //compute projected yields for other bkg using signal eff
    float OB_11 = do1015*totOB_11*id_cutS_11/allS_11;
    float OB_21 = do1015*totOB_21*id_cutS_21/allS_21;
    float OB_31 = do1015*totOB_31*id_cutS_31/allS_31;
    float OB_12 = do1520*totOB_12*id_cutS_12/allS_12;
    float OB_22 = do1520*totOB_22*id_cutS_22/allS_22;
    float OB_32 = do1520*totOB_32*id_cutS_32/allS_32;
    float OB_13 = do20up*totOB_13*id_cutS_13/allS_13;
    float OB_23 = do20up*totOB_23*id_cutS_23/allS_23;
    float OB_33 = do20up*totOB_33*id_cutS_33/allS_33;
    
    //compute projected yields for signal
    float S_11 = do1015*totS_11*id_cutS_11/allS_11;
    float S_21 = do1015*totS_21*id_cutS_21/allS_21;
    float S_31 = do1015*totS_31*id_cutS_31/allS_31;
    float S_12 = do1520*totS_12*id_cutS_12/allS_12;
    float S_22 = do1520*totS_22*id_cutS_22/allS_22;
    float S_32 = do1520*totS_32*id_cutS_32/allS_32;
    float S_13 = do20up*totS_13*id_cutS_13/allS_13;
    float S_23 = do20up*totS_23*id_cutS_23/allS_23;
    float S_33 = do20up*totS_33*id_cutS_33/allS_33;
    
    //cout << "id hww: " << Form("%4.2f",id_cut1/all1) << " - eg2010a: " << Form("%4.2f",id_cut0/all0) 
    //<< " - fom: " << (id_cut1/all1)*(id_cut1/all1)/(id_cut0/all0)  << endl;
    
    bool printIt = false;
    if (printIt) {
      cout << "EG2010 efficiencies" << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",id_cutB_11/allB_11,id_cutB_21/allB_21,id_cutB_31/allB_31) << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",id_cutB_12/allB_12,id_cutB_22/allB_22,id_cutB_32/allB_32) << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",id_cutB_13/allB_13,id_cutB_23/allB_23,id_cutB_33/allB_33) << endl;
      cout << "HWW130 efficiencies" << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",id_cutS_11/allS_11,id_cutS_21/allS_21,id_cutS_31/allS_31) << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",id_cutS_12/allS_12,id_cutS_22/allS_22,id_cutS_32/allS_32) << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",id_cutS_13/allS_13,id_cutS_23/allS_23,id_cutS_33/allS_33) << endl;
      
      cout << "Bkg from Fake Yields" << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",B_11,B_21,B_31) << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",B_12,B_22,B_32) << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",B_13,B_23,B_33) << endl;
      cout << "Non Fake Bkg Yields" << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",OB_11,OB_21,OB_31) << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",OB_12,OB_22,OB_32) << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",OB_13,OB_23,OB_33) << endl;
      cout << "HWW130 Yields" << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",S_11,S_21,S_31) << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",S_12,S_22,S_32) << endl;
      cout << Form("%3.2f - %3.2f - %3.2f",S_13,S_23,S_33) << endl;
    }
    
    float OB=OB_11+OB_21+OB_31+OB_12+OB_22+OB_32+OB_13+OB_23+OB_33;
    float FB=B_11+B_21+B_31+B_12+B_22+B_32+B_13+B_23+B_33;
    float B=OB+FB;
    float S=S_11+S_21+S_31+S_12+S_22+S_32+S_13+S_23+S_33;
    
    float fom1 = S/B;
    float fom2 = S/sqrt(S+B);
    float fom3 = S/sqrt(S+B+(0.35*B)*(0.35*B));

    names.push_back(cuts[id_it].GetName());
    fom1s.push_back(fom1);
    fom2s.push_back(fom2);
    fom3s.push_back(fom3);


    cout << Form("%2i) S: %5.2f FB: %5.2f OB: %5.2f fom1: %5.3f fom2: %5.3f fom3: %5.3f name: %10s",S,FB,OB,id_it,fom1,fom2,fom3,cuts[id_it].GetName()) << endl;
  }

  for (int i=0;i<fom1s.size();++i) cout << fom1s[i] << ",";
  cout << endl;
  for (int i=0;i<fom2s.size();++i) cout << fom2s[i] << ",";
  cout << endl;
  for (int i=0;i<fom3s.size();++i) cout << fom3s[i] << ",";
  cout << endl;
  for (int i=0;i<names.size();++i) cout << "\"" << names[i] << "\",";
  cout << endl;

  c1->SetBottomMargin(0.20);
  c1->SetGridy();

  TH1F* hfom1_d02d = new TH1F("fom1_d02d","FOM1: S/B",fom1s.size(),0,fom1s.size());
  hfom1_d02d->SetMarkerStyle(20);
  hfom1_d02d->SetMarkerColor(kBlue);
  hfom1_d02d->GetXaxis()->SetLabelColor(kBlue);
  TH1F* hfom1_d03d = new TH1F("fom1_d03d","FOM1: S/B",fom1s.size(),0,fom1s.size());
  hfom1_d03d->SetMarkerStyle(20);
  hfom1_d03d->SetMarkerColor(kRed);
  hfom1_d03d->GetXaxis()->SetLabelColor(kRed);
  TH1F* hfom1_d02ds = new TH1F("fom1_d02ds","FOM1: S/B",fom1s.size(),0,fom1s.size());
  hfom1_d02ds->SetMarkerStyle(21);
  hfom1_d02ds->SetMarkerColor(kBlue);
  hfom1_d02ds->GetXaxis()->SetLabelColor(kBlue);
  TH1F* hfom1_d03ds = new TH1F("fom1_d03ds","FOM1: S/B",fom1s.size(),0,fom1s.size());
  hfom1_d03ds->SetMarkerStyle(21);
  hfom1_d03ds->SetMarkerColor(kRed);
  hfom1_d03ds->GetXaxis()->SetLabelColor(kRed);
  for (int j=0;j<fom1s.size();++j) {
    if (names[j].Contains("d02d_")) {
      hfom1_d02d->SetBinContent(j+1,fom1s[j]);
      hfom1_d02d->GetXaxis()->SetBinLabel(j+1,"");
    }
    else if (names[j].Contains("d03d_")) {
      hfom1_d03d->SetBinContent(j+1,fom1s[j]);
      hfom1_d03d->GetXaxis()->SetBinLabel(j+1,"");
    }
    else if (names[j].Contains("d02ds_")) {
      hfom1_d02ds->SetBinContent(j+1,fom1s[j]);
      hfom1_d02ds->GetXaxis()->SetBinLabel(j+1,"");
    }
    else if (names[j].Contains("d03ds_")) {
      hfom1_d03ds->SetBinContent(j+1,fom1s[j]);
      hfom1_d03ds->GetXaxis()->SetBinLabel(j+1,"");
    }
  }
  float ymax = max(hfom1_d02d->GetMaximum(),max(hfom1_d03d->GetMaximum(),max(hfom1_d02ds->GetMaximum(),hfom1_d03ds->GetMaximum())))*1.2;
  hfom1_d02d->GetYaxis()->SetRangeUser(0,0.2);
  hfom1_d02d->Draw("P");
  Float_t x, y;
  y = -ymax/90.;
  TText t;
  t.SetTextAngle(90);
  t.SetTextSize(0.03);
  t.SetTextAlign(32);
  for (i=0;i<fom1s.size();i++) {
    x = hfom1_d02d->GetXaxis()->GetBinCenter(i+1);
    if (names[i].Contains("d02d")) t.SetTextColor(kBlue);
    else if (names[i].Contains("d03d")) t.SetTextColor(kRed);
    t.DrawText(x,y,names[i]);
  }
  hfom1_d03d->Draw("P,SAME");
  hfom1_d02ds->Draw("P,SAME");
  hfom1_d03ds->Draw("P,SAME");
  c1.SaveAs("ip_hww130_fom1_"+ptrange+".png");

  TH1F* hfom2_d02d = new TH1F("fom2_d02d","FOM2: S/sqrt(S+B)",fom2s.size(),0,fom2s.size());
  hfom2_d02d->SetMarkerStyle(20);
  hfom2_d02d->SetMarkerColor(kBlue);
  hfom2_d02d->GetXaxis()->SetLabelColor(kBlue);
  TH1F* hfom2_d03d = new TH1F("fom2_d03d","FOM2: S/sqrt(S+B)",fom2s.size(),0,fom2s.size());
  hfom2_d03d->SetMarkerStyle(20);
  hfom2_d03d->SetMarkerColor(kRed);
  hfom2_d03d->GetXaxis()->SetLabelColor(kRed);
  TH1F* hfom2_d02ds = new TH1F("fom2_d02ds","FOM2: S/sqrt(S+B)",fom2s.size(),0,fom2s.size());
  hfom2_d02ds->SetMarkerStyle(21);
  hfom2_d02ds->SetMarkerColor(kBlue);
  hfom2_d02ds->GetXaxis()->SetLabelColor(kBlue);
  TH1F* hfom2_d03ds = new TH1F("fom2_d03ds","FOM2: S/sqrt(S+B)",fom2s.size(),0,fom2s.size());
  hfom2_d03ds->SetMarkerStyle(21);
  hfom2_d03ds->SetMarkerColor(kRed);
  hfom2_d03ds->GetXaxis()->SetLabelColor(kRed);
  for (int j=0;j<fom2s.size();++j) {
    if (names[j].Contains("d02d_")) {
      hfom2_d02d->SetBinContent(j+1,fom2s[j]);
      hfom2_d02d->GetXaxis()->SetBinLabel(j+1,"");
    }
    else if (names[j].Contains("d03d_")) {
      hfom2_d03d->SetBinContent(j+1,fom2s[j]);
      hfom2_d03d->GetXaxis()->SetBinLabel(j+1,"");
    }
    else if (names[j].Contains("d02ds_")) {
      hfom2_d02ds->SetBinContent(j+1,fom2s[j]);
      hfom2_d02ds->GetXaxis()->SetBinLabel(j+1,"");
    }
    else if (names[j].Contains("d03ds_")) {
      hfom2_d03ds->SetBinContent(j+1,fom2s[j]);
      hfom2_d03ds->GetXaxis()->SetBinLabel(j+1,"");
    }
  }
  float ymax = max(hfom2_d02d->GetMaximum(),max(hfom2_d03d->GetMaximum(),max(hfom2_d02ds->GetMaximum(),hfom2_d03ds->GetMaximum())))*1.2;
  hfom2_d02d->GetYaxis()->SetRangeUser(0,1.0);
  hfom2_d02d->Draw("P");
  Float_t x, y;
  y = -ymax/90.;
  TText t;
  t.SetTextAngle(90);
  t.SetTextSize(0.03);
  t.SetTextAlign(32);
  for (i=0;i<fom2s.size();i++) {
    x = hfom2_d02d->GetXaxis()->GetBinCenter(i+1);
    if (names[i].Contains("d02d")) t.SetTextColor(kBlue);
    else if (names[i].Contains("d03d")) t.SetTextColor(kRed);
    t.DrawText(x,y,names[i]);
  }
  hfom2_d03d->Draw("P,SAME");
  hfom2_d02ds->Draw("P,SAME");
  hfom2_d03ds->Draw("P,SAME");
  c1.SaveAs("ip_hww130_fom2_"+ptrange+".png");

  TH1F* hfom3_d02d = new TH1F("fom3_d02d","FOM3: S/sqrt(S+B+(0.35*B)^2)",fom3s.size(),0,fom3s.size());
  hfom3_d02d->SetMarkerStyle(20);
  hfom3_d02d->SetMarkerColor(kBlue);
  hfom3_d02d->GetXaxis()->SetLabelColor(kBlue);
  TH1F* hfom3_d03d = new TH1F("fom3_d03d","FOM3: S/sqrt(S+B+(0.35*B)^2)",fom3s.size(),0,fom3s.size());
  hfom3_d03d->SetMarkerStyle(20);
  hfom3_d03d->SetMarkerColor(kRed);
  hfom3_d03d->GetXaxis()->SetLabelColor(kRed);
  TH1F* hfom3_d02ds = new TH1F("fom3_d02ds","FOM3: S/sqrt(S+B+(0.35*B)^2)",fom3s.size(),0,fom3s.size());
  hfom3_d02ds->SetMarkerStyle(21);
  hfom3_d02ds->SetMarkerColor(kBlue);
  hfom3_d02ds->GetXaxis()->SetLabelColor(kBlue);
  TH1F* hfom3_d03ds = new TH1F("fom3_d03ds","FOM3: S/sqrt(S+B+(0.35*B)^2)",fom3s.size(),0,fom3s.size());
  hfom3_d03ds->SetMarkerStyle(21);
  hfom3_d03ds->SetMarkerColor(kRed);
  hfom3_d03ds->GetXaxis()->SetLabelColor(kRed);
  for (int j=0;j<fom3s.size();++j) {
    if (names[j].Contains("d02d_")) {
      hfom3_d02d->SetBinContent(j+1,fom3s[j]);
      hfom3_d02d->GetXaxis()->SetBinLabel(j+1,"");
    }
    else if (names[j].Contains("d03d_")) {
      hfom3_d03d->SetBinContent(j+1,fom3s[j]);
      hfom3_d03d->GetXaxis()->SetBinLabel(j+1,"");
    }
    else if (names[j].Contains("d02ds_")) {
      hfom3_d02ds->SetBinContent(j+1,fom3s[j]);
      hfom3_d02ds->GetXaxis()->SetBinLabel(j+1,"");
    }
    else if (names[j].Contains("d03ds_")) {
      hfom3_d03ds->SetBinContent(j+1,fom3s[j]);
      hfom3_d03ds->GetXaxis()->SetBinLabel(j+1,"");
    }
  }
  float ymax = max(hfom3_d02d->GetMaximum(),max(hfom3_d03d->GetMaximum(),max(hfom3_d02ds->GetMaximum(),hfom3_d03ds->GetMaximum())))*1.2;
  hfom3_d02d->GetYaxis()->SetRangeUser(0,0.5);
  hfom3_d02d->Draw("P");
  Float_t x, y;
  y = -ymax/90.;
  TText t;
  t.SetTextAngle(90);
  t.SetTextSize(0.03);
  t.SetTextAlign(32);
  for (i=0;i<fom3s.size();i++) {
    x = hfom3_d02d->GetXaxis()->GetBinCenter(i+1);
    if (names[i].Contains("d02d")) t.SetTextColor(kBlue);
    else if (names[i].Contains("d03d")) t.SetTextColor(kRed);
    t.DrawText(x,y,names[i]);
  }
  hfom3_d03d->Draw("P,SAME");
  hfom3_d02ds->Draw("P,SAME");
  hfom3_d03ds->Draw("P,SAME");
  c1.SaveAs("ip_hww130_fom3_"+ptrange+".png");

  c1.Update();
}
