{

  gStyle->SetOptStat(0);

  TCanvas c1;

  TFile *_file_ref = TFile::Open("cards_def_oldEffSyst/125/hwwof_0j.input_8TeV.root");
  TFile *_file_new = TFile::Open("cards_def_newEffSyst/125/hwwof_0j.input_8TeV.root");

  TIter next(_file_ref->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1F *h_ref = (TH1F*)key->ReadObj();

    if (TString(h_ref->GetName()).Contains("ttH")) continue;

    TH1F* h_new2 = 0;
    h_new2 = (TH1F*) _file_new->Get(h_ref->GetName());
    if (h_new2 == 0) {
      //h_new = (TH1F*) _file_new->Get(TString(h_ref->GetName()).ReplaceAll("Bounding","Bounding_8TeV"));
      //if (h_new == 0) {
	cout << "Histo not present in test file: " << h_ref->GetName() << endl;
	continue;
	//}
    }

    //make sure binning is the same
    TH1F* h_new = (TH1F*) h_ref->Clone("h_new");
    h_new->Clear();
    for (int bin=1;bin<=h_new->GetNbinsX();++bin) {
      h_new->SetBinContent(bin,h_new2->GetBinContent(bin));
    }

    //check diff
    TH1F* diff = (TH1F*) h_ref->Clone("diff");
    diff->Add(h_new,-1.);
    for (int bin=1;bin<=diff->GetNbinsX();++bin) {
      diff->SetBinContent(bin,fabs(diff->GetBinContent(bin)));
    }
    float maxdiff = TMath::Max(diff->GetBinContent(diff->GetMaximumBin()),diff->GetBinContent(diff->GetMinimumBin()));    
    //check ratio
    TH1F* ratio = (TH1F*) diff->Clone("ratio");
    for (int bin=1;bin<=ratio->GetNbinsX();++bin) {
      if (fabs(h_ref->GetBinContent(bin)>0) ) {
	ratio->SetBinContent(bin,ratio->GetBinContent(bin)/h_ref->GetBinContent(bin));
      } else ratio->SetBinContent(bin,0);
    }
    float maxratio = TMath::Max(ratio->GetBinContent(ratio->GetMaximumBin()),ratio->GetBinContent(ratio->GetMinimumBin()));

    if (1/*maxdiff>1 && maxratio>0.1*/) {
      cout << "Maximum difference for " << h_ref->GetName() << " : " 
	   << " diff=" << maxdiff << " bin=" << diff->GetMaximumBin()  << " content=" << h_ref->GetBinContent(diff->GetMaximumBin())
	   << " ratio=" << maxratio << " bin=" << ratio->GetMaximumBin()  << " content=" << h_ref->GetBinContent(ratio->GetMaximumBin())
	   << endl;
    } else {
      continue;
    }

    c1.cd();

    h_ref->SetLineWidth(2);
    h_ref->SetTitle("");
    h_ref->GetXaxis()->SetTitle("unrolled bins");
    h_ref->Draw("HIST");

    h_new->SetLineColor(kRed);
    h_new->SetLineWidth(2);
    h_new->Draw("SAME");
    
    TLegend* leg = new TLegend(0.15,0.91,0.85,0.99);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0);
    leg->SetNColumns(2);
    leg->AddEntry(h_ref,"ref histo","l");
    leg->AddEntry(h_new,"new histo","l");
    leg->Draw();

    c1.SetLogy();
    
    //   TPaveText * labelcms  = new TPaveText(0.11,0.75,0.38,0.88,"NDCBR");
    //   labelcms->SetTextAlign(12);
    //   labelcms->SetTextSize(0.03);
    //   labelcms->SetFillColor(kWhite);
    //   labelcms->AddText("CMS, #sqrt{s} = 8 TeV, L_{int} = 5.1 fb^{-1}");
    //   labelcms->AddText("H#rightarrowWW#rightarrow2l2#nu mH=160 GeV");
    //   labelcms->AddText("0-jet,OF shape based analyis");
    //   labelcms->SetBorderSize(0);
    //   labelcms->SetTextFont(62);
    //   labelcms->SetLineWidth(2);
    //   labelcms->Draw();
    
    gSystem->Exec("mkdir -p plots2D");
    c1.SaveAs("plots2D/"+TString(h_ref->GetName())+".png");
  }
  
}
