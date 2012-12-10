#include "TH1F.h"
#include <iostream>

void printBins(TH1F* h) {
  cout << h->GetName() << endl;
  for (int bin=1;bin<=h->GetNbinsX();++bin) {
    cout << bin << " " << h->GetBinContent(bin) << endl;
  }
}

float ratioPoissErr(float nval, float nerr, float dval, float derr) {
  return sqrt( pow(nerr/dval,2) + pow(nval*derr/pow(dval,2),2) );
}

float efficiencyErr(float eff, float den) {
  if (eff<0 || eff>1 || den<=0) return 0;
  return sqrt( eff*(1-eff)/den );
}

void avoidNegativeBins(TH1F* h) {
  for (int bin=1;bin<=h->GetNbinsX();++bin) {
    if (h->GetBinContent(bin)<0) h->SetBinContent(bin,0);
  }
}

void fillDownMirrorUp(TH1F* central,TH1F* up,TH1F* down) {
  down->Add(up);
  down->Scale(-1);
  down->Add(central);
  down->Add(central);
  //need to avoid negative values...
  avoidNegativeBins(down);
}

void divideHisto(TH1F* num,TH1F* den){
  for (int bin=1;bin<=num->GetNbinsX();++bin) {
    if (fabs(den->GetBinContent(bin)>0) ) {
      //order is important here...
      num->SetBinError(bin,efficiencyErr(num->GetBinContent(bin)/den->GetBinContent(bin),den->GetBinContent(bin)));
      num->SetBinContent(bin,num->GetBinContent(bin)/den->GetBinContent(bin));
    } else num->SetBinContent(bin,0);
  }
}

void divideHistoProtected(TH1F* num,TH1F* den, float maxr=50.){
  for (int bin=1;bin<=num->GetNbinsX();++bin) {
    if (fabs(den->GetBinContent(bin)>0) ) {
      //fixme: avoid large flututaions
      if (num->GetBinContent(bin)/den->GetBinContent(bin)<maxr)
	num->SetBinContent(bin,num->GetBinContent(bin)/den->GetBinContent(bin));
      else {
	cout << "divideHisto - WARNING: setting bin to zero because dividing by a very small number - bin: " << bin <<  " - num " << num->GetName() << " " << num->GetBinContent(bin) 
	     << " den " << den->GetName() << " " << den->GetBinContent(bin) << endl;
	num->SetBinContent(bin,0);
      }
    } else num->SetBinContent(bin,0);
  }
}

void multiplyHisto(TH1F* num,TH1F* den){
  for (int bin=1;bin<=num->GetNbinsX();++bin) {
    num->SetBinContent(bin,num->GetBinContent(bin)*den->GetBinContent(bin));
  }
}

void scaleIntegral(TH1F* central,TH1F* other) {
  if (other->Integral()>0) other->Scale(central->Integral()/other->Integral());
}

void overFlowInLastBin(TH1F* h) {
  h->SetBinContent(h->GetNbinsX(),h->GetBinContent(h->GetNbinsX())+h->GetBinContent(h->GetNbinsX()+1));
  h->SetBinError(h->GetNbinsX(),sqrt(pow(h->GetBinError(h->GetNbinsX()),2)+pow(h->GetBinError(h->GetNbinsX()+1),2)));
  h->SetBinContent(h->GetNbinsX()+1,0.);
  h->SetBinError(h->GetNbinsX()+1,0.);
}

void makeVetoEfficHisto(TH1F* h, TH1F* i) {
  //inlcude overflow bin but not underflow
  float tot_integ = h->Integral(1,h->GetNbinsX()+1);
  if (tot_integ<=0) {
    cout << "INTEGRAL IS NOT POSITIVE!" << endl;
    return;
  }
  for (int bin=1;bin<=h->GetNbinsX();++bin) {
    float integ =  h->Integral(bin,h->GetNbinsX()+1);
    i->SetBinError(bin,efficiencyErr(integ/tot_integ,tot_integ));
    i->SetBinContent(bin,1.-integ/tot_integ);
  }
}
