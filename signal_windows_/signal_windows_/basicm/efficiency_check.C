#include "../Style/tdrstyle.C"
#include "../Style/CMS_lumi.C"

#include <TString.h>

void efficiency_check(){

 setTDRStyle();

 writeExtraText = true;
 extraText  = "Phase-2 simulation";

 TFile* histFile = new TFile("./Eff_nopu.root");

 TH1D* l1eg = (TH1D*)histFile->Get("l1eg_eta");
 TH1D* bpix1234 = (TH1D*)histFile->Get("bpix1234_eta");
 TH1D* bpix123d1 = (TH1D*)histFile->Get("bpix123d1_eta");
 TH1D* bpix12d12 = (TH1D*)histFile->Get("bpix12d12_eta");
 TH1D* bpix1d123 = (TH1D*)histFile->Get("bpix1d123_eta");
 TH1D* fpix1234 = (TH1D*)histFile->Get("fpix1234_eta");
 TH1D* fpix2345 = (TH1D*)histFile->Get("fpix2345_eta");

 TGraphAsymmErrors* bpix1234_eff = new TGraphAsymmErrors(bpix1234, l1eg,"B");
 TGraphAsymmErrors* bpix123d1_eff = new TGraphAsymmErrors(bpix123d1, l1eg,"B");
 TGraphAsymmErrors* bpix12d12_eff = new TGraphAsymmErrors(bpix12d12, l1eg,"B");
 TGraphAsymmErrors* bpix1d123_eff = new TGraphAsymmErrors(bpix1d123, l1eg,"B");
 TGraphAsymmErrors* fpix1234_eff = new TGraphAsymmErrors(fpix1234, l1eg,"B");
 TGraphAsymmErrors* fpix2345_eff = new TGraphAsymmErrors(fpix2345, l1eg,"B");

 TCanvas *c1 = new TCanvas("c1","c1",1800,1200);
 gStyle->SetOptStat(0);
 gStyle->SetPalette(1);
 c1->SetTopMargin(0.07);
 c1->SetRightMargin(0.05);
 c1->SetGridy();
 c1->SetGridx();
 TGaxis::SetMaxDigits(3);
 c1->cd(); 

 bpix1234_eff->GetXaxis()->SetTitle("#eta_{gen}");
 bpix1234_eff->GetYaxis()->SetTitle("Efficiency");
 bpix1234_eff->GetXaxis()->SetRangeUser(0.,3.);
 bpix1234_eff->SetMarkerColor(1);
 bpix1234_eff->SetLineColor(1);
 bpix1234_eff->SetLineWidth(1);
 bpix1234_eff->SetMarkerStyle(20);
 bpix1234_eff->SetMarkerSize(1.);
 bpix1234_eff->Draw("aplez");

 bpix123d1_eff->SetMarkerColor(2);
 bpix123d1_eff->SetLineColor(2);
 bpix123d1_eff->SetLineWidth(1);
 bpix123d1_eff->SetMarkerStyle(20);
 bpix123d1_eff->SetMarkerSize(1.);
 bpix123d1_eff->Draw("sameplez");

 bpix12d12_eff->SetMarkerColor(3);
 bpix12d12_eff->SetLineColor(3);
 bpix12d12_eff->SetLineWidth(1);
 bpix12d12_eff->SetMarkerStyle(20);
 bpix12d12_eff->SetMarkerSize(1.);
 bpix12d12_eff->Draw("sameplez");

 bpix1d123_eff->SetMarkerColor(4);
 bpix1d123_eff->SetLineColor(4);
 bpix1d123_eff->SetLineWidth(1);
 bpix1d123_eff->SetMarkerStyle(20);
 bpix1d123_eff->SetMarkerSize(1.);
 bpix1d123_eff->Draw("sameplez");

 fpix1234_eff->SetMarkerColor(5);
 fpix1234_eff->SetLineColor(5);
 fpix1234_eff->SetLineWidth(1);
 fpix1234_eff->SetMarkerStyle(20);
 fpix1234_eff->SetMarkerSize(1.);
 fpix1234_eff->Draw("sameplez");

 fpix2345_eff->SetMarkerColor(6);
 fpix2345_eff->SetLineColor(6);
 fpix2345_eff->SetLineWidth(1);
 fpix2345_eff->SetMarkerStyle(20);
 fpix2345_eff->SetMarkerSize(1.);
 fpix2345_eff->Draw("sameplez");

 TLegend *Lgd = new TLegend(0.2, 0.3, 0.4, 0.7);

 Lgd->SetFillColor(0);
 Lgd->SetTextFont(42);
 Lgd->SetTextSize(0.03);
 Lgd->SetBorderSize(0);
 Lgd->SetFillStyle(0);
 //Lgd->AddEntry(hEG,"Phase-2 L1 EG(barrel ECAL/ HGCAL)","lp");
 Lgd->AddEntry(bpix1234_eff,"Layer1234","lp");
 Lgd->AddEntry(bpix123d1_eff,"Layer123 Disk1","lp");
 Lgd->AddEntry(bpix12d12_eff,"Layer12 Disk12","lp");
 Lgd->AddEntry(bpix1d123_eff,"Layer1 Disk123","lp");
 Lgd->AddEntry(fpix1234_eff,"Disk1234","lp");
 Lgd->AddEntry(fpix2345_eff,"Disk2345","lp");
 Lgd->Draw();

 CMS_lumi( c1, 4, 0 );
 c1->Update();

 c1->SaveAs("effcheck.png");


}
