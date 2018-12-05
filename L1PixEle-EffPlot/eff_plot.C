#define eff_plot_cxx
#include "eff_plot.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void eff_plot::Loop()
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   
   int nbins = 60; float x1 = -3.0; float x2 = 3.0;

   TH1F* hEG_denom = new TH1F("hEG_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hEG_nom = new TH1F("hEG_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hEG_nom->Sumw2();
   hEG_denom->Sumw2();

   TH1F* hTrack_denom = new TH1F("hTrack_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hTrack_nom = new TH1F("hTrack_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hTrack_nom->Sumw2();
   hTrack_denom->Sumw2();

   TH1F* hTrackLoose_denom = new TH1F("hTrackLoose_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hTrackLoose_nom = new TH1F("hTrackLoose_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hTrackLoose_nom->Sumw2();
   hTrackLoose_denom->Sumw2();

   TH1F* hPix_denom = new TH1F("hPix_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hPix_nom = new TH1F("hPix_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hPix_nom->Sumw2();
   hPix_denom->Sumw2();
   
   TH1F* hPix_iso_denom = new TH1F("hPix_iso_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hPix_iso_nom = new TH1F("hPix_iso_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hPix_iso_nom->Sumw2();
   hPix_iso_denom->Sumw2();

   TH1F* hpixmatching_denom = new TH1F("hpixmatching_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hpixmatching_nom = new TH1F("hpixmatching_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hpixmatching_nom->Sumw2();
   hpixmatching_denom->Sumw2();

   TH1F* hTrk_denom = new TH1F("hTrk_denom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   TH1F* hTrk_nom = new TH1F("hTrk_nom",";#eta of generated electron; Efficiency",nbins,x1,x2);
   hTrk_nom->Sumw2();
   hTrk_denom->Sumw2();
   
   int width_bit[27] = {0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x80, 0x100, 0x200, 0x400, 0x800, 0x1000, 0x2000, 0x4000, 0x8000, 0x10000, 0x20000, 0x40000, 0x80000,
                        0x100000, 0x200000, 0x400000, 0x800000, 0x1000000, 0x2000000, 0x4000000};
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      if(nt_genPt < 20.) continue;
      if( fabs(nt_genEta) > 3.0 ) continue;
   
      int eg_size = ntEgEt->size();
      int gen_matched_eg = -1;
      float dr = 999.;


      // find the closest L1 egamma to gen electron
      for(int i = 0; i < eg_size; i++){
         float temp_dr = sqrt(pow(nt_genPhi-ntEgPhi->at(i),2)+pow(nt_genEta-ntEgEta->at(i),2));
         if( temp_dr < dr){
           dr = temp_dr;
           gen_matched_eg = i;
         }
      }// eg loop
   
      float pt_err = 0.; 
      if(gen_matched_eg != -1) {pt_err = fabs(nt_genPt-ntEgEt->at(gen_matched_eg))/nt_genPt; }
      if( gen_matched_eg != -1 && pt_err < 0.5){

        if(ntL1TkLooseEgEt->size() != 0){ 
           for(int i = 0; i < ntL1TkLooseEgEt->size(); i++){
              float temp_EgEt = ntL1TkLooseEgEt->at(i);
              float temp_EgEta = ntL1TkLooseEgEta->at(i);
              float temp_EgPhi = ntL1TkLooseEgPhi->at(i);
              if( ntEgEt->at(gen_matched_eg) == temp_EgEt && ntEgEta->at(gen_matched_eg) == temp_EgEta && ntEgPhi->at(gen_matched_eg) == temp_EgPhi) hTrackLoose_nom->Fill(nt_genEta, 1.); 
           }

        }
        if(ntL1TkEgEt->size() !=0){
           for(int i = 0; i < ntL1TkEgEt->size(); i++){
              float temp_EgEt = ntL1TkEgEt->at(i);
              float temp_EgEta = ntL1TkEgEta->at(i);
              float temp_EgPhi = ntL1TkEgPhi->at(i);
              if( ntEgEt->at(gen_matched_eg) == temp_EgEt && ntEgEta->at(gen_matched_eg) == temp_EgEta && ntEgPhi->at(gen_matched_eg) == temp_EgPhi) hTrack_nom->Fill(nt_genEta, 1.); 
           }
        }
        if((trigger_bit_width->at(gen_matched_eg)&width_bit[0])==width_bit[0]) hPix_nom->Fill(nt_genEta, 1.);
        if((trigger_bit_width->at(gen_matched_eg)&width_bit[0])==width_bit[0]) hpixmatching_nom->Fill(nt_genEta, 1.);
        //if( ntCl_iso_match->at(gen_matched_eg) ) hPix_iso_nom->Fill(nt_genEta, 1.);

        hEG_nom->Fill(nt_genEta, 1.);
        hpixmatching_denom->Fill(nt_genEta, 1.);
      }
   
      hEG_denom->Fill(nt_genEta, 1.);
      hPix_denom->Fill(nt_genEta, 1.);
      hTrack_denom->Fill(nt_genEta, 1.);
      hTrackLoose_denom->Fill(nt_genEta, 1.);
   
   } // event loop

   TGraphAsymmErrors* hEG = new TGraphAsymmErrors(hEG_nom, hEG_denom,"B");
   TGraphAsymmErrors* hPix = new TGraphAsymmErrors(hPix_nom, hPix_denom,"B");
   TGraphAsymmErrors* hTrack = new TGraphAsymmErrors(hTrack_nom, hTrack_denom,"B");
   TGraphAsymmErrors* hTrackLoose = new TGraphAsymmErrors(hTrackLoose_nom, hTrackLoose_denom,"B");
   TGraphAsymmErrors* hPixMatching = new TGraphAsymmErrors(hpixmatching_nom, hpixmatching_denom,"B");

   TGraphAsymmErrors* hPixEff = new TGraphAsymmErrors(hPix_nom, hEG_nom,"B");
   TGraphAsymmErrors* hPixEffIso = new TGraphAsymmErrors(hPix_iso_nom, hEG_nom,"B");
   TGraphAsymmErrors* hTrkEff = new TGraphAsymmErrors(hTrack_nom, hEG_nom,"B");
   TGraphAsymmErrors* hTrklooseEff = new TGraphAsymmErrors(hTrackLoose_nom, hEG_nom,"B");
   
   for(int pointNr=0;pointNr<hPixMatching->GetN();pointNr++){
     hPixMatching->SetPointEXhigh(pointNr,0);
     hPixMatching->SetPointEXlow(pointNr,0);

     hPixMatching->SetPointEYhigh(pointNr,0);
     hPixMatching->SetPointEYlow(pointNr,0);
   }

   //TCanvas *c1 = new TCanvas("c1","c1",800,700);
   TCanvas *c1 = new TCanvas("c1","c1",1200,800);
   gStyle->SetOptStat(0);
   gStyle->SetLineWidth(1); // axis width, default is 1
   c1->SetTopMargin(0.05);
   c1->SetBottomMargin(0.12);
   c1->SetRightMargin(0.03);
   c1->SetLeftMargin(0.15);
   c1->SetGrid();
   c1->SetTicky(1);
   c1->SetTickx(1);
   c1->cd();

   //hPix->SetTitle("CMSSW_10_1_0_pre3, Phase 2, <PU>=0");
   hPix->GetXaxis()->SetTitle("#eta_{gen} ");
   hPix->GetXaxis()->SetTitleOffset(1.);
   hPix->GetXaxis()->SetTitleSize(0.055);
   hPix->GetXaxis()->SetNdivisions(510);
   hPix->GetYaxis()->SetNdivisions(506);
   hPix->GetXaxis()->SetLabelSize(0.05);
   hPix->GetYaxis()->SetLabelSize(0.05);
   //hPix->GetXaxis()->SetRangeUser(-1.7, 1.7);
   hPix->GetXaxis()->SetRangeUser(-3.0, 3.0);
   //hPix->GetXaxis()->SetRangeUser(-2.5, 2.5);
   hPix->GetYaxis()->SetRangeUser(0.0, 1.1);
   //hPix->GetYaxis()->SetRangeUser(0.5, 1.1);
   hPix->GetYaxis()->SetTitle("Efficiency");
   hPix->GetYaxis()->SetTitleOffset(1.2);
   hPix->GetYaxis()->SetTitleSize(0.055);

   hPix->SetMarkerColor(4);
   hPix->SetLineColor(4);
   hPix->SetLineWidth(1);
   hPix->SetMarkerStyle(20);
   hPix->SetMarkerSize(1.);
   hPix->Draw("ape");

   hTrack->SetMarkerColor(1);
   hTrack->SetLineColor(1);
   hTrack->SetLineWidth(1);
   hTrack->SetMarkerStyle(20);
   hTrack->SetMarkerSize(1.);
   hTrack->Draw("pe same");

   hTrackLoose->SetMarkerColor(1);
   hTrackLoose->SetLineColor(1);
   hTrackLoose->SetLineWidth(1);
   hTrackLoose->SetMarkerStyle(24);
   hTrackLoose->SetMarkerSize(1.);
   hTrackLoose->Draw("pe same");

   hEG->SetMarkerColor(2);
   hEG->SetLineColor(2);
   hEG->SetLineWidth(1);
   hEG->SetMarkerStyle(20);
   hEG->SetMarkerSize(1.);
   hEG->Draw("pe same");
   
   hPixMatching->SetMarkerColor(4);
   hPixMatching->SetLineColor(4);
   hPixMatching->SetLineWidth(2);
   hPixMatching->SetLineStyle(2);
   hPixMatching->SetMarkerSize(0);
   hPixMatching->SetMarkerStyle(20);
   hPixMatching->Draw("l same");

   TLegend *Lgd = new TLegend(0.3, 0.15, 0.7, 0.45);

   //Lgd-> SetNColumns(2);
   Lgd->SetFillColor(0);
   Lgd->SetTextFont(42);
   Lgd->SetTextSize(0.035);
   Lgd->SetBorderSize(0);
   Lgd->SetFillStyle(0);
   //Lgd->AddEntry(hEG,"Phase-2 L1 EG(barrel ECAL/ HGCAL)","lp");
   Lgd->AddEntry(hEG,"Phase-2 L1 EG","lp");
   Lgd->AddEntry(hPix,"Phase-2 L1 EG + Pixel","lp");
   Lgd->AddEntry(hTrack,"Phase-2 L1 EG + L1 Track (Iso)","lp");
   Lgd->AddEntry(hTrackLoose,"Phase-2 L1 EG + L1 Track (Iso & loose)","lp");
   Lgd->AddEntry(hPixMatching,"Phase-2 Pixel","l");
   Lgd->Draw();

   float r_ = c1->GetRightMargin();
   float t_ = c1->GetTopMargin();

 //  TLatex t(1-r_,1-t_+0.2*t_,"CMSSW_10_1_7, Phase 2, <PU>=200");
   TLatex t(1-r_,1-t_+0.2*t_,"CMSSW_10_1_7, Phase 2, <PU>=0");
   t.SetNDC();
   t.SetTextFont(42);
   t.SetTextAlign(31);
   t.SetTextSize(0.6*t_);
   t.Draw();

   TLatex pt_cut(1.5,0.2,"p_{T}^{gen} > 20 GeV");
   pt_cut.SetTextSize(0.035);
   pt_cut.Draw();

   c1->Print("Eff.png");
    
   Double_t sum1 = 0.;
   Double_t sum2 = 0.;

   Double_t sum4 = 0.;

   for(Int_t i = 0; i <= 4; i++) {
      // 0 ~ 4 : eta 2.5 ~ 3.0
      // 5 ~ 14: eta 1.5 ~ 2.5
      // 15 ~ 29 : eta 0 ~ 1.5
      
      Double_t x1 = 0., y1 = 0.;
      Double_t x2 = 0., y2 = 0.;

      Double_t x4 = 0., y4 = 0.;

      hEG->GetPoint(i, x1, y1);
      hPix->GetPoint(i, x2, y2);

      hTrack->GetPoint(i, x4, y4);

      sum1 += y1;
      sum2 += y2;

      sum4 += y4;

   //   cout << "ith : " << i << endl;
   //   cout << " x1 : " << x1 << ", x2 : " << x2 << ", x3 : " << x3 << ", x4 : " << x4 << endl;
  //    cout << " y1 : " << y1 << ", y2 : " << y2 << ", y3 : " << y3 << ", y4 : " << y4 << endl;
  //    cout << endl;

   }

   sum1 /= 5.;
   sum2 /= 5.;

   sum4 /= 5.;

   cout << "Avg eff of L1 EG : " << sum1 << endl;
   cout << "Avg eff of L1 EG + pixel : " << sum2 << endl;
 
   cout << "Avg eff of L1 EG + L1 Track : " << sum4 << endl;
      

   /*
   TCanvas *c2 = new TCanvas("c2","c2",1200,800);
   c2->SetTopMargin(0.05);
   c2->SetBottomMargin(0.12);
   c2->SetRightMargin(0.03);
   c2->SetLeftMargin(0.15);
   c2->SetGrid();
   c2->SetTicky(1);
   c2->SetTickx(1);
   c2->cd();

   hPixEff->GetXaxis()->SetTitle("#eta_{gen} ");
   hPixEff->GetXaxis()->SetTitleOffset(1.);
   hPixEff->GetXaxis()->SetTitleSize(0.055);
   hPixEff->GetXaxis()->SetNdivisions(510);
   hPixEff->GetYaxis()->SetNdivisions(506);
   hPixEff->GetXaxis()->SetLabelSize(0.05);
   hPixEff->GetYaxis()->SetLabelSize(0.05);
   hPixEff->GetXaxis()->SetRangeUser(-3.0, 3.0);
   //hPixEff->GetXaxis()->SetRangeUser(-2.5, 2.5);
   hPixEff->GetYaxis()->SetRangeUser(0., 1.1);
   hPixEff->GetYaxis()->SetTitle("Efficiency");
   hPixEff->GetYaxis()->SetTitleOffset(1.2);
   hPixEff->GetYaxis()->SetTitleSize(0.055);

   hPixEff->SetMarkerColor(4);
   hPixEff->SetLineColor(4);
   hPixEff->SetLineWidth(1);
   hPixEff->SetMarkerStyle(20);
   hPixEff->SetMarkerSize(1.);
   hPixEff->Draw("ape");
   
   hPixEffIso->SetMarkerColor(8);
   hPixEffIso->SetLineColor(8);
   hPixEffIso->SetLineWidth(1);
   hPixEffIso->SetMarkerStyle(20);
   hPixEffIso->SetMarkerSize(1.);
   hPixEffIso->Draw("pe same");
   
   hTrkEff->SetMarkerColor(1);
   hTrkEff->SetLineColor(1);
   hTrkEff->SetLineWidth(1);
   hTrkEff->SetMarkerStyle(20);
   hTrkEff->SetMarkerSize(1.);
   hTrkEff->Draw("pe same");

   hTrklooseEff->SetMarkerColor(1);
   hTrklooseEff->SetLineColor(1);
   hTrklooseEff->SetLineWidth(1);
   hTrklooseEff->SetMarkerStyle(24);
   hTrklooseEff->SetMarkerSize(1.);
   hTrklooseEff->Draw("pe same");

   //TLegend *Lgd1 = new TLegend(0.25, 0.9, 0.95, 0.95);
   TLegend *Lgd1 = new TLegend(0.3, 0.15, 0.7, 0.45);

   //Lgd1->SetNColumns(4);
   Lgd1->SetFillColor(0);
   Lgd1->SetTextFont(42);
   Lgd1->SetTextSize(0.03);
   Lgd1->SetBorderSize(0);
   Lgd1->SetFillStyle(0);
   Lgd1->AddEntry(hPixEff,"Pixel matching","lp");
   Lgd1->AddEntry(hPixEffIso,"Pixel mathcing + Track isolation","lp");
   Lgd1->AddEntry(hTrkEff,"L1 Track","lp");
   Lgd1->AddEntry(hTrklooseEff,"L1 Track (loose matching)","lp");
   Lgd1->Draw();

   TLatex t1(-1.8,1.11,"CMS Preliminary Simulation, Phase 2, <PU>=200");
   t1.SetTextSize(0.035);
   t1.Draw();

   TLatex pt_cut1(2.,0.2,"p_{T}^{gen} > 20 GeV");
   pt_cut1.SetTextSize(0.035);
   pt_cut1.Draw();
 
   //c2->Print("Algorithm-Eff.png");
   //c2->Print("Algorithm-Eff.pdf");
   */  


} 
