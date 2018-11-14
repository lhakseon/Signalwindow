#define Make2Dplots_cxx
#include "Make2Dplots.h"
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGaxis.h>

#include <iostream>

//#include "../Style/tdrstyle.C"
#include "../Style/CMS_lumi.C"

double getMedian(const std::vector<float> &vec)
{
  double median=0.;

  //odd number, definate median
  if(vec.size() % 2 !=0) {
    int middleNr = (vec.size()+1)/2;
    median = vec[middleNr];
  }else{ //even number, take median as halfway between the two middle values
    int middleNr = (vec.size()+1)/2;
    median= vec[middleNr];
    if(middleNr+1 <(int) vec.size()) median+= vec[middleNr+1];
    median/=2.;
  }
  return median;
}

int getMedianIndex(const std::vector<float> &vec)
{

  //odd number, definate median
  if(vec.size() % 2 !=0) {
    int middleNr = (vec.size()+1)/2;
    return middleNr;
  }else{ //even number, take median as halfway between the two middle values
    int middleNr = (vec.size()+1)/2;
    return middleNr;
  }
}

void Make2Dplots::Loop(int eta_ = 1)
{
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  const TString eta_region[6] = {"|#eta|<0.8", "0.8<|#eta|<1.4", "1.4<|#eta|<1.7", "1.7<|#eta|<2.1", "2.1<|#eta|<2.7", "2.7<|#eta|<3."};
  
  const int etStep = 1;
  const int binSize = 90;
  float x[binSize], xErr[binSize];

  // eg-pixelpixel deta
  float pix12egDeta[binSize], pix12egDetaErr[binSize];
  float pix13egDeta[binSize], pix13egDetaErr[binSize];
  float pix14egDeta[binSize], pix14egDetaErr[binSize];
  float pix23egDeta[binSize], pix23egDetaErr[binSize];
  float pix24egDeta[binSize], pix24egDetaErr[binSize];
  float pix34egDeta[binSize], pix34egDetaErr[binSize];

  TH2F* pix12egDeta_dist = new TH2F("pix12egDeta_dist","pix12egDeta_dist", 90,10,100,100,-0.2,0.2);
  TH2F* pix13egDeta_dist = new TH2F("pix13egDeta_dist","pix13egDeta_dist", 90,10,100,100,-0.2,0.2);
  TH2F* pix14egDeta_dist = new TH2F("pix14egDeta_dist","pix14egDeta_dist", 90,10,100,100,-0.2,0.2);
  TH2F* pix23egDeta_dist = new TH2F("pix23egDeta_dist","pix23egDeta_dist", 90,10,100,100,-0.2,0.2);
  TH2F* pix24egDeta_dist = new TH2F("pix24egDeta_dist","pix24egDeta_dist", 90,10,100,100,-0.2,0.2);
  TH2F* pix34egDeta_dist = new TH2F("pix34egDeta_dist","pix34egDeta_dist", 90,10,100,100,-0.2,0.2);
  // pixel-pixel deta

  float pix123Deta[binSize], pix123DetaErr[binSize];
  float pix124Deta[binSize], pix124DetaErr[binSize];
  float pix134Deta[binSize], pix134DetaErr[binSize];
  float pix234Deta[binSize], pix234DetaErr[binSize];

  TH2F* pix123Deta_dist = new TH2F("pix123Deta_dist","pix123Deta_dist", 90,10,100,100,-0.02,0.02);
  TH2F* pix124Deta_dist = new TH2F("pix124Deta_dist","pix124Deta_dist", 90,10,100,100,-0.02,0.02);
  TH2F* pix134Deta_dist = new TH2F("pix134Deta_dist","pix134Deta_dist", 90,10,100,100,-0.02,0.02);
  TH2F* pix234Deta_dist = new TH2F("pix234Deta_dist","pix234Deta_dist", 90,10,100,100,-0.02,0.02);

      vector<float*> pixpixegDeta, pixpixegDetaErr;
  vector<float*> pixpixDeta, pixpixDetaErr;

  pixpixegDeta.push_back(pix12egDeta);
  pixpixegDeta.push_back(pix13egDeta);
  pixpixegDeta.push_back(pix14egDeta);
  pixpixegDeta.push_back(pix23egDeta);
  pixpixegDeta.push_back(pix24egDeta);
  pixpixegDeta.push_back(pix34egDeta);

  pixpixegDetaErr.push_back(pix12egDetaErr);
  pixpixegDetaErr.push_back(pix13egDetaErr);
  pixpixegDetaErr.push_back(pix14egDetaErr);
  pixpixegDetaErr.push_back(pix23egDetaErr);
  pixpixegDetaErr.push_back(pix24egDetaErr);
  pixpixegDetaErr.push_back(pix34egDetaErr);

  pixpixDeta.push_back(pix123Deta);
  pixpixDeta.push_back(pix124Deta);
  pixpixDeta.push_back(pix134Deta);
  pixpixDeta.push_back(pix234Deta);

  pixpixDetaErr.push_back(pix123DetaErr);
  pixpixDetaErr.push_back(pix124DetaErr);
  pixpixDetaErr.push_back(pix134DetaErr);
  pixpixDetaErr.push_back(pix234DetaErr);

  for(int nth = 0; nth < binSize; nth++){
  
     Long64_t nbytes = 0, nb = 0;

     vector<float> pix12egDeta_;
     vector<float> pix13egDeta_;
     vector<float> pix14egDeta_;
     vector<float> pix23egDeta_;
     vector<float> pix24egDeta_;
     vector<float> pix34egDeta_;

     vector<float> pix123Deta_;
     vector<float> pix124Deta_;
     vector<float> pix134Deta_;
     vector<float> pix234Deta_;

     for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if( eta_ == 1 && fabs(ntEgEta->at(0)) > 0.8 ) continue; // test with barrel only
        if( eta_ == 2 && (fabs(ntEgEta->at(0)) < 0.8 || fabs(ntEgEta->at(0)) > 1.4)) continue; // test with barrel only
        if( eta_ == 3 && (fabs(ntEgEta->at(0)) < 1.4 || fabs(ntEgEta->at(0)) > 1.7)) continue; // test with barrel only
        if( eta_ == 4 && (fabs(ntEgEta->at(0)) < 1.7 || fabs(ntEgEta->at(0)) > 2.1)) continue; // test with barrel only
        if( eta_ == 5 && (fabs(ntEgEta->at(0)) < 2.1 || fabs(ntEgEta->at(0)) > 2.7)) continue; // test with barrel only
        if( eta_ == 6 && (fabs(ntEgEta->at(0)) < 2.7 || fabs(ntEgEta->at(0)) > 3.)) continue; // test with barrel only
/*
        if( ntEgEt->at(0) > (float) (10 + nth) && ntEgEt->at(0) < (float) (11 + nth) ){

          if(ntPix12EGdeta->size() !=0) pix12egDeta_.push_back(ntPix12EGdeta->at(0));          
          if(ntPix13EGdeta->size() !=0) pix13egDeta_.push_back(ntPix13EGdeta->at(0));          
          if(ntPix14EGdeta->size() !=0) pix14egDeta_.push_back(ntPix14EGdeta->at(0));          
          if(ntPix23EGdeta->size() !=0) pix23egDeta_.push_back(ntPix23EGdeta->at(0));          
          if(ntPix24EGdeta->size() !=0) pix24egDeta_.push_back(ntPix24EGdeta->at(0));          
          if(ntPix34EGdeta->size() !=0) pix34egDeta_.push_back(ntPix34EGdeta->at(0));          

          if(ntPix123deta->size() !=0) pix123Deta_.push_back(ntPix123deta->at(0));          
          if(ntPix124deta->size() !=0) pix124Deta_.push_back(ntPix124deta->at(0));          
          if(ntPix134deta->size() !=0) pix134Deta_.push_back(ntPix134deta->at(0));          
          if(ntPix234deta->size() !=0) pix234Deta_.push_back(ntPix234deta->at(0));          
        }
*/

        if( ntEgEt->at(0) > (float) (10 + nth) && ntEgEt->at(0) < (float) (11 + nth) ){


          for(unsigned long i = 0; i < ntPix12EGdeta->size(); i++) { 
             if(fabs(ntPix12EGdeta->at(i)) > 0.2 ) continue;
             pix12egDeta_.push_back(ntPix12EGdeta->at(i));          
             pix12egDeta_dist->Fill(ntEgEt->at(0), ntPix12EGdeta->at(i));
          }
          for(unsigned long i = 0; i < ntPix13EGdeta->size(); i++) {
             if(fabs(ntPix13EGdeta->at(i)) > 0.2 ) continue;
             pix13egDeta_.push_back(ntPix13EGdeta->at(i));          
             pix13egDeta_dist->Fill(ntEgEt->at(0), ntPix13EGdeta->at(i));
          }
          for(unsigned long i = 0; i < ntPix14EGdeta->size(); i++) {
            if(fabs(ntPix14EGdeta->at(i)) > 0.2 ) continue;
             pix14egDeta_.push_back(ntPix14EGdeta->at(i));          
             pix14egDeta_dist->Fill(ntEgEt->at(0), ntPix14EGdeta->at(i));
          }
          for(unsigned long i = 0; i < ntPix23EGdeta->size(); i++) { 
            if(fabs(ntPix23EGdeta->at(i)) > 0.2 ) continue;
             pix23egDeta_.push_back(ntPix23EGdeta->at(i));          
             pix23egDeta_dist->Fill(ntEgEt->at(0), ntPix23EGdeta->at(i));
          }
          for(unsigned long i = 0; i < ntPix24EGdeta->size(); i++) { 
         if(fabs(ntPix24EGdeta->at(i)) > 0.2 ) continue;
             pix24egDeta_.push_back(ntPix24EGdeta->at(i));          
             pix24egDeta_dist->Fill(ntEgEt->at(0), ntPix24EGdeta->at(i));
          }
          for(unsigned long i = 0; i < ntPix34EGdeta->size(); i++) {
          if(fabs(ntPix34EGdeta->at(i)) > 0.2 ) continue;
             pix34egDeta_.push_back(ntPix34EGdeta->at(i));          
             pix34egDeta_dist->Fill(ntEgEt->at(0), ntPix34EGdeta->at(i));
          }
          for(unsigned long i = 0; i < ntPix123deta->size(); i++) {
             if(fabs(ntPix123deta->at(i)) > 0.015 ) continue;
             pix123Deta_.push_back(ntPix123deta->at(i));          
             pix123Deta_dist->Fill(ntEgEt->at(0), ntPix123deta->at(i));
          }
          for(unsigned long i = 0; i < ntPix124deta->size(); i++) {
             if(fabs(ntPix124deta->at(i)) > 0.015 ) continue;
            pix124Deta_.push_back(ntPix124deta->at(i));          
             pix124Deta_dist->Fill(ntEgEt->at(0), ntPix124deta->at(i));
          }
          for(unsigned long i = 0; i < ntPix134deta->size(); i++) {
             if(fabs(ntPix134deta->at(i)) > 0.015 ) continue;
             pix134Deta_.push_back(ntPix134deta->at(i));          
             pix134Deta_dist->Fill(ntEgEt->at(0), ntPix134deta->at(i));
          }
          for(unsigned long i = 0; i < ntPix234deta->size(); i++) {
             if(fabs(ntPix234deta->at(i)) > 0.015 ) continue;
             pix234Deta_.push_back(ntPix234deta->at(i));          
             pix234Deta_dist->Fill(ntEgEt->at(0), ntPix234deta->at(i));
//  cout<<"alp"<<endl;
	  }
        }
     }// event loop

     std::sort (pix12egDeta_.begin(), pix12egDeta_.end()); 
     std::sort (pix13egDeta_.begin(), pix13egDeta_.end()); 
     std::sort (pix14egDeta_.begin(), pix14egDeta_.end()); 
     std::sort (pix23egDeta_.begin(), pix23egDeta_.end()); 
     std::sort (pix24egDeta_.begin(), pix24egDeta_.end()); 
     std::sort (pix34egDeta_.begin(), pix34egDeta_.end()); 

     std::sort (pix123Deta_.begin(), pix123Deta_.end()); 
     std::sort (pix124Deta_.begin(), pix124Deta_.end()); 
     std::sort (pix134Deta_.begin(), pix134Deta_.end()); 
     std::sort (pix234Deta_.begin(), pix234Deta_.end()); 

     int pix12egDeta_size = pix12egDeta_.size();
     int pix13egDeta_size = pix13egDeta_.size();
     int pix14egDeta_size = pix14egDeta_.size();
     int pix23egDeta_size = pix23egDeta_.size();
     int pix24egDeta_size = pix24egDeta_.size();
     int pix34egDeta_size = pix34egDeta_.size();

     int pix123Deta_size = pix123Deta_.size();
     int pix124Deta_size = pix124Deta_.size();
     int pix134Deta_size = pix134Deta_.size();
     int pix234Deta_size = pix234Deta_.size();

     // get index for median

     int pix12egDeta_medianIndex = getMedianIndex(pix12egDeta_);
     int pix13egDeta_medianIndex = getMedianIndex(pix13egDeta_);
     int pix14egDeta_medianIndex = getMedianIndex(pix14egDeta_);
     int pix23egDeta_medianIndex = getMedianIndex(pix23egDeta_);
     int pix24egDeta_medianIndex = getMedianIndex(pix24egDeta_);
     int pix34egDeta_medianIndex = getMedianIndex(pix34egDeta_);

     int pix123Deta_medianIndex = getMedianIndex(pix123Deta_);
     int pix124Deta_medianIndex = getMedianIndex(pix124Deta_);
     int pix134Deta_medianIndex = getMedianIndex(pix134Deta_);
     int pix234Deta_medianIndex = getMedianIndex(pix234Deta_);

     // get median

     float pix12egDeta_median = getMedian(pix12egDeta_);
     float pix13egDeta_median = getMedian(pix13egDeta_);
     float pix14egDeta_median = getMedian(pix14egDeta_);
     float pix23egDeta_median = getMedian(pix23egDeta_);
     float pix24egDeta_median = getMedian(pix24egDeta_);
     float pix34egDeta_median = getMedian(pix34egDeta_);

     float pix123Deta_median = getMedian(pix123Deta_);
     float pix124Deta_median = getMedian(pix124Deta_);
     float pix134Deta_median = getMedian(pix134Deta_);
     float pix234Deta_median = getMedian(pix234Deta_);

     // calculate error

     int pix12egDeta_low = (int)(pix12egDeta_medianIndex - (0.668 * pix12egDeta_medianIndex));
     int pix13egDeta_low = (int)(pix13egDeta_medianIndex - (0.668 * pix13egDeta_medianIndex));
     int pix14egDeta_low = (int)(pix14egDeta_medianIndex - (0.668 * pix14egDeta_medianIndex));
     int pix23egDeta_low = (int)(pix23egDeta_medianIndex - (0.668 * pix23egDeta_medianIndex));
     int pix24egDeta_low = (int)(pix24egDeta_medianIndex - (0.668 * pix24egDeta_medianIndex));
     int pix34egDeta_low = (int)(pix34egDeta_medianIndex - (0.668 * pix34egDeta_medianIndex));

     int pix12egDeta_high =(int)(pix12egDeta_medianIndex + (0.668 * (pix12egDeta_size-1-pix12egDeta_medianIndex)));
     int pix13egDeta_high =(int)(pix13egDeta_medianIndex + (0.668 * (pix13egDeta_size-1-pix13egDeta_medianIndex)));
     int pix14egDeta_high =(int)(pix14egDeta_medianIndex + (0.668 * (pix14egDeta_size-1-pix14egDeta_medianIndex)));
     int pix23egDeta_high =(int)(pix23egDeta_medianIndex + (0.668 * (pix23egDeta_size-1-pix23egDeta_medianIndex)));
     int pix24egDeta_high =(int)(pix24egDeta_medianIndex + (0.668 * (pix24egDeta_size-1-pix24egDeta_medianIndex)));
     int pix34egDeta_high =(int)(pix34egDeta_medianIndex + (0.668 * (pix34egDeta_size-1-pix34egDeta_medianIndex)));

     int pix123Deta_low = (int)(pix123Deta_medianIndex - (0.668 * pix123Deta_medianIndex));
     int pix124Deta_low = (int)(pix124Deta_medianIndex - (0.668 * pix124Deta_medianIndex));
     int pix134Deta_low = (int)(pix134Deta_medianIndex - (0.668 * pix134Deta_medianIndex));
     int pix234Deta_low = (int)(pix234Deta_medianIndex - (0.668 * pix234Deta_medianIndex));

     int pix123Deta_high = (int)(pix123Deta_medianIndex + (0.668 * (pix123Deta_size-1-pix123Deta_medianIndex)));
     int pix124Deta_high = (int)(pix124Deta_medianIndex + (0.668 * (pix124Deta_size-1-pix124Deta_medianIndex)));
     int pix134Deta_high = (int)(pix134Deta_medianIndex + (0.668 * (pix134Deta_size-1-pix134Deta_medianIndex)));
     int pix234Deta_high = (int)(pix234Deta_medianIndex + (0.668 * (pix234Deta_size-1-pix234Deta_medianIndex)));

     float pix12egDeta_medianErr = pix12egDeta_.at(pix12egDeta_high) - pix12egDeta_.at(pix12egDeta_low);
     float pix13egDeta_medianErr = pix13egDeta_.at(pix13egDeta_high) - pix13egDeta_.at(pix13egDeta_low);
     float pix14egDeta_medianErr = pix14egDeta_.at(pix14egDeta_high) - pix14egDeta_.at(pix14egDeta_low);
     float pix23egDeta_medianErr = pix23egDeta_.at(pix23egDeta_high) - pix23egDeta_.at(pix23egDeta_low);
     float pix24egDeta_medianErr = pix24egDeta_.at(pix24egDeta_high) - pix24egDeta_.at(pix24egDeta_low);
     float pix34egDeta_medianErr = pix34egDeta_.at(pix34egDeta_high) - pix34egDeta_.at(pix34egDeta_low);

     float pix123Deta_medianErr = pix123Deta_.at(pix123Deta_high) - pix123Deta_.at(pix123Deta_low);
     float pix124Deta_medianErr = pix124Deta_.at(pix124Deta_high) - pix124Deta_.at(pix124Deta_low);
     float pix134Deta_medianErr = pix134Deta_.at(pix134Deta_high) - pix134Deta_.at(pix134Deta_low);
     float pix234Deta_medianErr = pix234Deta_.at(pix234Deta_high) - pix234Deta_.at(pix234Deta_low);

     x[nth] = 10.5 + float(nth);
     xErr[nth] = 0.;

     pix12egDeta[nth] = pix12egDeta_median;
     pix13egDeta[nth] = pix13egDeta_median;
     pix14egDeta[nth] = pix14egDeta_median;
     pix23egDeta[nth] = pix23egDeta_median;
     pix24egDeta[nth] = pix24egDeta_median;
     pix34egDeta[nth] = pix34egDeta_median;

     pix12egDetaErr[nth] = 0.;
     pix13egDetaErr[nth] = 0.;
     pix14egDetaErr[nth] = 0.;
     pix23egDetaErr[nth] = 0.;
     pix24egDetaErr[nth] = 0.;
     pix34egDetaErr[nth] = 0.;

     pix123Deta[nth] = pix123Deta_median;
     pix124Deta[nth] = pix124Deta_median;
     pix134Deta[nth] = pix134Deta_median;
     pix234Deta[nth] = pix234Deta_median;

     pix123DetaErr[nth] = 0.;
     pix124DetaErr[nth] = 0.;
     pix134DetaErr[nth] = 0.;
     pix234DetaErr[nth] = 0.;

  }// Et scanning loop


  const TString eg_pixpixSW[6] = {"12","13","14","23","24","34"};
  // draw eg-pixpix deta 
  for(int i = 0; i < 6; i++){

     setTDRStyle();
     writeExtraText = true;
     extraText  = "Phase-2 simulation";
     gROOT->SetBatch();
     TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
     gStyle->SetOptStat(0);
     gStyle->SetPalette(1);
     c1->SetTopMargin(0.07);
     c1->SetGridy();
     c1->SetGridx();
     TGaxis::SetMaxDigits(3);

//
     TGraphErrors* gr1 = new TGraphErrors(binSize,x,pixpixegDeta.at(i),xErr,pixpixegDetaErr.at(i)); //0.015 0.01
//     TGraphErrors* gr1 = new TGraphErrors(binSize,x,pixpixegDeta.at(i),xErr,pixpixegDetaErr.at(i));

     TH2F* scatter_plot = NULL;
     if(i==0) scatter_plot = pix12egDeta_dist;
     if(i==1) scatter_plot = pix13egDeta_dist;
     if(i==2) scatter_plot = pix14egDeta_dist;
     if(i==3) scatter_plot = pix23egDeta_dist;
     if(i==4) scatter_plot = pix24egDeta_dist;
     if(i==5) scatter_plot = pix34egDeta_dist;

     if(scatter_plot != NULL){
       scatter_plot->GetYaxis()->SetRangeUser(-0.2,0.2);
       scatter_plot->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
       scatter_plot->GetYaxis()->SetTitle("#Delta#eta [rad.]");
       scatter_plot->SetMarkerColor(kRed);
       //scatter_plot->SetMarkerSize(.25);
       scatter_plot->SetMarkerStyle(6);
       scatter_plot->Draw("scat=1.");

     }

     TString nthsw;
     nthsw.Form("%d", i+1);
//     gr1->SetTitle("#Delta#eta(L1 EG, pixel-pixel "+nthsw+")");
     gr1->SetTitle("#Delta#eta(L1 EG, pixel-pixel "+nthsw+")");
     gr1->GetYaxis()->SetRangeUser(-0.2,0.2);
     gr1->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
     gr1->GetYaxis()->SetTitle("#Delta#eta");
     gr1->SetLineColor(kBlue);
     gr1->SetMarkerColor(kBlue);
     gr1->SetMarkerSize(1.);
     gr1->SetMarkerStyle(20);
     gr1->SetFillColorAlpha(kBlue, 0.40);
     gr1->SetLineWidth(1);
     gr1->Draw("p");
//scatter_plot->Draw("p");
     TF1 *median_fitFunc = new TF1("func","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])+[4]) )", 10., 100.);
     median_fitFunc->SetLineColor(kGreen);
     median_fitFunc->SetLineStyle(1);
     median_fitFunc->SetLineWidth(2);

     gr1->Fit(median_fitFunc,"0");

     double ab[4]={0};
     median_fitFunc->GetParameters(ab);
     cout<<ab[0]*pow(2,0) + ab[1]*pow(2,ab[2])*exp(-pow(2,ab[3]))<<endl;


double x1[200],y1[200],x2[200],y2[200];
for ( int j = 10;j<120;j++){
    x1[j]=j;
    x2[j]=0.5;
    y1[j]=ab[0]*pow(j,0) + ab[1]*pow(j,ab[2])*exp(-pow(j,ab[3]));
    y2[j]=0.015;
}
    TGraphErrors* gr2 = new TGraphErrors(110,x1,y1,x2,y2); //0.0017 0.003
    
     gr2->SetLineColor(kBlue);
     gr2->SetFillColorAlpha(kBlue, 0.40);
//    
     gr2->Draw("sameE3");
     median_fitFunc->Draw("lsame");

     TLegend* leg = new TLegend(0.3,0.2,0.5,0.4);
     leg->AddEntry(scatter_plot, "Scattered #Delta#eta distribution as a function of L1 EG E_{T}","pl");
 
     leg->AddEntry(gr1, "Median point of #Delta#eta","lp");
     leg->AddEntry(median_fitFunc,"Fit of Median point", "l");
     leg->AddEntry(gr2,"Signal window", "f");
     leg->AddEntry((TObject*)0, eta_region[eta_-1], "");
     leg->SetTextFont(42);
     leg->SetTextSize(0.035);
     leg->SetBorderSize(0);
     leg->SetFillStyle(0);
     leg->Draw();

     CMS_lumi( c1, 4, 0 );
     c1->Update();

     c1->SaveAs("PixPixEG"+nthsw+"_"+eta_region[eta_-1]+"_eta.png");


     delete gr2;
     delete gr1;
     delete c1;
  }

  const TString pixpix1_SW[4] = {"12","12","13","23"};
  const TString pixpix2_SW[4] = {"23","24","34","34"};
  // draw pix-pixpix deta 
  for(int i = 0; i < 4; i++){

     setTDRStyle();
     writeExtraText = true;
     extraText  = "Phase-2 simulation";
     gROOT->SetBatch();
     TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
     gStyle->SetOptStat(0);
     gStyle->SetPalette(1);
     c1->SetTopMargin(0.07);
     c1->SetGridy();
     c1->SetGridx();
     TGaxis::SetMaxDigits(4);

     //  
     TGraphErrors* gr1 = new TGraphErrors(binSize,x,pixpixDeta.at(i),xErr,pixpixDetaErr.at(i)); //0.0017 0.003
     TH2F* scatter_plot = NULL;
     if(i==0) scatter_plot = pix123Deta_dist;
     if(i==1) scatter_plot = pix124Deta_dist;
     if(i==2) scatter_plot = pix134Deta_dist;
     if(i==3) scatter_plot = pix234Deta_dist;

     if(scatter_plot != NULL){

       scatter_plot->GetYaxis()->SetRangeUser(-0.015,0.015);
       scatter_plot->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
       scatter_plot->GetYaxis()->SetTitle("#Delta#eta [rad.]");
       scatter_plot->SetMarkerColor(kRed);
       //scatter_plot->SetMarkerSize(.25);
       scatter_plot->SetMarkerStyle(6);
       scatter_plot->Draw("scat=1.");
     }


     TString nthsw;
     nthsw.Form("%d", i+7);
     gr1->SetTitle("#Delta#eta(pixel-pixel, pixel-pixel "+nthsw+")");
     gr1->GetYaxis()->SetRangeUser(-0.015,0.015);
     gr1->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
     gr1->GetYaxis()->SetTitle("#Delta#eta");
     gr1->SetLineColor(kBlue);
     gr1->SetMarkerColor(kBlue);
     gr1->SetMarkerSize(1.);
     gr1->SetMarkerStyle(20);
 //    gr1->SetFillColorAlpha(kBlue, 0.20);
     gr1->SetLineWidth(1);
     gr1->Draw("p");

     TF1 *median_fitFunc = new TF1("func","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])) )", 10., 100.);

     median_fitFunc->SetParLimits(0, -0.05, 0.001);
     median_fitFunc->SetParLimits(1, -0.001, 5);
     median_fitFunc->SetParLimits(2, -2., 0.2);
     median_fitFunc->SetParLimits(3, -0.7, 0.5);

     median_fitFunc->SetLineColor(kGreen);
     median_fitFunc->SetLineStyle(1);
     median_fitFunc->SetLineWidth(2);

     gr1->Fit(median_fitFunc,"0");
     median_fitFunc->Draw("lsame");

     double ab[4]={0};
     median_fitFunc->GetParameters(ab);


double x1[200],y1[200],x2[200],y2[200];
for ( int j = 10;j<120;j++){
    x1[j]=j;
    x2[j]=0.5;
    y1[j]=ab[0]*pow(j,0) + ab[1]*pow(j,ab[2])*exp(-pow(j,ab[3]));
    y2[j]=0.003;
}
    TGraphErrors* gr2 = new TGraphErrors(110,x1,y1,x2,y2); //0.0017 0.003
    
     gr2->SetFillColorAlpha(kBlue, 0.40);
//    
     gr2->Draw("sameE3");
     TLegend* leg = new TLegend(0.2,0.2,0.4,0.4);
     leg->AddEntry(scatter_plot, "Scattered #Delta#eta distribution as a function of L1 EG E_{T}","pl");
     leg->AddEntry(gr1, "Median point of #Delta#eta","lp");
  //   leg->AddEntry(gr1, "Median point of #Delta#eta (pixel segment " + pixpix1_SW[i] + ", pixel segment " + pixpix2_SW[i] + ")", "pe");
     leg->AddEntry(median_fitFunc,"Fit of median point", "l");
     leg->AddEntry(gr2,"Signal window", "f");
     leg->AddEntry((TObject*)0, eta_region[eta_-1], "");
     leg->SetTextFont(42);
     leg->SetTextSize(0.035);
     leg->SetBorderSize(0);
     leg->SetFillStyle(0);
     leg->Draw();

     CMS_lumi( c1, 4, 0 );
     c1->Update();

     c1->SaveAs("PixPix"+nthsw+"_"+eta_region[eta_-1]+"_eta.png");

     delete gr2;
     delete gr1;
     delete c1; 
  }
 delete pix12egDeta_dist;
 delete pix13egDeta_dist;
 delete pix14egDeta_dist;
 delete pix23egDeta_dist;
 delete pix24egDeta_dist;
 delete pix34egDeta_dist;

 delete pix123Deta_dist;
 delete pix124Deta_dist;
 delete pix134Deta_dist;
 delete pix234Deta_dist;
 std::cout << "DONE" << std::endl;
}
