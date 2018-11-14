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
#include <fstream>

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

  const TString eta_region[6] = {"|#eta|<0.8", "0.8<|#eta|<1.4", "1.4<|#eta|<1.7", "1.7<|#eta|<2.1", "2.1<|#eta|<2.7", "2.7<|#eta|<3.0"};
  
  const int etStep = 1;
  const int binSize = 90;
  float x[binSize], xErr[binSize];

  // eg-pixel dphi
  float pix1egDphi[binSize], pix1egDphiErr[binSize];
  float pix2egDphi[binSize], pix2egDphiErr[binSize];
  float pix3egDphi[binSize], pix3egDphiErr[binSize];
  float pix4egDphi[binSize], pix4egDphiErr[binSize];

  TH2F* pix1egDphi_dist = new TH2F("pix1egDphi_dist","pix1egDphi_dist", 90,10,100,100,-0.2,0.2);
  TH2F* pix2egDphi_dist = new TH2F("pix2egDphi_dist","pix2egDphi_dist", 90,10,100,100,-0.2,0.2);
  TH2F* pix3egDphi_dist = new TH2F("pix3egDphi_dist","pix3egDphi_dist", 90,10,100,100,-0.2,0.2);
  TH2F* pix4egDphi_dist = new TH2F("pix4egDphi_dist","pix4egDphi_dist", 90,10,100,100,-0.2,0.2);

  // eg-pixelpixel dphi
  float pix12egDphi[binSize], pix12egDphiErr[binSize];
  float pix13egDphi[binSize], pix13egDphiErr[binSize];
  float pix14egDphi[binSize], pix14egDphiErr[binSize];
  float pix23egDphi[binSize], pix23egDphiErr[binSize];
  float pix24egDphi[binSize], pix24egDphiErr[binSize];
  float pix34egDphi[binSize], pix34egDphiErr[binSize];

  TH2F* pix12egDphi_dist = new TH2F("pix12egDphi_dist","pix12egDphi_dist", 90,10,100,100,-0.2,0.2);
  TH2F* pix13egDphi_dist = new TH2F("pix13egDphi_dist","pix13egDphi_dist", 90,10,100,100,-0.2,0.2);
  TH2F* pix14egDphi_dist = new TH2F("pix14egDphi_dist","pix14egDphi_dist", 90,10,100,100,-0.2,0.2);
  TH2F* pix23egDphi_dist = new TH2F("pix23egDphi_dist","pix23egDphi_dist", 90,10,100,100,-0.2,0.2);
  TH2F* pix24egDphi_dist = new TH2F("pix24egDphi_dist","pix24egDphi_dist", 90,10,100,100,-0.2,0.2);
  TH2F* pix34egDphi_dist = new TH2F("pix34egDphi_dist","pix34egDphi_dist", 90,10,100,100,-0.2,0.2);

  // pixel-pixel dphi
  float pix012Dphi[binSize], pix012DphiErr[binSize];
  float pix013Dphi[binSize], pix013DphiErr[binSize];
  float pix014Dphi[binSize], pix014DphiErr[binSize];
  float pix023Dphi[binSize], pix023DphiErr[binSize];
  float pix024Dphi[binSize], pix024DphiErr[binSize];
  float pix034Dphi[binSize], pix034DphiErr[binSize];
  float pix123Dphi[binSize], pix123DphiErr[binSize];
  float pix124Dphi[binSize], pix124DphiErr[binSize];
  float pix134Dphi[binSize], pix134DphiErr[binSize];
  float pix234Dphi[binSize], pix234DphiErr[binSize];

  TH2F* pix012Dphi_dist = new TH2F("pix012Dphi_dist","pix012Dphi_dist", 90,10,100,100,-0.02,0.02);
  TH2F* pix013Dphi_dist = new TH2F("pix013Dphi_dist","pix013Dphi_dist", 90,10,100,100,-0.02,0.02);
  TH2F* pix014Dphi_dist = new TH2F("pix014Dphi_dist","pix014Dphi_dist", 90,10,100,100,-0.02,0.02);
  TH2F* pix023Dphi_dist = new TH2F("pix023Dphi_dist","pix023Dphi_dist", 90,10,100,100,-0.02,0.02);
  TH2F* pix024Dphi_dist = new TH2F("pix024Dphi_dist","pix024Dphi_dist", 90,10,100,100,-0.02,0.02);
  TH2F* pix034Dphi_dist = new TH2F("pix034Dphi_dist","pix034Dphi_dist", 90,10,100,100,-0.02,0.02);
  TH2F* pix123Dphi_dist = new TH2F("pix123Dphi_dist","pix123Dphi_dist", 90,10,100,100,-0.02,0.02);
  TH2F* pix124Dphi_dist = new TH2F("pix124Dphi_dist","pix124Dphi_dist", 90,10,100,100,-0.02,0.02);
  TH2F* pix134Dphi_dist = new TH2F("pix134Dphi_dist","pix134Dphi_dist", 90,10,100,100,-0.02,0.02);
  TH2F* pix234Dphi_dist = new TH2F("pix234Dphi_dist","pix234Dphi_dist", 90,10,100,100,-0.02,0.02);

  vector<float*> pixegDphi, pixegDphiErr;
  vector<float*> pixpixegDphi, pixpixegDphiErr;
  vector<float*> pixpixDphi, pixpixDphiErr;

  pixegDphi.clear();
  pixegDphiErr.clear();
  pixpixegDphi.clear();
  pixpixegDphiErr.clear();
  pixpixDphi.clear(); 
  pixpixDphiErr.clear();

  pixegDphi.push_back(pix1egDphi);
  pixegDphi.push_back(pix2egDphi);
  pixegDphi.push_back(pix3egDphi);
  pixegDphi.push_back(pix4egDphi);

  pixegDphiErr.push_back(pix1egDphiErr);
  pixegDphiErr.push_back(pix2egDphiErr);
  pixegDphiErr.push_back(pix3egDphiErr);
  pixegDphiErr.push_back(pix4egDphiErr);

  pixpixegDphi.push_back(pix12egDphi);
  pixpixegDphi.push_back(pix13egDphi);
  pixpixegDphi.push_back(pix14egDphi);
  pixpixegDphi.push_back(pix23egDphi);
  pixpixegDphi.push_back(pix24egDphi);
  pixpixegDphi.push_back(pix34egDphi);

  pixpixegDphiErr.push_back(pix12egDphiErr);
  pixpixegDphiErr.push_back(pix13egDphiErr);
  pixpixegDphiErr.push_back(pix14egDphiErr);
  pixpixegDphiErr.push_back(pix23egDphiErr);
  pixpixegDphiErr.push_back(pix24egDphiErr);
  pixpixegDphiErr.push_back(pix34egDphiErr);

  pixpixDphi.push_back(pix012Dphi);
  pixpixDphi.push_back(pix013Dphi);
  pixpixDphi.push_back(pix014Dphi);
  pixpixDphi.push_back(pix023Dphi);
  pixpixDphi.push_back(pix024Dphi);
  pixpixDphi.push_back(pix034Dphi);
  pixpixDphi.push_back(pix123Dphi);
  pixpixDphi.push_back(pix124Dphi);
  pixpixDphi.push_back(pix134Dphi);
  pixpixDphi.push_back(pix234Dphi);

  pixpixDphiErr.push_back(pix012DphiErr);
  pixpixDphiErr.push_back(pix013DphiErr);
  pixpixDphiErr.push_back(pix014DphiErr);
  pixpixDphiErr.push_back(pix023DphiErr);
  pixpixDphiErr.push_back(pix024DphiErr);
  pixpixDphiErr.push_back(pix034DphiErr);
  pixpixDphiErr.push_back(pix123DphiErr);
  pixpixDphiErr.push_back(pix124DphiErr);
  pixpixDphiErr.push_back(pix134DphiErr);
  pixpixDphiErr.push_back(pix234DphiErr);

  for(int nth = 0; nth < binSize; nth++){
  
     Long64_t nbytes = 0, nb = 0;
     vector<float> pix1egDphi_;
     vector<float> pix2egDphi_;
     vector<float> pix3egDphi_;
     vector<float> pix4egDphi_;

     vector<float> pix12egDphi_;
     vector<float> pix13egDphi_;
     vector<float> pix14egDphi_;
     vector<float> pix23egDphi_;
     vector<float> pix24egDphi_;
     vector<float> pix34egDphi_;

     vector<float> pix012Dphi_;
     vector<float> pix013Dphi_;
     vector<float> pix014Dphi_;
     vector<float> pix023Dphi_;
     vector<float> pix024Dphi_;
     vector<float> pix034Dphi_;
     vector<float> pix123Dphi_;
     vector<float> pix124Dphi_;
     vector<float> pix134Dphi_;
     vector<float> pix234Dphi_;

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
        if( eta_ == 6 && (fabs(ntEgEta->at(0)) < 2.7 || fabs(ntEgEta->at(0)) > 3.0)) continue; // test with barrel only

        if( ntEgEt->at(0) > (float) (10 + nth) && ntEgEt->at(0) < (float) (11 + nth) ){
          for(unsigned long i = 0; i < ntPix1EGdphi->size(); i++) {
             pix1egDphi_.push_back(ntPix1EGdphi->at(i)); 
             pix1egDphi_dist->Fill(ntEgEt->at(0), ntPix1EGdphi->at(i));
          }          
          for(unsigned long i = 0; i < ntPix2EGdphi->size(); i++) {
             pix2egDphi_.push_back(ntPix2EGdphi->at(i));          
             pix2egDphi_dist->Fill(ntEgEt->at(0), ntPix2EGdphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix3EGdphi->size(); i++) { 
             pix3egDphi_.push_back(ntPix3EGdphi->at(i));          
             pix3egDphi_dist->Fill(ntEgEt->at(0), ntPix3EGdphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix4EGdphi->size(); i++) { 
             pix4egDphi_.push_back(ntPix4EGdphi->at(i));          
             pix4egDphi_dist->Fill(ntEgEt->at(0), ntPix4EGdphi->at(i));
          }

          for(unsigned long i = 0; i < ntPix12EGdphi->size(); i++) { 
             if(fabs(ntPix12EGdphi->at(i)) > 0.2 ) continue;
             pix12egDphi_.push_back(ntPix12EGdphi->at(i));          
             pix12egDphi_dist->Fill(ntEgEt->at(0), ntPix12EGdphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix13EGdphi->size(); i++) {
             if(fabs(ntPix13EGdphi->at(i)) > 0.2 ) continue;
             pix13egDphi_.push_back(ntPix13EGdphi->at(i));          
             pix13egDphi_dist->Fill(ntEgEt->at(0), ntPix13EGdphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix14EGdphi->size(); i++) {
             if(fabs(ntPix14EGdphi->at(i)) > 0.2 ) continue;
             pix14egDphi_.push_back(ntPix14EGdphi->at(i));          
             pix14egDphi_dist->Fill(ntEgEt->at(0), ntPix14EGdphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix23EGdphi->size(); i++) { 
             if(fabs(ntPix23EGdphi->at(i)) > 0.2 ) continue;
             pix23egDphi_.push_back(ntPix23EGdphi->at(i));          
             pix23egDphi_dist->Fill(ntEgEt->at(0), ntPix23EGdphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix24EGdphi->size(); i++) { 
             if(fabs(ntPix24EGdphi->at(i)) > 0.2 ) continue;
             pix24egDphi_.push_back(ntPix24EGdphi->at(i));          
             pix24egDphi_dist->Fill(ntEgEt->at(0), ntPix24EGdphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix34EGdphi->size(); i++) {
             if(fabs(ntPix34EGdphi->at(i)) > 0.2 ) continue;
             pix34egDphi_.push_back(ntPix34EGdphi->at(i));          
             pix34egDphi_dist->Fill(ntEgEt->at(0), ntPix34EGdphi->at(i));
          }

          for(unsigned long i = 0; i < ntPix012dphi->size(); i++) {
             if(fabs(ntPix012dphi->at(i)) > 0.015 ) continue;
             pix012Dphi_.push_back(ntPix012dphi->at(i));          
             pix012Dphi_dist->Fill(ntEgEt->at(0), ntPix012dphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix013dphi->size(); i++) {
             if(fabs(ntPix013dphi->at(i)) > 0.015 ) continue;
             pix013Dphi_.push_back(ntPix013dphi->at(i));          
             pix013Dphi_dist->Fill(ntEgEt->at(0), ntPix013dphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix014dphi->size(); i++) { 
             if(fabs(ntPix014dphi->at(i)) > 0.015 ) continue;
             pix014Dphi_.push_back(ntPix014dphi->at(i));          
             pix014Dphi_dist->Fill(ntEgEt->at(0), ntPix014dphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix023dphi->size(); i++) {
             if(fabs(ntPix023dphi->at(i)) > 0.015 ) continue;
             pix023Dphi_.push_back(ntPix023dphi->at(i));          
             pix023Dphi_dist->Fill(ntEgEt->at(0), ntPix023dphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix024dphi->size(); i++) {
             if(fabs(ntPix024dphi->at(i)) > 0.015 ) continue;
             pix024Dphi_.push_back(ntPix024dphi->at(i));          
             pix024Dphi_dist->Fill(ntEgEt->at(0), ntPix024dphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix034dphi->size(); i++) {
             if(fabs(ntPix034dphi->at(i)) > 0.015 ) continue;
             pix034Dphi_.push_back(ntPix034dphi->at(i));          
             pix034Dphi_dist->Fill(ntEgEt->at(0), ntPix034dphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix123dphi->size(); i++) {
             if(fabs(ntPix123dphi->at(i)) > 0.015 ) continue;
             pix123Dphi_.push_back(ntPix123dphi->at(i));          
             pix123Dphi_dist->Fill(ntEgEt->at(0), ntPix123dphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix124dphi->size(); i++) {
             if(fabs(ntPix124dphi->at(i)) > 0.015 ) continue;
             pix124Dphi_.push_back(ntPix124dphi->at(i));          
             pix124Dphi_dist->Fill(ntEgEt->at(0), ntPix124dphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix134dphi->size(); i++) {
             if(fabs(ntPix134dphi->at(i)) > 0.015 ) continue;
             pix134Dphi_.push_back(ntPix134dphi->at(i));          
             pix134Dphi_dist->Fill(ntEgEt->at(0), ntPix134dphi->at(i));
          }
          for(unsigned long i = 0; i < ntPix234dphi->size(); i++) {
             if(fabs(ntPix234dphi->at(i)) > 0.015 ) continue;
             pix234Dphi_.push_back(ntPix234dphi->at(i));          
             pix234Dphi_dist->Fill(ntEgEt->at(0), ntPix234dphi->at(i));
          }
        }
     }// event loop

     std::sort (pix1egDphi_.begin(), pix1egDphi_.end()); 
     std::sort (pix2egDphi_.begin(), pix2egDphi_.end()); 
     std::sort (pix3egDphi_.begin(), pix3egDphi_.end()); 
     std::sort (pix4egDphi_.begin(), pix4egDphi_.end()); 

     std::sort (pix12egDphi_.begin(), pix12egDphi_.end()); 
     std::sort (pix13egDphi_.begin(), pix13egDphi_.end()); 
     std::sort (pix14egDphi_.begin(), pix14egDphi_.end()); 
     std::sort (pix23egDphi_.begin(), pix23egDphi_.end()); 
     std::sort (pix24egDphi_.begin(), pix24egDphi_.end()); 
     std::sort (pix34egDphi_.begin(), pix34egDphi_.end()); 

     std::sort (pix012Dphi_.begin(), pix012Dphi_.end()); 
     std::sort (pix013Dphi_.begin(), pix013Dphi_.end()); 
     std::sort (pix014Dphi_.begin(), pix014Dphi_.end()); 
     std::sort (pix023Dphi_.begin(), pix023Dphi_.end()); 
     std::sort (pix024Dphi_.begin(), pix024Dphi_.end()); 
     std::sort (pix034Dphi_.begin(), pix034Dphi_.end()); 
     std::sort (pix123Dphi_.begin(), pix123Dphi_.end()); 
     std::sort (pix124Dphi_.begin(), pix124Dphi_.end()); 
     std::sort (pix134Dphi_.begin(), pix134Dphi_.end()); 
     std::sort (pix234Dphi_.begin(), pix234Dphi_.end()); 

     int pix1egDphi_size = pix1egDphi_.size();
     int pix2egDphi_size = pix2egDphi_.size();
     int pix3egDphi_size = pix3egDphi_.size();
     int pix4egDphi_size = pix4egDphi_.size();

     int pix12egDphi_size = pix12egDphi_.size();
     int pix13egDphi_size = pix13egDphi_.size();
     int pix14egDphi_size = pix14egDphi_.size();
     int pix23egDphi_size = pix23egDphi_.size();
     int pix24egDphi_size = pix24egDphi_.size();
     int pix34egDphi_size = pix34egDphi_.size();

     int pix012Dphi_size = pix012Dphi_.size();
     int pix013Dphi_size = pix013Dphi_.size();
     int pix014Dphi_size = pix014Dphi_.size();
     int pix023Dphi_size = pix023Dphi_.size();
     int pix024Dphi_size = pix024Dphi_.size();
     int pix034Dphi_size = pix034Dphi_.size();
     int pix123Dphi_size = pix123Dphi_.size();
     int pix124Dphi_size = pix124Dphi_.size();
     int pix134Dphi_size = pix134Dphi_.size();
     int pix234Dphi_size = pix234Dphi_.size();

     // get index for median
     int pix1egDphi_medianIndex = getMedianIndex(pix1egDphi_);
     int pix2egDphi_medianIndex = getMedianIndex(pix2egDphi_);
     int pix3egDphi_medianIndex = getMedianIndex(pix3egDphi_);
     int pix4egDphi_medianIndex = getMedianIndex(pix4egDphi_);

     int pix12egDphi_medianIndex = getMedianIndex(pix12egDphi_);
     int pix13egDphi_medianIndex = getMedianIndex(pix13egDphi_);
     int pix14egDphi_medianIndex = getMedianIndex(pix14egDphi_);
     int pix23egDphi_medianIndex = getMedianIndex(pix23egDphi_);
     int pix24egDphi_medianIndex = getMedianIndex(pix24egDphi_);
     int pix34egDphi_medianIndex = getMedianIndex(pix34egDphi_);

     int pix012Dphi_medianIndex = getMedianIndex(pix012Dphi_);
     int pix013Dphi_medianIndex = getMedianIndex(pix013Dphi_);
     int pix014Dphi_medianIndex = getMedianIndex(pix014Dphi_);
     int pix023Dphi_medianIndex = getMedianIndex(pix023Dphi_);
     int pix024Dphi_medianIndex = getMedianIndex(pix024Dphi_);
     int pix034Dphi_medianIndex = getMedianIndex(pix034Dphi_);
     int pix123Dphi_medianIndex = getMedianIndex(pix123Dphi_);
     int pix124Dphi_medianIndex = getMedianIndex(pix124Dphi_);
     int pix134Dphi_medianIndex = getMedianIndex(pix134Dphi_);
     int pix234Dphi_medianIndex = getMedianIndex(pix234Dphi_);

     // get median
     float pix1egDphi_median = getMedian(pix1egDphi_);
     float pix2egDphi_median = getMedian(pix2egDphi_);
     float pix3egDphi_median = getMedian(pix3egDphi_);
     float pix4egDphi_median = getMedian(pix4egDphi_);

     float pix12egDphi_median = getMedian(pix12egDphi_);
     float pix13egDphi_median = getMedian(pix13egDphi_);
     float pix14egDphi_median = getMedian(pix14egDphi_);
     float pix23egDphi_median = getMedian(pix23egDphi_);
     float pix24egDphi_median = getMedian(pix24egDphi_);
     float pix34egDphi_median = getMedian(pix34egDphi_);

     float pix012Dphi_median = getMedian(pix012Dphi_);
     float pix013Dphi_median = getMedian(pix013Dphi_);
     float pix014Dphi_median = getMedian(pix014Dphi_);
     float pix023Dphi_median = getMedian(pix023Dphi_);
     float pix024Dphi_median = getMedian(pix024Dphi_);
     float pix034Dphi_median = getMedian(pix034Dphi_);
     float pix123Dphi_median = getMedian(pix123Dphi_);
     float pix124Dphi_median = getMedian(pix124Dphi_);
     float pix134Dphi_median = getMedian(pix134Dphi_);
     float pix234Dphi_median = getMedian(pix234Dphi_);

     // calculate error
     int pix1egDphi_low = (int)(pix1egDphi_medianIndex - (0.668 * pix1egDphi_medianIndex));
     int pix2egDphi_low = (int)(pix2egDphi_medianIndex - (0.668 * pix2egDphi_medianIndex));
     int pix3egDphi_low = (int)(pix3egDphi_medianIndex - (0.668 * pix3egDphi_medianIndex));
     int pix4egDphi_low = (int)(pix4egDphi_medianIndex - (0.668 * pix4egDphi_medianIndex));

     int pix1egDphi_high =(int)(pix1egDphi_medianIndex + (0.668 * (pix1egDphi_size-1-pix1egDphi_medianIndex)));
     int pix2egDphi_high =(int)(pix2egDphi_medianIndex + (0.668 * (pix2egDphi_size-1-pix2egDphi_medianIndex)));
     int pix3egDphi_high =(int)(pix3egDphi_medianIndex + (0.668 * (pix3egDphi_size-1-pix3egDphi_medianIndex)));
     int pix4egDphi_high =(int)(pix4egDphi_medianIndex + (0.668 * (pix4egDphi_size-1-pix4egDphi_medianIndex)));

     int pix12egDphi_low = (int)(pix12egDphi_medianIndex - (0.668 * pix12egDphi_medianIndex));
     int pix13egDphi_low = (int)(pix13egDphi_medianIndex - (0.668 * pix13egDphi_medianIndex));
     int pix14egDphi_low = (int)(pix14egDphi_medianIndex - (0.668 * pix14egDphi_medianIndex));
     int pix23egDphi_low = (int)(pix23egDphi_medianIndex - (0.668 * pix23egDphi_medianIndex));
     int pix24egDphi_low = (int)(pix24egDphi_medianIndex - (0.668 * pix24egDphi_medianIndex));
     int pix34egDphi_low = (int)(pix34egDphi_medianIndex - (0.668 * pix34egDphi_medianIndex));

     int pix12egDphi_high =(int)(pix12egDphi_medianIndex + (0.668 * (pix12egDphi_size-1-pix12egDphi_medianIndex)));
     int pix13egDphi_high =(int)(pix13egDphi_medianIndex + (0.668 * (pix13egDphi_size-1-pix13egDphi_medianIndex)));
     int pix14egDphi_high =(int)(pix14egDphi_medianIndex + (0.668 * (pix14egDphi_size-1-pix14egDphi_medianIndex)));
     int pix23egDphi_high =(int)(pix23egDphi_medianIndex + (0.668 * (pix23egDphi_size-1-pix23egDphi_medianIndex)));
     int pix24egDphi_high =(int)(pix24egDphi_medianIndex + (0.668 * (pix24egDphi_size-1-pix24egDphi_medianIndex)));
     int pix34egDphi_high =(int)(pix34egDphi_medianIndex + (0.668 * (pix34egDphi_size-1-pix34egDphi_medianIndex)));

     int pix012Dphi_low = (int)(pix012Dphi_medianIndex - (0.668 * pix012Dphi_medianIndex));
     int pix013Dphi_low = (int)(pix013Dphi_medianIndex - (0.668 * pix013Dphi_medianIndex));
     int pix014Dphi_low = (int)(pix014Dphi_medianIndex - (0.668 * pix014Dphi_medianIndex));
     int pix023Dphi_low = (int)(pix023Dphi_medianIndex - (0.668 * pix023Dphi_medianIndex));
     int pix024Dphi_low = (int)(pix024Dphi_medianIndex - (0.668 * pix024Dphi_medianIndex));
     int pix034Dphi_low = (int)(pix034Dphi_medianIndex - (0.668 * pix034Dphi_medianIndex));
     int pix123Dphi_low = (int)(pix123Dphi_medianIndex - (0.668 * pix123Dphi_medianIndex));
     int pix124Dphi_low = (int)(pix124Dphi_medianIndex - (0.668 * pix124Dphi_medianIndex));
     int pix134Dphi_low = (int)(pix134Dphi_medianIndex - (0.668 * pix134Dphi_medianIndex));
     int pix234Dphi_low = (int)(pix234Dphi_medianIndex - (0.668 * pix234Dphi_medianIndex));

     int pix012Dphi_high = (int)(pix012Dphi_medianIndex + (0.668 * (pix012Dphi_size-1-pix012Dphi_medianIndex)));
     int pix013Dphi_high = (int)(pix013Dphi_medianIndex + (0.668 * (pix013Dphi_size-1-pix013Dphi_medianIndex)));
     int pix014Dphi_high = (int)(pix014Dphi_medianIndex + (0.668 * (pix014Dphi_size-1-pix014Dphi_medianIndex)));
     int pix023Dphi_high = (int)(pix023Dphi_medianIndex + (0.668 * (pix023Dphi_size-1-pix023Dphi_medianIndex)));
     int pix024Dphi_high = (int)(pix024Dphi_medianIndex + (0.668 * (pix024Dphi_size-1-pix024Dphi_medianIndex)));
     int pix034Dphi_high = (int)(pix034Dphi_medianIndex + (0.668 * (pix034Dphi_size-1-pix034Dphi_medianIndex)));
     int pix123Dphi_high = (int)(pix123Dphi_medianIndex + (0.668 * (pix123Dphi_size-1-pix123Dphi_medianIndex)));
     int pix124Dphi_high = (int)(pix124Dphi_medianIndex + (0.668 * (pix124Dphi_size-1-pix124Dphi_medianIndex)));
     int pix134Dphi_high = (int)(pix134Dphi_medianIndex + (0.668 * (pix134Dphi_size-1-pix134Dphi_medianIndex)));
     int pix234Dphi_high = (int)(pix234Dphi_medianIndex + (0.668 * (pix234Dphi_size-1-pix234Dphi_medianIndex)));

     float pix1egDphi_medianErr = 0.;
     float pix2egDphi_medianErr = 0.;
     float pix3egDphi_medianErr = 0.;
     float pix4egDphi_medianErr = 0.;

     float pix12egDphi_medianErr = 0.;
     float pix13egDphi_medianErr = 0.;
     float pix14egDphi_medianErr = 0.;
     float pix23egDphi_medianErr = 0.;
     float pix24egDphi_medianErr = 0.;
     float pix34egDphi_medianErr = 0.;

     float pix012Dphi_medianErr =0.;
     float pix013Dphi_medianErr =0.; 
     float pix014Dphi_medianErr =0.; 
     float pix023Dphi_medianErr =0.; 
     float pix024Dphi_medianErr =0.; 
     float pix034Dphi_medianErr =0.; 
     float pix123Dphi_medianErr =0.; 
     float pix124Dphi_medianErr =0.; 
     float pix134Dphi_medianErr =0.; 
     float pix234Dphi_medianErr =0.; 

     x[nth] = 10.5 + float(nth);
     xErr[nth] = 0.;

     pix1egDphi[nth] = pix1egDphi_median;
     pix2egDphi[nth] = pix2egDphi_median;
     pix3egDphi[nth] = pix3egDphi_median;
     pix4egDphi[nth] = pix4egDphi_median;

     pix1egDphiErr[nth] = pix1egDphi_medianErr/2.;
     pix2egDphiErr[nth] = pix2egDphi_medianErr/2.;
     pix3egDphiErr[nth] = pix3egDphi_medianErr/2.;
     pix4egDphiErr[nth] = pix4egDphi_medianErr/2.;

     pix12egDphi[nth] = pix12egDphi_median;
     pix13egDphi[nth] = pix13egDphi_median;
     pix14egDphi[nth] = pix14egDphi_median;
     pix23egDphi[nth] = pix23egDphi_median;
     pix24egDphi[nth] = pix24egDphi_median;
     pix34egDphi[nth] = pix34egDphi_median;

     pix12egDphiErr[nth] = pix12egDphi_medianErr/2.;
     pix13egDphiErr[nth] = pix13egDphi_medianErr/2.;
     pix14egDphiErr[nth] = pix14egDphi_medianErr/2.;
     pix23egDphiErr[nth] = pix23egDphi_medianErr/2.;
     pix24egDphiErr[nth] = pix24egDphi_medianErr/2.;
     pix34egDphiErr[nth] = pix34egDphi_medianErr/2.;

     pix012Dphi[nth] = pix012Dphi_median;
     pix013Dphi[nth] = pix013Dphi_median;
     pix014Dphi[nth] = pix014Dphi_median;
     pix023Dphi[nth] = pix023Dphi_median;
     pix024Dphi[nth] = pix024Dphi_median;
     pix034Dphi[nth] = pix034Dphi_median;
     pix123Dphi[nth] = pix123Dphi_median;
     pix124Dphi[nth] = pix124Dphi_median;
     pix134Dphi[nth] = pix134Dphi_median;
     pix234Dphi[nth] = pix234Dphi_median;

     pix012DphiErr[nth] = pix012Dphi_medianErr/2.;
     pix013DphiErr[nth] = pix013Dphi_medianErr/2.;
     pix014DphiErr[nth] = pix014Dphi_medianErr/2.;
     pix023DphiErr[nth] = pix023Dphi_medianErr/2.;
     pix024DphiErr[nth] = pix024Dphi_medianErr/2.;
     pix034DphiErr[nth] = pix034Dphi_medianErr/2.;
     pix123DphiErr[nth] = pix123Dphi_medianErr/2.;
     pix124DphiErr[nth] = pix124Dphi_medianErr/2.;
     pix134DphiErr[nth] = pix134Dphi_medianErr/2.;
     pix234DphiErr[nth] = pix234Dphi_medianErr/2.;

  }// Et scanning loop


  // save the fit parameters
  ofstream fit_result1;
 
  TString prefix1 = "ROI_"; 
  TString postfix1 = ".txt"; 
 
  //setTDRStyle();
  fit_result1.open(prefix1+eta_region[eta_-1]+postfix1);

  // draw eg-pix dphi 
  for(int i = 0; i < 4; i++){

     writeExtraText = true;
     extraText  = "Phase-2 simulation";
     gROOT->SetBatch();
     TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
     gStyle->SetOptStat(0);
     gStyle->SetPalette(1);
     c1->SetTopMargin(0.07);
     c1->SetRightMargin(0.05);
     c1->SetGridy();
     c1->SetGridx();
     TGaxis::SetMaxDigits(3);
     c1->cd();

     //

     TH2F* scatter_plot = NULL;
     if(i==0) scatter_plot = pix1egDphi_dist;
     if(i==1) scatter_plot = pix2egDphi_dist;
     if(i==2) scatter_plot = pix3egDphi_dist;
     if(i==3) scatter_plot = pix4egDphi_dist;

     if(scatter_plot != NULL){
       scatter_plot->GetYaxis()->SetRangeUser(-0.2,0.2);
       scatter_plot->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
       scatter_plot->GetYaxis()->SetTitle("#Delta#phi [rad.]");
       scatter_plot->SetMarkerColor(kRed);
       //scatter_plot->SetMarkerSize(.25);
       scatter_plot->SetMarkerStyle(6);
       scatter_plot->Draw("scat=1.");
     }
     if(scatter_plot == NULL) std::cout << "NULL" << std::endl;

     TGraphErrors* gr1 = new TGraphErrors(binSize,x,pixegDphi.at(i),xErr,pixegDphiErr.at(i));

     TString nthsw;
     nthsw.Form("%d", i+1);
     gr1->SetMarkerColor(kBlue);
     gr1->SetLineColor(kBlue);
     gr1->SetMarkerSize(.8);
     gr1->SetMarkerStyle(20);
     gr1->SetFillColorAlpha(kBlue, 0.4); // for error
     gr1->Draw("samepE3");

     TF1 *median_fitFunc = new TF1("func","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])) )", 10., 100.);
     median_fitFunc->SetLineColor(kGreen);
     median_fitFunc->SetLineStyle(1);
     median_fitFunc->SetLineWidth(2);

     gr1->Fit(median_fitFunc,"0");
     median_fitFunc->Draw("lsame");

     TLegend* leg = new TLegend(0.3,0.7,0.5,0.9);
     leg->AddEntry(scatter_plot, "Scattered #Delta#phi distribution as a function of L1 EG E_{T}","pl");
     leg->AddEntry(gr1, "Median point of #Delta#phi","fp");
     leg->AddEntry(median_fitFunc,"Fit of median point", "l");
     leg->AddEntry((TObject*)0, eta_region[eta_-1], "");
     leg->SetTextFont(42);
     leg->SetTextSize(0.035);
     leg->SetBorderSize(0);
     leg->SetFillStyle(0);
     leg->Draw();

     CMS_lumi( c1, 4, 0 );
     c1->Update();

     c1->SaveAs("PixEG"+nthsw+"_"+eta_region[eta_-1]+".png");

     fit_result1 << endl;
     fit_result1 << "if( region == " << eta_ << " && i == " << i << " ){" <<endl;
     for( int j=0; j < 4; j++){
         fit_result1 << "p[" << j << "] = " << median_fitFunc->GetParameter(j) << ";" << endl;
     }
     fit_result1 << "}" << endl;
     fit_result1 << endl;

     delete gr1;
     delete c1;
  }

  // save the fit parameters
  ofstream fit_result2;

  TString prefix2 = "EGmatching_";
  TString postfix2 = ".txt";

  fit_result2.open(prefix2+eta_region[eta_-1]+postfix2);

  const TString eg_pixpixSW[6] = {"12","13","14","23","24","34"};  
  // draw eg-pixpix dphi 
  for(int i = 0; i < 6; i++){

     writeExtraText = true;
     extraText  = "Phase-2 simulation";
     gROOT->SetBatch();
     TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
     gStyle->SetOptStat(0);
     gStyle->SetPalette(1);
     c1->SetTopMargin(0.07);
     c1->SetRightMargin(0.05);
     c1->SetGridy();
     c1->SetGridx();
     TGaxis::SetMaxDigits(3);
     c1->cd();

     //
     TGraphErrors* gr1 = new TGraphErrors(binSize,x,pixpixegDphi.at(i),xErr,pixpixegDphiErr.at(i));

     TH2F* scatter_plot = NULL;
     if(i==0) scatter_plot = pix12egDphi_dist;
     if(i==1) scatter_plot = pix13egDphi_dist;
     if(i==2) scatter_plot = pix14egDphi_dist;
     if(i==3) scatter_plot = pix23egDphi_dist;
     if(i==4) scatter_plot = pix24egDphi_dist;
     if(i==5) scatter_plot = pix34egDphi_dist;

     if(scatter_plot != NULL){
       scatter_plot->GetYaxis()->SetRangeUser(-0.2,0.2);
       scatter_plot->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
       scatter_plot->GetYaxis()->SetTitle("#Delta#phi [rad.]");
       scatter_plot->SetMarkerColor(kRed);
       //scatter_plot->SetMarkerSize(.25);
       scatter_plot->SetMarkerStyle(6);
       scatter_plot->Draw("scat=1.");

     }

     TString nthsw;
     nthsw.Form("%d", i+1);
     gr1->SetLineColor(kBlue);
     gr1->SetMarkerColor(kBlue);
     gr1->SetMarkerSize(.8);
     gr1->SetMarkerStyle(20);
     gr1->SetFillColorAlpha(kBlue, 0.4);
     gr1->SetLineWidth(1);
     gr1->Draw("samep");

     TF1 *median_fitFunc = new TF1("func","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])) )", 10., 100.);
     median_fitFunc->SetLineColor(kGreen);
     median_fitFunc->SetLineStyle(1);
     median_fitFunc->SetLineWidth(2);

     gr1->Fit(median_fitFunc,"0");
     median_fitFunc->Draw("lsame");

     double ab[4]={0};
     median_fitFunc->GetParameters(ab);
//MARK
double x1[200],y1[200],x2[200],y2[200];
for ( int j = 10;j<120;j++){
    x1[j]=j;
    x2[j]=0.5;
    y1[j]=ab[0]*pow(j,0) + ab[1]*pow(j,ab[2])*exp(-pow(j,ab[3]));
    y2[j]=0.03;
}
    TGraphErrors* gr2 = new TGraphErrors(110,x1,y1,x2,y2); //0.0017 0.003
     
     gr2->SetLineColor(kBlue);
     gr2->SetFillColorAlpha(kBlue, 0.40);
//    
     gr2->Draw("sameE3");
     median_fitFunc->Draw("lsame");
     TLegend* leg = new TLegend(0.3,0.2,0.5,0.4);
     leg->AddEntry(scatter_plot, "Scattered #Delta#phi distribution as a function of L1 EG E_{T}","pl");
     leg->AddEntry(gr1, "Median point of #Delta#phi","lp");
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

 //    c1->SaveAs("PixPixEG"+nthsw+"_"+eta_region[eta_-1]+".pdf");
     c1->SaveAs("PixPixEG"+nthsw+"_"+eta_region[eta_-1]+".png");

     fit_result2 << endl;
     fit_result2 << "if( region == " << eta_ << " && i == " << i << " ){" <<endl;
     for( int j=0; j < 4; j++){
         fit_result2 << "p[" << j << "] = " << median_fitFunc->GetParameter(j) << ";" << endl;
     }
     fit_result2 << "}" << endl;
     fit_result2 << endl;
delete gr2;
     delete gr1;
     delete c1;
  }

  // save the fit parameters
  ofstream fit_result3;

  TString prefix3 = "Pixelmatching_";
  TString postfix3 = ".txt";

  fit_result3.open(prefix3+eta_region[eta_-1]+postfix3);

  const TString pixpix1_SW[10] = {"01","01","01","02","02","03","12","12","13","23"};  
  const TString pixpix2_SW[10] = {"12","13","14","23","24","34","23","24","34","34"};  

  // draw pix-pixpix dphi 
  for(int i = 0; i < 10; i++){

     writeExtraText = true;
     extraText  = "Phase-2 simulation";
     gROOT->SetBatch();
     TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
     gStyle->SetOptStat(0);
     gStyle->SetPalette(1);
     c1->SetTopMargin(0.07);
     c1->SetRightMargin(0.05);
     c1->SetGridy();
     c1->SetGridx();
     TGaxis::SetMaxDigits(4);

     //  
     TGraphErrors* gr1 = new TGraphErrors(binSize,x,pixpixDphi.at(i),xErr,pixpixDphiErr.at(i));

     TH2F* scatter_plot = NULL;
     if(i==0) scatter_plot = pix012Dphi_dist;
     if(i==1) scatter_plot = pix013Dphi_dist;
     if(i==2) scatter_plot = pix014Dphi_dist;
     if(i==3) scatter_plot = pix023Dphi_dist;
     if(i==4) scatter_plot = pix024Dphi_dist;
     if(i==5) scatter_plot = pix034Dphi_dist;
     if(i==6) scatter_plot = pix123Dphi_dist;
     if(i==7) scatter_plot = pix124Dphi_dist;
     if(i==8) scatter_plot = pix134Dphi_dist;
     if(i==9) scatter_plot = pix234Dphi_dist;

     if(scatter_plot != NULL){

       scatter_plot->GetYaxis()->SetRangeUser(-0.015,0.015);
       scatter_plot->GetXaxis()->SetTitle("L1 E/gamma E_{T} [GeV]");
       scatter_plot->GetYaxis()->SetTitle("#Delta#phi [rad.]");
       scatter_plot->SetMarkerColor(kRed);
       //scatter_plot->SetMarkerSize(.25);
       scatter_plot->SetMarkerStyle(6);
       scatter_plot->Draw("scat=1.");
     }

     TString nthsw;
     nthsw.Form("%d", i+1);
     gr1->SetLineColor(kBlue);
     gr1->SetMarkerColor(kBlue);
     gr1->SetMarkerSize(.8);
     gr1->SetMarkerStyle(20);
     gr1->SetFillColorAlpha(kBlue, 0.4);
     gr1->SetLineWidth(1);
     gr1->Draw("samep");

     TF1 *median_fitFunc = new TF1("func","( [0]*pow(x,0) + [1]*pow(x,[2])*exp(-pow(x,[3])) )", 10., 100.);

     median_fitFunc->SetParLimits(0, -0.05, 0.001);
     median_fitFunc->SetParLimits(1, -0.001, 5);
     median_fitFunc->SetParLimits(2, -2., 0.2);
     median_fitFunc->SetParLimits(3, -0.7, 0.5);

     median_fitFunc->SetLineColor(kGreen);
     median_fitFunc->SetLineStyle(1);
     median_fitFunc->SetLineWidth(2);

     gr1->Fit(median_fitFunc,"0");

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
     
     gr2->SetLineColor(kBlue);
     gr2->SetFillColorAlpha(kBlue, 0.40);
//    
     gr2->Draw("sameE3");
     median_fitFunc->Draw("lsame");

     TLegend* leg = new TLegend(0.2,0.2,0.4,0.4);
   //  leg->AddEntry(scatter_plot, "Scattered #Delta#phi distribution as a function of L1 EG E_{T}","peZ");
     leg->AddEntry(scatter_plot, "Scattered #Delta#phi distribution as a function of L1 EG E_{T}","pl");
     leg->AddEntry(gr1, "Median point of #Delta#phi", "lp");
  //   leg->AddEntry(gr1, "Median point of #Delta#phi", "pe3");
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

     c1->SaveAs("PixPix"+nthsw+"_"+eta_region[eta_-1]+".png");

     fit_result3 << endl;
     fit_result3 << "if( region == " << eta_ << " && i == " << i << " ){" <<endl;
     for( int j=0; j < 4; j++){
         fit_result3 << "p[" << j << "] = " << median_fitFunc->GetParameter(j) << ";" << endl;
     }
     fit_result3 << "}" << endl;
     fit_result3 << endl;
delete gr2;
     delete gr1;
     delete c1; 
  }

 delete pix1egDphi_dist;
 delete pix2egDphi_dist;
 delete pix3egDphi_dist;
 delete pix4egDphi_dist;

 delete pix12egDphi_dist;
 delete pix13egDphi_dist;
 delete pix14egDphi_dist;
 delete pix23egDphi_dist;
 delete pix24egDphi_dist;
 delete pix34egDphi_dist;

 delete pix012Dphi_dist;
 delete pix013Dphi_dist;
 delete pix014Dphi_dist;
 delete pix023Dphi_dist;
 delete pix024Dphi_dist;
 delete pix034Dphi_dist;
 delete pix123Dphi_dist;
 delete pix124Dphi_dist;
 delete pix134Dphi_dist;
 delete pix234Dphi_dist;
 std::cout << "DONE" << std::endl;
}
