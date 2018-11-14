#define basic_study_cxx
#include "basic_study.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TLorentzVector.h>
void basic_study::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TFile *outfile;
   outfile = new TFile("./Eff_nopu.root","recreate");

   TH1F* dphi_genEle = new TH1F("dphi_genEle",";; ", 100, 0., 3.15);
   TH1F* deta_genEle = new TH1F("deta_genEle",";; ", 1000, 0., 100);
   dphi_genEle->Sumw2();
   deta_genEle->Sumw2();

   TH1F* genEle_eta = new TH1F("genEle_eta",";; ", 80, -4., 4.);
   TH1F* l1eg_eta = new TH1F("l1eg_eta",";; ", 80, -4., 4.);
   TH1F* l1eg_eta_be = new TH1F("l1eg_eta_be",";; ", 80, -4., 4.);
   TH1F* l1eg_eta_fe = new TH1F("l1eg_eta_fe",";; ", 80, -4., 4.);
   genEle_eta->Sumw2();
   l1eg_eta->Sumw2();
   l1eg_eta_be->Sumw2(); 
   l1eg_eta_fe->Sumw2();

   TH1F* bpix1_eta = new TH1F("bpix1_eta",";; ", 80, -4., 4.);
   TH1F* bpix2_eta = new TH1F("bpix2_eta",";; ", 80, -4., 4.);
   TH1F* bpix3_eta = new TH1F("bpix3_eta",";; ", 80, -4., 4.);
   TH1F* bpix4_eta = new TH1F("bpix4_eta",";; ", 80, -4., 4.);
   bpix1_eta->Sumw2(); 
   bpix2_eta->Sumw2();
   bpix3_eta->Sumw2();
   bpix4_eta->Sumw2();

   TH1F* fpix1_eta = new TH1F("fpix1_eta",";; ", 80, -4., 4.);
   TH1F* fpix2_eta = new TH1F("fpix2_eta",";; ", 80, -4., 4.);
   TH1F* fpix3_eta = new TH1F("fpix3_eta",";; ", 80, -4., 4.);
   TH1F* fpix4_eta = new TH1F("fpix4_eta",";; ", 80, -4., 4.);
   TH1F* fpix5_eta = new TH1F("fpix5_eta",";; ", 80, -4., 4.);
   fpix1_eta->Sumw2(); 
   fpix2_eta->Sumw2();
   fpix3_eta->Sumw2();
   fpix4_eta->Sumw2();
   fpix5_eta->Sumw2();

   TH1F* bpix1234_eta = new TH1F("bpix1234_eta",";; ", 80, -4., 4.);
   TH1F* bpix123d1_eta = new TH1F("bpix123d1_eta",";; ", 80, -4., 4.);
   TH1F* bpix12d12_eta = new TH1F("bpix12d12_eta",";; ", 80, -4., 4.);
   TH1F* bpix1d123_eta = new TH1F("bpix1d123_eta",";; ", 80, -4., 4.);
   TH1F* fpix1234_eta = new TH1F("fpix1234_eta",";; ", 80, -4., 4.);
   TH1F* fpix2345_eta = new TH1F("fpix2345_eta",";; ", 80, -4., 4.);
   bpix1234_eta->Sumw2();
   bpix123d1_eta->Sumw2();
   bpix12d12_eta->Sumw2();
   bpix1d123_eta->Sumw2();
   fpix1234_eta->Sumw2();
   fpix2345_eta->Sumw2();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      dphi_genEle->Fill(fabs(genPartPhi->at(0)-genPartPhi->at(1)),1.);
      deta_genEle->Fill(fabs(genPartEta->at(0)-genPartPhi->at(1)),1.);
      
      //find closest egamma object to the gen electron
      int EgN=egCrysClusterEt->size();
      int cl3d_N_ = cl3d_pt->size();

      float closest_dr = 9999.;
      int closest_eg = 0;
      int egCount = 0;

      // loop over barrel egamma objects
      for(int i=0; i < EgN;i++){

         float dPhi = deltaPhi(propgenElPartPhi->at(0), egCrysClusterPhi->at(i));

         float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-egCrysClusterEta->at(i),2));
         if(egCrysClusterEt->at(i) < 10) continue;
         if(current_dr < closest_dr){
           closest_dr = current_dr;
           closest_eg = i;
         }
      }// end of loop to find the closest egamma to gen electron 

      // HGCAL 3D cluster
      float closest_cl3d_dr = 9999.;
      int closest_cl3d = 0;
      int cl3d_Count = 0;

      for(int i=0; i < cl3d_N_;i++){

         if(cl3d_egid->at(i) != 1) continue; // egamma id

         float dPhi = deltaPhi(propgenElPartPhi->at(0), cl3d_phi->at(i));

         float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-cl3d_eta->at(i),2));
         if(cl3d_pt->at(i) < 10) continue;
         cl3d_Count++;
         if(current_dr < closest_cl3d_dr){
           closest_cl3d_dr = current_dr;
           closest_cl3d = i;
         }
      }// end of loop to find the closest egamma to gen electron 

      // check l1eg efficiency && pixel rechit efficiency
      float dr_cut = 0.1;
      if(propgenElPartPt->at(0) > 20){
        genEle_eta->Fill(propgenElPartEta->at(0), 1.);

        float EgPhi = 999.;
        if((closest_cl3d_dr < dr_cut)|| (closest_dr < dr_cut)){
          if( closest_dr < closest_cl3d_dr ){
            l1eg_eta->Fill(propgenElPartEta->at(0), 1.);
            EgPhi = egCrysClusterPhi->at(closest_eg);
          }
          else{
            l1eg_eta->Fill(propgenElPartEta->at(0), 1.);
            EgPhi = cl3d_phi->at(closest_cl3d);
          }
        }

        if((closest_cl3d_dr < dr_cut)){
          l1eg_eta_be->Fill(propgenElPartEta->at(0), 1.);
        }
        if((closest_dr < dr_cut)){
          l1eg_eta_fe->Fill(propgenElPartEta->at(0), 1.);
        }

        int  bpix1_ = 0;
        int  bpix2_ = 0;
        int  bpix3_ = 0;
        int  bpix4_ = 0;

        int  fpix1_ = 0;
        int  fpix2_ = 0;
        int  fpix3_ = 0;
        int  fpix4_ = 0;
        int  fpix5_ = 0;

        if( EgPhi != 999.){
          for(int a=0; a<bRecHitN; a++){
             TVector3 current_hit;
             current_hit.SetXYZ( bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a) );
             float Dphi = deltaPhi( current_hit.Phi(), EgPhi);
             if(fabs(Dphi) > 0.1) continue;
             if(bRecHitLayer->at(a)==1) bpix1_=1; 
             if(bRecHitLayer->at(a)==2) bpix2_=1; 
             if(bRecHitLayer->at(a)==3) bpix3_=1; 
             if(bRecHitLayer->at(a)==4) bpix4_=1; 
          }

          for(int a=0; a<fRecHitN; a++){
             TVector3 current_hit;
             current_hit.SetXYZ( fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) );
             float Dphi = deltaPhi( current_hit.Phi(), EgPhi);
             if(fabs(Dphi) > 0.1) continue;
             if(fRecHitDisk->at(a)==1) fpix1_=1;
             if(fRecHitDisk->at(a)==2) fpix2_=1;
             if(fRecHitDisk->at(a)==3) fpix3_=1;
             if(fRecHitDisk->at(a)==4) fpix4_=1;
             if(fRecHitDisk->at(a)==5) fpix5_=1;
          }
        }

        if(bpix1_) bpix1_eta->Fill(propgenElPartEta->at(0), 1.);
        if(bpix2_) bpix2_eta->Fill(propgenElPartEta->at(0), 1.);
        if(bpix3_) bpix3_eta->Fill(propgenElPartEta->at(0), 1.);
        if(bpix4_) bpix4_eta->Fill(propgenElPartEta->at(0), 1.);

        if(fpix1_) fpix1_eta->Fill(propgenElPartEta->at(0), 1.);
        if(fpix2_) fpix2_eta->Fill(propgenElPartEta->at(0), 1.);
        if(fpix3_) fpix3_eta->Fill(propgenElPartEta->at(0), 1.);
        if(fpix4_) fpix4_eta->Fill(propgenElPartEta->at(0), 1.);
        if(fpix5_) fpix5_eta->Fill(propgenElPartEta->at(0), 1.);

        if(bpix1_+bpix2_+bpix3_+bpix4_ >=3) bpix1234_eta->Fill(propgenElPartEta->at(0), 1.);
        if(bpix1_+bpix2_+bpix3_+fpix1_ >=3) bpix123d1_eta->Fill(propgenElPartEta->at(0), 1.);
        if(bpix1_+bpix2_+fpix1_+fpix2_ >=3) bpix12d12_eta->Fill(propgenElPartEta->at(0), 1.);
        if(bpix1_+fpix1_+fpix2_+fpix3_ >=3) bpix1d123_eta->Fill(propgenElPartEta->at(0), 1.);
        if(fpix1_+fpix2_+fpix3_+fpix4_ >=3) fpix1234_eta->Fill(propgenElPartEta->at(0), 1.);
        if(fpix2_+fpix3_+fpix4_+fpix5_ >=3) fpix2345_eta->Fill(propgenElPartEta->at(0), 1.);


      }// gen pt cut
   }// event loop
   outfile->Write();
}
