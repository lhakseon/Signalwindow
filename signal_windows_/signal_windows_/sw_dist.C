#define sw_dist_cxx
#include "sw_dist.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void sw_dist::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   float dr_cut = 0.1;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      ntEgEt.clear();
      ntEgEta.clear();
      ntEgPhi.clear();
      ntptErr.clear();

      ntPix1EGdphi.clear();
      ntPix2EGdphi.clear();
      ntPix3EGdphi.clear();
      ntPix4EGdphi.clear();

      ntPix12EGdphi.clear();
      ntPix13EGdphi.clear();
      ntPix14EGdphi.clear();
      ntPix23EGdphi.clear();
      ntPix24EGdphi.clear();
      ntPix34EGdphi.clear();

      ntPix12EGdeta.clear();
      ntPix13EGdeta.clear();
      ntPix14EGdeta.clear();
      ntPix23EGdeta.clear();
      ntPix24EGdeta.clear();
      ntPix34EGdeta.clear();

      ntPix012dphi.clear();
      ntPix013dphi.clear();
      ntPix014dphi.clear();
      ntPix023dphi.clear();
      ntPix024dphi.clear();
      ntPix034dphi.clear();
      ntPix123dphi.clear();
      ntPix124dphi.clear();
      ntPix134dphi.clear();
      ntPix234dphi.clear();
      ntPix123deta.clear();
      ntPix124deta.clear();
      ntPix134deta.clear();
      ntPix234deta.clear();

      float ptErr = 0.;
      float closest_dr = 9999.;
      int closest_eg = 0;
      int EgN = egCrysClusterEta->size();

      // loop over barrel egamma objects
      for(int i=0; i < EgN;i++){
         float dPhi = deltaPhi(propgenElPartPhi->at(0), egCrysClusterPhi->at(i));

         float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-egCrysClusterEta->at(i),2));
         if(egCrysClusterEt->at(i) < 5) continue;
         if(current_dr < closest_dr){
           closest_dr = current_dr;
           closest_eg = i;
         }
      }// end of loop to find the closest egamma to gen electron 

      // HGCAL 3D cluster
      float closest_cl3d_dr = 9999.;
      int closest_cl3d = 0;
      int cl3d_N_ = cl3d_eta->size();

      for(int i=0; i < cl3d_N_;i++){

         if(cl3d_egid->at(i) != 1) continue; // egamma id

         float dPhi = deltaPhi(propgenElPartPhi->at(0), cl3d_phi->at(i));

         float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-cl3d_eta->at(i),2));
         if(cl3d_pt->at(i) < 5) continue;
         if(current_dr < closest_cl3d_dr){
           closest_cl3d_dr = current_dr;
           closest_cl3d = i;
         }
      }// end of loop to find the closest egamma to gen electron 

     if((closest_cl3d_dr < dr_cut && closest_cl3d_dr != 9999.)|| (closest_dr < dr_cut && closest_dr != 9999.)){

      if( closest_dr < closest_cl3d_dr ){
        EgEt =egCrysClusterEt ->at(closest_eg);
        EgEta=egCrysClusterEta->at(closest_eg);
        EgPhi=egCrysClusterPhi->at(closest_eg);
        ptErr = (egCrysClusterEt ->at(closest_eg)-propgenElPartPt->at(0))/propgenElPartPt->at(0);

        float EgGx = egCrysClusterGx->at(closest_eg);
        float EgGy = egCrysClusterGy->at(closest_eg);
        float EgGz = egCrysClusterGz->at(closest_eg);
        emvector.SetXYZ(EgGx,EgGy,EgGz);
      }
      else{

          EgEt =cl3d_pt->at(closest_cl3d);
          EgEta=cl3d_eta->at(closest_cl3d);
          EgPhi=cl3d_phi->at(closest_cl3d);
          ptErr = (cl3d_pt->at(closest_eg)-propgenElPartPt->at(0))/propgenElPartPt->at(0);

          float EgGx = cl3d_x->at(closest_cl3d);
          float EgGy = cl3d_y->at(closest_cl3d);
          float EgGz = (float)cl3d_z->at(closest_cl3d);
          emvector.SetXYZ(EgGx,EgGy,EgGz);
      }
     }
      int eta_region = 0;
      if( fabs(EgEta) < 0.8 ) eta_region =1;
      if( fabs(EgEta) < 1.4 && fabs(EgEta) > 0.8 ) eta_region =2;
      if( fabs(EgEta) < 1.7 && fabs(EgEta) > 1.4 ) eta_region =3;
      if( fabs(EgEta) < 2.1 && fabs(EgEta) > 1.7 ) eta_region =4;
      if( fabs(EgEta) < 2.7 && fabs(EgEta) > 2.1 ) eta_region =5;
      if( fabs(EgEta) < 3.0 && fabs(EgEta) > 2.7 ) eta_region =6;
      if( fabs(EgEta) > 3. ) continue;

      layers[0] = 1; // beam spot
      layers[1] = 0; layers[2] = 0; layers[3] = 0; layers[4] = 0;

      first_layer_hits.clear();
      second_layer_hits.clear();
      third_layer_hits.clear();
      fourth_layer_hits.clear();
      // save 4 pixel rechits each from 4 different pixel layer/ disk according to eta of the egamma
      StorePixelHit(eta_region);

      ntEgEt.push_back(EgEt);
      ntEgEta.push_back(EgEta);
      ntEgPhi.push_back(EgPhi);
      ntptErr.push_back(ptErr);

      // pixel EG dphi
      for(unsigned long i = 0; i < first_layer_hits.size(); i++){
        ntPix1EGdphi.push_back(deltaPhi(first_layer_hits.at(i).Phi(), emvector.Phi()));
      }
      for(unsigned long i = 0; i < second_layer_hits.size(); i++){
        ntPix2EGdphi.push_back(deltaPhi(second_layer_hits.at(i).Phi(), emvector.Phi()));
      }
      for(unsigned long i = 0; i < third_layer_hits.size(); i ++){
        ntPix3EGdphi.push_back(deltaPhi(third_layer_hits.at(i).Phi(), emvector.Phi()));
      }
      for(unsigned long i = 0; i < fourth_layer_hits.size(); i++){
        ntPix4EGdphi.push_back(deltaPhi(fourth_layer_hits.at(i).Phi(), emvector.Phi()));
      }

      // pixel segment-EG dphi 
      // pixel1 pixel2
      for(unsigned long i = 0; i < first_layer_hits.size(); i++){
         for(unsigned long j = 0; j < second_layer_hits.size(); j++){
           TVector3 pixelVector = second_layer_hits.at(j) - first_layer_hits.at(i);
           TVector3 EM_pixelVector = emvector - second_layer_hits.at(j);
           ntPix12EGdphi.push_back(deltaPhi( EM_pixelVector.Phi(), pixelVector.Phi() ));
           ntPix12EGdeta.push_back(EM_pixelVector.Eta()-pixelVector.Eta());

           ntPix012dphi.push_back(deltaPhi(pixelVector.Phi(),first_layer_hits.at(i).Phi()));
         }
      }

      // pixel1 pixel3
      for(unsigned long i = 0; i < first_layer_hits.size(); i++){
         for(unsigned long j = 0; j < third_layer_hits.size(); j++){
           TVector3 pixelVector = third_layer_hits.at(j) - first_layer_hits.at(i);
           TVector3 EM_pixelVector = emvector - third_layer_hits.at(j);
           ntPix13EGdphi.push_back(deltaPhi( EM_pixelVector.Phi(), pixelVector.Phi() ));
           ntPix13EGdeta.push_back(EM_pixelVector.Eta()-pixelVector.Eta() );

           ntPix013dphi.push_back(deltaPhi(pixelVector.Phi(),first_layer_hits.at(i).Phi()));
         }
      }

      // pixel1 pixel4
      for(unsigned long i = 0; i < first_layer_hits.size(); i++){
         for(unsigned long j = 0; j < fourth_layer_hits.size(); j++){
           TVector3 pixelVector = fourth_layer_hits.at(j) - first_layer_hits.at(i);
           TVector3 EM_pixelVector = emvector - fourth_layer_hits.at(j);
           ntPix14EGdphi.push_back(deltaPhi( EM_pixelVector.Phi(), pixelVector.Phi() ));
           ntPix14EGdeta.push_back(EM_pixelVector.Eta()-pixelVector.Eta() );

           ntPix014dphi.push_back(deltaPhi(pixelVector.Phi(),first_layer_hits.at(i).Phi()));
         } 
      }

      // pixel2 pixel3
      for(unsigned long i = 0; i < second_layer_hits.size(); i++){
         for(unsigned long j = 0; j < third_layer_hits.size(); j++){
           TVector3 pixelVector = third_layer_hits.at(j) - second_layer_hits.at(i);
           TVector3 EM_pixelVector = emvector - third_layer_hits.at(j);
           ntPix23EGdphi.push_back(deltaPhi( EM_pixelVector.Phi(), pixelVector.Phi() ));
           ntPix23EGdeta.push_back(EM_pixelVector.Eta()-pixelVector.Eta() );

           ntPix023dphi.push_back(deltaPhi(pixelVector.Phi(),second_layer_hits.at(i).Phi()));
         }
      }

      // pixel2 pixel4
      for(unsigned long i = 0; i < second_layer_hits.size(); i++){
         for(unsigned long j = 0; j < fourth_layer_hits.size(); j++){
           TVector3 pixelVector = fourth_layer_hits.at(j) - second_layer_hits.at(i);
           TVector3 EM_pixelVector = emvector - fourth_layer_hits.at(j);
           ntPix24EGdphi.push_back(deltaPhi( EM_pixelVector.Phi(), pixelVector.Phi() ));
           ntPix24EGdeta.push_back(EM_pixelVector.Eta()- pixelVector.Eta() );

           ntPix024dphi.push_back(deltaPhi(pixelVector.Phi(),second_layer_hits.at(i).Phi()));
         }
      }

      // pixel3 pixel4
      for(unsigned long i = 0; i < third_layer_hits.size(); i++){
         for(unsigned long j = 0; j < fourth_layer_hits.size(); j++){
            TVector3 pixelVector = fourth_layer_hits.at(j) - third_layer_hits.at(i);
            TVector3 EM_pixelVector = emvector - fourth_layer_hits.at(j);
            ntPix34EGdphi.push_back(deltaPhi( EM_pixelVector.Phi(), pixelVector.Phi() ));
            ntPix34EGdeta.push_back(EM_pixelVector.Eta()-pixelVector.Eta() );

            ntPix034dphi.push_back(deltaPhi(pixelVector.Phi(),third_layer_hits.at(i).Phi()));
        }
      }

      // pixel1 pixel2 pixel3
      for(unsigned long i = 0; i < first_layer_hits.size(); i++){
      for(unsigned long j = 0; j < second_layer_hits.size(); j++){
      for(unsigned long k = 0; k < third_layer_hits.size(); k++){ 
        TVector3 pixelVector1 = second_layer_hits.at(j) - first_layer_hits.at(i);
        TVector3 pixelVector2 = third_layer_hits.at(k) - second_layer_hits.at(j);

        ntPix123dphi.push_back(deltaPhi(pixelVector2.Phi(),pixelVector1.Phi()));
        ntPix123deta.push_back(pixelVector2.Eta()-pixelVector1.Eta());
      }
      }
      }

      // pixel1 pixel2 pixel4
      for(unsigned long i = 0; i < first_layer_hits.size(); i++){
      for(unsigned long j = 0; j < second_layer_hits.size(); j++){
      for(unsigned long k = 0; k < fourth_layer_hits.size(); k++){ 
        TVector3 pixelVector1 = second_layer_hits.at(j) - first_layer_hits.at(i);
        TVector3 pixelVector2 = fourth_layer_hits.at(k) - second_layer_hits.at(j);

        ntPix124dphi.push_back(deltaPhi(pixelVector2.Phi(),pixelVector1.Phi()));
        ntPix124deta.push_back(pixelVector2.Eta()-pixelVector1.Eta());
      }
      }
      }

      // pixel1 pixel3 pixel4
      for(unsigned long i = 0; i < first_layer_hits.size(); i++){
      for(unsigned long j = 0; j < third_layer_hits.size(); j++){
      for(unsigned long k = 0; k < fourth_layer_hits.size(); k++){ 
        TVector3 pixelVector1 = third_layer_hits.at(j) - first_layer_hits.at(i);
        TVector3 pixelVector2 = fourth_layer_hits.at(k) - third_layer_hits.at(j);

        ntPix134dphi.push_back(deltaPhi(pixelVector2.Phi(),pixelVector1.Phi()));
        ntPix134deta.push_back(pixelVector2.Eta()-pixelVector1.Eta());
      }
      }
      }

      // pixel2 pixel3 pixel4
      for(unsigned long i = 0; i < second_layer_hits.size(); i++){
      for(unsigned long j = 0; j < third_layer_hits.size(); j++){
      for(unsigned long k = 0; k < fourth_layer_hits.size(); k++){ 
        TVector3 pixelVector1 = third_layer_hits.at(j) - second_layer_hits.at(i);
        TVector3 pixelVector2 = fourth_layer_hits.at(k) - third_layer_hits.at(j);

        ntPix234dphi.push_back(deltaPhi(pixelVector2.Phi(),pixelVector1.Phi()));
        ntPix234deta.push_back(pixelVector2.Eta()-pixelVector1.Eta());
      }
      }
      }


   sw_tree->Fill();
   }
   file->Write();
}
