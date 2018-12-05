#define test_cxx
#include "test.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <string>
#include <TLorentzVector.h>

using namespace std;

void test::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;
  
   float dr_cut = 0.1;

   bit1 = 0x1;

   debug = false;

   const double EM_PiX_dphi_width_[9] = {0.02, 0.03, 0.022, 0.024, 0.026, 0.028, 0.030, 0.032, 0.034}; // 
   const double EM_PiX_deta_width_[9] = {0.01, 0.015, 0.011, 0.012, 0.013, 0.014, 0.015, 0.16, 0.17}; // 

   const double PiX_PiX_dphi_width_[9] = {0.0017, 0.003, 0.0018, 0.0019, 0.0020, 0.0021, 0.0022, 0.0023, 0.0024};
   const double PiX_PiX_deta_width_[9] = {0.0017, 0.003, 0.0018, 0.0019, 0.0020, 0.0021, 0.0022, 0.0023, 0.0024};
   //nentries = 10;
   for (Long64_t jentry=0; jentry<nentries;jentry++) { //nentries
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if (!(jentry%1) ) cout << "Processing entry " << jentry << "/" << nentries << endl;
      FillCutFlow("NoCut", 1.);

      matchedEgEt  = -999.;
      matchedEgEta = -999.;
      matchedEgPhi = -999.;
      fired        = 0;;

      ntnEg2 = 0;
      ntEgEt.clear();
      ntEgEta.clear();
      ntEgPhi.clear();

      ntL1TkEgEt.clear();
      ntL1TkEgEta.clear();
      ntL1TkEgPhi.clear();

      ntL1TkLooseEgEt.clear();
      ntL1TkLooseEgEta.clear();
      ntL1TkLooseEgPhi.clear();

      PiXTRKbit.clear();
      trigger_bit_width.clear();
      pix_comb.clear();

      ntCl_match.clear();
      withoutEM_match.clear();
      withEM_match.clear();

      ntfirstPix.clear();
      ntsecondPix.clear();
      ntthirdPix.clear();
      ntfourthPix.clear();

      nt_EM12.clear();
      nt_EM13.clear();
      nt_EM23.clear();
      nt_EM24.clear();
      nt_EM34.clear();

      nt_SA012.clear();
      nt_SA013.clear();
      nt_SA023.clear();
      nt_SA123.clear();
      nt_SA124.clear();
      nt_SA134.clear();
      nt_SA234.clear();

      nt_genPhi = propgenElPartPhi->at(0);
      nt_genEta = propgenElPartEta->at(0);
      nt_genPt = propgenElPartPt->at(0);

      nt_lastSimtkpt = lastSimtkpt;
      nt_initialSimtkpt = initialSimtkpt;
 
      float tempDR = 999.;
      int   indx = -1;
      int   indx_closestEg = -1;

      //find closest egamma object to the gen electron
      EgN=egCrysClusterEt->size();
      int cl3d_N_ = cl3d_pt->size();
      int tower_N_ = tower_pt->size();

      float closest_dr = 9999.;
      int closest_eg = 0;
      int egCount = 0;

      // loop over barrel egamma objects
      for(int i=0; i < EgN;i++){

         float dPhi = deltaPhi(propgenElPartPhi->at(0), egCrysClusterPhi->at(i));

         float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-egCrysClusterEta->at(i),2));
         if(egCrysClusterEt->at(i) < 10) continue;
         egCount++;
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

         if(cl3d_egid->at(i) != 1) continue;

         float dPhi = deltaPhi(propgenElPartPhi->at(0), cl3d_phi->at(i));

         float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-cl3d_eta->at(i),2));
         if(cl3d_pt->at(i) < 10) continue;
         cl3d_Count++;
         if(current_dr < closest_cl3d_dr){
           closest_cl3d_dr = current_dr;
           closest_cl3d = i;
         }
      }// end of loop to find the closest egamma to gen electron 

      //for(int i=0; i < tower_N_;i++){

      //   float dPhi = deltaPhi(propgenElPartPhi->at(0), tower_phi->at(i));

      //   float current_dr = sqrt(pow(dPhi,2)+pow(propgenElPartEta->at(0)-tower_eta->at(i),2));
      //   if(tower_pt->at(i) < 10) continue;
      //   cl3d_Count++;
      //   if(current_dr < closest_cl3d_dr){
      //     closest_cl3d_dr = current_dr;
      //     closest_cl3d = i;
      //   }
      //}// end of loop to find the closest egamma to gen electron 

      pix_comb_ = 0x0;

      nPix123_segments = 0;
      nPix124_segments = 0;
      nPix134_segments = 0;
      nPix234_segments = 0;

      EM12 = false;
      EM13 = false;
      EM23 = false;
      EM24 = false;
      EM34 = false;

      SA012 = false;
      SA013 = false;
      SA023 = false;
      SA123 = false; 
      SA124 = false; 
      SA134 = false;
      SA234 = false;;

     // find egamma objects passing pixtrk signal windows
     if((closest_cl3d_dr < dr_cut && closest_cl3d_dr != 9999.)|| (closest_dr < dr_cut && closest_dr != 9999.)){


      //save L1TkElectron
      for(int i = 0; i < tkEGEt->size(); i++){
         ntL1TkEgEt.push_back(tkEGEt->at(i));
         ntL1TkEgEta.push_back(tkEGEta->at(i));
         ntL1TkEgPhi.push_back(tkEGPhi->at(i));
      }

      //save L1TkElectron
      for(int i = 0; i < tkEGLooseEt->size(); i++){
         ntL1TkLooseEgEt.push_back(tkEGLooseEt->at(i));
         ntL1TkLooseEgEta.push_back(tkEGLooseEta->at(i));
         ntL1TkLooseEgPhi.push_back(tkEGLoosePhi->at(i));
      }

      debug = false;

      indx++; //remove this variable

      if( closest_dr < closest_cl3d_dr ){
        EgEt =egCrysClusterEt ->at(closest_eg);
        EgEta=egCrysClusterEta->at(closest_eg);
        EgPhi=egCrysClusterPhi->at(closest_eg);

        float EgGx = egCrysClusterGx->at(closest_eg);
        float EgGy = egCrysClusterGy->at(closest_eg);
        float EgGz = egCrysClusterGz->at(closest_eg);
        emvector.SetXYZ(EgGx,EgGy,EgGz);

        tempDR = closest_dr; 
        matchedEgEta = EgEta;
        matchedEgPhi = EgPhi;
        matchedEgEt  = EgEt;
        indx_closestEg = indx;
      }
      else{

          EgEt =cl3d_pt->at(closest_cl3d);
          EgEta=cl3d_eta->at(closest_cl3d);
          EgPhi=cl3d_phi->at(closest_cl3d);

          float EgGx = cl3d_x->at(closest_cl3d);
          float EgGy = cl3d_y->at(closest_cl3d);
          float EgGz = (float)cl3d_z->at(closest_cl3d);
          emvector.SetXYZ(EgGx,EgGy,EgGz);

          tempDR = closest_cl3d_dr; 
          matchedEgEta = EgEta;
          matchedEgPhi = EgPhi;
          matchedEgEt  = EgEt;
          indx_closestEg = indx;

          //EgEt =tower_pt->at(closest_cl3d);
          //EgEta=tower_eta->at(closest_cl3d);
          //EgPhi=tower_phi->at(closest_cl3d);
        
          //float EgGx = tower_gx->at(closest_cl3d);
          //float EgGy = tower_gy->at(closest_cl3d);
          //float EgGz = (float)tower_gz->at(closest_cl3d);
          //emvector.SetXYZ(EgGx,EgGy,EgGz);

          //tempDR = closest_cl3d_dr;
          //matchedEgEta = EgEta;
          //matchedEgPhi = EgPhi;
          //matchedEgEt  = EgEt;
          //indx_closestEg = indx;
      }


      if( fabs(EgEta) <= 0.8 ) eta_region =1;
      if( fabs(EgEta) <= 1.4 && fabs(EgEta) > 0.8 ) eta_region =2;
      if( fabs(EgEta) <= 1.7 && fabs(EgEta) > 1.4 ) eta_region =3;
      if( fabs(EgEta) <= 2.1 && fabs(EgEta) > 1.7 ) eta_region =4;
      if( fabs(EgEta) <= 2.7 && fabs(EgEta) > 2.1 ) eta_region =5;
      if( fabs(EgEta) <= 3.0 && fabs(EgEta) > 2.7 ) eta_region =6;

      if( fabs(EgEta) > 3. ) continue;

      ntnEg2++;
      ntEgEt.push_back(EgEt);
      ntEgEta.push_back(EgEta);
      ntEgPhi.push_back(EgPhi);
      
      // set regin of interest
      SetROI(eta_region);

      // initialize pixel hit variables
      first_layer_hits.clear();
      second_layer_hits.clear();
      third_layer_hits.clear();
      fourth_layer_hits.clear();

      first_layer_hits_Ele_or_Pos.clear();
      second_layer_hits_Ele_or_Pos.clear();
      third_layer_hits_Ele_or_Pos.clear();
      fourth_layer_hits_Ele_or_Pos.clear();
      hitted_layers.clear();

      
      layers[0] = 1; // beam spot
      layers[1] = 0; layers[2] = 0; layers[3] = 0; layers[4] = 0;
      r = 0;

      StorePixelHit(eta_region); // save pixel hits in Region of Interest for the given eta region

      // check which pixel has hits
       for( int i=1; i < 5; i++){ 
          if( layers[i] != 0 ){ 
            hitted_layers.push_back(i); 
          }
       }

       trigger_bit_width_ = 0x0;
       // set pixtrk signal boundary
       for(int nth_eg_pix_deta = 0; nth_eg_pix_deta < 9; nth_eg_pix_deta++){
       if( nth_eg_pix_deta != 0) continue;
      
       if(eta_region == 1) SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta], EM_PiX_deta_width_[nth_eg_pix_deta], PiX_PiX_dphi_width_[nth_eg_pix_deta], PiX_PiX_deta_width_[nth_eg_pix_deta]);
       else SetSingalBoundary(eta_region, EM_PiX_dphi_width_[nth_eg_pix_deta+1], EM_PiX_deta_width_[nth_eg_pix_deta+1], PiX_PiX_dphi_width_[nth_eg_pix_deta+1], PiX_PiX_deta_width_[nth_eg_pix_deta+1]);

       //SetSingalBoundary(1);

       // PixTRK algorithm 
       PixTrkPassed = false;
       withoutEM_count_Ele = 0, withEM_count_Ele = 0;

       fourth_layer_missing = 0;
       third_layer_missing = 0;
       second_layer_missing = 0;
       first_layer_missing = 0;

       // loop over every 3 out of 4 pixel combination 
       for( std::vector<int>::iterator first_hit = hitted_layers.begin(); first_hit != hitted_layers.end(); first_hit++){
          for ( std::vector<int>::iterator second_hit = first_hit+1; second_hit != hitted_layers.end(); second_hit++){
              for ( std::vector<int>::iterator third_hit = second_hit+1; third_hit != hitted_layers.end(); third_hit++){

                 
                 // loop over every pixel hits in the given pixel combination
                 for( int k=0; k < layers[*first_hit]; k++){
                    for( int i=0; i < layers[*second_hit]; i++){
                        _pass_Ele = 0, _pass_Pos = 0;
                        L012_pass_Ele = 0, L012_pass_Pos = 0;
                        L013_pass_Ele = 0, L013_pass_Pos = 0;
                        L023_pass_Ele = 0, L023_pass_Pos = 0;

                        if( *first_hit == 1 && *second_hit == 2 )
                          TriggeringWith_1st2ndPixel(k,i);

                        if( *first_hit == 1 && *second_hit == 3 )
                          TriggeringWith_1st3rdPixel(k,i);

                        if( *first_hit == 2 && *second_hit == 3 )
                          TriggeringWith_2nd3rdPixel(k,i);

                        // skip only if both _pass_Ele and _pass_Pos are 0 i.e., both electron and positron signal window are not satisfied
                        if( !_pass_Ele && !_pass_Pos ) continue;

                        for( int j=0; j < layers[*third_hit]; j++){
                            all_cut_pass_Ele = 0, all_cut_pass_Pos = 0;
                            withoutEM_pass_Ele = 0, withEM_pass_Ele = 0;

                            L012_pass_Ele = 0, L012_pass_Pos = 0;
                            L013_pass_Ele = 0, L013_pass_Pos = 0;
                            L014_pass_Ele = 0, L014_pass_Pos = 0;
                            L023_pass_Ele = 0, L023_pass_Pos = 0;
                            L024_pass_Ele = 0, L024_pass_Pos = 0;
                            L034_pass_Ele = 0, L034_pass_Pos = 0;
                            L123_pass_Ele = 0, L123_pass_Pos = 0;
                            L124_pass_Ele = 0, L124_pass_Pos = 0;
                            L134_pass_Ele = 0, L134_pass_Pos = 0;
                            L234_pass_Ele = 0, L234_pass_Pos = 0;

                            L12_EM_Ele = 0, L12_EM_Pos = 0;
                            L13_EM_Ele = 0, L13_EM_Pos = 0;
                            L14_EM_Ele = 0, L14_EM_Pos = 0;
                            L23_EM_Ele = 0, L23_EM_Pos = 0;
                            L24_EM_Ele = 0, L24_EM_Pos = 0;
                            L34_EM_Ele = 0, L34_EM_Pos = 0;

            	            dPhi = StandaloneDPhi( *first_hit, *second_hit, *third_hit, k, i, j );
                            dEta = StandaloneDEta( *first_hit, *second_hit, *third_hit, k, i, j );

                              if( *first_hit == 1 && *second_hit == 2 && *third_hit == 3 ){ // for efficiency counting  !!caution of the position of this codition
                                // This is for the case that the first hit is in the first pixel layer and the second hit is in the second pixel layer and the third hit is in the third layer. 
                                TriggeringWithout_4thPixel(k, i, j);

                                if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 1 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[j] == 1 || third_layer_hits_Ele_or_Pos[j] ==3)){
                                      if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele && L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
                                      //if(L012_pass_Ele && L123_pass_Ele && L12_EM_Ele && L23_EM_Ele)
                                         all_cut_pass_Ele = 1; 

                                         if(L12_EM_Ele) EM12 = true;
                                         if(L23_EM_Ele) EM23 = true;
                                         if(L012_pass_Ele) SA012 = true;
                                         if(L123_pass_Ele )SA123 = true;

                                      if(L012_pass_Ele && L013_pass_Ele && L023_pass_Ele && L123_pass_Ele)
                                         withoutEM_pass_Ele = 1;
                                      if(L12_EM_Ele && L13_EM_Ele && L23_EM_Ele)
                                         withEM_pass_Ele = 1;
                                }
 
                                if( L012_pass_Pos && L013_pass_Pos && L023_pass_Pos && L123_pass_Pos && L12_EM_Pos && L13_EM_Pos && L23_EM_Pos &&
            			    (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[j] == 2 || third_layer_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

                                if( all_cut_pass_Ele == 1){
                                   pix_comb_ = pix_comb_ | (bit1 << 1);
                                   nPix123_segments++;
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 ){ // if pass exit loop
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
                                  fourth_layer_missing = 1;
                   		 }
                               }
                              }
                              if( *first_hit == 1 && *second_hit == 2 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                                // This is for the case that the first hit is in the first pixel layer and the second is in the second layer and the third hit is in the fourth layer.
                                TriggeringWithout_3rdPixel(k, i, j);

                                if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 1 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 1 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele && L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) all_cut_pass_Ele = 1;
                                    //if(L012_pass_Ele && L124_pass_Ele && L12_EM_Ele && L24_EM_Ele) all_cut_pass_Ele = 1;
                                    if(L012_pass_Ele && L014_pass_Ele && L024_pass_Ele && L124_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L12_EM_Ele && L14_EM_Ele && L24_EM_Ele) withEM_pass_Ele = 1;

                                         if(L12_EM_Ele) EM12 = true;
                                         if(L24_EM_Ele) EM24 = true;
                                         if(L012_pass_Ele) SA012 = true;
                                         if(L124_pass_Ele )SA124 = true;

                                } 

                                if( L012_pass_Pos && L014_pass_Pos && L024_pass_Pos && L124_pass_Pos && L12_EM_Pos && L14_EM_Pos && L24_EM_Pos &&
            			    (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (second_layer_hits_Ele_or_Pos[i] == 2 || second_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

                                if( all_cut_pass_Ele == 1){
                                   pix_comb_ = pix_comb_ | (bit1 << 2);
                                   nPix124_segments++;
                                }


                                if(skip){
                                if( all_cut_pass_Ele == 1 ){ // if pass exit the for loops
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
            	     	          third_layer_missing = 1;
                   		 }
                               }
                              }
                              if( *first_hit == 1 && *second_hit == 3 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                                TriggeringWithout_2ndPixel(k, i, j);

                                if( (first_layer_hits_Ele_or_Pos[k] == 1 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 1 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 1 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele && L13_EM_Ele && L14_EM_Ele && L34_EM_Ele)all_cut_pass_Ele = 1; 
                                    //if(L013_pass_Ele && L134_pass_Ele && L13_EM_Ele && L34_EM_Ele)all_cut_pass_Ele = 1; 
                                    if(L013_pass_Ele && L014_pass_Ele && L034_pass_Ele && L134_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L13_EM_Ele && L14_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;

                                         if(L13_EM_Ele) EM13 = true;
                                         if(L34_EM_Ele) EM34 = true;
                                         if(L013_pass_Ele) SA013 = true;
                                         if(L134_pass_Ele )SA134 = true;
                                }

                                if( L013_pass_Pos && L014_pass_Pos && L034_pass_Pos && L134_pass_Pos && L13_EM_Pos && L14_EM_Pos && L34_EM_Pos &&
            			    (first_layer_hits_Ele_or_Pos[k] == 2 || first_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 

                                if( all_cut_pass_Ele == 1){
                                   pix_comb_ = pix_comb_ | (bit1 << 3);
                                   nPix134_segments++;
                                }

                                if(skip){
                                if( all_cut_pass_Ele == 1 ){ // if pass exit the for loops
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
                                  second_layer_missing = 1; 
                   		 }
                               }
                              }
                              if( *first_hit == 2 && *second_hit == 3 && *third_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                                TriggeringWithout_1stPixel(k, i, j);

                                if( (second_layer_hits_Ele_or_Pos[k] == 1 || second_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 1 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 1 || fourth_layer_hits_Ele_or_Pos[j] ==3)){
                                    if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele && L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) all_cut_pass_Ele = 1;
                                    //if(L023_pass_Ele && L234_pass_Ele && L23_EM_Ele && L34_EM_Ele) all_cut_pass_Ele = 1;
                                    if(L023_pass_Ele && L024_pass_Ele && L034_pass_Ele && L234_pass_Ele) withoutEM_pass_Ele = 1;
                                    if(L23_EM_Ele && L24_EM_Ele && L34_EM_Ele) withEM_pass_Ele = 1;
                                         if(L23_EM_Ele) EM23 = true;
                                         if(L34_EM_Ele) EM34 = true;
                                         if(L023_pass_Ele) SA023 = true;
                                         if(L234_pass_Ele)SA234 = true;
                                } 

                                if( L023_pass_Pos && L024_pass_Pos && L034_pass_Pos && L234_pass_Pos && L23_EM_Pos && L24_EM_Pos && L34_EM_Pos &&
            			    (second_layer_hits_Ele_or_Pos[k] == 2 || second_layer_hits_Ele_or_Pos[k] ==3) &&
            			    (third_layer_hits_Ele_or_Pos[i] == 2 || third_layer_hits_Ele_or_Pos[i] ==3) &&
            			    (fourth_layer_hits_Ele_or_Pos[j] == 2 || fourth_layer_hits_Ele_or_Pos[j] ==3)) all_cut_pass_Pos = 1; 
                   
                                if( all_cut_pass_Ele == 1){
                                   pix_comb_ = pix_comb_ | (bit1 << 4);
                                   nPix234_segments++;
                                } 

                                if(skip){
                                if( all_cut_pass_Ele == 1 ){ // if pass exit the for loops
            		          k = layers[*first_hit];
                                  i = layers[*second_hit];
                                  j = layers[*third_hit]; 
                                  first_layer_missing = 1;
                   		 } 
                               }
                              }

                         if( all_cut_pass_Ele == 1 ) { PixTrkPassed = true;}
                         if( withoutEM_pass_Ele == 1 ) withoutEM_count_Ele = 1;
                         if( withEM_pass_Ele == 1 ) withEM_count_Ele = 1;
                       } // loop for third layer hits
                   } // loop for second layer hits       
                 } // loop for first layer hits

          }          
        }
      }
/*
     if( fabs(EgEta) <= 1.4 && fabs(EgEta) > 1.2 && PixTrkPassed == false) {
       for( std::vector<int>::iterator first_hit = hitted_layers.begin(); first_hit != hitted_layers.end(); first_hit++){
          for ( std::vector<int>::iterator second_hit = first_hit+1; second_hit != hitted_layers.end(); second_hit++){
                 
                 // loop over every pixel hits in the given pixel combination
                 for( int k=0; k < layers[*first_hit]; k++){
                    for( int i=0; i < layers[*second_hit]; i++){
                        _pass_Ele = 0, _pass_Pos = 0;

                        L012_pass_Ele = 0, L012_pass_Pos = 0;
                        L013_pass_Ele = 0, L013_pass_Pos = 0;
                        L014_pass_Ele = 0, L014_pass_Pos = 0;
                        L023_pass_Ele = 0, L023_pass_Pos = 0;
                        L024_pass_Ele = 0, L024_pass_Pos = 0;
                        L034_pass_Ele = 0, L034_pass_Pos = 0;
                        L123_pass_Ele = 0, L123_pass_Pos = 0;
                        L124_pass_Ele = 0, L124_pass_Pos = 0;
                        L134_pass_Ele = 0, L134_pass_Pos = 0;
                        L234_pass_Ele = 0, L234_pass_Pos = 0;

                        L12_EM_Ele = 0, L12_EM_Pos = 0;
                        L13_EM_Ele = 0, L13_EM_Pos = 0;
                        L14_EM_Ele = 0, L14_EM_Pos = 0;
                        L23_EM_Ele = 0, L23_EM_Pos = 0;
                        L24_EM_Ele = 0, L24_EM_Pos = 0;
                        L34_EM_Ele = 0, L34_EM_Pos = 0;

                        // skip only if both _pass_Ele and _pass_Pos are 0 i.e., both electron and positron signal window are not satisfied

                        if( *first_hit == 1 && *second_hit == 2 ){ // for efficiency counting  !!caution of the position of this codition

                          TriggeringWith_1st2ndPixel_v2(k,i); 

                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
            		    k = layers[*first_hit];
                            i = layers[*second_hit];
                   	   }
                         }
                        }

                        if( *first_hit == 1 && *second_hit == 3 ){ // for efficiency counting  !!caution of the position of this codition
                        
                          TriggeringWith_1st3rdPixel_v2(k,i);
                          
                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
                            k = layers[*first_hit];
                            i = layers[*second_hit];
                           }
                         } 
                        }

                        if( *first_hit == 1 && *second_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                        
                          TriggeringWith_1st4thPixel_v2(k,i);
                          
                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
                            k = layers[*first_hit];
                            i = layers[*second_hit];
                           }
                         } 
                        }

                        if( *first_hit == 2 && *second_hit == 3 ){ // for efficiency counting  !!caution of the position of this codition
                        
                          TriggeringWith_2nd3rdPixel_v2(k,i);
                          
                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
                            k = layers[*first_hit];
                            i = layers[*second_hit];
                           }
                         } 
                        }

                        if( *first_hit == 2 && *second_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition

                          TriggeringWith_2nd4thPixel_v2(k,i);

                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
                            k = layers[*first_hit];
                            i = layers[*second_hit];
                           }
                         }
                        }

                        if( *first_hit == 3 && *second_hit == 4 ){ // for efficiency counting  !!caution of the position of this codition
                        
                          TriggeringWith_3rd4thPixel_v2(k,i);
                          
                          if(skip){
                          if( _pass_Ele == 1 || _pass_Pos == 1 ){ // if pass exit loop
                            k = layers[*first_hit];
                            i = layers[*second_hit];
                           }
                         } 
                        }

                       if( _pass_Ele == 1 ) { PixTrkPassed = true;}
                   } // loop for second layer hits       
                 } // loop for first layer hits

        }
      }
     }
*/
      if( PixTrkPassed ){ 
          trigger_bit_width_ = trigger_bit_width_| (bit1 << nth_eg_pix_deta);
      }

      /////////////////////////////////
     }// end of signal window width loop
  } // end of egamma loop    

  nt_EM12.push_back(EM12); 
  nt_EM13.push_back(EM13); 
  nt_EM23.push_back(EM23); 
  nt_EM24.push_back(EM24); 
  nt_EM34.push_back(EM34); 
  
  nt_SA012.push_back(SA012);
  nt_SA013.push_back(SA013);
  nt_SA023.push_back(SA023);
  nt_SA123.push_back(SA123);
  nt_SA124.push_back(SA124);
  nt_SA134.push_back(SA134);
  nt_SA234.push_back(SA234);

  trigger_bit_width.push_back(trigger_bit_width_);
  pix_comb.push_back(pix_comb_);

  pixtrk_tree->Fill();
 } // end of entries loop 
   outfile->Write();
}

