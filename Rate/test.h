//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 17 20:40:59 2015 by ROOT version 5.34/19
// from TTree t/t
// found on file: SingleEle_ntuple_1.root
//////////////////////////////////////////////////////////

#ifndef test_h
#define test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
#include <iostream>
#include <fstream>
#include <string>

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include "../../../roi_v2.h"
#include "../../../withEM_v2.h"
#include "../../../withoutEM_v2.h"

using namespace std;
class test {
private:

map<TString, TH1*> maphist;

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           bunchN;
   vector<int>     *pileup;
   Float_t         beamSpotX0;
   Float_t         beamSpotY0;
   Float_t         beamSpotZ0;
   Float_t         beamSpotX0Error;
   Float_t         beamSpotY0Error;
   Float_t         beamSpotZ0Error;
   Float_t         beamWidthX;
   Float_t         beamWidthY;
   Float_t         beamSigmaZ;
   Float_t         beamWidthXError;
   Float_t         beamWidthYError;
   Float_t         beamSigmaZError;
   Int_t           genPartN;
   vector<float>   *genPartE;
   vector<float>   *genPartPt;
   vector<float>   *genPartEta;
   vector<float>   *genPartPhi;
   vector<int>     *genPartCharge;
   vector<int>     *genPartId;
   vector<float>   *propgenElPartE;
   vector<float>   *propgenElPartPt;
   vector<float>   *propgenElPartEta;
   vector<float>   *propgenElPartPhi;
   vector<int>     *propgenElPartCharge;
   vector<int>     *propgenElPartId;
   Int_t           simTrkN;
   vector<float>   *simTrkPt;
   vector<float>   *simTrkEta;
   vector<float>   *simTrkPhi;
   vector<int>     *simTrkId;
   vector<int>     *simTrkType;
   vector<float>   *simTrkVx;
   vector<float>   *simTrkVy;
   vector<float>   *simTrkVz;
   vector<float>   *simVx;
   vector<float>   *simVy;
   vector<float>   *simVz;
   Float_t         gammaBrem;
   Int_t           hardGamma;
   Int_t           egCrysN;
   vector<float>   *egCrysE;
   vector<float>   *egCrysEt;
   vector<float>   *egCrysEta;
   vector<float>   *egCrysPhi;
   vector<float>   *egCrysGx;
   vector<float>   *egCrysGy;
   vector<float>   *egCrysGz;
   vector<int>     *egCrysCharge;
   vector<bool>    *isTrackMatched;
   vector<float>   *isoConeNTrack;
   vector<float>   *isoConePtTrack;
   vector<float>   *trackHighestPt;
   vector<float>   *trackHighestPtEta;
   vector<float>   *trackHighestPtPhi;
   vector<float>   *trackHighestPtChi2;
   vector<float>   *trackHighestPtCutChi2;
   vector<float>   *trackHighestPtCutChi2Eta;
   vector<float>   *trackHighestPtCutChi2Phi;
   vector<float>   *trackHighestPtCutChi2Chi2;
   vector<float>   *trackmatchingdR;
   vector<bool>    *hgcal_isTrackMatched;
   vector<float>   *hgcal_isoConeNTrack;
   vector<float>   *hgcal_isoConePtTrack;
   vector<float>   *hgcal_trackHighestPt;
   vector<float>   *hgcal_trackHighestPtEta;
   vector<float>   *hgcal_trackHighestPtPhi;
   vector<float>   *hgcal_trackHighestPtChi2;
   vector<float>   *hgcal_trackHighestPtCutChi2;
   vector<float>   *hgcal_trackHighestPtCutChi2Eta;
   vector<float>   *hgcal_trackHighestPtCutChi2Phi;
   vector<float>   *hgcal_trackHighestPtCutChi2Chi2;
   vector<float>   *hgcal_trackmatchingdR;
   Int_t           cl3d_n;
   vector<float>   *cl3d_pt;
   vector<float>   *cl3d_egid;
   vector<float>   *cl3d_energy;
   vector<float>   *cl3d_eta;
   vector<float>   *cl3d_phi;
   vector<int>     *cl3d_nclu;
   vector<float>   *cl3d_x;
   vector<float>   *cl3d_y;
   vector<int>     *cl3d_z;
   vector<int>     *cl3d_coreshowerlength;
   vector<float>   *cl3d_srrtot;
   vector<int>     *cl3d_maxlayer;
   vector<int>     *cl3d_firstlayer;
   Int_t           l1tkegN;
   vector<float>   *l1tkegEt;
   vector<float>   *l1tkegEta;
   vector<float>   *l1tkegPhi;
   Int_t           l1tkegIsoN;
   vector<float>   *l1tkegIsoEt;
   vector<float>   *l1tkegIsoEta;
   vector<float>   *l1tkegIsoPhi;
   Int_t           egN;
   vector<float>   *egE;
   vector<float>   *egEt;
   vector<float>   *egEta;
   vector<float>   *egPhi;
   vector<float>   *egGx;
   vector<float>   *egGy;
   vector<float>   *egGz;
   vector<float>   *egclusterBS;
   vector<float>   *egclusterPt;
   vector<float>   *egclusterPtCore;
   vector<float>   *egclusterEta;
   vector<float>   *egclusterPhi;
   vector<double>  *egHoE;
   vector<double>  *clusterHoE;
   vector<double>  *seedHoE;
   vector<bool>    *clusterIsEG;
   vector<bool>    *seedIsEG;
   vector<bool>    *TriggerTower_IsBarrel;
   vector<double>  *seedHadEt;
   vector<double>  *seedEmEt;
   vector<double>  *HadEt;
   vector<double>  *EmEt;
   vector<int>     *bRecHitLayer;
   vector<int>     *bRecHitLadder;
   vector<int>     *bRecHitModule;
   vector<int>     *fRecHitDisk;
   vector<int>     *fRecHitBlade;
   vector<int>     *fRecHitSide;
   vector<int>     *fRecHitPanel;
   vector<int>     *fRecHitModule;
   Int_t           bRecHitN;
   Int_t           fRecHitN;
   vector<float>   *fRecHitGx;
   vector<float>   *fRecHitGy;
   vector<float>   *fRecHitGz;
   vector<float>   *fRhSize;
   vector<float>   *fRhSizeX;
   vector<float>   *fRhSizeY;
   vector<float>   *bRecHitGx;
   vector<float>   *bRecHitGy;
   vector<float>   *bRecHitGz;
   vector<float>   *bRhSize;
   vector<float>   *bRhSizeX;
   vector<float>   *bRhSizeY;

   UInt_t          ntkEG;
   vector<double>  *tkEGEt;
   vector<double>  *tkEGEta;
   vector<double>  *tkEGPhi;
   vector<int>     *tkEGBx;
   vector<double>  *tkEGTrkIso;
   vector<double>  *tkEGzVtx;
   vector<double>  *tkEGHwQual;
   vector<double>  *tkEGEGRefPt;
   vector<double>  *tkEGEGRefEta;
   vector<double>  *tkEGEGRefPhi;
   UInt_t          ntkEGLoose;
   vector<double>  *tkEGLooseEt;
   vector<double>  *tkEGLooseEta;
   vector<double>  *tkEGLoosePhi;
   vector<int>     *tkEGLooseBx;
   vector<double>  *tkEGLooseTrkIso;
   vector<double>  *tkEGLoosezVtx;
   vector<double>  *tkEGLooseHwQual;
   vector<double>  *tkEGLooseEGRefPt;
   vector<double>  *tkEGLooseEGRefEta;
   vector<double>  *tkEGLooseEGRefPhi;

   Int_t           tower_n;
   vector<float>   *tower_pt;
   vector<float>   *tower_energy;
   vector<float>   *tower_eta;
   vector<float>   *tower_phi;
   vector<float>   *tower_etEm;
   vector<float>   *tower_etHad;
   vector<float>   *tower_gx;
   vector<float>   *tower_gy;
   vector<float>   *tower_gz;
   vector<int>     *tower_iEta;
   vector<int>     *tower_iPhi;
   

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bunchN;   //!
   TBranch        *b_pileup;   //!
   TBranch        *b_beamSpotX0;   //!
   TBranch        *b_beamSpotY0;   //!
   TBranch        *b_beamSpotZ0;   //!
   TBranch        *b_beamSpotX0Error;   //!
   TBranch        *b_beamSpotY0Error;   //!
   TBranch        *b_beamSpotZ0Error;   //!
   TBranch        *b_beamWidthX;   //!
   TBranch        *b_beamWidthY;   //!
   TBranch        *b_beamSigmaZ;   //!
   TBranch        *b_beamWidthXError;   //!
   TBranch        *b_beamWidthYError;   //!
   TBranch        *b_beamSigmaZError;   //!
   TBranch        *b_genPartN;   //!
   TBranch        *b_genPartE;   //!
   TBranch        *b_genPartPt;   //!
   TBranch        *b_genPartEta;   //!
   TBranch        *b_genPartPhi;   //!
   TBranch        *b_genPartCharge;   //!
   TBranch        *b_genPartId;   //!
   TBranch        *b_propgenElPartE;   //!
   TBranch        *b_propgenElPartPt;   //!
   TBranch        *b_propgenElPartEta;   //!
   TBranch        *b_propgenElPartPhi;   //!
   TBranch        *b_propgenElPartCharge;   //!
   TBranch        *b_propgenElPartId;   //!
   TBranch        *b_simTrkN;   //!
   TBranch        *b_simTrkPt;   //!
   TBranch        *b_simTrkEta;   //!
   TBranch        *b_simTrkPhi;   //!
   TBranch        *b_simTrkId;   //!
   TBranch        *b_simTrkType;   //!
   TBranch        *b_simTrkVx;   //!
   TBranch        *b_simTrkVy;   //!
   TBranch        *b_simTrkVz;   //!
   TBranch        *b_simVx;   //!
   TBranch        *b_simVy;   //!
   TBranch        *b_simVz;   //!
   TBranch        *b_gammaBrem;   //!
   TBranch        *b_hardGamma;   //!
   TBranch        *b_egCrysN;   //!
   TBranch        *b_egCrysE;   //!
   TBranch        *b_egCrysEt;   //!
   TBranch        *b_egCrysEta;   //!
   TBranch        *b_egCrysPhi;   //!
   TBranch        *b_egCrysGx;   //!
   TBranch        *b_egCrysGy;   //!
   TBranch        *b_egCrysGz;   //!
   TBranch        *b_egCrysCharge;   //!
   TBranch        *b_isTrackMatched;   //!
   TBranch        *b_isoConeNTrack;   //!
   TBranch        *b_isoConePtTrack;   //!
   TBranch        *b_trackHighestPt;   //!
   TBranch        *b_trackHighestPtEta;   //!
   TBranch        *b_trackHighestPtPhi;   //!
   TBranch        *b_trackHighestPtChi2;   //!
   TBranch        *b_trackHighestPtCutChi2;   //!
   TBranch        *b_trackHighestPtCutChi2Eta;   //!
   TBranch        *b_trackHighestPtCutChi2Phi;   //!
   TBranch        *b_trackHighestPtCutChi2Chi2;   //!
   TBranch        *b_trackmatchingdR;
   TBranch        *b_hgcal_isTrackMatched;   //!
   TBranch        *b_hgcal_isoConeNTrack;   //!
   TBranch        *b_hgcal_isoConePtTrack;   //!
   TBranch        *b_hgcal_trackHighestPt;   //!
   TBranch        *b_hgcal_trackHighestPtEta;   //!
   TBranch        *b_hgcal_trackHighestPtPhi;   //!
   TBranch        *b_hgcal_trackHighestPtChi2;   //!
   TBranch        *b_hgcal_trackHighestPtCutChi2;   //!
   TBranch        *b_hgcal_trackHighestPtCutChi2Eta;   //!
   TBranch        *b_hgcal_trackHighestPtCutChi2Phi;   //!
   TBranch        *b_hgcal_trackHighestPtCutChi2Chi2;   //!
   TBranch        *b_hgcal_trackmatchingdR;
   TBranch        *b_cl3d_n;   //!
   TBranch        *b_cl3d_pt;   //!
   TBranch        *b_cl3d_egid;   //!
   TBranch        *b_cl3d_energy;   //!
   TBranch        *b_cl3d_eta;   //!
   TBranch        *b_cl3d_phi;   //!
   TBranch        *b_cl3d_nclu;   //!
   TBranch        *b_cl3d_x;   //!
   TBranch        *b_cl3d_y;   //!
   TBranch        *b_cl3d_z;   //!
   TBranch        *b_cl3d_coreshowerlength;
   TBranch        *b_cl3d_srrtot;
   TBranch        *b_cl3d_maxlayer;
   TBranch        *b_cl3d_firstlayer;
   TBranch        *b_l1tkegN;   //!
   TBranch        *b_l1tkegEt;   //!
   TBranch        *b_l1tkegEta;   //!
   TBranch        *b_l1tkegPhi;   //!
   TBranch        *b_l1tkegIsoN;   //!
   TBranch        *b_l1tkegIsoEt;   //!
   TBranch        *b_l1tkegIsoEta;   //!
   TBranch        *b_l1tkegIsoPhi;   //!
   TBranch        *b_egN;   //!
   TBranch        *b_egE;   //!
   TBranch        *b_egEt;   //!
   TBranch        *b_egEta;   //!
   TBranch        *b_egPhi;   //!
   TBranch        *b_egGx;   //!
   TBranch        *b_egGy;   //!
   TBranch        *b_egGz;   //!
   TBranch        *b_egclusterBS;   //!
   TBranch        *b_egclusterPt;   //!
   TBranch        *b_egclusterPtCore;   //!
   TBranch        *b_egclusterEta;   //!
   TBranch        *b_egclusterPhi;   //!
   TBranch        *b_egHoE;   //!
   TBranch        *b_clusterHoE;   //!
   TBranch        *b_seedHoE;   //!
   TBranch        *b_clusterIsEG;   //!
   TBranch        *b_seedIsEG;   //!
   TBranch        *b_TriggerTower_IsBarrel;   //!
   TBranch        *b_seedHadEt;   //!
   TBranch        *b_seedEmEt;   //!
   TBranch        *b_HadEt;   //!
   TBranch        *b_EmEt;   //!
   TBranch        *b_bRecHitLayer;   //!
   TBranch        *b_bRecHitLadder;   //!
   TBranch        *b_bRecHitModule;   //!
   TBranch        *b_fRecHitDisk;   //!
   TBranch        *b_fRecHitBlade;   //!
   TBranch        *b_fRecHitSide;   //!
   TBranch        *b_fRecHitPanel;   //!
   TBranch        *b_fRecHitModule;   //!
   TBranch        *b_bRecHitN;   //!
   TBranch        *b_fRecHitN;   //!
   TBranch        *b_fRecHitGx;   //!
   TBranch        *b_fRecHitGy;   //!
   TBranch        *b_fRecHitGz;   //!
   TBranch        *b_fRhSize;   //!
   TBranch        *b_fRhSizeX;   //!
   TBranch        *b_fRhSizeY;   //!
   TBranch        *b_bRecHitGx;   //!
   TBranch        *b_bRecHitGy;   //!
   TBranch        *b_bRecHitGz;   //!
   TBranch        *b_bRhSize;   //!
   TBranch        *b_bRhSizeX;   //!
   TBranch        *b_bRhSizeY;   //!
   TBranch        *b_tower_n;   //!
   TBranch        *b_tower_pt;   //!
   TBranch        *b_tower_energy;   //!
   TBranch        *b_tower_eta;   //!
   TBranch        *b_tower_phi;   //!
   TBranch        *b_tower_etEm;   //!
   TBranch        *b_tower_etHad;   //!
   TBranch        *b_tower_gx;   //!
   TBranch        *b_tower_gy;   //!
   TBranch        *b_tower_gz;   //!
   TBranch        *b_tower_iEta;   //!
   TBranch        *b_tower_iPhi;   //!
   TBranch        *b_ntkEG;   //!
   TBranch        *b_tkEGEt;   //!
   TBranch        *b_tkEGEta;   //!
   TBranch        *b_tkEGPhi;   //!
   TBranch        *b_tkEGBx;   //!
   TBranch        *b_tkEGTrkIso;   //!
   TBranch        *b_tkEGzVtx;   //!
   TBranch        *b_tkEGHwQual;   //!
   TBranch        *b_tkEGEGRefPt;   //!
   TBranch        *b_tkEGEGRefEta;   //!
   TBranch        *b_tkEGEGRefPhi;   //!
   TBranch        *b_ntkEGLoose;   //!
   TBranch        *b_tkEGLooseEt;   //!
   TBranch        *b_tkEGLooseEta;   //!
   TBranch        *b_tkEGLoosePhi;   //!
   TBranch        *b_tkEGLooseBx;   //!
   TBranch        *b_tkEGLooseTrkIso;   //!
   TBranch        *b_tkEGLoosezVtx;   //!
   TBranch        *b_tkEGLooseHwQual;   //!
   TBranch        *b_tkEGLooseEGRefPt;   //!
   TBranch        *b_tkEGLooseEGRefEta;   //!
   TBranch        *b_tkEGLooseEGRefPhi;   //!

   int trigger_bit_width_;
   int PiXTRKbit_4Hits_; 

   int Ele, Pos;
   int skip;
   float EgN;
   int eta_region;

   int withoutEM_count_Ele, withEM_count_Ele;
   double pass_count, pass_count_EleorPos, pass_count_Ele, pass_count_Pos;
   int pass_count_wo4thPix, pass_count_wo3thPix, pass_count_wo2thPix, pass_count_wo1thPix;
   int woEM_pass_count_wo4thPix, woEM_pass_count_wo3thPix, woEM_pass_count_wo2thPix, woEM_pass_count_wo1thPix;
   int wEM_pass_count_wo4thPix, wEM_pass_count_wo3thPix, wEM_pass_count_wo2thPix, wEM_pass_count_wo1thPix;

   double all_cut_pass_eg;
   int all_cut_pass_Ele, withoutEM_pass_Ele, withEM_pass_Ele, withoutEM_pass_Pos, withEM_pass_Pos;
   int all_cut_pass_Pos;
   int fourth_layer_missing;
   int third_layer_missing;
   int second_layer_missing;
   int first_layer_missing;

   double  L1_Dphi_cut1, L1_Dphi_cut2;
   double  L2_Dphi_cut1, L2_Dphi_cut2;
   double  L3_Dphi_cut1, L3_Dphi_cut2;
   double  L4_Dphi_cut1, L4_Dphi_cut2;
   double  D1_Dphi_cut1, D1_Dphi_cut2;
   double  D2_Dphi_cut1, D2_Dphi_cut2;
   double  D3_Dphi_cut1, D3_Dphi_cut2;

   double dPhi012;
   double dPhi013;
   double dPhi014;
   double dPhi023;
   double dPhi024;
   double dPhi034;

   double  L012_DPhi_cut1, L012_DPhi_cut2;
   double  L012_DEta_cut1, L012_DEta_cut2;
   double  L012_DR_cut1, L012_DR_cut2;
   
   double  L013_DPhi_cut1, L013_DPhi_cut2;
   double  L013_DEta_cut1, L013_DEta_cut2;
   double  L013_DR_cut1, L013_DR_cut2;
   
   double  L014_DPhi_cut1, L014_DPhi_cut2;
   double  L014_DEta_cut1, L014_DEta_cut2;
   double  L014_DR_cut1, L014_DR_cut2;
   
   double  L023_DPhi_cut1, L023_DPhi_cut2;
   double  L023_DEta_cut1, L023_DEta_cut2;
   double  L023_DR_cut1, L023_DR_cut2;
   
   double  L024_DPhi_cut1, L024_DPhi_cut2;
   double  L024_DEta_cut1, L024_DEta_cut2;
   double  L024_DR_cut1, L024_DR_cut2;
   
   double  L034_DPhi_cut1, L034_DPhi_cut2;
   double  L034_DEta_cut1, L034_DEta_cut2;
   double  L034_DR_cut1, L034_DR_cut2;
   
   double  L123_DPhi_cut1, L123_DPhi_cut2;
   double  L123_DEta_cut1, L123_DEta_cut2;
   double  L123_DR_cut1, L123_DR_cut2;
   
   double  L124_DPhi_cut1, L124_DPhi_cut2;
   double  L124_DEta_cut1, L124_DEta_cut2;
   double  L124_DR_cut1, L124_DR_cut2;
   
   double  L134_DPhi_cut1, L134_DPhi_cut2;
   double  L134_DEta_cut1, L134_DEta_cut2;
   double  L134_DR_cut1, L134_DR_cut2;
   
   double  L234_DPhi_cut1, L234_DPhi_cut2;
   double  L234_DEta_cut1, L234_DEta_cut2;
   double  L234_DR_cut1, L234_DR_cut2;

   double  L12_eta_upper, L13_eta_upper, L14_eta_upper, L23_eta_upper, L24_eta_upper, L34_eta_upper;
   double  L12_phi_upper, L13_phi_upper, L14_phi_upper, L23_phi_upper, L24_phi_upper, L34_phi_upper;
   double  L12_R_upper, L13_R_upper, L14_R_upper, L23_R_upper, L24_R_upper, L34_R_upper;
   double  L12_eta_bellow, L13_eta_bellow, L14_eta_bellow, L23_eta_bellow, L24_eta_bellow, L34_eta_bellow;
   double  L12_phi_bellow, L13_phi_bellow, L14_phi_bellow, L23_phi_bellow, L24_phi_bellow, L34_phi_bellow;
   double  L12_R_bellow, L13_R_bellow, L14_R_bellow, L23_R_bellow, L24_R_bellow, L34_R_bellow;

   double dPhi;
   double dEta;
   double dR;
   double dPhi_1, dPhi_2, dPhi_3;
   double dEta_1, dEta_2, dEta_3;
   double dR_1, dR_2, dR_3;
   TVector3 first_temp, second_temp;
   int _pass_Ele, _pass_Pos;

   int L012_pass_Ele, L012_pass_Pos; 
   int L013_pass_Ele, L013_pass_Pos;
   int L014_pass_Ele, L014_pass_Pos;
   int L023_pass_Ele, L023_pass_Pos;
   int L024_pass_Ele, L024_pass_Pos;
   int L034_pass_Ele, L034_pass_Pos;
   int L123_pass_Ele, L123_pass_Pos;
   int L124_pass_Ele, L124_pass_Pos;
   int L134_pass_Ele, L134_pass_Pos;
   int L234_pass_Ele, L234_pass_Pos;

   int L12_EM_Ele, L12_EM_Pos; 
   int L13_EM_Ele, L13_EM_Pos;
   int L14_EM_Ele, L14_EM_Pos;
   int L23_EM_Ele, L23_EM_Pos;
   int L24_EM_Ele, L24_EM_Pos;
   int L34_EM_Ele, L34_EM_Pos;

   std::vector<TVector3> first_layer_hits;
   std::vector<TVector3> second_layer_hits;
   std::vector<TVector3> third_layer_hits;
   std::vector<TVector3> fourth_layer_hits;

   double r; // r for radius of pixel tracker layer
   int layers[5];  // initialize as 0, layers contain # of hits on each pixel layer

   TVector3 emvector;
   float EgEt;
   float EgEta;
   float EgPhi;

   std::vector<int> first_layer_hits_Ele_or_Pos;
   std::vector<int> second_layer_hits_Ele_or_Pos;
   std::vector<int> third_layer_hits_Ele_or_Pos;
   std::vector<int> fourth_layer_hits_Ele_or_Pos;
   std::vector<int> hitted_layers;

   void MakeHistograms(TString hname, int nbins, float xmin, float xmax);
   TH1* GetHist(TString hname);
   void FillHist(TString histname, float value, float w, float xmin, float xmax, int nbins);

   void StorePixelHit( int region);
   double StandaloneDPhi( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit );
   double StandaloneDEta( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit );
   double EMmatchingDEta(TVector3& first_layer, TVector3& second_layer, TVector3& egvector);
   double EMmatchingDPhi(TVector3& first_layer, TVector3& second_layer, TVector3& egvector);
   int Signal_window_check( double upper, double value, double lower, int Ele_Pos);
   void FillCutFlow(TString cut, float weight);
   void SetROI(int region);
   void SetSingalBoundary(int region, double eg_dphi, double eg_deta, double sa_dphi, double sa_deta);
   void TriggeringWith_1st2ndPixel(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_1st3rdPixel(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_2nd3rdPixel(int nthFirstHit, int nthSecondHit);
   void TriggeringWithout_4thPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);
   void TriggeringWithout_3rdPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);
   void TriggeringWithout_2ndPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);
   void TriggeringWithout_1stPixel(int nthFirstHit, int nthSecondHit, int nthThirdHit);

   void TriggeringWith_1st2ndPixel_v2(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_1st3rdPixel_v2(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_1st4thPixel_v2(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_2nd3rdPixel_v2(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_2nd4thPixel_v2(int nthFirstHit, int nthSecondHit);
   void TriggeringWith_3rd4thPixel_v2(int nthFirstHit, int nthSecondHit);

   inline float deltaPhi(float phi1, float phi2){
     float result = phi1 - phi2;
     while(result > float(M_PI)) result -= float(2*M_PI);
     while(result <= -float(M_PI)) result += float(2*M_PI);
    return result;
   }

   test(TTree *tree=0);
   virtual ~test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   TFile *file3;
   TTree* pixtrk_tree;

   int count_Entry; 
   int pass_egobjects_check;
   int ntnEg2; 
   int event_denominator; 
   int event_nominator; 

   vector<float> ntEgEt; 
   vector<float> ntL1TkEleEt, ntL1TkEleIsoEt; 
   vector<float> ntEgEta; 
   vector<bool> ntclusterIsEG; 
   vector<float> ntL1TkEleEta, ntL1TkEleIsoEta; 
   vector<float> ntL1TkElePhi, ntL1TkEleIsoPhi;
   vector<bool> ntCl_match; 
   vector<bool> ntCl_iso_match; 
   vector<bool> withoutEM_match; 
   vector<bool> withEM_match; 
   vector<bool> isTrack_match;
   vector<float> chi2;
   vector<float> track_dr;

   vector<float> ntL1TkEgEt;
   vector<float> ntL1TkEgEta;
   vector<float> ntL1TkEgPhi;

   vector<float> ntL1TkLooseEgEt;
   vector<float> ntL1TkLooseEgEta;
   vector<float> ntL1TkLooseEgPhi;

   vector<bool> ntCl_match_wo4thPix, ntCl_match_wo3thPix, ntCl_match_wo2thPix, ntCl_match_wo1thPix;
   vector<int> Npass_woEM_wo4thPix, Npass_woEM_wo3thPix, Npass_woEM_wo2thPix, Npass_woEM_wo1thPix;
   vector<int> Npass_wEM_wo4thPix, Npass_wEM_wo3thPix, Npass_wEM_wo2thPix, Npass_wEM_wo1thPix;

   vector<float> dphi_L12EM_wo4thPix, dphi_L12EM_wo3thPix, dphi_L12EM_wo2thPix, dphi_L12EM_wo1thPix;
   vector<float> deta_L12EM_wo4thPix, deta_L12EM_wo3thPix, deta_L12EM_wo2thPix, deta_L12EM_wo1thPix;

   vector<int> pass_Ele; 
   vector<int> pass_Pos; 
   vector<int> pass_ElePos; 

   vector<int> ntfirstPix;
   vector<int> ntsecondPix;
   vector<int> ntthirdPix;
   vector<int> ntfourthPix;

   int L1TkEleN;
   int L1TkEleIsoN;

   vector<int> trigger_bit_width;
   vector<int> trigger_4Hits_bit_width;
   
   // For track isolation
   
   class track
   {
       public:
           float pos_x1; float pos_y1; float pos_z1;
           float pos_x2; float pos_y2; float pos_z2;
           float pos_x3; float pos_y3; float pos_z3;
           int index;

           track() 
           { 
               pos_x1 = 0; pos_y1 = 0; pos_z1 = 0; 
               pos_x2 = 0; pos_y2 = 0; pos_z2 = 0; 
               pos_x3 = 0; pos_y3 = 0; pos_z3 = 0; 
           }
           track(float a, float b, float c, float d, float e, float f, float g, float h, float i, int k) 
           { 
               pos_x3 = a; pos_x2 = b; pos_x1 = c; 
               pos_y3 = d; pos_y2 = e; pos_y1 = f; 
               pos_z3 = g; pos_z2 = h; pos_z1 = i; 
               index = k;
           }

           static bool comp3(const track &t1, const track &t2)
           {
               return ( t1.pos_x3 < t2.pos_x3 );
           }
           static bool uni3(const track &t1, const track &t2)
           {
               return( t1.pos_x3 == t2.pos_x3 );
           }

           static bool comp2(const track &t1, const track &t2)
           {
               return ( t1.pos_x2 < t2.pos_x2 );
           }
           static bool uni2(const track &t1, const track &t2)
           {
               return( t1.pos_x2 == t2.pos_x2 );
           }

           static bool comp1(const track &t1, const track &t2)
           {
               return ( t1.pos_x1 < t2.pos_x1 );
           }
           static bool uni1(const track &t1, const track &t2)
           {
               return( t1.pos_x1 == t2.pos_x1 );
           }
           
           static bool uni12(const track &t1, const track &t2)
           {
               return( t1.pos_x1 == t2.pos_x2 );
           }
           
           static bool uni21(const track &t1, const track &t2)
           {
               return( t1.pos_x2 == t2.pos_x1 );
           }
           
           static bool uni13(const track &t1, const track &t2)
           {
               return( t1.pos_x1 == t2.pos_x3 );
           }
           
           static bool uni31(const track &t1, const track &t2)
           {
               return( t1.pos_x3 == t2.pos_x1 );
           }
           
           static bool uni23(const track &t1, const track &t2)
           {
               return( t1.pos_x2 == t2.pos_x3 );
           }
           
           static bool uni32(const track &t1, const track &t2)
           {
               return( t1.pos_x3 == t2.pos_x2 );
           }
   };
   
   std::vector<TVector3> first_hits;
   std::vector<TVector3> second_hits;
   std::vector<TVector3> third_hits;
   std::vector<TVector3> fourth_hits;
  
   void StorePixelHitForIso(int region, float recoPV);
   void IsoWith_1st2nd3rd(int region, float recoPV);
   void IsoWith_1st2nd4th(int region, float recoPV);
   void IsoWith_1st3rd4th(int region, float recoPV);
   void IsoWith_2nd3rd4th(int region, float recoPV);

   std::vector<track> L123;
   std::vector<track> L124;
   std::vector<track> L134;
   std::vector<track> L234;

};

#endif

#ifdef test_cxx
test::test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TChain * chain = new TChain("l1PiXTRKTree/L1PiXTRKTree","");
      string line;
      ifstream myfile("txt_to_path"); // txt_to_path will be replaced by the name of txt file that contains the location of input files
      if (myfile.is_open())
        {
          while ( getline (myfile,line) )
          {
            if( line.length() != 0 ){
              cout << line << '\n';
              char file_path[300];
              strcpy(file_path, line.c_str());
              chain->Add(file_path);
            }
          }
          myfile.close();
        }
      tree = chain;
   }
   Init(tree);
 
   Ele = 1, Pos = 2;
   skip = 1;

   file3 = new TFile("../output_tmp/Tree_output","recreate");
   pixtrk_tree = new TTree("t","t");

   count_Entry = 1;
   pixtrk_tree->Branch("totalEvent", &count_Entry, "count_Entry/I");   pixtrk_tree->Branch("totalEgN", &EgN, "EgN/F");
   pixtrk_tree->Branch("ntnEg2", &ntnEg2, "ntnEg2/I");
   pixtrk_tree->Branch("event_with_selection_passed", &event_denominator, "event_denominator/I");
   pixtrk_tree->Branch("event_with_trigger_passed", &event_nominator, "event_nominator/I");

   pixtrk_tree->Branch("L1TkEleN", &L1TkEleN, "L1TkEleN/I");
   pixtrk_tree->Branch("L1TkEleIsoN", &L1TkEleIsoN, "L1TkEleIsoN/I");
   pixtrk_tree->Branch("ntEgEt",&ntEgEt);
   pixtrk_tree->Branch("ntL1TkEleEt",&ntL1TkEleEt);
   pixtrk_tree->Branch("ntL1TkEleIsoEt",&ntL1TkEleIsoEt);
   pixtrk_tree->Branch("ntEgEta",&ntEgEta);
   pixtrk_tree->Branch("ntclusterIsEG",&ntclusterIsEG);
   pixtrk_tree->Branch("ntL1TkEleEta",&ntL1TkEleEta);
   pixtrk_tree->Branch("ntL1TkEleIsoEta",&ntL1TkEleIsoEta);
   pixtrk_tree->Branch("ntL1TkElePhi",&ntL1TkElePhi);
   pixtrk_tree->Branch("ntL1TkEleIsoPhi",&ntL1TkEleIsoPhi);
   pixtrk_tree->Branch("ntCl_match",&ntCl_match);
   pixtrk_tree->Branch("ntCl_iso_match",&ntCl_iso_match);
   pixtrk_tree->Branch("withoutEM_match",&withoutEM_match);
   pixtrk_tree->Branch("withEM_match",&withEM_match);
   pixtrk_tree->Branch("isTrack_match",&isTrack_match);
   pixtrk_tree->Branch("chi2",&chi2);
   pixtrk_tree->Branch("track_dr",&track_dr);

   pixtrk_tree->Branch("ntCl_match_wo4thPix",&ntCl_match_wo4thPix);
   pixtrk_tree->Branch("ntCl_match_wo3thPix",&ntCl_match_wo3thPix);
   pixtrk_tree->Branch("ntCl_match_wo2thPix",&ntCl_match_wo2thPix);
   pixtrk_tree->Branch("ntCl_match_wo1thPix",&ntCl_match_wo1thPix);

   pixtrk_tree->Branch("Npass_woEM_wo4thPix", &Npass_woEM_wo4thPix);
   pixtrk_tree->Branch("Npass_woEM_wo3thPix", &Npass_woEM_wo3thPix);
   pixtrk_tree->Branch("Npass_woEM_wo2thPix", &Npass_woEM_wo2thPix);
   pixtrk_tree->Branch("Npass_woEM_wo1thPix", &Npass_woEM_wo1thPix);

   pixtrk_tree->Branch("Npass_wEM_wo4thPix", &Npass_wEM_wo4thPix);
   pixtrk_tree->Branch("Npass_wEM_wo3thPix", &Npass_wEM_wo3thPix);
   pixtrk_tree->Branch("Npass_wEM_wo2thPix", &Npass_wEM_wo2thPix);
   pixtrk_tree->Branch("Npass_wEM_wo1thPix", &Npass_wEM_wo1thPix);

   pixtrk_tree->Branch("deta_L12EM_wo4thPix", &deta_L12EM_wo4thPix);
   pixtrk_tree->Branch("deta_L12EM_wo3thPix", &deta_L12EM_wo3thPix);
   pixtrk_tree->Branch("deta_L12EM_wo2thPix", &deta_L12EM_wo2thPix);
   pixtrk_tree->Branch("deta_L12EM_wo1thPix", &deta_L12EM_wo1thPix);

   pixtrk_tree->Branch("dphi_L12EM_wo4thPix", &dphi_L12EM_wo4thPix);
   pixtrk_tree->Branch("dphi_L12EM_wo3thPix", &dphi_L12EM_wo3thPix);
   pixtrk_tree->Branch("dphi_L12EM_wo2thPix", &dphi_L12EM_wo2thPix);
   pixtrk_tree->Branch("dphi_L12EM_wo1thPix", &dphi_L12EM_wo1thPix);

   pixtrk_tree->Branch("pass_Ele", &pass_Ele);
   pixtrk_tree->Branch("pass_Pos", &pass_Pos);
   pixtrk_tree->Branch("pass_ElePos", &pass_ElePos);

   pixtrk_tree->Branch("ntfirstPix",&ntfirstPix);
   pixtrk_tree->Branch("ntsecondPix",&ntsecondPix);
   pixtrk_tree->Branch("ntthirdPix",&ntthirdPix);
   pixtrk_tree->Branch("ntfourthPix",&ntfourthPix);

   pixtrk_tree->Branch("trigger_bit_width",&trigger_bit_width);
   pixtrk_tree->Branch("trigger_4Hits_bit_width",&trigger_4Hits_bit_width);


   pixtrk_tree->Branch("ntL1TkEgEt",&ntL1TkEgEt);
   pixtrk_tree->Branch("ntL1TkEgEta",&ntL1TkEgEta);
   pixtrk_tree->Branch("ntL1TkEgPhi",&ntL1TkEgPhi);

   pixtrk_tree->Branch("ntL1TkLooseEgEt",&ntL1TkLooseEgEt);
   pixtrk_tree->Branch("ntL1TkLooseEgEta",&ntL1TkLooseEgEta);
   pixtrk_tree->Branch("ntL1TkLooseEgPhi",&ntL1TkLooseEgPhi);
}

test::~test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t test::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void test::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pileup = 0;
   genPartE = 0;
   genPartPt = 0;
   genPartEta = 0;
   genPartPhi = 0;
   genPartCharge = 0;
   genPartId = 0;
   propgenElPartE = 0;
   propgenElPartPt = 0;
   propgenElPartEta = 0;
   propgenElPartPhi = 0;
   propgenElPartCharge = 0;
   propgenElPartId = 0;
   simTrkPt = 0;
   simTrkEta = 0;
   simTrkPhi = 0;
   simTrkId = 0;
   simTrkType = 0;
   simTrkVx = 0;
   simTrkVy = 0;
   simTrkVz = 0;
   simVx = 0;
   simVy = 0;
   simVz = 0;
   egCrysE = 0;
   egCrysEt = 0;
   egCrysEta = 0;
   egCrysPhi = 0;
   egCrysGx = 0;
   egCrysGy = 0;
   egCrysGz = 0;
   egCrysCharge = 0;
   isTrackMatched = 0;
   isoConeNTrack = 0;
   isoConePtTrack = 0;
   trackHighestPt = 0;
   trackHighestPtEta = 0;
   trackHighestPtPhi = 0;
   trackHighestPtChi2 = 0;
   trackHighestPtCutChi2 = 0;
   trackHighestPtCutChi2Eta = 0;
   trackHighestPtCutChi2Phi = 0;
   trackHighestPtCutChi2Chi2 = 0;
   trackmatchingdR = 0;
   hgcal_isTrackMatched = 0;
   hgcal_isoConeNTrack = 0;
   hgcal_isoConePtTrack = 0;
   hgcal_trackHighestPt = 0;
   hgcal_trackHighestPtEta = 0;
   hgcal_trackHighestPtPhi = 0;
   hgcal_trackHighestPtChi2 = 0;
   hgcal_trackHighestPtCutChi2 = 0;
   hgcal_trackHighestPtCutChi2Eta = 0;
   hgcal_trackHighestPtCutChi2Phi = 0;
   hgcal_trackHighestPtCutChi2Chi2 = 0;
   hgcal_trackmatchingdR = 0;
   cl3d_pt = 0;
   cl3d_egid = 0;
   cl3d_energy = 0;
   cl3d_eta = 0;
   cl3d_phi = 0;
   cl3d_nclu = 0;
   cl3d_x = 0;
   cl3d_y = 0;
   cl3d_z = 0;
   cl3d_coreshowerlength = 0;
   cl3d_srrtot = 0;
   cl3d_maxlayer = 0;
   cl3d_firstlayer = 0;
   egE = 0;
   egEt = 0;
   egEta = 0;
   egPhi = 0;
   egGx = 0;
   egGy = 0;
   egGz = 0;
   egclusterBS = 0;
   egclusterPt = 0;
   egclusterPtCore = 0;
   egclusterEta = 0;
   egclusterPhi = 0;
   egHoE = 0;
   clusterHoE = 0;
   seedHoE = 0;
   clusterIsEG = 0;
   seedIsEG = 0;
   TriggerTower_IsBarrel = 0;
   seedHadEt = 0;
   seedEmEt = 0;
   HadEt = 0;
   EmEt = 0;
   bRecHitLayer = 0;
   bRecHitLadder = 0;
   bRecHitModule = 0;
   fRecHitDisk = 0;
   fRecHitBlade = 0;
   fRecHitSide = 0;
   fRecHitPanel = 0;
   fRecHitModule = 0;
   fRecHitGx = 0;
   fRecHitGy = 0;
   fRecHitGz = 0;
   fRhSize = 0;
   fRhSizeX = 0;
   fRhSizeY = 0;
   bRecHitGx = 0;
   bRecHitGy = 0;
   bRecHitGz = 0;
   bRhSize = 0;
   bRhSizeX = 0;
   bRhSizeY = 0;
   tower_pt = 0;
   tower_energy = 0;
   tower_eta = 0;
   tower_phi = 0;
   tower_etEm = 0;
   tower_etHad = 0;
   tower_gx = 0;
   tower_gy = 0;
   tower_gz = 0;
   tower_iEta = 0;
   tower_iPhi = 0;
   tkEGEt = 0;
   tkEGEta = 0;
   tkEGPhi = 0;
   tkEGBx = 0;
   tkEGTrkIso = 0;
   tkEGzVtx = 0;
   tkEGHwQual = 0;
   tkEGEGRefPt = 0;
   tkEGEGRefEta = 0;
   tkEGEGRefPhi = 0;
   tkEGLooseEt = 0;
   tkEGLooseEta = 0;
   tkEGLoosePhi = 0;
   tkEGLooseBx = 0;
   tkEGLooseTrkIso = 0;
   tkEGLoosezVtx = 0;
   tkEGLooseHwQual = 0;
   tkEGLooseEGRefPt = 0;
   tkEGLooseEGRefEta = 0;
   tkEGLooseEGRefPhi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
//   fChain->SetBranchAddress("bunchN", &bunchN, &b_bunchN);
   fChain->SetBranchAddress("pileup", &pileup, &b_pileup);
   fChain->SetBranchAddress("beamSpotX0", &beamSpotX0, &b_beamSpotX0);
   fChain->SetBranchAddress("beamSpotY0", &beamSpotY0, &b_beamSpotY0);
   fChain->SetBranchAddress("beamSpotZ0", &beamSpotZ0, &b_beamSpotZ0);
   fChain->SetBranchAddress("beamSpotX0Error", &beamSpotX0Error, &b_beamSpotX0Error);
   fChain->SetBranchAddress("beamSpotY0Error", &beamSpotY0Error, &b_beamSpotY0Error);
   fChain->SetBranchAddress("beamSpotZ0Error", &beamSpotZ0Error, &b_beamSpotZ0Error);
   fChain->SetBranchAddress("beamWidthX", &beamWidthX, &b_beamWidthX);
   fChain->SetBranchAddress("beamWidthY", &beamWidthY, &b_beamWidthY);
   fChain->SetBranchAddress("beamSigmaZ", &beamSigmaZ, &b_beamSigmaZ);
   fChain->SetBranchAddress("beamWidthXError", &beamWidthXError, &b_beamWidthXError);
   fChain->SetBranchAddress("beamWidthYError", &beamWidthYError, &b_beamWidthYError);
   fChain->SetBranchAddress("beamSigmaZError", &beamSigmaZError, &b_beamSigmaZError);
//   fChain->SetBranchAddress("genPartN", &genPartN, &b_genPartN);
   fChain->SetBranchAddress("genPartE", &genPartE, &b_genPartE);
   fChain->SetBranchAddress("genPartPt", &genPartPt, &b_genPartPt);
   fChain->SetBranchAddress("genPartEta", &genPartEta, &b_genPartEta);
   fChain->SetBranchAddress("genPartPhi", &genPartPhi, &b_genPartPhi);
   fChain->SetBranchAddress("genPartCharge", &genPartCharge, &b_genPartCharge);
   fChain->SetBranchAddress("genPartId", &genPartId, &b_genPartId);
   fChain->SetBranchAddress("propgenElPartE", &propgenElPartE, &b_propgenElPartE);
   fChain->SetBranchAddress("propgenElPartPt", &propgenElPartPt, &b_propgenElPartPt);
   fChain->SetBranchAddress("propgenElPartEta", &propgenElPartEta, &b_propgenElPartEta);
   fChain->SetBranchAddress("propgenElPartPhi", &propgenElPartPhi, &b_propgenElPartPhi);
   fChain->SetBranchAddress("propgenElPartCharge", &propgenElPartCharge, &b_propgenElPartCharge);
   fChain->SetBranchAddress("propgenElPartId", &propgenElPartId, &b_propgenElPartId);
   fChain->SetBranchAddress("simTrkN", &simTrkN, &b_simTrkN);
   fChain->SetBranchAddress("simTrkPt", &simTrkPt, &b_simTrkPt);
   fChain->SetBranchAddress("simTrkEta", &simTrkEta, &b_simTrkEta);
   fChain->SetBranchAddress("simTrkPhi", &simTrkPhi, &b_simTrkPhi);
   fChain->SetBranchAddress("simTrkId", &simTrkId, &b_simTrkId);
   fChain->SetBranchAddress("simTrkType", &simTrkType, &b_simTrkType);
   fChain->SetBranchAddress("simTrkVx", &simTrkVx, &b_simTrkVx);
   fChain->SetBranchAddress("simTrkVy", &simTrkVy, &b_simTrkVy);
   fChain->SetBranchAddress("simTrkVz", &simTrkVz, &b_simTrkVz);
   fChain->SetBranchAddress("simVx", &simVx, &b_simVx);
   fChain->SetBranchAddress("simVy", &simVy, &b_simVy);
   fChain->SetBranchAddress("simVz", &simVz, &b_simVz);
//   fChain->SetBranchAddress("gammaBrem", &gammaBrem, &b_gammaBrem);
//   fChain->SetBranchAddress("hardGamma", &hardGamma, &b_hardGamma);
   fChain->SetBranchAddress("egCrysN", &egCrysN, &b_egCrysN);
   fChain->SetBranchAddress("egCrysE", &egCrysE, &b_egCrysE);
   fChain->SetBranchAddress("egCrysEt", &egCrysEt, &b_egCrysEt);
   fChain->SetBranchAddress("egCrysEta", &egCrysEta, &b_egCrysEta);
   fChain->SetBranchAddress("egCrysPhi", &egCrysPhi, &b_egCrysPhi);
   fChain->SetBranchAddress("egCrysGx", &egCrysGx, &b_egCrysGx);
   fChain->SetBranchAddress("egCrysGy", &egCrysGy, &b_egCrysGy);
   fChain->SetBranchAddress("egCrysGz", &egCrysGz, &b_egCrysGz);
//   fChain->SetBranchAddress("egCrysCharge", &egCrysCharge, &b_egCrysCharge);
   fChain->SetBranchAddress("isTrackMatched", &isTrackMatched, &b_isTrackMatched);
   fChain->SetBranchAddress("isoConeNTrack", &isoConeNTrack, &b_isoConeNTrack);
   fChain->SetBranchAddress("isoConePtTrack", &isoConePtTrack, &b_isoConePtTrack);
   fChain->SetBranchAddress("trackHighestPt", &trackHighestPt, &b_trackHighestPt);
   fChain->SetBranchAddress("trackHighestPtEta", &trackHighestPtEta, &b_trackHighestPtEta);
   fChain->SetBranchAddress("trackHighestPtPhi", &trackHighestPtPhi, &b_trackHighestPtPhi);
   fChain->SetBranchAddress("trackHighestPtChi2", &trackHighestPtChi2, &b_trackHighestPtChi2);
   fChain->SetBranchAddress("trackHighestPtCutChi2", &trackHighestPtCutChi2, &b_trackHighestPtCutChi2);
   fChain->SetBranchAddress("trackHighestPtCutChi2Eta", &trackHighestPtCutChi2Eta, &b_trackHighestPtCutChi2Eta);
   fChain->SetBranchAddress("trackHighestPtCutChi2Phi", &trackHighestPtCutChi2Phi, &b_trackHighestPtCutChi2Phi);
   fChain->SetBranchAddress("trackHighestPtCutChi2Chi2", &trackHighestPtCutChi2Chi2, &b_trackHighestPtCutChi2Chi2);
   fChain->SetBranchAddress("trackmatchingdR", &trackmatchingdR, &b_trackmatchingdR);
   fChain->SetBranchAddress("hgcal_isTrackMatched", &hgcal_isTrackMatched, &b_hgcal_isTrackMatched);
   fChain->SetBranchAddress("hgcal_isoConeNTrack", &hgcal_isoConeNTrack, &b_hgcal_isoConeNTrack);
   fChain->SetBranchAddress("hgcal_isoConePtTrack", &hgcal_isoConePtTrack, &b_hgcal_isoConePtTrack);
   fChain->SetBranchAddress("hgcal_trackHighestPt", &hgcal_trackHighestPt, &b_hgcal_trackHighestPt);
   fChain->SetBranchAddress("hgcal_trackHighestPtEta", &hgcal_trackHighestPtEta, &b_hgcal_trackHighestPtEta);
   fChain->SetBranchAddress("hgcal_trackHighestPtPhi", &hgcal_trackHighestPtPhi, &b_hgcal_trackHighestPtPhi);
   fChain->SetBranchAddress("hgcal_trackHighestPtChi2", &hgcal_trackHighestPtChi2, &b_hgcal_trackHighestPtChi2);
   fChain->SetBranchAddress("hgcal_trackHighestPtCutChi2", &hgcal_trackHighestPtCutChi2, &b_hgcal_trackHighestPtCutChi2);
   fChain->SetBranchAddress("hgcal_trackHighestPtCutChi2Eta", &hgcal_trackHighestPtCutChi2Eta, &b_hgcal_trackHighestPtCutChi2Eta);
   fChain->SetBranchAddress("hgcal_trackHighestPtCutChi2Phi", &hgcal_trackHighestPtCutChi2Phi, &b_hgcal_trackHighestPtCutChi2Phi);
   fChain->SetBranchAddress("hgcal_trackHighestPtCutChi2Chi2", &hgcal_trackHighestPtCutChi2Chi2, &b_hgcal_trackHighestPtCutChi2Chi2);
   fChain->SetBranchAddress("hgcal_trackmatchingdR", &hgcal_trackmatchingdR, &b_hgcal_trackmatchingdR);
   fChain->SetBranchAddress("cl3d_n", &cl3d_n, &b_cl3d_n);
   fChain->SetBranchAddress("cl3d_pt", &cl3d_pt, &b_cl3d_pt);
   fChain->SetBranchAddress("cl3d_egid", &cl3d_egid, &b_cl3d_egid);
   fChain->SetBranchAddress("cl3d_energy", &cl3d_energy, &b_cl3d_energy);
   fChain->SetBranchAddress("cl3d_eta", &cl3d_eta, &b_cl3d_eta);
   fChain->SetBranchAddress("cl3d_phi", &cl3d_phi, &b_cl3d_phi);
   fChain->SetBranchAddress("cl3d_nclu", &cl3d_nclu, &b_cl3d_nclu);
   fChain->SetBranchAddress("cl3d_x", &cl3d_x, &b_cl3d_x);
   fChain->SetBranchAddress("cl3d_y", &cl3d_y, &b_cl3d_y);
   fChain->SetBranchAddress("cl3d_z", &cl3d_z, &b_cl3d_z);
   fChain->SetBranchAddress("cl3d_coreshowerlength", &cl3d_coreshowerlength, &b_cl3d_coreshowerlength);
   fChain->SetBranchAddress("cl3d_srrtot", &cl3d_srrtot, &b_cl3d_srrtot);
   fChain->SetBranchAddress("cl3d_maxlayer", &cl3d_maxlayer, &b_cl3d_maxlayer);
   fChain->SetBranchAddress("cl3d_firstlayer", &cl3d_firstlayer, &b_cl3d_firstlayer);
   fChain->SetBranchAddress("egHoE", &egHoE, &b_egHoE);
   fChain->SetBranchAddress("clusterHoE", &clusterHoE, &b_clusterHoE);
   fChain->SetBranchAddress("seedHoE", &seedHoE, &b_seedHoE);
   fChain->SetBranchAddress("clusterIsEG", &clusterIsEG, &b_clusterIsEG);
   fChain->SetBranchAddress("seedIsEG", &seedIsEG, &b_seedIsEG);
   fChain->SetBranchAddress("TriggerTower_IsBarrel", &TriggerTower_IsBarrel, &b_TriggerTower_IsBarrel);
   fChain->SetBranchAddress("seedHadEt", &seedHadEt, &b_seedHadEt);
   fChain->SetBranchAddress("seedEmEt", &seedEmEt, &b_seedEmEt);
   fChain->SetBranchAddress("HadEt", &HadEt, &b_HadEt);
   fChain->SetBranchAddress("EmEt", &EmEt, &b_EmEt);
   fChain->SetBranchAddress("l1tkegN", &l1tkegN, &b_l1tkegN);
   fChain->SetBranchAddress("l1tkegEta", &l1tkegEta, &b_l1tkegEta);
   fChain->SetBranchAddress("l1tkegEt", &l1tkegEt, &b_l1tkegEt);
   fChain->SetBranchAddress("l1tkegPhi", &l1tkegPhi, &b_l1tkegPhi);
   fChain->SetBranchAddress("l1tkegIsoN", &l1tkegIsoN, &b_l1tkegIsoN);
   fChain->SetBranchAddress("l1tkegIsoEta", &l1tkegIsoEta, &b_l1tkegIsoEta);
   fChain->SetBranchAddress("l1tkegIsoEt", &l1tkegIsoEt, &b_l1tkegIsoEt);
   fChain->SetBranchAddress("l1tkegIsoPhi", &l1tkegIsoPhi, &b_l1tkegIsoPhi);
   fChain->SetBranchAddress("egN", &egN, &b_egN);
   fChain->SetBranchAddress("egE", &egE, &b_egE);
   fChain->SetBranchAddress("egEt", &egEt, &b_egEt);
   fChain->SetBranchAddress("egEta", &egEta, &b_egEta);
   fChain->SetBranchAddress("egPhi", &egPhi, &b_egPhi);
   fChain->SetBranchAddress("egGx", &egGx, &b_egGx);
   fChain->SetBranchAddress("egGy", &egGy, &b_egGy);
   fChain->SetBranchAddress("egGz", &egGz, &b_egGz);
//   fChain->SetBranchAddress("egclusterBS", &egclusterBS, &b_egclusterBS);
//   fChain->SetBranchAddress("egclusterPt", &egclusterPt, &b_egclusterPt);
//   fChain->SetBranchAddress("egclusterPtCore", &egclusterPtCore, &b_egclusterPtCore);
//   fChain->SetBranchAddress("egclusterEta", &egclusterEta, &b_egclusterEta);
//   fChain->SetBranchAddress("egclusterPhi", &egclusterPhi, &b_egclusterPhi);
   fChain->SetBranchAddress("bRecHitLayer", &bRecHitLayer, &b_bRecHitLayer);
   fChain->SetBranchAddress("bRecHitLadder", &bRecHitLadder, &b_bRecHitLadder);
   fChain->SetBranchAddress("bRecHitModule", &bRecHitModule, &b_bRecHitModule);
   fChain->SetBranchAddress("fRecHitDisk", &fRecHitDisk, &b_fRecHitDisk);
   fChain->SetBranchAddress("fRecHitBlade", &fRecHitBlade, &b_fRecHitBlade);
   fChain->SetBranchAddress("fRecHitSide", &fRecHitSide, &b_fRecHitSide);
   fChain->SetBranchAddress("fRecHitPanel", &fRecHitPanel, &b_fRecHitPanel);
   fChain->SetBranchAddress("fRecHitModule", &fRecHitModule, &b_fRecHitModule);
   fChain->SetBranchAddress("bRecHitN", &bRecHitN, &b_bRecHitN);
   fChain->SetBranchAddress("fRecHitN", &fRecHitN, &b_fRecHitN);
   fChain->SetBranchAddress("fRecHitGx", &fRecHitGx, &b_fRecHitGx);
   fChain->SetBranchAddress("fRecHitGy", &fRecHitGy, &b_fRecHitGy);
   fChain->SetBranchAddress("fRecHitGz", &fRecHitGz, &b_fRecHitGz);
   fChain->SetBranchAddress("fRhSize", &fRhSize, &b_fRhSize);
   fChain->SetBranchAddress("fRhSizeX", &fRhSizeX, &b_fRhSizeX);
   fChain->SetBranchAddress("fRhSizeY", &fRhSizeY, &b_fRhSizeY);
   fChain->SetBranchAddress("bRecHitGx", &bRecHitGx, &b_bRecHitGx);
   fChain->SetBranchAddress("bRecHitGy", &bRecHitGy, &b_bRecHitGy);
   fChain->SetBranchAddress("bRecHitGz", &bRecHitGz, &b_bRecHitGz);
   fChain->SetBranchAddress("bRhSize", &bRhSize, &b_bRhSize);
   fChain->SetBranchAddress("bRhSizeX", &bRhSizeX, &b_bRhSizeX);
   fChain->SetBranchAddress("bRhSizeY", &bRhSizeY, &b_bRhSizeY);

   fChain->SetBranchAddress("tower_n", &tower_n, &b_tower_n);
   fChain->SetBranchAddress("tower_pt", &tower_pt, &b_tower_pt);
   fChain->SetBranchAddress("tower_energy", &tower_energy, &b_tower_energy);
   fChain->SetBranchAddress("tower_eta", &tower_eta, &b_tower_eta);
   fChain->SetBranchAddress("tower_phi", &tower_phi, &b_tower_phi);
   fChain->SetBranchAddress("tower_etEm", &tower_etEm, &b_tower_etEm);
   fChain->SetBranchAddress("tower_etHad", &tower_etHad, &b_tower_etHad);
   fChain->SetBranchAddress("tower_gx", &tower_gx, &b_tower_gx);
   fChain->SetBranchAddress("tower_gy", &tower_gy, &b_tower_gy);
   fChain->SetBranchAddress("tower_gz", &tower_gz, &b_tower_gz);
   fChain->SetBranchAddress("tower_iEta", &tower_iEta, &b_tower_iEta);
   fChain->SetBranchAddress("tower_iPhi", &tower_iPhi, &b_tower_iPhi);
   fChain->SetBranchAddress("ntkEG", &ntkEG, &b_ntkEG);
   fChain->SetBranchAddress("tkEGEt", &tkEGEt, &b_tkEGEt);
   fChain->SetBranchAddress("tkEGEta", &tkEGEta, &b_tkEGEta);
   fChain->SetBranchAddress("tkEGPhi", &tkEGPhi, &b_tkEGPhi);
   fChain->SetBranchAddress("tkEGBx", &tkEGBx, &b_tkEGBx);
   fChain->SetBranchAddress("tkEGTrkIso", &tkEGTrkIso, &b_tkEGTrkIso);
   fChain->SetBranchAddress("tkEGzVtx", &tkEGzVtx, &b_tkEGzVtx);
   fChain->SetBranchAddress("tkEGHwQual", &tkEGHwQual, &b_tkEGHwQual);
   fChain->SetBranchAddress("tkEGEGRefPt", &tkEGEGRefPt, &b_tkEGEGRefPt);
   fChain->SetBranchAddress("tkEGEGRefEta", &tkEGEGRefEta, &b_tkEGEGRefEta);
   fChain->SetBranchAddress("tkEGEGRefPhi", &tkEGEGRefPhi, &b_tkEGEGRefPhi);
   fChain->SetBranchAddress("ntkEGLoose", &ntkEGLoose, &b_ntkEGLoose);
   fChain->SetBranchAddress("tkEGLooseEt", &tkEGLooseEt, &b_tkEGLooseEt);
   fChain->SetBranchAddress("tkEGLooseEta", &tkEGLooseEta, &b_tkEGLooseEta);
   fChain->SetBranchAddress("tkEGLoosePhi", &tkEGLoosePhi, &b_tkEGLoosePhi);
   fChain->SetBranchAddress("tkEGLooseBx", &tkEGLooseBx, &b_tkEGLooseBx);
   fChain->SetBranchAddress("tkEGLooseTrkIso", &tkEGLooseTrkIso, &b_tkEGLooseTrkIso);
   fChain->SetBranchAddress("tkEGLoosezVtx", &tkEGLoosezVtx, &b_tkEGLoosezVtx);
   fChain->SetBranchAddress("tkEGLooseHwQual", &tkEGLooseHwQual, &b_tkEGLooseHwQual);
   fChain->SetBranchAddress("tkEGLooseEGRefPt", &tkEGLooseEGRefPt, &b_tkEGLooseEGRefPt);
   fChain->SetBranchAddress("tkEGLooseEGRefEta", &tkEGLooseEGRefEta, &b_tkEGLooseEGRefEta);
   fChain->SetBranchAddress("tkEGLooseEGRefPhi", &tkEGLooseEGRefPhi, &b_tkEGLooseEGRefPhi);

   Notify();
}

Bool_t test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void test::MakeHistograms(TString hname, int nbins, float xmin, float xmax){

 maphist[hname] =  new TH1F(hname.Data(),hname.Data(),nbins,xmin,xmax);

}

TH1* test::GetHist(TString hname){

 TH1* h = NULL;
 std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
 if(mapit != maphist.end()) return mapit->second;

 return h;

}

void test::FillHist(TString histname, float value, float w, float xmin, float xmax, int nbins){

 if(GetHist(histname)) GetHist(histname)->Fill(value, w);
 else{
//     cout << "Making histogram..." << endl;
     MakeHistograms(histname, nbins, xmin, xmax);
     if(GetHist(histname)) GetHist(histname)->Fill(value, w);
 }
}


double test::StandaloneDPhi( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit ){
  if( first_hit == 0 ){

     TVector3 temp;

     if( second_hit == 1 ) temp.SetXYZ( first_layer_hits[which_second_hit].X(), first_layer_hits[which_second_hit].Y(), first_layer_hits[which_second_hit].Z() );
     if( second_hit == 2 ) temp.SetXYZ( second_layer_hits[which_second_hit].X(), second_layer_hits[which_second_hit].Y(), second_layer_hits[which_second_hit].Z() );
     if( second_hit == 3 ) temp.SetXYZ( third_layer_hits[which_second_hit].X(), third_layer_hits[which_second_hit].Y(), third_layer_hits[which_second_hit].Z() );


     if( third_hit == 2 ) return deltaPhi( (second_layer_hits[which_third_hit] - temp).Phi(), temp.Phi());
     if( third_hit == 3 ) return deltaPhi( (third_layer_hits[which_third_hit] - temp).Phi(), temp.Phi());
     if( third_hit == 4 ) return deltaPhi( (fourth_layer_hits[which_third_hit] - temp).Phi(), temp.Phi());

  }
  if( first_hit != 0 ){
    TVector3 temp_first_layer;
    TVector3 temp_second_layer;
    TVector3 temp_third_layer;

    if( first_hit == 1 ) temp_first_layer = first_layer_hits[which_first_hit];
    if( first_hit == 2 ) temp_first_layer = second_layer_hits[which_first_hit];

    if( second_hit == 2 ) temp_second_layer = second_layer_hits[which_second_hit];
    if( second_hit == 3 ) temp_second_layer = third_layer_hits[which_second_hit];


    if( third_hit == 3 ) temp_third_layer = third_layer_hits[which_third_hit];
    if( third_hit == 4 ) temp_third_layer = fourth_layer_hits[which_third_hit];

    return deltaPhi( (temp_third_layer - temp_second_layer).Phi(), (temp_second_layer - temp_first_layer).Phi());
  }
  return 0.;
}

double test::StandaloneDEta( int first_hit, int second_hit, int third_hit, int which_first_hit, int which_second_hit, int which_third_hit ){

    TVector3 temp_first_layer;
    TVector3 temp_second_layer;
    TVector3 temp_third_layer;

    if( first_hit == 1 ) temp_first_layer = first_layer_hits[which_first_hit];
    if( first_hit == 2 ) temp_first_layer = second_layer_hits[which_first_hit];

    if( second_hit == 2 ) temp_second_layer = second_layer_hits[which_second_hit];
    if( second_hit == 3 ) temp_second_layer = third_layer_hits[which_second_hit];


    if( third_hit == 3 ) temp_third_layer = third_layer_hits[which_third_hit];
    if( third_hit == 4 ) temp_third_layer = fourth_layer_hits[which_third_hit];

    return (temp_third_layer - temp_second_layer).PseudoRapidity() - (temp_second_layer - temp_first_layer).PseudoRapidity();
}

double test::EMmatchingDEta(TVector3& first_layer, TVector3& second_layer, TVector3& egvector){

    TVector3 pixelVector = second_layer - first_layer;
    TVector3 EM_pixelVector = egvector - second_layer;
    return EM_pixelVector.Eta() - pixelVector.Eta();

}

double test::EMmatchingDPhi(TVector3& first_layer, TVector3& second_layer, TVector3& egvector){


    TVector3 pixelVector = second_layer - first_layer;
    TVector3 EM_pixelVector = egvector - second_layer;
    return deltaPhi( EM_pixelVector.Phi(), pixelVector.Phi());


}

int test::Signal_window_check( double upper, double value, double lower, int Ele_Pos){

  if( Ele_Pos == 1 ){ // 1 is Electron
    if( value <= upper && value >= lower){
      return true;
    }
    else
     return false;
  }
  if( Ele_Pos == 2 ){ // 2 is Positron
    if( value >= -upper && value <= -lower){
      return true;
    }
    else
     return false;
  }
 return 0;
}

void test::FillCutFlow(TString cut, float weight){


  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);

  }
  else{
    test::MakeHistograms("cutflow", 6,0.,6.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"MinEtCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"EtaCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"EvtCut");


  }
}

void test::StorePixelHit(int region){

        for(int a=0; a<bRecHitN; a++){
           int Dphi_Ele_pass = 0;
           int Dphi_Pos_pass = 0;
           double Dphi = 0.;
           int el_or_po = 0; // electron = 1, positron = 2, both = 3
           TVector3 current_hit;
           current_hit.SetXYZ( bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a) );
           Dphi = deltaPhi( current_hit.Phi(), EgPhi);

           if( region < 5 ){
             if( bRecHitLayer->at(a) == 1 ){ // First layer
               if( Dphi < L1_Dphi_cut1 && Dphi > L1_Dphi_cut2 ){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L1_Dphi_cut1 && Dphi < -L1_Dphi_cut2 ){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                  layers[1]++;
                  first_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
                  first_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
           }

           if( region == 1 || region == 2 || region == 3){
              if( bRecHitLayer->at(a) == 2 ){ // Second layer
                if( Dphi < L2_Dphi_cut1 && Dphi > L2_Dphi_cut2){
                   Dphi_Ele_pass = 1; el_or_po = 1;
                }
                if( Dphi > -L2_Dphi_cut1 && Dphi < -L2_Dphi_cut2){
                   Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
                }
                if( Dphi_Ele_pass || Dphi_Pos_pass ){
                  layers[2]++;
                  second_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
                  second_layer_hits_Ele_or_Pos.push_back(el_or_po);
                }
              }
           }
           if( region == 1 || region == 2){
             if( bRecHitLayer->at(a) == 3 ){ // Third layer 
               if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[3]++;
                 third_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
                 third_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
           }
           if( region == 1 ){
              if( bRecHitLayer->at(a) == 4 ){ // Fourth layer
                if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
                   Dphi_Ele_pass = 1; el_or_po = 1;
                }
                if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
                   Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
                }
                if( Dphi_Ele_pass || Dphi_Pos_pass ){
                  layers[4]++;
                  fourth_layer_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
                  fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
                }
              }
           }
       }
        for(int a=0; a<fRecHitN; a++){
           int Dphi_Ele_pass = 0;
           int Dphi_Pos_pass = 0;
           double Dphi = 0.;
           int el_or_po = 0; // electron = 1, positron = 2, both = 3
           TVector3 current_hit;
           current_hit.SetXYZ( fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) );
           Dphi = deltaPhi(current_hit.Phi(), EgPhi);

           if( region == 5 ){

             if( fRecHitDisk->at(a) == 1 ){ // second disk
               if( Dphi < L1_Dphi_cut1 && Dphi > L1_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L1_Dphi_cut1 && Dphi < -L1_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[1]++;
                 first_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 first_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }

             if( fRecHitDisk->at(a) == 2 ){ // third disk
               if( Dphi < L2_Dphi_cut1 && Dphi > L2_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L2_Dphi_cut1 && Dphi < -L2_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[2]++;
                 second_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 second_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
             if( fRecHitDisk->at(a) == 3 ){ // fourth disk
               if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[3]++;
                 third_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 third_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
             if( fRecHitDisk->at(a) == 4 ){ // fifth disk
               if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[4]++;
                 fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
           }
	 if( region == 6 ){

             if( fRecHitDisk->at(a) == 2 ){ // second disk
               if( Dphi < L1_Dphi_cut1 && Dphi > L1_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L1_Dphi_cut1 && Dphi < -L1_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[1]++;
                 first_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 first_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }

             if( fRecHitDisk->at(a) == 3 ){ // third disk
               if( Dphi < L2_Dphi_cut1 && Dphi > L2_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L2_Dphi_cut1 && Dphi < -L2_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[2]++;
                 second_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 second_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
             if( fRecHitDisk->at(a) == 4 ){ // fourth disk
               if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[3]++;
                 third_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 third_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }

             if( fRecHitDisk->at(a) == 5 ){ // fifth disk
               if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[4]++;
                 fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }

           }

           if( region == 4 ){

             if( fRecHitDisk->at(a) == 1 ){ // second disk
               if( Dphi < L1_Dphi_cut1 && Dphi > L1_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L1_Dphi_cut1 && Dphi < -L1_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[2]++;
                 second_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 second_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
             if( fRecHitDisk->at(a) == 2 ){ // third disk
               if( Dphi < L2_Dphi_cut1 && Dphi > L2_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L2_Dphi_cut1 && Dphi < -L2_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[3]++;
                 third_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 third_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
             if( fRecHitDisk->at(a) == 3 ){ // fourth disk
               if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[4]++;
                 fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
           }

           if( region == 3 ){
             if( fRecHitDisk->at(a) == 1 ){ // first disk
               if( Dphi < L3_Dphi_cut1 && Dphi > L3_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L3_Dphi_cut1 && Dphi < -L3_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[3]++;
                 third_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 third_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
             if( fRecHitDisk->at(a) == 2 ){ // second disk
               if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[4]++;
                 fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
           }

           if( region == 2 ){
             if( fRecHitDisk->at(a) == 1 ){ // second disk
               if( Dphi < L4_Dphi_cut1 && Dphi > L4_Dphi_cut2){
                  Dphi_Ele_pass = 1; el_or_po = 1;
               }
               if( Dphi > -L4_Dphi_cut1 && Dphi < -L4_Dphi_cut2){
                  Dphi_Pos_pass = 1; el_or_po = el_or_po + 2;
               }
               if( Dphi_Ele_pass || Dphi_Pos_pass ){
                 layers[4]++;
                 fourth_layer_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
                 fourth_layer_hits_Ele_or_Pos.push_back(el_or_po);
               }
             }
           }
        }
}

void test::SetROI(int region){


  float upper_width = 0.055;
  float lower_width = 0.055;

  L1_Dphi_cut1 = ROI_func(region, 0, EgEt);
  L1_Dphi_cut2 = ROI_func(region, 0, EgEt);

  L1_Dphi_cut1 = L1_Dphi_cut1 + upper_width;
  L1_Dphi_cut2 = L1_Dphi_cut2 - lower_width;

  L2_Dphi_cut1 = ROI_func(region, 1, EgEt);
  L2_Dphi_cut2 = ROI_func(region, 1, EgEt);

  L2_Dphi_cut1 = L2_Dphi_cut1 + upper_width;
  L2_Dphi_cut2 = L2_Dphi_cut2 - lower_width;

  L3_Dphi_cut1 = ROI_func(region, 2, EgEt);
  L3_Dphi_cut2 = ROI_func(region, 2, EgEt);

  L3_Dphi_cut1 = L3_Dphi_cut1 + upper_width;
  L3_Dphi_cut2 = L3_Dphi_cut2 - lower_width;

  L4_Dphi_cut1 = ROI_func(region, 3, EgEt);
  L4_Dphi_cut2 = ROI_func(region, 3, EgEt);

  L4_Dphi_cut1 = L4_Dphi_cut1 + upper_width;
  L4_Dphi_cut2 = L4_Dphi_cut2 - lower_width;

}

void test::SetSingalBoundary(int region, double eg_dphi, double eg_deta, double sa_dphi, double sa_deta){

      float EG_pixel_dphi_upper_width = eg_dphi;
      float EG_pixel_dphi_lower_width = eg_dphi;

      // pixel-EG 
      L12_phi_upper =  SW_func2_dphi_v2(region, 0, EgEt);
      L12_phi_bellow = SW_func2_dphi_v2(region, 0, EgEt);
 
      L12_phi_upper = L12_phi_upper   + EG_pixel_dphi_upper_width;
      L12_phi_bellow = L12_phi_bellow - EG_pixel_dphi_lower_width;

      L13_phi_upper =  SW_func2_dphi_v2(region, 1, EgEt);
      L13_phi_bellow = SW_func2_dphi_v2(region, 1, EgEt);

      L13_phi_upper =  L13_phi_upper  + EG_pixel_dphi_upper_width;
      L13_phi_bellow = L13_phi_bellow - EG_pixel_dphi_lower_width;

      L14_phi_upper =  SW_func2_dphi_v2(region, 2, EgEt);
      L14_phi_bellow = SW_func2_dphi_v2(region, 2, EgEt);

      L14_phi_upper =  L14_phi_upper  + EG_pixel_dphi_upper_width;
      L14_phi_bellow = L14_phi_bellow - EG_pixel_dphi_lower_width;

      L23_phi_upper =  SW_func2_dphi_v2(region, 3, EgEt);
      L23_phi_bellow = SW_func2_dphi_v2(region, 3, EgEt);

      L23_phi_upper =  L23_phi_upper  + EG_pixel_dphi_upper_width;
      L23_phi_bellow = L23_phi_bellow - EG_pixel_dphi_lower_width;

      L24_phi_upper =  SW_func2_dphi_v2(region, 4, EgEt);
      L24_phi_bellow = SW_func2_dphi_v2(region, 4, EgEt);

      L24_phi_upper =  L24_phi_upper  + EG_pixel_dphi_upper_width;
      L24_phi_bellow = L24_phi_bellow - EG_pixel_dphi_lower_width;

      L34_phi_upper =  SW_func2_dphi_v2(region, 5, EgEt);
      L34_phi_bellow = SW_func2_dphi_v2(region, 5, EgEt);

      L34_phi_upper = L34_phi_upper   + EG_pixel_dphi_upper_width;
      L34_phi_bellow = L34_phi_bellow - EG_pixel_dphi_lower_width;

      float EG_pixel_deta_upper_width = eg_deta;
      float EG_pixel_deta_lower_width = eg_deta;

      L12_eta_upper =  SW_func2_deta_v2(EgEt);
      L12_eta_bellow = SW_func2_deta_v2(EgEt);

      L12_eta_upper =  L12_eta_upper  + EG_pixel_deta_upper_width;
      L12_eta_bellow = L12_eta_bellow - EG_pixel_deta_lower_width;

      L13_eta_upper =  SW_func2_deta_v2(EgEt);
      L13_eta_bellow = SW_func2_deta_v2(EgEt);

      L13_eta_upper =  L13_eta_upper  + EG_pixel_deta_upper_width;
      L13_eta_bellow = L13_eta_bellow - EG_pixel_deta_lower_width;

      L14_eta_upper =  SW_func2_deta_v2(EgEt);
      L14_eta_bellow = SW_func2_deta_v2(EgEt);

      L14_eta_upper =  L14_eta_upper  + EG_pixel_deta_upper_width;
      L14_eta_bellow = L14_eta_bellow - EG_pixel_deta_lower_width;

      L23_eta_upper =  SW_func2_deta_v2(EgEt);
      L23_eta_bellow = SW_func2_deta_v2(EgEt);

      L23_eta_upper =  L23_eta_upper  + EG_pixel_deta_upper_width;
      L23_eta_bellow = L23_eta_bellow - EG_pixel_deta_lower_width;

      L24_eta_upper =  SW_func2_deta_v2(EgEt);
      L24_eta_bellow = SW_func2_deta_v2(EgEt);

      L24_eta_upper =  L24_eta_upper  + EG_pixel_deta_upper_width;
      L24_eta_bellow = L24_eta_bellow - EG_pixel_deta_lower_width;

      L34_eta_upper =  SW_func2_deta_v2(EgEt);
      L34_eta_bellow = SW_func2_deta_v2(EgEt);

      L34_eta_upper =  L34_eta_upper  + EG_pixel_deta_upper_width;
      L34_eta_bellow = L34_eta_bellow - EG_pixel_deta_lower_width;

      // pixel-pixel 
      float pixel_pixel_dphi_upper_width = sa_dphi;
      float pixel_pixel_dphi_lower_width = sa_dphi;

      L012_DPhi_cut1 = SW_func1_dphi_v2(region, 0, EgEt);
      L012_DPhi_cut2 = SW_func1_dphi_v2(region, 0, EgEt);

      L012_DPhi_cut1 = L012_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L012_DPhi_cut2 = L012_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L013_DPhi_cut1 = SW_func1_dphi_v2(region, 1, EgEt);
      L013_DPhi_cut2 = SW_func1_dphi_v2(region, 1, EgEt);

      L013_DPhi_cut1 = L013_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L013_DPhi_cut2 = L013_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L014_DPhi_cut1 = SW_func1_dphi_v2(region, 2, EgEt);
      L014_DPhi_cut2 = SW_func1_dphi_v2(region, 2, EgEt);

      L014_DPhi_cut1 = L014_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L014_DPhi_cut2 = L014_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L023_DPhi_cut1 = SW_func1_dphi_v2(region, 3, EgEt);
      L023_DPhi_cut2 = SW_func1_dphi_v2(region, 3, EgEt);

      L023_DPhi_cut1 = L023_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L023_DPhi_cut2 = L023_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L024_DPhi_cut1 = SW_func1_dphi_v2(region, 4, EgEt);
      L024_DPhi_cut2 = SW_func1_dphi_v2(region, 4, EgEt);

      L024_DPhi_cut1 = L024_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L024_DPhi_cut2 = L024_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L034_DPhi_cut1 = SW_func1_dphi_v2(region, 5, EgEt);
      L034_DPhi_cut2 = SW_func1_dphi_v2(region, 5, EgEt);

      L034_DPhi_cut1 = L034_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L034_DPhi_cut2 = L034_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L123_DPhi_cut1 = SW_func1_dphi_v2(region, 6, EgEt);
      L123_DPhi_cut2 = SW_func1_dphi_v2(region, 6, EgEt);

      L123_DPhi_cut1 = L123_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L123_DPhi_cut2 = L123_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L124_DPhi_cut1 = SW_func1_dphi_v2(region, 7, EgEt);
      L124_DPhi_cut2 = SW_func1_dphi_v2(region, 7 ,EgEt);

      L124_DPhi_cut1 = L124_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L124_DPhi_cut2 = L124_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L134_DPhi_cut1 = SW_func1_dphi_v2(region, 8, EgEt);
      L134_DPhi_cut2 = SW_func1_dphi_v2(region, 8, EgEt);

      L134_DPhi_cut1 = L134_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L134_DPhi_cut2 = L134_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      L234_DPhi_cut1 = SW_func1_dphi_v2(region, 9, EgEt);
      L234_DPhi_cut2 = SW_func1_dphi_v2(region, 9, EgEt);

      L234_DPhi_cut1 = L234_DPhi_cut1 + pixel_pixel_dphi_upper_width;
      L234_DPhi_cut2 = L234_DPhi_cut2 - pixel_pixel_dphi_lower_width;

      float pixel_pixel_deta_upper_width = sa_deta;
      float pixel_pixel_deta_lower_width = sa_deta;
 
      L123_DEta_cut1 = SW_func1_deta_v2(EgEt);
      L123_DEta_cut2 = SW_func1_deta_v2(EgEt);

      L123_DEta_cut1 = L123_DEta_cut1 + pixel_pixel_deta_upper_width;
      L123_DEta_cut2 = L123_DEta_cut2 - pixel_pixel_deta_lower_width;

      L124_DEta_cut1 = SW_func1_deta_v2(EgEt);
      L124_DEta_cut2 = SW_func1_deta_v2(EgEt);

      L124_DEta_cut1 = L124_DEta_cut1 + pixel_pixel_deta_upper_width;
      L124_DEta_cut2 = L124_DEta_cut2 - pixel_pixel_deta_lower_width;

      L134_DEta_cut1 = SW_func1_deta_v2(EgEt);
      L134_DEta_cut2 = SW_func1_deta_v2(EgEt);

      L134_DEta_cut1 = L134_DEta_cut1 + pixel_pixel_deta_upper_width;
      L134_DEta_cut2 = L134_DEta_cut2 - pixel_pixel_deta_lower_width;

      L234_DEta_cut1 = SW_func1_deta_v2(EgEt);
      L234_DEta_cut2 = SW_func1_deta_v2(EgEt);

      L234_DEta_cut1 = L234_DEta_cut1 + pixel_pixel_deta_upper_width;
      L234_DEta_cut2 = L234_DEta_cut2 - pixel_pixel_deta_lower_width;

}
void test::TriggeringWithout_4thPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit);
dPhi013 = StandaloneDPhi( 0, 1, 3, 0, nthFirstHit, nthThirdHit );
dPhi023 = StandaloneDPhi( 0, 2, 3, 0, nthSecondHit, nthThirdHit);

dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);
dR_1   = sqrt( pow(dPhi_1,2) + pow(dEta_1,2) );

dPhi_2 = EMmatchingDPhi(first_layer_hits[nthFirstHit], third_layer_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(first_layer_hits[nthFirstHit], third_layer_hits[nthThirdHit], emvector);
dR_2   = sqrt( pow(dPhi_2,2) + pow(dEta_2,2) );

dPhi_3 = EMmatchingDPhi(second_layer_hits[nthSecondHit], third_layer_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(second_layer_hits[nthSecondHit], third_layer_hits[nthThirdHit], emvector);
dR_3   = sqrt( pow(dPhi_3,2) + pow(dEta_3,2) );

L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );
L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );
L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );
/*
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele) && Signal_window_check(L12_R_upper, dR_1, L12_R_bellow, Ele) )
   L12_EM_Ele = 1;
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele) && Signal_window_check(L12_R_upper, dR_1, L12_R_bellow, Ele) )
   L12_EM_Pos = 1;

if( Signal_window_check(L13_phi_upper, dPhi_2, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_2, L13_eta_bellow, Ele) && Signal_window_check(L13_R_upper, dR_2, L13_R_bellow, Ele) )
   L13_EM_Ele = 1;
if( Signal_window_check(L13_phi_upper, dPhi_2, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_2, L13_eta_bellow, Ele) && Signal_window_check(L13_R_upper, dR_2, L13_R_bellow, Ele) )
   L13_EM_Pos = 1;

if( Signal_window_check(L23_phi_upper, dPhi_3, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_3, L23_eta_bellow, Ele) && Signal_window_check(L23_R_upper, dR_3, L23_R_bellow, Ele) )
   L23_EM_Ele = 1;
if( Signal_window_check(L23_phi_upper, dPhi_3, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_3, L23_eta_bellow, Ele) && Signal_window_check(L23_R_upper, dR_3, L23_R_bellow, Ele) )
   L23_EM_Pos = 1;

if( Signal_window_check(L123_DPhi_cut1, dPhi, L123_DPhi_cut2, Ele ) && Signal_window_check(L123_DEta_cut1, dEta, L123_DEta_cut2, Ele ) && Signal_window_check(L123_DR_cut1, dR, L123_DR_cut2, Ele ) )
   L123_pass_Ele = 1;
if( Signal_window_check(L123_DPhi_cut1, dPhi, L123_DPhi_cut2, Pos ) && Signal_window_check(L123_DEta_cut1, dEta, L123_DEta_cut2, Ele ) && Signal_window_check(L123_DR_cut1, dR, L123_DR_cut2, Ele ) ) 
   L123_pass_Pos = 1;
*/
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele) )
   L12_EM_Ele = 1;
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele) )
   L12_EM_Pos = 1;
   
if( Signal_window_check(L13_phi_upper, dPhi_2, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_2, L13_eta_bellow, Ele) )
   L13_EM_Ele = 1;
if( Signal_window_check(L13_phi_upper, dPhi_2, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_2, L13_eta_bellow, Ele) )
   L13_EM_Pos = 1;
   
if( Signal_window_check(L23_phi_upper, dPhi_3, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_3, L23_eta_bellow, Ele) )
   L23_EM_Ele = 1;
if( Signal_window_check(L23_phi_upper, dPhi_3, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_3, L23_eta_bellow, Ele) )
   L23_EM_Pos = 1;
   
if( Signal_window_check(L123_DPhi_cut1, dPhi, L123_DPhi_cut2, Ele ) && Signal_window_check(L123_DEta_cut1, dEta, L123_DEta_cut2, Ele ) )
   L123_pass_Ele = 1;
if( Signal_window_check(L123_DPhi_cut1, dPhi, L123_DPhi_cut2, Pos ) && Signal_window_check(L123_DEta_cut1, dEta, L123_DEta_cut2, Ele ) )
   L123_pass_Pos = 1;

}

void test::TriggeringWithout_3rdPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit);
dPhi014 = StandaloneDPhi( 0, 1, 4, 0, nthFirstHit, nthThirdHit );
dPhi024 = StandaloneDPhi( 0, 2, 4, 0, nthSecondHit, nthThirdHit);

dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);
dR_1   = sqrt( pow(dPhi_1,2) + pow(dEta_1,2) );

dPhi_2 = EMmatchingDPhi(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dR_2   = sqrt( pow(dPhi_2,2) + pow(dEta_2,2) );

dPhi_3 = EMmatchingDPhi(second_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(second_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dR_3   = sqrt( pow(dPhi_3,2) + pow(dEta_3,2) );

L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );
L014_pass_Ele = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Ele );
L014_pass_Pos = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Pos );
L024_pass_Ele = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Ele );
L024_pass_Pos = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Pos );

/*
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele) && Signal_window_check(L12_R_upper, dR_1, L12_R_bellow, Ele) ) 
   L12_EM_Ele = 1;
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele) && Signal_window_check(L12_R_upper, dR_1, L12_R_bellow, Ele) )
   L12_EM_Pos = 1;

if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Ele) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele) && Signal_window_check(L14_R_upper, dR_2, L14_R_bellow, Ele) )
   L14_EM_Ele = 1;
if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Pos) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele) && Signal_window_check(L14_R_upper, dR_2, L14_R_bellow, Ele) )
   L14_EM_Pos = 1;

if( Signal_window_check(L24_phi_upper, dPhi_3, L24_phi_bellow, Ele) && Signal_window_check(L24_eta_upper, dEta_3, L24_eta_bellow, Ele) && Signal_window_check(L24_R_upper, dR_3, L24_R_bellow, Ele) ) 
   L24_EM_Ele = 1;
if( Signal_window_check(L24_phi_upper, dPhi_3, L24_phi_bellow, Pos) && Signal_window_check(L24_eta_upper, dEta_3, L24_eta_bellow, Ele) && Signal_window_check(L24_R_upper, dR_3, L24_R_bellow, Ele) )
   L24_EM_Pos = 1;

if( Signal_window_check(L124_DPhi_cut1, dPhi, L124_DPhi_cut2, Ele ) && Signal_window_check(L124_DEta_cut1, dEta, L124_DEta_cut2, Ele ) && Signal_window_check(L124_DR_cut1, dR, L124_DR_cut2, Ele ) )
   L124_pass_Ele = 1;
if( Signal_window_check(L124_DPhi_cut1, dPhi, L124_DPhi_cut2, Pos ) && Signal_window_check(L124_DEta_cut1, dEta, L124_DEta_cut2, Ele ) && Signal_window_check(L124_DR_cut1, dR, L124_DR_cut2, Ele ) )
   L124_pass_Pos = 1;
*/

if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
   L12_EM_Ele = 1;
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
   L12_EM_Pos = 1;
   
if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Ele) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele)  )
   L14_EM_Ele = 1;
if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Pos) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele)  )
   L14_EM_Pos = 1;
   
if( Signal_window_check(L24_phi_upper, dPhi_3, L24_phi_bellow, Ele) && Signal_window_check(L24_eta_upper, dEta_3, L24_eta_bellow, Ele)  )
   L24_EM_Ele = 1;
if( Signal_window_check(L24_phi_upper, dPhi_3, L24_phi_bellow, Pos) && Signal_window_check(L24_eta_upper, dEta_3, L24_eta_bellow, Ele)  )
   L24_EM_Pos = 1;
   
if( Signal_window_check(L124_DPhi_cut1, dPhi, L124_DPhi_cut2, Ele ) && Signal_window_check(L124_DEta_cut1, dEta, L124_DEta_cut2, Ele )  )
   L124_pass_Ele = 1;
if( Signal_window_check(L124_DPhi_cut1, dPhi, L124_DPhi_cut2, Pos ) && Signal_window_check(L124_DEta_cut1, dEta, L124_DEta_cut2, Ele )  )
   L124_pass_Pos = 1;
}

void test::TriggeringWithout_2ndPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

dPhi013 = StandaloneDPhi( 0, 1, 3, 0, nthFirstHit, nthSecondHit);
dPhi014 = StandaloneDPhi( 0, 1, 4, 0, nthFirstHit, nthThirdHit );
dPhi034 = StandaloneDPhi( 0, 3, 4, 0, nthSecondHit, nthThirdHit);


dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
dR_1   = sqrt( pow(dPhi_1,2) + pow(dEta_1,2) );

dPhi_2 = EMmatchingDPhi(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dR_2   = sqrt( pow(dPhi_2,2) + pow(dEta_2,2) );

dPhi_3 = EMmatchingDPhi(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dR_3   = sqrt( pow(dPhi_3,2) + pow(dEta_3,2) );

L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );
L014_pass_Ele = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Ele );
L014_pass_Pos = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Pos );
L034_pass_Ele = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Ele );
L034_pass_Pos = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Pos );
/*
if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele) && Signal_window_check(L13_R_upper, dR_1, L13_R_bellow, Ele) ) 
   L13_EM_Ele = 1;
if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele) && Signal_window_check(L13_R_upper, dR_1, L13_R_bellow, Ele) )
   L13_EM_Pos = 1;

if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Ele) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele) && Signal_window_check(L14_R_upper, dR_2, L14_R_bellow, Ele) )
   L14_EM_Ele = 1;
if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Pos) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele) && Signal_window_check(L14_R_upper, dR_2, L14_R_bellow, Ele) )
   L14_EM_Pos = 1;

if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Ele) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele) && Signal_window_check(L34_R_upper, dR_3, L34_R_bellow, Ele) )
   L34_EM_Ele = 1;
if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Pos) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele) && Signal_window_check(L34_R_upper, dR_3, L34_R_bellow, Ele) ) 
   L34_EM_Pos = 1;

if( Signal_window_check(L134_DPhi_cut1, dPhi, L134_DPhi_cut2, Ele ) && Signal_window_check(L134_DEta_cut1, dEta, L134_DEta_cut2, Ele ) && Signal_window_check(L134_DR_cut1, dR, L134_DR_cut2, Ele ) )
  L134_pass_Ele = 1;
if( Signal_window_check(L134_DPhi_cut1, dPhi, L134_DPhi_cut2, Pos ) && Signal_window_check(L134_DEta_cut1, dEta, L134_DEta_cut2, Ele ) && Signal_window_check(L134_DR_cut1, dR, L134_DR_cut2, Ele ) ) 
  L134_pass_Pos = 1;
*/

if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele)  )
   L13_EM_Ele = 1;
if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele)  )
   L13_EM_Pos = 1;
   
if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Ele) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele)  )
   L14_EM_Ele = 1;
if( Signal_window_check(L14_phi_upper, dPhi_2, L14_phi_bellow, Pos) && Signal_window_check(L14_eta_upper, dEta_2, L14_eta_bellow, Ele)  )
   L14_EM_Pos = 1;
   
if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Ele) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele)  )
   L34_EM_Ele = 1;
if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Pos) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele)  )
   L34_EM_Pos = 1;
   
if( Signal_window_check(L134_DPhi_cut1, dPhi, L134_DPhi_cut2, Ele ) && Signal_window_check(L134_DEta_cut1, dEta, L134_DEta_cut2, Ele )  )
  L134_pass_Ele = 1;
if( Signal_window_check(L134_DPhi_cut1, dPhi, L134_DPhi_cut2, Pos ) && Signal_window_check(L134_DEta_cut1, dEta, L134_DEta_cut2, Ele )  )
  L134_pass_Pos = 1;
}

void test::TriggeringWithout_1stPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

dPhi023 = StandaloneDPhi( 0, 2, 3, 0, nthFirstHit, nthSecondHit);
dPhi024 = StandaloneDPhi( 0, 2, 4, 0, nthFirstHit, nthThirdHit );
dPhi034 = StandaloneDPhi( 0, 3, 4, 0, nthSecondHit, nthThirdHit);

dPhi_1 = EMmatchingDPhi(second_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(second_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
dR_1   = sqrt( pow(dPhi_1,2) + pow(dEta_1,2) );

dPhi_2 = EMmatchingDPhi(second_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(second_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dR_2   = sqrt( pow(dPhi_2,2) + pow(dEta_2,2) );

dPhi_3 = EMmatchingDPhi(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dR_3   = sqrt( pow(dPhi_3,2) + pow(dEta_3,2) );

L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );
L024_pass_Ele = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Ele );
L024_pass_Pos = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Pos );
L034_pass_Ele = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Ele );
L034_pass_Pos = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Pos );
/*
if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele) && Signal_window_check(L23_R_upper, dR_1, L23_R_bellow, Ele) )
   L23_EM_Ele = 1;
if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele) && Signal_window_check(L23_R_upper, dR_1, L23_R_bellow, Ele) )
   L23_EM_Pos = 1;

if( Signal_window_check(L24_phi_upper, dPhi_2, L24_phi_bellow, Ele) && Signal_window_check(L24_eta_upper, dEta_2, L24_eta_bellow, Ele) && Signal_window_check(L24_R_upper, dR_2, L24_R_bellow, Ele) )
   L24_EM_Ele = 1;
if( Signal_window_check(L24_phi_upper, dPhi_2, L24_phi_bellow, Pos) && Signal_window_check(L24_eta_upper, dEta_2, L24_eta_bellow, Ele) && Signal_window_check(L24_R_upper, dR_2, L24_R_bellow, Ele) )
   L24_EM_Pos = 1;

if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Ele) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele) && Signal_window_check(L34_R_upper, dR_3, L34_R_bellow, Ele) )
   L34_EM_Ele = 1;
if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Pos) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele) && Signal_window_check(L34_R_upper, dR_3, L34_R_bellow, Ele) )
   L34_EM_Pos = 1;

if( Signal_window_check(L234_DPhi_cut1, dPhi, L234_DPhi_cut2, Ele ) && Signal_window_check(L234_DEta_cut1, dEta, L234_DEta_cut2, Ele ) && Signal_window_check(L234_DR_cut1, dR, L234_DR_cut2, Ele ) )
   L234_pass_Ele = 1;
if( Signal_window_check(L234_DPhi_cut1, dPhi, L234_DPhi_cut2, Pos ) && Signal_window_check(L234_DEta_cut1, dEta, L234_DEta_cut2, Ele ) && Signal_window_check(L234_DR_cut1, dR, L234_DR_cut2, Ele ) )
   L234_pass_Pos = 1;
*/

if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele)  )
   L23_EM_Ele = 1;
if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele)  )
   L23_EM_Pos = 1;

if( Signal_window_check(L24_phi_upper, dPhi_2, L24_phi_bellow, Ele) && Signal_window_check(L24_eta_upper, dEta_2, L24_eta_bellow, Ele)  )
   L24_EM_Ele = 1;
if( Signal_window_check(L24_phi_upper, dPhi_2, L24_phi_bellow, Pos) && Signal_window_check(L24_eta_upper, dEta_2, L24_eta_bellow, Ele)  )
   L24_EM_Pos = 1;

if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Ele) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele)  )
   L34_EM_Ele = 1;
if( Signal_window_check(L34_phi_upper, dPhi_3, L34_phi_bellow, Pos) && Signal_window_check(L34_eta_upper, dEta_3, L34_eta_bellow, Ele)  )
   L34_EM_Pos = 1;

if( Signal_window_check(L234_DPhi_cut1, dPhi, L234_DPhi_cut2, Ele ) && Signal_window_check(L234_DEta_cut1, dEta, L234_DEta_cut2, Ele )  )
   L234_pass_Ele = 1;
if( Signal_window_check(L234_DPhi_cut1, dPhi, L234_DPhi_cut2, Pos ) && Signal_window_check(L234_DEta_cut1, dEta, L234_DEta_cut2, Ele )  )
   L234_pass_Pos = 1;
}
void test::TriggeringWith_1st2ndPixel(int nthFirstHit, int nthSecondHit){
dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );


if( L012_pass_Ele &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (second_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || second_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L012_pass_Pos &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (second_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || second_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;

}

void test::TriggeringWith_1st3rdPixel(int nthFirstHit, int nthSecondHit){
dPhi013 = StandaloneDPhi( 0, 1, 3, 0, nthFirstHit, nthSecondHit );
L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );

 if( L013_pass_Ele &&
     (first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
     (third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
 if( L013_pass_Pos &&
     (first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
     (third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;

}
void test::TriggeringWith_2nd3rdPixel(int nthFirstHit, int nthSecondHit){
dPhi023 = StandaloneDPhi( 0, 2, 3, 0, nthFirstHit, nthSecondHit );
L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );

 if( L023_pass_Ele &&
     (second_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || second_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
     (third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
 if( L023_pass_Pos &&
     (second_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || second_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
     (third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_1st2ndPixel_v2(int nthFirstHit, int nthSecondHit){
dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );

dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);

if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
   L12_EM_Ele = 1;
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
   L12_EM_Pos = 1;

if( L012_pass_Ele && L12_EM_Ele &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (second_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || second_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L012_pass_Pos && L12_EM_Pos &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (second_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || second_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_1st3rdPixel_v2(int nthFirstHit, int nthSecondHit){
dPhi013 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );

dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);

if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele)  )
   L13_EM_Ele = 1;
if( Signal_window_check(L13_phi_upper, dPhi_1, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_1, L13_eta_bellow, Ele)  )
   L13_EM_Pos = 1;

if( L013_pass_Ele && L13_EM_Ele &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L013_pass_Pos && L13_EM_Pos &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_1st4thPixel_v2(int nthFirstHit, int nthSecondHit){
dPhi014 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L014_pass_Ele = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Ele );
L014_pass_Pos = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Pos );

dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], fourth_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], fourth_layer_hits[nthSecondHit], emvector);

if( Signal_window_check(L14_phi_upper, dPhi_1, L14_phi_bellow, Ele) && Signal_window_check(L14_eta_upper, dEta_1, L14_eta_bellow, Ele)  )
   L14_EM_Ele = 1;
if( Signal_window_check(L14_phi_upper, dPhi_1, L14_phi_bellow, Pos) && Signal_window_check(L14_eta_upper, dEta_1, L14_eta_bellow, Ele)  )
   L14_EM_Pos = 1;

if( L014_pass_Ele && L14_EM_Ele &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (fourth_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || fourth_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L014_pass_Pos && L14_EM_Pos &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (fourth_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || fourth_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_2nd3rdPixel_v2(int nthFirstHit, int nthSecondHit){
dPhi023 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );

dPhi_1 = EMmatchingDPhi(second_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(second_layer_hits[nthFirstHit], third_layer_hits[nthSecondHit], emvector);

if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele)  )
   L23_EM_Ele = 1;
if( Signal_window_check(L23_phi_upper, dPhi_1, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_1, L23_eta_bellow, Ele)  )
   L23_EM_Pos = 1;

if( L023_pass_Ele && L23_EM_Ele &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L023_pass_Pos && L23_EM_Pos &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_2nd4thPixel_v2(int nthFirstHit, int nthSecondHit){
dPhi024 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L024_pass_Ele = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Ele );
L024_pass_Pos = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Pos );

dPhi_1 = EMmatchingDPhi(second_layer_hits[nthFirstHit], fourth_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(second_layer_hits[nthFirstHit], fourth_layer_hits[nthSecondHit], emvector);

if( Signal_window_check(L24_phi_upper, dPhi_1, L24_phi_bellow, Ele) && Signal_window_check(L24_eta_upper, dEta_1, L24_eta_bellow, Ele)  )
   L24_EM_Ele = 1;
if( Signal_window_check(L24_phi_upper, dPhi_1, L24_phi_bellow, Pos) && Signal_window_check(L24_eta_upper, dEta_1, L24_eta_bellow, Ele)  )
   L24_EM_Pos = 1;

if( L024_pass_Ele && L24_EM_Ele &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L024_pass_Pos && L24_EM_Pos &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::TriggeringWith_3rd4thPixel_v2(int nthFirstHit, int nthSecondHit){
dPhi034 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit );
L034_pass_Ele = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Ele );
L034_pass_Pos = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Pos );

dPhi_1 = EMmatchingDPhi(third_layer_hits[nthFirstHit], fourth_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(third_layer_hits[nthFirstHit], fourth_layer_hits[nthSecondHit], emvector);

if( Signal_window_check(L34_phi_upper, dPhi_1, L34_phi_bellow, Ele) && Signal_window_check(L34_eta_upper, dEta_1, L34_eta_bellow, Ele)  )
   L34_EM_Ele = 1;
if( Signal_window_check(L34_phi_upper, dPhi_1, L34_phi_bellow, Pos) && Signal_window_check(L34_eta_upper, dEta_1, L34_eta_bellow, Ele)  )
   L34_EM_Pos = 1;

if( L034_pass_Ele && L34_EM_Ele &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 1 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_layer_hits_Ele_or_Pos[nthSecondHit] == 1 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Ele = 1;
if( L034_pass_Pos && L34_EM_Pos &&
    (first_layer_hits_Ele_or_Pos[nthFirstHit] == 2 || first_layer_hits_Ele_or_Pos[nthFirstHit] ==3) &&
    (third_layer_hits_Ele_or_Pos[nthSecondHit] == 2 || third_layer_hits_Ele_or_Pos[nthSecondHit] ==3)) _pass_Pos = 1;
}

void test::StorePixelHitForIso(int region, float recoPV){

   // Store pixel clusters in vectors based on eta region
   // 0.0 < eta < 0.8 : region = 1  L1234
   // 0.8 < eta < 1.4 : region = 2  L123D1
   // 1.4 < eta < 1.7 : region = 3  L12D12
   // 1.7 < eta < 2.1 : region = 4  L1D123
   // 2.1 < eta < 2.7 : region = 5  D1234
   // 2.7 < eta < 3.0 : region = 6  D2345

        for(int a=0; a<bRecHitN; a++){
           TVector3 current_hit;
           current_hit.SetXYZ( bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a) - recoPV );
           Float_t Dphi = deltaPhi( current_hit.Phi(), EgPhi);
           Float_t Deta = current_hit.Eta() - EgEta;
           Float_t DR = sqrt(pow(Dphi,2)+pow(Deta,2));
           
           if( fabs(Dphi) > 0.3 ) continue;
           if( DR > 0.3 ) continue;

           if( region <= 4 ){
               if( bRecHitLayer->at(a) == 1 ){ // First layer
                   first_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
               }
           }
           
           if( region >= 1 && region <= 3 ) {
               if( bRecHitLayer->at(a) == 2 ){ // Second layer
                   second_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
               }
           }
           
           if( region >= 1 && region <= 2 ){
               if( bRecHitLayer->at(a) == 3 ){ // Third layer 
                   third_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
               }
           }
           
           if( region == 1 ){
               if( bRecHitLayer->at(a) == 4 ){ // Fourth layer
                   fourth_hits.push_back( TVector3(bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a)));
               }
           }
       } // Pixel layer loop
        
        for(int a=0; a<fRecHitN; a++){
           TVector3 current_hit;
           current_hit.SetXYZ( fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) - recoPV );
           Float_t Dphi = deltaPhi(current_hit.Phi(), EgPhi);
           if( fabs(Dphi) > 0.3 ) continue;
           
           Float_t Deta = current_hit.Eta() - EgEta;
           Float_t DR = sqrt(pow(Dphi,2)+pow(Deta,2));
           if( DR > 0.3 ) continue;

           if( region == 6 ){
               if( fRecHitDisk->at(a) == 2 ){ // second disk
                   first_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }
               if( fRecHitDisk->at(a) == 3 ){ // third disk
                   second_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }
               if( fRecHitDisk->at(a) == 4 ){ // fourth disk
                   third_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }
               if( fRecHitDisk->at(a) == 5 ){ // fifth disk
                   fourth_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }
           }

           if( region == 5 ){
               if( fRecHitDisk->at(a) == 1 ){ // first disk
                   first_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }
               if( fRecHitDisk->at(a) == 2 ){ // second disk
                   second_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }
               if( fRecHitDisk->at(a) == 3 ){ // third disk
                   third_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }
               if( fRecHitDisk->at(a) == 4 ){ // fourth disk
                   fourth_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }
           }
           if( region == 4 ){
               if( fRecHitDisk->at(a) == 1 ){ // first disk
                   second_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }

               if( fRecHitDisk->at(a) == 2 ){ // second disk
                   third_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }

               if( fRecHitDisk->at(a) == 3 ){ // third disk
                   fourth_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }
           }

           if( region == 3 ){
               if( fRecHitDisk->at(a) == 1 ){ // first disk
                   third_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }
               if( fRecHitDisk->at(a) == 2 ){ // second disk
                   fourth_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }
           }

           if( region == 2 ){
               if( fRecHitDisk->at(a) == 1 ){ // first disk
                   fourth_hits.push_back( TVector3(fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a)));
               }
           }

        } // pixel disk loop
}

void test::IsoWith_1st2nd3rd(int region, float recoPV)
{
    Float_t dzCut1 = 0., dzCut2 = 0., dzCut3 = 0.;
    Float_t dEtappCut1 = 0., dEtappCut2 = 0., dEtappCut3 = 0.;
    Float_t dEtaPvCut1 = 0., dEtaPvCut2 = 0.;
    Float_t ddPhi_left = 0., ddPhi_right = 0.;

    if( region == 1 ) {
        dzCut1 = 0.0305; dzCut2 = 0.0305; dzCut3 = 0.0325; 
        dEtappCut1 = 0.00225; dEtappCut2 = 0.00363; dEtappCut3 = 0.00181; 
        dEtaPvCut1 = 0.00103; dEtaPvCut2 = 0.00217;
        ddPhi_left = -0.0025; ddPhi_right = 0.0025;
    }
    if( region == 2 ) {
        dzCut1 = 0.0205; dzCut2 = 0.0185; dzCut3 = 0.0205; 
        dEtappCut1 = 0.00136; dEtappCut2 = 0.00226; dEtappCut3 = 0.00110; 
        dEtaPvCut1 = 0.00070; dEtaPvCut2 = 0.00136;
        ddPhi_left = -0.0025; ddPhi_right = 0.002;
    }
    if( region == 3 ) {
        dzCut1 = 0.0290; dzCut2 = 0.0280; dzCut3 = 0.0570; 
        dEtappCut1 = 0.00138; dEtappCut2 = 0.00230; dEtappCut3 = 0.00112; 
        dEtaPvCut1 = 0.00070; dEtaPvCut2 = 0.00092;
        ddPhi_left = -0.0035; ddPhi_right = 0.002;
    }
    if( region == 4 ) {
        dzCut1 = 0.0360; dzCut2 = 0.0280; dzCut3 = 0.0940; 
        dEtappCut1 = 0.00280; dEtappCut2 = 0.00560; dEtappCut3 = 0.00374; 
        dEtaPvCut1 = 0.00102; dEtaPvCut2 = 0.00092;
        ddPhi_left = -0.002; ddPhi_right = 0.0025;
    }
    if( region == 5 ) {
        dzCut1 = 0.184; dzCut2 = 0.163; dzCut3 = 0.189; 
        dEtappCut1 = 0.00702; dEtappCut2 = 0.01; dEtappCut3 = 0.00530; 
        dEtaPvCut1 = 0.00138; dEtaPvCut2 = 0.00168;
        ddPhi_left = -0.003; ddPhi_right = 0.003;
    }
    if( region == 6 ) {
        dzCut1 = 0.2; dzCut2 = 0.2; dzCut3 = 0.2; 
        dEtappCut1 = 0.00898; dEtappCut2 = 0.01; dEtappCut3 = 0.00720; 
        dEtaPvCut1 = 0.00176; dEtaPvCut2 = 0.00256;
        ddPhi_left = -0.004; ddPhi_right = 0.003;
    }
    
    for(std::vector<TVector3>::iterator a = first_hits.begin(); a != first_hits.end(); ++a)
    {
        TVector3 PVL1;
        PVL1.SetXYZ( (*a).X(), (*a).Y(), (*a).Z() - recoPV );
        Float_t PhiPVL1 = PVL1.Phi();
        Float_t EtaPVL1 = PVL1.Eta();
        
        Float_t R1 = sqrt(pow((*a).X(),2)+pow((*a).Y(),2));
        Float_t Z1 = (*a).Z();
        
        for(std::vector<TVector3>::iterator b = second_hits.begin(); b != second_hits.end(); ++b)
        {
            TVector3 L1L2;
            L1L2.SetXYZ( (*b).X()-(*a).X(), (*b).Y()-(*a).Y(), (*b).Z()-(*a).Z() );
            Float_t PhiL1L2 = L1L2.Phi();
            Float_t EtaL1L2 = L1L2.Eta();
            
            TVector3 PVL2;
            PVL2.SetXYZ( (*b).X(), (*b).Y(), (*b).Z() - recoPV );
            Float_t EtaPVL2 = PVL2.Eta();
            
            Float_t R2 = sqrt(pow((*b).X(),2)+pow((*b).Y(),2));
            Float_t Z2 = (*b).Z();
            
            // DeltaZ 1st Combination ( ex) L1-L2 )
            Float_t vertex1 = (R2*Z1-R1*Z2)/(R2-R1);
            Float_t deltaZ1 = fabs(vertex1-recoPV);

            if( deltaZ1 < dzCut1 )
            {
                for(std::vector<TVector3>::iterator c = third_hits.begin(); c != third_hits.end(); ++c)
                {
                    TVector3 L1L3;
                    L1L3.SetXYZ( (*c).X()-(*a).X(), (*c).Y()-(*a).Y(), (*c).Z()-(*a).Z() );
                    Float_t EtaL1L3 = L1L3.Eta();

                    TVector3 L2L3;
                    L2L3.SetXYZ( (*c).X()-(*b).X(), (*c).Y()-(*b).Y(), (*c).Z()-(*b).Z() );
                    Float_t PhiL2L3 = L2L3.Phi();
                    Float_t EtaL2L3 = L2L3.Eta();

                    TVector3 PVL3;
                    PVL3.SetXYZ( (*c).X(), (*c).Y(), (*c).Z() - recoPV );
                    Float_t EtaPVL3 = PVL3.Eta();

                    Float_t R3 = sqrt(pow((*c).X(),2)+pow((*c).Y(),2));
                    Float_t Z3 = (*c).Z();
                    
                    // DeltaZ 2nd Combination ( ex) L1-L3 )
                    Float_t vertex2 = (R3*Z1-R1*Z3)/(R3-R1);
                    Float_t deltaZ2 = fabs(vertex2-recoPV);

                    // DeltaZ 3rd Combination ( ex) L2-L3 )
                    Float_t vertex3 = (R3*Z2-R2*Z3)/(R3-R2);
                    Float_t deltaZ3 = fabs(vertex3-recoPV);
                    
                    // DeltaEta 1st Combination ( ex) L12-L13 )
                    Float_t dEtapp1 = fabs(EtaL1L2 - EtaL1L3);

                    // DeltaEta 2nd Combination ( ex) L12-L23 )
                    Float_t dEtapp2 = fabs(EtaL1L2 - EtaL2L3);

                    // DeltaEta 3rd Combination ( ex) L13-L23 )
                    Float_t dEtapp3 = fabs(EtaL1L3 - EtaL2L3);

                    // DeltaEta-PV 1st Combination ( ex) PVL3 - PVL2 )
                    Float_t dEtaPv1 = fabs(EtaPVL3 - EtaPVL2);

                    // DeltaEta-PV 2nd Combination ( ex) PVL3 - PVL1 )
                    Float_t dEtaPv2 = fabs(EtaPVL3 - EtaPVL1);
                    
                    // DeltaPhi difference
                    Float_t dphi1 = deltaPhi(PhiPVL1, PhiL1L2);
                    Float_t dphi2 = deltaPhi(PhiL1L2, PhiL2L3);
                    Float_t ddphi = dphi1 - dphi2;
                    
                    if( deltaZ2 < dzCut2 && deltaZ3 < dzCut3 )
                        if( dEtapp1 < dEtappCut1 && dEtapp2 < dEtappCut2 && dEtapp3 < dEtappCut3 )
                            if( dEtaPv1 < dEtaPvCut1 && dEtaPv2 < dEtaPvCut2 )
                                if( (dphi1 > 0 && dphi2 > 0) || (dphi1 < 0 && dphi2 < 0) )
                                    if( ddphi > ddPhi_left && ddphi < ddPhi_right )
                                        L123.push_back( track((*c).X(), (*b).X(), (*a).X(),
                                                              (*c).Y(), (*b).Y(), (*a).Y(),
                                                              (*c).Z(), (*b).Z(), (*a).Z(),
                                                              1) );
                } // third hit loop
            }
        } // second hit loop
    } // first hit loop
}

void test::IsoWith_1st2nd4th(int region, float recoPV)
{
    Float_t dzCut1 = 0., dzCut2 = 0., dzCut3 = 0.;
    Float_t dEtappCut1 = 0., dEtappCut2 = 0., dEtappCut3 = 0.;
    Float_t dEtaPvCut1 = 0., dEtaPvCut2 = 0.;
    Float_t ddPhi_left = 0., ddPhi_right = 0.;

    if( region == 1 ) {
        dzCut1 = 0.0315; dzCut2 = 0.0315; dzCut3 = 0.0315; 
        dEtappCut1 = 0.00223; dEtappCut2 = 0.00313; dEtappCut3 = 0.00091; 
        dEtaPvCut1 = 0.00101; dEtaPvCut2 = 0.00185;
        ddPhi_left = -0.004; ddPhi_right = 0.0015;
    }
    if( region == 2 ) {
        dzCut1 = 0.0185; dzCut2 = 0.0185; dzCut3 = 0.0185; 
        dEtappCut1 = 0.00110; dEtappCut2 = 0.00184; dEtappCut3 = 0.00108; 
        dEtaPvCut1 = 0.00068; dEtaPvCut2 = 0.00110;
        ddPhi_left = -0.004; ddPhi_right = 0.002;
    }
    if( region == 3 ) {
        dzCut1 = 0.0190; dzCut2 = 0.0180; dzCut3 = 0.0190; 
        dEtappCut1 = 0.00112; dEtappCut2 = 0.00230; dEtappCut3 = 0.00074; 
        dEtaPvCut1 = 0.00072; dEtaPvCut2 = 0.00092;
        ddPhi_left = -0.004; ddPhi_right = 0.002;
    }
    if( region == 4 ) {
        dzCut1 = 0.0360; dzCut2 = 0.0270; dzCut3 = 0.0560; 
        dEtappCut1 = 0.00232; dEtappCut2 = 0.00326; dEtappCut3 = 0.00138; 
        dEtaPvCut1 = 0.00088; dEtaPvCut2 = 0.00092;
        ddPhi_left = -0.0025; ddPhi_right = 0.0025;
    }
    if( region == 5 ) {
        dzCut1 = 0.190; dzCut2 = 0.111; dzCut3 = 0.138; 
        dEtappCut1 = 0.00698; dEtappCut2 = 0.00900; dEtappCut3 = 0.00248; 
        dEtaPvCut1 = 0.00130; dEtaPvCut2 = 0.00152;
        ddPhi_left = -0.003; ddPhi_right = 0.003;
    }
    if( region == 6 ) {
        dzCut1 = 0.2; dzCut2 = 0.149; dzCut3 = 0.189; 
        dEtappCut1 = 0.00894; dEtappCut2 = 0.01; dEtappCut3 = 0.00276; 
        dEtaPvCut1 = 0.00154; dEtaPvCut2 = 0.00182;
        ddPhi_left = -0.004; ddPhi_right = 0.004;
    }
    
    for(std::vector<TVector3>::iterator a = first_hits.begin(); a != first_hits.end(); ++a)
    {
        TVector3 PVL1;
        PVL1.SetXYZ( (*a).X(), (*a).Y(), (*a).Z() - recoPV );
        Float_t PhiPVL1 = PVL1.Phi();
        Float_t EtaPVL1 = PVL1.Eta();
        
        Float_t R1 = sqrt(pow((*a).X(),2)+pow((*a).Y(),2));
        Float_t Z1 = (*a).Z();
        
        for(std::vector<TVector3>::iterator b = second_hits.begin(); b != second_hits.end(); ++b)
        {
            TVector3 L1L2;
            L1L2.SetXYZ( (*b).X()-(*a).X(), (*b).Y()-(*a).Y(), (*b).Z()-(*a).Z() );
            Float_t PhiL1L2 = L1L2.Phi();
            Float_t EtaL1L2 = L1L2.Eta();
            
            TVector3 PVL2;
            PVL2.SetXYZ( (*b).X(), (*b).Y(), (*b).Z() - recoPV );
            Float_t EtaPVL2 = PVL2.Eta();
            
            Float_t R2 = sqrt(pow((*b).X(),2)+pow((*b).Y(),2));
            Float_t Z2 = (*b).Z();
            
            // DeltaZ 1st Combination ( ex) L1-L2 )
            Float_t vertex1 = (R2*Z1-R1*Z2)/(R2-R1);
            Float_t deltaZ1 = fabs(vertex1-recoPV);

            if( deltaZ1 < dzCut1 )
            {
                for(std::vector<TVector3>::iterator c = fourth_hits.begin(); c != fourth_hits.end(); ++c)
                {
                    TVector3 L1L4;
                    L1L4.SetXYZ( (*c).X()-(*a).X(), (*c).Y()-(*a).Y(), (*c).Z()-(*a).Z() );
                    Float_t EtaL1L4 = L1L4.Eta();

                    TVector3 L2L4;
                    L2L4.SetXYZ( (*c).X()-(*b).X(), (*c).Y()-(*b).Y(), (*c).Z()-(*b).Z() );
                    Float_t PhiL2L4 = L2L4.Phi();
                    Float_t EtaL2L4 = L2L4.Eta();

                    TVector3 PVL4;
                    PVL4.SetXYZ( (*c).X(), (*c).Y(), (*c).Z() - recoPV );
                    Float_t EtaPVL4 = PVL4.Eta();

                    Float_t R4 = sqrt(pow((*c).X(),2)+pow((*c).Y(),2));
                    Float_t Z4 = (*c).Z();
                    
                    // DeltaZ 2nd Combination ( ex) L1-L4 )
                    Float_t vertex2 = (R4*Z1-R1*Z4)/(R4-R1);
                    Float_t deltaZ2 = fabs(vertex2-recoPV);

                    // DeltaZ 3rd Combination ( ex) L2-L4 )
                    Float_t vertex3 = (R4*Z2-R2*Z4)/(R4-R2);
                    Float_t deltaZ3 = fabs(vertex3-recoPV);
                    
                    // DeltaEta 1st Combination ( ex) L12-L14 )
                    Float_t dEtapp1 = fabs(EtaL1L2 - EtaL1L4);

                    // DeltaEta 2nd Combination ( ex) L12-L24 )
                    Float_t dEtapp2 = fabs(EtaL1L2 - EtaL2L4);

                    // DeltaEta 3rd Combination ( ex) L14-L24 )
                    Float_t dEtapp3 = fabs(EtaL1L4 - EtaL2L4);

                    // DeltaEta-PV 1st Combination ( ex) PVL4 - PVL2 )
                    Float_t dEtaPv1 = fabs(EtaPVL4 - EtaPVL2);

                    // DeltaEta-PV 2nd Combination ( ex) PVL4 - PVL1 )
                    Float_t dEtaPv2 = fabs(EtaPVL4 - EtaPVL1);
                    
                    // DeltaPhi difference
                    Float_t dphi1 = deltaPhi(PhiPVL1, PhiL1L2);
                    Float_t dphi2 = deltaPhi(PhiL1L2, PhiL2L4);
                    Float_t ddphi = dphi1 - dphi2;
                    
                    if( deltaZ2 < dzCut2 && deltaZ3 < dzCut3 )
                        if( dEtapp1 < dEtappCut1 && dEtapp2 < dEtappCut2 && dEtapp3 < dEtappCut3 )
                            if( dEtaPv1 < dEtaPvCut1 && dEtaPv2 < dEtaPvCut2 )
                                if( (dphi1 > 0 && dphi2 > 0) || (dphi1 < 0 && dphi2 < 0) )
                                    if( ddphi > ddPhi_left && ddphi < ddPhi_right )
                                        L124.push_back( track((*c).X(), (*b).X(), (*a).X(),
                                                              (*c).Y(), (*b).Y(), (*a).Y(),
                                                              (*c).Z(), (*b).Z(), (*a).Z(),
                                                              2) );
                } // third hit loop
            }
        } // second hit loop
    } // first hit loop
}

void test::IsoWith_1st3rd4th(int region, float recoPV)
{
    Float_t dzCut1 = 0., dzCut2 = 0., dzCut3 = 0.;
    Float_t dEtappCut1 = 0., dEtappCut2 = 0., dEtappCut3 = 0.;
    Float_t dEtaPvCut1 = 0., dEtaPvCut2 = 0.;
    Float_t ddPhi_left = 0., ddPhi_right = 0.;

    if( region == 1 ) {
        dzCut1 = 0.0305; dzCut2 = 0.0315; dzCut3 = 0.0345; 
        dEtappCut1 = 0.00093; dEtappCut2 = 0.00315; dEtappCut3 = 0.00225; 
        dEtaPvCut1 = 0.00083; dEtaPvCut2 = 0.00185;
        ddPhi_left = -0.0025; ddPhi_right = 0.0025;
    }
    if( region == 2 ) {
        dzCut1 = 0.0305; dzCut2 = 0.0315; dzCut3 = 0.0315; 
        dEtappCut1 = 0.00066; dEtappCut2 = 0.00186; dEtappCut3 = 0.00112; 
        dEtaPvCut1 = 0.00044; dEtaPvCut2 = 0.00110;
        ddPhi_left = -0.002; ddPhi_right = 0.002;
    }
    if( region == 3 ) {
        dzCut1 = 0.0185; dzCut2 = 0.0180; dzCut3 = 0.0285; 
        dEtappCut1 = 0.00110; dEtappCut2 = 0.00278; dEtappCut3 = 0.00274; 
        dEtaPvCut1 = 0.00052; dEtaPvCut2 = 0.00092;
        ddPhi_left = -0.002; ddPhi_right = 0.002;
    }
    if( region == 4 ) {
        dzCut1 = 0.0280; dzCut2 = 0.0270; dzCut3 = 0.0970; 
        dEtappCut1 = 0.00186; dEtappCut2 = 0.00376; dEtappCut3 = 0.00328; 
        dEtaPvCut1 = 0.00072; dEtaPvCut2 = 0.00092;
        ddPhi_left = -0.002; ddPhi_right = 0.0025;
    }
    if( region == 5 ) {
        dzCut1 = 0.157; dzCut2 = 0.111; dzCut3 = 0.192; 
        dEtappCut1 = 0.00300; dEtappCut2 = 0.00698; dEtappCut3 = 0.00330; 
        dEtaPvCut1 = 0.00110; dEtaPvCut2 = 0.00152;
        ddPhi_left = -0.002; ddPhi_right = 0.003;
    }
    if( region == 6 ) {
        dzCut1 = 0.189; dzCut2 = 0.149; dzCut3 = 0.2; 
        dEtappCut1 = 0.00332; dEtappCut2 = 0.00904; dEtappCut3 = 0.00538; 
        dEtaPvCut1 = 0.00136; dEtaPvCut2 = 0.00182;
        ddPhi_left = -0.0025; ddPhi_right = 0.0035;
    }
    
    for(std::vector<TVector3>::iterator a = first_hits.begin(); a != first_hits.end(); ++a)
    {
        TVector3 PVL1;
        PVL1.SetXYZ( (*a).X(), (*a).Y(), (*a).Z() - recoPV );
        Float_t PhiPVL1 = PVL1.Phi();
        Float_t EtaPVL1 = PVL1.Eta();
        
        Float_t R1 = sqrt(pow((*a).X(),2)+pow((*a).Y(),2));
        Float_t Z1 = (*a).Z();
        
        for(std::vector<TVector3>::iterator b = third_hits.begin(); b != third_hits.end(); ++b)
        {
            TVector3 L1L3;
            L1L3.SetXYZ( (*b).X()-(*a).X(), (*b).Y()-(*a).Y(), (*b).Z()-(*a).Z() );
            Float_t PhiL1L3 = L1L3.Phi();
            Float_t EtaL1L3 = L1L3.Eta();
            
            TVector3 PVL3;
            PVL3.SetXYZ( (*b).X(), (*b).Y(), (*b).Z() - recoPV );
            Float_t EtaPVL3 = PVL3.Eta();
            
            Float_t R3 = sqrt(pow((*b).X(),2)+pow((*b).Y(),2));
            Float_t Z3 = (*b).Z();
            
            // DeltaZ 1st Combination ( ex) L1-L3 )
            Float_t vertex1 = (R3*Z1-R1*Z3)/(R3-R1);
            Float_t deltaZ1 = fabs(vertex1-recoPV);

            if( deltaZ1 < dzCut1 )
            {
                for(std::vector<TVector3>::iterator c = fourth_hits.begin(); c != fourth_hits.end(); ++c)
                {
                    TVector3 L1L4;
                    L1L4.SetXYZ( (*c).X()-(*a).X(), (*c).Y()-(*a).Y(), (*c).Z()-(*a).Z() );
                    Float_t EtaL1L4 = L1L4.Eta();

                    TVector3 L3L4;
                    L3L4.SetXYZ( (*c).X()-(*b).X(), (*c).Y()-(*b).Y(), (*c).Z()-(*b).Z() );
                    Float_t PhiL3L4 = L3L4.Phi();
                    Float_t EtaL3L4 = L3L4.Eta();

                    TVector3 PVL4;
                    PVL4.SetXYZ( (*c).X(), (*c).Y(), (*c).Z() - recoPV );
                    Float_t EtaPVL4 = PVL4.Eta();

                    Float_t R4 = sqrt(pow((*c).X(),2)+pow((*c).Y(),2));
                    Float_t Z4 = (*c).Z();
                    
                    // DeltaZ 2nd Combination ( ex) L1-L4 )
                    Float_t vertex2 = (R4*Z1-R1*Z4)/(R4-R1);
                    Float_t deltaZ2 = fabs(vertex2-recoPV);

                    // DeltaZ 3rd Combination ( ex) L3-L4 )
                    Float_t vertex3 = (R4*Z3-R3*Z4)/(R4-R3);
                    Float_t deltaZ3 = fabs(vertex3-recoPV);
                    
                    // DeltaEta 1st Combination ( ex) L13-L14 )
                    Float_t dEtapp1 = fabs(EtaL1L3 - EtaL1L4);

                    // DeltaEta 2nd Combination ( ex) L13-L34 )
                    Float_t dEtapp2 = fabs(EtaL1L3 - EtaL3L4);

                    // DeltaEta 3rd Combination ( ex) L14-L34 )
                    Float_t dEtapp3 = fabs(EtaL1L4 - EtaL3L4);

                    // DeltaEta-PV 1st Combination ( ex) PVL4 - PVL3 )
                    Float_t dEtaPv1 = fabs(EtaPVL4 - EtaPVL3);

                    // DeltaEta-PV 2nd Combination ( ex) PVL4 - PVL1 )
                    Float_t dEtaPv2 = fabs(EtaPVL4 - EtaPVL1);
                    
                    // DeltaPhi difference
                    Float_t dphi1 = deltaPhi(PhiPVL1, PhiL1L3);
                    Float_t dphi2 = deltaPhi(PhiL1L3, PhiL3L4);
                    Float_t ddphi = dphi1 - dphi2;
                    
                    if( deltaZ2 < dzCut2 && deltaZ3 < dzCut3 )
                        if( dEtapp1 < dEtappCut1 && dEtapp2 < dEtappCut2 && dEtapp3 < dEtappCut3 )
                            if( dEtaPv1 < dEtaPvCut1 && dEtaPv2 < dEtaPvCut2 )
                                if( (dphi1 > 0 && dphi2 > 0) || (dphi1 < 0 && dphi2 < 0) )
                                    if( ddphi > ddPhi_left && ddphi < ddPhi_right )
                                        L134.push_back( track((*c).X(), (*b).X(), (*a).X(),
                                                              (*c).Y(), (*b).Y(), (*a).Y(),
                                                              (*c).Z(), (*b).Z(), (*a).Z(),
                                                              3) );
                } // third hit loop
            }
        } // second hit loop
    } // first hit loop
}

void test::IsoWith_2nd3rd4th(int region, float recoPV)
{
    Float_t dzCut1 = 0., dzCut2 = 0., dzCut3 = 0.;
    Float_t dEtappCut1 = 0., dEtappCut2 = 0., dEtappCut3 = 0.;
    Float_t dEtaPvCut1 = 0., dEtaPvCut2 = 0.;
    Float_t ddPhi_left = 0., ddPhi_right = 0.;

    if( region == 1 ) {
        dzCut1 = 0.0325; dzCut2 = 0.0315; dzCut3 = 0.0345; 
        dEtappCut1 = 0.00181; dEtappCut2 = 0.00365; dEtappCut3 = 0.00225; 
        dEtaPvCut1 = 0.00083; dEtaPvCut2 = 0.00101;
        ddPhi_left = -0.0015; ddPhi_right = 0.003;
    }
    if( region == 2 ) {
        dzCut1 = 0.0315; dzCut2 = 0.0305; dzCut3 = 0.0315; 
        dEtappCut1 = 0.00090; dEtappCut2 = 0.00186; dEtappCut3 = 0.00138; 
        dEtaPvCut1 = 0.00044; dEtaPvCut2 = 0.00068;
        ddPhi_left = -0.002; ddPhi_right = 0.003;
    }
    if( region == 3 ) {
        dzCut1 = 0.0280; dzCut2 = 0.0190; dzCut3 = 0.0490; 
        dEtappCut1 = 0.00322; dEtappCut2 = 0.00184; dEtappCut3 = 0.00226; 
        dEtaPvCut1 = 0.00052; dEtaPvCut2 = 0.00072;
        ddPhi_left = -0.002; ddPhi_right = 0.003;
    }
    if( region == 4 ) {
        dzCut1 = 0.104; dzCut2 = 0.0730; dzCut3 = 0.106; 
        dEtappCut1 = 0.00318; dEtappCut2 = 0.00550; dEtappCut3 = 0.00268; 
        dEtaPvCut1 = 0.00072; dEtaPvCut2 = 0.00088;
        ddPhi_left = -0.002; ddPhi_right = 0.004;
    }
    if( region == 5 ) {
        dzCut1 = 0.189; dzCut2 = 0.165; dzCut3 = 0.169; 
        dEtappCut1 = 0.00530; dEtappCut2 = 0.00898; dEtappCut3 = 0.00330; 
        dEtaPvCut1 = 0.00110; dEtaPvCut2 = 0.00130;
        ddPhi_left = -0.003; ddPhi_right = 0.004;
    }
    if( region == 6 ) {
        dzCut1 = 0.2; dzCut2 = 0.190; dzCut3 = 0.2; 
        dEtappCut1 = 0.00724; dEtappCut2 = 0.01; dEtappCut3 = 0.00542; 
        dEtaPvCut1 = 0.00136; dEtaPvCut2 = 0.00154;
        ddPhi_left = -0.003; ddPhi_right = 0.004;
    }

    for(std::vector<TVector3>::iterator a = second_hits.begin(); a != second_hits.end(); ++a)
    {
        TVector3 PVL2;
        PVL2.SetXYZ( (*a).X(), (*a).Y(), (*a).Z() - recoPV );
        Float_t PhiPVL2 = PVL2.Phi();
        Float_t EtaPVL2 = PVL2.Eta();
        
        Float_t R2 = sqrt(pow((*a).X(),2)+pow((*a).Y(),2));
        Float_t Z2 = (*a).Z();
        
        for(std::vector<TVector3>::iterator b = third_hits.begin(); b != third_hits.end(); ++b)
        {
            TVector3 L2L3;
            L2L3.SetXYZ( (*b).X()-(*a).X(), (*b).Y()-(*a).Y(), (*b).Z()-(*a).Z() );
            Float_t PhiL2L3 = L2L3.Phi();
            Float_t EtaL2L3 = L2L3.Eta();
            
            TVector3 PVL3;
            PVL3.SetXYZ( (*b).X(), (*b).Y(), (*b).Z() - recoPV );
            Float_t EtaPVL3 = PVL3.Eta();
            
            Float_t R3 = sqrt(pow((*b).X(),2)+pow((*b).Y(),2));
            Float_t Z3 = (*b).Z();
            
            // DeltaZ 1st Combination ( ex) L2-L3 )
            Float_t vertex1 = (R3*Z2-R2*Z3)/(R3-R2);
            Float_t deltaZ1 = fabs(vertex1-recoPV);

            if( deltaZ1 < dzCut1 )
            {
                for(std::vector<TVector3>::iterator c = fourth_hits.begin(); c != fourth_hits.end(); ++c)
                {
                    TVector3 L2L4;
                    L2L4.SetXYZ( (*c).X()-(*a).X(), (*c).Y()-(*a).Y(), (*c).Z()-(*a).Z() );
                    Float_t EtaL2L4 = L2L4.Eta();

                    TVector3 L3L4;
                    L3L4.SetXYZ( (*c).X()-(*b).X(), (*c).Y()-(*b).Y(), (*c).Z()-(*b).Z() );
                    Float_t PhiL3L4 = L3L4.Phi();
                    Float_t EtaL3L4 = L3L4.Eta();

                    TVector3 PVL4;
                    PVL4.SetXYZ( (*c).X(), (*c).Y(), (*c).Z() - recoPV );
                    Float_t EtaPVL4 = PVL4.Eta();

                    Float_t R4 = sqrt(pow((*c).X(),2)+pow((*c).Y(),2));
                    Float_t Z4 = (*c).Z();
                    
                    // DeltaZ 2nd Combination ( ex) L2-L4 )
                    Float_t vertex2 = (R4*Z2-R2*Z4)/(R4-R2);
                    Float_t deltaZ2 = fabs(vertex2-recoPV);

                    // DeltaZ 3rd Combination ( ex) L3-L4 )
                    Float_t vertex3 = (R4*Z3-R3*Z4)/(R4-R3);
                    Float_t deltaZ3 = fabs(vertex3-recoPV);
                    
                    // DeltaEta 1st Combination ( ex) L23-L24 )
                    Float_t dEtapp1 = fabs(EtaL2L3 - EtaL2L4);

                    // DeltaEta 2nd Combination ( ex) L23-L34 )
                    Float_t dEtapp2 = fabs(EtaL2L3 - EtaL3L4);

                    // DeltaEta 3rd Combination ( ex) L24-L34 )
                    Float_t dEtapp3 = fabs(EtaL2L4 - EtaL3L4);

                    // DeltaEta-PV 1st Combination ( ex) PVL4 - PVL3 )
                    Float_t dEtaPv1 = fabs(EtaPVL4 - EtaPVL3);

                    // DeltaEta-PV 2nd Combination ( ex) PVL4 - PVL2 )
                    Float_t dEtaPv2 = fabs(EtaPVL4 - EtaPVL2);
                    
                    // DeltaPhi difference
                    Float_t dphi1 = deltaPhi(PhiPVL2, PhiL2L3);
                    Float_t dphi2 = deltaPhi(PhiL2L3, PhiL3L4);
                    Float_t ddphi = dphi1 - dphi2;
                    
                    if( deltaZ2 < dzCut2 && deltaZ3 < dzCut3 )
                        if( dEtapp1 < dEtappCut1 && dEtapp2 < dEtappCut2 && dEtapp3 < dEtappCut3 )
                            if( dEtaPv1 < dEtaPvCut1 && dEtaPv2 < dEtaPvCut2 )
                                if( (dphi1 > 0 && dphi2 > 0) || (dphi1 < 0 && dphi2 < 0) )
                                    if( ddphi > ddPhi_left && ddphi < ddPhi_right )
                                        L234.push_back( track((*c).X(), (*b).X(), (*a).X(),
                                                              (*c).Y(), (*b).Y(), (*a).Y(),
                                                              (*c).Z(), (*b).Z(), (*a).Z(),
                                                              4) );
                } // third hit loop
            }
        } // second hit loop
    } // first hit loop
}

#endif // #ifdef test_cxx

