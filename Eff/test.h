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

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nVtx;
   Int_t           nMeanPU;
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
   vector<float>   *propgenElPartx;
   vector<float>   *propgenElParty;
   vector<float>   *propgenElPartz;
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
   Float_t         maxBrempos;
   Float_t         lastSimtkpt;
   Float_t         initialSimtkpt;
   Int_t           bremflag;
   vector<float>   *Brempos_radius;
   vector<float>   *Brempos_x;
   vector<float>   *Brempos_y;
   vector<float>   *Brempos_z;
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
   vector<float>   *egCrysE;
   vector<float>   *egCrysEt;
   vector<float>   *egCrysEta;
   vector<float>   *egCrysPhi;
   vector<float>   *egCrysGx;
   vector<float>   *egCrysGy;
   vector<float>   *egCrysGz;
   vector<float>   *egCrysClusterE;
   vector<float>   *egCrysClusterEt;
   vector<float>   *egCrysClusterEta;
   vector<float>   *egCrysClusterPhi;
   vector<float>   *egCrysClusterGx;
   vector<float>   *egCrysClusterGy;
   vector<float>   *egCrysClusterGz;
   vector<float>   *egCrysClusterPGx;
   vector<float>   *egCrysClusterPGy;
   vector<float>   *egCrysClusterPGz;
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

   UShort_t        egN;
   vector<float>   *egEt;
   vector<float>   *egEta;
   vector<float>   *egPhi;
   vector<float>   *egGx;
   vector<float>   *egGy;
   vector<float>   *egGz;
   vector<short>   *egIEt;
   vector<short>   *egIEta;
   vector<short>   *egIPhi;
   vector<short>   *egIso;
   vector<short>   *egBx;
   vector<short>   *egTowerIPhi;
   vector<short>   *egTowerIEta;
   vector<short>   *egRawEt;
   vector<short>   *egIsoEt;
   vector<short>   *egFootprintEt;
   vector<short>   *egNTT;
   vector<short>   *egShape;
   vector<short>   *egTowerHoE;

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
   TBranch        *b_nVtx;   //!
   TBranch        *b_nMeanPU;   //!
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
   TBranch        *b_propgenElPartx;   //!
   TBranch        *b_propgenElParty;   //!
   TBranch        *b_propgenElPartz;   //!
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
   TBranch        *b_maxBrempos;   //!
   TBranch        *b_lastSimtkpt;   //!
   TBranch        *b_initialSimtkpt;   //!
   TBranch        *b_bremflag;   //!
   TBranch        *b_Brempos_radius;   //!
   TBranch        *b_Brempos_x;   //!
   TBranch        *b_Brempos_y;   //!
   TBranch        *b_Brempos_z;   //!
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
   TBranch        *b_egCrysE;   //!
   TBranch        *b_egCrysEt;   //!
   TBranch        *b_egCrysEta;   //!
   TBranch        *b_egCrysPhi;   //!
   TBranch        *b_egCrysGx;   //!
   TBranch        *b_egCrysGy;   //!
   TBranch        *b_egCrysGz;   //!
   TBranch        *b_egCrysClusterE;   //!
   TBranch        *b_egCrysClusterEt;   //!
   TBranch        *b_egCrysClusterEta;   //!
   TBranch        *b_egCrysClusterPhi;   //!
   TBranch        *b_egCrysClusterGx;   //!
   TBranch        *b_egCrysClusterGy;   //!
   TBranch        *b_egCrysClusterGz;   //!
   TBranch        *b_egCrysClusterPGx;   //!
   TBranch        *b_egCrysClusterPGy;   //!
   TBranch        *b_egCrysClusterPGz;   //!
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
   TBranch        *b_egN;   //!
   TBranch        *b_egEt;   //!
   TBranch        *b_egEta;   //!
   TBranch        *b_egPhi;   //!
   TBranch        *b_egGx;   //!
   TBranch        *b_egGy;   //!
   TBranch        *b_egGz;   //!
   TBranch        *b_egIEt;   //!
   TBranch        *b_egIEta;   //!
   TBranch        *b_egIPhi;   //!
   TBranch        *b_egIso;   //!
   TBranch        *b_egBx;   //!
   TBranch        *b_egTowerIPhi;   //!
   TBranch        *b_egTowerIEta;   //!
   TBranch        *b_egRawEt;   //!
   TBranch        *b_egIsoEt;   //!
   TBranch        *b_egFootprintEt;   //!
   TBranch        *b_egNTT;   //!
   TBranch        *b_egShape;   //!
   TBranch        *b_egTowerHoE;   //!
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


   int Ele, Pos;
   int skip;
   float EgN;
   int eta_region;

   int withoutEM_count_Ele, withEM_count_Ele;
   bool PixTrkPassed;
   int pass_count_wo4thPix, pass_count_wo3thPix, pass_count_wo2thPix, pass_count_wo1thPix;
   int woEM_pass_Ele_count_wo4thPix, woEM_pass_Ele_count_wo3thPix, woEM_pass_Ele_count_wo2thPix, woEM_pass_Ele_count_wo1thPix;
   int wEM_pass_Ele_count_wo4thPix, wEM_pass_Ele_count_wo3thPix, wEM_pass_Ele_count_wo2thPix, wEM_pass_Ele_count_wo1thPix;

   bool EM12, EM13, EM23, EM24,EM34;
   bool SA012, SA013, SA023, SA123, SA124, SA134, SA234;

   double all_cut_pass_eg;
   int all_cut_pass_Ele, withoutEM_pass_Ele, withEM_pass_Ele;
   int all_cut_pass_Pos;
   int fourth_layer_missing;
   int third_layer_missing;
   int second_layer_missing;
   int first_layer_missing;

   int bit1;
   int trigger_bit_width_;
   int pix_comb_;

   bool debug;

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
   
   double  L013_DPhi_cut1, L013_DPhi_cut2;
   
   double  L014_DPhi_cut1, L014_DPhi_cut2;
   
   double  L023_DPhi_cut1, L023_DPhi_cut2;
   
   double  L024_DPhi_cut1, L024_DPhi_cut2;
   
   double  L034_DPhi_cut1, L034_DPhi_cut2;
   
   double  L123_DPhi_cut1, L123_DPhi_cut2;
   double  L123_DEta_cut1, L123_DEta_cut2;
 
   double  L124_DPhi_cut1, L124_DPhi_cut2;
   double  L124_DEta_cut1, L124_DEta_cut2;
   
   double  L134_DPhi_cut1, L134_DPhi_cut2;
   double  L134_DEta_cut1, L134_DEta_cut2;
   
   double  L234_DPhi_cut1, L234_DPhi_cut2;
   double  L234_DEta_cut1, L234_DEta_cut2;

   double  L12_eta_upper, L13_eta_upper, L14_eta_upper, L23_eta_upper, L24_eta_upper, L34_eta_upper;
   double  L12_phi_upper, L13_phi_upper, L14_phi_upper, L23_phi_upper, L24_phi_upper, L34_phi_upper;
   double  L12_eta_bellow, L13_eta_bellow, L14_eta_bellow, L23_eta_bellow, L24_eta_bellow, L34_eta_bellow;
   double  L12_phi_bellow, L13_phi_bellow, L14_phi_bellow, L23_phi_bellow, L24_phi_bellow, L34_phi_bellow;
   double  L12_R_bellow, L13_R_bellow, L14_R_bellow, L23_R_bellow, L24_R_bellow, L34_R_bellow;

   double dPhi;
   double dEta;
   double dPhi_1, dPhi_2, dPhi_3;
   double dEta_1, dEta_2, dEta_3;
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


   inline float deltaPhi(float phi1, float phi2) { 
     float result = phi1 - phi2;
     while (result > float(M_PI)) result -= float(2*M_PI);
     while (result <= -float(M_PI)) result += float(2*M_PI);
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

   TFile *outfile;
   TTree* pixtrk_tree;

   int count_Entry; 
   int pass_egobjects_check;
   int ntnEg2; 
   int event_denominator; 
   int event_nominator; 

   int nPix123_segments;
   int nPix124_segments;
   int nPix134_segments;
   int nPix234_segments;

   vector<float> ntEgEt; 
   vector<float> ntEgEta; 
   vector<float> ntEgPhi; 

   vector<float> ntL1TkEgEt; 
   vector<float> ntL1TkEgEta; 
   vector<float> ntL1TkEgPhi; 

   vector<float> ntL1TkLooseEgEt; 
   vector<float> ntL1TkLooseEgEta; 
   vector<float> ntL1TkLooseEgPhi; 

   float matchedEgEt; 
   float matchedEgEta; 
   float matchedEgPhi; 
   int   fired;

   vector<int> PiXTRKbit;
   vector<int> pix_comb;
   vector<int> trigger_bit_width;

   vector<bool> ntCl_match; 
   vector<bool> withoutEM_match; 
   vector<bool> withEM_match; 

   vector<int> ntfirstPix;  
   vector<int> ntsecondPix; 
   vector<int> ntthirdPix;  
   vector<int> ntfourthPix; 

   float nt_lastSimtkpt;
   float nt_initialSimtkpt;

   float nt_genPhi;
   float nt_genEta;
   float nt_genPt;

   vector<bool> nt_EM12, nt_EM13, nt_EM23, nt_EM24, nt_EM34;
   vector<bool> nt_SA012, nt_SA013, nt_SA023, nt_SA123, nt_SA124, nt_SA134, nt_SA234;

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
   skip = 0;

   outfile = new TFile("../output_tmp/Tree_output","recreate");
   pixtrk_tree = new TTree("t","t");

   count_Entry = 1;
   pixtrk_tree->Branch("totalEvent", &count_Entry, "count_Entry/I");   pixtrk_tree->Branch("totalEgN", &EgN, "EgN/F");
   pixtrk_tree->Branch("ntnEg2", &ntnEg2, "ntnEg2/I");

   pixtrk_tree->Branch("ntEgEt",&ntEgEt);
   pixtrk_tree->Branch("ntEgEta",&ntEgEta);
   pixtrk_tree->Branch("ntEgPhi",&ntEgPhi);

   pixtrk_tree->Branch("ntL1TkEgEt",&ntL1TkEgEt);
   pixtrk_tree->Branch("ntL1TkEgEta",&ntL1TkEgEta);
   pixtrk_tree->Branch("ntL1TkEgPhi",&ntL1TkEgPhi);

   pixtrk_tree->Branch("ntL1TkLooseEgEt",&ntL1TkLooseEgEt);
   pixtrk_tree->Branch("ntL1TkLooseEgEta",&ntL1TkLooseEgEta);
   pixtrk_tree->Branch("ntL1TkLooseEgPhi",&ntL1TkLooseEgPhi);

   pixtrk_tree->Branch("PiXTRKbit",&PiXTRKbit);
   pixtrk_tree->Branch("trigger_bit_width",&trigger_bit_width);
   pixtrk_tree->Branch("pix_comb",&pix_comb);

   pixtrk_tree->Branch("nPix123_segments",&nPix123_segments,"nPix123_segments/I");
   pixtrk_tree->Branch("nPix124_segments",&nPix124_segments,"nPix124_segments/I");
   pixtrk_tree->Branch("nPix134_segments",&nPix134_segments,"nPix134_segments/I");
   pixtrk_tree->Branch("nPix234_segments",&nPix234_segments,"nPix234_segments/I");

   pixtrk_tree->Branch("ntCl_match",&ntCl_match);
   pixtrk_tree->Branch("withoutEM_match",&withoutEM_match);
   pixtrk_tree->Branch("withEM_match",&withEM_match);

   pixtrk_tree->Branch("ntfirstPix",&ntfirstPix);
   pixtrk_tree->Branch("ntsecondPix",&ntsecondPix);
   pixtrk_tree->Branch("ntthirdPix",&ntthirdPix);
   pixtrk_tree->Branch("ntfourthPix",&ntfourthPix);

   pixtrk_tree->Branch("nt_genPhi",&nt_genPhi,"nt_genPhi/F");
   pixtrk_tree->Branch("nt_genEta",&nt_genEta,"nt_genEta/F");
   pixtrk_tree->Branch("nt_genPt",&nt_genPt,"nt_genPt/F");

   pixtrk_tree->Branch("matchedEgEt",&matchedEgEt,"matchedEgEt/F");
   pixtrk_tree->Branch("matchedEgEta",&matchedEgEta,"matchedEgEta/F");
   pixtrk_tree->Branch("matchedEgPhi",&matchedEgPhi,"matchedEgPhi/F");

   pixtrk_tree->Branch("nt_lastSimtkpt",&nt_lastSimtkpt,"nt_lastSimtkpt/F");
   pixtrk_tree->Branch("nt_initialSimtkpt",&nt_initialSimtkpt,"nt_initialSimtkpt/F");

   pixtrk_tree->Branch("fired",&fired,"fired/I");

   pixtrk_tree->Branch("nt_EM12",&nt_EM12);
   pixtrk_tree->Branch("nt_EM13",&nt_EM13);
   pixtrk_tree->Branch("nt_EM23",&nt_EM23);
   pixtrk_tree->Branch("nt_EM24",&nt_EM24);
   pixtrk_tree->Branch("nt_EM34",&nt_EM34);
   pixtrk_tree->Branch("nt_SA012",&nt_SA012);
   pixtrk_tree->Branch("nt_SA013",&nt_SA013);
   pixtrk_tree->Branch("nt_SA023",&nt_SA023);
   pixtrk_tree->Branch("nt_SA123",&nt_SA123);
   pixtrk_tree->Branch("nt_SA124",&nt_SA124);
   pixtrk_tree->Branch("nt_SA134",&nt_SA134);
   pixtrk_tree->Branch("nt_SA234",&nt_SA234);

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
   propgenElPartx = 0;
   propgenElParty = 0;
   propgenElPartz = 0;
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
   Brempos_radius = 0;
   Brempos_x = 0;
   Brempos_y = 0;
   Brempos_z = 0;
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
   egCrysE = 0;
   egCrysEt = 0;
   egCrysEta = 0;
   egCrysPhi = 0;
   egCrysGx = 0;
   egCrysGy = 0;
   egCrysGz = 0;
   egCrysClusterE = 0;
   egCrysClusterEt = 0;
   egCrysClusterEta = 0;
   egCrysClusterPhi = 0;
   egCrysClusterGx = 0;
   egCrysClusterGy = 0;
   egCrysClusterGz = 0;
   egCrysClusterPGx = 0;
   egCrysClusterPGy = 0;
   egCrysClusterPGz = 0;
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
   egEt = 0;
   egEta = 0;
   egPhi = 0;
   egGx = 0;
   egGy = 0;
   egGz = 0;
   egIEt = 0;
   egIEta = 0;
   egIPhi = 0;
   egIso = 0;
   egBx = 0;
   egTowerIPhi = 0;
   egTowerIEta = 0;
   egRawEt = 0;
   egIsoEt = 0;
   egFootprintEt = 0;
   egNTT = 0;
   egShape = 0;
   egTowerHoE = 0;
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

   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nMeanPU", &nMeanPU, &b_nMeanPU);
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
   fChain->SetBranchAddress("propgenElPartx", &propgenElPartx, &b_propgenElPartx);
   fChain->SetBranchAddress("propgenElParty", &propgenElParty, &b_propgenElParty);
   fChain->SetBranchAddress("propgenElPartz", &propgenElPartz, &b_propgenElPartz);
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
   fChain->SetBranchAddress("maxBrempos", &maxBrempos, &b_maxBrempos);
   fChain->SetBranchAddress("lastSimtkpt", &lastSimtkpt, &b_lastSimtkpt);
   fChain->SetBranchAddress("initialSimtkpt", &initialSimtkpt, &b_initialSimtkpt);
   fChain->SetBranchAddress("bremflag", &bremflag, &b_bremflag);
   fChain->SetBranchAddress("Brempos_radius", &Brempos_radius, &b_Brempos_radius);
   fChain->SetBranchAddress("Brempos_x", &Brempos_x, &b_Brempos_x);
   fChain->SetBranchAddress("Brempos_y", &Brempos_y, &b_Brempos_y);
   fChain->SetBranchAddress("Brempos_z", &Brempos_z, &b_Brempos_z);
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
   fChain->SetBranchAddress("egCrysE", &egCrysE, &b_egCrysE);
   fChain->SetBranchAddress("egCrysEt", &egCrysEt, &b_egCrysEt);
   fChain->SetBranchAddress("egCrysEta", &egCrysEta, &b_egCrysEta);
   fChain->SetBranchAddress("egCrysPhi", &egCrysPhi, &b_egCrysPhi);
   fChain->SetBranchAddress("egCrysGx", &egCrysGx, &b_egCrysGx);
   fChain->SetBranchAddress("egCrysGy", &egCrysGy, &b_egCrysGy);
   fChain->SetBranchAddress("egCrysGz", &egCrysGz, &b_egCrysGz);
   fChain->SetBranchAddress("egCrysClusterE", &egCrysClusterE, &b_egCrysClusterE);
   fChain->SetBranchAddress("egCrysClusterEt", &egCrysClusterEt, &b_egCrysClusterEt);
   fChain->SetBranchAddress("egCrysClusterEta", &egCrysClusterEta, &b_egCrysClusterEta);
   fChain->SetBranchAddress("egCrysClusterPhi", &egCrysClusterPhi, &b_egCrysClusterPhi);
   fChain->SetBranchAddress("egCrysClusterGx", &egCrysClusterGx, &b_egCrysClusterGx);
   fChain->SetBranchAddress("egCrysClusterGy", &egCrysClusterGy, &b_egCrysClusterGy);
   fChain->SetBranchAddress("egCrysClusterGz", &egCrysClusterGz, &b_egCrysClusterGz);
   fChain->SetBranchAddress("egCrysClusterPGx", &egCrysClusterPGx, &b_egCrysClusterPGx);
   fChain->SetBranchAddress("egCrysClusterPGy", &egCrysClusterPGy, &b_egCrysClusterPGy);
   fChain->SetBranchAddress("egCrysClusterPGz", &egCrysClusterPGz, &b_egCrysClusterPGz);
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
   fChain->SetBranchAddress("egN", &egN, &b_egN);
   fChain->SetBranchAddress("egEt", &egEt, &b_egEt);
   fChain->SetBranchAddress("egEta", &egEta, &b_egEta);
   fChain->SetBranchAddress("egPhi", &egPhi, &b_egPhi);
   fChain->SetBranchAddress("egGx", &egGx, &b_egGx);
   fChain->SetBranchAddress("egGy", &egGy, &b_egGy);
   fChain->SetBranchAddress("egGz", &egGz, &b_egGz);
   fChain->SetBranchAddress("egIEt", &egIEt, &b_egIEt);
   fChain->SetBranchAddress("egIEta", &egIEta, &b_egIEta);
   fChain->SetBranchAddress("egIPhi", &egIPhi, &b_egIPhi);
   fChain->SetBranchAddress("egIso", &egIso, &b_egIso);
   fChain->SetBranchAddress("egBx", &egBx, &b_egBx);
   fChain->SetBranchAddress("egTowerIPhi", &egTowerIPhi, &b_egTowerIPhi);
   fChain->SetBranchAddress("egTowerIEta", &egTowerIEta, &b_egTowerIEta);
   fChain->SetBranchAddress("egRawEt", &egRawEt, &b_egRawEt);
   fChain->SetBranchAddress("egIsoEt", &egIsoEt, &b_egIsoEt);
   fChain->SetBranchAddress("egFootprintEt", &egFootprintEt, &b_egFootprintEt);
   fChain->SetBranchAddress("egNTT", &egNTT, &b_egNTT);
   fChain->SetBranchAddress("egShape", &egShape, &b_egShape);
   fChain->SetBranchAddress("egTowerHoE", &egTowerHoE, &b_egTowerHoE);
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
     if( third_hit == 3 ) return deltaPhi( (third_layer_hits[which_third_hit] - temp).Phi(),  temp.Phi());
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

    return deltaPhi( EM_pixelVector.Phi(), pixelVector.Phi() );
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
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"DRCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"PtErrCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"EvtCut");


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
           if( region == 1 || region == 2 || region == 3 ){
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
           if( region == 1){
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
             if( fRecHitDisk->at(a) == 2 ){ // first disk
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
             if( fRecHitDisk->at(a) == 3 ){ // second disk
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
             if( fRecHitDisk->at(a) == 1 ){ // first disk
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

dPhi_2 = EMmatchingDPhi(first_layer_hits[nthFirstHit], third_layer_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(first_layer_hits[nthFirstHit], third_layer_hits[nthThirdHit], emvector);

dPhi_3 = EMmatchingDPhi(second_layer_hits[nthSecondHit], third_layer_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(second_layer_hits[nthSecondHit], third_layer_hits[nthThirdHit], emvector);

L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );
L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );
L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );

if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Ele) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
   L12_EM_Ele = 1;
if( Signal_window_check(L12_phi_upper, dPhi_1, L12_phi_bellow, Pos) && Signal_window_check(L12_eta_upper, dEta_1, L12_eta_bellow, Ele)  )
   L12_EM_Pos = 1;

if( Signal_window_check(L13_phi_upper, dPhi_2, L13_phi_bellow, Ele) && Signal_window_check(L13_eta_upper, dEta_2, L13_eta_bellow, Ele)  )
   L13_EM_Ele = 1;
if( Signal_window_check(L13_phi_upper, dPhi_2, L13_phi_bellow, Pos) && Signal_window_check(L13_eta_upper, dEta_2, L13_eta_bellow, Ele)  )
   L13_EM_Pos = 1;

if( Signal_window_check(L23_phi_upper, dPhi_3, L23_phi_bellow, Ele) && Signal_window_check(L23_eta_upper, dEta_3, L23_eta_bellow, Ele)  )
   L23_EM_Ele = 1;
if( Signal_window_check(L23_phi_upper, dPhi_3, L23_phi_bellow, Pos) && Signal_window_check(L23_eta_upper, dEta_3, L23_eta_bellow, Ele)  )
   L23_EM_Pos = 1;

if( Signal_window_check(L123_DPhi_cut1, dPhi, L123_DPhi_cut2, Ele ) && Signal_window_check(L123_DEta_cut1, dEta, L123_DEta_cut2, Ele )  )
   L123_pass_Ele = 1;
if( Signal_window_check(L123_DPhi_cut1, dPhi, L123_DPhi_cut2, Pos ) && Signal_window_check(L123_DEta_cut1, dEta, L123_DEta_cut2, Ele )  ) 
   L123_pass_Pos = 1;

}

void test::TriggeringWithout_3rdPixel( int nthFirstHit, int nthSecondHit, int nthThirdHit){

dPhi012 = StandaloneDPhi( 0, 1, 2, 0, nthFirstHit, nthSecondHit);
dPhi014 = StandaloneDPhi( 0, 1, 4, 0, nthFirstHit, nthThirdHit );
dPhi024 = StandaloneDPhi( 0, 2, 4, 0, nthSecondHit, nthThirdHit);

dPhi_1 = EMmatchingDPhi(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);
dEta_1 = EMmatchingDEta(first_layer_hits[nthFirstHit], second_layer_hits[nthSecondHit], emvector);

dPhi_2 = EMmatchingDPhi(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);

dPhi_3 = EMmatchingDPhi(second_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(second_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);

L012_pass_Ele = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Ele );
L012_pass_Pos = Signal_window_check( L012_DPhi_cut1, dPhi012, L012_DPhi_cut2, Pos );
L014_pass_Ele = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Ele );
L014_pass_Pos = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Pos );
L024_pass_Ele = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Ele );
L024_pass_Pos = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Pos );

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

dPhi_2 = EMmatchingDPhi(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(first_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);

dPhi_3 = EMmatchingDPhi(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);

L013_pass_Ele = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Ele );
L013_pass_Pos = Signal_window_check( L013_DPhi_cut1, dPhi013, L013_DPhi_cut2, Pos );
L014_pass_Ele = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Ele );
L014_pass_Pos = Signal_window_check( L014_DPhi_cut1, dPhi014, L014_DPhi_cut2, Pos );
L034_pass_Ele = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Ele );
L034_pass_Pos = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Pos );

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

dPhi_2 = EMmatchingDPhi(second_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_2 = EMmatchingDEta(second_layer_hits[nthFirstHit], fourth_layer_hits[nthThirdHit], emvector);

dPhi_3 = EMmatchingDPhi(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);
dEta_3 = EMmatchingDEta(third_layer_hits[nthSecondHit], fourth_layer_hits[nthThirdHit], emvector);

L023_pass_Ele = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Ele );
L023_pass_Pos = Signal_window_check( L023_DPhi_cut1, dPhi023, L023_DPhi_cut2, Pos );
L024_pass_Ele = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Ele );
L024_pass_Pos = Signal_window_check( L024_DPhi_cut1, dPhi024, L024_DPhi_cut2, Pos );
L034_pass_Ele = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Ele );
L034_pass_Pos = Signal_window_check( L034_DPhi_cut1, dPhi034, L034_DPhi_cut2, Pos );

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


#endif // #ifdef test_cxx

