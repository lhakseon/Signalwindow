//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Oct  5 18:32:12 2018 by ROOT version 6.13/08
// from TTree L1PiXTRKTree/L1PiXTRKTree
// found on file: SingleElectron_NoPU.root
//////////////////////////////////////////////////////////

#ifndef sw_dist_h
#define sw_dist_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

class sw_dist {
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
   Float_t         lastSimtkpt;
   Float_t         initialSimtkpt;
   Int_t           bremflag;
   vector<float>   *Brempos_radius;
   vector<float>   *Brem_eLoss;
   vector<float>   *Brem_ptLoss;
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
   Int_t           bfastsimHitN;
   Int_t           ffastsimHitN;
   vector<int>     *bfastsimHitLayer;
   vector<float>   *bfastsimHitGx;
   vector<float>   *bfastsimHitGy;
   vector<float>   *bfastsimHitGz;
   vector<int>     *ffastsimHitLayer;
   vector<float>   *ffastsimHitGx;
   vector<float>   *ffastsimHitGy;
   vector<float>   *ffastsimHitGz;
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
   Int_t           cl3d_n;
   vector<float>   *cl3d_pt;
   vector<float>   *cl3d_energy;
   vector<float>   *cl3d_eta;
   vector<float>   *cl3d_phi;
   vector<int>     *cl3d_nclu;
   vector<float>   *cl3d_x;
   vector<float>   *cl3d_y;
   vector<int>     *cl3d_z;
   vector<float>   *cl3d_hovere;
   vector<int>     *cl3d_showerlength;
   vector<int>     *cl3d_coreshowerlength;
   vector<int>     *cl3d_firstlayer;
   vector<int>     *cl3d_maxlayer;
   vector<float>   *cl3d_seetot;
   vector<float>   *cl3d_seemax;
   vector<float>   *cl3d_spptot;
   vector<float>   *cl3d_sppmax;
   vector<float>   *cl3d_szz;
   vector<float>   *cl3d_srrtot;
   vector<float>   *cl3d_srrmax;
   vector<float>   *cl3d_srrmean;
   vector<float>   *cl3d_emaxe;
   vector<float>   *cl3d_bdteg;
   vector<float>   *cl3d_egid;
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
   TBranch        *b_lastSimtkpt;   //!
   TBranch        *b_initialSimtkpt;   //!
   TBranch        *b_bremflag;   //!
   TBranch        *b_Brempos_radius;   //!
   TBranch        *b_Brem_eLoss;   //!
   TBranch        *b_Brem_ptLoss;   //!
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
   TBranch        *b_bfastsimHitN;   //!
   TBranch        *b_ffastsimHitN;   //!
   TBranch        *b_bfastsimHitLayer;   //!
   TBranch        *b_bfastsimHitGx;   //!
   TBranch        *b_bfastsimHitGy;   //!
   TBranch        *b_bfastsimHitGz;   //!
   TBranch        *b_ffastsimHitLayer;   //!
   TBranch        *b_ffastsimHitGx;   //!
   TBranch        *b_ffastsimHitGy;   //!
   TBranch        *b_ffastsimHitGz;   //!
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
   TBranch        *b_cl3d_n;   //!
   TBranch        *b_cl3d_pt;   //!
   TBranch        *b_cl3d_energy;   //!
   TBranch        *b_cl3d_eta;   //!
   TBranch        *b_cl3d_phi;   //!
   TBranch        *b_cl3d_nclu;   //!
   TBranch        *b_cl3d_x;   //!
   TBranch        *b_cl3d_y;   //!
   TBranch        *b_cl3d_z;   //!
   TBranch        *b_cl3d_hovere;   //!
   TBranch        *b_cl3d_showerlength;   //!
   TBranch        *b_cl3d_coreshowerlength;   //!
   TBranch        *b_cl3d_firstlayer;   //!
   TBranch        *b_cl3d_maxlayer;   //!
   TBranch        *b_cl3d_seetot;   //!
   TBranch        *b_cl3d_seemax;   //!
   TBranch        *b_cl3d_spptot;   //!
   TBranch        *b_cl3d_sppmax;   //!
   TBranch        *b_cl3d_szz;   //!
   TBranch        *b_cl3d_srrtot;   //!
   TBranch        *b_cl3d_srrmax;   //!
   TBranch        *b_cl3d_srrmean;   //!
   TBranch        *b_cl3d_emaxe;   //!
   TBranch        *b_cl3d_bdteg;   //!
   TBranch        *b_cl3d_egid;   //!
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

   sw_dist(TTree *tree=0);
   virtual ~sw_dist();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   TVector3 emvector;
   float EgEt;
   float EgEta;
   float EgPhi;

   int layers[5];
   //
   std::vector<TVector3> first_layer_hits;
   std::vector<TVector3> second_layer_hits;
   std::vector<TVector3> third_layer_hits;
   std::vector<TVector3> fourth_layer_hits;

   void StorePixelHit(int region);

   inline float deltaPhi(float phi1, float phi2) {
     float result = phi1 - phi2;
     while (result > float(M_PI)) result -= float(2*M_PI);
     while (result <= -float(M_PI)) result += float(2*M_PI);
     return result;
   }


   TFile *file;
   TTree* sw_tree;

   // output
   vector<float> ntptErr;
   vector<float> ntEgEt;
   vector<float> ntEgEta;
   vector<float> ntEgPhi;

   // pixel-EG dphi
   vector<float> ntPix1EGdphi;
   vector<float> ntPix2EGdphi;
   vector<float> ntPix3EGdphi;
   vector<float> ntPix4EGdphi;

   // pixel segment-EG dphi
   vector<float> ntPix12EGdphi;
   vector<float> ntPix13EGdphi;
   vector<float> ntPix14EGdphi;
   vector<float> ntPix23EGdphi;
   vector<float> ntPix24EGdphi;
   vector<float> ntPix34EGdphi;

   // pixel segment-EG deta
   vector<float> ntPix12EGdeta;
   vector<float> ntPix13EGdeta;
   vector<float> ntPix14EGdeta;
   vector<float> ntPix23EGdeta;
   vector<float> ntPix24EGdeta;
   vector<float> ntPix34EGdeta;

   // pixel segment-pixel segment
   vector<float> ntPix012dphi;
   vector<float> ntPix013dphi;
   vector<float> ntPix014dphi;
   vector<float> ntPix023dphi;
   vector<float> ntPix024dphi;
   vector<float> ntPix034dphi;
   vector<float> ntPix123dphi;
   vector<float> ntPix124dphi;
   vector<float> ntPix134dphi;
   vector<float> ntPix234dphi;
   vector<float> ntPix123deta;
   vector<float> ntPix124deta;
   vector<float> ntPix134deta;
   vector<float> ntPix234deta;
};

#endif

#ifdef sw_dist_cxx
sw_dist::sw_dist(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("./SingleElectron_NoPU.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("./SingleElectron_NoPU.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("./SingleElectron_NoPU.root:/l1PiXTRKTree");
      dir->GetObject("L1PiXTRKTree",tree);

   }
   Init(tree);

   file = new TFile("sw.root","recreate");
   sw_tree = new TTree("t","t");

   sw_tree->Branch("ntEgEt",&ntEgEt);
   sw_tree->Branch("ntEgEta",&ntEgEta);
   sw_tree->Branch("ntEgPhi",&ntEgPhi);
   sw_tree->Branch("ntptErr",&ntptErr);

   sw_tree->Branch("ntPix1EGdphi",&ntPix1EGdphi);
   sw_tree->Branch("ntPix2EGdphi",&ntPix2EGdphi);
   sw_tree->Branch("ntPix3EGdphi",&ntPix3EGdphi);
   sw_tree->Branch("ntPix4EGdphi",&ntPix4EGdphi);

   sw_tree->Branch("ntPix12EGdphi",&ntPix12EGdphi);
   sw_tree->Branch("ntPix13EGdphi",&ntPix13EGdphi);
   sw_tree->Branch("ntPix14EGdphi",&ntPix14EGdphi);
   sw_tree->Branch("ntPix23EGdphi",&ntPix23EGdphi);
   sw_tree->Branch("ntPix24EGdphi",&ntPix24EGdphi);
   sw_tree->Branch("ntPix34EGdphi",&ntPix34EGdphi);

   sw_tree->Branch("ntPix12EGdeta",&ntPix12EGdeta);
   sw_tree->Branch("ntPix13EGdeta",&ntPix13EGdeta);
   sw_tree->Branch("ntPix14EGdeta",&ntPix14EGdeta);
   sw_tree->Branch("ntPix23EGdeta",&ntPix23EGdeta);
   sw_tree->Branch("ntPix24EGdeta",&ntPix24EGdeta);
   sw_tree->Branch("ntPix34EGdeta",&ntPix34EGdeta);

   sw_tree->Branch("ntPix012dphi",&ntPix012dphi);
   sw_tree->Branch("ntPix013dphi",&ntPix013dphi);
   sw_tree->Branch("ntPix014dphi",&ntPix014dphi);
   sw_tree->Branch("ntPix023dphi",&ntPix023dphi);
   sw_tree->Branch("ntPix024dphi",&ntPix024dphi);
   sw_tree->Branch("ntPix034dphi",&ntPix034dphi);
   sw_tree->Branch("ntPix123dphi",&ntPix123dphi);
   sw_tree->Branch("ntPix124dphi",&ntPix124dphi);
   sw_tree->Branch("ntPix134dphi",&ntPix134dphi);
   sw_tree->Branch("ntPix234dphi",&ntPix234dphi);
   sw_tree->Branch("ntPix123deta",&ntPix123deta);
   sw_tree->Branch("ntPix124deta",&ntPix124deta);
   sw_tree->Branch("ntPix134deta",&ntPix134deta);
   sw_tree->Branch("ntPix234deta",&ntPix234deta);
}

sw_dist::~sw_dist()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t sw_dist::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t sw_dist::LoadTree(Long64_t entry)
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

void sw_dist::Init(TTree *tree)
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
   Brem_eLoss = 0;
   Brem_ptLoss = 0;
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
   bfastsimHitLayer = 0;
   bfastsimHitGx = 0;
   bfastsimHitGy = 0;
   bfastsimHitGz = 0;
   ffastsimHitLayer = 0;
   ffastsimHitGx = 0;
   ffastsimHitGy = 0;
   ffastsimHitGz = 0;
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
   cl3d_pt = 0;
   cl3d_energy = 0;
   cl3d_eta = 0;
   cl3d_phi = 0;
   cl3d_nclu = 0;
   cl3d_x = 0;
   cl3d_y = 0;
   cl3d_z = 0;
   cl3d_hovere = 0;
   cl3d_showerlength = 0;
   cl3d_coreshowerlength = 0;
   cl3d_firstlayer = 0;
   cl3d_maxlayer = 0;
   cl3d_seetot = 0;
   cl3d_seemax = 0;
   cl3d_spptot = 0;
   cl3d_sppmax = 0;
   cl3d_szz = 0;
   cl3d_srrtot = 0;
   cl3d_srrmax = 0;
   cl3d_srrmean = 0;
   cl3d_emaxe = 0;
   cl3d_bdteg = 0;
   cl3d_egid = 0;
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
   fChain->SetBranchAddress("lastSimtkpt", &lastSimtkpt, &b_lastSimtkpt);
   fChain->SetBranchAddress("initialSimtkpt", &initialSimtkpt, &b_initialSimtkpt);
   fChain->SetBranchAddress("bremflag", &bremflag, &b_bremflag);
   fChain->SetBranchAddress("Brempos_radius", &Brempos_radius, &b_Brempos_radius);
   fChain->SetBranchAddress("Brem_eLoss", &Brem_eLoss, &b_Brem_eLoss);
   fChain->SetBranchAddress("Brem_ptLoss", &Brem_ptLoss, &b_Brem_ptLoss);
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
   fChain->SetBranchAddress("bfastsimHitN", &bfastsimHitN, &b_bfastsimHitN);
   fChain->SetBranchAddress("ffastsimHitN", &ffastsimHitN, &b_ffastsimHitN);
   fChain->SetBranchAddress("bfastsimHitLayer", &bfastsimHitLayer, &b_bfastsimHitLayer);
   fChain->SetBranchAddress("bfastsimHitGx", &bfastsimHitGx, &b_bfastsimHitGx);
   fChain->SetBranchAddress("bfastsimHitGy", &bfastsimHitGy, &b_bfastsimHitGy);
   fChain->SetBranchAddress("bfastsimHitGz", &bfastsimHitGz, &b_bfastsimHitGz);
   fChain->SetBranchAddress("ffastsimHitLayer", &ffastsimHitLayer, &b_ffastsimHitLayer);
   fChain->SetBranchAddress("ffastsimHitGx", &ffastsimHitGx, &b_ffastsimHitGx);
   fChain->SetBranchAddress("ffastsimHitGy", &ffastsimHitGy, &b_ffastsimHitGy);
   fChain->SetBranchAddress("ffastsimHitGz", &ffastsimHitGz, &b_ffastsimHitGz);
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
   fChain->SetBranchAddress("cl3d_n", &cl3d_n, &b_cl3d_n);
   fChain->SetBranchAddress("cl3d_pt", &cl3d_pt, &b_cl3d_pt);
   fChain->SetBranchAddress("cl3d_energy", &cl3d_energy, &b_cl3d_energy);
   fChain->SetBranchAddress("cl3d_eta", &cl3d_eta, &b_cl3d_eta);
   fChain->SetBranchAddress("cl3d_phi", &cl3d_phi, &b_cl3d_phi);
   fChain->SetBranchAddress("cl3d_nclu", &cl3d_nclu, &b_cl3d_nclu);
   fChain->SetBranchAddress("cl3d_x", &cl3d_x, &b_cl3d_x);
   fChain->SetBranchAddress("cl3d_y", &cl3d_y, &b_cl3d_y);
   fChain->SetBranchAddress("cl3d_z", &cl3d_z, &b_cl3d_z);
   fChain->SetBranchAddress("cl3d_hovere", &cl3d_hovere, &b_cl3d_hovere);
   fChain->SetBranchAddress("cl3d_showerlength", &cl3d_showerlength, &b_cl3d_showerlength);
   fChain->SetBranchAddress("cl3d_coreshowerlength", &cl3d_coreshowerlength, &b_cl3d_coreshowerlength);
   fChain->SetBranchAddress("cl3d_firstlayer", &cl3d_firstlayer, &b_cl3d_firstlayer);
   fChain->SetBranchAddress("cl3d_maxlayer", &cl3d_maxlayer, &b_cl3d_maxlayer);
   fChain->SetBranchAddress("cl3d_seetot", &cl3d_seetot, &b_cl3d_seetot);
   fChain->SetBranchAddress("cl3d_seemax", &cl3d_seemax, &b_cl3d_seemax);
   fChain->SetBranchAddress("cl3d_spptot", &cl3d_spptot, &b_cl3d_spptot);
   fChain->SetBranchAddress("cl3d_sppmax", &cl3d_sppmax, &b_cl3d_sppmax);
   fChain->SetBranchAddress("cl3d_szz", &cl3d_szz, &b_cl3d_szz);
   fChain->SetBranchAddress("cl3d_srrtot", &cl3d_srrtot, &b_cl3d_srrtot);
   fChain->SetBranchAddress("cl3d_srrmax", &cl3d_srrmax, &b_cl3d_srrmax);
   fChain->SetBranchAddress("cl3d_srrmean", &cl3d_srrmean, &b_cl3d_srrmean);
   fChain->SetBranchAddress("cl3d_emaxe", &cl3d_emaxe, &b_cl3d_emaxe);
   fChain->SetBranchAddress("cl3d_bdteg", &cl3d_bdteg, &b_cl3d_bdteg);
   fChain->SetBranchAddress("cl3d_egid", &cl3d_egid, &b_cl3d_egid);
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

void sw_dist::StorePixelHit(int region){

 // barrel pixels
 for(int a=0; a<bRecHitN; a++){

    double Dphi = 0.;
    TVector3 current_hit;
    current_hit.SetXYZ( bRecHitGx->at(a), bRecHitGy->at(a), bRecHitGz->at(a) );
    Dphi = deltaPhi( current_hit.Phi(), EgPhi);
    if(fabs(Dphi) > 0.1) continue;

    if( region == 1 ){
      if( bRecHitLayer->at(a) == 1 ){ // First layer
           layers[1]++;
           first_layer_hits.push_back(current_hit);
      }
      if( bRecHitLayer->at(a) == 2 ){ // First layer
           layers[2]++;
           second_layer_hits.push_back(current_hit);
      }
      if( bRecHitLayer->at(a) == 3 ){ // First layer
           layers[3]++;
           third_layer_hits.push_back(current_hit);
      }
      if( bRecHitLayer->at(a) == 4 ){ // First layer
           layers[4]++;
           fourth_layer_hits.push_back(current_hit);
      }
    } // |eta| < 0.8  

    if( region == 2 ){
      if( bRecHitLayer->at(a) == 1 ){ // First layer
           layers[1]++;
           first_layer_hits.push_back(current_hit);
      }   
      if( bRecHitLayer->at(a) == 2 ){ // First layer
           layers[2]++;
           second_layer_hits.push_back(current_hit);
      }   
      if( bRecHitLayer->at(a) == 3 ){ // First layer
           layers[3]++;
           third_layer_hits.push_back(current_hit);
      }   
    } // 0.8 < |eta| < 1.4 

    if( region == 3 ){
      if( bRecHitLayer->at(a) == 1 ){ // First layer
           layers[1]++;
           first_layer_hits.push_back(current_hit);
      }   
      if( bRecHitLayer->at(a) == 2 ){ // First layer
           layers[2]++;
           second_layer_hits.push_back(current_hit);
      }   
    } // 1.4 < |eta| < 1.7 

    if( region == 4 ){
      if( bRecHitLayer->at(a) == 1 ){ // First layer
           layers[1]++;
           first_layer_hits.push_back(current_hit);
      }   
    } // 1.7 < |eta| < 2.1
 }

 // disk pixels
 for(int a=0; a<fRecHitN; a++){

    double Dphi = 0.;
    TVector3 current_hit;
    current_hit.SetXYZ( fRecHitGx->at(a), fRecHitGy->at(a), fRecHitGz->at(a) );
    Dphi = deltaPhi( current_hit.Phi(), EgPhi);
    if(fabs(Dphi) > 0.1) continue;

    if( region == 2 ){
      if( fRecHitDisk->at(a) == 1 ){ // First disk
           layers[4]++;
           fourth_layer_hits.push_back(current_hit);
      }   
    } // 0.8 < |eta| < 1.4 

    if( region == 3 ){
      if( fRecHitDisk->at(a) == 1 ){ // First disk
           layers[3]++;
           third_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 2 ){ // First disk
           layers[4]++;
           fourth_layer_hits.push_back(current_hit);
      }
    } // 1.4 < |eta| < 1.7

    if( region == 4 ){
      if( fRecHitDisk->at(a) == 1 ){ // First disk
           layers[2]++;
           second_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 2 ){ // Second disk
           layers[3]++;
           third_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 3 ){ // Third disk
           layers[4]++;
           fourth_layer_hits.push_back(current_hit);
      }
    } // 1.7 < |eta| < 2.1

    if( region == 5 ){
      if( fRecHitDisk->at(a) == 1 ){ // First disk
           layers[1]++;
           first_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 2 ){ // Second disk
           layers[2]++;
           second_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 3 ){ // Third disk
           layers[3]++;
           third_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 4 ){ // Fourth disk
           layers[4]++;
           fourth_layer_hits.push_back(current_hit);
      }
    } // 2.7 < |eta| < 3.0

    if( region == 6 ){
      if( fRecHitDisk->at(a) == 2 ){ // First disk
           layers[1]++;
           first_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 3 ){ // Second disk
           layers[2]++;
           second_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 4 ){ // Third disk
           layers[3]++;
           third_layer_hits.push_back(current_hit);
      }
      if( fRecHitDisk->at(a) == 5 ){ // Fourth disk
           layers[4]++;
           fourth_layer_hits.push_back(current_hit);
      }
    } // 2.7 < |eta| < 3.0

 }

}

Bool_t sw_dist::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void sw_dist::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t sw_dist::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef sw_dist_cxx
