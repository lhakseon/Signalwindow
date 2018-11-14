//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Sep  8 11:11:07 2018 by ROOT version 6.13/08
// from TTree t/t
// found on file: sw.root
//////////////////////////////////////////////////////////

#ifndef Make2Dplots_h
#define Make2Dplots_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGaxis.h>

#include "../Style/tdrstyle.C"

class Make2Dplots {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<float>   *ntEgEt;
   vector<float>   *ntEgEta;
   vector<float>   *ntEgPhi;
   vector<float>   *ntptErr;
   vector<float>   *ntPix1EGdphi;
   vector<float>   *ntPix2EGdphi;
   vector<float>   *ntPix3EGdphi;
   vector<float>   *ntPix4EGdphi;
   vector<float>   *ntPix12EGdphi;
   vector<float>   *ntPix13EGdphi;
   vector<float>   *ntPix14EGdphi;
   vector<float>   *ntPix23EGdphi;
   vector<float>   *ntPix24EGdphi;
   vector<float>   *ntPix34EGdphi;
   vector<float>   *ntPix12EGdeta;
   vector<float>   *ntPix13EGdeta;
   vector<float>   *ntPix14EGdeta;
   vector<float>   *ntPix23EGdeta;
   vector<float>   *ntPix24EGdeta;
   vector<float>   *ntPix34EGdeta;
   vector<float>   *ntPix012dphi;
   vector<float>   *ntPix013dphi;
   vector<float>   *ntPix014dphi;
   vector<float>   *ntPix023dphi;
   vector<float>   *ntPix024dphi;
   vector<float>   *ntPix034dphi;
   vector<float>   *ntPix123dphi;
   vector<float>   *ntPix124dphi;
   vector<float>   *ntPix134dphi;
   vector<float>   *ntPix234dphi;
   vector<float>   *ntPix123deta;
   vector<float>   *ntPix124deta;
   vector<float>   *ntPix134deta;
   vector<float>   *ntPix234deta;

   // List of branches
   TBranch        *b_ntEgEt;   //!
   TBranch        *b_ntEgEta;   //!
   TBranch        *b_ntEgPhi;   //!
   TBranch        *b_ntptErr;   //!
   TBranch        *b_ntPix1EGdphi;   //!
   TBranch        *b_ntPix2EGdphi;   //!
   TBranch        *b_ntPix3EGdphi;   //!
   TBranch        *b_ntPix4EGdphi;   //!
   TBranch        *b_ntPix12EGdphi;   //!
   TBranch        *b_ntPix13EGdphi;   //!
   TBranch        *b_ntPix14EGdphi;   //!
   TBranch        *b_ntPix23EGdphi;   //!
   TBranch        *b_ntPix24EGdphi;   //!
   TBranch        *b_ntPix34EGdphi;   //!
   TBranch        *b_ntPix12EGdeta;   //!
   TBranch        *b_ntPix13EGdeta;   //!
   TBranch        *b_ntPix14EGdeta;   //!
   TBranch        *b_ntPix23EGdeta;   //!
   TBranch        *b_ntPix24EGdeta;   //!
   TBranch        *b_ntPix34EGdeta;   //!
   TBranch        *b_ntPix012dphi;   //!
   TBranch        *b_ntPix013dphi;   //!
   TBranch        *b_ntPix014dphi;   //!
   TBranch        *b_ntPix023dphi;   //!
   TBranch        *b_ntPix024dphi;   //!
   TBranch        *b_ntPix034dphi;   //!
   TBranch        *b_ntPix123dphi;   //!
   TBranch        *b_ntPix124dphi;   //!
   TBranch        *b_ntPix134dphi;   //!
   TBranch        *b_ntPix234dphi;   //!
   TBranch        *b_ntPix123deta;   //!
   TBranch        *b_ntPix124deta;   //!
   TBranch        *b_ntPix134deta;   //!
   TBranch        *b_ntPix234deta;   //!

   Make2Dplots(TTree *tree=0);
   virtual ~Make2Dplots();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int eta_);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Make2Dplots_cxx
Make2Dplots::Make2Dplots(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../sw.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../sw.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
}

Make2Dplots::~Make2Dplots()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Make2Dplots::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Make2Dplots::LoadTree(Long64_t entry)
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

void Make2Dplots::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ntEgEt = 0;
   ntEgEta = 0;
   ntEgPhi = 0;
   ntptErr = 0;
   ntPix1EGdphi = 0;
   ntPix2EGdphi = 0;
   ntPix3EGdphi = 0;
   ntPix4EGdphi = 0;
   ntPix12EGdphi = 0;
   ntPix13EGdphi = 0;
   ntPix14EGdphi = 0;
   ntPix23EGdphi = 0;
   ntPix24EGdphi = 0;
   ntPix34EGdphi = 0;
   ntPix12EGdeta = 0;
   ntPix13EGdeta = 0;
   ntPix14EGdeta = 0;
   ntPix23EGdeta = 0;
   ntPix24EGdeta = 0;
   ntPix34EGdeta = 0;
   ntPix012dphi = 0;
   ntPix013dphi = 0;
   ntPix014dphi = 0;
   ntPix023dphi = 0;
   ntPix024dphi = 0;
   ntPix034dphi = 0;
   ntPix123dphi = 0;
   ntPix124dphi = 0;
   ntPix134dphi = 0;
   ntPix234dphi = 0;
   ntPix123deta = 0;
   ntPix124deta = 0;
   ntPix134deta = 0;
   ntPix234deta = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ntEgEt", &ntEgEt, &b_ntEgEt);
   fChain->SetBranchAddress("ntEgEta", &ntEgEta, &b_ntEgEta);
   fChain->SetBranchAddress("ntEgPhi", &ntEgPhi, &b_ntEgPhi);
   fChain->SetBranchAddress("ntptErr", &ntptErr, &b_ntptErr);
   fChain->SetBranchAddress("ntPix1EGdphi", &ntPix1EGdphi, &b_ntPix1EGdphi);
   fChain->SetBranchAddress("ntPix2EGdphi", &ntPix2EGdphi, &b_ntPix2EGdphi);
   fChain->SetBranchAddress("ntPix3EGdphi", &ntPix3EGdphi, &b_ntPix3EGdphi);
   fChain->SetBranchAddress("ntPix4EGdphi", &ntPix4EGdphi, &b_ntPix4EGdphi);
   fChain->SetBranchAddress("ntPix12EGdphi", &ntPix12EGdphi, &b_ntPix12EGdphi);
   fChain->SetBranchAddress("ntPix13EGdphi", &ntPix13EGdphi, &b_ntPix13EGdphi);
   fChain->SetBranchAddress("ntPix14EGdphi", &ntPix14EGdphi, &b_ntPix14EGdphi);
   fChain->SetBranchAddress("ntPix23EGdphi", &ntPix23EGdphi, &b_ntPix23EGdphi);
   fChain->SetBranchAddress("ntPix24EGdphi", &ntPix24EGdphi, &b_ntPix24EGdphi);
   fChain->SetBranchAddress("ntPix34EGdphi", &ntPix34EGdphi, &b_ntPix34EGdphi);
   fChain->SetBranchAddress("ntPix12EGdeta", &ntPix12EGdeta, &b_ntPix12EGdeta);
   fChain->SetBranchAddress("ntPix13EGdeta", &ntPix13EGdeta, &b_ntPix13EGdeta);
   fChain->SetBranchAddress("ntPix14EGdeta", &ntPix14EGdeta, &b_ntPix14EGdeta);
   fChain->SetBranchAddress("ntPix23EGdeta", &ntPix23EGdeta, &b_ntPix23EGdeta);
   fChain->SetBranchAddress("ntPix24EGdeta", &ntPix24EGdeta, &b_ntPix24EGdeta);
   fChain->SetBranchAddress("ntPix34EGdeta", &ntPix34EGdeta, &b_ntPix34EGdeta);
   fChain->SetBranchAddress("ntPix012dphi", &ntPix012dphi, &b_ntPix012dphi);
   fChain->SetBranchAddress("ntPix013dphi", &ntPix013dphi, &b_ntPix013dphi);
   fChain->SetBranchAddress("ntPix014dphi", &ntPix014dphi, &b_ntPix014dphi);
   fChain->SetBranchAddress("ntPix023dphi", &ntPix023dphi, &b_ntPix023dphi);
   fChain->SetBranchAddress("ntPix024dphi", &ntPix024dphi, &b_ntPix024dphi);
   fChain->SetBranchAddress("ntPix034dphi", &ntPix034dphi, &b_ntPix034dphi);
   fChain->SetBranchAddress("ntPix123dphi", &ntPix123dphi, &b_ntPix123dphi);
   fChain->SetBranchAddress("ntPix124dphi", &ntPix124dphi, &b_ntPix124dphi);
   fChain->SetBranchAddress("ntPix134dphi", &ntPix134dphi, &b_ntPix134dphi);
   fChain->SetBranchAddress("ntPix234dphi", &ntPix234dphi, &b_ntPix234dphi);
   fChain->SetBranchAddress("ntPix123deta", &ntPix123deta, &b_ntPix123deta);
   fChain->SetBranchAddress("ntPix124deta", &ntPix124deta, &b_ntPix124deta);
   fChain->SetBranchAddress("ntPix134deta", &ntPix134deta, &b_ntPix134deta);
   fChain->SetBranchAddress("ntPix234deta", &ntPix234deta, &b_ntPix234deta);
   Notify();

   setTDRStyle();
}

Bool_t Make2Dplots::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Make2Dplots::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Make2Dplots::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Make2Dplots_cxx
