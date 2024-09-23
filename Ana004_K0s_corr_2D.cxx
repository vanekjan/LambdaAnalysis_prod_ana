//  Centrality binning:
//     Bin       Centrality (16)   Centrality (9)
//     0            75-80%            70-80%
//     1            70-75%            60-70%
//     2            65-70%            50-60%
//     3            60-65%            40-50%
//     4            55-60%            30-40%
//     5            50-55%            20-30%
//     6            45-50%            10-20%
//     7            40-45%             5-10%
//     8            35-40%             0- 5%
//     9            30-35%
//    10            25-30%
//    11            20-25%
//    12            15-20%
//    13            10-15%
//    14             5-10%
//    15             0- 5%
//centrality bins used (centrality 9):
// 0-10% (7+8), 10-40% (4+5+6), 40-80 (0+1+2+3), 0-80% (0-8)

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include"TH1.h"
#include"TH2.h"
#include"TH3.h"
#include"TF1.h"
#include"TF2.h"
#include"TF12.h"
#include"TMath.h"
#include"TCanvas.h"
#include"TFile.h"
#include"TLatex.h"
#include"TStyle.h"
#include"TPad.h"
#include"TLegend.h"
#include"TPaveText.h"
#include"TAxis.h"
#include"TTree.h"
#include"TFitResultPtr.h"
#include"TFitResult.h"
#include"TString.h"
#include"TLine.h"
#include"TChain.h"
#include"TLorentzVector.h"
#include"TGraphErrors.h"


using namespace std;

//const int nPtBins = 9;
//float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5., 7.};

const int nPtBins = 8;
float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};

const int nPtBins_corr = 2;
float const pT_bins_corr[nPtBins_corr+1] = { 0.5, 1.5, 5.};

const int nEtaBins = 3;
float const eta_bins[nEtaBins+1] = { -1, -0.2, 0.2, 1 };
//float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };

const double p_mass_PDG = 0.93827208816; //p mass on GeV/c^2 from latest PDG
const double pi_mass_PDG = 0.13957039; //p mass on GeV/c^2 from latest PDG
const double K0s_mass_PDG = 1.115683; //mass in GeV/c^2 from latest PDG

const float K0s_y_cut = 1;
//const int strictTOF_cut = 1; //0 - hybrid TOF for both daughters, 1 - strict TOF for pions, 2 - strict TOF for both pion and proton
//const float K0s_cos_theta_cut = 0.996;
//const float K0s_decayK0s_cut = 25;


Double_t LinFunc(Double_t *x, Double_t *par)
{
  if (x[0] > 1.795 && x[0] < 1.945)
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1]*x[0];
}

//calculate pz from pT and eta
double pz(double pt, double eta)
{
  return pt*TMath::SinH(eta);
}

//calculate rapidity from pT and eta
double rapidity(double pt, double eta, double mass)
{
  double p = pt*TMath::CosH(eta);
  double pz = pt*TMath::SinH(eta);
  double E = TMath::Sqrt(p*p + mass*mass);

  return 0.5*TMath::Log( (E + pz)/(E - pz) );
}

//calculate theta star for Lmabdda(-bar) pair
double LpairThetaStar(TLorentzVector *L1, TLorentzVector *p1, TLorentzVector *L2, TLorentzVector *p2)
{

  TLorentzVector const L1_reverse(-L1->Px(), -L1->Py(), -L1->Pz(), L1->E());
  p1->Boost(L1_reverse.BoostVector());


  TLorentzVector const L2_reverse(-L2->Px(), -L2->Py(), -L2->Pz(), L2->E());
  p2->Boost(L2_reverse.BoostVector());


  return p1->Angle(p2->Vect());
}

//calculate theta star for Lmabdda(-bar) pair
double LpairThetaStar(TLorentzVector L1, TLorentzVector p1, TLorentzVector L2, TLorentzVector p2)
{

  TLorentzVector const L1_reverse(-L1.Px(), -L1.Py(), -L1.Pz(), L1.E());
  p1.Boost(L1_reverse.BoostVector());


  TLorentzVector const L2_reverse(-L2.Px(), -L2.Py(), -L2.Pz(), L2.E());
  p2.Boost(L2_reverse.BoostVector());


  return p1.Angle(p2.Vect());
}

bool cuts(int year, float K0s_y, int hasTOFinfo, float pi1_DCA, float pi2_DCA, float pair_DCA, float dec_L, float cos_theta)
{

  if( !( TMath::Abs(K0s_y) < K0s_y_cut ) ) return false;
  if( year == 2017 && hasTOFinfo == 0 ) return false; //TOF matched pion 1, for Run17 only now
  //if( year == 2015 && hasTOFinfo == 0 ) return false; //TOF matched pion 1, for Run17 only now
  //if( strictTOF_cut == 2 && (pi1_hasTOFinfo == 0 || pi2_hasTOFinfo == 0) ) return false; //TOF matched pions and protons
  //if(cos(K0s_theta) < K0s_cos_theta_cut) return false;
  //if(K0s_decayL > K0s_decayK0s_cut) return false;

  if( pi1_DCA < 0.3 ) return false;
  if( pi2_DCA < 0.3 ) return false;

  if( pair_DCA > 1. ) return false;

  if( dec_L < 0.5 ) return false;

  if( cos_theta < 0.996 ) return false;

  return true;
}

//for tight topological cuts for sys. err.
bool cuts_topo_sys_err(float pi1_DCA, float pi2_DCA, float pair_DCA, float dec_L, float cos_theta)
{
  //if( pi1_DCA < 0.4 ) return false;
  //if( pi2_DCA < 0.4 ) return false;

  if( pair_DCA > 0.9 ) return false;

  if( dec_L < 0.6 ) return false;

  if( cos_theta < 0.997 ) return false;

  return true;
}

//for tight daughter pT cut for sys. err
bool cuts_pT_sys_err(float pi1_pT, float pi2_pT)
{
  if( pi1_pT < 0.17 ) return false;
  if( pi2_pT < 0.17 ) return false; //0.3 migt be too tight for pions

  return true;
}


//analyze invariant mass spectra
//arguments are vectors of bin numbers of the invariant mass peak
//the bins are determined via fit to the invariant mass spectra
//bool InvMass(TTree *K0s_tree, vector<int> &invMassBins_L, vector<int> &invMassBins_L0, vector<int> &invMassBins_L0bar, const int readMode)
bool DoAnalysis(TChain *K0s_tree, const int ReadMode = 0, const int cutType = 0, const int energy = 510, const int year = 2017)
{

  TFile *InvMassFile;

  if(cutType == 0) InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_K0s_2D_ana_cuts_work.root", year), "recreate");
  if(cutType == 1) InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_K0s_2D_tight_topo_cuts_work.root", year), "recreate");
  if(cutType == 2) InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_K0s_2D_tight_pT_cut_work.root", year), "recreate");

  //_______________________________________________________________________________________________________________________________________________


  Int_t charge;
  Float_t K0s_mass, K0s_pt, K0s_eta, K0s_phi, K0s_decayL, K0s_theta, K0s_DCAdaughters;

  Int_t pi1_InEventID, pi2_InEventID;
  Float_t pi1_pt, pi2_pt;
  Float_t pi1_eta, pi2_eta;
  Float_t pi1_phi, pi2_phi;
  Float_t pi1_dca, pi2_dca;
  Int_t pi1_ch;
  Int_t pi1_hasTOFinfo, pi2_hasTOFinfo;

  Int_t eventId;
  Float_t Vz;

  //---------------SET BARANCH ADDRESSES------------------------
  K0s_tree->SetBranchAddress("pair_charge", &charge);
  K0s_tree->SetBranchAddress("pair_mass", &K0s_mass);
  K0s_tree->SetBranchAddress("pair_pt", &K0s_pt);
  K0s_tree->SetBranchAddress("pair_eta", &K0s_eta);
  K0s_tree->SetBranchAddress("pair_phi", &K0s_phi);
  K0s_tree->SetBranchAddress("pair_decayL", &K0s_decayL);
  K0s_tree->SetBranchAddress("pair_theta", &K0s_theta);
  K0s_tree->SetBranchAddress("pair_DCAdaughters", &K0s_DCAdaughters);

  K0s_tree->SetBranchAddress("p1_InEventID", &pi1_InEventID),
  K0s_tree->SetBranchAddress("p1_pt", &pi1_pt);
  K0s_tree->SetBranchAddress("p1_eta", &pi1_eta);
  K0s_tree->SetBranchAddress("p1_phi", &pi1_phi);
  K0s_tree->SetBranchAddress("p1_dca", &pi1_dca);
  K0s_tree->SetBranchAddress("p1_ch", &pi1_ch); //charge of first pion in tha pair (can)
  K0s_tree->SetBranchAddress("p1_hasTOFinfo", &pi1_hasTOFinfo);

  K0s_tree->SetBranchAddress("p2_InEventID", &pi2_InEventID),
  K0s_tree->SetBranchAddress("p2_pt", &pi2_pt);
  K0s_tree->SetBranchAddress("p2_eta", &pi2_eta);
  K0s_tree->SetBranchAddress("p2_phi", &pi2_phi);
  K0s_tree->SetBranchAddress("p2_dca", &pi2_dca);
  K0s_tree->SetBranchAddress("p2_hasTOFinfo", &pi2_hasTOFinfo);

  K0s_tree->SetBranchAddress("eventId", &eventId);
  K0s_tree->SetBranchAddress("Vz", &Vz);

  //--------------------------------------------------------------------------


  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass_US[nPtBins_corr][nPtBins_corr];
  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //bacground when US is paired with LS
  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass_LS[nPtBins_corr][nPtBins_corr];
  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass[nPtBins_corr][nPtBins_corr];

  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass_US_proton_mass_one[nPtBins_corr][nPtBins_corr]; //for QA - how many K0s will fall into L mass
  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass_US_proton_mass_both[nPtBins_corr][nPtBins_corr]; //for QA - how many K0s will fall into L mass


  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2] = new TH2F(Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2),200, 0.45, 0.55, 200, 0.45, 0.55);
      K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2] = new TH2F(Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2),200, 0.45, 0.55, 200, 0.45, 0.55);
      K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2] = new TH2F(Form("K0s1_inv_mass_vs_K0s2_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_vs_K0s2_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2),200, 0.45, 0.55, 200, 0.45, 0.55);

      //one pair is corectly identified, other has pi misidentified as p
      K0s1_inv_mass_vs_K0s2_inv_mass_US_proton_mass_one[pTbin1][pTbin2] = new TH2F(Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_proton_mass_one_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_proton_mass_one_pT1_%i_pT2_%i", pTbin1, pTbin2), 200, 0.45, 0.55, 180, 1, 1.2);

      //both mothers have misidentified one pi as p
      K0s1_inv_mass_vs_K0s2_inv_mass_US_proton_mass_both[pTbin1][pTbin2] = new TH2F(Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_proton_mass_both_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_proton_mass_both_pT1_%i_pT2_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
    }
  }

  //daughter QA
  TH2F *K0s_K0s_pi_eta1_vs_eta2_US_hist = new TH2F("K0s_K0s_pi_eta1_vs_eta2_US_hist", "K0s_K0s_pi_eta1_vs_eta2_US_hist", 20, -1, 1, 20, -1, 1);
  TH2F *K0s_K0s_pi_phi1_vs_phi2_US_hist = new TH2F("K0s_K0s_pi_phi1_vs_phi2_US_hist", "K0s_K0s_pi_phi1_vs_phi2_US_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *K0s_K0s_pi_pT1_vs_pT2_US_hist = new TH2F("K0s_K0s_pi_pT1_vs_pT2_US_hist", "K0s_K0s_pi_pT1_vs_pT2_US_hist", 20, 0, 5, 20, 0, 5);


  TH1F *K0s_K0s_pi_Delta_eta = new TH1F("K0s_K0s_pi_Delta_eta", "K0s_K0s_pi_Delta_eta", 100, 0, 2);
  TH1F *K0s_K0s_pi_Delta_phi = new TH1F("K0s_K0s_pi_Delta_phi", "K0s_K0s_pi_Delta_phi", 100, 0, TMath::Pi());


  //________________________________________________________________________________________


  //to store Lorentz vectors for L pair analysis
  vector<TLorentzVector> K0s_vector;
  vector<TLorentzVector> K0s_vector_p_mass;
  vector<int> K0s_pT_bin_vector;
  vector<int> K0s_eta_bin_vector;

  vector<TLorentzVector> pi1_vector;
  vector<int> pi1_tag_vector;
  vector<int> pi1_charge;

  vector<TLorentzVector> pi2_vector;
  vector<int> pi2_tag_vector;


  vector<TLorentzVector> K0s_vector_background;
  vector<int> K0s_pT_bin_vector_background;
  vector<int> K0s_eta_bin_vector_background;

  vector<TLorentzVector> pi1_vector_background;
  vector<int> pi1_tag_vector_background;
  vector<int> pi1_charge_background;

  vector<TLorentzVector> pi2_vector_background;
  vector<int> pi2_tag_vector_background;

  //vectors for cos(theta*)
  vector<TLorentzVector> K0s_K0s_K0s1_mom;
  vector<TLorentzVector> K0s_K0s_K0s2_mom;
  vector<float> K0s_K0s_cos_theta;
  vector<float> K0s_K0s_Minv_K0s1;
  vector<float> K0s_K0s_Minv_K0s2;
  vector<int> K0s_K0s_pT_bin_K0s1;
  vector<int> K0s_K0s_pT_bin_K0s2;
  vector<int> K0s_K0s_eta_bin_K0s1;
  vector<int> K0s_K0s_eta_bin_K0s2;


  vector<TLorentzVector> K0s_K0s_K0s1_mom_back;
  vector<TLorentzVector> K0s_K0s_K0s2_mom_back;
  vector<float> K0s_K0s_cos_theta_back;
  vector<float> K0s_K0s_Minv_K0s1_back;
  vector<float> K0s_K0s_Minv_K0s2_back;
  vector<int> K0s_K0s_pT_bin_K0s1_back;
  vector<int> K0s_K0s_pT_bin_K0s2_back;
  vector<int> K0s_K0s_eta_bin_K0s1_back;
  vector<int> K0s_K0s_eta_bin_K0s2_back;

  //for reweight of ME - has to separate US-LS and LS-US
  //K0s1(US)-K0s2(LS) - 1, K0s1(LS)-K0s2(US) - 2
  //flag mathces index of histogram
  vector<int> K0s_K0s_US_LS_flag;

  //vectors for mixed-event
  //same event K0s-K0s
  //signal+ background
  vector<TLorentzVector> K0s_K0s_K0s1_vector_ME_SE;
  vector<int> K0s_K0s_K0s1_pT_bin_vector_ME_SE;
  vector<int> K0s_K0s_K0s1_eta_bin_vector_ME_SE;
  vector<TLorentzVector> K0s_K0s_pi1_K0s1_vector_ME_SE;
  vector<int> K0s_K0s_pi1_charge_K0s1_ME_SE;
  vector<TLorentzVector> K0s_K0s_pi2_K0s1_vector_ME_SE;

  vector<TLorentzVector> K0s_K0s_K0s2_vector_ME_SE;
  vector<int> K0s_K0s_K0s2_pT_bin_vector_ME_SE;
  vector<int> K0s_K0s_K0s2_eta_bin_vector_ME_SE;
  vector<TLorentzVector> K0s_K0s_pi1_K0s2_vector_ME_SE;
  vector<int> K0s_K0s_pi1_charge_K0s2_ME_SE;
  vector<TLorentzVector> K0s_K0s_pi2_K0s2_vector_ME_SE;

  //background
  vector<TLorentzVector> K0s_K0s_K0s1_vector_ME_SE_back;
  vector<int> K0s_K0s_K0s1_pT_bin_vector_ME_SE_back;
  vector<int> K0s_K0s_K0s1_eta_bin_vector_ME_SE_back;
  vector<TLorentzVector> K0s_K0s_pi1_K0s1_vector_ME_SE_back;
  vector<int> K0s_K0s_pi1_charge_K0s1_ME_SE_back;
  vector<TLorentzVector> K0s_K0s_pi2_K0s1_vector_ME_SE_back;

  vector<TLorentzVector> K0s_K0s_K0s2_vector_ME_SE_back;
  vector<int> K0s_K0s_K0s2_pT_bin_vector_ME_SE_back;
  vector<int> K0s_K0s_K0s2_eta_bin_vector_ME_SE_back;
  vector<TLorentzVector> K0s_K0s_pi1_K0s2_vector_ME_SE_back;
  vector<int> K0s_K0s_pi1_charge_K0s2_ME_SE_back;
  vector<TLorentzVector> K0s_K0s_pi2_K0s2_vector_ME_SE_back;

  //single K0s
  vector<TLorentzVector> K0s_vector_ME;
  vector<int> K0s_pT_bin_vector_ME;
  vector<int> K0s_eta_bin_vector_ME;
  vector<TLorentzVector> pi1_vector_ME;
  vector<int> pi1_charge_ME;
  vector<TLorentzVector> pi2_vector_ME;

  //LS pairs for mixed event (for background correction)
  vector<TLorentzVector> K0s_vector_ME_LS;
  vector<int> K0s_pT_bin_vector_ME_LS;
  vector<int> K0s_eta_bin_vector_ME_LS;
  vector<TLorentzVector> pi1_vector_ME_LS;
  vector<int> pi1_charge_ME_LS;
  vector<TLorentzVector> pi2_vector_ME_LS;


  int nFill_bckg = 0; //for fillibg of US-LS background

  int eventID_last = -1;


  Long64_t nEntries = 0; //total nEntries

  if(ReadMode == 0)
  {
    nEntries = K0s_tree->GetEntries();
    cout<<"nEntries = "<<nEntries<<endl;
  }

  for(Long64_t i = 0; i < nEntries; i++) //Read TTree only in ReadMode = 0
  {
    //if(ReadMode != 0) break;
    if(i%1000000 == 0)
    {
      cout<<i<<endl;
    }

    K0s_tree->GetEntry(i);

    //set Vz cut based on dataset
    //if( year == 2012 && fabs(Vz) > 30 ) continue; //for cross-check with original production with stricter Vz cut
    if( year == 2015 && fabs(Vz) > 30 ) continue; //double-check what Vz cut to use in Run15
    if( year == 2017 && fabs(Vz) > 30 ) continue; //Run17 Vz cut



    if(i == 0 ) //first iteration
    {
      eventID_last = eventId;
    }

    //------------------------------------------------------------------------------------------------------------------

    //cuts
    //one pi matched for Minv
    //if( !cuts(1, pi1_hasTOFinfo, pi2_hasTOFinfo, K0s_y) ) continue;

    //----------------------------------------------------------------------------------------------------------------

/*
    //fill all histograms for all pT and centrality bins
    int pT_bin = -1;

    //find pT bin of K0s
    for(int j = 0; j < nPtBins; j++) //loop over pT bins
    {
      if(K0s_pt > pT_bins[j] && K0s_pt <= pT_bins[j+1])
      {
        pT_bin = j;
        break; //stop after pT bin is found
      }
    }

    //if( pT_bin == -1 ) continue;
*/
    int pT_bin_corr = -1;

    //find pT bin of K0s
    for(int j = 0; j < nPtBins_corr; j++) //loop over pT bins
    {
      if(K0s_pt > pT_bins_corr[j] && K0s_pt <= pT_bins_corr[j+1])
      {
        pT_bin_corr = j;
        break; //stop after pT bin is found
      }
    }

    //if( pT_bin_corr == -1 ) continue;


    //fill all histograms for all eta and centrality bins
    int eta_bin = -1;

    //find eta bin of K0s
    for(int j = 0; j < nEtaBins; j++) //loop over eta bins
    {
      if(K0s_eta > eta_bins[j] && K0s_eta <= eta_bins[j+1])
      {
        eta_bin = j;
        break; //stop after eta bin is found
      }
    }

    //if( eta_bin == -1 ) continue;


    TLorentzVector K0s_Lorentz_vector(1,1,1,1);
    //K0s_Lorentz_vector.SetPtEtaPhiM(K0s_pt, K0s_eta, K0s_phi, K0s_mass);

    TLorentzVector pi1_Lorenz_vector(1,1,1,1);
    pi1_Lorenz_vector.SetPtEtaPhiM(pi1_pt, pi1_eta, pi1_phi, pi_mass_PDG);

    TLorentzVector pi2_Lorenz_vector(1,1,1,1);
    pi2_Lorenz_vector.SetPtEtaPhiM(pi2_pt, pi2_eta, pi2_phi, pi_mass_PDG);

    K0s_Lorentz_vector = pi1_Lorenz_vector + pi2_Lorenz_vector;


    double K0s_y = K0s_Lorentz_vector.Rapidity();
    //double K0s_y = rapidity(K0s_pt, K0s_eta, K0s_mass_PDG);

    //------------------------------------------

    TLorentzVector pi2_p_mass_Lorenz_vector(1,1,1,1);
    pi2_p_mass_Lorenz_vector.SetPtEtaPhiM(pi2_pt, pi2_eta, pi2_phi, p_mass_PDG);

    TLorentzVector K0s_Lorentz_vector_proton_mass = pi1_Lorenz_vector + pi2_p_mass_Lorenz_vector;


    if(eventId == eventID_last) //same event as in previous iteration and first event
    {
      if(charge == 0 ) //like-sign combinations
      {
        if( cuts(year, K0s_y, pi1_hasTOFinfo, pi1_dca, pi2_dca, K0s_DCAdaughters, K0s_decayL, cos(K0s_theta)) && pT_bin_corr != -1 && eta_bin != -1 &&
           ( cutType != 1 || cuts_topo_sys_err(pi1_dca, pi2_dca, K0s_DCAdaughters, K0s_decayL, cos(K0s_theta) ) ) && //cuts_topo_sys_err is evaluated only when cutType != 1
           ( cutType != 2 || cuts_pT_sys_err(pi1_pt, pi2_pt) ) )
        {
          K0s_vector.push_back(K0s_Lorentz_vector);
          K0s_vector_p_mass.push_back(K0s_Lorentz_vector_proton_mass);
          K0s_pT_bin_vector.push_back(pT_bin_corr);
          K0s_eta_bin_vector.push_back(eta_bin);

          pi1_vector.push_back(pi1_Lorenz_vector);
          pi1_tag_vector.push_back(pi1_InEventID);
          pi1_charge.push_back(pi1_ch);

          pi2_vector.push_back(pi2_Lorenz_vector);
          pi2_tag_vector.push_back(pi2_InEventID);
        }
      }
      else if( charge == 1 )//unlike-sign
      {
        if( cuts(year, K0s_y, pi1_hasTOFinfo, pi1_dca, pi2_dca, K0s_DCAdaughters, K0s_decayL, cos(K0s_theta)) && pT_bin_corr != -1 && eta_bin != -1 &&
           ( cutType != 1 || cuts_topo_sys_err(pi1_dca, pi2_dca, K0s_DCAdaughters, K0s_decayL, cos(K0s_theta) ) ) && //cuts_topo_sys_err is evaluated only when cutType != 1
           ( cutType != 2 || cuts_pT_sys_err(pi1_pt, pi2_pt) ) )
        {
          K0s_vector_background.push_back(K0s_Lorentz_vector);
          K0s_pT_bin_vector_background.push_back(pT_bin_corr);
          K0s_eta_bin_vector_background.push_back(eta_bin);

          pi1_vector_background.push_back(pi1_Lorenz_vector);
          pi1_tag_vector_background.push_back(pi1_InEventID);
          pi1_charge_background.push_back(pi1_ch);

          pi2_vector_background.push_back(pi2_Lorenz_vector);
          pi2_tag_vector_background.push_back(pi2_InEventID);
        }
      }

    }
    else //new event
    {
      eventID_last = eventId;

      //at least one L0-L0 pair in event
      //US-US charge combination
      if(K0s_vector.size() > 1 )
      {
        for(unsigned int iK0s1 = 0; iK0s1 < K0s_vector.size(); iK0s1++)
        {
          for(unsigned int iK0s2 = iK0s1+1; iK0s2 < K0s_vector.size(); iK0s2++)
          {
            //check that no daughters are shared by the K0s in the pair
            if(pi1_tag_vector.at(iK0s1) == pi1_tag_vector.at(iK0s2)) continue;
            if(pi2_tag_vector.at(iK0s1) == pi2_tag_vector.at(iK0s2)) continue;
            if(pi1_tag_vector.at(iK0s1) == pi2_tag_vector.at(iK0s2)) continue;
            if(pi2_tag_vector.at(iK0s1) == pi1_tag_vector.at(iK0s2)) continue;

            //if( (K0s_vector.at(iK0s1).Phi() - K0s_vector.at(iK0s2).Phi()) < 0.01 ) continue;

            K0s1_inv_mass_vs_K0s2_inv_mass_US[K0s_pT_bin_vector.at(iK0s1)][K0s_pT_bin_vector.at(iK0s2)]->Fill(K0s_vector.at(iK0s1).M(), K0s_vector.at(iK0s2).M());

            K0s1_inv_mass_vs_K0s2_inv_mass_US_proton_mass_one[K0s_pT_bin_vector.at(iK0s1)][K0s_pT_bin_vector.at(iK0s2)]->Fill(K0s_vector.at(iK0s1).M(), K0s_vector_p_mass.at(iK0s2).M());
            K0s1_inv_mass_vs_K0s2_inv_mass_US_proton_mass_both[K0s_pT_bin_vector.at(iK0s1)][K0s_pT_bin_vector.at(iK0s2)]->Fill(K0s_vector_p_mass.at(iK0s1).M(), K0s_vector_p_mass.at(iK0s2).M());

            //-----------------------------------------------

            K0s_K0s_pi_eta1_vs_eta2_US_hist->Fill(pi1_vector.at(iK0s1).Eta(), pi1_vector.at(iK0s2).Eta());
            K0s_K0s_pi_phi1_vs_phi2_US_hist->Fill(pi1_vector.at(iK0s1).Phi(), pi1_vector.at(iK0s2).Phi());
            K0s_K0s_pi_pT1_vs_pT2_US_hist->Fill(pi1_vector.at(iK0s1).Pt(), pi1_vector.at(iK0s2).Pt());


            float delta_phi = 0;

            //if( fabs(L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Phi() - L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
            if( fabs(pi1_vector.at(iK0s1).Phi() - pi1_vector.at(iK0s2).Phi()) <= TMath::Pi() )
            {
              //delta_phi = fabs(L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Phi() - L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).Phi());
              delta_phi = fabs(pi1_vector.at(iK0s1).Phi() - pi1_vector.at(iK0s2).Phi());
            }
            else
            {
              //delta_phi = 2*TMath::Pi() - fabs(L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Phi() - L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).Phi()) ;
              delta_phi = 2*TMath::Pi() - fabs(pi1_vector.at(iK0s1).Phi() - pi1_vector.at(iK0s2).Phi());
            }

            K0s_K0s_pi_Delta_phi->Fill( delta_phi );

            K0s_K0s_pi_Delta_eta->Fill(fabs(pi1_vector.at(iK0s1).Eta() - pi1_vector.at(iK0s2).Eta()));

            //match charges of pi used for calculation of theta*, pi1 is the reference
            double K0s_K0s_pairThetaStar = -99;

            if( pi1_charge.at(iK0s1) > 0 && pi1_charge.at(iK0s2) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi1_vector.at(iK0s1), K0s_vector.at(iK0s2), pi1_vector.at(iK0s2));
            if( pi1_charge.at(iK0s1) > 0 && pi1_charge.at(iK0s2) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi1_vector.at(iK0s1), K0s_vector.at(iK0s2), pi2_vector.at(iK0s2));
            if( pi1_charge.at(iK0s1) < 0 && pi1_charge.at(iK0s2) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi1_vector.at(iK0s1), K0s_vector.at(iK0s2), pi2_vector.at(iK0s2));
            if( pi1_charge.at(iK0s1) < 0 && pi1_charge.at(iK0s2) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi1_vector.at(iK0s1), K0s_vector.at(iK0s2), pi1_vector.at(iK0s2));

            //if( pi1_charge.at(iK0s1) < 0 && pi1_charge.at(iK0s2) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi2_vector.at(iK0s1), K0s_vector.at(iK0s2), pi1_vector.at(iK0s2));
            //if( pi1_charge.at(iK0s1) < 0 && pi1_charge.at(iK0s2) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi2_vector.at(iK0s1), K0s_vector.at(iK0s2), pi2_vector.at(iK0s2));

            //double K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi1_vector.at(iK0s1), K0s_vector.at(iK0s2), pi1_vector.at(iK0s2));

            K0s_K0s_K0s1_mom.push_back(K0s_vector.at(iK0s1));
            K0s_K0s_K0s2_mom.push_back(K0s_vector.at(iK0s2));

            K0s_K0s_cos_theta.push_back(TMath::Cos(K0s_K0s_pairThetaStar));

            K0s_K0s_Minv_K0s1.push_back(K0s_vector.at(iK0s1).M());
            K0s_K0s_Minv_K0s2.push_back(K0s_vector.at(iK0s2).M());

            K0s_K0s_pT_bin_K0s1.push_back(K0s_pT_bin_vector.at(iK0s1));
            K0s_K0s_pT_bin_K0s2.push_back(K0s_pT_bin_vector.at(iK0s2));

            K0s_K0s_eta_bin_K0s1.push_back(K0s_eta_bin_vector.at(iK0s1));
            K0s_K0s_eta_bin_K0s2.push_back(K0s_eta_bin_vector.at(iK0s2));

          }
        }
      }


      //US-LS charge combinations
      if(K0s_vector.size() > 0 && K0s_vector_background.size() > 0 )
      {
        for(unsigned int iK0s1 = 0; iK0s1 < K0s_vector.size(); iK0s1++)
        //for(unsigned int iK0s1 = 0; iK0s1 < 1; iK0s1++)
        {
          for(unsigned int iK0s2 = 0; iK0s2 < K0s_vector_background.size(); iK0s2++)
          //for(unsigned int iK0s2 = 0; iK0s2 < 1; iK0s2++)
          {
            //check auto-correlation
            if(pi1_tag_vector.at(iK0s1) == pi1_tag_vector_background.at(iK0s2)) continue;
            if(pi2_tag_vector.at(iK0s1) == pi2_tag_vector_background.at(iK0s2)) continue;
            if(pi1_tag_vector.at(iK0s1) == pi2_tag_vector_background.at(iK0s2)) continue;
            if(pi2_tag_vector.at(iK0s1) == pi1_tag_vector_background.at(iK0s2)) continue;

            //if(pi1_vector.at(iK0s1).Rapidity() == pi1_vector_background.at(iK0s2).Rapidity()) continue;
            //if(pi2_vector.at(iK0s1).Rapidity() == pi2_vector_background.at(iK0s2).Rapidity()) continue;
            //if(pi1_vector.at(iK0s1).Rapidity() == pi2_vector_background.at(iK0s2).Rapidity()) continue;
            //if(pi2_vector.at(iK0s1).Rapidity() == pi1_vector_background.at(iK0s2).Rapidity()) continue;

            //if(fabs(pi1_vector.at(iK0s1).Eta() - pi1_vector_background.at(iK0s2).Eta()) < 1e-5) continue;
            //if(fabs(pi2_vector.at(iK0s1).Eta() - pi2_vector_background.at(iK0s2).Eta()) < 1e-5) continue;
            //if(fabs(pi1_vector.at(iK0s1).Eta() - pi2_vector_background.at(iK0s2).Eta()) < 1e-5) continue;
            //if(fabs(pi2_vector.at(iK0s1).Eta() - pi1_vector_background.at(iK0s2).Eta()) < 1e-5) continue;

            if(nFill_bckg % 2 == 0)//fill US-LS
            {
              K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[K0s_pT_bin_vector.at(iK0s1)][K0s_pT_bin_vector_background.at(iK0s2)]->Fill(K0s_vector.at(iK0s1).M(), K0s_vector_background.at(iK0s2).M());

              //use pion 1 as reference particle
              //match charges of pi used for calculation of theta*, pi1 is the reference
              double K0s_K0s_pairThetaStar = -99;

              if( pi1_charge.at(iK0s1) > 0 && pi1_charge_background.at(iK0s2) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi1_vector.at(iK0s1), K0s_vector_background.at(iK0s2), pi1_vector_background.at(iK0s2));
              if( pi1_charge.at(iK0s1) > 0 && pi1_charge_background.at(iK0s2) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi1_vector.at(iK0s1), K0s_vector_background.at(iK0s2), pi2_vector_background.at(iK0s2));
              if( pi1_charge.at(iK0s1) < 0 && pi1_charge_background.at(iK0s2) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi1_vector.at(iK0s1), K0s_vector_background.at(iK0s2), pi2_vector_background.at(iK0s2));
              if( pi1_charge.at(iK0s1) < 0 && pi1_charge_background.at(iK0s2) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi1_vector.at(iK0s1), K0s_vector_background.at(iK0s2), pi1_vector_background.at(iK0s2));

              //double K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi1_vector.at(iK0s1), K0s_vector_background.at(iK0s2), pi1_vector_background.at(iK0s2));

              K0s_K0s_K0s1_mom_back.push_back(K0s_vector.at(iK0s1));
              K0s_K0s_K0s2_mom_back.push_back(K0s_vector_background.at(iK0s2));

              K0s_K0s_cos_theta_back.push_back(TMath::Cos(K0s_K0s_pairThetaStar));

              K0s_K0s_Minv_K0s1_back.push_back(K0s_vector.at(iK0s1).M());
              K0s_K0s_Minv_K0s2_back.push_back(K0s_vector_background.at(iK0s2).M());

              K0s_K0s_pT_bin_K0s1_back.push_back(K0s_pT_bin_vector.at(iK0s1));
              K0s_K0s_pT_bin_K0s2_back.push_back(K0s_pT_bin_vector_background.at(iK0s2));

              K0s_K0s_eta_bin_K0s1_back.push_back(K0s_eta_bin_vector.at(iK0s1));
              K0s_K0s_eta_bin_K0s2_back.push_back(K0s_eta_bin_vector_background.at(iK0s2));

              K0s_K0s_US_LS_flag.push_back(1);

            }
            else //fill LS-US
            {
              K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[K0s_pT_bin_vector_background.at(iK0s2)][K0s_pT_bin_vector.at(iK0s1)]->Fill(K0s_vector_background.at(iK0s2).M(), K0s_vector.at(iK0s1).M());

              //use pion 1 as reference particle
              //match charges of pi used for calculation of theta*, pi1 is the reference
              double K0s_K0s_pairThetaStar = -99;

              if( pi1_charge.at(iK0s1) > 0 && pi1_charge_background.at(iK0s2) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_background.at(iK0s2), pi1_vector_background.at(iK0s2), K0s_vector.at(iK0s1), pi1_vector.at(iK0s1));
              if( pi1_charge.at(iK0s1) > 0 && pi1_charge_background.at(iK0s2) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_background.at(iK0s2), pi2_vector_background.at(iK0s2), K0s_vector.at(iK0s1), pi1_vector.at(iK0s1));
              if( pi1_charge.at(iK0s1) < 0 && pi1_charge_background.at(iK0s2) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_background.at(iK0s2), pi2_vector_background.at(iK0s2), K0s_vector.at(iK0s1), pi1_vector.at(iK0s1));
              if( pi1_charge.at(iK0s1) < 0 && pi1_charge_background.at(iK0s2) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_background.at(iK0s2), pi1_vector_background.at(iK0s2), K0s_vector.at(iK0s1), pi1_vector.at(iK0s1));

              //double K0s_K0s_pairThetaStar = LpairThetaStar( K0s_vector_background.at(iK0s2), pi1_vector_background.at(iK0s2), K0s_vector.at(iK0s1), pi1_vector.at(iK0s1));

              K0s_K0s_K0s2_mom_back.push_back(K0s_vector.at(iK0s1));
              K0s_K0s_K0s1_mom_back.push_back(K0s_vector_background.at(iK0s2));

              K0s_K0s_cos_theta_back.push_back(TMath::Cos(K0s_K0s_pairThetaStar));

              K0s_K0s_Minv_K0s2_back.push_back(K0s_vector.at(iK0s1).M());
              K0s_K0s_Minv_K0s1_back.push_back(K0s_vector_background.at(iK0s2).M());

              K0s_K0s_pT_bin_K0s2_back.push_back(K0s_pT_bin_vector.at(iK0s1));
              K0s_K0s_pT_bin_K0s1_back.push_back(K0s_pT_bin_vector_background.at(iK0s2));

              K0s_K0s_eta_bin_K0s2_back.push_back(K0s_eta_bin_vector.at(iK0s1));
              K0s_K0s_eta_bin_K0s1_back.push_back(K0s_eta_bin_vector_background.at(iK0s2));

              K0s_K0s_US_LS_flag.push_back(2);
            }

            nFill_bckg++;

          }
        }
      }


      //at least one L0-L0 pair in event
      //LS-LS charge combination
      if(K0s_vector_background.size() > 1)
      {
        for(unsigned int iK0s1 = 0; iK0s1 < K0s_vector_background.size(); iK0s1++)
        //for(unsigned int iK0s1 = 0; iK0s1 < 2; iK0s1++)
        {
          for(unsigned int iK0s2 = iK0s1+1; iK0s2 < K0s_vector_background.size(); iK0s2++)
          //for(unsigned int iK0s2 = iK0s1+1; iK0s2 < 2; iK0s2++)
          {
            if(pi1_vector_background.at(iK0s1).Pt() == pi1_vector_background.at(iK0s2).Pt()) continue;
            if(pi2_vector_background.at(iK0s1).Pt() == pi2_vector_background.at(iK0s2).Pt()) continue;
            if(pi1_vector_background.at(iK0s1).Pt() == pi2_vector_background.at(iK0s2).Pt()) continue;
            if(pi2_vector_background.at(iK0s1).Pt() == pi1_vector_background.at(iK0s2).Pt()) continue;

            K0s1_inv_mass_vs_K0s2_inv_mass_LS[K0s_pT_bin_vector_background.at(iK0s1)][K0s_pT_bin_vector_background.at(iK0s2)]->Fill(K0s_vector_background.at(iK0s1).M(), K0s_vector_background.at(iK0s2).M());

          }
        }
      }
      //_____________________________________________________________________________________________

      //fill vectors for mixed event
      //selecting events with only one L or L-bar
      //need to select more - Minv cut will be applied later

      //same event K0s-K0s
      if( K0s_vector.size() > 1 && K0s_K0s_K0s1_vector_ME_SE.size() < 3e4)
      {
        //check auto-correlation for SE LL pair
        //no daughter can be shared
        if(pi1_tag_vector.at(0) != pi1_tag_vector.at(1) && pi2_tag_vector.at(0) != pi2_tag_vector.at(1) &&
           pi1_tag_vector.at(0) != pi2_tag_vector.at(1) && pi2_tag_vector.at(0) != pi1_tag_vector.at(1) )
        //if(pi1_vector.at(0).Rapidity() != pi1_vector.at(1).Rapidity() && pi2_vector.at(0).Rapidity() != pi2_vector.at(1).Rapidity() &&
        //   pi1_vector.at(0).Rapidity() != pi2_vector.at(1).Rapidity() && pi2_vector.at(0).Rapidity() != pi1_vector.at(1).Rapidity() )
        //if(fabs(pi1_vector.at(0).Eta() - pi1_vector.at(1).Eta()) > 1e-5 && fabs(pi2_vector.at(0).Eta() - pi2_vector.at(1).Eta()) > 1e-5 &&
        //   fabs(pi1_vector.at(0).Eta() - pi2_vector.at(1).Eta()) > 1e-5 && fabs(pi2_vector.at(0).Eta() - pi1_vector.at(1).Eta()) > 1e-5 )
        {
          //Minv pre-selection to boost usable ME statistics
          if( K0s_vector.at(0).M() > 0.475 && K0s_vector.at(0).M() < 0.525 && K0s_vector.at(1).M() > 0.475 && K0s_vector.at(1).M() < 0.525)
          {
            K0s_K0s_K0s1_vector_ME_SE.push_back(K0s_vector.at(0));
            K0s_K0s_pi1_K0s1_vector_ME_SE.push_back(pi1_vector.at(0));
            K0s_K0s_pi1_charge_K0s1_ME_SE.push_back(pi1_charge.at(0));
            K0s_K0s_pi2_K0s1_vector_ME_SE.push_back(pi2_vector.at(0));

            K0s_K0s_K0s1_pT_bin_vector_ME_SE.push_back(K0s_pT_bin_vector.at(0));
            K0s_K0s_K0s1_eta_bin_vector_ME_SE.push_back(K0s_eta_bin_vector.at(0));


            K0s_K0s_K0s2_vector_ME_SE.push_back(K0s_vector.at(1));
            K0s_K0s_pi1_K0s2_vector_ME_SE.push_back(pi1_vector.at(1));
            K0s_K0s_pi1_charge_K0s2_ME_SE.push_back(pi1_charge.at(1));
            K0s_K0s_pi2_K0s2_vector_ME_SE.push_back(pi2_vector.at(1));

            K0s_K0s_K0s2_pT_bin_vector_ME_SE.push_back(K0s_pT_bin_vector.at(1));
            K0s_K0s_K0s2_eta_bin_vector_ME_SE.push_back(K0s_eta_bin_vector.at(1));
          }
        }
      }

      if( K0s_vector.size() == 1 && K0s_vector_ME.size() < 1e5) //orig. 1e5
      {
        if( K0s_vector.at(0).M() > 0.475 && K0s_vector.at(0).M() < 0.525 )
        {
          K0s_vector_ME.push_back(K0s_vector.at(0));
          pi1_vector_ME.push_back(pi1_vector.at(0));
          pi1_charge_ME.push_back(pi1_charge.at(0));
          pi2_vector_ME.push_back(pi2_vector.at(0));

          K0s_pT_bin_vector_ME.push_back(K0s_pT_bin_vector.at(0));
          K0s_eta_bin_vector_ME.push_back(K0s_eta_bin_vector.at(0));
        }
      }

      //background
      if( K0s_vector.size() > 0 && K0s_vector_background.size() > 0 && K0s_K0s_K0s1_vector_ME_SE_back.size() < 3e4 )
      {
        //check auto-correlation for SE LL pair
        //no daughter can be shared
        if(pi1_tag_vector.at(0) != pi1_tag_vector_background.at(0) && pi2_tag_vector.at(0) != pi2_tag_vector_background.at(0) &&
           pi1_tag_vector.at(0) != pi2_tag_vector_background.at(0) && pi2_tag_vector.at(0) != pi1_tag_vector_background.at(0) )
        //if(pi1_vector.at(0).Rapidity() != pi1_vector_background.at(0).Rapidity() && pi2_vector.at(0).Rapidity() != pi2_vector_background.at(0).Rapidity() &&
        //   pi1_vector.at(0).Rapidity() != pi2_vector_background.at(0).Rapidity() && pi2_vector.at(0).Rapidity() != pi1_vector_background.at(0).Rapidity() )
        //if(fabs(pi1_vector.at(0).Eta() - pi1_vector_background.at(0).Eta()) > 1e-5 && fabs(pi2_vector.at(0).Eta() - pi2_vector_background.at(0).Eta()) > 1e-5 &&
        //   fabs(pi1_vector.at(0).Eta() - pi2_vector_background.at(0).Eta()) > 1e-5 && fabs(pi2_vector.at(0).Eta() - pi1_vector_background.at(0).Eta()) > 1e-5 )
        {
          //Minv pre-selection to boost usable ME statistics
          if( K0s_vector.at(0).M() > 0.475 && K0s_vector.at(0).M() < 0.525 && K0s_vector_background.at(0).M() > 0.475 && K0s_vector_background.at(0).M() < 0.525)
          {
            K0s_K0s_K0s1_vector_ME_SE_back.push_back(K0s_vector.at(0));
            K0s_K0s_pi1_K0s1_vector_ME_SE_back.push_back(pi1_vector.at(0));
            K0s_K0s_pi1_charge_K0s1_ME_SE_back.push_back(pi1_charge.at(0));
            K0s_K0s_pi2_K0s1_vector_ME_SE_back.push_back(pi2_vector.at(0));

            K0s_K0s_K0s1_pT_bin_vector_ME_SE_back.push_back(K0s_pT_bin_vector.at(0));
            K0s_K0s_K0s1_eta_bin_vector_ME_SE_back.push_back(K0s_eta_bin_vector.at(0));


            K0s_K0s_K0s2_vector_ME_SE_back.push_back(K0s_vector_background.at(0));
            K0s_K0s_pi1_K0s2_vector_ME_SE_back.push_back(pi1_vector_background.at(0));
            K0s_K0s_pi1_charge_K0s2_ME_SE_back.push_back(pi1_charge_background.at(0));
            K0s_K0s_pi2_K0s2_vector_ME_SE_back.push_back(pi2_vector_background.at(0));

            K0s_K0s_K0s2_pT_bin_vector_ME_SE_back.push_back(K0s_pT_bin_vector_background.at(0));
            K0s_K0s_K0s2_eta_bin_vector_ME_SE_back.push_back(K0s_eta_bin_vector_background.at(0));
          }
        }
      }

      if( K0s_vector_background.size() == 1 && K0s_vector_ME_LS.size() < 1e5) //orig 1e5, 1e4 for testing of SE
      {
        if( K0s_vector_background.at(0).M() > 0.475 && K0s_vector_background.at(0).M() < 0.525)
        {
          K0s_vector_ME_LS.push_back(K0s_vector_background.at(0));
          pi1_vector_ME_LS.push_back(pi1_vector_background.at(0));
          pi1_charge_ME_LS.push_back(pi1_charge_background.at(0));
          pi2_vector_ME_LS.push_back(pi2_vector_background.at(0));

          K0s_pT_bin_vector_ME_LS.push_back(K0s_pT_bin_vector_background.at(0));
          K0s_eta_bin_vector_ME_LS.push_back(K0s_eta_bin_vector_background.at(0));
        }
      }
      //_______________________________________________________________________________________________

      //reset number of K0ss and vectors
      K0s_vector.clear();
      K0s_vector_p_mass.clear();
      K0s_pT_bin_vector.clear();
      K0s_eta_bin_vector.clear();

      pi1_vector.clear();
      pi1_tag_vector.clear();
      pi1_charge.clear();

      pi2_vector.clear();
      pi2_tag_vector.clear();

      //check if we have good K0s in the new event
      if( charge == 0 && cuts(year, K0s_y, pi1_hasTOFinfo, pi1_dca, pi2_dca, K0s_DCAdaughters, K0s_decayL, cos(K0s_theta)) &&  pT_bin_corr != -1 && eta_bin != -1 &&
           ( cutType != 1 || cuts_topo_sys_err(pi1_dca, pi2_dca, K0s_DCAdaughters, K0s_decayL, cos(K0s_theta) ) ) && //cuts_topo_sys_err is evaluated only when cutType != 1
           ( cutType != 2 || cuts_pT_sys_err(pi1_pt, pi2_pt) ) )
      {
        K0s_vector.push_back(K0s_Lorentz_vector);
        K0s_vector_p_mass.push_back(K0s_Lorentz_vector_proton_mass);
        K0s_pT_bin_vector.push_back(pT_bin_corr);
        K0s_eta_bin_vector.push_back(eta_bin);

        pi1_vector.push_back(pi1_Lorenz_vector);
        pi1_tag_vector.push_back(pi1_InEventID);
        pi1_charge.push_back(pi1_ch);

        pi2_vector.push_back(pi2_Lorenz_vector);
        pi2_tag_vector.push_back(pi2_InEventID);
      }


      //reset number of K0ss and vectors
      K0s_vector_background.clear();
      K0s_pT_bin_vector_background.clear();
      K0s_eta_bin_vector_background.clear();

      pi1_vector_background.clear();
      pi1_tag_vector_background.clear();
      pi1_charge_background.clear();

      pi2_vector_background.clear();
      pi2_tag_vector_background.clear();

      if( charge != 0 && cuts(year, K0s_y, pi1_hasTOFinfo, pi1_dca, pi2_dca, K0s_DCAdaughters, K0s_decayL, cos(K0s_theta)) && pT_bin_corr != -1 && eta_bin != -1 &&
           ( cutType != 1 || cuts_topo_sys_err(pi1_dca, pi2_dca, K0s_DCAdaughters, K0s_decayL, cos(K0s_theta) ) ) && //cuts_topo_sys_err is evaluated only when cutType != 1
           ( cutType != 2 || cuts_pT_sys_err(pi1_pt, pi2_pt) ) )
      {
        K0s_vector_background.push_back(K0s_Lorentz_vector);
        K0s_pT_bin_vector_background.push_back(pT_bin_corr);
        K0s_eta_bin_vector_background.push_back(eta_bin);

        pi1_vector_background.push_back(pi1_Lorenz_vector);
        pi1_tag_vector_background.push_back(pi1_InEventID);
        pi1_charge_background.push_back(pi1_ch);

        pi2_vector_background.push_back(pi2_Lorenz_vector);
        pi2_tag_vector_background.push_back(pi2_InEventID);

      }
      //___________________________________________________________

    }//end else new event

  }//end loop over entries in NTuple


  //_____________________________________________________________________________________


  float sideBandScale_K0s_K0s[nPtBins_corr][nPtBins_corr];

  TFitResultPtr fit_res_gaus[nPtBins_corr][nPtBins_corr];

  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      //find bins for side band and for projection to x
      //projection bins just for testing, later can do projections with parameters of 2D fit
      int binLow = K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(0.49);
      int binHigh = K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(0.51);

      int binLow_sideBand = K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(0.52);
      int binHigh_sideBand = K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(0.54);

      //US-LS
      TF2 *doubleGauss_K0s_K0s = new TF2("doubleGauss_K0s_K0s", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 0.45, 0.55, 0.45, 0.55);
      //doubleGauss_K0s_K0s->SetParameters(1200, 0.499, 0.02, 0.499, 0.02);
      if( pTbin1 == 0 && pTbin2 == 0 ) doubleGauss_K0s_K0s->SetParameters(350, 0.499, 0.02, 0.499, 0.02);
      else if( fit_res_gaus[0][0]->IsValid() ) doubleGauss_K0s_K0s->SetParameters(fit_res_gaus[0][0]->Parameter(0), fit_res_gaus[0][0]->Parameter(1), fit_res_gaus[0][0]->Parameter(2), fit_res_gaus[0][0]->Parameter(3), fit_res_gaus[0][0]->Parameter(4));
      else
      {
        cout<<"Fit not valid for L-L pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
        return false;
      }

      TCanvas *K0s1_inv_mass_vs_K0s2_inv_mass_can = new TCanvas(Form("K0s1_inv_mass_vs_K0s2_inv_mass_can_%i_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_vs_K0s2_inv_mass_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      K0s1_inv_mass_vs_K0s2_inv_mass_can->cd();

      //scale LS to match US
      sideBandScale_K0s_K0s[pTbin1][pTbin2] = K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->Integral(binLow_sideBand, binHigh_sideBand, binLow, binHigh)/K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->Integral(binLow_sideBand, binHigh_sideBand, binLow, binHigh);
      //sideBandScale_K0s_K0s[pTbin1][pTbin2] = K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->Integral(binLow_sideBand, binHigh_sideBand, binLow, binHigh_sideBand)/K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow, binHigh_sideBand);

      K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->Sumw2();
      K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->Scale(sideBandScale_K0s_K0s[pTbin1][pTbin2]);

      K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2] = (TH2F*)K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->Clone(Form("K0s1_inv_mass_vs_K0s2_inv_mass_pT1_%i_pT2_%i", pTbin1, pTbin2));
      K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->Add(K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2], -1);
      fit_res_gaus[pTbin1][pTbin2] = K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_K0s_K0s, "s 0", "", 0.49, 0.51);
      //________________________________________________________________________________________________________

    }//end for eta

  }//end for pT

  InvMassFile->Write();
  InvMassFile->Close();


  //return true; //for peak fitting debugging, if needed
  //______________________________________________________________________________________________________________________


  //for extraction of reference polarization from K0s baseline
  const float L0_alpha = 0.732; //decay parameter of L0
  const float L0bar_alpha = -0.758; //decay paramteter of L0bar


  TFile *OutFile; //output file to store production plane histograms

  if(cutType == 0) //create production plane file from nTuple - run in this mode first
  {
    OutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_ana_cuts_work.root", year), "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT
  }
  else if(cutType == 1) //create production plane file from nTuple - run in this mode first
  {
    OutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_tight_topo_cuts_work.root", year), "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT
  }
  else if(cutType == 2) //create production plane file from nTuple - run in this mode first
  {
    OutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_tight_pT_cut_work.root", year), "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT
  }
  else  //read file to store wrong and good sign M_inv distributions and daughter pT
  {
    cout<<"Unable to open file with production plane histograms!"<<endl;

    return false;
  }

  OutFile->cd();

  //QA and ME reweight histograms

  //signal+background
  TH2F *K0s_K0s_eta1_vs_eta2_US_hist = new TH2F("K0s_K0s_eta1_vs_eta2_US_hist", "K0s_K0s_eta1_vs_eta2_US_hist", 20, -1, 1, 20, -1, 1);
  TH2F *K0s_K0s_phi1_vs_phi2_US_hist = new TH2F("K0s_K0s_phi1_vs_phi2_US_hist", "K0s_K0s_phi1_vs_phi2_US_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *K0s_K0s_pT1_vs_pT2_US_hist = new TH2F("K0s_K0s_pT1_vs_pT2_US_hist", "K0s_K0s_pT1_vs_pT2_US_hist", 20, 0, 5, 20, 0, 5);

  TH1F *K0s_K0s_Delta_eta = new TH1F("K0s_K0s_Delta_eta", "K0s_K0s_Delta_eta", 100, 0, 2);
  TH1F *K0s_K0s_Delta_phi = new TH1F("K0s_K0s_Delta_phi", "K0s_K0s_Delta_phi", 100, 0, TMath::Pi());

  //background
  //total - for QA
  TH2F *K0s_K0s_eta1_vs_eta2_US_LS_hist = new TH2F("K0s_K0s_eta1_vs_eta2_US_LS_hist", "K0s_K0s_eta1_vs_eta2_US_LS_hist", 20, -1, 1, 20, -1, 1);
  TH2F *K0s_K0s_phi1_vs_phi2_US_LS_hist = new TH2F("K0s_K0s_phi1_vs_phi2_US_LS_hist", "K0s_K0s_phi1_vs_phi2_US_LS_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *K0s_K0s_pT1_vs_pT2_US_LS_hist = new TH2F("K0s_K0s_pT1_vs_pT2_US_LS_hist", "K0s_K0s_pT1_vs_pT2_US_LS_hist", 20, 0, 5, 20, 0, 5);

  TH1F *K0s_K0s_Delta_LS_eta = new TH1F("K0s_K0s_Delta_LS_eta", "K0s_K0s_Delta_LS_eta", 100, 0, 2);
  TH1F *K0s_K0s_Delta_LS_phi = new TH1F("K0s_K0s_Delta_LS_phi", "K0s_K0s_Delta_LS_phi", 100, 0, TMath::Pi());


  //----------------------------------------------------------------------------------------

  //mixed-event histograms (only for reweight)
  //signal+background
  TH2F *K0s_K0s_eta1_vs_eta2_US_ME_hist = new TH2F("K0s_K0s_eta1_vs_eta2_US_ME_hist", "K0s_K0s_eta1_vs_eta2_US_ME_hist", 20, -1, 1, 20, -1, 1);
  TH2F *K0s_K0s_phi1_vs_phi2_US_ME_hist = new TH2F("K0s_K0s_phi1_vs_phi2_US_ME_hist", "K0s_K0s_phi1_vs_phi2_US_ME_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *K0s_K0s_pT1_vs_pT2_US_ME_hist = new TH2F("K0s_K0s_pT1_vs_pT2_US_ME_hist", "K0s_K0s_pT1_vs_pT2_US_ME_hist", 20, 0, 5, 20, 0, 5);

  TH1F *K0s_K0s_Delta_ME_eta = new TH1F("K0s_K0s_Delta_ME_eta", "K0s_K0s_Delta_ME_eta", 100, 0, 2);
  TH1F *K0s_K0s_Delta_ME_phi = new TH1F("K0s_K0s_Delta_ME_phi", "K0s_K0s_Delta_ME_phi", 100, 0, TMath::Pi());


  TH2F *K0s_K0s_eta1_vs_eta2_US_ME_2_hist = new TH2F("K0s_K0s_eta1_vs_eta2_US_ME_2_hist", "K0s_K0s_eta1_vs_eta2_US_ME_2_hist", 20, -1, 1, 20, -1, 1);
  TH2F *K0s_K0s_phi1_vs_phi2_US_ME_2_hist = new TH2F("K0s_K0s_phi1_vs_phi2_US_ME_2_hist", "K0s_K0s_phi1_vs_phi2_US_ME_2_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *K0s_K0s_pT1_vs_pT2_US_ME_2_hist = new TH2F("K0s_K0s_pT1_vs_pT2_US_ME_2_hist", "K0s_K0s_pT1_vs_pT2_US_ME_2_hist", 20, 0, 5, 20, 0, 5);


  //background
  TH2F *K0s_K0s_eta1_vs_eta2_US_LS_ME_hist = new TH2F("K0s_K0s_eta1_vs_eta2_US_LS_ME_hist", "K0s_K0s_eta1_vs_eta2_US_LS_ME_hist", 20, -1, 1, 20, -1, 1);
  TH2F *K0s_K0s_phi1_vs_phi2_US_LS_ME_hist = new TH2F("K0s_K0s_phi1_vs_phi2_US_LS_ME_hist", "K0s_K0s_phi1_vs_phi2_US_LS_ME_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *K0s_K0s_pT1_vs_pT2_US_LS_ME_hist = new TH2F("K0s_K0s_pT1_vs_pT2_US_LS_ME_hist", "K0s_K0s_pT1_vs_pT2_US_LS_ME_hist", 20, 0, 5, 20, 0, 5);

  TH1F *K0s_K0s_Delta_LS_ME_eta = new TH1F("K0s_K0s_Delta_LS_ME_eta", "K0s_K0s_Delta_LS_ME_eta", 100, 0, 2);
  TH1F *K0s_K0s_Delta_LS_ME_phi = new TH1F("K0s_K0s_Delta_LS_ME_phi", "K0s_K0s_Delta_LS_ME_phi", 100, 0, TMath::Pi());



  //_____________________________________________________________________________

  //data histograms

  TH1F *K0s_K0s_cosThetaProdPlane_US_hist;
  TH1F *K0s_K0s_cosThetaProdPlane_LS_hist;
  TH1F *K0s_K0s_cosThetaProdPlane_ME_hist;
  TH1F *K0s_K0s_cosThetaProdPlane_ME_LS_hist;


  TH1F *K0s_K0s_cosThetaProdPlane_pT_US_hist[nPtBins_corr][nPtBins_corr];
  TH1F *K0s_K0s_cosThetaProdPlane_pT_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1F *K0s_K0s_cosThetaProdPlane_pT_ME_hist[nPtBins_corr][nPtBins_corr];
  TH1F *K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[nPtBins_corr][nPtBins_corr];

  TH1F *K0s_K0s_cosThetaProdPlane_eta_US_hist[nEtaBins][nEtaBins];
  TH1F *K0s_K0s_cosThetaProdPlane_eta_LS_hist[nEtaBins][nEtaBins];
  TH1F *K0s_K0s_cosThetaProdPlane_eta_ME_hist[nEtaBins][nEtaBins];
  TH1F *K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[nEtaBins][nEtaBins];



  K0s_K0s_cosThetaProdPlane_US_hist = new TH1F("K0s_K0s_cosThetaProdPlane_US_hist", "K0s_K0s_cosThetaProdPlane_US_hist", 10, -1, 1);
  K0s_K0s_cosThetaProdPlane_LS_hist = new TH1F("K0s_K0s_cosThetaProdPlane_LS_hist", "K0s_K0s_cosThetaProdPlane_LS_hist", 10, -1, 1);
  K0s_K0s_cosThetaProdPlane_ME_hist = new TH1F("K0s_K0s_cosThetaProdPlane_ME_hist", "K0s_K0s_cosThetaProdPlane_ME_hist", 10, -1, 1);
  K0s_K0s_cosThetaProdPlane_ME_LS_hist = new TH1F("K0s_K0s_cosThetaProdPlane_ME_LS_hist", "K0s_K0s_cosThetaProdPlane_ME_LS_hist", 10, -1, 1);

  //Delta phi bins
  TH2F *K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_hist = new TH2F("K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_hist", "K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_hist", 10, -1, 1, 60, 0, TMath::Pi());
  TH2F *K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist = new TH2F("K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist", "K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist", 10, -1, 1, 60, 0, TMath::Pi());
  TH2F *K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist = new TH2F("K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist", "K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist", 10, -1, 1, 60, 0, TMath::Pi());
  TH2F *K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist = new TH2F("K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist", "K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist", 10, -1, 1, 60, 0, TMath::Pi());

  //Delta y bins
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_US_hist = new TH2F("K0s_K0s_cos_theta_star_vs_delta_eta_US_hist", "K0s_K0s_cos_theta_star_vs_delta_eta_US_hist", 10, -1, 1, 2, 0, 2);
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_US_LS_hist = new TH2F("K0s_K0s_cos_theta_star_vs_delta_eta_US_LS_hist", "K0s_K0s_cos_theta_star_vs_delta_eta_US_LS_hist", 10, -1, 1, 2, 0, 2);
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_US_ME_hist = new TH2F("K0s_K0s_cos_theta_star_vs_delta_eta_US_ME_hist", "K0s_K0s_cos_theta_star_vs_delta_eta_US_ME_hist", 10, -1, 1, 2, 0, 2);
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_US_LS_ME_hist = new TH2F("K0s_K0s_cos_theta_star_vs_delta_eta_US_LS_ME_hist", "K0s_K0s_cos_theta_star_vs_delta_eta_US_LS_ME_hist", 10, -1, 1, 2, 0, 2);

  //Delta y and Delta phi bins
  //Bin 1 is "in cone", bin 2 is out of cone
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_hist = new TH2F("K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_hist", "K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_hist", 10, -1, 1, 2, 0, 2);
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist = new TH2F("K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist", "K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist", 10, -1, 1, 2, 0, 2);
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist = new TH2F("K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist", "K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist", 10, -1, 1, 2, 0, 2);
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist = new TH2F("K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist", "K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist", 10, -1, 1, 2, 0, 2);

  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = new TH1F(Form("K0s_K0s_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = new TH1F(Form("K0s_K0s_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = new TH1F(Form("K0s_K0s_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[pTbin1][pTbin2] = new TH1F(Form("K0s_K0s_cosThetaProdPlane_ME_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_ME_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = new TH1F(Form("K0s_K0s_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = new TH1F(Form("K0s_K0s_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = new TH1F(Form("K0s_K0s_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[etaBin1][etaBin2] = new TH1F(Form("K0s_K0s_cosThetaProdPlane_ME_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_ME_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
    }
  }



  //fill cos(theta*) histograms within the Minv window determined above

  //US combinations (signal + background)
  for(unsigned int K0s_K0s_index = 0; K0s_K0s_index < K0s_K0s_cos_theta.size(); K0s_K0s_index++)
  {
    //if( fabs(K0s_K0s_K0s1_mom.at(K0s_K0s_index).Phi() - K0s_K0s_K0s2_mom.at(K0s_K0s_index).Phi()) < 0.1 ) continue;

    float K0s1_peak_mean = fit_res_gaus[K0s_K0s_pT_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2.at(K0s_K0s_index)]->Parameter(1);
    float K0s1_peak_sigma = fabs(fit_res_gaus[K0s_K0s_pT_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2.at(K0s_K0s_index)]->Parameter(2));

    float K0s2_peak_mean = fit_res_gaus[K0s_K0s_pT_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2.at(K0s_K0s_index)]->Parameter(3);
    float K0s2_peak_sigma = fabs(fit_res_gaus[K0s_K0s_pT_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2.at(K0s_K0s_index)]->Parameter(4));


    if( K0s_K0s_Minv_K0s1.at(K0s_K0s_index) > K0s1_peak_mean-2*K0s1_peak_sigma && K0s_K0s_Minv_K0s1.at(K0s_K0s_index) < K0s1_peak_mean+2*K0s1_peak_sigma &&
        K0s_K0s_Minv_K0s2.at(K0s_K0s_index) > K0s2_peak_mean-2*K0s2_peak_sigma && K0s_K0s_Minv_K0s2.at(K0s_K0s_index) < K0s2_peak_mean+2*K0s2_peak_sigma )
    {
      K0s_K0s_cosThetaProdPlane_US_hist->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index));
      K0s_K0s_cosThetaProdPlane_pT_US_hist[K0s_K0s_pT_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2.at(K0s_K0s_index)]->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index));
      K0s_K0s_cosThetaProdPlane_eta_US_hist[K0s_K0s_eta_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_eta_bin_K0s2.at(K0s_K0s_index)]->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index));

      //fill kinematics histograms for QA and for ME reweight
      K0s_K0s_eta1_vs_eta2_US_hist->Fill(K0s_K0s_K0s1_mom.at(K0s_K0s_index).Eta(), K0s_K0s_K0s2_mom.at(K0s_K0s_index).Eta());
      K0s_K0s_phi1_vs_phi2_US_hist->Fill(K0s_K0s_K0s1_mom.at(K0s_K0s_index).Phi(), K0s_K0s_K0s2_mom.at(K0s_K0s_index).Phi());
      K0s_K0s_pT1_vs_pT2_US_hist->Fill(K0s_K0s_K0s1_mom.at(K0s_K0s_index).Pt(), K0s_K0s_K0s2_mom.at(K0s_K0s_index).Pt());


      float delta_phi = 0;

      if( fabs(K0s_K0s_K0s1_mom.at(K0s_K0s_index).Phi() - K0s_K0s_K0s2_mom.at(K0s_K0s_index).Phi()) <= TMath::Pi() )
      {
        delta_phi = fabs(K0s_K0s_K0s1_mom.at(K0s_K0s_index).Phi() - K0s_K0s_K0s2_mom.at(K0s_K0s_index).Phi());
      }
      else
      {
        delta_phi = 2*TMath::Pi() - fabs(K0s_K0s_K0s1_mom.at(K0s_K0s_index).Phi() - K0s_K0s_K0s2_mom.at(K0s_K0s_index).Phi());
      }

      K0s_K0s_Delta_phi->Fill( delta_phi );

      K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_hist->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index), delta_phi);

      //-------------

      float delta_eta = fabs(fabs(K0s_K0s_K0s1_mom.at(K0s_K0s_index).Rapidity() - K0s_K0s_K0s2_mom.at(K0s_K0s_index).Rapidity()));

      K0s_K0s_Delta_eta->Fill( delta_eta );

      if( delta_eta < 0.5 )
      {
        K0s_K0s_cos_theta_star_vs_delta_eta_US_hist->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index), 0.5);
      }
      else
      {
        K0s_K0s_cos_theta_star_vs_delta_eta_US_hist->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index), 1.5);
      }


      if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
      {
        //fill in-cone
        K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_hist->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index), 0.5);
      }
      else
      {
        //fill out-of-cone
        K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_hist->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index), 1.5);
      }


    }

  }

  //background - US paired with LS
  for(unsigned int K0s_K0s_index = 0; K0s_K0s_index < K0s_K0s_cos_theta_back.size(); K0s_K0s_index++)
  {
    float K0s1_peak_mean = fit_res_gaus[K0s_K0s_pT_bin_K0s1_back.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2_back.at(K0s_K0s_index)]->Parameter(1);
    float K0s1_peak_sigma = fabs(fit_res_gaus[K0s_K0s_pT_bin_K0s1_back.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2_back.at(K0s_K0s_index)]->Parameter(2));

    float K0s2_peak_mean = fit_res_gaus[K0s_K0s_pT_bin_K0s1_back.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2_back.at(K0s_K0s_index)]->Parameter(3);
    float K0s2_peak_sigma = fabs(fit_res_gaus[K0s_K0s_pT_bin_K0s1_back.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2_back.at(K0s_K0s_index)]->Parameter(4));


    //weight for LS - need to scale LS to match realistic sig/bckg ratio known from Minv
    float LS_weight = sideBandScale_K0s_K0s[K0s_K0s_pT_bin_K0s1_back.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2_back.at(K0s_K0s_index)];


    if( K0s_K0s_Minv_K0s1_back.at(K0s_K0s_index) > K0s1_peak_mean-2*K0s1_peak_sigma && K0s_K0s_Minv_K0s1_back.at(K0s_K0s_index) < K0s1_peak_mean+2*K0s1_peak_sigma &&
        K0s_K0s_Minv_K0s2_back.at(K0s_K0s_index) > K0s2_peak_mean-2*K0s2_peak_sigma && K0s_K0s_Minv_K0s2_back.at(K0s_K0s_index) < K0s2_peak_mean+2*K0s2_peak_sigma)
    {
      K0s_K0s_cosThetaProdPlane_LS_hist->Fill(K0s_K0s_cos_theta_back.at(K0s_K0s_index), LS_weight);
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[K0s_K0s_pT_bin_K0s1_back.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2_back.at(K0s_K0s_index)]->Fill(K0s_K0s_cos_theta_back.at(K0s_K0s_index), LS_weight);
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[K0s_K0s_eta_bin_K0s1_back.at(K0s_K0s_index)][K0s_K0s_eta_bin_K0s2_back.at(K0s_K0s_index)]->Fill(K0s_K0s_cos_theta_back.at(K0s_K0s_index), LS_weight);

      //fill kinematics histograms for QA and for ME reweight
      K0s_K0s_eta1_vs_eta2_US_LS_hist->Fill(K0s_K0s_K0s1_mom_back.at(K0s_K0s_index).Eta(), K0s_K0s_K0s2_mom_back.at(K0s_K0s_index).Eta(), LS_weight);
      K0s_K0s_phi1_vs_phi2_US_LS_hist->Fill(K0s_K0s_K0s1_mom_back.at(K0s_K0s_index).Phi(), K0s_K0s_K0s2_mom_back.at(K0s_K0s_index).Phi(), LS_weight);
      K0s_K0s_pT1_vs_pT2_US_LS_hist->Fill(K0s_K0s_K0s1_mom_back.at(K0s_K0s_index).Pt(), K0s_K0s_K0s2_mom_back.at(K0s_K0s_index).Pt(), LS_weight);


      float delta_phi = 0;

      if( fabs(K0s_K0s_K0s1_mom_back.at(K0s_K0s_index).Phi() - K0s_K0s_K0s2_mom_back.at(K0s_K0s_index).Phi()) <= TMath::Pi() )
      {
        delta_phi = fabs(K0s_K0s_K0s1_mom_back.at(K0s_K0s_index).Phi() - K0s_K0s_K0s2_mom_back.at(K0s_K0s_index).Phi());
      }
      else
      {
        delta_phi = 2*TMath::Pi() - fabs(K0s_K0s_K0s1_mom_back.at(K0s_K0s_index).Phi() - K0s_K0s_K0s2_mom_back.at(K0s_K0s_index).Phi());
      }

      K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist->Fill(K0s_K0s_cos_theta_back.at(K0s_K0s_index), delta_phi, LS_weight);


      float delta_eta = fabs(fabs(K0s_K0s_K0s1_mom_back.at(K0s_K0s_index).Rapidity() - K0s_K0s_K0s2_mom_back.at(K0s_K0s_index).Rapidity()));

      K0s_K0s_Delta_LS_eta->Fill( delta_eta );
      K0s_K0s_Delta_LS_phi->Fill( delta_phi );

      if( delta_eta < 0.5 )
      {
        K0s_K0s_cos_theta_star_vs_delta_eta_US_LS_hist->Fill(K0s_K0s_cos_theta_back.at(K0s_K0s_index), 0.5, LS_weight);
      }
      else
      {
        K0s_K0s_cos_theta_star_vs_delta_eta_US_LS_hist->Fill(K0s_K0s_cos_theta_back.at(K0s_K0s_index), 1.5, LS_weight);
      }


      if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
      {
        //fill in-cone
        K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist->Fill(K0s_K0s_cos_theta_back.at(K0s_K0s_index), 0.5, LS_weight);
      }
      else
      {
        //fill out-of-cone
        K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist->Fill(K0s_K0s_cos_theta_back.at(K0s_K0s_index), 1.5, LS_weight);
      }

    }

  }

  //__________________________________________________________________________________________________________________________________

  //mixed event
  vector<float> n_fill_K0sK0s_weight;

  for(unsigned int iK0s_ME_1 = 0; iK0s_ME_1 < K0s_K0s_K0s1_vector_ME_SE.size(); iK0s_ME_1++)
  {
    int fill_ME = 0;

    for(unsigned int iK0s_ME_2 = 0; iK0s_ME_2 < K0s_vector_ME.size(); iK0s_ME_2++)
    {
      float delta_phi_ME = 0;

      //if( fabs(L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Phi() - L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
      if( fabs( K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Phi() - K0s_vector_ME.at(iK0s_ME_2).Phi() ) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Phi() - K0s_vector_ME.at(iK0s_ME_2).Phi() );
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Phi() - K0s_vector_ME.at(iK0s_ME_2).Phi() );
      }

      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Eta() - K0s_vector_ME.at(iK0s_ME_2).Eta() ) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Pt() - K0s_vector_ME.at(iK0s_ME_2).Pt() ) > 0.1 ) continue;

      fill_ME++;

    }

    if( fill_ME > 0) n_fill_K0sK0s_weight.push_back(1./fill_ME);
    else n_fill_K0sK0s_weight.push_back(0);

  }


  int nFill_ME = 0;

  for(unsigned int iK0s_ME_1 = 0; iK0s_ME_1 < K0s_K0s_K0s1_vector_ME_SE.size(); iK0s_ME_1++)
  {
    if( n_fill_K0sK0s_weight.at(iK0s_ME_1) == 0 ) continue;

    for(unsigned int iK0s_ME_2 = 0; iK0s_ME_2 < K0s_vector_ME.size(); iK0s_ME_2++)
    {

      float delta_phi_ME = 0;

      //if( fabs(L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Phi() - L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
      if( fabs( K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Phi() - K0s_vector_ME.at(iK0s_ME_2).Phi() ) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Phi() - K0s_vector_ME.at(iK0s_ME_2).Phi() );
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Phi() - K0s_vector_ME.at(iK0s_ME_2).Phi() );
      }

      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Eta() - K0s_vector_ME.at(iK0s_ME_2).Eta() ) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Pt() - K0s_vector_ME.at(iK0s_ME_2).Pt() ) > 0.1 ) continue;

      //if(nFill_ME % 2 == 0)
      //{
        float K0s1_peak_mean = fit_res_gaus[K0s_K0s_K0s1_pT_bin_vector_ME_SE.at(iK0s_ME_1)][K0s_pT_bin_vector_ME.at(iK0s_ME_2)]->Parameter(1);
        float K0s1_peak_sigma = fabs(fit_res_gaus[K0s_K0s_K0s1_pT_bin_vector_ME_SE.at(iK0s_ME_1)][K0s_pT_bin_vector_ME.at(iK0s_ME_2)]->Parameter(2));

        float K0s2_peak_mean = fit_res_gaus[K0s_K0s_K0s1_pT_bin_vector_ME_SE.at(iK0s_ME_1)][K0s_pT_bin_vector_ME.at(iK0s_ME_2)]->Parameter(3);
        float K0s2_peak_sigma = fabs(fit_res_gaus[K0s_K0s_K0s1_pT_bin_vector_ME_SE.at(iK0s_ME_1)][K0s_pT_bin_vector_ME.at(iK0s_ME_2)]->Parameter(4));

        if( K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).M() < K0s1_peak_mean-2*K0s1_peak_sigma || K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).M() > K0s1_peak_mean+2*K0s1_peak_sigma ) continue;
        if( K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).M() < K0s2_peak_mean-2*K0s2_peak_sigma || K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).M() > K0s2_peak_mean+2*K0s2_peak_sigma ) continue;
        if( K0s_vector_ME.at(iK0s_ME_2).M() < K0s2_peak_mean-2*K0s2_peak_sigma || K0s_vector_ME.at(iK0s_ME_2).M() > K0s2_peak_mean+2*K0s2_peak_sigma ) continue;


        //check what are the charges of the pi in the pair (+- vs. -+)
        //use pi from ME partner with the same charge as the first pi from K0s2 from the SE pair
        double K0s_K0s_pairThetaStar = -99;

        if( K0s_K0s_pi1_charge_K0s1_ME_SE.at(iK0s_ME_1) > 0 && pi1_charge_ME.at(iK0s_ME_2) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_K0s_pi1_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_vector_ME.at(iK0s_ME_2), pi1_vector_ME.at(iK0s_ME_2));
        if( K0s_K0s_pi1_charge_K0s1_ME_SE.at(iK0s_ME_1) > 0 && pi1_charge_ME.at(iK0s_ME_2) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_K0s_pi1_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_vector_ME.at(iK0s_ME_2), pi2_vector_ME.at(iK0s_ME_2));
        if( K0s_K0s_pi1_charge_K0s1_ME_SE.at(iK0s_ME_1) < 0 && pi1_charge_ME.at(iK0s_ME_2) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_K0s_pi1_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_vector_ME.at(iK0s_ME_2), pi2_vector_ME.at(iK0s_ME_2));
        if( K0s_K0s_pi1_charge_K0s1_ME_SE.at(iK0s_ME_1) < 0 && pi1_charge_ME.at(iK0s_ME_2) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_K0s_pi1_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_vector_ME.at(iK0s_ME_2), pi1_vector_ME.at(iK0s_ME_2));

        //double K0s_K0s_pairThetaStar = LpairThetaStar(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_K0s_pi1_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_vector_ME.at(iK0s_ME_2), pi1_vector_ME.at(iK0s_ME_2));


        K0s_K0s_cosThetaProdPlane_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), n_fill_K0sK0s_weight.at(iK0s_ME_1));
        K0s_K0s_cosThetaProdPlane_pT_ME_hist[K0s_K0s_K0s1_pT_bin_vector_ME_SE.at(iK0s_ME_1)][K0s_pT_bin_vector_ME.at(iK0s_ME_2)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar), n_fill_K0sK0s_weight.at(iK0s_ME_1));
        K0s_K0s_cosThetaProdPlane_eta_ME_hist[K0s_K0s_K0s1_eta_bin_vector_ME_SE.at(iK0s_ME_1)][K0s_eta_bin_vector_ME.at(iK0s_ME_2)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar), n_fill_K0sK0s_weight.at(iK0s_ME_1));

        //final ME pair
        K0s_K0s_eta1_vs_eta2_US_ME_hist->Fill(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Eta(), K0s_vector_ME.at(iK0s_ME_2).Eta(), n_fill_K0sK0s_weight.at(iK0s_ME_1));
        K0s_K0s_phi1_vs_phi2_US_ME_hist->Fill(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Phi(), K0s_vector_ME.at(iK0s_ME_2).Phi(), n_fill_K0sK0s_weight.at(iK0s_ME_1));
        K0s_K0s_pT1_vs_pT2_US_ME_hist->Fill(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Pt(), K0s_vector_ME.at(iK0s_ME_2).Pt(), n_fill_K0sK0s_weight.at(iK0s_ME_1));

        //the SE pair used to create ME pair
        K0s_K0s_eta1_vs_eta2_US_ME_2_hist->Fill(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Eta(), K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Eta(), n_fill_K0sK0s_weight.at(iK0s_ME_1));
        K0s_K0s_phi1_vs_phi2_US_ME_2_hist->Fill(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Phi(), K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Phi(), n_fill_K0sK0s_weight.at(iK0s_ME_1));
        K0s_K0s_pT1_vs_pT2_US_ME_2_hist->Fill(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Pt(), K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Pt(), n_fill_K0sK0s_weight.at(iK0s_ME_1));


        float delta_phi = 0;

        if( fabs(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Phi() - K0s_vector_ME.at(iK0s_ME_2).Phi()) <= TMath::Pi() )
        {
          delta_phi = fabs(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Phi() - K0s_vector_ME.at(iK0s_ME_2).Phi());
        }
        else
        {
          delta_phi = 2*TMath::Pi() - fabs(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Phi() - K0s_vector_ME.at(iK0s_ME_2).Phi());
        }

        K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), delta_phi, n_fill_K0sK0s_weight.at(iK0s_ME_1));


        float delta_eta = fabs(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Rapidity() - K0s_vector_ME.at(iK0s_ME_2).Rapidity());

        if( delta_eta < 0.5 )
        {
          K0s_K0s_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 0.5, n_fill_K0sK0s_weight.at(iK0s_ME_1));
        }
        else
        {
          K0s_K0s_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 1.5, n_fill_K0sK0s_weight.at(iK0s_ME_1));
        }

        K0s_K0s_Delta_ME_eta->Fill( delta_eta );
        K0s_K0s_Delta_ME_phi->Fill( delta_phi );


        if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
        {
          //fill in-cone
          K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 0.5, n_fill_K0sK0s_weight.at(iK0s_ME_1));
        }
        else
        {
          //fill out-of-cone
          K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 1.5, n_fill_K0sK0s_weight.at(iK0s_ME_1));
        }


        //fill_ME = 1;
/*
        nFill_ME++;

      }
      else
      {
        float K0s1_peak_mean = fit_res_gaus[K0s_pT_bin_vector_ME.at(iK0s_ME_2)][K0s_K0s_K0s1_pT_bin_vector_ME_SE.at(iK0s_ME_1)]->Parameter(1);
        float K0s1_peak_sigma = fabs(fit_res_gaus[K0s_pT_bin_vector_ME.at(iK0s_ME_2)][K0s_K0s_K0s1_pT_bin_vector_ME_SE.at(iK0s_ME_1)]->Parameter(2));

        float K0s2_peak_mean = fit_res_gaus[K0s_pT_bin_vector_ME.at(iK0s_ME_2)][K0s_K0s_K0s1_pT_bin_vector_ME_SE.at(iK0s_ME_1)]->Parameter(3);
        float K0s2_peak_sigma = fabs(fit_res_gaus[K0s_pT_bin_vector_ME.at(iK0s_ME_2)][K0s_K0s_K0s1_pT_bin_vector_ME_SE.at(iK0s_ME_1)]->Parameter(4));

        if( K0s_vector_ME.at(iK0s_ME_2).M() < K0s1_peak_mean-2*K0s1_peak_sigma || K0s_vector_ME.at(iK0s_ME_2).M() > K0s1_peak_mean+2*K0s1_peak_sigma ) continue;
        if( K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).M() < K0s1_peak_mean-2*K0s1_peak_sigma || K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).M() > K0s1_peak_mean+2*K0s1_peak_sigma ) continue;
        if( K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).M() < K0s2_peak_mean-2*K0s2_peak_sigma || K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).M() > K0s2_peak_mean+2*K0s2_peak_sigma ) continue;


        double K0s_K0s_pairThetaStar = -99;

        if( K0s_K0s_pi1_charge_K0s1_ME_SE.at(iK0s_ME_1) > 0 && pi1_charge_ME.at(iK0s_ME_2) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_ME.at(iK0s_ME_2), pi1_vector_ME.at(iK0s_ME_2), K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_K0s_pi1_K0s1_vector_ME_SE.at(iK0s_ME_1));
        if( K0s_K0s_pi1_charge_K0s1_ME_SE.at(iK0s_ME_1) > 0 && pi1_charge_ME.at(iK0s_ME_2) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_ME.at(iK0s_ME_2), pi1_vector_ME.at(iK0s_ME_2), K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_K0s_pi2_K0s1_vector_ME_SE.at(iK0s_ME_1));
        if( K0s_K0s_pi1_charge_K0s1_ME_SE.at(iK0s_ME_1) < 0 && pi1_charge_ME.at(iK0s_ME_2) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_ME.at(iK0s_ME_2), pi1_vector_ME.at(iK0s_ME_2), K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_K0s_pi2_K0s1_vector_ME_SE.at(iK0s_ME_1));
        if( K0s_K0s_pi1_charge_K0s1_ME_SE.at(iK0s_ME_1) < 0 && pi1_charge_ME.at(iK0s_ME_2) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_ME.at(iK0s_ME_2), pi1_vector_ME.at(iK0s_ME_2), K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_K0s_pi1_K0s1_vector_ME_SE.at(iK0s_ME_1));

        //double K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_ME.at(iK0s_ME_2), pi1_vector_ME.at(iK0s_ME_2), K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1), K0s_K0s_pi1_K0s1_vector_ME_SE.at(iK0s_ME_1));

        K0s_K0s_cosThetaProdPlane_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), n_fill_K0sK0s_weight.at(iK0s_ME_1));
        K0s_K0s_cosThetaProdPlane_pT_ME_hist[K0s_pT_bin_vector_ME.at(iK0s_ME_2)][K0s_K0s_K0s1_pT_bin_vector_ME_SE.at(iK0s_ME_1)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar), n_fill_K0sK0s_weight.at(iK0s_ME_1));
        K0s_K0s_cosThetaProdPlane_eta_ME_hist[K0s_eta_bin_vector_ME.at(iK0s_ME_2)][K0s_K0s_K0s1_eta_bin_vector_ME_SE.at(iK0s_ME_1)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar), n_fill_K0sK0s_weight.at(iK0s_ME_1));

        //final ME pair
        K0s_K0s_eta1_vs_eta2_US_ME_hist->Fill(K0s_vector_ME.at(iK0s_ME_2).Eta(), K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Eta(), n_fill_K0sK0s_weight.at(iK0s_ME_1));
        K0s_K0s_phi1_vs_phi2_US_ME_hist->Fill(K0s_vector_ME.at(iK0s_ME_2).Phi(), K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Phi(), n_fill_K0sK0s_weight.at(iK0s_ME_1));
        K0s_K0s_pT1_vs_pT2_US_ME_hist->Fill(K0s_vector_ME.at(iK0s_ME_2).Pt(), K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Pt(), n_fill_K0sK0s_weight.at(iK0s_ME_1));

        //the SE pair used to create ME pair
        K0s_K0s_eta1_vs_eta2_US_ME_2_hist->Fill(K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Eta(), K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Eta(), n_fill_K0sK0s_weight.at(iK0s_ME_1));
        K0s_K0s_phi1_vs_phi2_US_ME_2_hist->Fill(K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Phi(), K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Phi(), n_fill_K0sK0s_weight.at(iK0s_ME_1));
        K0s_K0s_pT1_vs_pT2_US_ME_2_hist->Fill(K0s_K0s_K0s2_vector_ME_SE.at(iK0s_ME_1).Pt(), K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Phi(), n_fill_K0sK0s_weight.at(iK0s_ME_1));


        float delta_phi = 0;

        if( fabs(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Phi() - K0s_vector_ME.at(iK0s_ME_2).Phi()) <= TMath::Pi() )
        {
          delta_phi = fabs(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Phi() - K0s_vector_ME.at(iK0s_ME_2).Phi());
        }
        else
        {
          delta_phi = 2*TMath::Pi() - fabs(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Phi() - K0s_vector_ME.at(iK0s_ME_2).Phi());
        }

        K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), delta_phi, n_fill_K0sK0s_weight.at(iK0s_ME_1));


        float delta_eta = fabs(K0s_K0s_K0s1_vector_ME_SE.at(iK0s_ME_1).Rapidity() - K0s_vector_ME.at(iK0s_ME_2).Rapidity());

        if( delta_eta < 0.5 )
        {
          K0s_K0s_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 0.5, n_fill_K0sK0s_weight.at(iK0s_ME_1));
        }
        else
        {
          K0s_K0s_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 1.5, n_fill_K0sK0s_weight.at(iK0s_ME_1));
        }

        K0s_K0s_Delta_ME_eta->Fill( delta_eta );
        K0s_K0s_Delta_ME_phi->Fill( delta_phi );

        if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
        {
          //fill in-cone
          K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 0.5, n_fill_K0sK0s_weight.at(iK0s_ME_1));
        }
        else
        {
          //fill out-of-cone
          K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 1.5, n_fill_K0sK0s_weight.at(iK0s_ME_1));
        }

        //fill_ME = 1;

        nFill_ME++;

      }
*/
      //nFill_ME++;

    }
  }

  //-------------------------------------------------------------------------------------------------------------------

  //ME for background - US paired with LS
  vector<float> n_fill_K0sK0s_back_weight;

  for(unsigned int iK0s_ME = 0; iK0s_ME < K0s_K0s_K0s1_vector_ME_SE_back.size(); iK0s_ME++)
  {
    int fill_ME = 0;

    for(unsigned int iK0s_ME_LS = 0; iK0s_ME_LS < K0s_vector_ME_LS.size(); iK0s_ME_LS++)
    {
      float delta_phi_ME = 0;

      if( fabs( K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).Phi() - K0s_vector_ME_LS.at(iK0s_ME_LS).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).Phi() - K0s_vector_ME_LS.at(iK0s_ME_LS).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).Phi() - K0s_vector_ME_LS.at(iK0s_ME_LS).Phi());
      }

      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).Eta() - K0s_vector_ME_LS.at(iK0s_ME_LS).Eta()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).Pt() - K0s_vector_ME_LS.at(iK0s_ME_LS).Pt()) > 0.1 ) continue;

      fill_ME++;

    }

    if( fill_ME > 0) n_fill_K0sK0s_back_weight.push_back(1./fill_ME);
    else n_fill_K0sK0s_back_weight.push_back(0);

  }


  int nFill_ME_LS = 0;

  for(unsigned int iK0s_ME = 0; iK0s_ME < K0s_K0s_K0s1_vector_ME_SE_back.size(); iK0s_ME++)
  {
    if(n_fill_K0sK0s_back_weight.at(iK0s_ME) == 0) continue;

    for(unsigned int iK0s_ME_LS = 0; iK0s_ME_LS < K0s_vector_ME_LS.size(); iK0s_ME_LS++)
    {
      float delta_phi_ME = 0;

      if( fabs( K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).Phi() - K0s_vector_ME_LS.at(iK0s_ME_LS).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).Phi() - K0s_vector_ME_LS.at(iK0s_ME_LS).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).Phi() - K0s_vector_ME_LS.at(iK0s_ME_LS).Phi());
      }

      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).Eta() - K0s_vector_ME_LS.at(iK0s_ME_LS).Eta()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).Pt() - K0s_vector_ME_LS.at(iK0s_ME_LS).Pt()) > 0.1 ) continue;

      //fill US-LS and LS-US combinations
      //if( nFill_ME_LS % 2 == 0 )
      //{
        //for US
        float K0s1_peak_mean = fit_res_gaus[K0s_K0s_K0s1_pT_bin_vector_ME_SE_back.at(iK0s_ME)][K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)]->Parameter(1);
        float K0s1_peak_sigma = fabs(fit_res_gaus[K0s_K0s_K0s1_pT_bin_vector_ME_SE_back.at(iK0s_ME)][K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)]->Parameter(2));

        //for LS
        float K0s2_peak_mean = fit_res_gaus[K0s_K0s_K0s1_pT_bin_vector_ME_SE_back.at(iK0s_ME)][K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)]->Parameter(3);
        float K0s2_peak_sigma = fabs(fit_res_gaus[K0s_K0s_K0s1_pT_bin_vector_ME_SE_back.at(iK0s_ME)][K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)]->Parameter(4));

        //US
        if( K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).M() < K0s1_peak_mean-2*K0s1_peak_sigma || K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).M() > K0s1_peak_mean+2*K0s1_peak_sigma ) continue;
        if( K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).M() < K0s2_peak_mean-2*K0s2_peak_sigma || K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).M() > K0s2_peak_mean+2*K0s2_peak_sigma ) continue;
        //LS
        if( K0s_vector_ME_LS.at(iK0s_ME_LS).M() < K0s2_peak_mean-2*K0s2_peak_sigma || K0s_vector_ME_LS.at(iK0s_ME_LS).M() > K0s2_peak_mean+2*K0s2_peak_sigma ) continue;

        //use pi from ME partner with the same charge as the first pi from K0s2 from the SE pair
        double K0s_K0s_pairThetaStar = -99;

        if( K0s_K0s_pi1_charge_K0s1_ME_SE_back.at(iK0s_ME) > 0 && pi1_charge_ME_LS.at(iK0s_ME_LS) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_K0s_pi1_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_vector_ME_LS.at(iK0s_ME_LS), pi1_vector_ME_LS.at(iK0s_ME_LS));
        if( K0s_K0s_pi1_charge_K0s1_ME_SE_back.at(iK0s_ME) > 0 && pi1_charge_ME_LS.at(iK0s_ME_LS) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_K0s_pi1_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_vector_ME_LS.at(iK0s_ME_LS), pi2_vector_ME_LS.at(iK0s_ME_LS));
        if( K0s_K0s_pi1_charge_K0s1_ME_SE_back.at(iK0s_ME) < 0 && pi1_charge_ME_LS.at(iK0s_ME_LS) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_K0s_pi1_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_vector_ME_LS.at(iK0s_ME_LS), pi2_vector_ME_LS.at(iK0s_ME_LS));
        if( K0s_K0s_pi1_charge_K0s1_ME_SE_back.at(iK0s_ME) < 0 && pi1_charge_ME_LS.at(iK0s_ME_LS) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_K0s_pi1_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_vector_ME_LS.at(iK0s_ME_LS), pi1_vector_ME_LS.at(iK0s_ME_LS));

        //double K0s_K0s_pairThetaStar = LpairThetaStar(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_K0s_pi1_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_vector_ME_LS.at(iK0s_ME_LS), pi1_vector_ME_LS.at(iK0s_ME_LS));

        K0s_K0s_cosThetaProdPlane_ME_LS_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), n_fill_K0sK0s_back_weight.at(iK0s_ME));
        K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[K0s_K0s_K0s1_pT_bin_vector_ME_SE_back.at(iK0s_ME)][K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar), n_fill_K0sK0s_back_weight.at(iK0s_ME));
        K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[K0s_K0s_K0s1_eta_bin_vector_ME_SE_back.at(iK0s_ME)][K0s_eta_bin_vector_ME_LS.at(iK0s_ME_LS)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar), n_fill_K0sK0s_back_weight.at(iK0s_ME));


        K0s_K0s_eta1_vs_eta2_US_LS_ME_hist->Fill(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Eta(), K0s_vector_ME_LS.at(iK0s_ME_LS).Eta(), n_fill_K0sK0s_back_weight.at(iK0s_ME));
        K0s_K0s_phi1_vs_phi2_US_LS_ME_hist->Fill(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Phi(), K0s_vector_ME_LS.at(iK0s_ME_LS).Phi(), n_fill_K0sK0s_back_weight.at(iK0s_ME));
        K0s_K0s_pT1_vs_pT2_US_LS_ME_hist->Fill(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Pt(), K0s_vector_ME_LS.at(iK0s_ME_LS).Pt(), n_fill_K0sK0s_back_weight.at(iK0s_ME));


        float delta_phi = 0;

        if( fabs(fabs(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Phi() - K0s_vector_ME_LS.at(iK0s_ME_LS).Phi())) <= TMath::Pi() )
        {
          delta_phi = fabs(fabs(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Phi() - K0s_vector_ME_LS.at(iK0s_ME_LS).Phi()));
        }
        else
        {
          delta_phi = 2*TMath::Pi() - fabs(fabs(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Phi() - K0s_vector_ME_LS.at(iK0s_ME_LS).Phi()));
        }

        K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), delta_phi, n_fill_K0sK0s_back_weight.at(iK0s_ME));


        float delta_eta = fabs(fabs(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Rapidity() - K0s_vector_ME_LS.at(iK0s_ME_LS).Rapidity()));

        //K0s_K0s_Delta_eta->Fill( delta_eta );

        if( delta_eta < 0.5 )
        {
          K0s_K0s_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 0.5, n_fill_K0sK0s_back_weight.at(iK0s_ME));
        }
        else
        {
          K0s_K0s_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 1.5, n_fill_K0sK0s_back_weight.at(iK0s_ME));
        }

        K0s_K0s_Delta_LS_ME_eta->Fill( delta_eta );
        K0s_K0s_Delta_LS_ME_phi->Fill( delta_phi );

        if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
        {
          //fill in-cone
          K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 0.5, n_fill_K0sK0s_back_weight.at(iK0s_ME));
        }
        else
        {
          //fill out-of-cone
          K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 1.5, n_fill_K0sK0s_back_weight.at(iK0s_ME));
        }
/*
        nFill_ME_LS++;

      }
      else //background matched with US (opposite order)
      {
        //for LS
        float K0s1_peak_mean = fit_res_gaus[K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)][K0s_K0s_K0s1_pT_bin_vector_ME_SE_back.at(iK0s_ME)]->Parameter(1);
        float K0s1_peak_sigma = fabs(fit_res_gaus[K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)][K0s_K0s_K0s1_pT_bin_vector_ME_SE_back.at(iK0s_ME)]->Parameter(2));

        //for US
        float K0s2_peak_mean = fit_res_gaus[K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)][K0s_K0s_K0s1_pT_bin_vector_ME_SE_back.at(iK0s_ME)]->Parameter(3);
        float K0s2_peak_sigma = fabs(fit_res_gaus[K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)][K0s_K0s_K0s1_pT_bin_vector_ME_SE_back.at(iK0s_ME)]->Parameter(4));

        //US
        if( K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).M() < K0s2_peak_mean-2*K0s2_peak_sigma || K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).M() > K0s2_peak_mean+2*K0s2_peak_sigma ) continue;
        if( K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).M() < K0s1_peak_mean-2*K0s1_peak_sigma || K0s_K0s_K0s2_vector_ME_SE_back.at(iK0s_ME).M() > K0s1_peak_mean+2*K0s1_peak_sigma ) continue;
        //LS
        if( K0s_vector_ME_LS.at(iK0s_ME_LS).M() < K0s1_peak_mean-2*K0s1_peak_sigma || K0s_vector_ME_LS.at(iK0s_ME_LS).M() > K0s1_peak_mean+2*K0s1_peak_sigma ) continue;

        double K0s_K0s_pairThetaStar = -99;

        if( K0s_K0s_pi1_charge_K0s1_ME_SE_back.at(iK0s_ME) > 0 && pi1_charge_ME_LS.at(iK0s_ME_LS) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_ME_LS.at(iK0s_ME_LS), pi1_vector_ME_LS.at(iK0s_ME_LS), K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_K0s_pi1_K0s1_vector_ME_SE_back.at(iK0s_ME));
        if( K0s_K0s_pi1_charge_K0s1_ME_SE_back.at(iK0s_ME) > 0 && pi1_charge_ME_LS.at(iK0s_ME_LS) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_ME_LS.at(iK0s_ME_LS), pi2_vector_ME_LS.at(iK0s_ME_LS), K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_K0s_pi1_K0s1_vector_ME_SE_back.at(iK0s_ME));
        if( K0s_K0s_pi1_charge_K0s1_ME_SE_back.at(iK0s_ME) < 0 && pi1_charge_ME_LS.at(iK0s_ME_LS) > 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_ME_LS.at(iK0s_ME_LS), pi2_vector_ME_LS.at(iK0s_ME_LS), K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_K0s_pi1_K0s1_vector_ME_SE_back.at(iK0s_ME));
        if( K0s_K0s_pi1_charge_K0s1_ME_SE_back.at(iK0s_ME) < 0 && pi1_charge_ME_LS.at(iK0s_ME_LS) < 0 ) K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_ME_LS.at(iK0s_ME_LS), pi1_vector_ME_LS.at(iK0s_ME_LS), K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_K0s_pi1_K0s1_vector_ME_SE_back.at(iK0s_ME));

        //double K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_ME_LS.at(iK0s_ME_LS), pi1_vector_ME_LS.at(iK0s_ME_LS), K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME), K0s_K0s_pi1_K0s1_vector_ME_SE_back.at(iK0s_ME) );

        K0s_K0s_cosThetaProdPlane_ME_LS_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), n_fill_K0sK0s_back_weight.at(iK0s_ME));
        K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)][K0s_K0s_K0s1_pT_bin_vector_ME_SE_back.at(iK0s_ME)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar), n_fill_K0sK0s_back_weight.at(iK0s_ME));
        K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[K0s_eta_bin_vector_ME_LS.at(iK0s_ME_LS)][K0s_K0s_K0s1_eta_bin_vector_ME_SE_back.at(iK0s_ME)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar), n_fill_K0sK0s_back_weight.at(iK0s_ME));

        K0s_K0s_eta1_vs_eta2_US_LS_ME_hist->Fill(K0s_vector_ME_LS.at(iK0s_ME_LS).Eta(), K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Eta(), n_fill_K0sK0s_back_weight.at(iK0s_ME));
        K0s_K0s_phi1_vs_phi2_US_LS_ME_hist->Fill(K0s_vector_ME_LS.at(iK0s_ME_LS).Phi(), K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Phi(), n_fill_K0sK0s_back_weight.at(iK0s_ME));
        K0s_K0s_pT1_vs_pT2_US_LS_ME_hist->Fill(K0s_vector_ME_LS.at(iK0s_ME_LS).Pt(), K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Pt(), n_fill_K0sK0s_back_weight.at(iK0s_ME));


        float delta_phi = 0;

        if( fabs(fabs(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Phi() - K0s_vector_ME_LS.at(iK0s_ME_LS).Phi())) <= TMath::Pi() )
        {
          delta_phi = fabs(fabs(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Phi() - K0s_vector_ME_LS.at(iK0s_ME_LS).Phi()));
        }
        else
        {
          delta_phi = 2*TMath::Pi() - fabs(fabs(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Phi() - K0s_vector_ME_LS.at(iK0s_ME_LS).Phi()));
        }

        K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), delta_phi, n_fill_K0sK0s_back_weight.at(iK0s_ME));



        float delta_eta = fabs(fabs(K0s_K0s_K0s1_vector_ME_SE_back.at(iK0s_ME).Rapidity() - K0s_vector_ME_LS.at(iK0s_ME_LS).Rapidity()));

        //K0s_K0s_Delta_eta->Fill( delta_eta );

        if( delta_eta < 0.5 )
        {
          K0s_K0s_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 0.5, n_fill_K0sK0s_back_weight.at(iK0s_ME));
        }
        else
        {
          K0s_K0s_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 1.5, n_fill_K0sK0s_back_weight.at(iK0s_ME));
        }

        K0s_K0s_Delta_LS_ME_eta->Fill( delta_eta );
        K0s_K0s_Delta_LS_ME_phi->Fill( delta_phi );


        if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
        {
          //fill in-cone
          K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 0.5, n_fill_K0sK0s_back_weight.at(iK0s_ME));
        }
        else
        {
          //fill out-of-cone
          K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar), 1.5, n_fill_K0sK0s_back_weight.at(iK0s_ME));
        }

        nFill_ME_LS++;
      }
*/
      //nFill_ME_LS++;

    }
  }

  OutFile->Write();
  OutFile->Close();

  //________________________________________________________________________________________


  return true;

}
//__________________________________________________________________________________________________________________________

//ReadMode = 0 - read TTree, ReadMode = 1 - read histograms - First run in ReadMode = 0 to save relevant histograms, then can run in ReadMode = 1 to read just histograms and save time
//cutType = 0 - analysis cuts, 1 - tight topological cuts, 2 - tight daughter pT cut
//energy - collision energy in GeV
void Ana004_K0s_corr_2D(const int ReadMode = 0, const int cutType = 0, const int energy = 510, const int year = 2017)
{
  ifstream fileList;

  if(energy == 510) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run17_MB_noTOF/fileList.list");
  //else if(energy == 200) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF/fileList.list");
  else if(energy == 200)
  {
    //if(year == 2012) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF/fileList.list");
    //if(year == 2012) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF_open_cuts/fileList.list");
    //if(year == 2012) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF_open_cuts_new/fileList.list");
    if(year == 2012) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF_open_cuts_daughter_tag_helix/fileList.list");
    if(year == 2015) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run15_MB/fileList.list");
    if(year == 2016)
    {
      fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run16_AuAu/fileList.list");
    }
  }

  else
  {
    cout<<"Not a valid colliison energy! Abborting!"<<endl;
    return;
  }


  TH1F * hEventStat1;

  TChain *myChain = new TChain("ntp_K0s");

  string fileFromList;

  int iteration = 0;

  while(getline(fileList, fileFromList))
  {
    myChain->Add(fileFromList.c_str());

    //for event stat. histogram
    TFile *inFile = new TFile(fileFromList.c_str(), "read");

    TList *histoList = (TList*)inFile->Get("picoLambdaAnaMaker");

    TH1F * hEventStat1_work = static_cast<TH1F*>(histoList->FindObject("hEventStat1"));

    if(iteration == 0)
    {
      hEventStat1 = (TH1F*)hEventStat1_work->Clone();
    }
    else
    {
      hEventStat1->Add(hEventStat1_work);
    }

    iteration++;

  }
  //____________________________________________________________________


  TCanvas *testCan = new TCanvas("testCan", "testCan", 1200, 1000);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  testCan->cd();
  hEventStat1->GetXaxis()->SetRange(1,6);
  hEventStat1->GetYaxis()->SetMaxDigits(3);
  hEventStat1->SetMinimum(0);
  hEventStat1->Draw();

  testCan->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/EventStat.png");


  bool AnaFinish = DoAnalysis(myChain, ReadMode, cutType, energy , year);

  if(!AnaFinish)
  {
    cout<<"Analysis of invariant spectra ended abnormally. Aborting!"<<endl;

    return;
  }

  cout<<endl;
  cout<<"Nubmer of accepted events: "<<hEventStat1->GetBinContent(6)<<endl;

  fileList.close();

  return;
}
