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
#include"TLine.h"


using namespace std;

//const int nPtBins = 9;
//float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5., 7.};

//const int nPtBins = 8;
//float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};

const int nPtBins = 8;
float const pT_bins[nPtBins+1] = { 0.5, 0.75, 1.,1.5, 2., 2.5, 3., 4., 5.};

const int nPtBins_corr = 2;
float const pT_bins_corr[nPtBins_corr+1] = { 0.5, 1.5, 5.};

const int nEtaBins = 3;
//float const eta_bins[nEtaBins+1] = { -1, -0.2, 0.2, 1 };
float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };

const double pi_mass_PDG = 0.13957039; //p mass on GeV/c^2 from latest PDG
const double p_mass_PDG = 0.93827208816; //p mass on GeV/c^2 from latest PDG
const double L_mass_PDG = 1.115683; //mass in GeV/c^2 from latest PDG

const float L_y_cut = 1;
//0 - hybrid TOF for both daughters, 1 - strict TOF for pions, 2 - strict TOF for both pion and proton
//also check candidate tree - may have TOF requirement in production
//const int strictTOF_cut = 0;
const float L_cos_theta_cut = 0.996;
const float L_decayL_cut = 25;


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

  //TLorentzVector const L1_reverse(-L1->Px(), -L1->Py(), -L1->Pz(), L1->E());
  //p1->Boost(L1_reverse.BoostVector());


  //TLorentzVector const L2_reverse(-L2->Px(), -L2->Py(), -L2->Pz(), L2->E());
  //p2->Boost(L2_reverse.BoostVector());

  p1->Boost(-L1->BoostVector());

  p2->Boost(-L2->BoostVector());


  return p1->Angle(p2->Vect());
}

//calculate theta star for Lmabdda(-bar) pair
double LpairThetaStar(TLorentzVector L1, TLorentzVector p1, TLorentzVector L2, TLorentzVector p2)
{

  //TLorentzVector const L1_reverse(-L1.Px(), -L1.Py(), -L1.Pz(), L1.E());
  //p1.Boost(L1_reverse.BoostVector());


  //TLorentzVector const L2_reverse(-L2.Px(), -L2.Py(), -L2.Pz(), L2.E());
  //p2.Boost(L2_reverse.BoostVector());


  p1.Boost(-L1.BoostVector());

  p2.Boost(-L2.BoostVector());


  return p1.Angle(p2.Vect());
}

double BoostTest(TLorentzVector L, TLorentzVector p, TLorentzVector pi)
{
  //p.Boost(-L.Vect());
  p.Boost(-L.BoostVector());

  //pi.Boost(-L.Vect());
  pi.Boost(-L.BoostVector());

  return (p + pi).Vect().Mag();
}

//bool cuts(int strictTOF_cut, int pi_hasTOFinfo, int p_hasTOFinfo, float L_y)
bool cuts(int year, float L_y, int hasTOFinfo, float p_DCA, float pi_DCA, float pair_DCA, float dec_L, float cos_theta)
{

  if( !( TMath::Abs(L_y) < L_y_cut ) ) return false;
  if( year == 2017 && hasTOFinfo == 0 ) return false; //TOF matched daughter (choose when calling cuts())
  //if( year == 2016 && hasTOFinfo == 0 ) return false; //TOF matched daughter (choose when calling cuts())
  if( year == 2015 && hasTOFinfo == 0 ) return false; //TOF matched daughter (choose when calling cuts())
  //if( year == 2012 && hasTOFinfo == 0 ) return false; //TOF matched daughter (choose when calling cuts())
  //if( strictTOF_cut == 2 && (pi_hasTOFinfo == 0 || p_hasTOFinfo == 0) ) return false; //TOF matched pions and protons
  //if(cos(L_theta) < L_cos_theta_cut) return false;
  //if(L_decayL > L_decayL_cut) return false;

  //if( p_DCA < 0.1 ) return false;
  //if( pi_DCA < 0.3 ) return false;

  //if( pair_DCA > 1. ) return false;

  //if( dec_L < 2 ) return false;

  //if( cos_theta < 0.996 ) return false;


  if(year == 2016)
  {
    //if(hasTOFinfo == 0) return false;
    //if( pair_DCA > 0.01 ) return false;
  }

  return true;

}


//for tight topological cuts for sys. err.
bool cuts_topo_sys_err(float p_DCA, float pi_DCA, float pair_DCA, float dec_L, float cos_theta)
{
  //if( p_DCA < 0.2 ) return false;
  //if( pi_DCA < 0.4 ) return false;

  if( pair_DCA > 0.9 ) return false;

  if( dec_L < 3 ) return false;

  if( cos_theta < 0.997 ) return false;

  return true;
}

//for tight daughter pT cut for sys. err
bool cuts_pT_sys_err(float p_pT, float pi_pT)
{
  //if( p_pT < 0.3 ) return false;
  if( pi_pT < 0.17 ) return false; //0.2 or 0.3 might be too tight for pions

  return true;
}

bool DoAnalysis(TChain *L_tree, const int ReadMode, const int cut_type = 0, const int energy = 510, const int year = 2017)
{
  TFile *InvMassFile; //output file to store invariant mass histograms


  if(cut_type == 0) InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts_work.root", year), "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT
  if(cut_type == 1) InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_tight_topo_cuts_work.root", year), "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT
  if(cut_type == 2) InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_tight_pT_cut_work.root", year), "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT


  //_______________________________________________________________________________________________________________________________________________

  InvMassFile->cd();

  //without any bins
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_all_US = new TH2F("L0_inv_mass_vs_L0bar_inv_mass_all_US", "L0_inv_mass_vs_L0bar_inv_mass_all_US", 180, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_all_US_LS = new TH2F("L0_inv_mass_vs_L0bar_inv_mass_all_US_LS", "L0_inv_mass_vs_L0bar_inv_mass_all_US_LS", 180, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs;
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_all_LS = new TH2F("L0_inv_mass_vs_L0bar_inv_mass_all_LS", "L0_inv_mass_vs_L0bar_inv_mass_all_LS", 180, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs;
  //TH2F *L0_inv_mass_vs_L0bar_inv_mass_all;

  //adding pT bins together
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_all_2_US = new TH2F("L0_inv_mass_vs_L0bar_inv_mass_all_2_US", "L0_inv_mass_vs_L0bar_inv_mass_all_2_US", 180, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs;
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_all_2_US_LS = new TH2F("L0_inv_mass_vs_L0bar_inv_mass_all_2_US_LS", "L0_inv_mass_vs_L0bar_inv_mass_all_2_US_LS", 180, 1, 1.2, 180, 1, 1.2); //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_all_2_LS = new TH2F("L0_inv_mass_vs_L0bar_inv_mass_all_2_LS", "L0_inv_mass_vs_L0bar_inv_mass_all_2_LS", 180, 1, 1.2, 180, 1, 1.2); //for LS-LS Lambda pairs
  //TH2F *L0_inv_mass_vs_L0bar_inv_mass_all_2;

  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_LS[nPtBins_corr][nPtBins_corr]; //for LS-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass[nPtBins_corr][nPtBins_corr];

  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_zoom[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_zoom[nPtBins_corr][nPtBins_corr];

  //------------------------------

  //without any bins
  TH2F *L0_inv_mass_vs_L0_inv_mass_all_US = new TH2F("L0_inv_mass_vs_L0_inv_mass_all_US", "L0_inv_mass_vs_L0_inv_mass_all_US", 180, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_all_US_LS = new TH2F("L0_inv_mass_vs_L0_inv_mass_all_US_LS", "L0_inv_mass_vs_L0_inv_mass_all_US_LS", 180, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs;
  TH2F *L0_inv_mass_vs_L0_inv_mass_all_LS = new TH2F("L0_inv_mass_vs_L0_inv_mass_all_LS", "L0_inv_mass_vs_L0_inv_mass_all_LS", 180, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs;
  TH2F *L0_inv_mass_vs_L0_inv_mass_all;

  //adding pT bins together
  TH2F *L0_inv_mass_vs_L0_inv_mass_all_2_US = new TH2F("L0_inv_mass_vs_L0_inv_mass_all_2_US", "L0_inv_mass_vs_L0_inv_mass_all_2_US", 180, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs;
  TH2F *L0_inv_mass_vs_L0_inv_mass_all_2_US_LS = new TH2F("L0_inv_mass_vs_L0_inv_mass_all_2_US_LS", "L0_inv_mass_vs_L0_inv_mass_all_2_US_LS", 180, 1, 1.2, 180, 1, 1.2); //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_all_2_LS = new TH2F("L0_inv_mass_vs_L0_inv_mass_all_2_LS", "L0_inv_mass_vs_L0_inv_mass_all_2_LS", 180, 1, 1.2, 180, 1, 1.2); //for LS-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_all_2;

  TH2F *L0_inv_mass_vs_L0_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_LS[nPtBins_corr][nPtBins_corr]; //for LS-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass[nPtBins_corr][nPtBins_corr];

  TH2F *L0_inv_mass_vs_L0_inv_mass_US_zoom[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_US_LS_zoom[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_zoom[nPtBins_corr][nPtBins_corr];

  //------------------------------

  //without any bins
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_all_US = new TH2F("L0bar_inv_mass_vs_L0bar_inv_mass_all_US", "L0bar_inv_mass_vs_L0bar_inv_mass_all_US", 180, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS = new TH2F("L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS", "L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS", 180, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs;
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_all_LS = new TH2F("L0bar_inv_mass_vs_L0bar_inv_mass_all_LS", "L0bar_inv_mass_vs_L0bar_inv_mass_all_LS", 180, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs;
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_all;

  //adding pT bins together
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_all_2_US = new TH2F("L0bar_inv_mass_vs_L0bar_inv_mass_all_2_US", "L0bar_inv_mass_vs_L0bar_inv_mass_all_2_US", 180, 1, 1.2, 180, 1, 1.2); //for US-US Lambda pairs;
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_all_2_US_LS = new TH2F("L0bar_inv_mass_vs_L0bar_inv_mass_all_2_US_LS", "L0bar_inv_mass_vs_L0bar_inv_mass_all_2_US_LS", 180, 1, 1.2, 180, 1, 1.2); //for US-LS Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_all_2_LS = new TH2F("L0bar_inv_mass_vs_L0bar_inv_mass_all_2_LS", "L0bar_inv_mass_vs_L0bar_inv_mass_all_2_LS", 180, 1, 1.2, 180, 1, 1.2); //for LS-LS Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_all_2;

  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_LS[nPtBins_corr][nPtBins_corr]; //for LS-LS Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass[nPtBins_corr][nPtBins_corr];

  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_zoom[nPtBins_corr][nPtBins_corr];


  //------------------------------------------------------------------------------

  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_L[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_Lbar[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_both[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs

  TH2F *L0_inv_mass_vs_L0_inv_mass_US_pi_mass_one[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_US_pi_mass_both[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs

  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_pi_mass_one[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_pi_mass_both[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs


  //QA histograms

  //number of L pairs

  TH1F *n_LLbar_per_event_hist = new TH1F("n_LLbar_per_event_hist", "n_LLbar_per_event_hist", 5, 0, 5);
  TH1F *n_L_per_event_hist = new TH1F("n_L_per_event_hist", "n_L_per_event_hist", 5, 0, 5);
  TH1F *n_Lbar_per_event_hist = new TH1F("n_Lbar_per_event_hist", "n_Lbar_per_event_hist", 5, 0, 5);

  //-----------------------------------------------------------
  //for comparison with PYTHIA or embedding
  TH1F *L_pT_US = new TH1F("L_pT_US", "L_pT_US", 100, 0, 5);
  TH1F *L_pT_LS = new TH1F("L_pT_LS", "L_pT_LS", 100, 0, 5);

  TH1F *L_y_US = new TH1F("L_y_US", "L_y_US", 100, -2, 2);
  TH1F *L_y_LS = new TH1F("L_y_LS", "L_y_LS", 100, -2, 2);

  TH1F *L_phi_US = new TH1F("L_phi_US", "L_phi_US", 180, -TMath::Pi(), TMath::Pi());
  TH1F *L_phi_LS = new TH1F("L_phi_LS", "L_phi_LS", 180, -TMath::Pi(), TMath::Pi());


  TH1F *Lbar_pT_US = new TH1F("Lbar_pT_US", "Lbar_pT_US", 100, 0, 5);
  TH1F *Lbar_pT_LS = new TH1F("Lbar_pT_LS", "Lbar_pT_LS", 100, 0, 5);

  TH1F *Lbar_y_US = new TH1F("Lbar_y_US", "Lbar_y_US", 100, -2, 2);
  TH1F *Lbar_y_LS = new TH1F("Lbar_y_LS", "Lbar_y_LS", 100, -2, 2);

  TH1F *Lbar_phi_US = new TH1F("Lbar_phi_US", "Lbar_phi_US", 180, -TMath::Pi(), TMath::Pi());
  TH1F *Lbar_phi_LS = new TH1F("Lbar_phi_LS", "Lbar_phi_LS", 180, -TMath::Pi(), TMath::Pi());

  //-------------------------------------------------------------------------------------

  TH1F *L_p_pT_US = new TH1F("L_p_pT_US", "L_p_pT_US", 100, 0,5);
  TH1F *L_pi_pT_US = new TH1F("L_pi_pT_US", "L_pi_pT_US", 100, 0,2);

  TH1F *L_p_eta_US = new TH1F("L_p_eta_US", "L_p_eta_US",  99, -0.1, 0.1);
  TH1F *L_pi_eta_US = new TH1F("L_pi_eta_US", "L_pi_eta_US",  99, -0.1, 0.1);

  TH1F *L_p_phi_US = new TH1F("L_p_phi_US", "L_p_phi_US", 20, -TMath::Pi(), TMath::Pi());
  TH1F *L_pi_phi_US = new TH1F("L_pi_phi_US", "L_pi_phi_US", 20, -TMath::Pi(), TMath::Pi());


  TH1F *L_p_pT_LS = new TH1F("L_p_pT_LS", "L_p_pT_LS", 100, 0,5);
  TH1F *L_pi_pT_LS = new TH1F("L_pi_pT_LS", "L_pi_pT_LS", 100, 0,2);

  TH1F *L_p_eta_LS = new TH1F("L_p_eta_LS", "L_p_eta_LS",  99, -0.1, 0.1);
  TH1F *L_pi_eta_LS = new TH1F("L_pi_eta_LS", "L_pi_eta_LS",  99, -0.1, 0.1);

  TH1F *L_p_phi_LS = new TH1F("L_p_phi_LS", "L_p_phi_LS", 20, -TMath::Pi(), TMath::Pi());
  TH1F *L_pi_phi_LS = new TH1F("L_pi_phi_LS", "L_pi_phi_LS", 20, -TMath::Pi(), TMath::Pi());

  //-------------------------------------------------------------------------------------

  TH1F *Lbar_p_pT_US = new TH1F("Lbar_p_pT_US", "Lbar_p_pT_US", 100, 0,5);
  TH1F *Lbar_pi_pT_US = new TH1F("Lbar_pi_pT_US", "Lbar_pi_pT_US", 100, 0,2);

  TH1F *Lbar_p_eta_US = new TH1F("Lbar_p_eta_US", "Lbar_p_eta_US",  99, -0.1, 0.1);
  TH1F *Lbar_pi_eta_US = new TH1F("Lbar_pi_eta_US", "Lbar_pi_eta_US",  99, -0.1, 0.1);

  TH1F *Lbar_p_phi_US = new TH1F("Lbar_p_phi_US", "Lbar_p_phi_US", 20, -TMath::Pi(), TMath::Pi());
  TH1F *Lbar_pi_phi_US = new TH1F("Lbar_pi_phi_US", "Lbar_pi_phi_US", 20, -TMath::Pi(), TMath::Pi());


  TH1F *Lbar_p_pT_LS = new TH1F("Lbar_p_pT_LS", "Lbar_p_pT_LS", 100, 0,5);
  TH1F *Lbar_pi_pT_LS = new TH1F("Lbar_pi_pT_LS", "Lbar_pi_pT_LS", 100, 0,2);

  TH1F *Lbar_p_eta_LS = new TH1F("Lbar_p_eta_LS", "Lbar_p_eta_LS",  99, -0.1, 0.1);
  TH1F *Lbar_pi_eta_LS = new TH1F("Lbar_pi_eta_LS", "Lbar_pi_eta_LS",  99, -0.1, 0.1);

  TH1F *Lbar_p_phi_LS = new TH1F("Lbar_p_phi_LS", "Lbar_p_phi_LS", 20, -TMath::Pi(), TMath::Pi());
  TH1F *Lbar_pi_phi_LS = new TH1F("Lbar_pi_phi_LS", "Lbar_pi_phi_LS", 20, -TMath::Pi(), TMath::Pi());

  //----------------------------------------------------------------------------


  TH2F *L_pT_vs_pi_pT_US = new TH2F("L_pT_vs_pi_pT_US", "L_pT_vs_pi_pT_US", 100, 0, 5, 100, 0, 1.5);
  TH2F *L_pT_vs_pi_pT_LS = new TH2F("L_pT_vs_pi_pT_LS", "L_pT_vs_pi_pT_LS", 100, 0, 5, 100, 0, 1.5);

  TH2F *Lbar_pT_vs_pi_pT_US = new TH2F("Lbar_pT_vs_pi_pT_US", "Lbar_pT_vs_pi_pT_US", 100, 0, 5, 100, 0, 1.5);
  TH2F *Lbar_pT_vs_pi_pT_LS = new TH2F("Lbar_pT_vs_pi_pT_LS", "Lbar_pT_vs_pi_pT_LS", 100, 0, 5, 100, 0, 1.5);

  //----------------------------------------------------------------------------

  TH2F *L_phi_Pico_vs_local_US = new TH2F("L_phi_Pico_vs_local_US", "L_phi_Pico_vs_local_US", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L_phi_Pico_vs_local_LS = new TH2F("L_phi_Pico_vs_local_LS", "L_phi_Pico_vs_local_LS", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());

  TH2F *Lbar_phi_Pico_vs_local_US = new TH2F("Lbar_phi_Pico_vs_local_US", "Lbar_phi_Pico_vs_local_US", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *Lbar_phi_Pico_vs_local_LS = new TH2F("Lbar_phi_Pico_vs_local_LS", "Lbar_phi_Pico_vs_local_LS", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());

  //----------------------------------------------------------------------------

  TH1F *L_Minv_US[nPtBins];
  TH1F *L_Minv_LS[nPtBins];

  TH1F *Lbar_Minv_US[nPtBins];
  TH1F *Lbar_Minv_LS[nPtBins];

  //----------------------------------------------------------------------------

  //boost test histograms
  TH1F *BoostTest_L_hist = new TH1F("BoostTest_L_hist", "BoostTest_L_hist", 20, 0, 1);
  TH1F *BoostTest_Lbar_hist = new TH1F("BoostTest_Lbar_hist", "BoostTest_Lbar_hist", 20, 0, 1);


  for(unsigned int pTbin = 0; pTbin < nPtBins; pTbin++)
  {
    L_Minv_US[pTbin] = new TH1F(Form("L_Minv_US_pT_%i", pTbin), Form("L_Minv_US_pT_%i", pTbin), 180, 1, 1.2);
    L_Minv_LS[pTbin] = new TH1F(Form("L_Minv_LS_pT_%i", pTbin), Form("L_Minv_LS_pT_%i", pTbin), 180, 1, 1.2);

    Lbar_Minv_US[pTbin] = new TH1F(Form("Lbar_Minv_US_pT_%i", pTbin), Form("Lbar_Minv_US_pT_%i", pTbin), 180, 1, 1.2);
    Lbar_Minv_LS[pTbin] = new TH1F(Form("Lbar_Minv_LS_pT_%i", pTbin), Form("Lbar_Minv_LS_pT_%i", pTbin), 180, 1, 1.2);
  }


  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      //invariant mass histograms
      //old bins 200
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);

      L0_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2), 60, 1.1, 1.13, 60, 1.1, 1.13);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2), 60, 1.1, 1.13, 60, 1.1, 1.13);


      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
      L0_inv_mass_vs_L0_inv_mass_LS[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);

      L0_inv_mass_vs_L0_inv_mass_US_zoom[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2), 60, 1.1, 1.13, 60, 1.1, 1.13);
      L0_inv_mass_vs_L0_inv_mass_US_LS_zoom[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_LS_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2), 60, 1.1, 1.13, 60, 1.1, 1.13);


      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);

      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2), 60, 1.1, 1.13, 60, 1.1, 1.13);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2), 60, 1.1, 1.13, 60, 1.1, 1.13);

      //______________________________________________________________________________________________________________________________


      L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_L[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_L_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_L_pT1_%i_pT2_%i", pTbin1, pTbin2), 200, 0.45, 0.55, 180, 1, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_Lbar[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_Lbar_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_Lbar_pT1_%i_pT2_%i", pTbin1, pTbin2), 180, 1, 1.2, 200, 0.45, 0.55);
      L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_both[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_both_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_both_pT1_%i_pT2_%i", pTbin1, pTbin2), 200, 0.45, 0.55, 200, 0.45, 0.55);

      L0_inv_mass_vs_L0_inv_mass_US_pi_mass_one[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_pi_mass_one_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_pi_mass_one_pT1_%i_pT2_%i", pTbin1, pTbin2), 180, 1, 1.2, 200, 0.45, 0.55);
      L0_inv_mass_vs_L0_inv_mass_US_pi_mass_both[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_pi_mass_both_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_pi_mass_both_pT1_%i_pT2_%i", pTbin1, pTbin2), 200, 0.45, 0.55, 200, 0.45, 0.55);

      L0bar_inv_mass_vs_L0bar_inv_mass_US_pi_mass_one[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_pi_mass_one_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_pi_mass_one_pT1_%i_pT2_%i", pTbin1, pTbin2), 180, 1, 1.2, 200, 0.45, 0.55);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_pi_mass_both[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_pi_mass_both_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_pi_mass_both_pT1_%i_pT2_%i", pTbin1, pTbin2), 200, 0.45, 0.55, 200, 0.45, 0.55);

    }
  }

  //________________________________________________________________________________________


  int eventID_last = -1; //to store event ID from last L candidate
  int eventID_last_background = -1; //to store event ID from last L candidate

  Long64_t nEntries = 0; //total nEntries

  Int_t charge;
  Float_t L_mass, L_pt, L_eta, L_phi, L_decayL, L_theta, L_DCAdaughters;

  Int_t pi_InEventID, p_InEventID;
  Float_t pi_pt, p_pt;
  Float_t pi_eta, p_eta;
  Float_t pi_phi, p_phi;
  Int_t p_ch;
  Float_t pi_dca, p_dca;
  Int_t pi_hasTOFinfo, p_hasTOFinfo;

  //Float_t thetaProdPlane;

  Float_t Vz;
  Int_t eventId;

  if(ReadMode == 0)
  {
    //new variable names
    //---------------SET BARANCH ADDRESSES------------------------
    L_tree->SetBranchAddress("pair_charge", &charge);
    L_tree->SetBranchAddress("pair_mass", &L_mass);
    L_tree->SetBranchAddress("pair_pt", &L_pt);
    L_tree->SetBranchAddress("pair_eta", &L_eta);
    L_tree->SetBranchAddress("pair_phi", &L_phi);
    L_tree->SetBranchAddress("pair_decayL", &L_decayL);
    L_tree->SetBranchAddress("pair_theta", &L_theta);
    L_tree->SetBranchAddress("pair_DCAdaughters", &L_DCAdaughters);

    //proton is particle 1 in the pair inside the TTree
    L_tree->SetBranchAddress("p1_InEventID", &p_InEventID),
    L_tree->SetBranchAddress("p1_pt", &p_pt);
    L_tree->SetBranchAddress("p1_eta", &p_eta);
    L_tree->SetBranchAddress("p1_phi", &p_phi);
    L_tree->SetBranchAddress("p1_dca", &p_dca);
    L_tree->SetBranchAddress("p1_ch", &p_ch);
    L_tree->SetBranchAddress("p1_hasTOFinfo", &p_hasTOFinfo);

    //pion is particle 2 in the pair inside the TTree
    L_tree->SetBranchAddress("p2_InEventID", &pi_InEventID),
    L_tree->SetBranchAddress("p2_pt", &pi_pt);
    L_tree->SetBranchAddress("p2_eta", &pi_eta);
    L_tree->SetBranchAddress("p2_phi", &pi_phi);
    L_tree->SetBranchAddress("p2_dca", &pi_dca);
    L_tree->SetBranchAddress("p2_hasTOFinfo", &pi_hasTOFinfo);

    //L_tree->SetBranchAddress("thetaProdPlane", &thetaProdPlane);

    L_tree->SetBranchAddress("Vz", &Vz);
    L_tree->SetBranchAddress("eventId", &eventId);

    //--------------------------------------------------------------------------


    nEntries = L_tree->GetEntries();
    cout<<"nEntries = "<<nEntries<<endl;
  }

  //return false; //for testing

  //to store Lorentz vectors for L pair analysis
  //p vectors to calculate cos(theta*) and for auto-correlation check of L-L and Lbar-Lbar
  //pi vectors for auto-correlation check of L-L and Lbar-Lbar

  //unlike-sign (signal+background)
  vector<TLorentzVector> L_vector;
  vector<TLorentzVector> L_vector_pi_mass;
  vector<int> L_pT_bin_vector;
  vector<int> L_eta_bin_vector;

  vector<TLorentzVector> pi_vector;
  vector<int> pi_tag_vector;
  vector<TLorentzVector> p_vector;
  vector<int> p_tag_vector;


  vector<TLorentzVector> Lbar_vector;
  vector<TLorentzVector> Lbar_vector_pi_mass;
  vector<int> Lbar_pT_bin_vector;
  vector<int> Lbar_eta_bin_vector;

  vector<TLorentzVector> piBar_vector;
  vector<int> piBar_tag_vector;
  vector<TLorentzVector> pBar_vector;
  vector<int> pBar_tag_vector;

  //-------------------------------------------------

  //vectors to store cos(theta*) info
  //analyzed later, after Minv window is found
  vector<TLorentzVector> L_Lbar_L_mom;
  vector<TLorentzVector> L_Lbar_Lbar_mom;
  vector<TVector3> L_Lbar_p_mom;
  vector<TVector3> L_Lbar_pi_mom;
  vector<TVector3> L_Lbar_pBar_mom;
  vector<TVector3> L_Lbar_piBar_mom;

  vector<float> L_Lbar_cos_theta;
  vector<float> L_Lbar_Minv_L;
  vector<float> L_Lbar_Minv_Lbar;
  vector<int> L_Lbar_pT_bin_L;
  vector<int> L_Lbar_pT_bin_Lbar;
  vector<int> L_Lbar_eta_bin_L;
  vector<int> L_Lbar_eta_bin_Lbar;

  //-------------------------------

  vector<TLorentzVector> L_L_L1_mom;
  vector<TLorentzVector> L_L_L2_mom;

  vector<TVector3> L_L_p1_mom;
  vector<TVector3> L_L_pi1_mom;
  vector<TVector3> L_L_p2_mom;
  vector<TVector3> L_L_pi2_mom;

  vector<float> L_L_cos_theta;
  vector<float> L_L_Minv_L1;
  vector<float> L_L_Minv_L2;
  vector<int> L_L_pT_bin_L2;
  vector<int> L_L_pT_bin_L1;
  vector<int> L_L_eta_bin_L1;
  vector<int> L_L_eta_bin_L2;

  //-------------------------------

  vector<TLorentzVector> Lbar_Lbar_Lbar1_mom;
  vector<TLorentzVector> Lbar_Lbar_Lbar2_mom;

  vector<TVector3> Lbar_Lbar_pBar1_mom;
  vector<TVector3> Lbar_Lbar_pBar2_mom;
  vector<TVector3> Lbar_Lbar_piBar1_mom;
  vector<TVector3> Lbar_Lbar_piBar2_mom;

  vector<float> Lbar_Lbar_cos_theta;
  vector<float> Lbar_Lbar_Minv_Lbar1;
  vector<float> Lbar_Lbar_Minv_Lbar2;
  vector<int> Lbar_Lbar_pT_bin_Lbar1;
  vector<int> Lbar_Lbar_pT_bin_Lbar2;
  vector<int> Lbar_Lbar_eta_bin_Lbar1;
  vector<int> Lbar_Lbar_eta_bin_Lbar2;

  //______________________________________________________

  //for combinatorial background
  vector<TLorentzVector> L_vector_background;
  vector<int> L_pT_bin_vector_background;
  vector<int> L_eta_bin_vector_background;

  vector<TLorentzVector> pi_vector_background;
  vector<int> pi_tag_vector_background;

  vector<TLorentzVector> p_vector_background;
  vector<int> p_tag_vector_background;


  vector<TLorentzVector> Lbar_vector_background;
  vector<int> Lbar_pT_bin_vector_background;
  vector<int> Lbar_eta_bin_vector_background;

  vector<TLorentzVector> piBar_vector_background;
  vector<int> piBar_tag_vector_background;

  vector<TLorentzVector> pBar_vector_background;
  vector<int> pBar_tag_vector_background;


  //pair vectors
  //L-Lbar
  //all background pairs
  vector<TLorentzVector> L_Lbar_L_mom_back;
  vector<TLorentzVector> L_Lbar_Lbar_mom_back;

  vector<TVector3> L_Lbar_p_mom_back;
  vector<TVector3> L_Lbar_pi_mom_back;
  vector<TVector3> L_Lbar_pBar_mom_back;
  vector<TVector3> L_Lbar_piBar_mom_back;

  vector<float> L_Lbar_cos_theta_back;
  vector<float> L_Lbar_Minv_L_back;
  vector<float> L_Lbar_Minv_Lbar_back;
  vector<int> L_Lbar_pT_bin_L_back;
  vector<int> L_Lbar_pT_bin_Lbar_back;
  vector<int> L_Lbar_eta_bin_L_back;
  vector<int> L_Lbar_eta_bin_Lbar_back;

  //for reweight of ME - has to separate US-LS and LS-US
  //L(US)-Lbar(LS) - 1, L(LS)-Lbar(US) - 2
  //flag mathces index of histogram
  vector<int> L_Lbar_US_LS_flag;

  //----------------------------------------------

  //L-L
  //all background pairs
  vector<TLorentzVector> L_L_L1_mom_back;
  vector<TLorentzVector> L_L_L2_mom_back;

  vector<TVector3> L_L_p1_mom_back;
  vector<TVector3> L_L_pi1_mom_back;
  vector<TVector3> L_L_p2_mom_back;
  vector<TVector3> L_L_pi2_mom_back;

  vector<float> L_L_cos_theta_back;
  vector<float> L_L_Minv_L1_back;
  vector<float> L_L_Minv_L2_back;
  vector<int> L_L_pT_bin_L1_back;
  vector<int> L_L_pT_bin_L2_back;
  vector<int> L_L_eta_bin_L1_back;
  vector<int> L_L_eta_bin_L2_back;

  //for reweight of ME - has to separate US-LS and LS-US
  //L1(US)-L2(LS) - 1, L1(LS)-L2(US) - 2
  //flag mathces index of histogram
  vector<int> L_L_US_LS_flag;

  //----------------------------------------------

  //Lbar-Lbar
  //all background pairs
  vector<TLorentzVector> Lbar_Lbar_Lbar1_mom_back;
  vector<TLorentzVector> Lbar_Lbar_Lbar2_mom_back;

  vector<TVector3> Lbar_Lbar_pBar1_mom_back;
  vector<TVector3> Lbar_Lbar_piBar1_mom_back;
  vector<TVector3> Lbar_Lbar_pBar2_mom_back;
  vector<TVector3> Lbar_Lbar_piBar2_mom_back;

  vector<float> Lbar_Lbar_cos_theta_back;
  vector<float> Lbar_Lbar_Minv_Lbar1_back;
  vector<float> Lbar_Lbar_Minv_Lbar2_back;
  vector<int> Lbar_Lbar_pT_bin_Lbar1_back;
  vector<int> Lbar_Lbar_pT_bin_Lbar2_back;
  vector<int> Lbar_Lbar_eta_bin_Lbar1_back;
  vector<int> Lbar_Lbar_eta_bin_Lbar2_back;

  //for reweight of ME - has to separate US-LS and LS-US
  //Lbar1(US)-Lbar2(LS) - 1, Lbar1(LS)-Lbar2(US) - 2
  //flag mathces index of histogram
  vector<int> Lbar_Lbar_US_LS_flag;

  //______________________________________________________

  //vectors for mixed-event

  //same event L-Lbar - only L or Lbar from SE pair will be used in ME pairs
  //signal+ background
  vector<TLorentzVector> L_Lbar_L_vector_ME_SE;
  vector<int> L_Lbar_L_pT_bin_vector_ME_SE;
  vector<int> L_Lbar_L_eta_bin_vector_ME_SE;
  vector<TLorentzVector> L_Lbar_p_vector_ME_SE;
  vector<TLorentzVector> L_Lbar_pi_vector_ME_SE;

  vector<TLorentzVector> L_Lbar_Lbar_vector_ME_SE;
  vector<int> L_Lbar_Lbar_pT_bin_vector_ME_SE;
  vector<int> L_Lbar_Lbar_eta_bin_vector_ME_SE;
  vector<TLorentzVector> L_Lbar_pBar_vector_ME_SE;
  vector<TLorentzVector> L_Lbar_piBar_vector_ME_SE;


  //background
  //L(US)-Lbar(LS)
  vector<TLorentzVector> L_Lbar_L_vector_ME_SE_back_US_LS;
  vector<int> L_Lbar_L_pT_bin_vector_ME_SE_back_US_LS;
  vector<int> L_Lbar_L_eta_bin_vector_ME_SE_back_US_LS;
  vector<TLorentzVector> L_Lbar_p_vector_ME_SE_back_US_LS;
  vector<TLorentzVector> L_Lbar_pi_vector_ME_SE_back_US_LS;

  vector<TLorentzVector> L_Lbar_Lbar_vector_ME_SE_back_US_LS;
  vector<int> L_Lbar_Lbar_pT_bin_vector_ME_SE_back_US_LS;
  vector<int> L_Lbar_Lbar_eta_bin_vector_ME_SE_back_US_LS;
  vector<TLorentzVector> L_Lbar_pBar_vector_ME_SE_back_US_LS;
  vector<TLorentzVector> L_Lbar_piBar_vector_ME_SE_back_US_LS;

  //L(LS)-Lbar(US)
  vector<TLorentzVector> L_Lbar_L_vector_ME_SE_back_LS_US;
  vector<int> L_Lbar_L_pT_bin_vector_ME_SE_back_LS_US;
  vector<int> L_Lbar_L_eta_bin_vector_ME_SE_back_LS_US;
  vector<TLorentzVector> L_Lbar_p_vector_ME_SE_back_LS_US;
  vector<TLorentzVector> L_Lbar_pi_vector_ME_SE_back_LS_US;

  vector<TLorentzVector> L_Lbar_Lbar_vector_ME_SE_back_LS_US;
  vector<int> L_Lbar_Lbar_pT_bin_vector_ME_SE_back_LS_US;
  vector<int> L_Lbar_Lbar_eta_bin_vector_ME_SE_back_LS_US;
  vector<TLorentzVector> L_Lbar_pBar_vector_ME_SE_back_LS_US;
  vector<TLorentzVector> L_Lbar_piBar_vector_ME_SE_back_LS_US;


  //----------------------------------------------------


  //same event L-L - only L or Lbar from SE pair will be used in ME pairs
  //signal+ background
  vector<TLorentzVector> L_L_L1_vector_ME_SE;
  vector<int> L_L_L1_pT_bin_vector_ME_SE;
  vector<int> L_L_L1_eta_bin_vector_ME_SE;
  vector<TLorentzVector> L_L_p1_vector_ME_SE;
  vector<TLorentzVector> L_L_pi1_vector_ME_SE;

  vector<TLorentzVector> L_L_L2_vector_ME_SE;
  vector<int> L_L_L2_pT_bin_vector_ME_SE;
  vector<int> L_L_L2_eta_bin_vector_ME_SE;
  vector<TLorentzVector> L_L_p2_vector_ME_SE;
  vector<TLorentzVector> L_L_pi2_vector_ME_SE;

  //background
  vector<TLorentzVector> L_L_L1_vector_ME_SE_back;
  vector<int> L_L_L1_pT_bin_vector_ME_SE_back;
  vector<int> L_L_L1_eta_bin_vector_ME_SE_back;
  vector<TLorentzVector> L_L_p1_vector_ME_SE_back;
  vector<TLorentzVector> L_L_pi1_vector_ME_SE_back;

  vector<TLorentzVector> L_L_L2_vector_ME_SE_back;
  vector<int> L_L_L2_pT_bin_vector_ME_SE_back;
  vector<int> L_L_L2_eta_bin_vector_ME_SE_back;
  vector<TLorentzVector> L_L_p2_vector_ME_SE_back;
  vector<TLorentzVector> L_L_pi2_vector_ME_SE_back;

  //----------------------------------------------------


  //same event Lbar-Lbar - onLbary Lbar or Lbarbar from SE pair wiLbarLbar be used in ME pairs
  //signal+ background
  vector<TLorentzVector> Lbar_Lbar_Lbar1_vector_ME_SE;
  vector<int> Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE;
  vector<int> Lbar_Lbar_Lbar1_eta_bin_vector_ME_SE;
  vector<TLorentzVector> Lbar_Lbar_pBar1_vector_ME_SE;
  vector<TLorentzVector> Lbar_Lbar_piBar1_vector_ME_SE;

  vector<TLorentzVector> Lbar_Lbar_Lbar2_vector_ME_SE;
  vector<int> Lbar_Lbar_Lbar2_pT_bin_vector_ME_SE;
  vector<int> Lbar_Lbar_Lbar2_eta_bin_vector_ME_SE;
  vector<TLorentzVector> Lbar_Lbar_pBar2_vector_ME_SE;
  vector<TLorentzVector> Lbar_Lbar_piBar2_vector_ME_SE;

  //background
  vector<TLorentzVector> Lbar_Lbar_Lbar1_vector_ME_SE_back;
  vector<int> Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back;
  vector<int> Lbar_Lbar_Lbar1_eta_bin_vector_ME_SE_back;
  vector<TLorentzVector> Lbar_Lbar_pBar1_vector_ME_SE_back;
  vector<TLorentzVector> Lbar_Lbar_piBar1_vector_ME_SE_back;

  vector<TLorentzVector> Lbar_Lbar_Lbar2_vector_ME_SE_back;
  vector<int> Lbar_Lbar_Lbar2_pT_bin_vector_ME_SE_back;
  vector<int> Lbar_Lbar_Lbar2_eta_bin_vector_ME_SE_back;
  vector<TLorentzVector> Lbar_Lbar_pBar2_vector_ME_SE_back;
  vector<TLorentzVector> Lbar_Lbar_piBar2_vector_ME_SE_back;


  //----------------------------------------------------------------------------------------------------

  //single L and Lbat to be mixed with same event
  vector<TLorentzVector> L_vector_ME;
  vector<int> L_pT_bin_vector_ME;
  vector<int> L_eta_bin_vector_ME;
  vector<TLorentzVector> p_vector_ME;
  vector<TLorentzVector> pi_vector_ME;

  vector<TLorentzVector> Lbar_vector_ME;
  vector<int> Lbar_pT_bin_vector_ME;
  vector<int> Lbar_eta_bin_vector_ME;
  vector<TLorentzVector> pBar_vector_ME;
  vector<TLorentzVector> piBar_vector_ME;

  //background ME vectors
  vector<TLorentzVector> L_vector_ME_LS;
  vector<int> L_pT_bin_vector_ME_LS;
  vector<int> L_eta_bin_vector_ME_LS;
  vector<TLorentzVector> p_vector_ME_LS;
  vector<TLorentzVector> pi_vector_ME_LS;

  vector<TLorentzVector> Lbar_vector_ME_LS;
  vector<int> Lbar_pT_bin_vector_ME_LS;
  vector<int> Lbar_eta_bin_vector_ME_LS;
  vector<TLorentzVector> pBar_vector_ME_LS;
  vector<TLorentzVector> piBar_vector_ME_LS;


  //TLorentzVector *L_Lorentz_vector = new TLorentzVector(1,1,1,1);
  //TLorentzVector *p_Lorenz_vector = new TLorentzVector(1,1,1,1);

  //for filling of US-LS L-L and Lbar-Lbar pair background
  //to be able to fill background Minv histograms correctly
  int nFillLLbar_sig = 0; //for testing

  int nFillLL_sig = 0; //for testing
  int nFillLL_bckg = 0;

  int nFillLbarLbar_sig = 0; //for testing
  int nFillLbarLbar_backg = 0;

  Long64_t nEvents = 0;

  for(Long64_t i = 0; i < nEntries; i++) //Read TTree only in ReadMode = 0
  {
    L_tree->GetEntry(i);

     //if(ReadMode != 0) break;
    if(i%1000000 == 0)
    //if(i%100000 == 0)
    {
      cout<<i<<endl;

      //cout<<L_Lbar_cos_theta.size()<<endl;
    }

    if(i == 0 ) //first iteration
    {
      eventID_last = eventId;
    }

    //set Vz cut based on dataset
    //if( year == 2012 && fabs(Vz) > 30 ) continue; //for cross-check with original production with stricter Vz cut
    //if( year == 2015 && fabs(Vz) > 30 ) continue; //double-check what Vz cut to use in Run15
    if( year == 2015 && fabs(Vz) > 6 ) continue; //double-check what Vz cut to use in Run15
    if( year == 2016 && fabs(Vz) > 5 ) continue; //double-check what Vz cut to use in Run15
    if( year == 2017 && fabs(Vz) > 30 ) continue; //Run17 Vz cut

    //testing of pion pT cut
    //if( pi_pt < 0.21 ) continue;

    //if(p_dca > 1.5 ) continue;
    //if(pi_dca > 1.5 ) continue;

    //double L_xF = fabs(pz(L_pt, L_eta))/energy/2.; //energy is in CMS, need energz of one proton

    //calculate Lambda rapidity y
    double L_y = rapidity(L_pt, L_eta, L_mass_PDG);

    //------------------------------------------------------------------------------------------------------------------

    TLorentzVector L_Lorentz_vector(1,1,1,1);
    //L_Lorentz_vector.SetPtEtaPhiM(L_pt, L_eta, L_phi, L_mass);

    TLorentzVector p_Lorenz_vector(1,1,1,1);
    p_Lorenz_vector.SetPtEtaPhiM(p_pt, p_eta, p_phi, p_mass_PDG);

    TLorentzVector pi_Lorenz_vector(1,1,1,1);
    pi_Lorenz_vector.SetPtEtaPhiM(pi_pt, pi_eta, pi_phi, pi_mass_PDG);

    L_Lorentz_vector = p_Lorenz_vector + pi_Lorenz_vector;
    TLorentzVector L_Lorentz_vector_local = p_Lorenz_vector + pi_Lorenz_vector;

    //---------------------------------------------------------------------------------------

    //find bins
    int pT_bin = -1;

    //find pT bin of Lambda
    for(int j = 0; j < nPtBins; j++) //loop over pT bins
    {
      if(L_pt > pT_bins[j] && L_pt <= pT_bins[j+1])
      {
        pT_bin = j;
        break; //stop after pT bin is found
      }
    }

    int pT_bin_corr = -1;

    //find pT bin of Lambda
    for(int j = 0; j < nPtBins_corr; j++) //loop over pT bins
    {
      if(L_pt > pT_bins_corr[j] && L_pt <= pT_bins_corr[j+1])
      {
        pT_bin_corr = j;
        break; //stop after pT bin is found
      }
    }

    //check of pT bin done later
    //if( pT_bin_corr == -1 ) continue;


    //fill all histograms for all eta and centrality bins
    int eta_bin = -1;

    //find eta bin of Lambda
    for(int j = 0; j < nEtaBins; j++) //loop over eta bins
    {
      //if(L_eta > eta_bins[j] && L_eta <= eta_bins[j+1])
      if(L_y > eta_bins[j] && L_y <= eta_bins[j+1])
      {
        eta_bin = j;
        break; //stop after eta bin is found
      }
    }

    //check of eta bin done later
    //if( eta_bin == -1 ) continue;

    //for comparison with MC
    //daughter kinematics

    if(charge == 0)//US
    {
      if(p_ch == 1)//L
      {
        if( L_mass > 1.112 && L_mass < 1.12)
        {
          L_p_pT_US->Fill(p_pt);
          L_pi_pT_US->Fill(pi_pt);

          L_p_eta_US->Fill(p_eta);
          L_pi_eta_US->Fill(pi_eta);

          L_p_phi_US->Fill(p_phi);
          L_pi_phi_US->Fill(pi_phi);

          L_pT_US->Fill(L_pt);
          L_y_US->Fill(L_y);
          L_phi_US->Fill(L_phi);

          //-------------------------

          L_phi_Pico_vs_local_US->Fill(L_Lorentz_vector.Phi(), L_Lorentz_vector_local.Phi());
        }

        //L_Minv_US->Fill(L_mass);
      }

      if(p_ch == -1)//Lbar
      {
        if( L_mass > 1.112 && L_mass < 1.12)
        {
          Lbar_p_pT_US->Fill(p_pt);
          Lbar_pi_pT_US->Fill(pi_pt);

          Lbar_p_eta_US->Fill(p_eta);
          Lbar_pi_eta_US->Fill(pi_eta);

          Lbar_p_phi_US->Fill(p_phi);
          Lbar_pi_phi_US->Fill(pi_phi);


          Lbar_pT_US->Fill(L_pt);
          Lbar_y_US->Fill(L_y);
          Lbar_phi_US->Fill(L_phi);

          //-------------------------

          Lbar_phi_Pico_vs_local_US->Fill(L_Lorentz_vector.Phi(), L_Lorentz_vector_local.Phi());
        }

        //Lbar_Minv_US->Fill(L_mass);
      }
    }
    if(charge != 0)//LS
    {
      if(p_ch == 1)//L
      {
        if( L_mass > 1.112 && L_mass < 1.12)
        {
          L_p_pT_LS->Fill(p_pt);
          L_pi_pT_LS->Fill(pi_pt);

          L_p_eta_LS->Fill(p_eta);
          L_pi_eta_LS->Fill(pi_eta);

          L_p_phi_LS->Fill(p_phi);
          L_pi_phi_LS->Fill(pi_phi);


          L_pT_LS->Fill(L_pt);
          L_y_LS->Fill(L_y);
          L_phi_LS->Fill(L_phi);

          //-------------------------

          L_phi_Pico_vs_local_LS->Fill(L_Lorentz_vector.Phi(), L_Lorentz_vector_local.Phi());
        }

        //L_Minv_LS->Fill(L_mass);
      }

      if(p_ch == -1)//Lbar
      {
        if( L_mass > 1.112 && L_mass < 1.12)
        {
          Lbar_p_pT_LS->Fill(p_pt);
          Lbar_pi_pT_LS->Fill(pi_pt);

          Lbar_p_eta_LS->Fill(p_eta);
          Lbar_pi_eta_LS->Fill(pi_eta);

          Lbar_p_phi_LS->Fill(p_phi);
          Lbar_pi_phi_LS->Fill(pi_phi);


          Lbar_pT_LS->Fill(L_pt);
          Lbar_y_LS->Fill(L_y);
          Lbar_phi_LS->Fill(L_phi);

          //-------------------------

          Lbar_phi_Pico_vs_local_LS->Fill(L_Lorentz_vector.Phi(), L_Lorentz_vector_local.Phi());
        }

        //Lbar_Minv_LS->Fill(L_mass);
      }
    }


/*
    TLorentzVector *L_Lorentz_vector = new TLorentzVector(1,1,1,1);
    L_Lorentz_vector->SetPtEtaPhiM(L_pt, L_eta, L_phi, L_mass);

    TLorentzVector *p_Lorenz_vector = new TLorentzVector(1,1,1,1);
    p_Lorenz_vector->SetPtEtaPhiM(p_pt, p_eta, p_phi, p_mass_PDG);
*/


    //------------------------------------------------------------------------

    TLorentzVector p_Lorenz_vector_pi_mass(1,1,1,1);
    p_Lorenz_vector_pi_mass.SetPtEtaPhiM(p_pt, p_eta, p_phi, pi_mass_PDG);

    TLorentzVector L_Lorentz_vector_pi_mass = p_Lorenz_vector_pi_mass + pi_Lorenz_vector;

    //------------------------------------------------------------------------

    if(eventId == eventID_last)
    {
      //US charge combinations
      if(charge == 0)
      {
        //cuts
        //first line is evaluated everz time
        if( cuts(year, L_y, pi_hasTOFinfo, p_dca, pi_dca, L_DCAdaughters, L_decayL, cos(L_theta)) && pT_bin_corr != -1 && eta_bin != -1 &&
            ( cut_type != 1 || cuts_topo_sys_err(p_dca, pi_dca, L_DCAdaughters, L_decayL, cos(L_theta) ) ) && //cuts_topo_sys_err is evaluated only when cut_type != 1
            ( cut_type != 2 || cuts_pT_sys_err(p_pt, pi_pt) ) )                                               //cuts_pT_sys_err is evaluated only when cut_type != 2
        //if( cuts(year, L_y, p_hasTOFinfo) && pT_bin_corr != -1 && eta_bin != -1 )
        {
          if( p_ch == 1 )
          {
            //QA histograms
            L_Minv_US[pT_bin]->Fill(L_mass);

            L_pT_vs_pi_pT_US->Fill(L_Lorentz_vector.Pt(), pi_Lorenz_vector.Pt());

            BoostTest_L_hist->Fill(BoostTest(L_Lorentz_vector, p_Lorenz_vector, pi_Lorenz_vector));

            //fill vectors
            L_vector.push_back(L_Lorentz_vector);
            L_vector_pi_mass.push_back(L_Lorentz_vector_pi_mass);
            L_pT_bin_vector.push_back(pT_bin_corr);
            L_eta_bin_vector.push_back(eta_bin);

            pi_vector.push_back(pi_Lorenz_vector);
            pi_tag_vector.push_back(pi_InEventID);

            p_vector.push_back(p_Lorenz_vector);
            p_tag_vector.push_back(p_InEventID);

          }
          else if( p_ch == -1 )
          {
            //QA histograms
            Lbar_Minv_US[pT_bin]->Fill(L_mass);

            Lbar_pT_vs_pi_pT_US->Fill(L_Lorentz_vector.Pt(), pi_Lorenz_vector.Pt());

            BoostTest_Lbar_hist->Fill(BoostTest(L_Lorentz_vector, p_Lorenz_vector, pi_Lorenz_vector));

            //fill vectors
            Lbar_vector.push_back(L_Lorentz_vector);
            Lbar_vector_pi_mass.push_back(L_Lorentz_vector_pi_mass);
            Lbar_pT_bin_vector.push_back(pT_bin_corr);
            Lbar_eta_bin_vector.push_back(eta_bin);

            piBar_vector.push_back(pi_Lorenz_vector);
            piBar_tag_vector.push_back(pi_InEventID);

            pBar_vector.push_back(p_Lorenz_vector);
            pBar_tag_vector.push_back(p_InEventID);

          }
        }
      }

      //LS charge combinations
      if(charge != 0) //same event as in previous iteration
      {
        //cuts
        if( cuts(year, L_y, pi_hasTOFinfo, p_dca, pi_dca, L_DCAdaughters, L_decayL, cos(L_theta)) && pT_bin_corr != -1 && eta_bin != -1 &&
            ( cut_type != 1 || cuts_topo_sys_err(p_dca, pi_dca, L_DCAdaughters, L_decayL, cos(L_theta) ) ) &&
            ( cut_type != 2 || cuts_pT_sys_err(p_pt, pi_pt) ) )
        //if( cuts(year, L_y, p_hasTOFinfo) && pT_bin_corr != -1 && eta_bin != -1)
        {
          if( p_ch == 1 )
          {
            //QA histograms
            L_Minv_LS[pT_bin]->Fill(L_mass);

            L_pT_vs_pi_pT_LS->Fill(L_Lorentz_vector.Pt(), pi_Lorenz_vector.Pt());

            //fill vectors
            L_vector_background.push_back(L_Lorentz_vector);
            L_pT_bin_vector_background.push_back(pT_bin_corr);
            L_eta_bin_vector_background.push_back(eta_bin);

            pi_vector_background.push_back(pi_Lorenz_vector);
            pi_tag_vector_background.push_back(pi_InEventID);

            p_vector_background.push_back(p_Lorenz_vector);
            p_tag_vector_background.push_back(p_InEventID);

          }
          else if( p_ch == -1 )
          {
            //QA histograms
            Lbar_Minv_LS[pT_bin]->Fill(L_mass);

            Lbar_pT_vs_pi_pT_LS->Fill(L_Lorentz_vector.Pt(), pi_Lorenz_vector.Pt());

            //fill vectors
            Lbar_vector_background.push_back(L_Lorentz_vector);
            Lbar_pT_bin_vector_background.push_back(pT_bin_corr);
            Lbar_eta_bin_vector_background.push_back(eta_bin);

            piBar_vector_background.push_back(pi_Lorenz_vector);
            piBar_tag_vector_background.push_back(pi_InEventID);

            pBar_vector_background.push_back(p_Lorenz_vector);
            pBar_tag_vector_background.push_back(p_InEventID);

          }
        }
      }

    }
    else if(eventId != eventID_last)//new event
    {
      nEvents++;

      eventID_last = eventId; //store new event ID

      //reset nFill flags
      nFillLLbar_sig = 0;
      nFillLL_sig = 0;
      nFillLbarLbar_sig = 0;

      // US-US L pairs - signal + background
      //at least one L-Lbar pair in event
      if(L_vector.size() > 0 && Lbar_vector.size() > 0)
      {
        //cout<<eventId<<endl;
        n_LLbar_per_event_hist->Fill(L_vector.size()*Lbar_vector.size());

        for(unsigned int iLambda = 0; iLambda < L_vector.size(); iLambda++)
        {
          for(unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector.size(); iLambdaBar++)
          {
            //check that L and Lbar don't share decay daughters
            //here possible only with mis-PID
            if(pi_tag_vector.at(iLambda) == pBar_tag_vector.at(iLambdaBar)) continue;
            if(p_tag_vector.at(iLambda) == piBar_tag_vector.at(iLambdaBar)) continue;

            //if(pi_vector.at(iLambda).Eta() == pBar_vector.at(iLambdaBar).Eta()) continue;
            //if(p_vector.at(iLambda).Eta() == piBar_vector.at(iLambdaBar).Eta()) continue;

            //if(fabs(pi_vector.at(iLambda).Eta() - pBar_vector.at(iLambdaBar).Eta()) < 1e-5) continue;
            //if(fabs(p_vector.at(iLambda).Eta() - piBar_vector.at(iLambdaBar).Eta()) < 1e-5) continue;

            //fill Minv histogram

            //if(nFillLLbar_sig > 1) continue;

            L0_inv_mass_vs_L0bar_inv_mass_all_US->Fill(L_vector.at(iLambda).M(), Lbar_vector.at(iLambdaBar).M());

            L0_inv_mass_vs_L0bar_inv_mass_US[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(L_vector.at(iLambda).M(), Lbar_vector.at(iLambdaBar).M());
            L0_inv_mass_vs_L0bar_inv_mass_US_zoom[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(L_vector.at(iLambda).M(), Lbar_vector.at(iLambdaBar).M());

            //------------------------------------------------------

            L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_Lbar[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(L_vector.at(iLambda).M(), Lbar_vector_pi_mass.at(iLambdaBar).M());
            L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_L[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(L_vector_pi_mass.at(iLambda).M(), Lbar_vector.at(iLambdaBar).M());
            L0_inv_mass_vs_L0bar_inv_mass_US_pi_mass_both[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(L_vector_pi_mass.at(iLambda).M(), Lbar_vector_pi_mass.at(iLambdaBar).M());

            //---------------------------------



            //fill new vectors here with Minv1 and Minv2, and with cos(theta*)
            //analze stored pairs later, just within selected Minv window
            //store also pT bin of both L in the pair (maybe eta)

            //save info for cos(theta*) to be analyzed later
            double L_Lbar_pairThetaStar = LpairThetaStar(L_vector.at(iLambda), p_vector.at(iLambda), Lbar_vector.at(iLambdaBar), pBar_vector.at(iLambdaBar));

            L_Lbar_L_mom.push_back(L_vector.at(iLambda));
            L_Lbar_Lbar_mom.push_back(Lbar_vector.at(iLambdaBar));

            L_Lbar_p_mom.push_back(p_vector.at(iLambda).Vect());
            L_Lbar_pi_mom.push_back(pi_vector.at(iLambda).Vect());

            L_Lbar_pBar_mom.push_back(pBar_vector.at(iLambdaBar).Vect());
            L_Lbar_piBar_mom.push_back(piBar_vector.at(iLambdaBar).Vect());

            L_Lbar_cos_theta.push_back(TMath::Cos(L_Lbar_pairThetaStar));

            L_Lbar_Minv_L.push_back(L_vector.at(iLambda).M());
            L_Lbar_Minv_Lbar.push_back(Lbar_vector.at(iLambdaBar).M());

            L_Lbar_pT_bin_L.push_back(L_pT_bin_vector.at(iLambda));
            L_Lbar_pT_bin_Lbar.push_back(Lbar_pT_bin_vector.at(iLambdaBar));

            L_Lbar_eta_bin_L.push_back(L_eta_bin_vector.at(iLambda));
            L_Lbar_eta_bin_Lbar.push_back(Lbar_eta_bin_vector.at(iLambdaBar));

            //nFillLLbar_sig++;//for testing

          }
        }
      }

      //at least one L0-L0 pair in event
      if(L_vector.size() > 1)
      {
        n_L_per_event_hist->Fill(L_vector.size());

        for(unsigned int iLambda1 = 0; iLambda1 < L_vector.size(); iLambda1++)
        {
          for(unsigned int iLambda2 = iLambda1+1; iLambda2 < L_vector.size(); iLambda2++)
          {
            //check that no daughters are shared by the two L in the pair
            if(pi_tag_vector.at(iLambda1) == pi_tag_vector.at(iLambda2)) continue;
            if(p_tag_vector.at(iLambda1) == p_tag_vector.at(iLambda2)) continue;

            L0_inv_mass_vs_L0_inv_mass_all_US->Fill(L_vector.at(iLambda1).M(), L_vector.at(iLambda2).M());

            L0_inv_mass_vs_L0_inv_mass_US[L_pT_bin_vector.at(iLambda1)][L_pT_bin_vector.at(iLambda2)]->Fill(L_vector.at(iLambda1).M(), L_vector.at(iLambda2).M());
            L0_inv_mass_vs_L0_inv_mass_US_zoom[L_pT_bin_vector.at(iLambda1)][L_pT_bin_vector.at(iLambda2)]->Fill(L_vector.at(iLambda1).M(), L_vector.at(iLambda2).M());


            L0_inv_mass_vs_L0_inv_mass_US_pi_mass_one[L_pT_bin_vector.at(iLambda1)][L_pT_bin_vector.at(iLambda2)]->Fill(L_vector.at(iLambda1).M(), L_vector_pi_mass.at(iLambda2).M());
            L0_inv_mass_vs_L0_inv_mass_US_pi_mass_both[L_pT_bin_vector.at(iLambda1)][L_pT_bin_vector.at(iLambda2)]->Fill(L_vector_pi_mass.at(iLambda1).M(), L_vector_pi_mass.at(iLambda2).M());

            //---------------------------------



            //save info for cos(theta*) to be analyzed later
            double L_L_pairThetaStar = LpairThetaStar(L_vector.at(iLambda1), p_vector.at(iLambda1), L_vector.at(iLambda2), p_vector.at(iLambda2));

            L_L_L1_mom.push_back(L_vector.at(iLambda1));
            L_L_L2_mom.push_back(L_vector.at(iLambda2));

            L_L_p1_mom.push_back(p_vector.at(iLambda1).Vect());
            L_L_pi1_mom.push_back(pi_vector.at(iLambda1).Vect());

            L_L_p2_mom.push_back(p_vector.at(iLambda2).Vect());
            L_L_pi2_mom.push_back(pi_vector.at(iLambda2).Vect());

            L_L_cos_theta.push_back(TMath::Cos(L_L_pairThetaStar));

            L_L_Minv_L1.push_back(L_vector.at(iLambda1).M());
            L_L_Minv_L2.push_back(L_vector.at(iLambda2).M());

            L_L_pT_bin_L1.push_back(L_pT_bin_vector.at(iLambda1));
            L_L_pT_bin_L2.push_back(L_pT_bin_vector.at(iLambda2));

            L_L_eta_bin_L1.push_back(L_eta_bin_vector.at(iLambda1));
            L_L_eta_bin_L2.push_back(L_eta_bin_vector.at(iLambda2));

            //nFillLL_sig++;//for testing

          }
        }
      }

      //at least one L0bar-L0bar pair in event
      if(Lbar_vector.size() > 1)
      {
        n_Lbar_per_event_hist->Fill(Lbar_vector.size());

        //for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < Lbar_vector.size(); iLambdaBar1++)
        for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < Lbar_vector.size(); iLambdaBar1++)
        {
          //for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < Lbar_vector.size(); iLambdaBar2++)
          for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < Lbar_vector.size(); iLambdaBar2++)
          {
            //check that no daughters are shared by the two L in the pair
            if(piBar_tag_vector.at(iLambdaBar1) == piBar_tag_vector.at(iLambdaBar2)) continue;
            if(pBar_tag_vector.at(iLambdaBar1) == pBar_tag_vector.at(iLambdaBar2)) continue;

            //if(piBar_vector.at(iLambdaBar1).Eta() == piBar_vector.at(iLambdaBar2).Eta()) continue;
            //if(pBar_vector.at(iLambdaBar1).Eta() == pBar_vector.at(iLambdaBar2).Eta()) continue;

            //if(fabs(piBar_vector.at(iLambdaBar1).Eta() - piBar_vector.at(iLambdaBar2).Eta()) < 1e-5) continue;
            //if(fabs(pBar_vector.at(iLambdaBar1).Eta() - pBar_vector.at(iLambdaBar2).Eta()) < 1e-5) continue;

            //L0bar_inv_mass_vs_L0bar_inv_mass_US[Lbar_pT_bin_vector.at(iLambdaBar1)][Lbar_pT_bin_vector.at(iLambdaBar2)]->Fill(Lbar_vector.at(iLambdaBar1).M(), Lbar_vector.at(iLambdaBar2).M());

            //if(nFillLbarLbar_sig > 1) continue;

            L0bar_inv_mass_vs_L0bar_inv_mass_all_US->Fill(Lbar_vector.at(iLambdaBar1).M(), Lbar_vector.at(iLambdaBar2).M());

            L0bar_inv_mass_vs_L0bar_inv_mass_US[Lbar_pT_bin_vector.at(iLambdaBar1)][Lbar_pT_bin_vector.at(iLambdaBar2)]->Fill(Lbar_vector.at(iLambdaBar1).M(), Lbar_vector.at(iLambdaBar2).M());
            L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[Lbar_pT_bin_vector.at(iLambdaBar1)][Lbar_pT_bin_vector.at(iLambdaBar2)]->Fill(Lbar_vector.at(iLambdaBar1).M(), Lbar_vector.at(iLambdaBar2).M());


            L0bar_inv_mass_vs_L0bar_inv_mass_US_pi_mass_one[Lbar_pT_bin_vector.at(iLambdaBar1)][Lbar_pT_bin_vector.at(iLambdaBar2)]->Fill(Lbar_vector.at(iLambdaBar1).M(), Lbar_vector_pi_mass.at(iLambdaBar2).M());
            L0bar_inv_mass_vs_L0bar_inv_mass_US_pi_mass_both[Lbar_pT_bin_vector.at(iLambdaBar1)][Lbar_pT_bin_vector.at(iLambdaBar2)]->Fill(Lbar_vector_pi_mass.at(iLambdaBar1).M(), Lbar_vector_pi_mass.at(iLambdaBar2).M());

            //---------------------------------


            //save info for cos(theta*) to be analyzed later
            double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_vector.at(iLambdaBar1), pBar_vector.at(iLambdaBar1), Lbar_vector.at(iLambdaBar2), pBar_vector.at(iLambdaBar2));

            Lbar_Lbar_Lbar1_mom.push_back(Lbar_vector.at(iLambdaBar1));
            Lbar_Lbar_Lbar2_mom.push_back(Lbar_vector.at(iLambdaBar2));

            Lbar_Lbar_pBar1_mom.push_back(pBar_vector.at(iLambdaBar1).Vect());
            Lbar_Lbar_piBar1_mom.push_back(piBar_vector.at(iLambdaBar1).Vect());

            Lbar_Lbar_pBar2_mom.push_back(pBar_vector.at(iLambdaBar2).Vect());
            Lbar_Lbar_piBar2_mom.push_back(piBar_vector.at(iLambdaBar2).Vect());

            Lbar_Lbar_cos_theta.push_back(TMath::Cos(Lbar_Lbar_pairThetaStar));

            Lbar_Lbar_Minv_Lbar1.push_back(Lbar_vector.at(iLambdaBar1).M());
            Lbar_Lbar_Minv_Lbar2.push_back(Lbar_vector.at(iLambdaBar2).M());

            Lbar_Lbar_pT_bin_Lbar1.push_back(Lbar_pT_bin_vector.at(iLambdaBar1));
            Lbar_Lbar_pT_bin_Lbar2.push_back(Lbar_pT_bin_vector.at(iLambdaBar2));

            Lbar_Lbar_eta_bin_Lbar1.push_back(Lbar_eta_bin_vector.at(iLambdaBar1));
            Lbar_Lbar_eta_bin_Lbar2.push_back(Lbar_eta_bin_vector.at(iLambdaBar2));

            //nFillLbarLbar_sig++; //for testing
          }
        }
      }
      //_________________________________________________________________________________________________________




      // LS-LS L pairs - background - just continuum

      if(L_vector_background.size() > 0 && Lbar_vector_background.size() > 0)
      {
        for(unsigned int iLambda = 0; iLambda < L_vector_background.size(); iLambda++)
        {
          for(unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector_background.size(); iLambdaBar++)
          {
            //if( fabs(p_vector_background.at(iLambda).Rapidity()) < 1e-5  ) continue;

            L0_inv_mass_vs_L0bar_inv_mass_all_LS->Fill(L_vector_background.at(iLambda).M(), Lbar_vector_background.at(iLambdaBar).M());

            L0_inv_mass_vs_L0bar_inv_mass_LS[L_pT_bin_vector_background.at(iLambda)][Lbar_pT_bin_vector_background.at(iLambdaBar)]->Fill(L_vector_background.at(iLambda).M(), Lbar_vector_background.at(iLambdaBar).M());

          }
        }
      }


      //at least one L0-L0 pair in event
      if(L_vector_background.size() > 1)
      {
        for(unsigned int iLambda1 = 0; iLambda1 < L_vector_background.size(); iLambda1++)
        {
          for(unsigned int iLambda2 = iLambda1+1; iLambda2 < L_vector_background.size(); iLambda2++)
          {
            //to prevent auto-correlation
            if(pi_vector_background.at(iLambda1).Rapidity() == pi_vector_background.at(iLambda2).Rapidity()) continue;
            if(p_vector_background.at(iLambda1).Rapidity() == p_vector_background.at(iLambda2).Rapidity()) continue;

            L0_inv_mass_vs_L0_inv_mass_all_LS->Fill(L_vector_background.at(iLambda1).M(), L_vector_background.at(iLambda2).M());

            L0_inv_mass_vs_L0_inv_mass_LS[L_pT_bin_vector_background.at(iLambda1)][L_pT_bin_vector_background.at(iLambda2)]->Fill(L_vector_background.at(iLambda1).M(), L_vector_background.at(iLambda2).M());

          }
        }
      }

      //at least one L0bar-L0bar pair in event
      if(Lbar_vector_background.size() > 1)
      {
        for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < Lbar_vector_background.size(); iLambdaBar1++)
        {
          for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < Lbar_vector_background.size(); iLambdaBar2++)
          {
            if(piBar_vector_background.at(iLambdaBar1).Rapidity() == piBar_vector_background.at(iLambdaBar2).Rapidity()) continue;
            if(pBar_vector_background.at(iLambdaBar1).Rapidity() == pBar_vector_background.at(iLambdaBar2).Rapidity()) continue;

            L0bar_inv_mass_vs_L0bar_inv_mass_all_LS->Fill(Lbar_vector_background.at(iLambdaBar1).M(), Lbar_vector_background.at(iLambdaBar2).M());

            L0bar_inv_mass_vs_L0bar_inv_mass_LS[Lbar_pT_bin_vector_background.at(iLambdaBar1)][Lbar_pT_bin_vector_background.at(iLambdaBar2)]->Fill(Lbar_vector_background.at(iLambdaBar1).M(), Lbar_vector_background.at(iLambdaBar2).M());

          }
        }
      }

      //___________________________________________________________________________________________________


      // US-LS and LS-US L pairs - background - case when US L (can be signal or background) is paired with background L

      //Background for L-Lbar
      // L - US, Lbar - LS
      if(L_vector.size() > 0 && Lbar_vector_background.size() > 0)
      {
        //cout<<eventId<<endl;

        for(unsigned int iLambda = 0; iLambda < L_vector.size(); iLambda++)
        {
          for(unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector_background.size(); iLambdaBar++)
          {
            //check that no daughters are shared, including mis-PID particles
            //L(p pi-)-Lbar_LS(pbar pi-)
            if(pi_tag_vector.at(iLambda) == piBar_tag_vector_background.at(iLambdaBar)) continue;
            if(pi_tag_vector.at(iLambda) == pBar_tag_vector_background.at(iLambdaBar)) continue;

            //if(pi_vector.at(iLambda).Eta() == piBar_vector_background.at(iLambdaBar).Eta()) continue;
            //if(p_vector.at(iLambda).Eta() == pBar_vector_background.at(iLambdaBar).Eta()) continue;

            //L(p pi-) Lbar_LS(pBar pi-)
            //if(fabs(pi_vector.at(iLambda).Eta() - pBar_vector_background.at(iLambdaBar).Eta()) < 1e-5) continue;
            //if(fabs(pi_vector.at(iLambda).Eta() - piBar_vector_background.at(iLambdaBar).Eta()) < 1e-5) continue;

            //fill Minv histogram
            L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->Fill(L_vector.at(iLambda).M(), Lbar_vector_background.at(iLambdaBar).M());

            L0_inv_mass_vs_L0bar_inv_mass_US_LS[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector_background.at(iLambdaBar)]->Fill(L_vector.at(iLambda).M(), Lbar_vector_background.at(iLambdaBar).M());
            L0_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector_background.at(iLambdaBar)]->Fill(L_vector.at(iLambda).M(), Lbar_vector_background.at(iLambdaBar).M());


            //save info for cos(theta*) to be analyzed later
            double L_Lbar_pairThetaStar = LpairThetaStar(L_vector.at(iLambda), p_vector.at(iLambda), Lbar_vector_background.at(iLambdaBar), pBar_vector_background.at(iLambdaBar));

            //all pairs
            L_Lbar_L_mom_back.push_back(L_vector.at(iLambda));
            L_Lbar_Lbar_mom_back.push_back(Lbar_vector_background.at(iLambdaBar));

            L_Lbar_p_mom_back.push_back(p_vector.at(iLambda).Vect());
            L_Lbar_pi_mom_back.push_back(pi_vector.at(iLambda).Vect());

            L_Lbar_pBar_mom_back.push_back(pBar_vector_background.at(iLambdaBar).Vect());
            L_Lbar_piBar_mom_back.push_back(piBar_vector_background.at(iLambdaBar).Vect());

            L_Lbar_cos_theta_back.push_back(TMath::Cos(L_Lbar_pairThetaStar));

            L_Lbar_Minv_L_back.push_back(L_vector.at(iLambda).M());
            L_Lbar_Minv_Lbar_back.push_back(Lbar_vector_background.at(iLambdaBar).M());

            L_Lbar_pT_bin_L_back.push_back(L_pT_bin_vector.at(iLambda));
            L_Lbar_pT_bin_Lbar_back.push_back(Lbar_pT_bin_vector_background.at(iLambdaBar));

            L_Lbar_eta_bin_L_back.push_back(L_eta_bin_vector.at(iLambda));
            L_Lbar_eta_bin_Lbar_back.push_back(Lbar_eta_bin_vector_background.at(iLambdaBar));

            //L(US)-Lbar(LS) for ME reweight
            L_Lbar_US_LS_flag.push_back(1);

          }
        }
      }

      // L - LS, Lbar - US
      if(L_vector_background.size() > 0 && Lbar_vector.size() > 0)
      {
        for(unsigned int iLambda = 0; iLambda < L_vector_background.size(); iLambda++)
        {
          for(unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector.size(); iLambdaBar++)
          {
            //check that no daughters are shared, including mis-PID particles
            //L_LS(p pi+)-Lbar_LS(pbar pi+)
            if(pi_tag_vector_background.at(iLambda) == piBar_tag_vector.at(iLambdaBar)) continue;
            if(p_tag_vector_background.at(iLambda) == piBar_tag_vector.at(iLambdaBar)) continue;

            //if(pi_vector_background.at(iLambda).Eta() == piBar_vector.at(iLambdaBar).Eta()) continue;
            //if(p_vector_background.at(iLambda).Eta() == pBar_vector.at(iLambdaBar).Eta()) continue;

            //L_LS(p pi+) Lbar(pBar pi+)
            //if(fabs(pi_vector_background.at(iLambda).Eta() - piBar_vector.at(iLambdaBar).Eta()) < 1e-5) continue;
            //if(fabs(p_vector_background.at(iLambda).Eta() - piBar_vector.at(iLambdaBar).Eta()) < 1e-5) continue;

            //fill Minv histogram
            L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->Fill(L_vector_background.at(iLambda).M(), Lbar_vector.at(iLambdaBar).M());

            L0_inv_mass_vs_L0bar_inv_mass_US_LS[L_pT_bin_vector_background.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(L_vector_background.at(iLambda).M(), Lbar_vector.at(iLambdaBar).M());
            L0_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[L_pT_bin_vector_background.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(L_vector_background.at(iLambda).M(), Lbar_vector.at(iLambdaBar).M());

            //save info for cos(theta*) to be analyzed later
            double L_Lbar_pairThetaStar = LpairThetaStar(L_vector_background.at(iLambda), p_vector_background.at(iLambda), Lbar_vector.at(iLambdaBar), pBar_vector.at(iLambdaBar));

            //all pairs
            L_Lbar_L_mom_back.push_back(L_vector_background.at(iLambda));
            L_Lbar_Lbar_mom_back.push_back(Lbar_vector.at(iLambdaBar));

            L_Lbar_p_mom_back.push_back(p_vector_background.at(iLambda).Vect());
            L_Lbar_pi_mom_back.push_back(pi_vector_background.at(iLambda).Vect());

            L_Lbar_pBar_mom_back.push_back(pBar_vector.at(iLambdaBar).Vect());
            L_Lbar_piBar_mom_back.push_back(piBar_vector.at(iLambdaBar).Vect());

            L_Lbar_cos_theta_back.push_back(TMath::Cos(L_Lbar_pairThetaStar));

            L_Lbar_Minv_L_back.push_back(L_vector_background.at(iLambda).M());
            L_Lbar_Minv_Lbar_back.push_back(Lbar_vector.at(iLambdaBar).M());

            L_Lbar_pT_bin_L_back.push_back(L_pT_bin_vector_background.at(iLambda));
            L_Lbar_pT_bin_Lbar_back.push_back(Lbar_pT_bin_vector.at(iLambdaBar));

            L_Lbar_eta_bin_L_back.push_back(L_eta_bin_vector_background.at(iLambda));
            L_Lbar_eta_bin_Lbar_back.push_back(Lbar_eta_bin_vector.at(iLambdaBar));

            //L(LS)-Lbar(US) for ME reweight
            L_Lbar_US_LS_flag.push_back(2);

          }
        }
      }

      //background for L-L
      // L - US, L - LS
      if(L_vector.size() > 0 && L_vector_background.size() > 0)
      {
        for(unsigned int iLambda = 0; iLambda < L_vector.size(); iLambda++)
        {
          for(unsigned int iLambdaBckg = 0; iLambdaBckg < L_vector_background.size(); iLambdaBckg++)
          {
            //check that no daughters are shared, including mis-PID particles
            //L(p pi-)-L_LS(p pi+)
            if(p_tag_vector.at(iLambda) == p_tag_vector_background.at(iLambdaBckg) ) continue;
            if(p_tag_vector.at(iLambda) == pi_tag_vector_background.at(iLambdaBckg) ) continue;

            //if(p_vector.at(iLambda).Eta() == p_vector_background.at(iLambdaBckg).Eta() ) continue;
            //if(pi_vector.at(iLambda).Eta() == pi_vector_background.at(iLambdaBckg).Eta() ) continue;

            //L(p pi-) L_LS(p pi+)
            //if(fabs(p_vector.at(iLambda).Eta() - p_vector_background.at(iLambdaBckg).Eta() ) < 1e-5) continue;
            //if(fabs(p_vector.at(iLambda).Eta() - pi_vector_background.at(iLambdaBckg).Eta() ) < 1e-5) continue;

            if( nFillLL_bckg % 2 == 0 ) //fill LL pair wit US at x axis and LS at y axis
            {
              L0_inv_mass_vs_L0_inv_mass_all_US_LS->Fill(L_vector.at(iLambda).M(), L_vector_background.at(iLambdaBckg).M());

              L0_inv_mass_vs_L0_inv_mass_US_LS[L_pT_bin_vector.at(iLambda)][L_pT_bin_vector_background.at(iLambdaBckg)]->Fill(L_vector.at(iLambda).M(), L_vector_background.at(iLambdaBckg).M());
              L0_inv_mass_vs_L0_inv_mass_US_LS_zoom[L_pT_bin_vector.at(iLambda)][L_pT_bin_vector_background.at(iLambdaBckg)]->Fill(L_vector.at(iLambda).M(), L_vector_background.at(iLambdaBckg).M());

              //save info for cos(theta*) to be analyzed later
              double L_L_pairThetaStar = LpairThetaStar(L_vector.at(iLambda), p_vector.at(iLambda), L_vector_background.at(iLambdaBckg), p_vector_background.at(iLambdaBckg));

              L_L_L1_mom_back.push_back(L_vector.at(iLambda));
              L_L_L2_mom_back.push_back(L_vector_background.at(iLambdaBckg));

              L_L_p1_mom_back.push_back(p_vector.at(iLambda).Vect());
              L_L_pi1_mom_back.push_back(pi_vector.at(iLambda).Vect());

              L_L_p2_mom_back.push_back(p_vector_background.at(iLambdaBckg).Vect());
              L_L_pi2_mom_back.push_back(pi_vector_background.at(iLambdaBckg).Vect());

              L_L_cos_theta_back.push_back(TMath::Cos(L_L_pairThetaStar));

              L_L_Minv_L1_back.push_back(L_vector.at(iLambda).M());
              L_L_Minv_L2_back.push_back(L_vector_background.at(iLambdaBckg).M());

              L_L_pT_bin_L1_back.push_back(L_pT_bin_vector.at(iLambda));
              L_L_pT_bin_L2_back.push_back(L_pT_bin_vector_background.at(iLambdaBckg));

              L_L_eta_bin_L1_back.push_back(L_eta_bin_vector.at(iLambda));
              L_L_eta_bin_L2_back.push_back(L_eta_bin_vector_background.at(iLambdaBckg));

              //L(US)-L(LS) for ME reweight
              L_L_US_LS_flag.push_back(1);

            }
            else //fill LL pair wit US at y axis and LS at x axis, i.e. opposite order than above
            {
              L0_inv_mass_vs_L0_inv_mass_all_US_LS->Fill(L_vector_background.at(iLambdaBckg).M(), L_vector.at(iLambda).M());

              L0_inv_mass_vs_L0_inv_mass_US_LS[L_pT_bin_vector_background.at(iLambdaBckg)][L_pT_bin_vector.at(iLambda)]->Fill(L_vector_background.at(iLambdaBckg).M(), L_vector.at(iLambda).M());
              L0_inv_mass_vs_L0_inv_mass_US_LS_zoom[L_pT_bin_vector_background.at(iLambdaBckg)][L_pT_bin_vector.at(iLambda)]->Fill(L_vector_background.at(iLambdaBckg).M(), L_vector.at(iLambda).M());

              //save info for cos(theta*) to be analyzed later
              double L_L_pairThetaStar = LpairThetaStar(L_vector_background.at(iLambdaBckg), p_vector_background.at(iLambdaBckg), L_vector.at(iLambda), p_vector.at(iLambda));

              L_L_L1_mom_back.push_back(L_vector_background.at(iLambdaBckg));
              L_L_L2_mom_back.push_back(L_vector.at(iLambda));

              L_L_p1_mom_back.push_back(p_vector_background.at(iLambdaBckg).Vect());
              L_L_pi1_mom_back.push_back(pi_vector_background.at(iLambdaBckg).Vect());

              L_L_p2_mom_back.push_back(p_vector.at(iLambda).Vect());
              L_L_pi2_mom_back.push_back(pi_vector.at(iLambda).Vect());

              L_L_cos_theta_back.push_back(TMath::Cos(L_L_pairThetaStar));

              L_L_Minv_L1_back.push_back(L_vector_background.at(iLambdaBckg).M());
              L_L_Minv_L2_back.push_back(L_vector.at(iLambda).M());

              L_L_pT_bin_L1_back.push_back(L_pT_bin_vector_background.at(iLambdaBckg));
              L_L_pT_bin_L2_back.push_back(L_pT_bin_vector.at(iLambda));

              L_L_eta_bin_L1_back.push_back(L_eta_bin_vector_background.at(iLambdaBckg));
              L_L_eta_bin_L2_back.push_back(L_eta_bin_vector.at(iLambda));

              //L(LS)-L(US) for ME reweight
              L_L_US_LS_flag.push_back(2);

            }

            nFillLL_bckg++;


          }
        }

      }


      //background for Lbar-Lbar
      // Lbar - US, Lbar - LS
      if(Lbar_vector.size() > 0 && Lbar_vector_background.size() > 0)
      {
        for(unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector.size(); iLambdaBar++)
        {
          for(unsigned int iLambdaBarBckg = 0; iLambdaBarBckg < Lbar_vector_background.size(); iLambdaBarBckg++)
          {
            //check that no daughters are shared, including mis-PID particles
            //Lbar(pbar pi+)-Lbar_LS(pbar pi-)
            if(pBar_tag_vector.at(iLambdaBar) == pBar_tag_vector_background.at(iLambdaBarBckg) ) continue;
            if(pBar_tag_vector.at(iLambdaBar) == piBar_tag_vector_background.at(iLambdaBarBckg) ) continue;

            //if(pBar_vector.at(iLambdaBar).Eta() == pBar_vector_background.at(iLambdaBarBckg).Eta() ) continue;
            //if(piBar_vector.at(iLambdaBar).Eta() == piBar_vector_background.at(iLambdaBarBckg).Eta() ) continue;

            //Lbar(pBar pi+) Lbar_LS(pBar pi-)
            //if(fabs(pBar_vector.at(iLambdaBar).Eta() - pBar_vector_background.at(iLambdaBarBckg).Eta() ) < 1e-5) continue;
            //if(fabs(pBar_vector.at(iLambdaBar).Eta() - piBar_vector_background.at(iLambdaBarBckg).Eta() ) < 1e-5) continue;

            if( nFillLbarLbar_backg % 2 == 0 ) //fill LL pair wit US at x axis and LS at y axis
            {
              L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->Fill(Lbar_vector.at(iLambdaBar).M(), Lbar_vector_background.at(iLambdaBarBckg).M());

              L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[Lbar_pT_bin_vector.at(iLambdaBar)][Lbar_pT_bin_vector_background.at(iLambdaBarBckg)]->Fill(Lbar_vector.at(iLambdaBar).M(), Lbar_vector_background.at(iLambdaBarBckg).M());
              L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[Lbar_pT_bin_vector.at(iLambdaBar)][Lbar_pT_bin_vector_background.at(iLambdaBarBckg)]->Fill(Lbar_vector.at(iLambdaBar).M(), Lbar_vector_background.at(iLambdaBarBckg).M());

              //save info for cos(theta*) to be analyzed later
              double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_vector.at(iLambdaBar), pBar_vector.at(iLambdaBar), Lbar_vector_background.at(iLambdaBarBckg), pBar_vector_background.at(iLambdaBarBckg));

              Lbar_Lbar_Lbar1_mom_back.push_back(Lbar_vector.at(iLambdaBar));
              Lbar_Lbar_Lbar2_mom_back.push_back(Lbar_vector_background.at(iLambdaBarBckg));

              Lbar_Lbar_pBar1_mom_back.push_back(pBar_vector.at(iLambdaBar).Vect());
              Lbar_Lbar_piBar1_mom_back.push_back(piBar_vector.at(iLambdaBar).Vect());

              Lbar_Lbar_pBar2_mom_back.push_back(pBar_vector_background.at(iLambdaBarBckg).Vect());
              Lbar_Lbar_piBar2_mom_back.push_back(piBar_vector_background.at(iLambdaBarBckg).Vect());

              Lbar_Lbar_cos_theta_back.push_back(TMath::Cos(Lbar_Lbar_pairThetaStar));

              Lbar_Lbar_Minv_Lbar1_back.push_back(Lbar_vector.at(iLambdaBar).M());
              Lbar_Lbar_Minv_Lbar2_back.push_back(Lbar_vector_background.at(iLambdaBarBckg).M());

              Lbar_Lbar_pT_bin_Lbar1_back.push_back(Lbar_pT_bin_vector.at(iLambdaBar));
              Lbar_Lbar_pT_bin_Lbar2_back.push_back(Lbar_pT_bin_vector_background.at(iLambdaBarBckg));

              Lbar_Lbar_eta_bin_Lbar1_back.push_back(Lbar_eta_bin_vector.at(iLambdaBar));
              Lbar_Lbar_eta_bin_Lbar2_back.push_back(Lbar_eta_bin_vector_background.at(iLambdaBarBckg));

              //Lbar(US)-Lbar(LS) for ME reweight
              Lbar_Lbar_US_LS_flag.push_back(1);

            }
            else //fill LL pair wit US at y axis and LS at x axis, i.e. opposite order than above
            {
              L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->Fill(Lbar_vector_background.at(iLambdaBarBckg).M(), Lbar_vector.at(iLambdaBar).M());

              L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[Lbar_pT_bin_vector_background.at(iLambdaBarBckg)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(Lbar_vector_background.at(iLambdaBarBckg).M(), Lbar_vector.at(iLambdaBar).M());
              L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[Lbar_pT_bin_vector_background.at(iLambdaBarBckg)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(Lbar_vector_background.at(iLambdaBarBckg).M(), Lbar_vector.at(iLambdaBar).M());

              //save info for cos(theta*) to be analyzed later
              double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_vector_background.at(iLambdaBarBckg), pBar_vector_background.at(iLambdaBarBckg), Lbar_vector.at(iLambdaBar), pBar_vector.at(iLambdaBar));

              Lbar_Lbar_Lbar1_mom_back.push_back(Lbar_vector_background.at(iLambdaBarBckg));
              Lbar_Lbar_Lbar2_mom_back.push_back(Lbar_vector.at(iLambdaBar));

              Lbar_Lbar_pBar1_mom_back.push_back(pBar_vector_background.at(iLambdaBarBckg).Vect());
              Lbar_Lbar_piBar1_mom_back.push_back(piBar_vector_background.at(iLambdaBarBckg).Vect());

              Lbar_Lbar_pBar2_mom_back.push_back(pBar_vector.at(iLambdaBar).Vect());
              Lbar_Lbar_piBar2_mom_back.push_back(piBar_vector.at(iLambdaBar).Vect());

              Lbar_Lbar_cos_theta_back.push_back(TMath::Cos(Lbar_Lbar_pairThetaStar));

              Lbar_Lbar_Minv_Lbar1_back.push_back(Lbar_vector_background.at(iLambdaBarBckg).M());
              Lbar_Lbar_Minv_Lbar2_back.push_back(Lbar_vector.at(iLambdaBar).M());

              Lbar_Lbar_pT_bin_Lbar1_back.push_back(Lbar_pT_bin_vector_background.at(iLambdaBarBckg));
              Lbar_Lbar_pT_bin_Lbar2_back.push_back(Lbar_pT_bin_vector.at(iLambdaBar));

              Lbar_Lbar_eta_bin_Lbar1_back.push_back(Lbar_eta_bin_vector_background.at(iLambdaBarBckg));
              Lbar_Lbar_eta_bin_Lbar2_back.push_back(Lbar_eta_bin_vector.at(iLambdaBar));

              //Lbar(LS)-Lbar(US) for ME reweight
              Lbar_Lbar_US_LS_flag.push_back(2);

            }

            nFillLbarLbar_backg++;

          }
        }

      }


      //___________________________________________________________________________________________________

      //fill vectors for mixed event
      //need to select more - Minv cut will be applied later

      //signal+background
      //same event L-Lbar for ME
      if( L_vector.size() > 0 && Lbar_vector.size() > 0 && L_Lbar_L_vector_ME_SE.size() < 3e4 )
      {
        //if( fabs(p_vector.at(0).Eta() - piBar_vector.at(0).Eta()) > 1e-5 && fabs(pBar_vector.at(0).Eta() - pi_vector.at(0).Eta()) > 1e-5 )
        if( p_tag_vector.at(0) != piBar_tag_vector.at(0) && pBar_tag_vector.at(0) != pi_tag_vector.at(0) )
        {
          //Minv pre-selection to boost usable ME statistics
          if( L_vector.at(0).M() > 1.105 && L_vector.at(0).M() < 1.13 && Lbar_vector.at(0).M() > 1.105 && Lbar_vector.at(0).M() < 1.13)
          {
            L_Lbar_L_vector_ME_SE.push_back(L_vector.at(0));
            L_Lbar_p_vector_ME_SE.push_back(p_vector.at(0));
            L_Lbar_pi_vector_ME_SE.push_back(pi_vector.at(0));

            L_Lbar_L_pT_bin_vector_ME_SE.push_back(L_pT_bin_vector.at(0));
            L_Lbar_L_eta_bin_vector_ME_SE.push_back(L_eta_bin_vector.at(0));


            L_Lbar_Lbar_vector_ME_SE.push_back(Lbar_vector.at(0));
            L_Lbar_pBar_vector_ME_SE.push_back(pBar_vector.at(0));
            L_Lbar_piBar_vector_ME_SE.push_back(piBar_vector.at(0));

            L_Lbar_Lbar_pT_bin_vector_ME_SE.push_back(Lbar_pT_bin_vector.at(0));
            L_Lbar_Lbar_eta_bin_vector_ME_SE.push_back(Lbar_eta_bin_vector.at(0));
          }
        }
      }

      //same event L-L
      if( L_vector.size() > 1 && L_L_L1_vector_ME_SE.size() < 3e4)
      {
        //check auto-correlation for SE LL pair
        //no daughter can be shared
        //if(fabs(pi_vector.at(0).Eta() - pi_vector.at(1).Eta()) > 1e-5 && fabs(p_vector.at(0).Eta() - p_vector.at(1).Eta()) > 1e-5 )
        if( pi_tag_vector.at(0) != pi_tag_vector.at(1) && p_tag_vector.at(0) != p_tag_vector.at(1) )
        {
          //Minv pre-selection to boost usable ME statistics
          if( L_vector.at(0).M() > 1.105 && L_vector.at(0).M() < 1.13 && L_vector.at(1).M() > 1.105 && L_vector.at(1).M() < 1.13)
          {
            L_L_L1_vector_ME_SE.push_back(L_vector.at(0));
            L_L_p1_vector_ME_SE.push_back(p_vector.at(0));
            L_L_pi1_vector_ME_SE.push_back(pi_vector.at(0));

            L_L_L1_pT_bin_vector_ME_SE.push_back(L_pT_bin_vector.at(0));
            L_L_L1_eta_bin_vector_ME_SE.push_back(L_eta_bin_vector.at(0));


            L_L_L2_vector_ME_SE.push_back(L_vector.at(1));
            L_L_p2_vector_ME_SE.push_back(p_vector.at(1));
            L_L_pi2_vector_ME_SE.push_back(pi_vector.at(1));

            L_L_L2_pT_bin_vector_ME_SE.push_back(L_pT_bin_vector.at(1));
            L_L_L2_eta_bin_vector_ME_SE.push_back(L_eta_bin_vector.at(1));
          }
        }


      }

      //same event Lbar-Lbar
      if( Lbar_vector.size() > 1 && Lbar_Lbar_Lbar1_vector_ME_SE.size() < 3e4)
      {
        //if(fabs(piBar_vector.at(0).Eta() - piBar_vector.at(1).Eta()) > 1e-5 && fabs(pBar_vector.at(0).Eta() - pBar_vector.at(1).Eta()) > 1e-5 )
        if( piBar_tag_vector.at(0) != piBar_tag_vector.at(1) && pBar_tag_vector.at(0) != pBar_tag_vector.at(1) )
        {
          //Minv pre-selection to boost usable ME statistics
          if( Lbar_vector.at(0).M() > 1.105 && Lbar_vector.at(0).M() < 1.13 && Lbar_vector.at(1).M() > 1.105 && Lbar_vector.at(1).M() < 1.13)
          {
            Lbar_Lbar_Lbar1_vector_ME_SE.push_back(Lbar_vector.at(0));
            Lbar_Lbar_pBar1_vector_ME_SE.push_back(pBar_vector.at(0));
            Lbar_Lbar_piBar1_vector_ME_SE.push_back(piBar_vector.at(0));

            Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.push_back(Lbar_pT_bin_vector.at(0));
            Lbar_Lbar_Lbar1_eta_bin_vector_ME_SE.push_back(Lbar_eta_bin_vector.at(0));


            Lbar_Lbar_Lbar2_vector_ME_SE.push_back(Lbar_vector.at(1));
            Lbar_Lbar_pBar2_vector_ME_SE.push_back(pBar_vector.at(1));
            Lbar_Lbar_piBar2_vector_ME_SE.push_back(piBar_vector.at(1));

            Lbar_Lbar_Lbar2_pT_bin_vector_ME_SE.push_back(Lbar_pT_bin_vector.at(1));
            Lbar_Lbar_Lbar2_eta_bin_vector_ME_SE.push_back(Lbar_eta_bin_vector.at(1));
          }
        }
      }

      //--------------------------------------------------------------------------------------------------------

      //single L and Lbar to be mixed with

      if( L_vector.size() == 1 && Lbar_vector.size() == 0 && L_vector_ME.size() < 1e5)
      {
        //Minv pre-selection to boost usable ME statistics
        if( L_vector.at(0).M() > 1.105 && L_vector.at(0).M() < 1.13 )
        {
          L_vector_ME.push_back(L_vector.at(0));
          p_vector_ME.push_back(p_vector.at(0));
          pi_vector_ME.push_back(pi_vector.at(0));

          L_pT_bin_vector_ME.push_back(L_pT_bin_vector.at(0));
          L_eta_bin_vector_ME.push_back(L_eta_bin_vector.at(0));
        }
      }

      if( L_vector.size() == 0 && Lbar_vector.size() == 1 && Lbar_vector_ME.size() < 1e5)
      {
        //Minv pre-selection to boost usable ME statistics
        if( Lbar_vector.at(0).M() > 1.105 && Lbar_vector.at(0).M() < 1.13 )
        {
          Lbar_vector_ME.push_back(Lbar_vector.at(0));
          pBar_vector_ME.push_back(pBar_vector.at(0));
          piBar_vector_ME.push_back(piBar_vector.at(0));

          Lbar_pT_bin_vector_ME.push_back(Lbar_pT_bin_vector.at(0));
          Lbar_eta_bin_vector_ME.push_back(Lbar_eta_bin_vector.at(0));
        }
      }

      //_________________________________________________________________________________________________________


      //background
      //L-Lbar
      //L(US)-Lbar(LS)
      if( L_vector.size() > 0 && Lbar_vector_background.size() > 0 && L_Lbar_L_vector_ME_SE_back_US_LS.size() < 3e4)
      {
        //L(p pi-) Lbar_LS(pBar pi-)
        //if(fabs(pi_vector.at(0).Eta() - piBar_vector_background.at(0).Eta()) > 1e-5 && fabs(pi_vector.at(0).Eta() - pBar_vector_background.at(0).Eta()) > 1e-5 )
        if( pi_tag_vector.at(0) != piBar_tag_vector_background.at(0) && pi_tag_vector.at(0) != pBar_tag_vector_background.at(0) )
        {
          //Minv pre-selection to boost usable ME statistics
          if( L_vector.at(0).M() > 1.105 && L_vector.at(0).M() < 1.13 && Lbar_vector_background.at(0).M() > 1.105 && Lbar_vector_background.at(0).M() < 1.13)
          {

            L_Lbar_L_vector_ME_SE_back_US_LS.push_back(L_vector.at(0));
            L_Lbar_p_vector_ME_SE_back_US_LS.push_back(p_vector.at(0));
            L_Lbar_pi_vector_ME_SE_back_US_LS.push_back(pi_vector.at(0));

            L_Lbar_L_pT_bin_vector_ME_SE_back_US_LS.push_back(L_pT_bin_vector.at(0));
            L_Lbar_L_eta_bin_vector_ME_SE_back_US_LS.push_back(L_eta_bin_vector.at(0));


            L_Lbar_Lbar_vector_ME_SE_back_US_LS.push_back(Lbar_vector_background.at(0));
            L_Lbar_pBar_vector_ME_SE_back_US_LS.push_back(pBar_vector_background.at(0));
            L_Lbar_piBar_vector_ME_SE_back_US_LS.push_back(piBar_vector_background.at(0));

            L_Lbar_Lbar_pT_bin_vector_ME_SE_back_US_LS.push_back(Lbar_pT_bin_vector_background.at(0));
            L_Lbar_Lbar_eta_bin_vector_ME_SE_back_US_LS.push_back(Lbar_eta_bin_vector_background.at(0));

          }
        }

      }


      //L(LS)-Lbar(US)
      if( L_vector_background.size() > 0 && Lbar_vector.size() > 0 && L_Lbar_L_vector_ME_SE_back_LS_US.size() < 3e4)
      {
        //L_LS(p pi+) Lbar(pBar pi+)
        //if(fabs(pi_vector_background.at(0).Eta() - piBar_vector.at(0).Eta()) > 1e-5 && fabs(p_vector_background.at(0).Eta() - piBar_vector.at(0).Eta()) > 1e-5 )
        if( pi_tag_vector_background.at(0) != piBar_tag_vector.at(0) && p_tag_vector_background.at(0) != piBar_tag_vector.at(0) )
        {
          //Minv pre-selection to boost usable ME statistics
          if( L_vector_background.at(0).M() > 1.105 && L_vector_background.at(0).M() < 1.13 && Lbar_vector.at(0).M() > 1.105 && Lbar_vector.at(0).M() < 1.13)
          {
            L_Lbar_L_vector_ME_SE_back_LS_US.push_back(L_vector_background.at(0));
            L_Lbar_p_vector_ME_SE_back_LS_US.push_back(p_vector_background.at(0));
            L_Lbar_pi_vector_ME_SE_back_LS_US.push_back(pi_vector_background.at(0));

            L_Lbar_L_pT_bin_vector_ME_SE_back_LS_US.push_back(L_pT_bin_vector_background.at(0));
            L_Lbar_L_eta_bin_vector_ME_SE_back_LS_US.push_back(L_eta_bin_vector_background.at(0));


            L_Lbar_Lbar_vector_ME_SE_back_LS_US.push_back(Lbar_vector.at(0));
            L_Lbar_pBar_vector_ME_SE_back_LS_US.push_back(pBar_vector.at(0));
            L_Lbar_piBar_vector_ME_SE_back_LS_US.push_back(piBar_vector.at(0));

            L_Lbar_Lbar_pT_bin_vector_ME_SE_back_LS_US.push_back(Lbar_pT_bin_vector.at(0));
            L_Lbar_Lbar_eta_bin_vector_ME_SE_back_LS_US.push_back(Lbar_eta_bin_vector.at(0));
          }
        }
      }


      //----------------------------------------------------------------------------------

      //same event L-L
      if( L_vector.size() > 0 && L_vector_background.size() > 0 && L_L_L1_vector_ME_SE_back.size() < 3e4)
      {
        //if(fabs(p_vector.at(0).Eta() - p_vector_background.at(0).Eta()) > 1e-5 && fabs(p_vector.at(0).Eta() - pi_vector_background.at(0).Eta()) > 1e-5 )
        if( p_tag_vector.at(0) != p_tag_vector_background.at(0) && p_tag_vector.at(0) != pi_tag_vector_background.at(0) )
        {
          //Minv pre-selection to boost usable ME statistics
          if( L_vector.at(0).M() > 1.105 && L_vector.at(0).M() < 1.13 && L_vector_background.at(0).M() > 1.105 && L_vector_background.at(0).M() < 1.13)
          {
            L_L_L1_vector_ME_SE_back.push_back(L_vector.at(0));
            L_L_p1_vector_ME_SE_back.push_back(p_vector.at(0));
            L_L_pi1_vector_ME_SE_back.push_back(pi_vector.at(0));

            L_L_L1_pT_bin_vector_ME_SE_back.push_back(L_pT_bin_vector.at(0));
            L_L_L1_eta_bin_vector_ME_SE_back.push_back(L_eta_bin_vector.at(0));


            L_L_L2_vector_ME_SE_back.push_back(L_vector_background.at(0));
            L_L_p2_vector_ME_SE_back.push_back(p_vector_background.at(0));
            L_L_pi2_vector_ME_SE_back.push_back(pi_vector_background.at(0));

            L_L_L2_pT_bin_vector_ME_SE_back.push_back(L_pT_bin_vector_background.at(0));
            L_L_L2_eta_bin_vector_ME_SE_back.push_back(L_eta_bin_vector_background.at(0));
          }
        }
      }

      //same event Lbar-Lbar
      if(Lbar_vector.size() > 0 && Lbar_vector_background.size() > 0 && Lbar_Lbar_Lbar1_vector_ME_SE_back.size() < 3e4)
      {
        //if(fabs(pBar_vector.at(0).Eta() - pBar_vector_background.at(0).Eta()) > 1e-5 && fabs(pBar_vector.at(0).Eta() - piBar_vector_background.at(0).Eta()) > 1e-5 )
        if( pBar_tag_vector.at(0) != pBar_tag_vector_background.at(0) && pBar_tag_vector.at(0) != piBar_tag_vector_background.at(0) )
        {
          //Minv pre-selection to boost usable ME statistics
          if( Lbar_vector.at(0).M() > 1.105 && Lbar_vector.at(0).M() < 1.13 && Lbar_vector_background.at(0).M() > 1.105 && Lbar_vector_background.at(0).M() < 1.13)
          {
            Lbar_Lbar_Lbar1_vector_ME_SE_back.push_back(Lbar_vector.at(0));
            Lbar_Lbar_pBar1_vector_ME_SE_back.push_back(pBar_vector.at(0));
            Lbar_Lbar_piBar1_vector_ME_SE_back.push_back(piBar_vector.at(0));

            Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.push_back(Lbar_pT_bin_vector.at(0));
            Lbar_Lbar_Lbar1_eta_bin_vector_ME_SE_back.push_back(Lbar_eta_bin_vector.at(0));


            Lbar_Lbar_Lbar2_vector_ME_SE_back.push_back(Lbar_vector_background.at(0));
            Lbar_Lbar_pBar2_vector_ME_SE_back.push_back(pBar_vector_background.at(0));
            Lbar_Lbar_piBar2_vector_ME_SE_back.push_back(piBar_vector_background.at(0));

            Lbar_Lbar_Lbar2_pT_bin_vector_ME_SE_back.push_back(Lbar_pT_bin_vector_background.at(0));
            Lbar_Lbar_Lbar2_eta_bin_vector_ME_SE_back.push_back(Lbar_eta_bin_vector_background.at(0));
          }
        }
      }


      //these LS pairs will be matched with SE US pairs to get ME for background
      if( L_vector_background.size() == 1 && Lbar_vector_background.size() == 0 && L_vector_ME_LS.size() < 1e5)
      {
        if( L_vector_background.at(0).M() > 1.105 && L_vector_background.at(0).M() < 1.13)
        {
          L_vector_ME_LS.push_back(L_vector_background.at(0));
          p_vector_ME_LS.push_back(p_vector_background.at(0));
          pi_vector_ME_LS.push_back(pi_vector_background.at(0));

          L_pT_bin_vector_ME_LS.push_back(L_pT_bin_vector_background.at(0));
          L_eta_bin_vector_ME_LS.push_back(L_eta_bin_vector_background.at(0));
        }
      }

      if( L_vector_background.size() == 0 && Lbar_vector_background.size() == 1 && Lbar_vector_ME_LS.size() < 1e5)
      {
        if( Lbar_vector_background.at(0).M() > 1.105 && Lbar_vector_background.at(0).M() < 1.13)
        {
          Lbar_vector_ME_LS.push_back(Lbar_vector_background.at(0));
          pBar_vector_ME_LS.push_back(pBar_vector_background.at(0));
          piBar_vector_ME_LS.push_back(piBar_vector_background.at(0));

          Lbar_pT_bin_vector_ME_LS.push_back(Lbar_pT_bin_vector_background.at(0));
          Lbar_eta_bin_vector_ME_LS.push_back(Lbar_eta_bin_vector_background.at(0));
        }
      }
      //_________________________________________________________________________________________________________



      //clear all vectors and fill new event

      //clear US vectors
      L_vector.clear();
      L_vector_pi_mass.clear();
      L_pT_bin_vector.clear();
      L_eta_bin_vector.clear();

      pi_vector.clear();
      pi_tag_vector.clear();

      p_vector.clear();
      p_tag_vector.clear();

      //----------------------------------

      Lbar_vector.clear();
      Lbar_vector_pi_mass.clear();
      Lbar_pT_bin_vector.clear();
      Lbar_eta_bin_vector.clear();

      piBar_vector.clear();
      piBar_tag_vector.clear();

      pBar_vector.clear();
      pBar_tag_vector.clear();

      //fill appropriate vectors from the new event
      if(charge == 0) //US charage combinations
      {
        //fill L or Lbar from new event
        //need to check cuts again in the new event
        if( cuts(year, L_y, pi_hasTOFinfo, p_dca, pi_dca, L_DCAdaughters, L_decayL, cos(L_theta)) && pT_bin_corr != -1 && eta_bin != -1 &&
            ( cut_type != 1 || cuts_topo_sys_err(p_dca, pi_dca, L_DCAdaughters, L_decayL, cos(L_theta) ) ) &&
            ( cut_type != 2 || cuts_pT_sys_err(p_pt, pi_pt) ) )
        {
          if( p_ch == 1 )
          {
            //QA histograms
            L_Minv_US[pT_bin]->Fill(L_mass);

            L_pT_vs_pi_pT_US->Fill(L_Lorentz_vector.Pt(), pi_Lorenz_vector.Pt());

            //fill vectors
            L_vector.push_back(L_Lorentz_vector);
            L_vector_pi_mass.push_back(L_Lorentz_vector_pi_mass);
            L_pT_bin_vector.push_back(pT_bin_corr);
            L_eta_bin_vector.push_back(eta_bin);

            pi_vector.push_back(pi_Lorenz_vector);
            pi_tag_vector.push_back(pi_InEventID);

            p_vector.push_back(p_Lorenz_vector);
            p_tag_vector.push_back(p_InEventID);

          }
          else if( p_ch == -1 )
          {
            //QA histograms
            Lbar_Minv_US[pT_bin]->Fill(L_mass);

            Lbar_pT_vs_pi_pT_US->Fill(L_Lorentz_vector.Pt(), pi_Lorenz_vector.Pt());

            //fill vectors
            Lbar_vector.push_back(L_Lorentz_vector);
            Lbar_vector_pi_mass.push_back(L_Lorentz_vector_pi_mass);
            Lbar_pT_bin_vector.push_back(pT_bin_corr);
            Lbar_eta_bin_vector.push_back(eta_bin);

            piBar_vector.push_back(pi_Lorenz_vector);
            piBar_tag_vector.push_back(pi_InEventID);

            pBar_vector.push_back(p_Lorenz_vector);
            pBar_tag_vector.push_back(p_InEventID);

          }

        } //end if cuts

      }//end if charge US

      //____________________________________________________________________________________________________________

      //clear LS vectors
      L_vector_background.clear();
      L_pT_bin_vector_background.clear();
      L_eta_bin_vector_background.clear();

      pi_vector_background.clear();
      pi_tag_vector_background.clear();

      p_vector_background.clear();
      p_tag_vector_background.clear();


      Lbar_vector_background.clear();
      Lbar_pT_bin_vector_background.clear();
      Lbar_eta_bin_vector_background.clear();

      piBar_vector_background.clear();
      piBar_tag_vector_background.clear();

      pBar_vector_background.clear();
      pBar_tag_vector_background.clear();


      if(charge != 0 ) //LS
      {
        //cuts
        if( cuts(year, L_y, pi_hasTOFinfo, p_dca, pi_dca, L_DCAdaughters, L_decayL, cos(L_theta)) && pT_bin_corr != -1 && eta_bin != -1 &&
          ( cut_type != 1 || cuts_topo_sys_err(p_dca, pi_dca, L_DCAdaughters, L_decayL, cos(L_theta) ) ) &&
          ( cut_type != 2 || cuts_pT_sys_err(p_pt, pi_pt) ) )
        {
          if( p_ch == 1 )
          {
            //QA histograms
            L_Minv_LS[pT_bin]->Fill(L_mass);

            L_pT_vs_pi_pT_LS->Fill(L_Lorentz_vector.Pt(), pi_Lorenz_vector.Pt());

            //fill vectors
            L_vector_background.push_back(L_Lorentz_vector);
            L_pT_bin_vector_background.push_back(pT_bin_corr);
            L_eta_bin_vector_background.push_back(eta_bin);

            pi_vector_background.push_back(pi_Lorenz_vector);
            pi_tag_vector_background.push_back(pi_InEventID);

            p_vector_background.push_back(p_Lorenz_vector);
            p_tag_vector_background.push_back(p_InEventID);

          }
          else if( p_ch == -1 )
          {
            //QA histograms
            Lbar_Minv_LS[pT_bin]->Fill(L_mass);

            Lbar_pT_vs_pi_pT_LS->Fill(L_Lorentz_vector.Pt(), pi_Lorenz_vector.Pt());

            //fill vectors
            Lbar_vector_background.push_back(L_Lorentz_vector);
            Lbar_pT_bin_vector_background.push_back(pT_bin_corr);
            Lbar_eta_bin_vector_background.push_back(eta_bin);

            piBar_vector_background.push_back(pi_Lorenz_vector);
            piBar_tag_vector_background.push_back(pi_InEventID);

            pBar_vector_background.push_back(p_Lorenz_vector);
            pBar_tag_vector_background.push_back(p_InEventID);

          }
        }
      }//end if LS

    }//end else for new event

  }//end loop over entries in NTuple

  cout<<"nEvenrs = "<<nEvents<<endl;


  //________________________________________________________________________________________________________

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(42);

  TFitResultPtr fit_res_gaus_L0_L0Bar[nPtBins_corr][nPtBins_corr];
  TFitResultPtr fit_res_gaus_L0_L0[nPtBins_corr][nPtBins_corr];
  TFitResultPtr fit_res_gaus_L0Bar_L0Bar[nPtBins_corr][nPtBins_corr];

  //for scaling of background to match (US-US) background levels
  float sideBandScale_L_Lbar[nPtBins_corr][nPtBins_corr];
  float sideBandScale_L_L[nPtBins_corr][nPtBins_corr];
  float sideBandScale_Lbar_Lbar[nPtBins_corr][nPtBins_corr];

  //to save fit parameters
  //L-Lbar
  float L_peak_mean_fit[nPtBins_corr][nPtBins_corr];
  float L_peak_sigma_fit[nPtBins_corr][nPtBins_corr];

  float Lbar_peak_mean_fit[nPtBins_corr][nPtBins_corr];
  float Lbar_peak_sigma_fit[nPtBins_corr][nPtBins_corr];

  //L-L
  float L1_peak_mean_fit[nPtBins_corr][nPtBins_corr];
  float L1_peak_sigma_fit[nPtBins_corr][nPtBins_corr];

  float L2_peak_mean_fit[nPtBins_corr][nPtBins_corr];
  float L2_peak_sigma_fit[nPtBins_corr][nPtBins_corr];

  //Lbar-Lbar
  float Lbar1_peak_mean_fit[nPtBins_corr][nPtBins_corr];
  float Lbar1_peak_sigma_fit[nPtBins_corr][nPtBins_corr];

  float Lbar2_peak_mean_fit[nPtBins_corr][nPtBins_corr];
  float Lbar2_peak_sigma_fit[nPtBins_corr][nPtBins_corr];

  //gives option to override fitting and set Minv range manually
  const int fit_Minv = 1;

  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      if( pTbin1 == 0 && pTbin2 == 0 )
      {
        L0_inv_mass_vs_L0bar_inv_mass_all_2_US = (TH2F*)L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Clone("L0_inv_mass_vs_L0bar_inv_mass_all_2_US");
        L0_inv_mass_vs_L0bar_inv_mass_all_2_LS = (TH2F*)L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Clone("L0_inv_mass_vs_L0bar_inv_mass_all_2_LS");
        L0_inv_mass_vs_L0bar_inv_mass_all_2_US_LS = (TH2F*)L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Clone("L0_inv_mass_vs_L0bar_inv_mass_all_2_US_LS");

        L0_inv_mass_vs_L0_inv_mass_all_2_US = (TH2F*)L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->Clone("L0_inv_mass_vs_L0_inv_mass_all_2_US");
        L0_inv_mass_vs_L0_inv_mass_all_2_LS = (TH2F*)L0_inv_mass_vs_L0_inv_mass_LS[pTbin1][pTbin2]->Clone("L0_inv_mass_vs_L0_inv_mass_all_2_LS");
        L0_inv_mass_vs_L0_inv_mass_all_2_US_LS = (TH2F*)L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->Clone("L0_inv_mass_vs_L0_inv_mass_all_2_US_LS");

        L0bar_inv_mass_vs_L0bar_inv_mass_all_2_US = (TH2F*)L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Clone("L0bar_inv_mass_vs_L0bar_inv_mass_all_2_US");
        L0bar_inv_mass_vs_L0bar_inv_mass_all_2_LS = (TH2F*)L0bar_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Clone("L0bar_inv_mass_vs_L0bar_inv_mass_all_2_LS");
        L0bar_inv_mass_vs_L0bar_inv_mass_all_2_US_LS = (TH2F*)L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Clone("L0bar_inv_mass_vs_L0bar_inv_mass_all_2_US_LS");
      }
      else
      {
        L0_inv_mass_vs_L0bar_inv_mass_all_2_US->Add(L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2], 1);
        L0_inv_mass_vs_L0bar_inv_mass_all_2_LS->Add(L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2], 1);
        L0_inv_mass_vs_L0bar_inv_mass_all_2_US_LS->Add(L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2], 1);

        L0_inv_mass_vs_L0_inv_mass_all_2_US->Add(L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2], 1);
        L0_inv_mass_vs_L0_inv_mass_all_2_LS->Add(L0_inv_mass_vs_L0_inv_mass_LS[pTbin1][pTbin2], 1);
        L0_inv_mass_vs_L0_inv_mass_all_2_US_LS->Add(L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2], 1);

        L0bar_inv_mass_vs_L0bar_inv_mass_all_2_US->Add(L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2], 1);
        L0bar_inv_mass_vs_L0bar_inv_mass_all_2_LS->Add(L0bar_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2], 1);
        L0bar_inv_mass_vs_L0bar_inv_mass_all_2_US_LS->Add(L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2], 1);
      }

      //find bins for side band and for projection to x
      //projection bins just for testing, later can do projections with parameters of 2D fit
      int binLow = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(1.11);
      int binHigh = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(1.12);

      int binLow_sideBand = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(1.14);
      int binHigh_sideBand = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(1.18);

      //scale background to match US-US
      sideBandScale_L_Lbar[pTbin1][pTbin2] = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh)/L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh);

      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Sumw2();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Scale(sideBandScale_L_Lbar[pTbin1][pTbin2]);

      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2] = (TH2F*)L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Clone(Form("L0_inv_mass_vs_L0bar_inv_mass_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Add(L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2], -1);

      //------------------------------------------------

      //scale LS to match US
      sideBandScale_L_L[pTbin1][pTbin2] = L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh)/L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh);

      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->Sumw2();
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->Scale(sideBandScale_L_L[pTbin1][pTbin2]);

      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2] = (TH2F*)L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->Clone(Form("L0_inv_mass_vs_L0_inv_mass_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->Add(L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2], -1);

      //------------------------------------------------

      //scale LS to match US
      sideBandScale_Lbar_Lbar[pTbin1][pTbin2] = L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Integral(binLow,binHigh, binLow_sideBand,binHigh_sideBand)/L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Integral( binLow,binHigh, binLow_sideBand,binHigh_sideBand);

      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Sumw2();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Scale(sideBandScale_Lbar_Lbar[pTbin1][pTbin2]);

      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2] = (TH2F*)L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Clone(Form("L0bar_inv_mass_vs_L0bar_inv_mass_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Add(L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2], -1);


      if(fit_Minv == 1) //gives option to override fitting and set Minv range manually
      {

        //L-Lbar
        TF2 *doubleGauss_L_Lbar = new TF2("doubleGauss_L_Lbar", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 1.07, 1.2, 1.07, 1.2);

        if( pTbin1 == 0 && pTbin2 == 0 )
        {
          if( energy == 200 ) doubleGauss_L_Lbar->SetParameters(1200, 1.116, 0.002, 1.116, 0.002);
          if( energy == 510 ) doubleGauss_L_Lbar->SetParameters(3000, 1.116, 0.002, 1.116, 0.002);
        }
        else if( pTbin1 > 0 && fit_res_gaus_L0_L0Bar[pTbin1-1][0]->IsValid() ) doubleGauss_L_Lbar->SetParameters(fit_res_gaus_L0_L0Bar[pTbin1-1][0]->Parameter(0)/10, fit_res_gaus_L0_L0Bar[pTbin1-1][0]->Parameter(1), fit_res_gaus_L0_L0Bar[pTbin1-1][0]->Parameter(2), fit_res_gaus_L0_L0Bar[pTbin1-1][0]->Parameter(3), fit_res_gaus_L0_L0Bar[pTbin1-1][0]->Parameter(4));
        else if( pTbin2 > 0 && fit_res_gaus_L0_L0Bar[0][pTbin2-1]->IsValid() ) doubleGauss_L_Lbar->SetParameters(fit_res_gaus_L0_L0Bar[0][pTbin2-1]->Parameter(0)/10, fit_res_gaus_L0_L0Bar[0][pTbin2-1]->Parameter(1), fit_res_gaus_L0_L0Bar[0][pTbin2-1]->Parameter(2), fit_res_gaus_L0_L0Bar[0][pTbin2-1]->Parameter(3), fit_res_gaus_L0_L0Bar[0][pTbin2-1]->Parameter(4));
        else
        {
          cout<<"Fit not valid for L-Lbar pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
          return false;
        }


        fit_res_gaus_L0_L0Bar[pTbin1][pTbin2] = L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_L_Lbar, "s 0", "", 1.11, 1.125);

        if( !fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->IsValid() )
        {
          cout<<"Fit not valid for L-Lbar pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
          return false;
        }


        //________________________________________________________________________________________________________

        //L-L
        TF2 *doubleGauss_L_L = new TF2("doubleGauss_L_L", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 1.07, 1.2, 1.07, 1.2);

        if( pTbin1 == 0 && pTbin2 == 0 )
        {
          if( energy == 200 ) doubleGauss_L_L->SetParameters(1000, 1.116, 0.002, 1.116, 0.002);
          if( energy == 510 ) doubleGauss_L_L->SetParameters(1000, 1.116, 0.002, 1.116, 0.002);
        }
        else if( pTbin1 > 0 && fit_res_gaus_L0_L0[pTbin1-1][0]->IsValid() ) doubleGauss_L_L->SetParameters(fit_res_gaus_L0_L0[pTbin1-1][0]->Parameter(0)/10, fit_res_gaus_L0_L0[pTbin1-1][0]->Parameter(1), fit_res_gaus_L0_L0[pTbin1-1][0]->Parameter(2), fit_res_gaus_L0_L0[pTbin1-1][0]->Parameter(3), fit_res_gaus_L0_L0[pTbin1-1][0]->Parameter(4));
        else if( pTbin2 > 0 && fit_res_gaus_L0_L0[0][pTbin2-1]->IsValid() ) doubleGauss_L_L->SetParameters(fit_res_gaus_L0_L0[0][pTbin2-1]->Parameter(0)/10, fit_res_gaus_L0_L0[0][pTbin2-1]->Parameter(1), fit_res_gaus_L0_L0[0][pTbin2-1]->Parameter(2), fit_res_gaus_L0_L0[0][pTbin2-1]->Parameter(3), fit_res_gaus_L0_L0[0][pTbin2-1]->Parameter(4));
        else
        {
          cout<<"Fit not valid for L-L pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
          return false;
        }

        fit_res_gaus_L0_L0[pTbin1][pTbin2] = L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_L_L, "s 0", "", 1.11, 1.125);

        if( !fit_res_gaus_L0_L0[pTbin1][pTbin2]->IsValid() )
        {
          cout<<"Fit not valid for L-L pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
          return false;
        }

        //________________________________________________________________________________________________________

        //Lbar-Lbar
        TF2 *doubleGauss_Lbar_Lbar = new TF2("doubleGauss_Lbar_Lbar", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 1.07, 1.2, 1.07, 1.2);

        if( pTbin1 == 0 && pTbin2 == 0 )
        {
          if( energy == 200 && year == 2012 )
          {
            if(cut_type == 2) doubleGauss_Lbar_Lbar->SetParameters(200, 1.116, 0.002, 1.116, 0.002);
            else doubleGauss_Lbar_Lbar->SetParameters(1000, 1.116, 0.002, 1.116, 0.002);
          }
          if( energy == 200 && year == 2015 ) doubleGauss_Lbar_Lbar->SetParameters(200, 1.116, 0.002, 1.116, 0.002);
          if( energy == 200 && year == 2016 ) doubleGauss_Lbar_Lbar->SetParameters(600, 1.116, 0.002, 1.116, 0.002);
          if( energy == 510 ) doubleGauss_Lbar_Lbar->SetParameters(600, 1.116, 0.002, 1.116, 0.002);
        }
        else if( pTbin1 > 0 && fit_res_gaus_L0Bar_L0Bar[pTbin1-1][0]->IsValid() ) doubleGauss_Lbar_Lbar->SetParameters(fit_res_gaus_L0Bar_L0Bar[pTbin1-1][0]->Parameter(0)/10, fit_res_gaus_L0Bar_L0Bar[pTbin1-1][0]->Parameter(1), fit_res_gaus_L0Bar_L0Bar[pTbin1-1][0]->Parameter(2), fit_res_gaus_L0Bar_L0Bar[pTbin1-1][0]->Parameter(3), fit_res_gaus_L0Bar_L0Bar[pTbin1-1][0]->Parameter(4));
        else if( pTbin2 > 0 && fit_res_gaus_L0Bar_L0Bar[0][pTbin2-1]->IsValid() ) doubleGauss_Lbar_Lbar->SetParameters(fit_res_gaus_L0Bar_L0Bar[0][pTbin2-1]->Parameter(0)/10, fit_res_gaus_L0Bar_L0Bar[0][pTbin2-1]->Parameter(1), fit_res_gaus_L0Bar_L0Bar[0][pTbin2-1]->Parameter(2), fit_res_gaus_L0Bar_L0Bar[0][pTbin2-1]->Parameter(3), fit_res_gaus_L0Bar_L0Bar[0][pTbin2-1]->Parameter(4));
        else
        {
          cout<<"Fit not valid for Lbar-Lbar pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
          return false;
        }

        fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2] = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_Lbar_Lbar, "s 0", "", 1.113, 1.125);

        if( !fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->IsValid() )
        {
          cout<<"Fit not valid for Lbar-Lbar pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
          return false;
        }

        //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        //save fit parameters

        //L-Lbar
        L_peak_mean_fit[pTbin1][pTbin2] = fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(1);
        L_peak_sigma_fit[pTbin1][pTbin2] = fabs(fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(2));

        Lbar_peak_mean_fit[pTbin1][pTbin2] = fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(3);
        Lbar_peak_sigma_fit[pTbin1][pTbin2] = fabs(fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(4));

        //L-L
        L1_peak_mean_fit[pTbin1][pTbin2] = fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(1);
        L1_peak_sigma_fit[pTbin1][pTbin2] = fabs(fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(2));

        L2_peak_mean_fit[pTbin1][pTbin2] = fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(3);
        L2_peak_sigma_fit[pTbin1][pTbin2] = fabs(fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(4));

        //Lbar-Lbar
        Lbar1_peak_mean_fit[pTbin1][pTbin2] = fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(1);
        Lbar1_peak_sigma_fit[pTbin1][pTbin2] = fabs(fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(2));

        Lbar2_peak_mean_fit[pTbin1][pTbin2] = fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(3);
        Lbar2_peak_sigma_fit[pTbin1][pTbin2] = fabs(fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(4));

      } //end if fit_Minv
      else //set parameters manually
      {
        //L-Lbar
        L_peak_mean_fit[pTbin1][pTbin2] = L_mass_PDG;
        L_peak_sigma_fit[pTbin1][pTbin2] = 2e-3;

        Lbar_peak_mean_fit[pTbin1][pTbin2] = L_mass_PDG;
        Lbar_peak_sigma_fit[pTbin1][pTbin2] = 2e-3;

        //L-L
        L1_peak_mean_fit[pTbin1][pTbin2] = L_mass_PDG;
        L1_peak_sigma_fit[pTbin1][pTbin2] = 2e-3;

        L2_peak_mean_fit[pTbin1][pTbin2] = L_mass_PDG;
        L2_peak_sigma_fit[pTbin1][pTbin2] = 2e-3;

        //Lbar-Lbar
        Lbar1_peak_mean_fit[pTbin1][pTbin2] = L_mass_PDG;
        Lbar1_peak_sigma_fit[pTbin1][pTbin2] = 2e-3;

        Lbar2_peak_mean_fit[pTbin1][pTbin2] = L_mass_PDG;
        Lbar2_peak_sigma_fit[pTbin1][pTbin2] = 2e-3;
      }


    }//end for pT bin 2
  }//end for pT bin 1


  InvMassFile->Write();
  InvMassFile->Close();

  return false; //just for debugging of fitting


  //________________________________________________________________________________________________________________________________________________________________________________

  //analyze stored Lambda pairs and save cos(theta*) histograms

  //TFile *Delta_phi_weight_file = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Weight/output_Lambda_pp_200_MB_1B_events_hists_100_SE_in_ME.root", "read");

  //TH1F *Delta_phi_reweight_hist = (TH1F*)Delta_phi_weight_file->Get("Delta_phi_reweight");

  //------------------------------------------

  TFile *LLbarOutFile; //output file to store production plane histograms;

  if(cut_type == 0)
  {
    LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_work.root", year), "recreate"); //analysis cuts
  }
  else if(cut_type == 1)
  {
    LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_topo_tight_cuts_work.root", year), "recreate"); //analysis cuts
  }
  else if(cut_type == 2)
  {
    LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_pT_tight_cuts_work.root", year), "recreate"); //analysis cuts
  }
  else
  {
    cout<<"Wrong cut type"<<endl;

    return false;
  }

  LLbarOutFile->cd();

  //QA and ME reweight histograms

  //signal+background
  TH2F *L0_L0bar_eta1_vs_eta2_US_hist = new TH2F("L0_L0bar_eta1_vs_eta2_US_hist", "L0_L0bar_eta1_vs_eta2_US_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_US_hist = new TH2F("L0_L0bar_phi1_vs_phi2_US_hist", "L0_L0bar_phi1_vs_phi2_US_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_US_hist = new TH2F("L0_L0bar_pT1_vs_pT2_US_hist", "L0_L0bar_pT1_vs_pT2_US_hist", 20, 0, 5, 20, 0, 5);

  TH1F *L0_L0bar_delta_phi_US_hist = new TH1F("L0_L0bar_delta_phi_US_hist", "L0_L0bar_delta_phi_US_hist", 100, 0, TMath::Pi());
  TH1F *L0_L0bar_delta_eta_US_hist = new TH1F("L0_L0bar_delta_eta_US_hist", "L0_L0bar_delta_eta_US_hist", 100, 0, 2);
  TH1F *L0_L0bar_delta_R_US_hist = new TH1F("L0_L0bar_delta_R_US_hist", "L0_L0bar_delta_R_US_hist", 100, 0, 4);
  TH1F *L0_L0bar_Q_US_hist = new TH1F("L0_L0bar_Q_US_hist", "L0_L0bar_Q_US_hist", 100, 0, 10);

  TH2F* L0_L0bar_delta_R_vs_Q_US_hist = new TH2F("L0_L0bar_delta_R_vs_Q_US_hist", "L0_L0bar_delta_R_vs_Q_US_hist", 100, 0, 4, 100, 0, 10 );


  TH2F *L0_L0bar_delta_phi_vs_delta_eta_US_hist = new TH2F("L0_L0bar_delta_phi_vs_delta_eta_US_hist", "L0_L0bar_delta_phi_vs_delta_eta_US_hist", 100, 0, TMath::Pi(), 100, 0, 2);

  //for J/Psi Minv
  TH2F *L0_L0bar_Minv_vs_delta_phi_US_hist = new TH2F("L0_L0bar_Minv_vs_delta_phi_US_hist", "L0_L0bar_Minv_vs_delta_phi_US_hist", 20, 2.5, 3.5, 3, 0, TMath::Pi() ); //|Delta phi| < cut
  TH2F *L0_L0bar_Minv_vs_delta_phi_2_US_hist = new TH2F("L0_L0bar_Minv_vs_delta_phi_2_US_hist", "L0_L0bar_Minv_vs_delta_phi_2_US_hist", 20, 2.5, 3.5, 3, 0, TMath::Pi() ); //|Delta phi| > cut

  TH2F *L0_L0bar_Minv_vs_delta_eta_US_hist = new TH2F("L0_L0bar_Minv_vs_delta_eta_US_hist", "L0_L0bar_Minv_vs_delta_eta_US_hist", 20, 2.5, 3.5, 2, 0, 2 ); //|Delta eta|

  //background
  //total - for QA
  TH2F *L0_L0bar_eta1_vs_eta2_US_LS_hist = new TH2F("L0_L0bar_eta1_vs_eta2_US_LS_hist", "L0_L0bar_eta1_vs_eta2_US_LS_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_US_LS_hist = new TH2F("L0_L0bar_phi1_vs_phi2_US_LS_hist", "L0_L0bar_phi1_vs_phi2_US_LS_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_US_LS_hist = new TH2F("L0_L0bar_pT1_vs_pT2_US_LS_hist", "L0_L0bar_pT1_vs_pT2_US_LS_hist", 20, 0, 5, 20, 0, 5);

  TH1F *L0_L0bar_delta_phi_US_LS_hist = new TH1F("L0_L0bar_delta_phi_US_LS_hist", "L0_L0bar_delta_phi_US_LS_hist", 100, 0, TMath::Pi());
  TH1F *L0_L0bar_delta_eta_US_LS_hist = new TH1F("L0_L0bar_delta_eta_US_LS_hist", "L0_L0bar_delta_eta_US_LS_hist", 100, 0, 2);
  TH1F *L0_L0bar_delta_R_US_LS_hist = new TH1F("L0_L0bar_delta_R_US_LS_hist", "L0_L0bar_delta_R_US_LS_hist", 100, 0, 4);
  TH1F *L0_L0bar_Q_US_LS_hist = new TH1F("L0_L0bar_Q_US_LS_hist", "L0_L0bar_Q_US_LS_hist", 100, 0, 10);

  TH2F* L0_L0bar_delta_R_vs_Q_US_LS_hist = new TH2F("L0_L0bar_delta_R_vs_Q_US_LS_hist", "L0_L0bar_delta_R_vs_Q_US_LS_hist", 100, 0, 4, 100, 0, 10 );

  TH2F *L0_L0bar_delta_phi_vs_delta_eta_US_LS_hist = new TH2F("L0_L0bar_delta_phi_vs_delta_eta_US_LS_hist", "L0_L0bar_delta_phi_vs_delta_eta_US_LS_hist", 100, 0, TMath::Pi(), 100, 0, 2);

  //for J/Psi Minv
  TH2F *L0_L0bar_Minv_vs_delta_phi_US_LS_hist = new TH2F("L0_L0bar_Minv_vs_delta_phi_US_LS_hist", "L0_L0bar_Minv_vs_delta_phi_US_LS_hist", 20, 2.5, 3.5, 3, 0, TMath::Pi() ); //|Delta phi| < cut
  TH2F *L0_L0bar_Minv_vs_delta_phi_2_US_LS_hist = new TH2F("L0_L0bar_Minv_vs_delta_phi_2_US_LS_hist", "L0_L0bar_Minv_vs_delta_phi_2_US_LS_hist", 20, 2.5, 3.5, 3, 0, TMath::Pi() ); //|Delta phi| > cut

  TH2F *L0_L0bar_Minv_vs_delta_eta_US_LS_hist = new TH2F("L0_L0bar_Minv_vs_delta_eta_US_LS_hist", "L0_L0bar_Minv_vs_delta_eta_US_LS_hist", 20, 2.5, 3.5, 2, 0, 2 ); //|Delta eta|

  //for ME reweight
  //L(US)-Lbar(LS)
  TH2F *L0_L0bar_eta1_vs_eta2_US_LS_1_hist = new TH2F("L0_L0bar_eta1_vs_eta2_US_LS_1_hist", "L0_L0bar_eta1_vs_eta2_US_LS_1_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_US_LS_1_hist = new TH2F("L0_L0bar_phi1_vs_phi2_US_LS_1_hist", "L0_L0bar_phi1_vs_phi2_US_LS_1_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_US_LS_1_hist = new TH2F("L0_L0bar_pT1_vs_pT2_US_LS_1_hist", "L0_L0bar_pT1_vs_pT2_US_LS_1_hist", 20, 0, 5, 20, 0, 5);

  //L(LS)-Lbar(US)
  TH2F *L0_L0bar_eta1_vs_eta2_US_LS_2_hist = new TH2F("L0_L0bar_eta1_vs_eta2_US_LS_2_hist", "L0_L0bar_eta1_vs_eta2_US_LS_2_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_US_LS_2_hist = new TH2F("L0_L0bar_phi1_vs_phi2_US_LS_2_hist", "L0_L0bar_phi1_vs_phi2_US_LS_2_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_US_LS_2_hist = new TH2F("L0_L0bar_pT1_vs_pT2_US_LS_2_hist", "L0_L0bar_pT1_vs_pT2_US_LS_2_hist", 20, 0, 5, 20, 0, 5);

  //------------------------------------------------------------------

  //daughter kinematics

  //US-US
  TH2F *L0_L0bar_p1_pT1_vs_p2_pT2_US_hist = new TH2F("L0_L0bar_p1_pT1_vs_p2_pT2_US_hist", "L0_L0bar_p1_pT1_vs_p2_pT2_US_hist", 100, 0, 5, 100, 0, 5);
  TH2F *L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist = new TH2F("L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist", "L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist", 100, 0, 2, 100, 0, 2);

  TH2F *L0_L0bar_p1_eta1_vs_p2_eta2_US_hist = new TH2F("L0_L0bar_p1_eta1_vs_p2_eta2_US_hist", "L0_L0bar_p1_eta1_vs_p2_eta2_US_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist = new TH2F("L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist", "L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist", 20, -1, 1, 20, -1, 1);

  TH2F *L0_L0bar_p1_phi1_vs_p2_phi2_US_hist = new TH2F("L0_L0bar_p1_phi1_vs_p2_phi2_US_hist", "L0_L0bar_p1_phi1_vs_p2_phi2_US_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist = new TH2F("L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist", "L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());

  //US-LS
  TH2F *L0_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist = new TH2F("L0_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist", "L0_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist", 100, 0, 5, 100, 0, 5);
  TH2F *L0_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist = new TH2F("L0_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist", "L0_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist", 100, 0, 2, 100, 0, 2);

  TH2F *L0_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist = new TH2F("L0_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist", "L0_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist = new TH2F("L0_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist", "L0_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist", 20, -1, 1, 20, -1, 1);

  TH2F *L0_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist = new TH2F("L0_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist", "L0_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist = new TH2F("L0_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist", "L0_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());

  //pi pT vs. cos(theta*)
  TH2F *L0_L0bar_pi1_pT_vs_cos_theta_star_US_hist = new TH2F("L0_L0bar_pi1_pT_vs_cos_theta_star_US_hist", "L0_L0bar_pi1_pT_vs_cos_theta_star_US_hist", 50, 0, 1, 10, -1, 1);
  TH2F *L0_L0bar_pi2_pT_vs_cos_theta_star_US_hist = new TH2F("L0_L0bar_pi2_pT_vs_cos_theta_star_US_hist", "L0_L0bar_pi2_pT_vs_cos_theta_star_US_hist", 50, 0, 1, 10, -1, 1);

  TH2F *L0_L0bar_pi1_pT_vs_cos_theta_star_US_LS_hist = new TH2F("L0_L0bar_pi1_pT_vs_cos_theta_star_US_LS_hist", "L0_L0bar_pi1_pT_vs_cos_theta_star_US_LS_hist", 50, 0, 1, 10, -1, 1);
  TH2F *L0_L0bar_pi2_pT_vs_cos_theta_star_US_LS_hist = new TH2F("L0_L0bar_pi2_pT_vs_cos_theta_star_US_LS_hist", "L0_L0bar_pi2_pT_vs_cos_theta_star_US_LS_hist", 50, 0, 1, 10, -1, 1);

  //------------------------------------------------------------------------------------------------------------------------------

  //signal+background
  TH2F *L0_L0_eta1_vs_eta2_US_hist = new TH2F("L0_L0_eta1_vs_eta2_US_hist", "L0_L0_eta1_vs_eta2_US_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0_phi1_vs_phi2_US_hist = new TH2F("L0_L0_phi1_vs_phi2_US_hist", "L0_L0_phi1_vs_phi2_US_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0_pT1_vs_pT2_US_hist = new TH2F("L0_L0_pT1_vs_pT2_US_hist", "L0_L0_pT1_vs_pT2_US_hist", 20, 0, 5, 20, 0, 5);

  TH1F *L0_L0_delta_phi_US_hist = new TH1F("L0_L0_delta_phi_US_hist", "L0_L0_delta_phi_US_hist", 100, 0, TMath::Pi());
  TH1F *L0_L0_delta_eta_US_hist = new TH1F("L0_L0_delta_eta_US_hist", "L0_L0_delta_eta_US_hist", 100, 0, 2);
  TH1F *L0_L0_delta_R_US_hist = new TH1F("L0_L0_delta_R_US_hist", "L0_L0_delta_R_US_hist", 100, 0, 4);
  TH1F *L0_L0_Q_US_hist = new TH1F("L0_L0_Q_US_hist", "L0_L0_Q_US_hist", 100, 0, 10);

  TH2F* L0_L0_delta_R_vs_Q_US_hist = new TH2F("L0_L0_delta_R_vs_Q_US_hist", "L0_L0_delta_R_vs_Q_US_hist", 100, 0, 4, 100, 0, 10 );

  TH2F *L0_L0_delta_phi_vs_delta_eta_US_hist = new TH2F("L0_L0_delta_phi_vs_delta_eta_US_hist", "L0_L0_delta_phi_vs_delta_eta_US_hist", 100, 0, TMath::Pi(), 100, 0, 2);

  //TH1F *L0_L0_cos_theta_star_vs_delta_phi_US_hist = new TH1F("L0_L0_cos_theta_star_vs_delta_phi_US_hist", "L0_L0_cos_theta_star_vs_delta_phi_US_hist", 10, -1, 1, 3, 0, TMath::Pi());

  //background
  //total - for QA
  TH2F *L0_L0_eta1_vs_eta2_US_LS_hist = new TH2F("L0_L0_eta1_vs_eta2_US_LS_hist", "L0_L0_eta1_vs_eta2_US_LS_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0_phi1_vs_phi2_US_LS_hist = new TH2F("L0_L0_phi1_vs_phi2_US_LS_hist", "L0_L0_phi1_vs_phi2_US_LS_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0_pT1_vs_pT2_US_LS_hist = new TH2F("L0_L0_pT1_vs_pT2_US_LS_hist", "L0_L0_pT1_vs_pT2_US_LS_hist", 20, 0, 5, 20, 0, 5);

  TH1F *L0_L0_delta_phi_US_LS_hist = new TH1F("L0_L0_delta_phi_US_LS_hist", "L0_L0_delta_phi_US_LS_hist", 100, 0, TMath::Pi());
  TH1F *L0_L0_delta_eta_US_LS_hist = new TH1F("L0_L0_delta_eta_US_LS_hist", "L0_L0_delta_eta_US_LS_hist", 100, 0, 2);
  TH1F *L0_L0_delta_R_US_LS_hist = new TH1F("L0_L0_delta_R_US_LS_hist", "L0_L0_delta_R_US_LS_hist", 100, 0, 4);
  TH1F *L0_L0_Q_US_LS_hist = new TH1F("L0_L0_Q_US_LS_hist", "L0_L0_Q_US_LS_hist", 100, 0, 10);

  TH2F* L0_L0_delta_R_vs_Q_US_LS_hist = new TH2F("L0_L0_delta_R_vs_Q_US_LS_hist", "L0_L0_delta_R_vs_Q_US_LS_hist", 100, 0, 4, 100, 0, 10 );

  TH2F *L0_L0_delta_phi_vs_delta_eta_US_LS_hist = new TH2F("L0_L0_delta_phi_vs_delta_eta_US_LS_hist", "L0_L0_delta_phi_vs_delta_eta_US_LS_hist", 100, 0, TMath::Pi(), 100, 0, 2);

  //TH1F *L0_L0_cos_theta_star_vs_delta_phi_US_LS_hist = new TH1F("L0_L0_cos_theta_star_vs_delta_phi_US_LS_hist", "L0_L0_cos_theta_star_vs_delta_phi_US_LS_hist", 10, -1, 1, 3, 0, TMath::Pi());

  //for ME reweight
  //L1(US)-L2(LS)
  TH2F *L0_L0_eta1_vs_eta2_US_LS_1_hist = new TH2F("L0_L0_eta1_vs_eta2_US_LS_1_hist", "L0_L0_eta1_vs_eta2_US_LS_1_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0_phi1_vs_phi2_US_LS_1_hist = new TH2F("L0_L0_phi1_vs_phi2_US_LS_1_hist", "L0_L0_phi1_vs_phi2_US_LS_1_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0_pT1_vs_pT2_US_LS_1_hist = new TH2F("L0_L0_pT1_vs_pT2_US_LS_1_hist", "L0_L0_pT1_vs_pT2_US_LS_1_hist", 20, 0, 5, 20, 0, 5);

  //L1(LS)-L2(US)
  TH2F *L0_L0_eta1_vs_eta2_US_LS_2_hist = new TH2F("L0_L0_eta1_vs_eta2_US_LS_2_hist", "L0_L0_eta1_vs_eta2_US_LS_2_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0_phi1_vs_phi2_US_LS_2_hist = new TH2F("L0_L0_phi1_vs_phi2_US_LS_2_hist", "L0_L0_phi1_vs_phi2_US_LS_2_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0_pT1_vs_pT2_US_LS_2_hist = new TH2F("L0_L0_pT1_vs_pT2_US_LS_2_hist", "L0_L0_pT1_vs_pT2_US_LS_2_hist", 20, 0, 5, 20, 0, 5);

  //----------------------------------------------------------------------------------------

  //daughter kinematics
  //US-US
  TH2F *L0_L0_p1_pT1_vs_p2_pT2_US_hist = new TH2F("L0_L0_p1_pT1_vs_p2_pT2_US_hist", "L0_L0_p1_pT1_vs_p2_pT2_US_hist", 100, 0, 5, 100, 0, 5);
  TH2F *L0_L0_pi1_pT1_vs_pi2_pT2_US_hist = new TH2F("L0_L0_pi1_pT1_vs_pi2_pT2_US_hist", "L0_L0_pi1_pT1_vs_pi2_pT2_US_hist", 100, 0, 2, 100, 0, 2);

  TH2F *L0_L0_p1_eta1_vs_p2_eta2_US_hist = new TH2F("L0_L0_p1_eta1_vs_p2_eta2_US_hist", "L0_L0_p1_eta1_vs_p2_eta2_US_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0_pi1_eta1_vs_pi2_eta2_US_hist = new TH2F("L0_L0_pi1_eta1_vs_pi2_eta2_US_hist", "L0_L0_pi1_eta1_vs_pi2_eta2_US_hist", 20, -1, 1, 20, -1, 1);

  TH2F *L0_L0_p1_phi1_vs_p2_phi2_US_hist = new TH2F("L0_L0_p1_phi1_vs_p2_phi2_US_hist", "L0_L0_p1_phi1_vs_p2_phi2_US_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0_pi1_phi1_vs_pi2_phi2_US_hist = new TH2F("L0_L0_pi1_phi1_vs_pi2_phi2_US_hist", "L0_L0_pi1_phi1_vs_pi2_phi2_US_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());

  //US-LS
  TH2F *L0_L0_p1_pT1_vs_p2_pT2_US_LS_hist = new TH2F("L0_L0_p1_pT1_vs_p2_pT2_US_LS_hist", "L0_L0_p1_pT1_vs_p2_pT2_US_LS_hist", 100, 0, 5, 100, 0, 5);
  TH2F *L0_L0_pi1_pT1_vs_pi2_pT2_US_LS_hist = new TH2F("L0_L0_pi1_pT1_vs_pi2_pT2_US_LS_hist", "L0_L0_pi1_pT1_vs_pi2_pT2_US_LS_hist", 100, 0, 2, 100, 0, 2);

  TH2F *L0_L0_p1_eta1_vs_p2_eta2_US_LS_hist = new TH2F("L0_L0_p1_eta1_vs_p2_eta2_US_LS_hist", "L0_L0_p1_eta1_vs_p2_eta2_US_LS_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0_pi1_eta1_vs_pi2_eta2_US_LS_hist = new TH2F("L0_L0_pi1_eta1_vs_pi2_eta2_US_LS_hist", "L0_L0_pi1_eta1_vs_pi2_eta2_US_LS_hist", 20, -1, 1, 20, -1, 1);

  TH2F *L0_L0_p1_phi1_vs_p2_phi2_US_LS_hist = new TH2F("L0_L0_p1_phi1_vs_p2_phi2_US_LS_hist", "L0_L0_p1_phi1_vs_p2_phi2_US_LS_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0_pi1_phi1_vs_pi2_phi2_US_LS_hist = new TH2F("L0_L0_pi1_phi1_vs_pi2_phi2_US_LS_hist", "L0_L0_pi1_phi1_vs_pi2_phi2_US_LS_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());

  //------------------------------------------------------------------------------------------------------------------------------

  //signal+background
  TH2F *L0bar_L0bar_eta1_vs_eta2_US_hist = new TH2F("L0bar_L0bar_eta1_vs_eta2_US_hist", "L0bar_L0bar_eta1_vs_eta2_US_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0bar_L0bar_phi1_vs_phi2_US_hist = new TH2F("L0bar_L0bar_phi1_vs_phi2_US_hist", "L0bar_L0bar_phi1_vs_phi2_US_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0bar_L0bar_pT1_vs_pT2_US_hist = new TH2F("L0bar_L0bar_pT1_vs_pT2_US_hist", "L0bar_L0bar_pT1_vs_pT2_US_hist", 20, 0, 5, 20, 0, 5);

  TH1F *L0bar_L0bar_delta_phi_US_hist = new TH1F("L0bar_L0bar_delta_phi_US_hist", "L0bar_L0bar_delta_phi_US_hist", 100, 0, TMath::Pi());
  TH1F *L0bar_L0bar_delta_eta_US_hist = new TH1F("L0bar_L0bar_delta_eta_US_hist", "L0bar_L0bar_delta_eta_US_hist", 100, 0, 2);
  TH1F *L0bar_L0bar_delta_R_US_hist = new TH1F("L0bar_L0bar_delta_R_US_hist", "L0bar_L0bar_delta_R_US_hist", 100, 0, 4);
  TH1F *L0bar_L0bar_Q_US_hist = new TH1F("L0bar_L0bar_Q_US_hist", "L0bar_L0bar_Q_US_hist", 100, 0, 10);

  TH2F* L0bar_L0bar_delta_R_vs_Q_US_hist = new TH2F("L0bar_L0bar_delta_R_vs_Q_US_hist", "L0bar_L0bar_delta_R_vs_Q_US_hist", 100, 0, 4, 100, 0, 10 );


  TH2F *L0bar_L0bar_delta_phi_vs_delta_eta_US_hist = new TH2F("L0bar_L0bar_delta_phi_vs_delta_eta_US_hist", "L0bar_L0bar_delta_phi_vs_delta_eta_US_hist", 100, 0, TMath::Pi(), 100, 0, 2);

  //TH1F *L0bar_L0bar_cos_theta_star_vs_delta_phi_US_hist = new TH1F("L0bar_L0bar_cos_theta_star_vs_delta_phi_US_hist", "L0bar_L0bar_cos_theta_star_vs_delta_phi_US_hist", 10, -1, 1, 3, 0, TMath::Pi());

  //background
  //total - for QA
  TH2F *L0bar_L0bar_eta1_vs_eta2_US_LS_hist = new TH2F("L0bar_L0bar_eta1_vs_eta2_US_LS_hist", "L0bar_L0bar_eta1_vs_eta2_US_LS_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0bar_L0bar_phi1_vs_phi2_US_LS_hist = new TH2F("L0bar_L0bar_phi1_vs_phi2_US_LS_hist", "L0bar_L0bar_phi1_vs_phi2_US_LS_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0bar_L0bar_pT1_vs_pT2_US_LS_hist = new TH2F("L0bar_L0bar_pT1_vs_pT2_US_LS_hist", "L0bar_L0bar_pT1_vs_pT2_US_LS_hist", 20, 0, 5, 20, 0, 5);

  TH1F *L0bar_L0bar_delta_phi_US_LS_hist = new TH1F("L0bar_L0bar_delta_phi_US_LS_hist", "L0bar_L0bar_delta_phi_US_LS_hist", 100, 0, TMath::Pi());
  TH1F *L0bar_L0bar_delta_eta_US_LS_hist = new TH1F("L0bar_L0bar_delta_eta_US_LS_hist", "L0bar_L0bar_delta_eta_US_LS_hist", 100, 0, 2);
  TH1F *L0bar_L0bar_delta_R_US_LS_hist = new TH1F("L0bar_L0bar_delta_R_US_LS_hist", "L0bar_L0bar_delta_R_US_LS_hist", 100, 0, 4);
  TH1F *L0bar_L0bar_Q_US_LS_hist = new TH1F("L0bar_L0bar_Q_US_LS_hist", "L0bar_L0bar_Q_US_LS_hist", 100, 0, 10);

  TH2F* L0bar_L0bar_delta_R_vs_Q_US_LS_hist = new TH2F("L0bar_L0bar_delta_R_vs_Q_US_LS_hist", "L0bar_L0bar_delta_R_vs_Q_US_LS_hist", 100, 0, 4, 100, 0, 10 );

  TH2F *L0bar_L0bar_delta_phi_vs_delta_eta_US_LS_hist = new TH2F("L0bar_L0bar_delta_phi_vs_delta_eta_US_LS_hist", "L0bar_L0bar_delta_phi_vs_delta_eta_US_LS_hist", 100, 0, TMath::Pi(), 100, 0, 2);

  //TH1F *L0bar_L0bar_cos_theta_star_vs_delta_phi_US_LS_hist = new TH1F("L0bar_L0bar_cos_theta_star_vs_delta_phi_US_LS_hist", "L0bar_L0bar_cos_theta_star_vs_delta_phi_US_LS_hist", 10, -1, 1, 3, 0, TMath::Pi());

  //for ME reweight
  //Lbar1(US)-Lbar2(LS)
  TH2F *L0bar_L0bar_eta1_vs_eta2_US_LS_1_hist = new TH2F("L0bar_L0bar_eta1_vs_eta2_US_LS_1_hist", "L0bar_L0bar_eta1_vs_eta2_US_LS_1_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0bar_L0bar_phi1_vs_phi2_US_LS_1_hist = new TH2F("L0bar_L0bar_phi1_vs_phi2_US_LS_1_hist", "L0bar_L0bar_phi1_vs_phi2_US_LS_1_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0bar_L0bar_pT1_vs_pT2_US_LS_1_hist = new TH2F("L0bar_L0bar_pT1_vs_pT2_US_LS_1_hist", "L0bar_L0bar_pT1_vs_pT2_US_LS_1_hist", 20, 0, 5, 20, 0, 5);

  //Lbar1(US)-Lbar2(LS)
  TH2F *L0bar_L0bar_eta1_vs_eta2_US_LS_2_hist = new TH2F("L0bar_L0bar_eta1_vs_eta2_US_LS_2_hist", "L0bar_L0bar_eta1_vs_eta2_US_LS_2_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0bar_L0bar_phi1_vs_phi2_US_LS_2_hist = new TH2F("L0bar_L0bar_phi1_vs_phi2_US_LS_2_hist", "L0bar_L0bar_phi1_vs_phi2_US_LS_2_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0bar_L0bar_pT1_vs_pT2_US_LS_2_hist = new TH2F("L0bar_L0bar_pT1_vs_pT2_US_LS_2_hist", "L0bar_L0bar_pT1_vs_pT2_US_LS_2_hist", 20, 0, 5, 20, 0, 5);

  //----------------------------------------------------------------------------------------

  //daughter kinematics
  //US-US
  TH2F *L0bar_L0bar_p1_pT1_vs_p2_pT2_US_hist = new TH2F("L0bar_L0bar_p1_pT1_vs_p2_pT2_US_hist", "L0bar_L0bar_p1_pT1_vs_p2_pT2_US_hist", 100, 0, 5, 100, 0, 5);
  TH2F *L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_hist = new TH2F("L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_hist", "L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_hist", 100, 0, 2, 100, 0, 2);

  TH2F *L0bar_L0bar_p1_eta1_vs_p2_eta2_US_hist = new TH2F("L0bar_L0bar_p1_eta1_vs_p2_eta2_US_hist", "L0bar_L0bar_p1_eta1_vs_p2_eta2_US_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_hist = new TH2F("L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_hist", "L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_hist", 20, -1, 1, 20, -1, 1);

  TH2F *L0bar_L0bar_p1_phi1_vs_p2_phi2_US_hist = new TH2F("L0bar_L0bar_p1_phi1_vs_p2_phi2_US_hist", "L0bar_L0bar_p1_phi1_vs_p2_phi2_US_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_hist = new TH2F("L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_hist", "L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());

  //US-LS
  TH2F *L0bar_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist = new TH2F("L0bar_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist", "L0bar_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist", 100, 0, 5, 100, 0, 5);
  TH2F *L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist = new TH2F("L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist", "L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist", 100, 0, 2, 100, 0, 2);

  TH2F *L0bar_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist = new TH2F("L0bar_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist", "L0bar_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist = new TH2F("L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist", "L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist", 20, -1, 1, 20, -1, 1);

  TH2F *L0bar_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist = new TH2F("L0bar_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist", "L0bar_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist = new TH2F("L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist", "L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());

  //______________________________________________________________________________________________________________________________________________________________________________

  //mixed-event histograms
  //L-Lbar
  //signal+background
  TH2F *L0_L0bar_eta1_vs_eta2_US_ME_hist = new TH2F("L0_L0bar_eta1_vs_eta2_US_ME_hist", "L0_L0bar_eta1_vs_eta2_US_ME_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_US_ME_hist = new TH2F("L0_L0bar_phi1_vs_phi2_US_ME_hist", "L0_L0bar_phi1_vs_phi2_US_ME_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_US_ME_hist = new TH2F("L0_L0bar_pT1_vs_pT2_US_ME_hist", "L0_L0bar_pT1_vs_pT2_US_ME_hist", 20, 0, 5, 20, 0, 5);

  TH1F *L0_L0bar_delta_phi_US_ME_hist = new TH1F("L0_L0bar_delta_phi_US_ME_hist", "L0_L0bar_delta_phi_US_ME_hist", 100, 0, TMath::Pi());
  TH1F *L0_L0bar_delta_eta_US_ME_hist = new TH1F("L0_L0bar_delta_eta_US_ME_hist", "L0_L0bar_delta_eta_US_ME_hist", 100, 0, 2);

  TH2F *L0_L0bar_delta_phi_vs_delta_eta_US_ME_hist = new TH2F("L0_L0bar_delta_phi_vs_delta_eta_US_ME_hist", "L0_L0bar_delta_phi_vs_delta_eta_US_ME_hist", 100, 0, TMath::Pi(), 100, 0, 2);

  TH1F *L0_L0bar_delta_phi_US_ME_nFill_hist = new TH1F("L0_L0bar_delta_phi_US_ME_nFill_hist", "L0_L0bar_delta_phi_US_ME_nFill_hist", 100, 0, 100);
  TH1F *L0_L0bar_delta_phi_weighted_US_ME_nFill_hist = new TH1F("L0_L0bar_delta_phi_weighted_US_ME_nFill_hist", "L0_L0bar_delta_phi_weighted_US_ME_nFill_hist", 100, 0, 100);

  TH2F *L0_L0bar_delta_phi_vs_nFill_US_ME_hist = new TH2F("L0_L0bar_delta_phi_vs_nFill_US_ME_hist", "L0_L0bar_delta_phi_vs_nFill_US_ME_hist", 100, 0, TMath::Pi(), 100, 0, 100);

  //background
  //L(US)-Lbar(LS)
  TH2F *L0_L0bar_eta1_vs_eta2_US_LS_ME_1_hist = new TH2F("L0_L0bar_eta1_vs_eta2_US_LS_ME_1_hist", "L0_L0bar_eta1_vs_eta2_US_LS_ME_1_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_US_LS_ME_1_hist = new TH2F("L0_L0bar_phi1_vs_phi2_US_LS_ME_1_hist", "L0_L0bar_phi1_vs_phi2_US_LS_ME_1_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_US_LS_ME_1_hist = new TH2F("L0_L0bar_pT1_vs_pT2_US_LS_ME_1_hist", "L0_L0bar_pT1_vs_pT2_US_LS_ME_1_hist", 20, 0, 5, 20, 0, 5);

  //L(LS)-Lbar(US)
  TH2F *L0_L0bar_eta1_vs_eta2_US_LS_ME_2_hist = new TH2F("L0_L0bar_eta1_vs_eta2_US_LS_ME_2_hist", "L0_L0bar_eta1_vs_eta2_US_LS_ME_2_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0bar_phi1_vs_phi2_US_LS_ME_2_hist = new TH2F("L0_L0bar_phi1_vs_phi2_US_LS_ME_2_hist", "L0_L0bar_phi1_vs_phi2_US_LS_ME_2_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0bar_pT1_vs_pT2_US_LS_ME_2_hist = new TH2F("L0_L0bar_pT1_vs_pT2_US_LS_ME_2_hist", "L0_L0bar_pT1_vs_pT2_US_LS_ME_2_hist", 20, 0, 5, 20, 0, 5);

  TH1F *L0_L0bar_delta_phi_US_LS_ME_hist = new TH1F("L0_L0bar_delta_phi_US_LS_ME_hist", "L0_L0bar_delta_phi_US_LS_ME_hist", 100, 0, TMath::Pi());
  TH1F *L0_L0bar_delta_eta_US_LS_ME_hist = new TH1F("L0_L0bar_delta_eta_US_LS_ME_hist", "L0_L0bar_delta_eta_US_LS_ME_hist", 100, 0, 2);

  TH2F *L0_L0bar_delta_phi_vs_delta_eta_US_LS_ME_hist = new TH2F("L0_L0bar_delta_phi_vs_delta_eta_US_LS_ME_hist", "L0_L0bar_delta_phi_vs_delta_eta_US_LS_ME_hist", 100, 0, TMath::Pi(), 100, 0, 2);

  //------------------------------------------------------------------------------------------------------------------------------

  //pion kinematics
  TH2F *L0_L0bar_pi1_pT1_vs_pi2_pT2_US_ME_hist = new TH2F("L0_L0bar_pi1_pT1_vs_pi2_pT2_US_ME_hist", "L0_L0bar_pi1_pT1_vs_pi2_pT2_US_ME_hist", 100, 0, 2, 100, 0, 2);

  TH2F *L0_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_ME_hist = new TH2F("L0_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_ME_hist", "L0_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_ME_hist", 100, 0, 2, 100, 0, 2);


  //pi pT vs. cos(theta*)
  TH2F *L0_L0bar_pi1_pT_vs_cos_theta_star_US_ME_hist = new TH2F("L0_L0bar_pi1_pT_vs_cos_theta_star_US_ME_hist", "L0_L0bar_pi1_pT_vs_cos_theta_star_US_ME_hist", 50, 0, 1, 10, -1, 1);
  TH2F *L0_L0bar_pi2_pT_vs_cos_theta_star_US_ME_hist = new TH2F("L0_L0bar_pi2_pT_vs_cos_theta_star_US_ME_hist", "L0_L0bar_pi2_pT_vs_cos_theta_star_US_ME_hist", 50, 0, 1, 10, -1, 1);

  TH2F *L0_L0bar_pi1_pT_vs_cos_theta_star_US_LS_ME_hist = new TH2F("L0_L0bar_pi1_pT_vs_cos_theta_star_US_LS_ME_hist", "L0_L0bar_pi1_pT_vs_cos_theta_star_US_LS_ME_hist", 50, 0, 1, 10, -1, 1);
  TH2F *L0_L0bar_pi2_pT_vs_cos_theta_star_US_LS_ME_hist = new TH2F("L0_L0bar_pi2_pT_vs_cos_theta_star_US_LS_ME_hist", "L0_L0bar_pi2_pT_vs_cos_theta_star_US_LS_ME_hist", 50, 0, 1, 10, -1, 1);

  //------------------------------------------------------------------------------------------------------------------------------

  //L-L
  //signal+background
  TH2F *L0_L0_eta1_vs_eta2_US_ME_hist = new TH2F("L0_L0_eta1_vs_eta2_US_ME_hist", "L0_L0_eta1_vs_eta2_US_ME_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0_phi1_vs_phi2_US_ME_hist = new TH2F("L0_L0_phi1_vs_phi2_US_ME_hist", "L0_L0_phi1_vs_phi2_US_ME_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0_pT1_vs_pT2_US_ME_hist = new TH2F("L0_L0_pT1_vs_pT2_US_ME_hist", "L0_L0_pT1_vs_pT2_US_ME_hist", 20, 0, 5, 20, 0, 5);

  TH1F *L0_L0_delta_phi_US_ME_hist = new TH1F("L0_L0_delta_phi_US_ME_hist", "L0_L0_delta_phi_US_ME_hist", 100, 0, TMath::Pi());
  //TH1F *L0_L0_cos_theta_star_vs_delta_phi_US_ME_hist = new TH1F("L0_L0_cos_theta_star_vs_delta_phi_US_ME_hist", "L0_L0_cos_theta_star_vs_delta_phi_US_ME_hist", 10, -1, 1, 3, 0, TMath::Pi());

  TH1F *L0_L0_delta_phi_US_ME_nFill_hist = new TH1F("L0_L0_delta_phi_US_ME_nFill_hist", "L0_L0_delta_phi_US_ME_nFill_hist", 100, 0, 100);
  TH1F *L0_L0_delta_phi_weighted_US_ME_nFill_hist = new TH1F("L0_L0_delta_phi_weighted_US_ME_nFill_hist", "L0_L0_delta_phi_weighted_US_ME_nFill_hist", 100, 0, 100);


  //background
  //L1(US)-L2(LS)
  TH2F *L0_L0_eta1_vs_eta2_US_LS_ME_1_hist = new TH2F("L0_L0_eta1_vs_eta2_US_LS_ME_1_hist", "L0_L0_eta1_vs_eta2_US_LS_ME_1_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0_phi1_vs_phi2_US_LS_ME_1_hist = new TH2F("L0_L0_phi1_vs_phi2_US_LS_ME_1_hist", "L0_L0_phi1_vs_phi2_US_LS_ME_1_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0_pT1_vs_pT2_US_LS_ME_1_hist = new TH2F("L0_L0_pT1_vs_pT2_US_LS_ME_1_hist", "L0_L0_pT1_vs_pT2_US_LS_ME_1_hist", 20, 0, 5, 20, 0, 5);

  //L1(LS)-L2(US)
  TH2F *L0_L0_eta1_vs_eta2_US_LS_ME_2_hist = new TH2F("L0_L0_eta1_vs_eta2_US_LS_ME_2_hist", "L0_L0_eta1_vs_eta2_US_LS_ME_2_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0_L0_phi1_vs_phi2_US_LS_ME_2_hist = new TH2F("L0_L0_phi1_vs_phi2_US_LS_ME_2_hist", "L0_L0_phi1_vs_phi2_US_LS_ME_2_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0_L0_pT1_vs_pT2_US_LS_ME_2_hist = new TH2F("L0_L0_pT1_vs_pT2_US_LS_ME_2_hist", "L0_L0_pT1_vs_pT2_US_LS_ME_2_hist", 20, 0, 5, 20, 0, 5);

  TH1F *L0_L0_delta_phi_US_LS_ME_hist = new TH1F("L0_L0_delta_phi_US_LS_ME_hist", "L0_L0_delta_phi_US_LS_ME_hist", 100, 0, TMath::Pi());
  //TH1F *L0_L0_cos_theta_star_vs_delta_phi_US_LS_ME_hist = new TH1F("L0_L0_cos_theta_star_vs_delta_phi_US_LS_ME_hist", "L0_L0_cos_theta_star_vs_delta_phi_US_LS_ME_hist", 10, -1, 1, 3, 0, TMath::Pi());


  //----------------------------------------------------------------------------------------

  //Lbar-Lbar
  //signal+background
  TH2F *L0bar_L0bar_eta1_vs_eta2_US_ME_hist = new TH2F("L0bar_L0bar_eta1_vs_eta2_US_ME_hist", "L0bar_L0bar_eta1_vs_eta2_US_ME_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0bar_L0bar_phi1_vs_phi2_US_ME_hist = new TH2F("L0bar_L0bar_phi1_vs_phi2_US_ME_hist", "L0bar_L0bar_phi1_vs_phi2_US_ME_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0bar_L0bar_pT1_vs_pT2_US_ME_hist = new TH2F("L0bar_L0bar_pT1_vs_pT2_US_ME_hist", "L0bar_L0bar_pT1_vs_pT2_US_ME_hist", 20, 0, 5, 20, 0, 5);

  TH1F *L0bar_L0bar_delta_phi_US_ME_hist = new TH1F("L0bar_L0bar_delta_phi_US_ME_hist", "L0bar_L0bar_delta_phi_US_ME_hist", 100, 0, TMath::Pi());
  //TH1F *L0bar_L0bar_cos_theta_star_vs_delta_phi_US_ME_hist = new TH1F("L0bar_L0bar_cos_theta_star_vs_delta_phi_US_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_phi_US_ME_hist", 10, -1, 1, 3, 0, TMath::Pi());

  TH1F *L0bar_L0bar_delta_phi_US_ME_nFill_hist = new TH1F("L0bar_L0bar_delta_phi_US_ME_nFill_hist", "L0bar_L0bar_delta_phi_US_ME_nFill_hist", 100, 0, 100);
  TH1F *L0bar_L0bar_delta_phi_weighted_US_ME_nFill_hist = new TH1F("L0bar_L0bar_delta_phi_weighted_US_ME_nFill_hist", "L0bar_L0bar_delta_phi_weighted_US_ME_nFill_hist", 100, 0, 100);


  //background
  //Lbar1(US)-Lbar2(LS)
  TH2F *L0bar_L0bar_eta1_vs_eta2_US_LS_ME_1_hist = new TH2F("L0bar_L0bar_eta1_vs_eta2_US_LS_ME_1_hist", "L0bar_L0bar_eta1_vs_eta2_US_LS_ME_1_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0bar_L0bar_phi1_vs_phi2_US_LS_ME_1_hist = new TH2F("L0bar_L0bar_phi1_vs_phi2_US_LS_ME_1_hist", "L0bar_L0bar_phi1_vs_phi2_US_LS_ME_1_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0bar_L0bar_pT1_vs_pT2_US_LS_ME_1_hist = new TH2F("L0bar_L0bar_pT1_vs_pT2_US_LS_ME_1_hist", "L0bar_L0bar_pT1_vs_pT2_US_LS_ME_1_hist", 20, 0, 5, 20, 0, 5);

  //Lbar1(LS)-Lbar2(US)
  TH2F *L0bar_L0bar_eta1_vs_eta2_US_LS_ME_2_hist = new TH2F("L0bar_L0bar_eta1_vs_eta2_US_LS_ME_2_hist", "L0bar_L0bar_eta1_vs_eta2_US_LS_ME_2_hist", 20, -1, 1, 20, -1, 1);
  TH2F *L0bar_L0bar_phi1_vs_phi2_US_LS_ME_2_hist = new TH2F("L0bar_L0bar_phi1_vs_phi2_US_LS_ME_2_hist", "L0bar_L0bar_phi1_vs_phi2_US_LS_ME_2_hist", 20, -TMath::Pi(), TMath::Pi(), 20, -TMath::Pi(), TMath::Pi());
  TH2F *L0bar_L0bar_pT1_vs_pT2_US_LS_ME_2_hist = new TH2F("L0bar_L0bar_pT1_vs_pT2_US_LS_ME_2_hist", "L0bar_L0bar_pT1_vs_pT2_US_LS_ME_2_hist", 20, 0, 5, 20, 0, 5);

  TH1F *L0bar_L0bar_delta_phi_US_LS_ME_hist = new TH1F("L0bar_L0bar_delta_phi_US_LS_ME_hist", "L0bar_L0bar_delta_phi_US_LS_ME_hist", 100, 0, TMath::Pi());
  //TH1F *L0bar_L0bar_cos_theta_star_vs_delta_phi_US_LS_ME_hist = new TH1F("L0bar_L0bar_cos_theta_star_vs_delta_phi_US_LS_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_phi_US_LS_ME_hist", 10, -1, 1, 3, 0, TMath::Pi());

  //______________________________________________________________________________________________________________________________________________________________________________________

  //ME Minv histograms for QA

  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_ME[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  //TH2F *L0_inv_mass_vs_L0bar_inv_mass_ME[nPtBins_corr][nPtBins_corr];


  TH2F *L0_inv_mass_vs_L0_inv_mass_US_ME[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_US_LS_ME[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  //TH2F *L0_inv_mass_vs_L0_inv_mass_ME[nPtBins_corr][nPtBins_corr];


  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  //TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_ME[nPtBins_corr][nPtBins_corr];

  //______________________________________________________________________________________________________________________________________________________________________________________

  //LLbarOutFile->cd();

  //data histograms
  TH1F *L0_L0bar_cosThetaProdPlane_US_hist = new TH1F("L0_L0bar_cosThetaProdPlane_US_hist", "L0_L0bar_cosThetaProdPlane_US_hist", 10, -1, 1);
  TH1F *L0_L0bar_cosThetaProdPlane_LS_hist = new TH1F("L0_L0bar_cosThetaProdPlane_LS_hist", "L0_L0bar_cosThetaProdPlane_LS_hist", 10, -1, 1);

  TH1F *L0_L0bar_cosThetaProdPlane_ME_hist = new TH1F("L0_L0bar_cosThetaProdPlane_ME_hist", "L0_L0bar_cosThetaProdPlane_ME_hist", 10, -1, 1);
  TH1F *L0_L0bar_cosThetaProdPlane_ME_LS_hist = new TH1F("L0_L0bar_cosThetaProdPlane_ME_LS_hist", "L0_L0bar_cosThetaProdPlane_ME_LS_hist", 10, -1, 1);

  //TH1F *L0_L0bar_cosThetaProdPlane_US_side_band_hist = new TH1F("L0_L0bar_cosThetaProdPlane_US_side_band_hist", "L0_L0bar_cosThetaProdPlane_US_side_band_hist", 10, -1, 1);

  //---------------

  //Delta phi bins
  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_hist", "L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_hist", 10, -1, 1, 60, 0, TMath::Pi());
  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist", "L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist", 10, -1, 1, 60, 0, TMath::Pi());

  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist", "L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist", 10, -1, 1, 60, 0, TMath::Pi());
  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist", "L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist", 10, -1, 1, 60, 0, TMath::Pi());

  TH3F *L0_L0bar_cos_theta_star_vs_delta_phi_vs_nFill_US_ME_hist = new TH3F("L0_L0bar_cos_theta_star_vs_delta_phi_vs_nFill_US_ME_hist", "L0_L0bar_cos_theta_star_vs_delta_phi_vs_nFill_US_ME_hist", 10, -1, 1, 60, 0, TMath::Pi(), 10, 0, 100);

  //---------------

  //Delta y bins
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_US_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_eta_US_hist", "L0_L0bar_cos_theta_star_vs_delta_eta_US_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_US_LS_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_eta_US_LS_hist", "L0_L0bar_cos_theta_star_vs_delta_eta_US_LS_hist", 10, -1, 1, 2, 0, 2);

  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist", "L0_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist", "L0_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist", 10, -1, 1, 2, 0, 2);

  //Delta y bins scan
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_hist", "L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_hist", 10, -1, 1, 20, 0, 2);
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_hist", "L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_hist", 10, -1, 1, 20, 0, 2);

  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_ME_hist", "L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_ME_hist", 10, -1, 1, 20, 0, 2);
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist", "L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist", 10, -1, 1, 20, 0, 2);

  //---------------

  //Delta y and Delta phi bins
  //Bin 1 is "in cone", bin 2 is out of cone
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist", "L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist", "L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist", "L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist", "L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist", 10, -1, 1, 2, 0, 2);

  //---------------

  //Delta R
  //Bin 1 is "in cone", bin 2 is out of cone
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_US_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_R_US_hist", "L0_L0bar_cos_theta_star_vs_delta_R_US_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_US_LS_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_R_US_LS_hist", "L0_L0bar_cos_theta_star_vs_delta_R_US_LS_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_US_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_R_US_ME_hist", "L0_L0bar_cos_theta_star_vs_delta_R_US_ME_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist", "L0_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist", 10, -1, 1, 2, 0, 2);

  //Delta R for scan
  //Bin 1 is "in cone", bin 2 is out of cone
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_scan_US_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_R_scan_US_hist", "L0_L0bar_cos_theta_star_vs_delta_R_scan_US_hist", 10, -1, 1, 80, 0, 4);
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_hist", "L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_hist", 10, -1, 1, 80, 0, 4);
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist", "L0_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist", 10, -1, 1, 80, 0, 4);
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist", "L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist", 10, -1, 1, 80, 0, 4);

  //---------------

  //Q bins (see ALEPH paper)
  //Bin 1 is "in cone", bin 2 is out of cone
  TH2F *L0_L0bar_cos_theta_star_vs_Q_US_hist = new TH2F("L0_L0bar_cos_theta_star_vs_Q_US_hist", "L0_L0bar_cos_theta_star_vs_Q_US_hist", 10, -1, 1, 20, 0, 10);
  TH2F *L0_L0bar_cos_theta_star_vs_Q_US_LS_hist = new TH2F("L0_L0bar_cos_theta_star_vs_Q_US_LS_hist", "L0_L0bar_cos_theta_star_vs_Q_US_LS_hist", 10, -1, 1, 20, 0, 10);
  TH2F *L0_L0bar_cos_theta_star_vs_Q_US_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_Q_US_ME_hist", "L0_L0bar_cos_theta_star_vs_Q_US_ME_hist", 10, -1, 1, 20, 0, 10);
  TH2F *L0_L0bar_cos_theta_star_vs_Q_US_LS_ME_hist = new TH2F("L0_L0bar_cos_theta_star_vs_Q_US_LS_ME_hist", "L0_L0bar_cos_theta_star_vs_Q_US_LS_ME_hist", 10, -1, 1, 20, 0, 10);

  //---------------------------------------------------------------------------

  TH1F *L0_L0bar_cosThetaProdPlane_pT_US_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_pT_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_pT_ME_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_pT_ME_LS_hist[nPtBins_corr][nPtBins_corr];

  //---------------------------------------------------------------------------

  TH1F *L0_L0bar_cosThetaProdPlane_eta_US_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_LS_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_ME_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist[nEtaBins][nEtaBins];

  //___________________________________________________________________________________________________________________

  TH1F *L0_L0_cosThetaProdPlane_US_hist = new TH1F("L0_L0_cosThetaProdPlane_US_hist", "L0_L0_cosThetaProdPlane_US_hist", 10, -1, 1);
  TH1F *L0_L0_cosThetaProdPlane_LS_hist = new TH1F("L0_L0_cosThetaProdPlane_LS_hist", "L0_L0_cosThetaProdPlane_LS_hist", 10, -1, 1);

  TH1F *L0_L0_cosThetaProdPlane_ME_hist = new TH1F("L0_L0_cosThetaProdPlane_ME_hist", "L0_L0_cosThetaProdPlane_ME_hist", 10, -1, 1);
  TH1F *L0_L0_cosThetaProdPlane_ME_LS_hist = new TH1F("L0_L0_cosThetaProdPlane_ME_LS_hist", "L0_L0_cosThetaProdPlane_ME_LS_hist", 10, -1, 1);

  //TH1F *L0_L0_cosThetaProdPlane_US_side_band_hist = new TH1F("L0_L0_cosThetaProdPlane_US_side_band_hist", "L0_L0_cosThetaProdPlane_US_side_band_hist", 10, -1, 1);

  //Delta phi bins
  TH2F *L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_hist", "L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_hist", 10, -1, 1, 60, 0, TMath::Pi());
  TH2F *L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist", "L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist", 10, -1, 1, 60, 0, TMath::Pi());
  TH2F *L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist", "L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist", 10, -1, 1, 60, 0, TMath::Pi());
  TH2F *L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist", "L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist", 10, -1, 1, 60, 0, TMath::Pi());

  //Delta y bins
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_US_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_eta_US_hist", "L0_L0_cos_theta_star_vs_delta_eta_US_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_US_LS_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_eta_US_LS_hist", "L0_L0_cos_theta_star_vs_delta_eta_US_LS_hist", 10, -1, 1, 2, 0, 2);

  TH2F *L0_L0_cos_theta_star_vs_delta_eta_US_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_eta_US_ME_hist", "L0_L0_cos_theta_star_vs_delta_eta_US_ME_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_US_LS_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_eta_US_LS_ME_hist", "L0_L0_cos_theta_star_vs_delta_eta_US_LS_ME_hist", 10, -1, 1, 2, 0, 2);


  //Delta y bins scan
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_scan_US_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_eta_scan_US_hist", "L0_L0_cos_theta_star_vs_delta_eta_scan_US_hist", 10, -1, 1, 20, 0, 2);
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_scan_US_LS_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_eta_scan_US_LS_hist", "L0_L0_cos_theta_star_vs_delta_eta_scan_US_LS_hist", 10, -1, 1, 20, 0, 2);

  TH2F *L0_L0_cos_theta_star_vs_delta_eta_scan_US_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_eta_scan_US_ME_hist", "L0_L0_cos_theta_star_vs_delta_eta_scan_US_ME_hist", 10, -1, 1, 20, 0, 2);
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist", "L0_L0_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist", 10, -1, 1, 20, 0, 2);


  //Delta y and Delta phi bins
  //Bin 1 is "in cone", bin 2 is out of cone
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_hist", "L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist", "L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist", "L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist", "L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist", 10, -1, 1, 2, 0, 2);

  //Delta R
  //Bin 1 is "in cone", bin 2 is out of cone
  TH2F *L0_L0_cos_theta_star_vs_delta_R_US_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_R_US_hist", "L0_L0_cos_theta_star_vs_delta_R_US_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0_cos_theta_star_vs_delta_R_US_LS_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_R_US_LS_hist", "L0_L0_cos_theta_star_vs_delta_R_US_LS_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0_cos_theta_star_vs_delta_R_US_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_R_US_ME_hist", "L0_L0_cos_theta_star_vs_delta_R_US_ME_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0_L0_cos_theta_star_vs_delta_R_US_LS_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_R_US_LS_ME_hist", "L0_L0_cos_theta_star_vs_delta_R_US_LS_ME_hist", 10, -1, 1, 2, 0, 2);

  //Delta R for scan
  //Bin 1 is "in cone", bin 2 is out of cone
  TH2F *L0_L0_cos_theta_star_vs_delta_R_scan_US_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_R_scan_US_hist", "L0_L0_cos_theta_star_vs_delta_R_scan_US_hist", 10, -1, 1, 80, 0, 4);
  TH2F *L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_hist", "L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_hist", 10, -1, 1, 80, 0, 4);
  TH2F *L0_L0_cos_theta_star_vs_delta_R_scan_US_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_R_scan_US_ME_hist", "L0_L0_cos_theta_star_vs_delta_R_scan_US_ME_hist", 10, -1, 1, 80, 0, 4);
  TH2F *L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist", "L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist", 10, -1, 1, 80, 0, 4);


  //Q bins (see ALEPH paper)
  //Bin 1 is "in cone", bin 2 is out of cone
  TH2F *L0_L0_cos_theta_star_vs_Q_US_hist = new TH2F("L0_L0_cos_theta_star_vs_Q_US_hist", "L0_L0_cos_theta_star_vs_Q_US_hist", 10, -1, 1, 20, 0, 10);
  TH2F *L0_L0_cos_theta_star_vs_Q_US_LS_hist = new TH2F("L0_L0_cos_theta_star_vs_Q_US_LS_hist", "L0_L0_cos_theta_star_vs_Q_US_LS_hist", 10, -1, 1, 20, 0, 10);
  TH2F *L0_L0_cos_theta_star_vs_Q_US_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_Q_US_ME_hist", "L0_L0_cos_theta_star_vs_Q_US_ME_hist", 10, -1, 1, 20, 0, 10);
  TH2F *L0_L0_cos_theta_star_vs_Q_US_LS_ME_hist = new TH2F("L0_L0_cos_theta_star_vs_Q_US_LS_ME_hist", "L0_L0_cos_theta_star_vs_Q_US_LS_ME_hist", 10, -1, 1, 20, 0, 10);


  TH1F *L0_L0_cosThetaProdPlane_pT_US_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_pT_US_side_band_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_pT_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_pT_ME_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_pT_ME_LS_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_eta_US_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0_cosThetaProdPlane_eta_US_side_band_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0_cosThetaProdPlane_eta_LS_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0_cosThetaProdPlane_eta_ME_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0_cosThetaProdPlane_eta_ME_LS_hist[nEtaBins][nEtaBins];

  //_________________________________________________

  TH1F *L0bar_L0bar_cosThetaProdPlane_US_hist = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_hist", "L0bar_L0bar_cosThetaProdPlane_US_hist", 10, -1, 1);
  TH1F *L0bar_L0bar_cosThetaProdPlane_LS_hist = new TH1F("L0bar_L0bar_cosThetaProdPlane_LS_hist", "L0bar_L0bar_cosThetaProdPlane_LS_hist", 10, -1, 1);

  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_hist = new TH1F("L0bar_L0bar_cosThetaProdPlane_ME_hist", "L0bar_L0bar_cosThetaProdPlane_ME_hist", 10, -1, 1);
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_LS_hist = new TH1F("L0bar_L0bar_cosThetaProdPlane_ME_LS_hist", "L0bar_L0bar_cosThetaProdPlane_ME_LS_hist", 10, -1, 1);

  //TH1F *L0bar_L0bar_cosThetaProdPlane_US_side_band_hist = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_side_band_hist", "L0bar_L0bar_cosThetaProdPlane_US_side_band_hist", 10, -1, 1);

  //Delta phi bins
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_hist", "L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_hist", 10, -1, 1, 60, 0, TMath::Pi());
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist", "L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist", 10, -1, 1, 60, 0, TMath::Pi());
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist", 10, -1, 1, 60, 0, TMath::Pi());
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist", 10, -1, 1, 60, 0, TMath::Pi());

  //Delta y bins
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_US_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_eta_US_hist", "L0bar_L0bar_cos_theta_star_vs_delta_eta_US_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_US_LS_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_eta_US_LS_hist", "L0bar_L0bar_cos_theta_star_vs_delta_eta_US_LS_hist", 10, -1, 1, 2, 0, 2);

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist", 10, -1, 1, 2, 0, 2);

  //Delta y bins scan
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_hist", "L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_hist", 10, -1, 1, 20, 0, 2);
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_hist", "L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_hist", 10, -1, 1, 20, 0, 2);

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_ME_hist", 10, -1, 1, 20, 0, 2);
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist", 10, -1, 1, 20, 0, 2);


  //Delta y and Delta phi bins
  //Bin 1 is "in cone", bin 2 is out of cone
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist", "L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist", "L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist", 10, -1, 1, 2, 0, 2);

  //Delta R
  //Bin 1 is "in cone", bin 2 is out of cone
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_US_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_R_US_hist", "L0bar_L0bar_cos_theta_star_vs_delta_R_US_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_hist", "L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_US_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_R_US_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_R_US_ME_hist", 10, -1, 1, 2, 0, 2);
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist", 10, -1, 1, 2, 0, 2);

  //Delta R for scan
  //Bin 1 is "in cone", bin 2 is out of cone
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_hist", "L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_hist", 10, -1, 1, 80, 0, 4);
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_hist", "L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_hist", 10, -1, 1, 80, 0, 4);
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist", 10, -1, 1, 80, 0, 4);
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist", "L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist", 10, -1, 1, 80, 0, 4);


  //Q bins (see ALEPH paper)
  //Bin 1 is "in cone", bin 2 is out of cone
  TH2F *L0bar_L0bar_cos_theta_star_vs_Q_US_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_Q_US_hist", "L0bar_L0bar_cos_theta_star_vs_Q_US_hist", 10, -1, 1, 20, 0, 10);
  TH2F *L0bar_L0bar_cos_theta_star_vs_Q_US_LS_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_Q_US_LS_hist", "L0bar_L0bar_cos_theta_star_vs_Q_US_LS_hist", 10, -1, 1, 20, 0, 10);
  TH2F *L0bar_L0bar_cos_theta_star_vs_Q_US_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_Q_US_ME_hist", "L0bar_L0bar_cos_theta_star_vs_Q_US_ME_hist", 10, -1, 1, 20, 0, 10);
  TH2F *L0bar_L0bar_cos_theta_star_vs_Q_US_LS_ME_hist = new TH2F("L0bar_L0bar_cos_theta_star_vs_Q_US_LS_ME_hist", "L0bar_L0bar_cos_theta_star_vs_Q_US_LS_ME_hist", 10, -1, 1, 20, 0, 10);


  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_US_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_ME_LS_hist[nPtBins_corr][nPtBins_corr];


  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_US_hist[nEtaBins][nEtaBins];
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[nEtaBins][nEtaBins];
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[nEtaBins][nEtaBins];
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[nEtaBins][nEtaBins];
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist[nEtaBins][nEtaBins];
  //___________________________________________________________________________________________

  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_side_band_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_side_band_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_pT_ME_LS_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_side_band_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_side_band_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_pT_ME_LS_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_side_band_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_side_band_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_pT_ME_LS_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      //------------------------------------------------------

      //alternative Minv binning:  60, 1.1, 1.13, 60, 1.1, 1.13
      //full range Minv binning:  180, 1, 1.2, 180, 1, 1.2

      L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_ME_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_ME_pT1_%i_pT2_%i", pTbin1, pTbin2), 60, 1.1, 1.13, 60, 1.1, 1.13);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME_pT1_%i_pT2_%i", pTbin1, pTbin2), 60, 1.1, 1.13, 60, 1.1, 1.13);


      L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_ME_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_ME_pT1_%i_pT2_%i", pTbin1, pTbin2), 60, 1.1, 1.13, 60, 1.1, 1.13);
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_ME_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_LS_ME_pT1_%i_pT2_%i", pTbin1, pTbin2), 60, 1.1, 1.13, 60, 1.1, 1.13);


      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_ME_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_ME_pT1_%i_pT2_%i", pTbin1, pTbin2), 60, 1.1, 1.13, 60, 1.1, 1.13);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME_pT1_%i_pT2_%i", pTbin1, pTbin2), 60, 1.1, 1.13, 60, 1.1, 1.13);
    }
  }



  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_side_band_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_side_band_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_side_band_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_side_band_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_eta_ME_LS_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_side_band_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_side_band_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
    }
  }

  //________________________________________________________________________________________


  cout<<"Analyzing SE signal pairs"<<endl;

  //Unlike-Sign
  //L-Lbar
  for(unsigned int L_Lbar_index = 0; L_Lbar_index < L_Lbar_cos_theta.size(); L_Lbar_index++)
  {
    //float L_peak_mean = fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)]->Parameter(1);
    //float L_peak_sigma = fabs(fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)]->Parameter(2));

    //float Lbar_peak_mean = fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)]->Parameter(3);
    //float Lbar_peak_sigma = fabs(fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)]->Parameter(4));

    float L_peak_mean = L_peak_mean_fit[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)];
    float L_peak_sigma = L_peak_sigma_fit[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)];

    float Lbar_peak_mean = Lbar_peak_mean_fit[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)];
    float Lbar_peak_sigma = Lbar_peak_sigma_fit[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)];

    if( L_Lbar_Minv_L.at(L_Lbar_index) > L_peak_mean-2*L_peak_sigma && L_Lbar_Minv_L.at(L_Lbar_index) < L_peak_mean+2*L_peak_sigma &&
        L_Lbar_Minv_Lbar.at(L_Lbar_index) > Lbar_peak_mean-2*Lbar_peak_sigma && L_Lbar_Minv_Lbar.at(L_Lbar_index) < Lbar_peak_mean+2*Lbar_peak_sigma)
    {


      //-------------------------

      L0_L0bar_cosThetaProdPlane_US_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index));
      L0_L0bar_cosThetaProdPlane_pT_US_hist[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)]->Fill(L_Lbar_cos_theta.at(L_Lbar_index));
      L0_L0bar_cosThetaProdPlane_eta_US_hist[L_Lbar_eta_bin_L.at(L_Lbar_index)][L_Lbar_eta_bin_Lbar.at(L_Lbar_index)]->Fill(L_Lbar_cos_theta.at(L_Lbar_index));

      //fill kinematics histograms for QA and for ME reweight
      L0_L0bar_eta1_vs_eta2_US_hist->Fill(L_Lbar_L_mom.at(L_Lbar_index).Rapidity(), L_Lbar_Lbar_mom.at(L_Lbar_index).Rapidity());
      L0_L0bar_phi1_vs_phi2_US_hist->Fill(L_Lbar_L_mom.at(L_Lbar_index).Phi(), L_Lbar_Lbar_mom.at(L_Lbar_index).Phi());
      L0_L0bar_pT1_vs_pT2_US_hist->Fill(L_Lbar_L_mom.at(L_Lbar_index).Pt(), L_Lbar_Lbar_mom.at(L_Lbar_index).Pt());

      //For J/psi Minv
      TLorentzVector L0_L0bar_4mom_US = L_Lbar_L_mom.at(L_Lbar_index) + L_Lbar_Lbar_mom.at(L_Lbar_index);

      float delta_phi = 0;

      if( fabs(L_Lbar_L_mom.at(L_Lbar_index).Phi() - L_Lbar_Lbar_mom.at(L_Lbar_index).Phi()) <= TMath::Pi() )
      {
        delta_phi = fabs(L_Lbar_L_mom.at(L_Lbar_index).Phi() - L_Lbar_Lbar_mom.at(L_Lbar_index).Phi());
      }
      else
      {
        delta_phi = 2*TMath::Pi() - fabs(L_Lbar_L_mom.at(L_Lbar_index).Phi() - L_Lbar_Lbar_mom.at(L_Lbar_index).Phi()) ;
      }

      L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index), delta_phi);

      L0_L0bar_delta_phi_US_hist->Fill( delta_phi );


      float delta_eta = fabs(L_Lbar_L_mom.at(L_Lbar_index).Rapidity() - L_Lbar_Lbar_mom.at(L_Lbar_index).Rapidity());

      L0_L0bar_delta_eta_US_hist->Fill( delta_eta );
      L0_L0bar_delta_phi_vs_delta_eta_US_hist->Fill( delta_phi, delta_eta );

      L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index), delta_eta);

      if( delta_eta < 0.5 )
      {
        L0_L0bar_cos_theta_star_vs_delta_eta_US_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index), 0.5);

        L0_L0bar_Minv_vs_delta_eta_US_hist->Fill(L0_L0bar_4mom_US.M(), 0.5);
      }
      else
      {
        L0_L0bar_cos_theta_star_vs_delta_eta_US_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index), 1.5);

        L0_L0bar_Minv_vs_delta_eta_US_hist->Fill(L0_L0bar_4mom_US.M(), 1.5);
      }


      if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
      {
        //fill in-cone
        L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index), 0.5);
      }
      else
      {
        //fill out-of-cone
        L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index), 1.5);
      }


      float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

      L0_L0bar_delta_R_US_hist->Fill( delta_R );

      L0_L0bar_cos_theta_star_vs_delta_R_scan_US_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index), delta_R);

      //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
      if( delta_R < 0.93 )
      {
        //fill in-cone
        L0_L0bar_cos_theta_star_vs_delta_R_US_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index), 0.5);
      }
      else
      {
        //fill out-of-cone
        L0_L0bar_cos_theta_star_vs_delta_R_US_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index), 1.5);
      }


      //Q bins

      TLorentzVector LLbar_fourMom_diff = L_Lbar_L_mom.at(L_Lbar_index) - L_Lbar_Lbar_mom.at(L_Lbar_index);

      float LLbar_Q = sqrt(-LLbar_fourMom_diff.Dot(LLbar_fourMom_diff));

      L0_L0bar_Q_US_hist->Fill(LLbar_Q);
      L0_L0bar_delta_R_vs_Q_US_hist->Fill(delta_R, LLbar_Q);

      L0_L0bar_cos_theta_star_vs_Q_US_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index), LLbar_Q);




      L0_L0bar_p1_eta1_vs_p2_eta2_US_hist->Fill(L_Lbar_p_mom.at(L_Lbar_index).Eta(), L_Lbar_pBar_mom.at(L_Lbar_index).Eta());
      L0_L0bar_p1_phi1_vs_p2_phi2_US_hist->Fill(L_Lbar_p_mom.at(L_Lbar_index).Phi(), L_Lbar_pBar_mom.at(L_Lbar_index).Phi());
      L0_L0bar_p1_pT1_vs_p2_pT2_US_hist->Fill(L_Lbar_p_mom.at(L_Lbar_index).Pt(), L_Lbar_pBar_mom.at(L_Lbar_index).Pt());

      L0_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->Fill(L_Lbar_pi_mom.at(L_Lbar_index).Eta(), L_Lbar_piBar_mom.at(L_Lbar_index).Eta());
      L0_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->Fill(L_Lbar_pi_mom.at(L_Lbar_index).Phi(), L_Lbar_piBar_mom.at(L_Lbar_index).Phi());
      L0_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->Fill(L_Lbar_pi_mom.at(L_Lbar_index).Pt(), L_Lbar_piBar_mom.at(L_Lbar_index).Pt());

      //pT vs. cos theta*
      L0_L0bar_pi1_pT_vs_cos_theta_star_US_hist->Fill(L_Lbar_pi_mom.at(L_Lbar_index).Pt(), L_Lbar_cos_theta.at(L_Lbar_index));
      L0_L0bar_pi2_pT_vs_cos_theta_star_US_hist->Fill(L_Lbar_piBar_mom.at(L_Lbar_index).Pt(), L_Lbar_cos_theta.at(L_Lbar_index));

    }
/*
    //continuum only - probably not usefull
    //if( L_Lbar_Minv_L.at(L_Lbar_index) > L_peak_mean+2*L_peak_sigma &&  L_Lbar_Minv_Lbar.at(L_Lbar_index) > Lbar_peak_mean+2*Lbar_peak_sigma)
    //shoud select ridges
    if( (L_Lbar_Minv_L.at(L_Lbar_index) > L_peak_mean+2*L_peak_sigma && L_Lbar_Minv_Lbar.at(L_Lbar_index) > Lbar_peak_mean-2*Lbar_peak_sigma && L_Lbar_Minv_Lbar.at(L_Lbar_index) < Lbar_peak_mean+2*Lbar_peak_sigma) ||
        (L_Lbar_Minv_Lbar.at(L_Lbar_index) > Lbar_peak_mean+2*Lbar_peak_sigma && L_Lbar_Minv_L.at(L_Lbar_index) > L_peak_mean-2*L_peak_sigma && L_Lbar_Minv_L.at(L_Lbar_index) < L_peak_mean+2*L_peak_sigma) )
    {
      L0_L0bar_cosThetaProdPlane_US_side_band_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index));
      L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)]->Fill(L_Lbar_cos_theta.at(L_Lbar_index));
      L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[L_Lbar_eta_bin_L.at(L_Lbar_index)][L_Lbar_eta_bin_Lbar.at(L_Lbar_index)]->Fill(L_Lbar_cos_theta.at(L_Lbar_index));
    }
*/
  }

  //------------------------------------------------------------

  //L-L
  for(unsigned int L_L_index = 0; L_L_index < L_L_cos_theta.size(); L_L_index++)
  {
    //float L1_peak_mean = fit_res_gaus_L0_L0[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)]->Parameter(1);
    //float L1_peak_sigma = fabs(fit_res_gaus_L0_L0[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)]->Parameter(2));

    //float L2_peak_mean = fit_res_gaus_L0_L0[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)]->Parameter(3);
    //float L2_peak_sigma = fabs(fit_res_gaus_L0_L0[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)]->Parameter(4));

    float L1_peak_mean = L1_peak_mean_fit[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)];
    float L1_peak_sigma = L1_peak_sigma_fit[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)];

    float L2_peak_mean = L2_peak_mean_fit[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)];
    float L2_peak_sigma = L2_peak_sigma_fit[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)];


    if( L_L_Minv_L1.at(L_L_index) > L1_peak_mean-2*L1_peak_sigma && L_L_Minv_L1.at(L_L_index) < L1_peak_mean+2*L1_peak_sigma &&
        L_L_Minv_L2.at(L_L_index) > L2_peak_mean-2*L2_peak_sigma && L_L_Minv_L2.at(L_L_index) < L2_peak_mean+2*L2_peak_sigma )
    {
      //-------------------------


      L0_L0_cosThetaProdPlane_US_hist->Fill(L_L_cos_theta.at(L_L_index));
      L0_L0_cosThetaProdPlane_pT_US_hist[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)]->Fill(L_L_cos_theta.at(L_L_index));
      L0_L0_cosThetaProdPlane_eta_US_hist[L_L_eta_bin_L1.at(L_L_index)][L_L_eta_bin_L2.at(L_L_index)]->Fill(L_L_cos_theta.at(L_L_index));

      //fill kinematics histograms for QA and for ME reweight
      L0_L0_eta1_vs_eta2_US_hist->Fill(L_L_L1_mom.at(L_L_index).Rapidity(), L_L_L2_mom.at(L_L_index).Rapidity());
      L0_L0_phi1_vs_phi2_US_hist->Fill(L_L_L1_mom.at(L_L_index).Phi(), L_L_L2_mom.at(L_L_index).Phi());
      L0_L0_pT1_vs_pT2_US_hist->Fill(L_L_L1_mom.at(L_L_index).Pt(), L_L_L2_mom.at(L_L_index).Pt());


      float delta_phi = 0;

      if( fabs(L_L_L1_mom.at(L_L_index).Phi() - L_L_L2_mom.at(L_L_index).Phi()) <= TMath::Pi() )
      {
        delta_phi = fabs(L_L_L1_mom.at(L_L_index).Phi() - L_L_L2_mom.at(L_L_index).Phi());
      }
      else
      {
        delta_phi = 2*TMath::Pi() - fabs(L_L_L1_mom.at(L_L_index).Phi() - L_L_L2_mom.at(L_L_index).Phi());
      }

      L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_hist->Fill(L_L_cos_theta.at(L_L_index), delta_phi);

      L0_L0_delta_phi_US_hist->Fill( delta_phi );


      float delta_eta = fabs(L_L_L1_mom.at(L_L_index).Rapidity() - L_L_L2_mom.at(L_L_index).Rapidity());

      L0_L0_delta_eta_US_hist->Fill( delta_eta );
      L0_L0_delta_phi_vs_delta_eta_US_hist->Fill( delta_phi, delta_eta );

      L0_L0_cos_theta_star_vs_delta_eta_scan_US_hist->Fill(L_L_cos_theta.at(L_L_index), delta_eta);

      if( delta_eta < 0.5 )
      {
        L0_L0_cos_theta_star_vs_delta_eta_US_hist->Fill(L_L_cos_theta.at(L_L_index), 0.5);
      }
      else
      {
        L0_L0_cos_theta_star_vs_delta_eta_US_hist->Fill(L_L_cos_theta.at(L_L_index), 1.5);
      }

      if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
      {
        //fill in-cone
        L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_hist->Fill(L_L_cos_theta.at(L_L_index), 0.5);
      }
      else
      {
        //fill out-of-cone
        L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_hist->Fill(L_L_cos_theta.at(L_L_index), 1.5);
      }


      float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

      L0_L0_delta_R_US_hist->Fill( delta_R );

      L0_L0_cos_theta_star_vs_delta_R_scan_US_hist->Fill(L_L_cos_theta.at(L_L_index), delta_R);

      //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
      if( delta_R < 0.93 )
      {
        //fill in-cone
        L0_L0_cos_theta_star_vs_delta_R_US_hist->Fill(L_L_cos_theta.at(L_L_index), 0.5);
      }
      else
      {
        //fill out-of-cone
        L0_L0_cos_theta_star_vs_delta_R_US_hist->Fill(L_L_cos_theta.at(L_L_index), 1.5);
      }



      //Q bins

      TLorentzVector LL_fourMom_diff = L_L_L1_mom.at(L_L_index) - L_L_L2_mom.at(L_L_index);

      float LL_Q = sqrt(-LL_fourMom_diff.Dot(LL_fourMom_diff));

      L0_L0_Q_US_hist->Fill(LL_Q);
      L0_L0_delta_R_vs_Q_US_hist->Fill(delta_R, LL_Q);

      L0_L0_cos_theta_star_vs_Q_US_hist->Fill(L_L_cos_theta.at(L_L_index), LL_Q);


      //---------------------------------------------------

      L0_L0_p1_eta1_vs_p2_eta2_US_hist->Fill(L_L_p1_mom.at(L_L_index).Eta(), L_L_p2_mom.at(L_L_index).Eta());
      L0_L0_p1_phi1_vs_p2_phi2_US_hist->Fill(L_L_p1_mom.at(L_L_index).Phi(), L_L_p2_mom.at(L_L_index).Phi());
      L0_L0_p1_pT1_vs_p2_pT2_US_hist->Fill(L_L_p1_mom.at(L_L_index).Pt(), L_L_p2_mom.at(L_L_index).Pt());

      L0_L0_pi1_eta1_vs_pi2_eta2_US_hist->Fill(L_L_pi1_mom.at(L_L_index).Eta(), L_L_pi2_mom.at(L_L_index).Eta());
      L0_L0_pi1_phi1_vs_pi2_phi2_US_hist->Fill(L_L_pi1_mom.at(L_L_index).Phi(), L_L_pi2_mom.at(L_L_index).Phi());
      L0_L0_pi1_pT1_vs_pi2_pT2_US_hist->Fill(L_L_pi1_mom.at(L_L_index).Pt(), L_L_pi2_mom.at(L_L_index).Pt());

    }
/*
    //continuum only - probably not usefull
    //if(  L_L_Minv_L1.at(L_L_index) > L1_peak_mean+2*L1_peak_sigma && L_L_Minv_L2.at(L_L_index) > L2_peak_mean+2*L2_peak_sigma )
    //shoud select ridges
    if( (L_L_Minv_L1.at(L_L_index) > L1_peak_mean+2*L1_peak_sigma && L_L_Minv_L2.at(L_L_index) > L2_peak_mean-2*L2_peak_sigma && L_L_Minv_L2.at(L_L_index) < L2_peak_mean+2*L2_peak_sigma) ||
        (L_L_Minv_L2.at(L_L_index) > L2_peak_mean+2*L2_peak_sigma && L_L_Minv_L1.at(L_L_index) > L1_peak_mean-2*L1_peak_sigma && L_L_Minv_L1.at(L_L_index) < L1_peak_mean+2*L1_peak_sigma) )
    {
      L0_L0_cosThetaProdPlane_US_side_band_hist->Fill(L_L_cos_theta.at(L_L_index));
      L0_L0_cosThetaProdPlane_pT_US_side_band_hist[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)]->Fill(L_L_cos_theta.at(L_L_index));
      L0_L0_cosThetaProdPlane_eta_US_side_band_hist[L_L_eta_bin_L1.at(L_L_index)][L_L_eta_bin_L2.at(L_L_index)]->Fill(L_L_cos_theta.at(L_L_index));
    }*/
  }

  //------------------------------------------------------------

  //Lbar-Lbar
  for(unsigned int Lbar_Lbar_index = 0; Lbar_Lbar_index < Lbar_Lbar_cos_theta.size(); Lbar_Lbar_index++)
  {
    //float Lbar1_peak_mean = fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)]->Parameter(1);
    //float Lbar1_peak_sigma = fabs(fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)]->Parameter(2));

    //float Lbar2_peak_mean = fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)]->Parameter(3);
    //float Lbar2_peak_sigma = fabs(fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)]->Parameter(4));

    float Lbar1_peak_mean = Lbar1_peak_mean_fit[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)];
    float Lbar1_peak_sigma = Lbar1_peak_sigma_fit[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)];

    float Lbar2_peak_mean = Lbar2_peak_mean_fit[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)];
    float Lbar2_peak_sigma = Lbar2_peak_sigma_fit[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)];


    if( Lbar_Lbar_Minv_Lbar1.at(Lbar_Lbar_index) > Lbar1_peak_mean-2*Lbar1_peak_sigma && Lbar_Lbar_Minv_Lbar1.at(Lbar_Lbar_index) < Lbar1_peak_mean+2*Lbar1_peak_sigma &&
        Lbar_Lbar_Minv_Lbar2.at(Lbar_Lbar_index) > Lbar2_peak_mean-2*Lbar2_peak_sigma && Lbar_Lbar_Minv_Lbar2.at(Lbar_Lbar_index) < Lbar2_peak_mean+2*Lbar2_peak_sigma)
    {
      //-------------------------

      L0bar_L0bar_cosThetaProdPlane_US_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index));
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)]->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index));
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[Lbar_Lbar_eta_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_eta_bin_Lbar2.at(Lbar_Lbar_index)]->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index));

      //fill kinematics histograms for QA and for ME reweight
      L0bar_L0bar_eta1_vs_eta2_US_hist->Fill(Lbar_Lbar_Lbar1_mom.at(Lbar_Lbar_index).Rapidity(), Lbar_Lbar_Lbar2_mom.at(Lbar_Lbar_index).Rapidity());
      L0bar_L0bar_phi1_vs_phi2_US_hist->Fill(Lbar_Lbar_Lbar1_mom.at(Lbar_Lbar_index).Phi(), Lbar_Lbar_Lbar2_mom.at(Lbar_Lbar_index).Phi());
      L0bar_L0bar_pT1_vs_pT2_US_hist->Fill(Lbar_Lbar_Lbar1_mom.at(Lbar_Lbar_index).Pt(), Lbar_Lbar_Lbar2_mom.at(Lbar_Lbar_index).Pt());


      float delta_phi = 0;

      if( fabs(Lbar_Lbar_Lbar1_mom.at(Lbar_Lbar_index).Phi() - Lbar_Lbar_Lbar2_mom.at(Lbar_Lbar_index).Phi()) <= TMath::Pi() )
      {
        delta_phi = fabs(Lbar_Lbar_Lbar1_mom.at(Lbar_Lbar_index).Phi() - Lbar_Lbar_Lbar2_mom.at(Lbar_Lbar_index).Phi());
      }
      else
      {
        delta_phi = 2*TMath::Pi() - fabs(Lbar_Lbar_Lbar1_mom.at(Lbar_Lbar_index).Phi() - Lbar_Lbar_Lbar2_mom.at(Lbar_Lbar_index).Phi());
      }

      L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index), delta_phi);

      L0bar_L0bar_delta_phi_US_hist->Fill( delta_phi );


      float delta_eta = fabs(Lbar_Lbar_Lbar1_mom.at(Lbar_Lbar_index).Rapidity() - Lbar_Lbar_Lbar2_mom.at(Lbar_Lbar_index).Rapidity());

      L0bar_L0bar_delta_eta_US_hist->Fill( delta_eta );
      L0bar_L0bar_delta_phi_vs_delta_eta_US_hist->Fill( delta_phi, delta_eta );

      L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index), delta_eta);

      if( delta_eta < 0.5 )
      {
        L0bar_L0bar_cos_theta_star_vs_delta_eta_US_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index), 0.5);
      }
      else
      {
        L0bar_L0bar_cos_theta_star_vs_delta_eta_US_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index), 1.5);
      }


      if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
      {
        //fill in-cone
        L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index), 0.5);
      }
      else
      {
        //fill out-of-cone
        L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index), 1.5);
      }


      float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

      L0bar_L0bar_delta_R_US_hist->Fill( delta_R );

      L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index), delta_R);

      //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
      if( delta_R < 0.93 )
      {
        //fill in-cone
        L0bar_L0bar_cos_theta_star_vs_delta_R_US_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index), 0.5);
      }
      else
      {
        //fill out-of-cone
        L0bar_L0bar_cos_theta_star_vs_delta_R_US_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index), 1.5);
      }

      //Q bins

      TLorentzVector LbarLbar_fourMom_diff = Lbar_Lbar_Lbar1_mom.at(Lbar_Lbar_index) - Lbar_Lbar_Lbar2_mom.at(Lbar_Lbar_index);

      float LbarLbar_Q = sqrt(-LbarLbar_fourMom_diff.Dot(LbarLbar_fourMom_diff));

      L0bar_L0bar_Q_US_hist->Fill(LbarLbar_Q);
      L0bar_L0bar_delta_R_vs_Q_US_hist->Fill(delta_R, LbarLbar_Q);

      L0bar_L0bar_cos_theta_star_vs_Q_US_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index), LbarLbar_Q);

      //---------------------------------------------------

      L0bar_L0bar_p1_eta1_vs_p2_eta2_US_hist->Fill(Lbar_Lbar_pBar1_mom.at(Lbar_Lbar_index).Eta(), Lbar_Lbar_pBar2_mom.at(Lbar_Lbar_index).Eta());
      L0bar_L0bar_p1_phi1_vs_p2_phi2_US_hist->Fill(Lbar_Lbar_pBar1_mom.at(Lbar_Lbar_index).Phi(), Lbar_Lbar_pBar2_mom.at(Lbar_Lbar_index).Phi());
      L0bar_L0bar_p1_pT1_vs_p2_pT2_US_hist->Fill(Lbar_Lbar_pBar1_mom.at(Lbar_Lbar_index).Pt(), Lbar_Lbar_pBar2_mom.at(Lbar_Lbar_index).Pt());

      L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_hist->Fill(Lbar_Lbar_piBar1_mom.at(Lbar_Lbar_index).Eta(), Lbar_Lbar_piBar2_mom.at(Lbar_Lbar_index).Eta());
      L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_hist->Fill(Lbar_Lbar_piBar1_mom.at(Lbar_Lbar_index).Phi(), Lbar_Lbar_piBar2_mom.at(Lbar_Lbar_index).Phi());
      L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_hist->Fill(Lbar_Lbar_piBar1_mom.at(Lbar_Lbar_index).Pt(), Lbar_Lbar_piBar2_mom.at(Lbar_Lbar_index).Pt());

    }
/*
    //continuum only - probably not usefull
    //if( Lbar_Lbar_Minv_Lbar1.at(Lbar_Lbar_index) > Lbar1_peak_mean+2*Lbar1_peak_sigma && Lbar_Lbar_Minv_Lbar2.at(Lbar_Lbar_index) > Lbar2_peak_mean+2*Lbar2_peak_sigma)
    //shoud select ridges
    if( (Lbar_Lbar_Minv_Lbar1.at(Lbar_Lbar_index) > Lbar1_peak_mean+2*Lbar1_peak_sigma && Lbar_Lbar_Minv_Lbar2.at(Lbar_Lbar_index) > Lbar2_peak_mean-2*Lbar2_peak_sigma && Lbar_Lbar_Minv_Lbar2.at(Lbar_Lbar_index) < Lbar2_peak_mean+2*Lbar2_peak_sigma) ||
        (Lbar_Lbar_Minv_Lbar2.at(Lbar_Lbar_index) > Lbar2_peak_mean+2*Lbar2_peak_sigma && Lbar_Lbar_Minv_Lbar1.at(Lbar_Lbar_index) > Lbar1_peak_mean-2*Lbar1_peak_sigma && Lbar_Lbar_Minv_Lbar1.at(Lbar_Lbar_index) < Lbar1_peak_mean+2*Lbar1_peak_sigma) )
    {
      L0bar_L0bar_cosThetaProdPlane_US_side_band_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index));
      L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)]->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index));
      L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[Lbar_Lbar_eta_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_eta_bin_Lbar2.at(Lbar_Lbar_index)]->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index));
    }
*/
  }

  //__________________________________________________________________________________________________


  cout<<"Analyzing SE background pairs"<<endl;

  //background
  //L-Lbar (total)
  for(unsigned int L_Lbar_index = 0; L_Lbar_index < L_Lbar_cos_theta_back.size(); L_Lbar_index++)
  {
    //float L_peak_mean = fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)]->Parameter(1);
    //float L_peak_sigma = fabs(fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)]->Parameter(2));

    //float Lbar_peak_mean = fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)]->Parameter(3);
    //float Lbar_peak_sigma = fabs(fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)]->Parameter(4));

    float L_peak_mean = L_peak_mean_fit[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)];
    float L_peak_sigma = L_peak_sigma_fit[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)];

    float Lbar_peak_mean = Lbar_peak_mean_fit[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)];
    float Lbar_peak_sigma = Lbar_peak_sigma_fit[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)];


    //weight for LS - need to scale LS to match realistic sig/bckg ratio known from Minv
    float LS_weight = sideBandScale_L_Lbar[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)];
    //float LS_weight = 1;

    //signal region LS
    if( L_Lbar_Minv_L_back.at(L_Lbar_index) > L_peak_mean-2*L_peak_sigma && L_Lbar_Minv_L_back.at(L_Lbar_index) < L_peak_mean+2*L_peak_sigma &&
        L_Lbar_Minv_Lbar_back.at(L_Lbar_index) > Lbar_peak_mean-2*Lbar_peak_sigma && L_Lbar_Minv_Lbar_back.at(L_Lbar_index) < Lbar_peak_mean+2*Lbar_peak_sigma )
    {

      //-------------------------

      L0_L0bar_cosThetaProdPlane_LS_hist->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), LS_weight);
      L0_L0bar_cosThetaProdPlane_pT_LS_hist[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)]->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), LS_weight);
      L0_L0bar_cosThetaProdPlane_eta_LS_hist[L_Lbar_eta_bin_L_back.at(L_Lbar_index)][L_Lbar_eta_bin_Lbar_back.at(L_Lbar_index)]->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), LS_weight);

      //fill kinematics histograms for QA and for ME reweight
      L0_L0bar_eta1_vs_eta2_US_LS_hist->Fill(L_Lbar_L_mom_back.at(L_Lbar_index).Rapidity(), L_Lbar_Lbar_mom_back.at(L_Lbar_index).Rapidity(), LS_weight);
      L0_L0bar_phi1_vs_phi2_US_LS_hist->Fill(L_Lbar_L_mom_back.at(L_Lbar_index).Phi(), L_Lbar_Lbar_mom_back.at(L_Lbar_index).Phi(), LS_weight);
      L0_L0bar_pT1_vs_pT2_US_LS_hist->Fill(L_Lbar_L_mom_back.at(L_Lbar_index).Pt(), L_Lbar_Lbar_mom_back.at(L_Lbar_index).Pt(), LS_weight);

      //For J/psi Minv
      TLorentzVector L0_L0bar_4mom_US_LS = L_Lbar_L_mom_back.at(L_Lbar_index) + L_Lbar_Lbar_mom_back.at(L_Lbar_index);


      float delta_phi = 0;

      if( fabs(L_Lbar_L_mom_back.at(L_Lbar_index).Phi() - L_Lbar_Lbar_mom_back.at(L_Lbar_index).Phi()) <= TMath::Pi() )
      {
        delta_phi = fabs(L_Lbar_L_mom_back.at(L_Lbar_index).Phi() - L_Lbar_Lbar_mom_back.at(L_Lbar_index).Phi());
      }
      else
      {
        delta_phi = 2*TMath::Pi() - fabs(L_Lbar_L_mom_back.at(L_Lbar_index).Phi() - L_Lbar_Lbar_mom_back.at(L_Lbar_index).Phi()) ;
      }

      L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), delta_phi, LS_weight);

      L0_L0bar_delta_phi_US_LS_hist->Fill( delta_phi, LS_weight );


      float delta_eta = fabs(L_Lbar_L_mom_back.at(L_Lbar_index).Rapidity() - L_Lbar_Lbar_mom_back.at(L_Lbar_index).Rapidity());

      L0_L0bar_delta_eta_US_LS_hist->Fill( delta_eta, LS_weight );
      L0_L0bar_delta_phi_vs_delta_eta_US_LS_hist->Fill( delta_phi, delta_eta, LS_weight );

      L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_hist->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), delta_eta, LS_weight);

      if( delta_eta < 0.5 )
      {
        L0_L0bar_cos_theta_star_vs_delta_eta_US_LS_hist->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), 0.5, LS_weight);

        L0_L0bar_Minv_vs_delta_eta_US_LS_hist->Fill(L0_L0bar_4mom_US_LS.M(), 0.5, LS_weight);
      }
      else
      {
        L0_L0bar_cos_theta_star_vs_delta_eta_US_LS_hist->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), 1.5, LS_weight);

        L0_L0bar_Minv_vs_delta_eta_US_LS_hist->Fill(L0_L0bar_4mom_US_LS.M(), 1.5, LS_weight);
      }


      if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
      {
        //fill in-cone
        L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), 0.5, LS_weight);
      }
      else
      {
        //fill out-of-cone
        L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), 1.5, LS_weight);
      }


      float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

      L0_L0bar_delta_R_US_LS_hist->Fill( delta_R, LS_weight );

      L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_hist->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), delta_R, LS_weight);

      //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
      if( delta_R < 0.93 )
      {
        //fill in-cone
        L0_L0bar_cos_theta_star_vs_delta_R_US_LS_hist->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), 0.5, LS_weight);
      }
      else
      {
        //fill out-of-cone
        L0_L0bar_cos_theta_star_vs_delta_R_US_LS_hist->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), 1.5, LS_weight);
      }


      //Q bins

      TLorentzVector LLbar_fourMom_diff = L_Lbar_L_mom_back.at(L_Lbar_index) - L_Lbar_Lbar_mom_back.at(L_Lbar_index);

      float LLbar_Q = sqrt(-LLbar_fourMom_diff.Dot(LLbar_fourMom_diff));

      L0_L0bar_Q_US_LS_hist->Fill(LLbar_Q, LS_weight);
      L0_L0bar_delta_R_vs_Q_US_LS_hist->Fill(delta_R, LLbar_Q, LS_weight);

      L0_L0bar_cos_theta_star_vs_Q_US_LS_hist->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), LLbar_Q, LS_weight);

      //-----------------------------------------------------------------------------------------------


      //L(US)-Lbar(LS) for ME reweight
      if(L_Lbar_US_LS_flag.at(L_Lbar_index) == 1)
      {
        L0_L0bar_eta1_vs_eta2_US_LS_1_hist->Fill(L_Lbar_L_mom_back.at(L_Lbar_index).Rapidity(), L_Lbar_Lbar_mom_back.at(L_Lbar_index).Rapidity(), LS_weight);
        L0_L0bar_phi1_vs_phi2_US_LS_1_hist->Fill(L_Lbar_L_mom_back.at(L_Lbar_index).Phi(), L_Lbar_Lbar_mom_back.at(L_Lbar_index).Phi(), LS_weight);
        L0_L0bar_pT1_vs_pT2_US_LS_1_hist->Fill(L_Lbar_L_mom_back.at(L_Lbar_index).Pt(), L_Lbar_Lbar_mom_back.at(L_Lbar_index).Pt(), LS_weight);
      }

      //L(LS)-Lbar(US) for ME reweight
      if(L_Lbar_US_LS_flag.at(L_Lbar_index) == 2)
      {
        L0_L0bar_eta1_vs_eta2_US_LS_2_hist->Fill(L_Lbar_L_mom_back.at(L_Lbar_index).Rapidity(), L_Lbar_Lbar_mom_back.at(L_Lbar_index).Rapidity(), LS_weight);
        L0_L0bar_phi1_vs_phi2_US_LS_2_hist->Fill(L_Lbar_L_mom_back.at(L_Lbar_index).Phi(), L_Lbar_Lbar_mom_back.at(L_Lbar_index).Phi(), LS_weight);
        L0_L0bar_pT1_vs_pT2_US_LS_2_hist->Fill(L_Lbar_L_mom_back.at(L_Lbar_index).Pt(), L_Lbar_Lbar_mom_back.at(L_Lbar_index).Pt(), LS_weight);
      }

      L0_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist->Fill(L_Lbar_p_mom_back.at(L_Lbar_index).Eta(), L_Lbar_pBar_mom_back.at(L_Lbar_index).Eta(), LS_weight);
      L0_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist->Fill(L_Lbar_p_mom_back.at(L_Lbar_index).Phi(), L_Lbar_pBar_mom_back.at(L_Lbar_index).Phi(), LS_weight);
      L0_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist->Fill(L_Lbar_p_mom_back.at(L_Lbar_index).Pt(), L_Lbar_pBar_mom_back.at(L_Lbar_index).Pt(), LS_weight);

      L0_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist->Fill(L_Lbar_pi_mom_back.at(L_Lbar_index).Eta(), L_Lbar_piBar_mom_back.at(L_Lbar_index).Eta(), LS_weight);
      L0_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist->Fill(L_Lbar_pi_mom_back.at(L_Lbar_index).Phi(), L_Lbar_piBar_mom_back.at(L_Lbar_index).Phi(), LS_weight);
      L0_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist->Fill(L_Lbar_pi_mom_back.at(L_Lbar_index).Pt(), L_Lbar_piBar_mom_back.at(L_Lbar_index).Pt(), LS_weight);

      //pT vs. cos theta*
      L0_L0bar_pi1_pT_vs_cos_theta_star_US_LS_hist->Fill(L_Lbar_pi_mom_back.at(L_Lbar_index).Pt(), L_Lbar_cos_theta_back.at(L_Lbar_index), LS_weight);
      L0_L0bar_pi2_pT_vs_cos_theta_star_US_LS_hist->Fill(L_Lbar_piBar_mom_back.at(L_Lbar_index).Pt(), L_Lbar_cos_theta_back.at(L_Lbar_index), LS_weight);

    }

  }


  //------------------------------------------------------------------------------------------------------------

  //L-L
  for(unsigned int L_L_index = 0; L_L_index < L_L_cos_theta_back.size(); L_L_index++)
  {
    //float L1_peak_mean = fit_res_gaus_L0_L0[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)]->Parameter(1);
    //float L1_peak_sigma = fabs(fit_res_gaus_L0_L0[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)]->Parameter(2));

    //float L2_peak_mean = fit_res_gaus_L0_L0[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)]->Parameter(3);
    //float L2_peak_sigma = fabs(fit_res_gaus_L0_L0[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)]->Parameter(4));

    float L1_peak_mean = L1_peak_mean_fit[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)];
    float L1_peak_sigma = L1_peak_sigma_fit[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)];

    float L2_peak_mean = L2_peak_mean_fit[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)];
    float L2_peak_sigma = L2_peak_sigma_fit[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)];


    float LS_weight = sideBandScale_L_L[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)];
    //float LS_weight = 1.;


    if( L_L_Minv_L1_back.at(L_L_index) > L1_peak_mean-2*L1_peak_sigma && L_L_Minv_L1_back.at(L_L_index) < L1_peak_mean+2*L1_peak_sigma &&
        L_L_Minv_L2_back.at(L_L_index) > L2_peak_mean-2*L2_peak_sigma && L_L_Minv_L2_back.at(L_L_index) < L2_peak_mean+2*L2_peak_sigma)
    {
      //-------------------------


      L0_L0_cosThetaProdPlane_LS_hist->Fill(L_L_cos_theta_back.at(L_L_index), LS_weight);
      L0_L0_cosThetaProdPlane_pT_LS_hist[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)]->Fill(L_L_cos_theta_back.at(L_L_index), LS_weight);
      L0_L0_cosThetaProdPlane_eta_LS_hist[L_L_eta_bin_L1_back.at(L_L_index)][L_L_eta_bin_L2_back.at(L_L_index)]->Fill(L_L_cos_theta_back.at(L_L_index), LS_weight);

      //fill kinematics histograms for QA and for ME reweight
      L0_L0_eta1_vs_eta2_US_LS_hist->Fill(L_L_L1_mom_back.at(L_L_index).Rapidity(), L_L_L2_mom_back.at(L_L_index).Rapidity(), LS_weight);
      L0_L0_phi1_vs_phi2_US_LS_hist->Fill(L_L_L1_mom_back.at(L_L_index).Phi(), L_L_L2_mom_back.at(L_L_index).Phi(), LS_weight);
      L0_L0_pT1_vs_pT2_US_LS_hist->Fill(L_L_L1_mom_back.at(L_L_index).Pt(), L_L_L2_mom_back.at(L_L_index).Pt(), LS_weight);


      float delta_phi = 0;

      if( fabs(L_L_L1_mom_back.at(L_L_index).Phi() - L_L_L2_mom_back.at(L_L_index).Phi() ) <= TMath::Pi() )
      {
        delta_phi = fabs(L_L_L1_mom_back.at(L_L_index).Phi() - L_L_L2_mom_back.at(L_L_index).Phi() );
      }
      else
      {
        delta_phi = 2*TMath::Pi() - fabs(L_L_L1_mom_back.at(L_L_index).Phi() - L_L_L2_mom_back.at(L_L_index).Phi() ) ;
      }

      L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist->Fill(L_L_cos_theta_back.at(L_L_index), delta_phi, LS_weight);

      L0_L0_delta_phi_US_LS_hist->Fill( delta_phi, LS_weight );


      float delta_eta = fabs(L_L_L1_mom_back.at(L_L_index).Rapidity() - L_L_L2_mom_back.at(L_L_index).Rapidity() );

      L0_L0_delta_eta_US_LS_hist->Fill( delta_eta, LS_weight );
      L0_L0_delta_phi_vs_delta_eta_US_LS_hist->Fill( delta_phi, delta_eta, LS_weight );

      L0_L0_cos_theta_star_vs_delta_eta_scan_US_LS_hist->Fill(L_L_cos_theta_back.at(L_L_index), delta_eta, LS_weight);

      if( delta_eta < 0.5 )
      {
        L0_L0_cos_theta_star_vs_delta_eta_US_LS_hist->Fill(L_L_cos_theta_back.at(L_L_index), 0.5, LS_weight);
      }
      else
      {
        L0_L0_cos_theta_star_vs_delta_eta_US_LS_hist->Fill(L_L_cos_theta_back.at(L_L_index), 1.5, LS_weight);
      }


      if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
      {
        //fill in-cone
        L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist->Fill(L_L_cos_theta_back.at(L_L_index), 0.5, LS_weight);
      }
      else
      {
        //fill out-of-cone
        L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist->Fill(L_L_cos_theta_back.at(L_L_index), 1.5, LS_weight);
      }


      float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

      L0_L0_delta_R_US_LS_hist->Fill( delta_R, LS_weight );

      L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_hist->Fill(L_L_cos_theta_back.at(L_L_index), delta_R, LS_weight);

      //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
      if( delta_R < 0.93 )
      {
        //fill in-cone
        L0_L0_cos_theta_star_vs_delta_R_US_LS_hist->Fill(L_L_cos_theta_back.at(L_L_index), 0.5, LS_weight);
      }
      else
      {
        //fill out-of-cone
        L0_L0_cos_theta_star_vs_delta_R_US_LS_hist->Fill(L_L_cos_theta_back.at(L_L_index), 1.5, LS_weight);
      }


      //Q bins

      TLorentzVector LL_fourMom_diff = L_L_L1_mom_back.at(L_L_index) - L_L_L2_mom_back.at(L_L_index);

      float LL_Q = sqrt(-LL_fourMom_diff.Dot(LL_fourMom_diff));

      L0_L0_Q_US_LS_hist->Fill(LL_Q, LS_weight);
      L0_L0_delta_R_vs_Q_US_LS_hist->Fill(delta_R, LL_Q, LS_weight);

      L0_L0_cos_theta_star_vs_Q_US_LS_hist->Fill(L_L_cos_theta_back.at(L_L_index), LL_Q, LS_weight);


      //L(US)-L(LS) for ME reweight
      if(L_L_US_LS_flag.at(L_L_index) == 1)
      {
        L0_L0_eta1_vs_eta2_US_LS_1_hist->Fill(L_L_L1_mom_back.at(L_L_index).Rapidity(), L_L_L2_mom_back.at(L_L_index).Rapidity(), LS_weight);
        L0_L0_phi1_vs_phi2_US_LS_1_hist->Fill(L_L_L1_mom_back.at(L_L_index).Phi(), L_L_L2_mom_back.at(L_L_index).Phi(), LS_weight);
        L0_L0_pT1_vs_pT2_US_LS_1_hist->Fill(L_L_L1_mom_back.at(L_L_index).Pt(), L_L_L2_mom_back.at(L_L_index).Pt(), LS_weight);
      }

      //L(LS)-L(US) for ME reweight
      if(L_L_US_LS_flag.at(L_L_index) == 2)
      {
        L0_L0_eta1_vs_eta2_US_LS_2_hist->Fill(L_L_L1_mom_back.at(L_L_index).Rapidity(), L_L_L2_mom_back.at(L_L_index).Rapidity(), LS_weight);
        L0_L0_phi1_vs_phi2_US_LS_2_hist->Fill(L_L_L1_mom_back.at(L_L_index).Phi(), L_L_L2_mom_back.at(L_L_index).Phi(), LS_weight);
        L0_L0_pT1_vs_pT2_US_LS_2_hist->Fill(L_L_L1_mom_back.at(L_L_index).Pt(), L_L_L2_mom_back.at(L_L_index).Pt(), LS_weight);
      }


      L0_L0_p1_eta1_vs_p2_eta2_US_LS_hist->Fill(L_L_p1_mom_back.at(L_L_index).Eta(), L_L_p2_mom_back.at(L_L_index).Eta(), LS_weight);
      L0_L0_p1_phi1_vs_p2_phi2_US_LS_hist->Fill(L_L_p1_mom_back.at(L_L_index).Phi(), L_L_p2_mom_back.at(L_L_index).Phi(), LS_weight);
      L0_L0_p1_pT1_vs_p2_pT2_US_LS_hist->Fill(L_L_p1_mom_back.at(L_L_index).Pt(), L_L_p2_mom_back.at(L_L_index).Pt(), LS_weight);

      L0_L0_pi1_eta1_vs_pi2_eta2_US_LS_hist->Fill(L_L_pi1_mom_back.at(L_L_index).Eta(), L_L_pi2_mom_back.at(L_L_index).Eta(), LS_weight);
      L0_L0_pi1_phi1_vs_pi2_phi2_US_LS_hist->Fill(L_L_pi1_mom_back.at(L_L_index).Phi(), L_L_pi2_mom_back.at(L_L_index).Phi(), LS_weight);
      L0_L0_pi1_pT1_vs_pi2_pT2_US_LS_hist->Fill(L_L_pi1_mom_back.at(L_L_index).Pt(), L_L_pi2_mom_back.at(L_L_index).Pt(), LS_weight);

    }

  }

  //------------------------------------------------------------------------------------------------------------

  //Lbar-Lbar
  for(unsigned int Lbar_Lbar_index = 0; Lbar_Lbar_index < Lbar_Lbar_cos_theta_back.size(); Lbar_Lbar_index++)
  {
    //float Lbar1_peak_mean = fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)]->Parameter(1);
    //float Lbar1_peak_sigma = fabs(fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)]->Parameter(2));

    //float Lbar2_peak_mean = fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)]->Parameter(3);
    //float Lbar2_peak_sigma = fabs(fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)]->Parameter(4));

    float Lbar1_peak_mean = Lbar1_peak_mean_fit[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)];
    float Lbar1_peak_sigma = Lbar1_peak_sigma_fit[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)];

    float Lbar2_peak_mean = Lbar2_peak_mean_fit[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)];
    float Lbar2_peak_sigma = Lbar2_peak_sigma_fit[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)];


    float LS_weight = sideBandScale_Lbar_Lbar[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)];
    //float LS_weight = 1.;

    if( Lbar_Lbar_Minv_Lbar1_back.at(Lbar_Lbar_index) > Lbar1_peak_mean-2*Lbar1_peak_sigma && Lbar_Lbar_Minv_Lbar1_back.at(Lbar_Lbar_index) < Lbar1_peak_mean+2*Lbar1_peak_sigma &&
        Lbar_Lbar_Minv_Lbar2_back.at(Lbar_Lbar_index) > Lbar2_peak_mean-2*Lbar2_peak_sigma && Lbar_Lbar_Minv_Lbar2_back.at(Lbar_Lbar_index) < Lbar2_peak_mean+2*Lbar2_peak_sigma)
    {
      //-------------------------

      L0bar_L0bar_cosThetaProdPlane_LS_hist->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), LS_weight);
      L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)]->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), LS_weight);
      L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[Lbar_Lbar_eta_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_eta_bin_Lbar2_back.at(Lbar_Lbar_index)]->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), LS_weight);

      //fill kinematics histograms for QA and for ME reweight
      L0bar_L0bar_eta1_vs_eta2_US_LS_hist->Fill(Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index).Rapidity(), Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index).Rapidity(), LS_weight);
      L0bar_L0bar_phi1_vs_phi2_US_LS_hist->Fill(Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index).Phi(), Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index).Phi(), LS_weight);
      L0bar_L0bar_pT1_vs_pT2_US_LS_hist->Fill(Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index).Pt(), Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index).Pt(), LS_weight);


      float delta_phi = 0;

      if( fabs(Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index).Phi() - Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index).Phi() ) <= TMath::Pi() )
      {
        delta_phi = fabs(Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index).Phi() - Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index).Phi() );
      }
      else
      {
        delta_phi = 2*TMath::Pi() - fabs(Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index).Phi() - Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index).Phi() ) ;
      }

      L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), delta_phi, LS_weight);

      L0bar_L0bar_delta_phi_US_LS_hist->Fill( delta_phi, LS_weight );


      float delta_eta = fabs(Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index).Rapidity() - Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index).Rapidity() );

      L0bar_L0bar_delta_eta_US_LS_hist->Fill( delta_eta, LS_weight );
      L0bar_L0bar_delta_phi_vs_delta_eta_US_LS_hist->Fill( delta_phi, delta_eta, LS_weight );

      L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_hist->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), delta_eta, LS_weight);

      if( delta_eta < 0.5 )
      {
        L0bar_L0bar_cos_theta_star_vs_delta_eta_US_LS_hist->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), 0.5, LS_weight);
      }
      else
      {
        L0bar_L0bar_cos_theta_star_vs_delta_eta_US_LS_hist->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), 1.5, LS_weight);
      }


      if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
      {
        //fill in-cone
        L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), 0.5, LS_weight);
      }
      else
      {
        //fill out-of-cone
        L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), 1.5, LS_weight);
      }


      float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

      L0bar_L0bar_delta_R_US_LS_hist->Fill( delta_R, LS_weight );

      L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_hist->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), delta_R, LS_weight);

      //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
      if( delta_R < 0.93 )
      {
        //fill in-cone
        L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_hist->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), 0.5, LS_weight);
      }
      else
      {
        //fill out-of-cone
        L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_hist->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), 1.5, LS_weight);
      }


      //Q bins

      TLorentzVector LbarLbar_fourMom_diff = Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index) - Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index);

      float LbarLbar_Q = sqrt(-LbarLbar_fourMom_diff.Dot(LbarLbar_fourMom_diff));

      L0bar_L0bar_Q_US_LS_hist->Fill(LbarLbar_Q, LS_weight);
      L0bar_L0bar_delta_R_vs_Q_US_LS_hist->Fill(delta_R, LbarLbar_Q, LS_weight);

      L0bar_L0bar_cos_theta_star_vs_Q_US_LS_hist->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), LbarLbar_Q, LS_weight);


      //-------------------------------------------------------------------------------------------


      //Lbar(US)-Lbar(LS) for ME reweight
      if(Lbar_Lbar_US_LS_flag.at(Lbar_Lbar_index) == 1)
      {
        L0bar_L0bar_eta1_vs_eta2_US_LS_1_hist->Fill(Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index).Rapidity(), Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index).Rapidity(), LS_weight);
        L0bar_L0bar_phi1_vs_phi2_US_LS_1_hist->Fill(Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index).Phi(), Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index).Phi(), LS_weight);
        L0bar_L0bar_pT1_vs_pT2_US_LS_1_hist->Fill(Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index).Pt(), Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index).Pt(), LS_weight);
      }

      //Lbar(LS)-Lbar(US) for ME reweight
      if(Lbar_Lbar_US_LS_flag.at(Lbar_Lbar_index) == 2)
      {
        L0bar_L0bar_eta1_vs_eta2_US_LS_2_hist->Fill(Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index).Rapidity(), Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index).Rapidity(), LS_weight);
        L0bar_L0bar_phi1_vs_phi2_US_LS_2_hist->Fill(Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index).Phi(), Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index).Phi(), LS_weight);
        L0bar_L0bar_pT1_vs_pT2_US_LS_2_hist->Fill(Lbar_Lbar_Lbar1_mom_back.at(Lbar_Lbar_index).Pt(), Lbar_Lbar_Lbar2_mom_back.at(Lbar_Lbar_index).Pt(), LS_weight);
      }

      L0bar_L0bar_p1_eta1_vs_p2_eta2_US_LS_hist->Fill(Lbar_Lbar_pBar1_mom_back.at(Lbar_Lbar_index).Eta(), Lbar_Lbar_pBar2_mom_back.at(Lbar_Lbar_index).Eta(), LS_weight);
      L0bar_L0bar_p1_phi1_vs_p2_phi2_US_LS_hist->Fill(Lbar_Lbar_pBar1_mom_back.at(Lbar_Lbar_index).Phi(), Lbar_Lbar_pBar2_mom_back.at(Lbar_Lbar_index).Phi(), LS_weight);
      L0bar_L0bar_p1_pT1_vs_p2_pT2_US_LS_hist->Fill(Lbar_Lbar_pBar1_mom_back.at(Lbar_Lbar_index).Pt(), Lbar_Lbar_pBar2_mom_back.at(Lbar_Lbar_index).Pt(), LS_weight);

      L0bar_L0bar_pi1_eta1_vs_pi2_eta2_US_LS_hist->Fill(Lbar_Lbar_piBar1_mom_back.at(Lbar_Lbar_index).Eta(), Lbar_Lbar_piBar2_mom_back.at(Lbar_Lbar_index).Eta(), LS_weight);
      L0bar_L0bar_pi1_phi1_vs_pi2_phi2_US_LS_hist->Fill(Lbar_Lbar_piBar1_mom_back.at(Lbar_Lbar_index).Phi(), Lbar_Lbar_piBar2_mom_back.at(Lbar_Lbar_index).Phi(), LS_weight);
      L0bar_L0bar_pi1_pT1_vs_pi2_pT2_US_LS_hist->Fill(Lbar_Lbar_piBar1_mom_back.at(Lbar_Lbar_index).Pt(), Lbar_Lbar_piBar2_mom_back.at(Lbar_Lbar_index).Pt(), LS_weight);

    }

  }

  //________________________________________________________________________________________________________________________________________________________________________________________________


  //analyze mixed-event

  cout<<"Analyzing ME signal pairs"<<endl;

  //signal + background
  //L-Lbar
  //L is from SE, Lbar is from ME

  vector<float> n_fill_LLbar_weight;

  for(unsigned int iLambda_ME = 0; iLambda_ME < L_Lbar_L_vector_ME_SE.size(); iLambda_ME++)
  {
    int fill_LLbar_ME = 0; //flag indicating if given SE pair was already used - want to use just once

    for(unsigned int iLambdaBar_ME = 0; iLambdaBar_ME < Lbar_vector_ME.size(); iLambdaBar_ME++)
    {
      if( (iLambda_ME % 2) !=0 ) break; //choose just even LLbar pairs. Odd will be used later, for L(ME)-Lbar(SE)

      float delta_phi_ME = 0;

      if( fabs( L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).Phi() - Lbar_vector_ME.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).Phi() - Lbar_vector_ME.at(iLambdaBar_ME).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).Phi() - Lbar_vector_ME.at(iLambdaBar_ME).Phi());
      }

      //limit kinematics of ME Lbar based on kinematics of same event Lbar - double check precision
      if( fabs( L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).Rapidity() - Lbar_vector_ME.at(iLambdaBar_ME).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).Pt() - Lbar_vector_ME.at(iLambdaBar_ME).Pt()) > 0.1 ) continue;

      float L_peak_mean = L_peak_mean_fit[L_Lbar_L_pT_bin_vector_ME_SE.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)];
      float L_peak_sigma = L_peak_sigma_fit[L_Lbar_L_pT_bin_vector_ME_SE.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)];

      float Lbar_peak_mean = Lbar_peak_mean_fit[L_Lbar_L_pT_bin_vector_ME_SE.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)];
      float Lbar_peak_sigma = Lbar_peak_sigma_fit[L_Lbar_L_pT_bin_vector_ME_SE.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)];

      if( L_Lbar_L_vector_ME_SE.at(iLambda_ME).M() < L_peak_mean-2*L_peak_sigma || L_Lbar_L_vector_ME_SE.at(iLambda_ME).M() > L_peak_mean+2*L_peak_sigma ) continue;
      if( L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).M() < Lbar_peak_mean-2*Lbar_peak_sigma || L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).M() > Lbar_peak_mean+2*Lbar_peak_sigma ) continue;
      if( Lbar_vector_ME.at(iLambdaBar_ME).M() < Lbar_peak_mean-2*Lbar_peak_sigma || Lbar_vector_ME.at(iLambdaBar_ME).M() > Lbar_peak_mean+2*Lbar_peak_sigma ) continue;

      fill_LLbar_ME++;

    }

    //if( fill_LLbar_ME > 0) n_fill_LLbar_weight.push_back(1./fill_LLbar_ME/2.);
    if( fill_LLbar_ME > 0) n_fill_LLbar_weight.push_back(1./fill_LLbar_ME);
    else n_fill_LLbar_weight.push_back(0);

    L0_L0bar_delta_phi_US_ME_nFill_hist->Fill(fill_LLbar_ME);
    L0_L0bar_delta_phi_weighted_US_ME_nFill_hist->Fill(fill_LLbar_ME, fill_LLbar_ME);

  }

  for(unsigned int iLambda_ME = 0; iLambda_ME < L_Lbar_L_vector_ME_SE.size(); iLambda_ME++)
  {
    if(n_fill_LLbar_weight.at(iLambda_ME) == 0) continue;

    for(unsigned int iLambdaBar_ME = 0; iLambdaBar_ME < Lbar_vector_ME.size(); iLambdaBar_ME++)
    {
      float delta_phi_ME = 0;

      if( fabs( L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).Phi() - Lbar_vector_ME.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).Phi() - Lbar_vector_ME.at(iLambdaBar_ME).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).Phi() - Lbar_vector_ME.at(iLambdaBar_ME).Phi());
      }

      //limit kinematics of ME Lbar based on kinematics of same event Lbar - double check precision
      if( fabs( L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).Rapidity() - Lbar_vector_ME.at(iLambdaBar_ME).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).Pt() - Lbar_vector_ME.at(iLambdaBar_ME).Pt()) > 0.1 ) continue;

      //ME Minv
      //L0_inv_mass_vs_L0bar_inv_mass_US_ME[L_Lbar_L_pT_bin_vector_ME_SE.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)]->Fill(L_Lbar_L_vector_ME_SE.at(iLambda_ME).M(), Lbar_vector_ME.at(iLambdaBar_ME).M(), n_fill_LLbar_weight.at(iLambda_ME));

      float L_peak_mean = L_peak_mean_fit[L_Lbar_L_pT_bin_vector_ME_SE.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)];
      float L_peak_sigma = L_peak_sigma_fit[L_Lbar_L_pT_bin_vector_ME_SE.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)];

      float Lbar_peak_mean = Lbar_peak_mean_fit[L_Lbar_L_pT_bin_vector_ME_SE.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)];
      float Lbar_peak_sigma = Lbar_peak_sigma_fit[L_Lbar_L_pT_bin_vector_ME_SE.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)];

      if( L_Lbar_L_vector_ME_SE.at(iLambda_ME).M() < L_peak_mean-2*L_peak_sigma || L_Lbar_L_vector_ME_SE.at(iLambda_ME).M() > L_peak_mean+2*L_peak_sigma ) continue;
      if( L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).M() < Lbar_peak_mean-2*Lbar_peak_sigma || L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).M() > Lbar_peak_mean+2*Lbar_peak_sigma ) continue;
      if( Lbar_vector_ME.at(iLambdaBar_ME).M() < Lbar_peak_mean-2*Lbar_peak_sigma || Lbar_vector_ME.at(iLambdaBar_ME).M() > Lbar_peak_mean+2*Lbar_peak_sigma ) continue;

      L0_inv_mass_vs_L0bar_inv_mass_US_ME[L_Lbar_L_pT_bin_vector_ME_SE.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)]->Fill(L_Lbar_L_vector_ME_SE.at(iLambda_ME).M(), Lbar_vector_ME.at(iLambdaBar_ME).M(), n_fill_LLbar_weight.at(iLambda_ME));

      //pi kinematics QA
      L0_L0bar_pi1_pT1_vs_pi2_pT2_US_ME_hist->Fill(L_Lbar_pi_vector_ME_SE.at(iLambda_ME).Pt(), piBar_vector_ME.at(iLambdaBar_ME).Pt(), n_fill_LLbar_weight.at(iLambda_ME));

      //-------------------------

      L0_L0bar_eta1_vs_eta2_US_ME_hist->Fill(L_Lbar_L_vector_ME_SE.at(iLambda_ME).Rapidity(), Lbar_vector_ME.at(iLambdaBar_ME).Rapidity(), n_fill_LLbar_weight.at(iLambda_ME));
      L0_L0bar_phi1_vs_phi2_US_ME_hist->Fill(L_Lbar_L_vector_ME_SE.at(iLambda_ME).Phi(), Lbar_vector_ME.at(iLambdaBar_ME).Phi(), n_fill_LLbar_weight.at(iLambda_ME));
      L0_L0bar_pT1_vs_pT2_US_ME_hist->Fill(L_Lbar_L_vector_ME_SE.at(iLambda_ME).Pt(), Lbar_vector_ME.at(iLambdaBar_ME).Pt(), n_fill_LLbar_weight.at(iLambda_ME));

      //------------------------------

      double L_Lbar_pairThetaStar = LpairThetaStar(L_Lbar_L_vector_ME_SE.at(iLambda_ME), L_Lbar_p_vector_ME_SE.at(iLambda_ME), Lbar_vector_ME.at(iLambdaBar_ME), pBar_vector_ME.at(iLambdaBar_ME));

      L0_L0bar_cosThetaProdPlane_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_weight.at(iLambda_ME));
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[L_Lbar_L_pT_bin_vector_ME_SE.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)]->Fill(TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_weight.at(iLambda_ME));
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[L_Lbar_L_eta_bin_vector_ME_SE.at(iLambda_ME)][Lbar_eta_bin_vector_ME.at(iLambdaBar_ME)]->Fill(TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_weight.at(iLambda_ME));

      //pi pT vs. cos(theta*)
      L0_L0bar_pi1_pT_vs_cos_theta_star_US_ME_hist->Fill(L_Lbar_pi_vector_ME_SE.at(iLambda_ME).Pt(), TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_weight.at(iLambda_ME));
      L0_L0bar_pi2_pT_vs_cos_theta_star_US_ME_hist->Fill(piBar_vector_ME.at(iLambdaBar_ME).Pt(), TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_weight.at(iLambda_ME));

      float delta_phi = 0;

      //if( fabs(L_Lbar_L_vector_ME_SE.at(iLambda_ME).Phi() - L_Lbar_Lbar_vector_ME_SE.at(iLambda_ME).Phi()) <= TMath::Pi() )
      if( fabs(L_Lbar_L_vector_ME_SE.at(iLambda_ME).Phi() - Lbar_vector_ME.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
      {
        delta_phi = fabs(L_Lbar_L_vector_ME_SE.at(iLambda_ME).Phi() - Lbar_vector_ME.at(iLambdaBar_ME).Phi());
      }
      else
      {
        delta_phi = 2*TMath::Pi() - fabs(L_Lbar_L_vector_ME_SE.at(iLambda_ME).Phi() - Lbar_vector_ME.at(iLambdaBar_ME).Phi()) ;
      }


      float delta_eta = fabs(L_Lbar_L_vector_ME_SE.at(iLambda_ME).Rapidity() - Lbar_vector_ME.at(iLambdaBar_ME).Rapidity());

      L0_L0bar_delta_phi_US_ME_hist->Fill( delta_phi, n_fill_LLbar_weight.at(iLambda_ME) );

      L0_L0bar_delta_phi_vs_nFill_US_ME_hist->Fill( delta_phi, 1./n_fill_LLbar_weight.at(iLambda_ME) ); // 1./n_fill_LLbar_weight.at(iLambda_ME) shound be N_fill

      L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), delta_eta, n_fill_LLbar_weight.at(iLambda_ME));

      if( delta_eta < 0.5 )
      {
        L0_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 0.5, n_fill_LLbar_weight.at(iLambda_ME));
      }
      else
      {
        L0_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 1.5, n_fill_LLbar_weight.at(iLambda_ME));
      }

      L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), delta_phi, n_fill_LLbar_weight.at(iLambda_ME));

      L0_L0bar_cos_theta_star_vs_delta_phi_vs_nFill_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), delta_phi, 1./n_fill_LLbar_weight.at(iLambda_ME));

      if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
      {
        //fill in-cone
        L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 0.5, n_fill_LLbar_weight.at(iLambda_ME));
      }
      else
      {
        //fill out-of-cone
        L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 1.5, n_fill_LLbar_weight.at(iLambda_ME));
      }


      L0_L0bar_delta_eta_US_ME_hist->Fill( delta_eta, n_fill_LLbar_weight.at(iLambda_ME) );
      L0_L0bar_delta_phi_vs_delta_eta_US_ME_hist->Fill( delta_phi, delta_eta, n_fill_LLbar_weight.at(iLambda_ME) );




      float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

      L0_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), delta_R, n_fill_LLbar_weight.at(iLambda_ME));

      //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
      if( delta_R < 0.93 )
      {
        //fill in-cone
        L0_L0bar_cos_theta_star_vs_delta_R_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 0.5, n_fill_LLbar_weight.at(iLambda_ME));
      }
      else
      {
        //fill out-of-cone
        L0_L0bar_cos_theta_star_vs_delta_R_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 1.5, n_fill_LLbar_weight.at(iLambda_ME));
      }


      //Q bins

      TLorentzVector LLbar_fourMom_diff = L_Lbar_L_vector_ME_SE.at(iLambda_ME) - Lbar_vector_ME.at(iLambdaBar_ME);

      float LLbar_Q = sqrt(-LLbar_fourMom_diff.Dot(LLbar_fourMom_diff));

      //L0_L0bar_Q_US_ME_hist->Fill(LLbar_Q);

      L0_L0bar_cos_theta_star_vs_Q_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), LLbar_Q, n_fill_LLbar_weight.at(iLambda_ME));


    }

  }

  //----------------

  n_fill_LLbar_weight.clear();

  //L is from ME and L bar is from SE
  for(unsigned int iLambdaBar_ME = 0; iLambdaBar_ME < L_Lbar_Lbar_vector_ME_SE.size(); iLambdaBar_ME++)
  {
    int fill_LLbar_ME = 0; //flag indicating if given SE pair was already used - want to use just once

    for(unsigned int iLambda_ME = 0; iLambda_ME < L_vector_ME.size(); iLambda_ME++)
    {
      if( (iLambdaBar_ME % 2) == 0 ) break; //choose just odd LLbar pairs. Even were used earlier, for L(SE)-Lbar(ME)

      float delta_phi_ME = 0;

      if( fabs( L_vector_ME.at(iLambda_ME).Phi() - L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( L_vector_ME.at(iLambda_ME).Phi() - L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( L_vector_ME.at(iLambda_ME).Phi() - L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Phi());
      }

      //limit kinematics of ME L based on kinematics of same event L - double check precision
      if( fabs( L_vector_ME.at(iLambda_ME).Rapidity() - L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( L_vector_ME.at(iLambda_ME).Pt() -  L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Pt()) > 0.1 ) continue;

      float L_peak_mean = L_peak_mean_fit[L_pT_bin_vector_ME.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE.at(iLambdaBar_ME)];
      float L_peak_sigma = L_peak_sigma_fit[L_pT_bin_vector_ME.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE.at(iLambdaBar_ME)];

      float Lbar_peak_mean = Lbar_peak_mean_fit[L_pT_bin_vector_ME.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE.at(iLambdaBar_ME)];
      float Lbar_peak_sigma = Lbar_peak_sigma_fit[L_pT_bin_vector_ME.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE.at(iLambdaBar_ME)];

      if( L_vector_ME.at(iLambda_ME).M() < L_peak_mean-2*L_peak_sigma || L_vector_ME.at(iLambda_ME).M() > L_peak_mean+2*L_peak_sigma ) continue;
      if( L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).M() < L_peak_mean-2*L_peak_sigma || L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).M() > L_peak_mean+2*L_peak_sigma ) continue;
      if( L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).M() < Lbar_peak_mean-2*Lbar_peak_sigma || L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).M() > Lbar_peak_mean+2*Lbar_peak_sigma ) continue;

      fill_LLbar_ME++;

    }

    //if( fill_LLbar_ME > 0) n_fill_LLbar_weight.push_back(1./fill_LLbar_ME/2.);
    if( fill_LLbar_ME > 0) n_fill_LLbar_weight.push_back(1./fill_LLbar_ME);
    else n_fill_LLbar_weight.push_back(0);

    L0_L0bar_delta_phi_US_ME_nFill_hist->Fill(fill_LLbar_ME);
    L0_L0bar_delta_phi_weighted_US_ME_nFill_hist->Fill(fill_LLbar_ME, fill_LLbar_ME);

  }

  //L is from ME and L bar is from SE
  for(unsigned int iLambdaBar_ME = 0; iLambdaBar_ME < L_Lbar_Lbar_vector_ME_SE.size(); iLambdaBar_ME++)
  {
    if(n_fill_LLbar_weight.at(iLambdaBar_ME) == 0) continue;

    for(unsigned int iLambda_ME = 0; iLambda_ME < L_vector_ME.size(); iLambda_ME++)
    {
      float delta_phi_ME = 0;

      if( fabs( L_vector_ME.at(iLambda_ME).Phi() - L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( L_vector_ME.at(iLambda_ME).Phi() - L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( L_vector_ME.at(iLambda_ME).Phi() - L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Phi());
      }

      //limit kinematics of ME L based on kinematics of same event L - double check precision
      if( fabs( L_vector_ME.at(iLambda_ME).Rapidity() - L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( L_vector_ME.at(iLambda_ME).Pt() -  L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).Pt()) > 0.1 ) continue;

      //L0_inv_mass_vs_L0bar_inv_mass_US_ME[L_pT_bin_vector_ME.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE.at(iLambdaBar_ME)]->Fill(L_vector_ME.at(iLambda_ME).M(), L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).M(), n_fill_LLbar_weight.at(iLambdaBar_ME));

      float L_peak_mean = L_peak_mean_fit[L_pT_bin_vector_ME.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE.at(iLambdaBar_ME)];
      float L_peak_sigma = L_peak_sigma_fit[L_pT_bin_vector_ME.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE.at(iLambdaBar_ME)];

      float Lbar_peak_mean = Lbar_peak_mean_fit[L_pT_bin_vector_ME.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE.at(iLambdaBar_ME)];
      float Lbar_peak_sigma = Lbar_peak_sigma_fit[L_pT_bin_vector_ME.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE.at(iLambdaBar_ME)];

      if( L_vector_ME.at(iLambda_ME).M() < L_peak_mean-2*L_peak_sigma || L_vector_ME.at(iLambda_ME).M() > L_peak_mean+2*L_peak_sigma ) continue;
      if( L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).M() < L_peak_mean-2*L_peak_sigma || L_Lbar_L_vector_ME_SE.at(iLambdaBar_ME).M() > L_peak_mean+2*L_peak_sigma ) continue;
      if( L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).M() < Lbar_peak_mean-2*Lbar_peak_sigma || L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).M() > Lbar_peak_mean+2*Lbar_peak_sigma ) continue;

      L0_inv_mass_vs_L0bar_inv_mass_US_ME[L_pT_bin_vector_ME.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE.at(iLambdaBar_ME)]->Fill(L_vector_ME.at(iLambda_ME).M(), L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).M(), n_fill_LLbar_weight.at(iLambdaBar_ME));

      //-------------------------

      //pi kinematics QA
      L0_L0bar_pi1_pT1_vs_pi2_pT2_US_ME_hist->Fill(pi_vector_ME.at(iLambda_ME).Pt(), L_Lbar_piBar_vector_ME_SE.at(iLambdaBar_ME).Pt(), n_fill_LLbar_weight.at(iLambdaBar_ME));


      L0_L0bar_eta1_vs_eta2_US_ME_hist->Fill(L_vector_ME.at(iLambda_ME).Rapidity(), L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).Rapidity(), n_fill_LLbar_weight.at(iLambdaBar_ME));
      L0_L0bar_phi1_vs_phi2_US_ME_hist->Fill(L_vector_ME.at(iLambda_ME).Phi(), L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).Phi(), n_fill_LLbar_weight.at(iLambdaBar_ME));
      L0_L0bar_pT1_vs_pT2_US_ME_hist->Fill(L_vector_ME.at(iLambda_ME).Pt(), L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).Pt(), n_fill_LLbar_weight.at(iLambdaBar_ME));

      //----------

      double L_Lbar_pairThetaStar = LpairThetaStar(L_vector_ME.at(iLambda_ME), p_vector_ME.at(iLambda_ME), L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME), L_Lbar_pBar_vector_ME_SE.at(iLambdaBar_ME));

      L0_L0bar_cosThetaProdPlane_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_weight.at(iLambdaBar_ME));
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[L_pT_bin_vector_ME.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE.at(iLambdaBar_ME)]->Fill(TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_weight.at(iLambdaBar_ME));
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[L_eta_bin_vector_ME.at(iLambda_ME)][L_Lbar_Lbar_eta_bin_vector_ME_SE.at(iLambdaBar_ME)]->Fill(TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_weight.at(iLambdaBar_ME));

      //pi pT vs. cos(theta*)
      L0_L0bar_pi1_pT_vs_cos_theta_star_US_ME_hist->Fill(pi_vector_ME.at(iLambda_ME).Pt(), TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_weight.at(iLambdaBar_ME));
      L0_L0bar_pi2_pT_vs_cos_theta_star_US_ME_hist->Fill(L_Lbar_piBar_vector_ME_SE.at(iLambdaBar_ME).Pt(), TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_weight.at(iLambdaBar_ME));


      float delta_phi = 0;

      if( fabs(L_vector_ME.at(iLambda_ME).Phi() - L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
      {
        delta_phi = fabs(L_vector_ME.at(iLambda_ME).Phi() - L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).Phi());
      }
      else
      {
        delta_phi = 2*TMath::Pi() - fabs(L_vector_ME.at(iLambda_ME).Phi() - L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).Phi()) ;
      }


      float delta_eta = fabs(L_vector_ME.at(iLambda_ME).Rapidity() - L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME).Rapidity());


      L0_L0bar_delta_phi_US_ME_hist->Fill( delta_phi, n_fill_LLbar_weight.at(iLambdaBar_ME) );

      L0_L0bar_delta_eta_US_ME_hist->Fill( delta_eta, n_fill_LLbar_weight.at(iLambdaBar_ME));
      L0_L0bar_delta_phi_vs_delta_eta_US_ME_hist->Fill( delta_phi, delta_eta, n_fill_LLbar_weight.at(iLambdaBar_ME) );

      L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), delta_eta, n_fill_LLbar_weight.at(iLambdaBar_ME));

      if( delta_eta < 0.5 )
      {
        L0_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 0.5, n_fill_LLbar_weight.at(iLambdaBar_ME));
      }
      else
      {
        L0_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 1.5, n_fill_LLbar_weight.at(iLambdaBar_ME));
      }

      L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), delta_phi, n_fill_LLbar_weight.at(iLambdaBar_ME));

      if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
      {
        //fill in-cone
        L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 0.5, n_fill_LLbar_weight.at(iLambdaBar_ME));
      }
      else
      {
        //fill out-of-cone
        L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 1.5, n_fill_LLbar_weight.at(iLambdaBar_ME));
      }




      float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

      L0_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), delta_R, n_fill_LLbar_weight.at(iLambdaBar_ME));

      //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
      if( delta_R < 0.93 )
      {
        //fill in-cone
        L0_L0bar_cos_theta_star_vs_delta_R_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 0.5, n_fill_LLbar_weight.at(iLambdaBar_ME));
      }
      else
      {
        //fill out-of-cone
        L0_L0bar_cos_theta_star_vs_delta_R_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 1.5, n_fill_LLbar_weight.at(iLambdaBar_ME));
      }


      //Q bins

      TLorentzVector LLbar_fourMom_diff = L_vector_ME.at(iLambda_ME) - L_Lbar_Lbar_vector_ME_SE.at(iLambdaBar_ME);

      float LLbar_Q = sqrt(-LLbar_fourMom_diff.Dot(LLbar_fourMom_diff));

      //L0_L0bar_Q_US_ME_hist->Fill(LLbar_Q);

      L0_L0bar_cos_theta_star_vs_Q_US_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), LLbar_Q, n_fill_LLbar_weight.at(iLambdaBar_ME));

    }

  }

  //------------------------------------------------------------

  //L-L

  vector<float> n_fill_LL_weight;

  //L1 is from SE, L2 is from ME
  for(unsigned int iLambda_ME_1 = 0; iLambda_ME_1 < L_L_L1_vector_ME_SE.size(); iLambda_ME_1++)
  {
    int fill_LL_ME = 0;

    for(unsigned int iLambda_ME_2 = 0; iLambda_ME_2 < L_vector_ME.size(); iLambda_ME_2++)
    {
      float delta_phi_ME = 0;

      if( fabs( L_L_L2_vector_ME_SE.at(iLambda_ME_1).Phi() - L_vector_ME.at(iLambda_ME_2).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( L_L_L2_vector_ME_SE.at(iLambda_ME_1).Phi() - L_vector_ME.at(iLambda_ME_2).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( L_L_L2_vector_ME_SE.at(iLambda_ME_1).Phi() - L_vector_ME.at(iLambda_ME_2).Phi());
      }

      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( L_L_L2_vector_ME_SE.at(iLambda_ME_1).Rapidity() - L_vector_ME.at(iLambda_ME_2).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( L_L_L2_vector_ME_SE.at(iLambda_ME_1).Pt() - L_vector_ME.at(iLambda_ME_2).Pt()) > 0.1 ) continue;

      float L1_peak_mean = L1_peak_mean_fit[L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)];
      float L1_peak_sigma = L1_peak_sigma_fit[L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)];

      float L2_peak_mean = L2_peak_mean_fit[L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)];
      float L2_peak_sigma = L2_peak_sigma_fit[L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)];

      if( L_L_L1_vector_ME_SE.at(iLambda_ME_1).M() < L1_peak_mean-2*L1_peak_sigma || L_L_L1_vector_ME_SE.at(iLambda_ME_1).M() > L1_peak_mean+2*L1_peak_sigma ) continue;
      if( L_L_L2_vector_ME_SE.at(iLambda_ME_1).M() < L2_peak_mean-2*L2_peak_sigma || L_L_L2_vector_ME_SE.at(iLambda_ME_1).M() > L2_peak_mean+2*L2_peak_sigma ) continue;
      if( L_vector_ME.at(iLambda_ME_2).M() < L2_peak_mean-2*L2_peak_sigma || L_vector_ME.at(iLambda_ME_2).M() > L2_peak_mean+2*L2_peak_sigma ) continue;

      fill_LL_ME++;

    }

    if( fill_LL_ME > 0) n_fill_LL_weight.push_back(1./fill_LL_ME);
    else n_fill_LL_weight.push_back(0);

    L0_L0_delta_phi_US_ME_nFill_hist->Fill(fill_LL_ME);
    L0_L0_delta_phi_weighted_US_ME_nFill_hist->Fill(fill_LL_ME, fill_LL_ME);

  }


  int nFillLL_ME = 0;

  //L1 is from SE, L2 is from ME
  for(unsigned int iLambda_ME_1 = 0; iLambda_ME_1 < L_L_L1_vector_ME_SE.size(); iLambda_ME_1++)
  {
    if(n_fill_LL_weight.at(iLambda_ME_1) == 0) continue;

    for(unsigned int iLambda_ME_2 = 0; iLambda_ME_2 < L_vector_ME.size(); iLambda_ME_2++)
    {
      float delta_phi_ME = 0;

      if( fabs( L_L_L2_vector_ME_SE.at(iLambda_ME_1).Phi() - L_vector_ME.at(iLambda_ME_2).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( L_L_L2_vector_ME_SE.at(iLambda_ME_1).Phi() - L_vector_ME.at(iLambda_ME_2).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( L_L_L2_vector_ME_SE.at(iLambda_ME_1).Phi() - L_vector_ME.at(iLambda_ME_2).Phi());
      }

      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( L_L_L2_vector_ME_SE.at(iLambda_ME_1).Rapidity() - L_vector_ME.at(iLambda_ME_2).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( L_L_L2_vector_ME_SE.at(iLambda_ME_1).Pt() - L_vector_ME.at(iLambda_ME_2).Pt()) > 0.1 ) continue;

      //if( nFillLL_ME % 2 == 0 )//L1 is from SE, L2 from ME
      //{
        L0_inv_mass_vs_L0_inv_mass_US_ME[L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)]->Fill(L_L_L1_vector_ME_SE.at(iLambda_ME_1).M(), L_vector_ME.at(iLambda_ME_2).M(), n_fill_LL_weight.at(iLambda_ME_1));

        float L1_peak_mean = L1_peak_mean_fit[L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)];
        float L1_peak_sigma = L1_peak_sigma_fit[L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)];

        float L2_peak_mean = L2_peak_mean_fit[L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)];
        float L2_peak_sigma = L2_peak_sigma_fit[L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)];

        if( L_L_L1_vector_ME_SE.at(iLambda_ME_1).M() < L1_peak_mean-2*L1_peak_sigma || L_L_L1_vector_ME_SE.at(iLambda_ME_1).M() > L1_peak_mean+2*L1_peak_sigma ) continue;
        if( L_L_L2_vector_ME_SE.at(iLambda_ME_1).M() < L2_peak_mean-2*L2_peak_sigma || L_L_L2_vector_ME_SE.at(iLambda_ME_1).M() > L2_peak_mean+2*L2_peak_sigma ) continue;
        if( L_vector_ME.at(iLambda_ME_2).M() < L2_peak_mean-2*L2_peak_sigma || L_vector_ME.at(iLambda_ME_2).M() > L2_peak_mean+2*L2_peak_sigma ) continue;

        //-------------------------

        double L_L_pairThetaStar = LpairThetaStar(L_L_L1_vector_ME_SE.at(iLambda_ME_1), L_L_p1_vector_ME_SE.at(iLambda_ME_1), L_vector_ME.at(iLambda_ME_2), p_vector_ME.at(iLambda_ME_2));

        L0_L0_cosThetaProdPlane_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), n_fill_LL_weight.at(iLambda_ME_1));
        L0_L0_cosThetaProdPlane_pT_ME_hist[L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)]->Fill(TMath::Cos(L_L_pairThetaStar), n_fill_LL_weight.at(iLambda_ME_1));
        L0_L0_cosThetaProdPlane_eta_ME_hist[L_L_L1_eta_bin_vector_ME_SE.at(iLambda_ME_1)][L_eta_bin_vector_ME.at(iLambda_ME_2)]->Fill(TMath::Cos(L_L_pairThetaStar), n_fill_LL_weight.at(iLambda_ME_1));


        L0_L0_eta1_vs_eta2_US_ME_hist->Fill(L_L_L1_vector_ME_SE.at(iLambda_ME_1).Rapidity(), L_vector_ME.at(iLambda_ME_2).Rapidity(), n_fill_LL_weight.at(iLambda_ME_1));
        L0_L0_phi1_vs_phi2_US_ME_hist->Fill(L_L_L1_vector_ME_SE.at(iLambda_ME_1).Phi(), L_vector_ME.at(iLambda_ME_2).Phi(), n_fill_LL_weight.at(iLambda_ME_1));
        L0_L0_pT1_vs_pT2_US_ME_hist->Fill(L_L_L1_vector_ME_SE.at(iLambda_ME_1).Pt(), L_vector_ME.at(iLambda_ME_2).Pt(), n_fill_LL_weight.at(iLambda_ME_1));


        float delta_phi = 0;

        //if( fabs(L_L_L_vector_ME_SE.at(iLambda_ME).Phi() - L_L_Lbar_vector_ME_SE.at(iLambda_ME).Phi()) <= TMath::Pi() )
        if( fabs(L_L_L1_vector_ME_SE.at(iLambda_ME_1).Phi() - L_vector_ME.at(iLambda_ME_2).Phi()) <= TMath::Pi() )
        {
          delta_phi = fabs(L_L_L1_vector_ME_SE.at(iLambda_ME_1).Phi() - L_vector_ME.at(iLambda_ME_2).Phi());
        }
        else
        {
          delta_phi = 2*TMath::Pi() - fabs(L_L_L1_vector_ME_SE.at(iLambda_ME_1).Phi() - L_vector_ME.at(iLambda_ME_2).Phi()) ;
        }


        L0_L0_delta_phi_US_ME_hist->Fill( delta_phi, n_fill_LL_weight.at(iLambda_ME_1) );


        float delta_eta = fabs(L_L_L1_vector_ME_SE.at(iLambda_ME_1).Rapidity() - L_vector_ME.at(iLambda_ME_2).Rapidity());
        //float delta_eta_SE = fabs(L_L_L1_vector_ME_SE.at(iLambda_ME_1).Rapidity() - L_L_L2_vector_ME_SE.at(iLambda_ME_1).Rapidity());


        L0_L0_cos_theta_star_vs_delta_eta_scan_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), delta_eta, n_fill_LL_weight.at(iLambda_ME_1));

        if( delta_eta < 0.5 )
        {
          L0_L0_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 0.5, n_fill_LL_weight.at(iLambda_ME_1));
        }
        else
        {
          L0_L0_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 1.5, n_fill_LL_weight.at(iLambda_ME_1));
        }

        L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), delta_phi, n_fill_LL_weight.at(iLambda_ME_1));

        if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
        {
          //fill in-cone
          L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 0.5, n_fill_LL_weight.at(iLambda_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 1.5, n_fill_LL_weight.at(iLambda_ME_1));
        }



        float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

        L0_L0_cos_theta_star_vs_delta_R_scan_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), delta_R, n_fill_LL_weight.at(iLambda_ME_1));

        //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
        if( delta_R < 0.93 )
        {
          //fill in-cone
          L0_L0_cos_theta_star_vs_delta_R_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 0.5, n_fill_LL_weight.at(iLambda_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0_L0_cos_theta_star_vs_delta_R_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 1.5, n_fill_LL_weight.at(iLambda_ME_1));
        }

        //Q bins
        TLorentzVector LL_fourMom_diff = L_L_L1_vector_ME_SE.at(iLambda_ME_1) - L_vector_ME.at(iLambda_ME_2);

        float LL_Q = sqrt(-LL_fourMom_diff.Dot(LL_fourMom_diff));

        //L0_L0bar_Q_US_ME_hist->Fill(LL_Q);

        L0_L0_cos_theta_star_vs_Q_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), LL_Q, n_fill_LL_weight.at(iLambda_ME_1));
/*
        nFillLL_ME++;

      }
     else//L1 is from ME, L2 from SE (here just fill in reverse order)
      {
        L0_inv_mass_vs_L0_inv_mass_US_ME[L_pT_bin_vector_ME.at(iLambda_ME_2)][L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)]->Fill(L_vector_ME.at(iLambda_ME_2).M(), L_L_L1_vector_ME_SE.at(iLambda_ME_1).M(), n_fill_LL_weight.at(iLambda_ME_1));

        float L1_peak_mean = L1_peak_mean_fit[L_pT_bin_vector_ME.at(iLambda_ME_2)][L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)];
        float L1_peak_sigma = L1_peak_sigma_fit[L_pT_bin_vector_ME.at(iLambda_ME_2)][L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)];

        float L2_peak_mean = L2_peak_mean_fit[L_pT_bin_vector_ME.at(iLambda_ME_2)][L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)];
        float L2_peak_sigma = L2_peak_sigma_fit[L_pT_bin_vector_ME.at(iLambda_ME_2)][L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)];

        if( L_vector_ME.at(iLambda_ME_2).M() < L1_peak_mean-2*L1_peak_sigma || L_vector_ME.at(iLambda_ME_2).M() > L1_peak_mean+2*L1_peak_sigma ) continue;
        if( L_L_L2_vector_ME_SE.at(iLambda_ME_1).M() < L1_peak_mean-2*L1_peak_sigma || L_L_L2_vector_ME_SE.at(iLambda_ME_1).M() > L1_peak_mean+2*L1_peak_sigma ) continue;
        if( L_L_L1_vector_ME_SE.at(iLambda_ME_1).M() < L2_peak_mean-2*L2_peak_sigma || L_L_L1_vector_ME_SE.at(iLambda_ME_1).M() > L2_peak_mean+2*L2_peak_sigma ) continue;

        //-------------------------

        double L_L_pairThetaStar = LpairThetaStar(L_vector_ME.at(iLambda_ME_2), p_vector_ME.at(iLambda_ME_2), L_L_L1_vector_ME_SE.at(iLambda_ME_1), L_L_p1_vector_ME_SE.at(iLambda_ME_1));

        L0_L0_cosThetaProdPlane_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), n_fill_LL_weight.at(iLambda_ME_1));
        L0_L0_cosThetaProdPlane_pT_ME_hist[L_pT_bin_vector_ME.at(iLambda_ME_2)][L_L_L1_pT_bin_vector_ME_SE.at(iLambda_ME_1)]->Fill(TMath::Cos(L_L_pairThetaStar), n_fill_LL_weight.at(iLambda_ME_1));
        L0_L0_cosThetaProdPlane_eta_ME_hist[L_eta_bin_vector_ME.at(iLambda_ME_2)][L_L_L1_eta_bin_vector_ME_SE.at(iLambda_ME_1)]->Fill(TMath::Cos(L_L_pairThetaStar), n_fill_LL_weight.at(iLambda_ME_1));


        L0_L0_eta1_vs_eta2_US_ME_hist->Fill(L_vector_ME.at(iLambda_ME_2).Rapidity(), L_L_L1_vector_ME_SE.at(iLambda_ME_1).Rapidity(), n_fill_LL_weight.at(iLambda_ME_1));
        L0_L0_phi1_vs_phi2_US_ME_hist->Fill(L_vector_ME.at(iLambda_ME_2).Phi(), L_L_L1_vector_ME_SE.at(iLambda_ME_1).Phi(), n_fill_LL_weight.at(iLambda_ME_1));
        L0_L0_pT1_vs_pT2_US_ME_hist->Fill(L_vector_ME.at(iLambda_ME_2).Pt(), L_L_L1_vector_ME_SE.at(iLambda_ME_1).Pt(), n_fill_LL_weight.at(iLambda_ME_1));


        float delta_phi = 0;

        //if( fabs(L_L_L_vector_ME_SE.at(iLambda_ME).Phi() - L_L_Lbar_vector_ME_SE.at(iLambda_ME).Phi()) <= TMath::Pi() )
        if( fabs(L_L_L1_vector_ME_SE.at(iLambda_ME_1).Phi() - L_vector_ME.at(iLambda_ME_2).Phi()) <= TMath::Pi() )
        {
          delta_phi = fabs(L_L_L1_vector_ME_SE.at(iLambda_ME_1).Phi() - L_vector_ME.at(iLambda_ME_2).Phi());
        }
        else
        {
          delta_phi = 2*TMath::Pi() - fabs(L_L_L1_vector_ME_SE.at(iLambda_ME_1).Phi() - L_vector_ME.at(iLambda_ME_2).Phi()) ;
        }

        L0_L0_delta_phi_US_ME_hist->Fill( delta_phi, n_fill_LL_weight.at(iLambda_ME_1) );


        float delta_eta = fabs(L_L_L1_vector_ME_SE.at(iLambda_ME_1).Rapidity() - L_vector_ME.at(iLambda_ME_2).Rapidity());


        L0_L0_cos_theta_star_vs_delta_eta_scan_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), delta_eta, n_fill_LL_weight.at(iLambda_ME_1));

        if( delta_eta < 0.5 )
        {
          L0_L0_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 0.5, n_fill_LL_weight.at(iLambda_ME_1));
        }
        else
        {
          L0_L0_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 1.5, n_fill_LL_weight.at(iLambda_ME_1));
        }
        L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), delta_phi, n_fill_LL_weight.at(iLambda_ME_1));

        if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
        {
          //fill in-cone
          L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 0.5, n_fill_LL_weight.at(iLambda_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 1.5, n_fill_LL_weight.at(iLambda_ME_1));
        }




        float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

        L0_L0_cos_theta_star_vs_delta_R_scan_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), delta_R, n_fill_LL_weight.at(iLambda_ME_1));

        //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
        if( delta_R < 0.93 )
        {
          //fill in-cone
          L0_L0_cos_theta_star_vs_delta_R_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 0.5, n_fill_LL_weight.at(iLambda_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0_L0_cos_theta_star_vs_delta_R_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 1.5, n_fill_LL_weight.at(iLambda_ME_1));
        }

        //Q bins
        TLorentzVector LL_fourMom_diff = L_L_L1_vector_ME_SE.at(iLambda_ME_1) - L_vector_ME.at(iLambda_ME_2);

        float LL_Q = sqrt(-LL_fourMom_diff.Dot(LL_fourMom_diff));

        //L0_L0bar_Q_US_ME_hist->Fill(LL_Q);

        L0_L0_cos_theta_star_vs_Q_US_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), LL_Q, n_fill_LL_weight.at(iLambda_ME_1));

        nFillLL_ME++;

      }
*/


    }
  }

  //------------------------------------------------------------

  //Lbar-Lbar

  vector<float> n_fill_LbarLbar_weight;

  //Lbar1 is from SE, Lbar2 is from ME
  for(unsigned int iLambdaBar_ME_1 = 0; iLambdaBar_ME_1 < Lbar_Lbar_Lbar1_vector_ME_SE.size(); iLambdaBar_ME_1++)
  {
    int fill_LbarLbar_ME = 0;

    for(unsigned int iLambdaBar_ME_2 = 0; iLambdaBar_ME_2 < Lbar_vector_ME.size(); iLambdaBar_ME_2++)
    {
      float delta_phi_ME = 0;

      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME.at(iLambdaBar_ME_2).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME.at(iLambdaBar_ME_2).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME.at(iLambdaBar_ME_2).Phi());
      }

      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).Rapidity() - Lbar_vector_ME.at(iLambdaBar_ME_2).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).Pt() - Lbar_vector_ME.at(iLambdaBar_ME_2).Pt()) > 0.1 ) continue;

      float Lbar1_peak_mean = Lbar1_peak_mean_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)];
      float Lbar1_peak_sigma = Lbar1_peak_sigma_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)];

      float Lbar2_peak_mean = Lbar2_peak_mean_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)];
      float Lbar2_peak_sigma = Lbar2_peak_sigma_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)];

      if( Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).M() < Lbar1_peak_mean-2*Lbar1_peak_sigma || Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).M() > Lbar1_peak_mean+2*Lbar1_peak_sigma ) continue;
      if( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).M() < Lbar2_peak_mean-2*Lbar2_peak_sigma || Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).M() > Lbar2_peak_mean+2*Lbar2_peak_sigma ) continue;
      if( Lbar_vector_ME.at(iLambdaBar_ME_2).M() < Lbar2_peak_mean-2*Lbar2_peak_sigma || Lbar_vector_ME.at(iLambdaBar_ME_2).M() > Lbar2_peak_mean+2*Lbar2_peak_sigma ) continue;


      fill_LbarLbar_ME++;
    }

    if( fill_LbarLbar_ME > 0) n_fill_LbarLbar_weight.push_back(1./fill_LbarLbar_ME);
    else n_fill_LbarLbar_weight.push_back(0);

    L0bar_L0bar_delta_phi_US_ME_nFill_hist->Fill(fill_LbarLbar_ME);
    L0bar_L0bar_delta_phi_weighted_US_ME_nFill_hist->Fill(fill_LbarLbar_ME, fill_LbarLbar_ME);

  }


  int nFillLbarLbar_ME = 0;

  //Lbar1 is from SE, Lbar2 is from ME
  for(unsigned int iLambdaBar_ME_1 = 0; iLambdaBar_ME_1 < Lbar_Lbar_Lbar1_vector_ME_SE.size(); iLambdaBar_ME_1++)
  {
    if(n_fill_LbarLbar_weight.at(iLambdaBar_ME_1) == 0) continue;

    for(unsigned int iLambdaBar_ME_2 = 0; iLambdaBar_ME_2 < Lbar_vector_ME.size(); iLambdaBar_ME_2++)
    {
      float delta_phi_ME = 0;

      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME.at(iLambdaBar_ME_2).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME.at(iLambdaBar_ME_2).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME.at(iLambdaBar_ME_2).Phi());
      }

      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).Rapidity() - Lbar_vector_ME.at(iLambdaBar_ME_2).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).Pt() - Lbar_vector_ME.at(iLambdaBar_ME_2).Pt()) > 0.1 ) continue;

      //if(nFillLbarLbar_ME % 2 == 0)
      //{
        L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)]->Fill(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).M(), Lbar_vector_ME.at(iLambdaBar_ME_2).M(), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));

        float Lbar1_peak_mean = Lbar1_peak_mean_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)];
        float Lbar1_peak_sigma = Lbar1_peak_sigma_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)];

        float Lbar2_peak_mean = Lbar2_peak_mean_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)];
        float Lbar2_peak_sigma = Lbar2_peak_sigma_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)];

        if( Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).M() < Lbar1_peak_mean-2*Lbar1_peak_sigma || Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).M() > Lbar1_peak_mean+2*Lbar1_peak_sigma ) continue;
        if( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).M() < Lbar2_peak_mean-2*Lbar2_peak_sigma || Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).M() > Lbar2_peak_mean+2*Lbar2_peak_sigma ) continue;
        if( Lbar_vector_ME.at(iLambdaBar_ME_2).M() < Lbar2_peak_mean-2*Lbar2_peak_sigma || Lbar_vector_ME.at(iLambdaBar_ME_2).M() > Lbar2_peak_mean+2*Lbar2_peak_sigma ) continue;

        //-------------------------

        double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1), Lbar_Lbar_pBar1_vector_ME_SE.at(iLambdaBar_ME_1), Lbar_vector_ME.at(iLambdaBar_ME_2), pBar_vector_ME.at(iLambdaBar_ME_2));

        L0bar_L0bar_cosThetaProdPlane_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[Lbar_Lbar_Lbar1_eta_bin_vector_ME_SE.at(iLambdaBar_ME_1)][Lbar_eta_bin_vector_ME.at(iLambdaBar_ME_2)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));


        L0bar_L0bar_eta1_vs_eta2_US_ME_hist->Fill(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Rapidity(), Lbar_vector_ME.at(iLambdaBar_ME_2).Rapidity(), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_phi1_vs_phi2_US_ME_hist->Fill(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Phi(), Lbar_vector_ME.at(iLambdaBar_ME_2).Phi(), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_pT1_vs_pT2_US_ME_hist->Fill(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Pt(), Lbar_vector_ME.at(iLambdaBar_ME_2).Pt(), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));


        float delta_phi = 0;

        //if( fabs(L_L_L_vector_ME_SE.at(iLambda_ME).Phi() - L_L_Lbar_vector_ME_SE.at(iLambda_ME).Phi()) <= TMath::Pi() )
        if( fabs(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME.at(iLambdaBar_ME_2).Phi()) <= TMath::Pi() )
        {
          delta_phi = fabs(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME.at(iLambdaBar_ME_2).Phi());
        }
        else
        {
          delta_phi = 2*TMath::Pi() - fabs(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME.at(iLambdaBar_ME_2).Phi());
        }


        L0bar_L0bar_delta_phi_US_ME_hist->Fill( delta_phi, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1) );


        float delta_eta = fabs(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Rapidity() - Lbar_vector_ME.at(iLambdaBar_ME_2).Rapidity());


        L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), delta_eta, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));

        if( delta_eta < 0.5 )
        {
          L0bar_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 0.5, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        }
        else
        {
          L0bar_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 1.5, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        }

        L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), delta_phi, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));

        if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
        {
          //fill in-cone
          L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 0.5, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 1.5, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        }




        float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

        L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), delta_R, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));

        //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
        if( delta_R < 0.93 )
        {
          //fill in-cone
          L0bar_L0bar_cos_theta_star_vs_delta_R_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 0.5, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0bar_L0bar_cos_theta_star_vs_delta_R_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 1.5, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        }

        //Q bins
        TLorentzVector LbarLbar_fourMom_diff = Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1) - Lbar_vector_ME.at(iLambdaBar_ME_2);

        float LbarLbar_Q = sqrt(-LbarLbar_fourMom_diff.Dot(LbarLbar_fourMom_diff));

        //L0_L0bar_Q_US_ME_hist->Fill(LbarLbar_Q);

        L0bar_L0bar_cos_theta_star_vs_Q_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), LbarLbar_Q, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
/*
        nFillLbarLbar_ME++;

      }
      else
      {
        L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)]->Fill(Lbar_vector_ME.at(iLambdaBar_ME_2).M(), Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).M(), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));

        float Lbar1_peak_mean = Lbar1_peak_mean_fit[Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)];
        float Lbar1_peak_sigma = Lbar1_peak_sigma_fit[Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)];

        float Lbar2_peak_mean = Lbar2_peak_mean_fit[Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)];
        float Lbar2_peak_sigma = Lbar2_peak_sigma_fit[Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)];

        if( Lbar_vector_ME.at(iLambdaBar_ME_2).M() < Lbar1_peak_mean-2*Lbar1_peak_sigma || Lbar_vector_ME.at(iLambdaBar_ME_2).M() > Lbar1_peak_mean+2*Lbar1_peak_sigma ) continue;
        if( Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).M() < Lbar1_peak_mean-2*Lbar1_peak_sigma || Lbar_Lbar_Lbar2_vector_ME_SE.at(iLambdaBar_ME_1).M() > Lbar1_peak_mean+2*Lbar1_peak_sigma ) continue;
        if( Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).M() < Lbar2_peak_mean-2*Lbar2_peak_sigma || Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).M() > Lbar2_peak_mean+2*Lbar2_peak_sigma ) continue;

        //-------------------------

        double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_vector_ME.at(iLambdaBar_ME_2), pBar_vector_ME.at(iLambdaBar_ME_2), Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1), Lbar_Lbar_pBar1_vector_ME_SE.at(iLambdaBar_ME_1));

        L0bar_L0bar_cosThetaProdPlane_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE.at(iLambdaBar_ME_1)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[Lbar_eta_bin_vector_ME.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_eta_bin_vector_ME_SE.at(iLambdaBar_ME_1)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));


        L0bar_L0bar_eta1_vs_eta2_US_ME_hist->Fill(Lbar_vector_ME.at(iLambdaBar_ME_2).Rapidity(), Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Rapidity(), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_phi1_vs_phi2_US_ME_hist->Fill(Lbar_vector_ME.at(iLambdaBar_ME_2).Phi(), Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Phi(), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_pT1_vs_pT2_US_ME_hist->Fill(Lbar_vector_ME.at(iLambdaBar_ME_2).Pt(), Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Pt(), n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));


        float delta_phi = 0;

        //if( fabs(L_L_L_vector_ME_SE.at(iLambda_ME).Phi() - L_L_Lbar_vector_ME_SE.at(iLambda_ME).Phi()) <= TMath::Pi() )
        if( fabs(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME.at(iLambdaBar_ME_2).Phi()) <= TMath::Pi() )
        {
          delta_phi = fabs(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME.at(iLambdaBar_ME_2).Phi());
        }
        else
        {
          delta_phi = 2*TMath::Pi() - fabs(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME.at(iLambdaBar_ME_2).Phi());
        }

        L0bar_L0bar_delta_phi_US_ME_hist->Fill( delta_phi, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1) );


        float delta_eta = fabs(Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1).Rapidity() - Lbar_vector_ME.at(iLambdaBar_ME_2).Rapidity());


        L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), delta_eta, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));

        if( delta_eta < 0.5 )
        {
          L0bar_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 0.5, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        }
        else
        {
          L0bar_L0bar_cos_theta_star_vs_delta_eta_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 1.5, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        }

        L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), delta_phi, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));

        if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
        {
          //fill in-cone
          L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 0.5, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 1.5, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        }



        float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

        L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), delta_R, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));

        //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
        if( delta_R < 0.93 )
        {
          //fill in-cone
          L0bar_L0bar_cos_theta_star_vs_delta_R_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 0.5, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0bar_L0bar_cos_theta_star_vs_delta_R_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 1.5, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));
        }


        //Q bins
        TLorentzVector LbarLbar_fourMom_diff = Lbar_Lbar_Lbar1_vector_ME_SE.at(iLambdaBar_ME_1) - Lbar_vector_ME.at(iLambdaBar_ME_2);

        float LbarLbar_Q = sqrt(-LbarLbar_fourMom_diff.Dot(LbarLbar_fourMom_diff));

        //L0_L0bar_Q_US_ME_hist->Fill(LbarLbar_Q);

        L0bar_L0bar_cos_theta_star_vs_Q_US_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), LbarLbar_Q, n_fill_LbarLbar_weight.at(iLambdaBar_ME_1));

        //fill_LbarLbar_ME = 1;

        nFillLbarLbar_ME++;
      }
*/


    }
  }

  //________________________________________________________________________________________________________________________________________________________________________________________________________

  //background

  cout<<"Analyzing ME background pairs"<<endl;

  //US paired with LS pairs
  //L-Lbar

  vector<float> n_fill_LLbar_back_weight;

  //L(US) is from SE, Lbar(LS) is from ME
  for(unsigned int iLambda_ME = 0; iLambda_ME < L_Lbar_L_vector_ME_SE_back_US_LS.size(); iLambda_ME++)
  {
    int fill_LLbar_ME = 0;

    for(unsigned int iLambdaBar_ME = 0; iLambdaBar_ME < Lbar_vector_ME_LS.size(); iLambdaBar_ME++)
    {
      //if( (iLambda_ME % 2) !=0 ) break; //choose just even LLbar pairs. Odd will be used later, for L(ME)-Lbar(SE)

      float delta_phi_ME = 0;

      if( fabs( L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME).Phi());
      }

      //limit kinematics of ME Lbar based on kinematics of same event Lbar - double check precision
      if( fabs( L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).Rapidity() - Lbar_vector_ME_LS.at(iLambdaBar_ME).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).Pt() - Lbar_vector_ME_LS.at(iLambdaBar_ME).Pt()) > 0.1 ) continue;

      float L_peak_mean = L_peak_mean_fit[L_Lbar_L_pT_bin_vector_ME_SE_back_US_LS.at(iLambda_ME)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME)];
      float L_peak_sigma = L_peak_sigma_fit[L_Lbar_L_pT_bin_vector_ME_SE_back_US_LS.at(iLambda_ME)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME)];

      float Lbar_peak_mean = Lbar_peak_mean_fit[L_Lbar_L_pT_bin_vector_ME_SE_back_US_LS.at(iLambda_ME)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME)];
      float Lbar_peak_sigma = Lbar_peak_sigma_fit[L_Lbar_L_pT_bin_vector_ME_SE_back_US_LS.at(iLambda_ME)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME)];

      if( L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME).M() < L_peak_mean-2*L_peak_sigma || L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME).M() > L_peak_mean+2*L_peak_sigma ) continue;
      if( L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).M() < Lbar_peak_mean-2*Lbar_peak_sigma || L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).M() > Lbar_peak_mean+2*Lbar_peak_sigma ) continue;
      if( Lbar_vector_ME_LS.at(iLambdaBar_ME).M() < Lbar_peak_mean-2*Lbar_peak_sigma || Lbar_vector_ME_LS.at(iLambdaBar_ME).M() > Lbar_peak_mean+2*Lbar_peak_sigma ) continue;


      fill_LLbar_ME++;

    }

    if( fill_LLbar_ME > 0) n_fill_LLbar_back_weight.push_back(1./fill_LLbar_ME);
    else n_fill_LLbar_back_weight.push_back(0);

  }

  //L(US) is from SE, Lbar(LS) is from ME
  for(unsigned int iLambda_ME = 0; iLambda_ME < L_Lbar_L_vector_ME_SE_back_US_LS.size(); iLambda_ME++)
  {
    if(n_fill_LLbar_back_weight.at(iLambda_ME) == 0) continue;

    for(unsigned int iLambdaBar_ME = 0; iLambdaBar_ME < Lbar_vector_ME_LS.size(); iLambdaBar_ME++)
    {
      float delta_phi_ME = 0;

      if( fabs( L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME).Phi());
      }

      //limit kinematics of ME Lbar based on kinematics of same event Lbar - double check precision
      if( fabs( L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).Rapidity() - Lbar_vector_ME_LS.at(iLambdaBar_ME).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).Pt() - Lbar_vector_ME_LS.at(iLambdaBar_ME).Pt()) > 0.1 ) continue;

      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[L_Lbar_L_pT_bin_vector_ME_SE_back_US_LS.at(iLambda_ME)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME)]->Fill(L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME).M(), Lbar_vector_ME_LS.at(iLambdaBar_ME).M(), n_fill_LLbar_back_weight.at(iLambda_ME));

      float L_peak_mean = L_peak_mean_fit[L_Lbar_L_pT_bin_vector_ME_SE_back_US_LS.at(iLambda_ME)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME)];
      float L_peak_sigma = L_peak_sigma_fit[L_Lbar_L_pT_bin_vector_ME_SE_back_US_LS.at(iLambda_ME)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME)];

      float Lbar_peak_mean = Lbar_peak_mean_fit[L_Lbar_L_pT_bin_vector_ME_SE_back_US_LS.at(iLambda_ME)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME)];
      float Lbar_peak_sigma = Lbar_peak_sigma_fit[L_Lbar_L_pT_bin_vector_ME_SE_back_US_LS.at(iLambda_ME)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME)];

      if( L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME).M() < L_peak_mean-2*L_peak_sigma || L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME).M() > L_peak_mean+2*L_peak_sigma ) continue;
      if( L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).M() < Lbar_peak_mean-2*Lbar_peak_sigma || L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).M() > Lbar_peak_mean+2*Lbar_peak_sigma ) continue;
      if( Lbar_vector_ME_LS.at(iLambdaBar_ME).M() < Lbar_peak_mean-2*Lbar_peak_sigma || Lbar_vector_ME_LS.at(iLambdaBar_ME).M() > Lbar_peak_mean+2*Lbar_peak_sigma ) continue;

      //-------------------------

      double L_Lbar_pairThetaStar = LpairThetaStar(L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME), L_Lbar_p_vector_ME_SE_back_US_LS.at(iLambda_ME), Lbar_vector_ME_LS.at(iLambdaBar_ME), pBar_vector_ME_LS.at(iLambdaBar_ME));

      L0_L0bar_cosThetaProdPlane_ME_LS_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_back_weight.at(iLambda_ME));
      L0_L0bar_cosThetaProdPlane_pT_ME_LS_hist[L_Lbar_L_pT_bin_vector_ME_SE_back_US_LS.at(iLambda_ME)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME)]->Fill(TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_back_weight.at(iLambda_ME));
      L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist[L_Lbar_L_eta_bin_vector_ME_SE_back_US_LS.at(iLambda_ME)][Lbar_eta_bin_vector_ME_LS.at(iLambdaBar_ME)]->Fill(TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_back_weight.at(iLambda_ME));

      //pi kinematics QA
      L0_L0bar_pi1_pT1_vs_pi2_pT2_US_ME_hist->Fill(L_Lbar_pi_vector_ME_SE_back_US_LS.at(iLambda_ME).Pt(), piBar_vector_ME.at(iLambdaBar_ME).Pt(), n_fill_LLbar_back_weight.at(iLambda_ME));

      //pi pT vs. cos(theta*)
      L0_L0bar_pi1_pT_vs_cos_theta_star_US_ME_hist->Fill(L_Lbar_pi_vector_ME_SE_back_US_LS.at(iLambda_ME).Pt(), TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_back_weight.at(iLambda_ME));
      L0_L0bar_pi2_pT_vs_cos_theta_star_US_ME_hist->Fill(piBar_vector_ME_LS.at(iLambdaBar_ME).Pt(), TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_back_weight.at(iLambda_ME));

      L0_L0bar_eta1_vs_eta2_US_LS_ME_1_hist->Fill(L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME).Rapidity(), Lbar_vector_ME_LS.at(iLambdaBar_ME).Rapidity(), n_fill_LLbar_back_weight.at(iLambda_ME));
      L0_L0bar_phi1_vs_phi2_US_LS_ME_1_hist->Fill(L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME).Phi(), Lbar_vector_ME_LS.at(iLambdaBar_ME).Phi(), n_fill_LLbar_back_weight.at(iLambda_ME));
      L0_L0bar_pT1_vs_pT2_US_LS_ME_1_hist->Fill(L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME).Pt(), Lbar_vector_ME_LS.at(iLambdaBar_ME).Pt(), n_fill_LLbar_back_weight.at(iLambda_ME));


      float delta_phi = 0;

      if( fabs(L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME).Phi() - L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).Phi()) <= TMath::Pi() )
      {
        delta_phi = fabs(L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME).Phi() - L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).Phi());
      }
      else
      {
        delta_phi = 2*TMath::Pi() - fabs(L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME).Phi() - L_Lbar_Lbar_vector_ME_SE_back_US_LS.at(iLambda_ME).Phi()) ;
      }


      float delta_eta = fabs(L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME).Rapidity() - Lbar_vector_ME_LS.at(iLambdaBar_ME).Rapidity());

      L0_L0bar_delta_eta_US_LS_ME_hist->Fill( delta_eta, n_fill_LLbar_back_weight.at(iLambda_ME) );
      L0_L0bar_delta_phi_vs_delta_eta_US_LS_ME_hist->Fill( delta_phi, delta_eta, n_fill_LLbar_back_weight.at(iLambda_ME) );


      L0_L0bar_delta_phi_US_LS_ME_hist->Fill( delta_phi, n_fill_LLbar_back_weight.at(iLambda_ME) );

      L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), delta_eta, n_fill_LLbar_back_weight.at(iLambda_ME));

      if( delta_eta < 0.5 )
      {
        L0_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 0.5, n_fill_LLbar_back_weight.at(iLambda_ME));
      }
      else
      {
        L0_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 1.5, n_fill_LLbar_back_weight.at(iLambda_ME));
      }

      L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), delta_phi, n_fill_LLbar_back_weight.at(iLambda_ME));

      if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
      {
        //fill in-cone
        L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 0.5, n_fill_LLbar_back_weight.at(iLambda_ME));
      }
      else
      {
        //fill out-of-cone
        L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 1.5, n_fill_LLbar_back_weight.at(iLambda_ME));
      }



      float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

      L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), delta_R, n_fill_LLbar_back_weight.at(iLambda_ME));

      //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
      if( delta_R < 0.93 )
      {
        //fill in-cone
        L0_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 0.5, n_fill_LLbar_back_weight.at(iLambda_ME));
      }
      else
      {
        //fill out-of-cone
        L0_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 1.5, n_fill_LLbar_back_weight.at(iLambda_ME));
      }

      //Q bins
      TLorentzVector LLbar_fourMom_diff = L_Lbar_L_vector_ME_SE_back_US_LS.at(iLambda_ME) - Lbar_vector_ME_LS.at(iLambdaBar_ME);

      float LLbar_Q = sqrt(-LLbar_fourMom_diff.Dot(LLbar_fourMom_diff));

      //L0_L0bar_Q_US_ME_hist->Fill(LLbar_Q);

      L0_L0bar_cos_theta_star_vs_Q_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), LLbar_Q, n_fill_LLbar_back_weight.at(iLambda_ME));

      //fill_LLbar_ME = 1;
    }
  }

  //----------

  n_fill_LLbar_back_weight.clear();

  //L(LS) is from ME and Lbar(US) bar is from SE
  for(unsigned int iLambdaBar_ME = 0; iLambdaBar_ME < L_Lbar_Lbar_vector_ME_SE_back_LS_US.size(); iLambdaBar_ME++)
  {
    int fill_LLbar_ME = 0;

    for(unsigned int iLambda_ME = 0; iLambda_ME < L_vector_ME_LS.size(); iLambda_ME++)
    {
      //if( (iLambdaBar_ME % 2) == 0 ) break; //choose just odd LLbar pairs. Even were used earlier, for L(SE)-Lbar(ME)

      float delta_phi_ME = 0;

      if( fabs( L_vector_ME_LS.at(iLambda_ME).Phi() - L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( L_vector_ME_LS.at(iLambda_ME).Phi() - L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( L_vector_ME_LS.at(iLambda_ME).Phi() - L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Phi());
      }

      //limit kinematics of ME L based on kinematics of same event L - double check precision
      if( fabs( L_vector_ME_LS.at(iLambda_ME).Rapidity() - L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( L_vector_ME_LS.at(iLambda_ME).Pt() -  L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Pt()) > 0.1 ) continue;

      float L_peak_mean = L_peak_mean_fit[L_pT_bin_vector_ME_LS.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE_back_LS_US.at(iLambdaBar_ME)];
      float L_peak_sigma = L_peak_sigma_fit[L_pT_bin_vector_ME_LS.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE_back_LS_US.at(iLambdaBar_ME)];

      float Lbar_peak_mean = Lbar_peak_mean_fit[L_pT_bin_vector_ME_LS.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE_back_LS_US.at(iLambdaBar_ME)];
      float Lbar_peak_sigma = Lbar_peak_sigma_fit[L_pT_bin_vector_ME_LS.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE_back_LS_US.at(iLambdaBar_ME)];

      if( L_vector_ME_LS.at(iLambda_ME).M() < L_peak_mean-2*L_peak_sigma || L_vector_ME_LS.at(iLambda_ME).M() > L_peak_mean+2*L_peak_sigma ) continue;
      if( L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).M() < L_peak_mean-2*L_peak_sigma || L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).M() > L_peak_mean+2*L_peak_sigma ) continue;
      if( L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).M() < Lbar_peak_mean-2*Lbar_peak_sigma || L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).M() > Lbar_peak_mean+2*Lbar_peak_sigma ) continue;


      fill_LLbar_ME++;

    }

    if( fill_LLbar_ME > 0) n_fill_LLbar_back_weight.push_back(1./fill_LLbar_ME);
    else n_fill_LLbar_back_weight.push_back(0);

  }


  //L(LS) is from ME and Lbar(US) bar is from SE
  for(unsigned int iLambdaBar_ME = 0; iLambdaBar_ME < L_Lbar_Lbar_vector_ME_SE_back_LS_US.size(); iLambdaBar_ME++)
  {
    if(n_fill_LLbar_back_weight.at(iLambdaBar_ME) == 0) continue;

    for(unsigned int iLambda_ME = 0; iLambda_ME < L_vector_ME_LS.size(); iLambda_ME++)
    {
      float delta_phi_ME = 0;

      if( fabs( L_vector_ME_LS.at(iLambda_ME).Phi() - L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( L_vector_ME_LS.at(iLambda_ME).Phi() - L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( L_vector_ME_LS.at(iLambda_ME).Phi() - L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Phi());
      }

      //limit kinematics of ME L based on kinematics of same event L - double check precision
      if( fabs( L_vector_ME_LS.at(iLambda_ME).Rapidity() - L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( L_vector_ME_LS.at(iLambda_ME).Pt() -  L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Pt()) > 0.1 ) continue;

      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[L_pT_bin_vector_ME_LS.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE_back_LS_US.at(iLambdaBar_ME)]->Fill(L_vector_ME_LS.at(iLambda_ME).M(), L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).M(), n_fill_LLbar_back_weight.at(iLambdaBar_ME));

      float L_peak_mean = L_peak_mean_fit[L_pT_bin_vector_ME_LS.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE_back_LS_US.at(iLambdaBar_ME)];
      float L_peak_sigma = L_peak_sigma_fit[L_pT_bin_vector_ME_LS.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE_back_LS_US.at(iLambdaBar_ME)];

      float Lbar_peak_mean = Lbar_peak_mean_fit[L_pT_bin_vector_ME_LS.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE_back_LS_US.at(iLambdaBar_ME)];
      float Lbar_peak_sigma = Lbar_peak_sigma_fit[L_pT_bin_vector_ME_LS.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE_back_LS_US.at(iLambdaBar_ME)];

      if( L_vector_ME_LS.at(iLambda_ME).M() < L_peak_mean-2*L_peak_sigma || L_vector_ME_LS.at(iLambda_ME).M() > L_peak_mean+2*L_peak_sigma ) continue;
      if( L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).M() < L_peak_mean-2*L_peak_sigma || L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).M() > L_peak_mean+2*L_peak_sigma ) continue;
      if( L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).M() < Lbar_peak_mean-2*Lbar_peak_sigma || L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).M() > Lbar_peak_mean+2*Lbar_peak_sigma ) continue;

      //-------------------------

      double L_Lbar_pairThetaStar = LpairThetaStar(L_vector_ME_LS.at(iLambda_ME), p_vector_ME_LS.at(iLambda_ME), L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME), L_Lbar_pBar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME));


      L0_L0bar_cosThetaProdPlane_ME_LS_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_back_weight.at(iLambdaBar_ME));
      L0_L0bar_cosThetaProdPlane_pT_ME_LS_hist[L_pT_bin_vector_ME_LS.at(iLambda_ME)][L_Lbar_Lbar_pT_bin_vector_ME_SE_back_LS_US.at(iLambdaBar_ME)]->Fill(TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_back_weight.at(iLambdaBar_ME));
      L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist[L_eta_bin_vector_ME_LS.at(iLambda_ME)][L_Lbar_Lbar_eta_bin_vector_ME_SE_back_LS_US.at(iLambdaBar_ME)]->Fill(TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_back_weight.at(iLambdaBar_ME));

      //pi kinematics QA
      L0_L0bar_pi1_pT1_vs_pi2_pT2_US_ME_hist->Fill(pi_vector_ME_LS.at(iLambda_ME).Pt(), L_Lbar_piBar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Pt(), n_fill_LLbar_back_weight.at(iLambdaBar_ME));

      //pi pT vs. cos(theta*)
      L0_L0bar_pi1_pT_vs_cos_theta_star_US_ME_hist->Fill(pi_vector_ME_LS.at(iLambda_ME).Pt(), TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_back_weight.at(iLambdaBar_ME));
      L0_L0bar_pi2_pT_vs_cos_theta_star_US_ME_hist->Fill(L_Lbar_piBar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Pt(), TMath::Cos(L_Lbar_pairThetaStar), n_fill_LLbar_back_weight.at(iLambdaBar_ME));

      L0_L0bar_eta1_vs_eta2_US_LS_ME_2_hist->Fill(L_vector_ME_LS.at(iLambda_ME).Rapidity(), L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Rapidity(), n_fill_LLbar_back_weight.at(iLambdaBar_ME));
      L0_L0bar_phi1_vs_phi2_US_LS_ME_2_hist->Fill(L_vector_ME_LS.at(iLambda_ME).Phi(), L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Phi(), n_fill_LLbar_back_weight.at(iLambdaBar_ME));
      L0_L0bar_pT1_vs_pT2_US_LS_ME_2_hist->Fill(L_vector_ME_LS.at(iLambda_ME).Pt(), L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Pt(), n_fill_LLbar_back_weight.at(iLambdaBar_ME));


      float delta_phi = 0;

      if( fabs(L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Phi() - L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Phi()) <= TMath::Pi() )
      {
        delta_phi = fabs(L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Phi() - L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Phi());
      }
      else
      {
        delta_phi = 2*TMath::Pi() - fabs(L_Lbar_L_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Phi() - L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Phi()) ;
      }


      L0_L0bar_delta_phi_US_LS_ME_hist->Fill( delta_phi, n_fill_LLbar_back_weight.at(iLambdaBar_ME) );


      float delta_eta = fabs(L_vector_ME_LS.at(iLambda_ME).Rapidity() - L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME).Rapidity());

      L0_L0bar_delta_eta_US_LS_ME_hist->Fill( delta_eta, n_fill_LLbar_back_weight.at(iLambdaBar_ME) );
      L0_L0bar_delta_phi_vs_delta_eta_US_LS_ME_hist->Fill( delta_phi, delta_eta, n_fill_LLbar_back_weight.at(iLambdaBar_ME) );


      L0_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), delta_eta, n_fill_LLbar_back_weight.at(iLambdaBar_ME));

      if( delta_eta < 0.5 )
      {
        L0_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 0.5, n_fill_LLbar_back_weight.at(iLambdaBar_ME));
      }
      else
      {
        L0_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 1.5, n_fill_LLbar_back_weight.at(iLambdaBar_ME));
      }

      L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), delta_phi, n_fill_LLbar_back_weight.at(iLambdaBar_ME));

      if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
      {
        //fill in-cone
        L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 0.5, n_fill_LLbar_back_weight.at(iLambdaBar_ME));
      }
      else
      {
        //fill out-of-cone
        L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 1.5, n_fill_LLbar_back_weight.at(iLambdaBar_ME));
      }




      float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

      L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), delta_R, n_fill_LLbar_back_weight.at(iLambdaBar_ME));

      //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
      if( delta_R < 0.93 )
      {
        //fill in-cone
        L0_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 0.5, n_fill_LLbar_back_weight.at(iLambdaBar_ME));
      }
      else
      {
        //fill out-of-cone
        L0_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), 1.5, n_fill_LLbar_back_weight.at(iLambdaBar_ME));
      }

      //Q bins
      TLorentzVector LLbar_fourMom_diff = L_vector_ME_LS.at(iLambda_ME) - L_Lbar_Lbar_vector_ME_SE_back_LS_US.at(iLambdaBar_ME);

      float LLbar_Q = sqrt(-LLbar_fourMom_diff.Dot(LLbar_fourMom_diff));

      //L0_L0bar_Q_US_ME_hist->Fill(LLbar_Q);

      L0_L0bar_cos_theta_star_vs_Q_US_LS_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar), LLbar_Q, n_fill_LLbar_back_weight.at(iLambdaBar_ME));

    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  //L-L

  vector<float> n_fill_LL_back_weight;

  //in SE pair, firtst L is US, second is LS
  for(unsigned int iLambda_ME_1 = 0; iLambda_ME_1 < L_L_L1_vector_ME_SE_back.size(); iLambda_ME_1++)
  {
    int fill_LL_ME = 0;

    for(unsigned int iLambda_ME_2 = 0; iLambda_ME_2 < L_vector_ME_LS.size(); iLambda_ME_2++)
    {
      float delta_phi_ME = 0;

      if( fabs( L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).Phi() - L_vector_ME_LS.at(iLambda_ME_2).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).Phi() - L_vector_ME_LS.at(iLambda_ME_2).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).Phi() - L_vector_ME_LS.at(iLambda_ME_2).Phi());
      }

      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).Rapidity() - L_vector_ME_LS.at(iLambda_ME_2).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).Pt() - L_vector_ME_LS.at(iLambda_ME_2).Pt()) > 0.1 ) continue;

      float L1_peak_mean = L1_peak_mean_fit[L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)][L_pT_bin_vector_ME_LS.at(iLambda_ME_2)];
      float L1_peak_sigma = L1_peak_sigma_fit[L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)][L_pT_bin_vector_ME_LS.at(iLambda_ME_2)];

      float L2_peak_mean = L2_peak_mean_fit[L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)][L_pT_bin_vector_ME_LS.at(iLambda_ME_2)];
      float L2_peak_sigma = L2_peak_sigma_fit[L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)][L_pT_bin_vector_ME_LS.at(iLambda_ME_2)];

      if( L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).M() < L1_peak_mean-2*L1_peak_sigma || L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).M() > L1_peak_mean+2*L1_peak_sigma ) continue;
      if( L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).M() < L2_peak_mean-2*L2_peak_sigma || L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).M() > L2_peak_mean+2*L2_peak_sigma ) continue;
      if( L_vector_ME_LS.at(iLambda_ME_2).M() < L2_peak_mean-2*L2_peak_sigma || L_vector_ME_LS.at(iLambda_ME_2).M() > L2_peak_mean+2*L2_peak_sigma ) continue;

      fill_LL_ME++;
    }

    if( fill_LL_ME > 0) n_fill_LL_back_weight.push_back(1./fill_LL_ME);
    else n_fill_LL_back_weight.push_back(0);

  }

  int nFillLL_ME_bckg = 0;

  //in SE pair, firtst L is US, second is LS
  for(unsigned int iLambda_ME_1 = 0; iLambda_ME_1 < L_L_L1_vector_ME_SE_back.size(); iLambda_ME_1++)
  {
    if(n_fill_LL_back_weight.at(iLambda_ME_1) == 0) continue;

    for(unsigned int iLambda_ME_2 = 0; iLambda_ME_2 < L_vector_ME_LS.size(); iLambda_ME_2++)
    {
      float delta_phi_ME = 0;

      if( fabs( L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).Phi() - L_vector_ME_LS.at(iLambda_ME_2).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).Phi() - L_vector_ME_LS.at(iLambda_ME_2).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).Phi() - L_vector_ME_LS.at(iLambda_ME_2).Phi());
      }

      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).Rapidity() - L_vector_ME_LS.at(iLambda_ME_2).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).Pt() - L_vector_ME_LS.at(iLambda_ME_2).Pt()) > 0.1 ) continue;

      //if( nFillLL_ME_bckg % 2 == 0 ) //L1(US) is from SE, L2(LS) is from ME
      //{
        L0_inv_mass_vs_L0_inv_mass_US_LS_ME[L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)][L_pT_bin_vector_ME_LS.at(iLambda_ME_2)]->Fill(L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).M(), L_vector_ME_LS.at(iLambda_ME_2).M(), n_fill_LL_back_weight.at(iLambda_ME_1));

        float L1_peak_mean = L1_peak_mean_fit[L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)][L_pT_bin_vector_ME_LS.at(iLambda_ME_2)];
        float L1_peak_sigma = L1_peak_sigma_fit[L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)][L_pT_bin_vector_ME_LS.at(iLambda_ME_2)];

        float L2_peak_mean = L2_peak_mean_fit[L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)][L_pT_bin_vector_ME_LS.at(iLambda_ME_2)];
        float L2_peak_sigma = L2_peak_sigma_fit[L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)][L_pT_bin_vector_ME_LS.at(iLambda_ME_2)];

        if( L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).M() < L1_peak_mean-2*L1_peak_sigma || L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).M() > L1_peak_mean+2*L1_peak_sigma ) continue;
        if( L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).M() < L2_peak_mean-2*L2_peak_sigma || L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).M() > L2_peak_mean+2*L2_peak_sigma ) continue;
        if( L_vector_ME_LS.at(iLambda_ME_2).M() < L2_peak_mean-2*L2_peak_sigma || L_vector_ME_LS.at(iLambda_ME_2).M() > L2_peak_mean+2*L2_peak_sigma ) continue;

        //-------------------------

        double L_L_pairThetaStar = LpairThetaStar(L_L_L1_vector_ME_SE_back.at(iLambda_ME_1), L_L_p1_vector_ME_SE_back.at(iLambda_ME_1), L_vector_ME_LS.at(iLambda_ME_2), p_vector_ME_LS.at(iLambda_ME_2));

        L0_L0_cosThetaProdPlane_ME_LS_hist->Fill(TMath::Cos(L_L_pairThetaStar), n_fill_LL_back_weight.at(iLambda_ME_1));
        L0_L0_cosThetaProdPlane_pT_ME_LS_hist[L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)][L_pT_bin_vector_ME_LS.at(iLambda_ME_2)]->Fill(TMath::Cos(L_L_pairThetaStar), n_fill_LL_back_weight.at(iLambda_ME_1));
        L0_L0_cosThetaProdPlane_eta_ME_LS_hist[L_L_L1_eta_bin_vector_ME_SE_back.at(iLambda_ME_1)][L_eta_bin_vector_ME_LS.at(iLambda_ME_2)]->Fill(TMath::Cos(L_L_pairThetaStar), n_fill_LL_back_weight.at(iLambda_ME_1));


        L0_L0_eta1_vs_eta2_US_LS_ME_1_hist->Fill(L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Rapidity(), L_vector_ME_LS.at(iLambda_ME_2).Rapidity(), n_fill_LL_back_weight.at(iLambda_ME_1));
        L0_L0_phi1_vs_phi2_US_LS_ME_1_hist->Fill(L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Phi(), L_vector_ME_LS.at(iLambda_ME_2).Phi(), n_fill_LL_back_weight.at(iLambda_ME_1));
        L0_L0_pT1_vs_pT2_US_LS_ME_1_hist->Fill(L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Pt(), L_vector_ME_LS.at(iLambda_ME_2).Pt(), n_fill_LL_back_weight.at(iLambda_ME_1));


        float delta_phi = 0;

        if( fabs(L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Phi() - L_vector_ME_LS.at(iLambda_ME_2).Phi()) <= TMath::Pi() )
        {
          delta_phi = fabs(L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Phi() - L_vector_ME_LS.at(iLambda_ME_2).Phi());
        }
        else
        {
          delta_phi = 2*TMath::Pi() - fabs(L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Phi() - L_vector_ME_LS.at(iLambda_ME_2).Phi());
        }

        L0_L0_delta_phi_US_LS_ME_hist->Fill( delta_phi, n_fill_LL_back_weight.at(iLambda_ME_1) );


        float delta_eta = fabs(L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Rapidity() - L_vector_ME_LS.at(iLambda_ME_2).Rapidity());


        L0_L0_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), delta_eta, n_fill_LL_back_weight.at(iLambda_ME_1));

        if( delta_eta < 0.5 )
        {
          L0_L0_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 0.5, n_fill_LL_back_weight.at(iLambda_ME_1));
        }
        else
        {
          L0_L0_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 1.5, n_fill_LL_back_weight.at(iLambda_ME_1));
        }

        L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), delta_phi, n_fill_LL_back_weight.at(iLambda_ME_1));

        if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
        {
          //fill in-cone
          L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 0.5, n_fill_LL_back_weight.at(iLambda_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 1.5, n_fill_LL_back_weight.at(iLambda_ME_1));
        }


        float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

        L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), delta_R, n_fill_LL_back_weight.at(iLambda_ME_1));

        //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
        if( delta_R < 0.93 )
        {
          //fill in-cone
          L0_L0_cos_theta_star_vs_delta_R_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 0.5, n_fill_LL_back_weight.at(iLambda_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0_L0_cos_theta_star_vs_delta_R_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 1.5, n_fill_LL_back_weight.at(iLambda_ME_1));
        }

        //Q bins
        TLorentzVector LL_fourMom_diff = L_L_L1_vector_ME_SE_back.at(iLambda_ME_1) - L_vector_ME_LS.at(iLambda_ME_2);

        float LL_Q = sqrt(-LL_fourMom_diff.Dot(LL_fourMom_diff));

        //L0_L0bar_Q_US_ME_hist->Fill(LL_Q);

        L0_L0_cos_theta_star_vs_Q_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), LL_Q, n_fill_LL_back_weight.at(iLambda_ME_1));

        //fill_LL_ME = 1;
/*
        nFillLL_ME_bckg++;

      }
      else //L1(LS) is from ME, L2(US) is from SE
      {
        L0_inv_mass_vs_L0_inv_mass_US_LS_ME[L_pT_bin_vector_ME_LS.at(iLambda_ME_2)][L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)]->Fill(L_vector_ME_LS.at(iLambda_ME_2).M(), L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).M(), n_fill_LL_back_weight.at(iLambda_ME_1));

        float L1_peak_mean = L1_peak_mean_fit[L_pT_bin_vector_ME_LS.at(iLambda_ME_2)][L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)];
        float L1_peak_sigma = L1_peak_sigma_fit[L_pT_bin_vector_ME_LS.at(iLambda_ME_2)][L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)];

        float L2_peak_mean = L2_peak_mean_fit[L_pT_bin_vector_ME_LS.at(iLambda_ME_2)][L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)];
        float L2_peak_sigma = L2_peak_sigma_fit[L_pT_bin_vector_ME_LS.at(iLambda_ME_2)][L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)];

        if( L_vector_ME_LS.at(iLambda_ME_2).M() < L1_peak_mean-2*L1_peak_sigma || L_vector_ME_LS.at(iLambda_ME_2).M() > L1_peak_mean+2*L1_peak_sigma ) continue;
        if( L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).M() < L1_peak_mean-2*L1_peak_sigma || L_L_L2_vector_ME_SE_back.at(iLambda_ME_1).M() > L1_peak_mean+2*L1_peak_sigma ) continue;
        if( L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).M() < L2_peak_mean-2*L2_peak_sigma || L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).M() > L2_peak_mean+2*L2_peak_sigma ) continue;

        //-------------------------

        double L_L_pairThetaStar = LpairThetaStar(L_vector_ME_LS.at(iLambda_ME_2), p_vector_ME_LS.at(iLambda_ME_2), L_L_L1_vector_ME_SE_back.at(iLambda_ME_1), L_L_p1_vector_ME_SE_back.at(iLambda_ME_1));

        L0_L0_cosThetaProdPlane_ME_LS_hist->Fill(TMath::Cos(L_L_pairThetaStar), n_fill_LL_back_weight.at(iLambda_ME_1));
        L0_L0_cosThetaProdPlane_pT_ME_LS_hist[L_pT_bin_vector_ME_LS.at(iLambda_ME_2)][L_L_L1_pT_bin_vector_ME_SE_back.at(iLambda_ME_1)]->Fill(TMath::Cos(L_L_pairThetaStar), n_fill_LL_back_weight.at(iLambda_ME_1));
        L0_L0_cosThetaProdPlane_eta_ME_LS_hist[L_eta_bin_vector_ME_LS.at(iLambda_ME_2)][L_L_L1_eta_bin_vector_ME_SE_back.at(iLambda_ME_1)]->Fill(TMath::Cos(L_L_pairThetaStar), n_fill_LL_back_weight.at(iLambda_ME_1));


        L0_L0_eta1_vs_eta2_US_LS_ME_2_hist->Fill(L_vector_ME_LS.at(iLambda_ME_2).Rapidity(), L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Rapidity(), n_fill_LL_back_weight.at(iLambda_ME_1));
        L0_L0_phi1_vs_phi2_US_LS_ME_2_hist->Fill(L_vector_ME_LS.at(iLambda_ME_2).Phi(), L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Phi(), n_fill_LL_back_weight.at(iLambda_ME_1));
        L0_L0_pT1_vs_pT2_US_LS_ME_2_hist->Fill(L_vector_ME_LS.at(iLambda_ME_2).Pt(), L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Pt(), n_fill_LL_back_weight.at(iLambda_ME_1));


        float delta_phi = 0;

        if( fabs(L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Phi() - L_vector_ME_LS.at(iLambda_ME_2).Phi()) <= TMath::Pi() )
        {
          delta_phi = fabs(L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Phi() - L_vector_ME_LS.at(iLambda_ME_2).Phi());
        }
        else
        {
          delta_phi = 2*TMath::Pi() - fabs(L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Phi() - L_vector_ME_LS.at(iLambda_ME_2).Phi());
        }


        L0_L0_delta_phi_US_LS_ME_hist->Fill( delta_phi, n_fill_LL_back_weight.at(iLambda_ME_1) );


        float delta_eta = fabs(L_L_L1_vector_ME_SE_back.at(iLambda_ME_1).Rapidity() - L_vector_ME_LS.at(iLambda_ME_2).Rapidity());


        L0_L0_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), delta_eta, n_fill_LL_back_weight.at(iLambda_ME_1));

        if( delta_eta < 0.5 )
        {
          L0_L0_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 0.5, n_fill_LL_back_weight.at(iLambda_ME_1));
        }
        else
        {
          L0_L0_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 1.5, n_fill_LL_back_weight.at(iLambda_ME_1));
        }

        L0_L0_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), delta_phi, n_fill_LL_back_weight.at(iLambda_ME_1));

        if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
        {
          //fill in-cone
          L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 0.5, n_fill_LL_back_weight.at(iLambda_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 1.5, n_fill_LL_back_weight.at(iLambda_ME_1));
        }




        float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

        L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), delta_R, n_fill_LL_back_weight.at(iLambda_ME_1));

        //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
        if( delta_R < 0.93 )
        {
          //fill in-cone
          L0_L0_cos_theta_star_vs_delta_R_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 0.5, n_fill_LL_back_weight.at(iLambda_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0_L0_cos_theta_star_vs_delta_R_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), 1.5, n_fill_LL_back_weight.at(iLambda_ME_1));
        }


        //Q bins
        TLorentzVector LL_fourMom_diff = L_L_L1_vector_ME_SE_back.at(iLambda_ME_1) - L_vector_ME_LS.at(iLambda_ME_2);

        float LL_Q = sqrt(-LL_fourMom_diff.Dot(LL_fourMom_diff));

        //L0_L0bar_Q_US_ME_hist->Fill(LL_Q);

        L0_L0_cos_theta_star_vs_Q_US_LS_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar), LL_Q, n_fill_LL_back_weight.at(iLambda_ME_1));

        nFillLL_ME_bckg++;

      }
*/
    }
  }

  //------------------------------------------------------------

  vector<float> n_fill_LbarLbar_back_weight;

  //in SE pair, firtst Lbar is US, second is LS
  for(unsigned int iLambdaBar_ME_1 = 0; iLambdaBar_ME_1 < Lbar_Lbar_Lbar1_vector_ME_SE_back.size(); iLambdaBar_ME_1++)
  {
    int fill_LbarLbar_ME = 0;

    for(unsigned int iLambdaBar_ME_2 = 0; iLambdaBar_ME_2 < Lbar_vector_ME_LS.size(); iLambdaBar_ME_2++)
    {
      float delta_phi_ME = 0;

      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi());
      }

      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).Rapidity() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).Pt() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Pt()) > 0.1 ) continue;

      float Lbar1_peak_mean = Lbar1_peak_mean_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)];
      float Lbar1_peak_sigma = Lbar1_peak_sigma_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)];

      float Lbar2_peak_mean = Lbar2_peak_mean_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)];
      float Lbar2_peak_sigma = Lbar2_peak_sigma_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)];

      if( Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).M() < Lbar1_peak_mean-2*Lbar1_peak_sigma || Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).M() > Lbar1_peak_mean+2*Lbar1_peak_sigma ) continue;
      if( Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).M() < Lbar2_peak_mean-2*Lbar2_peak_sigma || Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).M() > Lbar2_peak_mean+2*Lbar2_peak_sigma ) continue;
      if( Lbar_vector_ME_LS.at(iLambdaBar_ME_2).M() < Lbar2_peak_mean-2*Lbar2_peak_sigma || Lbar_vector_ME_LS.at(iLambdaBar_ME_2).M() > Lbar2_peak_mean+2*Lbar2_peak_sigma ) continue;


      fill_LbarLbar_ME++;

    }

    if( fill_LbarLbar_ME > 0) n_fill_LbarLbar_back_weight.push_back(1./fill_LbarLbar_ME);
    else n_fill_LbarLbar_back_weight.push_back(0);

  }


  int nFillLbarLbar_ME_backg = 0;

  //in SE pair, firtst Lbar is US, second is LS
  for(unsigned int iLambdaBar_ME_1 = 0; iLambdaBar_ME_1 < Lbar_Lbar_Lbar1_vector_ME_SE_back.size(); iLambdaBar_ME_1++)
  {
    if(n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1) == 0) continue;

    for(unsigned int iLambdaBar_ME_2 = 0; iLambdaBar_ME_2 < Lbar_vector_ME_LS.size(); iLambdaBar_ME_2++)
    {
      float delta_phi_ME = 0;

      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi()) <= TMath::Pi() )
      {
        delta_phi_ME = fabs( Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi());
      }
      else
      {
        delta_phi_ME = 2*TMath::Pi() - fabs( Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi());
      }

      //limit kinematics of ME L2 based on kinematics of same event L2 - double check precision
      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).Rapidity() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Rapidity()) > 0.1 ) continue;
      if( delta_phi_ME > 0.1 ) continue;
      if( fabs( Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).Pt() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Pt()) > 0.1 ) continue;

      //if( nFillLbarLbar_ME_backg % 2 == 0 ) //Lbar1(US) is from SE, Lbar2(LS) is from ME
      //{
        L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)]->Fill(Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).M(), Lbar_vector_ME_LS.at(iLambdaBar_ME_2).M(), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));

        float Lbar1_peak_mean = Lbar1_peak_mean_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)];
        float Lbar1_peak_sigma = Lbar1_peak_sigma_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)];

        float Lbar2_peak_mean = Lbar2_peak_mean_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)];
        float Lbar2_peak_sigma = Lbar2_peak_sigma_fit[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)];

        if( Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).M() < Lbar1_peak_mean-2*Lbar1_peak_sigma || Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).M() > Lbar1_peak_mean+2*Lbar1_peak_sigma ) continue;
        if( Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).M() < Lbar2_peak_mean-2*Lbar2_peak_sigma || Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).M() > Lbar2_peak_mean+2*Lbar2_peak_sigma ) continue;
        if( Lbar_vector_ME_LS.at(iLambdaBar_ME_2).M() < Lbar2_peak_mean-2*Lbar2_peak_sigma || Lbar_vector_ME_LS.at(iLambdaBar_ME_2).M() > Lbar2_peak_mean+2*Lbar2_peak_sigma ) continue;

        //-------------------------

        double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1), Lbar_Lbar_pBar1_vector_ME_SE_back.at(iLambdaBar_ME_1), Lbar_vector_ME_LS.at(iLambdaBar_ME_2), pBar_vector_ME_LS.at(iLambdaBar_ME_2));

        L0bar_L0bar_cosThetaProdPlane_ME_LS_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_cosThetaProdPlane_pT_ME_LS_hist[Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist[Lbar_Lbar_Lbar1_eta_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)][Lbar_eta_bin_vector_ME_LS.at(iLambdaBar_ME_2)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));


        L0bar_L0bar_eta1_vs_eta2_US_LS_ME_1_hist->Fill(Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Rapidity(), Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Rapidity(), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_phi1_vs_phi2_US_LS_ME_1_hist->Fill(Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi(), Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi(), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_pT1_vs_pT2_US_LS_ME_1_hist->Fill(Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Pt(), Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Pt(), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));


        float delta_phi = 0;

        if( fabs(Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi()) <= TMath::Pi() )
        {
          delta_phi = fabs(Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi());
        }
        else
        {
          delta_phi = 2*TMath::Pi() - fabs(Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi());
        }


        L0bar_L0bar_delta_phi_US_LS_ME_hist->Fill( delta_phi, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1) );


        float delta_eta = fabs(Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Rapidity() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Rapidity());

        L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), delta_eta, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));

        if( delta_eta < 0.5 )
        {
          L0bar_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 0.5, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        }
        else
        {
          L0bar_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 1.5, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        }

        L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), delta_phi, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));

        if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
        {
          //fill in-cone
          L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 0.5, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 1.5, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        }




        float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

        L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), delta_R, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));

        //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
        if( delta_R < 0.93 )
        {
          //fill in-cone
          L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 0.5, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 1.5, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        }


        //Q bins
        TLorentzVector LbarLbar_fourMom_diff = Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1) - Lbar_vector_ME_LS.at(iLambdaBar_ME_2);

        float LbarLbar_Q = sqrt(-LbarLbar_fourMom_diff.Dot(LbarLbar_fourMom_diff));

        //L0bar_L0barbar_Q_US_ME_hist->Fill(LbarLbar_Q);

        L0bar_L0bar_cos_theta_star_vs_Q_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), LbarLbar_Q, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
/*
        nFillLbarLbar_ME_backg++;

      }
      else //Lbar1(LS) is from ME, Lbar2(US) is from SE
      {
        L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)]->Fill(Lbar_vector_ME_LS.at(iLambdaBar_ME_2).M(), Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).M(), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));

        float Lbar1_peak_mean = Lbar1_peak_mean_fit[Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)];
        float Lbar1_peak_sigma = Lbar1_peak_sigma_fit[Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)];

        float Lbar2_peak_mean = Lbar2_peak_mean_fit[Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)];
        float Lbar2_peak_sigma = Lbar2_peak_sigma_fit[Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)];

        if( Lbar_vector_ME_LS.at(iLambdaBar_ME_2).M() < Lbar1_peak_mean-2*Lbar1_peak_sigma || Lbar_vector_ME_LS.at(iLambdaBar_ME_2).M() > Lbar1_peak_mean+2*Lbar1_peak_sigma ) continue;
        if( Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).M() < Lbar1_peak_mean-2*Lbar1_peak_sigma || Lbar_Lbar_Lbar2_vector_ME_SE_back.at(iLambdaBar_ME_1).M() > Lbar1_peak_mean+2*Lbar1_peak_sigma ) continue;
        if( Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).M() < Lbar2_peak_mean-2*Lbar2_peak_sigma || Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).M() > Lbar2_peak_mean+2*Lbar2_peak_sigma ) continue;

        //-------------------------

        double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_vector_ME_LS.at(iLambdaBar_ME_2), pBar_vector_ME_LS.at(iLambdaBar_ME_2), Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1), Lbar_Lbar_pBar1_vector_ME_SE_back.at(iLambdaBar_ME_1));

        L0bar_L0bar_cosThetaProdPlane_ME_LS_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_cosThetaProdPlane_pT_ME_LS_hist[Lbar_pT_bin_vector_ME_LS.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_pT_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist[Lbar_eta_bin_vector_ME_LS.at(iLambdaBar_ME_2)][Lbar_Lbar_Lbar1_eta_bin_vector_ME_SE_back.at(iLambdaBar_ME_1)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));


        L0bar_L0bar_eta1_vs_eta2_US_LS_ME_2_hist->Fill(Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Rapidity(), Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Rapidity(), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_phi1_vs_phi2_US_LS_ME_2_hist->Fill(Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi(), Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi(), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        L0bar_L0bar_pT1_vs_pT2_US_LS_ME_2_hist->Fill(Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Pt(), Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Pt(), n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));


        float delta_phi = 0;

        if( fabs(Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi()) <= TMath::Pi() )
        {
          delta_phi = fabs(Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi());
        }
        else
        {
          delta_phi = 2*TMath::Pi() - fabs(Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Phi() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Phi());
        }

        L0bar_L0bar_delta_phi_US_LS_ME_hist->Fill( delta_phi, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1) );


        float delta_eta = fabs(Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1).Rapidity() - Lbar_vector_ME_LS.at(iLambdaBar_ME_2).Rapidity());


        L0bar_L0bar_cos_theta_star_vs_delta_eta_scan_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), delta_eta, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));

        if( delta_eta < 0.5 )
        {
          L0bar_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 0.5, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        }
        else
        {
          L0bar_L0bar_cos_theta_star_vs_delta_eta_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 1.5, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        }

        L0bar_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), delta_phi, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));

        if( delta_eta < 0.5 && delta_phi < TMath::Pi()/3. )
        {
          //fill in-cone
          L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 0.5, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 1.5, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        }


        float delta_R = sqrt( delta_eta*delta_eta + delta_phi*delta_phi );

        L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), delta_R, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));

        //Delta R = 0.93 is about the same as delta_eta < 0.5 && delta_phi < TMath::Pi()/3.
        if( delta_R < 0.93 )
        {
          //fill in-cone
          L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 0.5, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        }
        else
        {
          //fill out-of-cone
          L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), 1.5, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));
        }

        //Q bins
        TLorentzVector LbarLbar_fourMom_diff = Lbar_Lbar_Lbar1_vector_ME_SE_back.at(iLambdaBar_ME_1) - Lbar_vector_ME_LS.at(iLambdaBar_ME_2);

        float LbarLbar_Q = sqrt(-LbarLbar_fourMom_diff.Dot(LbarLbar_fourMom_diff));

        //L0bar_L0barbar_Q_US_ME_hist->Fill(LbarLbar_Q);

        L0bar_L0bar_cos_theta_star_vs_Q_US_LS_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar), LbarLbar_Q, n_fill_LbarLbar_back_weight.at(iLambdaBar_ME_1));

        nFillLbarLbar_ME_backg++;

      }
*/
    }

  }


  LLbarOutFile->Write();
  LLbarOutFile->Close();

  //________________________________________________________________________________________________________

  return true;
}
//__________________________________________________________________________________________________________________________

//ReadMode = 0 - read TTree, ReadMode = 1 - read histograms - First run in ReadMode = 0 to save relevant histograms, then can run in ReadMode = 1 to read just histograms and save time
//cut_type = 0 - analysis cuts, 1 - tight topological cuts, 2 - tight daughter pT cut
//energy - collision energy in GeV
void Ana003_Lambda_corr_2D(const int ReadMode = 0, const int cut_type = 0, const int energy = 510, const int year = 2017)
{
  ifstream fileList;

  if(energy == 510)
  {
    //MB without strict TOF matching
    fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run17_MB_noTOF/fileList.list");

  }
  else if(energy == 200)
  {
    if(year == 2012)
    {
      //MB without strict TOF matching
      //fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF_open_cuts/fileList.list");

      //fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF_open_cuts_new/fileList.list");
      //fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF_open_cuts_daughter_tag/fileList.list");
      fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF_open_cuts_daughter_tag_helix/fileList.list");
      //fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF_open_cuts_daughter_tag_straight_line/fileList.list");

      //fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF/fileList.list");
      //fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF_new/fileList.list");
    }

    if(year == 2015)
    {
      //fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run15_MB/fileList.list");
      fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run15_MB_TOF_match_pi/fileList.list");
    }

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

  TChain *myChain = new TChain("ntp_Lambda");

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

  cout<<"nEvents from hist = "<<hEventStat1->GetBinContent(6)<<endl;


  TCanvas *testCan = new TCanvas("testCan", "testCan", 1200, 1000);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  testCan->cd();
  hEventStat1->GetXaxis()->SetRange(1,6);
  hEventStat1->GetYaxis()->SetMaxDigits(3);
  hEventStat1->SetMinimum(0);
  hEventStat1->Draw();

  testCan->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/EventStat_test.png");

  //return; //for reploting of event stats

  bool AnaFinish = DoAnalysis(myChain, ReadMode, cut_type, energy , year);

  if(!AnaFinish)
  {
    cout<<"Analysis ended abnormally. Aborting!"<<endl;

    return;
  }


  cout<<endl;
  cout<<"Nubmer of accepted events: "<<hEventStat1->GetBinContent(6)<<endl;


  return;
}
