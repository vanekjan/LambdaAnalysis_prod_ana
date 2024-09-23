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

//bool cuts(int strictTOF_cut, int pi_hasTOFinfo, int p_hasTOFinfo, float L_y)
bool cuts(int energy, float L_y, int hasTOFinfo)
{

  if( !( TMath::Abs(L_y) < L_y_cut ) ) return false;
  if( energy == 510 && hasTOFinfo == 0 ) return false; //TOF matched daughter (choose when calling cuts())
  //if( strictTOF_cut == 2 && (pi_hasTOFinfo == 0 || p_hasTOFinfo == 0) ) return false; //TOF matched pions and protons
  //if(cos(L_theta) < L_cos_theta_cut) return false;
  //if(L_decayL > L_decayL_cut) return false;

  return true;

}


//for tight topological cuts for sys. err.
bool cuts_topo_sys_err(float p_DCA, float pi_DCA, float pair_DCA, float dec_L, float cos_theta)
{
  if( p_DCA < 0.4 ) return false;
  if( pi_DCA < 0.2 ) return false;

  if( pair_DCA > 0.9 ) return false;

  if( dec_L < 3 ) return false;

  if( cos_theta < 0.997 ) return false;

  return true;
}

//for tight daughter pT cut for sys. err
bool cuts_pT_sys_err(float p_pT, float pi_pT)
{
  if( p_pT < 0.3 ) return false;
  if( pi_pT < 0.2 ) return false; //0.3 migt be too tight for pions

  return true;
}

bool Ana003_Lambda_corr_2D_get_corr_Delta_R(const int cut_type = 0, const int energy = 510, const int year = 2017)
{
  //analyze stored Lambda pairs and save cos(theta*) histograms

  const float L0_alpha = 0.732; //decay parameter of L0
  const float L0_alpha_relat_err = 0.014/L0_alpha; //relative error of decay parameter

  const float L0bar_alpha = -0.758; //decay paramteter of L0bar
  const float L0bar_alpha_relat_err = 0.012/fabs(L0bar_alpha); //relative error of decay paramteter


  //_______________________________________________________________________________________________________________________________________________

  //systematic uncertainties
  //residual effect from PYTHIA closure test
  TFile *SysErrSlopeFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErrSlope_nocorr.root", year), "read");

  TH1F *SysErrSlope_delta_R_hist[2];

  for( unsigned int delta_R_bin = 0; delta_R_bin < 2; delta_R_bin++ )
  {
    //SysErrSlope_delta_R_hist[delta_R_bin] = (TH1F*)SysErrSlopeFile->Get(Form("SysErrSlope_delta_R_hist_%i", delta_R_bin));
    SysErrSlope_delta_R_hist[delta_R_bin] = (TH1F*)SysErrSlopeFile->Get(Form("SysErrSlope_delta_eta_delta_phi_hist_%i", delta_R_bin));
  }

  //systematic uncertainty histograms and values
  //alpha
  float sysErr_alpha_L0_L0bar = sqrt(L0_alpha_relat_err*L0_alpha_relat_err + L0bar_alpha_relat_err*L0bar_alpha_relat_err);
  float sysErr_alpha_L0_L0 = sqrt(L0_alpha_relat_err*L0_alpha_relat_err + L0_alpha_relat_err*L0_alpha_relat_err);
  float sysErr_alpha_L0bar_L0bar = sqrt(L0bar_alpha_relat_err*L0bar_alpha_relat_err + L0bar_alpha_relat_err*L0bar_alpha_relat_err);


  //cuts variation
  //have to run this code with cuts_type = 1 and 2 first
  TFile *SysErrCutsTopo;
  TFile *SysErrCutsPt;

  if(cut_type == 0 )
  {
    //SysErrCutsTopo = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_topo_cuts_delta_R.root", year), "read");
    SysErrCutsTopo = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_topo_cuts_Delta_y_Delta_phi_work.root", year), "read");

    //SysErrCutsPt = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_pT_cuts_delta_R.root", year), "read");
    SysErrCutsPt = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_pT_cuts_Delta_y_Delta_phi_work.root", year), "read");
  }
  else
  {
    if( cut_type == 1 ) SysErrCutsTopo = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_topo_cuts_delta_R_work.root", year), "recreate");
    if( cut_type == 2 ) SysErrCutsPt = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_pT_cuts_delta_R_work.root", year), "recreate");
  }

  //tight cuts
  TF1 *fitL0_L0bar_tight_topo_cuts_delta_R[2];
  TF1 *fitL0_L0_tight_topo_cuts_delta_R[2];
  TF1 *fitL0bar_L0bar_tight_topo_cuts_delta_R[2];

  TF1 *fitL0_L0bar_tight_pT_cuts_delta_R[2];
  TF1 *fitL0_L0_tight_pT_cuts_delta_R[2];
  TF1 *fitL0bar_L0bar_tight_pT_cuts_delta_R[2];

  for( unsigned int delta_R_bin = 0; delta_R_bin < 2; delta_R_bin++)
  {
    if(cut_type == 0 )
    {
      /*
      fitL0_L0bar_tight_topo_cuts_delta_R[delta_R_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0_L0bar_tight_topo_cuts_delta_R_%i", delta_R_bin));
      fitL0_L0_tight_topo_cuts_delta_R[delta_R_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0_L0_tight_topo_cuts_delta_R_%i", delta_R_bin));
      fitL0bar_L0bar_tight_topo_cuts_delta_R[delta_R_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0bar_L0bar_tight_topo_cuts_delta_R_%i", delta_R_bin));

      fitL0_L0bar_tight_pT_cuts_delta_R[delta_R_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0_L0bar_tight_pT_cuts_delta_R_%i", delta_R_bin));
      fitL0_L0_tight_pT_cuts_delta_R[delta_R_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0_L0_tight_pT_cuts_delta_R_%i", delta_R_bin));
      fitL0bar_L0bar_tight_pT_cuts_delta_R[delta_R_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0bar_L0bar_tight_pT_cuts_delta_R_%i", delta_R_bin));
      */
      fitL0_L0bar_tight_topo_cuts_delta_R[delta_R_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_%i", delta_R_bin));
      fitL0_L0_tight_topo_cuts_delta_R[delta_R_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_%i", delta_R_bin));
      fitL0bar_L0bar_tight_topo_cuts_delta_R[delta_R_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_%i", delta_R_bin));

      fitL0_L0bar_tight_pT_cuts_delta_R[delta_R_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_%i", delta_R_bin));
      fitL0_L0_tight_pT_cuts_delta_R[delta_R_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_%i", delta_R_bin));
      fitL0bar_L0bar_tight_pT_cuts_delta_R[delta_R_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_%i", delta_R_bin));

    }
    else
    {
      fitL0_L0bar_tight_topo_cuts_delta_R[delta_R_bin] = new TF1(Form("fitL0_L0bar_tight_topo_cuts_delta_R_%i", delta_R_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_tight_topo_cuts_delta_R[delta_R_bin]->SetParameters(100, 0.5);

      fitL0_L0_tight_topo_cuts_delta_R[delta_R_bin] = new TF1(Form("fitL0_L0_tight_topo_cuts_delta_R_%i", delta_R_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_tight_topo_cuts_delta_R[delta_R_bin]->SetParameters(100, 0.5);

      fitL0bar_L0bar_tight_topo_cuts_delta_R[delta_R_bin] = new TF1(Form("fitL0bar_L0bar_tight_topo_cuts_delta_R_%i", delta_R_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_tight_topo_cuts_delta_R[delta_R_bin]->SetParameters(100, 0.5);


      fitL0_L0bar_tight_pT_cuts_delta_R[delta_R_bin] = new TF1(Form("fitL0_L0bar_tight_pT_cuts_delta_R_%i", delta_R_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_tight_pT_cuts_delta_R[delta_R_bin]->SetParameters(100, 0.5);

      fitL0_L0_tight_pT_cuts_delta_R[delta_R_bin] = new TF1(Form("fitL0_L0_tight_pT_cuts_delta_R_%i", delta_R_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_tight_pT_cuts_delta_R[delta_R_bin]->SetParameters(100, 0.5);

      fitL0bar_L0bar_tight_pT_cuts_delta_R[delta_R_bin] = new TF1(Form("fitL0bar_L0bar_tight_pT_cuts_delta_R_%i", delta_R_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_tight_pT_cuts_delta_R[delta_R_bin]->SetParameters(100, 0.5);
    }

  }


  //----------------------------------------------------------------------

  //output file with polarization graphs
  TFile *out_file = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/Polarization/%i/Polarization_delta_R.root", year), "recreate");


  TFile *LLbarOutFile; //input file with data histograms

  if(cut_type == 0)
  {
    if(year == 2012)
    {
      //Run12
      //work file with any recent updates
      //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_work.root", year), "read"); //analysis cuts

      //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/2024_06_new_ME/ProdPlane_Lambda_ana_cuts_new_ME.root", year), "read"); //analysis cuts

      //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/20240729_tests_n_sigma_signal/ProdPlane_Lambda_ana_cuts_2_sigma_signal_alt_ME_Delta_pi_third_new_scan_hists_fixed_bckg.root", year), "read"); //analysis cuts
      //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/20240729_tests_n_sigma_signal/ProdPlane_Lambda_ana_cuts_3_sigma_signal_alt_ME_Delta_pi_third_new_scan_hists_fixed_bckg.root", year), "read"); //analysis cuts

      //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/cuts_tests/ProdPlane_Lambda_ana_cuts_Delta_phi_pi_5.root", year), "read"); //analysis cuts
      //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/cuts_tests/ProdPlane_Lambda_ana_cuts_Delta_y_0.4_Delta_R_0.85.root", year), "read"); //analysis cuts
      //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/cuts_tests/ProdPlane_Lambda_ana_cuts_Delta_y_0.6_Delta_R_1.5.root", year), "read"); //analysis cuts
      //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/cuts_tests/ProdPlane_Lambda_ana_cuts_ana_cuts_R_scan_hist.root", year), "read"); //analysis cuts

      LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/Delta_R_tests/ProdPlane_Lambda_ana_cuts_Delta_R_0.93.root", year), "read"); //analysis cuts
      //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/Delta_R_tests/ProdPlane_Lambda_ana_cuts_Delta_R_0.75.root", year), "read"); //analysis cuts
      //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/Delta_R_tests/ProdPlane_Lambda_ana_cuts_Delta_R_0.5.root", year), "read"); //analysis cuts


    }
    else if(year == 2015)
    {
      LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_work.root", year), "read"); //analysis cuts

      //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_new_no_TOF.root", year), "read"); //analysis cuts
    }
    else
    {
      //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts.root", year), "read"); //analysis cuts
      LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_work.root", year), "read"); //analysis cuts
    }

    //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  }
  else if(cut_type == 1)
  {
    LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_topo_tight_cuts_work.root", year), "read");

    //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_topo_tight_cuts.root", year), "read");

    //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_topo_tight_cuts_new_prod.root", year), "read");
  }
  else if(cut_type == 2)
  {
    LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_pT_tight_cuts_work.root", year), "read");

    //LLbarOutFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/tests_Delta_phi/ProdPlane_Lambda_pT_tight_cuts_delta_phi_pi_pT_170.root", year), "read");
  }
  else
  {
    cout<<"Wrong cut type"<<endl;

    return false;
  }


  //data histograms
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_US_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_R_US_hist");

  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_US_LS_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_R_US_LS_hist");

  //--------------------------------------------------------

  TH2F *L0_L0_cos_theta_star_vs_delta_R_US_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_R_US_hist");

  TH2F *L0_L0_cos_theta_star_vs_delta_R_US_LS_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_R_US_LS_hist");

  //--------------------------------------------------------

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_US_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_R_US_hist");

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_hist");

  //--------------------------------------------------------

  //mixed event
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_US_ME_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_R_US_ME_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist");

  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_US_ME_alt_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_R_US_ME_alt_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_alt_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_alt_hist");

  //--------------------------------------------------------

  TH2F *L0_L0_cos_theta_star_vs_delta_R_US_ME_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_R_US_ME_hist");
  TH2F *L0_L0_cos_theta_star_vs_delta_R_US_LS_ME_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_R_US_LS_ME_hist");

  TH2F *L0_L0_cos_theta_star_vs_delta_R_US_ME_alt_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_R_US_ME_alt_hist");
  TH2F *L0_L0_cos_theta_star_vs_delta_R_US_LS_ME_alt_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_R_US_LS_ME_alt_hist");

  //--------------------------------------------------------

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_US_ME_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_R_US_ME_hist");
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist");

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_US_ME_alt_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_R_US_ME_alt_hist");
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_alt_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_alt_hist");

  //________________________________________________________________________________________



  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(42);

  //polarization graph for Delta y
  //first graph is for |Delta y| < 0.5, second is for 0.5 < |Delta y| < 2.0
  TGraphErrors *PolarizationGraph_delta_R[2];
  TGraphErrors *PolarizationGraph_delta_R_sys_err[2];


  TCanvas *PolarizationGraph_delta_R_can = new TCanvas("PolarizationGraph_delta_R_can", "PolarizationGraph_delta_R_can", 2000, 1200);
  PolarizationGraph_delta_R_can->Divide(2,1);

  TH1F *DefaultHist = new TH1F("DefaultHist", "DefaultHist", 10, -0.5, 0.5);
  DefaultHist->GetXaxis()->SetTitle("P_{#Lambda_{1}#Lambda_{2}}");
  DefaultHist->GetXaxis()->CenterTitle();
  DefaultHist->GetXaxis()->SetRangeUser(-0.42, 0.42);
  //DefaultHist->GetYaxis()->SetRangeUser(0.01, 4.99); //with K0s
  DefaultHist->GetYaxis()->SetRangeUser(0.01, 3.99); //without K0s
  //DefaultHist->GetYaxis()->SetNdivisions(505);
  //DefaultHist->GetYaxis()->SetNdivisions(06); //with K0s
  DefaultHist->GetYaxis()->SetNdivisions(06); //without K0s
  //DefaultHist->GetYaxis()->SetNdivisions(-5);
  DefaultHist->GetYaxis()->ChangeLabel(1,-1,-1,-1,-1,-1,"#Lambda#bar{#Lambda}");
  DefaultHist->GetYaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"#Lambda#Lambda");
  DefaultHist->GetYaxis()->ChangeLabel(3,-1,-1,-1,-1,-1,"#bar{#Lambda}#bar{#Lambda}");
  //DefaultHist->GetYaxis()->ChangeLabel(4,-1,-1,-1,-1,-1,"K_{s}^{0}K_{s}^{0}");

  TLine *ZeroLine_eta = new TLine(0,0.01,0,3.99);
  ZeroLine_eta->SetLineStyle(9);
  ZeroLine_eta->SetLineColor(1);
  //ZeroLine_eta->Draw("same");

  //TPaveText *Polarization_text = new TPaveText(0.12, 0.58, 0.49, 0.89, "NDC"); //with K0s
  TPaveText *Polarization_text = new TPaveText(0.55, 0.75, 0.89, 0.89, "NDC"); //without K0s
  Polarization_text->SetTextFont(42);
  //Polarization_text->AddText("STAR");
  //Polarization_text->AddText("STAR preliminary");
  //((TText*)Polarization_text->GetListOfLines()->Last())->SetTextColor(2);
  Polarization_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  Polarization_text->AddText("Minimum bias");
  Polarization_text->AddText("|#it{y}| < 1");
  Polarization_text->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
  Polarization_text->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
  Polarization_text->SetFillColorAlpha(0, 0.01);
  //Polarization_text->Draw("same");

  TPaveText *bin_text[2];

  out_file->cd();

  for( unsigned int delta_R_bin = 1; delta_R_bin < 3; delta_R_bin++)
  {
    //create polarization graphs, defined earlier

    PolarizationGraph_delta_R[delta_R_bin-1] = new TGraphErrors(3);
    PolarizationGraph_delta_R_sys_err[delta_R_bin-1] = new TGraphErrors(3);

    //-----------------------------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_eta_no_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_eta_no_corr_can_%i", delta_R_bin), Form("L0_L0bar_cosThetaProdPlane_eta_no_corr_can_%i", delta_R_bin), 1200, 1000);

    L0_L0bar_cosThetaProdPlane_eta_no_corr_can->cd();

    TH1D *L0_L0bar_cosThetaProdPlane_eta_US_hist = L0_L0bar_cos_theta_star_vs_delta_R_US_hist->ProjectionX( Form("proj_US_%i", delta_R_bin), delta_R_bin, delta_R_bin);

    L0_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->CenterTitle();
    //L0_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTextSizePixels(30);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->CenterTitle();
    //L0_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTextSizePixels(30);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetMaxDigits(3);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetMarkerStyle(20);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetMarkerColor(kRed);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetLineColor(kRed);
    double nLLbar = L0_L0bar_cosThetaProdPlane_eta_US_hist->Integral();
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Sumw2();
    //L0_L0bar_cosThetaProdPlane_eta_US_hist->Divide(L0_L0bar_cosThetaProdPlane_eta_eff);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1));
    //L0_L0bar_cosThetaProdPlane_eta_US_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_US_hist->Integral());
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //ME here just for plotting, used loser
    TH1D *L0_L0bar_cosThetaProdPlane_eta_ME_hist = L0_L0bar_cos_theta_star_vs_delta_R_US_ME_hist->ProjectionX( Form("proj_US_ME_%i", delta_R_bin), delta_R_bin, delta_R_bin);

    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerStyle(24);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerColor(1);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetLineColor(1);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Sumw2();
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Scale(nLLbar/L0_L0bar_cosThetaProdPlane_eta_ME_hist->Integral()); //scale ME to US
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_ME_hist->GetXaxis()->GetBinWidth(1)); //US is scaled by bin width, ME needs to be scaled as well
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Draw("p e same");

    TF1 *fitL0_L0bar_US_ThetaStar_no_corr_ME = new TF1("fitL0_L0bar_US_ThetaStar_no_corr_ME", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0bar_US_ThetaStar_no_corr_ME->SetParameters(100, 0.5);

    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Fit(fitL0_L0bar_US_ThetaStar_no_corr_ME, "s i 0 r");

    float P_L0_L0bar_no_corr_ME = fitL0_L0bar_US_ThetaStar_no_corr_ME->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_no_corr_ME_err = fitL0_L0bar_US_ThetaStar_no_corr_ME->GetParError(1)/(L0_alpha*L0bar_alpha);


    TH1D *L0_L0bar_cosThetaProdPlane_eta_LS_hist = L0_L0bar_cos_theta_star_vs_delta_R_US_LS_hist->ProjectionX( Form("proj_LS_%i", delta_R_bin), delta_R_bin, delta_R_bin);

    L0_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerStyle(21);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerColor(kBlue);
    double nLLbar_back = L0_L0bar_cosThetaProdPlane_eta_LS_hist->Integral();
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Sumw2();
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1));
    //L0_L0bar_cosThetaProdPlane_eta_LS_hist->Divide(L0_L0bar_cosThetaProdPlane_eta_eff);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Draw("p e same");


    TH1D *L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist = L0_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist->ProjectionX( Form("proj_LS_ME_%i", delta_R_bin), delta_R_bin, delta_R_bin);

    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerStyle(25);
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerColor(kMagenta+1);
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetLineColor(kMagenta+1);
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Sumw2();
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Scale(nLLbar_back/L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Integral()); //scale ME_LS to background
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->GetXaxis()->GetBinWidth(1)); //background is scaled by bin width, ME_LS needs to be scaled as well
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Draw("p e same");


    TF1 *fitL0_L0bar_US_ThetaStar_no_corr = new TF1("fitL0_L0bar_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0bar_US_ThetaStar_no_corr->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Fit(fitL0_L0bar_US_ThetaStar_no_corr, "s i 0 r");

    float P_L0_L0bar_no_corr = fitL0_L0bar_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_no_corr_err = fitL0_L0bar_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

    fitL0_L0bar_US_ThetaStar_no_corr->SetLineColor(1);
    //fitL0_L0bar_US_ThetaStar_no_corr->Draw("same");

    TLegend *L0_L0bar_leg = new TLegend(0.15, 0.45, 0.45, 0.69);
    L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_eta_US_hist, "(US-US) p#pi");
    L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_eta_ME_hist, "(US-US) ME");
    L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_eta_LS_hist, "Combinatorial bckg.");
    L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist, "Bckg. ME");
    //L0_L0bar_leg->AddEntry(fitL0_L0bar_US_ThetaStar_no_corr, "Linear fit to US");
    L0_L0bar_leg->SetBorderSize(0);
    L0_L0bar_leg->SetFillColorAlpha(0, 0.01);
    L0_L0bar_leg->Draw("same");

    TPaveText *L0_L0bar_text_no_corr = new TPaveText(0.5, 0.4, 0.85, 0.75, "NDC");
    L0_L0bar_text_no_corr->SetTextFont(42);
    //L0_L0bar_text_no_corr->AddText("STAR Internal");
    //L0_L0bar_text_no_corr->AddText("STAR preliminary");
    //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
    L0_L0bar_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0_L0bar_text_no_corr->AddText("Minimum bias, no correction");
    L0_L0bar_text_no_corr->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
    L0_L0bar_text_no_corr->AddText("|#it{y}| < 1");
    L0_L0bar_text_no_corr->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
    L0_L0bar_text_no_corr->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
    L0_L0bar_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar_no_corr, fabs(P_L0_L0bar_no_corr_err)));
    L0_L0bar_text_no_corr->AddText(Form("P_{ME} = %.3f #pm %.3f", P_L0_L0bar_no_corr_ME, fabs(P_L0_L0bar_no_corr_ME_err)));
    L0_L0bar_text_no_corr->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_no_corr->Draw("same");

    L0_L0bar_cosThetaProdPlane_eta_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations_delta_R/L0_L0bar_cosThetaProdPlane_eta_no_corr_less_delta_R_%i.png", delta_R_bin));

    //----------------------------------------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_eta_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_eta_can_%i", delta_R_bin), Form("L0_L0bar_cosThetaProdPlane_eta_can_%i", delta_R_bin), 1200, 1000);

    L0_L0bar_cosThetaProdPlane_eta_can->cd();

    //ME histogram higher

    L0_L0bar_cosThetaProdPlane_eta_US_hist->Divide(L0_L0bar_cosThetaProdPlane_eta_ME_hist); // correct US using ME
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Scale(nLLbar/L0_L0bar_cosThetaProdPlane_eta_US_hist->Integral()); //scale back to raw US
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //L0_L0bar_cosThetaProdPlane_eta_ME_hist->Draw("same p e");

    TH1D *L0_L0bar_cosThetaProdPlane_eta_LS_hist_clone = (TH1D*)L0_L0bar_cosThetaProdPlane_eta_LS_hist->Clone("L0_L0bar_cosThetaProdPlane_eta_LS_hist_clone");
    L0_L0bar_cosThetaProdPlane_eta_LS_hist_clone->Divide(L0_L0bar_cosThetaProdPlane_eta_ME_hist); // correct US using ME

    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Divide(L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist); //correct background using ME
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Scale(nLLbar_back/L0_L0bar_cosThetaProdPlane_eta_LS_hist->Integral()); //scale back to raw background
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Draw("p e same");

    //fit dN/dcos(theta*)
    //signal + bacground
    TF1 *fitL0_L0bar_US_ThetaStar = new TF1("fitL0_L0bar_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0bar_US_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Fit(fitL0_L0bar_US_ThetaStar, "s i 0 r");

    float P_L0_L0bar = fitL0_L0bar_US_ThetaStar->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_err = fitL0_L0bar_US_ThetaStar->GetParError(1)/(L0_alpha*L0bar_alpha);

    fitL0_L0bar_US_ThetaStar->SetLineColor(1);
    fitL0_L0bar_US_ThetaStar->Draw("same");

    //background
    TF1 *fitL0_L0bar_US_LS_ThetaStar = new TF1("fitL0_L0bar_US_LS_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0bar_US_LS_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaUS_LS_wrong_sign = L_inv_mass_US_LS->Fit(fitGaUS_LSsBack, "s i 0", "", 1.07, 1.4);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Fit(fitL0_L0bar_US_LS_ThetaStar, "s i 0 r");

    float P_L0_L0bar_back = fitL0_L0bar_US_LS_ThetaStar->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_back_err = fitL0_L0bar_US_LS_ThetaStar->GetParError(1)/(L0_alpha*L0bar_alpha);

    fitL0_L0bar_US_LS_ThetaStar->SetLineColor(1);
    fitL0_L0bar_US_LS_ThetaStar->Draw("same");


    TPaveText *L0_L0bar_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
    L0_L0bar_text->SetTextFont(42);
    //L0_L0bar_text->AddText("STAR Internal");
    //L0_L0bar_text->AddText("STAR preliminary");
    //((TText*)L0_L0bar_text->GetListOfLines()->Last())->SetTextColor(2);
    L0_L0bar_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0_L0bar_text->AddText("Minimum bias");
    L0_L0bar_text->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
    L0_L0bar_text->AddText("|#it{y}| < 1");
    L0_L0bar_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
    L0_L0bar_text->AddText(Form("P_{tot} = %.2f #pm %.2f", P_L0_L0bar, fabs(P_L0_L0bar_err)));
    L0_L0bar_text->AddText(Form("P_{bckg} = %.2f #pm %.2f", P_L0_L0bar_back, fabs(P_L0_L0bar_back_err)));
    L0_L0bar_text->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text->Draw("same");

    L0_L0bar_leg->Draw("same");

    L0_L0bar_cosThetaProdPlane_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations_delta_R/L0_L0bar_cosThetaProdPlane_eta_less_delta_R_%i.png", delta_R_bin));

    //----------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_eta_can_2 = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_eta_can_2_%i", delta_R_bin), Form("L0_L0bar_cosThetaProdPlane_eta_can_2_%i", delta_R_bin), 1200, 1000);

    L0_L0bar_cosThetaProdPlane_eta_can_2->cd();

    L0_L0bar_cosThetaProdPlane_eta_US_hist->Add(L0_L0bar_cosThetaProdPlane_eta_LS_hist, -1);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Write(Form("L0_L0bar_cosThetaProdPlane_delta_R_%i", delta_R_bin-1));
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");


    TF1 *fitL0_L0bar_US_ThetaStar_2 = new TF1("fitL0_L0bar_US_ThetaStar_2", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0bar_US_ThetaStar_2->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Fit(fitL0_L0bar_US_ThetaStar_2, "s i 0 r");


    //store fit result for systematic errors
    //tight topo cuts
    if(cut_type == 1)
    {
      SysErrCutsTopo->cd();

      fitL0_L0bar_tight_topo_cuts_delta_R[delta_R_bin-1]->SetParameters(fitL0_L0bar_US_ThetaStar_2->GetParameter(0), fitL0_L0bar_US_ThetaStar_2->GetParameter(1));
      fitL0_L0bar_tight_topo_cuts_delta_R[delta_R_bin-1]->SetParError(0, fitL0_L0bar_US_ThetaStar_2->GetParError(0));
      fitL0_L0bar_tight_topo_cuts_delta_R[delta_R_bin-1]->SetParError(1, fitL0_L0bar_US_ThetaStar_2->GetParError(1));
      fitL0_L0bar_tight_topo_cuts_delta_R[delta_R_bin-1]->Write();
    }

    //tight pT cuts
    if(cut_type == 2)
    {
      SysErrCutsPt->cd();

      fitL0_L0bar_tight_pT_cuts_delta_R[delta_R_bin-1]->SetParameters(fitL0_L0bar_US_ThetaStar_2->GetParameter(0), fitL0_L0bar_US_ThetaStar_2->GetParameter(1));
      fitL0_L0bar_tight_pT_cuts_delta_R[delta_R_bin-1]->SetParError(0, fitL0_L0bar_US_ThetaStar_2->GetParError(0));
      fitL0_L0bar_tight_pT_cuts_delta_R[delta_R_bin-1]->SetParError(1, fitL0_L0bar_US_ThetaStar_2->GetParError(1));
      fitL0_L0bar_tight_pT_cuts_delta_R[delta_R_bin-1]->Write();
    }

    float P_L0_L0bar_2 = fitL0_L0bar_US_ThetaStar_2->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_err_2 = fitL0_L0bar_US_ThetaStar_2->GetParError(1)/(L0_alpha*L0bar_alpha);

    fitL0_L0bar_US_ThetaStar_2->SetLineColor(1);
    fitL0_L0bar_US_ThetaStar_2->Draw("same");


    //calculate total systematic uncertainty for L-Lbar

    //slope
    float SysErrSlope_L0_L0bar = SysErrSlope_delta_R_hist[delta_R_bin-1]->GetBinContent(1);

    //----------------------------

    //background subtraction
    float P_L0_L0bar_from_fits = P_L0_L0bar - nLLbar_back/nLLbar*P_L0_L0bar_back;
    float P_L0_L0bar_from_fits_err = sqrt( P_L0_L0bar_err*P_L0_L0bar_err + nLLbar_back/nLLbar*nLLbar_back/nLLbar*P_L0_L0bar_back_err*P_L0_L0bar_back_err );

    //statistical error correction for systematic error
    float SysErrBackground_L0_L0bar_corr = sqrt( fabs( P_L0_L0bar_err_2*P_L0_L0bar_err_2 - P_L0_L0bar_from_fits_err*P_L0_L0bar_from_fits_err ) );

    float SysErrBackground_L0_L0bar_work = ( fabs( P_L0_L0bar_2 - P_L0_L0bar_from_fits) - SysErrBackground_L0_L0bar_corr )/fabs(P_L0_L0bar_2);

    float SysErrBackground_L0_L0bar = 0;

    if( SysErrBackground_L0_L0bar_work > 0 ) SysErrBackground_L0_L0bar = SysErrBackground_L0_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations

    //----------------------------

    //cuts variation
    float SysErr_tight_cuts_L0_L0bar = 0;

    //float P_L0_L0bar_tight_topo_cuts = 0;
    //float P_L0_L0bar_tight_pT_cuts = 0;


    if(cut_type == 0)
    {
      float P_L0_L0bar_tight_topo_cuts = fitL0_L0bar_tight_topo_cuts_delta_R[delta_R_bin-1]->GetParameter(1)/(L0_alpha*L0bar_alpha);
      float P_L0_L0bar_tight_topo_cuts_err = fitL0_L0bar_tight_topo_cuts_delta_R[delta_R_bin-1]->GetParError(1)/(L0_alpha*L0bar_alpha);

            cout<<"hovno"<<endl;
      //statistical error correction for systematic error
      float SysErr_tight_topo_cuts_L0_L0bar_corr = sqrt( fabs( P_L0_L0bar_err_2*P_L0_L0bar_err_2 - P_L0_L0bar_tight_topo_cuts_err*P_L0_L0bar_tight_topo_cuts_err ) );

      float SysErr_tight_topo_cuts_L0_L0bar_work = ( fabs( P_L0_L0bar_2 - P_L0_L0bar_tight_topo_cuts) - SysErr_tight_topo_cuts_L0_L0bar_corr )/fabs(P_L0_L0bar_2);

      float SysErr_tight_topo_cuts_L0_L0bar = 0;

      if( SysErr_tight_topo_cuts_L0_L0bar_work > 0 ) SysErr_tight_topo_cuts_L0_L0bar = SysErr_tight_topo_cuts_L0_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations


      float P_L0_L0bar_tight_pT_cuts = fitL0_L0bar_tight_pT_cuts_delta_R[delta_R_bin-1]->GetParameter(1)/(L0_alpha*L0bar_alpha);
      float P_L0_L0bar_tight_pT_cuts_err = fitL0_L0bar_tight_pT_cuts_delta_R[delta_R_bin-1]->GetParError(1)/(L0_alpha*L0bar_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_pT_cuts_L0_L0bar_corr = sqrt( fabs( P_L0_L0bar_err_2*P_L0_L0bar_err_2 - P_L0_L0bar_tight_pT_cuts_err*P_L0_L0bar_tight_pT_cuts_err ) );

      float SysErr_tight_pT_cuts_L0_L0bar_work = ( fabs( P_L0_L0bar_2 - P_L0_L0bar_tight_pT_cuts) - SysErr_tight_pT_cuts_L0_L0bar_corr )/fabs(P_L0_L0bar_2);

      float SysErr_tight_pT_cuts_L0_L0bar = 0;

      if( SysErr_tight_pT_cuts_L0_L0bar_work > 0 ) SysErr_tight_pT_cuts_L0_L0bar = SysErr_tight_pT_cuts_L0_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_cuts_L0_L0bar = sqrt(SysErr_tight_topo_cuts_L0_L0bar*SysErr_tight_topo_cuts_L0_L0bar + SysErr_tight_pT_cuts_L0_L0bar*SysErr_tight_pT_cuts_L0_L0bar);
      //SysErr_tight_cuts_L0_L0bar = SysErr_tight_topo_cuts_L0_L0bar;
      //SysErr_tight_cuts_L0_L0bar = SysErr_tight_pT_cuts_L0_L0bar;
    }

    //total
    //float SysErrTot_L0_L0bar = sqrt( sysErr_alpha_L0_L0bar*sysErr_alpha_L0_L0bar + SysErrSlope_L0_L0bar*SysErrSlope_L0_L0bar + SysErrBackground_L0_L0bar*SysErrBackground_L0_L0bar  );
    float SysErrTot_L0_L0bar = sqrt( sysErr_alpha_L0_L0bar*sysErr_alpha_L0_L0bar + SysErrSlope_L0_L0bar*SysErrSlope_L0_L0bar + SysErrBackground_L0_L0bar*SysErrBackground_L0_L0bar + SysErr_tight_cuts_L0_L0bar*SysErr_tight_cuts_L0_L0bar  );


    //TPaveText *L0_L0bar_text_2 = new TPaveText(0.47, 0.15, 0.85, 0.53, "NDC"); //with polarization value
    TPaveText *L0_L0bar_text_2 = new TPaveText(0.47, 0.15, 0.8, 0.53, "NDC"); //without polarization value
    L0_L0bar_text_2->SetTextFont(42);
    //L0_L0bar_text_2->SetTextSize(15);
    //L0_L0bar_text_2->AddText("STAR");
    //L0_L0bar_text_2->AddText("STAR preliminary");
    //((TText*)L0_L0bar_text_2->GetListOfLines()->Last())->SetTextColor(2);
    L0_L0bar_text_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0_L0bar_text_2->AddText("Minimum bias");
    L0_L0bar_text_2->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
    L0_L0bar_text_2->AddText("|#it{y}| < 1");
    L0_L0bar_text_2->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
    L0_L0bar_text_2->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
    L0_L0bar_text_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2}} = %.3f #pm %.3f #pm %.3f", P_L0_L0bar_2, fabs(P_L0_L0bar_err_2), fabs(SysErrTot_L0_L0bar*P_L0_L0bar_2)));
    //L0_L0bar_text_2->AddText(Form("P_{topo} = %.2f", P_L0_L0bar_tight_topo_cuts));
    //L0_L0bar_text_2->AddText(Form("P_{pT} = %.2f", P_L0_L0bar_tight_pT_cuts));
    L0_L0bar_text_2->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_2->Draw("same");


    TLegend *L0_L0bar_2_leg = new TLegend(0.15, 0.3, 0.4, 0.49);
    //L0_L0bar_2_leg->SetTextSizePixels(15);
    L0_L0bar_2_leg->AddEntry(L0_L0bar_cosThetaProdPlane_eta_US_hist, "(US-US)-Bckg.");
    L0_L0bar_2_leg->AddEntry(fitL0_L0bar_US_ThetaStar_2, "Fit", "l");
    L0_L0bar_2_leg->SetBorderSize(0);
    L0_L0bar_2_leg->SetFillColorAlpha(0, 0.01);
    L0_L0bar_2_leg->Draw("same");

    L0_L0bar_cosThetaProdPlane_eta_can_2->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations_delta_R/L0_L0bar_cosThetaProdPlane_eta_subtract_less_delta_R_%i.png", delta_R_bin));

    //-----------------------


    PolarizationGraph_delta_R[delta_R_bin-1]->SetPoint(1, P_L0_L0bar_2, 1);
    PolarizationGraph_delta_R[delta_R_bin-1]->SetPointError(1, fabs(P_L0_L0bar_err_2), 0);


    PolarizationGraph_delta_R_sys_err[delta_R_bin-1]->SetPoint(1, P_L0_L0bar_2, 1);
    PolarizationGraph_delta_R_sys_err[delta_R_bin-1]->SetPointError(1, fabs(SysErrTot_L0_L0bar*P_L0_L0bar_2), 0.045);

    //____________________________________________________________________________________________________________________________________________________________________________________________________________

    TCanvas *L0_L0_cosThetaProdPlane_eta_no_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_eta_no_corr_can_%i", delta_R_bin), Form("L0_L0_cosThetaProdPlane_eta_no_corr_can_%i", delta_R_bin), 1200, 1000);

    L0_L0_cosThetaProdPlane_eta_no_corr_can->cd();

    TH1D *L0_L0_cosThetaProdPlane_eta_US_hist = L0_L0_cos_theta_star_vs_delta_R_US_hist->ProjectionX( Form("proj_US_%i", delta_R_bin), delta_R_bin, delta_R_bin);

    L0_L0_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0_cosThetaProdPlane_eta_US_hist->GetXaxis()->CenterTitle();
    //L0_L0_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTextSizePixels(30);
    L0_L0_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0_cosThetaProdPlane_eta_US_hist->GetYaxis()->CenterTitle();
    //L0_L0_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTextSizePixels(30);
    L0_L0_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetMaxDigits(3);
    L0_L0_cosThetaProdPlane_eta_US_hist->SetMarkerSize(1.5);
    L0_L0_cosThetaProdPlane_eta_US_hist->SetMarkerStyle(20);
    L0_L0_cosThetaProdPlane_eta_US_hist->SetMarkerColor(kRed);
    L0_L0_cosThetaProdPlane_eta_US_hist->SetLineColor(kRed);
    double nLL = L0_L0_cosThetaProdPlane_eta_US_hist->Integral();
    L0_L0_cosThetaProdPlane_eta_US_hist->Sumw2();
    //L0_L0_cosThetaProdPlane_eta_US_hist->Divide(L0_L0_cosThetaProdPlane_eta_eff);
    L0_L0_cosThetaProdPlane_eta_US_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1));
    //L0_L0_cosThetaProdPlane_eta_US_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_US_hist->Integral());
    L0_L0_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0_L0_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //ME here just for plotting, used loser
    TH1D *L0_L0_cosThetaProdPlane_eta_ME_hist = L0_L0_cos_theta_star_vs_delta_R_US_ME_hist->ProjectionX( Form("proj_US_ME_%i", delta_R_bin), delta_R_bin, delta_R_bin);

    L0_L0_cosThetaProdPlane_eta_ME_hist->SetMarkerSize(1.5);
    L0_L0_cosThetaProdPlane_eta_ME_hist->SetMarkerStyle(24);
    L0_L0_cosThetaProdPlane_eta_ME_hist->SetMarkerColor(1);
    L0_L0_cosThetaProdPlane_eta_ME_hist->SetLineColor(1);
    L0_L0_cosThetaProdPlane_eta_ME_hist->Sumw2();
    L0_L0_cosThetaProdPlane_eta_ME_hist->Scale(nLL/L0_L0_cosThetaProdPlane_eta_ME_hist->Integral()); //scale ME to US
    L0_L0_cosThetaProdPlane_eta_ME_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_ME_hist->GetXaxis()->GetBinWidth(1)); //US is scaled by bin width, ME needs to be scaled as well
    L0_L0_cosThetaProdPlane_eta_ME_hist->SetMinimum(0);
    L0_L0_cosThetaProdPlane_eta_ME_hist->Draw("p e same");

    TF1 *fitL0_L0_US_ThetaStar_no_corr_ME = new TF1("fitL0_L0_US_ThetaStar_no_corr_ME", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0_US_ThetaStar_no_corr_ME->SetParameters(100, 0.5);

    L0_L0_cosThetaProdPlane_eta_ME_hist->Fit(fitL0_L0_US_ThetaStar_no_corr_ME, "s i 0 r");

    float P_L0_L0_no_corr_ME = fitL0_L0_US_ThetaStar_no_corr_ME->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_L0_L0_no_corr_ME_err = fitL0_L0_US_ThetaStar_no_corr_ME->GetParError(1)/(L0_alpha*L0_alpha);


    TH1D *L0_L0_cosThetaProdPlane_eta_LS_hist = L0_L0_cos_theta_star_vs_delta_R_US_LS_hist->ProjectionX( Form("proj_LS_%i", delta_R_bin), delta_R_bin, delta_R_bin);

    L0_L0_cosThetaProdPlane_eta_LS_hist->SetMarkerSize(1.5);
    L0_L0_cosThetaProdPlane_eta_LS_hist->SetMarkerStyle(21);
    L0_L0_cosThetaProdPlane_eta_LS_hist->SetMarkerColor(kBlue);
    double nLL_back = L0_L0_cosThetaProdPlane_eta_LS_hist->Integral();
    L0_L0_cosThetaProdPlane_eta_LS_hist->Sumw2();
    L0_L0_cosThetaProdPlane_eta_LS_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1));
    //L0_L0_cosThetaProdPlane_eta_LS_hist->Divide(L0_L0_cosThetaProdPlane_eta_eff);
    L0_L0_cosThetaProdPlane_eta_LS_hist->Draw("p e same");


    TH1D *L0_L0_cosThetaProdPlane_eta_ME_LS_hist = L0_L0_cos_theta_star_vs_delta_R_US_LS_ME_hist->ProjectionX( Form("proj_LS_ME_%i", delta_R_bin), delta_R_bin, delta_R_bin);

    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerSize(1.5);
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerStyle(25);
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerColor(kMagenta+1);
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->SetLineColor(kMagenta+1);
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->Sumw2();
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->Scale(nLL_back/L0_L0_cosThetaProdPlane_eta_ME_LS_hist->Integral()); //scale ME_LS to background
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_ME_LS_hist->GetXaxis()->GetBinWidth(1)); //background is scaled by bin width, ME_LS needs to be scaled as well
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->SetMinimum(0);
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->Draw("p e same");


    TF1 *fitL0_L0_US_ThetaStar_no_corr = new TF1("fitL0_L0_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0_US_ThetaStar_no_corr->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0_L0_cosThetaProdPlane_eta_US_hist->Fit(fitL0_L0_US_ThetaStar_no_corr, "s i 0 r");

    float P_L0_L0_no_corr = fitL0_L0_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_L0_L0_no_corr_err = fitL0_L0_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0_alpha);

    fitL0_L0_US_ThetaStar_no_corr->SetLineColor(1);
    //fitL0_L0_US_ThetaStar_no_corr->Draw("same");

    TLegend *L0_L0_leg = new TLegend(0.15, 0.45, 0.45, 0.69);
    L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_eta_US_hist, "(US-US) p#pi");
    L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_eta_ME_hist, "(US-US) ME");
    L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_eta_LS_hist, "Combinatorial bckg.");
    L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_eta_ME_LS_hist, "Bckg. ME");
    //L0_L0_leg->AddEntry(fitL0_L0_US_ThetaStar_no_corr, "Linear fit to US");
    L0_L0_leg->SetBorderSize(0);
    L0_L0_leg->SetFillColorAlpha(0, 0.01);
    L0_L0_leg->Draw("same");

    TPaveText *L0_L0_text_no_corr = new TPaveText(0.5, 0.4, 0.85, 0.75, "NDC");
    L0_L0_text_no_corr->SetTextFont(42);
    //L0_L0_text_no_corr->AddText("STAR Internal");
    //L0_L0_text_no_corr->AddText("STAR preliminary");
    //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
    L0_L0_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0_L0_text_no_corr->AddText("Minimum bias, no correction");
    L0_L0_text_no_corr->AddText("#Lambda^{0}-#Lambda^{0}");
    L0_L0_text_no_corr->AddText("|#it{y}| < 1");
    L0_L0_text_no_corr->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
    L0_L0_text_no_corr->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
    L0_L0_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0_no_corr, fabs(P_L0_L0_no_corr_err)));
    L0_L0_text_no_corr->AddText(Form("P_{ME} = %.3f #pm %.3f", P_L0_L0_no_corr_ME, fabs(P_L0_L0_no_corr_ME_err)));
    L0_L0_text_no_corr->SetFillColorAlpha(0, 0.01);
    L0_L0_text_no_corr->Draw("same");

    L0_L0_cosThetaProdPlane_eta_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations_delta_R/L0_L0_cosThetaProdPlane_eta_no_corr_less_delta_R_%i.png", delta_R_bin));

    //----------------------------------------------------------------------------------------------------

    TCanvas *L0_L0_cosThetaProdPlane_eta_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_eta_can_%i", delta_R_bin), Form("L0_L0_cosThetaProdPlane_eta_can_%i", delta_R_bin), 1200, 1000);

    L0_L0_cosThetaProdPlane_eta_can->cd();

    //ME histogram higher

    L0_L0_cosThetaProdPlane_eta_US_hist->Divide(L0_L0_cosThetaProdPlane_eta_ME_hist); // correct US using ME
    L0_L0_cosThetaProdPlane_eta_US_hist->Scale(nLL/L0_L0_cosThetaProdPlane_eta_US_hist->Integral()); //scale back to raw US
    L0_L0_cosThetaProdPlane_eta_US_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    L0_L0_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0_L0_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //L0_L0_cosThetaProdPlane_eta_ME_hist->Draw("same p e");

    L0_L0_cosThetaProdPlane_eta_LS_hist->Divide(L0_L0_cosThetaProdPlane_eta_ME_LS_hist); //correct background using ME
    L0_L0_cosThetaProdPlane_eta_LS_hist->Scale(nLL_back/L0_L0_cosThetaProdPlane_eta_LS_hist->Integral()); //scale back to raw background
    L0_L0_cosThetaProdPlane_eta_LS_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    L0_L0_cosThetaProdPlane_eta_LS_hist->SetMinimum(0);
    L0_L0_cosThetaProdPlane_eta_LS_hist->Draw("p e same");

    //fit dN/dcos(theta*)
    //signal + bacground
    TF1 *fitL0_L0_US_ThetaStar = new TF1("fitL0_L0_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0_US_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0_L0_cosThetaProdPlane_eta_US_hist->Fit(fitL0_L0_US_ThetaStar, "s i 0 r");

    float P_L0_L0 = fitL0_L0_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_L0_L0_err = fitL0_L0_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

    fitL0_L0_US_ThetaStar->SetLineColor(1);
    fitL0_L0_US_ThetaStar->Draw("same");

    //background
    TF1 *fitL0_L0_US_LS_ThetaStar = new TF1("fitL0_L0_US_LS_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0_US_LS_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaUS_LS_wrong_sign = L_inv_mass_US_LS->Fit(fitGaUS_LSsBack, "s i 0", "", 1.07, 1.4);
    L0_L0_cosThetaProdPlane_eta_LS_hist->Fit(fitL0_L0_US_LS_ThetaStar, "s i 0 r");

    float P_L0_L0_back = fitL0_L0_US_LS_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_L0_L0_back_err = fitL0_L0_US_LS_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

    fitL0_L0_US_LS_ThetaStar->SetLineColor(1);
    fitL0_L0_US_LS_ThetaStar->Draw("same");


    TPaveText *L0_L0_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
    L0_L0_text->SetTextFont(42);
    //L0_L0_text->AddText("STAR Internal");
    //L0_L0_text->AddText("STAR preliminary");
    //((TText*)L0_L0_text->GetListOfLines()->Last())->SetTextColor(2);
    L0_L0_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0_L0_text->AddText("Minimum bias");
    L0_L0_text->AddText("#Lambda^{0}-#Lambda^{0}");
    L0_L0_text->AddText("|#it{y}| < 1");
    L0_L0_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
    L0_L0_text->AddText(Form("P_{tot} = %.2f #pm %.2f", P_L0_L0, fabs(P_L0_L0_err)));
    L0_L0_text->AddText(Form("P_{bckg} = %.2f #pm %.2f", P_L0_L0_back, fabs(P_L0_L0_back_err)));
    L0_L0_text->SetFillColorAlpha(0, 0.01);
    L0_L0_text->Draw("same");

    L0_L0_leg->Draw("same");

    L0_L0_cosThetaProdPlane_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations_delta_R/L0_L0_cosThetaProdPlane_eta_less_delta_R_%i.png", delta_R_bin));

    //----------------------------------------------------------------------

    TCanvas *L0_L0_cosThetaProdPlane_eta_can_2 = new TCanvas(Form("L0_L0_cosThetaProdPlane_eta_can_2_%i", delta_R_bin), Form("L0_L0_cosThetaProdPlane_eta_can_2_%i", delta_R_bin), 1200, 1000);

    L0_L0_cosThetaProdPlane_eta_can_2->cd();

    L0_L0_cosThetaProdPlane_eta_US_hist->Add(L0_L0_cosThetaProdPlane_eta_LS_hist, -1);
    L0_L0_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0_L0_cosThetaProdPlane_eta_US_hist->Write(Form("L0_L0_cosThetaProdPlane_delta_R_%i", delta_R_bin-1));
    L0_L0_cosThetaProdPlane_eta_US_hist->Draw("p e");


    TF1 *fitL0_L0_US_ThetaStar_2 = new TF1("fitL0_L0_US_ThetaStar_2", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0_US_ThetaStar_2->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0_L0_cosThetaProdPlane_eta_US_hist->Fit(fitL0_L0_US_ThetaStar_2, "s i 0 r");

    //store fit result for systematic errors
    //tight topo cuts
    if(cut_type == 1)
    {
      SysErrCutsTopo->cd();

      fitL0_L0_tight_topo_cuts_delta_R[delta_R_bin-1]->SetParameters(fitL0_L0_US_ThetaStar_2->GetParameter(0), fitL0_L0_US_ThetaStar_2->GetParameter(1));
      fitL0_L0_tight_topo_cuts_delta_R[delta_R_bin-1]->SetParError(0, fitL0_L0_US_ThetaStar_2->GetParError(0));
      fitL0_L0_tight_topo_cuts_delta_R[delta_R_bin-1]->SetParError(1, fitL0_L0_US_ThetaStar_2->GetParError(1));
      fitL0_L0_tight_topo_cuts_delta_R[delta_R_bin-1]->Write();
    }

    //tight pT cuts
    if(cut_type == 2)
    {
      SysErrCutsPt->cd();

      fitL0_L0_tight_pT_cuts_delta_R[delta_R_bin-1]->SetParameters(fitL0_L0_US_ThetaStar_2->GetParameter(0), fitL0_L0_US_ThetaStar_2->GetParameter(1));
      fitL0_L0_tight_pT_cuts_delta_R[delta_R_bin-1]->SetParError(0, fitL0_L0_US_ThetaStar_2->GetParError(0));
      fitL0_L0_tight_pT_cuts_delta_R[delta_R_bin-1]->SetParError(1, fitL0_L0_US_ThetaStar_2->GetParError(1));
      fitL0_L0_tight_pT_cuts_delta_R[delta_R_bin-1]->Write();
    }

    float P_L0_L0_2 = fitL0_L0_US_ThetaStar_2->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_L0_L0_err_2 = fitL0_L0_US_ThetaStar_2->GetParError(1)/(L0_alpha*L0_alpha);

    fitL0_L0_US_ThetaStar_2->SetLineColor(1);
    fitL0_L0_US_ThetaStar_2->Draw("same");


    //TPaveText *L0_L0_text_2 = new TPaveText(0.47, 0.15, 0.85, 0.53, "NDC"); //with polarization value
    TPaveText *L0_L0_text_2 = new TPaveText(0.47, 0.15, 0.8, 0.53, "NDC"); //without polarization value
    L0_L0_text_2->SetTextFont(42);
    //L0_L0_text_2->SetTextSize(15);
    //L0_L0_text_2->AddText("STAR");
    //L0_L0_text_2->AddText("STAR preliminary");
    //((TText*)L0_L0_text_2->GetListOfLines()->Last())->SetTextColor(2);
    L0_L0_text_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0_L0_text_2->AddText("Minimum bias");
    L0_L0_text_2->AddText("#Lambda^{0}-#Lambda^{0}");
    L0_L0_text_2->AddText("|#it{y}| < 1");
    L0_L0_text_2->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
    L0_L0_text_2->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
    L0_L0_text_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2}} = %.3f #pm %.3f", P_L0_L0_2, fabs(P_L0_L0_err_2)));
    //L0_L0_text_2->AddText(Form("P_{topo} = %.2f", P_L0_L0_tight_topo_cuts));
    //L0_L0_text_2->AddText(Form("P_{pT} = %.2f", P_L0_L0_tight_pT_cuts));
    L0_L0_text_2->SetFillColorAlpha(0, 0.01);
    L0_L0_text_2->Draw("same");


    TLegend *L0_L0_2_leg = new TLegend(0.15, 0.3, 0.4, 0.49);
    //L0_L0_2_leg->SetTextSizePixels(15);
    L0_L0_2_leg->AddEntry(L0_L0_cosThetaProdPlane_eta_US_hist, "(US-US)-Bckg.");
    L0_L0_2_leg->AddEntry(fitL0_L0_US_ThetaStar_2, "Fit", "l");
    L0_L0_2_leg->SetBorderSize(0);
    L0_L0_2_leg->SetFillColorAlpha(0, 0.01);
    L0_L0_2_leg->Draw("same");

    L0_L0_cosThetaProdPlane_eta_can_2->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations_delta_R/L0_L0_cosThetaProdPlane_eta_subtract_less_delta_R_%i.png", delta_R_bin));

    PolarizationGraph_delta_R[delta_R_bin-1]->SetPoint(2, P_L0_L0_2, 2);
    PolarizationGraph_delta_R[delta_R_bin-1]->SetPointError(2, fabs(P_L0_L0_err_2), 0);

    //calculate total systematic uncertainty for L-Lbar

    //slope
    float SysErrSlope_L0_L0 = SysErrSlope_delta_R_hist[delta_R_bin-1]->GetBinContent(2);

    //----------------------------

    //background subtraction
    float P_L0_L0_from_fits = P_L0_L0 - nLL_back/nLL*P_L0_L0_back;
    float P_L0_L0_from_fits_err = sqrt( P_L0_L0_err*P_L0_L0_err + nLL_back/nLL*nLL_back/nLL*P_L0_L0_back_err*P_L0_L0_back_err );

    //statistical error correction for systematic error
    float SysErrBackground_L0_L0_corr = sqrt( fabs( P_L0_L0_err_2*P_L0_L0_err_2 - P_L0_L0_from_fits_err*P_L0_L0_from_fits_err ) );

    float SysErrBackground_L0_L0_work = ( fabs( P_L0_L0_2 - P_L0_L0_from_fits) - SysErrBackground_L0_L0_corr )/fabs(P_L0_L0_2);

    float SysErrBackground_L0_L0 = 0;

    if( SysErrBackground_L0_L0_work > 0 ) SysErrBackground_L0_L0 = SysErrBackground_L0_L0_work; //store sys. err. only if it is larger than statistical fluctuations

    //----------------------------

    //cuts variation
    float SysErr_tight_cuts_L0_L0 = 0;

    //float P_L0_L0_tight_topo_cuts = 0;
    //float P_L0_L0_tight_pT_cuts = 0;

    if(cut_type == 0)
    {
      float P_L0_L0_tight_topo_cuts = fitL0_L0_tight_topo_cuts_delta_R[delta_R_bin-1]->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_L0_L0_tight_topo_cuts_err = fitL0_L0_tight_topo_cuts_delta_R[delta_R_bin-1]->GetParError(1)/(L0_alpha*L0_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_topo_cuts_L0_L0_corr = sqrt( fabs( P_L0_L0_err_2*P_L0_L0_err_2 - P_L0_L0_tight_topo_cuts_err*P_L0_L0_tight_topo_cuts_err ) );

      float SysErr_tight_topo_cuts_L0_L0_work = ( fabs( P_L0_L0_2 - P_L0_L0_tight_topo_cuts) - SysErr_tight_topo_cuts_L0_L0_corr )/fabs(P_L0_L0_2);

      float SysErr_tight_topo_cuts_L0_L0 = 0;

      if( SysErr_tight_topo_cuts_L0_L0_work > 0 ) SysErr_tight_topo_cuts_L0_L0 = SysErr_tight_topo_cuts_L0_L0_work; //store sys. err. only if it is larger than statistical fluctuations


      float P_L0_L0_tight_pT_cuts = fitL0_L0_tight_pT_cuts_delta_R[delta_R_bin-1]->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_L0_L0_tight_pT_cuts_err = fitL0_L0_tight_pT_cuts_delta_R[delta_R_bin-1]->GetParError(1)/(L0_alpha*L0_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_pT_cuts_L0_L0_corr = sqrt( fabs( P_L0_L0_err_2*P_L0_L0_err_2 - P_L0_L0_tight_pT_cuts_err*P_L0_L0_tight_pT_cuts_err ) );

      float SysErr_tight_pT_cuts_L0_L0_work = ( fabs( P_L0_L0_2 - P_L0_L0_tight_pT_cuts) - SysErr_tight_pT_cuts_L0_L0_corr )/fabs(P_L0_L0_2);

      float SysErr_tight_pT_cuts_L0_L0 = 0;

      if( SysErr_tight_pT_cuts_L0_L0_work > 0 ) SysErr_tight_pT_cuts_L0_L0 = SysErr_tight_pT_cuts_L0_L0_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_cuts_L0_L0 = sqrt(SysErr_tight_topo_cuts_L0_L0*SysErr_tight_topo_cuts_L0_L0 + SysErr_tight_pT_cuts_L0_L0*SysErr_tight_pT_cuts_L0_L0);
      //SysErr_tight_cuts_L0_L0 = SysErr_tight_topo_cuts_L0_L0;
      //SysErr_tight_cuts_L0_L0 = SysErr_tight_pT_cuts_L0_L0;
    }

    //total
    //float SysErrTot_L0_L0 = sqrt( sysErr_alpha_L0_L0*sysErr_alpha_L0_L0 + SysErrSlope_L0_L0*SysErrSlope_L0_L0 + SysErrBackground_L0_L0*SysErrBackground_L0_L0  );
    float SysErrTot_L0_L0 = sqrt( sysErr_alpha_L0_L0*sysErr_alpha_L0_L0 + SysErrSlope_L0_L0*SysErrSlope_L0_L0 + SysErrBackground_L0_L0*SysErrBackground_L0_L0 + SysErr_tight_cuts_L0_L0*SysErr_tight_cuts_L0_L0  );

    PolarizationGraph_delta_R_sys_err[delta_R_bin-1]->SetPoint(2, P_L0_L0_2, 2);
    PolarizationGraph_delta_R_sys_err[delta_R_bin-1]->SetPointError(2, fabs(SysErrTot_L0_L0*P_L0_L0_2), 0.045);

    //____________________________________________________________________________________________________________________________________________________________________________________________________________

    TCanvas *L0bar_L0bar_cosThetaProdPlane_eta_no_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_eta_no_corr_can_%i", delta_R_bin), Form("L0bar_L0bar_cosThetaProdPlane_eta_no_corr_can_%i", delta_R_bin), 1200, 1000);

    L0bar_L0bar_cosThetaProdPlane_eta_no_corr_can->cd();

    TH1D *L0bar_L0bar_cosThetaProdPlane_eta_US_hist = L0bar_L0bar_cos_theta_star_vs_delta_R_US_hist->ProjectionX( Form("proj_US_%i", delta_R_bin), delta_R_bin, delta_R_bin);

    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->CenterTitle();
    //L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTextSizePixels(30);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->CenterTitle();
    //L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTextSizePixels(30);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetMaxDigits(3);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetMarkerSize(1.5);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetMarkerStyle(20);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetMarkerColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetLineColor(kRed);
    double nLbarLbar = L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Integral();
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Sumw2();
    //L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Divide(L0bar_L0bar_cosThetaProdPlane_eta_eff);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1));
    //L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Integral());
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //ME here just for plotting, used loser
    TH1D *L0bar_L0bar_cosThetaProdPlane_eta_ME_hist = L0bar_L0bar_cos_theta_star_vs_delta_R_US_ME_hist->ProjectionX( Form("proj_US_ME_%i", delta_R_bin), delta_R_bin, delta_R_bin);

    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerSize(1.5);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerStyle(24);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerColor(1);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->SetLineColor(1);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Scale(nLbarLbar/L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Integral()); //scale ME to US
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->GetXaxis()->GetBinWidth(1)); //US is scaled by bin width, ME needs to be scaled as well
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Draw("p e same");

    TF1 *fitL0bar_L0bar_US_ThetaStar_no_corr_ME = new TF1("fitL0bar_L0bar_US_ThetaStar_no_corr_ME", "[0]*(1 + [1]*x)", -1, 1);
    fitL0bar_L0bar_US_ThetaStar_no_corr_ME->SetParameters(100, 0.5);

    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Fit(fitL0bar_L0bar_US_ThetaStar_no_corr_ME, "s i 0 r");

    float P_L0bar_L0bar_no_corr_ME = fitL0bar_L0bar_US_ThetaStar_no_corr_ME->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float P_L0bar_L0bar_no_corr_ME_err = fitL0bar_L0bar_US_ThetaStar_no_corr_ME->GetParError(1)/(L0bar_alpha*L0bar_alpha);


    TH1D *L0bar_L0bar_cosThetaProdPlane_eta_LS_hist = L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_hist->ProjectionX( Form("proj_LS_%i", delta_R_bin), delta_R_bin, delta_R_bin);

    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerSize(1.5);
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerStyle(21);
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerColor(kBlue);
    double nLbarLbar_back = L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Integral();
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1));
    //L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Divide(L0bar_L0bar_cosThetaProdPlane_eta_eff);
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Draw("p e same");


    TH1D *L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist = L0bar_L0bar_cos_theta_star_vs_delta_R_US_LS_ME_hist->ProjectionX( Form("proj_LS_ME_%i", delta_R_bin), delta_R_bin, delta_R_bin);

    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerSize(1.5);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerStyle(25);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerColor(kMagenta+1);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetLineColor(kMagenta+1);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Scale(nLbarLbar_back/L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Integral()); //scale ME_LS to background
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->GetXaxis()->GetBinWidth(1)); //background is scaled by bin width, ME_LS needs to be scaled as well
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Draw("p e same");


    TF1 *fitL0bar_L0bar_US_ThetaStar_no_corr = new TF1("fitL0bar_L0bar_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
    fitL0bar_L0bar_US_ThetaStar_no_corr->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Fit(fitL0bar_L0bar_US_ThetaStar_no_corr, "s i 0 r");

    float P_L0bar_L0bar_no_corr = fitL0bar_L0bar_US_ThetaStar_no_corr->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float P_L0bar_L0bar_no_corr_err = fitL0bar_L0bar_US_ThetaStar_no_corr->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    fitL0bar_L0bar_US_ThetaStar_no_corr->SetLineColor(1);
    //fitL0bar_L0bar_US_ThetaStar_no_corr->Draw("same");

    TLegend *L0bar_L0bar_leg = new TLegend(0.15, 0.45, 0.45, 0.69);
    L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_eta_US_hist, "(US-US) p#pi");
    L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_eta_ME_hist, "(US-US) ME");
    L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_eta_LS_hist, "Combinatorial bckg.");
    L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist, "Bckg. ME");
    //L0bar_L0bar_leg->AddEntry(fitL0bar_L0bar_US_ThetaStar_no_corr, "Linear fit to US");
    L0bar_L0bar_leg->SetBorderSize(0);
    L0bar_L0bar_leg->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_leg->Draw("same");

    TPaveText *L0bar_L0bar_text_no_corr = new TPaveText(0.5, 0.4, 0.85, 0.75, "NDC");
    L0bar_L0bar_text_no_corr->SetTextFont(42);
    //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
    //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
    //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
    L0bar_L0bar_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0bar_L0bar_text_no_corr->AddText("Minimum bias, no correction");
    L0bar_L0bar_text_no_corr->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
    L0bar_L0bar_text_no_corr->AddText("|#it{y}| < 1");
    L0bar_L0bar_text_no_corr->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
    L0bar_L0bar_text_no_corr->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
    L0bar_L0bar_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar_no_corr, fabs(P_L0bar_L0bar_no_corr_err)));
    L0bar_L0bar_text_no_corr->AddText(Form("P_{ME} = %.3f #pm %.3f", P_L0bar_L0bar_no_corr_ME, fabs(P_L0bar_L0bar_no_corr_ME_err)));
    L0bar_L0bar_text_no_corr->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_no_corr->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_eta_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations_delta_R/L0bar_L0bar_cosThetaProdPlane_eta_no_corr_less_delta_R_%i.png", delta_R_bin));

    //----------------------------------------------------------------------------------------------------

    TCanvas *L0bar_L0bar_cosThetaProdPlane_eta_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_eta_can_%i", delta_R_bin), Form("L0bar_L0bar_cosThetaProdPlane_eta_can_%i", delta_R_bin), 1200, 1000);

    L0bar_L0bar_cosThetaProdPlane_eta_can->cd();

    //ME histogram higher

    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Divide(L0bar_L0bar_cosThetaProdPlane_eta_ME_hist); // correct US using ME
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Scale(nLbarLbar/L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Integral()); //scale back to raw US
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Draw("same p e");

    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Divide(L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist); //correct background using ME
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Scale(nLbarLbar_back/L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Integral()); //scale back to raw background
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Draw("p e same");

    //fit dN/dcos(theta*)
    //signal + bacground
    TF1 *fitL0bar_L0bar_US_ThetaStar = new TF1("fitL0bar_L0bar_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitL0bar_L0bar_US_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Fit(fitL0bar_L0bar_US_ThetaStar, "s i 0 r");

    float P_L0bar_L0bar = fitL0bar_L0bar_US_ThetaStar->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float P_L0bar_L0bar_err = fitL0bar_L0bar_US_ThetaStar->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    fitL0bar_L0bar_US_ThetaStar->SetLineColor(1);
    fitL0bar_L0bar_US_ThetaStar->Draw("same");

    //background
    TF1 *fitL0bar_L0bar_US_LS_ThetaStar = new TF1("fitL0bar_L0bar_US_LS_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitL0bar_L0bar_US_LS_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaUS_LS_wrong_sign = L_inv_mass_US_LS->Fit(fitGaUS_LSsBack, "s i 0", "", 1.07, 1.4);
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Fit(fitL0bar_L0bar_US_LS_ThetaStar, "s i 0 r");

    float P_L0bar_L0bar_back = fitL0bar_L0bar_US_LS_ThetaStar->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float P_L0bar_L0bar_back_err = fitL0bar_L0bar_US_LS_ThetaStar->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    fitL0bar_L0bar_US_LS_ThetaStar->SetLineColor(1);
    fitL0bar_L0bar_US_LS_ThetaStar->Draw("same");


    TPaveText *L0bar_L0bar_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
    L0bar_L0bar_text->SetTextFont(42);
    //L0bar_L0bar_text->AddText("STAR Internal");
    //L0bar_L0bar_text->AddText("STAR preliminary");
    //((TText*)L0bar_L0bar_text->GetListOfLines()->Last())->SetTextColor(2);
    L0bar_L0bar_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0bar_L0bar_text->AddText("Minimum bias");
    L0bar_L0bar_text->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
    L0bar_L0bar_text->AddText("|#it{y}| < 1");
    L0bar_L0bar_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
    L0bar_L0bar_text->AddText(Form("P_{tot} = %.2f #pm %.2f", P_L0bar_L0bar, fabs(P_L0bar_L0bar_err)));
    L0bar_L0bar_text->AddText(Form("P_{bckg} = %.2f #pm %.2f", P_L0bar_L0bar_back, fabs(P_L0bar_L0bar_back_err)));
    L0bar_L0bar_text->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text->Draw("same");

    L0bar_L0bar_leg->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations_delta_R/L0bar_L0bar_cosThetaProdPlane_eta_less_delta_R_%i.png", delta_R_bin));

    //----------------------------------------------------------------------

    TCanvas *L0bar_L0bar_cosThetaProdPlane_eta_can_2 = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_eta_can_2_%i", delta_R_bin), Form("L0bar_L0bar_cosThetaProdPlane_eta_can_2_%i", delta_R_bin), 1200, 1000);

    L0bar_L0bar_cosThetaProdPlane_eta_can_2->cd();

    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Add(L0bar_L0bar_cosThetaProdPlane_eta_LS_hist, -1);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Write(Form("L0bar_L0bar_cosThetaProdPlane_delta_R_%i", delta_R_bin-1));
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");


    TF1 *fitL0bar_L0bar_US_ThetaStar_2 = new TF1("fitL0bar_L0bar_US_ThetaStar_2", "[0]*(1 + [1]*x)", -1, 1);
    fitL0bar_L0bar_US_ThetaStar_2->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Fit(fitL0bar_L0bar_US_ThetaStar_2, "s i 0 r");

    //store fit result for systematic errors
    //tight topo cuts
    if(cut_type == 1)
    {
      SysErrCutsTopo->cd();

      fitL0bar_L0bar_tight_topo_cuts_delta_R[delta_R_bin-1]->SetParameters(fitL0bar_L0bar_US_ThetaStar_2->GetParameter(0), fitL0bar_L0bar_US_ThetaStar_2->GetParameter(1));
      fitL0bar_L0bar_tight_topo_cuts_delta_R[delta_R_bin-1]->SetParError(0, fitL0bar_L0bar_US_ThetaStar_2->GetParError(0));
      fitL0bar_L0bar_tight_topo_cuts_delta_R[delta_R_bin-1]->SetParError(1, fitL0bar_L0bar_US_ThetaStar_2->GetParError(1));
      fitL0bar_L0bar_tight_topo_cuts_delta_R[delta_R_bin-1]->Write();
    }

    //tight pT cuts
    if(cut_type == 2)
    {
      SysErrCutsPt->cd();

      fitL0bar_L0bar_tight_pT_cuts_delta_R[delta_R_bin-1]->SetParameters(fitL0bar_L0bar_US_ThetaStar_2->GetParameter(0), fitL0bar_L0bar_US_ThetaStar_2->GetParameter(1));
      fitL0bar_L0bar_tight_pT_cuts_delta_R[delta_R_bin-1]->SetParError(0, fitL0bar_L0bar_US_ThetaStar_2->GetParError(0));
      fitL0bar_L0bar_tight_pT_cuts_delta_R[delta_R_bin-1]->SetParError(1, fitL0bar_L0bar_US_ThetaStar_2->GetParError(1));
      fitL0bar_L0bar_tight_pT_cuts_delta_R[delta_R_bin-1]->Write();
    }

    float P_L0bar_L0bar_2 = fitL0bar_L0bar_US_ThetaStar_2->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float P_L0bar_L0bar_err_2 = fitL0bar_L0bar_US_ThetaStar_2->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    fitL0bar_L0bar_US_ThetaStar_2->SetLineColor(1);
    fitL0bar_L0bar_US_ThetaStar_2->Draw("same");


    //TPaveText *L0bar_L0bar_text_2 = new TPaveText(0.47, 0.15, 0.85, 0.53, "NDC"); //with polarization value
    TPaveText *L0bar_L0bar_text_2 = new TPaveText(0.47, 0.15, 0.8, 0.53, "NDC"); //without polarization value
    L0bar_L0bar_text_2->SetTextFont(42);
    //L0bar_L0bar_text_2->SetTextSize(15);
    //L0bar_L0bar_text_2->AddText("STAR");
    //L0bar_L0bar_text_2->AddText("STAR preliminary");
    //((TText*)L0bar_L0bar_text_2->GetListOfLines()->Last())->SetTextColor(2);
    L0bar_L0bar_text_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0bar_L0bar_text_2->AddText("Minimum bias");
    L0bar_L0bar_text_2->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
    L0bar_L0bar_text_2->AddText("|#it{y}| < 1");
    L0bar_L0bar_text_2->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
    L0bar_L0bar_text_2->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
    L0bar_L0bar_text_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2}} = %.3f #pm %.3f", P_L0bar_L0bar_2, fabs(P_L0bar_L0bar_err_2)));
    //L0bar_L0bar_text_2->AddText(Form("P_{topo} = %.2f", P_L0bar_L0bar_tight_topo_cuts));
    //L0bar_L0bar_text_2->AddText(Form("P_{pT} = %.2f", P_L0bar_L0bar_tight_pT_cuts));
    L0bar_L0bar_text_2->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_2->Draw("same");


    TLegend *L0bar_L0bar_2_leg = new TLegend(0.15, 0.3, 0.4, 0.49);
    //L0bar_L0bar_2_leg->SetTextSizePixels(15);
    L0bar_L0bar_2_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_eta_US_hist, "(US-US)-Bckg.");
    L0bar_L0bar_2_leg->AddEntry(fitL0bar_L0bar_US_ThetaStar_2, "Fit", "l");
    L0bar_L0bar_2_leg->SetBorderSize(0);
    L0bar_L0bar_2_leg->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_2_leg->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_eta_can_2->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations_delta_R/L0bar_L0bar_cosThetaProdPlane_eta_subtract_less_delta_R_%i.png", delta_R_bin));

    PolarizationGraph_delta_R[delta_R_bin-1]->SetPoint(3, P_L0bar_L0bar_2, 3);
    PolarizationGraph_delta_R[delta_R_bin-1]->SetPointError(3, fabs(P_L0bar_L0bar_err_2), 0);

    //calculate total systematic uncertainty for L-Lbar

    //slope
    float SysErrSlope_L0bar_L0bar = SysErrSlope_delta_R_hist[delta_R_bin-1]->GetBinContent(3);

    //----------------------------

    //background subtraction
    float P_L0bar_L0bar_from_fits = P_L0bar_L0bar - nLbarLbar_back/nLbarLbar*P_L0bar_L0bar_back;
    float P_L0bar_L0bar_from_fits_err = sqrt( P_L0bar_L0bar_err*P_L0bar_L0bar_err + nLbarLbar_back/nLbarLbar*nLbarLbar_back/nLbarLbar*P_L0bar_L0bar_back_err*P_L0bar_L0bar_back_err );

    //statistical error correction for systematic error
    float SysErrBackground_L0bar_L0bar_corr = sqrt( fabs( P_L0bar_L0bar_err_2*P_L0bar_L0bar_err_2 - P_L0bar_L0bar_from_fits_err*P_L0bar_L0bar_from_fits_err ) );

    float SysErrBackground_L0bar_L0bar_work = ( fabs( P_L0bar_L0bar_2 - P_L0bar_L0bar_from_fits) - SysErrBackground_L0bar_L0bar_corr )/fabs(P_L0bar_L0bar_2);

    float SysErrBackground_L0bar_L0bar = 0;

    if( SysErrBackground_L0bar_L0bar_work > 0 ) SysErrBackground_L0bar_L0bar = SysErrBackground_L0bar_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations

    //----------------------------

    //cuts variation
    float SysErr_tight_cuts_L0bar_L0bar = 0;

    //float P_L0bar_L0bar_tight_topo_cuts = 0;
    //float P_L0bar_L0bar_tight_pT_cuts = 0;

    if(cut_type == 0)
    {
      float P_L0bar_L0bar_tight_topo_cuts = fitL0bar_L0bar_tight_topo_cuts_delta_R[delta_R_bin-1]->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
      float P_L0bar_L0bar_tight_topo_cuts_err = fitL0bar_L0bar_tight_topo_cuts_delta_R[delta_R_bin-1]->GetParError(1)/(L0bar_alpha*L0bar_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_topo_cuts_L0bar_L0barbar_corr = sqrt( fabs( P_L0bar_L0bar_err_2*P_L0bar_L0bar_err_2 - P_L0bar_L0bar_tight_topo_cuts_err*P_L0bar_L0bar_tight_topo_cuts_err ) );

      float SysErr_tight_topo_cuts_L0bar_L0barbar_work = ( fabs( P_L0bar_L0bar_2 - P_L0bar_L0bar_tight_topo_cuts) - SysErr_tight_topo_cuts_L0bar_L0barbar_corr )/fabs(P_L0bar_L0bar_2);

      float SysErr_tight_topo_cuts_L0bar_L0bar = 0;

      if( SysErr_tight_topo_cuts_L0bar_L0barbar_work > 0 ) SysErr_tight_topo_cuts_L0bar_L0bar = SysErr_tight_topo_cuts_L0bar_L0barbar_work; //store sys. err. only if it is larger than statistical fluctuations

      //SysErr_tight_topo_cuts_L0bar_L0bar = fabs( P_L0bar_L0bar_2 - P_L0bar_L0bar_tight_topo_cuts)/fabs(P_L0bar_L0bar_2);

      float P_L0bar_L0bar_tight_pT_cuts = fitL0bar_L0bar_tight_pT_cuts_delta_R[delta_R_bin-1]->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
      float P_L0bar_L0bar_tight_pT_cuts_err = fitL0bar_L0bar_tight_pT_cuts_delta_R[delta_R_bin-1]->GetParError(1)/(L0bar_alpha*L0bar_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_pT_cuts_L0bar_L0barbar_corr = sqrt( fabs( P_L0bar_L0bar_err_2*P_L0bar_L0bar_err_2 - P_L0bar_L0bar_tight_pT_cuts_err*P_L0bar_L0bar_tight_pT_cuts_err ) );

      float SysErr_tight_pT_cuts_L0bar_L0barbar_work = ( fabs( P_L0bar_L0bar_2 - P_L0bar_L0bar_tight_pT_cuts) - SysErr_tight_pT_cuts_L0bar_L0barbar_corr )/fabs(P_L0bar_L0bar_2);

      float SysErr_tight_pT_cuts_L0bar_L0bar = 0;

      if( SysErr_tight_pT_cuts_L0bar_L0barbar_work > 0 ) SysErr_tight_pT_cuts_L0bar_L0bar = SysErr_tight_pT_cuts_L0bar_L0barbar_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_cuts_L0bar_L0bar = sqrt(SysErr_tight_topo_cuts_L0bar_L0bar*SysErr_tight_topo_cuts_L0bar_L0bar + SysErr_tight_pT_cuts_L0bar_L0bar*SysErr_tight_pT_cuts_L0bar_L0bar);
      //SysErr_tight_cuts_L0bar_L0bar = SysErr_tight_topo_cuts_L0bar_L0bar;
      //SysErr_tight_cuts_L0bar_L0bar = SysErr_tight_pT_cuts_L0bar_L0bar;
    }

    //total
    //float SysErrTot_L0bar_L0bar = sqrt( sysErr_alpha_L0bar_L0bar*sysErr_alpha_L0bar_L0bar + SysErrSlope_L0bar_L0bar*SysErrSlope_L0bar_L0bar + SysErrBackground_L0bar_L0bar*SysErrBackground_L0bar_L0bar  );
    float SysErrTot_L0bar_L0bar = sqrt( sysErr_alpha_L0bar_L0bar*sysErr_alpha_L0bar_L0bar + SysErrSlope_L0bar_L0bar*SysErrSlope_L0bar_L0bar + SysErrBackground_L0bar_L0bar*SysErrBackground_L0bar_L0bar + SysErr_tight_cuts_L0bar_L0bar*SysErr_tight_cuts_L0bar_L0bar );

    PolarizationGraph_delta_R_sys_err[delta_R_bin-1]->SetPoint(3, P_L0bar_L0bar_2, 3);
    PolarizationGraph_delta_R_sys_err[delta_R_bin-1]->SetPointError(3, fabs(SysErrTot_L0bar_L0bar*P_L0bar_L0bar_2), 0.045);

    //____________________________________________________________________________________________________________________________________________________________________________________________________________

    //save polarization graphs
    out_file->cd();

    PolarizationGraph_delta_R[delta_R_bin-1]->Write(Form("PolarizationGraph_delta_R_%i", delta_R_bin-1));
    PolarizationGraph_delta_R_sys_err[delta_R_bin-1]->Write(Form("PolarizationGraph_delta_R_sys_err_%i", delta_R_bin-1));

    //____________________________________________________________________________________________________________________________________________________________________________________________________________

    //plot polarization graphs in bins
    PolarizationGraph_delta_R_can->cd(delta_R_bin);

    //gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.05);

    DefaultHist->Draw();

    PolarizationGraph_delta_R[delta_R_bin-1]->SetMarkerStyle(20);
    PolarizationGraph_delta_R[delta_R_bin-1]->SetMarkerSize(2);
    PolarizationGraph_delta_R[delta_R_bin-1]->SetMarkerColor(kRed);
    PolarizationGraph_delta_R[delta_R_bin-1]->SetLineColor(kRed);
    PolarizationGraph_delta_R[delta_R_bin-1]->Draw("p e same");

    PolarizationGraph_delta_R_sys_err[delta_R_bin-1]->SetMarkerSize(2);
    PolarizationGraph_delta_R_sys_err[delta_R_bin-1]->SetFillColorAlpha(kRed, 0.25);
    PolarizationGraph_delta_R_sys_err[delta_R_bin-1]->SetLineColor(kRed);
    PolarizationGraph_delta_R_sys_err[delta_R_bin-1]->Draw("2 same");

    ZeroLine_eta->Draw("same");

    bin_text[delta_R_bin-1] = new TPaveText(0.12, 0.85, 0.49, 0.89, "NDC");
    bin_text[delta_R_bin-1]->SetTextFont(42);
    if(delta_R_bin == 1) bin_text[delta_R_bin-1]->AddText("In-cone");
    if(delta_R_bin == 2) bin_text[delta_R_bin-1]->AddText("Out-of-cone");
    bin_text[delta_R_bin-1]->SetFillColorAlpha(0, 0.01);
    bin_text[delta_R_bin-1]->Draw("same");

    if(delta_R_bin == 1 ) Polarization_text->Draw("same");

    //----------------------------------L pair stats---------------------------------
    cout<<endl;
    cout<<"Bin "<<delta_R_bin<<endl;
    cout<<endl;
    cout<<"N L-Lbar pairs from hist: "<<nLLbar<<endl;
    cout<<"N L-Lbar background pairs from hist: "<<nLLbar_back<<endl;
    cout<<endl;
    cout<<"N L-L pairs from hist: "<<nLL<<endl;
    cout<<"N L-L background pairs from hist: "<<nLL_back<<endl;
    cout<<endl;
    cout<<"N Lbar-Lbar pairs from hist: "<<nLbarLbar<<endl;
    cout<<"N Lbar-Lbar background pairs from hist: "<<nLbarLbar_back<<endl;
    cout<<endl;

  }


  PolarizationGraph_delta_R_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_polarization/L_polarization_delta_R.png");

  //__________________________________________________________________________________________________


  LLbarOutFile->Close();

  return true;
}
