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

bool Ana004_K0s_corr_2D_get_corr_Delta_phi_new_old(const int cut_type = 0, const int energy = 510, const int year = 2017)
{
  //analyze stored Lambda pairs and save cos(theta*) histograms

  const float L0_alpha = 0.732; //decay parameter of L0
  const float L0_alpha_relat_err = 0.014/L0_alpha; //relative error of decay parameter

  const float L0bar_alpha = -0.758; //decay paramteter of L0bar
  const float L0bar_alpha_relat_err = 0.012/fabs(L0bar_alpha); //relative error of decay paramteter


  //_______________________________________________________________________________________________________________________________________________

  //systematic uncertainties
  //residual effect from PYTHIA closure test
  //residual effect from PYTHIA closure test
  TFile *SysErrSlopeFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErrSlope.root", year), "read");

  TH1F *SysErrSlope_delta_eta_hist[2];

  for( unsigned int delta_eta_bin = 0; delta_eta_bin < 2; delta_eta_bin++ )
  {
    SysErrSlope_delta_eta_hist[delta_eta_bin] = (TH1F*)SysErrSlopeFile->Get(Form("SysErrSlope_delta_eta_hist_%i", delta_eta_bin));
  }

/*
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
    //SysErrCutsTopo = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_topo_cuts_Delta_phi.root", year), "read");
    SysErrCutsTopo = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_topo_cuts_Delta_phi_work.root", year), "read");

    //SysErrCutsPt = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_pT_cuts_Delta_phi.root", year), "read");
    SysErrCutsPt = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_pT_cuts_Delta_phi_work.root", year), "read");
  }
  else
  {
    if( cut_type == 1 ) SysErrCutsTopo = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_topo_cuts_Delta_phi_work.root", year), "recreate");
    if( cut_type == 2 ) SysErrCutsPt = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_pT_cuts_Delta_phi_work.root", year), "recreate");
  }

  //tight cuts
  TF1 *fitL0_L0_tight_topo_cuts_Delta_phi[2];

  TF1 *fitL0_L0_tight_pT_cuts_Delta_phi[2];

  for( unsigned int delta_eta_bin = 0; delta_eta_bin < 2; delta_eta_bin++)
  {
    if(cut_type == 0 )
    {
      fitL0_L0_tight_topo_cuts_Delta_phi[delta_eta_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0_L0_tight_topo_cuts_Delta_phi_%i", delta_eta_bin));

      fitL0_L0_tight_pT_cuts_Delta_phi[delta_eta_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0_L0_tight_pT_cuts_Delta_phi_%i", delta_eta_bin));
    }
    else
    {
      fitL0_L0bar_tight_topo_cuts_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0_L0bar_tight_topo_cuts_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_tight_topo_cuts_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0_L0_tight_topo_cuts_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0_L0_tight_topo_cuts_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_tight_topo_cuts_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0bar_L0bar_tight_topo_cuts_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0bar_L0bar_tight_topo_cuts_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_tight_topo_cuts_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);


      fitL0_L0bar_tight_pT_cuts_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0_L0bar_tight_pT_cuts_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_tight_pT_cuts_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0_L0_tight_pT_cuts_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0_L0_tight_pT_cuts_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_tight_pT_cuts_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0bar_L0bar_tight_pT_cuts_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0bar_L0bar_tight_pT_cuts_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_tight_pT_cuts_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);
    }

  }
*/
  //_________________________________________________________

    //output file with polarization graphs
  TFile *out_file = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/Polarization/%i/Polarization_K0s_Delta_phi.root", year), "recreate");

  //systematic uncertainty histograms and values
  //alpha
  float sysErr_alpha_K0s_K0sbar = sqrt(L0_alpha_relat_err*L0_alpha_relat_err + L0bar_alpha_relat_err*L0bar_alpha_relat_err);
  float sysErr_alpha_K0s_K0s = sqrt(L0_alpha_relat_err*L0_alpha_relat_err + L0_alpha_relat_err*L0_alpha_relat_err);
  float sysErr_alpha_K0sbar_K0sbar = sqrt(L0bar_alpha_relat_err*L0bar_alpha_relat_err + L0bar_alpha_relat_err*L0bar_alpha_relat_err);

  //----------------------------------------------------------------------

  TFile *inFile; //output file to store production plane histograms;

  if(cut_type == 0)
  {
    if(year == 2012)
    {
      inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_ana_cuts_work.root", year), "read");
      //inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_ana_cuts_orig_full_prod.root", year), "read");
      //inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_ana_cuts_new_prod.root", year), "read");
    }
    else
    {
      inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_ana_cuts.root", year), "read");
    }

    //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  }
  else if(cut_type == 1) //create production plane file from nTuple - run in this mode first
  {
    inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_tight_topo_cuts.root", year), "read");
  }
  else if(cut_type == 2) //create production plane file from nTuple - run in this mode first
  {
    inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_tight_pT_cut.root", year), "read");
  }
  else  //read file to store wrong and good sign M_inv distributions and daughter pT
  {
    cout<<"Wrong cut type!"<<endl;
    return false;
  }


  //data histograms

  TH2F *K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_hist");

  TH2F *K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist");

  //--------------------------------------------------------

  //mixed event

  TH2F *K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist");

  TH2F *K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist");


  //________________________________________________________________________________________



  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(42);

  //polarization graph for Delta y
  //first graph is for |Delta phi| < pi/4, second is for |Delta phi| > pi/4
  TGraphErrors *PolarizationGraph_delta_phi[2];
  TGraphErrors *PolarizationGraph_delta_phi_sys_err[2];

  for( unsigned int delta_eta_bin = 1; delta_eta_bin < 3; delta_eta_bin++)
  {
    //create polarization graphs, defined earlier

    PolarizationGraph_delta_phi[delta_eta_bin-1] = new TGraphErrors(1);
    PolarizationGraph_delta_phi_sys_err[delta_eta_bin-1] = new TGraphErrors(1);
  }



  out_file->cd();

  //----------------------------------------

  TCanvas *K0s_K0s_cosThetaProdPlane_no_corr_can = new TCanvas("K0s_K0s_cosThetaProdPlane_no_corr_can", "K0s_K0s_cosThetaProdPlane_no_corr_can", 1200, 1000);

  K0s_K0s_cosThetaProdPlane_no_corr_can->cd();

  TH1D *K0s_K0s_cosThetaProdPlane_US_hist = K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_hist->ProjectionX( "proj_LL_US", 1, 5); // |Delta_phi| < pi/4

  K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->CenterTitle();
  //K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->SetTextSizePixels(30);
  K0s_K0s_cosThetaProdPlane_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_US_hist->GetYaxis()->CenterTitle();
  //K0s_K0s_cosThetaProdPlane_US_hist->GetYaxis()->SetTextSizePixels(30);
  K0s_K0s_cosThetaProdPlane_US_hist->GetYaxis()->SetMaxDigits(3);
  K0s_K0s_cosThetaProdPlane_US_hist->SetMarkerSize(1.5);
  K0s_K0s_cosThetaProdPlane_US_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_US_hist->SetMarkerColor(kRed);
  K0s_K0s_cosThetaProdPlane_US_hist->SetLineColor(kRed);
  double nLL = K0s_K0s_cosThetaProdPlane_US_hist->Integral();
  K0s_K0s_cosThetaProdPlane_US_hist->Sumw2();
  //K0s_K0s_cosThetaProdPlane_US_hist->Divide(K0s_K0s_cosThetaProdPlane_eff);
  K0s_K0s_cosThetaProdPlane_US_hist->Scale(1./K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->GetBinWidth(1));
  //K0s_K0s_cosThetaProdPlane_US_hist->Scale(1./K0s_K0s_cosThetaProdPlane_US_hist->Integral());
  K0s_K0s_cosThetaProdPlane_US_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_US_hist->Draw("p e");

  //ME here just for plotting, used loser
  TH1D *K0s_K0s_cosThetaProdPlane_ME_hist = K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist->ProjectionX( "proj_LL_US_ME", 1, 5); // |Delta_phi| < pi/4

  K0s_K0s_cosThetaProdPlane_ME_hist->SetMarkerSize(1.5);
  K0s_K0s_cosThetaProdPlane_ME_hist->SetMarkerStyle(24);
  K0s_K0s_cosThetaProdPlane_ME_hist->SetMarkerColor(1);
  K0s_K0s_cosThetaProdPlane_ME_hist->SetLineColor(1);
  K0s_K0s_cosThetaProdPlane_ME_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_ME_hist->Scale(nLL/K0s_K0s_cosThetaProdPlane_ME_hist->Integral()); //scale ME to US
  K0s_K0s_cosThetaProdPlane_ME_hist->Scale(1./K0s_K0s_cosThetaProdPlane_ME_hist->GetXaxis()->GetBinWidth(1)); //US is scaled by bin width, ME needs to be scaled as well
  K0s_K0s_cosThetaProdPlane_ME_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_ME_hist->Draw("p e same");

  TF1 *fitK0s_K0s_US_ThetaStar_no_corr_ME = new TF1("fitK0s_K0s_US_ThetaStar_no_corr_ME", "[0]*(1 + [1]*x)", -1, 1);
  fitK0s_K0s_US_ThetaStar_no_corr_ME->SetParameters(100, 0.5);

  K0s_K0s_cosThetaProdPlane_ME_hist->Fit(fitK0s_K0s_US_ThetaStar_no_corr_ME, "s i 0 r");

  float P_K0s_K0s_no_corr_ME = fitK0s_K0s_US_ThetaStar_no_corr_ME->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_K0s_K0s_no_corr_ME_err = fitK0s_K0s_US_ThetaStar_no_corr_ME->GetParError(1)/(L0_alpha*L0_alpha);


  TH1D *K0s_K0s_cosThetaProdPlane_LS_hist = K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist->ProjectionX( "proj_LL_LS_%i", 1, 5); // |Delta_phi| < pi/4

  K0s_K0s_cosThetaProdPlane_LS_hist->SetMarkerSize(1.5);
  K0s_K0s_cosThetaProdPlane_LS_hist->SetMarkerStyle(21);
  K0s_K0s_cosThetaProdPlane_LS_hist->SetMarkerColor(kBlue);
  double nLL_back = K0s_K0s_cosThetaProdPlane_LS_hist->Integral();
  K0s_K0s_cosThetaProdPlane_LS_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_LS_hist->Scale(1./K0s_K0s_cosThetaProdPlane_LS_hist->GetXaxis()->GetBinWidth(1));
  //K0s_K0s_cosThetaProdPlane_LS_hist->Divide(K0s_K0s_cosThetaProdPlane_eff);
  K0s_K0s_cosThetaProdPlane_LS_hist->Draw("p e same");


  TH1D *K0s_K0s_cosThetaProdPlane_ME_LS_hist = K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist->ProjectionX( "proj_LL_LS_ME", 1, 5); // |Delta_phi| < pi/4

  K0s_K0s_cosThetaProdPlane_ME_LS_hist->SetMarkerSize(1.5);
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->SetMarkerStyle(25);
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->SetMarkerColor(kMagenta+1);
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->SetLineColor(kMagenta+1);
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->Scale(nLL_back/K0s_K0s_cosThetaProdPlane_ME_LS_hist->Integral()); //scale ME_LS to background
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->Scale(1./K0s_K0s_cosThetaProdPlane_ME_LS_hist->GetXaxis()->GetBinWidth(1)); //background is scaled by bin width, ME_LS needs to be scaled as well
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->Draw("p e same");


  TF1 *fitK0s_K0s_US_ThetaStar_no_corr = new TF1("fitK0s_K0s_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
  fitK0s_K0s_US_ThetaStar_no_corr->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  K0s_K0s_cosThetaProdPlane_US_hist->Fit(fitK0s_K0s_US_ThetaStar_no_corr, "s i 0 r");

  float P_K0s_K0s_no_corr = fitK0s_K0s_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_K0s_K0s_no_corr_err = fitK0s_K0s_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0_alpha);

  fitK0s_K0s_US_ThetaStar_no_corr->SetLineColor(1);
  //fitK0s_K0s_US_ThetaStar_no_corr->Draw("same");

  TLegend *K0s_K0s_leg = new TLegend(0.15, 0.3, 0.45, 0.54);
  K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_US_hist, "(US-US) p#pi");
  K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_ME_hist, "(US-US) ME");
  K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_LS_hist, "Combinatorial bckg.");
  K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_ME_LS_hist, "Bckg. ME");
  //K0s_K0s_leg->AddEntry(fitK0s_K0s_US_ThetaStar_no_corr, "Linear fit to US");
  K0s_K0s_leg->SetBorderSize(0);
  K0s_K0s_leg->SetFillColorAlpha(0, 0.01);
  K0s_K0s_leg->Draw("same");

  TPaveText *K0s_K0s_text_no_corr = new TPaveText(0.5, 0.3, 0.85, 0.65, "NDC");
  K0s_K0s_text_no_corr->SetTextFont(42);
  //K0s_K0s_text_no_corr->AddText("STAR Internal");
  //K0s_K0s_text_no_corr->AddText("STAR preliminary");
  //((TText*)K0s_K0s_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  K0s_K0s_text_no_corr->AddText("Minimum bias, no correction");
  K0s_K0s_text_no_corr->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_text_no_corr->AddText("|#it{y}| < 1");
  K0s_K0s_text_no_corr->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
  K0s_K0s_text_no_corr->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
  K0s_K0s_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_K0s_K0s_no_corr, fabs(P_K0s_K0s_no_corr_err)));
  K0s_K0s_text_no_corr->AddText(Form("P_{ME} = %.3f #pm %.3f", P_K0s_K0s_no_corr_ME, fabs(P_K0s_K0s_no_corr_ME_err)));
  K0s_K0s_text_no_corr->SetFillColorAlpha(0, 0.01);
  K0s_K0s_text_no_corr->Draw("same");

  K0s_K0s_cosThetaProdPlane_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations_delta_phi_new/K0s_K0s_cosThetaProdPlane_no_corr_delta_phi_new_%i.png", 0));

  //----------------------------------------------------------------------------------------------------

  TCanvas *K0s_K0s_cosThetaProdPlane_can = new TCanvas("K0s_K0s_cosThetaProdPlane_can", "K0s_K0s_cosThetaProdPlane_can", 1200, 1000);

  K0s_K0s_cosThetaProdPlane_can->cd();

  //ME histogram higher

  K0s_K0s_cosThetaProdPlane_US_hist->Divide(K0s_K0s_cosThetaProdPlane_ME_hist); // correct US using ME
  K0s_K0s_cosThetaProdPlane_US_hist->Scale(nLL/K0s_K0s_cosThetaProdPlane_US_hist->Integral()); //scale back to raw US
  K0s_K0s_cosThetaProdPlane_US_hist->Scale(1./K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
  K0s_K0s_cosThetaProdPlane_US_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_US_hist->Draw("p e");

  //K0s_K0s_cosThetaProdPlane_ME_hist->Draw("same p e");

  K0s_K0s_cosThetaProdPlane_LS_hist->Divide(K0s_K0s_cosThetaProdPlane_ME_LS_hist); //correct background using ME
  K0s_K0s_cosThetaProdPlane_LS_hist->Scale(nLL_back/K0s_K0s_cosThetaProdPlane_LS_hist->Integral()); //scale back to raw background
  K0s_K0s_cosThetaProdPlane_LS_hist->Scale(1./K0s_K0s_cosThetaProdPlane_LS_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
  K0s_K0s_cosThetaProdPlane_LS_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_LS_hist->Draw("p e same");

  //fit dN/dcos(theta*)
  //signal + bacground
  TF1 *fitK0s_K0s_US_ThetaStar = new TF1("fitK0s_K0s_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fitK0s_K0s_US_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  K0s_K0s_cosThetaProdPlane_US_hist->Fit(fitK0s_K0s_US_ThetaStar, "s i 0 r");

  float P_K0s_K0s = fitK0s_K0s_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_K0s_K0s_err = fitK0s_K0s_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

  fitK0s_K0s_US_ThetaStar->SetLineColor(1);
  fitK0s_K0s_US_ThetaStar->Draw("same");

  //background
  TF1 *fitK0s_K0s_US_LS_ThetaStar = new TF1("fitK0s_K0s_US_LS_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fitK0s_K0s_US_LS_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaUS_LS_wrong_sign = L_inv_mass_US_LS->Fit(fitGaUS_LSsBack, "s i 0", "", 1.07, 1.4);
  K0s_K0s_cosThetaProdPlane_LS_hist->Fit(fitK0s_K0s_US_LS_ThetaStar, "s i 0 r");

  float P_K0s_K0s_back = fitK0s_K0s_US_LS_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_K0s_K0s_back_err = fitK0s_K0s_US_LS_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

  fitK0s_K0s_US_LS_ThetaStar->SetLineColor(1);
  fitK0s_K0s_US_LS_ThetaStar->Draw("same");


  TPaveText *K0s_K0s_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
  K0s_K0s_text->SetTextFont(42);
  //K0s_K0s_text->AddText("STAR Internal");
  //K0s_K0s_text->AddText("STAR preliminary");
  //((TText*)K0s_K0s_text->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  K0s_K0s_text->AddText("Minimum bias");
  K0s_K0s_text->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_text->AddText("|#it{y}| < 1");
  K0s_K0s_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  K0s_K0s_text->AddText(Form("P_{tot} = %.2f #pm %.2f", P_K0s_K0s, fabs(P_K0s_K0s_err)));
  K0s_K0s_text->AddText(Form("P_{bckg} = %.2f #pm %.2f", P_K0s_K0s_back, fabs(P_K0s_K0s_back_err)));
  K0s_K0s_text->SetFillColorAlpha(0, 0.01);
  K0s_K0s_text->Draw("same");

  K0s_K0s_leg->Draw("same");

  K0s_K0s_cosThetaProdPlane_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations_delta_phi_new/K0s_K0s_cosThetaProdPlane_delta_phi_new_%i.png", 0));

  //----------------------------------------------------------------------

  TCanvas *K0s_K0s_cosThetaProdPlane_can_2 = new TCanvas("K0s_K0s_cosThetaProdPlane_can_2", "K0s_K0s_cosThetaProdPlane_can_2", 1200, 1000);

  K0s_K0s_cosThetaProdPlane_can_2->cd();

  K0s_K0s_cosThetaProdPlane_US_hist->Add(K0s_K0s_cosThetaProdPlane_LS_hist, -1);
  K0s_K0s_cosThetaProdPlane_US_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_US_hist->Write(Form("K0s_K0s_cosThetaProdPlane_delta_phi_%i", 0));
  K0s_K0s_cosThetaProdPlane_US_hist->Draw("p e");


  TF1 *fitK0s_K0s_US_ThetaStar_2 = new TF1("fitK0s_K0s_US_ThetaStar_2", "[0]*(1 + [1]*x)", -1, 1);
  fitK0s_K0s_US_ThetaStar_2->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  K0s_K0s_cosThetaProdPlane_US_hist->Fit(fitK0s_K0s_US_ThetaStar_2, "s i 0 r");
/*
  //store fit result for systematic errors
  //tight topo cuts
  if(cut_type == 1)
  {
    SysErrCutsTopo->cd();

    fitK0s_K0s_tight_topo_cuts_Delta_phi[0]->SetParameters(fitK0s_K0s_US_ThetaStar_2->GetParameter(0), fitK0s_K0s_US_ThetaStar_2->GetParameter(1));
    fitK0s_K0s_tight_topo_cuts_Delta_phi[0]->SetParError(0, fitK0s_K0s_US_ThetaStar_2->GetParError(0));
    fitK0s_K0s_tight_topo_cuts_Delta_phi[0]->SetParError(1, fitK0s_K0s_US_ThetaStar_2->GetParError(1));
    fitK0s_K0s_tight_topo_cuts_Delta_phi[0]->Write();
  }

  //tight pT cuts
  if(cut_type == 2)
  {
    SysErrCutsPt->cd();

    fitK0s_K0s_tight_pT_cuts_Delta_phi[0]->SetParameters(fitK0s_K0s_US_ThetaStar_2->GetParameter(0), fitK0s_K0s_US_ThetaStar_2->GetParameter(1));
    fitK0s_K0s_tight_pT_cuts_Delta_phi[0]->SetParError(0, fitK0s_K0s_US_ThetaStar_2->GetParError(0));
    fitK0s_K0s_tight_pT_cuts_Delta_phi[0]->SetParError(1, fitK0s_K0s_US_ThetaStar_2->GetParError(1));
    fitK0s_K0s_tight_pT_cuts_Delta_phi[0]->Write();
  }
*/
  float P_K0s_K0s_2 = fitK0s_K0s_US_ThetaStar_2->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_K0s_K0s_err_2 = fitK0s_K0s_US_ThetaStar_2->GetParError(1)/(L0_alpha*L0_alpha);

  fitK0s_K0s_US_ThetaStar_2->SetLineColor(1);
  fitK0s_K0s_US_ThetaStar_2->Draw("same");


  //calculate total systematic uncertainty for L-Lbar

  //slope
  //float SysErrSlope_K0s_K0s = SysErrSlope_delta_phi_less_hist[0]->GetBinContent(1);
  float SysErrSlope_K0s_K0s = SysErrSlope_delta_eta_hist[0]->GetBinContent(2);

  //----------------------------

  //background subtraction
  float P_K0s_K0s_from_fits = P_K0s_K0s - nLL_back/nLL*P_K0s_K0s_back;
  float P_K0s_K0s_from_fits_err = sqrt( P_K0s_K0s_err*P_K0s_K0s_err + nLL_back/nLL*nLL_back/nLL*P_K0s_K0s_back_err*P_K0s_K0s_back_err );

  //statistical error correction for systematic error
  float SysErrBackground_K0s_K0s_corr = sqrt( fabs( P_K0s_K0s_err_2*P_K0s_K0s_err_2 - P_K0s_K0s_from_fits_err*P_K0s_K0s_from_fits_err ) );

  float SysErrBackground_K0s_K0s_work = ( fabs( P_K0s_K0s_2 - P_K0s_K0s_from_fits) - SysErrBackground_K0s_K0s_corr )/fabs(P_K0s_K0s_2);

  float SysErrBackground_K0s_K0s = 0;

  if( SysErrBackground_K0s_K0s_work > 0 ) SysErrBackground_K0s_K0s = SysErrBackground_K0s_K0s_work; //store sys. err. only if it is larger than statistical fluctuations

  //----------------------------

  //cuts variation
  float SysErr_tight_cuts_K0s_K0s = 0;

  //float P_K0s_K0s_tight_topo_cuts = 0;
  //float P_K0s_K0s_tight_pT_cuts = 0;
/*
  if(cut_type == 0)
  {
    float P_K0s_K0s_tight_topo_cuts = fitK0s_K0s_tight_topo_cuts_Delta_phi[0]->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_K0s_K0s_tight_topo_cuts_err = fitK0s_K0s_tight_topo_cuts_Delta_phi[0]->GetParError(1)/(L0_alpha*L0_alpha);

    //statistical error correction for systematic error
    float SysErr_tight_topo_cuts_K0s_K0s_corr = sqrt( fabs( P_K0s_K0s_err_2*P_K0s_K0s_err_2 - P_K0s_K0s_tight_topo_cuts_err*P_K0s_K0s_tight_topo_cuts_err ) );

    float SysErr_tight_topo_cuts_K0s_K0s_work = ( fabs( P_K0s_K0s_2 - P_K0s_K0s_tight_topo_cuts) - SysErr_tight_topo_cuts_K0s_K0s_corr )/fabs(P_K0s_K0s_2);

    float SysErr_tight_topo_cuts_K0s_K0s = 0;

    if( SysErr_tight_topo_cuts_K0s_K0s_work > 0 ) SysErr_tight_topo_cuts_K0s_K0s = SysErr_tight_topo_cuts_K0s_K0s_work; //store sys. err. only if it is larger than statistical fluctuations


    float P_K0s_K0s_tight_pT_cuts = fitK0s_K0s_tight_pT_cuts_Delta_phi[0]->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_K0s_K0s_tight_pT_cuts_err = fitK0s_K0s_tight_pT_cuts_Delta_phi[0]->GetParError(1)/(L0_alpha*L0_alpha);

    //statistical error correction for systematic error
    float SysErr_tight_pT_cuts_K0s_K0s_corr = sqrt( fabs( P_K0s_K0s_err_2*P_K0s_K0s_err_2 - P_K0s_K0s_tight_pT_cuts_err*P_K0s_K0s_tight_pT_cuts_err ) );

    float SysErr_tight_pT_cuts_K0s_K0s_work = ( fabs( P_K0s_K0s_2 - P_K0s_K0s_tight_pT_cuts) - SysErr_tight_pT_cuts_K0s_K0s_corr )/fabs(P_K0s_K0s_2);

    float SysErr_tight_pT_cuts_K0s_K0s = 0;

    if( SysErr_tight_pT_cuts_K0s_K0s_work > 0 ) SysErr_tight_pT_cuts_K0s_K0s = SysErr_tight_pT_cuts_K0s_K0s_work; //store sys. err. only if it is larger than statistical fluctuations

    SysErr_tight_cuts_K0s_K0s = sqrt(SysErr_tight_topo_cuts_K0s_K0s*SysErr_tight_topo_cuts_K0s_K0s + SysErr_tight_pT_cuts_K0s_K0s*SysErr_tight_pT_cuts_K0s_K0s);
    //SysErr_tight_cuts_K0s_K0s = SysErr_tight_topo_cuts_K0s_K0s;
    //SysErr_tight_cuts_K0s_K0s = SysErr_tight_pT_cuts_K0s_K0s;
  }
*/
  //total
  float SysErrTot_K0s_K0s = sqrt( sysErr_alpha_K0s_K0s*sysErr_alpha_K0s_K0s + SysErrSlope_K0s_K0s*SysErrSlope_K0s_K0s + SysErrBackground_K0s_K0s*SysErrBackground_K0s_K0s  );
  //float SysErrTot_K0s_K0s = sqrt( sysErr_alpha_K0s_K0s*sysErr_alpha_K0s_K0s + SysErrSlope_K0s_K0s*SysErrSlope_K0s_K0s + SysErrBackground_K0s_K0s*SysErrBackground_K0s_K0s + SysErr_tight_cuts_K0s_K0s*SysErr_tight_cuts_K0s_K0s  );




  //TPaveText *K0s_K0s_text_2 = new TPaveText(0.47, 0.15, 0.85, 0.53, "NDC"); //with polarization value
  TPaveText *K0s_K0s_text_2 = new TPaveText(0.47, 0.15, 0.8, 0.53, "NDC"); //without polarization value
  K0s_K0s_text_2->SetTextFont(42);
  //K0s_K0s_text_2->SetTextSize(15);
  //K0s_K0s_text_2->AddText("STAR");
  //K0s_K0s_text_2->AddText("STAR preliminary");
  //((TText*)K0s_K0s_text_2->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_text_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  K0s_K0s_text_2->AddText("Minimum bias");
  K0s_K0s_text_2->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_text_2->AddText("|#it{y}| < 1");
  K0s_K0s_text_2->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
  K0s_K0s_text_2->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
  K0s_K0s_text_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2}} = %.3f #pm %.3f", P_K0s_K0s_2, fabs(P_K0s_K0s_err_2)));
  //K0s_K0s_text_2->AddText(Form("P_{topo} = %.2f", P_K0s_K0s_tight_topo_cuts));
  //K0s_K0s_text_2->AddText(Form("P_{pT} = %.2f", P_K0s_K0s_tight_pT_cuts));
  K0s_K0s_text_2->SetFillColorAlpha(0, 0.01);
  K0s_K0s_text_2->Draw("same");


  TLegend *K0s_K0s_2_leg = new TLegend(0.15, 0.3, 0.4, 0.49);
  //K0s_K0s_2_leg->SetTextSizePixels(15);
  K0s_K0s_2_leg->AddEntry(K0s_K0s_cosThetaProdPlane_US_hist, "(US-US)-Bckg.");
  K0s_K0s_2_leg->AddEntry(fitK0s_K0s_US_ThetaStar_2, "Fit", "l");
  K0s_K0s_2_leg->SetBorderSize(0);
  K0s_K0s_2_leg->SetFillColorAlpha(0, 0.01);
  K0s_K0s_2_leg->Draw("same");

  K0s_K0s_cosThetaProdPlane_can_2->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations_delta_phi_new/K0s_K0s_cosThetaProdPlane_subtract_delta_phi_new_%i.png", 0));


  PolarizationGraph_delta_phi[0]->SetPoint(1, P_K0s_K0s_2, 4);
  PolarizationGraph_delta_phi[0]->SetPointError(1, fabs(P_K0s_K0s_err_2), 0);

  PolarizationGraph_delta_phi_sys_err[0]->SetPoint(1, P_K0s_K0s_2, 4);
  PolarizationGraph_delta_phi_sys_err[0]->SetPointError(1, fabs(SysErrTot_K0s_K0s*P_K0s_K0s_2), 0.045);

  //____________________________________________________________________________________________________________________________________________________________________________________________________________


  TCanvas *K0s_K0s_cosThetaProdPlane_no_corr_2_can = new TCanvas("K0s_K0s_cosThetaProdPlane_no_corr_2_can", "K0s_K0s_cosThetaProdPlane_no_corr_2_can", 1200, 1000);

  K0s_K0s_cosThetaProdPlane_no_corr_2_can->cd();

  TH1D *K0s_K0s_cosThetaProdPlane_US_2_hist = K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_hist->ProjectionX( "proj_LL_US_2", 6, 20);

  K0s_K0s_cosThetaProdPlane_US_2_hist->GetXaxis()->SetTitle("cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_US_2_hist->GetXaxis()->CenterTitle();
  //K0s_K0s_cosThetaProdPlane_US_2_hist->GetXaxis()->SetTextSizePixels(30);
  K0s_K0s_cosThetaProdPlane_US_2_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_US_2_hist->GetYaxis()->CenterTitle();
  //K0s_K0s_cosThetaProdPlane_US_2_hist->GetYaxis()->SetTextSizePixels(30);
  K0s_K0s_cosThetaProdPlane_US_2_hist->GetYaxis()->SetMaxDigits(3);
  K0s_K0s_cosThetaProdPlane_US_2_hist->SetMarkerSize(1.5);
  K0s_K0s_cosThetaProdPlane_US_2_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_US_2_hist->SetMarkerColor(kRed);
  K0s_K0s_cosThetaProdPlane_US_2_hist->SetLineColor(kRed);
  double nLL_2 = K0s_K0s_cosThetaProdPlane_US_2_hist->Integral();
  K0s_K0s_cosThetaProdPlane_US_2_hist->Sumw2();
  //K0s_K0s_cosThetaProdPlane_US_2_hist->Divide(K0s_K0s_cosThetaProdPlane_eff);
  K0s_K0s_cosThetaProdPlane_US_2_hist->Scale(1./K0s_K0s_cosThetaProdPlane_US_2_hist->GetXaxis()->GetBinWidth(1));
  //K0s_K0s_cosThetaProdPlane_US_2_hist->Scale(1./K0s_K0s_cosThetaProdPlane_US_2_hist->Integral());
  K0s_K0s_cosThetaProdPlane_US_2_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_US_2_hist->Draw("p e");

  //ME here just for plotting, used loser
  TH1D *K0s_K0s_cosThetaProdPlane_ME_2_hist = K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist->ProjectionX( "proj_LL_US_ME_2", 6, 20);

  K0s_K0s_cosThetaProdPlane_ME_2_hist->SetMarkerSize(1.5);
  K0s_K0s_cosThetaProdPlane_ME_2_hist->SetMarkerStyle(24);
  K0s_K0s_cosThetaProdPlane_ME_2_hist->SetMarkerColor(1);
  K0s_K0s_cosThetaProdPlane_ME_2_hist->SetLineColor(1);
  K0s_K0s_cosThetaProdPlane_ME_2_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_ME_2_hist->Scale(nLL_2/K0s_K0s_cosThetaProdPlane_ME_2_hist->Integral()); //scale ME to US
  K0s_K0s_cosThetaProdPlane_ME_2_hist->Scale(1./K0s_K0s_cosThetaProdPlane_ME_2_hist->GetXaxis()->GetBinWidth(1)); //US is scaled by bin width, ME needs to be scaled as well
  K0s_K0s_cosThetaProdPlane_ME_2_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_ME_2_hist->Draw("p e same");

  TF1 *fit_2_K0s_K0s_US_ThetaStar_no_corr_ME = new TF1("fit_2_K0s_K0s_US_ThetaStar_no_corr_ME", "[0]*(1 + [1]*x)", -1, 1);
  fit_2_K0s_K0s_US_ThetaStar_no_corr_ME->SetParameters(100, 0.5);

  K0s_K0s_cosThetaProdPlane_ME_2_hist->Fit(fit_2_K0s_K0s_US_ThetaStar_no_corr_ME, "s i 0 r");

  float P_2_K0s_K0s_no_corr_ME = fit_2_K0s_K0s_US_ThetaStar_no_corr_ME->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_2_K0s_K0s_no_corr_ME_err = fit_2_K0s_K0s_US_ThetaStar_no_corr_ME->GetParError(1)/(L0_alpha*L0_alpha);


  TH1D *K0s_K0s_cosThetaProdPlane_LS_2_hist = K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist->ProjectionX( "proj_LL_LS_2", 6, 20);

  K0s_K0s_cosThetaProdPlane_LS_2_hist->SetMarkerSize(1.5);
  K0s_K0s_cosThetaProdPlane_LS_2_hist->SetMarkerStyle(21);
  K0s_K0s_cosThetaProdPlane_LS_2_hist->SetMarkerColor(kBlue);
  double nLL_2_back = K0s_K0s_cosThetaProdPlane_LS_2_hist->Integral();
  K0s_K0s_cosThetaProdPlane_LS_2_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_LS_2_hist->Scale(1./K0s_K0s_cosThetaProdPlane_LS_2_hist->GetXaxis()->GetBinWidth(1));
  //K0s_K0s_cosThetaProdPlane_LS_2_hist->Divide(K0s_K0s_cosThetaProdPlane_eff);
  K0s_K0s_cosThetaProdPlane_LS_2_hist->Draw("p e same");


  TH1D *K0s_K0s_cosThetaProdPlane_ME_LS_2_hist = K0s_K0s_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist->ProjectionX( "proj_LL_LS_ME_2", 6, 20);

  K0s_K0s_cosThetaProdPlane_ME_LS_2_hist->SetMarkerSize(1.5);
  K0s_K0s_cosThetaProdPlane_ME_LS_2_hist->SetMarkerStyle(25);
  K0s_K0s_cosThetaProdPlane_ME_LS_2_hist->SetMarkerColor(kMagenta+1);
  K0s_K0s_cosThetaProdPlane_ME_LS_2_hist->SetLineColor(kMagenta+1);
  K0s_K0s_cosThetaProdPlane_ME_LS_2_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_ME_LS_2_hist->Scale(nLL_2_back/K0s_K0s_cosThetaProdPlane_ME_LS_2_hist->Integral()); //scale ME_LS to background
  K0s_K0s_cosThetaProdPlane_ME_LS_2_hist->Scale(1./K0s_K0s_cosThetaProdPlane_ME_LS_2_hist->GetXaxis()->GetBinWidth(1)); //background is scaled by bin width, ME_LS needs to be scaled as well
  K0s_K0s_cosThetaProdPlane_ME_LS_2_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_ME_LS_2_hist->Draw("p e same");


  TF1 *fit_2_K0s_K0s_US_ThetaStar_no_corr = new TF1("fit_2_K0s_K0s_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
  fit_2_K0s_K0s_US_ThetaStar_no_corr->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  K0s_K0s_cosThetaProdPlane_US_2_hist->Fit(fit_2_K0s_K0s_US_ThetaStar_no_corr, "s i 0 r");

  float P_2_K0s_K0s_no_corr = fit_2_K0s_K0s_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_2_K0s_K0s_no_corr_err = fit_2_K0s_K0s_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0_alpha);

  fit_2_K0s_K0s_US_ThetaStar_no_corr->SetLineColor(1);
  //fit_2_K0s_K0s_US_ThetaStar_no_corr->Draw("same");

  TLegend *K0s_K0s_more_leg = new TLegend(0.55, 0.3, 0.85, 0.54);
  K0s_K0s_more_leg->AddEntry(K0s_K0s_cosThetaProdPlane_US_hist, "(US-US) p#pi");
  K0s_K0s_more_leg->AddEntry(K0s_K0s_cosThetaProdPlane_ME_hist, "(US-US) ME");
  K0s_K0s_more_leg->AddEntry(K0s_K0s_cosThetaProdPlane_LS_hist, "Combinatorial bckg.");
  K0s_K0s_more_leg->AddEntry(K0s_K0s_cosThetaProdPlane_ME_LS_hist, "Bckg. ME");
  //K0s_K0s_more_leg->AddEntry(fitK0s_K0s_US_ThetaStar_no_corr, "Linear fit to US");
  K0s_K0s_more_leg->SetBorderSize(0);
  K0s_K0s_more_leg->SetFillColorAlpha(0, 0.01);
  K0s_K0s_more_leg->Draw("same");

  TPaveText *K0s_K0s_text_no_corr_2 = new TPaveText(0.15, 0.3, 0.55, 0.65, "NDC");
  K0s_K0s_text_no_corr_2->SetTextFont(42);
  //K0s_K0s_text_no_corr_2->AddText("STAR Internal");
  //K0s_K0s_text_no_corr_2->AddText("STAR preliminary");
  //((TText*)K0s_K0s_text_no_corr_2->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_text_no_corr_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  K0s_K0s_text_no_corr_2->AddText("Minimum bias, no correction");
  K0s_K0s_text_no_corr_2->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_text_no_corr_2->AddText("|#it{y}| < 1");
  K0s_K0s_text_no_corr_2->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
  K0s_K0s_text_no_corr_2->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
  K0s_K0s_text_no_corr_2->AddText(Form("P = %.2f #pm %.2f", P_2_K0s_K0s_no_corr, fabs(P_2_K0s_K0s_no_corr_err)));
  K0s_K0s_text_no_corr_2->AddText(Form("P_{ME} = %.3f #pm %.3f", P_2_K0s_K0s_no_corr_ME, fabs(P_2_K0s_K0s_no_corr_ME_err)));
  K0s_K0s_text_no_corr_2->SetFillColorAlpha(0, 0.01);
  K0s_K0s_text_no_corr_2->Draw("same");
  K0s_K0s_text_no_corr_2->Draw("same");

  K0s_K0s_cosThetaProdPlane_no_corr_2_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations_delta_phi_new/K0s_K0s_cosThetaProdPlane_no_corr_delta_phi_new_%i.png", 1));

  //----------------------------------------------------------------------------------------------------

  TCanvas *K0s_K0s_cosThetaProdPlane_2_can = new TCanvas("K0s_K0s_cosThetaProdPlane_2_can", "K0s_K0s_cosThetaProdPlane_2_can", 1200, 1000);

  K0s_K0s_cosThetaProdPlane_2_can->cd();

  //ME histogram higher

  K0s_K0s_cosThetaProdPlane_US_2_hist->Divide(K0s_K0s_cosThetaProdPlane_ME_2_hist); // correct US using ME
  K0s_K0s_cosThetaProdPlane_US_2_hist->Scale(nLL_2/K0s_K0s_cosThetaProdPlane_US_2_hist->Integral()); //scale back to raw US
  K0s_K0s_cosThetaProdPlane_US_2_hist->Scale(1./K0s_K0s_cosThetaProdPlane_US_2_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
  K0s_K0s_cosThetaProdPlane_US_2_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_US_2_hist->Draw("p e");

  //K0s_K0s_cosThetaProdPlane_ME_2_hist->Draw("same p e");

  K0s_K0s_cosThetaProdPlane_LS_2_hist->Divide(K0s_K0s_cosThetaProdPlane_ME_LS_2_hist); //correct background using ME
  K0s_K0s_cosThetaProdPlane_LS_2_hist->Scale(nLL_2_back/K0s_K0s_cosThetaProdPlane_LS_2_hist->Integral()); //scale back to raw background
  K0s_K0s_cosThetaProdPlane_LS_2_hist->Scale(1./K0s_K0s_cosThetaProdPlane_LS_2_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
  K0s_K0s_cosThetaProdPlane_LS_2_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_LS_2_hist->Draw("p e same");

  //fit dN/dcos(theta*)
  //signal + bacground
  TF1 *fit_2_K0s_K0s_US_ThetaStar = new TF1("fit_2_K0s_K0s_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fit_2_K0s_K0s_US_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  K0s_K0s_cosThetaProdPlane_US_2_hist->Fit(fit_2_K0s_K0s_US_ThetaStar, "s i 0 r");

  float P_2_K0s_K0s = fit_2_K0s_K0s_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_2_K0s_K0s_err = fit_2_K0s_K0s_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

  fit_2_K0s_K0s_US_ThetaStar->SetLineColor(1);
  fit_2_K0s_K0s_US_ThetaStar->Draw("same");

  //background
  TF1 *fit_2_K0s_K0s_US_LS_ThetaStar = new TF1("fit_2_K0s_K0s_US_LS_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fit_2_K0s_K0s_US_LS_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaUS_LS_wrong_sign = L_inv_mass_US_LS->Fit(fitGaUS_LSsBack, "s i 0", "", 1.07, 1.4);
  K0s_K0s_cosThetaProdPlane_LS_2_hist->Fit(fit_2_K0s_K0s_US_LS_ThetaStar, "s i 0 r");

  float P_2_K0s_K0s_back = fit_2_K0s_K0s_US_LS_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_2_K0s_K0s_back_err = fit_2_K0s_K0s_US_LS_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

  fit_2_K0s_K0s_US_LS_ThetaStar->SetLineColor(1);
  fit_2_K0s_K0s_US_LS_ThetaStar->Draw("same");

  TPaveText *K0s_K0s_text_more = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
  K0s_K0s_text_more->SetTextFont(42);
  //K0s_K0s_text_more->AddText("STAR Internal");
  //K0s_K0s_text_more->AddText("STAR preliminary");
  //((TText*)K0s_K0s_text_more->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_text_more->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  K0s_K0s_text_more->AddText("Minimum bias");
  K0s_K0s_text_more->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_text_more->AddText("|#it{y}| < 1");
  K0s_K0s_text_more->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  K0s_K0s_text_more->AddText(Form("P_{tot} = %.2f #pm %.2f", P_2_K0s_K0s, fabs(P_2_K0s_K0s_err)));
  K0s_K0s_text_more->AddText(Form("P_{bckg} = %.2f #pm %.2f", P_2_K0s_K0s_back, fabs(P_2_K0s_K0s_back_err)));
  K0s_K0s_text_more->SetFillColorAlpha(0, 0.01);
  K0s_K0s_text_more->Draw("same");


  K0s_K0s_more_leg->Draw("same");

  K0s_K0s_cosThetaProdPlane_2_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations_delta_phi_new/K0s_K0s_cosThetaProdPlane_delta_phi_new_%i.png", 1));

  //----------------------------------------------------------------------

  TCanvas *K0s_K0s_cosThetaProdPlane_2_can_2 = new TCanvas("K0s_K0s_cosThetaProdPlane_2_can_2", "K0s_K0s_cosThetaProdPlane_2_can_2", 1200, 1000);

  K0s_K0s_cosThetaProdPlane_2_can_2->cd();

  K0s_K0s_cosThetaProdPlane_US_2_hist->Add(K0s_K0s_cosThetaProdPlane_LS_2_hist, -1);
  K0s_K0s_cosThetaProdPlane_US_2_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_US_2_hist->Write(Form("K0s_K0s_cosThetaProdPlane_delta_phi_%i", 1));
  K0s_K0s_cosThetaProdPlane_US_2_hist->Draw("p e");


  TF1 *fit_2_K0s_K0s_US_ThetaStar_2 = new TF1("fit_2_K0s_K0s_US_ThetaStar_2", "[0]*(1 + [1]*x)", -1, 1);
  fit_2_K0s_K0s_US_ThetaStar_2->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  K0s_K0s_cosThetaProdPlane_US_2_hist->Fit(fit_2_K0s_K0s_US_ThetaStar_2, "s i 0 r");
/*
  //store fit result for systematic errors
  //tight topo cuts
  if(cut_type == 1)
  {
    SysErrCutsTopo->cd();

    fitK0s_K0s_tight_topo_cuts_Delta_phi[1]->SetParameters(fit_2_K0s_K0s_US_ThetaStar_2->GetParameter(0), fitK0s_K0s_US_ThetaStar_2->GetParameter(1));
    fitK0s_K0s_tight_topo_cuts_Delta_phi[1]->SetParError(0, fit_2_K0s_K0s_US_ThetaStar_2->GetParError(0));
    fitK0s_K0s_tight_topo_cuts_Delta_phi[1]->SetParError(1, fit_2_K0s_K0s_US_ThetaStar_2->GetParError(1));
    fitK0s_K0s_tight_topo_cuts_Delta_phi[1]->Write();
  }

  //tight pT cuts
  if(cut_type == 2)
  {
    SysErrCutsPt->cd();

    fitK0s_K0s_tight_pT_cuts_Delta_phi[1]->SetParameters(fit_2_K0s_K0s_US_ThetaStar_2->GetParameter(0), fitK0s_K0s_US_ThetaStar_2->GetParameter(1));
    fitK0s_K0s_tight_pT_cuts_Delta_phi[1]->SetParError(0, fit_2_K0s_K0s_US_ThetaStar_2->GetParError(0));
    fitK0s_K0s_tight_pT_cuts_Delta_phi[1]->SetParError(1, fit_2_K0s_K0s_US_ThetaStar_2->GetParError(1));
    fitK0s_K0s_tight_pT_cuts_Delta_phi[1]->Write();
  }
*/
  float P_2_K0s_K0s_2 = fit_2_K0s_K0s_US_ThetaStar_2->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_2_K0s_K0s_err_2 = fit_2_K0s_K0s_US_ThetaStar_2->GetParError(1)/(L0_alpha*L0_alpha);

  fit_2_K0s_K0s_US_ThetaStar_2->SetLineColor(1);
  fit_2_K0s_K0s_US_ThetaStar_2->Draw("same");


  //calculate total systematic uncertainty for L-Lbar

  //slope
  //float SysErrSlope_K0s_K0s_2 = SysErrSlope_delta_phi_more_hist[0]->GetBinContent(1);
  float SysErrSlope_K0s_K0s_2 = SysErrSlope_delta_eta_hist[0]->GetBinContent(2);

  //----------------------------

  //background subtraction
  float P_2_K0s_K0s_from_fits = P_2_K0s_K0s - nLL_2_back/nLL_2*P_K0s_K0s_back;
  float P_2_K0s_K0s_from_fits_err = sqrt( P_2_K0s_K0s_err*P_2_K0s_K0s_err + nLL_2_back/nLL_2*nLL_2_back/nLL_2*P_2_K0s_K0s_back_err*P_2_K0s_K0s_back_err );

  //statistical error correction for systematic error
  float SysErrBackground_K0s_K0s_corr_2 = sqrt( fabs( P_2_K0s_K0s_err_2*P_2_K0s_K0s_err_2 - P_2_K0s_K0s_from_fits_err*P_2_K0s_K0s_from_fits_err ) );

  float SysErrBackground_K0s_K0s_work_2 = ( fabs( P_2_K0s_K0s_2 - P_2_K0s_K0s_from_fits) - SysErrBackground_K0s_K0s_corr_2 )/fabs(P_2_K0s_K0s_2);

  float SysErrBackground_K0s_K0s_2 = 0;

  if( SysErrBackground_K0s_K0s_work_2 > 0 ) SysErrBackground_K0s_K0s_2 = SysErrBackground_K0s_K0s_work_2; //store sys. err. only if it is larger than statistical fluctuations

  //----------------------------

  //cuts variation
  float SysErr_tight_cuts_K0s_K0s_2 = 0;

  //float P_K0s_K0s_tight_topo_cuts = 0;
  //float P_K0s_K0s_tight_pT_cuts = 0;
/*
  if(cut_type == 0)
  {
    float P_K0s_K0s_tight_topo_cuts = fitK0s_K0s_tight_topo_cuts_Delta_phi[1]->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_K0s_K0s_tight_topo_cuts_err = fitK0s_K0s_tight_topo_cuts_Delta_phi[1]->GetParError(1)/(L0_alpha*L0_alpha);

    //statistical error correction for systematic error
    float SysErr_tight_topo_cuts_K0s_K0s_corr = sqrt( fabs( P_K0s_K0s_err_2*P_K0s_K0s_err_2 - P_K0s_K0s_tight_topo_cuts_err*P_K0s_K0s_tight_topo_cuts_err ) );

    float SysErr_tight_topo_cuts_K0s_K0s_work = ( fabs( P_K0s_K0s_2 - P_K0s_K0s_tight_topo_cuts) - SysErr_tight_topo_cuts_K0s_K0s_corr )/fabs(P_K0s_K0s_2);

    float SysErr_tight_topo_cuts_K0s_K0s = 0;

    if( SysErr_tight_topo_cuts_K0s_K0s_work > 0 ) SysErr_tight_topo_cuts_K0s_K0s = SysErr_tight_topo_cuts_K0s_K0s_work; //store sys. err. only if it is larger than statistical fluctuations


    float P_K0s_K0s_tight_pT_cuts = fitK0s_K0s_tight_pT_cuts_Delta_phi[1]->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_K0s_K0s_tight_pT_cuts_err = fitK0s_K0s_tight_pT_cuts_Delta_phi[1]->GetParError(1)/(L0_alpha*L0_alpha);

    //statistical error correction for systematic error
    float SysErr_tight_pT_cuts_K0s_K0s_corr = sqrt( fabs( P_K0s_K0s_err_2*P_K0s_K0s_err_2 - P_K0s_K0s_tight_pT_cuts_err*P_K0s_K0s_tight_pT_cuts_err ) );

    float SysErr_tight_pT_cuts_K0s_K0s_work = ( fabs( P_K0s_K0s_2 - P_K0s_K0s_tight_pT_cuts) - SysErr_tight_pT_cuts_K0s_K0s_corr )/fabs(P_K0s_K0s_2);

    float SysErr_tight_pT_cuts_K0s_K0s = 0;

    if( SysErr_tight_pT_cuts_K0s_K0s_work > 0 ) SysErr_tight_pT_cuts_K0s_K0s = SysErr_tight_pT_cuts_K0s_K0s_work; //store sys. err. only if it is larger than statistical fluctuations

    SysErr_tight_cuts_K0s_K0s = sqrt(SysErr_tight_topo_cuts_K0s_K0s*SysErr_tight_topo_cuts_K0s_K0s + SysErr_tight_pT_cuts_K0s_K0s*SysErr_tight_pT_cuts_K0s_K0s);
    //SysErr_tight_cuts_K0s_K0s = SysErr_tight_topo_cuts_K0s_K0s;
    //SysErr_tight_cuts_K0s_K0s = SysErr_tight_pT_cuts_K0s_K0s;
  }
*/
  //total
  float SysErrTot_K0s_K0s_2 = sqrt( sysErr_alpha_K0s_K0s*sysErr_alpha_K0s_K0s + SysErrSlope_K0s_K0s_2*SysErrSlope_K0s_K0s_2 + SysErrBackground_K0s_K0s_2*SysErrBackground_K0s_K0s_2  );
  //float SysErrTot_K0s_K0s_2 = sqrt( sysErr_alpha_K0s_K0s*sysErr_alpha_K0s_K0s + SysErrSlope_K0s_K0s_2*SysErrSlope_K0s_K0s_2 + SysErrBackground_K0s_K0s_2*SysErrBackground_K0s_K0s_2 + SysErr_tight_cuts_K0s_K0s_2*SysErr_tight_cuts_K0s_K0s_2 );




  TPaveText *K0s_K0s_text_more_2 = new TPaveText(0.47, 0.15, 0.8, 0.53, "NDC"); //without polarization value
  K0s_K0s_text_more_2->SetTextFont(42);
  //K0s_K0s_text_more_2->SetTextSize(15);
  //K0s_K0s_text_more_2->AddText("STAR");
  //K0s_K0s_text_more_2->AddText("STAR preliminary");
  //((TText*)K0s_K0s_text_more_2->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_text_more_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  K0s_K0s_text_more_2->AddText("Minimum bias");
  K0s_K0s_text_more_2->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_text_more_2->AddText("|#it{y}| < 1");
  K0s_K0s_text_more_2->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
  K0s_K0s_text_more_2->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
  K0s_K0s_text_more_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2}} = %.3f #pm %.3f", P_2_K0s_K0s_2, fabs(P_2_K0s_K0s_err_2)));
  //K0s_K0s_text_more_2->AddText(Form("P_{topo} = %.2f", P_K0s_K0s_tight_topo_cuts));
  //K0s_K0s_text_more_2->AddText(Form("P_{pT} = %.2f", P_K0s_K0s_tight_pT_cuts));
  K0s_K0s_text_more_2->SetFillColorAlpha(0, 0.01);
  K0s_K0s_text_more_2->Draw("same");

  K0s_K0s_2_leg->Draw("same");

  K0s_K0s_cosThetaProdPlane_2_can_2->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations_delta_phi_new/K0s_K0s_cosThetaProdPlane_subtract_delta_phi_new_%i.png", 1));

  PolarizationGraph_delta_phi[1]->SetPoint(1, P_2_K0s_K0s_2, 4);
  PolarizationGraph_delta_phi[1]->SetPointError(1, fabs(P_2_K0s_K0s_err_2), 0);

  PolarizationGraph_delta_phi_sys_err[1]->SetPoint(1, P_2_K0s_K0s_2, 4);
  PolarizationGraph_delta_phi_sys_err[1]->SetPointError(1, fabs(SysErrTot_K0s_K0s_2*P_2_K0s_K0s_2), 0.045);

  //____________________________________________________________________________________________________________________________________________________________________________________________________________


  //plot polarization graphs in bins

  TCanvas *PolarizationGraph_delta_phi_can = new TCanvas("PolarizationGraph_delta_phi_can", "PolarizationGraph_delta_phi_can", 2000, 1200);
  PolarizationGraph_delta_phi_can->Divide(2,1);

  TH1F *DefaultHist = new TH1F("DefaultHist", "DefaultHist", 10, -0.5, 0.5);
  DefaultHist->GetXaxis()->SetTitle("P_{#Lambda_{1}#Lambda_{2}}");
  DefaultHist->GetXaxis()->CenterTitle();
  DefaultHist->GetXaxis()->SetRangeUser(-0.22, 0.22);
  //DefaultHist->GetYaxis()->SetRangeUser(0.01, 4.99); //with K0s
  DefaultHist->GetYaxis()->SetRangeUser(3.01, 4.99); //without K0s
  //DefaultHist->GetYaxis()->SetNdivisions(505);
  //DefaultHist->GetYaxis()->SetNdivisions(06); //with K0s
  DefaultHist->GetYaxis()->SetNdivisions(00); //without K0s
  //DefaultHist->GetYaxis()->SetNdivisions(-5);
  DefaultHist->GetYaxis()->ChangeLabel(1,-1,-1,-1,-1,-1,"K_{s}^{0}K_{s}^{0}");
  //DefaultHist->GetYaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"#Lambda#Lambda");
  //DefaultHist->GetYaxis()->ChangeLabel(3,-1,-1,-1,-1,-1,"#bar{#Lambda}#bar{#Lambda}");
  //DefaultHist->GetYaxis()->ChangeLabel(4,-1,-1,-1,-1,-1,"K_{s}^{0}K_{s}^{0}");

  TLine *ZeroLine_eta = new TLine(0,3.01,0,4.99);
  ZeroLine_eta->SetLineStyle(9);
  ZeroLine_eta->SetLineColor(1);
  //ZeroLine_eta->Draw("same");

  //TPaveText *Polarization_text = new TPaveText(0.12, 0.58, 0.49, 0.89, "NDC"); //with K0s
  TPaveText *Polarization_text = new TPaveText(0.12, 0.55, 0.49, 0.89, "NDC"); //without K0s
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

  for(unsigned int bin = 1; bin < 3; bin++)
  {

    //plot polarization graphs in bins
    PolarizationGraph_delta_phi_can->cd(bin);

    //gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.05);

    DefaultHist->Draw();

    PolarizationGraph_delta_phi[bin-1]->SetMarkerStyle(20);
    PolarizationGraph_delta_phi[bin-1]->SetMarkerSize(2);
    PolarizationGraph_delta_phi[bin-1]->SetMarkerColor(kRed);
    PolarizationGraph_delta_phi[bin-1]->SetLineColor(kRed);
    PolarizationGraph_delta_phi[bin-1]->Write(Form("PolarizationGraph_delta_phi_%i", bin-1));
    PolarizationGraph_delta_phi[bin-1]->Draw("p e same");

    PolarizationGraph_delta_phi_sys_err[bin-1]->SetMarkerSize(2);
    PolarizationGraph_delta_phi_sys_err[bin-1]->SetFillColorAlpha(kRed, 0.25);
    PolarizationGraph_delta_phi_sys_err[bin-1]->SetLineColor(kRed);
    PolarizationGraph_delta_phi_sys_err[bin-1]->Write(Form("PolarizationGraph_delta_phi_sys_err_%i", bin-1));
    PolarizationGraph_delta_phi_sys_err[bin-1]->Draw("2 same");

    ZeroLine_eta->Draw("same");

    bin_text[bin-1] = new TPaveText(0.12, 0.85, 0.49, 0.89, "NDC");
    bin_text[bin-1]->SetTextFont(42);
    if(bin == 1) bin_text[bin-1]->AddText("|#Delta#phi| < #pi/4");
    if(bin == 2) bin_text[bin-1]->AddText("|#Delta#phi| > #pi/4");
    bin_text[bin-1]->SetFillColorAlpha(0, 0.01);
    bin_text[bin-1]->Draw("same");

    if(bin == 1 ) Polarization_text->Draw("same");
  }

  PolarizationGraph_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_polarization/K0s_polarization_delta_phi_new.png");


  //__________________________________________________________________________________________________


  //LLbarOutFile->Write();
  inFile->Close();
  out_file->Close();

  return true;
}
