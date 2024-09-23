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

bool Ana004_K0s_corr_2D_get_corr_Delta_y_Delta_phi(const int cut_type = 0, const int energy = 510, const int year = 2017)
{
  //analyze stored Lambda pairs and save cos(theta*) histograms

  const float L0_alpha = 0.732; //decay parameter of K0s
  const float L0_alpha_relat_err = 0.014/L0_alpha; //relative error of decay parameter

  const float L0bar_alpha = -0.758; //decay paramteter of K0sbar
  const float L0bar_alpha_relat_err = 0.012/fabs(L0bar_alpha); //relative error of decay paramteter


  //_______________________________________________________________________________________________________________________________________________

  //systematic uncertainties
  //residual effect from PYTHIA closure test
  TFile *SysErrSlopeFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErrSlope_K0s_nocorr.root", year), "read");

  TH1F *SysErrSlope_delta_eta_delta_phi_hist[2];
  TH1F *ResidualPolarization_delta_eta_delta_phi_hist[2];

  for( unsigned int delta_eta_bin = 0; delta_eta_bin < 2; delta_eta_bin++ )
  {
    SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin] = (TH1F*)SysErrSlopeFile->Get(Form("SysErrSlope_delta_eta_delta_phi_hist_%i", delta_eta_bin));
    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin] = (TH1F*)SysErrSlopeFile->Get(Form("ResidualPolarization_K0s_delta_eta_delta_phi_hist_%i", delta_eta_bin));
  }

  //----------------------------------------------------------------------

  //systematic uncertainty histograms and values
  //alpha
  float sysErr_alpha_K0s_K0sbar = sqrt(L0_alpha_relat_err*L0_alpha_relat_err + L0bar_alpha_relat_err*L0bar_alpha_relat_err);
  float sysErr_alpha_K0s_K0s = sqrt(L0_alpha_relat_err*L0_alpha_relat_err + L0_alpha_relat_err*L0_alpha_relat_err);
  float sysErr_alpha_K0sbar_K0sbar = sqrt(L0bar_alpha_relat_err*L0bar_alpha_relat_err + L0bar_alpha_relat_err*L0bar_alpha_relat_err);

  //----------------------------------------------------------------------

  //cuts variation
  //have to run this code with cuts_type = 1 and 2 first
  TFile *SysErrCutsTopo;
  TFile *SysErrCutsPt;

  if(cut_type == 0 )
  {
    SysErrCutsTopo = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_K0s_tight_topo_cuts_Delta_y_Delta_phi_work.root", year), "read");

    SysErrCutsPt = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_K0s_tight_pT_cuts_Delta_y_Delta_phi_work.root", year), "read");
  }
  else
  {
    if( cut_type == 1 ) SysErrCutsTopo = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_K0s_tight_topo_cuts_Delta_y_Delta_phi_work.root", year), "recreate");
    if( cut_type == 2 ) SysErrCutsPt = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_K0s_tight_pT_cuts_Delta_y_Delta_phi_work.root", year), "recreate");
  }


  //tight cuts
  //topological cuts
  TF1 *fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US[2];
  TF1 *fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS[2];
  TF1 *fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi[2];


  //daughter pT cuts
  TF1 *fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US[2];
  TF1 *fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS[2];
  TF1 *fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi[2];


  for( unsigned int delta_eta_bin = 0; delta_eta_bin < 2; delta_eta_bin++)
  {
    if(cut_type == 0 )
    {
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin));
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin));
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_%i", delta_eta_bin));


      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = (TF1*)SysErrCutsPt->Get(Form("fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin));
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = (TF1*)SysErrCutsPt->Get(Form("fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin));
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin] = (TF1*)SysErrCutsPt->Get(Form("fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_%i", delta_eta_bin));
    }
    else
    {
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = new TF1(Form("fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin]->SetParameters(100, 0.5);

      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = new TF1(Form("fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin]->SetParameters(100, 0.5);

      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin] = new TF1(Form("fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);

      //-----------------------------------------------------------------------------------------

      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = new TF1(Form("fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin]->SetParameters(100, 0.5);

      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = new TF1(Form("fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin]->SetParameters(100, 0.5);

      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin] = new TF1(Form("fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);
    }

  }

  //------------------------------------------------------------------------------------------------------------------

  //for weighted average of uncorrected errors (absolute errors)

  float SysErr_sum_from_fits_topo_cuts = 0;
  float SysErr_sum_from_fits_of_w_topo_cuts = 0;

  float SysErr_sum_from_fits_pT_cut = 0;
  float SysErr_sum_from_fits_of_w_pT_cut = 0;

  //----------------------------------

  float SysErr_sum_topo_cuts = 0;
  float SysErr_sum_of_w_topo_cuts = 0;

  float SysErr_sum_pT_cut = 0;
  float SysErr_sum_of_w_pT_cut = 0;

  float SysErr_sum_background = 0;
  float SysErr_sum_of_w_background = 0;

  float SysErr_sum_ME = 0;
  float SysErr_sum_of_w_ME = 0;

  //-----------------------------------------------------------------------------------------------------------------

  //output file with polarization graphs
  TFile *out_file = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/Polarization/%i/Polarization_K0s_Delta_y_Delta_phi.root", year), "recreate");

  //----------------------------------------------------------------------------

  TFile *inFile; //output file to store production plane histograms;

  if(cut_type == 0)
  {
    if(year == 2012)
    {
      //inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_ana_cuts_work.root", year), "read");

      inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/2024_08_final/ProdPlane_K0s_2D_ana_cuts_2_sigma_Delta_phi_third.root", year), "read");

      //inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/2024_08_fixed_ME_LLbar/ProdPlane_K0s_2D_ana_cuts_3_sigma_Delta_phi_third.root", year), "read");

      //inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/2024_06_new_ME/ProdPlane_K0s_2D_ana_cuts_new_ME.root", year), "read");

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
    //inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_tight_topo_cuts.root", year), "read");

    inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/2024_08_final/ProdPlane_K0s_2D_tight_topo_cuts_2_sigma_Delta_phi_third.root", year), "read");
  }
  else if(cut_type == 2) //create production plane file from nTuple - run in this mode first
  {
    //inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_tight_pT_cut.root", year), "read");

    inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/2024_08_final/ProdPlane_K0s_2D_tight_pT_cut_2_sigma_Delta_phi_third.root", year), "read");
  }
  else  //read file to store wrong and good sign M_inv distributions and daughter pT
  {
    cout<<"Wrong cut type!"<<endl;
    return false;
  }



  //data histograms

  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_hist");

  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist");

  //--------------------------------------------------------

  //mixed event
  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist");

  TH2F *K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist = (TH2F*)inFile->Get("K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist");

  //________________________________________________________________________________________



  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(42);

  //polarization graph for Delta y
  //first graph is for |Delta y| < 0.5, second is for 0.5 < |Delta y| < 2.0
  TGraphErrors *PolarizationGraph_delta_eta_delta_phi[2];
  TGraphErrors *PolarizationGraph_delta_eta_delta_phi_sys_err[2];
  TGraphErrors *PolarizationGraph_delta_eta_delta_phi_sys_err_average[2];


  TCanvas *PolarizationGraph_delta_eta_delta_phi_can = new TCanvas("PolarizationGraph_delta_eta_delta_phi_can", "PolarizationGraph_delta_eta_delta_phi_can", 2000, 1200);
  PolarizationGraph_delta_eta_delta_phi_can->Divide(2,1);

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

  TPaveText *Delta_eta_text_1 = new TPaveText(0.15, 0.11, 0.35, 0.2, "NDC");
  Delta_eta_text_1->SetTextFont(43);
  Delta_eta_text_1->SetTextSize(33);
  Delta_eta_text_1->SetFillColorAlpha(0, 0.01);
  Delta_eta_text_1->AddText("|#Delta#it{y}| < 0.5, |#Delta#phi| < #pi/3");

  TPaveText *Delta_eta_text_2 = new TPaveText(0.15, 0.11, 0.45, 0.2, "NDC");
  Delta_eta_text_2->SetTextFont(43);
  Delta_eta_text_2->SetTextSize(33);
  Delta_eta_text_2->SetFillColorAlpha(0, 0.01);
  Delta_eta_text_2->AddText("0.5 < |#Delta#it{y}| < 2.0 or |#Delta#phi| > #pi/3");

  out_file->cd();

  for( unsigned int delta_eta_bin = 1; delta_eta_bin < 3; delta_eta_bin++)
  {
    //create polarization graphs, defined earlier

    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1] = new TGraphErrors(1);
    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1] = new TGraphErrors(1);
    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1] = new TGraphErrors(1);

    //____________________________________________________________________________________________________________________________________________________________________________________________________________

    TCanvas *K0s_K0s_cosThetaProdPlane_eta_no_corr_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_eta_no_corr_can_%i", delta_eta_bin), Form("K0s_K0s_cosThetaProdPlane_eta_no_corr_can_%i", delta_eta_bin), 1200, 1000);

    K0s_K0s_cosThetaProdPlane_eta_no_corr_can->cd();

    TH1D *K0s_K0s_cosThetaProdPlane_eta_US_hist = K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_hist->ProjectionX( Form("proj_US_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);

    K0s_K0s_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_eta_US_hist->GetXaxis()->CenterTitle();
    //K0s_K0s_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTextSizePixels(30);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    K0s_K0s_cosThetaProdPlane_eta_US_hist->GetYaxis()->CenterTitle();
    //K0s_K0s_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTextSizePixels(30);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetMaxDigits(3);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->SetMarkerSize(1.5);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->SetMarkerStyle(20);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->SetMarkerColor(kRed);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->SetLineColor(kRed);
    double nK0sK0s = K0s_K0s_cosThetaProdPlane_eta_US_hist->Integral();
    K0s_K0s_cosThetaProdPlane_eta_US_hist->Sumw2();
    //K0s_K0s_cosThetaProdPlane_eta_US_hist->Divide(K0s_K0s_cosThetaProdPlane_eta_eff);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->Scale(1./K0s_K0s_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1));
    //K0s_K0s_cosThetaProdPlane_eta_US_hist->Scale(1./K0s_K0s_cosThetaProdPlane_eta_US_hist->Integral());
    K0s_K0s_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //ME here just for plotting, used loser
    TH1D *K0s_K0s_cosThetaProdPlane_eta_ME_hist = K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->ProjectionX( Form("proj_US_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);

    K0s_K0s_cosThetaProdPlane_eta_ME_hist->SetMarkerSize(1.5);
    K0s_K0s_cosThetaProdPlane_eta_ME_hist->SetMarkerStyle(24);
    K0s_K0s_cosThetaProdPlane_eta_ME_hist->SetMarkerColor(1);
    K0s_K0s_cosThetaProdPlane_eta_ME_hist->SetLineColor(1);
    K0s_K0s_cosThetaProdPlane_eta_ME_hist->Sumw2();
    K0s_K0s_cosThetaProdPlane_eta_ME_hist->Scale(nK0sK0s/K0s_K0s_cosThetaProdPlane_eta_ME_hist->Integral()); //scale ME to US
    K0s_K0s_cosThetaProdPlane_eta_ME_hist->Scale(1./K0s_K0s_cosThetaProdPlane_eta_ME_hist->GetXaxis()->GetBinWidth(1)); //US is scaled by bin width, ME needs to be scaled as well
    K0s_K0s_cosThetaProdPlane_eta_ME_hist->SetMinimum(0);
    K0s_K0s_cosThetaProdPlane_eta_ME_hist->Draw("p e same");

    TF1 *fitK0s_K0s_US_ThetaStar_no_corr_ME = new TF1("fitK0s_K0s_US_ThetaStar_no_corr_ME", "[0]*(1 + [1]*x)", -1, 1);
    fitK0s_K0s_US_ThetaStar_no_corr_ME->SetParameters(100, 0.5);

    K0s_K0s_cosThetaProdPlane_eta_ME_hist->Fit(fitK0s_K0s_US_ThetaStar_no_corr_ME, "s i 0 r");

    float P_K0s_K0s_no_corr_ME = fitK0s_K0s_US_ThetaStar_no_corr_ME->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_K0s_K0s_no_corr_ME_err = fitK0s_K0s_US_ThetaStar_no_corr_ME->GetParError(1)/(L0_alpha*L0_alpha);


    TH1D *K0s_K0s_cosThetaProdPlane_eta_LS_hist = K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist->ProjectionX( Form("proj_LS_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);

    K0s_K0s_cosThetaProdPlane_eta_LS_hist->SetMarkerSize(1.5);
    K0s_K0s_cosThetaProdPlane_eta_LS_hist->SetMarkerStyle(21);
    K0s_K0s_cosThetaProdPlane_eta_LS_hist->SetMarkerColor(kBlue);
    double nK0sK0s_back = K0s_K0s_cosThetaProdPlane_eta_LS_hist->Integral();
    K0s_K0s_cosThetaProdPlane_eta_LS_hist->Sumw2();
    K0s_K0s_cosThetaProdPlane_eta_LS_hist->Scale(1./K0s_K0s_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1));
    //K0s_K0s_cosThetaProdPlane_eta_LS_hist->Divide(K0s_K0s_cosThetaProdPlane_eta_eff);
    K0s_K0s_cosThetaProdPlane_eta_LS_hist->Draw("p e same");


    TH1D *K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist = K0s_K0s_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->ProjectionX( Form("proj_LS_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);

    K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerSize(1.5);
    K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerStyle(25);
    K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerColor(kMagenta+1);
    K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist->SetLineColor(kMagenta+1);
    K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist->Sumw2();
    K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist->Scale(nK0sK0s_back/K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist->Integral()); //scale ME_LS to background
    K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist->Scale(1./K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist->GetXaxis()->GetBinWidth(1)); //background is scaled by bin width, ME_LS needs to be scaled as well
    K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist->SetMinimum(0);
    K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist->Draw("p e same");


    TF1 *fitK0s_K0s_US_ThetaStar_no_corr = new TF1("fitK0s_K0s_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
    fitK0s_K0s_US_ThetaStar_no_corr->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->Fit(fitK0s_K0s_US_ThetaStar_no_corr, "s i 0 r");

    float P_K0s_K0s_no_corr = fitK0s_K0s_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_K0s_K0s_no_corr_err = fitK0s_K0s_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0_alpha);

    fitK0s_K0s_US_ThetaStar_no_corr->SetLineColor(1);
    //fitK0s_K0s_US_ThetaStar_no_corr->Draw("same");

    TLegend *K0s_K0s_leg = new TLegend(0.15, 0.45, 0.45, 0.69);
    K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_eta_US_hist, "(US-US) p#pi");
    K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_eta_ME_hist, "(US-US) ME");
    K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_eta_LS_hist, "Combinatorial bckg.");
    K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist, "Bckg. ME");
    //K0s_K0s_leg->AddEntry(fitK0s_K0s_US_ThetaStar_no_corr, "Linear fit to US");
    K0s_K0s_leg->SetBorderSize(0);
    K0s_K0s_leg->SetFillColorAlpha(0, 0.01);
    K0s_K0s_leg->Draw("same");

    TPaveText *K0s_K0s_text_no_corr = new TPaveText(0.5, 0.4, 0.85, 0.75, "NDC");
    K0s_K0s_text_no_corr->SetTextFont(42);
    //K0s_K0s_text_no_corr->AddText("STAR Internal");
    //K0s_K0s_text_no_corr->AddText("STAR preliminary");
    //((TText*)K0s_K0s_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
    K0s_K0s_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    K0s_K0s_text_no_corr->AddText("Minimum bias, no correction");
    K0s_K0s_text_no_corr->AddText("K_{s}^{0}K_{s}^{0}");
    K0s_K0s_text_no_corr->AddText("|#it{y}| < 1");
    K0s_K0s_text_no_corr->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
    K0s_K0s_text_no_corr->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
    K0s_K0s_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_K0s_K0s_no_corr, fabs(P_K0s_K0s_no_corr_err)));
    K0s_K0s_text_no_corr->AddText(Form("P_{ME} = %.3f #pm %.3f", P_K0s_K0s_no_corr_ME, fabs(P_K0s_K0s_no_corr_ME_err)));
    K0s_K0s_text_no_corr->SetFillColorAlpha(0, 0.01);
    K0s_K0s_text_no_corr->Draw("same");

    if(delta_eta_bin == 1) Delta_eta_text_1->Draw("same");
    if(delta_eta_bin == 2) Delta_eta_text_2->Draw("same");

    K0s_K0s_cosThetaProdPlane_eta_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations_delta_y_delta_phi/K0s_K0s_cosThetaProdPlane_eta_no_corr_less_delta_eta_delta_phi_%i.png", delta_eta_bin));

    //----------------------------------------------------------------------------------------------------

    TCanvas *K0s_K0s_cosThetaProdPlane_eta_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_eta_can_%i", delta_eta_bin), Form("K0s_K0s_cosThetaProdPlane_eta_can_%i", delta_eta_bin), 1200, 1000);

    K0s_K0s_cosThetaProdPlane_eta_can->cd();

    //ME histogram higher

    K0s_K0s_cosThetaProdPlane_eta_US_hist->Divide(K0s_K0s_cosThetaProdPlane_eta_ME_hist); // correct US using ME
    K0s_K0s_cosThetaProdPlane_eta_US_hist->Scale(nK0sK0s/K0s_K0s_cosThetaProdPlane_eta_US_hist->Integral()); //scale back to raw US
    K0s_K0s_cosThetaProdPlane_eta_US_hist->Scale(1./K0s_K0s_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    K0s_K0s_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //K0s_K0s_cosThetaProdPlane_eta_ME_hist->Draw("same p e");

    K0s_K0s_cosThetaProdPlane_eta_LS_hist->Divide(K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist); //correct background using ME
    K0s_K0s_cosThetaProdPlane_eta_LS_hist->Scale(nK0sK0s_back/K0s_K0s_cosThetaProdPlane_eta_LS_hist->Integral()); //scale back to raw background
    K0s_K0s_cosThetaProdPlane_eta_LS_hist->Scale(1./K0s_K0s_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    K0s_K0s_cosThetaProdPlane_eta_LS_hist->SetMinimum(0);
    K0s_K0s_cosThetaProdPlane_eta_LS_hist->Draw("p e same");

    //fit dN/dcos(theta*)
    //signal + bacground
    TF1 *fitK0s_K0s_US_ThetaStar = new TF1("fitK0s_K0s_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitK0s_K0s_US_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->Fit(fitK0s_K0s_US_ThetaStar, "s i 0 r");

    float P_K0s_K0s = fitK0s_K0s_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_K0s_K0s_err = fitK0s_K0s_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

    fitK0s_K0s_US_ThetaStar->SetLineColor(1);
    fitK0s_K0s_US_ThetaStar->Draw("same");

    //background
    TF1 *fitK0s_K0s_US_LS_ThetaStar = new TF1("fitK0s_K0s_US_LS_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitK0s_K0s_US_LS_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaUS_LS_wrong_sign = L_inv_mass_US_LS->Fit(fitGaUS_LSsBack, "s i 0", "", 1.07, 1.4);
    K0s_K0s_cosThetaProdPlane_eta_LS_hist->Fit(fitK0s_K0s_US_LS_ThetaStar, "s i 0 r");

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
    K0s_K0s_text->AddText("K_{s}^{0}K_{s}^{0}");
    K0s_K0s_text->AddText("|#it{y}| < 1");
    K0s_K0s_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
    K0s_K0s_text->AddText(Form("P_{tot} = %.2f #pm %.2f", P_K0s_K0s, fabs(P_K0s_K0s_err)));
    K0s_K0s_text->AddText(Form("P_{bckg} = %.2f #pm %.2f", P_K0s_K0s_back, fabs(P_K0s_K0s_back_err)));
    K0s_K0s_text->SetFillColorAlpha(0, 0.01);
    K0s_K0s_text->Draw("same");

    K0s_K0s_leg->Draw("same");

    if(delta_eta_bin == 1) Delta_eta_text_1->Draw("same");
    if(delta_eta_bin == 2) Delta_eta_text_2->Draw("same");

    K0s_K0s_cosThetaProdPlane_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations_delta_y_delta_phi/K0s_K0s_cosThetaProdPlane_eta_less_delta_eta_delta_phi_%i.png", delta_eta_bin));

    //----------------------------------------------------------------------

    TCanvas *K0s_K0s_cosThetaProdPlane_eta_can_2 = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_eta_can_2_%i", delta_eta_bin), Form("K0s_K0s_cosThetaProdPlane_eta_can_2_%i", delta_eta_bin), 1200, 1000);

    K0s_K0s_cosThetaProdPlane_eta_can_2->cd();

    K0s_K0s_cosThetaProdPlane_eta_US_hist->Add(K0s_K0s_cosThetaProdPlane_eta_LS_hist, -1);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->Write(Form("K0s_K0s_cosThetaProdPlane_delta_y_%i", delta_eta_bin-1));
    K0s_K0s_cosThetaProdPlane_eta_US_hist->Draw("p e");


    TF1 *fitK0s_K0s_US_ThetaStar_2 = new TF1("fitK0s_K0s_US_ThetaStar_2", "[0]*(1 + [1]*x)", -1, 1);
    fitK0s_K0s_US_ThetaStar_2->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    K0s_K0s_cosThetaProdPlane_eta_US_hist->Fit(fitK0s_K0s_US_ThetaStar_2, "s i 0 r");


    //store fit result for systematic errors
    //tight topo cuts
    if(cut_type == 1)
    {
      SysErrCutsTopo->cd();

      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParameters(fitK0s_K0s_US_ThetaStar->GetParameter(0), fitK0s_K0s_US_ThetaStar->GetParameter(1));
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(0, fitK0s_K0s_US_ThetaStar->GetParError(0));
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(1, fitK0s_K0s_US_ThetaStar->GetParError(1));
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->Write();


      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParameters(fitK0s_K0s_US_LS_ThetaStar->GetParameter(0), fitK0s_K0s_US_LS_ThetaStar->GetParameter(1));
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(0, fitK0s_K0s_US_LS_ThetaStar->GetParError(0));
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(1, fitK0s_K0s_US_LS_ThetaStar->GetParError(1));
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->Write();

      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParameters(fitK0s_K0s_US_ThetaStar_2->GetParameter(0), fitK0s_K0s_US_ThetaStar_2->GetParameter(1));
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(0, fitK0s_K0s_US_ThetaStar_2->GetParError(0));
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(1, fitK0s_K0s_US_ThetaStar_2->GetParError(1));
      fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->Write();
    }

    //tight pT cuts
    if(cut_type == 2)
    {
      SysErrCutsPt->cd();

      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParameters(fitK0s_K0s_US_ThetaStar->GetParameter(0), fitK0s_K0s_US_ThetaStar->GetParameter(1));
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(0, fitK0s_K0s_US_ThetaStar->GetParError(0));
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(1, fitK0s_K0s_US_ThetaStar->GetParError(1));
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->Write();


      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParameters(fitK0s_K0s_US_LS_ThetaStar->GetParameter(0), fitK0s_K0s_US_LS_ThetaStar->GetParameter(1));
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(0, fitK0s_K0s_US_LS_ThetaStar->GetParError(0));
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(1, fitK0s_K0s_US_LS_ThetaStar->GetParError(1));
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->Write();

      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParameters(fitK0s_K0s_US_ThetaStar_2->GetParameter(0), fitK0s_K0s_US_ThetaStar_2->GetParameter(1));
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(0, fitK0s_K0s_US_ThetaStar_2->GetParError(0));
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(1, fitK0s_K0s_US_ThetaStar_2->GetParError(1));
      fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->Write();
    }



    float P_K0s_K0s_2 = fitK0s_K0s_US_ThetaStar_2->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_K0s_K0s_err_2 = fitK0s_K0s_US_ThetaStar_2->GetParError(1)/(L0_alpha*L0_alpha);

    fitK0s_K0s_US_ThetaStar_2->SetLineColor(1);
    fitK0s_K0s_US_ThetaStar_2->Draw("same");


    //TPaveText *K0s_K0s_text_2 = new TPaveText(0.47, 0.15, 0.85, 0.53, "NDC"); //with polarization value
    TPaveText *K0s_K0s_text_2 = new TPaveText(0.47, 0.15, 0.8, 0.53, "NDC"); //without polarization value
    K0s_K0s_text_2->SetTextFont(42);
    //K0s_K0s_text_2->SetTextSize(15);
    //K0s_K0s_text_2->AddText("STAR");
    //K0s_K0s_text_2->AddText("STAR preliminary");
    //((TText*)K0s_K0s_text_2->GetListOfLines()->Last())->SetTextColor(2);
    K0s_K0s_text_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    K0s_K0s_text_2->AddText("Minimum bias");
    K0s_K0s_text_2->AddText("K_{s}^{0}K_{s}^{0}");
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
    K0s_K0s_2_leg->AddEntry(K0s_K0s_cosThetaProdPlane_eta_US_hist, "(US-US)-Bckg.");
    K0s_K0s_2_leg->AddEntry(fitK0s_K0s_US_ThetaStar_2, "Fit", "l");
    K0s_K0s_2_leg->SetBorderSize(0);
    K0s_K0s_2_leg->SetFillColorAlpha(0, 0.01);
    K0s_K0s_2_leg->Draw("same");

    if(delta_eta_bin == 1) Delta_eta_text_1->Draw("same");
    if(delta_eta_bin == 2) Delta_eta_text_2->Draw("same");

    K0s_K0s_cosThetaProdPlane_eta_can_2->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations_delta_y_delta_phi/K0s_K0s_cosThetaProdPlane_eta_subtract_less_delta_eta_delta_phi_%i.png", delta_eta_bin));


    //calculate total systematic uncertainty for L-Lbar

    //slope
    //float SysErrSlope_K0s_K0s = SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinContent(2);
    //float SysErrSlope_K0s_K0s = SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinContent(2);

    float SysErrSlope_K0s_K0s = fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinContent(1))/fabs(P_K0s_K0s_2);

    SysErr_sum_ME += fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinContent(1))/fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinError(1));
    SysErr_sum_of_w_ME += 1./fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinError(1));


    //----------------------------

    //background subtraction
    //float P_K0s_K0s_from_fits = P_K0s_K0s - nK0sK0s_back/nK0sK0s*P_K0s_K0s_back;
    //float P_K0s_K0s_from_fits_err = sqrt( P_K0s_K0s_err*P_K0s_K0s_err + nK0sK0s_back/nK0sK0s*nK0sK0s_back/nK0sK0s*P_K0s_K0s_back_err*P_K0s_K0s_back_err );

    float nK0sK0s_fit = fitK0s_K0s_US_ThetaStar->GetParameter(0);
    float nK0sK0s_fit_err = fitK0s_K0s_US_ThetaStar->GetParError(0);

    float nK0sK0s_back_fit = fitK0s_K0s_US_LS_ThetaStar->GetParameter(0);
    float nK0sK0s_back_fit_err = fitK0s_K0s_US_LS_ThetaStar->GetParError(0);

    float P_K0s_K0s_from_fits = P_K0s_K0s + nK0sK0s_back_fit/(nK0sK0s_fit-nK0sK0s_back_fit)*(P_K0s_K0s - P_K0s_K0s_back);

    float P_K0s_K0s_from_fits_err = sqrt( P_K0s_K0s_err*P_K0s_K0s_err + nK0sK0s_back_fit/(nK0sK0s_fit-nK0sK0s_back_fit)*nK0sK0s_back_fit/(nK0sK0s_fit-nK0sK0s_back_fit)*P_K0s_K0s_err*P_K0s_K0s_err +
                                           nK0sK0s_back_fit/(nK0sK0s_fit-nK0sK0s_back_fit)*nK0sK0s_back_fit/(nK0sK0s_fit-nK0sK0s_back_fit)*P_K0s_K0s_back_err*P_K0s_K0s_back_err +
                                           nK0sK0s_back_fit*nK0sK0s_back_fit/(nK0sK0s_fit-nK0sK0s_back_fit)/(nK0sK0s_fit-nK0sK0s_back_fit)/(nK0sK0s_fit-nK0sK0s_back_fit)/(nK0sK0s_fit-nK0sK0s_back_fit)*(P_K0s_K0s - P_K0s_K0s_back)*(P_K0s_K0s - P_K0s_K0s_back)*nK0sK0s_fit_err*nK0sK0s_fit_err +
                                           (P_K0s_K0s - P_K0s_K0s_back)*(P_K0s_K0s - P_K0s_K0s_back)/(nK0sK0s_fit-nK0sK0s_back_fit)/(nK0sK0s_fit-nK0sK0s_back_fit)/(nK0sK0s_fit-nK0sK0s_back_fit)/(nK0sK0s_fit-nK0sK0s_back_fit)*nK0sK0s_back_fit_err*nK0sK0s_back_fit_err );



    cout<<"Polarization from fit: "<<P_K0s_K0s_from_fits<<" +- "<<P_K0s_K0s_from_fits_err<<endl;
    cout<<endl;

    //statistical error correction for systematic error
    float SysErrBackground_K0s_K0s_corr = sqrt( fabs( P_K0s_K0s_err_2*P_K0s_K0s_err_2 - P_K0s_K0s_from_fits_err*P_K0s_K0s_from_fits_err ) );

    float SysErrBackground_K0s_K0s_work = ( fabs( P_K0s_K0s_2 - P_K0s_K0s_from_fits) - SysErrBackground_K0s_K0s_corr )/fabs(P_K0s_K0s_2);

    float SysErrBackground_K0s_K0s = 0;

    if( SysErrBackground_K0s_K0s_work > 0 ) SysErrBackground_K0s_K0s = SysErrBackground_K0s_K0s_work; //store sys. err. only if it is larger than statistical fluctuations

    //SysErr_sum_background += SysErrBackground_K0s_K0s*fabs(P_K0s_K0s_2)/SysErrBackground_K0s_K0s_corr;
    SysErr_sum_background += SysErrBackground_K0s_K0s*fabs(P_K0s_K0s_from_fits)/SysErrBackground_K0s_K0s_corr;
    SysErr_sum_of_w_background += 1./SysErrBackground_K0s_K0s_corr;

    //----------------------------------------------------------------------------------------------------------

    //cuts variation
    float SysErr_tight_cuts_K0s_K0s_from_fits = 0;

    float SysErr_tight_topo_cuts_K0s_K0s_from_fits = 0;
    float SysErr_tight_pT_cuts_K0s_K0s_from_fits = 0;


    float SysErr_tight_cuts_K0s_K0s_from_fits_no_corr = 0;

    float SysErr_tight_topo_cuts_K0s_K0s_from_fits_no_corr = 0;
    float SysErr_tight_pT_cuts_K0s_K0s_from_fits_no_corr = 0;

    //------------------------------------------------------

    float SysErr_tight_cuts_K0s_K0s = 0;

    float SysErr_tight_topo_cuts_K0s_K0s = 0;
    float SysErr_tight_pT_cuts_K0s_K0s = 0;


    float SysErr_tight_cuts_K0s_K0s_no_corr = 0;

    float SysErr_tight_topo_cuts_K0s_K0s_no_corr = 0;
    float SysErr_tight_pT_cuts_K0s_K0s_no_corr = 0;

    if(cut_type == 0)
    {
      //independent fit of US-US and background

      float P_K0s_K0s_tight_topo_cuts_from_fits = fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_K0s_K0s_tight_topo_cuts_from_fits_err = fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->GetParError(1)/(L0_alpha*L0_alpha);

      float nK0sK0s_tight_topo_cuts_fit = fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->GetParameter(0);
      float nK0sK0s_tight_topo_cuts_fit_err = fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->GetParError(0);


      float P_K0s_K0s_tight_topo_cuts_from_fits_back = fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_K0s_K0s_tight_topo_cuts_from_fits_back_err = fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->GetParError(1)/(L0_alpha*L0_alpha);

      float nK0sK0s_tight_topo_cuts_back_fit = fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->GetParameter(0);
      float nK0sK0s_tight_topo_cuts_back_fit_err = fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->GetParError(0);


      float P_K0s_K0s_tight_topo_cuts_from_fits_signal = P_K0s_K0s_tight_topo_cuts_from_fits + nK0sK0s_tight_topo_cuts_back_fit/(nK0sK0s_tight_topo_cuts_fit-nK0sK0s_tight_topo_cuts_back_fit)*(P_K0s_K0s_tight_topo_cuts_from_fits - P_K0s_K0s_tight_topo_cuts_from_fits_back);

      float P_K0s_K0s_tight_topo_cuts_from_fits_signal_err = sqrt( P_K0s_K0s_tight_topo_cuts_from_fits_err*P_K0s_K0s_tight_topo_cuts_from_fits_err + nK0sK0s_tight_topo_cuts_back_fit/(nK0sK0s_tight_topo_cuts_fit-nK0sK0s_tight_topo_cuts_back_fit)*nK0sK0s_tight_topo_cuts_back_fit/(nK0sK0s_tight_topo_cuts_fit-nK0sK0s_tight_topo_cuts_back_fit)*P_K0s_K0s_tight_topo_cuts_from_fits_err*P_K0s_K0s_tight_topo_cuts_from_fits_err +
                                             nK0sK0s_tight_topo_cuts_back_fit/(nK0sK0s_tight_topo_cuts_fit-nK0sK0s_tight_topo_cuts_back_fit)*nK0sK0s_tight_topo_cuts_back_fit/(nK0sK0s_tight_topo_cuts_fit-nK0sK0s_tight_topo_cuts_back_fit)*P_K0s_K0s_tight_topo_cuts_from_fits_back_err*P_K0s_K0s_tight_topo_cuts_from_fits_back_err +
                                             nK0sK0s_tight_topo_cuts_back_fit*nK0sK0s_tight_topo_cuts_back_fit/(nK0sK0s_tight_topo_cuts_fit-nK0sK0s_tight_topo_cuts_back_fit)/(nK0sK0s_tight_topo_cuts_fit-nK0sK0s_tight_topo_cuts_back_fit)/(nK0sK0s_tight_topo_cuts_fit-nK0sK0s_tight_topo_cuts_back_fit)/(nK0sK0s_tight_topo_cuts_fit-nK0sK0s_tight_topo_cuts_back_fit)*(P_K0s_K0s_tight_topo_cuts_from_fits - P_K0s_K0s_tight_topo_cuts_from_fits_back)*(P_K0s_K0s_tight_topo_cuts_from_fits - P_K0s_K0s_tight_topo_cuts_from_fits_back)*nK0sK0s_tight_topo_cuts_fit_err*nK0sK0s_tight_topo_cuts_fit_err +
                                             (P_K0s_K0s_tight_topo_cuts_from_fits - P_K0s_K0s_tight_topo_cuts_from_fits_back)*(P_K0s_K0s_tight_topo_cuts_from_fits - P_K0s_K0s_tight_topo_cuts_from_fits_back)/(nK0sK0s_tight_topo_cuts_fit-nK0sK0s_tight_topo_cuts_back_fit)/(nK0sK0s_tight_topo_cuts_fit-nK0sK0s_tight_topo_cuts_back_fit)/(nK0sK0s_tight_topo_cuts_fit-nK0sK0s_tight_topo_cuts_back_fit)/(nK0sK0s_tight_topo_cuts_fit-nK0sK0s_tight_topo_cuts_back_fit)*nK0sK0s_tight_topo_cuts_back_fit_err*nK0sK0s_tight_topo_cuts_back_fit_err );

      //statistical error correction for systematic error
      float SysErr_tight_topo_cuts_K0s_K0s_from_fits_corr = sqrt( fabs( P_K0s_K0s_from_fits_err*P_K0s_K0s_from_fits_err - P_K0s_K0s_tight_topo_cuts_from_fits_signal_err*P_K0s_K0s_tight_topo_cuts_from_fits_signal_err ) );

      float SysErr_tight_topo_cuts_K0s_K0s_from_fits_work = ( fabs( P_K0s_K0s_from_fits - P_K0s_K0s_tight_topo_cuts_from_fits_signal) - SysErr_tight_topo_cuts_K0s_K0s_from_fits_corr )/fabs(P_K0s_K0s_from_fits);

      if( SysErr_tight_topo_cuts_K0s_K0s_from_fits_work > 0 ) SysErr_tight_topo_cuts_K0s_K0s_from_fits = SysErr_tight_topo_cuts_K0s_K0s_from_fits_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_topo_cuts_K0s_K0s_from_fits_no_corr = fabs( P_K0s_K0s_from_fits - P_K0s_K0s_tight_topo_cuts_from_fits_signal)/fabs(P_K0s_K0s_from_fits);


      SysErr_sum_from_fits_topo_cuts += SysErr_tight_topo_cuts_K0s_K0s_from_fits*fabs(P_K0s_K0s_from_fits)/SysErr_tight_topo_cuts_K0s_K0s_from_fits_corr;
      SysErr_sum_from_fits_of_w_topo_cuts += 1./SysErr_tight_topo_cuts_K0s_K0s_from_fits_corr;

      //----------------------------------------------------------------

      //fit after background subtraction

      float P_K0s_K0s_tight_topo_cuts = fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_K0s_K0s_tight_topo_cuts_err = fitK0s_K0s_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParError(1)/(L0_alpha*L0_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_topo_cuts_K0s_K0s_corr = sqrt( fabs( P_K0s_K0s_err_2*P_K0s_K0s_err_2 - P_K0s_K0s_tight_topo_cuts_err*P_K0s_K0s_tight_topo_cuts_err ) );

      float SysErr_tight_topo_cuts_K0s_K0s_work = ( fabs( P_K0s_K0s_2 - P_K0s_K0s_tight_topo_cuts) - SysErr_tight_topo_cuts_K0s_K0s_corr )/fabs(P_K0s_K0s_2);

      if( SysErr_tight_topo_cuts_K0s_K0s_work > 0 ) SysErr_tight_topo_cuts_K0s_K0s = SysErr_tight_topo_cuts_K0s_K0s_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_topo_cuts_K0s_K0s_no_corr = fabs( P_K0s_K0s_2 - P_K0s_K0s_tight_topo_cuts)/fabs(P_K0s_K0s_2);


      SysErr_sum_topo_cuts += SysErr_tight_topo_cuts_K0s_K0s*fabs(P_K0s_K0s_2)/SysErr_tight_topo_cuts_K0s_K0s_corr;
      SysErr_sum_of_w_topo_cuts += 1./SysErr_tight_topo_cuts_K0s_K0s_corr;


      //___________________________________________________________________________________________________________________________________

      //independent fit of US-US and background

      float P_K0s_K0s_tight_pT_cuts_from_fits = fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_K0s_K0s_tight_pT_cuts_from_fits_err = fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->GetParError(1)/(L0_alpha*L0_alpha);

      float nK0sK0s_tight_pT_cuts_fit = fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->GetParameter(0);
      float nK0sK0s_tight_pT_cuts_fit_err = fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->GetParError(0);


      float P_K0s_K0s_tight_pT_cuts_from_fits_back = fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_K0s_K0s_tight_pT_cuts_from_fits_back_err = fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->GetParError(1)/(L0_alpha*L0_alpha);

      float nK0sK0s_tight_pT_cuts_back_fit = fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->GetParameter(0);
      float nK0sK0s_tight_pT_cuts_back_fit_err = fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->GetParError(0);


      float P_K0s_K0s_tight_pT_cuts_from_fits_signal = P_K0s_K0s_tight_pT_cuts_from_fits + nK0sK0s_tight_pT_cuts_back_fit/(nK0sK0s_tight_pT_cuts_fit-nK0sK0s_tight_pT_cuts_back_fit)*(P_K0s_K0s_tight_pT_cuts_from_fits - P_K0s_K0s_tight_pT_cuts_from_fits_back);

      float P_K0s_K0s_tight_pT_cuts_from_fits_signal_err = sqrt( P_K0s_K0s_tight_pT_cuts_from_fits_err*P_K0s_K0s_tight_pT_cuts_from_fits_err + nK0sK0s_tight_pT_cuts_back_fit/(nK0sK0s_tight_pT_cuts_fit-nK0sK0s_tight_pT_cuts_back_fit)*nK0sK0s_tight_pT_cuts_back_fit/(nK0sK0s_tight_pT_cuts_fit-nK0sK0s_tight_pT_cuts_back_fit)*P_K0s_K0s_tight_pT_cuts_from_fits_err*P_K0s_K0s_tight_pT_cuts_from_fits_err +
                                             nK0sK0s_tight_pT_cuts_back_fit/(nK0sK0s_tight_pT_cuts_fit-nK0sK0s_tight_pT_cuts_back_fit)*nK0sK0s_tight_pT_cuts_back_fit/(nK0sK0s_tight_pT_cuts_fit-nK0sK0s_tight_pT_cuts_back_fit)*P_K0s_K0s_tight_pT_cuts_from_fits_back_err*P_K0s_K0s_tight_pT_cuts_from_fits_back_err +
                                             nK0sK0s_tight_pT_cuts_back_fit*nK0sK0s_tight_pT_cuts_back_fit/(nK0sK0s_tight_pT_cuts_fit-nK0sK0s_tight_pT_cuts_back_fit)/(nK0sK0s_tight_pT_cuts_fit-nK0sK0s_tight_pT_cuts_back_fit)/(nK0sK0s_tight_pT_cuts_fit-nK0sK0s_tight_pT_cuts_back_fit)/(nK0sK0s_tight_pT_cuts_fit-nK0sK0s_tight_pT_cuts_back_fit)*(P_K0s_K0s_tight_pT_cuts_from_fits - P_K0s_K0s_tight_pT_cuts_from_fits_back)*(P_K0s_K0s_tight_pT_cuts_from_fits - P_K0s_K0s_tight_pT_cuts_from_fits_back)*nK0sK0s_tight_pT_cuts_fit_err*nK0sK0s_tight_pT_cuts_fit_err +
                                             (P_K0s_K0s_tight_pT_cuts_from_fits - P_K0s_K0s_tight_pT_cuts_from_fits_back)*(P_K0s_K0s_tight_pT_cuts_from_fits - P_K0s_K0s_tight_pT_cuts_from_fits_back)/(nK0sK0s_tight_pT_cuts_fit-nK0sK0s_tight_pT_cuts_back_fit)/(nK0sK0s_tight_pT_cuts_fit-nK0sK0s_tight_pT_cuts_back_fit)/(nK0sK0s_tight_pT_cuts_fit-nK0sK0s_tight_pT_cuts_back_fit)/(nK0sK0s_tight_pT_cuts_fit-nK0sK0s_tight_pT_cuts_back_fit)*nK0sK0s_tight_pT_cuts_back_fit_err*nK0sK0s_tight_pT_cuts_back_fit_err );

      //statistical error correction for systematic error
      float SysErr_tight_pT_cuts_K0s_K0s_from_fits_corr = sqrt( fabs( P_K0s_K0s_from_fits_err*P_K0s_K0s_from_fits_err - P_K0s_K0s_tight_pT_cuts_from_fits_signal_err*P_K0s_K0s_tight_pT_cuts_from_fits_signal_err ) );

      float SysErr_tight_pT_cuts_K0s_K0s_from_fits_work = ( fabs( P_K0s_K0s_from_fits - P_K0s_K0s_tight_pT_cuts_from_fits_signal) - SysErr_tight_pT_cuts_K0s_K0s_from_fits_corr )/fabs(P_K0s_K0s_from_fits);

      if( SysErr_tight_pT_cuts_K0s_K0s_from_fits_work > 0 ) SysErr_tight_pT_cuts_K0s_K0s_from_fits = SysErr_tight_pT_cuts_K0s_K0s_from_fits_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_pT_cuts_K0s_K0s_from_fits_no_corr = fabs( P_K0s_K0s_from_fits - P_K0s_K0s_tight_pT_cuts_from_fits_signal)/fabs(P_K0s_K0s_from_fits);


      SysErr_sum_from_fits_pT_cut += SysErr_tight_pT_cuts_K0s_K0s_from_fits*fabs(P_K0s_K0s_from_fits)/SysErr_tight_pT_cuts_K0s_K0s_from_fits_corr;
      SysErr_sum_from_fits_of_w_pT_cut += 1./SysErr_tight_pT_cuts_K0s_K0s_from_fits_corr;

      //----------------------------------------------------------------


      float P_K0s_K0s_tight_pT_cuts = fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_K0s_K0s_tight_pT_cuts_err = fitK0s_K0s_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParError(1)/(L0_alpha*L0_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_pT_cuts_K0s_K0s_corr = sqrt( fabs( P_K0s_K0s_err_2*P_K0s_K0s_err_2 - P_K0s_K0s_tight_pT_cuts_err*P_K0s_K0s_tight_pT_cuts_err ) );

      float SysErr_tight_pT_cuts_K0s_K0s_work = ( fabs( P_K0s_K0s_2 - P_K0s_K0s_tight_pT_cuts) - SysErr_tight_pT_cuts_K0s_K0s_corr )/fabs(P_K0s_K0s_2);

      if( SysErr_tight_pT_cuts_K0s_K0s_work > 0 ) SysErr_tight_pT_cuts_K0s_K0s = SysErr_tight_pT_cuts_K0s_K0s_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_pT_cuts_K0s_K0s_no_corr = fabs( P_K0s_K0s_2 - P_K0s_K0s_tight_pT_cuts)/fabs(P_K0s_K0s_2);


      SysErr_sum_pT_cut += SysErr_tight_pT_cuts_K0s_K0s*fabs(P_K0s_K0s_2)/SysErr_tight_pT_cuts_K0s_K0s_corr;
      SysErr_sum_of_w_pT_cut += 1./SysErr_tight_pT_cuts_K0s_K0s_corr;

      //___________________________________________________________________________________________________________________________________

      SysErr_tight_cuts_K0s_K0s = sqrt(SysErr_tight_topo_cuts_K0s_K0s*SysErr_tight_topo_cuts_K0s_K0s + SysErr_tight_pT_cuts_K0s_K0s*SysErr_tight_pT_cuts_K0s_K0s);

      SysErr_tight_cuts_K0s_K0s_no_corr = sqrt(SysErr_tight_topo_cuts_K0s_K0s_no_corr*SysErr_tight_topo_cuts_K0s_K0s_no_corr + SysErr_tight_pT_cuts_K0s_K0s_no_corr*SysErr_tight_pT_cuts_K0s_K0s_no_corr);

      //----------------------------

      SysErr_tight_cuts_K0s_K0s_from_fits = sqrt(SysErr_tight_topo_cuts_K0s_K0s_from_fits*SysErr_tight_topo_cuts_K0s_K0s_from_fits + SysErr_tight_pT_cuts_K0s_K0s_from_fits*SysErr_tight_pT_cuts_K0s_K0s_from_fits);

      SysErr_tight_cuts_K0s_K0s_from_fits_no_corr = sqrt(SysErr_tight_topo_cuts_K0s_K0s_from_fits_no_corr*SysErr_tight_topo_cuts_K0s_K0s_from_fits_no_corr + SysErr_tight_pT_cuts_K0s_K0s_from_fits_no_corr*SysErr_tight_pT_cuts_K0s_K0s_from_fits_no_corr);


    }

    //--------------------------------------------------------------------------------------------------------------------------

    //total
    //float SysErrTot_K0s_K0s = sqrt( sysErr_alpha_K0s_K0s*sysErr_alpha_K0s_K0s + SysErrSlope_K0s_K0s*SysErrSlope_K0s_K0s + SysErrBackground_K0s_K0s*SysErrBackground_K0s_K0s  );
    //float SysErrTot_K0s_K0s = sqrt( sysErr_alpha_K0s_K0s*sysErr_alpha_K0s_K0s + SysErrSlope_K0s_K0s*SysErrSlope_K0s_K0s + SysErrBackground_K0s_K0s*SysErrBackground_K0s_K0s + SysErr_tight_cuts_K0s_K0s*SysErr_tight_cuts_K0s_K0s);
    float SysErrTot_K0s_K0s = sqrt( sysErr_alpha_K0s_K0s*sysErr_alpha_K0s_K0s + SysErrSlope_K0s_K0s*SysErrSlope_K0s_K0s + SysErrBackground_K0s_K0s*SysErrBackground_K0s_K0s + SysErr_tight_cuts_K0s_K0s_from_fits*SysErr_tight_cuts_K0s_K0s_from_fits);

    //-----------------------------------------------------------------------------------------------------------

    //PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPoint(1, P_K0s_K0s_2, 4);
    //PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPointError(1, fabs(P_K0s_K0s_err_2), 0);

    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPoint(1, P_K0s_K0s_from_fits, 4);
    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPointError(1, fabs(P_K0s_K0s_from_fits_err), 0);


    //PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPoint(1, P_K0s_K0s_2, 4);
    //PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPointError(1, fabs(SysErrTot_K0s_K0s*P_K0s_K0s_2), 0.045);

    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPoint(1, P_K0s_K0s_from_fits, 4);
    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPointError(1, fabs(SysErrTot_K0s_K0s*P_K0s_K0s_from_fits), 0.045);


    //____________________________________________________________________________________________________________________________________________________________________________________________________________

    //save polarization graphs
    out_file->cd();

    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->Write(Form("PolarizationGraph_delta_y_delta_phi_%i", delta_eta_bin-1));
    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->Write(Form("PolarizationGraph_delta_y_delta_phi_sys_err_%i", delta_eta_bin-1));

    //____________________________________________________________________________________________________________________________________________________________________________________________________________

  }

  float SysErr_topo_weights = SysErr_sum_topo_cuts/SysErr_sum_of_w_topo_cuts;
  float SysErr_pT_weights = SysErr_sum_pT_cut/SysErr_sum_of_w_pT_cut;
  float SysErr_bckg_weights = SysErr_sum_background/SysErr_sum_of_w_background;
  float SysErr_ME_weights = SysErr_sum_ME/SysErr_sum_of_w_ME;

  cout<<"Systematic uncetainty weighted averages by source:"<<endl;
  cout<<"s_topo | s_pT | s_bckg | s_ME"<<endl;
  cout<<SysErr_topo_weights<<" | "<<SysErr_pT_weights<<" | "<<SysErr_bckg_weights<<" | "<<SysErr_ME_weights<<endl;
  cout<<endl;

  float SysErrTot_weights_work = sqrt(SysErr_topo_weights*SysErr_topo_weights + SysErr_pT_weights*SysErr_pT_weights + SysErr_bckg_weights*SysErr_bckg_weights + SysErr_ME_weights*SysErr_ME_weights);

  cout<<"Total uncertainty from weighter average:"<<endl;
  cout<<SysErrTot_weights_work<<endl;
  cout<<endl;

  //-----------------------------------------------------------------------------------------------------------

  float SysErr_topo_weights_from_fits = SysErr_sum_from_fits_topo_cuts/SysErr_sum_from_fits_of_w_topo_cuts;
  float SysErr_pT_weights_from_fits = SysErr_sum_from_fits_pT_cut/SysErr_sum_from_fits_of_w_pT_cut;

  cout<<"Systematic uncetainty weighted averages by source:"<<endl;
  cout<<"s_topo | s_pT | s_bckg | s_ME"<<endl;
  cout<<SysErr_topo_weights_from_fits<<" | "<<SysErr_pT_weights_from_fits<<" | "<<SysErr_bckg_weights<<" | "<<SysErr_ME_weights<<endl;
  cout<<endl;

  float SysErrTot_weights_from_fits_work = sqrt(SysErr_topo_weights_from_fits*SysErr_topo_weights_from_fits + SysErr_pT_weights_from_fits*SysErr_pT_weights_from_fits + SysErr_bckg_weights*SysErr_bckg_weights + SysErr_ME_weights*SysErr_ME_weights);

  cout<<"Total uncertainty from weighter average:"<<endl;
  cout<<SysErrTot_weights_from_fits_work<<endl;
  cout<<endl;


  for( unsigned int delta_eta_bin = 1; delta_eta_bin < 3; delta_eta_bin++)
  {

    //float K0sK0s_SysErrTot_weights_work_relat = SysErrTot_weights_work/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(1);
    float K0sK0s_SysErrTot_weights_work_relat = SysErrTot_weights_from_fits_work/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(1);

    float K0sK0s_SysErrTot_weights = sqrt( K0sK0s_SysErrTot_weights_work_relat*K0sK0s_SysErrTot_weights_work_relat + sysErr_alpha_K0s_K0s*sysErr_alpha_K0s_K0s);

    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetPoint(1, PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(1), 4);
    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetPointError(1, K0sK0s_SysErrTot_weights*PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(1), 0.045);

    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->Write(Form("PolarizationGraph_delta_y_delta_phi_sys_err_average_%i", delta_eta_bin-1));

    //-----------------------------------------------------------------------------


    //plot polarization graphs in bins
    PolarizationGraph_delta_eta_delta_phi_can->cd(delta_eta_bin);

    //gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.05);

    DefaultHist->Draw();

    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetMarkerStyle(20);
    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetMarkerSize(2);
    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetMarkerColor(kRed);
    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetLineColor(kRed);
    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->Draw("p e same");

    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetMarkerSize(2);
    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetFillColorAlpha(kRed, 0.25);
    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetLineColor(kRed);
    //PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->Draw("2 same");

    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetMarkerSize(2);
    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetFillColorAlpha(kRed, 0.25);
    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetLineColor(kRed);
    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->Draw("2 same");

    ZeroLine_eta->Draw("same");

    if(delta_eta_bin == 1 ) Polarization_text->Draw("same");

  }


  PolarizationGraph_delta_eta_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_polarization/K0s_polarization_delta_eta_delta_phi.png");

  //__________________________________________________________________________________________________


  inFile->Close();

  return true;
}
