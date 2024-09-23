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

//const double p_mass_PDG = 0.93827208816; //p mass on GeV/c^2 from latest PDG
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

bool cuts(int energy, float K0s_y, int hasTOFinfo)
{

  if( !( TMath::Abs(K0s_y) < K0s_y_cut ) ) return false;
  if( energy == 510 && hasTOFinfo == 0 ) return false; //TOF matched pion 1, for Run17 only now
  //if( strictTOF_cut == 2 && (pi1_hasTOFinfo == 0 || pi2_hasTOFinfo == 0) ) return false; //TOF matched pions and protons
  //if(cos(K0s_theta) < K0s_cos_theta_cut) return false;
  //if(K0s_decayL > K0s_decayK0s_cut) return false;

  return true;
}


//analyze invariant mass spectra
//arguments are vectors of bin numbers of the invariant mass peak
//the bins are determined via fit to the invariant mass spectra
//bool InvMass(TTree *K0s_tree, vector<int> &invMassBins_L, vector<int> &invMassBins_L0, vector<int> &invMassBins_L0bar, const int readMode)
bool Ana004_K0s_corr_2D_plot_Minv(const int cutType, const int energy = 510, const int year = 2017)
{

  TFile *InvMassFile;


  if(cutType == 0)
  {
    if(year == 2012)
    {
      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_K0s_2D_ana_cuts_orig_full_prod.root", year), "read");

      InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/2024_06_new_ME/InvMass_K0s_2D_ana_cuts_new_ME.root", year), "read");

      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_K0s_2D_ana_cuts_new_prod.root", year), "read");
    }
    else
    {
      InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_K0s_2D_ana_cuts.root", year), "read");
    }
  }

  if(cutType == 1)InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_K0s_2D_tight_topo_cuts.root", year), "read");
  if(cutType == 2)InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_K0s_2D_tight_pT_cut.root", year), "read");

  if( !(InvMassFile->IsOpen()) )
  {
    cout<<"Unable to open file with invariant mass histograms!"<<endl;
    return false;
  }

  //__________________________________________________________________________________________

  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass_US[nPtBins_corr][nPtBins_corr];
  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //background when US is paired with LS
  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass_LS[nPtBins_corr][nPtBins_corr];
  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass[nPtBins_corr][nPtBins_corr];


  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2));
      K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2));
      K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("K0s1_inv_mass_vs_K0s2_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2));
    }
  }


  //_____________________________________________________________________________________

  //projection Minv histograms
  TH1F *K0s1_inv_mass_K0sK0s_projection_US[nPtBins_corr][nPtBins_corr];
  TH1F *K0s1_inv_mass_K0sK0s_projection_LS[nPtBins_corr][nPtBins_corr];
  TH1F *K0s1_inv_mass_K0sK0s_projection[nPtBins_corr][nPtBins_corr];
  TH1F *K0s1_inv_mass_K0sK0s_projection_full[nPtBins_corr][nPtBins_corr];

  TH1F *K0s2_inv_mass_K0sK0s_projection_US[nPtBins_corr][nPtBins_corr];
  TH1F *K0s2_inv_mass_K0sK0s_projection_LS[nPtBins_corr][nPtBins_corr];
  TH1F *K0s2_inv_mass_K0sK0s_projection[nPtBins_corr][nPtBins_corr];
  TH1F *K0s2_inv_mass_K0sK0s_projection_full[nPtBins_corr][nPtBins_corr];


  float sideBandScale_K0s_K0s[nPtBins_corr][nPtBins_corr];

  float nK0s_pairs_US = 0;
  float nK0s_pairs_US_LS = 0;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


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


      //US
      TCanvas *K0s1_inv_mass_vs_K0s2_inv_mass_US_can = new TCanvas(Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_can_%i_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      K0s1_inv_mass_vs_K0s2_inv_mass_US_can->cd();

      K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      //K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      //K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *K0s_K0s_text_Minv_2D_US = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      K0s_K0s_text_Minv_2D_US->SetTextFont(42);
      //K0s_K0s_text_no_corr->AddText("STAR Internal");
      //K0s_K0s_text_no_corr->AddText("STAR preliminary");
      //((TText*)K0s_K0s_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      K0s_K0s_text_Minv_2D_US->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      K0s_K0s_text_Minv_2D_US->AddText("Minimum bias");
      K0s_K0s_text_Minv_2D_US->AddText("K^{0}_{s}-K^{0}_{s}");
      K0s_K0s_text_Minv_2D_US->AddText("US vs. US");
      K0s_K0s_text_Minv_2D_US->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_Minv_2D_US->Draw("same");

      TPaveText *K0s_K0s_text_Minv_2D_kine = new TPaveText(0.75, 0.83, 0.95, 0.98, "NDC");
      K0s_K0s_text_Minv_2D_kine->SetTextFont(42);
      K0s_K0s_text_Minv_2D_kine->AddText("|#it{y}| < 1");
      K0s_K0s_text_Minv_2D_kine->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      K0s_K0s_text_Minv_2D_kine->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      K0s_K0s_text_Minv_2D_kine->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_Minv_2D_kine->Draw("same");

      K0s1_inv_mass_vs_K0s2_inv_mass_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/Minv/US/K0s1_inv_mass_vs_K0s2_inv_mass_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      // US paired with LS for background estimation
      TCanvas *K0s1_inv_mass_vs_K0s2_inv_mass_US_LS_can = new TCanvas(Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      K0s1_inv_mass_vs_K0s2_inv_mass_US_LS_can->cd();

      K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      //K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      //K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *K0s_K0s_text_Minv_2D_US_LS = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      K0s_K0s_text_Minv_2D_US_LS->SetTextFont(42);
      //K0s_K0s_text_no_corr->AddText("STAR Internal");
      //K0s_K0s_text_no_corr->AddText("STAR preliminary");
      //((TText*)K0s_K0s_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      K0s_K0s_text_Minv_2D_US_LS->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      K0s_K0s_text_Minv_2D_US_LS->AddText("Minimum bias");
      K0s_K0s_text_Minv_2D_US_LS->AddText("K^{0}_{s}-K^{0}_{s}");
      K0s_K0s_text_Minv_2D_US_LS->AddText("US vs. LS");
      K0s_K0s_text_Minv_2D_US_LS->AddText("LS vs. US");
      K0s_K0s_text_Minv_2D_US_LS->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_Minv_2D_US_LS->Draw("same");

      K0s_K0s_text_Minv_2D_kine->Draw("same");

      K0s1_inv_mass_vs_K0s2_inv_mass_US_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/Minv/US/K0s1_inv_mass_vs_K0s2_inv_mass_US_LS_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //LS - pure continuum background
      TCanvas *K0s1_inv_mass_vs_K0s2_inv_mass_LS_can = new TCanvas(Form("K0s1_inv_mass_vs_K0s2_inv_mass_LS_can_%i_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_vs_K0s2_inv_mass_LS_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      K0s1_inv_mass_vs_K0s2_inv_mass_LS_can->cd();

      K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      //K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->GetXaxis()->SetRangeLSer(1.07, 1.2);
      K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      //K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->GetYaxis()->SetRangeLSer(1.07,1.2);
      K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *K0s_K0s_text_Minv_2D_LS = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      K0s_K0s_text_Minv_2D_LS->SetTextFont(42);
      //K0s_K0s_text_no_corr->AddText("STAR Internal");
      //K0s_K0s_text_no_corr->AddText("STAR preliminary");
      //((TText*)K0s_K0s_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      K0s_K0s_text_Minv_2D_LS->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      K0s_K0s_text_Minv_2D_LS->AddText("Minimum bias");
      K0s_K0s_text_Minv_2D_LS->AddText("K^{0}_{s}-K^{0}_{s}");
      K0s_K0s_text_Minv_2D_LS->AddText("LS vs. LS");
      K0s_K0s_text_Minv_2D_LS->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_Minv_2D_LS->Draw("same");

      K0s_K0s_text_Minv_2D_kine->Draw("same");

      K0s1_inv_mass_vs_K0s2_inv_mass_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/Minv/US/K0s1_inv_mass_vs_K0s2_inv_mass_LS_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *K0s1_inv_mass_K0sK0s_projection_US_can = new TCanvas(Form("K0s1_inv_mass_K0sK0s_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_K0sK0s_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      K0s1_inv_mass_K0sK0s_projection_US_can->cd();

      K0s1_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2] = (TH1F*)K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->ProjectionX(Form("K0s1_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      K0s1_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s1_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      K0s1_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      K0s1_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      K0s1_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s1_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->SetMarkerSize(2);
      K0s1_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->SetMarkerColor(kRed);
      K0s1_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->SetLineColor(kRed);
      K0s1_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->Draw("p e");

      //K0s1_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2] = (TH1F*)K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->ProjectionX(Form("K0s1_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      K0s1_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2] = (TH1F*)K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->ProjectionX(Form("K0s1_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      K0s1_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2]->SetMarkerStyle(24);
      K0s1_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2]->SetMarkerSize(2);
      K0s1_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      K0s1_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2]->SetLineColor(kBlue);
      K0s1_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2]->Draw("p e same");

      TPaveText *K0s_K0s_text_Minv_K0s1 = new TPaveText(0.6, 0.55, 0.85, 0.9, "NDC");
      K0s_K0s_text_Minv_K0s1->SetTextFont(42);
      //K0s_K0s_text_no_corr->AddText("STAR Internal");
      //K0s_K0s_text_no_corr->AddText("STAR preliminary");
      //((TText*)K0s_K0s_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      K0s_K0s_text_Minv_K0s1->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      K0s_K0s_text_Minv_K0s1->AddText("Minimum bias");
      K0s_K0s_text_Minv_K0s1->AddText("Projection to K^{0}_{s}(1)");
      //K0s_K0s_text_Minv_K0s1->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      K0s_K0s_text_Minv_K0s1->AddText("|#it{y}| < 1");
      K0s_K0s_text_Minv_K0s1->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      K0s_K0s_text_Minv_K0s1->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      K0s_K0s_text_Minv_K0s1->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_Minv_K0s1->Draw("same");

      TLegend *K0s_K0s_leg_Minv_K0s1 = new TLegend(0.6, 0.35, 0.85, 0.54);
      K0s_K0s_leg_Minv_K0s1->AddEntry(K0s1_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2], "US");
      K0s_K0s_leg_Minv_K0s1->AddEntry(K0s1_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2], "Conbinatorial bckg.");
      K0s_K0s_leg_Minv_K0s1->SetBorderSize(0);
      K0s_K0s_leg_Minv_K0s1->Draw("same");

      K0s1_inv_mass_K0sK0s_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/Minv/US/K0s1_inv_mass_K0sK0s_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      TCanvas *K0s2_inv_mass_K0sK0s_projection_US_can = new TCanvas(Form("K0s2_inv_mass_K0sK0s_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s2_inv_mass_K0sK0s_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      K0s2_inv_mass_K0sK0s_projection_US_can->cd();

      K0s2_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2] = (TH1F*)K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->ProjectionY(Form("K0s2_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      K0s2_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s2_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      K0s2_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      K0s2_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      K0s2_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s2_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->SetMarkerSize(2);
      K0s2_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->SetMarkerColor(kRed);
      K0s2_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->SetLineColor(kRed);
      K0s2_inv_mass_K0sK0s_projection_US[pTbin1][pTbin2]->Draw("p e");

      //K0s2_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2] = (TH1F*)K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->ProjectionY(Form("K0s2_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      K0s2_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2] = (TH1F*)K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->ProjectionY(Form("K0s2_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      K0s2_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2]->SetMarkerStyle(24);
      K0s2_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2]->SetMarkerSize(2);
      K0s2_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      K0s2_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2]->SetLineColor(kBlue);
      K0s2_inv_mass_K0sK0s_projection_LS[pTbin1][pTbin2]->Draw("p e same");

      TPaveText *K0s_K0s_text_Minv_K0s2 = new TPaveText(0.6, 0.55, 0.85, 0.9, "NDC");
      K0s_K0s_text_Minv_K0s2->SetTextFont(42);
      //K0s_K0s_text_no_corr->AddText("STAR Internal");
      //K0s_K0s_text_no_corr->AddText("STAR preliminary");
      //((TText*)K0s_K0s_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      K0s_K0s_text_Minv_K0s2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      K0s_K0s_text_Minv_K0s2->AddText("Minimum bias");
      K0s_K0s_text_Minv_K0s2->AddText("Projection to K^{0}_{s}(2)");
      //K0s_K0s_text_Minv_K0s2->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      K0s_K0s_text_Minv_K0s2->AddText("|#it{y}| < 1");
      K0s_K0s_text_Minv_K0s2->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      K0s_K0s_text_Minv_K0s2->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      K0s_K0s_text_Minv_K0s2->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_Minv_K0s2->Draw("same");

      K0s_K0s_leg_Minv_K0s1->Draw("same");

      K0s2_inv_mass_K0sK0s_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/Minv/US/K0s2_inv_mass_K0sK0s_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      //US-LS
      TF2 *doubleGauss_K0s_K0s = new TF2("doubleGauss_K0s_K0s", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 0.45, 0.55, 0.45, 0.55);
      //doubleGauss_K0s_K0s->SetParameters(1200, 0.499, 0.02, 0.499, 0.02);
      if( pTbin1 == 0 && pTbin2 == 0 ) doubleGauss_K0s_K0s->SetParameters(150, 0.499, 0.02, 0.499, 0.02);
      else if( fit_res_gaus[0][0]->IsValid() ) doubleGauss_K0s_K0s->SetParameters(fit_res_gaus[0][0]->Parameter(0), fit_res_gaus[0][0]->Parameter(1), fit_res_gaus[0][0]->Parameter(2), fit_res_gaus[0][0]->Parameter(3), fit_res_gaus[0][0]->Parameter(4));
      else
      {
        cout<<"Fit not valid for L-L pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
        return false;
      }

      TCanvas *K0s1_inv_mass_vs_K0s2_inv_mass_can = new TCanvas(Form("K0s1_inv_mass_vs_K0s2_inv_mass_can_%i_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_vs_K0s2_inv_mass_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      K0s1_inv_mass_vs_K0s2_inv_mass_can->cd();
/*
      //background already scaled from main macro
      //scale LS to match US
      sideBandScale_K0s_K0s[pTbin1][pTbin2] = K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh)/K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh);

      K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->Sumw2();
      K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->Scale(sideBandScale_K0s_K0s[pTbin1][pTbin2]);
*/
      K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2] = (TH2F*)K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->Clone(Form("K0s1_inv_mass_vs_K0s2_inv_mass_pT1_%i_pT2_%i", pTbin1, pTbin2));
      K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      //K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07,1.2);
      K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      //K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      //K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->Add(K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2], -1);
      K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->Add(K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2], -1);
      fit_res_gaus[pTbin1][pTbin2] = K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_K0s_K0s, "s 0", "", 0.49, 0.51);
      K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *K0s_K0s_text_Minv_2D_sig = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      K0s_K0s_text_Minv_2D_sig->SetTextFont(42);
      //K0s_K0s_text_no_corr->AddText("STAR Internal");
      //K0s_K0s_text_no_corr->AddText("STAR preliminary");
      //((TText*)K0s_K0s_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      K0s_K0s_text_Minv_2D_sig->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      K0s_K0s_text_Minv_2D_sig->AddText("Minimum bias");
      K0s_K0s_text_Minv_2D_sig->AddText("K^{0}_{s}-K^{0}_{s}");
      K0s_K0s_text_Minv_2D_sig->AddText("US - Bckg.");
      K0s_K0s_text_Minv_2D_sig->SetFillColorAlpha(0, 0.01);
      K0s_K0s_text_Minv_2D_sig->Draw("same");

      K0s_K0s_text_Minv_2D_kine->Draw("same");

      K0s1_inv_mass_vs_K0s2_inv_mass_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/Minv/US-LS/K0s1_inv_mass_vs_K0s2_inv_mass_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      if( !fit_res_gaus[pTbin1][pTbin2]->IsValid() )
      {
        cout<<"Fit not valid for pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;

        return false;
      }


      TCanvas *K0s1_inv_mass_K0sK0s_projection_can = new TCanvas(Form("K0s1_inv_mass_K0sK0s_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_K0sK0s_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      K0s1_inv_mass_K0sK0s_projection_can->cd();

      int K0sK0s_projectionBin_L1 = K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus[pTbin1][pTbin2]->Parameter(1));

      //K0s_inv_mass_K0sK0s_projection[pTbin1][pTbin2] = (TH1F*)K0s_inv_mass_vs_K0sbar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("K0s_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      K0s1_inv_mass_K0sK0s_projection[pTbin1][pTbin2] = (TH1F*)K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("K0s1_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), K0sK0s_projectionBin_L1, K0sK0s_projectionBin_L1);
      K0s1_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s1_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      K0s1_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      K0s1_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      K0s1_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s1_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      K0s1_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      K0s1_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      K0s1_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->Draw("p e");


      TF12 *K0sK0s_Gaus_K0s = new TF12("K0sK0s_Gaus_K0s", doubleGauss_K0s_K0s, fit_res_gaus[pTbin1][pTbin2]->Parameter(1), "x");
      K0sK0s_Gaus_K0s->SetLineColor(1);
      K0sK0s_Gaus_K0s->Draw("same");

      K0s_K0s_text_Minv_K0s1->Draw("same");

      TLegend *K0s_K0s_leg_Minv_L_sig = new TLegend(0.6, 0.45, 0.85, 0.54);
      K0s_K0s_leg_Minv_L_sig->AddEntry(K0s1_inv_mass_K0sK0s_projection[pTbin1][pTbin2], "US - Bckg.");
      K0s_K0s_leg_Minv_L_sig->AddEntry(K0sK0s_Gaus_K0s, "Gaussian fit");
      K0s_K0s_leg_Minv_L_sig->SetBorderSize(0);
      K0s_K0s_leg_Minv_L_sig->Draw("same");

      K0s1_inv_mass_K0sK0s_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/Minv/US-LS/K0s1_inv_mass_K0sK0s_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *K0s2_inv_mass_K0sK0s_projection_can = new TCanvas(Form("K0s2_inv_mass_K0sK0s_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s2_inv_mass_K0sK0s_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      K0s2_inv_mass_K0sK0s_projection_can->cd();

      int K0sK0s_projectionBin_L2 = K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus[pTbin1][pTbin2]->Parameter(3));

      //K0s2_inv_mass_K0sK0s_projection[pTbin1][pTbin2] = (TH1F*)K0s2_inv_mass_vs_K0s2bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("K0s2_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      K0s2_inv_mass_K0sK0s_projection[pTbin1][pTbin2] = (TH1F*)K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("K0s2_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), K0sK0s_projectionBin_L2, K0sK0s_projectionBin_L2);
      K0s2_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s2_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      K0s2_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      K0s2_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      K0s2_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s2_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      K0s2_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      K0s2_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      K0s2_inv_mass_K0sK0s_projection[pTbin1][pTbin2]->Draw("p e");


      TF12 *K0sK0s_Gaus_K0s2 = new TF12("K0sK0s_Gaus_K0s2", doubleGauss_K0s_K0s, fit_res_gaus[pTbin1][pTbin2]->Parameter(3), "y");
      K0sK0s_Gaus_K0s2->SetLineColor(1);
      K0sK0s_Gaus_K0s2->Draw("same");

      K0s_K0s_text_Minv_K0s2->Draw("same");

      K0s_K0s_leg_Minv_L_sig->Draw("same");

      K0s2_inv_mass_K0sK0s_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/Minv/US-LS/K0s2_inv_mass_K0sK0s_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //________________________________________________________________________________________________________

      TCanvas *K0s1_inv_mass_K0sK0s_projection_full_can = new TCanvas(Form("K0s1_inv_mass_K0sK0s_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_K0sK0s_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      K0s1_inv_mass_K0sK0s_projection_full_can->cd();

      int K0sK0s_projection_fullBin_L1 = K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus[pTbin1][pTbin2]->Parameter(1));

      K0s1_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2] = (TH1F*)K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("K0s_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      //K0s1_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2] = (TH1F*)K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("K0s1_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      K0s1_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s1_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      K0s1_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      K0s1_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      K0s1_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s1_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->SetMarkerSize(2);
      K0s1_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->SetMarkerColor(kRed);
      K0s1_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->SetLineColor(kRed);
      K0s1_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->Draw("p e");



      K0s_K0s_text_Minv_K0s1->Draw("same");

      TLegend *K0s_K0s_leg_Minv_L_sig_full = new TLegend(0.6, 0.45, 0.85, 0.54);
      K0s_K0s_leg_Minv_L_sig_full->AddEntry(K0s1_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2], "US - Bckg.");
      K0s_K0s_leg_Minv_L_sig_full->SetBorderSize(0);
      K0s_K0s_leg_Minv_L_sig_full->Draw("same");

      K0s1_inv_mass_K0sK0s_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/Minv/US-LS/K0s1_inv_mass_K0sK0s_projection_full_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *K0s2_inv_mass_K0sK0s_projection_full_can = new TCanvas(Form("K0s2_inv_mass_K0sK0s_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s2_inv_mass_K0sK0s_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      K0s2_inv_mass_K0sK0s_projection_full_can->cd();

      int K0sK0s_projection_fullBin_L2 = K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus[pTbin1][pTbin2]->Parameter(3));

      K0s2_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2] = (TH1F*)K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("K0s2_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      //K0s2_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2] = (TH1F*)K0s1_inv_mass_vs_K0s2_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("K0s2_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      K0s2_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#pi^{-}}");
      K0s2_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      K0s2_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      K0s2_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      K0s2_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s2_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->SetMarkerSize(2);
      K0s2_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->SetMarkerColor(kRed);
      K0s2_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->SetLineColor(kRed);
      K0s2_inv_mass_K0sK0s_projection_full[pTbin1][pTbin2]->Draw("p e");

      K0s_K0s_text_Minv_K0s2->Draw("same");

      K0s_K0s_leg_Minv_L_sig_full->Draw("same");

      K0s2_inv_mass_K0sK0s_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/Minv/US-LS/K0s2_inv_mass_K0sK0s_projection_full_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //________________________________________________________________________________________________________

      float K0s1_peak_mean = fit_res_gaus[pTbin1][pTbin2]->Parameter(1);
      float K0s1_peak_sigma = fabs(fit_res_gaus[pTbin1][pTbin2]->Parameter(2));

      float K0s2_peak_mean = fit_res_gaus[pTbin1][pTbin2]->Parameter(3);
      float K0s2_peak_sigma = fabs(fit_res_gaus[pTbin1][pTbin2]->Parameter(4));

      int bin_K0s_1_low = K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(K0s1_peak_mean-3*K0s1_peak_sigma);
      int bin_K0s_1_high = K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(K0s1_peak_mean+3*K0s1_peak_sigma);

      int bin_K0s_2_low = K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->FindBin(K0s2_peak_mean-3*K0s2_peak_sigma);
      int bin_K0s_2_high = K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->FindBin(K0s2_peak_mean+3*K0s2_peak_sigma);


      nK0s_pairs_US += K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->Integral(bin_K0s_1_low, bin_K0s_1_high, bin_K0s_2_low, bin_K0s_2_high);
      nK0s_pairs_US_LS += K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->Integral(bin_K0s_1_low, bin_K0s_1_high, bin_K0s_2_low, bin_K0s_2_high);




    }//end for eta

  }//end for pT

  cout<<"Number US K0s pairs from Minv "<<nK0s_pairs_US<<endl;
  cout<<"Number US_LS K0s pairs from Minv "<<nK0s_pairs_US_LS<<endl;


  InvMassFile->Close();



  return true;

}
