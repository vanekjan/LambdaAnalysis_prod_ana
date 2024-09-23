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
bool Ana004_K0s_corr_2D_get_corr(const int cutType, const int energy = 510, const int year = 2017)
{

  //for extraction of reference polarization from K0s baseline
  const float L0_alpha = 0.732; //decay parameter of L0
  const float L0_alpha_relat_err = 0.014/L0_alpha; //relative error of decay parameter

  const float L0bar_alpha = -0.758; //decay paramteter of L0bar
  const float L0bar_alpha_relat_err = 0.012/fabs(L0bar_alpha); //relative error of decay paramteter

  //_________________________________________________________________________

  //systematic uncertainties
  //residual effect from PYTHIA closure test
  //TFile *SysErrSlopeFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErrSlope.root", year), "read");

  //cuts variation
  //have to run this code with cuts_type = 1 and 2 first
  TFile *SysErrCutsTopo;
  TFile *SysErrCutsPt;

  if(cutType == 0 )
  {
    SysErrCutsTopo = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_topo_cuts_K0s.root", year), "read");
    SysErrCutsPt = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_pT_cuts_K0s.root", year), "read");
  }
  else
  {
    if( cutType == 1 ) SysErrCutsTopo = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_topo_cuts_K0s.root", year), "recreate");
    if( cutType == 2 ) SysErrCutsPt = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/SysErr/%i/SysErr_tight_pT_cuts_K0s.root", year), "recreate");
  }

  //_______________________________________________________________________________________________________________________________________________

  //systematic uncertainty histograms and values
  //alpha
  //float sysErr_alpha_L0_L0bar = sqrt(L0_alpha_relat_err*L0_alpha_relat_err + L0bar_alpha_relat_err*L0bar_alpha_relat_err);
  float sysErr_alpha_K0s_K0s = sqrt(L0_alpha_relat_err*L0_alpha_relat_err + L0_alpha_relat_err*L0_alpha_relat_err);
  //float sysErr_alpha_L0bar_L0bar = sqrt(L0bar_alpha_relat_err*L0bar_alpha_relat_err + L0bar_alpha_relat_err*L0bar_alpha_relat_err);

  //slope
  //TH1F *SysErrSlope_hist = (TH1F*)SysErrSlopeFile->Get("SysErrSlope_hist");

  //tight cuts
  TF1 *fitK0s_K0s_tight_topo_cuts;

  TF1 *fitK0s_K0s_tight_pT_cuts;

  if(cutType == 0 )
  {
    fitK0s_K0s_tight_topo_cuts = (TF1*)SysErrCutsTopo->Get("fitK0s_K0s_tight_topo_cuts");

    fitK0s_K0s_tight_pT_cuts = (TF1*)SysErrCutsPt->Get("fitK0s_K0s_tight_pT_cuts");
  }
  else
  {
    fitK0s_K0s_tight_topo_cuts = new TF1("fitK0s_K0s_tight_topo_cuts", "[0]*(1 + [1]*x)", -1, 1);
    fitK0s_K0s_tight_topo_cuts->SetParameters(100, 0.5);


    fitK0s_K0s_tight_pT_cuts = new TF1("fitK0s_K0s_tight_pT_cuts", "[0]*(1 + [1]*x)", -1, 1);
    fitK0s_K0s_tight_pT_cuts->SetParameters(100, 0.5);
  }

  //______________________________________________________________________________________________________

  TFile *OutFile_pol = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/Polarization_K0s.root", year), "recreate");


  TFile *inFile; //output file to store production plane histograms

  if(cutType == 0) //create production plane file from nTuple - run in this mode first
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

  }
  else if(cutType == 1) //create production plane file from nTuple - run in this mode first
  {
    inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_tight_topo_cuts.root", year), "read");
  }
  else if(cutType == 2) //create production plane file from nTuple - run in this mode first
  {
    inFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_K0s_2D_tight_pT_cut.root", year), "read");
  }
  else  //read file to store wrong and good sign M_inv distributions and daughter pT
  {
    cout<<"Wrong cut type!"<<endl;
    return false;
  }

  if( !(inFile->IsOpen()) )
  {
    cout<<"Unable to open file with production plane histograms!"<<endl;
    return false;
  }

  inFile->cd();

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



  K0s_K0s_cosThetaProdPlane_US_hist = (TH1F*)inFile->Get("K0s_K0s_cosThetaProdPlane_US_hist");
  K0s_K0s_cosThetaProdPlane_LS_hist = (TH1F*)inFile->Get("K0s_K0s_cosThetaProdPlane_LS_hist");
  K0s_K0s_cosThetaProdPlane_ME_hist = (TH1F*)inFile->Get("K0s_K0s_cosThetaProdPlane_ME_hist");
  K0s_K0s_cosThetaProdPlane_ME_LS_hist = (TH1F*)inFile->Get("K0s_K0s_cosThetaProdPlane_ME_LS_hist");

  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      K0s_K0s_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[pTbin1][pTbin2] = (TH1F*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_ME_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      K0s_K0s_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[etaBin1][etaBin2] = (TH1F*)inFile->Get(Form("K0s_K0s_cosThetaProdPlane_ME_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
    }
  }


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TCanvas *K0s_K0s_cosThetaProdPlane_no_corr_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_no_corr_can"), Form("K0s_K0s_cosThetaProdPlane_no_corr_can"), 1200, 1000);

  K0s_K0s_cosThetaProdPlane_no_corr_can->cd();

  //signal + background
  K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_US_hist->GetYaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_US_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_US_hist->SetMarkerColor(kRed);
  K0s_K0s_cosThetaProdPlane_US_hist->SetLineColor(kRed);
  double nK0sK0s = K0s_K0s_cosThetaProdPlane_US_hist->Integral();
  K0s_K0s_cosThetaProdPlane_US_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_US_hist->Scale(1./K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->GetBinWidth(1)); //bin width
  K0s_K0s_cosThetaProdPlane_US_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_US_hist->Draw("p e");

  K0s_K0s_cosThetaProdPlane_ME_hist->SetMarkerStyle(24);
  K0s_K0s_cosThetaProdPlane_ME_hist->SetMarkerColor(1);
  K0s_K0s_cosThetaProdPlane_ME_hist->SetLineColor(1);
  K0s_K0s_cosThetaProdPlane_ME_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_ME_hist->Scale(nK0sK0s/K0s_K0s_cosThetaProdPlane_ME_hist->Integral()); //scale ME to match signal+background
  K0s_K0s_cosThetaProdPlane_ME_hist->Scale(1./K0s_K0s_cosThetaProdPlane_ME_hist->GetXaxis()->GetBinWidth(1));
  K0s_K0s_cosThetaProdPlane_ME_hist->Draw("same p e");

  //background
  K0s_K0s_cosThetaProdPlane_LS_hist->SetMarkerStyle(21);
  K0s_K0s_cosThetaProdPlane_LS_hist->SetMarkerColor(kBlue);
  double nK0sK0s_back = K0s_K0s_cosThetaProdPlane_LS_hist->Integral();
  K0s_K0s_cosThetaProdPlane_LS_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_LS_hist->Scale(1./K0s_K0s_cosThetaProdPlane_LS_hist->GetXaxis()->GetBinWidth(1)); //bin width
  K0s_K0s_cosThetaProdPlane_LS_hist->Draw("p e same");

  K0s_K0s_cosThetaProdPlane_ME_LS_hist->SetMarkerStyle(25);
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->SetMarkerColor(kMagenta+1);
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->SetLineColor(kMagenta+1);
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->Scale(nK0sK0s_back/K0s_K0s_cosThetaProdPlane_ME_LS_hist->Integral()); //scale ME to match background
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->Scale(1./K0s_K0s_cosThetaProdPlane_ME_LS_hist->GetXaxis()->GetBinWidth(1));
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->Draw("p e same");

  TPaveText *K0s_K0s_text_no_corr = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_text_no_corr->SetTextFont(42);
  K0s_K0s_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  K0s_K0s_text_no_corr->AddText("Minimum bias, no correction");
  K0s_K0s_text_no_corr->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_text_no_corr->AddText("|#eta| < 1");
  K0s_K0s_text_no_corr->AddText("p_{T} integrated");
  //K0s_K0s_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_K0s_K0s, P_K0s_K0s_err));
  K0s_K0s_text_no_corr->SetFillColorAlpha(0, 0.01);
  K0s_K0s_text_no_corr->Draw("same");

  TLegend *K0s_K0s_leg_no_corr = new TLegend(0.15, 0.5, 0.39, 0.69);
  K0s_K0s_leg_no_corr->AddEntry(K0s_K0s_cosThetaProdPlane_US_hist, "Unlike-sign");
  K0s_K0s_leg_no_corr->AddEntry(K0s_K0s_cosThetaProdPlane_LS_hist, "Combinatorial bckg.");
  K0s_K0s_leg_no_corr->SetBorderSize(0);
  K0s_K0s_leg_no_corr->SetFillColorAlpha(0, 0.01);
  K0s_K0s_leg_no_corr->Draw("same");

  K0s_K0s_cosThetaProdPlane_no_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations/K0s_K0s_cosThetaProdPlane_no_corr.png");

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  TCanvas *K0s_K0s_cosThetaProdPlane_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_can"), Form("K0s_K0s_cosThetaProdPlane_can"), 1200, 1000);

  K0s_K0s_cosThetaProdPlane_can->cd();

  //signal + background


  K0s_K0s_cosThetaProdPlane_US_hist->Divide(K0s_K0s_cosThetaProdPlane_ME_hist); //correction of signal+bacground using ME
  K0s_K0s_cosThetaProdPlane_US_hist->Scale(nK0sK0s/K0s_K0s_cosThetaProdPlane_US_hist->Integral()); //scale back to US integral
  K0s_K0s_cosThetaProdPlane_US_hist->Scale(1./K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->GetBinWidth(1)); //bin width
  K0s_K0s_cosThetaProdPlane_US_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_US_hist->Draw("p e");

  //K0s_K0s_cosThetaProdPlane_ME_hist->Draw("p e same"); //for plotting ME with uncorrected distributions


  //background
  K0s_K0s_cosThetaProdPlane_LS_hist->Divide(K0s_K0s_cosThetaProdPlane_ME_LS_hist); //correction of background using ME
  K0s_K0s_cosThetaProdPlane_LS_hist->Scale(nK0sK0s_back/K0s_K0s_cosThetaProdPlane_LS_hist->Integral()); //scale back to raw background integral
  K0s_K0s_cosThetaProdPlane_LS_hist->Scale(1./K0s_K0s_cosThetaProdPlane_LS_hist->GetXaxis()->GetBinWidth(1)); //bin width
  K0s_K0s_cosThetaProdPlane_LS_hist->Draw("p e same");

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

  //fit_res_gaus_wrong_sign = L_inv_mass_US_LS->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  K0s_K0s_cosThetaProdPlane_LS_hist->Fit(fitK0s_K0s_US_LS_ThetaStar, "s i 0 r");

  float P_K0s_K0s_back = fitK0s_K0s_US_LS_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_K0s_K0s_back_err = fitK0s_K0s_US_LS_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

  fitK0s_K0s_US_LS_ThetaStar->SetLineColor(1);
  fitK0s_K0s_US_LS_ThetaStar->Draw("same");


  TPaveText *K0s_K0s_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_text->SetTextFont(42);
  K0s_K0s_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  K0s_K0s_text->AddText("Minimum bias");
  K0s_K0s_text->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_text->AddText("|#eta| < 1");
  K0s_K0s_text->AddText("p_{T} integrated");
  K0s_K0s_text->AddText(Form("P = %.2f #pm %.2f", P_K0s_K0s, P_K0s_K0s_err));
  K0s_K0s_text->SetFillColorAlpha(0, 0.01);
  K0s_K0s_text->Draw("same");


  TLegend *K0s_K0s_leg = new TLegend(0.15, 0.5, 0.39, 0.69);
  K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_US_hist, "Unlike-sign");
  //K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_US_side_band_hist, "Side-band p#pi");
  K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_LS_hist, "Combinatorial bckg.");
  //K0s_K0s_leg->AddEntry(K0s_K0s_cosThetaProdPlane_ME_hist, "Mixed event");
  K0s_K0s_leg->AddEntry(fitK0s_K0s_US_ThetaStar, "Linear fit to US");
  K0s_K0s_leg->SetBorderSize(0);
  K0s_K0s_leg->SetFillColorAlpha(0, 0.01);
  K0s_K0s_leg->Draw("same");

  K0s_K0s_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations/K0s_K0s_cosThetaProdPlane.png");

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  TCanvas *K0s_K0s_cosThetaProdPlane_2_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_2_can"), Form("K0s_K0s_cosThetaProdPlane_2_can"), 1200, 1000);

  K0s_K0s_cosThetaProdPlane_2_can->cd();

  //for background subtraction
  K0s_K0s_cosThetaProdPlane_US_hist->Add(K0s_K0s_cosThetaProdPlane_LS_hist, -1); //subtract background from signal+background
  K0s_K0s_cosThetaProdPlane_US_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_US_hist->Draw("p e");

  TF1 *fitK0s_K0s_US_ThetaStar_2 = new TF1("fitK0s_K0s_US_ThetaStar_2", "[0]*(1 + [1]*x)", -1, 1);
  fitK0s_K0s_US_ThetaStar_2->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  K0s_K0s_cosThetaProdPlane_US_hist->Fit(fitK0s_K0s_US_ThetaStar_2, "s i 0 r");

  //store fit result for systematic errors
  //tight topo cuts
  if(cutType == 1)
  {
    SysErrCutsTopo->cd();

    fitK0s_K0s_tight_topo_cuts->SetParameters(fitK0s_K0s_US_ThetaStar_2->GetParameter(0), fitK0s_K0s_US_ThetaStar_2->GetParameter(1));
    fitK0s_K0s_tight_topo_cuts->Write();
  }

  //tight pT cuts
  if(cutType == 2)
  {
    SysErrCutsPt->cd();

    fitK0s_K0s_tight_pT_cuts->SetParameters(fitK0s_K0s_US_ThetaStar_2->GetParameter(0), fitK0s_K0s_US_ThetaStar_2->GetParameter(1));
    fitK0s_K0s_tight_pT_cuts->Write();
  }

  float P_K0s_K0s_2 = fitK0s_K0s_US_ThetaStar_2->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_K0s_K0s_err_2 = fitK0s_K0s_US_ThetaStar_2->GetParError(1)/(L0_alpha*L0_alpha);

  fitK0s_K0s_US_ThetaStar_2->SetLineColor(1);
  fitK0s_K0s_US_ThetaStar_2->Draw("same");

/*
  //calculate total systematic uncertainty for L-Lbar

  //slope
  //float SysErrSlope_K0s_K0s = SysErrSlope_hist->GetBinContent(2);
  float SysErrSlope_K0s_K0s = 0;

  //background subtraction
  float P_K0s_K0s_from_fits = P_K0s_K0s - nK0sK0s_back/(nK0sK0s-nK0sK0s_back)*P_K0s_K0s_back;
  float P_K0s_K0s_from_fits_err = sqrt( P_K0s_K0s_err*P_K0s_K0s_err + nK0sK0s_back/(nK0sK0s-nK0sK0s_back)*nK0sK0s_back/(nK0sK0s-nK0sK0s_back)*P_K0s_K0s_back_err*P_K0s_K0s_back_err );

  //statistical error correction for systematic error
  float SysErrBackground_K0s_K0s_corr = sqrt( fabs( P_K0s_K0s_err_2*P_K0s_K0s_err_2 - P_K0s_K0s_from_fits_err*P_K0s_K0s_from_fits_err ) );

  float SysErrBackground_K0s_K0s_work = ( fabs( P_K0s_K0s_2 - P_K0s_K0s_from_fits) - SysErrBackground_K0s_K0s_corr )/fabs(P_K0s_K0s_2);

  float SysErrBackground_K0s_K0s = 0;

  if( SysErrBackground_K0s_K0s_work > 0 ) SysErrBackground_K0s_K0s = SysErrBackground_K0s_K0s_work; //store sys. err. only if it is larger than statistical fluctuations


  //cuts variation
  float SysErr_tight_cuts_K0s_K0s = 0;

  if(cutType == 0)
  {
    float P_K0s_K0s_tight_topo_cuts = fitK0s_K0s_tight_topo_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_K0s_K0s_tight_topo_cuts_err = fitK0s_K0s_tight_topo_cuts->GetParError(1)/(L0_alpha*L0_alpha);

    //statistical error correction for systematic error
    float SysErr_tight_topo_cuts_K0s_K0s_corr = sqrt( fabs( P_K0s_K0s_err_2*P_K0s_K0s_err_2 - P_K0s_K0s_tight_topo_cuts_err*P_K0s_K0s_tight_topo_cuts_err ) );

    float SysErr_tight_topo_cuts_K0s_K0s_work = ( fabs( P_K0s_K0s_2 - P_K0s_K0s_tight_topo_cuts) - SysErr_tight_topo_cuts_K0s_K0s_corr )/fabs(P_K0s_K0s_2);

    float SysErr_tight_topo_cuts_K0s_K0s = 0;

    if( SysErr_tight_topo_cuts_K0s_K0s_work > 0 ) SysErr_tight_topo_cuts_K0s_K0s = SysErr_tight_topo_cuts_K0s_K0s_work; //store sys. err. only if it is larger than statistical fluctuations


    float P_K0s_K0s_tight_pT_cuts = fitK0s_K0s_tight_pT_cuts->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_K0s_K0s_tight_pT_cuts_err = fitK0s_K0s_tight_pT_cuts->GetParError(1)/(L0_alpha*L0_alpha);

    //statistical error correction for systematic error
    float SysErr_tight_pT_cuts_K0s_K0s_corr = sqrt( fabs( P_K0s_K0s_err_2*P_K0s_K0s_err_2 - P_K0s_K0s_tight_pT_cuts_err*P_K0s_K0s_tight_pT_cuts_err ) );

    float SysErr_tight_pT_cuts_K0s_K0s_work = ( fabs( P_K0s_K0s_2 - P_K0s_K0s_tight_pT_cuts) - SysErr_tight_pT_cuts_K0s_K0s_corr )/fabs(P_K0s_K0s_2);

    float SysErr_tight_pT_cuts_K0s_K0s = 0;

    if( SysErr_tight_pT_cuts_K0s_K0s_work > 0 ) SysErr_tight_pT_cuts_K0s_K0s = SysErr_tight_pT_cuts_K0s_K0s_work; //store sys. err. only if it is larger than statistical fluctuations


    SysErr_tight_cuts_K0s_K0s = sqrt(SysErr_tight_topo_cuts_K0s_K0s*SysErr_tight_topo_cuts_K0s_K0s + SysErr_tight_pT_cuts_K0s_K0s+SysErr_tight_pT_cuts_K0s_K0s);

  }

  //total
  //float SysErrTot_K0s_K0s = sqrt( sysErr_alpha_K0s_K0s*sysErr_alpha_K0s_K0s + SysErrSlope_K0s_K0s*SysErrSlope_K0s_K0s + SysErrBackground_K0s_K0s*SysErrBackground_K0s_K0s );
  float SysErrTot_K0s_K0s = sqrt( sysErr_alpha_K0s_K0s*sysErr_alpha_K0s_K0s + SysErrSlope_K0s_K0s*SysErrSlope_K0s_K0s + SysErrBackground_K0s_K0s*SysErrBackground_K0s_K0s + SysErr_tight_cuts_K0s_K0s*SysErr_tight_cuts_K0s_K0s);

  //--------------------------------------------------------------------------
*/
  TPaveText *K0s_K0s_text_2 = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
  //cent_text_2->AddText("STAR preliminary");
  //((TText*)cent_text_2->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_text_2->SetTextFont(42);
  K0s_K0s_text_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  K0s_K0s_text_2->AddText("Minimum bias");
  K0s_K0s_text_2->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_text_2->AddText("|#eta| < 1");
  K0s_K0s_text_2->AddText("p_{T} integrated");
  K0s_K0s_text_2->AddText(Form("P = %.2f #pm %.2f", P_K0s_K0s_2, P_K0s_K0s_err_2));
  K0s_K0s_text_2->SetFillColorAlpha(0, 0.01);
  K0s_K0s_text_2->Draw("same");

  TLegend *K0s_K0s_2_leg = new TLegend(0.15, 0.5, 0.39, 0.69);
  K0s_K0s_2_leg->AddEntry(K0s_K0s_cosThetaProdPlane_US_hist, "US-Bckg.");
  //K0s_K0s_2_leg->AddEntry(K0s_K0s_cosThetaProdPlane_ME_hist, "Mixed event");
  K0s_K0s_2_leg->AddEntry(fitK0s_K0s_US_ThetaStar_2, "Linear fit to US-Bckg.");
  K0s_K0s_2_leg->SetBorderSize(0);
  K0s_K0s_2_leg->SetFillColorAlpha(0, 0.01);
  K0s_K0s_2_leg->Draw("same");

  K0s_K0s_cosThetaProdPlane_2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations/K0s_K0s_cosThetaProdPlane_subtract.png");

  //__________________________________________________________________________________________________________________________________________________________________________

  OutFile_pol->cd();

  TGraphErrors *PolarizationGraph_K0s = new TGraphErrors(1);
  PolarizationGraph_K0s->SetNameTitle("PolarizationGraph_K0s","PolarizationGraph_K0s");

  PolarizationGraph_K0s->SetPoint(1, P_K0s_K0s_2, 1);
  PolarizationGraph_K0s->SetPointError(1, fabs(P_K0s_K0s_err_2), 0);

  PolarizationGraph_K0s->Write();

  //__________________________________________________________________________________________________________________________________________________________________________



  //----------------------------------Histograms, pT binnig, main analyisis---------------------------------

  cout<<endl;
  //cout<<"N events with K0s pair: "<<nEventsWithK0sPair_US_hist->GetBinContent(2)<<endl;
  cout<<"N K0s pairs from hist: "<<nK0sK0s<<endl;
  cout<<"N K0s background pairs from hist: "<<nK0sK0s_back<<endl;
/*
  cout<<endl;
  cout<<"Systematic errors"<<endl;
  cout<<"K0s-K0s "<<SysErrSlope_K0s_K0s<<" "<<SysErrBackground_K0s_K0s<<" "<<sysErr_alpha_K0s_K0s<<" "<<SysErr_tight_cuts_K0s_K0s<<endl;
*/
  inFile->Close();
  OutFile_pol->Close();


  return true;

}
