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
bool cuts(float L_y)
{

  if( !( TMath::Abs(L_y) < L_y_cut ) ) return false;
  //if( strictTOF_cut == 1 && pi_hasTOFinfo == 0 ) return false; //TOF matched pions
  //if( strictTOF_cut == 2 && (pi_hasTOFinfo == 0 || p_hasTOFinfo == 0) ) return false; //TOF matched pions and protons
  //if(cos(L_theta) < L_cos_theta_cut) return false;
  //if(L_decayL > L_decayL_cut) return false;

  return true;

}


bool DoAnalysis(TChain *L_tree, const int ReadMode, const int energy = 510, const int year = 2017)
{


  TFile *InvMassFile; //output file to store invariant mass histograms

  if(ReadMode == 0) //create invariant mass file from nTuple - run in this mode first
  {
    InvMassFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/InvMass_Lambda_work.root", "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT
  }
  else  //read file to store wrong and good sign M_inv distributions and daughter pT
  {

    InvMassFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/InvMass_Lambda_work.root", "read"); //old, non-optimized cuts

    if( !(InvMassFile->IsOpen()) )
    {
      cout<<"Unable to open file with invariant mass histograms!"<<endl;
      return false;
    }
  }
  //_______________________________________________________________________________________________________________________________________________

  //need to change binning
  // 2D - pT1 vs. pT2
  // 2D - eta1 vs. eta2 (maybe laterr)
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_LS[nPtBins_corr][nPtBins_corr]; //for LS-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass[nPtBins_corr][nPtBins_corr];

  TH2F *L0_inv_mass_vs_L0_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_LS[nPtBins_corr][nPtBins_corr]; //for LS-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass[nPtBins_corr][nPtBins_corr];

  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_LS[nPtBins_corr][nPtBins_corr]; //for LS-LS Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass[nPtBins_corr][nPtBins_corr];



  if(ReadMode == 0) //create histograms to be saved into file
  {
    for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
    {
      for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
      {

        //invariant mass histograms
        //old bins 200
        L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_pT_%i_eta_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_pT_%i_eta_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
        L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_pT_%i_eta_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_pT_%i_eta_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
        L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0bar_inv_mass_LS_pT_%i_eta_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_LS_pT_%i_eta_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);

        L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_pT_%i_eta_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_pT_%i_eta_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
        L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_pT_%i_eta_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_LS_pT_%i_eta_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
        L0_inv_mass_vs_L0_inv_mass_LS[pTbin1][pTbin2] = new TH2F(Form("L0_inv_mass_vs_L0_inv_mass_LS_pT_%i_eta_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_LS_pT_%i_eta_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);

        L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_pT_%i_eta_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_pT_%i_eta_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
        L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_pT_%i_eta_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_pT_%i_eta_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
        L0bar_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2] = new TH2F(Form("L0bar_inv_mass_vs_L0bar_inv_mass_LS_pT_%i_eta_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_LS_pT_%i_eta_%i", pTbin1, pTbin2), 180, 1, 1.2, 180, 1, 1.2);
        //______________________________________________________________________________________________________________________________

      }
    }
  }
  else //load histograms from file
  {
    for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
    {
      for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
      {
        L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_US_pT_%i_eta_%i", pTbin1, pTbin2));
        L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_pT_%i_eta_%i", pTbin1, pTbin2));
        L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_LS_pT_%i_eta_%i", pTbin1, pTbin2));

        L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_US_pT_%i_eta_%i", pTbin1, pTbin2));
        L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_pT_%i_eta_%i", pTbin1, pTbin2));
        L0_inv_mass_vs_L0_inv_mass_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_LS_pT_%i_eta_%i", pTbin1, pTbin2));

        L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_pT_%i_eta_%i", pTbin1, pTbin2));
        L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_pT_%i_eta_%i", pTbin1, pTbin2));
        L0bar_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_LS_pT_%i_eta_%i", pTbin1, pTbin2));
        //______________________________________________________________________________________________________________________________


      }
    }
  }
  //________________________________________________________________________________________


  int eventID_last = -1; //to store event ID from last L candidate
  int eventID_last_background = -1; //to store event ID from last L candidate

  Long64_t nEntries = 0; //total nEntries

  Int_t charge;
  Float_t L_mass, L_pt, L_eta, L_phi;//, L_decayL, L_theta, L_DCAdaughters;

  Float_t pi_pt, p_pt;
  Float_t pi_eta, p_eta;
  Float_t pi_phi, p_phi;
  Int_t p_ch;
  //Float_t pi_dca, p_dca;
  Int_t pi_hasTOFinfo, p_hasTOFinfo;

  //Float_t thetaProdPlane;

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
    //L_tree->SetBranchAddress("pair_decayL", &L_decayL);
    //L_tree->SetBranchAddress("pair_theta", &L_theta);
    //L_tree->SetBranchAddress("pair_DCAdaughters", &L_DCAdaughters);

    //pion is particle 2 in the pair niside the TTree
    L_tree->SetBranchAddress("p2_pt", &pi_pt);
    L_tree->SetBranchAddress("p2_eta", &pi_eta);
    L_tree->SetBranchAddress("p2_phi", &pi_phi);
    //L_tree->SetBranchAddress("p2_dca", &pi_dca);
    L_tree->SetBranchAddress("p2_hasTOFinfo", &pi_hasTOFinfo);

    //proton is particle 1 in the pair niside the TTree
    L_tree->SetBranchAddress("p1_pt", &p_pt);
    L_tree->SetBranchAddress("p1_eta", &p_eta);
    L_tree->SetBranchAddress("p1_phi", &p_phi);
    //L_tree->SetBranchAddress("p1_dca", &p_dca);
    L_tree->SetBranchAddress("p1_ch", &p_ch);
    L_tree->SetBranchAddress("p1_hasTOFinfo", &p_hasTOFinfo);

    //L_tree->SetBranchAddress("thetaProdPlane", &thetaProdPlane);

    L_tree->SetBranchAddress("eventId", &eventId);

    //--------------------------------------------------------------------------


    nEntries = L_tree->GetEntries();
    cout<<"nEntries = "<<nEntries<<endl;
  }

  //to store Lorentz vectors for L pair analysis
  //p vectors to calculate cos(theta*) and for auto-correlation check of L-L and Lbar-Lbar
  //pi vectors for auto-correlation check of L-L and Lbar-Lbar

  //unlike-sign (signal+background)
  vector<TLorentzVector> L_vector;
  vector<int> L_pT_bin_vector;
  vector<int> L_eta_bin_vector;
  vector<TLorentzVector> pi_vector;
  vector<TLorentzVector> p_vector;

  vector<TLorentzVector> Lbar_vector;
  vector<int> Lbar_pT_bin_vector;
  vector<int> Lbar_eta_bin_vector;
  vector<TLorentzVector> piBar_vector;
  vector<TLorentzVector> pBar_vector;

  //vectors to store cos(theta*) info
  //analyzed later, after Minv window is found
  vector<float> L_Lbar_cos_theta;
  vector<float> L_Lbar_Minv_L;
  vector<float> L_Lbar_Minv_Lbar;
  vector<int> L_Lbar_pT_bin_L;
  vector<int> L_Lbar_pT_bin_Lbar;
  vector<int> L_Lbar_eta_bin_L;
  vector<int> L_Lbar_eta_bin_Lbar;

  vector<float> L_L_cos_theta;
  vector<float> L_L_Minv_L1;
  vector<float> L_L_Minv_L2;
  vector<int> L_L_pT_bin_L1;
  vector<int> L_L_pT_bin_L2;
  vector<int> L_L_eta_bin_L1;
  vector<int> L_L_eta_bin_L2;

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
  vector<TLorentzVector> p_vector_background;

  vector<TLorentzVector> Lbar_vector_background;
  vector<int> Lbar_pT_bin_vector_background;
  vector<int> Lbar_eta_bin_vector_background;
  vector<TLorentzVector> piBar_vector_background;
  vector<TLorentzVector> pBar_vector_background;


  vector<float> L_Lbar_cos_theta_back;
  vector<float> L_Lbar_Minv_L_back;
  vector<float> L_Lbar_Minv_Lbar_back;
  vector<int> L_Lbar_pT_bin_L_back;
  vector<int> L_Lbar_pT_bin_Lbar_back;
  vector<int> L_Lbar_eta_bin_L_back;
  vector<int> L_Lbar_eta_bin_Lbar_back;

  vector<float> L_L_cos_theta_back;
  vector<float> L_L_Minv_L1_back;
  vector<float> L_L_Minv_L2_back;
  vector<int> L_L_pT_bin_L1_back;
  vector<int> L_L_pT_bin_L2_back;
  vector<int> L_L_eta_bin_L1_back;
  vector<int> L_L_eta_bin_L2_back;

  vector<float> Lbar_Lbar_cos_theta_back;
  vector<float> Lbar_Lbar_Minv_Lbar1_back;
  vector<float> Lbar_Lbar_Minv_Lbar2_back;
  vector<int> Lbar_Lbar_pT_bin_Lbar1_back;
  vector<int> Lbar_Lbar_pT_bin_Lbar2_back;
  vector<int> Lbar_Lbar_eta_bin_Lbar1_back;
  vector<int> Lbar_Lbar_eta_bin_Lbar2_back;


  //vectors for mixed-event
  vector<TLorentzVector> L_vector_ME;
  vector<int> L_pT_bin_vector_ME;
  vector<int> L_eta_bin_vector_ME;
  vector<TLorentzVector> p_vector_ME;

  vector<TLorentzVector> Lbar_vector_ME;
  vector<int> Lbar_pT_bin_vector_ME;
  vector<int> Lbar_eta_bin_vector_ME;
  vector<TLorentzVector> pBar_vector_ME;


  //TLorentzVector *L_Lorentz_vector = new TLorentzVector(1,1,1,1);
  //TLorentzVector *p_Lorenz_vector = new TLorentzVector(1,1,1,1);

  //for filling of US-LS L-L and Lbar-Lbar pair background
  //to be able to fill bacground Minv histograms correctly
  int nFillLL_bckg = 0;
  int nFillLbarLbar_backg = 0;

  for(Long64_t i = 0; i < nEntries; i++) //Read TTree only in ReadMode = 0
  {
    L_tree->GetEntry(i);

     //if(ReadMode != 0) break;
    if(i%1000000 == 0)
    {
      cout<<i<<endl;
    }

    //double L_xF = fabs(pz(L_pt, L_eta))/energy/2.; //energy is in CMS, need energz of one proton

    //calculate Lambda rapidity y
    double L_y = rapidity(L_pt, L_eta, L_mass_PDG);

    if(i == 0 ) //first iteration
    {
      eventID_last = eventId;
    }

    //------------------------------------------------------------------------------------------------------------------

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

    if( pT_bin_corr == -1 ) continue;


    //fill all histograms for all eta and centrality bins
    int eta_bin = -1;

    //find eta bin of Lambda
    for(int j = 0; j < nEtaBins; j++) //loop over eta bins
    {
      if(L_eta > eta_bins[j] && L_eta <= eta_bins[j+1])
      {
        eta_bin = j;
        break; //stop after eta bin is found
      }
    }

    if( eta_bin == -1 ) continue;

/*
    TLorentzVector *L_Lorentz_vector = new TLorentzVector(1,1,1,1);
    L_Lorentz_vector->SetPtEtaPhiM(L_pt, L_eta, L_phi, L_mass);

    TLorentzVector *p_Lorenz_vector = new TLorentzVector(1,1,1,1);
    p_Lorenz_vector->SetPtEtaPhiM(p_pt, p_eta, p_phi, p_mass_PDG);
*/
    TLorentzVector L_Lorentz_vector(1,1,1,1);
    L_Lorentz_vector.SetPtEtaPhiM(L_pt, L_eta, L_phi, L_mass);

    TLorentzVector p_Lorenz_vector(1,1,1,1);
    p_Lorenz_vector.SetPtEtaPhiM(p_pt, p_eta, p_phi, p_mass_PDG);

    TLorentzVector pi_Lorenz_vector(1,1,1,1);
    pi_Lorenz_vector.SetPtEtaPhiM(pi_pt, pi_eta, pi_phi, pi_mass_PDG);

    if(eventId == eventID_last)
    //if(charge == 0 ) //like-sign combinations
    {
      //US charge combinations
      if(charge == 0)
      {
        //cuts
        //if( cuts(requireTOF_Lcorr, pi_hasTOFinfo, p_hasTOFinfo, L_y) && pT_bin != -1 && eta_bin != -1 && pT_bin_corr != -1)
        if( cuts(L_y) && pT_bin_corr != -1)
        {
          if( p_ch == 1 )
          {
            L_vector.push_back(L_Lorentz_vector);
            L_pT_bin_vector.push_back(pT_bin_corr);
            L_eta_bin_vector.push_back(eta_bin);

            pi_vector.push_back(pi_Lorenz_vector);
            p_vector.push_back(p_Lorenz_vector);

          }
          else if( p_ch == -1 )
          {
            Lbar_vector.push_back(L_Lorentz_vector);
            Lbar_pT_bin_vector.push_back(pT_bin_corr);
            Lbar_eta_bin_vector.push_back(eta_bin);

            piBar_vector.push_back(pi_Lorenz_vector);
            pBar_vector.push_back(p_Lorenz_vector);

          }
        }
      }

      //LS charge combinations
      if(charge != 0) //same event as in previous iteration
      {
        //cuts
        if( cuts(L_y) && pT_bin_corr != -1)
        {
          if( p_ch == 1 )
          {
            L_vector_background.push_back(L_Lorentz_vector);
            L_pT_bin_vector_background.push_back(pT_bin_corr);
            L_eta_bin_vector_background.push_back(eta_bin);

            pi_vector_background.push_back(pi_Lorenz_vector);
            p_vector_background.push_back(p_Lorenz_vector);
          }
          else if( p_ch == -1 )
          {
            Lbar_vector_background.push_back(L_Lorentz_vector);
            Lbar_pT_bin_vector_background.push_back(pT_bin_corr);
            Lbar_eta_bin_vector_background.push_back(eta_bin);

            piBar_vector_background.push_back(pi_Lorenz_vector);
            pBar_vector_background.push_back(p_Lorenz_vector);

          }
        }
      }

    }
    //else //unlike-sign combinations
    else if(eventId != eventID_last)//new event
    {
      eventID_last = eventId; //store new event ID

      // US-US L pairs - signal + bacground
      //at least one L-Lbar pair in event
      if(L_vector.size() > 0 && Lbar_vector.size() > 0)
      {
        //cout<<eventId<<endl;

        for(unsigned int iLambda = 0; iLambda < L_vector.size(); iLambda++)
        {
          for(unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector.size(); iLambdaBar++)
          {
            //fill Minv histogram

            L0_inv_mass_vs_L0bar_inv_mass_US[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(L_vector.at(iLambda).M(), Lbar_vector.at(iLambdaBar).M());

            //fill new vectors here with Minv1 and Minv2, and with cos(theta*)
            //analze stored pairs later, just within selected Minv window
            //store also pT bin of both L in the pair (maybe eta)

            //save info for cos(theta*) to be analyzed later
            double L_Lbar_pairThetaStar = LpairThetaStar(L_vector.at(iLambda), p_vector.at(iLambda), Lbar_vector.at(iLambdaBar), pBar_vector.at(iLambdaBar));

            L_Lbar_cos_theta.push_back(TMath::Cos(L_Lbar_pairThetaStar));

            L_Lbar_Minv_L.push_back(L_vector.at(iLambda).M());
            L_Lbar_Minv_Lbar.push_back(Lbar_vector.at(iLambdaBar).M());

            L_Lbar_pT_bin_L.push_back(L_pT_bin_vector.at(iLambda));
            L_Lbar_pT_bin_Lbar.push_back(Lbar_pT_bin_vector.at(iLambdaBar));

            L_Lbar_eta_bin_L.push_back(L_eta_bin_vector.at(iLambda));
            L_Lbar_eta_bin_Lbar.push_back(Lbar_eta_bin_vector.at(iLambdaBar));

          }
        }
      }

      //at least one L0-L0 pair in event
      if(L_vector.size() > 1)
      {
        for(unsigned int iLambda1 = 0; iLambda1 < L_vector.size(); iLambda1++)
        {
          for(unsigned int iLambda2 = iLambda1+1; iLambda2 < L_vector.size(); iLambda2++)
          {
            //if( fabs(p_vector.at(iLambda1).Rapidity()) < 0.001 || fabs(p_vector.at(iLambda2).Rapidity()) < 0.001 ) continue;

            if(pi_vector.at(iLambda1).Rapidity() == pi_vector.at(iLambda2).Rapidity()) continue;
            if(p_vector.at(iLambda1).Rapidity() == p_vector.at(iLambda2).Rapidity()) continue;

            L0_inv_mass_vs_L0_inv_mass_US[L_pT_bin_vector.at(iLambda1)][L_pT_bin_vector.at(iLambda2)]->Fill(L_vector.at(iLambda1).M(), L_vector.at(iLambda2).M());

            //save info for cos(theta*) to be analyzed later
            double L_L_pairThetaStar = LpairThetaStar(L_vector.at(iLambda1), p_vector.at(iLambda1), L_vector.at(iLambda2), p_vector.at(iLambda2));

            L_L_cos_theta.push_back(TMath::Cos(L_L_pairThetaStar));

            L_L_Minv_L1.push_back(L_vector.at(iLambda1).M());
            L_L_Minv_L2.push_back(L_vector.at(iLambda2).M());

            L_L_pT_bin_L1.push_back(L_pT_bin_vector.at(iLambda1));
            L_L_pT_bin_L2.push_back(L_pT_bin_vector.at(iLambda2));

            L_L_eta_bin_L1.push_back(L_eta_bin_vector.at(iLambda1));
            L_L_eta_bin_L2.push_back(L_eta_bin_vector.at(iLambda2));

          }
        }
      }

      //at least one L0bar-L0bar pair in event
      if(Lbar_vector.size() > 1)
      {
        for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < Lbar_vector.size(); iLambdaBar1++)
        {
          for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < Lbar_vector.size(); iLambdaBar2++)
          {
            if(piBar_vector.at(iLambdaBar1).Rapidity() == piBar_vector.at(iLambdaBar2).Rapidity()) continue;
            if(pBar_vector.at(iLambdaBar1).Rapidity() == pBar_vector.at(iLambdaBar2).Rapidity()) continue;

            //L0bar_inv_mass_vs_L0bar_inv_mass_US[Lbar_pT_bin_vector.at(iLambdaBar1)][Lbar_pT_bin_vector.at(iLambdaBar2)]->Fill(Lbar_vector.at(iLambdaBar1).M(), Lbar_vector.at(iLambdaBar2).M());

            //try swithc axes - for testing of the fit
            L0bar_inv_mass_vs_L0bar_inv_mass_US[Lbar_pT_bin_vector.at(iLambdaBar2)][Lbar_pT_bin_vector.at(iLambdaBar1)]->Fill(Lbar_vector.at(iLambdaBar2).M(), Lbar_vector.at(iLambdaBar1).M());

            //save info for cos(theta*) to be analyzed later
            double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_vector.at(iLambdaBar1), pBar_vector.at(iLambdaBar1), Lbar_vector.at(iLambdaBar2), pBar_vector.at(iLambdaBar2));

            Lbar_Lbar_cos_theta.push_back(TMath::Cos(Lbar_Lbar_pairThetaStar));

            Lbar_Lbar_Minv_Lbar1.push_back(Lbar_vector.at(iLambdaBar1).M());
            Lbar_Lbar_Minv_Lbar2.push_back(Lbar_vector.at(iLambdaBar2).M());

            Lbar_Lbar_pT_bin_Lbar1.push_back(Lbar_pT_bin_vector.at(iLambdaBar1));
            Lbar_Lbar_pT_bin_Lbar2.push_back(Lbar_pT_bin_vector.at(iLambdaBar2));

            Lbar_Lbar_eta_bin_Lbar1.push_back(Lbar_eta_bin_vector.at(iLambdaBar1));
            Lbar_Lbar_eta_bin_Lbar2.push_back(Lbar_eta_bin_vector.at(iLambdaBar2));


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
            //if( fabs(p_vector_background.at(iLambda).Rapidity()) < 0.001  ) continue;

            L0_inv_mass_vs_L0bar_inv_mass_LS[L_pT_bin_vector_background.at(iLambda)][Lbar_pT_bin_vector_background.at(iLambdaBar)]->Fill(L_vector_background.at(iLambda).M(), Lbar_vector_background.at(iLambdaBar).M());

            /*
            //save info for cos(theta*) to be analyzed later
            double L_Lbar_pairThetaStar = LpairThetaStar(L_vector_background.at(iLambda), p_vector_background.at(iLambda), Lbar_vector_background.at(iLambdaBar), pBar_vector_background.at(iLambdaBar));

            L_Lbar_cos_theta_back.push_back(TMath::Cos(L_Lbar_pairThetaStar));

            L_Lbar_Minv_L_back.push_back(L_vector_background.at(iLambda).M());
            L_Lbar_Minv_Lbar_back.push_back(Lbar_vector_background.at(iLambdaBar).M());

            L_Lbar_pT_bin_L_back.push_back(L_pT_bin_vector_background.at(iLambda));
            L_Lbar_pT_bin_Lbar_back.push_back(Lbar_pT_bin_vector_background.at(iLambdaBar));

            L_Lbar_eta_bin_L_back.push_back(L_eta_bin_vector_background.at(iLambda));
            L_Lbar_eta_bin_Lbar_back.push_back(Lbar_eta_bin_vector_background.at(iLambdaBar));
            */

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

            L0_inv_mass_vs_L0_inv_mass_LS[L_pT_bin_vector_background.at(iLambda1)][L_pT_bin_vector_background.at(iLambda2)]->Fill(L_vector_background.at(iLambda1).M(), L_vector_background.at(iLambda2).M());

            /*
            //save info for cos(theta*) to be analyzed later
            double L_L_pairThetaStar = LpairThetaStar(L_vector_background.at(iLambda1), p_vector_background.at(iLambda1), L_vector_background.at(iLambda2), p_vector_background.at(iLambda2));

            L_L_cos_theta_back.push_back(TMath::Cos(L_L_pairThetaStar));

            L_L_Minv_L1_back.push_back(L_vector_background.at(iLambda1).M());
            L_L_Minv_L2_back.push_back(L_vector_background.at(iLambda2).M());

            L_L_pT_bin_L1_back.push_back(L_pT_bin_vector_background.at(iLambda1));
            L_L_pT_bin_L2_back.push_back(L_pT_bin_vector_background.at(iLambda2));

            L_L_eta_bin_L1_back.push_back(L_eta_bin_vector_background.at(iLambda1));
            L_L_eta_bin_L2_back.push_back(L_eta_bin_vector_background.at(iLambda2));
            */

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

            L0bar_inv_mass_vs_L0bar_inv_mass_LS[Lbar_pT_bin_vector_background.at(iLambdaBar1)][Lbar_pT_bin_vector_background.at(iLambdaBar2)]->Fill(Lbar_vector_background.at(iLambdaBar1).M(), Lbar_vector_background.at(iLambdaBar2).M());

            /*
             //save info for cos(theta*) to be analyzed later
            double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_vector_background.at(iLambdaBar1), pBar_vector_background.at(iLambdaBar1), Lbar_vector_background.at(iLambdaBar2), pBar_vector_background.at(iLambdaBar2));

            Lbar_Lbar_cos_theta_back.push_back(TMath::Cos(Lbar_Lbar_pairThetaStar));

            Lbar_Lbar_Minv_Lbar1_back.push_back(Lbar_vector_background.at(iLambdaBar1).M());
            Lbar_Lbar_Minv_Lbar2_back.push_back(Lbar_vector_background.at(iLambdaBar2).M());

            Lbar_Lbar_pT_bin_Lbar1_back.push_back(Lbar_pT_bin_vector_background.at(iLambdaBar1));
            Lbar_Lbar_pT_bin_Lbar2_back.push_back(Lbar_pT_bin_vector_background.at(iLambdaBar2));

            Lbar_Lbar_eta_bin_Lbar1_back.push_back(Lbar_eta_bin_vector_background.at(iLambdaBar1));
            Lbar_Lbar_eta_bin_Lbar2_back.push_back(Lbar_eta_bin_vector_background.at(iLambdaBar2));
            */

          }
        }
      }

      //___________________________________________________________________________________________________


      // US-LS and LS-US L pairs - background - case when US L (can be signal or background) is paired with bacground L (+ continuum?)

      //Bacground for L-Lbar
      // L - US, Lbar - LS
      if(L_vector.size() > 0 && Lbar_vector_background.size() > 0)
      {
        //cout<<eventId<<endl;

        for(unsigned int iLambda = 0; iLambda < L_vector.size(); iLambda++)
        {
          for(unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector_background.size(); iLambdaBar++)
          {
            //fill Minv histogram
            if(pi_vector.at(iLambda).Rapidity() == piBar_vector_background.at(iLambdaBar).Rapidity()) continue;

            L0_inv_mass_vs_L0bar_inv_mass_US_LS[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector_background.at(iLambdaBar)]->Fill(L_vector.at(iLambda).M(), Lbar_vector_background.at(iLambdaBar).M());

            //save info for cos(theta*) to be analyzed later
            double L_Lbar_pairThetaStar = LpairThetaStar(L_vector.at(iLambda), p_vector.at(iLambda), Lbar_vector_background.at(iLambdaBar), pBar_vector_background.at(iLambdaBar));

            L_Lbar_cos_theta_back.push_back(TMath::Cos(L_Lbar_pairThetaStar));

            L_Lbar_Minv_L_back.push_back(L_vector.at(iLambda).M());
            L_Lbar_Minv_Lbar_back.push_back(Lbar_vector_background.at(iLambdaBar).M());

            L_Lbar_pT_bin_L_back.push_back(L_pT_bin_vector.at(iLambda));
            L_Lbar_pT_bin_Lbar_back.push_back(Lbar_pT_bin_vector_background.at(iLambdaBar));

            L_Lbar_eta_bin_L_back.push_back(L_eta_bin_vector.at(iLambda));
            L_Lbar_eta_bin_Lbar_back.push_back(Lbar_eta_bin_vector_background.at(iLambdaBar));

          }
        }
      }
      // L - LS, Lbar - US
      if(L_vector_background.size() > 0 && Lbar_vector.size() > 0)
      {
        //cout<<eventId<<endl;

        for(unsigned int iLambda = 0; iLambda < L_vector_background.size(); iLambda++)
        {
          for(unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector.size(); iLambdaBar++)
          {
            //fill Minv histogram
            if(pi_vector_background.at(iLambda).Rapidity() == piBar_vector.at(iLambdaBar).Rapidity()) continue;

            L0_inv_mass_vs_L0bar_inv_mass_US_LS[L_pT_bin_vector_background.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(L_vector_background.at(iLambda).M(), Lbar_vector.at(iLambdaBar).M());


            //save info for cos(theta*) to be analyzed later
            double L_Lbar_pairThetaStar = LpairThetaStar(L_vector_background.at(iLambda), p_vector_background.at(iLambda), Lbar_vector.at(iLambdaBar), pBar_vector.at(iLambdaBar));

            L_Lbar_cos_theta_back.push_back(TMath::Cos(L_Lbar_pairThetaStar));

            L_Lbar_Minv_L_back.push_back(L_vector_background.at(iLambda).M());
            L_Lbar_Minv_Lbar_back.push_back(Lbar_vector.at(iLambdaBar).M());

            L_Lbar_pT_bin_L_back.push_back(L_pT_bin_vector_background.at(iLambda));
            L_Lbar_pT_bin_Lbar_back.push_back(Lbar_pT_bin_vector.at(iLambdaBar));

            L_Lbar_eta_bin_L_back.push_back(L_eta_bin_vector_background.at(iLambda));
            L_Lbar_eta_bin_Lbar_back.push_back(Lbar_eta_bin_vector.at(iLambdaBar));


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
            if(p_vector.at(iLambda).Rapidity() == p_vector_background.at(iLambdaBckg).Rapidity() ) continue;

            if( nFillLL_bckg % 2 == 0 ) //fill LL pair wit US at x axis and LS at y axis
            {
              L0_inv_mass_vs_L0_inv_mass_US_LS[L_pT_bin_vector.at(iLambda)][L_pT_bin_vector_background.at(iLambdaBckg)]->Fill(L_vector.at(iLambda).M(), L_vector_background.at(iLambdaBckg).M());

              //save info for cos(theta*) to be analyzed later
              double L_L_pairThetaStar = LpairThetaStar(L_vector.at(iLambda), p_vector.at(iLambda), L_vector_background.at(iLambdaBckg), p_vector_background.at(iLambdaBckg));

              L_L_cos_theta_back.push_back(TMath::Cos(L_L_pairThetaStar));

              L_L_Minv_L1_back.push_back(L_vector.at(iLambda).M());
              L_L_Minv_L2_back.push_back(L_vector_background.at(iLambdaBckg).M());

              L_L_pT_bin_L1_back.push_back(L_pT_bin_vector.at(iLambda));
              L_L_pT_bin_L2_back.push_back(L_pT_bin_vector_background.at(iLambdaBckg));

              L_L_eta_bin_L1_back.push_back(L_eta_bin_vector.at(iLambda));
              L_L_eta_bin_L2_back.push_back(L_eta_bin_vector_background.at(iLambdaBckg));

            }
            else //fill LL pair wit US at y axis and LS at x axis, i.e. opposite order than above
            {
              L0_inv_mass_vs_L0_inv_mass_US_LS[L_pT_bin_vector_background.at(iLambdaBckg)][L_pT_bin_vector.at(iLambda)]->Fill(L_vector_background.at(iLambdaBckg).M(), L_vector.at(iLambda).M());


              //save info for cos(theta*) to be analyzed later
              double L_L_pairThetaStar = LpairThetaStar(L_vector_background.at(iLambdaBckg), p_vector_background.at(iLambdaBckg), L_vector.at(iLambda), p_vector.at(iLambda));

              L_L_cos_theta_back.push_back(TMath::Cos(L_L_pairThetaStar));

              L_L_Minv_L1_back.push_back(L_vector_background.at(iLambdaBckg).M());
              L_L_Minv_L2_back.push_back(L_vector.at(iLambda).M());

              L_L_pT_bin_L1_back.push_back(L_pT_bin_vector_background.at(iLambdaBckg));
              L_L_pT_bin_L2_back.push_back(L_pT_bin_vector.at(iLambda));

              L_L_eta_bin_L1_back.push_back(L_eta_bin_vector_background.at(iLambdaBckg));
              L_L_eta_bin_L2_back.push_back(L_eta_bin_vector.at(iLambda));

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
            if(pBar_vector.at(iLambdaBar).Rapidity() == pBar_vector_background.at(iLambdaBarBckg).Rapidity() ) continue;

            if( nFillLbarLbar_backg % 2 == 0 ) //fill LL pair wit US at x axis and LS at y axis
            {
              L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[Lbar_pT_bin_vector.at(iLambdaBar)][Lbar_pT_bin_vector_background.at(iLambdaBarBckg)]->Fill(Lbar_vector.at(iLambdaBar).M(), Lbar_vector_background.at(iLambdaBarBckg).M());

              //save info for cos(theta*) to be analyzed later
              double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_vector.at(iLambdaBar), pBar_vector.at(iLambdaBar), Lbar_vector_background.at(iLambdaBarBckg), pBar_vector_background.at(iLambdaBarBckg));

              Lbar_Lbar_cos_theta_back.push_back(TMath::Cos(Lbar_Lbar_pairThetaStar));

              Lbar_Lbar_Minv_Lbar1_back.push_back(Lbar_vector.at(iLambdaBar).M());
              Lbar_Lbar_Minv_Lbar2_back.push_back(Lbar_vector_background.at(iLambdaBarBckg).M());

              Lbar_Lbar_pT_bin_Lbar1_back.push_back(Lbar_pT_bin_vector.at(iLambdaBar));
              Lbar_Lbar_pT_bin_Lbar2_back.push_back(Lbar_pT_bin_vector_background.at(iLambdaBarBckg));

              Lbar_Lbar_eta_bin_Lbar1_back.push_back(Lbar_eta_bin_vector.at(iLambdaBar));
              Lbar_Lbar_eta_bin_Lbar2_back.push_back(Lbar_eta_bin_vector_background.at(iLambdaBarBckg));

            }
            else //fill LL pair wit US at y axis and LS at x axis, i.e. opposite order than above
            {
              L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[Lbar_pT_bin_vector_background.at(iLambdaBarBckg)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(Lbar_vector_background.at(iLambdaBarBckg).M(), Lbar_vector.at(iLambdaBar).M());

              //save info for cos(theta*) to be analyzed later
              double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_vector_background.at(iLambdaBarBckg), pBar_vector_background.at(iLambdaBarBckg), Lbar_vector.at(iLambdaBar), pBar_vector.at(iLambdaBar));

              Lbar_Lbar_cos_theta_back.push_back(TMath::Cos(Lbar_Lbar_pairThetaStar));

              Lbar_Lbar_Minv_Lbar1_back.push_back(Lbar_vector_background.at(iLambdaBarBckg).M());
              Lbar_Lbar_Minv_Lbar2_back.push_back(Lbar_vector.at(iLambdaBar).M());

              Lbar_Lbar_pT_bin_Lbar1_back.push_back(Lbar_pT_bin_vector_background.at(iLambdaBarBckg));
              Lbar_Lbar_pT_bin_Lbar2_back.push_back(Lbar_pT_bin_vector.at(iLambdaBar));

              Lbar_Lbar_eta_bin_Lbar1_back.push_back(Lbar_eta_bin_vector_background.at(iLambdaBarBckg));
              Lbar_Lbar_eta_bin_Lbar2_back.push_back(Lbar_eta_bin_vector.at(iLambdaBar));

            }

            nFillLbarLbar_backg++;

          }
        }

      }


      //___________________________________________________________________________________________________

      //fill vectors for mixed event
      //selecting events with only one L or L-bar
      //need to select more - Minv cut will be applied later

      if( L_vector.size() == 1 && Lbar_vector.size() == 0 && L_vector_ME.size() < 1e4)
      {
        L_vector_ME.push_back(L_vector.at(0));
        p_vector_ME.push_back(p_vector.at(0));

        L_pT_bin_vector_ME.push_back(L_pT_bin_vector.at(0));
        L_eta_bin_vector_ME.push_back(L_eta_bin_vector.at(0));

      }

      if( L_vector.size() == 0 && Lbar_vector.size() == 1 && Lbar_vector_ME.size() < 1e4)
      {
        Lbar_vector_ME.push_back(Lbar_vector.at(0));
        pBar_vector_ME.push_back(pBar_vector.at(0));

        Lbar_pT_bin_vector_ME.push_back(Lbar_pT_bin_vector.at(0));
        Lbar_eta_bin_vector_ME.push_back(Lbar_eta_bin_vector.at(0));

      }
      //_________________________________________________________________________________________________________



      //clear all vectors and fill new event


      //clear US vectors
      L_vector.clear();
      L_pT_bin_vector.clear();
      L_eta_bin_vector.clear();

      pi_vector.clear();
      p_vector.clear();


      Lbar_vector.clear();
      Lbar_pT_bin_vector.clear();
      Lbar_eta_bin_vector.clear();

      piBar_vector.clear();
      pBar_vector.clear();

      //fill appropriate vectors from the new event
      if(charge == 0) //US charage combinations
      {
        //fill L or Lbar from new event
        //need to check cuts again in the new event
        if( cuts(L_y) && pT_bin_corr != -1)
        {
          if( p_ch == 1 )
          {

            L_vector.push_back(L_Lorentz_vector);
            L_pT_bin_vector.push_back(pT_bin_corr);
            L_eta_bin_vector.push_back(eta_bin);

            pi_vector.push_back(pi_Lorenz_vector);
            p_vector.push_back(p_Lorenz_vector);
          }
          else if( p_ch == -1 )
          {

            Lbar_vector.push_back(L_Lorentz_vector);
            Lbar_pT_bin_vector.push_back(pT_bin_corr);
            Lbar_eta_bin_vector.push_back(eta_bin);

            piBar_vector.push_back(pi_Lorenz_vector);
            pBar_vector.push_back(p_Lorenz_vector);

          }

        } //end if cuts

      }//end if charge US

      //____________________________________________________________________________________________________________

      //clear LS vectors
      L_vector_background.clear();
      L_pT_bin_vector_background.clear();
      L_eta_bin_vector_background.clear();

      pi_vector_background.clear();
      p_vector_background.clear();


      Lbar_vector_background.clear();
      Lbar_pT_bin_vector_background.clear();
      Lbar_eta_bin_vector_background.clear();

      piBar_vector_background.clear();
      pBar_vector_background.clear();

      if(charge != 0 ) //new event
      {

        //cuts
        if( cuts(L_y) && pT_bin_corr != -1)
        {
          if( p_ch == 1 )
          {
            L_vector_background.push_back(L_Lorentz_vector);
            L_pT_bin_vector_background.push_back(pT_bin_corr);
            L_eta_bin_vector_background.push_back(eta_bin);

            pi_vector_background.push_back(pi_Lorenz_vector);
            p_vector_background.push_back(p_Lorenz_vector);

          }
          else if( p_ch == -1 )
          {
            Lbar_vector_background.push_back(L_Lorentz_vector);
            Lbar_pT_bin_vector_background.push_back(pT_bin_corr);
            Lbar_eta_bin_vector_background.push_back(eta_bin);

            piBar_vector_background.push_back(pi_Lorenz_vector);
            pBar_vector_background.push_back(p_Lorenz_vector);

          }
        }
      }//end if charge LS

    }//end else for new event

  }//end loop over entries in NTuple

  //________________________________________________________________________________________________________

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TH1F *L0_inv_mass_LLbar_projection_US[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_inv_mass_LLbar_projection_US[nPtBins_corr][nPtBins_corr];
  TH1F *L0_inv_mass_LLbar_projection_LS[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_inv_mass_LLbar_projection_LS[nPtBins_corr][nPtBins_corr];
  TH1F *L0_inv_mass_LLbar_projection[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_inv_mass_LLbar_projection[nPtBins_corr][nPtBins_corr];

  TH1F *L01_inv_mass_LL_projection_US[nPtBins_corr][nPtBins_corr];
  TH1F *L02_inv_mass_LL_projection_US[nPtBins_corr][nPtBins_corr];
  TH1F *L01_inv_mass_LL_projection_LS[nPtBins_corr][nPtBins_corr];
  TH1F *L02_inv_mass_LL_projection_LS[nPtBins_corr][nPtBins_corr];
  TH1F *L01_inv_mass_LL_projection[nPtBins_corr][nPtBins_corr];
  TH1F *L02_inv_mass_LL_projection[nPtBins_corr][nPtBins_corr];

  TH1F *L0bar1_inv_mass_LbarLbar_projection_US[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar2_inv_mass_LbarLbar_projection_US[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar1_inv_mass_LbarLbar_projection_LS[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar2_inv_mass_LbarLbar_projection_LS[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar1_inv_mass_LbarLbar_projection[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar2_inv_mass_LbarLbar_projection[nPtBins_corr][nPtBins_corr];

  TFitResultPtr fit_res_gaus_L0_L0Bar[nPtBins_corr][nPtBins_corr];
  TFitResultPtr fit_res_gaus_L0_L0[nPtBins_corr][nPtBins_corr];
  TFitResultPtr fit_res_gaus_L0Bar_L0Bar[nPtBins_corr][nPtBins_corr];

  float sideBandScale_L_Lbar[nPtBins_corr][nPtBins_corr];
  float sideBandScale_L_L[nPtBins_corr][nPtBins_corr];
  float sideBandScale_Lbar_Lbar[nPtBins_corr][nPtBins_corr];

  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      //find bins for side band and for projection to x
      //projection bins just for testing, later can do projections with parameters of 2D fit
      int binLow = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(1.11);
      int binHigh = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(1.125);

      int binLow_sideBand = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(1.14);
      int binHigh_sideBand = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(1.18);

      //L-Lbar

      //US and LS
      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_US_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_US_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_US_can->cd();

      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Draw("surf1");

      L0_inv_mass_vs_L0bar_inv_mass_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_LS_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_LS_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_LS_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_LS_can->cd();

      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Draw("surf1");

      L0_inv_mass_vs_L0bar_inv_mass_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_LS_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_US_LS_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_can->cd();

      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Draw("surf1");

      L0_inv_mass_vs_L0bar_inv_mass_US_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L0_inv_mass_LLbar_projection_US_can = new TCanvas(Form("L0_inv_mass_LLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_LLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_LLbar_projection_US_can->cd();

      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      //L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binHigh, -1); //projection in side-band for testing
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->Draw("p e");

      //scale LS to match US
      sideBandScale_L_Lbar[pTbin1][pTbin2] = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh)/L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh);

      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Sumw2();
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Scale(sideBandScale_L_Lbar[pTbin1][pTbin2]);

      //L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->SetMarkerStyle(24);
      L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->SetLineColor(kBlue);
      L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->Draw("p e same");

      L0_inv_mass_LLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_LLbar_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      TCanvas *L0bar_inv_mass_LLbar_projection_US_can = new TCanvas(Form("L0bar_inv_mass_LLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_LLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_LLbar_projection_US_can->cd();

      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->ProjectionY(Form("L0bar_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->Draw("p e");


      //L0bar_inv_mass_LLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->ProjectionY(Form("L0bar_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar_inv_mass_LLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->ProjectionY(Form("L0bar_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->SetMarkerStyle(24);
      L0bar_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0bar_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->SetLineColor(kBlue);

      L0bar_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->Draw("p e same");

      L0bar_inv_mass_LLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0bar_inv_mass_LLbar_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      //US-LS
      TF2 *doubleGauss_L_Lbar = new TF2("doubleGauss_L_Lbar", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 1.07, 1.2, 1.07, 1.2);
      if( pTbin1 == 0 && pTbin2 == 0 ) doubleGauss_L_Lbar->SetParameters(1200, 1.116, 0.002, 1.116, 0.002);
      else if( fit_res_gaus_L0_L0Bar[0][0]->IsValid() ) doubleGauss_L_Lbar->SetParameters(fit_res_gaus_L0_L0Bar[0][0]->Parameter(0), fit_res_gaus_L0_L0Bar[0][0]->Parameter(1), fit_res_gaus_L0_L0Bar[0][0]->Parameter(2), fit_res_gaus_L0_L0Bar[0][0]->Parameter(3), fit_res_gaus_L0_L0Bar[0][0]->Parameter(4));
      else
      {
        cout<<"Fit not valid for L-Lbar pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
        return false;
      }
/*
      //scale LS to match US
      sideBandScale_L_Lbar[pTbin1][pTbin2] = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh)/L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh);

      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Sumw2();
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Scale(sideBandScale_L_Lbar[pTbin1][pTbin2]);
*/
      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_can->cd();

      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2] = (TH2F*)L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Clone(Form("L0_inv_mass_vs_L0bar_inv_mass_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      //L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Add(L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2], -1);
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Add(L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2], -1);
      fit_res_gaus_L0_L0Bar[pTbin1][pTbin2] = L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_L_Lbar, "s 0", "", 1.11, 1.125);
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Draw("surf1");

      //doubleGauss_L_Lbar->Draw("surf1");

      L0_inv_mass_vs_L0bar_inv_mass_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0_inv_mass_vs_L0bar_inv_mass_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L0_inv_mass_LLbar_projection_can = new TCanvas(Form("L0_inv_mass_LLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_LLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_LLbar_projection_can->cd();

      int LLbar_projectionBin_L = L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(1));

      //L0_inv_mass_LLbar_projection[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LLbar_projectionBin_L, LLbar_projectionBin_L);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->Draw("p e");


      TF12 *LLbar_Gaus_L0 = new TF12("LLbar_Gaus_L0", doubleGauss_L_Lbar, fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(1), "x");
      LLbar_Gaus_L0->Draw("same");

      L0_inv_mass_LLbar_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0_inv_mass_LLbar_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L0bar_inv_mass_LLbar_projection_can = new TCanvas(Form("L0bar_inv_mass_LLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_LLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_LLbar_projection_can->cd();

      int LLbar_projectionBin_Lbar = L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(3));

      //L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0barbar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0bar_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L0bar_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LLbar_projectionBin_Lbar, LLbar_projectionBin_Lbar);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->Draw("p e");


      TF12 *LLbar_Gaus_L0bar = new TF12("LLbar_Gaus_L0bar", doubleGauss_L_Lbar, fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(3), "y");
      LLbar_Gaus_L0bar->Draw("same");

      L0bar_inv_mass_LLbar_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0bar_inv_mass_LLbar_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      //________________________________________________________________________________________________________

      //L-L

      //US
      TCanvas *L0_inv_mass_vs_L0_inv_mass_US_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_US_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_US_can->cd();

      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->Draw("surf1");

      L0_inv_mass_vs_L0_inv_mass_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L0_inv_mass_vs_L0_inv_mass_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      TCanvas *L0_inv_mass_vs_L0_inv_mass_US_LS_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_US_LS_can->cd();

      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->Draw("surf1");

      L0_inv_mass_vs_L0_inv_mass_US_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L0_inv_mass_vs_L0_inv_mass_US_LS_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L01_inv_mass_LL_projection_US_can = new TCanvas(Form("L01_inv_mass_LL_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L01_inv_mass_LL_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L01_inv_mass_LL_projection_US_can->cd();

      L01_inv_mass_LL_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->ProjectionX(Form("L01_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L01_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L01_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L01_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L01_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L01_inv_mass_LL_projection_US[pTbin1][pTbin2]->SetMarkerStyle(20);
      L01_inv_mass_LL_projection_US[pTbin1][pTbin2]->SetMarkerSize(2);
      L01_inv_mass_LL_projection_US[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L01_inv_mass_LL_projection_US[pTbin1][pTbin2]->SetLineColor(kRed);
      L01_inv_mass_LL_projection_US[pTbin1][pTbin2]->Draw("p e");

      //L01_inv_mass_LL_projection_LS[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass_LS[pTbin1][pTbin2]->ProjectionX(Form("L01_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L01_inv_mass_LL_projection_LS[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->ProjectionX(Form("L01_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L01_inv_mass_LL_projection_LS[pTbin1][pTbin2]->SetMarkerStyle(24);
      L01_inv_mass_LL_projection_LS[pTbin1][pTbin2]->SetMarkerSize(2);
      L01_inv_mass_LL_projection_LS[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L01_inv_mass_LL_projection_LS[pTbin1][pTbin2]->SetLineColor(kBlue);
      L01_inv_mass_LL_projection_LS[pTbin1][pTbin2]->Draw("p e same");

      L01_inv_mass_LL_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L01_inv_mass_LL_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      TCanvas *L02_inv_mass_LL_projection_US_can = new TCanvas(Form("L02_inv_mass_LL_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L02_inv_mass_LL_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L02_inv_mass_LL_projection_US_can->cd();

      L02_inv_mass_LL_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->ProjectionY(Form("L02_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L02_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L02_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L02_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L02_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L02_inv_mass_LL_projection_US[pTbin1][pTbin2]->SetMarkerStyle(20);
      L02_inv_mass_LL_projection_US[pTbin1][pTbin2]->SetMarkerSize(2);
      L02_inv_mass_LL_projection_US[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L02_inv_mass_LL_projection_US[pTbin1][pTbin2]->SetLineColor(kRed);
      L02_inv_mass_LL_projection_US[pTbin1][pTbin2]->Draw("p e");

      //L02_inv_mass_LL_projection_LS[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass_LS[pTbin1][pTbin2]->ProjectionY(Form("L02_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L02_inv_mass_LL_projection_LS[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->ProjectionY(Form("L02_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L02_inv_mass_LL_projection_LS[pTbin1][pTbin2]->SetMarkerStyle(24);
      L02_inv_mass_LL_projection_LS[pTbin1][pTbin2]->SetMarkerSize(2);
      L02_inv_mass_LL_projection_LS[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L02_inv_mass_LL_projection_LS[pTbin1][pTbin2]->SetLineColor(kBlue);
      L02_inv_mass_LL_projection_LS[pTbin1][pTbin2]->Draw("p e same");

      L02_inv_mass_LL_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L02_inv_mass_LL_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      //US-LS
      TF2 *doubleGauss_L_L = new TF2("doubleGauss_L_L", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 1.07, 1.2, 1.07, 1.2);
      doubleGauss_L_L->SetParameters(1200, 1.116, 0.002, 1.116, 0.002);
      if( pTbin1 == 0 && pTbin2 == 0 ) doubleGauss_L_L->SetParameters(1200, 1.116, 0.002, 1.116, 0.002);
      else if( fit_res_gaus_L0_L0[0][0]->IsValid() ) doubleGauss_L_L->SetParameters(fit_res_gaus_L0_L0[0][0]->Parameter(0), fit_res_gaus_L0_L0[0][0]->Parameter(1), fit_res_gaus_L0_L0[0][0]->Parameter(2), fit_res_gaus_L0_L0[0][0]->Parameter(3), fit_res_gaus_L0_L0[0][0]->Parameter(4));
      else
      {
        cout<<"Fit not valid for L-L pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
        return false;
      }

      TCanvas *L0_inv_mass_vs_L0_inv_mass_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_can->cd();

      //scale LS to match US
      sideBandScale_L_L[pTbin1][pTbin2] = L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh)/L0_inv_mass_vs_L0_inv_mass_LS[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh);

      L0_inv_mass_vs_L0_inv_mass_LS[pTbin1][pTbin2]->Sumw2();
      L0_inv_mass_vs_L0_inv_mass_LS[pTbin1][pTbin2]->Scale(sideBandScale_L_L[pTbin1][pTbin2]);

      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2] = (TH2F*)L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->Clone(Form("L0_inv_mass_vs_L0_inv_mass_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      //L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->Add(L0_inv_mass_vs_L0_inv_mass_LS[pTbin1][pTbin2], -1);
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->Add(L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2], -1);
      fit_res_gaus_L0_L0[pTbin1][pTbin2] = L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_L_L, "s 0", "", 1.11, 1.125);
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->Draw("surf1");

      L0_inv_mass_vs_L0_inv_mass_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L0_inv_mass_vs_L0_inv_mass_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      TCanvas *L01_inv_mass_LL_projection_can = new TCanvas(Form("L01_inv_mass_LL_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L01_inv_mass_LL_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L01_inv_mass_LL_projection_can->cd();

      int LL_projectionBin_L1 = L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(1));

      //L0_inv_mass_LL_projection[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L01_inv_mass_LL_projection[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L01_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LL_projectionBin_L1, LL_projectionBin_L1);
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->Draw("p e");


      TF12 *LL_Gaus_L0 = new TF12("LL_Gaus_L0", doubleGauss_L_L, fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(1), "x");
      LL_Gaus_L0->Draw("same");

      L01_inv_mass_LL_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L01_inv_mass_LL_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L02_inv_mass_LL_projection_can = new TCanvas(Form("L02_inv_mass_LL_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L02_inv_mass_LL_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L02_inv_mass_LL_projection_can->cd();

      int LL_projectionBin_L2 = L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(3));

      //L02_inv_mass_LL_projection[pTbin1][pTbin2] = (TH1F*)L02_inv_mass_vs_L02bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L02_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L02_inv_mass_LL_projection[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L02_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LL_projectionBin_L2, LL_projectionBin_L2);
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->Draw("p e");


      TF12 *LL_Gaus_L02 = new TF12("LL_Gaus_L02", doubleGauss_L_L, fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(3), "y");
      LL_Gaus_L02->Draw("same");

      L02_inv_mass_LL_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L02_inv_mass_LL_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //________________________________________________________________________________________________________

      //Lbar-Lbar

      //US
      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_US_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_can->cd();

      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Draw("surf1");

      L0bar_inv_mass_vs_L0bar_inv_mass_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));


      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can->cd();

      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Draw("surf1");

      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L0bar1_inv_mass_LbarLbar_projection_US_can = new TCanvas(Form("L0bar1_inv_mass_LbarLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar1_inv_mass_LbarLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar1_inv_mass_LbarLbar_projection_US_can->cd();

      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->ProjectionX(Form("L0bar1_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->Draw("p e");

      //L0bar1_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0bar1_inv_mass_vs_L0bar1_inv_mass_LS[pTbin1][pTbin2]->ProjectionX(Form("L0bar1_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar1_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->ProjectionX(Form("L0bar1_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar1_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->SetMarkerStyle(24);
      L0bar1_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar1_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0bar1_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->SetLineColor(kBlue);
      L0bar1_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->Draw("p e same");

      L0bar1_inv_mass_LbarLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar1_inv_mass_LbarLbar_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L0bar2_inv_mass_LbarLbar_projection_US_can = new TCanvas(Form("L0bar2_inv_mass_LbarLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar2_inv_mass_LbarLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar2_inv_mass_LbarLbar_projection_US_can->cd();

      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->ProjectionX(Form("L0bar2_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->Draw("p e");

      //L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0bar2_inv_mass_vs_L0bar2_inv_mass_LS[pTbin1][pTbin2]->ProjectionX(Form("L0bar2_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->ProjectionX(Form("L0bar2_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->SetMarkerStyle(24);
      L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->SetLineColor(kBlue);
      L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->Draw("p e same");

      L0bar2_inv_mass_LbarLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar2_inv_mass_LbarLbar_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //US-LS
      TF2 *doubleGauss_Lbar_Lbar = new TF2("doubleGauss_Lbar_Lbar", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 1.07, 1.2, 1.07, 1.2);
      doubleGauss_Lbar_Lbar->SetParameters(1200, 1.116, 0.002, 1.116, 0.002);
      if( pTbin1 == 0 && pTbin2 == 0 ) doubleGauss_Lbar_Lbar->SetParameters(1200, 1.116, 0.002, 1.116, 0.002);
      else if( fit_res_gaus_L0Bar_L0Bar[0][0]->IsValid() ) doubleGauss_Lbar_Lbar->SetParameters(200, fit_res_gaus_L0Bar_L0Bar[0][0]->Parameter(1), fit_res_gaus_L0Bar_L0Bar[0][0]->Parameter(2), fit_res_gaus_L0Bar_L0Bar[0][0]->Parameter(3), fit_res_gaus_L0Bar_L0Bar[0][0]->Parameter(4));
      else
      {
        cout<<"Fit not valid for Lbar-Lbar pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
        return false;
      }

      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_can->cd();

      //scale LS to match US
      if( pTbin1 == 0 && pTbin2 == 1 ) sideBandScale_Lbar_Lbar[pTbin1][pTbin2] = L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh)/L0bar_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Integral( binLow_sideBand,binHigh_sideBand, binLow,binHigh);
      else sideBandScale_Lbar_Lbar[pTbin1][pTbin2] = L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Integral(binLow,binHigh, binLow_sideBand,binHigh_sideBand)/L0bar_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Integral( binLow,binHigh, binLow_sideBand,binHigh_sideBand);

      L0bar_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Sumw2();
      L0bar_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Scale(sideBandScale_Lbar_Lbar[pTbin1][pTbin2]);


      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2] = (TH2F*)L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Clone(Form("L0bar_inv_mass_vs_L0bar_inv_mass_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      //L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Add(L0bar_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2], -1);
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Add(L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2], -1);
      //if( pTbin1 == 0 && pTbin2 == 1 ) fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2] = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_Lbar_Lbar, "s 0", "", 1.1, 1.12);
      //else if( pTbin1 == 1 && pTbin2 == 1 ) fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2] = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_Lbar_Lbar, "s 0", "", 1.1, 1.12);
      //else fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2] = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_Lbar_Lbar, "s 0", "", 1.11, 1.125);
      fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2] = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_Lbar_Lbar, "s 0", "", 1.114, 1.125);
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Draw("surf1");

      L0bar_inv_mass_vs_L0bar_inv_mass_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar_inv_mass_vs_L0bar_inv_mass_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L0bar1_inv_mass_LbarLbar_projection_can = new TCanvas(Form("L0bar1_inv_mass_LbarLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar1_inv_mass_LbarLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar1_inv_mass_LbarLbar_projection_can->cd();

      int LbarLbar_projectionBin_Lbar1 = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(1));

      //L0bar_inv_mass_LbarLbar_projection[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0barbar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0bar_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0bar1_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LbarLbar_projectionBin_Lbar1, LbarLbar_projectionBin_Lbar1);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->Draw("p e");


      TF12 *LbarLbar_Gaus_L0bar = new TF12("LbarLbar_Gaus_L0bar", doubleGauss_Lbar_Lbar, fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(1), "x");
      LbarLbar_Gaus_L0bar->Draw("same");

      L0bar1_inv_mass_LbarLbar_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar1_inv_mass_LbarLbar_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));



      TCanvas *L0bar2_inv_mass_LbarLbar_projection_can = new TCanvas(Form("L0bar2_inv_mass_LbarLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar2_inv_mass_LbarLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar2_inv_mass_LbarLbar_projection_can->cd();

      int LbarLbar_projectionBin_Lbar2 = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(3));

      //L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L0bar2_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LbarLbar_projectionBin_Lbar2, LbarLbar_projectionBin_Lbar2);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L0bar2_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LbarLbar_projectionBin_Lbar1, LbarLbar_projectionBin_Lbar1);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->Draw("p e");


      TF12 *LbarLbar_Gaus_L0bar2 = new TF12("LbarLbar_Gaus_L0bar2", doubleGauss_Lbar_Lbar, fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(3), "y");
      LbarLbar_Gaus_L0bar2->Draw("same");

      L0bar2_inv_mass_LbarLbar_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar2_inv_mass_LbarLbar_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //________________________________________________________________________________________________________

    }
  }

  InvMassFile->Close();


  //________________________________________________________________________________________________________________________________________________________________________________

  //analyze stored Lambda pairs and save cos(theta*) histograms

  const float L0_alpha = 0.732; //decay parameter of L0
  const float L0bar_alpha = -0.758; //decay paramteter of L0bar

  TFile *LLbarOutFile; //output file to store production plane histograms

  if(ReadMode == 0) //create production plane file from nTuple - run in this mode first
  {
    LLbarOutFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/ProdPlane_Lambda_work.root", "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT
  }
  else  //read file to store wrong and good sign M_inv distributions and daughter pT
  {

    LLbarOutFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/ProdPlane_Lambda_work.root", "read"); //old, non-optimized cuts

    if( !(LLbarOutFile->IsOpen()) )
    {
      cout<<"Unable to open file with production plane histograms!"<<endl;
      return false;
    }
  }


  //update efficiency files
  TFile *EffFile;

  if(energy == 510) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run17_1B_ME.root", "read");
  else if(energy == 200) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run12_1B_ME.root", "read");
  //if(energy == 510) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run17_1B_ME_tight_eta.root", "read");
  //else if(energy == 200) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run12_1B_ME_tight_eta.root", "read");
  else
  {
    cout<<"Not a valid colliison energy! Abborting!"<<endl;
    return false;
  }

  TFile *EffFileEmbedd;

  if(energy == 510) EffFileEmbedd = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run12_embedd.root", "read");
  else if(energy == 200) EffFileEmbedd = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run12_embedd.root", "read");
  //if(energy == 510) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run17_1B_ME_tight_eta.root", "read");
  //else if(energy == 200) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run12_1B_ME_tight_eta.root", "read");
  else
  {
    cout<<"Not a valid colliison energy! Abborting!"<<endl;
    return false;
  }

  //_______________________________________________________________________________________________________________________________________________

  //efficiency histograms

  TH1F *L0_L0bar_cosThetaProdPlane_eff = (TH1F*)EffFileEmbedd->Get("L0_L0bar_cosThetaProdPlane_RC_hist");
  //TH1F *L0_L0bar_cosThetaProdPlane_eff = (TH1F*)EffFile->Get("L0_L0bar_cosThetaProdPlane_eff");
  TH1F *L0_L0_cosThetaProdPlane_eff = (TH1F*)EffFile->Get("L0_L0_cosThetaProdPlane_eff");
  TH1F *L0bar_L0bar_cosThetaProdPlane_eff = (TH1F*)EffFile->Get("L0bar_L0bar_cosThetaProdPlane_eff");

  TH1F *L0_L0bar_cosThetaProdPlane_pT_eff[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_eff[nEtaBins][nEtaBins];

  TH1F *L0_L0_cosThetaProdPlane_pT_eff[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_eta_eff[nEtaBins][nEtaBins];

  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_eff[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_eff[nEtaBins][nEtaBins];

  //_________________________________________________________________________________________________

  TH1F *L0_L0bar_cosThetaProdPlane_ME_eff = (TH1F*)EffFile->Get("L0_L0bar_cosThetaProdPlane_ME_eff");
  TH1F *L0_L0_cosThetaProdPlane_ME_eff = (TH1F*)EffFile->Get("L0_L0_cosThetaProdPlane_ME_eff");
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_eff = (TH1F*)EffFile->Get("L0bar_L0bar_cosThetaProdPlane_ME_eff");

  TH1F *L0_L0bar_cosThetaProdPlane_ME_pT_eff[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_ME_eta_eff[nEtaBins][nEtaBins];

  TH1F *L0_L0_cosThetaProdPlane_ME_pT_eff[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_ME_eta_eff[nEtaBins][nEtaBins];

  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_pT_eff[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_eta_eff[nEtaBins][nEtaBins];

  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      L0_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2] = (TH1F*)EffFileEmbedd->Get(Form("L0_L0bar_cosThetaProdPlane_RC_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      //L0_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2] = (TH1F*)EffFile->Get(Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0_L0_cosThetaProdPlane_pT_eff[pTbin1][pTbin2] = (TH1F*)EffFile->Get(Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0bar_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2] = (TH1F*)EffFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));

      L0_L0bar_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2] = (TH1F*)EffFile->Get(Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0_L0_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2] = (TH1F*)EffFile->Get(Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0bar_L0bar_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2] = (TH1F*)EffFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      L0_L0bar_cosThetaProdPlane_eta_eff[etaBin1][etaBin2] = (TH1F*)EffFile->Get(Form("L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0_L0_cosThetaProdPlane_eta_eff[etaBin1][etaBin2] = (TH1F*)EffFile->Get(Form("L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0bar_L0bar_cosThetaProdPlane_eta_eff[etaBin1][etaBin2] = (TH1F*)EffFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));

      L0_L0bar_cosThetaProdPlane_ME_eta_eff[etaBin1][etaBin2] = (TH1F*)EffFile->Get(Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0_L0_cosThetaProdPlane_ME_eta_eff[etaBin1][etaBin2] = (TH1F*)EffFile->Get(Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_eff[etaBin1][etaBin2] = (TH1F*)EffFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
    }
  }
  //_____________________________________________________________________________


  //data histograms
  TH1F *L0_L0bar_cosThetaProdPlane_US_hist;
  TH1F *L0_L0bar_cosThetaProdPlane_US_side_band_hist;
  TH1F *L0_L0bar_cosThetaProdPlane_LS_hist;
  TH1F *L0_L0bar_cosThetaProdPlane_ME_hist; //mixed event

  TH1F *L0_L0bar_cosThetaProdPlane_pT_US_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_pT_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0bar_cosThetaProdPlane_pT_ME_hist[nPtBins_corr][nPtBins_corr];


  TH1F *L0_L0bar_cosThetaProdPlane_eta_US_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_LS_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0bar_cosThetaProdPlane_eta_ME_hist[nEtaBins][nEtaBins];

  //_________________________________________

  TH1F *L0_L0_cosThetaProdPlane_US_hist;
  TH1F *L0_L0_cosThetaProdPlane_US_side_band_hist;
  TH1F *L0_L0_cosThetaProdPlane_LS_hist;
  TH1F *L0_L0_cosThetaProdPlane_ME_hist;

  TH1F *L0_L0_cosThetaProdPlane_pT_US_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_pT_US_side_band_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_pT_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0_L0_cosThetaProdPlane_pT_ME_hist[nPtBins_corr][nPtBins_corr];

  TH1F *L0_L0_cosThetaProdPlane_eta_US_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0_cosThetaProdPlane_eta_US_side_band_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0_cosThetaProdPlane_eta_LS_hist[nEtaBins][nEtaBins];
  TH1F *L0_L0_cosThetaProdPlane_eta_ME_hist[nEtaBins][nEtaBins];

  //_________________________________________________

  TH1F *L0bar_L0bar_cosThetaProdPlane_US_hist;
  TH1F *L0bar_L0bar_cosThetaProdPlane_US_side_band_hist;
  TH1F *L0bar_L0bar_cosThetaProdPlane_LS_hist;
  TH1F *L0bar_L0bar_cosThetaProdPlane_ME_hist;

  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_US_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[nPtBins_corr][nPtBins_corr];


  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_US_hist[nEtaBins][nEtaBins];
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[nEtaBins][nEtaBins];
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[nEtaBins][nEtaBins];
  TH1F *L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[nEtaBins][nEtaBins];
  //________________________________________________________________




  L0_L0bar_cosThetaProdPlane_US_hist = new TH1F("L0_L0bar_cosThetaProdPlane_US_hist", "L0_L0bar_cosThetaProdPlane_US_hist", 10, -1, 1);
  L0_L0bar_cosThetaProdPlane_US_side_band_hist = new TH1F("L0_L0bar_cosThetaProdPlane_US_side_band_hist", "L0_L0bar_cosThetaProdPlane_US_side_band_hist", 10, -1, 1);
  L0_L0bar_cosThetaProdPlane_LS_hist = new TH1F("L0_L0bar_cosThetaProdPlane_LS_hist", "L0_L0bar_cosThetaProdPlane_LS_hist", 10, -1, 1);
  L0_L0bar_cosThetaProdPlane_ME_hist = new TH1F("L0_L0bar_cosThetaProdPlane_ME_hist", "L0_L0bar_cosThetaProdPlane_ME_hist", 10, -1, 1);

  L0_L0_cosThetaProdPlane_US_hist = new TH1F("L0_L0_cosThetaProdPlane_US_hist", "L0_L0_cosThetaProdPlane_US_hist", 10, -1, 1);
  L0_L0_cosThetaProdPlane_US_side_band_hist = new TH1F("L0_L0_cosThetaProdPlane_US_side_band_hist", "L0_L0_cosThetaProdPlane_US_side_band_hist", 10, -1, 1);
  L0_L0_cosThetaProdPlane_LS_hist = new TH1F("L0_L0_cosThetaProdPlane_LS_hist", "L0_L0_cosThetaProdPlane_LS_hist", 10, -1, 1);
  L0_L0_cosThetaProdPlane_ME_hist = new TH1F("L0_L0_cosThetaProdPlane_ME_hist", "L0_L0_cosThetaProdPlane_ME_hist", 10, -1, 1);

  L0bar_L0bar_cosThetaProdPlane_US_hist = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_hist", "L0bar_L0bar_cosThetaProdPlane_US_hist", 10, -1, 1);
  L0bar_L0bar_cosThetaProdPlane_US_side_band_hist = new TH1F("L0bar_L0bar_cosThetaProdPlane_US_side_band_hist", "L0bar_L0bar_cosThetaProdPlane_US_side_band_hist", 10, -1, 1);
  L0bar_L0bar_cosThetaProdPlane_LS_hist = new TH1F("L0bar_L0bar_cosThetaProdPlane_LS_hist", "L0bar_L0bar_cosThetaProdPlane_LS_hist", 10, -1, 1);
  L0bar_L0bar_cosThetaProdPlane_ME_hist = new TH1F("L0bar_L0bar_cosThetaProdPlane_ME_hist", "L0bar_L0bar_cosThetaProdPlane_ME_hist", 10, -1, 1);

  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_US_side_band_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_side_band_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_side_band_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_side_band_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_side_band_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_side_band_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
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

      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_US_side_band_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_side_band_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = new TH1F(Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_US_side_band_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_side_band_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = new TH1F(Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
    }
  }

  //________________________________________________________________________________________




  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //Unlike-Sign

  //L-Lbar
  for(unsigned int L_Lbar_index = 0; L_Lbar_index < L_Lbar_cos_theta.size(); L_Lbar_index++)
  {
    float L_peak_mean = fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)]->Parameter(1);
    float L_peak_sigma = fabs(fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)]->Parameter(2));

    float Lbar_peak_mean = fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)]->Parameter(3);
    float Lbar_peak_sigma = fabs(fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)]->Parameter(4));

    if( L_Lbar_Minv_L.at(L_Lbar_index) > L_peak_mean-3*L_peak_sigma && L_Lbar_Minv_L.at(L_Lbar_index) < L_peak_mean+3*L_peak_sigma &&
        L_Lbar_Minv_Lbar.at(L_Lbar_index) > Lbar_peak_mean-3*Lbar_peak_sigma && L_Lbar_Minv_Lbar.at(L_Lbar_index) < Lbar_peak_mean+3*Lbar_peak_sigma)
    {
      L0_L0bar_cosThetaProdPlane_US_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index));

      L0_L0bar_cosThetaProdPlane_pT_US_hist[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)]->Fill(L_Lbar_cos_theta.at(L_Lbar_index));

      L0_L0bar_cosThetaProdPlane_eta_US_hist[L_Lbar_eta_bin_L.at(L_Lbar_index)][L_Lbar_eta_bin_Lbar.at(L_Lbar_index)]->Fill(L_Lbar_cos_theta.at(L_Lbar_index));
    }


    if( L_Lbar_Minv_L.at(L_Lbar_index) > L_peak_mean+5*L_peak_sigma &&  L_Lbar_Minv_Lbar.at(L_Lbar_index) > Lbar_peak_mean+5*Lbar_peak_sigma)
    {
      L0_L0bar_cosThetaProdPlane_US_side_band_hist->Fill(L_Lbar_cos_theta.at(L_Lbar_index));

      L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[L_Lbar_pT_bin_L.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar.at(L_Lbar_index)]->Fill(L_Lbar_cos_theta.at(L_Lbar_index));

      L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[L_Lbar_eta_bin_L.at(L_Lbar_index)][L_Lbar_eta_bin_Lbar.at(L_Lbar_index)]->Fill(L_Lbar_cos_theta.at(L_Lbar_index));
    }


  }



  //L-L
  for(unsigned int L_L_index = 0; L_L_index < L_L_cos_theta.size(); L_L_index++)
  {
    float L1_peak_mean = fit_res_gaus_L0_L0[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)]->Parameter(1);
    float L1_peak_sigma = fabs(fit_res_gaus_L0_L0[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)]->Parameter(2));

    float L2_peak_mean = fit_res_gaus_L0_L0[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)]->Parameter(3);
    float L2_peak_sigma = fabs(fit_res_gaus_L0_L0[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)]->Parameter(4));


    if( L_L_Minv_L1.at(L_L_index) > L1_peak_mean-3*L1_peak_sigma && L_L_Minv_L1.at(L_L_index) < L1_peak_mean+3*L1_peak_sigma &&
        L_L_Minv_L2.at(L_L_index) > L2_peak_mean-3*L2_peak_sigma && L_L_Minv_L2.at(L_L_index) < L2_peak_mean+3*L2_peak_sigma)
    {
      L0_L0_cosThetaProdPlane_US_hist->Fill(L_L_cos_theta.at(L_L_index));

      L0_L0_cosThetaProdPlane_pT_US_hist[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)]->Fill(L_L_cos_theta.at(L_L_index));

      L0_L0_cosThetaProdPlane_eta_US_hist[L_L_eta_bin_L1.at(L_L_index)][L_L_eta_bin_L2.at(L_L_index)]->Fill(L_L_cos_theta.at(L_L_index));
    }


    if(  L_L_Minv_L1.at(L_L_index) > L1_peak_mean+5*L1_peak_sigma && L_L_Minv_L2.at(L_L_index) > L2_peak_mean+5*L2_peak_sigma )
    {
      L0_L0_cosThetaProdPlane_US_side_band_hist->Fill(L_L_cos_theta.at(L_L_index));

      L0_L0_cosThetaProdPlane_pT_US_side_band_hist[L_L_pT_bin_L1.at(L_L_index)][L_L_pT_bin_L2.at(L_L_index)]->Fill(L_L_cos_theta.at(L_L_index));

      L0_L0_cosThetaProdPlane_eta_US_side_band_hist[L_L_eta_bin_L1.at(L_L_index)][L_L_eta_bin_L2.at(L_L_index)]->Fill(L_L_cos_theta.at(L_L_index));
    }


  }



  //Lbar-Lbar
  for(unsigned int Lbar_Lbar_index = 0; Lbar_Lbar_index < Lbar_Lbar_cos_theta.size(); Lbar_Lbar_index++)
  {
    float Lbar1_peak_mean = fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)]->Parameter(1);
    float Lbar1_peak_sigma = fabs(fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)]->Parameter(2));

    float Lbar2_peak_mean = fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)]->Parameter(3);
    float Lbar2_peak_sigma = fabs(fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)]->Parameter(4));


    if( Lbar_Lbar_Minv_Lbar1.at(Lbar_Lbar_index) > Lbar1_peak_mean-3*Lbar1_peak_sigma && Lbar_Lbar_Minv_Lbar1.at(Lbar_Lbar_index) < Lbar1_peak_mean+3*Lbar1_peak_sigma &&
        Lbar_Lbar_Minv_Lbar2.at(Lbar_Lbar_index) > Lbar2_peak_mean-3*Lbar2_peak_sigma && Lbar_Lbar_Minv_Lbar2.at(Lbar_Lbar_index) < Lbar2_peak_mean+3*Lbar2_peak_sigma)
    {
      L0bar_L0bar_cosThetaProdPlane_US_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index));

      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)]->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index));

      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[Lbar_Lbar_eta_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_eta_bin_Lbar2.at(Lbar_Lbar_index)]->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index));
    }


    if( Lbar_Lbar_Minv_Lbar1.at(Lbar_Lbar_index) > Lbar1_peak_mean+5*Lbar1_peak_sigma && Lbar_Lbar_Minv_Lbar2.at(Lbar_Lbar_index) > Lbar2_peak_mean+5*Lbar2_peak_sigma)
    {
      L0bar_L0bar_cosThetaProdPlane_US_side_band_hist->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index));

      L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[Lbar_Lbar_pT_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2.at(Lbar_Lbar_index)]->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index));

      L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[Lbar_Lbar_eta_bin_Lbar1.at(Lbar_Lbar_index)][Lbar_Lbar_eta_bin_Lbar2.at(Lbar_Lbar_index)]->Fill(Lbar_Lbar_cos_theta.at(Lbar_Lbar_index));
    }

  }

  //__________________________________________________________________________________________________



  //background

  //L-Lbar
  for(unsigned int L_Lbar_index = 0; L_Lbar_index < L_Lbar_cos_theta_back.size(); L_Lbar_index++)
  {
    float L_peak_mean = fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)]->Parameter(1);
    float L_peak_sigma = fabs(fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)]->Parameter(2));

    float Lbar_peak_mean = fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)]->Parameter(3);
    float Lbar_peak_sigma = fabs(fit_res_gaus_L0_L0Bar[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)]->Parameter(4));


    //weight for LS - need to scale LS to match realistic sig/bckg ratio known from Minv
    //float LS_weight = sideBandScale_L_Lbar[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)];
    float LS_weight = 1;

    //signal region LS
    if( L_Lbar_Minv_L_back.at(L_Lbar_index) > L_peak_mean-3*L_peak_sigma && L_Lbar_Minv_L_back.at(L_Lbar_index) < L_peak_mean+3*L_peak_sigma &&
        L_Lbar_Minv_Lbar_back.at(L_Lbar_index) > Lbar_peak_mean-3*Lbar_peak_sigma && L_Lbar_Minv_Lbar_back.at(L_Lbar_index) < Lbar_peak_mean+3*Lbar_peak_sigma )
    {
      L0_L0bar_cosThetaProdPlane_LS_hist->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), LS_weight);

      L0_L0bar_cosThetaProdPlane_pT_LS_hist[L_Lbar_pT_bin_L_back.at(L_Lbar_index)][L_Lbar_pT_bin_Lbar_back.at(L_Lbar_index)]->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), LS_weight);

      L0_L0bar_cosThetaProdPlane_eta_LS_hist[L_Lbar_eta_bin_L_back.at(L_Lbar_index)][L_Lbar_eta_bin_Lbar_back.at(L_Lbar_index)]->Fill(L_Lbar_cos_theta_back.at(L_Lbar_index), LS_weight);
    }

  }




  //L-L
  for(unsigned int L_L_index = 0; L_L_index < L_L_cos_theta_back.size(); L_L_index++)
  {
    float L1_peak_mean = fit_res_gaus_L0_L0[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)]->Parameter(1);
    float L1_peak_sigma = fabs(fit_res_gaus_L0_L0[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)]->Parameter(2));

    float L2_peak_mean = fit_res_gaus_L0_L0[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)]->Parameter(3);
    float L2_peak_sigma = fabs(fit_res_gaus_L0_L0[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)]->Parameter(4));


    //float LS_weight = sideBandScale_L_L[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)];
    float LS_weight = 1.;


    if( L_L_Minv_L1_back.at(L_L_index) > L1_peak_mean-3*L1_peak_sigma && L_L_Minv_L1_back.at(L_L_index) < L1_peak_mean+3*L1_peak_sigma &&
        L_L_Minv_L2_back.at(L_L_index) > L2_peak_mean-3*L2_peak_sigma && L_L_Minv_L2_back.at(L_L_index) < L2_peak_mean+3*L2_peak_sigma)
    {
      L0_L0_cosThetaProdPlane_LS_hist->Fill(L_L_cos_theta_back.at(L_L_index), LS_weight);

      L0_L0_cosThetaProdPlane_pT_LS_hist[L_L_pT_bin_L1_back.at(L_L_index)][L_L_pT_bin_L2_back.at(L_L_index)]->Fill(L_L_cos_theta_back.at(L_L_index), LS_weight);

      L0_L0_cosThetaProdPlane_eta_LS_hist[L_L_eta_bin_L1_back.at(L_L_index)][L_L_eta_bin_L2_back.at(L_L_index)]->Fill(L_L_cos_theta_back.at(L_L_index), LS_weight);
    }

  }




  //Lbar-Lbar
  for(unsigned int Lbar_Lbar_index = 0; Lbar_Lbar_index < Lbar_Lbar_cos_theta_back.size(); Lbar_Lbar_index++)
  {
    float Lbar1_peak_mean = fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)]->Parameter(1);
    float Lbar1_peak_sigma = fabs(fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)]->Parameter(2));

    float Lbar2_peak_mean = fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)]->Parameter(3);
    float Lbar2_peak_sigma = fabs(fit_res_gaus_L0Bar_L0Bar[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)]->Parameter(4));


    //float LS_weight = sideBandScale_Lbar_Lbar[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)];
    float LS_weight = 1.;

    if( Lbar_Lbar_Minv_Lbar1_back.at(Lbar_Lbar_index) > Lbar1_peak_mean-3*Lbar1_peak_sigma && Lbar_Lbar_Minv_Lbar1_back.at(Lbar_Lbar_index) < Lbar1_peak_mean+3*Lbar1_peak_sigma &&
        Lbar_Lbar_Minv_Lbar2_back.at(Lbar_Lbar_index) > Lbar2_peak_mean-3*Lbar2_peak_sigma && Lbar_Lbar_Minv_Lbar2_back.at(Lbar_Lbar_index) < Lbar2_peak_mean+3*Lbar2_peak_sigma)
    {
      L0bar_L0bar_cosThetaProdPlane_LS_hist->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), LS_weight);

      L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[Lbar_Lbar_pT_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_pT_bin_Lbar2_back.at(Lbar_Lbar_index)]->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), LS_weight);

      L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[Lbar_Lbar_eta_bin_Lbar1_back.at(Lbar_Lbar_index)][Lbar_Lbar_eta_bin_Lbar2_back.at(Lbar_Lbar_index)]->Fill(Lbar_Lbar_cos_theta_back.at(Lbar_Lbar_index), LS_weight);
    }

  }

  //__________________________________________________________________________________________________



  //analyze mixed-event

  //L-Lbar
  for(unsigned int iLambda_ME = 0; iLambda_ME < L_vector_ME.size(); iLambda_ME++)
  {

    for(unsigned int iLambdaBar_ME = 0; iLambdaBar_ME < Lbar_vector_ME.size(); iLambdaBar_ME++)
    {

      float L_peak_mean = fit_res_gaus_L0_L0Bar[L_pT_bin_vector_ME.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)]->Parameter(1);
      float L_peak_sigma = fabs(fit_res_gaus_L0_L0Bar[L_pT_bin_vector_ME.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)]->Parameter(2));

      float Lbar_peak_mean = fit_res_gaus_L0_L0Bar[L_pT_bin_vector_ME.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)]->Parameter(3);
      float Lbar_peak_sigma = fabs(fit_res_gaus_L0_L0Bar[L_pT_bin_vector_ME.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)]->Parameter(4));

      if( L_vector_ME.at(iLambda_ME).M() < L_peak_mean-3*L_peak_sigma || L_vector_ME.at(iLambda_ME).M() > L_peak_mean+3*L_peak_sigma ) continue;
      if( Lbar_vector_ME.at(iLambda_ME).M() < Lbar_peak_mean-3*Lbar_peak_sigma || Lbar_vector_ME.at(iLambda_ME).M() > Lbar_peak_mean+3*Lbar_peak_sigma ) continue;


      double L_Lbar_pairThetaStar = LpairThetaStar(L_vector_ME.at(iLambda_ME), p_vector_ME.at(iLambda_ME), Lbar_vector_ME.at(iLambdaBar_ME), pBar_vector_ME.at(iLambdaBar_ME));

      L0_L0bar_cosThetaProdPlane_ME_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar));
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[L_pT_bin_vector_ME.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)]->Fill(TMath::Cos(L_Lbar_pairThetaStar));
      //if( fabs(L_vector_ME.at(iLambda_ME).Rapidity()) < 0.2 && fabs(Lbar_vector_ME.at(iLambdaBar_ME).Rapidity()) < 0.2 ) L0_L0bar_cosThetaProdPlane_pT_ME_hist[L_pT_bin_vector_ME.at(iLambda_ME)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME)]->Fill(TMath::Cos(L_Lbar_pairThetaStar));
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[L_eta_bin_vector_ME.at(iLambda_ME)][Lbar_eta_bin_vector_ME.at(iLambdaBar_ME)]->Fill(TMath::Cos(L_Lbar_pairThetaStar));
    }
  }

  //L-L
  for(unsigned int iLambda_ME_1 = 0; iLambda_ME_1 < L_vector_ME.size(); iLambda_ME_1++)
  {
    for(unsigned int iLambda_ME_2 = iLambda_ME_1+1; iLambda_ME_2 < L_vector_ME.size(); iLambda_ME_2++)
    {

      float L1_peak_mean = fit_res_gaus_L0_L0[L_pT_bin_vector_ME.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)]->Parameter(1);
      float L1_peak_sigma = fabs(fit_res_gaus_L0_L0[L_pT_bin_vector_ME.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)]->Parameter(2));

      float L2_peak_mean = fit_res_gaus_L0_L0[L_pT_bin_vector_ME.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)]->Parameter(3);
      float L2_peak_sigma = fabs(fit_res_gaus_L0_L0[L_pT_bin_vector_ME.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)]->Parameter(4));

      if( L_vector_ME.at(iLambda_ME_1).M() < L1_peak_mean-3*L1_peak_sigma || L_vector_ME.at(iLambda_ME_1).M() > L1_peak_mean+3*L1_peak_sigma ) continue;
      if( L_vector_ME.at(iLambda_ME_2).M() < L2_peak_mean-3*L2_peak_sigma || L_vector_ME.at(iLambda_ME_2).M() > L2_peak_mean+3*L2_peak_sigma ) continue;


      double L_L_pairThetaStar = LpairThetaStar(L_vector_ME.at(iLambda_ME_1), p_vector_ME.at(iLambda_ME_1), L_vector_ME.at(iLambda_ME_2), p_vector_ME.at(iLambda_ME_2));

      L0_L0_cosThetaProdPlane_ME_hist->Fill(TMath::Cos(L_L_pairThetaStar));
      L0_L0_cosThetaProdPlane_pT_ME_hist[L_pT_bin_vector_ME.at(iLambda_ME_1)][L_pT_bin_vector_ME.at(iLambda_ME_2)]->Fill(TMath::Cos(L_L_pairThetaStar));
      L0_L0_cosThetaProdPlane_eta_ME_hist[L_eta_bin_vector_ME.at(iLambda_ME_1)][L_eta_bin_vector_ME.at(iLambda_ME_2)]->Fill(TMath::Cos(L_L_pairThetaStar));
    }
  }


  //Lbar-Lbar
  for(unsigned int iLambdaBar_ME_1 = 0; iLambdaBar_ME_1 < Lbar_vector_ME.size(); iLambdaBar_ME_1++)
  {
    for(unsigned int iLambdaBar_ME_2 = iLambdaBar_ME_1+1; iLambdaBar_ME_2 < Lbar_vector_ME.size(); iLambdaBar_ME_2++)
    {

      float Lbar1_peak_mean = fit_res_gaus_L0Bar_L0Bar[Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)]->Parameter(1);
      float Lbar1_peak_sigma = fabs(fit_res_gaus_L0Bar_L0Bar[Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)]->Parameter(2));

      float Lbar2_peak_mean = fit_res_gaus_L0Bar_L0Bar[Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)]->Parameter(3);
      float Lbar2_peak_sigma = fabs(fit_res_gaus_L0Bar_L0Bar[Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)]->Parameter(4));

      if( Lbar_vector_ME.at(iLambdaBar_ME_1).M() < Lbar1_peak_mean-3*Lbar1_peak_sigma || Lbar_vector_ME.at(iLambdaBar_ME_1).M() > Lbar1_peak_mean+3*Lbar1_peak_sigma ) continue;
      if( Lbar_vector_ME.at(iLambdaBar_ME_2).M() < Lbar2_peak_mean-3*Lbar2_peak_sigma || Lbar_vector_ME.at(iLambdaBar_ME_2).M() > Lbar2_peak_mean+3*Lbar2_peak_sigma ) continue;

      double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_vector_ME.at(iLambdaBar_ME_1), pBar_vector_ME.at(iLambdaBar_ME_1), Lbar_vector_ME.at(iLambdaBar_ME_2), pBar_vector_ME.at(iLambdaBar_ME_2));

      L0bar_L0bar_cosThetaProdPlane_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar));
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar));
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[Lbar_eta_bin_vector_ME.at(iLambdaBar_ME_1)][Lbar_eta_bin_vector_ME.at(iLambdaBar_ME_2)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar));
    }
  }

  //________________________________________________________________________________________________________





  TCanvas *L0_L0bar_cosThetaProdPlane_no_corr_can = new TCanvas("L0_L0bar_cosThetaProdPlane_no_corr_can", "L0_L0bar_cosThetaProdPlane_no_corr_can", 1200, 1000);

  L0_L0bar_cosThetaProdPlane_no_corr_can->cd();

  L0_L0bar_cosThetaProdPlane_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_US_hist->GetYaxis()->SetMaxDigits(3);
  L0_L0bar_cosThetaProdPlane_US_hist->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_US_hist->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_US_hist->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane_US_hist->SetLineColor(kRed);
  double nLLbar = L0_L0bar_cosThetaProdPlane_US_hist->Integral();
  L0_L0bar_cosThetaProdPlane_US_hist->Sumw2();
  //L0_L0bar_cosThetaProdPlane_US_hist->Divide(L0_L0bar_cosThetaProdPlane_eff);
  L0_L0bar_cosThetaProdPlane_US_hist->Scale(1./L0_L0bar_cosThetaProdPlane_US_hist->GetXaxis()->GetBinWidth(1));
  //L0_L0bar_cosThetaProdPlane_US_hist->Scale(1./L0_L0bar_cosThetaProdPlane_US_hist->Integral());
  L0_L0bar_cosThetaProdPlane_US_hist->SetMinimum(0);
  L0_L0bar_cosThetaProdPlane_US_hist->Draw("p e");

  L0_L0bar_cosThetaProdPlane_US_side_band_hist->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_US_side_band_hist->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_US_side_band_hist->SetMarkerColor(kBlue);
  L0_L0bar_cosThetaProdPlane_US_side_band_hist->Sumw2();
  L0_L0bar_cosThetaProdPlane_US_side_band_hist->Scale(L0_L0bar_cosThetaProdPlane_LS_hist->Integral()/L0_L0bar_cosThetaProdPlane_US_side_band_hist->Integral()); //scale side band backgroun to match combinatorial bckg. levels
  L0_L0bar_cosThetaProdPlane_US_side_band_hist->Scale(1./L0_L0bar_cosThetaProdPlane_US_side_band_hist->GetXaxis()->GetBinWidth(1));
  //L0_L0bar_cosThetaProdPlane_US_side_band_hist->Draw("p e same");

  L0_L0bar_cosThetaProdPlane_LS_hist->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_LS_hist->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_LS_hist->SetMarkerColor(kBlue);
  double nLLbar_back = L0_L0bar_cosThetaProdPlane_LS_hist->Integral();
  L0_L0bar_cosThetaProdPlane_LS_hist->Sumw2();
  L0_L0bar_cosThetaProdPlane_LS_hist->Scale(1./L0_L0bar_cosThetaProdPlane_LS_hist->GetXaxis()->GetBinWidth(1));
  //L0_L0bar_cosThetaProdPlane_LS_hist->Divide(L0_L0bar_cosThetaProdPlane_eff);
  L0_L0bar_cosThetaProdPlane_LS_hist->Draw("p e same");

  L0_L0bar_cosThetaProdPlane_ME_hist->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_ME_hist->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_ME_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_ME_hist->GetYaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_ME_hist->GetYaxis()->SetMaxDigits(3);
  L0_L0bar_cosThetaProdPlane_ME_hist->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_ME_hist->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_ME_hist->SetMarkerColor(kBlue);
  L0_L0bar_cosThetaProdPlane_ME_hist->SetLineColor(kBlue);
  L0_L0bar_cosThetaProdPlane_ME_hist->Sumw2();
  //L0_L0bar_cosThetaProdPlane_ME_hist->Divide(L0_L0bar_cosThetaProdPlane_ME_eff);
  L0_L0bar_cosThetaProdPlane_ME_hist->Scale(nLLbar_back/L0_L0bar_cosThetaProdPlane_ME_hist->Integral()); //scale ME to expected background levels using LS
  L0_L0bar_cosThetaProdPlane_ME_hist->Scale(1./L0_L0bar_cosThetaProdPlane_ME_hist->GetXaxis()->GetBinWidth(1));
  L0_L0bar_cosThetaProdPlane_ME_hist->SetMinimum(0);
  //L0_L0bar_cosThetaProdPlane_ME_hist->Draw("same p e");


  TF1 *fitL0_L0bar_US_ThetaStar_no_corr = new TF1("fitL0_L0bar_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
  fitL0_L0bar_US_ThetaStar_no_corr->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0_L0bar_cosThetaProdPlane_US_hist->Fit(fitL0_L0bar_US_ThetaStar_no_corr, "s i 0 r");

  float P_L0_L0bar_no_corr = fitL0_L0bar_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
  float P_L0_L0bar_no_corr_err = fitL0_L0bar_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

  fitL0_L0bar_US_ThetaStar_no_corr->SetLineColor(1);
  fitL0_L0bar_US_ThetaStar_no_corr->Draw("same");

  TLegend *L0_L0bar_leg = new TLegend(0.15, 0.4, 0.39, 0.69);
  L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_US_hist, "Unlike-sign p#pi");
  //L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_US_side_band_hist, "Side-band p#pi");
  L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_LS_hist, "Combinatorial bckg.");
  //L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_ME_hist, "Mixed event");
  L0_L0bar_leg->AddEntry(fitL0_L0bar_US_ThetaStar_no_corr, "Linear fit to US");
  L0_L0bar_leg->SetBorderSize(0);
  L0_L0bar_leg->SetFillColorAlpha(0, 0.01);
  L0_L0bar_leg->Draw("same");

  TPaveText *L0_L0bar_text_no_corr = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
  L0_L0bar_text_no_corr->SetTextFont(42);
  //L0_L0bar_text_no_corr->AddText("STAR Internal");
  //L0_L0bar_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0bar_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0_L0bar_text_no_corr->AddText("Minimum bias, no correction");
  L0_L0bar_text_no_corr->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
  L0_L0bar_text_no_corr->AddText("|#it{y}| < 1");
  L0_L0bar_text_no_corr->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L0_L0bar_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar_no_corr, P_L0_L0bar_no_corr_err));
  L0_L0bar_text_no_corr->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_no_corr->Draw("same");

  L0_L0bar_cosThetaProdPlane_no_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0_L0bar_cosThetaProdPlane_no_corr.png");
  //_____________________________________________________________________________________

  TCanvas *L0_L0bar_cosThetaProdPlane_can = new TCanvas("L0_L0bar_cosThetaProdPlane_can", "L0_L0bar_cosThetaProdPlane_can", 1200, 1000);

  L0_L0bar_cosThetaProdPlane_can->cd();

  L0_L0bar_cosThetaProdPlane_US_hist->Divide(L0_L0bar_cosThetaProdPlane_eff);
  //L0_L0bar_cosThetaProdPlane_US_hist->Scale(1./L0_L0bar_cosThetaProdPlane_US_hist->Integral());
  L0_L0bar_cosThetaProdPlane_US_hist->SetMinimum(0);
  L0_L0bar_cosThetaProdPlane_US_hist->Draw("p e");

  L0_L0bar_cosThetaProdPlane_ME_hist->Divide(L0_L0bar_cosThetaProdPlane_ME_eff);
  //L0_L0bar_cosThetaProdPlane_ME_hist->Scale(1./L0_L0bar_cosThetaProdPlane_ME_hist->Integral());
  //L0_L0bar_cosThetaProdPlane_ME_hist->Draw("same p e");

  TF1 *fitL0_L0bar_US_ThetaStar = new TF1("fitL0_L0bar_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fitL0_L0bar_US_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0_L0bar_cosThetaProdPlane_US_hist->Fit(fitL0_L0bar_US_ThetaStar, "s i 0 r");

  float P_L0_L0bar = fitL0_L0bar_US_ThetaStar->GetParameter(1)/(L0_alpha*L0bar_alpha);
  float P_L0_L0bar_err = fitL0_L0bar_US_ThetaStar->GetParError(1)/(L0_alpha*L0bar_alpha);

  fitL0_L0bar_US_ThetaStar->SetLineColor(1);
  fitL0_L0bar_US_ThetaStar->Draw("same");


  TPaveText *L0_L0bar_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
  L0_L0bar_text->SetTextFont(42);
  //L0_L0bar_text->AddText("STAR Internal");
  //L0_L0bar_text->AddText("STAR preliminary");
  //((TText*)L0_L0bar_text->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0bar_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0_L0bar_text->AddText("Minimum bias");
  L0_L0bar_text->AddText("JP2");
  L0_L0bar_text->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
  L0_L0bar_text->AddText("|#it{y}| < 1");
  L0_L0bar_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L0_L0bar_text->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar, P_L0_L0bar_err));
  L0_L0bar_text->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text->Draw("same");

  L0_L0bar_leg->Draw("same");

  L0_L0bar_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0_L0bar_cosThetaProdPlane.png");
  //____________________________________________________________________


  TCanvas *L0_L0_cosThetaProdPlane_no_corr_can = new TCanvas("L0_L0_cosThetaProdPlane_no_corr_can", "L0_L0_cosThetaProdPlane_no_corr_can", 1200, 1000);

  L0_L0_cosThetaProdPlane_no_corr_can->cd();

  L0_L0_cosThetaProdPlane_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_US_hist->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0_cosThetaProdPlane_US_hist->GetYaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_US_hist->SetMarkerSize(1.5);
  L0_L0_cosThetaProdPlane_US_hist->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_US_hist->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane_US_hist->SetLineColor(kRed);
  double nLL = L0_L0_cosThetaProdPlane_US_hist->Integral();
  L0_L0_cosThetaProdPlane_US_hist->Sumw2();
  //L0_L0_cosThetaProdPlane_US_hist->Divide(L0_L0_cosThetaProdPlane_eff);
  L0_L0_cosThetaProdPlane_US_hist->Scale(1./L0_L0_cosThetaProdPlane_US_hist->GetXaxis()->GetBinWidth(1));
  //L0_L0_cosThetaProdPlane_US_hist->Scale(1./L0_L0_cosThetaProdPlane_US_hist->Integral());
  L0_L0_cosThetaProdPlane_US_hist->SetMinimum(0);
  L0_L0_cosThetaProdPlane_US_hist->Draw("p e");

  L0_L0_cosThetaProdPlane_US_side_band_hist->SetMarkerSize(1.5);
  L0_L0_cosThetaProdPlane_US_side_band_hist->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_US_side_band_hist->SetMarkerColor(kBlue);
  L0_L0_cosThetaProdPlane_US_side_band_hist->Sumw2();
  L0_L0_cosThetaProdPlane_US_side_band_hist->Scale(L0_L0_cosThetaProdPlane_LS_hist->Integral()/L0_L0_cosThetaProdPlane_US_side_band_hist->Integral());
  L0_L0_cosThetaProdPlane_US_side_band_hist->Scale(1./L0_L0_cosThetaProdPlane_US_side_band_hist->GetXaxis()->GetBinWidth(1));
  //L0_L0_cosThetaProdPlane_US_side_band_hist->Draw("p e same");

  L0_L0_cosThetaProdPlane_LS_hist->SetMarkerSize(1.5);
  L0_L0_cosThetaProdPlane_LS_hist->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_LS_hist->SetMarkerColor(kBlue);
  double nLL_back = L0_L0_cosThetaProdPlane_LS_hist->Integral();
  L0_L0_cosThetaProdPlane_LS_hist->Sumw2();
  L0_L0_cosThetaProdPlane_LS_hist->Scale(1./L0_L0_cosThetaProdPlane_LS_hist->GetXaxis()->GetBinWidth(1));
  //L0_L0_cosThetaProdPlane_LS_hist->Divide(L0_L0_cosThetaProdPlane_eff);
  L0_L0_cosThetaProdPlane_LS_hist->Draw("p e same");

  L0_L0_cosThetaProdPlane_ME_hist->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_ME_hist->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_ME_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0_L0_cosThetaProdPlane_ME_hist->GetYaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_ME_hist->GetYaxis()->SetMaxDigits(3);
  L0_L0_cosThetaProdPlane_ME_hist->SetMarkerSize(1.5);
  L0_L0_cosThetaProdPlane_ME_hist->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_ME_hist->SetMarkerColor(kBlue);
  L0_L0_cosThetaProdPlane_ME_hist->SetLineColor(kBlue);
  L0_L0_cosThetaProdPlane_ME_hist->Sumw2();
  //L0_L0_cosThetaProdPlane_ME_hist->Divide(L0_L0_cosThetaProdPlane_ME_eff);
  L0_L0_cosThetaProdPlane_ME_hist->Scale(nLL_back/L0_L0_cosThetaProdPlane_ME_hist->Integral());
  L0_L0_cosThetaProdPlane_ME_hist->Scale(1./L0_L0_cosThetaProdPlane_ME_hist->GetXaxis()->GetBinWidth(1));
  L0_L0_cosThetaProdPlane_ME_hist->SetMinimum(0);
  //L0_L0_cosThetaProdPlane_ME_hist->Draw("same p e");


  TF1 *fitL0_L0_US_ThetaStar_no_corr = new TF1("fitL0_L0_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
  fitL0_L0_US_ThetaStar_no_corr->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0_L0_cosThetaProdPlane_US_hist->Fit(fitL0_L0_US_ThetaStar_no_corr, "s i 0 r");

  float P_L0_L0_no_corr = fitL0_L0_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_L0_L0_no_corr_err = fitL0_L0_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0_alpha);

  fitL0_L0_US_ThetaStar_no_corr->SetLineColor(1);
  fitL0_L0_US_ThetaStar_no_corr->Draw("same");

  TLegend *L0_L0_leg = new TLegend(0.15, 0.4, 0.39, 0.69);
  L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_US_hist, "Unlike-sign p#pi");
  //L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_US_side_band_hist, "Side-band p#pi");
  L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_LS_hist, "Combinatorial bckg.");
  //L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_ME_hist, "Mixed event");
  L0_L0_leg->AddEntry(fitL0_L0_US_ThetaStar_no_corr, "Linear fit to US");
  L0_L0_leg->SetBorderSize(0);
  L0_L0_leg->SetFillColorAlpha(0, 0.01);
  L0_L0_leg->Draw("same");

  TPaveText *L0_L0_text_no_corr = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
  L0_L0_text_no_corr->SetTextFont(42);
  //L0_L0_text_no_corr->AddText("STAR Internal");
  //L0_L0_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0_L0_text_no_corr->AddText("Minimum bias, no correction");
  L0_L0_text_no_corr->AddText("#Lambda^{0}-#Lambda^{0}");
  L0_L0_text_no_corr->AddText("|#it{y}| < 1");
  L0_L0_text_no_corr->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L0_L0_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0_no_corr, P_L0_L0_no_corr_err));
  L0_L0_text_no_corr->SetFillColorAlpha(0, 0.01);
  L0_L0_text_no_corr->Draw("same");

  L0_L0_cosThetaProdPlane_no_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0_L0_cosThetaProdPlane_no_corr.png");
  //________________________________________________________________________

  TCanvas *L0_L0_cosThetaProdPlane_can = new TCanvas("L0_L0_cosThetaProdPlane_can", "L0_L0_cosThetaProdPlane_can", 1200, 1000);

  L0_L0_cosThetaProdPlane_can->cd();

  L0_L0_cosThetaProdPlane_US_hist->Divide(L0_L0_cosThetaProdPlane_eff);
  //L0_L0_cosThetaProdPlane_US_hist->Scale(1./L0_L0_cosThetaProdPlane_US_hist->Integral());
  L0_L0_cosThetaProdPlane_US_hist->SetMinimum(0);
  L0_L0_cosThetaProdPlane_US_hist->Draw("p e");

  L0_L0_cosThetaProdPlane_ME_hist->Divide(L0_L0_cosThetaProdPlane_ME_eff);
  //L0_L0_cosThetaProdPlane_ME_hist->Scale(1./L0_L0_cosThetaProdPlane_ME_hist->Integral());
  //L0_L0_cosThetaProdPlane_ME_hist->Draw("same p e");

  TF1 *fitL0_L0_US_ThetaStar = new TF1("fitL0_L0_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fitL0_L0_US_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0_L0_cosThetaProdPlane_US_hist->Fit(fitL0_L0_US_ThetaStar, "s i 0 r");

  float P_L0_L0 = fitL0_L0_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_L0_L0_err = fitL0_L0_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

  fitL0_L0_US_ThetaStar->SetLineColor(1);
  fitL0_L0_US_ThetaStar->Draw("same");

  TPaveText *L0_L0_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
  L0_L0_text->SetTextFont(42);
  //L0_L0_text->AddText("STAR Internal");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0_L0_text->AddText("Minimum bias");
  L0_L0_text->AddText("JP2");
  L0_L0_text->AddText("#Lambda^{0}-#Lambda^{0}");
  L0_L0_text->AddText("|#it{y}| < 1");
  L0_L0_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L0_L0_text->AddText(Form("P = %.2f #pm %.2f", P_L0_L0, P_L0_L0_err));
  L0_L0_text->SetFillColorAlpha(0, 0.01);
  L0_L0_text->Draw("same");

  L0_L0_leg->Draw("same");

  L0_L0_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0_L0_cosThetaProdPlane.png");
  //____________________________________________________________________


  TCanvas *L0bar_L0bar_cosThetaProdPlane_no_corr_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_no_corr_can", "L0bar_L0bar_cosThetaProdPlane_no_corr_can", 1200, 1000);

  L0bar_L0bar_cosThetaProdPlane_no_corr_can->cd();

  L0bar_L0bar_cosThetaProdPlane_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_US_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_US_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_US_hist->SetMarkerSize(1.5);
  L0bar_L0bar_cosThetaProdPlane_US_hist->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_US_hist->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_US_hist->SetLineColor(kRed);
  double nLbarLbar = L0bar_L0bar_cosThetaProdPlane_US_hist->Integral();
  L0bar_L0bar_cosThetaProdPlane_US_hist->Sumw2();
  //L0bar_L0bar_cosThetaProdPlane_US_hist->Divide(L0bar_L0bar_cosThetaProdPlane_ME_eff);
  L0bar_L0bar_cosThetaProdPlane_US_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_hist->GetXaxis()->GetBinWidth(1));
  //L0bar_L0bar_cosThetaProdPlane_US_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_hist->Integral());
  L0bar_L0bar_cosThetaProdPlane_US_hist->SetMinimum(0);
  L0bar_L0bar_cosThetaProdPlane_US_hist->Draw("p e");

  L0bar_L0bar_cosThetaProdPlane_US_side_band_hist->SetMarkerSize(1.5);
  L0bar_L0bar_cosThetaProdPlane_US_side_band_hist->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_US_side_band_hist->SetMarkerColor(kBlue);
  L0bar_L0bar_cosThetaProdPlane_US_side_band_hist->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_US_side_band_hist->Scale(L0bar_L0bar_cosThetaProdPlane_LS_hist->Integral()/L0bar_L0bar_cosThetaProdPlane_US_side_band_hist->Integral());
  L0bar_L0bar_cosThetaProdPlane_US_side_band_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_side_band_hist->GetXaxis()->GetBinWidth(1));
  //L0bar_L0bar_cosThetaProdPlane_US_side_band_hist->Draw("p e same");

  L0bar_L0bar_cosThetaProdPlane_LS_hist->SetMarkerSize(1.5);
  L0bar_L0bar_cosThetaProdPlane_LS_hist->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_LS_hist->SetMarkerColor(kBlue);
  double nLbarLbar_back = L0bar_L0bar_cosThetaProdPlane_LS_hist->Integral();
  L0bar_L0bar_cosThetaProdPlane_LS_hist->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_LS_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_LS_hist->GetXaxis()->GetBinWidth(1));
  //L0bar_L0bar_cosThetaProdPlane_LS_hist->Divide(L0bar_L0bar_cosThetaProdPlane_eff);
  L0bar_L0bar_cosThetaProdPlane_LS_hist->Draw("p e same");

  L0bar_L0bar_cosThetaProdPlane_ME_hist->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_ME_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_ME_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_ME_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_ME_hist->GetYaxis()->SetMaxDigits(3);
  L0bar_L0bar_cosThetaProdPlane_ME_hist->SetMarkerSize(1.5);
  L0bar_L0bar_cosThetaProdPlane_ME_hist->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_ME_hist->SetMarkerColor(kBlue);
  L0bar_L0bar_cosThetaProdPlane_ME_hist->SetLineColor(kBlue);
  L0bar_L0bar_cosThetaProdPlane_ME_hist->Sumw2();
  //L0bar_L0bar_cosThetaProdPlane_ME_hist->Divide(L0bar_L0bar_cosThetaProdPlane_ME_eff);
  L0bar_L0bar_cosThetaProdPlane_ME_hist->Scale(nLbarLbar_back/L0bar_L0bar_cosThetaProdPlane_ME_hist->Integral());
  L0bar_L0bar_cosThetaProdPlane_ME_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_ME_hist->GetXaxis()->GetBinWidth(1));
  L0bar_L0bar_cosThetaProdPlane_ME_hist->SetMinimum(0);
  //L0bar_L0bar_cosThetaProdPlane_ME_hist->Draw("same p e");


  TF1 *fitL0bar_L0bar_US_ThetaStar_no_corr = new TF1("fitL0bar_L0bar_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
  fitL0bar_L0bar_US_ThetaStar_no_corr->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0bar_L0bar_cosThetaProdPlane_US_hist->Fit(fitL0bar_L0bar_US_ThetaStar_no_corr, "s i 0 r");

  float P_L0bar_L0bar_no_corr = fitL0bar_L0bar_US_ThetaStar_no_corr->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
  float P_L0bar_L0bar_no_corr_err = fitL0bar_L0bar_US_ThetaStar_no_corr->GetParError(1)/(L0bar_alpha*L0bar_alpha);

  fitL0bar_L0bar_US_ThetaStar_no_corr->SetLineColor(1);
  fitL0bar_L0bar_US_ThetaStar_no_corr->Draw("same");


  TLegend *L0bar_L0bar_leg = new TLegend(0.15, 0.4, 0.39, 0.69);
  L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_US_hist, "Unlike-sign p#pi");
  //L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_US_side_band_hist, "Side-band p#pi");
  L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_LS_hist, "Combinatorial bckg.");
  //L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_ME_hist, "Mixed event");
  L0bar_L0bar_leg->AddEntry(fitL0bar_L0bar_US_ThetaStar_no_corr, "Linear fit to US");
  L0bar_L0bar_leg->SetBorderSize(0);
  L0bar_L0bar_leg->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_leg->Draw("same");


  TPaveText *L0bar_L0bar_text_no_corr = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
  L0bar_L0bar_text_no_corr->SetTextFont(42);
  //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
  //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0bar_L0bar_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0bar_L0bar_text_no_corr->AddText("Minimum bias, no correction");
  L0bar_L0bar_text_no_corr->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
  L0bar_L0bar_text_no_corr->AddText("|#it{y}| < 1");
  L0bar_L0bar_text_no_corr->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L0bar_L0bar_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar_no_corr, P_L0bar_L0bar_no_corr_err));
  L0bar_L0bar_text_no_corr->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_no_corr->Draw("same");


  L0bar_L0bar_cosThetaProdPlane_no_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0bar_L0bar_cosThetaProdPlane_no_corr.png");
  //_______________________________________________________________________________________

  TCanvas *L0bar_L0bar_cosThetaProdPlane_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_can", "L0bar_L0bar_cosThetaProdPlane_can", 1200, 1000);

  L0bar_L0bar_cosThetaProdPlane_can->cd();

  L0bar_L0bar_cosThetaProdPlane_US_hist->Divide(L0bar_L0bar_cosThetaProdPlane_eff);
  //L0bar_L0bar_cosThetaProdPlane_US_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_hist->Integral());
  L0bar_L0bar_cosThetaProdPlane_US_hist->SetMinimum(0);
  L0bar_L0bar_cosThetaProdPlane_US_hist->Draw("p e");

  L0bar_L0bar_cosThetaProdPlane_ME_hist->Divide(L0bar_L0bar_cosThetaProdPlane_ME_eff);
  //L0bar_L0bar_cosThetaProdPlane_ME_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_ME_hist->Integral());
  //L0bar_L0bar_cosThetaProdPlane_ME_hist->Draw("same p e");

    TF1 *fitL0bar_L0bar_US_ThetaStar = new TF1("fitL0bar_L0bar_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fitL0bar_L0bar_US_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0bar_L0bar_cosThetaProdPlane_US_hist->Fit(fitL0bar_L0bar_US_ThetaStar, "s i 0 r");

  float P_L0bar_L0bar = fitL0bar_L0bar_US_ThetaStar->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
  float P_L0bar_L0bar_err = fitL0bar_L0bar_US_ThetaStar->GetParError(1)/(L0bar_alpha*L0bar_alpha);

  fitL0bar_L0bar_US_ThetaStar->SetLineColor(1);
  fitL0bar_L0bar_US_ThetaStar->Draw("same");

  TPaveText *L0bar_L0bar_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
  L0bar_L0bar_text->SetTextFont(42);
  //L0bar_L0bar_text->AddText("STAR Internal");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  L0bar_L0bar_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0bar_L0bar_text->AddText("Minimum bias");
  L0bar_L0bar_text->AddText("JP2");
  L0bar_L0bar_text->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
  L0bar_L0bar_text->AddText("|#it{y}| < 1");
  L0bar_L0bar_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L0bar_L0bar_text->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar, P_L0bar_L0bar_err));
  L0bar_L0bar_text->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text->Draw("same");


  L0bar_L0bar_leg->Draw("same");


  L0bar_L0bar_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0bar_L0bar_cosThetaProdPlane.png");
  //____________________________________________________________________

  double nLLbar_pT[nPtBins_corr][nPtBins_corr];
  double nLL_pT[nPtBins_corr][nPtBins_corr];
  double nLbarLbar_pT[nPtBins_corr][nPtBins_corr];

  double nLLbar_pT_back[nPtBins_corr][nPtBins_corr];
  double nLL_pT_back[nPtBins_corr][nPtBins_corr];
  double nLbarLbar_pT_back[nPtBins_corr][nPtBins_corr];


  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      TCanvas *L0_L0bar_cosThetaProdPlane_pT_no_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0_L0bar_cosThetaProdPlane_pT_no_corr_can->cd();

      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      nLLbar_pT[pTbin1][pTbin2] = L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral();
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Sumw2();
      //L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      //L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral());
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Draw("p e");

      L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      nLLbar_pT_back[pTbin1][pTbin2] = L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Integral();
      L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Scale(L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral()/L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Integral());
      L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Draw("p e same");

      L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      nLLbar_pT_back[pTbin1][pTbin2] = L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral();
      L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Sumw2();
      //L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      //L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      //L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Draw("p e same");

      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetLineColor(kBlue);
      //double nLLbar = L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral();
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Sumw2();
      //L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2]);
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral()/L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral());
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMinimum(0);
      //L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Draw("same p e");


      TF1 *fitL0_L0bar_pT_US_ThetaStar_no_corr = new TF1("fitL0_L0bar_pT_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_pT_US_ThetaStar_no_corr->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitL0_L0bar_pT_US_ThetaStar_no_corr, "s i 0 r");

      float P_L0_L0bar_pT_no_corr = fitL0_L0bar_pT_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
      float P_L0_L0bar_pT_no_corr_err = fitL0_L0bar_pT_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

      fitL0_L0bar_pT_US_ThetaStar_no_corr->SetLineColor(1);
      fitL0_L0bar_pT_US_ThetaStar_no_corr->Draw("same");

      TPaveText *L0_L0bar_pT_text_no_corr = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      L0_L0bar_pT_text_no_corr->SetTextFont(42);
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0bar_pT_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_pT_text_no_corr->AddText("Minimum bias, no correction");
      L0_L0bar_pT_text_no_corr->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      //L0_L0bar_pT_text_no_corr->AddText("|#it{y}| < 0.2");
      L0_L0bar_pT_text_no_corr->AddText("|#it{y}| < 1");
      L0_L0bar_pT_text_no_corr->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0bar_pT_text_no_corr->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0bar_pT_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar_pT_no_corr, P_L0_L0bar_pT_no_corr_err));
      L0_L0bar_pT_text_no_corr->SetFillColorAlpha(0, 0.01);
      L0_L0bar_pT_text_no_corr->Draw("same");

      L0_L0bar_leg->Draw("same");


      L0_L0bar_cosThetaProdPlane_pT_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0_L0bar_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //_______________________________________________



      TCanvas *L0_L0bar_cosThetaProdPlane_pT_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0_L0bar_cosThetaProdPlane_pT_can->cd();

      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      //L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral());
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Draw("p e");

      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2]);
      //L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral());
      //L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Draw("same p e");

       TF1 *fitL0_L0bar_pT_US_ThetaStar = new TF1("fitL0_L0bar_pT_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_pT_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitL0_L0bar_pT_US_ThetaStar, "s i 0 r");

      float P_L0_L0bar_pT = fitL0_L0bar_pT_US_ThetaStar->GetParameter(1)/(L0_alpha*L0bar_alpha);
      float P_L0_L0bar_pT_err = fitL0_L0bar_pT_US_ThetaStar->GetParError(1)/(L0_alpha*L0bar_alpha);

      fitL0_L0bar_pT_US_ThetaStar->SetLineColor(1);
      fitL0_L0bar_pT_US_ThetaStar->Draw("same");

      TPaveText *L0_L0bar_pT_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      L0_L0bar_pT_text->SetTextFont(42);
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0bar_pT_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_pT_text->AddText("Minimum bias");
      L0_L0bar_pT_text->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_pT_text->AddText("|#it{y}| < 1");
      L0_L0bar_pT_text->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0bar_pT_text->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0bar_pT_text->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar_pT, P_L0_L0bar_pT_err));
      L0_L0bar_pT_text->SetFillColorAlpha(0, 0.01);
      L0_L0bar_pT_text->Draw("same");

      L0_L0bar_leg->Draw("same");

      L0_L0bar_cosThetaProdPlane_pT_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //____________________________________________________________________________________________________


      TCanvas *L0_L0_cosThetaProdPlane_pT_no_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0_L0_cosThetaProdPlane_pT_no_corr_can->cd();

      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      nLL_pT[pTbin1][pTbin2] = L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral();
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Sumw2();
      //L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      //L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral());
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Draw("p e");

      L0_L0_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0_L0_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_L0_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Sumw2();
      L0_L0_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Scale(L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral()/L0_L0_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Integral());
      L0_L0_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0_L0_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Draw("p e same");

      L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      nLL_pT_back[pTbin1][pTbin2] = L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral();
      L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Sumw2();
      //L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      //L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Draw("p e same");

      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetLineColor(kBlue);
      //double nLLbar = L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral();
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Sumw2();
      //L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2]);
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral()/L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral());
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMinimum(0);
      //L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Draw("same p e");

      TF1 *fitL0_L0_pT_US_ThetaStar_no_corr = new TF1("fitL0_L0_pT_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_pT_US_ThetaStar_no_corr->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitL0_L0_pT_US_ThetaStar_no_corr, "s i 0 r");

      float P_L0_L0_pT_no_corr = fitL0_L0_pT_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_L0_L0_pT_no_corr_err = fitL0_L0_pT_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0_alpha);

      fitL0_L0_pT_US_ThetaStar_no_corr->SetLineColor(1);
      fitL0_L0_pT_US_ThetaStar_no_corr->Draw("same");

      TPaveText *L0_L0_pT_text_no_corr = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      L0_L0_pT_text_no_corr->SetTextFont(42);
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_pT_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_pT_text_no_corr->AddText("Minimum bias, no correction");
      L0_L0_pT_text_no_corr->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_pT_text_no_corr->AddText("|#it{y}| < 1");
      L0_L0_pT_text_no_corr->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0_pT_text_no_corr->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0_pT_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0_pT_no_corr, P_L0_L0_pT_no_corr_err));
      L0_L0_pT_text_no_corr->SetFillColorAlpha(0, 0.01);
      L0_L0_pT_text_no_corr->Draw("same");

      L0_L0_leg->Draw("same");

      L0_L0_cosThetaProdPlane_pT_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0_L0_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //_______________________________________________________________________

      TCanvas *L0_L0_cosThetaProdPlane_pT_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0_L0_cosThetaProdPlane_pT_can->cd();


      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      //L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral());
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Draw("p e");

      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2]);
      //L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral());
      //L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Draw("same p e");


      TF1 *fitL0_L0_pT_US_ThetaStar = new TF1("fitL0_L0_pT_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_pT_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitL0_L0_pT_US_ThetaStar, "s i 0 r");

      float P_L0_L0_pT = fitL0_L0_pT_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_L0_L0_pT_err = fitL0_L0_pT_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

      fitL0_L0_pT_US_ThetaStar->SetLineColor(1);
      fitL0_L0_pT_US_ThetaStar->Draw("same");


      TPaveText *L0_L0_pT_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      L0_L0_pT_text->SetTextFont(42);
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_pT_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_pT_text->AddText("Minimum bias");
      L0_L0_pT_text->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_pT_text->AddText("|#it{y}| < 1");
      L0_L0_pT_text->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0_pT_text->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0_pT_text->AddText(Form("P = %.2f #pm %.2f", P_L0_L0_pT, P_L0_L0_pT_err));
      L0_L0_pT_text->SetFillColorAlpha(0, 0.01);
      L0_L0_pT_text->Draw("same");

      L0_L0_leg->Draw("same");


      L0_L0_cosThetaProdPlane_pT_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //____________________________________________________________________________________________________


      TCanvas *L0bar_L0bar_cosThetaProdPlane_pT_no_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0bar_L0bar_cosThetaProdPlane_pT_no_corr_can->cd();

      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      nLbarLbar_pT[pTbin1][pTbin2] = L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral();
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      //L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Draw("p e");

      L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Scale(L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral()/L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0bar_L0bar_cosThetaProdPlane_pT_US_side_band_hist[pTbin1][pTbin2]->Draw("p e same");

      L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      nLbarLbar_pT_back[pTbin1][pTbin2] = L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral();
      L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Sumw2();
      //L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      //L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Draw("p e same");

      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetLineColor(kBlue);
      //double nLLbar = L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral();
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Sumw2();
      //L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2]);
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral()/L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMinimum(0);
      //L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Draw("same p e");



      TF1 *fitL0bar_L0bar_pT_US_ThetaStar_no_corr = new TF1("fitL0bar_L0bar_pT_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_pT_US_ThetaStar_no_corr->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitL0bar_L0bar_pT_US_ThetaStar_no_corr, "s i 0 r");

      float P_L0bar_L0bar_pT_no_corr = fitL0bar_L0bar_pT_US_ThetaStar_no_corr->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
      float P_L0bar_L0bar_pT_no_corr_err = fitL0bar_L0bar_pT_US_ThetaStar_no_corr->GetParError(1)/(L0bar_alpha*L0bar_alpha);

      fitL0bar_L0bar_pT_US_ThetaStar_no_corr->SetLineColor(1);
      fitL0bar_L0bar_pT_US_ThetaStar_no_corr->Draw("same");

      TPaveText *L0bar_L0bar_pT_text_no_corr = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      L0bar_L0bar_pT_text_no_corr->SetTextFont(42);
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_pT_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_pT_text_no_corr->AddText("Minimum bias, no correction");
      L0bar_L0bar_pT_text_no_corr->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_pT_text_no_corr->AddText("|#it{y}| < 1");
      L0bar_L0bar_pT_text_no_corr->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0bar_L0bar_pT_text_no_corr->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0bar_L0bar_pT_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar_pT_no_corr, P_L0bar_L0bar_pT_no_corr_err));
      L0bar_L0bar_pT_text_no_corr->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_pT_text_no_corr->Draw("same");

      L0bar_L0bar_leg->Draw("same");

      L0bar_L0bar_cosThetaProdPlane_pT_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0bar_L0bar_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //_____________________________________________________________________________

      TCanvas *L0bar_L0bar_cosThetaProdPlane_pT_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0bar_L0bar_cosThetaProdPlane_pT_can->cd();

      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      //L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Draw("p e");

      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2]);
      //L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral());
      //L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Draw("same p e");

      TF1 *fitL0bar_L0bar_pT_US_ThetaStar = new TF1("fitL0bar_L0bar_pT_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_pT_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitL0bar_L0bar_pT_US_ThetaStar, "s i 0 r");

      float P_L0bar_L0bar_pT = fitL0bar_L0bar_pT_US_ThetaStar->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
      float P_L0bar_L0bar_pT_err = fitL0bar_L0bar_pT_US_ThetaStar->GetParError(1)/(L0bar_alpha*L0bar_alpha);

      fitL0bar_L0bar_pT_US_ThetaStar->SetLineColor(1);
      fitL0bar_L0bar_pT_US_ThetaStar->Draw("same");

      TPaveText *L0bar_L0bar_pT_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      L0bar_L0bar_pT_text->SetTextFont(42);
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_pT_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_pT_text->AddText("Minimum bias");
      L0bar_L0bar_pT_text->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_pT_text->AddText("|#it{y}| < 1");
      L0bar_L0bar_pT_text->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0bar_L0bar_pT_text->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0bar_L0bar_pT_text->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar_pT, P_L0bar_L0bar_pT_err));
      L0bar_L0bar_pT_text->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_pT_text->Draw("same");

      L0bar_L0bar_leg->Draw("same");

      L0bar_L0bar_cosThetaProdPlane_pT_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

    }
  }

  //___________________________________________________________________________________________________________________________________________________________________

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      TCanvas *L0_L0bar_cosThetaProdPlane_eta_no_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i_can", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i_can", etaBin1, etaBin2), 1200, 1000);

      L0_L0bar_cosThetaProdPlane_eta_no_corr_can->cd();

      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMarkerSize(2);
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      //double nLLbar = L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Integral();
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Sumw2();
      //L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Divide(L0_L0bar_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Scale(1./L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      //L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Scale(1./L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Integral());
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Draw("p e");

      L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->SetMarkerSize(2);
      L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Scale(L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral()/L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Integral());
      L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Scale(1./L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      L0_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Draw("p e same");

      L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerSize(2);
      //double nLLbar_back = L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral();
      L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Sumw2();
      //L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Scale(1./L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      //L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Divide(L0_L0bar_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      //L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Draw("p e same");

      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMarkerSize(2);
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetLineColor(kBlue);
      //double nLLbar = L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Integral();
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Sumw2();
      //L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Divide(L0_L0bar_cosThetaProdPlane_ME_eta_eff[etaBin1][etaBin2]);
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Scale(L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral()/L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Integral());
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Scale(1./L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMinimum(0);
      //L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Draw("same p e");


      TF1 *fitL0_L0bar_eta_US_ThetaStar_no_corr = new TF1("fitL0_L0bar_eta_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_eta_US_ThetaStar_no_corr->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_eta_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Fit(fitL0_L0bar_eta_US_ThetaStar_no_corr, "s i 0 r");

      float P_L0_L0bar_eta_no_corr = fitL0_L0bar_eta_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
      float P_L0_L0bar_eta_no_corr_err = fitL0_L0bar_eta_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

      fitL0_L0bar_eta_US_ThetaStar_no_corr->SetLineColor(1);
      fitL0_L0bar_eta_US_ThetaStar_no_corr->Draw("same");

      TPaveText *L0_L0bar_eta_text_no_corr = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0bar_eta_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_eta_text_no_corr->AddText("Minimum bias, no correction");
      L0_L0bar_eta_text_no_corr->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_eta_text_no_corr->AddText("|#it{y}| < 1");
      L0_L0bar_eta_text_no_corr->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0bar_eta_text_no_corr->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0bar_eta_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar_eta_no_corr, P_L0_L0bar_eta_no_corr_err));
      L0_L0bar_eta_text_no_corr->SetFillColorAlpha(0, 0.01);
      L0_L0bar_eta_text_no_corr->Draw("same");

      L0_L0bar_leg->Draw("same");

      L0_L0bar_cosThetaProdPlane_eta_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0_L0bar_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i.png", etaBin1, etaBin2));
      //_____________________________________________________________________

      TCanvas *L0_L0bar_cosThetaProdPlane_eta_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_can", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_can", etaBin1, etaBin2), 1200, 1000);

      L0_L0bar_cosThetaProdPlane_eta_can->cd();

      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Divide(L0_L0bar_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      //L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Scale(1./L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Integral());
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Draw("p e");

      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Divide(L0_L0bar_cosThetaProdPlane_ME_eta_eff[etaBin1][etaBin2]);
      //L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Scale(1./L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Integral());
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Draw("same p e");

      TF1 *fitL0_L0bar_eta_US_ThetaStar = new TF1("fitL0_L0bar_eta_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_eta_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_eta_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Fit(fitL0_L0bar_eta_US_ThetaStar, "s i 0 r");

      float P_L0_L0bar_eta = fitL0_L0bar_eta_US_ThetaStar->GetParameter(1)/(L0_alpha*L0bar_alpha);
      float P_L0_L0bar_eta_err = fitL0_L0bar_eta_US_ThetaStar->GetParError(1)/(L0_alpha*L0bar_alpha);

      fitL0_L0bar_eta_US_ThetaStar->SetLineColor(1);
      fitL0_L0bar_eta_US_ThetaStar->Draw("same");

      TPaveText *L0_L0bar_eta_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0bar_eta_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_eta_text->AddText("Minimum bias");
      L0_L0bar_eta_text->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_eta_text->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0bar_eta_text->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0bar_eta_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
      L0_L0bar_eta_text->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar_eta, P_L0_L0bar_eta_err));
      L0_L0bar_eta_text->SetFillColorAlpha(0, 0.01);
      L0_L0bar_eta_text->Draw("same");

      L0_L0bar_leg->Draw("same");

      L0_L0bar_cosThetaProdPlane_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i.png", etaBin1, etaBin2));
      //____________________________________________________________________________________________________



      TCanvas *L0_L0_cosThetaProdPlane_eta_no_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i_can", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i_can", etaBin1, etaBin2), 1200, 1000);

      L0_L0_cosThetaProdPlane_eta_no_corr_can->cd();

      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMarkerSize(2);
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      //double nLL = L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Integral();
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Sumw2();
      //L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Divide(L0_L0_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Scale(1./L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      //L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Scale(1./L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Integral());
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Draw("p e");

      L0_L0_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      L0_L0_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->SetMarkerSize(2);
      L0_L0_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Sumw2();
      L0_L0_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Scale(L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral()/L0_L0_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Integral());
      L0_L0_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Scale(1./L0_L0_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      L0_L0_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Draw("p e same");

      L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerSize(2);
      //double nLL_back = L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral();
      L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Sumw2();
      //L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Scale(1./L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      //L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Divide(L0_L0_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      //L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Draw("p e same");

      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMarkerSize(2);
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetLineColor(kBlue);
      //double nLLbar = L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Integral();
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Sumw2();
      //L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Divide(L0_L0_cosThetaProdPlane_ME_eta_eff[etaBin1][etaBin2]);
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Scale(L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral()/L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Integral());
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Scale(1./L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMinimum(0);
      //L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Draw("same p e");


      TF1 *fitL0_L0_eta_US_ThetaStar_no_corr = new TF1("fitL0_L0_eta_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_eta_US_ThetaStar_no_corr->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_eta_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Fit(fitL0_L0_eta_US_ThetaStar_no_corr, "s i 0 r");

      float P_L0_L0_eta_no_corr = fitL0_L0_eta_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_L0_L0_eta_no_corr_err = fitL0_L0_eta_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0_alpha);

      fitL0_L0_eta_US_ThetaStar_no_corr->SetLineColor(1);
      fitL0_L0_eta_US_ThetaStar_no_corr->Draw("same");

      TPaveText *L0_L0_eta_text_no_corr = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_eta_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_eta_text_no_corr->AddText("Minimum bias, no correction");
      L0_L0_eta_text_no_corr->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_eta_text_no_corr->AddText("|#it{y}| < 1");
      L0_L0_eta_text_no_corr->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0_eta_text_no_corr->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0_eta_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0_eta_no_corr, P_L0_L0_eta_no_corr_err));
      L0_L0_eta_text_no_corr->SetFillColorAlpha(0, 0.01);
      L0_L0_eta_text_no_corr->Draw("same");

      L0_L0_leg->Draw("same");


      L0_L0_cosThetaProdPlane_eta_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0_L0_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i.png", etaBin1, etaBin2));
      //____________________________________________________________________

      TCanvas *L0_L0_cosThetaProdPlane_eta_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_can", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_can", etaBin1, etaBin2), 1200, 1000);

      L0_L0_cosThetaProdPlane_eta_can->cd();

      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Divide(L0_L0_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      //L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Scale(1./L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Integral());
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Draw("p e");

      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Divide(L0_L0_cosThetaProdPlane_ME_eta_eff[etaBin1][etaBin2]);
      //L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Scale(1./L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Integral());
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Draw("same p e");

      TF1 *fitL0_L0_eta_US_ThetaStar = new TF1("fitL0_L0_eta_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_eta_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_eta_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Fit(fitL0_L0_eta_US_ThetaStar, "s i 0 r");

      float P_L0_L0_eta = fitL0_L0_eta_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_L0_L0_eta_err = fitL0_L0_eta_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

      fitL0_L0_eta_US_ThetaStar->SetLineColor(1);
      fitL0_L0_eta_US_ThetaStar->Draw("same");

      TPaveText *L0_L0_eta_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_eta_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_eta_text->AddText("Minimum bias");
      L0_L0_eta_text->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_eta_text->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0_eta_text->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0_eta_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
      L0_L0_eta_text->AddText(Form("P = %.2f #pm %.2f", P_L0_L0_eta, P_L0_L0_eta_err));
      L0_L0_eta_text->SetFillColorAlpha(0, 0.01);
      L0_L0_eta_text->Draw("same");

      L0_L0_leg->Draw("same");

      L0_L0_cosThetaProdPlane_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i.png", etaBin1, etaBin2));
      //____________________________________________________________________________________________________



      TCanvas *L0bar_L0bar_cosThetaProdPlane__eta_no_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i_can", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i_can", etaBin1, etaBin2), 1200, 1000);

      L0bar_L0bar_cosThetaProdPlane__eta_no_corr_can->cd();

      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      //double nLbarLbar = L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Integral();
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Divide(L0bar_L0bar_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      //L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Draw("p e");

      L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Sumw2();
      L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Scale(L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral()/L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      L0bar_L0bar_cosThetaProdPlane_eta_US_side_band_hist[etaBin1][etaBin2]->Draw("p e same");

      L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      //double nLbarLbar_back = L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral();
      L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Sumw2();
      //L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      //L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Divide(L0bar_L0bar_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      //L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Draw("p e same");

      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetLineColor(kBlue);
      //double nLLbar = L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Integral();
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Sumw2();
      //L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Divide(L0bar_L0bar_cosThetaProdPlane_ME_eta_eff[etaBin1][etaBin2]);
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Scale(L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral()/L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMinimum(0);
      //L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Draw("same p e");

      TF1 *fitL0bar_L0bar_eta_US_ThetaStar_no_corr = new TF1("fitL0bar_L0bar_eta_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_eta_US_ThetaStar_no_corr->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_eta_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Fit(fitL0bar_L0bar_eta_US_ThetaStar_no_corr, "s i 0 r");

      float P_L0bar_L0bar_eta_no_corr = fitL0bar_L0bar_eta_US_ThetaStar_no_corr->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
      float P_L0bar_L0bar_eta_no_corr_err = fitL0bar_L0bar_eta_US_ThetaStar_no_corr->GetParError(1)/(L0bar_alpha*L0bar_alpha);

      fitL0bar_L0bar_eta_US_ThetaStar_no_corr->SetLineColor(1);
      fitL0bar_L0bar_eta_US_ThetaStar_no_corr->Draw("same");

      TPaveText *L0bar_L0bar_eta_text_no_corr = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_eta_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_eta_text_no_corr->AddText("Minimum bias, no correction");
      L0bar_L0bar_eta_text_no_corr->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_eta_text_no_corr->AddText("|#it{y}| < 1");
      L0bar_L0bar_eta_text_no_corr->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0bar_L0bar_eta_text_no_corr->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0bar_L0bar_eta_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar_eta_no_corr, P_L0bar_L0bar_eta_no_corr_err));
      L0bar_L0bar_eta_text_no_corr->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_eta_text_no_corr->Draw("same");

      L0bar_L0bar_leg->Draw("same");


      L0bar_L0bar_cosThetaProdPlane__eta_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0bar_L0bar_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i.png", etaBin1, etaBin2));
      //_________________________________________________________

      TCanvas *L0bar_L0bar_cosThetaProdPlane_eta_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_can", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_can", etaBin1, etaBin2), 1200, 1000);

      L0bar_L0bar_cosThetaProdPlane_eta_can->cd();

      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Divide(L0bar_L0bar_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      //L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Draw("p e");

      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Divide(L0bar_L0bar_cosThetaProdPlane_ME_eta_eff[etaBin1][etaBin2]);
      //L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Draw("p e same");

      TF1 *fitL0bar_L0bar_eta_US_ThetaStar = new TF1("fitL0bar_L0bar_eta_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_eta_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_eta_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Fit(fitL0bar_L0bar_eta_US_ThetaStar, "s i 0 r");

      float P_L0bar_L0bar_eta = fitL0bar_L0bar_eta_US_ThetaStar->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
      float P_L0bar_L0bar_eta_err = fitL0bar_L0bar_eta_US_ThetaStar->GetParError(1)/(L0bar_alpha*L0bar_alpha);

      fitL0bar_L0bar_eta_US_ThetaStar->SetLineColor(1);
      fitL0bar_L0bar_eta_US_ThetaStar->Draw("same");


      TPaveText *L0bar_L0bar_eta_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_eta_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_eta_text->AddText("Minimum bias");
      L0bar_L0bar_eta_text->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_eta_text->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0bar_L0bar_eta_text->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0bar_L0bar_eta_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
      L0bar_L0bar_eta_text->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar_eta, P_L0bar_L0bar_eta_err));
      L0bar_L0bar_eta_text->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_eta_text->Draw("same");

      L0bar_L0bar_leg->Draw("same");


      L0bar_L0bar_cosThetaProdPlane_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/L_correlations/L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i.png", etaBin1, etaBin2));
      //____________________________________________________________________________________________________
    }
  }


  //----------------------------------L pair stats---------------------------------

  cout<<endl;
  cout<<"N L-Lbar pairs from hist: "<<nLLbar<<endl;
  cout<<"N L-Lbar background pairs from hist: "<<nLLbar_back<<endl;
  cout<<endl;
  cout<<"N L-L pairs from hist: "<<nLL<<endl;
  cout<<"N L-L background pairs from hist: "<<nLL_back<<endl;
  cout<<endl;
  cout<<"N Lbar-Lbar pairs from hist: "<<nLbarLbar<<endl;
  cout<<"N Lbar-Lbar background pairs from hist: "<<nLbarLbar_back<<endl;

  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      cout<<endl;
      cout<<"pT1: "<<pTbin1<<", pT2: "<<pTbin2<<endl;
      cout<<"N L-Lbar pairs from hist: "<<nLLbar_pT[pTbin1][pTbin2]<<endl;
      cout<<"N L-Lbar background pairs from hist: "<<nLLbar_pT_back[pTbin1][pTbin2]<<endl;
      cout<<endl;
      cout<<"N L-L pairs from hist: "<<nLL_pT[pTbin1][pTbin2]<<endl;
      cout<<"N L-L background pairs from hist: "<<nLL_pT_back[pTbin1][pTbin2]<<endl;
      cout<<endl;
      cout<<"N Lbar-Lbar pairs from hist: "<<nLbarLbar_pT[pTbin1][pTbin2]<<endl;
      cout<<"N Lbar-Lbar background pairs from hist: "<<nLbarLbar_pT_back[pTbin1][pTbin2]<<endl;
    }
  }




  LLbarOutFile->Close();


  return true;
}
//__________________________________________________________________________________________________________________________

//ReadMode = 0 - read TTree, ReadMode = 1 - read histograms - First run in ReadMode = 0 to save relevant histograms, then can run in ReadMode = 1 to read just histograms and save time
//energy - collision energy in GeV
void Ana003_Lambda_corr_2D(const int ReadMode = 0, const int energy = 510, const int year = 2017)
{
  ifstream fileList;

  if(energy == 510)
  {
    //MB without strict TOF matching
    fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run17_MB_noTOF/fileList.list");

  }
  else if(energy == 200)
  {
    //MB without strict TOF matching
    fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF/fileList.list");

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


  TCanvas *testCan = new TCanvas("testCan", "testCan", 1200, 1000);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  testCan->cd();
  hEventStat1->GetXaxis()->SetRange(1,6);
  hEventStat1->GetYaxis()->SetMaxDigits(3);
  hEventStat1->SetMinimum(0);
  hEventStat1->Draw();

  testCan->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/EventStat_test.png");


  bool AnaFinish = DoAnalysis(myChain, ReadMode, energy , year);

  if(!AnaFinish)
  {
    cout<<"Analysis ended abnormally. Aborting!"<<endl;

    return;
  }


  cout<<endl;
  cout<<"Nubmer of accepted events: "<<hEventStat1->GetBinContent(6)<<endl;


  return;
}
