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

bool cuts(float K0s_y)
{

  if( !( TMath::Abs(K0s_y) < K0s_y_cut ) ) return false;
  //if( strictTOF_cut == 1 && pi1_hasTOFinfo == 0 ) return false; //TOF matched pions
  //if( strictTOF_cut == 2 && (pi1_hasTOFinfo == 0 || pi2_hasTOFinfo == 0) ) return false; //TOF matched pions and protons
  //if(cos(K0s_theta) < K0s_cos_theta_cut) return false;
  //if(K0s_decayL > K0s_decayK0s_cut) return false;

  return true;

}

//analyze invariant mass spectra
//arguments are vectors of bin numbers of the invariant mass peak
//the bins are determined via fit to the invariant mass spectra
//bool InvMass(TTree *K0s_tree, vector<int> &invMassBins_L, vector<int> &invMassBins_L0, vector<int> &invMassBins_L0bar, const int readMode)
bool DoAnalysis(TChain *K0s_tree, const int ReadMode, const int energy = 510, const int year = 2017)
{

  TFile *InvMassFile;

  if(ReadMode == 0) //create invariant mass file Form nTuple - run in this mode first
  {
    InvMassFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/InvMass_K0s_2D_work.root", "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT
  }
  else  //read file to store wrong and good sign M_inv distributions and daughter pT
  {

    InvMassFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/InvMass_K0s_2D_work.root", "read"); //old, non-optimized cuts

    if( !(InvMassFile->IsOpen()) )
    {
      cout<<"Unable to open file with invariant mass histograms!"<<endl;
      return false;
    }
  }



  //_______________________________________________________________________________________________________________________________________________


  Int_t charge;
  Float_t K0s_mass, K0s_pt, K0s_eta, K0s_phi;//, K0s_decayL, K0s_theta, K0s_DCAdaughters;

  Float_t pi1_pt, pi2_pt;
  Float_t pi1_eta, pi2_eta;
  Float_t pi1_phi, pi2_phi;
  //Float_t pi1_ch;
  Int_t pi1_hasTOFinfo, pi2_hasTOFinfo;

  Int_t eventId;

  //---------------SET BARANCH ADDRESSES------------------------
  K0s_tree->SetBranchAddress("pair_charge", &charge);
  K0s_tree->SetBranchAddress("pair_mass", &K0s_mass);
  K0s_tree->SetBranchAddress("pair_pt", &K0s_pt);
  K0s_tree->SetBranchAddress("pair_eta", &K0s_eta);
  K0s_tree->SetBranchAddress("pair_phi", &K0s_phi);
  //K0s_tree->SetBranchAddress("pair_decayL", &K0s_decayL);
  //K0s_tree->SetBranchAddress("pair_theta", &K0s_theta);
  //K0s_tree->SetBranchAddress("pair_DCAdaughters", &K0s_DCAdaughters);

  K0s_tree->SetBranchAddress("p1_pt", &pi1_pt);
  K0s_tree->SetBranchAddress("p1_eta", &pi1_eta);
  K0s_tree->SetBranchAddress("p1_phi", &pi1_phi);
  //K0s_tree->SetBranchAddress("p1_ch", &pi1_ch);
  K0s_tree->SetBranchAddress("p1_hasTOFinfo", &pi1_hasTOFinfo);

  K0s_tree->SetBranchAddress("p2_pt", &pi2_pt);
  K0s_tree->SetBranchAddress("p2_eta", &pi2_eta);
  K0s_tree->SetBranchAddress("p2_phi", &pi2_phi);
  K0s_tree->SetBranchAddress("p2_hasTOFinfo", &pi2_hasTOFinfo);

  K0s_tree->SetBranchAddress("eventId", &eventId);

  //--------------------------------------------------------------------------


  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass_US[nPtBins_corr][nPtBins_corr];
  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //bacground when US is paired with LS
  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass_LS[nPtBins_corr][nPtBins_corr];
  TH2F *K0s1_inv_mass_vs_K0s2_inv_mass[nPtBins_corr][nPtBins_corr];


  if(ReadMode == 0) //create histograms to be saved into file
  {
    for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
    {
      for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
      {
        K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2] = new TH2F(Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2),200, 0.45, 0.55, 200, 0.45, 0.55);
        K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2] = new TH2F(Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2),200, 0.45, 0.55, 200, 0.45, 0.55);
        K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2] = new TH2F(Form("K0s1_inv_mass_vs_K0s2_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("K0s1_inv_mass_vs_K0s2_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2),200, 0.45, 0.55, 200, 0.45, 0.55);
      }
    }
  }
  else //load histograms Form file
  {
    for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
    {
      for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
      {
        K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2));
        K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("K0s1_inv_mass_vs_K0s2_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2));
        K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("K0s1_inv_mass_vs_K0s2_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2));
      }
    }
  }
  //________________________________________________________________________________________


  //to store Lorentz vectors for L pair analysis
  vector<TLorentzVector> K0s_vector;
  vector<int> K0s_pT_bin_vector;
  vector<int> K0s_eta_bin_vector;
  vector<TLorentzVector> pi1_vector;
  vector<TLorentzVector> pi2_vector;

  vector<TLorentzVector> K0s_vector_background;
  vector<int> K0s_pT_bin_vector_background;
  vector<int> K0s_eta_bin_vector_background;
  vector<TLorentzVector> pi1_vector_background;
  vector<TLorentzVector> pi2_vector_background;

  //vectors for cos(theta*)
  vector<float> K0s_K0s_cos_theta;
  vector<float> K0s_K0s_Minv_K0s1;
  vector<float> K0s_K0s_Minv_K0s2;
  vector<int> K0s_K0s_pT_bin_K0s1;
  vector<int> K0s_K0s_pT_bin_K0s2;
  vector<int> K0s_K0s_eta_bin_K0s1;
  vector<int> K0s_K0s_eta_bin_K0s2;


  vector<float> K0s_K0s_cos_theta_back;
  vector<float> K0s_K0s_Minv_K0s1_back;
  vector<float> K0s_K0s_Minv_K0s2_back;
  vector<int> K0s_K0s_pT_bin_K0s1_back;
  vector<int> K0s_K0s_pT_bin_K0s2_back;
  vector<int> K0s_K0s_eta_bin_K0s1_back;
  vector<int> K0s_K0s_eta_bin_K0s2_back;

  //vectors for mixed-event
  vector<TLorentzVector> K0s_vector_ME;
  vector<int> K0s_pT_bin_vector_ME;
  vector<int> K0s_eta_bin_vector_ME;
  vector<TLorentzVector> pi_vector_ME;

  //LS pairs for mixed event (for background correction)
  vector<TLorentzVector> K0s_vector_ME_LS;
  vector<int> K0s_pT_bin_vector_ME_LS;
  vector<int> K0s_eta_bin_vector_ME_LS;
  vector<TLorentzVector> pi_vector_ME_LS;


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


    //calculate K0s rapidity y
    double K0s_y = rapidity(K0s_pt, K0s_eta, K0s_mass_PDG);

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
    K0s_Lorentz_vector.SetPtEtaPhiM(K0s_pt, K0s_eta, K0s_phi, K0s_mass);

    TLorentzVector pi1_Lorenz_vector(1,1,1,1);
    pi1_Lorenz_vector.SetPtEtaPhiM(pi1_pt, pi1_eta, pi1_phi, pi_mass_PDG);

    TLorentzVector pi2_Lorenz_vector(1,1,1,1);
    pi2_Lorenz_vector.SetPtEtaPhiM(pi2_pt, pi2_eta, pi2_phi, pi_mass_PDG);


    if(eventId == eventID_last) //same event as in previous iteration and first event
    {
      if(charge == 0 ) //like-sign combinations
      {
        if( cuts(K0s_y) && pT_bin_corr != -1 && eta_bin != -1)
        {
          K0s_vector.push_back(K0s_Lorentz_vector);
          K0s_pT_bin_vector.push_back(pT_bin_corr);
          K0s_eta_bin_vector.push_back(eta_bin);

          pi1_vector.push_back(pi1_Lorenz_vector);
          pi2_vector.push_back(pi2_Lorenz_vector);
        }
      }
      else if( charge == 1 )//unlike-sign
      {
        if( cuts(K0s_y) && pT_bin_corr != -1 && eta_bin != -1)
        {
          K0s_vector_background.push_back(K0s_Lorentz_vector);
          K0s_pT_bin_vector_background.push_back(pT_bin_corr);
          K0s_eta_bin_vector_background.push_back(eta_bin);

          pi1_vector_background.push_back(pi1_Lorenz_vector);
          pi2_vector_background.push_back(pi2_Lorenz_vector);
        }
      }

    }
    else //new event
    {
      eventID_last = eventId;

      //at least one L0-L0 pair in event
      //US-US charge combination
      if(K0s_vector.size() > 1)
      {
        for(unsigned int iK0s1 = 0; iK0s1 < K0s_vector.size(); iK0s1++)
        {
          for(unsigned int iK0s2 = iK0s1+1; iK0s2 < K0s_vector.size(); iK0s2++)
          {
            //check auto-correlation
            if(pi1_vector.at(iK0s1).Phi() == pi1_vector.at(iK0s2).Phi()) continue;
            if(pi2_vector.at(iK0s1).Phi() == pi2_vector.at(iK0s2).Phi()) continue;
            if(pi1_vector.at(iK0s1).Phi() == pi2_vector.at(iK0s2).Phi()) continue;
            if(pi2_vector.at(iK0s1).Phi() == pi1_vector.at(iK0s2).Phi()) continue;

            K0s1_inv_mass_vs_K0s2_inv_mass_US[K0s_pT_bin_vector.at(iK0s1)][K0s_pT_bin_vector.at(iK0s2)]->Fill(K0s_vector.at(iK0s1).M(), K0s_vector.at(iK0s2).M());

            //use pion 1 as reference particle
            double K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi1_vector.at(iK0s1), K0s_vector.at(iK0s2), pi1_vector.at(iK0s2));

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
        {
          for(unsigned int iK0s2 = 0; iK0s2 < K0s_vector_background.size(); iK0s2++)
          {
            //check auto-correlation
            if(pi1_vector.at(iK0s1).Phi() == pi1_vector_background.at(iK0s2).Phi()) continue;
            if(pi2_vector.at(iK0s1).Phi() == pi2_vector_background.at(iK0s2).Phi()) continue;
            if(pi1_vector.at(iK0s1).Phi() == pi2_vector_background.at(iK0s2).Phi()) continue;
            if(pi2_vector.at(iK0s1).Phi() == pi1_vector_background.at(iK0s2).Phi()) continue;

            if(nFill_bckg % 2 == 0)//fill US-LS
            {
              K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[K0s_pT_bin_vector.at(iK0s1)][K0s_pT_bin_vector_background.at(iK0s2)]->Fill(K0s_vector.at(iK0s1).M(), K0s_vector_background.at(iK0s2).M());

              //use pion 1 as reference particle
              double K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi1_vector.at(iK0s1), K0s_vector_background.at(iK0s2), pi1_vector_background.at(iK0s2));

              K0s_K0s_cos_theta_back.push_back(TMath::Cos(K0s_K0s_pairThetaStar));

              K0s_K0s_Minv_K0s1_back.push_back(K0s_vector.at(iK0s1).M());
              K0s_K0s_Minv_K0s2_back.push_back(K0s_vector_background.at(iK0s2).M());

              K0s_K0s_pT_bin_K0s1_back.push_back(K0s_pT_bin_vector.at(iK0s1));
              K0s_K0s_pT_bin_K0s2_back.push_back(K0s_pT_bin_vector_background.at(iK0s2));

              K0s_K0s_eta_bin_K0s1_back.push_back(K0s_eta_bin_vector.at(iK0s1));
              K0s_K0s_eta_bin_K0s2_back.push_back(K0s_eta_bin_vector_background.at(iK0s2));

            }
            else //fill LS-US
            {
              K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[K0s_pT_bin_vector_background.at(iK0s2)][K0s_pT_bin_vector.at(iK0s1)]->Fill(K0s_vector_background.at(iK0s2).M(), K0s_vector.at(iK0s1).M());

              //use pion 1 as reference particle
              double K0s_K0s_pairThetaStar = LpairThetaStar( K0s_vector_background.at(iK0s2), pi1_vector_background.at(iK0s2), K0s_vector.at(iK0s1), pi1_vector.at(iK0s1));

              K0s_K0s_cos_theta_back.push_back(TMath::Cos(K0s_K0s_pairThetaStar));

              K0s_K0s_Minv_K0s2_back.push_back(K0s_vector.at(iK0s1).M());
              K0s_K0s_Minv_K0s1_back.push_back(K0s_vector_background.at(iK0s2).M());

              K0s_K0s_pT_bin_K0s2_back.push_back(K0s_pT_bin_vector.at(iK0s1));
              K0s_K0s_pT_bin_K0s1_back.push_back(K0s_pT_bin_vector_background.at(iK0s2));

              K0s_K0s_eta_bin_K0s2_back.push_back(K0s_eta_bin_vector.at(iK0s1));
              K0s_K0s_eta_bin_K0s1_back.push_back(K0s_eta_bin_vector_background.at(iK0s2));
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
        {
          for(unsigned int iK0s2 = iK0s1+1; iK0s2 < K0s_vector_background.size(); iK0s2++)
          {
            if(pi1_vector_background.at(iK0s1).Phi() == pi1_vector_background.at(iK0s2).Phi()) continue;
            if(pi2_vector_background.at(iK0s1).Phi() == pi2_vector_background.at(iK0s2).Phi()) continue;
            if(pi1_vector_background.at(iK0s1).Phi() == pi2_vector_background.at(iK0s2).Phi()) continue;
            if(pi2_vector_background.at(iK0s1).Phi() == pi1_vector_background.at(iK0s2).Phi()) continue;

            K0s1_inv_mass_vs_K0s2_inv_mass_LS[K0s_pT_bin_vector_background.at(iK0s1)][K0s_pT_bin_vector_background.at(iK0s2)]->Fill(K0s_vector_background.at(iK0s1).M(), K0s_vector_background.at(iK0s2).M());

          }
        }
      }
      //_____________________________________________________________________________________________


      //fill vectors for mixed event
      //selecting events with only one L or L-bar
      //need to select more - Minv cut will be applied later

      if( K0s_vector.size() == 1 && K0s_vector_ME.size() < 1e4)
      {
        K0s_vector_ME.push_back(K0s_vector.at(0));
        pi_vector_ME.push_back(pi1_vector.at(0));

        K0s_pT_bin_vector_ME.push_back(K0s_pT_bin_vector.at(0));
        K0s_eta_bin_vector_ME.push_back(K0s_eta_bin_vector.at(0));

      }

      if( K0s_vector_background.size() == 1 && K0s_vector_ME_LS.size() < 1e4)
      {
        K0s_vector_ME_LS.push_back(K0s_vector_background.at(0));
        pi_vector_ME_LS.push_back(pi1_vector_background.at(0));

        K0s_pT_bin_vector_ME_LS.push_back(K0s_pT_bin_vector_background.at(0));
        K0s_eta_bin_vector_ME_LS.push_back(K0s_eta_bin_vector_background.at(0));
      }
      //_______________________________________________________________________________________________

      //reset number of K0ss and vectors
      K0s_vector.clear();
      K0s_pT_bin_vector.clear();
      K0s_eta_bin_vector.clear();
      pi1_vector.clear();
      pi2_vector.clear();

      //check if we have good K0s in the new event
      if( charge == 0 && cuts(K0s_y) &&  pT_bin_corr != -1 && eta_bin != -1)
      {
        K0s_vector.push_back(K0s_Lorentz_vector);
        K0s_pT_bin_vector.push_back(pT_bin_corr);
        K0s_eta_bin_vector.push_back(eta_bin);

        pi1_vector.push_back(pi1_Lorenz_vector);
        pi2_vector.push_back(pi2_Lorenz_vector);
      }


      //reset number of K0ss and vectors
      K0s_vector_background.clear();
      K0s_pT_bin_vector_background.clear();
      K0s_eta_bin_vector_background.clear();

      pi1_vector_background.clear();
      pi2_vector_background.clear();

      if( charge != 0 && cuts(K0s_y) && pT_bin_corr != -1 && eta_bin != -1)
      {
        K0s_vector_background.push_back(K0s_Lorentz_vector);
        K0s_pT_bin_vector_background.push_back(pT_bin_corr);
        K0s_eta_bin_vector_background.push_back(eta_bin);

        pi1_vector_background.push_back(pi1_Lorenz_vector);
        pi2_vector_background.push_back(pi2_Lorenz_vector);

      }
      //___________________________________________________________

    }//end else new event

  }//end loop over entries in NTuple

  //_________________________________________________________________________________________________________


  if(ReadMode == 0) //if in ReadMode = 0, save all histograms to output file
  {
    InvMassFile->cd();

    //hEventStat1->Write();

    for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
    {
      for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
      {
        K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->Write();
        K0s1_inv_mass_vs_K0s2_inv_mass_US_LS[pTbin1][pTbin2]->Write();
        K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->Write();

      }
    }

  } //end if RreadMode = 0



  //_____________________________________________________________________________________

  //projection Minv histograms
  TH1F *K0s1_inv_mass_K0sK0s_projection_US[nPtBins_corr][nPtBins_corr];
  TH1F *K0s1_inv_mass_K0sK0s_projection_LS[nPtBins_corr][nPtBins_corr];
  TH1F *K0s1_inv_mass_K0sK0s_projection[nPtBins_corr][nPtBins_corr];

  TH1F *K0s2_inv_mass_K0sK0s_projection_US[nPtBins_corr][nPtBins_corr];
  TH1F *K0s2_inv_mass_K0sK0s_projection_LS[nPtBins_corr][nPtBins_corr];
  TH1F *K0s2_inv_mass_K0sK0s_projection[nPtBins_corr][nPtBins_corr];


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

      //scale LS to match US
      sideBandScale_K0s_K0s[pTbin1][pTbin2] = K0s1_inv_mass_vs_K0s2_inv_mass_US[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh)/K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->Integral(binLow_sideBand,binHigh_sideBand, binLow,binHigh);

      K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->Sumw2();
      K0s1_inv_mass_vs_K0s2_inv_mass_LS[pTbin1][pTbin2]->Scale(sideBandScale_K0s_K0s[pTbin1][pTbin2]);

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

    }//end for eta

  }//end for pT


  InvMassFile->Close();

  //______________________________________________________________________________________________________________________


  //for extraction of reference polarization from K0s baseline
  const float L0_alpha = 0.732; //decay parameter of L0
  const float L0bar_alpha = -0.758; //decay paramteter of L0bar


  TFile *OutFile; //output file to store production plane histograms

  if(ReadMode == 0) //create production plane file from nTuple - run in this mode first
  {
    OutFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/ProdPlane_K0s_2D_work.root", "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT
  }
  else  //read file to store wrong and good sign M_inv distributions and daughter pT
  {

    OutFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/ProdPlane_K0s_2D_work.root", "read"); //old, non-optimized cuts

    if( !(OutFile->IsOpen()) )
    {
      cout<<"Unable to open file with production plane histograms!"<<endl;
      return false;
    }
  }


  //update efficiency files
  TFile *EffFile;

  if(energy == 510) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/K0s_cosThetaStar_eff_Run17_1B.root", "read");
  else if(energy == 200) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/K0s_cosThetaStar_eff_Run12_1B.root", "read");
  else
  {
    cout<<"Not a valid colliison energy! Abborting!"<<endl;
    return false;
  }

  //_______________________________________________________________________________________________________________________________________________



  //efficiency histograms
  TH1D* K0s_K0s_costhetaProdPlane_eff = (TH1D*)EffFile->Get("K0s_K0s_cosThetaProdPlane_eff");

  TH1D *K0s_K0s_cosThetaProdPlane_pT_eff[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_eta_eff[nEtaBins][nEtaBins];


  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      K0s_K0s_cosThetaProdPlane_pT_eff[pTbin1][pTbin2] = (TH1D*)EffFile->Get(Form("K0s_K0s_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      K0s_K0s_cosThetaProdPlane_eta_eff[etaBin1][etaBin2] = (TH1D*)EffFile->Get(Form("K0s_K0s_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
    }
  }
  //_____________________________________________________________________________

  //data histograms

  TH1D *K0s_K0s_cosThetaProdPlane_US_hist;
  TH1D *K0s_K0s_cosThetaProdPlane_LS_hist;
  TH1D *K0s_K0s_cosThetaProdPlane_ME_hist;
  TH1D *K0s_K0s_cosThetaProdPlane_ME_LS_hist;


  TH1D *K0s_K0s_cosThetaProdPlane_pT_US_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_pT_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_pT_ME_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[nPtBins_corr][nPtBins_corr];

  TH1D *K0s_K0s_cosThetaProdPlane_eta_US_hist[nEtaBins][nEtaBins];
  TH1D *K0s_K0s_cosThetaProdPlane_eta_LS_hist[nEtaBins][nEtaBins];
  TH1D *K0s_K0s_cosThetaProdPlane_eta_ME_hist[nEtaBins][nEtaBins];
  TH1D *K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[nEtaBins][nEtaBins];



  K0s_K0s_cosThetaProdPlane_US_hist = new TH1D("K0s_K0s_cosThetaProdPlane_US_hist", "K0s_K0s_cosThetaProdPlane_US_hist", 10, -1, 1);
  K0s_K0s_cosThetaProdPlane_LS_hist = new TH1D("K0s_K0s_cosThetaProdPlane_LS_hist", "K0s_K0s_cosThetaProdPlane_LS_hist", 10, -1, 1);
  K0s_K0s_cosThetaProdPlane_ME_hist = new TH1D("K0s_K0s_cosThetaProdPlane_ME_hist", "K0s_K0s_cosThetaProdPlane_ME_hist", 10, -1, 1);
  K0s_K0s_cosThetaProdPlane_ME_LS_hist = new TH1D("K0s_K0s_cosThetaProdPlane_ME_LS_hist", "K0s_K0s_cosThetaProdPlane_ME_LS_hist", 10, -1, 1);

  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_ME_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_ME_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_ME_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_ME_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
    }
  }



  //fill cos(theta*) histograms within the Minv window determined above

  //US combinations (signal + background)
  for(unsigned int K0s_K0s_index = 0; K0s_K0s_index < K0s_K0s_cos_theta.size(); K0s_K0s_index++)
  {
    float K0s1_peak_mean = fit_res_gaus[K0s_K0s_pT_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2.at(K0s_K0s_index)]->Parameter(1);
    float K0s1_peak_sigma = fabs(fit_res_gaus[K0s_K0s_pT_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2.at(K0s_K0s_index)]->Parameter(2));

    float K0s2_peak_mean = fit_res_gaus[K0s_K0s_pT_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2.at(K0s_K0s_index)]->Parameter(3);
    float K0s2_peak_sigma = fabs(fit_res_gaus[K0s_K0s_pT_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2.at(K0s_K0s_index)]->Parameter(4));


    if( K0s_K0s_Minv_K0s1.at(K0s_K0s_index) > K0s1_peak_mean-3*K0s1_peak_sigma && K0s_K0s_Minv_K0s1.at(K0s_K0s_index) < K0s1_peak_mean+3*K0s1_peak_sigma &&
        K0s_K0s_Minv_K0s2.at(K0s_K0s_index) > K0s2_peak_mean-3*K0s2_peak_sigma && K0s_K0s_Minv_K0s2.at(K0s_K0s_index) < K0s2_peak_mean+3*K0s2_peak_sigma)
    {
      K0s_K0s_cosThetaProdPlane_US_hist->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index));

      K0s_K0s_cosThetaProdPlane_pT_US_hist[K0s_K0s_pT_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2.at(K0s_K0s_index)]->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index));

      K0s_K0s_cosThetaProdPlane_eta_US_hist[K0s_K0s_eta_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_eta_bin_K0s2.at(K0s_K0s_index)]->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index));
    }

/*
    if(  K0s_K0s_Minv_K0s1.at(K0s_K0s_index) > K0s1_peak_mean+5*K0s1_peak_sigma && K0s_K0s_Minv_K0s2.at(K0s_K0s_index) > K0s2_peak_mean+5*K0s2_peak_sigma )
    {
      K0s_K0s_cosThetaProdPlane_US_side_band_hist->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index));

      K0s_K0s_cosThetaProdPlane_pT_US_side_band_hist[K0s_K0s_pT_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2.at(K0s_K0s_index)]->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index));

      K0s_K0s_cosThetaProdPlane_eta_US_side_band_hist[K0s_K0s_eta_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_eta_bin_K0s2.at(K0s_K0s_index)]->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index));
    }
*/

  }


  //background - US paired with LS
  for(unsigned int K0s_K0s_index = 0; K0s_K0s_index < K0s_K0s_cos_theta_back.size(); K0s_K0s_index++)
  {
    float K0s1_peak_mean = fit_res_gaus[K0s_K0s_pT_bin_K0s1_back.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2_back.at(K0s_K0s_index)]->Parameter(1);
    float K0s1_peak_sigma = fabs(fit_res_gaus[K0s_K0s_pT_bin_K0s1_back.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2_back.at(K0s_K0s_index)]->Parameter(2));

    float K0s2_peak_mean = fit_res_gaus[K0s_K0s_pT_bin_K0s1_back.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2_back.at(K0s_K0s_index)]->Parameter(3);
    float K0s2_peak_sigma = fabs(fit_res_gaus[K0s_K0s_pT_bin_K0s1_back.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2_back.at(K0s_K0s_index)]->Parameter(4));


    if( K0s_K0s_Minv_K0s1_back.at(K0s_K0s_index) > K0s1_peak_mean-3*K0s1_peak_sigma && K0s_K0s_Minv_K0s1_back.at(K0s_K0s_index) < K0s1_peak_mean+3*K0s1_peak_sigma &&
        K0s_K0s_Minv_K0s2_back.at(K0s_K0s_index) > K0s2_peak_mean-3*K0s2_peak_sigma && K0s_K0s_Minv_K0s2_back.at(K0s_K0s_index) < K0s2_peak_mean+3*K0s2_peak_sigma)
    {
      K0s_K0s_cosThetaProdPlane_LS_hist->Fill(K0s_K0s_cos_theta_back.at(K0s_K0s_index));

      K0s_K0s_cosThetaProdPlane_pT_LS_hist[K0s_K0s_pT_bin_K0s1_back.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2_back.at(K0s_K0s_index)]->Fill(K0s_K0s_cos_theta_back.at(K0s_K0s_index));

      K0s_K0s_cosThetaProdPlane_eta_LS_hist[K0s_K0s_eta_bin_K0s1_back.at(K0s_K0s_index)][K0s_K0s_eta_bin_K0s2_back.at(K0s_K0s_index)]->Fill(K0s_K0s_cos_theta_back.at(K0s_K0s_index));
    }

/*
    if(  K0s_K0s_Minv_K0s1.at(K0s_K0s_index) > K0s1_peak_mean+5*K0s1_peak_sigma && K0s_K0s_Minv_K0s2.at(K0s_K0s_index) > K0s2_peak_mean+5*K0s2_peak_sigma )
    {
      K0s_K0s_cosThetaProdPlane_US_side_band_hist->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index));

      K0s_K0s_cosThetaProdPlane_pT_US_side_band_hist[K0s_K0s_pT_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_pT_bin_K0s2.at(K0s_K0s_index)]->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index));

      K0s_K0s_cosThetaProdPlane_eta_US_side_band_hist[K0s_K0s_eta_bin_K0s1.at(K0s_K0s_index)][K0s_K0s_eta_bin_K0s2.at(K0s_K0s_index)]->Fill(K0s_K0s_cos_theta.at(K0s_K0s_index));
    }
*/

  }

  //mixed event
  for(unsigned int iK0s_ME_1 = 0; iK0s_ME_1 < K0s_vector_ME.size(); iK0s_ME_1++)
  {
    for(unsigned int iK0s_ME_2 = iK0s_ME_1+1; iK0s_ME_2 < K0s_vector_ME.size(); iK0s_ME_2++)
    {
      float K0s1_peak_mean = fit_res_gaus[K0s_pT_bin_vector_ME.at(iK0s_ME_1)][K0s_pT_bin_vector_ME.at(iK0s_ME_2)]->Parameter(1);
      float K0s1_peak_sigma = fabs(fit_res_gaus[K0s_pT_bin_vector_ME.at(iK0s_ME_1)][K0s_pT_bin_vector_ME.at(iK0s_ME_2)]->Parameter(2));

      float K0s2_peak_mean = fit_res_gaus[K0s_pT_bin_vector_ME.at(iK0s_ME_1)][K0s_pT_bin_vector_ME.at(iK0s_ME_2)]->Parameter(3);
      float K0s2_peak_sigma = fabs(fit_res_gaus[K0s_pT_bin_vector_ME.at(iK0s_ME_1)][K0s_pT_bin_vector_ME.at(iK0s_ME_2)]->Parameter(4));

      if( K0s_vector_ME.at(iK0s_ME_1).M() < K0s1_peak_mean-3*K0s1_peak_sigma || K0s_vector_ME.at(iK0s_ME_1).M() > K0s1_peak_mean+3*K0s1_peak_sigma ) continue;
      if( K0s_vector_ME.at(iK0s_ME_2).M() < K0s2_peak_mean-3*K0s2_peak_sigma || K0s_vector_ME.at(iK0s_ME_2).M() > K0s2_peak_mean+3*K0s2_peak_sigma ) continue;


      double K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_ME.at(iK0s_ME_1), pi_vector_ME.at(iK0s_ME_1), K0s_vector_ME.at(iK0s_ME_2), pi_vector_ME.at(iK0s_ME_2));

      K0s_K0s_cosThetaProdPlane_ME_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar));
      K0s_K0s_cosThetaProdPlane_pT_ME_hist[K0s_pT_bin_vector_ME.at(iK0s_ME_1)][K0s_pT_bin_vector_ME.at(iK0s_ME_2)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar));
      K0s_K0s_cosThetaProdPlane_eta_ME_hist[K0s_eta_bin_vector_ME.at(iK0s_ME_1)][K0s_eta_bin_vector_ME.at(iK0s_ME_2)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar));
    }
  }

  //ME for bacground - US paired with LS
  int nFill_ME_LS = 0;

  for(unsigned int iK0s_ME = 0; iK0s_ME < K0s_vector_ME.size(); iK0s_ME++)
  {
    for(unsigned int iK0s_ME_LS = 0; iK0s_ME_LS < K0s_vector_ME_LS.size(); iK0s_ME_LS++)
    {

      //fill US-LS and LS-US combinations
      if( nFill_ME_LS % 2 == 0 )
      {
        //for US
        float K0s1_peak_mean = fit_res_gaus[K0s_pT_bin_vector_ME.at(iK0s_ME)][K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)]->Parameter(1);
        float K0s1_peak_sigma = fabs(fit_res_gaus[K0s_pT_bin_vector_ME.at(iK0s_ME)][K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)]->Parameter(2));

        //for LS
        float K0s2_peak_mean = fit_res_gaus[K0s_pT_bin_vector_ME.at(iK0s_ME)][K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)]->Parameter(3);
        float K0s2_peak_sigma = fabs(fit_res_gaus[K0s_pT_bin_vector_ME.at(iK0s_ME)][K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)]->Parameter(4));

        //US
        if( K0s_vector_ME.at(iK0s_ME).M() < K0s1_peak_mean-3*K0s1_peak_sigma || K0s_vector_ME.at(iK0s_ME).M() > K0s1_peak_mean+3*K0s1_peak_sigma ) continue;
        //LS
        if( K0s_vector_ME_LS.at(iK0s_ME_LS).M() < K0s2_peak_mean-3*K0s2_peak_sigma || K0s_vector_ME_LS.at(iK0s_ME_LS).M() > K0s2_peak_mean+3*K0s2_peak_sigma ) continue;

        double K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_ME.at(iK0s_ME), pi_vector_ME.at(iK0s_ME), K0s_vector_ME_LS.at(iK0s_ME_LS), pi_vector_ME_LS.at(iK0s_ME_LS));

        K0s_K0s_cosThetaProdPlane_ME_LS_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar));
        K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[K0s_pT_bin_vector_ME.at(iK0s_ME)][K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar));
        K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[K0s_eta_bin_vector_ME.at(iK0s_ME)][K0s_eta_bin_vector_ME_LS.at(iK0s_ME_LS)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar));
      }
      else
      {
        //for LS
        float K0s1_peak_mean = fit_res_gaus[K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)][K0s_pT_bin_vector_ME.at(iK0s_ME)]->Parameter(1);
        float K0s1_peak_sigma = fabs(fit_res_gaus[K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)][K0s_pT_bin_vector_ME.at(iK0s_ME)]->Parameter(2));

        //for US
        float K0s2_peak_mean = fit_res_gaus[K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)][K0s_pT_bin_vector_ME.at(iK0s_ME)]->Parameter(3);
        float K0s2_peak_sigma = fabs(fit_res_gaus[K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)][K0s_pT_bin_vector_ME.at(iK0s_ME)]->Parameter(4));

        //US
        if( K0s_vector_ME.at(iK0s_ME).M() < K0s2_peak_mean-3*K0s2_peak_sigma || K0s_vector_ME.at(iK0s_ME).M() > K0s2_peak_mean+3*K0s2_peak_sigma ) continue;
        //LS
        if( K0s_vector_ME_LS.at(iK0s_ME_LS).M() < K0s1_peak_mean-3*K0s1_peak_sigma || K0s_vector_ME_LS.at(iK0s_ME_LS).M() > K0s1_peak_mean+3*K0s1_peak_sigma ) continue;

        double K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_ME_LS.at(iK0s_ME_LS), pi_vector_ME_LS.at(iK0s_ME_LS), K0s_vector_ME.at(iK0s_ME), pi_vector_ME.at(iK0s_ME) );

        K0s_K0s_cosThetaProdPlane_ME_LS_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar));
        K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[K0s_pT_bin_vector_ME_LS.at(iK0s_ME_LS)][K0s_pT_bin_vector_ME.at(iK0s_ME)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar));
        K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[K0s_eta_bin_vector_ME_LS.at(iK0s_ME_LS)][K0s_eta_bin_vector_ME.at(iK0s_ME)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar));
      }

      nFill_ME_LS++;

    }
  }



  //________________________________________________________________________________________

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  TCanvas *K0s_K0s_cosThetaProdPlane_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_can"), Form("K0s_K0s_cosThetaProdPlane_can"), 1200, 1000);

  K0s_K0s_cosThetaProdPlane_can->cd();

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

  K0s_K0s_cosThetaProdPlane_ME_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_ME_hist->SetMarkerColor(kOrange);
  K0s_K0s_cosThetaProdPlane_ME_hist->SetLineColor(kOrange);
  K0s_K0s_cosThetaProdPlane_ME_hist->Sumw2();
  //K0s_K0s_cosThetaProdPlane_ME_hist->Scale(nK0sK0s_back/K0s_K0s_cosThetaProdPlane_ME_hist->Integral()); //scale ME to match background
  K0s_K0s_cosThetaProdPlane_ME_hist->Scale(nK0sK0s/K0s_K0s_cosThetaProdPlane_ME_hist->Integral()); //scale ME to match signal+background
  K0s_K0s_cosThetaProdPlane_ME_hist->Scale(1./K0s_K0s_cosThetaProdPlane_ME_hist->GetXaxis()->GetBinWidth(1));

  //K0s_K0s_cosThetaProdPlane_US_hist->Divide(K0s_K0s_cosThetaProdPlane_ME_hist); //correction of signal+bacground using ME
  //K0s_K0s_cosThetaProdPlane_US_hist->Scale(nK0sK0s/K0s_K0s_cosThetaProdPlane_US_hist->Integral()); //scale back to US integral
  K0s_K0s_cosThetaProdPlane_US_hist->Scale(1./K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->GetBinWidth(1)); //bin width
  K0s_K0s_cosThetaProdPlane_US_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_US_hist->Draw("p e");

  K0s_K0s_cosThetaProdPlane_ME_hist->Draw("p e same"); //for plotting ME with uncorrected distributions


  //background
  K0s_K0s_cosThetaProdPlane_LS_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_LS_hist->SetMarkerColor(kBlue);
  double nK0sK0s_back = K0s_K0s_cosThetaProdPlane_LS_hist->Integral();
  K0s_K0s_cosThetaProdPlane_LS_hist->Sumw2();
  //K0s_K0s_cosThetaProdPlane_LS_hist->Scale(1./;
  //K0s_K0s_cosThetaProdPlane_LS_hist->Divide(K0s_K0s_costhetaProdPlane_eff);
  //K0s_K0s_cosThetaProdPlane_LS_hist->Draw("p e same");

  K0s_K0s_cosThetaProdPlane_ME_LS_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->SetMarkerColor(kOrange);
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->SetLineColor(kOrange);
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->Scale(nK0sK0s_back/K0s_K0s_cosThetaProdPlane_ME_LS_hist->Integral()); //scale ME to match background
  //K0s_K0s_cosThetaProdPlane_ME_LS_hist->Scale(nK0sK0s/K0s_K0s_cosThetaProdPlane_ME_LS_hist->Integral()); //scale ME to match signal+background
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->Scale(1./K0s_K0s_cosThetaProdPlane_ME_LS_hist->GetXaxis()->GetBinWidth(1));
  K0s_K0s_cosThetaProdPlane_ME_LS_hist->Draw("p e same");

  //K0s_K0s_cosThetaProdPlane_LS_hist->Divide(K0s_K0s_cosThetaProdPlane_ME_LS_hist); //correction of background using ME
  //K0s_K0s_cosThetaProdPlane_LS_hist->Scale(nK0sK0s_back/K0s_K0s_cosThetaProdPlane_LS_hist->Integral()); //scale back to raw background integral
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

  //__________________________________________________________________________________________________________________________________________________________________________

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

  float P_K0s_K0s_2 = fitK0s_K0s_US_ThetaStar_2->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_K0s_K0s_err_2 = fitK0s_K0s_US_ThetaStar_2->GetParError(1)/(L0_alpha*L0_alpha);

  fitK0s_K0s_US_ThetaStar_2->SetLineColor(1);
  fitK0s_K0s_US_ThetaStar_2->Draw("same");

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


  float nK0sK0s_pT[nPtBins_corr][nPtBins_corr];
  float nK0sK0s_back_pT[nPtBins_corr][nPtBins_corr];

  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      TCanvas *K0s_K0s_cosThetaProdPlane_pT_can = new TCanvas("K0s_K0s_cosThetaProdPlane_pT_can", "K0s_K0s_cosThetaProdPlane_pT_can", 1200, 1000);

      K0s_K0s_cosThetaProdPlane_pT_can->cd();

      //signal + background
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      nK0sK0s_pT[pTbin1][pTbin2] = K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral();
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Sumw2();
      //K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Add(K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2], -1);

      K0s_K0s_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      K0s_K0s_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(nK0sK0s_pT[pTbin1][pTbin2]/K0s_K0s_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral()); //scale ME to match US
      //K0s_K0s_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Draw("p e same");

      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Divide(K0s_K0s_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]); //correct US using ME
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(nK0sK0s_pT[pTbin1][pTbin2]/K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral()); //scale corrected back to US integral
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1)); //bin width
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMinimum(0);
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Draw("p e");


      //background
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      nK0sK0s_back_pT[pTbin1][pTbin2] = K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral();
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Sumw2();

      K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[pTbin1][pTbin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[pTbin1][pTbin2]->Scale(nK0sK0s_pT[pTbin1][pTbin2]/K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[pTbin1][pTbin2]->Integral()); //scale ME to match background
      //K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[pTbin1][pTbin2]->Draw("p e same");

      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Divide(K0s_K0s_cosThetaProdPlane_pT_ME_LS_hist[pTbin1][pTbin2]); //correct background using ME
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Scale(nK0sK0s_back_pT[pTbin1][pTbin2]/K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral()); //scale corrected back to background integral + bin width
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Scale(1./K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1)); //scale corrected back to background integral + bin width
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Draw("p e same");


      TF1 *fitK0s_K0s_pT_US_ThetaStar = new TF1("fitK0s_K0s_pT_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitK0s_K0s_pT_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitK0s_K0s_pT_US_ThetaStar, "s i 0 r");

      float P_K0s_K0s_pT = fitK0s_K0s_pT_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_K0s_K0s_pT_err = fitK0s_K0s_pT_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

      fitK0s_K0s_pT_US_ThetaStar->SetLineColor(1);
      fitK0s_K0s_pT_US_ThetaStar->Draw("same");


      TPaveText *K0s_K0s_pT_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      K0s_K0s_pT_text->SetTextFont(42);
      K0s_K0s_pT_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      K0s_K0s_pT_text->AddText("Minimum bias");
      K0s_K0s_pT_text->AddText("K_{s}^{0}-K_{s}^{0}");
      K0s_K0s_pT_text->AddText("|#eta| < 1");
      K0s_K0s_pT_text->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      K0s_K0s_pT_text->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      K0s_K0s_pT_text->AddText(Form("P = %.2f #pm %.2f", P_K0s_K0s_pT, P_K0s_K0s_pT_err));
      K0s_K0s_pT_text->SetFillColorAlpha(0, 0.01);
      K0s_K0s_pT_text->Draw("same");

      K0s_K0s_leg->Draw("same");

      K0s_K0s_cosThetaProdPlane_pT_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations/K0s_K0s_cosThetaProdPlane_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

      //________________________________________________________________________________________________________________________________________________________

      TCanvas *K0s_K0s_cosThetaProdPlane_pT_2_can = new TCanvas("K0s_K0s_cosThetaProdPlane_pT_2_can", "K0s_K0s_cosThetaProdPlane_pT_2_can", 1200, 1000);

      K0s_K0s_cosThetaProdPlane_pT_2_can->cd();

      //for background subtraction
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Add(K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2], -1);
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMinimum(0);
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Draw("p e");


      TF1 *fitK0s_K0s_pT_US_ThetaStar_2 = new TF1("fitK0s_K0s_pT_US_ThetaStar_2", "[0]*(1 + [1]*x)", -1, 1);
      fitK0s_K0s_pT_US_ThetaStar_2->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitK0s_K0s_pT_US_ThetaStar_2, "s i 0 r");

      float P_K0s_K0s_pT_2 = fitK0s_K0s_pT_US_ThetaStar_2->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_K0s_K0s_pT_err_2 = fitK0s_K0s_pT_US_ThetaStar_2->GetParError(1)/(L0_alpha*L0_alpha);

      fitK0s_K0s_pT_US_ThetaStar_2->SetLineColor(1);
      fitK0s_K0s_pT_US_ThetaStar_2->Draw("same");

      TPaveText *K0s_K0s_pT_text_2 = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      //cent_text_2->AddText("STAR preliminary");
      //((TText*)cent_text_2->GetListOfLines()->Last())->SetTextColor(2);
      K0s_K0s_pT_text_2->SetTextFont(42);
      K0s_K0s_pT_text_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      K0s_K0s_pT_text_2->AddText("Minimum bias");
      K0s_K0s_pT_text_2->AddText("K_{s}^{0}-K_{s}^{0}");
      K0s_K0s_pT_text_2->AddText("|#eta| < 1");
      K0s_K0s_pT_text_2->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      K0s_K0s_pT_text_2->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      K0s_K0s_pT_text_2->AddText(Form("P = %.2f #pm %.2f", P_K0s_K0s_pT_2, P_K0s_K0s_pT_err_2));
      K0s_K0s_pT_text_2->SetFillColorAlpha(0, 0.01);
      K0s_K0s_pT_text_2->Draw("same");

      K0s_K0s_2_leg->Draw("same");

      K0s_K0s_cosThetaProdPlane_pT_2_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations/K0s_K0s_cosThetaProdPlane_subtract_pT1_%i_pT2_%i.png", pTbin1, pTbin2));

    }
  }


  float nK0sK0s_eta[nEtaBins][nEtaBins];
  float nK0sK0s_back_eta[nEtaBins][nEtaBins];

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      TCanvas *K0s_K0s_cosThetaProdPlane_eta_can = new TCanvas("K0s_K0s_cosThetaProdPlane_eta_can", "K0s_K0s_cosThetaProdPlane_eta_can", 1200, 1000);

      K0s_K0s_cosThetaProdPlane_eta_can->cd();

      //signal + background
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      nK0sK0s_eta[etaBin1][etaBin2] = K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Integral();
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Sumw2();

      K0s_K0s_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      K0s_K0s_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Scale(nK0sK0s_eta[etaBin1][etaBin2]/K0s_K0s_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Integral());
      //K0s_K0s_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Draw("p e same");

      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Divide(K0s_K0s_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]);
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Scale(nK0sK0s_eta[etaBin1][etaBin2]/K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Integral()); //scale back to match raw US
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Scale(1./K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1)); //bin width
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMinimum(0);
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Draw("p e");


      //background
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      nK0sK0s_back_eta[etaBin1][etaBin2] = K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral();
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Sumw2();

      K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[etaBin1][etaBin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[etaBin1][etaBin2]->Scale(nK0sK0s_back_eta[etaBin1][etaBin2]/K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[etaBin1][etaBin2]->Integral());
      //K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[etaBin1][etaBin2]->Draw("p e same");

      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Divide(K0s_K0s_cosThetaProdPlane_eta_ME_LS_hist[etaBin1][etaBin2]);
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Scale(nK0sK0s_back_eta[etaBin1][etaBin2]/K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral());//scale back to match raw background
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Scale(1./K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1)); //bin width
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Draw("p e same");


      TF1 *fitK0s_K0s_eta_US_ThetaStar = new TF1("fitK0s_K0s_eta_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitK0s_K0s_eta_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_eta_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Fit(fitK0s_K0s_eta_US_ThetaStar, "s i 0 r");

      float P_K0s_K0s_eta = fitK0s_K0s_eta_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_K0s_K0s_eta_err = fitK0s_K0s_eta_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

      fitK0s_K0s_eta_US_ThetaStar->SetLineColor(1);
      fitK0s_K0s_eta_US_ThetaStar->Draw("same");

      TPaveText *K0s_K0s_eta_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      K0s_K0s_eta_text->SetTextFont(42);
      K0s_K0s_eta_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      K0s_K0s_eta_text->AddText("Minimum bias");
      K0s_K0s_eta_text->AddText("K_{s}^{0}-K_{s}^{0}");
      K0s_K0s_eta_text->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      K0s_K0s_eta_text->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      K0s_K0s_eta_text->AddText("p_{T} integrated");
      K0s_K0s_eta_text->AddText(Form("P = %.2f #pm %.2f", P_K0s_K0s_eta, P_K0s_K0s_eta_err));
      K0s_K0s_eta_text->SetFillColorAlpha(0, 0.01);
      K0s_K0s_eta_text->Draw("same");

      K0s_K0s_leg->Draw("same");

      K0s_K0s_cosThetaProdPlane_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations/K0s_K0s_cosThetaProdPlane_eta1_%i_eta2_%i.png", etaBin1, etaBin2));

      //_______________________________________________________________________________________________________________________________________________________

      TCanvas *K0s_K0s_cosThetaProdPlane_eta_2_can = new TCanvas("K0s_K0s_cosThetaProdPlane_eta_2_can", "K0s_K0s_cosThetaProdPlane_eta_2_can", 1200, 1000);

      K0s_K0s_cosThetaProdPlane_eta_2_can->cd();

      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Add(K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2], -1);
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMinimum(0);
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Draw("p e");

      TF1 *fitK0s_K0s_eta_US_ThetaStar_2 = new TF1("fitK0s_K0s_eta_US_ThetaStar_2", "[0]*(1 + [1]*x)", -1, 1);
      fitK0s_K0s_eta_US_ThetaStar_2->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_eta_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Fit(fitK0s_K0s_eta_US_ThetaStar_2, "s i 0 r");

      float P_K0s_K0s_eta_2 = fitK0s_K0s_eta_US_ThetaStar_2->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_K0s_K0s_eta_err_2 = fitK0s_K0s_eta_US_ThetaStar_2->GetParError(1)/(L0_alpha*L0_alpha);

      fitK0s_K0s_eta_US_ThetaStar_2->SetLineColor(1);
      fitK0s_K0s_eta_US_ThetaStar_2->Draw("same");

      TPaveText *K0s_K0s_eta_text_2 = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
      //cent_text_2->AddText("STAR preliminary");
      //((TText*)cent_text_2->GetListOfLines()->Last())->SetTextColor(2);
      K0s_K0s_eta_text_2->SetTextFont(42);
      K0s_K0s_eta_text_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      K0s_K0s_eta_text_2->AddText("Minimum bias");
      K0s_K0s_eta_text_2->AddText("K_{s}^{0}-K_{s}^{0}");
      K0s_K0s_eta_text_2->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      K0s_K0s_eta_text_2->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      K0s_K0s_eta_text_2->AddText("p_{T} integrated");
      K0s_K0s_eta_text_2->AddText(Form("P = %.2f #pm %.2f", P_K0s_K0s_eta_2, P_K0s_K0s_eta_err_2));
      K0s_K0s_eta_text_2->SetFillColorAlpha(0, 0.01);
      K0s_K0s_eta_text_2->Draw("same");

      K0s_K0s_2_leg->Draw("same");

      K0s_K0s_cosThetaProdPlane_eta_2_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/2D_analysis/K0s_correlations/K0s_K0s_cosThetaProdPlane_subtract_eta1_%i_eta2_%i.png", etaBin1, etaBin2));


    }
  }


  //----------------------------------Histograms, pT binnig, main analyisis---------------------------------

  cout<<endl;
  //cout<<"N events with K0s pair: "<<nEventsWithK0sPair_US_hist->GetBinContent(2)<<endl;
  cout<<"N K0s pairs from hist: "<<nK0sK0s<<endl;
  cout<<"N K0s background pairs from hist: "<<nK0sK0s_back<<endl;

  OutFile->Close();


  return true;

}
//__________________________________________________________________________________________________________________________

//ReadMode = 0 - read TTree, ReadMode = 1 - read histograms - First run in ReadMode = 0 to save relevant histograms, then can run in ReadMode = 1 to read just histograms and save time
//energy - collision energy in GeV
void Ana004_K0s_corr_2D(const int ReadMode = 0, const int energy = 510, const int year = 2017)
{
  ifstream fileList;

  if(energy == 510) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run17_MB_noTOF/fileList.list");
  else if(energy == 200) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF/fileList.list");
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


  bool AnaFinish = DoAnalysis(myChain, ReadMode, energy , year);

  if(!AnaFinish)
  {
    cout<<"Analysis of invariant spectra ended abnormally. Abborting!"<<endl;

    return;
  }

  cout<<endl;
  cout<<"Nubmer of accepted events: "<<hEventStat1->GetBinContent(6)<<endl;

  fileList.close();

  return;
}
