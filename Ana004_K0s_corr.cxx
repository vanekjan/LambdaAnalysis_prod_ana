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

bool cuts(int strictTOF_cut, int pi1_hasTOFinfo, int pi2_hasTOFinfo, float K0s_y)
{

  if( !( TMath::Abs(K0s_y) < K0s_y_cut ) ) return false;
  if( strictTOF_cut == 1 && pi1_hasTOFinfo == 0 ) return false; //TOF matched pions
  if( strictTOF_cut == 2 && (pi1_hasTOFinfo == 0 || pi2_hasTOFinfo == 0) ) return false; //TOF matched pions and protons
  //if(cos(K0s_theta) < K0s_cos_theta_cut) return false;
  //if(K0s_decayL > K0s_decayK0s_cut) return false;

  return true;

}

//analyze invariant mass spectra
//arguments are vectors of bin numbers of the invariant mass peak
//the bins are determined via fit to the invariant mass spectra
//bool InvMass(TTree *K0s_tree, vector<int> &invMassBins_L, vector<int> &invMassBins_L0, vector<int> &invMassBins_L0bar, const int readMode)
bool InvMass(TChain *K0s_tree, double (&invMassRange)[2][nPtBins+1][nEtaBins+1], double (&kappa)[nPtBins+1][nEtaBins+1], const int ReadMode, const int energy = 510, const int year = 2017)
{

  TFile *InvMassFile;

  if(ReadMode == 0) //create invariant mass file Form nTuple - run in this mode first
  {
    InvMassFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/InvMass_K0s_work.root", "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT
  }
  else  //read file to store wrong and good sign M_inv distributions and daughter pT
  {

    InvMassFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/InvMass_K0s_work.root", "read"); //old, non-optimized cuts

    if( !(InvMassFile->IsOpen()) )
    {
      cout<<"Unable to open file with invariant mass histograms!"<<endl;
      return false;
    }
  }



  //_______________________________________________________________________________________________________________________________________________


  Int_t charge;
  Float_t K0s_mass, K0s_pt, K0s_eta;//, K0s_decayL, K0s_theta, K0s_DCAdaughters;


  Float_t pi1_pt, pi2_pt;
  Float_t pi1_eta, pi2_eta;
  Float_t pi1_phi, pi2_phi;
  //Float_t pi1_ch;
  Int_t pi1_hasTOFinfo, pi2_hasTOFinfo;


  //Float_t thetaProdPlane;


  //---------------SET BARANCH ADDRESSES------------------------
  K0s_tree->SetBranchAddress("pair_charge", &charge);
  K0s_tree->SetBranchAddress("pair_mass", &K0s_mass);
  K0s_tree->SetBranchAddress("pair_pt", &K0s_pt);
  K0s_tree->SetBranchAddress("pair_eta", &K0s_eta);
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

  //K0s_tree->SetBranchAddress("thetaProdPlane", &thetaProdPlane);

  //--------------------------------------------------------------------------


  TH1D *K0s_inv_mass_US[nPtBins+1][nEtaBins+1];
  TH1D *K0s_inv_mass_LS[nPtBins+1][nEtaBins+1];


  if(ReadMode == 0) //create histograms to be saved into file
  {


    for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
    {
      for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
      {
        K0s_inv_mass_US[pTbin][etaBin] = new TH1D(Form("K0s_inv_mass_US_pT_%i_eta_%i", pTbin, etaBin), Form("K0s_inv_mass_US_pT_%i_eta_%i", pTbin, etaBin), 200, 0.45, 0.55);
        K0s_inv_mass_LS[pTbin][etaBin] = new TH1D(Form("K0s_inv_mass_LS_pT_%i_eta_%i", pTbin, etaBin), Form("K0s_inv_mass_LS_pT_%i_eta_%i", pTbin, etaBin), 200, 0.45, 0.55);
      }
    }
  }
  else //load histograms Form file
  {

    for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
    {
      for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
      {
        K0s_inv_mass_US[pTbin][etaBin] = (TH1D*)InvMassFile->Get(Form("K0s_inv_mass_US_pT_%i_eta_%i", pTbin, etaBin));
        K0s_inv_mass_LS[pTbin][etaBin] = (TH1D*)InvMassFile->Get(Form("K0s_inv_mass_LS_pT_%i_eta_%i", pTbin, etaBin));

      }
    }
  }
  //________________________________________________________________________________________


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



    double K0s_y = rapidity(K0s_pt, K0s_eta, K0s_mass_PDG);

    //------------------------------------------------------------------------------------------------------------------

    //cuts
    //one pi matched for Minv
    if( !cuts(1, pi1_hasTOFinfo, pi2_hasTOFinfo, K0s_y) ) continue;

    //----------------------------------------------------------------------------------------------------------------



    //fill all histograms for all pT and centrality bins
    int pT_bin = -1;

    //find pT bin of Dpm
    for(int j = 0; j < nPtBins; j++) //loop over pT bins
    {
      if(K0s_pt > pT_bins[j] && K0s_pt <= pT_bins[j+1])
      {
        pT_bin = j;
        break; //stop after pT bin is found
      }
    }

    if( pT_bin == -1 ) continue;


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

    if( eta_bin == -1 ) continue;



    if(charge == 0 ) //unlike-sign combinations
    {
      K0s_inv_mass_US[pT_bin][eta_bin]->Fill(K0s_mass);
      K0s_inv_mass_US[nPtBins][eta_bin]->Fill(K0s_mass);
      K0s_inv_mass_US[pT_bin][nEtaBins]->Fill(K0s_mass);
      K0s_inv_mass_US[nPtBins][nEtaBins]->Fill(K0s_mass); //pT integrated

    }
    else //like-sign combinations
    {
      K0s_inv_mass_LS[pT_bin][eta_bin]->Fill(K0s_mass);
      K0s_inv_mass_LS[nPtBins][eta_bin]->Fill(K0s_mass);
      K0s_inv_mass_LS[pT_bin][nEtaBins]->Fill(K0s_mass);
      K0s_inv_mass_LS[nPtBins][nEtaBins]->Fill(K0s_mass);

    }//end else for flag_int check
  }//end loop over entries in NTuple

  //_________________________________________________________________________________________________________


  if(ReadMode == 0) //if in ReadMode = 0, save all histograms to output file
  {
    InvMassFile->cd();

    //hEventStat1->Write();

    for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
    {
      for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
      {
        K0s_inv_mass_US[pTbin][etaBin]->Write();

        K0s_inv_mass_LS[pTbin][etaBin]->Write();
      }
    }

  } //end if RreadMode = 0

  //_____________________________________________________________________________________



  TFitResultPtr fit_res_gaus;

  for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
  {
    for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
    {


      TF1 *fitGauss = new TF1("fitGauss", "gaus(0)", 0.45, 0.55);

      if( pTbin == 0 ) fitGauss->SetParameters(10000, 1.116, 0.002);
      else
      {
        //cout<<"Valid"<<endl;
        fitGauss->SetParameters(fit_res_gaus->Parameter(0), fit_res_gaus->Parameter(1), fit_res_gaus->Parameter(2));
      }


      //fit_res_gaus_wrong_sign = K0s_inv_mass_US[pTbin][etaBin]->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      fit_res_gaus = K0s_inv_mass_US[pTbin][etaBin]->Fit(fitGauss, "s i l 0", "", 0.49, 0.51);


      if(!fit_res_gaus->IsValid())
      {
        cout<<"Fit not valid for pT "<<pTbin<<" eta "<<etaBin<<endl;
        return false;
      }




      float peak_range_min = fit_res_gaus->Parameter(1) - 3*fit_res_gaus->Parameter(2);
      float peak_range_max = fit_res_gaus->Parameter(1) + 3*fit_res_gaus->Parameter(2);


      invMassRange[0][pTbin][etaBin] = peak_range_min;
      invMassRange[1][pTbin][etaBin] = peak_range_max;


      int peak_range_min_bin = K0s_inv_mass_US[pTbin][etaBin]->GetXaxis()->FindBin(peak_range_min);
      int peak_range_max_bin = K0s_inv_mass_US[pTbin][etaBin]->GetXaxis()->FindBin(peak_range_max);

      //cout<<"hovno"<<endl;

      kappa[pTbin][etaBin] = K0s_inv_mass_US[pTbin][etaBin]->Integral(peak_range_min_bin, peak_range_max_bin)/K0s_inv_mass_LS[pTbin][etaBin]->Integral(peak_range_min_bin, peak_range_max_bin);


    }//end for eta

  }//end for pT

  return true;

}
//_______________________________________________________________________________________________________________________________________________________________________________________


//for analysis of K0s polarization
void K0sK0sSpinCorr(TChain *K0s_tree, double invMassRange[2][nPtBins+1][nEtaBins+1], double kappa[nPtBins+1][nEtaBins+1], const int ReadMode, const int energy = 510, const int year = 2017)
{
  //for extraction of reference polarization from K0s baseline
  const float L0_alpha = 0.732; //decay parameter of L0
  const float L0bar_alpha = -0.758; //decay paramteter of L0bar


  TFile *OutFile; //output file to store production plane histograms

  if(ReadMode == 0) //create production plane file from nTuple - run in this mode first
  {
    OutFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/ProdPlane_K0s_work.root", "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT
  }
  else  //read file to store wrong and good sign M_inv distributions and daughter pT
  {

    OutFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/ProdPlane_K0s_work.root", "read"); //old, non-optimized cuts

    if( !(OutFile->IsOpen()) )
    {
      cout<<"Unable to open file with production plane histograms!"<<endl;
      return;
    }
  }


  //update efficiency files
  TFile *EffFile;

  if(energy == 510) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/K0s_cosThetaStar_eff_Run17_1B.root", "read");
  else if(energy == 200) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/K0s_cosThetaStar_eff_Run12_1B.root", "read");
  else
  {
    cout<<"Not a valid colliison energy! Abborting!"<<endl;
    return;
  }

  //_______________________________________________________________________________________________________________________________________________


  //now not using binnig due to statistics

  TH1D *K0s_K0s_cosThetaProdPlane_US_hist;
  TH1D *K0s_K0s_cosThetaProdPlane_LS_hist;


  TH1D *K0s_K0s_cosThetaProdPlane_pT_US_hist[nPtBins_corr][nPtBins_corr];
  TH1D *K0s_K0s_cosThetaProdPlane_pT_LS_hist[nPtBins_corr][nPtBins_corr];

  TH1D *K0s_K0s_cosThetaProdPlane_eta_US_hist[nEtaBins][nEtaBins];
  TH1D *K0s_K0s_cosThetaProdPlane_eta_LS_hist[nEtaBins][nEtaBins];

  //________________________________________________________________

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

  //______________________________________________________________


  //TH1D *nK0sPairsInEvent_US_hist = new TH1D("nK0sPairsInEvent_US_hist", "nK0sPairsInEvent_US_hist", 10, 0, 10);
  TH1D *nK0ssInEvent_US_hist = new TH1D("nK0ssInEvent_US_hist", "nK0ssInEvent_US_hist", 10, 0, 10);

  TH1D *nEventsWithK0sPair_US_hist = new TH1D("nEventsWithK0sPair_US_hist", "nEventsWithK0sPair_US_hist", 2, 0, 2);

  //TH1D *nK0sPairsInEvent_LS_hist = new TH1D("nK0sPairsInEvent_LS_hist", "nK0sPairsInEvent_LS_hist", 10, 0, 10);
  TH1D *nK0ssInEvent_LS_hist = new TH1D("nK0ssInEvent_LS_hist", "nK0ssInEvent_LS_hist", 10, 0, 10);

  TH1D *nEventsWithK0sPair_LS_hist = new TH1D("nEventsWithK0sPair_LS_hist", "nEventsWithK0sPair_LS_hist", 2, 0, 2);


  if(ReadMode == 0) //create histograms to be saved into file
  {
    K0s_K0s_cosThetaProdPlane_US_hist = new TH1D("K0s_K0s_cosThetaProdPlane_US_hist", "K0s_K0s_cosThetaProdPlane_US_hist", 10, -1, 1);
    K0s_K0s_cosThetaProdPlane_LS_hist = new TH1D("K0s_K0s_cosThetaProdPlane_LS_hist", "K0s_K0s_cosThetaProdPlane_LS_hist", 10, -1, 1);

    for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
    {
      for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
      {
        K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
        K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("K0s_K0s_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      }
    }

    for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
    {
      for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
      {
        K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
        K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = new TH1D(Form("K0s_K0s_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("K0s_K0s_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      }
    }

  }
  else //load histograms from file
  {
    K0s_K0s_cosThetaProdPlane_US_hist = (TH1D*)OutFile->Get("K0s_K0s_cosThetaProdPlane_US_hist");
    K0s_K0s_cosThetaProdPlane_LS_hist = (TH1D*)OutFile->Get("K0s_K0s_cosThetaProdPlane_LS_hist");

    for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
    {
      for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
      {
        K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = (TH1D*)OutFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
        K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = (TH1D*)OutFile->Get(Form("K0s_K0s_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      }
    }

    for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
    {
      for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
      {
        K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = (TH1D*)OutFile->Get(Form("K0s_K0s_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
        K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = (TH1D*)OutFile->Get(Form("K0s_K0s_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      }
    }

  }
  //________________________________________________________________________________________

  //Labmda and K0s_bar stats

  int nK0s = -1; //default value

  int nK0s_background = -1; //default value

  int eventID_last = -1; //to store event ID from last L candidate
  int eventID_last_background = -1; //to store event ID from last L candidate

  Long64_t nEntries = 0; //total nEntries

  Int_t charge;
  Float_t K0s_mass, K0s_pt, K0s_eta, K0s_phi;//, K0s_decayL, K0s_theta, K0s_DCAdaughters;

  Float_t pi2_pt, pi1_pt;
  Float_t pi2_eta, pi1_eta;
  Float_t pi2_phi, pi1_phi;
  //Float_t pi2_ch;
  Float_t pi2_dca, pi1_dca;
  Int_t pi2_hasTOFinfo, pi1_hasTOFinfo;

  //Float_t thetaProdPlane;

  Int_t eventId;

  if(ReadMode == 0)
  {
    //new variable names
    //---------------SET BARANCH ADDRESSES------------------------
    K0s_tree->SetBranchAddress("pair_charge", &charge);
    K0s_tree->SetBranchAddress("pair_mass", &K0s_mass);
    K0s_tree->SetBranchAddress("pair_pt", &K0s_pt);
    K0s_tree->SetBranchAddress("pair_eta", &K0s_eta);
    K0s_tree->SetBranchAddress("pair_phi", &K0s_phi);
    //K0s_tree->SetBranchAddress("pair_decayL", &K0s_decayL);
    //K0s_tree->SetBranchAddress("pair_theta", &K0s_theta);
    //K0s_tree->SetBranchAddress("pair_DCAdaughters", &K0s_DCAdaughters);

    //pion 1
    K0s_tree->SetBranchAddress("p1_pt", &pi1_pt);
    K0s_tree->SetBranchAddress("p1_eta", &pi1_eta);
    K0s_tree->SetBranchAddress("p1_phi", &pi1_phi);
    //K0s_tree->SetBranchAddress("p1_dca", &pi1_dca);
    K0s_tree->SetBranchAddress("p1_hasTOFinfo", &pi1_hasTOFinfo);

    //pion 2
    K0s_tree->SetBranchAddress("p2_pt", &pi2_pt);
    K0s_tree->SetBranchAddress("p2_eta", &pi2_eta);
    K0s_tree->SetBranchAddress("p2_phi", &pi2_phi);
    //K0s_tree->SetBranchAddress("p2_ch", &pi2_ch);
    //K0s_tree->SetBranchAddress("p2_dca", &pi2_dca);
    K0s_tree->SetBranchAddress("p2_hasTOFinfo", &pi2_hasTOFinfo);

    //K0s_tree->SetBranchAddress("thetaProdPlane", &thetaProdPlane);

    K0s_tree->SetBranchAddress("eventId", &eventId);

    //--------------------------------------------------------------------------


    nEntries = K0s_tree->GetEntries();
    cout<<"nEntries = "<<nEntries<<endl;
  }

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



  for(Long64_t i = 0; i < nEntries; i++) //Read TTree only in ReadMode = 0
  {
    K0s_tree->GetEntry(i);

     //if(ReadMode != 0) break;
    if(i%1000000 == 0)
    {
      cout<<i<<endl;
    }

    //double K0s_xF = fabs(pz(K0s_pt, K0s_eta))/energy/2.; //energy is in CMS, need energz of one proton

    //calculate K0s rapidity y
    double K0s_y = rapidity(K0s_pt, K0s_eta, K0s_mass_PDG);

    if(nK0s == -1 ) //first iteration
    {
      eventID_last = eventId;
      eventID_last_background = eventId;

      nK0s = 0;

      nK0s_background = 0;
    }

    //------------------------------------------------------------------------------------------------------------------





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


    if(charge == 0 ) //like-sign combinations
    {

      if(eventId == eventID_last) //same event as in previous iteration and first event
      {
        if( cuts(0, pi1_hasTOFinfo, pi2_hasTOFinfo, K0s_y) && pT_bin != -1 && pT_bin_corr != -1 && eta_bin != -1)
        {
          if( K0s_mass > invMassRange[0][pT_bin][eta_bin] && K0s_mass < invMassRange[1][pT_bin][eta_bin] )
          {
            K0s_vector.push_back(K0s_Lorentz_vector);
            K0s_pT_bin_vector.push_back(pT_bin_corr);
            K0s_eta_bin_vector.push_back(eta_bin);

            pi1_vector.push_back(pi1_Lorenz_vector);
            pi2_vector.push_back(pi2_Lorenz_vector);

            nK0s++;
          }
        }

      }
      else if(eventId != eventID_last) //new event
      {
        //cout<<"new event"<<endl;

        //at least one L0-L0 pair in event
        if(nK0s > 1)
        {
          for(unsigned int iK0s1 = 0; iK0s1 < K0s_vector.size(); iK0s1++)
          {
            for(unsigned int iK0s2 = iK0s1+1; iK0s2 < K0s_vector.size(); iK0s2++)
            {
              //check auto-correlation
              if(pi1_vector.at(iK0s1).Phi() == pi1_vector.at(iK0s2).Phi()) continue;
              if(pi2_vector.at(iK0s1).Phi() == pi2_vector.at(iK0s2).Phi()) continue;
              if(pi1_vector.at(iK0s1).Phi() == pi2_vector.at(iK0s2).Phi()) continue;
              if(pi1_vector.at(iK0s2).Phi() == pi2_vector.at(iK0s1).Phi()) continue;

              //use pion 1 as reference particle
              double K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector.at(iK0s1), pi1_vector.at(iK0s1), K0s_vector.at(iK0s2), pi1_vector.at(iK0s2));

              K0s_K0s_cosThetaProdPlane_US_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar));

              K0s_K0s_cosThetaProdPlane_pT_US_hist[K0s_pT_bin_vector.at(iK0s1)][K0s_pT_bin_vector.at(iK0s2)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar));
              K0s_K0s_cosThetaProdPlane_eta_US_hist[K0s_eta_bin_vector.at(iK0s1)][K0s_eta_bin_vector.at(iK0s2)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar));

              nEventsWithK0sPair_US_hist->Fill(1.5);
            }
          }
        }
        else
        {
          nEventsWithK0sPair_US_hist->Fill(0.5);
        }

        nK0ssInEvent_US_hist->Fill(nK0s);

        //reset number of K0ss and vectors
        K0s_vector.clear();
        K0s_pT_bin_vector.clear();
        K0s_eta_bin_vector.clear();
        pi1_vector.clear();
        pi2_vector.clear();

        eventID_last = eventId;

        //check if we have good K0s in the new event
        if( cuts(0, pi1_hasTOFinfo, pi2_hasTOFinfo, K0s_y) && pT_bin != -1 && pT_bin_corr != -1 && eta_bin != -1)
        {
          if( K0s_mass > invMassRange[0][pT_bin][eta_bin] && K0s_mass < invMassRange[1][pT_bin][eta_bin] )
          {
            K0s_vector.push_back(K0s_Lorentz_vector);
            K0s_pT_bin_vector.push_back(pT_bin_corr);
            K0s_eta_bin_vector.push_back(eta_bin);

            pi1_vector.push_back(pi1_Lorenz_vector);
            pi2_vector.push_back(pi2_Lorenz_vector);

            nK0s = 1;
          }
          else
          {
            nK0s = 0;
          }
        }
        else
        {
          nK0s = 0;
        }

      }
    }
    else //unlike-sign combinations
    {
      if(eventId == eventID_last_background) //same event as in previous iteration
      {
        if( cuts(0, pi1_hasTOFinfo, pi2_hasTOFinfo, K0s_y) && pT_bin != -1 && pT_bin_corr != -1 && eta_bin != -1)
        {
          if( K0s_mass > invMassRange[0][pT_bin][eta_bin] && K0s_mass < invMassRange[1][pT_bin][eta_bin] )
          {
            K0s_vector_background.push_back(K0s_Lorentz_vector);
            K0s_pT_bin_vector_background.push_back(pT_bin_corr);
            K0s_eta_bin_vector_background.push_back(eta_bin);

            pi1_vector_background.push_back(pi1_Lorenz_vector);
            pi2_vector_background.push_back(pi2_Lorenz_vector);

            nK0s_background++;
          }
        }
      }
      else if(eventId != eventID_last_background) //new event
      {
        //at least one L0-L0 pair in event
        if(nK0s_background > 1)
        {
          for(unsigned int iK0s1 = 0; iK0s1 < K0s_vector_background.size(); iK0s1++)
          {
            for(unsigned int iK0s2 = iK0s1+1; iK0s2 < K0s_vector_background.size(); iK0s2++)
            {
              if(pi1_vector_background.at(iK0s1).Phi() == pi1_vector_background.at(iK0s2).Phi()) continue;
              if(pi2_vector_background.at(iK0s1).Phi() == pi2_vector_background.at(iK0s2).Phi()) continue;
              if(pi1_vector_background.at(iK0s1).Phi() == pi2_vector_background.at(iK0s2).Phi()) continue;
              if(pi1_vector_background.at(iK0s2).Phi() == pi2_vector_background.at(iK0s1).Phi()) continue;

              double K0s_K0s_pairThetaStar = LpairThetaStar(K0s_vector_background.at(iK0s1), pi1_vector_background.at(iK0s1), K0s_vector_background.at(iK0s2), pi1_vector_background.at(iK0s2));

              K0s_K0s_cosThetaProdPlane_LS_hist->Fill(TMath::Cos(K0s_K0s_pairThetaStar));

              K0s_K0s_cosThetaProdPlane_pT_LS_hist[K0s_pT_bin_vector_background.at(iK0s1)][K0s_pT_bin_vector_background.at(iK0s2)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar));
              K0s_K0s_cosThetaProdPlane_eta_LS_hist[K0s_eta_bin_vector_background.at(iK0s1)][K0s_eta_bin_vector_background.at(iK0s2)]->Fill(TMath::Cos(K0s_K0s_pairThetaStar));

              nEventsWithK0sPair_LS_hist->Fill(1.5);
            }
          }
        }
        else
        {
          nEventsWithK0sPair_LS_hist->Fill(0.5);
        }

        nK0ssInEvent_LS_hist->Fill(nK0s_background);

        //reset number of K0ss and vectors
        K0s_vector_background.clear();
        K0s_pT_bin_vector_background.clear();
        K0s_eta_bin_vector_background.clear();

        pi1_vector_background.clear();
        pi2_vector_background.clear();

        eventID_last_background = eventId;

        if( cuts(0, pi1_hasTOFinfo, pi2_hasTOFinfo, K0s_y) && pT_bin != -1 && pT_bin_corr != -1 && eta_bin != -1)
        {
          if( K0s_mass > invMassRange[0][pT_bin][eta_bin] && K0s_mass < invMassRange[1][pT_bin][eta_bin] )
          {
            K0s_vector_background.push_back(K0s_Lorentz_vector);
            K0s_pT_bin_vector_background.push_back(pT_bin_corr);
            K0s_eta_bin_vector_background.push_back(eta_bin);

            pi1_vector_background.push_back(pi1_Lorenz_vector);
            pi2_vector_background.push_back(pi2_Lorenz_vector);

            nK0s_background = 1;
          }
          else
          {
            nK0s_background = 0;
          }
        }
        else
        {
          nK0s_background = 0;
        }

      }

    }//end else for flag_int check
  }//end loop over entries in NTuple

  //________________________________________________________________________________________________________

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  TCanvas *K0s_K0s_cosThetaProdPlane_can = new TCanvas(Form("K0s_K0s_cosThetaProdPlane_can"), Form("K0s_K0s_cosThetaProdPlane_can"), 1200, 1000);

  K0s_K0s_cosThetaProdPlane_can->cd();

  K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  K0s_K0s_cosThetaProdPlane_US_hist->GetYaxis()->CenterTitle();
  K0s_K0s_cosThetaProdPlane_US_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_US_hist->SetMarkerColor(kRed);
  K0s_K0s_cosThetaProdPlane_US_hist->SetLineColor(kRed);
  double nK0sK0s = K0s_K0s_cosThetaProdPlane_US_hist->Integral();
  K0s_K0s_cosThetaProdPlane_US_hist->Sumw2();
  //K0s_K0s_cosThetaProdPlane_US_hist->Add(K0s_K0s_cosThetaProdPlane_LS_hist, -1);
  K0s_K0s_cosThetaProdPlane_US_hist->Scale(1./K0s_K0s_cosThetaProdPlane_US_hist->GetXaxis()->GetBinWidth(1));
  K0s_K0s_cosThetaProdPlane_US_hist->Divide(K0s_K0s_costhetaProdPlane_eff);
  //K0s_K0s_cosThetaProdPlane_US_hist->Scale(1./K0s_K0s_cosThetaProdPlane_US_hist->Integral());
  K0s_K0s_cosThetaProdPlane_US_hist->SetMinimum(0);
  K0s_K0s_cosThetaProdPlane_US_hist->Draw("p e");

  TF1 *fitK0s_K0s_US_ThetaStar = new TF1("fitK0s_K0s_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fitK0s_K0s_US_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  K0s_K0s_cosThetaProdPlane_US_hist->Fit(fitK0s_K0s_US_ThetaStar, "s i 0 r");

  float P_K0s_K0s = fitK0s_K0s_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_K0s_K0s_err = fitK0s_K0s_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

  fitK0s_K0s_US_ThetaStar->SetLineColor(1);
  fitK0s_K0s_US_ThetaStar->Draw("same");

  K0s_K0s_cosThetaProdPlane_LS_hist->SetMarkerStyle(20);
  K0s_K0s_cosThetaProdPlane_LS_hist->SetMarkerColor(kBlue);
  double nK0sK0s_back = K0s_K0s_cosThetaProdPlane_LS_hist->Integral();
  K0s_K0s_cosThetaProdPlane_LS_hist->Sumw2();
  K0s_K0s_cosThetaProdPlane_LS_hist->Scale(1./K0s_K0s_cosThetaProdPlane_LS_hist->GetXaxis()->GetBinWidth(1));
  //K0s_K0s_cosThetaProdPlane_LS_hist->Divide(K0s_K0s_costhetaProdPlane_eff);
  K0s_K0s_cosThetaProdPlane_LS_hist->Scale(1./K0s_K0s_cosThetaProdPlane_LS_hist->Integral());
  //K0s_K0s_cosThetaProdPlane_LS_hist->Draw("p e same");


  TPaveText *K0s_K0s_text = new TPaveText(0.35, 0.3, 0.75, 0.5, "NDC");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  K0s_K0s_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  K0s_K0s_text->AddText("Minimum bias");
  K0s_K0s_text->AddText("K_{s}^{0}-K_{s}^{0}");
  K0s_K0s_text->AddText("|#eta| < 1");
  K0s_K0s_text->AddText("p_{T} integrated");
  K0s_K0s_text->AddText(Form("P = %.2f #pm %.2f", P_K0s_K0s, P_K0s_K0s_err));
  K0s_K0s_text->SetFillColorAlpha(0, 0.01);
  K0s_K0s_text->Draw("same");

  K0s_K0s_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/K0s_correlations/K0s_K0s_cosThetaProdPlane.png");


  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      TCanvas *K0s_K0s_cosThetaProdPlane_pT_can = new TCanvas("K0s_K0s_cosThetaProdPlane_pT_can", "K0s_K0s_cosThetaProdPlane_pT_can", 1200, 1000);

      K0s_K0s_cosThetaProdPlane_pT_can->cd();

      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMarkerColor(kRed);
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetLineColor(kRed);
      //double nLL = K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral();
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Sumw2();
      //K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Add(K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2], -1);
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Divide(K0s_K0s_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      //K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral());
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMinimum(0);
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Draw("p e");

      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      //double nLL_back = K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral();
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Scale(1./K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      //K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Divide(K0s_K0s_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Scale(1./K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral());
      //K0s_K0s_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Draw("p e same");


      TF1 *fitK0s_K0s_pT_US_ThetaStar = new TF1("fitK0s_K0s_pT_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitK0s_K0s_pT_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      K0s_K0s_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitK0s_K0s_pT_US_ThetaStar, "s i 0 r");

      float P_K0s_K0s_pT = fitK0s_K0s_pT_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_K0s_K0s_pT_err = fitK0s_K0s_pT_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

      fitK0s_K0s_pT_US_ThetaStar->SetLineColor(1);
      fitK0s_K0s_pT_US_ThetaStar->Draw("same");


      TPaveText *K0s_K0s_pT_text = new TPaveText(0.35, 0.2, 0.75, 0.4, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      K0s_K0s_pT_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      K0s_K0s_pT_text->AddText("Minimum bias");
      K0s_K0s_pT_text->AddText("K_{s}^{0}-K_{s}^{0}");
      K0s_K0s_pT_text->AddText("|#eta| < 1");
      K0s_K0s_pT_text->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      K0s_K0s_pT_text->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      K0s_K0s_pT_text->AddText(Form("P = %.2f #pm %.2f", P_K0s_K0s_pT, P_K0s_K0s_pT_err));
      K0s_K0s_pT_text->SetFillColorAlpha(0, 0.01);
      K0s_K0s_pT_text->Draw("same");

      K0s_K0s_cosThetaProdPlane_pT_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/K0s_correlations/K0s_K0s_cosThetaProdPlane_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      TCanvas *K0s_K0s_cosThetaProdPlane_eta_can = new TCanvas("K0s_K0s_cosThetaProdPlane_eta_can", "K0s_K0s_cosThetaProdPlane_eta_can", 1200, 1000);

      K0s_K0s_cosThetaProdPlane_eta_can->cd();

      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMarkerColor(kRed);
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetLineColor(kRed);
      //double nLL = K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Integral();
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Sumw2();
      //K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Add(K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2], -1);
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Scale(1./K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Divide(K0s_K0s_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      //K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Scale(1./K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Integral());
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->SetMinimum(0);
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Draw("p e");

      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      //double nLL_back = K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral();
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Sumw2();
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Scale(1./K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      //K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Divide(K0s_K0s_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Scale(1./K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral());
      //K0s_K0s_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Draw("p e same");


      TF1 *fitK0s_K0s_eta_US_ThetaStar = new TF1("fitK0s_K0s_eta_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitK0s_K0s_eta_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_eta_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      K0s_K0s_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Fit(fitK0s_K0s_eta_US_ThetaStar, "s i 0 r");

      float P_K0s_K0s_eta = fitK0s_K0s_eta_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_K0s_K0s_eta_err = fitK0s_K0s_eta_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

      fitK0s_K0s_eta_US_ThetaStar->SetLineColor(1);
      fitK0s_K0s_eta_US_ThetaStar->Draw("same");

      TPaveText *K0s_K0s_eta_text = new TPaveText(0.35, 0.2, 0.75, 0.4, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      K0s_K0s_eta_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      K0s_K0s_eta_text->AddText("Minimum bias");
      K0s_K0s_eta_text->AddText("K_{s}^{0}-K_{s}^{0}");
      K0s_K0s_eta_text->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      K0s_K0s_eta_text->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      K0s_K0s_eta_text->AddText("p_{T} integrated");
      K0s_K0s_eta_text->AddText(Form("P = %.2f #pm %.2f", P_K0s_K0s_eta, P_K0s_K0s_eta_err));
      K0s_K0s_eta_text->SetFillColorAlpha(0, 0.01);
      K0s_K0s_eta_text->Draw("same");

      K0s_K0s_cosThetaProdPlane_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/K0s_correlations/K0s_K0s_cosThetaProdPlane_eta1_%i_eta2_%i.png", etaBin1, etaBin2));
    }
  }


  //----------------------------------Histograms, pT binnig, main analyisis---------------------------------

  cout<<endl;
  cout<<"N events with K0s pair: "<<nEventsWithK0sPair_US_hist->GetBinContent(2)<<endl;
  cout<<"N K0s pairs from hist: "<<nK0sK0s<<endl;
  cout<<"N K0s background pairs from hist: "<<nK0sK0s_back<<endl;

  OutFile->Close();

  return;

}
//__________________________________________________________________________________________________________________________

//ReadMode = 0 - read TTree, ReadMode = 1 - read histograms - First run in ReadMode = 0 to save relevant histograms, then can run in ReadMode = 1 to read just histograms and save time
//energy - collision energy in GeV
void Ana004_K0s_corr(const int ReadMode = 0, const int energy = 510, const int year = 2017)
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


  //arrays to store invariant mass peak ranges
  //2 is for lower and upper range
  double invMassRange[2][nPtBins+1][nEtaBins+1];


  //K0s signal to signal+background ratio
  //for bacground subtraction for polarization measurement
  double kappa[nPtBins+1][nEtaBins+1];



  bool invMassFinish = InvMass(myChain, invMassRange, kappa, ReadMode, energy , year);

  if(!invMassFinish)
  {
    cout<<"Analysis of invariant spectra ended abnormally. Abborting!"<<endl;

    return;
  }

  K0sK0sSpinCorr(myChain, invMassRange, kappa, ReadMode, energy , year);

  cout<<endl;
  cout<<"Nubmer of accepted events: "<<hEventStat1->GetBinContent(6)<<endl;

  fileList.close();

  return;
}
