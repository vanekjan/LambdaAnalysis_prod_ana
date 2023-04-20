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

bool cuts(int strictTOF_cut, int pi_hasTOFinfo, int p_hasTOFinfo, float L_y)
{

  if( !( TMath::Abs(L_y) < L_y_cut ) ) return false;
  if( strictTOF_cut == 1 && pi_hasTOFinfo == 0 ) return false; //TOF matched pions
  if( strictTOF_cut == 2 && (pi_hasTOFinfo == 0 || p_hasTOFinfo == 0) ) return false; //TOF matched pions and protons
  //if(cos(L_theta) < L_cos_theta_cut) return false;
  //if(L_decayL > L_decayL_cut) return false;

  return true;

}

//analyze invariant mass spectra
//arguments are vectors of bin numbers of the invariant mass peak
//the bins are determined via fit to the invariant mass spectra
//bool InvMass(TTree *L_tree, vector<int> &invMassBins_L, vector<int> &invMassBins_L0, vector<int> &invMassBins_L0bar, const int readMode)
bool InvMass(TChain *L_tree, double (&invMassRange_L)[2][nPtBins+1][nEtaBins+1], double (&invMassRange_L0)[2][nPtBins+1][nEtaBins+1], double (&invMassRange_L0bar)[2][nPtBins+1][nEtaBins+1], double (&L0_kappa)[nPtBins+1][nEtaBins+1], double (&L0bar_kappa)[nPtBins+1][nEtaBins+1], const int ReadMode, const int trigger = 0,const int energy = 510, const int year = 2017)
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

  TH1D *L_inv_mass_US[nPtBins+1][nEtaBins+1];
  TH1D *L_inv_mass_LS[nPtBins+1][nEtaBins+1];
  TH1D *L_inv_mass[nPtBins+1][nEtaBins+1];

  TH1D *L0_inv_mass_US[nPtBins+1][nEtaBins+1];
  TH1D *L0_inv_mass_LS[nPtBins+1][nEtaBins+1];
  TH1D *L0_inv_mass[nPtBins+1][nEtaBins+1];

  TH1D *L0bar_inv_mass_US[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_inv_mass_LS[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_inv_mass[nPtBins+1][nEtaBins+1];



  if(ReadMode == 0) //create histograms to be saved into file
  {
    for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
    {
      for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
      {

        //invariant mass histograms
        L_inv_mass_US[pTbin][etaBin] = new TH1D(Form("L_inv_mass_US_pT_%i_eta_%i", pTbin, etaBin), Form("L_inv_mass_US_pT_%i_eta_%i", pTbin, etaBin), 200, 1, 1.2);
        L_inv_mass_LS[pTbin][etaBin] = new TH1D(Form("L_inv_mass_LS_pT_%i_eta_%i", pTbin, etaBin), Form("L_inv_mass_LS_pT_%i_eta_%i", pTbin, etaBin), 200, 1, 1.2);

        L0_inv_mass_US[pTbin][etaBin] = new TH1D(Form("L0_inv_mass_US_pT_%i_eta_%i", pTbin, etaBin), Form("L0_inv_mass_US_pT_%i_eta_%i", pTbin, etaBin), 200, 1, 1.2);
        L0_inv_mass_LS[pTbin][etaBin] = new TH1D(Form("L0_inv_mass_LS_pT_%i_eta_%i", pTbin, etaBin), Form("L0_inv_mass_LS_pT_%i_eta_%i", pTbin, etaBin), 200, 1, 1.2);

        L0bar_inv_mass_US[pTbin][etaBin] = new TH1D(Form("L0bar_inv_mass_US_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_inv_mass_US_pT_%i_eta_%i", pTbin, etaBin), 200, 1, 1.2);
        L0bar_inv_mass_LS[pTbin][etaBin] = new TH1D(Form("L0bar_inv_mass_LS_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_inv_mass_LS_pT_%i_eta_%i", pTbin, etaBin), 200, 1, 1.2);
        //______________________________________________________________________________________________________________________________

      }
    }
  }
  else //load histograms from file
  {
    for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
    {
      for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
      {
        L_inv_mass_US[pTbin][etaBin] = (TH1D*)InvMassFile->Get(Form("L_inv_mass_US_pT_%i_eta_%i", pTbin, etaBin));
        L_inv_mass_LS[pTbin][etaBin] = (TH1D*)InvMassFile->Get(Form("L_inv_mass_LS_pT_%i_eta_%i", pTbin, etaBin));

        L0_inv_mass_US[pTbin][etaBin] = (TH1D*)InvMassFile->Get(Form("L0_inv_mass_US_pT_%i_eta_%i", pTbin, etaBin));
        L0_inv_mass_LS[pTbin][etaBin] = (TH1D*)InvMassFile->Get(Form("L0_inv_mass_LS_pT_%i_eta_%i", pTbin, etaBin));

        L0bar_inv_mass_US[pTbin][etaBin] = (TH1D*)InvMassFile->Get(Form("L0bar_inv_mass_US_pT_%i_eta_%i", pTbin, etaBin));
        L0bar_inv_mass_LS[pTbin][etaBin] = (TH1D*)InvMassFile->Get(Form("L0bar_inv_mass_LS_pT_%i_eta_%i", pTbin, etaBin));

      }
    }
  }
  //________________________________________________________________________________________



  Long64_t nEntries = 0; //total nEntries


  Int_t charge;
  Float_t L_mass, L_pt, L_eta; //, L_decayL, L_theta, L_DCAdaughters;

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



  for(Long64_t i = 0; i < nEntries; i++) //Read TTree only in ReadMode = 0
  {
    L_tree->GetEntry(i);

      //if(ReadMode != 0) break;
    if(i%1000000 == 0)
    {
      cout<<i<<endl;
    }

    //calculate Lambda rapidity y

    //double L_pz = pz(L_pt, L_eta);

    double L_y = rapidity(L_pt, L_eta, L_mass_PDG);



    //cuts
    //always require at least pion matched to TOF for imvariant mass
    //need to suppress pileup to nice peaks and determine M_inv window for each bin
    if( !cuts(1, pi_hasTOFinfo, p_hasTOFinfo, L_y) ) continue;



    //fill all histograms for all pT and centrality bins
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

    if( pT_bin == -1 ) continue;
    //----------------------------------------------------------

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
    //-----------------------------------------------------------

    if(charge == 0 ) //like-sign combinations
    {
      L_inv_mass_US[pT_bin][eta_bin]->Fill(L_mass);
      L_inv_mass_US[nPtBins][eta_bin]->Fill(L_mass);
      L_inv_mass_US[pT_bin][nEtaBins]->Fill(L_mass);
      L_inv_mass_US[nPtBins][nEtaBins]->Fill(L_mass); //pT integrated

      if( p_ch == 1)
      {
        L0_inv_mass_US[pT_bin][eta_bin]->Fill(L_mass);
        L0_inv_mass_US[nPtBins][eta_bin]->Fill(L_mass);
        L0_inv_mass_US[pT_bin][nEtaBins]->Fill(L_mass);
        L0_inv_mass_US[nPtBins][nEtaBins]->Fill(L_mass);
      }

      if(p_ch == -1)
      {
        L0bar_inv_mass_US[pT_bin][eta_bin]->Fill(L_mass);
        L0bar_inv_mass_US[nPtBins][eta_bin]->Fill(L_mass);
        L0bar_inv_mass_US[pT_bin][nEtaBins]->Fill(L_mass);
        L0bar_inv_mass_US[nPtBins][nEtaBins]->Fill(L_mass);
      }


      }
      else //unlike-sign combinations
      {


        L_inv_mass_LS[pT_bin][eta_bin]->Fill(L_mass);
        L_inv_mass_LS[nPtBins][eta_bin]->Fill(L_mass);
        L_inv_mass_LS[pT_bin][nEtaBins]->Fill(L_mass);
        L_inv_mass_LS[nPtBins][nEtaBins]->Fill(L_mass);



        if( p_ch == 1)
        {
          L0_inv_mass_LS[pT_bin][eta_bin]->Fill(L_mass);
          L0_inv_mass_LS[nPtBins][eta_bin]->Fill(L_mass);
          L0_inv_mass_LS[pT_bin][nEtaBins]->Fill(L_mass);
          L0_inv_mass_LS[nPtBins][nEtaBins]->Fill(L_mass);
        }

        if(p_ch == -1)
        {
          L0bar_inv_mass_LS[pT_bin][eta_bin]->Fill(L_mass);
          L0bar_inv_mass_LS[nPtBins][eta_bin]->Fill(L_mass);
          L0bar_inv_mass_LS[pT_bin][nEtaBins]->Fill(L_mass);
          L0bar_inv_mass_LS[nPtBins][nEtaBins]->Fill(L_mass);
        }

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
          L_inv_mass_US[pTbin][etaBin]->Write();
          L_inv_mass_US[pTbin][etaBin]->Write();

          L_inv_mass_LS[pTbin][etaBin]->Write();
          L_inv_mass_LS[pTbin][etaBin]->Write();
        }
      }

      //InvMassFile->Close();

    } //end if RreadMode = 0


    //______________________________________________________________________________


    TCanvas *InvMassCan_US_LS[nPtBins+1][nEtaBins+1];
    TCanvas *InvMassCan[nPtBins+1][nEtaBins+1];

    TCanvas *InvMassCanL0_US_LS[nPtBins+1][nEtaBins+1];
    TCanvas *InvMassCanL0[nPtBins+1][nEtaBins+1];

    TCanvas *InvMassCanL0bar_US_LS[nPtBins+1][nEtaBins+1];
    TCanvas *InvMassCanL0bar[nPtBins+1][nEtaBins+1];


    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TFitResultPtr fit_res_gaus_L;

    TFitResultPtr fit_res_gaus_L0;

    TFitResultPtr fit_res_gaus_L0bar;


    for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
    {

      for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
      {

        //cout<< nEtaBins+1<<endl;

        InvMassCan_US_LS[pTbin][etaBin] = new TCanvas(Form("InvMassCan_US_LS_pT_%i_eta_%i", pTbin, etaBin), Form("InvMassCan_US_LS_pT_%i_eta_%i", pTbin, etaBin), 1300, 1000);
        InvMassCan[pTbin][etaBin] = new TCanvas(Form("InvMassCan_pT_%i_eta_%i", pTbin, etaBin), Form("InvMassCan_pT_%i_eta_%i", pTbin, etaBin), 1300, 1000);

        InvMassCanL0_US_LS[pTbin][etaBin] = new TCanvas(Form("L0_InvMassCan_US_LS_pT_%i_eta_%i", pTbin, etaBin), Form("L0_InvMassCan_US_LS_pT_%i_eta_%i", pTbin, etaBin), 1300, 1000);
        InvMassCanL0[pTbin][etaBin] = new TCanvas(Form("L0_InvMassCan_pT_%i_eta_%i", pTbin, etaBin), Form("L0_InvMassCan_pT_%i_eta_%i", pTbin, etaBin), 1300, 1000);

        InvMassCanL0bar_US_LS[pTbin][etaBin] = new TCanvas(Form("L0bar_InvMassCan_US_LS_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_InvMassCan_US_LS_pT_%i_eta_%i", pTbin, etaBin), 1300, 1000);
        InvMassCanL0bar[pTbin][etaBin] = new TCanvas(Form("L0bar_InvMassCan_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_InvMassCan_pT_%i_eta_%i", pTbin, etaBin), 1300, 1000);


        TString *pT_range = new TString();
        if(pTbin < nPtBins) pT_range->Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c}", pT_bins[pTbin], pT_bins[pTbin+1]);
        else pT_range->Form("0.5 < p_{T} < 5.0 GeV/#it{c}");


        TString *eta_range = new TString();
        if(etaBin < nEtaBins) eta_range->Form("%0.1f < #eta < %0.1f", eta_bins[etaBin], eta_bins[etaBin+1]);
        else eta_range->Form("-1 < #eta < 1");



        TF1 *fitGauss_L = new TF1("fitGauss_L", "gaus(0)", 1.07, 2.);

        TF1 *fitGauss_L0 = new TF1("fitGaussL0", "gaus(0)", 1.07, 2.);

        TF1 *fitGauss_L0bar = new TF1("fitGauss_L0bar", "gaus(0)", 1.07, 2.);


        InvMassCan_US_LS[pTbin][etaBin]->cd();

        gPad->SetLeftMargin(0.1);

        L_inv_mass_US[pTbin][etaBin]->SetMarkerStyle(20);
        L_inv_mass_US[pTbin][etaBin]->SetMarkerColor(kRed);
        L_inv_mass_US[pTbin][etaBin]->SetLineColor(kRed);
        L_inv_mass_US[pTbin][etaBin]->GetXaxis()->SetTitle("M_{inv}^{p#pi} (GeV/#it{c}^{2})");
        L_inv_mass_US[pTbin][etaBin]->GetXaxis()->CenterTitle();
        L_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
        L_inv_mass_US[pTbin][etaBin]->GetYaxis()->CenterTitle();
        //L_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 7000);
        //L_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 3000);
        L_inv_mass_US[pTbin][etaBin]->Draw("p e");

        L_inv_mass_LS[pTbin][etaBin]->SetMarkerStyle(24);
        L_inv_mass_LS[pTbin][etaBin]->SetMarkerColor(kBlue);
        L_inv_mass_LS[pTbin][etaBin]->SetLineColor(kBlue);
        L_inv_mass_LS[pTbin][etaBin]->Draw("p e same");

        TLegend *SigAndBckgLeg = new TLegend(0.15, 0.3, 0.4, 0.45 );
        SigAndBckgLeg->SetTextFont(42);
        SigAndBckgLeg->AddEntry(L_inv_mass_US[pTbin][etaBin], "Unlike-sign");
        SigAndBckgLeg->AddEntry(L_inv_mass_LS[pTbin][etaBin], "Like-sign");
        SigAndBckgLeg->SetBorderSize(0);
        SigAndBckgLeg->Draw("same");

        TPaveText *cent_text = new TPaveText(0.15, 0.5, 0.45, 0.8, "NDC");
        cent_text->SetTextFont(42);
        //cent_text->AddText("STAR preliminary");
        //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
        cent_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
        if(trigger == 0 || trigger == 2) cent_text->AddText("Minimum bias");
        if(trigger == 1) cent_text->AddText("JP2");
        cent_text->AddText("#Lambda^{0} and #bar{#Lambda^{0}}");
        cent_text->AddText(eta_range->Data());
        cent_text->AddText(pT_range->Data());
        cent_text->SetFillColorAlpha(0, 0.01);
        cent_text->Draw("same");

        InvMassCan_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_inv_mass_US_LS_pT_%i_eta_%i.png", pTbin, etaBin));
        //InvMassCan_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_inv_mass_US_LS_pT_%i_eta_%i.pdf", pTbin, etaBin));


        InvMassCan[pTbin][etaBin]->cd();

        gPad->SetLeftMargin(0.1);

        L_inv_mass[pTbin][etaBin] = (TH1D*)L_inv_mass_US[pTbin][etaBin]->Clone();
        L_inv_mass[pTbin][etaBin]->Sumw2();
        L_inv_mass[pTbin][etaBin]->Add(L_inv_mass_LS[pTbin][etaBin], -1);
        L_inv_mass[pTbin][etaBin]->GetXaxis()->SetTitle("M_{inv}^{p#pi} (GeV/#it{c}^{2})");
        L_inv_mass[pTbin][etaBin]->GetXaxis()->CenterTitle();
        L_inv_mass[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
        L_inv_mass[pTbin][etaBin]->GetYaxis()->CenterTitle();
        L_inv_mass[pTbin][etaBin]->SetMarkerStyle(20);
        L_inv_mass[pTbin][etaBin]->SetMarkerColor(kRed);
        L_inv_mass[pTbin][etaBin]->SetLineColor(kRed);
        L_inv_mass[pTbin][etaBin]->Draw("p e");


        if(pTbin > 0)//skip first pT bin - no peak visisble
        {
          if(pTbin == 1) fitGauss_L->SetParameters(2000, 1.116, 0.002);
          else fitGauss_L->SetParameters(fit_res_gaus_L->Parameter(0), fit_res_gaus_L->Parameter(1), fit_res_gaus_L->Parameter(2));


          fit_res_gaus_L = L_inv_mass[pTbin][etaBin]->Fit(fitGauss_L, "s i l 0", "", 1.11, 1.12);

          if( !fit_res_gaus_L->IsValid())
          {
            cout<<"Fit not valid for L bin pT = "<<pTbin<<" eta = "<<etaBin<<endl;
            return false;
          }

          fitGauss_L->SetLineColor(1);
          fitGauss_L->Draw("same");


          invMassRange_L[0][pTbin][etaBin] = fit_res_gaus_L->Parameter(1) - 3*fit_res_gaus_L->Parameter(2);
          invMassRange_L[1][pTbin][etaBin] = fit_res_gaus_L->Parameter(1) + 3*fit_res_gaus_L->Parameter(2);

        }
        else
        {
          invMassRange_L[0][pTbin][etaBin] = 0; //min for pT bin == 0
          invMassRange_L[1][pTbin][etaBin] = 0; //max for pT bin == 0
        }



        TF1 *zeroLine = new TF1("zeroLine", "[0]", 1, 1.2);
        zeroLine->SetParameter(0, 0);
        zeroLine->SetLineColor(1);
        zeroLine->SetLineStyle(9);
        zeroLine->Draw("same");

        TLegend *SigMinusBckgLeg = new TLegend(0.15, 0.3, 0.4, 0.45 );
        SigMinusBckgLeg->AddEntry(L_inv_mass[pTbin][etaBin], "US-LS");
        SigMinusBckgLeg->AddEntry(fitGauss_L, "Fit: Gauss");
        SigMinusBckgLeg->SetBorderSize(0);
        SigMinusBckgLeg->Draw("same");

        cent_text->Draw("same");

        InvMassCan[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_inv_mass_signal_pT_%i_eta_%i.png", pTbin, etaBin));
        //InvMassCan[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_inv_mass_signal_pT_%i_eta_%i.pdf", pTbin, etaBin));

        //------------------------------------------------------------------------------------------------------------



        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);


        InvMassCanL0_US_LS[pTbin][etaBin]->cd();

        gPad->SetLeftMargin(0.1);

        L0_inv_mass_US[pTbin][etaBin]->SetMarkerStyle(20);
        L0_inv_mass_US[pTbin][etaBin]->SetMarkerColor(kRed);
        L0_inv_mass_US[pTbin][etaBin]->SetLineColor(kRed);
        L0_inv_mass_US[pTbin][etaBin]->GetXaxis()->SetTitle("M_{inv}^{p#pi^{-}} (GeV/#it{c}^{2})");
        L0_inv_mass_US[pTbin][etaBin]->GetXaxis()->CenterTitle();
        L0_inv_mass_US[pTbin][etaBin]->GetXaxis()->SetTitleSize(0.04);
        L0_inv_mass_US[pTbin][etaBin]->GetXaxis()->SetLabelSize(0.04);
        L0_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
        L0_inv_mass_US[pTbin][etaBin]->GetYaxis()->CenterTitle();
        L0_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetTitleSize(0.04);
        L0_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetLabelSize(0.04);
        L0_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetMaxDigits(3);
        //L0_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 7000);
        //L0_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 3000);
        L0_inv_mass_US[pTbin][etaBin]->Draw("p e");

        L0_inv_mass_LS[pTbin][etaBin]->SetMarkerStyle(24);
        L0_inv_mass_LS[pTbin][etaBin]->SetMarkerColor(kBlue);
        L0_inv_mass_LS[pTbin][etaBin]->SetLineColor(kBlue);
        L0_inv_mass_LS[pTbin][etaBin]->Draw("p e same");


        SigAndBckgLeg->Draw("same");

        TPaveText *cent_text_2 = new TPaveText(0.15, 0.5, 0.45, 0.8, "NDC");
        cent_text_2->SetTextFont(42);
        //cent_text_2->AddText("STAR Internal");
        //cent_text->AddText("STAR preliminary");
        //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
        cent_text_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
        if(trigger == 0|| trigger == 2) cent_text_2->AddText("Minimum bias");
        if(trigger == 1) cent_text_2->AddText("JP2");
        cent_text_2->AddText("#Lambda^{0}");
        cent_text_2->AddText(eta_range->Data());
        cent_text_2->AddText(pT_range->Data());
        cent_text_2->SetFillColorAlpha(0, 0.01);
        cent_text_2->Draw("same");

        InvMassCanL0_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0_inv_mass_US_LS_pT_%i_eta_%i.png", pTbin, etaBin));
        //InvMassCanL0_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0_inv_mass_US_LS_pT_%i_eta_%i.pdf", pTbin, etaBin));


        InvMassCanL0[pTbin][etaBin]->cd();

        gPad->SetLeftMargin(0.1);

        L0_inv_mass[pTbin][etaBin] = (TH1D*)L0_inv_mass_US[pTbin][etaBin]->Clone();
        L0_inv_mass[pTbin][etaBin]->Sumw2();
        L0_inv_mass[pTbin][etaBin]->Add(L0_inv_mass_LS[pTbin][etaBin], -1);
        L0_inv_mass[pTbin][etaBin]->GetXaxis()->SetTitle("M_{inv}^{p#pi} (GeV/#it{c}^{2})");
        L0_inv_mass[pTbin][etaBin]->GetXaxis()->CenterTitle();
        L0_inv_mass[pTbin][etaBin]->GetXaxis()->SetTitleSize(0.04);
        L0_inv_mass[pTbin][etaBin]->GetXaxis()->SetLabelSize(0.04);
        L0_inv_mass[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
        L0_inv_mass[pTbin][etaBin]->GetYaxis()->CenterTitle();
        L0_inv_mass[pTbin][etaBin]->GetYaxis()->SetTitleSize(0.04);
        L0_inv_mass[pTbin][etaBin]->GetYaxis()->SetLabelSize(0.04);
        L0_inv_mass[pTbin][etaBin]->SetMarkerStyle(20);
        L0_inv_mass[pTbin][etaBin]->SetMarkerColor(kRed);
        L0_inv_mass[pTbin][etaBin]->SetLineColor(kRed);
        L0_inv_mass[pTbin][etaBin]->Draw("p e");


        if( ( (energy == 510 ) && pTbin > 0 ) || energy == 200 )//skip first pT bin for 510 GeV - no peak visisble
        {
          //if( pTbin == 1 ) fitGauss_L0->SetParameters(2000, 1.116, 0.002);
          if( ( ( energy == 510 ) && pTbin == 1) ) fitGauss_L0->SetParameters(2000, 1.116, 0.002);
          else if( energy == 200) fitGauss_L0->SetParameters(50, 1.116, 0.002);
          else fitGauss_L0->SetParameters(fit_res_gaus_L0->Parameter(0), fit_res_gaus_L0->Parameter(1), fit_res_gaus_L0->Parameter(2));


          fit_res_gaus_L0 = L0_inv_mass[pTbin][etaBin]->Fit(fitGauss_L0, "s i l 0", "", 1.11, 1.12);

          if( /*!fit_res_gaus_L0 ||*/ !fit_res_gaus_L0->IsValid())
          {
            cout<<"Fit not valid for L0 bin pT = "<<pTbin<<" eta = "<<etaBin<<endl;
            return false;
          }

          fitGauss_L0->SetLineColor(1);
          fitGauss_L0->Draw("same");

          float peak_range_min = fit_res_gaus_L0->Parameter(1) - 3*fit_res_gaus_L0->Parameter(2);
          float peak_range_max = fit_res_gaus_L0->Parameter(1) + 3*fit_res_gaus_L0->Parameter(2);

          invMassRange_L0[0][pTbin][etaBin] = peak_range_min;
          invMassRange_L0[1][pTbin][etaBin] = peak_range_max;

          int peak_range_min_bin = L0_inv_mass[pTbin][etaBin]->GetXaxis()->FindBin(peak_range_min);
          int peak_range_max_bin = L0_inv_mass[pTbin][etaBin]->GetXaxis()->FindBin(peak_range_max);

          L0_kappa[pTbin][etaBin] = L0_inv_mass[pTbin][etaBin]->Integral(peak_range_min_bin, peak_range_max_bin)/L0_inv_mass_US[pTbin][etaBin]->Integral(peak_range_min_bin, peak_range_max_bin);


        }
        else
        {
          invMassRange_L0[0][pTbin][etaBin] = 0; //min for pT bin == 0
          invMassRange_L0[1][pTbin][etaBin] = 0; //max for pT bin == 0
        }


        zeroLine->Draw("same");

        SigMinusBckgLeg->Draw("same");

        cent_text_2->Draw("same");

        InvMassCanL0[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0_inv_mass_signal_pT_%i_eta_%i.png", pTbin, etaBin));
        //InvMassCanL0[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0_inv_mass_signal_pT_%i_eta_%i.pdf", pTbin, etaBin));

        //------------------------------------------------------------------------------------------------------------





        gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);


        InvMassCanL0bar_US_LS[pTbin][etaBin]->cd();

        gPad->SetLeftMargin(0.1);

        L0bar_inv_mass_US[pTbin][etaBin]->SetMarkerStyle(20);
        L0bar_inv_mass_US[pTbin][etaBin]->SetMarkerColor(kRed);
        L0bar_inv_mass_US[pTbin][etaBin]->SetLineColor(kRed);
        L0bar_inv_mass_US[pTbin][etaBin]->GetXaxis()->SetTitle("M_{inv}^{#bar{p}#pi^{+}} (GeV/#it{c}^{2})");
        L0bar_inv_mass_US[pTbin][etaBin]->GetXaxis()->CenterTitle();
        L0bar_inv_mass_US[pTbin][etaBin]->GetXaxis()->SetTitleSize(0.04);
        L0bar_inv_mass_US[pTbin][etaBin]->GetXaxis()->SetLabelSize(0.04);
        L0bar_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
        L0bar_inv_mass_US[pTbin][etaBin]->GetYaxis()->CenterTitle();
        L0bar_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetTitleSize(0.04);
        L0bar_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetLabelSize(0.04);
        //L0bar_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 7000);
        //L0bar_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 3000);
        L0bar_inv_mass_US[pTbin][etaBin]->Draw("p e");

        L0bar_inv_mass_LS[pTbin][etaBin]->SetMarkerStyle(24);
        L0bar_inv_mass_LS[pTbin][etaBin]->SetMarkerColor(kBlue);
        L0bar_inv_mass_LS[pTbin][etaBin]->SetLineColor(kBlue);
        L0bar_inv_mass_LS[pTbin][etaBin]->Draw("p e same");


        SigAndBckgLeg->Draw("same");

        TPaveText *cent_text_3 = new TPaveText(0.15, 0.5, 0.45, 0.8, "NDC");
        cent_text_3->SetTextFont(42);
        //cent_text->AddText("STAR preliminary");
        //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
        cent_text_3->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
        if(trigger == 0 || trigger == 2) cent_text_3->AddText("Minimum bias");
        if(trigger == 1) cent_text_3->AddText("JP2");
        cent_text_3->AddText("#bar{#Lambda^{0}}");
        cent_text_3->AddText(eta_range->Data());
        cent_text_3->AddText(pT_range->Data());
        cent_text_3->SetFillColorAlpha(0, 0.01);
        cent_text_3->Draw("same");

        InvMassCanL0bar_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0bar_inv_mass_US_LS_pT_%i_eta_%i.png", pTbin, etaBin));
        //InvMassCanL0bar_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0bar_inv_mass_US_LS_pT_%i_eta_%i.pdf", pTbin, etaBin));


        InvMassCanL0bar[pTbin][etaBin]->cd();

        gPad->SetLeftMargin(0.1);

        L0bar_inv_mass[pTbin][etaBin] = (TH1D*)L0bar_inv_mass_US[pTbin][etaBin]->Clone();
        L0bar_inv_mass[pTbin][etaBin]->Sumw2();
        L0bar_inv_mass[pTbin][etaBin]->Add(L0bar_inv_mass_LS[pTbin][etaBin], -1);
        L0bar_inv_mass[pTbin][etaBin]->GetXaxis()->SetTitle("M_{inv}^{#bar{p}#pi^{+}} (GeV/#it{c}^{2})");
        L0bar_inv_mass[pTbin][etaBin]->GetXaxis()->CenterTitle();
        L0bar_inv_mass[pTbin][etaBin]->GetXaxis()->SetTitleSize(0.04);
        L0bar_inv_mass[pTbin][etaBin]->GetXaxis()->SetLabelSize(0.04);
        L0bar_inv_mass[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
        L0bar_inv_mass[pTbin][etaBin]->GetYaxis()->CenterTitle();
        L0bar_inv_mass[pTbin][etaBin]->GetYaxis()->SetTitleSize(0.04);
        L0bar_inv_mass[pTbin][etaBin]->GetYaxis()->SetLabelSize(0.04);
        L0bar_inv_mass[pTbin][etaBin]->SetMarkerStyle(20);
        L0bar_inv_mass[pTbin][etaBin]->SetMarkerColor(kRed);
        L0bar_inv_mass[pTbin][etaBin]->SetLineColor(kRed);
        L0bar_inv_mass[pTbin][etaBin]->Draw("p e");


        if( ( (energy == 510 ) && pTbin > 0 ) || energy == 200 )//skip first pT bin for 510 GeV - no peak visisble
        {
          //if(   pTbin == 1 ) fitGauss_L0bar->SetParameters(2000, 1.116, 0.002);
          if( ( (energy == 510 ) && pTbin == 1) ) fitGauss_L0bar->SetParameters(2000, 1.116, 0.002);
          else if( energy == 200 ) fitGauss_L0bar->SetParameters(50, 1.116, 0.002);
          else fitGauss_L0bar->SetParameters(fit_res_gaus_L0bar->Parameter(0), fit_res_gaus_L0bar->Parameter(1), fit_res_gaus_L0bar->Parameter(2));


          fit_res_gaus_L0bar = L0bar_inv_mass[pTbin][etaBin]->Fit(fitGauss_L0bar, "s i l 0", "", 1.11, 1.12);

          if(/*!fit_res_gaus_L0bar ||*/ !fit_res_gaus_L0bar->IsValid())
          {
            cout<<"Fit not valid for L0bar bin pT = "<<pTbin<<" eta = "<<etaBin<<endl;
            return false;
          }

          fitGauss_L0bar->SetLineColor(1);
          fitGauss_L0bar->Draw("same");

          float peak_range_min = fit_res_gaus_L0bar->Parameter(1) - 3*fit_res_gaus_L0bar->Parameter(2);
          float peak_range_max = fit_res_gaus_L0bar->Parameter(1) + 3*fit_res_gaus_L0bar->Parameter(2);

          invMassRange_L0bar[0][pTbin][etaBin] = peak_range_min;
          invMassRange_L0bar[1][pTbin][etaBin] = peak_range_max;

          int peak_range_min_bin = L0bar_inv_mass[pTbin][etaBin]->GetXaxis()->FindBin(peak_range_min);
          int peak_range_max_bin = L0bar_inv_mass[pTbin][etaBin]->GetXaxis()->FindBin(peak_range_max);

          L0bar_kappa[pTbin][etaBin] = L0bar_inv_mass[pTbin][etaBin]->Integral(peak_range_min_bin, peak_range_max_bin)/L0bar_inv_mass_US[pTbin][etaBin]->Integral(peak_range_min_bin, peak_range_max_bin);

        }
        else
        {
          invMassRange_L0bar[0][pTbin][etaBin] = 0; //min for pT bin == 0
          invMassRange_L0bar[1][pTbin][etaBin] = 0; //max for pT bin == 0
        }

        zeroLine->Draw("same");

        SigMinusBckgLeg->Draw("same");

        cent_text_3->Draw("same");

        InvMassCanL0bar[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0bar_inv_mass_signal_pT_%i_eta_%i.png", pTbin, etaBin));
        //InvMassCanL0bar[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0bar_inv_mass_signal_pT_%i_eta_%i.pdf", pTbin, etaBin));


        //save histograms to multi-page PDF
        if(pTbin == 0 && etaBin == 0)
        {
          InvMassCan_US_LS[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_inv_mass_US_LS.pdf(", "pdf");
          InvMassCanL0_US_LS[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/InvMassCanL0_US_LS.pdf(", "pdf");
          InvMassCanL0bar_US_LS[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/InvMassCanL0bar_US_LS.pdf(", "pdf");

          InvMassCan[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_inv_mass.pdf(", "pdf");
          InvMassCanL0[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/InvMassCanL0.pdf(", "pdf");
          InvMassCanL0bar[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/InvMassCanL0bar.pdf(", "pdf");

        }
        else if( pTbin == nPtBins && etaBin == nEtaBins)
        {
          InvMassCan_US_LS[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_inv_mass_US_LS.pdf)", "pdf");
          InvMassCanL0_US_LS[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/InvMassCanL0_US_LS.pdf)", "pdf");
          InvMassCanL0bar_US_LS[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/InvMassCanL0bar_US_LS.pdf)", "pdf");

          InvMassCan[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_inv_mass.pdf)", "pdf");
          InvMassCanL0[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/InvMassCanL0.pdf)", "pdf");
          InvMassCanL0bar[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/InvMassCanL0bar.pdf)", "pdf");
        }
        else
        {
          InvMassCan_US_LS[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_inv_mass_US_LS.pdf", "pdf");
          InvMassCanL0_US_LS[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/InvMassCanL0_US_LS.pdf", "pdf");
          InvMassCanL0bar_US_LS[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/InvMassCanL0bar_US_LS.pdf", "pdf");

          InvMassCan[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_inv_mass.pdf", "pdf");
          InvMassCanL0[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/InvMassCanL0.pdf", "pdf");
          InvMassCanL0bar[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/InvMassCanL0bar.pdf", "pdf");
        }

      }
    }

    //------------------------------------------------------------------------------------------------------------

    if(ReadMode !=0 ) InvMassFile->Close();

    InvMassFile->Close();

    return true;

}
//_______________________________________________________________________________________________________________________________________________________________________________________


//for analysis of Lambda polarization
void LambdaLambdaBarSpinCorr(TChain *L_tree, double invMassRange_L0[2][nPtBins+1][nEtaBins+1], double invMassRange_L0bar[2][nPtBins+1][nEtaBins+1], double L0_kappa[nPtBins+1][nEtaBins+1], double L0bar_kappa[nPtBins+1][nEtaBins+1], const int ReadMode, const int trigger = 0,const int energy = 510, const int year = 2017)
{
  const float L0_alpha = 0.732; //decay parameter of L0
  const float L0bar_alpha = -0.758; //decay paramteter of L0bar

  //set TOF requirement based on trigger
  int requireTOF_Lcorr = 0; //dafault is no TOF requirement

  if(trigger == 0 || trigger == 1) requireTOF_Lcorr = 1;

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
      return;
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
    return;
  }

  TFile *EffFileEmbedd;

  if(energy == 510) EffFileEmbedd = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run12_embedd.root", "read");
  else if(energy == 200) EffFileEmbedd = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run12_embedd.root", "read");
  //if(energy == 510) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run17_1B_ME_tight_eta.root", "read");
  //else if(energy == 200) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run12_1B_ME_tight_eta.root", "read");
  else
  {
    cout<<"Not a valid colliison energy! Abborting!"<<endl;
    return;
  }

  //_______________________________________________________________________________________________________________________________________________

  //efficiency histograms

  //TH1D *L0_L0bar_cosThetaProdPlane_eff = (TH1D*)EffFileEmbedd->Get("L0_L0bar_cosThetaProdPlane_RC_hist");
  TH1D *L0_L0bar_cosThetaProdPlane_eff = (TH1D*)EffFile->Get("L0_L0bar_cosThetaProdPlane_eff");
  TH1D *L0_L0_cosThetaProdPlane_eff = (TH1D*)EffFile->Get("L0_L0_cosThetaProdPlane_eff");
  TH1D *L0bar_L0bar_cosThetaProdPlane_eff = (TH1D*)EffFile->Get("L0bar_L0bar_cosThetaProdPlane_eff");

  TH1D *L0_L0bar_cosThetaProdPlane_pT_eff[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0bar_cosThetaProdPlane_eta_eff[nEtaBins][nEtaBins];

  TH1D *L0_L0_cosThetaProdPlane_pT_eff[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0_cosThetaProdPlane_eta_eff[nEtaBins][nEtaBins];

  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_eff[nPtBins_corr][nPtBins_corr];
  TH1D *L0bar_L0bar_cosThetaProdPlane_eta_eff[nEtaBins][nEtaBins];


  TH1D *L0_L0bar_cosThetaProdPlane_ME_eff = (TH1D*)EffFile->Get("L0_L0bar_cosThetaProdPlane_ME_eff");
  TH1D *L0_L0_cosThetaProdPlane_ME_eff = (TH1D*)EffFile->Get("L0_L0_cosThetaProdPlane_ME_eff");
  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_eff = (TH1D*)EffFile->Get("L0bar_L0bar_cosThetaProdPlane_ME_eff");

  TH1D *L0_L0bar_cosThetaProdPlane_ME_pT_eff[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0bar_cosThetaProdPlane_ME_eta_eff[nEtaBins][nEtaBins];

  TH1D *L0_L0_cosThetaProdPlane_ME_pT_eff[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0_cosThetaProdPlane_ME_eta_eff[nEtaBins][nEtaBins];

  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_pT_eff[nPtBins_corr][nPtBins_corr];
  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_eta_eff[nEtaBins][nEtaBins];

  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      //L0_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2] = (TH1D*)EffFileEmbedd->Get(Form("L0_L0bar_cosThetaProdPlane_RC_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      L0_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2] = (TH1D*)EffFile->Get(Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0_L0_cosThetaProdPlane_pT_eff[pTbin1][pTbin2] = (TH1D*)EffFile->Get(Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0bar_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2] = (TH1D*)EffFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));

      L0_L0bar_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2] = (TH1D*)EffFile->Get(Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0_L0_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2] = (TH1D*)EffFile->Get(Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
      L0bar_L0bar_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2] = (TH1D*)EffFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_eff", pTbin1, pTbin2));
    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      L0_L0bar_cosThetaProdPlane_eta_eff[etaBin1][etaBin2] = (TH1D*)EffFile->Get(Form("L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0_L0_cosThetaProdPlane_eta_eff[etaBin1][etaBin2] = (TH1D*)EffFile->Get(Form("L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0bar_L0bar_cosThetaProdPlane_eta_eff[etaBin1][etaBin2] = (TH1D*)EffFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));

      L0_L0bar_cosThetaProdPlane_ME_eta_eff[etaBin1][etaBin2] = (TH1D*)EffFile->Get(Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0_L0_cosThetaProdPlane_ME_eta_eff[etaBin1][etaBin2] = (TH1D*)EffFile->Get(Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
      L0bar_L0bar_cosThetaProdPlane_ME_eta_eff[etaBin1][etaBin2] = (TH1D*)EffFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_eff", etaBin1, etaBin2));
    }
  }
  //_____________________________________________________________________________

  //stat. err. projection histograms

  TH1D *L0_L0bar_cosThetaProdPlane_stat_err = new TH1D("L0_L0bar_cosThetaProdPlane_stat_err", "L0_L0bar_cosThetaProdPlane_stat_err", 10, -1, 1);
  TH1D *L0_L0_cosThetaProdPlane_stat_err = new TH1D("L0_L0_cosThetaProdPlane_stat_err",  "L0_L0_cosThetaProdPlane_stat_err",10, -1, 1);
  TH1D *L0bar_L0bar_cosThetaProdPlane_stat_err = new TH1D("L0bar_L0bar_cosThetaProdPlane_stat_err", "L0bar_L0bar_cosThetaProdPlane_stat_err", 10, -1, 1);

  TH1D *L0_L0bar_cosThetaProdPlane_pT_stat_err[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0bar_cosThetaProdPlane_eta_stat_err[nEtaBins][nEtaBins];

  TH1D *L0_L0_cosThetaProdPlane_pT_stat_err[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0_cosThetaProdPlane_eta_stat_err[nEtaBins][nEtaBins];

  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_stat_err[nPtBins_corr][nPtBins_corr];
  TH1D *L0bar_L0bar_cosThetaProdPlane_eta_stat_err[nEtaBins][nEtaBins];


  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {

      L0_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_stat_err", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_stat_err", pTbin1, pTbin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_stat_err", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_stat_err", pTbin1, pTbin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_stat_err", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_stat_err", pTbin1, pTbin2), 10, -1, 1);

    }
  }

  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      L0_L0bar_cosThetaProdPlane_eta_stat_err[etaBin1][etaBin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_stat_err", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_stat_err", etaBin1, etaBin2), 10, -1, 1);
      L0_L0_cosThetaProdPlane_eta_stat_err[etaBin1][etaBin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_stat_err", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i_stat_err", etaBin1, etaBin2), 10, -1, 1);
      L0bar_L0bar_cosThetaProdPlane_eta_stat_err[etaBin1][etaBin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_stat_err", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i_stat_err", etaBin1, etaBin2), 10, -1, 1);

    }
  }
  //_____________________________________________________________________________

  //data histograms
  TH1D *L0_L0bar_cosThetaProdPlane_US_hist;
  TH1D *L0_L0bar_cosThetaProdPlane_LS_hist;
  TH1D *L0_L0bar_cosThetaProdPlane_ME_hist; //mixed event

  TH1D *L0_L0bar_cosThetaProdPlane_pT_US_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0bar_cosThetaProdPlane_pT_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0bar_cosThetaProdPlane_pT_ME_hist[nPtBins_corr][nPtBins_corr];

  //|eta| < 0.2
  TH1D *L0_L0bar_cosThetaProdPlane_pT_tight_eta_US_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0bar_cosThetaProdPlane_pT_tight_eta_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0bar_cosThetaProdPlane_pT_tight_eta_ME_hist[nPtBins_corr][nPtBins_corr];

  //low pT is fixed, upper is different for each bin
  TH1D *L0_L0bar_cosThetaProdPlane_pUp_US_hist[nPtBins_corr][nPtBins-1];
  TH1D *L0_L0bar_cosThetaProdPlane_pUp_LS_hist[nPtBins_corr][nPtBins-1];
  TH1D *L0_L0bar_cosThetaProdPlane_pUp_ME_hist[nPtBins_corr][nPtBins-1];

  TH1D *L0_L0bar_cosThetaProdPlane_eta_US_hist[nEtaBins][nEtaBins];
  TH1D *L0_L0bar_cosThetaProdPlane_eta_LS_hist[nEtaBins][nEtaBins];
  TH1D *L0_L0bar_cosThetaProdPlane_eta_ME_hist[nEtaBins][nEtaBins];

  TH2D *L0_L0bar_delta_eta_vs_delta_phi_US_hist;
  TH2D *L0_L0bar_delta_eta_vs_delta_phi_LS_hist;

  TH2D *L0_L0bar_y1_vs_y2_US_hist;
  TH2D *L0_L0bar_y1_vs_y2_LS_hist;

  TH2D *L0_L0bar_pT1_vs_pT2_US_hist;
  TH2D *L0_L0bar_pT1_vs_pT2_LS_hist;

  TH2D *L0_L0bar_phi1_vs_phi2_US_hist;
  TH2D *L0_L0bar_phi1_vs_phi2_LS_hist;


  TH1D *L0_L0_cosThetaProdPlane_US_hist;
  TH1D *L0_L0_cosThetaProdPlane_LS_hist;
  TH1D *L0_L0_cosThetaProdPlane_ME_hist;

  TH1D *L0_L0_cosThetaProdPlane_pT_US_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0_cosThetaProdPlane_pT_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0_cosThetaProdPlane_pT_ME_hist[nPtBins_corr][nPtBins_corr];

  //|eta| < 0.2
  TH1D *L0_L0_cosThetaProdPlane_pT_tight_eta_US_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0_cosThetaProdPlane_pT_tight_eta_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0_L0_cosThetaProdPlane_pT_tight_eta_ME_hist[nPtBins_corr][nPtBins_corr];

  //low pT is fixed, upper is different for each bin
  TH1D *L0_L0_cosThetaProdPlane_pUp_US_hist[nPtBins_corr][nPtBins-1];
  TH1D *L0_L0_cosThetaProdPlane_pUp_LS_hist[nPtBins_corr][nPtBins-1];
  TH1D *L0_L0_cosThetaProdPlane_pUp_ME_hist[nPtBins_corr][nPtBins-1];

  TH1D *L0_L0_cosThetaProdPlane_eta_US_hist[nEtaBins][nEtaBins];
  TH1D *L0_L0_cosThetaProdPlane_eta_LS_hist[nEtaBins][nEtaBins];
  TH1D *L0_L0_cosThetaProdPlane_eta_ME_hist[nEtaBins][nEtaBins];

  TH2D *L0_L0_delta_eta_vs_delta_phi_US_hist;
  TH2D *L0_L0_delta_eta_vs_delta_phi_LS_hist;

  TH2D *L0_L0_delta_eta_vs_delta_phi_US_zoom_hist;
  TH2D *L0_L0_delta_eta_vs_delta_phi_LS_zoom_hist;

  TH2D *L0_L0_y1_vs_y2_US_hist;
  TH2D *L0_L0_y1_vs_y2_LS_hist;

  TH2D *L0_p_L0_p_y1_vs_y2_US_hist;
  TH2D *L0_p_L0_p_y1_vs_y2_LS_hist;

  TH2D *L0_pi_L0_pi_y1_vs_y2_US_hist;
  TH2D *L0_pi_L0_pi_y1_vs_y2_LS_hist;

  TH2D *L0_L0_pT1_vs_pT2_US_hist;
  TH2D *L0_L0_pT1_vs_pT2_LS_hist;

  TH2D *L0_L0_phi1_vs_phi2_US_hist;
  TH2D *L0_L0_phi1_vs_phi2_LS_hist;

  TH2D *L0_p_L0_p_phi1_vs_phi2_US_hist;
  TH2D *L0_p_L0_p_phi1_vs_phi2_LS_hist;

  TH2D *L0_pi_L0_pi_phi1_vs_phi2_US_hist;
  TH2D *L0_pi_L0_pi_phi1_vs_phi2_LS_hist;


  TH1D *L0bar_L0bar_cosThetaProdPlane_US_hist;
  TH1D *L0bar_L0bar_cosThetaProdPlane_LS_hist;
  TH1D *L0bar_L0bar_cosThetaProdPlane_ME_hist;

  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_US_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[nPtBins_corr][nPtBins_corr];

  //|eta| < 0.2
  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_tight_eta_US_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_tight_eta_LS_hist[nPtBins_corr][nPtBins_corr];
  TH1D *L0bar_L0bar_cosThetaProdPlane_pT_tight_eta_ME_hist[nPtBins_corr][nPtBins_corr];

  //low pT is fixed, upper is different for each bin
  TH1D *L0bar_L0bar_cosThetaProdPlane_pUp_US_hist[nPtBins_corr][nPtBins-1];
  TH1D *L0bar_L0bar_cosThetaProdPlane_pUp_LS_hist[nPtBins_corr][nPtBins-1];
  TH1D *L0bar_L0bar_cosThetaProdPlane_pUp_ME_hist[nPtBins_corr][nPtBins-1];

  TH1D *L0bar_L0bar_cosThetaProdPlane_eta_US_hist[nEtaBins][nEtaBins];
  TH1D *L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[nEtaBins][nEtaBins];
  TH1D *L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[nEtaBins][nEtaBins];

  TH2D *L0bar_L0bar_delta_eta_vs_delta_phi_US_hist;
  TH2D *L0bar_L0bar_delta_eta_vs_delta_phi_LS_hist;

  TH2D *L0bar_L0bar_delta_eta_vs_delta_phi_US_zoom_hist;
  TH2D *L0bar_L0bar_delta_eta_vs_delta_phi_LS_zoom_hist;

  TH2D *L0bar_L0bar_y1_vs_y2_US_hist;
  TH2D *L0bar_L0bar_y1_vs_y2_LS_hist;

  TH2D *L0bar_p_L0bar_p_y1_vs_y2_US_hist;
  TH2D *L0bar_p_L0bar_p_y1_vs_y2_LS_hist;

  TH2D *L0bar_pi_L0bar_pi_y1_vs_y2_US_hist;
  TH2D *L0bar_pi_L0bar_pi_y1_vs_y2_LS_hist;

  TH2D *L0bar_L0bar_pT1_vs_pT2_US_hist;
  TH2D *L0bar_L0bar_pT1_vs_pT2_LS_hist;

  TH2D *L0bar_L0bar_phi1_vs_phi2_US_hist;
  TH2D *L0bar_L0bar_phi1_vs_phi2_LS_hist;

  TH2D *L0bar_p_L0bar_p_phi1_vs_phi2_US_hist;
  TH2D *L0bar_p_L0bar_p_phi1_vs_phi2_LS_hist;

  TH2D *L0bar_pi_L0bar_pi_phi1_vs_phi2_US_hist;
  TH2D *L0bar_pi_L0bar_pi_phi1_vs_phi2_LS_hist;


  TH1D *L_pair_cosThetaProdPlane_US_hist;
  TH1D *L_pair_cosThetaProdPlane_LS_hist;


  //polarization histograms

  //one bin for each lower pT selection
  //histograms are as function of upper pT
  TH1D *L0_L0bar_polarization_pT_hist[nPtBins_corr];
  TH1D *L0_L0_polarization_pT_hist[nPtBins_corr];
  TH1D *L0bar_L0bar_polarization_pT_hist[nPtBins_corr];


  //invariant mass histigrams for L from pairs
  TH1D *L_inv_mass_US_from_L_Lbar = new TH1D("L_inv_mass_US_from_L_Lbar", "L_inv_mass_US_from_L_Lbar", 200, 1, 1.2);
  TH1D *L_inv_mass_LS_from_L_Lbar = new TH1D("L_inv_mass_LS_from_L_Lbar", "L_inv_mass_LS_from_L_Lbar", 200, 1, 1.2);

  TH1D *Lbar_inv_mass_US_from_L_Lbar = new TH1D("Lbar_inv_mass_US_from_L_Lbar", "Lbar_inv_mass_US_from_L_Lbar", 200, 1, 1.2);
  TH1D *Lbar_inv_mass_LS_from_L_Lbar = new TH1D("Lbar_inv_mass_LS_from_L_Lbar", "Lbar_inv_mass_LS_from_L_Lbar", 200, 1, 1.2);

  TH2D *L_vs_Lbar_mass_US = new TH2D("L_vs_Lbar_mass_US", "L_vs_Lbar_mass_US", 200, 1, 1.2, 200, 1, 1.2);
  TH2D *L_vs_Lbar_mass_LS = new TH2D("L_vs_Lbar_mass_LS", "L_vs_Lbar_mass_LS", 200, 1, 1.2, 200, 1, 1.2);

  TH1D *L_inv_mass_US_from_L_L = new TH1D("L_inv_mass_US_from_L_L", "L_inv_mass_US_from_L_L", 200, 1, 1.2);
  TH1D *L_inv_mass_LS_from_L_L = new TH1D("L_inv_mass_LS_from_L_L", "L_inv_mass_LS_from_L_L", 200, 1, 1.2);

  TH2D *L_vs_L_mass_US = new TH2D("L_vs_L_mass_US", "L_vs_L_mass_US", 200, 1, 1.2, 200, 1, 1.2);
  TH2D *L_vs_L_mass_LS = new TH2D("L_vs_L_mass_LS", "L_vs_L_mass_LS", 200, 1, 1.2, 200, 1, 1.2);

  TH1D *Lbar_inv_mass_US_from_Lbar_Lbar = new TH1D("Lbar_inv_mass_US_from_Lbar_Lbar", "Lbar_inv_mass_US_from_Lbar_Lbar", 200, 1, 1.2);
  TH1D *Lbar_inv_mass_LS_from_Lbar_Lbar = new TH1D("Lbar_inv_mass_LS_from_Lbar_Lbar", "Lbar_inv_mass_LS_from_Lbar_Lbar", 200, 1, 1.2);

  TH2D *Lbar_vs_Lbar_mass_US = new TH2D("Lbar_vs_Lbar_mass_US", "Lbar_vs_Lbar_mass_US", 200, 1, 1.2, 200, 1, 1.2);
  TH2D *Lbar_vs_Lbar_mass_LS = new TH2D("Lbar_vs_Lbar_mass_LS", "Lbar_vs_Lbar_mass_LS", 200, 1, 1.2, 200, 1, 1.2);



  //________________________________________________________________

  //efficiency histograms
  //need to produce efficiency histograms for L-Lbar correlations

  //______________________________________________________________


  //TH1D *nLambdaPairsInEvent_US_hist = new TH1D("nLambdaPairsInEvent_US_hist", "nLambdaPairsInEvent_US_hist", 10, 0, 10);
  TH1D *nLambdasInEvent_US_hist = new TH1D("nLambdasInEvent_US_hist", "nLambdasInEvent_US_hist", 10, 0, 10);
  TH1D *nLambdaBarsInEvent_US_hist = new TH1D("nLambdaBarsInEvent_US_hist", "nLambdaBarsInEvent_US_hist", 10, 0, 10);

  TH1D *nEventsWithLambdaPair_US_hist = new TH1D("nEventsWithLambdaPair_US_hist", "nEventsWithLambdaPair_US_hist", 2, 0, 2);

  //TH1D *nLambdaPairsInEvent_LS_hist = new TH1D("nLambdaPairsInEvent_LS_hist", "nLambdaPairsInEvent_LS_hist", 10, 0, 10);
  TH1D *nLambdasInEvent_LS_hist = new TH1D("nLambdasInEvent_LS_hist", "nLambdasInEvent_LS_hist", 10, 0, 10);
  TH1D *nLambdaBarsInEvent_LS_hist = new TH1D("nLambdaBarsInEvent_LS_hist", "nLambdaBarsInEvent_LS_hist", 10, 0, 10);

  TH1D *nEventsWithLambdaPair_LS_hist = new TH1D("nEventsWithLambdaPair_LS_hist", "nEventsWithLambdaPair_LS_hist", 2, 0, 2);


  if(ReadMode == 0) //create histograms to be saved into file
  {

    for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
    {
      for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
      {
        L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
        L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
        L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

        L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
        L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
        L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);

        L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
        L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
        L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2), 10, -1, 1);
      }
    }

    for(unsigned int pTupbin = 0; pTupbin < nPtBins_corr; pTupbin++)
    {
      for(unsigned int pTbin = 0; pTbin < nPtBins-1; pTbin++)
      {
        L0_L0bar_cosThetaProdPlane_pUp_US_hist[pTupbin][pTbin] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_pUp_US_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), Form("L0_L0bar_cosThetaProdPlane_pUp_US_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), 10, -1, 1);
        L0_L0bar_cosThetaProdPlane_pUp_LS_hist[pTupbin][pTbin] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_pUp_LS_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), Form("L0_L0bar_cosThetaProdPlane_pUp_LS_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), 10, -1, 1);
        L0_L0bar_cosThetaProdPlane_pUp_ME_hist[pTupbin][pTbin] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_pUp_ME_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), Form("L0_L0bar_cosThetaProdPlane_pUp_ME_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), 10, -1, 1);

        L0_L0_cosThetaProdPlane_pUp_US_hist[pTupbin][pTbin] = new TH1D(Form("L0_L0_cosThetaProdPlane_pUp_US_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), Form("L0_L0_cosThetaProdPlane_pUp_US_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), 10, -1, 1);
        L0_L0_cosThetaProdPlane_pUp_LS_hist[pTupbin][pTbin] = new TH1D(Form("L0_L0_cosThetaProdPlane_pUp_LS_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), Form("L0_L0_cosThetaProdPlane_pUp_LS_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), 10, -1, 1);
        L0_L0_cosThetaProdPlane_pUp_ME_hist[pTupbin][pTbin] = new TH1D(Form("L0_L0_cosThetaProdPlane_pUp_ME_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), Form("L0_L0_cosThetaProdPlane_pUp_ME_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), 10, -1, 1);

        L0bar_L0bar_cosThetaProdPlane_pUp_US_hist[pTupbin][pTbin] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_pUp_US_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), Form("L0bar_L0bar_cosThetaProdPlane_pUp_US_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), 10, -1, 1);
        L0bar_L0bar_cosThetaProdPlane_pUp_LS_hist[pTupbin][pTbin] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_pUp_LS_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), Form("L0bar_L0bar_cosThetaProdPlane_pUp_LS_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), 10, -1, 1);
        L0bar_L0bar_cosThetaProdPlane_pUp_ME_hist[pTupbin][pTbin] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_pUp_ME_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), Form("L0bar_L0bar_cosThetaProdPlane_pUp_ME_pUp1_%i_pUp2_%i_hist", pTupbin, pTbin), 10, -1, 1);
      }

      //one bin for each lower pT selection
      //histograms are as function of upper pT
      L0_L0bar_polarization_pT_hist[pTupbin] = new TH1D(Form("L0_L0bar_polarization_pT_pTup_%i", pTupbin), Form("L0_L0bar_polarization_pT_pTup_%i", pTupbin), nPtBins, pT_bins);
      L0_L0_polarization_pT_hist[pTupbin] = new TH1D(Form("L0_L0_polarization_pT_pTup_%i", pTupbin), Form("L0_L0_polarization_pT_pTup_%i", pTupbin), nPtBins, pT_bins);
      L0bar_L0bar_polarization_pT_hist[pTupbin] = new TH1D(Form("L0bar_L0bar_polarization_pT_pTup_%i", pTupbin), Form("L0bar_L0bar_polarization_pT_pTup_%i", pTupbin), nPtBins, pT_bins);
    }

    for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
    {
      for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
      {
        L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
        L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
        L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

        L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
        L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
        L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = new TH1D(Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);

        L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
        L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
        L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = new TH1D(Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2), 10, -1, 1);
      }
    }

    L0_L0bar_cosThetaProdPlane_US_hist = new TH1D("L0_L0bar_cosThetaProdPlane_US_hist", "L0_L0bar_cosThetaProdPlane_US_hist", 10, -1, 1);
    L0_L0bar_cosThetaProdPlane_LS_hist = new TH1D("L0_L0bar_cosThetaProdPlane_LS_hist", "L0_L0bar_cosThetaProdPlane_LS_hist", 10, -1, 1);
    L0_L0bar_cosThetaProdPlane_ME_hist = new TH1D("L0_L0bar_cosThetaProdPlane_ME_hist", "L0_L0bar_cosThetaProdPlane_ME_hist", 10, -1, 1);

    L0_L0bar_delta_eta_vs_delta_phi_US_hist = new TH2D("L0_L0bar_delta_eta_vs_delta_phi_US_hist", "L0_L0bar_delta_eta_vs_delta_phi_US_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
    L0_L0bar_delta_eta_vs_delta_phi_LS_hist = new TH2D("L0_L0bar_delta_eta_vs_delta_phi_LS_hist", "L0_L0bar_delta_eta_vs_delta_phi_LS_hist", 100, 0, 2, 100, 0, TMath::TwoPi());

    L0_L0bar_y1_vs_y2_US_hist = new TH2D("L0_L0bar_y1_vs_y2_US_hist", "L0_L0bar_y1_vs_y2_US_hist", 100, -1, 1, 100, -1, 1);
    L0_L0bar_y1_vs_y2_LS_hist = new TH2D("L0_L0bar_y1_vs_y2_LS_hist", "L0_L0bar_y1_vs_y2_LS_hist", 100, -1, 1, 100, -1, 1);

    L0_L0bar_pT1_vs_pT2_US_hist = new TH2D("L0_L0bar_pT1_vs_pT2_US_hist", "L0_L0bar_pT1_vs_pT2_US_hist", 50, 0, 5, 50, 0, 5);
    L0_L0bar_pT1_vs_pT2_LS_hist = new TH2D("L0_L0bar_pT1_vs_pT2_LS_hist", "L0_L0bar_pT1_vs_pT2_LS_hist", 50, 0, 5, 50, 0, 5);

    L0_L0bar_phi1_vs_phi2_US_hist = new TH2D("L0_L0bar_phi1_vs_phi2_US_hist", "L0_L0bar_phi1_vs_phi2_US_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
    L0_L0bar_phi1_vs_phi2_LS_hist = new TH2D("L0_L0bar_phi1_vs_phi2_LS_hist", "L0_L0bar_phi1_vs_phi2_LS_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());


    L0_L0_cosThetaProdPlane_US_hist = new TH1D("L0_L0_cosThetaProdPlane_US_hist", "L0_L0_cosThetaProdPlane_US_hist", 10, -1, 1);
    L0_L0_cosThetaProdPlane_LS_hist = new TH1D("L0_L0_cosThetaProdPlane_LS_hist", "L0_L0_cosThetaProdPlane_LS_hist", 10, -1, 1);
    L0_L0_cosThetaProdPlane_ME_hist = new TH1D("L0_L0_cosThetaProdPlane_ME_hist", "L0_L0_cosThetaProdPlane_ME_hist", 10, -1, 1);

    L0_L0_delta_eta_vs_delta_phi_US_hist = new TH2D("L0_L0_delta_eta_vs_delta_phi_US_hist", "L0_L0_delta_eta_vs_delta_phi_US_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
    L0_L0_delta_eta_vs_delta_phi_LS_hist = new TH2D("L0_L0_delta_eta_vs_delta_phi_LS_hist", "L0_L0_delta_eta_vs_delta_phi_LS_hist", 100, 0, 2, 100, 0, TMath::TwoPi());

    L0_L0_delta_eta_vs_delta_phi_US_zoom_hist = new TH2D("L0_L0_delta_eta_vs_delta_phi_US_zoom_hist", "L0_L0_delta_eta_vs_delta_phi_US_zoom_hist", 100, 0, 0.2, 100, 0, 0.2);
    L0_L0_delta_eta_vs_delta_phi_LS_zoom_hist = new TH2D("L0_L0_delta_eta_vs_delta_phi_LS_zoom_hist", "L0_L0_delta_eta_vs_delta_phi_LS_zoom_hist", 100, 0, 0.2, 100, 0, 0.2);

    L0_L0_y1_vs_y2_US_hist = new TH2D("L0_L0_y1_vs_y2_US_hist", "L0_L0_y1_vs_y2_US_hist", 100, -1, 1, 100, -1, 1);
    L0_L0_y1_vs_y2_LS_hist = new TH2D("L0_L0_y1_vs_y2_LS_hist", "L0_L0_y1_vs_y2_LS_hist", 100, -1, 1, 100, -1, 1);

    L0_p_L0_p_y1_vs_y2_US_hist = new TH2D("L0_p_L0_p_y1_vs_y2_US_hist", "L0_p_L0_p_y1_vs_y2_US_hist", 100, -1, 1, 100, -1, 1);
    L0_p_L0_p_y1_vs_y2_LS_hist = new TH2D("L0_p_L0_p_y1_vs_y2_LS_hist", "L0_p_L0_p_y1_vs_y2_LS_hist", 100, -1, 1, 100, -1, 1);

    L0_pi_L0_pi_y1_vs_y2_US_hist = new TH2D("L0_pi_L0_pi_y1_vs_y2_US_hist", "L0_pi_L0_pi_y1_vs_y2_US_hist", 100, -1, 1, 100, -1, 1);
    L0_pi_L0_pi_y1_vs_y2_LS_hist = new TH2D("L0_pi_L0_pi_y1_vs_y2_LS_hist", "L0_pi_L0_pi_y1_vs_y2_LS_hist", 100, -1, 1, 100, -1, 1);

    L0_L0_pT1_vs_pT2_US_hist = new TH2D("L0_L0_pT1_vs_pT2_US_hist", "L0_L0_pT1_vs_pT2_US_hist", 50, 0, 5, 50, 0, 5);
    L0_L0_pT1_vs_pT2_LS_hist = new TH2D("L0_L0_pT1_vs_pT2_LS_hist", "L0_L0_pT1_vs_pT2_LS_hist", 50, 0, 5, 50, 0, 5);

    L0_L0_phi1_vs_phi2_US_hist = new TH2D("L0_L0_phi1_vs_phi2_US_hist", "L0_L0_phi1_vs_phi2_US_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
    L0_L0_phi1_vs_phi2_LS_hist = new TH2D("L0_L0_phi1_vs_phi2_LS_hist", "L0_L0_phi1_vs_phi2_LS_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());

    L0_p_L0_p_phi1_vs_phi2_US_hist = new TH2D("L0_p_L0_p_phi1_vs_phi2_US_hist", "L0_p_L0_p_phi1_vs_phi2_US_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
    L0_p_L0_p_phi1_vs_phi2_LS_hist = new TH2D("L0_p_L0_p_phi1_vs_phi2_LS_hist", "L0_p_L0_p_phi1_vs_phi2_LS_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());

    L0_pi_L0_pi_phi1_vs_phi2_US_hist = new TH2D("L0_pi_L0_pi_phi1_vs_phi2_US_hist", "L0_pi_L0_pi_phi1_vs_phi2_US_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
    L0_pi_L0_pi_phi1_vs_phi2_LS_hist = new TH2D("L0_pi_L0_pi_phi1_vs_phi2_LS_hist", "L0_pi_L0_pi_phi1_vs_phi2_LS_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());


    L0bar_L0bar_cosThetaProdPlane_US_hist = new TH1D("L0bar_L0bar_cosThetaProdPlane_US_hist", "L0bar_L0bar_cosThetaProdPlane_US_hist", 10, -1, 1);
    L0bar_L0bar_cosThetaProdPlane_LS_hist = new TH1D("L0bar_L0bar_cosThetaProdPlane_LS_hist", "L0bar_L0bar_cosThetaProdPlane_LS_hist", 10, -1, 1);
    L0bar_L0bar_cosThetaProdPlane_ME_hist = new TH1D("L0bar_L0bar_cosThetaProdPlane_ME_hist", "L0bar_L0bar_cosThetaProdPlane_ME_hist", 10, -1, 1);

    L0bar_L0bar_delta_eta_vs_delta_phi_US_hist = new TH2D("L0bar_L0bar_delta_eta_vs_delta_phi_US_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_US_hist", 100, 0, 2, 100, 0, TMath::TwoPi());
    L0bar_L0bar_delta_eta_vs_delta_phi_LS_hist = new TH2D("L0bar_L0bar_delta_eta_vs_delta_phi_LS_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_LS_hist", 100, 0, 2, 100, 0, TMath::TwoPi());

    L0bar_L0bar_delta_eta_vs_delta_phi_US_zoom_hist = new TH2D("L0bar_L0bar_delta_eta_vs_delta_phi_US_zoom_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_US_zoom_hist", 100, 0, 0.2, 100, 0, 0.2);
    L0bar_L0bar_delta_eta_vs_delta_phi_LS_zoom_hist = new TH2D("L0bar_L0bar_delta_eta_vs_delta_phi_LS_zoom_hist", "L0bar_L0bar_delta_eta_vs_delta_phi_LS_zoom_hist", 100, 0, 0.2, 100, 0, 0.2);

    L0bar_L0bar_y1_vs_y2_US_hist = new TH2D("L0bar_L0bar_y1_vs_y2_US_hist", "L0bar_L0bar_y1_vs_y2_US_hist", 100, -1, 1, 100, -1, 1);
    L0bar_L0bar_y1_vs_y2_LS_hist = new TH2D("L0bar_L0bar_y1_vs_y2_LS_hist", "L0bar_L0bar_y1_vs_y2_LS_hist", 100, -1, 1, 100, -1, 1);

    L0bar_p_L0bar_p_y1_vs_y2_US_hist = new TH2D("L0bar_p_L0bar_p_y1_vs_y2_US_hist", "L0bar_p_L0bar_p_y1_vs_y2_US_hist", 100, -1, 1, 100, -1, 1);
    L0bar_p_L0bar_p_y1_vs_y2_LS_hist = new TH2D("L0bar_p_L0bar_p_y1_vs_y2_LS_hist", "L0bar_p_L0bar_p_y1_vs_y2_LS_hist", 100, -1, 1, 100, -1, 1);

    L0bar_pi_L0bar_pi_y1_vs_y2_US_hist = new TH2D("L0bar_pi_L0bar_pi_y1_vs_y2_US_hist", "L0bar_pi_L0bar_pi_y1_vs_y2_US_hist", 100, -1, 1, 100, -1, 1);
    L0bar_pi_L0bar_pi_y1_vs_y2_LS_hist = new TH2D("L0bar_pi_L0bar_pi_y1_vs_y2_LS_hist", "L0bar_pi_L0bar_pi_y1_vs_y2_LS_hist", 100, -1, 1, 100, -1, 1);

    L0bar_L0bar_pT1_vs_pT2_US_hist = new TH2D("L0bar_L0bar_pT1_vs_pT2_US_hist", "L0bar_L0bar_pT1_vs_pT2_US_hist", 50, 0, 5, 50, 0, 5);
    L0bar_L0bar_pT1_vs_pT2_LS_hist = new TH2D("L0bar_L0bar_pT1_vs_pT2_LS_hist", "L0bar_L0bar_pT1_vs_pT2_LS_hist", 50, 0, 5, 50, 0, 5);

    L0bar_L0bar_phi1_vs_phi2_US_hist = new TH2D("L0bar_L0bar_phi1_vs_phi2_US_hist", "L0bar_L0bar_phi1_vs_phi2_US_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
    L0bar_L0bar_phi1_vs_phi2_LS_hist = new TH2D("L0bar_L0bar_phi1_vs_phi2_LS_hist", "L0bar_L0bar_phi1_vs_phi2_LS_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());

    L0bar_p_L0bar_p_phi1_vs_phi2_US_hist = new TH2D("L0bar_p_L0bar_p_phi1_vs_phi2_US_hist", "L0bar_p_L0bar_p_phi1_vs_phi2_US_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
    L0bar_p_L0bar_p_phi1_vs_phi2_LS_hist = new TH2D("L0bar_p_L0bar_p_phi1_vs_phi2_LS_hist", "L0bar_p_L0bar_p_phi1_vs_phi2_LS_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());

    L0bar_pi_L0bar_pi_phi1_vs_phi2_US_hist = new TH2D("L0bar_pi_L0bar_pi_phi1_vs_phi2_US_hist", "L0bar_pi_L0bar_pi_phi1_vs_phi2_US_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
    L0bar_pi_L0bar_pi_phi1_vs_phi2_LS_hist = new TH2D("L0bar_pi_L0bar_pi_phi1_vs_phi2_LS_hist", "L0bar_pi_L0bar_pi_phi1_vs_phi2_LS_hist", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());

    L_pair_cosThetaProdPlane_US_hist = new TH1D("L_pair_cosThetaProdPlane_US_hist", "L_pair_cosThetaProdPlane_US_hist", 10, -1, 1);
    L_pair_cosThetaProdPlane_LS_hist = new TH1D("L_pair_cosThetaProdPlane_LS_hist", "L_pair_cosThetaProdPlane_LS_hist", 10, -1, 1);


  }
  else //load histograms from file
  {

    for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
    {
      for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
      {
        L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = (TH1D*)LLbarOutFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
        L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = (TH1D*)LLbarOutFile->Get(Form("L0_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
        L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = (TH1D*)LLbarOutFile->Get(Form("L0_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

        L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = (TH1D*)LLbarOutFile->Get(Form("L0_L0_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
        L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = (TH1D*)LLbarOutFile->Get(Form("L0_L0_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
        L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = (TH1D*)LLbarOutFile->Get(Form("L0_L0_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));

        L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2] = (TH1D*)LLbarOutFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
        L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2] = (TH1D*)LLbarOutFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_LS_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
        L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2] = (TH1D*)LLbarOutFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_ME_pT1_%i_pT2_%i_hist", pTbin1, pTbin2));
      }
    }

    for(unsigned int pTupbin = 0; pTupbin < nPtBins_corr; pTupbin++)
    {
      for(unsigned int pTbin = 0; pTbin < nPtBins-1; pTbin++)
      {
        L0_L0bar_cosThetaProdPlane_pT_US_hist[pTupbin][pTbin] = (TH1D*)LLbarOutFile->Get(Form("L0_L0bar_cosThetaProdPlane_pT_US_pT1_%i_pT2_%i_hist", pTupbin, pTbin));
        L0_L0bar_cosThetaProdPlane_pT_LS_hist[pTupbin][pTbin] = (TH1D*)LLbarOutFile->Get(Form("L0_L0bar_cosThetaProdPlane_pT_LS_pT1_%i_pT2_%i_hist", pTupbin, pTbin));
        L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTupbin][pTbin] = (TH1D*)LLbarOutFile->Get(Form("L0_L0bar_cosThetaProdPlane_pT_ME_pT1_%i_pT2_%i_hist", pTupbin, pTbin));

        L0_L0_cosThetaProdPlane_pT_US_hist[pTupbin][pTbin] = (TH1D*)LLbarOutFile->Get(Form("L0_L0_cosThetaProdPlane_pT_US_pT1_%i_pT2_%i_hist", pTupbin, pTbin));
        L0_L0_cosThetaProdPlane_pT_LS_hist[pTupbin][pTbin] = (TH1D*)LLbarOutFile->Get(Form("L0_L0_cosThetaProdPlane_pT_LS_pT1_%i_pT2_%i_hist", pTupbin, pTbin));
        L0_L0_cosThetaProdPlane_pT_ME_hist[pTupbin][pTbin] = (TH1D*)LLbarOutFile->Get(Form("L0_L0_cosThetaProdPlane_pT_ME_pT1_%i_pT2_%i_hist", pTupbin, pTbin));

        L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTupbin][pTbin] = (TH1D*)LLbarOutFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_pT_US_pT1_%i_pT2_%i_hist", pTupbin, pTbin));
        L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTupbin][pTbin] = (TH1D*)LLbarOutFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_pT_LS_pT1_%i_pT2_%i_hist", pTupbin, pTbin));
        L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTupbin][pTbin] = (TH1D*)LLbarOutFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_pT_ME_pT1_%i_pT2_%i_hist", pTupbin, pTbin));
      }

      //one bin for each lower pT selection
      //histograms are as function of upper pT
      L0_L0bar_polarization_pT_hist[pTupbin] = (TH1D*)LLbarOutFile->Get(Form("L0_L0bar_polarization_pT_pTup_%i", pTupbin));
      L0_L0_polarization_pT_hist[pTupbin] = (TH1D*)LLbarOutFile->Get(Form("L0_L0_polarization_pT_pTup_%i", pTupbin));
      L0bar_L0bar_polarization_pT_hist[pTupbin] = (TH1D*)LLbarOutFile->Get(Form("L0bar_L0bar_polarization_pT_pTup_%i", pTupbin));
    }

    for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
    {
      for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
      {
        L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = (TH1D*)LLbarOutFile->Get(Form("L0_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
        L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = (TH1D*)LLbarOutFile->Get(Form("L0_L0bar_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
        L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = (TH1D*)LLbarOutFile->Get(Form("L0_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

        L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = (TH1D*)LLbarOutFile->Get(Form("L0_L0_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
        L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = (TH1D*)LLbarOutFile->Get(Form("L0_L0_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
        L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = (TH1D*)LLbarOutFile->Get(Form("L0_L0_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));

        L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2] = (TH1D*)LLbarOutFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_US_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
        L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2] = (TH1D*)LLbarOutFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_LS_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
        L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2] = (TH1D*)LLbarOutFile->Get(Form("L0bar_L0bar_cosThetaProdPlane_ME_eta1_%i_eta2_%i_hist", etaBin1, etaBin2));
      }
    }

    L0_L0bar_cosThetaProdPlane_US_hist = (TH1D*)LLbarOutFile->Get("L0_L0bar_cosThetaProdPlane_US_hist");
    L0_L0bar_cosThetaProdPlane_LS_hist = (TH1D*)LLbarOutFile->Get("L0_L0bar_cosThetaProdPlane_LS_hist");
    L0_L0bar_cosThetaProdPlane_ME_hist = (TH1D*)LLbarOutFile->Get("L0_L0bar_cosThetaProdPlane_ME_hist");

    L0_L0bar_delta_eta_vs_delta_phi_US_hist = (TH2D*)LLbarOutFile->Get("L0_L0bar_delta_eta_vs_delta_phi_US_hist");
    L0_L0bar_delta_eta_vs_delta_phi_LS_hist = (TH2D*)LLbarOutFile->Get("L0_L0bar_delta_eta_vs_delta_phi_LS_hist");

    L0_L0bar_y1_vs_y2_US_hist = (TH2D*)LLbarOutFile->Get("L0_L0bar_y1_vs_y2_US_hist");
    L0_L0bar_y1_vs_y2_LS_hist = (TH2D*)LLbarOutFile->Get("L0_L0bar_y1_vs_y2_LS_hist");

    L0_L0bar_pT1_vs_pT2_US_hist = (TH2D*)LLbarOutFile->Get("L0_L0bar_pT1_vs_pT2_US_hist");
    L0_L0bar_pT1_vs_pT2_LS_hist = (TH2D*)LLbarOutFile->Get("L0_L0bar_pT1_vs_pT2_LS_hist");

    L0_L0bar_phi1_vs_phi2_US_hist = (TH2D*)LLbarOutFile->Get("L0_L0bar_phi1_vs_phi2_US_hist");
    L0_L0bar_phi1_vs_phi2_LS_hist = (TH2D*)LLbarOutFile->Get("L0_L0bar_phi1_vs_phi2_LS_hist");


    L0_L0_cosThetaProdPlane_US_hist = (TH1D*)LLbarOutFile->Get("L0_L0_cosThetaProdPlane_US_hist");
    L0_L0_cosThetaProdPlane_LS_hist = (TH1D*)LLbarOutFile->Get("L0_L0_cosThetaProdPlane_LS_hist");
    L0_L0_cosThetaProdPlane_ME_hist = (TH1D*)LLbarOutFile->Get("L0_L0_cosThetaProdPlane_ME_hist");

    L0_L0_delta_eta_vs_delta_phi_US_hist = (TH2D*)LLbarOutFile->Get("L0_L0_delta_eta_vs_delta_phi_US_hist");
    L0_L0_delta_eta_vs_delta_phi_LS_hist = (TH2D*)LLbarOutFile->Get("L0_L0_delta_eta_vs_delta_phi_LS_hist");

    L0_L0_delta_eta_vs_delta_phi_US_zoom_hist = (TH2D*)LLbarOutFile->Get("L0_L0_delta_eta_vs_delta_phi_US_zoom_hist");
    L0_L0_delta_eta_vs_delta_phi_LS_zoom_hist = (TH2D*)LLbarOutFile->Get("L0_L0_delta_eta_vs_delta_phi_LS_zoom_hist");

    L0_L0_y1_vs_y2_US_hist = (TH2D*)LLbarOutFile->Get("L0_L0_y1_vs_y2_US_hist");
    L0_L0_y1_vs_y2_LS_hist = (TH2D*)LLbarOutFile->Get("L0_L0_y1_vs_y2_LS_hist");

    L0_p_L0_p_y1_vs_y2_US_hist = (TH2D*)LLbarOutFile->Get("L0_p_L0_p_y1_vs_y2_US_hist");
    L0_p_L0_p_y1_vs_y2_LS_hist = (TH2D*)LLbarOutFile->Get("L0_p_L0_p_y1_vs_y2_LS_hist");

    L0_pi_L0_pi_y1_vs_y2_US_hist = (TH2D*)LLbarOutFile->Get("L0_pi_L0_pi_y1_vs_y2_US_hist");
    L0_pi_L0_pi_y1_vs_y2_LS_hist = (TH2D*)LLbarOutFile->Get("L0_pi_L0_pi_y1_vs_y2_LS_hist");

    L0_L0_pT1_vs_pT2_US_hist = (TH2D*)LLbarOutFile->Get("L0_L0_pT1_vs_pT2_US_hist");
    L0_L0_pT1_vs_pT2_LS_hist = (TH2D*)LLbarOutFile->Get("L0_L0_pT1_vs_pT2_LS_hist");

    L0_L0_phi1_vs_phi2_US_hist = (TH2D*)LLbarOutFile->Get("L0_L0_phi1_vs_phi2_US_hist");
    L0_L0_phi1_vs_phi2_LS_hist = (TH2D*)LLbarOutFile->Get("L0_L0_phi1_vs_phi2_LS_hist");

    L0_p_L0_p_phi1_vs_phi2_US_hist = (TH2D*)LLbarOutFile->Get("L0_p_L0_p_phi1_vs_phi2_US_hist");
    L0_p_L0_p_phi1_vs_phi2_LS_hist = (TH2D*)LLbarOutFile->Get("L0_p_L0_p_phi1_vs_phi2_LS_hist");

    L0_pi_L0_pi_phi1_vs_phi2_US_hist = (TH2D*)LLbarOutFile->Get("L0_pi_L0_pi_phi1_vs_phi2_US_hist");
    L0_pi_L0_pi_phi1_vs_phi2_LS_hist = (TH2D*)LLbarOutFile->Get("L0_pi_L0_pi_phi1_vs_phi2_LS_hist");


    L0bar_L0bar_cosThetaProdPlane_US_hist = (TH1D*)LLbarOutFile->Get("L0bar_L0bar_cosThetaProdPlane_US_hist");
    L0bar_L0bar_cosThetaProdPlane_LS_hist = (TH1D*)LLbarOutFile->Get("L0bar_L0bar_cosThetaProdPlane_LS_hist");
    L0bar_L0bar_cosThetaProdPlane_ME_hist = (TH1D*)LLbarOutFile->Get("L0bar_L0bar_cosThetaProdPlane_ME_hist");

    L0bar_L0bar_delta_eta_vs_delta_phi_US_hist = (TH2D*)LLbarOutFile->Get("L0bar_L0bar_delta_eta_vs_delta_phi_US_hist");
    L0bar_L0bar_delta_eta_vs_delta_phi_LS_hist = (TH2D*)LLbarOutFile->Get("L0bar_L0bar_delta_eta_vs_delta_phi_LS_hist");

    L0bar_L0bar_delta_eta_vs_delta_phi_US_zoom_hist = (TH2D*)LLbarOutFile->Get("L0bar_L0bar_delta_eta_vs_delta_phi_US_zoom_hist");
    L0bar_L0bar_delta_eta_vs_delta_phi_LS_zoom_hist = (TH2D*)LLbarOutFile->Get("L0bar_L0bar_delta_eta_vs_delta_phi_LS_zoom_hist");

    L0bar_L0bar_y1_vs_y2_US_hist = (TH2D*)LLbarOutFile->Get("L0bar_L0bar_y1_vs_y2_US_hist");
    L0bar_L0bar_y1_vs_y2_LS_hist = (TH2D*)LLbarOutFile->Get("L0bar_L0bar_y1_vs_y2_LS_hist");

    L0bar_p_L0bar_p_y1_vs_y2_US_hist = (TH2D*)LLbarOutFile->Get("L0bar_p_L0bar_p_y1_vs_y2_US_hist");
    L0bar_p_L0bar_p_y1_vs_y2_LS_hist = (TH2D*)LLbarOutFile->Get("L0bar_p_L0bar_p_y1_vs_y2_LS_hist");

    L0bar_pi_L0bar_pi_y1_vs_y2_US_hist = (TH2D*)LLbarOutFile->Get("L0bar_pi_L0bar_pi_y1_vs_y2_US_hist");
    L0bar_pi_L0bar_pi_y1_vs_y2_LS_hist = (TH2D*)LLbarOutFile->Get("L0bar_pi_L0bar_pi_y1_vs_y2_LS_hist");

    L0bar_L0bar_pT1_vs_pT2_US_hist = (TH2D*)LLbarOutFile->Get("L0bar_L0bar_pT1_vs_pT2_US_hist");
    L0bar_L0bar_pT1_vs_pT2_LS_hist = (TH2D*)LLbarOutFile->Get("L0bar_L0bar_pT1_vs_pT2_LS_hist");

    L0bar_L0bar_phi1_vs_phi2_US_hist = (TH2D*)LLbarOutFile->Get("L0bar_L0bar_phi1_vs_phi2_US_hist");
    L0bar_L0bar_phi1_vs_phi2_LS_hist = (TH2D*)LLbarOutFile->Get("L0bar_L0bar_phi1_vs_phi2_LS_hist");

    L0bar_p_L0bar_p_phi1_vs_phi2_US_hist = (TH2D*)LLbarOutFile->Get("L0bar_p_L0bar_p_phi1_vs_phi2_US_hist");
    L0bar_p_L0bar_p_phi1_vs_phi2_LS_hist = (TH2D*)LLbarOutFile->Get("L0bar_p_L0bar_p_phi1_vs_phi2_LS_hist");

    L0bar_pi_L0bar_pi_phi1_vs_phi2_US_hist = (TH2D*)LLbarOutFile->Get("L0bar_pi_L0bar_pi_phi1_vs_phi2_US_hist");
    L0bar_pi_L0bar_pi_phi1_vs_phi2_LS_hist = (TH2D*)LLbarOutFile->Get("L0bar_pi_L0bar_pi_phi1_vs_phi2_LS_hist");


  }

  //________________________________________________________________________________________

  //Labmda and Lambda_bar stats

  int nLambda = -1; //default value
  int nLambdaBar = -1;

  int nLambda_ME = -1;
  int nLambdaBar_ME = -1;


  int nLambda_background = -1; //default value
  int nLambdaBar_background = -1;


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
  vector<TLorentzVector> L_vector;
  vector<int> L_pT_bin_vector;
  vector<int> L_eta_bin_vector;
  vector<TLorentzVector> p_vector;
  vector<TLorentzVector> pi_vector;

  vector<TLorentzVector> Lbar_vector;
  vector<int> Lbar_pT_bin_vector;
  vector<int> Lbar_eta_bin_vector;
  vector<TLorentzVector> pbar_vector;
  vector<TLorentzVector> pibar_vector;


  vector<TLorentzVector> L_vector_background;
  vector<int> L_pT_bin_vector_background;
  vector<int> L_eta_bin_vector_background;
  vector<TLorentzVector> p_vector_background;
  vector<TLorentzVector> pi_vector_background;

  vector<TLorentzVector> Lbar_vector_background;
  vector<int> Lbar_pT_bin_vector_background;
  vector<int> Lbar_eta_bin_vector_background;
  vector<TLorentzVector> pbar_vector_background;
  vector<TLorentzVector> pibar_vector_background;

  //vectors for mixed-event
  vector<TLorentzVector> L_vector_ME;
  vector<int> L_pT_bin_vector_ME;
  vector<int> L_eta_bin_vector_ME;
  vector<TLorentzVector> p_vector_ME;

  vector<TLorentzVector> Lbar_vector_ME;
  vector<int> Lbar_pT_bin_vector_ME;
  vector<int> Lbar_eta_bin_vector_ME;
  vector<TLorentzVector> pbar_vector_ME;

  //TLorentzVector *L_Lorentz_vector = new TLorentzVector(1,1,1,1);
  //TLorentzVector *p_Lorenz_vector = new TLorentzVector(1,1,1,1);

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



    if(nLambda == -1 ) //first iteration
    {
      eventID_last = eventId;
      eventID_last_background = eventId;

      //cout<<eventId<<endl;

      nLambda = 0;
      nLambdaBar = 0;

      nLambda_background = 0;
      nLambdaBar_background = 0;
    }

    //------------------------------------------------------------------------------------------------------------------



    //----------------------------------------------------------------------------------------------------------------



    //fill all histograms for all pT and centrality bins
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

    //if( pT_bin == -1 ) continue;


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

    //if( eta_bin == -1 ) continue;

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


    if(charge == 0 ) //like-sign combinations
    {

      if(eventId == eventID_last) //same event as in previous iteration and first event
      {
        //cuts
        if( cuts(requireTOF_Lcorr, pi_hasTOFinfo, p_hasTOFinfo, L_y) && pT_bin != -1 && eta_bin != -1 && pT_bin_corr != -1)
        {
          if( p_ch == 1 && L_mass > invMassRange_L0[0][pT_bin][eta_bin] && L_mass < invMassRange_L0[1][pT_bin][eta_bin])
          {
            L_vector.push_back(L_Lorentz_vector);
            L_pT_bin_vector.push_back(pT_bin_corr);
            L_eta_bin_vector.push_back(eta_bin);
            p_vector.push_back(p_Lorenz_vector);
            pi_vector.push_back(pi_Lorenz_vector);

            nLambda++;
          }
          else if( p_ch == -1 && L_mass > invMassRange_L0bar[0][pT_bin][eta_bin] && L_mass < invMassRange_L0bar[1][pT_bin][eta_bin])
          {
            Lbar_vector.push_back(L_Lorentz_vector);
            Lbar_pT_bin_vector.push_back(pT_bin_corr);
            Lbar_eta_bin_vector.push_back(eta_bin);
            pbar_vector.push_back(p_Lorenz_vector);
            pibar_vector.push_back(pi_Lorenz_vector);

            nLambdaBar++;
          }

        }

      }
      else if(eventId != eventID_last) //new event
      {


        //at least one L-Lbar pair in event
        if(nLambda > 0 && nLambdaBar > 0)
        {
          //cout<<eventId<<endl;

          for(unsigned int iLambda = 0; iLambda < L_vector.size(); iLambda++)
          {
            for(unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector.size(); iLambdaBar++)
            {
              double L_Lbar_pairThetaStar = LpairThetaStar(L_vector.at(iLambda), p_vector.at(iLambda), Lbar_vector.at(iLambdaBar), pbar_vector.at(iLambdaBar));

              //if( fabs(p_vector.at(iLambda).Rapidity()) < 0.001  ) continue;

              L0_L0bar_cosThetaProdPlane_US_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar));
              L0_L0bar_cosThetaProdPlane_pT_US_hist[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(TMath::Cos(L_Lbar_pairThetaStar));
              //if( fabs(L_vector.at(iLambda).Rapidity()) < 0.2 && fabs(Lbar_vector.at(iLambdaBar).Rapidity()) < 0.2 ) L0_L0bar_cosThetaProdPlane_pT_US_hist[L_pT_bin_vector.at(iLambda)][Lbar_pT_bin_vector.at(iLambdaBar)]->Fill(TMath::Cos(L_Lbar_pairThetaStar));
              L0_L0bar_cosThetaProdPlane_eta_US_hist[L_eta_bin_vector.at(iLambda)][Lbar_eta_bin_vector.at(iLambdaBar)]->Fill(TMath::Cos(L_Lbar_pairThetaStar));
              //_____________________________________________________________________________

              double delta_eta = fabs( L_vector.at(iLambda).Eta() - Lbar_vector.at(iLambdaBar).Eta() );
              double delta_phi = fabs( L_vector.at(iLambda).Phi() - Lbar_vector.at(iLambdaBar).Phi() );

              L0_L0bar_delta_eta_vs_delta_phi_US_hist->Fill(delta_eta, delta_phi);
              //___________________________________________________________________________________

              L0_L0bar_y1_vs_y2_US_hist->Fill( L_vector.at(iLambda).Rapidity(), Lbar_vector.at(iLambdaBar).Rapidity() );

              L0_L0bar_pT1_vs_pT2_US_hist->Fill( L_vector.at(iLambda).Pt(), Lbar_vector.at(iLambdaBar).Pt() );

              L0_L0bar_phi1_vs_phi2_US_hist->Fill( L_vector.at(iLambda).Phi(), Lbar_vector.at(iLambdaBar).Phi() );

              //invariant mass
              L_inv_mass_US_from_L_Lbar->Fill(L_vector.at(iLambda).M());
              Lbar_inv_mass_US_from_L_Lbar->Fill(Lbar_vector.at(iLambdaBar).M());

              L_vs_Lbar_mass_US->Fill(L_vector.at(iLambda).M(),Lbar_vector.at(iLambdaBar).M());


              nEventsWithLambdaPair_US_hist->Fill(1.5);
            }
          }

        }
        else
        {
          nEventsWithLambdaPair_US_hist->Fill(0.5);
        }

        //at least one L0-L0 pair in event
        if(nLambda > 1)
        {
          for(unsigned int iLambda1 = 0; iLambda1 < L_vector.size(); iLambda1++)
          {
            for(unsigned int iLambda2 = iLambda1+1; iLambda2 < L_vector.size(); iLambda2++)
            {
              //if( fabs(p_vector.at(iLambda1).Rapidity()) < 0.001 || fabs(p_vector.at(iLambda2).Rapidity()) < 0.001 ) continue;

              if(p_vector.at(iLambda1).Rapidity() == p_vector.at(iLambda2).Rapidity()) continue;
              if(pi_vector.at(iLambda1).Rapidity() == pi_vector.at(iLambda2).Rapidity()) continue;

              double L_L_pairThetaStar = LpairThetaStar(L_vector.at(iLambda1), p_vector.at(iLambda1), L_vector.at(iLambda2), p_vector.at(iLambda2));

              L0_L0_cosThetaProdPlane_US_hist->Fill(TMath::Cos(L_L_pairThetaStar));
              L0_L0_cosThetaProdPlane_pT_US_hist[L_pT_bin_vector.at(iLambda1)][L_pT_bin_vector.at(iLambda2)]->Fill(TMath::Cos(L_L_pairThetaStar));
              L0_L0_cosThetaProdPlane_eta_US_hist[L_eta_bin_vector.at(iLambda1)][L_eta_bin_vector.at(iLambda2)]->Fill(TMath::Cos(L_L_pairThetaStar));
              //_____________________________________________________________________________

              double delta_eta = fabs( L_vector.at(iLambda1).Eta() - L_vector.at(iLambda2).Eta() );
              double delta_phi = fabs( L_vector.at(iLambda1).Phi() - L_vector.at(iLambda2).Phi() );

              L0_L0_delta_eta_vs_delta_phi_US_hist->Fill(delta_eta, delta_phi);
              L0_L0_delta_eta_vs_delta_phi_US_zoom_hist->Fill(delta_eta, delta_phi);
              //___________________________________________________________________________________

              L0_L0_y1_vs_y2_US_hist->Fill( L_vector.at(iLambda1).Rapidity(), L_vector.at(iLambda2).Rapidity() );
              L0_p_L0_p_y1_vs_y2_US_hist->Fill( p_vector.at(iLambda1).Rapidity(), p_vector.at(iLambda2).Rapidity() );
              L0_pi_L0_pi_y1_vs_y2_US_hist->Fill( pi_vector.at(iLambda1).Rapidity(), pi_vector.at(iLambda2).Rapidity() );

              L0_L0_pT1_vs_pT2_US_hist->Fill( L_vector.at(iLambda1).Pt(), L_vector.at(iLambda2).Pt() );

              L0_L0_phi1_vs_phi2_US_hist->Fill( L_vector.at(iLambda1).Phi(), L_vector.at(iLambda2).Phi() );
              L0_p_L0_p_phi1_vs_phi2_US_hist->Fill( p_vector.at(iLambda1).Phi(), p_vector.at(iLambda2).Phi() );
              L0_pi_L0_pi_phi1_vs_phi2_US_hist->Fill( pi_vector.at(iLambda1).Phi(), pi_vector.at(iLambda2).Phi() );

              //invariant mass
              L_inv_mass_US_from_L_L->Fill(L_vector.at(iLambda1).M());
              L_inv_mass_US_from_L_L->Fill(L_vector.at(iLambda2).M());

              L_vs_L_mass_US->Fill(L_vector.at(iLambda1).M(),L_vector.at(iLambda2).M());

            }
          }
        }

        //at least one L0bar-L0bar pair in event
        if(nLambdaBar > 1)
        {
          for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < Lbar_vector.size(); iLambdaBar1++)
          {
            for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < Lbar_vector.size(); iLambdaBar2++)
            {
              if(pbar_vector.at(iLambdaBar1).Rapidity() == pbar_vector.at(iLambdaBar2).Rapidity()) continue;
              if(pibar_vector.at(iLambdaBar1).Rapidity() == pibar_vector.at(iLambdaBar2).Rapidity()) continue;

              double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_vector.at(iLambdaBar1), pbar_vector.at(iLambdaBar1), Lbar_vector.at(iLambdaBar2), pbar_vector.at(iLambdaBar2));

              L0bar_L0bar_cosThetaProdPlane_US_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar));
              L0bar_L0bar_cosThetaProdPlane_pT_US_hist[Lbar_pT_bin_vector.at(iLambdaBar1)][Lbar_pT_bin_vector.at(iLambdaBar2)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar));
              L0bar_L0bar_cosThetaProdPlane_eta_US_hist[Lbar_eta_bin_vector.at(iLambdaBar1)][Lbar_eta_bin_vector.at(iLambdaBar2)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar));
              //_____________________________________________________________________________

              double delta_eta = fabs( Lbar_vector.at(iLambdaBar1).Eta() - Lbar_vector.at(iLambdaBar2).Eta() );
              double delta_phi = fabs( Lbar_vector.at(iLambdaBar1).Phi() - Lbar_vector.at(iLambdaBar2).Phi() );

              L0bar_L0bar_delta_eta_vs_delta_phi_US_hist->Fill(delta_eta, delta_phi);
              L0bar_L0bar_delta_eta_vs_delta_phi_US_zoom_hist->Fill(delta_eta, delta_phi);
              //___________________________________________________________________________________

              L0bar_L0bar_y1_vs_y2_US_hist->Fill( Lbar_vector.at(iLambdaBar1).Rapidity(), Lbar_vector.at(iLambdaBar2).Rapidity() );
              L0bar_p_L0bar_p_y1_vs_y2_US_hist->Fill( pbar_vector.at(iLambdaBar1).Rapidity(), pbar_vector.at(iLambdaBar2).Rapidity() );
              L0bar_pi_L0bar_pi_y1_vs_y2_US_hist->Fill( pibar_vector.at(iLambdaBar1).Rapidity(), pibar_vector.at(iLambdaBar2).Rapidity() );

              L0bar_L0bar_pT1_vs_pT2_US_hist->Fill( Lbar_vector.at(iLambdaBar1).Pt(), Lbar_vector.at(iLambdaBar2).Pt() );

              L0bar_L0bar_phi1_vs_phi2_US_hist->Fill( Lbar_vector.at(iLambdaBar1).Phi(), Lbar_vector.at(iLambdaBar2).Phi() );
              L0bar_p_L0bar_p_phi1_vs_phi2_US_hist->Fill( pbar_vector.at(iLambdaBar1).Phi(), pbar_vector.at(iLambdaBar2).Phi() );
              L0bar_pi_L0bar_pi_phi1_vs_phi2_US_hist->Fill( pibar_vector.at(iLambdaBar1).Phi(), pibar_vector.at(iLambdaBar2).Phi() );

              //invariant mass
              Lbar_inv_mass_US_from_Lbar_Lbar->Fill(Lbar_vector.at(iLambdaBar1).M());
              Lbar_inv_mass_US_from_Lbar_Lbar->Fill(Lbar_vector.at(iLambdaBar2).M());

              Lbar_vs_Lbar_mass_US->Fill(Lbar_vector.at(iLambdaBar1).M(),Lbar_vector.at(iLambdaBar2).M());

            }
          }
        }
        //_________________________________________________________________________________________________________

        //fill vectors for mixed event
        //selecting events with only one L or L-bar
        //cout<<"nLambda = "<<nLambda<<endl;
        //cout<<"L vector size = "<<L_vector.size()<<endl;

        //cout<<"nLambdaBar = "<<nLambdaBar<<endl;
        //cout<<"Lbar vector size = "<<Lbar_vector.size()<<endl;

        if( nLambda == 1 && nLambdaBar == 0 && L_vector_ME.size() < 1e4)
        {
          L_vector_ME.push_back(L_vector.at(0));
          p_vector_ME.push_back(p_vector.at(0));

          L_pT_bin_vector_ME.push_back(L_pT_bin_vector.at(0));
          L_eta_bin_vector_ME.push_back(L_eta_bin_vector.at(0));

        }

        if( nLambda == 0 && nLambdaBar == 1 && Lbar_vector_ME.size() < 1e4)
        {
          Lbar_vector_ME.push_back(Lbar_vector.at(0));
          pbar_vector_ME.push_back(pbar_vector.at(0));

          Lbar_pT_bin_vector_ME.push_back(Lbar_pT_bin_vector.at(0));
          Lbar_eta_bin_vector_ME.push_back(Lbar_eta_bin_vector.at(0));

        }
        //_________________________________________________________________________________________________________



        nLambdasInEvent_US_hist->Fill(nLambda);
        nLambdaBarsInEvent_US_hist->Fill(nLambdaBar);


        //reset number of lambdas and vectors
        L_vector.clear();
        L_pT_bin_vector.clear();
        L_eta_bin_vector.clear();
        p_vector.clear();
        pi_vector.clear();

        Lbar_vector.clear();
        Lbar_pT_bin_vector.clear();
        Lbar_eta_bin_vector.clear();
        pbar_vector.clear();
        pibar_vector.clear();

        eventID_last = eventId;


        //fill L or Lbar from new event
        //need to check cuts again in the new event
        if( cuts(requireTOF_Lcorr,pi_hasTOFinfo, p_hasTOFinfo, L_y) && pT_bin != -1 && eta_bin != -1 && pT_bin_corr != -1)
        {
          if( p_ch == 1 && L_mass > invMassRange_L0[0][pT_bin][eta_bin] && L_mass < invMassRange_L0[1][pT_bin][eta_bin])
          {

            L_vector.push_back(L_Lorentz_vector);
            L_pT_bin_vector.push_back(pT_bin_corr);
            L_eta_bin_vector.push_back(eta_bin);
            p_vector.push_back(p_Lorenz_vector);
            pi_vector.push_back(pi_Lorenz_vector);

            nLambda = 1;
            nLambdaBar = 0;

          }
          else if( p_ch == -1 && L_mass > invMassRange_L0bar[0][pT_bin][eta_bin] && L_mass < invMassRange_L0bar[1][pT_bin][eta_bin])
          {

            Lbar_vector.push_back(L_Lorentz_vector);
            Lbar_pT_bin_vector.push_back(pT_bin_corr);
            Lbar_eta_bin_vector.push_back(eta_bin);
            pbar_vector.push_back(p_Lorenz_vector);
            pibar_vector.push_back(pi_Lorenz_vector);

            nLambda = 0;
            nLambdaBar = 1;
          }
          else
          {
            nLambda = 0;
            nLambdaBar = 0;
          }

        }
        else
        {
          nLambda = 0;
          nLambdaBar = 0;
        }

      }

    }
    else //unlike-sign combinations
    {

      if(eventId == eventID_last_background) //same event as in previous iteration
      {
        //cuts
        if( cuts(requireTOF_Lcorr,pi_hasTOFinfo, p_hasTOFinfo, L_y) && pT_bin != -1 && eta_bin != -1 && pT_bin_corr != -1)
        {
          if( p_ch == 1 && L_mass > invMassRange_L0[0][pT_bin][eta_bin] && L_mass < invMassRange_L0[1][pT_bin][eta_bin])
          {
            L_vector_background.push_back(L_Lorentz_vector);
            L_pT_bin_vector_background.push_back(pT_bin_corr);
            L_eta_bin_vector_background.push_back(eta_bin);
            p_vector_background.push_back(p_Lorenz_vector);
            pi_vector_background.push_back(pi_Lorenz_vector);

            nLambda_background++;
          }
          else if( p_ch == -1 && L_mass > invMassRange_L0bar[0][pT_bin][eta_bin] && L_mass < invMassRange_L0bar[1][pT_bin][eta_bin])
          {
            Lbar_vector_background.push_back(L_Lorentz_vector);
            Lbar_pT_bin_vector_background.push_back(pT_bin_corr);
            Lbar_eta_bin_vector_background.push_back(eta_bin);
            pbar_vector_background.push_back(p_Lorenz_vector);
            pibar_vector_background.push_back(pi_Lorenz_vector);

            nLambdaBar_background++;
          }
        }
      }
      else if(eventId != eventID_last_background) //new event
      {
        eventID_last_background = eventId;

        if(nLambda_background > 0 && nLambdaBar_background > 0)
        {
          //nEventsWithLambdaPair_LS_hist->Fill(1.5);

          for(unsigned int iLambda = 0; iLambda < L_vector_background.size(); iLambda++)
          {
            for(unsigned int iLambdaBar = 0; iLambdaBar < Lbar_vector_background.size(); iLambdaBar++)
            {
              //if( fabs(p_vector_background.at(iLambda).Rapidity()) < 0.001  ) continue;

              double L_Lbar_pairThetaStar = LpairThetaStar(L_vector_background.at(iLambda), p_vector_background.at(iLambda), Lbar_vector_background.at(iLambdaBar), pbar_vector_background.at(iLambdaBar));

              L0_L0bar_cosThetaProdPlane_LS_hist->Fill(TMath::Cos(L_Lbar_pairThetaStar));
              L0_L0bar_cosThetaProdPlane_pT_LS_hist[L_pT_bin_vector_background.at(iLambda)][Lbar_pT_bin_vector_background.at(iLambdaBar)]->Fill(TMath::Cos(L_Lbar_pairThetaStar));
              //if( L_vector_background.at(iLambda).Rapidity() < 0.2 && Lbar_vector_background.at(iLambdaBar).Rapidity() < 0.2 ) L0_L0bar_cosThetaProdPlane_pT_LS_hist[L_pT_bin_vector_background.at(iLambda)][Lbar_pT_bin_vector_background.at(iLambdaBar)]->Fill(TMath::Cos(L_Lbar_pairThetaStar));
              L0_L0bar_cosThetaProdPlane_eta_LS_hist[L_eta_bin_vector_background.at(iLambda)][Lbar_eta_bin_vector_background.at(iLambdaBar)]->Fill(TMath::Cos(L_Lbar_pairThetaStar));

              double delta_eta = fabs( L_vector_background.at(iLambda).Eta() - Lbar_vector_background.at(iLambdaBar).Eta() );
              double delta_phi = fabs( L_vector_background.at(iLambda).Phi() - Lbar_vector_background.at(iLambdaBar).Phi() );

              L0_L0bar_delta_eta_vs_delta_phi_LS_hist->Fill(delta_eta, delta_phi);
              //___________________________________________________________________________________

              L0_L0bar_y1_vs_y2_LS_hist->Fill( L_vector_background.at(iLambda).Rapidity(), Lbar_vector_background.at(iLambdaBar).Rapidity() );

              L0_L0bar_pT1_vs_pT2_LS_hist->Fill( L_vector_background.at(iLambda).Pt(), Lbar_vector_background.at(iLambdaBar).Pt() );

              L0_L0bar_phi1_vs_phi2_LS_hist->Fill( L_vector_background.at(iLambda).Phi(), Lbar_vector_background.at(iLambdaBar).Phi() );

              //invariant mass
              L_inv_mass_LS_from_L_Lbar->Fill(L_vector_background.at(iLambda).M());
              Lbar_inv_mass_LS_from_L_Lbar->Fill(Lbar_vector_background.at(iLambdaBar).M());

              L_vs_Lbar_mass_LS->Fill(L_vector_background.at(iLambda).M(),Lbar_vector_background.at(iLambdaBar).M());

              nEventsWithLambdaPair_LS_hist->Fill(1.5);
            }
          }
        }
        else
        {
          nEventsWithLambdaPair_LS_hist->Fill(0.5);
        }


        //at least one L0-L0 pair in event
        if(nLambda_background > 1)
        {
          for(unsigned int iLambda1 = 0; iLambda1 < L_vector_background.size(); iLambda1++)
          {
            for(unsigned int iLambda2 = iLambda1+1; iLambda2 < L_vector_background.size(); iLambda2++)
            {
              //if( fabs(p_vector_background.at(iLambda1).Rapidity()) < 0.001 || fabs(p_vector_background.at(iLambda2).Rapidity()) < 0.001 ) continue;

              if(p_vector_background.at(iLambda1).Rapidity() == p_vector_background.at(iLambda2).Rapidity()) continue;
              if(pi_vector_background.at(iLambda1).Rapidity() == pi_vector_background.at(iLambda2).Rapidity()) continue;

              double L_L_pairThetaStar = LpairThetaStar(L_vector_background.at(iLambda1), p_vector_background.at(iLambda1), L_vector_background.at(iLambda2), p_vector_background.at(iLambda2));

              L0_L0_cosThetaProdPlane_LS_hist->Fill(TMath::Cos(L_L_pairThetaStar));
              L0_L0_cosThetaProdPlane_pT_LS_hist[L_pT_bin_vector_background.at(iLambda1)][L_pT_bin_vector_background.at(iLambda2)]->Fill(TMath::Cos(L_L_pairThetaStar));
              L0_L0_cosThetaProdPlane_eta_LS_hist[L_eta_bin_vector_background.at(iLambda1)][L_eta_bin_vector_background.at(iLambda2)]->Fill(TMath::Cos(L_L_pairThetaStar));
              //_____________________________________________________________________________

              double delta_eta = fabs( L_vector_background.at(iLambda1).Eta() - L_vector_background.at(iLambda2).Eta() );
              double delta_phi = fabs( L_vector_background.at(iLambda1).Phi() - L_vector_background.at(iLambda2).Phi() );

              L0_L0_delta_eta_vs_delta_phi_LS_hist->Fill(delta_eta, delta_phi);
              L0_L0_delta_eta_vs_delta_phi_LS_zoom_hist->Fill(delta_eta, delta_phi);
              //___________________________________________________________________________________

              L0_L0_y1_vs_y2_LS_hist->Fill( L_vector_background.at(iLambda1).Rapidity(), L_vector_background.at(iLambda2).Rapidity() );
              L0_p_L0_p_y1_vs_y2_LS_hist->Fill( p_vector_background.at(iLambda1).Rapidity(), p_vector_background.at(iLambda2).Rapidity() );
              L0_pi_L0_pi_y1_vs_y2_LS_hist->Fill( pi_vector_background.at(iLambda1).Rapidity(), pi_vector_background.at(iLambda2).Rapidity() );

              L0_L0_pT1_vs_pT2_LS_hist->Fill( L_vector_background.at(iLambda1).Pt(), L_vector_background.at(iLambda2).Pt() );

              L0_L0_phi1_vs_phi2_LS_hist->Fill( L_vector_background.at(iLambda1).Phi(), L_vector_background.at(iLambda2).Phi() );
              L0_p_L0_p_phi1_vs_phi2_LS_hist->Fill( p_vector_background.at(iLambda1).Phi(), p_vector_background.at(iLambda2).Phi() );
              L0_pi_L0_pi_phi1_vs_phi2_LS_hist->Fill( pi_vector_background.at(iLambda1).Phi(), pi_vector_background.at(iLambda2).Phi() );

              //invariant mass
              L_inv_mass_LS_from_L_L->Fill(L_vector_background.at(iLambda1).M());
              L_inv_mass_LS_from_L_L->Fill(L_vector_background.at(iLambda2).M());

              L_vs_L_mass_LS->Fill(L_vector_background.at(iLambda1).M(), L_vector_background.at(iLambda2).M());

            }
          }
        }

        //at least one L0bar-L0bar pair in event
        if(nLambdaBar_background > 1)
        {
          for(unsigned int iLambdaBar1 = 0; iLambdaBar1 < Lbar_vector_background.size(); iLambdaBar1++)
          {
            for(unsigned int iLambdaBar2 = iLambdaBar1+1; iLambdaBar2 < Lbar_vector_background.size(); iLambdaBar2++)
            {
              if(pbar_vector_background.at(iLambdaBar1).Rapidity() == pbar_vector_background.at(iLambdaBar2).Rapidity()) continue;
              if(pibar_vector_background.at(iLambdaBar1).Rapidity() == pibar_vector_background.at(iLambdaBar2).Rapidity()) continue;

              double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_vector_background.at(iLambdaBar1), pbar_vector_background.at(iLambdaBar1), Lbar_vector_background.at(iLambdaBar2), pbar_vector_background.at(iLambdaBar2));

              L0bar_L0bar_cosThetaProdPlane_LS_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar));
              L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[Lbar_pT_bin_vector_background.at(iLambdaBar1)][Lbar_pT_bin_vector_background.at(iLambdaBar2)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar));
              L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[Lbar_eta_bin_vector_background.at(iLambdaBar1)][Lbar_eta_bin_vector_background.at(iLambdaBar2)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar));
              //_____________________________________________________________________________

              double delta_eta = fabs( Lbar_vector_background.at(iLambdaBar1).Eta() - Lbar_vector_background.at(iLambdaBar2).Eta() );
              double delta_phi = fabs( Lbar_vector_background.at(iLambdaBar1).Phi() - Lbar_vector_background.at(iLambdaBar2).Phi() );

              L0bar_L0bar_delta_eta_vs_delta_phi_LS_hist->Fill(delta_eta, delta_phi);
              L0bar_L0bar_delta_eta_vs_delta_phi_LS_zoom_hist->Fill(delta_eta, delta_phi);
              //___________________________________________________________________________________

              L0bar_L0bar_y1_vs_y2_LS_hist->Fill( Lbar_vector_background.at(iLambdaBar1).Rapidity(), Lbar_vector_background.at(iLambdaBar2).Rapidity() );
              L0bar_p_L0bar_p_y1_vs_y2_LS_hist->Fill( pbar_vector_background.at(iLambdaBar1).Rapidity(), pbar_vector_background.at(iLambdaBar2).Rapidity() );
              L0bar_pi_L0bar_pi_y1_vs_y2_LS_hist->Fill( pibar_vector_background.at(iLambdaBar1).Rapidity(), pibar_vector_background.at(iLambdaBar2).Rapidity() );

              L0bar_L0bar_pT1_vs_pT2_LS_hist->Fill( Lbar_vector_background.at(iLambdaBar1).Pt(), Lbar_vector_background.at(iLambdaBar2).Pt() );

              L0bar_L0bar_phi1_vs_phi2_LS_hist->Fill( Lbar_vector_background.at(iLambdaBar1).Phi(), Lbar_vector_background.at(iLambdaBar2).Phi() );
              L0bar_p_L0bar_p_phi1_vs_phi2_LS_hist->Fill( pbar_vector_background.at(iLambdaBar1).Phi(), pbar_vector_background.at(iLambdaBar2).Phi() );
              L0bar_pi_L0bar_pi_phi1_vs_phi2_LS_hist->Fill( pibar_vector_background.at(iLambdaBar1).Phi(), pibar_vector_background.at(iLambdaBar2).Phi() );

              //invariant mass
              Lbar_inv_mass_LS_from_Lbar_Lbar->Fill(Lbar_vector_background.at(iLambdaBar1).M());
              Lbar_inv_mass_LS_from_Lbar_Lbar->Fill(Lbar_vector_background.at(iLambdaBar2).M());

              Lbar_vs_Lbar_mass_LS->Fill(Lbar_vector_background.at(iLambdaBar1).M(),Lbar_vector_background.at(iLambdaBar2).M());

            }
          }
        }

        nLambdasInEvent_LS_hist->Fill(nLambda_background);
        nLambdaBarsInEvent_LS_hist->Fill(nLambdaBar_background);


        //reset number of lambdas and vectors
        L_vector_background.clear();
        L_pT_bin_vector_background.clear();
        L_eta_bin_vector_background.clear();
        p_vector_background.clear();
        pi_vector_background.clear();

        Lbar_vector_background.clear();
        Lbar_pT_bin_vector_background.clear();
        Lbar_eta_bin_vector_background.clear();
        pbar_vector_background.clear();
        pibar_vector_background.clear();

        //cuts
        if( cuts(requireTOF_Lcorr,pi_hasTOFinfo, p_hasTOFinfo, L_y) && pT_bin != -1 && eta_bin != -1 && pT_bin_corr != -1)
        {
          if( p_ch == 1 && L_mass > invMassRange_L0[0][pT_bin][eta_bin] && L_mass < invMassRange_L0[1][pT_bin][eta_bin])
          {
            L_vector_background.push_back(L_Lorentz_vector);
            L_pT_bin_vector_background.push_back(pT_bin_corr);
            L_eta_bin_vector_background.push_back(eta_bin);
            p_vector_background.push_back(p_Lorenz_vector);
            pi_vector_background.push_back(pi_Lorenz_vector);

            nLambda_background = 1;
            nLambdaBar_background = 0;
          }
          else if( p_ch == -1 && L_mass > invMassRange_L0bar[0][pT_bin][eta_bin] && L_mass < invMassRange_L0bar[1][pT_bin][eta_bin])
          {
            Lbar_vector_background.push_back(L_Lorentz_vector);
            Lbar_pT_bin_vector_background.push_back(pT_bin_corr);
            Lbar_eta_bin_vector_background.push_back(eta_bin);
            pbar_vector_background.push_back(p_Lorenz_vector);
            pibar_vector_background.push_back(pi_Lorenz_vector);

            nLambda_background = 0;
            nLambdaBar_background = 1;
          }
          else
          {
            nLambda_background = 0;
            nLambdaBar_background = 0;
          }

        }
        else
        {
          nLambda_background = 0;
          nLambdaBar_background = 0;
        }

      }


    }//end else for flag_int check
  }//end loop over entries in NTuple

  //________________________________________________________________________________________________________

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  //analyze mixed-event

  //L-Lbar
  for(unsigned int iLambda_ME = 0; iLambda_ME < L_vector_ME.size(); iLambda_ME++)
  {
    for(unsigned int iLambdaBar_ME = 0; iLambdaBar_ME < Lbar_vector_ME.size(); iLambdaBar_ME++)
    {
      double L_Lbar_pairThetaStar = LpairThetaStar(L_vector_ME.at(iLambda_ME), p_vector_ME.at(iLambda_ME), Lbar_vector_ME.at(iLambdaBar_ME), pbar_vector_ME.at(iLambdaBar_ME));

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
      double Lbar_Lbar_pairThetaStar = LpairThetaStar(Lbar_vector_ME.at(iLambdaBar_ME_1), pbar_vector_ME.at(iLambdaBar_ME_1), Lbar_vector_ME.at(iLambdaBar_ME_2), pbar_vector_ME.at(iLambdaBar_ME_2));

      L0bar_L0bar_cosThetaProdPlane_ME_hist->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar));
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_1)][Lbar_pT_bin_vector_ME.at(iLambdaBar_ME_2)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar));
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[Lbar_eta_bin_vector_ME.at(iLambdaBar_ME_1)][Lbar_eta_bin_vector_ME.at(iLambdaBar_ME_2)]->Fill(TMath::Cos(Lbar_Lbar_pairThetaStar));
    }
  }

  //________________________________________________________________________________________________________

  TCanvas *L0_L0bar_cosThetaProdPlane_stat_err_can = new TCanvas("L0_L0bar_cosThetaProdPlane_stat_err_can", "L0_L0bar_cosThetaProdPlane_stat_err_can", 1200, 1000);

  L0_L0bar_cosThetaProdPlane_stat_err_can->cd();

  for(unsigned int iBin = 0; iBin < 10; iBin++)
  {
    L0_L0bar_cosThetaProdPlane_stat_err->SetBinContent(iBin+1, L0_L0bar_cosThetaProdPlane_US_hist->GetBinError(iBin+1)/L0_L0bar_cosThetaProdPlane_US_hist->GetBinContent(iBin+1)*100 );
    L0_L0bar_cosThetaProdPlane_stat_err->SetBinError(iBin+1, 0 );
  }
  L0_L0bar_cosThetaProdPlane_stat_err->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_stat_err->GetXaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_stat_err->GetYaxis()->SetTitle("Relative error (%)");
  L0_L0bar_cosThetaProdPlane_stat_err->GetYaxis()->CenterTitle();
  L0_L0bar_cosThetaProdPlane_stat_err->SetMinimum(0);
  L0_L0bar_cosThetaProdPlane_stat_err->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_stat_err->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_stat_err->SetMarkerColor(kRed);
  L0_L0bar_cosThetaProdPlane_stat_err->SetLineColor(kRed);
  L0_L0bar_cosThetaProdPlane_stat_err->Draw("p e");

  TPaveText *L0_L0bar_text = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
  L0_L0bar_text->SetTextFont(42);
  //L0_L0bar_text->AddText("STAR Internal");
  //L0_L0bar_text->AddText("STAR preliminary");
  //((TText*)L0_L0bar_text->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0bar_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  if(trigger == 0 || trigger == 2) L0_L0bar_text->AddText("Minimum bias");
  if(trigger == 1) L0_L0bar_text->AddText("JP2");
  L0_L0bar_text->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
  L0_L0bar_text->AddText("|#it{y}| < 1");
  L0_L0bar_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L0_L0bar_text->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text->Draw("same");


  L0_L0bar_cosThetaProdPlane_stat_err_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0bar_cosThetaProdPlane_stat_err.png");
  //________________________________________________________________


  TCanvas *L0_L0bar_cosThetaProdPlane_no_corr_can = new TCanvas("L0_L0bar_cosThetaProdPlane_no_corr_can", "L0_L0bar_cosThetaProdPlane_no_corr_can", 1200, 1000);

  L0_L0bar_cosThetaProdPlane_no_corr_can->cd();

  L0_L0bar_cosThetaProdPlane_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_US_hist->GetXaxis()->CenterTitle();
  //L0_L0bar_cosThetaProdPlane_US_hist->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
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

  L0_L0bar_cosThetaProdPlane_LS_hist->SetMarkerSize(1.5);
  L0_L0bar_cosThetaProdPlane_LS_hist->SetMarkerStyle(20);
  L0_L0bar_cosThetaProdPlane_LS_hist->SetMarkerColor(kBlue);
  double nLLbar_back = L0_L0bar_cosThetaProdPlane_LS_hist->Integral();
  L0_L0bar_cosThetaProdPlane_LS_hist->Sumw2();
  L0_L0bar_cosThetaProdPlane_LS_hist->Scale(1./L0_L0bar_cosThetaProdPlane_LS_hist->GetXaxis()->GetBinWidth(1));
  //L0_L0bar_cosThetaProdPlane_LS_hist->Divide(L0_L0bar_cosThetaProdPlane_eff);
  //L0_L0bar_cosThetaProdPlane_LS_hist->Draw("p e same");

  L0_L0bar_cosThetaProdPlane_ME_hist->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0bar_cosThetaProdPlane_ME_hist->GetXaxis()->CenterTitle();
  //L0_L0bar_cosThetaProdPlane_ME_hist->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
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
  L0_L0bar_cosThetaProdPlane_ME_hist->Draw("same p e");


  TF1 *fitL0_L0bar_US_ThetaStar_no_corr = new TF1("fitL0_L0bar_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
  fitL0_L0bar_US_ThetaStar_no_corr->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0_L0bar_cosThetaProdPlane_US_hist->Fit(fitL0_L0bar_US_ThetaStar_no_corr, "s i 0 r");

  float P_L0_L0bar_no_corr = fitL0_L0bar_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
  float P_L0_L0bar_no_corr_err = fitL0_L0bar_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

  fitL0_L0bar_US_ThetaStar_no_corr->SetLineColor(1);
  fitL0_L0bar_US_ThetaStar_no_corr->Draw("same");

  TLegend *L0_L0bar_leg = new TLegend(0.15, 0.2, 0.39, 0.49);
  L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_US_hist, "Unlike-sign p#pi");
  //L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_LS_hist, "Like-sign p#pi");
  L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_LS_hist, "Mixed event");
  L0_L0bar_leg->AddEntry(fitL0_L0bar_US_ThetaStar_no_corr, "Linear fit to US");
  L0_L0bar_leg->SetBorderSize(0);
  L0_L0bar_leg->SetFillColorAlpha(0, 0.01);
  L0_L0bar_leg->Draw("same");

  TPaveText *L0_L0bar_text_no_corr = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
  L0_L0bar_text_no_corr->SetTextFont(42);
  //L0_L0bar_text_no_corr->AddText("STAR Internal");
  //L0_L0bar_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0bar_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  if(trigger == 0 || trigger == 2) L0_L0bar_text_no_corr->AddText("Minimum bias, no correction");
  if(trigger == 1) L0_L0bar_text_no_corr->AddText("JP2, no correction");
  L0_L0bar_text_no_corr->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
  L0_L0bar_text_no_corr->AddText("|#it{y}| < 1");
  L0_L0bar_text_no_corr->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L0_L0bar_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar_no_corr, P_L0_L0bar_no_corr_err));
  L0_L0bar_text_no_corr->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_no_corr->Draw("same");

  L0_L0bar_cosThetaProdPlane_no_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0bar_cosThetaProdPlane_no_corr.png");
  //_____________________________________________________________________________________

  TCanvas *L0_L0bar_cosThetaProdPlane_can = new TCanvas("L0_L0bar_cosThetaProdPlane_can", "L0_L0bar_cosThetaProdPlane_can", 1200, 1000);

  L0_L0bar_cosThetaProdPlane_can->cd();

  L0_L0bar_cosThetaProdPlane_US_hist->Divide(L0_L0bar_cosThetaProdPlane_eff);
  //L0_L0bar_cosThetaProdPlane_US_hist->Scale(1./L0_L0bar_cosThetaProdPlane_US_hist->Integral());
  L0_L0bar_cosThetaProdPlane_US_hist->SetMinimum(0);
  L0_L0bar_cosThetaProdPlane_US_hist->Draw("p e");

  L0_L0bar_cosThetaProdPlane_ME_hist->Divide(L0_L0bar_cosThetaProdPlane_ME_eff);
  //L0_L0bar_cosThetaProdPlane_ME_hist->Scale(1./L0_L0bar_cosThetaProdPlane_ME_hist->Integral());
  L0_L0bar_cosThetaProdPlane_ME_hist->Draw("same p e");

  TF1 *fitL0_L0bar_US_ThetaStar = new TF1("fitL0_L0bar_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fitL0_L0bar_US_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0_L0bar_cosThetaProdPlane_US_hist->Fit(fitL0_L0bar_US_ThetaStar, "s i 0 r");

  float P_L0_L0bar = fitL0_L0bar_US_ThetaStar->GetParameter(1)/(L0_alpha*L0bar_alpha);
  float P_L0_L0bar_err = fitL0_L0bar_US_ThetaStar->GetParError(1)/(L0_alpha*L0bar_alpha);

  fitL0_L0bar_US_ThetaStar->SetLineColor(1);
  fitL0_L0bar_US_ThetaStar->Draw("same");


  L0_L0bar_text->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar, P_L0_L0bar_err));
  L0_L0bar_text->Draw("same");

  L0_L0bar_leg->Draw("same");

  L0_L0bar_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0bar_cosThetaProdPlane.png");
  //____________________________________________________________________


  TCanvas *L0_L0_cosThetaProdPlane_stat_err_can = new TCanvas("L0_L0_cosThetaProdPlane_stat_err_can", "L0_L0_cosThetaProdPlane_stat_err_can", 1200, 1000);

  L0_L0_cosThetaProdPlane_stat_err_can->cd();

  for(unsigned int iBin = 0; iBin < 10; iBin++)
  {
    L0_L0_cosThetaProdPlane_stat_err->SetBinContent(iBin+1, L0_L0_cosThetaProdPlane_US_hist->GetBinError(iBin+1)/L0_L0_cosThetaProdPlane_US_hist->GetBinContent(iBin+1)*100 );
    L0_L0_cosThetaProdPlane_stat_err->SetBinError(iBin+1, 0 );
  }
  L0_L0_cosThetaProdPlane_stat_err->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_stat_err->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_stat_err->GetYaxis()->SetTitle("Relative error (%)");
  L0_L0_cosThetaProdPlane_stat_err->GetYaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_stat_err->SetMinimum(0);
  L0_L0_cosThetaProdPlane_stat_err->SetMarkerSize(1.5);
  L0_L0_cosThetaProdPlane_stat_err->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_stat_err->SetMarkerColor(kRed);
  L0_L0_cosThetaProdPlane_stat_err->SetLineColor(kRed);
  L0_L0_cosThetaProdPlane_stat_err->Draw("p e");


  TPaveText *L0_L0_text = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
  L0_L0_text->SetTextFont(42);
  //L0_L0_text->AddText("STAR Internal");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  if(trigger == 0 || trigger == 2) L0_L0_text->AddText("Minimum bias");
  if(trigger == 1) L0_L0_text->AddText("JP2");
  L0_L0_text->AddText("#Lambda^{0}-#Lambda^{0}");
  L0_L0_text->AddText("|#it{y}| < 1");
  L0_L0_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L0_L0_text->SetFillColorAlpha(0, 0.01);
  L0_L0_text->Draw("same");


  L0_L0_cosThetaProdPlane_stat_err_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0_cosThetaProdPlane_stat_err.png");
  //________________________________


  TCanvas *L0_L0_cosThetaProdPlane_no_corr_can = new TCanvas("L0_L0_cosThetaProdPlane_no_corr_can", "L0_L0_cosThetaProdPlane_no_corr_can", 1200, 1000);

  L0_L0_cosThetaProdPlane_no_corr_can->cd();

  L0_L0_cosThetaProdPlane_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_US_hist->GetXaxis()->CenterTitle();
  //L0_L0_cosThetaProdPlane_US_hist->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
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

  L0_L0_cosThetaProdPlane_LS_hist->SetMarkerSize(1.5);
  L0_L0_cosThetaProdPlane_LS_hist->SetMarkerStyle(20);
  L0_L0_cosThetaProdPlane_LS_hist->SetMarkerColor(kBlue);
  double nLL_back = L0_L0_cosThetaProdPlane_LS_hist->Integral();
  L0_L0_cosThetaProdPlane_LS_hist->Sumw2();
  L0_L0_cosThetaProdPlane_LS_hist->Scale(1./L0_L0_cosThetaProdPlane_LS_hist->GetXaxis()->GetBinWidth(1));
  L0_L0_cosThetaProdPlane_LS_hist->Divide(L0_L0_cosThetaProdPlane_eff);
  //L0_L0_cosThetaProdPlane_LS_hist->Draw("p e same");

  L0_L0_cosThetaProdPlane_ME_hist->GetXaxis()->SetTitle("cos(#theta*)");
  L0_L0_cosThetaProdPlane_ME_hist->GetXaxis()->CenterTitle();
  L0_L0_cosThetaProdPlane_ME_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  //L0_L0_cosThetaProdPlane_ME_hist->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
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
  L0_L0_cosThetaProdPlane_ME_hist->Draw("same p e");


  TF1 *fitL0_L0_US_ThetaStar_no_corr = new TF1("fitL0_L0_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
  fitL0_L0_US_ThetaStar_no_corr->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0_L0_cosThetaProdPlane_US_hist->Fit(fitL0_L0_US_ThetaStar_no_corr, "s i 0 r");

  float P_L0_L0_no_corr = fitL0_L0_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_L0_L0_no_corr_err = fitL0_L0_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0_alpha);

  fitL0_L0_US_ThetaStar_no_corr->SetLineColor(1);
  fitL0_L0_US_ThetaStar_no_corr->Draw("same");

  TLegend *L0_L0_leg = new TLegend(0.15, 0.2, 0.39, 0.49);
  L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_US_hist, "Unlike-sign p#pi");
  //L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_LS_hist, "Like-sign p#pi");
  L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_LS_hist, "Mixed event");
  L0_L0_leg->AddEntry(fitL0_L0_US_ThetaStar_no_corr, "Linear fit to US");
  L0_L0_leg->SetBorderSize(0);
  L0_L0_leg->SetFillColorAlpha(0, 0.01);
  L0_L0_leg->Draw("same");

  TPaveText *L0_L0_text_no_corr = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
  L0_L0_text_no_corr->SetTextFont(42);
  //L0_L0_text_no_corr->AddText("STAR Internal");
  //L0_L0_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  if(trigger == 0 || trigger == 2) L0_L0_text_no_corr->AddText("Minimum bias, no correction");
  if(trigger == 1) L0_L0_text_no_corr->AddText("JP2, no correction");
  L0_L0_text_no_corr->AddText("#Lambda^{0}-#Lambda^{0}");
  L0_L0_text_no_corr->AddText("|#it{y}| < 1");
  L0_L0_text_no_corr->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L0_L0_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0_no_corr, P_L0_L0_no_corr_err));
  L0_L0_text_no_corr->SetFillColorAlpha(0, 0.01);
  L0_L0_text_no_corr->Draw("same");

  L0_L0_cosThetaProdPlane_no_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0_cosThetaProdPlane_no_corr.png");
  //________________________________________________________________________

  TCanvas *L0_L0_cosThetaProdPlane_can = new TCanvas("L0_L0_cosThetaProdPlane_can", "L0_L0_cosThetaProdPlane_can", 1200, 1000);

  L0_L0_cosThetaProdPlane_can->cd();

  L0_L0_cosThetaProdPlane_US_hist->Divide(L0_L0_cosThetaProdPlane_eff);
  //L0_L0_cosThetaProdPlane_US_hist->Scale(1./L0_L0_cosThetaProdPlane_US_hist->Integral());
  L0_L0_cosThetaProdPlane_US_hist->SetMinimum(0);
  L0_L0_cosThetaProdPlane_US_hist->Draw("p e");

  L0_L0_cosThetaProdPlane_ME_hist->Divide(L0_L0_cosThetaProdPlane_ME_eff);
  //L0_L0_cosThetaProdPlane_ME_hist->Scale(1./L0_L0_cosThetaProdPlane_ME_hist->Integral());
  L0_L0_cosThetaProdPlane_ME_hist->Draw("same p e");

  TF1 *fitL0_L0_US_ThetaStar = new TF1("fitL0_L0_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fitL0_L0_US_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0_L0_cosThetaProdPlane_US_hist->Fit(fitL0_L0_US_ThetaStar, "s i 0 r");

  float P_L0_L0 = fitL0_L0_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
  float P_L0_L0_err = fitL0_L0_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

  fitL0_L0_US_ThetaStar->SetLineColor(1);
  fitL0_L0_US_ThetaStar->Draw("same");

  L0_L0_text->AddText(Form("P = %.2f #pm %.2f", P_L0_L0, P_L0_L0_err));
  L0_L0_text->Draw("same");

  L0_L0_leg->Draw("same");

  L0_L0_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0_cosThetaProdPlane.png");
  //____________________________________________________________________

  TCanvas *L0bar_L0bar_cosThetaProdPlane_stat_err_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_stat_err_can", "L0bar_L0bar_cosThetaProdPlane_stat_err_can", 1200, 1000);

  L0bar_L0bar_cosThetaProdPlane_stat_err_can->cd();

  for(unsigned int iBin = 0; iBin < 10; iBin++)
  {
    L0bar_L0bar_cosThetaProdPlane_stat_err->SetBinContent(iBin+1, L0bar_L0bar_cosThetaProdPlane_US_hist->GetBinError(iBin+1)/L0bar_L0bar_cosThetaProdPlane_US_hist->GetBinContent(iBin+1)*100 );
    L0bar_L0bar_cosThetaProdPlane_stat_err->SetBinError(iBin+1, 0 );
  }
  L0bar_L0bar_cosThetaProdPlane_stat_err->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_stat_err->GetXaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_stat_err->GetYaxis()->SetTitle("Relative error (%)");
  L0bar_L0bar_cosThetaProdPlane_stat_err->GetYaxis()->CenterTitle();
  L0bar_L0bar_cosThetaProdPlane_stat_err->SetMinimum(0);
  L0bar_L0bar_cosThetaProdPlane_stat_err->SetMarkerSize(1.5);
  L0bar_L0bar_cosThetaProdPlane_stat_err->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_stat_err->SetMarkerColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_stat_err->SetLineColor(kRed);
  L0bar_L0bar_cosThetaProdPlane_stat_err->Draw("p e");

  TPaveText *L0bar_L0bar_text = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
  L0bar_L0bar_text->SetTextFont(42);
  //L0bar_L0bar_text->AddText("STAR Internal");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  L0bar_L0bar_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  if(trigger == 0 || trigger == 2) L0bar_L0bar_text->AddText("Minimum bias");
  if(trigger == 1) L0bar_L0bar_text->AddText("JP2");
  L0bar_L0bar_text->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
  L0bar_L0bar_text->AddText("|#it{y}| < 1");
  L0bar_L0bar_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L0bar_L0bar_text->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text->Draw("same");


  L0bar_L0bar_cosThetaProdPlane_stat_err_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_L0bar_cosThetaProdPlane_stat_err.png");
  //________________________________


  TCanvas *L0bar_L0bar_cosThetaProdPlane_no_corr_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_no_corr_can", "L0bar_L0bar_cosThetaProdPlane_no_corr_can", 1200, 1000);

  L0bar_L0bar_cosThetaProdPlane_no_corr_can->cd();

  L0bar_L0bar_cosThetaProdPlane_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_US_hist->GetXaxis()->CenterTitle();
  //L0bar_L0bar_cosThetaProdPlane_US_hist->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
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

  L0bar_L0bar_cosThetaProdPlane_LS_hist->SetMarkerSize(1.5);
  L0bar_L0bar_cosThetaProdPlane_LS_hist->SetMarkerStyle(20);
  L0bar_L0bar_cosThetaProdPlane_LS_hist->SetMarkerColor(kBlue);
  double nLbarLbar_back = L0bar_L0bar_cosThetaProdPlane_LS_hist->Integral();
  L0bar_L0bar_cosThetaProdPlane_LS_hist->Sumw2();
  L0bar_L0bar_cosThetaProdPlane_LS_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_LS_hist->GetXaxis()->GetBinWidth(1));
  //L0bar_L0bar_cosThetaProdPlane_LS_hist->Divide(L0bar_L0bar_cosThetaProdPlane_eff);
  //L0bar_L0bar_cosThetaProdPlane_LS_hist->Draw("p e same");

  L0bar_L0bar_cosThetaProdPlane_ME_hist->GetXaxis()->SetTitle("cos(#theta*)");
  L0bar_L0bar_cosThetaProdPlane_ME_hist->GetXaxis()->CenterTitle();
  //L0bar_L0bar_cosThetaProdPlane_ME_hist->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
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
  L0bar_L0bar_cosThetaProdPlane_ME_hist->Draw("same p e");


  TF1 *fitL0bar_L0bar_US_ThetaStar_no_corr = new TF1("fitL0bar_L0bar_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
  fitL0bar_L0bar_US_ThetaStar_no_corr->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0bar_L0bar_cosThetaProdPlane_US_hist->Fit(fitL0bar_L0bar_US_ThetaStar_no_corr, "s i 0 r");

  float P_L0bar_L0bar_no_corr = fitL0bar_L0bar_US_ThetaStar_no_corr->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
  float P_L0bar_L0bar_no_corr_err = fitL0bar_L0bar_US_ThetaStar_no_corr->GetParError(1)/(L0bar_alpha*L0bar_alpha);

  fitL0bar_L0bar_US_ThetaStar_no_corr->SetLineColor(1);
  fitL0bar_L0bar_US_ThetaStar_no_corr->Draw("same");


  TLegend *L0bar_L0bar_leg = new TLegend(0.15, 0.2, 0.39, 0.49);
  L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_US_hist, "Unlike-sign p#pi");
  //L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_LS_hist, "Like-sign p#pi");
  L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_LS_hist, "Mixed event");
  L0bar_L0bar_leg->AddEntry(fitL0bar_L0bar_US_ThetaStar_no_corr, "Linear fit to US");
  L0bar_L0bar_leg->SetBorderSize(0);
  L0bar_L0bar_leg->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_leg->Draw("same");


  TPaveText *L0bar_L0bar_text_no_corr = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
  L0bar_L0bar_text_no_corr->SetTextFont(42);
  //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
  //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0bar_L0bar_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  if(trigger == 0 || trigger == 2) L0bar_L0bar_text_no_corr->AddText("Minimum bias, no correction");
  if(trigger == 1) L0bar_L0bar_text_no_corr->AddText("JP2, no correction");
  L0bar_L0bar_text_no_corr->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
  L0bar_L0bar_text_no_corr->AddText("|#it{y}| < 1");
  L0bar_L0bar_text_no_corr->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L0bar_L0bar_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar_no_corr, P_L0bar_L0bar_no_corr_err));
  L0bar_L0bar_text_no_corr->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_no_corr->Draw("same");


  L0bar_L0bar_cosThetaProdPlane_no_corr_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_L0bar_cosThetaProdPlane_no_corr.png");
  //_______________________________________________________________________________________

  TCanvas *L0bar_L0bar_cosThetaProdPlane_can = new TCanvas("L0bar_L0bar_cosThetaProdPlane_can", "L0bar_L0bar_cosThetaProdPlane_can", 1200, 1000);

  L0bar_L0bar_cosThetaProdPlane_can->cd();

  L0bar_L0bar_cosThetaProdPlane_US_hist->Divide(L0bar_L0bar_cosThetaProdPlane_eff);
  //L0bar_L0bar_cosThetaProdPlane_US_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_US_hist->Integral());
  L0bar_L0bar_cosThetaProdPlane_US_hist->SetMinimum(0);
  L0bar_L0bar_cosThetaProdPlane_US_hist->Draw("p e");

  L0bar_L0bar_cosThetaProdPlane_ME_hist->Divide(L0bar_L0bar_cosThetaProdPlane_ME_eff);
  //L0bar_L0bar_cosThetaProdPlane_ME_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_ME_hist->Integral());
  L0bar_L0bar_cosThetaProdPlane_ME_hist->Draw("same p e");

    TF1 *fitL0bar_L0bar_US_ThetaStar = new TF1("fitL0bar_L0bar_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
  fitL0bar_L0bar_US_ThetaStar->SetParameters(100, 0.5);

  //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
  L0bar_L0bar_cosThetaProdPlane_US_hist->Fit(fitL0bar_L0bar_US_ThetaStar, "s i 0 r");

  float P_L0bar_L0bar = fitL0bar_L0bar_US_ThetaStar->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
  float P_L0bar_L0bar_err = fitL0bar_L0bar_US_ThetaStar->GetParError(1)/(L0bar_alpha*L0bar_alpha);

  fitL0bar_L0bar_US_ThetaStar->SetLineColor(1);
  fitL0bar_L0bar_US_ThetaStar->Draw("same");

  L0bar_L0bar_text->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar, P_L0bar_L0bar_err));
  L0bar_L0bar_text->Draw("same");


  L0bar_L0bar_leg->Draw("same");


  L0bar_L0bar_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_L0bar_cosThetaProdPlane.png");
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

      TCanvas *L0_L0bar_cosThetaProdPlane_pT_stat_err_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_stat_err_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_stat_err_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0_L0bar_cosThetaProdPlane_pT_stat_err_can->cd();

      for(unsigned int iBin = 0; iBin < 10; iBin++)
      {
        L0_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetBinContent(iBin+1, L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetBinError(iBin+1)/L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetBinContent(iBin+1)*100 );
        L0_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetBinError(iBin+1, 0 );
      }
      L0_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->GetYaxis()->SetTitle("Relative error (%)");
      L0_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetMarkerSize(1.5);
      L0_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->Draw("p e");

      TPaveText *L0_L0bar_pT_text = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
      L0_L0bar_pT_text->SetTextFont(42);
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0bar_pT_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      if(trigger == 0 || trigger == 2) L0_L0bar_pT_text->AddText("Minimum bias");
      if(trigger == 1) L0_L0bar_pT_text->AddText("JP2");
      L0_L0bar_pT_text->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_pT_text->AddText("|#it{y}| < 1");
      L0_L0bar_pT_text->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0bar_pT_text->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0bar_pT_text->SetFillColorAlpha(0, 0.01);
      L0_L0bar_pT_text->Draw("same");


      L0_L0bar_cosThetaProdPlane_pT_stat_err_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0bar_cosThetaProdPlane_stat_err_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //________________________________________________________________

      TCanvas *L0_L0bar_cosThetaProdPlane_pT_no_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0_L0bar_cosThetaProdPlane_pT_no_corr_can->cd();

      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      //L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
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
      //L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetLineColor(kBlue);
      //double nLLbar = L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral();
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Sumw2();
      //L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2]);
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(nLLbar_pT_back[pTbin1][pTbin2]/L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral());
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Draw("same p e");


      TF1 *fitL0_L0bar_pT_US_ThetaStar_no_corr = new TF1("fitL0_L0bar_pT_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_pT_US_ThetaStar_no_corr->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitL0_L0bar_pT_US_ThetaStar_no_corr, "s i 0 r");

      float P_L0_L0bar_pT_no_corr = fitL0_L0bar_pT_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
      float P_L0_L0bar_pT_no_corr_err = fitL0_L0bar_pT_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

      fitL0_L0bar_pT_US_ThetaStar_no_corr->SetLineColor(1);
      fitL0_L0bar_pT_US_ThetaStar_no_corr->Draw("same");

      TPaveText *L0_L0bar_pT_text_no_corr = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
      L0_L0bar_pT_text_no_corr->SetTextFont(42);
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0bar_pT_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      if(trigger == 0 || trigger == 2) L0_L0bar_pT_text_no_corr->AddText("Minimum bias, no correction");
      if(trigger == 1) L0_L0bar_pT_text_no_corr->AddText("JP2, no correction");
      L0_L0bar_pT_text_no_corr->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      //L0_L0bar_pT_text_no_corr->AddText("|#it{y}| < 0.2");
      L0_L0bar_pT_text_no_corr->AddText("|#it{y}| < 1");
      L0_L0bar_pT_text_no_corr->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0bar_pT_text_no_corr->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0bar_pT_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar_pT_no_corr, P_L0_L0bar_pT_no_corr_err));
      L0_L0bar_pT_text_no_corr->SetFillColorAlpha(0, 0.01);
      L0_L0bar_pT_text_no_corr->Draw("same");

      L0_L0bar_leg->Draw("same");


      L0_L0bar_cosThetaProdPlane_pT_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0bar_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //_______________________________________________



      TCanvas *L0_L0bar_cosThetaProdPlane_pT_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0_L0bar_cosThetaProdPlane_pT_can->cd();

      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      //L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral());
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Draw("p e");

      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Divide(L0_L0bar_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2]);
      //L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(1./L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral());
      L0_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Draw("same p e");

       TF1 *fitL0_L0bar_pT_US_ThetaStar = new TF1("fitL0_L0bar_pT_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_pT_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitL0_L0bar_pT_US_ThetaStar, "s i 0 r");

      float P_L0_L0bar_pT = fitL0_L0bar_pT_US_ThetaStar->GetParameter(1)/(L0_alpha*L0bar_alpha);
      float P_L0_L0bar_pT_err = fitL0_L0bar_pT_US_ThetaStar->GetParError(1)/(L0_alpha*L0bar_alpha);

      fitL0_L0bar_pT_US_ThetaStar->SetLineColor(1);
      fitL0_L0bar_pT_US_ThetaStar->Draw("same");


      L0_L0bar_pT_text->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar_pT, P_L0_L0bar_pT_err));
      L0_L0bar_pT_text->Draw("same");

      L0_L0bar_leg->Draw("same");

      L0_L0bar_cosThetaProdPlane_pT_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //____________________________________________________________________________________________________


      TCanvas *L0_L0_cosThetaProdPlane_pT_stat_err_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_stat_err_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_stat_err_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0_L0_cosThetaProdPlane_pT_stat_err_can->cd();

      for(unsigned int iBin = 0; iBin < 10; iBin++)
      {
        L0_L0_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetBinContent(iBin+1, L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetBinError(iBin+1)/L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetBinContent(iBin+1)*100 );
        L0_L0_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetBinError(iBin+1, 0 );
      }
      L0_L0_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->GetYaxis()->SetTitle("Relative error (%)");
      L0_L0_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetMarkerSize(1.5);
      L0_L0_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_L0_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_L0_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->Draw("p e");


      TPaveText *L0_L0_pT_text = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
      L0_L0_pT_text->SetTextFont(42);
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_pT_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      if(trigger == 0 || trigger == 2) L0_L0_pT_text->AddText("Minimum bias");
      if(trigger == 1) L0_L0_pT_text->AddText("JP2");
      L0_L0_pT_text->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_pT_text->AddText("|#it{y}| < 1");
      L0_L0_pT_text->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0_pT_text->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0_pT_text->SetFillColorAlpha(0, 0.01);
      L0_L0_pT_text->Draw("same");


      L0_L0_cosThetaProdPlane_pT_stat_err_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0_cosThetaProdPlane_stat_err_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //________________________________________________________________


      TCanvas *L0_L0_cosThetaProdPlane_pT_no_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0_L0_cosThetaProdPlane_pT_no_corr_can->cd();

      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      //L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
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

      L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      nLL_pT_back[pTbin1][pTbin2] = L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral();
      L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Sumw2();
      //L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      //L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      //L0_L0_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Draw("p e same");

      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      //L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetLineColor(kBlue);
      //double nLLbar = L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral();
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Sumw2();
      //L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2]);
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(nLL_pT_back[pTbin1][pTbin2]/L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral());
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Draw("same p e");

      TF1 *fitL0_L0_pT_US_ThetaStar_no_corr = new TF1("fitL0_L0_pT_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_pT_US_ThetaStar_no_corr->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitL0_L0_pT_US_ThetaStar_no_corr, "s i 0 r");

      float P_L0_L0_pT_no_corr = fitL0_L0_pT_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_L0_L0_pT_no_corr_err = fitL0_L0_pT_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0_alpha);

      fitL0_L0_pT_US_ThetaStar_no_corr->SetLineColor(1);
      fitL0_L0_pT_US_ThetaStar_no_corr->Draw("same");

      TPaveText *L0_L0_pT_text_no_corr = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
      L0_L0_pT_text_no_corr->SetTextFont(42);
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_pT_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      if(trigger == 0 || trigger == 2) L0_L0_pT_text_no_corr->AddText("Minimum bias, no correction");
      if(trigger == 1) L0_L0_pT_text_no_corr->AddText("JP2, no correction");
      L0_L0_pT_text_no_corr->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_pT_text_no_corr->AddText("|#it{y}| < 1");
      L0_L0_pT_text_no_corr->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0_pT_text_no_corr->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0_pT_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0_pT_no_corr, P_L0_L0_pT_no_corr_err));
      L0_L0_pT_text_no_corr->SetFillColorAlpha(0, 0.01);
      L0_L0_pT_text_no_corr->Draw("same");

      L0_L0_leg->Draw("same");

      L0_L0_cosThetaProdPlane_pT_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //_______________________________________________________________________

      TCanvas *L0_L0_cosThetaProdPlane_pT_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0_L0_cosThetaProdPlane_pT_can->cd();


      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      //L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral());
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Draw("p e");

      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Divide(L0_L0_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2]);
      //L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(1./L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral());
      L0_L0_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Draw("same p e");


      TF1 *fitL0_L0_pT_US_ThetaStar = new TF1("fitL0_L0_pT_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_pT_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitL0_L0_pT_US_ThetaStar, "s i 0 r");

      float P_L0_L0_pT = fitL0_L0_pT_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_L0_L0_pT_err = fitL0_L0_pT_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

      fitL0_L0_pT_US_ThetaStar->SetLineColor(1);
      fitL0_L0_pT_US_ThetaStar->Draw("same");


      L0_L0_pT_text->AddText(Form("P = %.2f #pm %.2f", P_L0_L0_pT, P_L0_L0_pT_err));
      L0_L0_pT_text->Draw("same");

      L0_L0_leg->Draw("same");


      L0_L0_cosThetaProdPlane_pT_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0_cosThetaProdPlane_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //____________________________________________________________________________________________________

      TCanvas *L0bar_L0bar_cosThetaProdPlane_pT_stat_err_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_stat_err_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_stat_err_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0bar_L0bar_cosThetaProdPlane_pT_stat_err_can->cd();

      for(unsigned int iBin = 0; iBin < 10; iBin++)
      {
        L0bar_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetBinContent(iBin+1, L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetBinError(iBin+1)/L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetBinContent(iBin+1)*100 );
        L0bar_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetBinError(iBin+1, 0 );
      }
      L0bar_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->GetYaxis()->SetTitle("Relative error (%)");
      L0bar_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetMarkerSize(1.5);
      L0bar_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar_L0bar_cosThetaProdPlane_pT_stat_err[pTbin1][pTbin2]->Draw("p e");

      TPaveText *L0bar_L0bar_pT_text = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
      L0bar_L0bar_pT_text->SetTextFont(42);
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_pT_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      if(trigger == 0 || trigger == 2) L0bar_L0bar_pT_text->AddText("Minimum bias");
      if(trigger == 1) L0bar_L0bar_pT_text->AddText("JP2");
      L0bar_L0bar_pT_text->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_pT_text->AddText("|#it{y}| < 1");
      L0bar_L0bar_pT_text->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0bar_L0bar_pT_text->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0bar_L0bar_pT_text->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_pT_text->Draw("same");


      L0bar_L0bar_cosThetaProdPlane_pT_stat_err_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_L0bar_cosThetaProdPlane_stat_err_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //________________________________________________________________


      TCanvas *L0bar_L0bar_cosThetaProdPlane_pT_no_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0bar_L0bar_cosThetaProdPlane_pT_no_corr_can->cd();

      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      //L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
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

      L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      nLbarLbar_pT_back[pTbin1][pTbin2] = L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Integral();
      L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Sumw2();
      //L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      //L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      //L0bar_L0bar_cosThetaProdPlane_pT_LS_hist[pTbin1][pTbin2]->Draw("p e same");

      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      //L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetLineColor(kBlue);
      //double nLLbar = L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral();
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Sumw2();
      //L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2]);
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(nLbarLbar_pT_back[pTbin1][pTbin2]/L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->GetXaxis()->GetBinWidth(1));
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Draw("same p e");



      TF1 *fitL0bar_L0bar_pT_US_ThetaStar_no_corr = new TF1("fitL0bar_L0bar_pT_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_pT_US_ThetaStar_no_corr->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitL0bar_L0bar_pT_US_ThetaStar_no_corr, "s i 0 r");

      float P_L0bar_L0bar_pT_no_corr = fitL0bar_L0bar_pT_US_ThetaStar_no_corr->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
      float P_L0bar_L0bar_pT_no_corr_err = fitL0bar_L0bar_pT_US_ThetaStar_no_corr->GetParError(1)/(L0bar_alpha*L0bar_alpha);

      fitL0bar_L0bar_pT_US_ThetaStar_no_corr->SetLineColor(1);
      fitL0bar_L0bar_pT_US_ThetaStar_no_corr->Draw("same");

      TPaveText *L0bar_L0bar_pT_text_no_corr = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
      L0bar_L0bar_pT_text_no_corr->SetTextFont(42);
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_pT_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      if(trigger == 0 || trigger == 2) L0bar_L0bar_pT_text_no_corr->AddText("Minimum bias, no correction");
      if(trigger == 1) L0bar_L0bar_pT_text_no_corr->AddText("JP2, no correction");
      L0bar_L0bar_pT_text_no_corr->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_pT_text_no_corr->AddText("|#it{y}| < 1");
      L0bar_L0bar_pT_text_no_corr->AddText(Form("%.2f < p_{T}^{1} < %.2f", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0bar_L0bar_pT_text_no_corr->AddText(Form("%.2f < p_{T}^{2} < %.2f", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0bar_L0bar_pT_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar_pT_no_corr, P_L0bar_L0bar_pT_no_corr_err));
      L0bar_L0bar_pT_text_no_corr->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_pT_text_no_corr->Draw("same");

      L0bar_L0bar_leg->Draw("same");

      L0bar_L0bar_cosThetaProdPlane_pT_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_L0bar_cosThetaProdPlane_no_corr_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //_____________________________________________________________________________

      TCanvas *L0bar_L0bar_cosThetaProdPlane_pT_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_can", pTbin1, pTbin2), Form("L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i_can", pTbin1, pTbin2), 1200, 1000);

      L0bar_L0bar_cosThetaProdPlane_pT_can->cd();

      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_pT_eff[pTbin1][pTbin2]);
      //L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Draw("p e");

      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Divide(L0bar_L0bar_cosThetaProdPlane_ME_pT_eff[pTbin1][pTbin2]);
      //L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_pT_ME_hist[pTbin1][pTbin2]->Draw("same p e");

      TF1 *fitL0bar_L0bar_pT_US_ThetaStar = new TF1("fitL0bar_L0bar_pT_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_pT_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_pT_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0bar_L0bar_cosThetaProdPlane_pT_US_hist[pTbin1][pTbin2]->Fit(fitL0bar_L0bar_pT_US_ThetaStar, "s i 0 r");

      float P_L0bar_L0bar_pT = fitL0bar_L0bar_pT_US_ThetaStar->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
      float P_L0bar_L0bar_pT_err = fitL0bar_L0bar_pT_US_ThetaStar->GetParError(1)/(L0bar_alpha*L0bar_alpha);

      fitL0bar_L0bar_pT_US_ThetaStar->SetLineColor(1);
      fitL0bar_L0bar_pT_US_ThetaStar->Draw("same");

      L0bar_L0bar_pT_text->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar_pT, P_L0bar_L0bar_pT_err));
      L0bar_L0bar_pT_text->Draw("same");

      L0bar_L0bar_leg->Draw("same");

      L0bar_L0bar_cosThetaProdPlane_pT_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_L0bar_cosThetaProdPlane_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      //____________________________________________________________________________________________________
    }
  }


  for(unsigned int etaBin1 = 0; etaBin1 < nEtaBins; etaBin1++)
  {
    for(unsigned int etaBin2 = 0; etaBin2 < nEtaBins; etaBin2++)
    {
      TCanvas *L0_L0bar_cosThetaProdPlane_eta_no_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i_can", etaBin1, etaBin2), Form("L0_L0bar_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i_can", etaBin1, etaBin2), 1200, 1000);

      L0_L0bar_cosThetaProdPlane_eta_no_corr_can->cd();

      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      //L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
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

      L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerSize(2);
      //double nLLbar_back = L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral();
      L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Sumw2();
      //L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Scale(1./L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      //L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Divide(L0_L0bar_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      //L0_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Draw("p e same");

      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      //L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
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
      L0_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Draw("same p e");


      TF1 *fitL0_L0bar_eta_US_ThetaStar_no_corr = new TF1("fitL0_L0bar_eta_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_eta_US_ThetaStar_no_corr->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_eta_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Fit(fitL0_L0bar_eta_US_ThetaStar_no_corr, "s i 0 r");

      float P_L0_L0bar_eta_no_corr = fitL0_L0bar_eta_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
      float P_L0_L0bar_eta_no_corr_err = fitL0_L0bar_eta_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

      fitL0_L0bar_eta_US_ThetaStar_no_corr->SetLineColor(1);
      fitL0_L0bar_eta_US_ThetaStar_no_corr->Draw("same");

      TPaveText *L0_L0bar_eta_text_no_corr = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0bar_eta_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      if(trigger == 0 || trigger == 2) L0_L0bar_eta_text_no_corr->AddText("Minimum bias, no correction");
      if(trigger == 1) L0_L0bar_eta_text_no_corr->AddText("JP2, no correction");
      L0_L0bar_eta_text_no_corr->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_eta_text_no_corr->AddText("|#it{y}| < 1");
      L0_L0bar_eta_text_no_corr->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0bar_eta_text_no_corr->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0bar_eta_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar_eta_no_corr, P_L0_L0bar_eta_no_corr_err));
      L0_L0bar_eta_text_no_corr->SetFillColorAlpha(0, 0.01);
      L0_L0bar_eta_text_no_corr->Draw("same");

      L0_L0bar_leg->Draw("same");

      L0_L0bar_cosThetaProdPlane_eta_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0bar_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i.png", etaBin1, etaBin2));
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

      TPaveText *L0_L0bar_eta_text = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0bar_eta_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      if(trigger == 0 || trigger == 2) L0_L0bar_eta_text->AddText("Minimum bias");
      if(trigger == 1) L0_L0bar_eta_text->AddText("JP2");
      L0_L0bar_eta_text->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_eta_text->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0bar_eta_text->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0bar_eta_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
      L0_L0bar_eta_text->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar_eta, P_L0_L0bar_eta_err));
      L0_L0bar_eta_text->SetFillColorAlpha(0, 0.01);
      L0_L0bar_eta_text->Draw("same");

      L0_L0bar_leg->Draw("same");

      L0_L0bar_cosThetaProdPlane_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i.png", etaBin1, etaBin2));
      //____________________________________________________________________________________________________



      TCanvas *L0_L0_cosThetaProdPlane_eta_no_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i_can", etaBin1, etaBin2), Form("L0_L0_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i_can", etaBin1, etaBin2), 1200, 1000);

      L0_L0_cosThetaProdPlane_eta_no_corr_can->cd();

      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      //L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
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

      L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerSize(2);
      //double nLL_back = L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral();
      L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Sumw2();
      //L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Scale(1./L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      //L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Divide(L0_L0_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      //L0_L0_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Draw("p e same");

      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      //L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
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
      L0_L0_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Draw("same p e");


      TF1 *fitL0_L0_eta_US_ThetaStar_no_corr = new TF1("fitL0_L0_eta_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_eta_US_ThetaStar_no_corr->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_eta_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0_L0_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Fit(fitL0_L0_eta_US_ThetaStar_no_corr, "s i 0 r");

      float P_L0_L0_eta_no_corr = fitL0_L0_eta_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_L0_L0_eta_no_corr_err = fitL0_L0_eta_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0_alpha);

      fitL0_L0_eta_US_ThetaStar_no_corr->SetLineColor(1);
      fitL0_L0_eta_US_ThetaStar_no_corr->Draw("same");

      TPaveText *L0_L0_eta_text_no_corr = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_eta_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      if(trigger == 0 || trigger == 2) L0_L0_eta_text_no_corr->AddText("Minimum bias, no correction");
      if(trigger == 1) L0_L0_eta_text_no_corr->AddText("JP2, no correction");
      L0_L0_eta_text_no_corr->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_eta_text_no_corr->AddText("|#it{y}| < 1");
      L0_L0_eta_text_no_corr->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0_eta_text_no_corr->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0_eta_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0_eta_no_corr, P_L0_L0_eta_no_corr_err));
      L0_L0_eta_text_no_corr->SetFillColorAlpha(0, 0.01);
      L0_L0_eta_text_no_corr->Draw("same");

      L0_L0_leg->Draw("same");


      L0_L0_cosThetaProdPlane_eta_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i.png", etaBin1, etaBin2));
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

      TPaveText *L0_L0_eta_text = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_eta_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      if(trigger == 0 || trigger == 2) L0_L0_eta_text->AddText("Minimum bias");
      if(trigger == 1) L0_L0_eta_text->AddText("JP2");
      L0_L0_eta_text->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_eta_text->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0_L0_eta_text->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0_L0_eta_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
      L0_L0_eta_text->AddText(Form("P = %.2f #pm %.2f", P_L0_L0_eta, P_L0_L0_eta_err));
      L0_L0_eta_text->SetFillColorAlpha(0, 0.01);
      L0_L0_eta_text->Draw("same");

      L0_L0_leg->Draw("same");

      L0_L0_cosThetaProdPlane_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0_cosThetaProdPlane_eta1_%i_eta2_%i.png", etaBin1, etaBin2));
      //____________________________________________________________________________________________________



      TCanvas *L0bar_L0bar_cosThetaProdPlane__eta_no_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i_can", etaBin1, etaBin2), Form("L0bar_L0bar_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i_can", etaBin1, etaBin2), 1200, 1000);

      L0bar_L0bar_cosThetaProdPlane__eta_no_corr_can->cd();

      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      //L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
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

      L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      //double nLbarLbar_back = L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral();
      L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Sumw2();
      //L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      //L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Divide(L0bar_L0bar_cosThetaProdPlane_eta_eff[etaBin1][etaBin2]);
      //L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Draw("p e same");

      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetXaxis()->CenterTitle();
      //L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetYaxis()->CenterTitle();
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMarkerStyle(20);
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMarkerColor(kBlue);
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetLineColor(kBlue);
      //double nLLbar = L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Integral();
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Sumw2();
      //L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Divide(L0bar_L0bar_cosThetaProdPlane_ME_eta_eff[etaBin1][etaBin2]);
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Scale(L0bar_L0bar_cosThetaProdPlane_eta_LS_hist[etaBin1][etaBin2]->Integral()/L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Integral());
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->GetXaxis()->GetBinWidth(1));
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->SetMinimum(0);
      L0bar_L0bar_cosThetaProdPlane_eta_ME_hist[etaBin1][etaBin2]->Draw("same p e");

      TF1 *fitL0bar_L0bar_eta_US_ThetaStar_no_corr = new TF1("fitL0bar_L0bar_eta_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_eta_US_ThetaStar_no_corr->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_eta_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      L0bar_L0bar_cosThetaProdPlane_eta_US_hist[etaBin1][etaBin2]->Fit(fitL0bar_L0bar_eta_US_ThetaStar_no_corr, "s i 0 r");

      float P_L0bar_L0bar_eta_no_corr = fitL0bar_L0bar_eta_US_ThetaStar_no_corr->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
      float P_L0bar_L0bar_eta_no_corr_err = fitL0bar_L0bar_eta_US_ThetaStar_no_corr->GetParError(1)/(L0bar_alpha*L0bar_alpha);

      fitL0bar_L0bar_eta_US_ThetaStar_no_corr->SetLineColor(1);
      fitL0bar_L0bar_eta_US_ThetaStar_no_corr->Draw("same");

      TPaveText *L0bar_L0bar_eta_text_no_corr = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_eta_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      if(trigger == 0 || trigger == 2) L0bar_L0bar_eta_text_no_corr->AddText("Minimum bias, no correction");
      if(trigger == 1) L0bar_L0bar_eta_text_no_corr->AddText("JP2, no correction");
      L0bar_L0bar_eta_text_no_corr->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_eta_text_no_corr->AddText("|#it{y}| < 1");
      L0bar_L0bar_eta_text_no_corr->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0bar_L0bar_eta_text_no_corr->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0bar_L0bar_eta_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar_eta_no_corr, P_L0bar_L0bar_eta_no_corr_err));
      L0bar_L0bar_eta_text_no_corr->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_eta_text_no_corr->Draw("same");

      L0bar_L0bar_leg->Draw("same");


      L0bar_L0bar_cosThetaProdPlane__eta_no_corr_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_L0bar_cosThetaProdPlane_no_corr_eta1_%i_eta2_%i.png", etaBin1, etaBin2));
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


      TPaveText *L0bar_L0bar_eta_text = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_eta_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      if(trigger == 0 || trigger == 2) L0bar_L0bar_eta_text->AddText("Minimum bias");
      if(trigger == 1) L0bar_L0bar_eta_text->AddText("JP2");
      L0bar_L0bar_eta_text->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_eta_text->AddText(Form("%.2f < #eta^{1} < %.2f", eta_bins[etaBin1], eta_bins[etaBin1+1]));
      L0bar_L0bar_eta_text->AddText(Form("%.2f < #eta^{2} < %.2f", eta_bins[etaBin2], eta_bins[etaBin2+1]));
      L0bar_L0bar_eta_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
      L0bar_L0bar_eta_text->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar_eta, P_L0bar_L0bar_eta_err));
      L0bar_L0bar_eta_text->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_eta_text->Draw("same");

      L0bar_L0bar_leg->Draw("same");


      L0bar_L0bar_cosThetaProdPlane_eta_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_L0bar_cosThetaProdPlane_eta1_%i_eta2_%i.png", etaBin1, etaBin2));
      //____________________________________________________________________________________________________
    }
  }


  //_______________________________________________________________________________________________




  //QA histograms
  TCanvas *L0_L0bar_delta_eta_vs_delta_phi_can = new TCanvas("L0_L0bar_delta_eta_vs_delta_phi_can", "L0_L0bar_delta_eta_vs_delta_phi_can", 1200, 1000);

  L0_L0bar_delta_eta_vs_delta_phi_can->cd();

  L0_L0bar_delta_eta_vs_delta_phi_US_hist->Add(L0_L0bar_delta_eta_vs_delta_phi_LS_hist, -1);
  L0_L0bar_delta_eta_vs_delta_phi_US_hist->GetXaxis()->SetTitle("#Delta#eta");
  L0_L0bar_delta_eta_vs_delta_phi_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_delta_eta_vs_delta_phi_US_hist->GetYaxis()->SetTitle("#Delta#phi");
  L0_L0bar_delta_eta_vs_delta_phi_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_delta_eta_vs_delta_phi_US_hist->Draw("colz");

  L0_L0bar_delta_eta_vs_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0bar_delta_eta_vs_delta_phi.png");


  TCanvas *L0_L0bar_y1_vs_y2_can = new TCanvas("L0_L0bar_y1_vs_y2_can", "L0_L0bar_y1_vs_y2_can", 1200, 1000);

  L0_L0bar_y1_vs_y2_can->cd();

  L0_L0bar_y1_vs_y2_US_hist->Add(L0_L0bar_y1_vs_y2_LS_hist, -1);
  L0_L0bar_y1_vs_y2_US_hist->GetXaxis()->SetTitle("y_{#Lambda}");
  L0_L0bar_y1_vs_y2_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_y1_vs_y2_US_hist->GetYaxis()->SetTitle("y_{#bar{#Lambda}}");
  L0_L0bar_y1_vs_y2_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_y1_vs_y2_US_hist->Draw("colz");

  L0_L0bar_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0bar_y1_vs_y2.png");


  TCanvas *L0_L0bar_pT1_vs_pT2_can = new TCanvas("L0_L0bar_pT1_vs_pT2_can", "L0_L0bar_pT1_vs_pT2_can", 1200, 1000);

  L0_L0bar_pT1_vs_pT2_can->cd();

  L0_L0bar_pT1_vs_pT2_US_hist->Add(L0_L0bar_pT1_vs_pT2_LS_hist, -1);
  L0_L0bar_pT1_vs_pT2_US_hist->GetXaxis()->SetTitle("p_{T}(#Lambda) (GeV/#it{c})");
  L0_L0bar_pT1_vs_pT2_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_US_hist->GetYaxis()->SetTitle("p_{T}(#bar{#Lambda}) (GeV/#it{c})");
  L0_L0bar_pT1_vs_pT2_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_pT1_vs_pT2_US_hist->Draw("colz");

  L0_L0bar_pT1_vs_pT2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0bar_pT1_vs_pT2.png");


  TCanvas *L0_L0bar_phi1_vs_phi2_can = new TCanvas("L0_L0bar_phi1_vs_phi2_can", "L0_L0bar_phi1_vs_phi2_can", 1200, 1000);

  L0_L0bar_phi1_vs_phi2_can->cd();

  L0_L0bar_phi1_vs_phi2_US_hist->Add(L0_L0bar_phi1_vs_phi2_LS_hist, -1);
  L0_L0bar_phi1_vs_phi2_US_hist->GetXaxis()->SetTitle("#phi_{#Lambda}");
  L0_L0bar_phi1_vs_phi2_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_US_hist->GetYaxis()->SetTitle("#phi_{#bar{#Lambda}}");
  L0_L0bar_phi1_vs_phi2_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_phi1_vs_phi2_US_hist->Draw("colz");

  L0_L0bar_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0bar_phi1_vs_phi2.png");

  //_______________________________________________________________________________________________________________________________




  TCanvas *L0_L0_delta_eta_vs_delta_phi_can = new TCanvas("L0_L0_delta_eta_vs_delta_phi_can", "L0_L0_delta_eta_vs_delta_phi_can", 1200, 1000);

  L0_L0_delta_eta_vs_delta_phi_can->cd();

  L0_L0_delta_eta_vs_delta_phi_US_hist->Add(L0_L0_delta_eta_vs_delta_phi_LS_hist, -1);
  L0_L0_delta_eta_vs_delta_phi_US_hist->GetXaxis()->SetTitle("#Delta#eta");
  L0_L0_delta_eta_vs_delta_phi_US_hist->GetXaxis()->CenterTitle();
  L0_L0_delta_eta_vs_delta_phi_US_hist->GetYaxis()->SetTitle("#Delta#phi");
  L0_L0_delta_eta_vs_delta_phi_US_hist->GetYaxis()->CenterTitle();
  L0_L0_delta_eta_vs_delta_phi_US_hist->Draw("colz");

  L0_L0_delta_eta_vs_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0_delta_eta_vs_delta_phi.png");


  TCanvas *L0_L0_delta_eta_vs_delta_phi_zoom_can = new TCanvas("L0_L0_delta_eta_vs_delta_phi_zoom_can", "L0_L0_delta_eta_vs_delta_phi_zoom_can", 1200, 1000);

  L0_L0_delta_eta_vs_delta_phi_zoom_can->cd();

  L0_L0_delta_eta_vs_delta_phi_US_zoom_hist->Add(L0_L0_delta_eta_vs_delta_phi_LS_zoom_hist, -1);
  L0_L0_delta_eta_vs_delta_phi_US_zoom_hist->GetXaxis()->SetTitle("#Delta#eta");
  L0_L0_delta_eta_vs_delta_phi_US_zoom_hist->GetXaxis()->CenterTitle();
  L0_L0_delta_eta_vs_delta_phi_US_zoom_hist->GetYaxis()->SetTitle("#Delta#phi");
  L0_L0_delta_eta_vs_delta_phi_US_zoom_hist->GetYaxis()->CenterTitle();
  L0_L0_delta_eta_vs_delta_phi_US_zoom_hist->Draw("colz");

  L0_L0_delta_eta_vs_delta_phi_zoom_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0_delta_eta_vs_delta_phi_zoom.png");


  TCanvas *L0_L0_y1_vs_y2_can = new TCanvas("L0_L0_y1_vs_y2_can", "L0_L0_y1_vs_y2_can", 1200, 1000);

  L0_L0_y1_vs_y2_can->cd();

  L0_L0_y1_vs_y2_US_hist->Add(L0_L0_y1_vs_y2_LS_hist, -1);
  L0_L0_y1_vs_y2_US_hist->GetXaxis()->SetTitle("y_{#Lambda}");
  L0_L0_y1_vs_y2_US_hist->GetXaxis()->CenterTitle();
  L0_L0_y1_vs_y2_US_hist->GetYaxis()->SetTitle("y_{#Lambda}");
  L0_L0_y1_vs_y2_US_hist->GetYaxis()->CenterTitle();
  L0_L0_y1_vs_y2_US_hist->Draw("colz");

  L0_L0_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0_y1_vs_y2.png");


  TCanvas *L0_p_L0_p_y1_vs_y2_can = new TCanvas("L0_p_L0_p_y1_vs_y2_can", "L0_p_L0_p_y1_vs_y2_can", 1200, 1000);

  L0_p_L0_p_y1_vs_y2_can->cd();

  L0_p_L0_p_y1_vs_y2_US_hist->Add(L0_p_L0_p_y1_vs_y2_LS_hist, -1);
  L0_p_L0_p_y1_vs_y2_US_hist->GetXaxis()->SetTitle("y_{p}");
  L0_p_L0_p_y1_vs_y2_US_hist->GetXaxis()->CenterTitle();
  L0_p_L0_p_y1_vs_y2_US_hist->GetYaxis()->SetTitle("y_{p}");
  L0_p_L0_p_y1_vs_y2_US_hist->GetYaxis()->CenterTitle();
  L0_p_L0_p_y1_vs_y2_US_hist->Draw("colz");

  L0_p_L0_p_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_p_L0_p_y1_vs_y2.png");


  TCanvas *L0_pi_L0_pi_y1_vs_y2_can = new TCanvas("L0_pi_L0_pi_y1_vs_y2_can", "L0_pi_L0_pi_y1_vs_y2_can", 1200, 1000);

  L0_pi_L0_pi_y1_vs_y2_can->cd();

  L0_pi_L0_pi_y1_vs_y2_US_hist->Add(L0_pi_L0_pi_y1_vs_y2_LS_hist, -1);
  L0_pi_L0_pi_y1_vs_y2_US_hist->GetXaxis()->SetTitle("y_{#pi^{-}}");
  L0_pi_L0_pi_y1_vs_y2_US_hist->GetXaxis()->CenterTitle();
  L0_pi_L0_pi_y1_vs_y2_US_hist->GetYaxis()->SetTitle("y_{#pi^{-}}");
  L0_pi_L0_pi_y1_vs_y2_US_hist->GetYaxis()->CenterTitle();
  L0_pi_L0_pi_y1_vs_y2_US_hist->Draw("colz");

  L0_pi_L0_pi_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_pi_L0_pi_y1_vs_y2.png");


  TCanvas *L0_L0_pT1_vs_pT2_can = new TCanvas("L0_L0_pT1_vs_pT2_can", "L0_L0_pT1_vs_pT2_can", 1200, 1000);

  L0_L0_pT1_vs_pT2_can->cd();

  L0_L0_pT1_vs_pT2_US_hist->Add(L0_L0_pT1_vs_pT2_LS_hist, -1);
  L0_L0_pT1_vs_pT2_US_hist->GetXaxis()->SetTitle("p_{T}(#Lambda) (GeV/#it{c})");
  L0_L0_pT1_vs_pT2_US_hist->GetXaxis()->CenterTitle();
  L0_L0_pT1_vs_pT2_US_hist->GetYaxis()->SetTitle("p_{T}(#Lambda) (GeV/#it{c})");
  L0_L0_pT1_vs_pT2_US_hist->GetYaxis()->CenterTitle();
  L0_L0_pT1_vs_pT2_US_hist->Draw("colz");

  L0_L0_pT1_vs_pT2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0_pT1_vs_pT2.png");


  TCanvas *L0_L0_phi1_vs_phi2_can = new TCanvas("L0_L0_phi1_vs_phi2_can", "L0_L0_phi1_vs_phi2_can", 1200, 1000);

  L0_L0_phi1_vs_phi2_can->cd();

  L0_L0_phi1_vs_phi2_US_hist->Add(L0_L0_phi1_vs_phi2_LS_hist, -1);
  L0_L0_phi1_vs_phi2_US_hist->GetXaxis()->SetTitle("#phi_{#Lambda}");
  L0_L0_phi1_vs_phi2_US_hist->GetXaxis()->CenterTitle();
  L0_L0_phi1_vs_phi2_US_hist->GetYaxis()->SetTitle("#phi_{#Lambda}");
  L0_L0_phi1_vs_phi2_US_hist->GetYaxis()->CenterTitle();
  L0_L0_phi1_vs_phi2_US_hist->Draw("colz");

  L0_L0_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_L0_phi1_vs_phi2.png");


  TCanvas *L0_p_L0_p_phi1_vs_phi2_can = new TCanvas("L0_p_L0_p_phi1_vs_phi2_can", "L0_p_L0_p_phi1_vs_phi2_can", 1200, 1000);

  L0_p_L0_p_phi1_vs_phi2_can->cd();

  L0_p_L0_p_phi1_vs_phi2_US_hist->Add(L0_p_L0_p_phi1_vs_phi2_LS_hist, -1);
  L0_p_L0_p_phi1_vs_phi2_US_hist->GetXaxis()->SetTitle("#phi_{p}");
  L0_p_L0_p_phi1_vs_phi2_US_hist->GetXaxis()->CenterTitle();
  L0_p_L0_p_phi1_vs_phi2_US_hist->GetYaxis()->SetTitle("#phi_{p}");
  L0_p_L0_p_phi1_vs_phi2_US_hist->GetYaxis()->CenterTitle();
  L0_p_L0_p_phi1_vs_phi2_US_hist->Draw("colz");

  L0_p_L0_p_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_p_L0_p_phi1_vs_phi2.png");


  TCanvas *L0_pi_L0_pi_phi1_vs_phi2_can = new TCanvas("L0_pi_L0_pi_phi1_vs_phi2_can", "L0_pi_L0_pi_phi1_vs_phi2_can", 1200, 1000);

  L0_pi_L0_pi_phi1_vs_phi2_can->cd();

  L0_pi_L0_pi_phi1_vs_phi2_US_hist->Add(L0_pi_L0_pi_phi1_vs_phi2_LS_hist, -1);
  L0_pi_L0_pi_phi1_vs_phi2_US_hist->GetXaxis()->SetTitle("#phi_{#pi^{-}}");
  L0_pi_L0_pi_phi1_vs_phi2_US_hist->GetXaxis()->CenterTitle();
  L0_pi_L0_pi_phi1_vs_phi2_US_hist->GetYaxis()->SetTitle("#phi_{#pi^{-}}");
  L0_pi_L0_pi_phi1_vs_phi2_US_hist->GetYaxis()->CenterTitle();
  L0_pi_L0_pi_phi1_vs_phi2_US_hist->Draw("colz");

  L0_pi_L0_pi_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0_pi_L0_pi_phi1_vs_phi2.png");

  //_______________________________________________________________________________________________________________________________


  TCanvas *L0bar_L0bar_delta_eta_vs_delta_phi_can = new TCanvas("L0bar_L0bar_delta_eta_vs_delta_phi_can", "L0bar_L0bar_delta_eta_vs_delta_phi_can", 1200, 1000);

  L0bar_L0bar_delta_eta_vs_delta_phi_can->cd();

  L0bar_L0bar_delta_eta_vs_delta_phi_US_hist->Add(L0bar_L0bar_delta_eta_vs_delta_phi_LS_hist, -1);
  L0bar_L0bar_delta_eta_vs_delta_phi_US_hist->GetXaxis()->SetTitle("#Delta#eta");
  L0bar_L0bar_delta_eta_vs_delta_phi_US_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_delta_eta_vs_delta_phi_US_hist->GetYaxis()->SetTitle("#Delta#phi");
  L0bar_L0bar_delta_eta_vs_delta_phi_US_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_delta_eta_vs_delta_phi_US_hist->Draw("colz");

  L0bar_L0bar_delta_eta_vs_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_L0bar_delta_eta_vs_delta_phi.png");


  TCanvas *L0bar_L0bar_delta_eta_vs_delta_phi_zoom_can = new TCanvas("L0bar_L0bar_delta_eta_vs_delta_phi_zoom_can", "L0bar_L0bar_delta_eta_vs_delta_phi_zoom_can", 1200, 1000);

  L0bar_L0bar_delta_eta_vs_delta_phi_zoom_can->cd();

  L0bar_L0bar_delta_eta_vs_delta_phi_US_zoom_hist->Add(L0bar_L0bar_delta_eta_vs_delta_phi_LS_zoom_hist, -1);
  L0bar_L0bar_delta_eta_vs_delta_phi_US_zoom_hist->GetXaxis()->SetTitle("#Delta#eta");
  L0bar_L0bar_delta_eta_vs_delta_phi_US_zoom_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_delta_eta_vs_delta_phi_US_zoom_hist->GetYaxis()->SetTitle("#Delta#phi");
  L0bar_L0bar_delta_eta_vs_delta_phi_US_zoom_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_delta_eta_vs_delta_phi_US_zoom_hist->Draw("colz");

  L0bar_L0bar_delta_eta_vs_delta_phi_zoom_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_L0bar_delta_eta_vs_delta_phi_zoom.png");


  TCanvas *L0bar_L0bar_y1_vs_y2_can = new TCanvas("L0bar_L0bar_y1_vs_y2_can", "L0bar_L0bar_y1_vs_y2_can", 1200, 1000);

  L0bar_L0bar_y1_vs_y2_can->cd();

  L0bar_L0bar_y1_vs_y2_US_hist->Add(L0bar_L0bar_y1_vs_y2_LS_hist, -1);
  L0bar_L0bar_y1_vs_y2_US_hist->GetXaxis()->SetTitle("y_{#bar{#Lambda}}");
  L0bar_L0bar_y1_vs_y2_US_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_y1_vs_y2_US_hist->GetYaxis()->SetTitle("y_{#bar{#Lambda}}");
  L0bar_L0bar_y1_vs_y2_US_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_y1_vs_y2_US_hist->Draw("colz");

  L0bar_L0bar_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_L0bar_y1_vs_y2.png");


  TCanvas *L0bar_p_L0bar_p_y1_vs_y2_can = new TCanvas("L0bar_p_L0bar_p_y1_vs_y2_can", "L0bar_p_L0bar_p_y1_vs_y2_can", 1200, 1000);

  L0bar_p_L0bar_p_y1_vs_y2_can->cd();

  L0bar_p_L0bar_p_y1_vs_y2_US_hist->Add(L0bar_p_L0bar_p_y1_vs_y2_LS_hist, -1);
  L0bar_p_L0bar_p_y1_vs_y2_US_hist->GetXaxis()->SetTitle("y_{#bar{p}}");
  L0bar_p_L0bar_p_y1_vs_y2_US_hist->GetXaxis()->CenterTitle();
  L0bar_p_L0bar_p_y1_vs_y2_US_hist->GetYaxis()->SetTitle("y_{#bar{p}}");
  L0bar_p_L0bar_p_y1_vs_y2_US_hist->GetYaxis()->CenterTitle();
  L0bar_p_L0bar_p_y1_vs_y2_US_hist->Draw("colz");

  L0bar_p_L0bar_p_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_p_L0bar_p_y1_vs_y2.png");


  TCanvas *L0bar_pi_L0bar_pi_y1_vs_y2_can = new TCanvas("L0bar_pi_L0bar_pi_y1_vs_y2_can", "L0bar_pi_L0bar_pi_y1_vs_y2_can", 1200, 1000);

  L0bar_pi_L0bar_pi_y1_vs_y2_can->cd();

  L0bar_pi_L0bar_pi_y1_vs_y2_US_hist->Add(L0bar_pi_L0bar_pi_y1_vs_y2_LS_hist, -1);
  L0bar_pi_L0bar_pi_y1_vs_y2_US_hist->GetXaxis()->SetTitle("y_{#pi^{+}}");
  L0bar_pi_L0bar_pi_y1_vs_y2_US_hist->GetXaxis()->CenterTitle();
  L0bar_pi_L0bar_pi_y1_vs_y2_US_hist->GetYaxis()->SetTitle("y_{#pi^{+}}");
  L0bar_pi_L0bar_pi_y1_vs_y2_US_hist->GetYaxis()->CenterTitle();
  L0bar_pi_L0bar_pi_y1_vs_y2_US_hist->Draw("colz");

  L0bar_pi_L0bar_pi_y1_vs_y2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_pi_L0bar_pi_y1_vs_y2.png");


  TCanvas *L0bar_L0bar_pT1_vs_pT2_can = new TCanvas("L0bar_L0bar_pT1_vs_pT2_can", "L0bar_L0bar_pT1_vs_pT2_can", 1200, 1000);

  L0bar_L0bar_pT1_vs_pT2_can->cd();

  L0bar_L0bar_pT1_vs_pT2_US_hist->Add(L0bar_L0bar_pT1_vs_pT2_LS_hist, -1);
  L0bar_L0bar_pT1_vs_pT2_US_hist->GetXaxis()->SetTitle("p_{T}(#bar{#Lambda}) (GeV/#it{c})");
  L0bar_L0bar_pT1_vs_pT2_US_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_pT1_vs_pT2_US_hist->GetYaxis()->SetTitle("p_{T}(#bar{#Lambda}) (GeV/#it{c})");
  L0bar_L0bar_pT1_vs_pT2_US_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_pT1_vs_pT2_US_hist->Draw("colz");

  L0bar_L0bar_pT1_vs_pT2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_L0bar_pT1_vs_pT2.png");


  TCanvas *L0bar_L0bar_phi1_vs_phi2_can = new TCanvas("L0bar_L0bar_phi1_vs_phi2_can", "L0bar_L0bar_phi1_vs_phi2_can", 1200, 1000);

  L0bar_L0bar_phi1_vs_phi2_can->cd();

  L0bar_L0bar_phi1_vs_phi2_US_hist->Add(L0bar_L0bar_phi1_vs_phi2_LS_hist, -1);
  L0bar_L0bar_phi1_vs_phi2_US_hist->GetXaxis()->SetTitle("#phi_{#bar{#Lambda}}");
  L0bar_L0bar_phi1_vs_phi2_US_hist->GetXaxis()->CenterTitle();
  L0bar_L0bar_phi1_vs_phi2_US_hist->GetYaxis()->SetTitle("#phi_{#bar{#Lambda}}");
  L0bar_L0bar_phi1_vs_phi2_US_hist->GetYaxis()->CenterTitle();
  L0bar_L0bar_phi1_vs_phi2_US_hist->Draw("colz");

  L0bar_L0bar_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_L0bar_phi1_vs_phi2.png");


  TCanvas *L0bar_p_L0bar_p_phi1_vs_phi2_can = new TCanvas("L0bar_p_L0bar_p_phi1_vs_phi2_can", "L0bar_p_L0bar_p_phi1_vs_phi2_can", 1200, 1000);

  L0bar_p_L0bar_p_phi1_vs_phi2_can->cd();

  L0bar_p_L0bar_p_phi1_vs_phi2_US_hist->Add(L0bar_p_L0bar_p_phi1_vs_phi2_LS_hist, -1);
  L0bar_p_L0bar_p_phi1_vs_phi2_US_hist->GetXaxis()->SetTitle("#phi_{#bar{p}}");
  L0bar_p_L0bar_p_phi1_vs_phi2_US_hist->GetXaxis()->CenterTitle();
  L0bar_p_L0bar_p_phi1_vs_phi2_US_hist->GetYaxis()->SetTitle("#phi_{#bar{p}}");
  L0bar_p_L0bar_p_phi1_vs_phi2_US_hist->GetYaxis()->CenterTitle();
  L0bar_p_L0bar_p_phi1_vs_phi2_US_hist->Draw("colz");

  L0bar_p_L0bar_p_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_p_L0bar_p_phi1_vs_phi2.png");


  TCanvas *L0bar_pi_L0bar_pi_phi1_vs_phi2_can = new TCanvas("L0bar_pi_L0bar_pi_phi1_vs_phi2_can", "L0bar_pi_L0bar_pi_phi1_vs_phi2_can", 1200, 1000);

  L0bar_pi_L0bar_pi_phi1_vs_phi2_can->cd();

  L0bar_pi_L0bar_pi_phi1_vs_phi2_US_hist->Add(L0bar_pi_L0bar_pi_phi1_vs_phi2_LS_hist, -1);
  L0bar_pi_L0bar_pi_phi1_vs_phi2_US_hist->GetXaxis()->SetTitle("#phi_{#pi^{+}}");
  L0bar_pi_L0bar_pi_phi1_vs_phi2_US_hist->GetXaxis()->CenterTitle();
  L0bar_pi_L0bar_pi_phi1_vs_phi2_US_hist->GetYaxis()->SetTitle("#phi_{#pi^{+}}");
  L0bar_pi_L0bar_pi_phi1_vs_phi2_US_hist->GetYaxis()->CenterTitle();
  L0bar_pi_L0bar_pi_phi1_vs_phi2_US_hist->Draw("colz");

  L0bar_pi_L0bar_pi_phi1_vs_phi2_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L0bar_pi_L0bar_pi_phi1_vs_phi2.png");

  //_______________________________________________________________________________________________________________________________



  TCanvas *L_pair_cosThetaProdPlane_can = new TCanvas("L_pair_cosThetaProdPlane_can", "L_pair_cosThetaProdPlane_can", 1200, 1000);

  L_pair_cosThetaProdPlane_can->cd();

  L_pair_cosThetaProdPlane_US_hist = (TH1D*)L0_L0_cosThetaProdPlane_US_hist->Clone("L_pair_cosThetaProdPlane_US_hist");
  L_pair_cosThetaProdPlane_US_hist->Add(L0bar_L0bar_cosThetaProdPlane_US_hist);
  L_pair_cosThetaProdPlane_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
  L_pair_cosThetaProdPlane_US_hist->GetXaxis()->CenterTitle();
  //L_pair_cosThetaProdPlane_US_hist->GetYaxis()->SetTitle("1/N_{pair} d#it{N}/d cos(#theta*)");
  L_pair_cosThetaProdPlane_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
  L_pair_cosThetaProdPlane_US_hist->GetYaxis()->CenterTitle();
  L_pair_cosThetaProdPlane_US_hist->SetMarkerStyle(20);
  L_pair_cosThetaProdPlane_US_hist->SetMarkerColor(kRed);
  L_pair_cosThetaProdPlane_US_hist->SetLineColor(kRed);
  L_pair_cosThetaProdPlane_US_hist->Sumw2();
  L_pair_cosThetaProdPlane_US_hist->Scale(1./L_pair_cosThetaProdPlane_US_hist->GetXaxis()->GetBinWidth(1));
  L_pair_cosThetaProdPlane_US_hist->SetMinimum(0);
  L_pair_cosThetaProdPlane_US_hist->Draw("p e");

  L_pair_cosThetaProdPlane_LS_hist = (TH1D*)L0_L0_cosThetaProdPlane_LS_hist->Clone("L_pair_cosThetaProdPlane_LS_hist");
  L_pair_cosThetaProdPlane_LS_hist->Add(L0bar_L0bar_cosThetaProdPlane_LS_hist);
  L_pair_cosThetaProdPlane_LS_hist->SetMarkerStyle(20);
  L_pair_cosThetaProdPlane_LS_hist->SetMarkerColor(kBlue);
  L_pair_cosThetaProdPlane_LS_hist->Sumw2();
  L_pair_cosThetaProdPlane_LS_hist->Scale(1./L_pair_cosThetaProdPlane_LS_hist->GetXaxis()->GetBinWidth(1));
  L_pair_cosThetaProdPlane_LS_hist->Draw("p e same");

  TPaveText *LL_LbarLbar_text = new TPaveText(0.6, 0.2, 0.85, 0.49, "NDC");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  LL_LbarLbar_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  if(trigger == 0) LL_LbarLbar_text->AddText("Minimum bias");
  if(trigger == 1) LL_LbarLbar_text->AddText("JP2");
  LL_LbarLbar_text->AddText("(#bar{#Lambda^{0}}-#bar{#Lambda^{0}} + #Lambda^{0}-#Lambda^{0})");
  LL_LbarLbar_text->AddText("|#it{y}| < 1");
  LL_LbarLbar_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  LL_LbarLbar_text->SetFillColorAlpha(0, 0.01);
  LL_LbarLbar_text->Draw("same");

  L_pair_cosThetaProdPlane_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L_pair_cosThetaProdPlane_can.png");

  //________________________________________________________________________________________________________________________________________________________________

  //L and Lbar invariant mass from pairs with combinatorial background

  TCanvas *L_inv_mass_from_L_Lbar_can = new TCanvas("L_inv_mass_from_L_Lbar_can", "L_inv_mass_from_L_Lbar_can", 1200, 1000);
  L_inv_mass_from_L_Lbar_can->cd();

  L_inv_mass_US_from_L_Lbar->SetMarkerStyle(20);
  L_inv_mass_US_from_L_Lbar->SetMarkerColor(kRed);
  L_inv_mass_US_from_L_Lbar->SetLineColor(kRed);
  L_inv_mass_US_from_L_Lbar->GetXaxis()->SetTitle("M_{inv}^{p#pi^{-}} (GeV/#it{c}^{2})");
  L_inv_mass_US_from_L_Lbar->GetXaxis()->CenterTitle();
  L_inv_mass_US_from_L_Lbar->GetXaxis()->SetRangeUser(1.1, 1.14);
  L_inv_mass_US_from_L_Lbar->GetYaxis()->SetTitle("Counts");
  L_inv_mass_US_from_L_Lbar->GetYaxis()->CenterTitle();
  //L_inv_mass_US_from_L_Lbar->GetYaxis()->SetRangeUser(0, 3000);
  L_inv_mass_US_from_L_Lbar->Draw("p e");

  L_inv_mass_LS_from_L_Lbar->SetMarkerStyle(20);
  L_inv_mass_LS_from_L_Lbar->SetMarkerColor(kBlue);
  L_inv_mass_LS_from_L_Lbar->SetLineColor(kBlue);
  L_inv_mass_LS_from_L_Lbar->Draw("p e same");


  TLegend *SigAndBckgLeg = new TLegend(0.2, 0.65, 0.4, 0.75 );
  SigAndBckgLeg->SetTextFont(42);
  SigAndBckgLeg->AddEntry(L_inv_mass_US_from_L_Lbar, "Unlike-sign");
  SigAndBckgLeg->AddEntry(L_inv_mass_LS_from_L_Lbar, "Like-sign");
  SigAndBckgLeg->SetBorderSize(0);
  SigAndBckgLeg->Draw("same");

  TPaveText *L_Minv_text_pair = new TPaveText(0.6, 0.6, 0.89, 0.75, "NDC");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  L_Minv_text_pair->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  if(trigger == 0) L_Minv_text_pair->AddText("Minimum bias");
  if(trigger == 1) L_Minv_text_pair->AddText("JP2");
  L_Minv_text_pair->AddText("#Lambda^{0} from #Lambda^{0}-#bar{#Lambda^{0}}");
  L_Minv_text_pair->AddText("|#it{y}| < 1");
  L_Minv_text_pair->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L_Minv_text_pair->SetFillColorAlpha(0, 0.01);
  L_Minv_text_pair->Draw("same");


  L_inv_mass_from_L_Lbar_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L_inv_mass_from_L_Lbar.png");


  TCanvas *Lbar_inv_mass_from_L_Lbar_can = new TCanvas("Lbar_inv_mass_from_L_Lbar_can", "Lbar_inv_mass_from_L_Lbar_can", 1200, 1000);
  Lbar_inv_mass_from_L_Lbar_can->cd();

  Lbar_inv_mass_US_from_L_Lbar->SetMarkerStyle(20);
  Lbar_inv_mass_US_from_L_Lbar->SetMarkerColor(kRed);
  Lbar_inv_mass_US_from_L_Lbar->SetLineColor(kRed);
  Lbar_inv_mass_US_from_L_Lbar->GetXaxis()->SetTitle("M_{inv}^{#bar{p}#pi^{+}} (GeV/#it{c}^{2})");
  Lbar_inv_mass_US_from_L_Lbar->GetXaxis()->CenterTitle();
  Lbar_inv_mass_US_from_L_Lbar->GetXaxis()->SetRangeUser(1.1, 1.14);
  Lbar_inv_mass_US_from_L_Lbar->GetYaxis()->SetTitle("Counts");
  Lbar_inv_mass_US_from_L_Lbar->GetYaxis()->CenterTitle();
  //Lbar_inv_mass_US_from_L_Lbar->GetYaxis()->SetRangeUser(0, 7000);
  //Lbar_inv_mass_US_from_L_Lbar->GetYaxis()->SetRangeUser(0, 3000);
  Lbar_inv_mass_US_from_L_Lbar->Draw("p e");

  Lbar_inv_mass_LS_from_L_Lbar->SetMarkerStyle(20);
  Lbar_inv_mass_LS_from_L_Lbar->SetMarkerColor(kBlue);
  Lbar_inv_mass_LS_from_L_Lbar->SetLineColor(kBlue);
  Lbar_inv_mass_LS_from_L_Lbar->Draw("p e same");

  SigAndBckgLeg->Draw("same");

  TPaveText *Lbar_Minv_text_pair = new TPaveText(0.6, 0.6, 0.89, 0.75, "NDC");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  Lbar_Minv_text_pair->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  if(trigger == 0) Lbar_Minv_text_pair->AddText("Minimum bias");
  if(trigger == 1) Lbar_Minv_text_pair->AddText("JP2");
  Lbar_Minv_text_pair->AddText("#bar{#Lambda^{0}} from #Lambda^{0}-#bar{#Lambda^{0}}");
  Lbar_Minv_text_pair->AddText("|#it{y}| < 1");
  Lbar_Minv_text_pair->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  Lbar_Minv_text_pair->SetFillColorAlpha(0, 0.01);
  Lbar_Minv_text_pair->Draw("same");

  Lbar_inv_mass_from_L_Lbar_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/Lbar_inv_mass_from_L_Lbar.png");


  TCanvas *L_vs_Lbar_mass_can = new TCanvas("L_vs_Lbar_mass_can", "L_vs_Lbar_mass_can", 1200, 1200);
  L_vs_Lbar_mass_can->cd();

  L_vs_Lbar_mass_US->Add(L_vs_Lbar_mass_LS, -1);
  L_vs_Lbar_mass_US->GetXaxis()->SetTitle("M_{inv}^{p#pi^{-}} (GeV/#it{c}^{2})");
  L_vs_Lbar_mass_US->GetXaxis()->CenterTitle();
  L_vs_Lbar_mass_US->GetXaxis()->SetRangeUser(1.1, 1.14);
  L_vs_Lbar_mass_US->GetYaxis()->SetTitle("M_{inv}^{#bar{p}#pi^{+}} (GeV/#it{c}^{2})");
  L_vs_Lbar_mass_US->GetYaxis()->CenterTitle();
  L_vs_Lbar_mass_US->GetYaxis()->SetRangeUser(1.1, 1.14);
  L_vs_Lbar_mass_US->Draw("colz");

  L_vs_Lbar_mass_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L_vs_Lbar_mass.png");


  //_____________________________________________________________________________________

  TCanvas *L_inv_mass_from_L_L_can = new TCanvas("L_inv_mass_from_L_L_can", "L_inv_mass_from_L_L_can", 1200, 1000);
  L_inv_mass_from_L_L_can->cd();

  L_inv_mass_US_from_L_L->SetMarkerStyle(20);
  L_inv_mass_US_from_L_L->SetMarkerColor(kRed);
  L_inv_mass_US_from_L_L->SetLineColor(kRed);
  L_inv_mass_US_from_L_L->GetXaxis()->SetTitle("M_{inv}^{p#pi^{-}} (GeV/#it{c}^{2})");
  L_inv_mass_US_from_L_L->GetXaxis()->CenterTitle();
  L_inv_mass_US_from_L_L->GetXaxis()->SetRangeUser(1.1, 1.14);
  L_inv_mass_US_from_L_L->GetYaxis()->SetTitle("Counts");
  L_inv_mass_US_from_L_L->GetYaxis()->CenterTitle();
  //L_inv_mass_US_from_L_L->GetYaxis()->SetRangeUser(0, 7000);
  //L_inv_mass_US_from_L_L->GetYaxis()->SetRangeUser(0, 3000);
  L_inv_mass_US_from_L_L->Draw("p e");

  L_inv_mass_LS_from_L_L->SetMarkerStyle(20);
  L_inv_mass_LS_from_L_L->SetMarkerColor(kBlue);
  L_inv_mass_LS_from_L_L->SetLineColor(kBlue);
  L_inv_mass_LS_from_L_L->Draw("p e same");


  SigAndBckgLeg->Draw("same");

  TPaveText *L_Minv_text_pair_2 = new TPaveText(0.6, 0.6, 0.89, 0.75, "NDC");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  L_Minv_text_pair_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  if(trigger == 0) L_Minv_text_pair_2->AddText("Minimum bias");
  if(trigger == 1) L_Minv_text_pair_2->AddText("JP2");
  L_Minv_text_pair_2->AddText("#Lambda^{0} from #Lambda^{0}-#Lambda^{0}");
  L_Minv_text_pair_2->AddText("|#it{y}| < 1");
  L_Minv_text_pair_2->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  L_Minv_text_pair_2->SetFillColorAlpha(0, 0.01);
  L_Minv_text_pair_2->Draw("same");

  L_inv_mass_from_L_L_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L_inv_mass_from_L_L.png");


  TCanvas *L_vs_L_mass_can = new TCanvas("L_vs_L_mass_can", "L_vs_L_mass_can", 1200, 1200);
  L_vs_L_mass_can->cd();

  L_vs_L_mass_US->Add(L_vs_L_mass_LS, -1);
  L_vs_L_mass_US->GetXaxis()->SetTitle("M_{inv}^{p#pi^{-}} (GeV/#it{c}^{2})");
  L_vs_L_mass_US->GetXaxis()->CenterTitle();
  L_vs_L_mass_US->GetXaxis()->SetRangeUser(1.1, 1.14);
  L_vs_L_mass_US->GetYaxis()->SetTitle("M_{inv}^{#bar{p}#pi^{+}} (GeV/#it{c}^{2})");
  L_vs_L_mass_US->GetYaxis()->CenterTitle();
  L_vs_L_mass_US->GetYaxis()->SetRangeUser(1.1, 1.14);
  L_vs_L_mass_US->Draw("colz");

  L_vs_L_mass_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/L_vs_L_mass.png");



  //_____________________________________________________________________________________

  TCanvas *Lbar_inv_mass_from_Lbar_Lbar_can = new TCanvas("Lbar_inv_mass_from_Lbar_Lbar_can", "Lbar_inv_mass_from_Lbar_Lbar_can", 1200, 1000);
  Lbar_inv_mass_from_Lbar_Lbar_can->cd();

  Lbar_inv_mass_US_from_Lbar_Lbar->SetMarkerStyle(20);
  Lbar_inv_mass_US_from_Lbar_Lbar->SetMarkerColor(kRed);
  Lbar_inv_mass_US_from_Lbar_Lbar->SetLineColor(kRed);
  Lbar_inv_mass_US_from_Lbar_Lbar->GetXaxis()->SetTitle("M_{inv}^{p#pi^{-}} (GeV/#it{c}^{2})");
  Lbar_inv_mass_US_from_Lbar_Lbar->GetXaxis()->CenterTitle();
  Lbar_inv_mass_US_from_Lbar_Lbar->GetXaxis()->SetRangeUser(1.1, 1.14);
  Lbar_inv_mass_US_from_Lbar_Lbar->GetYaxis()->SetTitle("Counts");
  Lbar_inv_mass_US_from_Lbar_Lbar->GetYaxis()->CenterTitle();
  //Lbar_inv_mass_US_from_Lbar_Lbar->GetYaxis()->SetRangeUser(0, 7000);
  //Lbar_inv_mass_US_from_Lbar_Lbar->GetYaxis()->SetRangeUser(0, 3000);
  Lbar_inv_mass_US_from_Lbar_Lbar->Draw("p e");

  Lbar_inv_mass_LS_from_Lbar_Lbar->SetMarkerStyle(20);
  Lbar_inv_mass_LS_from_Lbar_Lbar->SetMarkerColor(kBlue);
  Lbar_inv_mass_LS_from_Lbar_Lbar->SetLineColor(kBlue);
  Lbar_inv_mass_LS_from_Lbar_Lbar->Draw("p e same");


  SigAndBckgLeg->Draw("same");

  TPaveText *Lbar_Minv_text_pair_2 = new TPaveText(0.6, 0.6, 0.89, 0.75, "NDC");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  Lbar_Minv_text_pair_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  if(trigger == 0) Lbar_Minv_text_pair_2->AddText("Minimum bias");
  if(trigger == 1) Lbar_Minv_text_pair_2->AddText("JP2");
  Lbar_Minv_text_pair_2->AddText("#bar{#Lambda^{0}} from #bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
  Lbar_Minv_text_pair_2->AddText("|#it{y}| < 1");
  Lbar_Minv_text_pair_2->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
  Lbar_Minv_text_pair_2->SetFillColorAlpha(0, 0.01);
  Lbar_Minv_text_pair_2->Draw("same");

  Lbar_inv_mass_from_Lbar_Lbar_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/Lbar_inv_mass_from_Lbar_Lbar.png");


  TCanvas *Lbar_vs_Lbar_mass_can = new TCanvas("Lbar_vs_Lbar_mass_can", "Lbar_vs_Lbar_mass_can", 1200, 1200);
  Lbar_vs_Lbar_mass_can->cd();

  Lbar_vs_Lbar_mass_US->Add(Lbar_vs_Lbar_mass_LS, -1);
  Lbar_vs_Lbar_mass_US->GetXaxis()->SetTitle("M_{inv}^{p#pi^{-}} (GeV/#it{c}^{2})");
  Lbar_vs_Lbar_mass_US->GetXaxis()->CenterTitle();
  Lbar_vs_Lbar_mass_US->GetXaxis()->SetRangeUser(1.1, 1.14);
  Lbar_vs_Lbar_mass_US->GetYaxis()->SetTitle("M_{inv}^{#bar{p}#pi^{+}} (GeV/#it{c}^{2})");
  Lbar_vs_Lbar_mass_US->GetYaxis()->CenterTitle();
  Lbar_vs_Lbar_mass_US->GetYaxis()->SetRangeUser(1.1, 1.14);
  Lbar_vs_Lbar_mass_US->Draw("colz");

  Lbar_vs_Lbar_mass_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_correlations/Lbar_vs_Lbar_mass.png");



  //----------------------------------L pair stats---------------------------------

  cout<<endl;
  cout<<"N events with L-Lbar pair: "<<nEventsWithLambdaPair_US_hist->GetBinContent(2)<<endl;
  cout<<"N L-Lbar pairs from hist: "<<nLLbar<<endl;
  cout<<"N L-Lbar background pairs from hist: "<<nLLbar_back<<endl;
  cout<<endl;
  cout<<"N events with 2+ L in event: "<<nLambdasInEvent_US_hist->Integral(3, 10)<<endl;
  cout<<"N L-L pairs from hist: "<<nLL<<endl;
  cout<<"N L-L background pairs from hist: "<<nLL_back<<endl;
  cout<<endl;
  cout<<"N events with 2+ Lbar in event: "<<nLambdaBarsInEvent_US_hist->Integral(3, 10)<<endl;
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


  return;

}
//__________________________________________________________________________________________________________________________

//ReadMode = 0 - read TTree, ReadMode = 1 - read histograms - First run in ReadMode = 0 to save relevant histograms, then can run in ReadMode = 1 to read just histograms and save time
//trigger = 0 - MB pi TOF, 1 - JP2, 2 - MB no TOF
//energy - collision energy in GeV
void Ana003_Lambda_corr(const int ReadMode = 0, const int trigger = 0, const int energy = 510, const int year = 2017)
{
  ifstream fileList;

  if(energy == 510)
  {
    //MB
    //if(trigger == 0) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run17/fileList.list");

    //JP2
    if(trigger == 1) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run17_JP2/fileList.list");

    //MB without strict TOF matching
    if(trigger == 0 || trigger == 2) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run17_MB_noTOF/fileList.list");

  }
  else if(energy == 200)
  {
    //MB
    //if(trigger == 0) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12/fileList.list");

    //JP2
    if(trigger == 1) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_JP2/fileList.list");

    //MB without strict TOF matching
    if(trigger == 0 ||trigger == 2) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12_MB_noTOF/fileList.list");

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


  //arrays to store invariant mass peak ranges
  //2 is for lower and upper range
  double invMassRange_L[2][nPtBins+1][nEtaBins+1];
  double invMassRange_L0[2][nPtBins+1][nEtaBins+1];
  double invMassRange_L0bar[2][nPtBins+1][nEtaBins+1];

  //Lambda signal to signal+background ratio
  //for bacground subtraction for polarization measurement
  double L0_kappa[nPtBins+1][nEtaBins+1];
  double L0bar_kappa[nPtBins+1][nEtaBins+1];


  bool invMassFinish = InvMass(myChain, invMassRange_L, invMassRange_L0, invMassRange_L0bar, L0_kappa, L0bar_kappa, ReadMode, trigger, energy , year);

  if(!invMassFinish)
  {
    cout<<"Analysis of invariant spectra ended abnormally. Abborting!"<<endl;

    return;
  }

  LambdaLambdaBarSpinCorr(myChain, invMassRange_L0, invMassRange_L0bar, L0_kappa, L0bar_kappa, ReadMode, trigger, energy , year);

  cout<<endl;
  cout<<"Nubmer of accepted events: "<<hEventStat1->GetBinContent(6)<<endl;


  return;
}
