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


using namespace std;

//const int nPtBins = 9;
//float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5., 7.};

const int nPtBins = 8;
float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};

const int nEtaBins = 3;
float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };

const double L_mass_PDG = 1.115683; //mass in GeV/c^2 from latest PDG

const float L_y_cut = 1;
const int strictTOF_cut = 2; //0 - hybrid TOF for both daughters, 1 - strict TOF for pions, 2 - strict TOF for both pion and proton
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

bool cuts(int pi_hasTOFinfo, int p_hasTOFinfo, float L_y, float L_theta, float L_decayL)
{

  if( !( TMath::Abs(L_y) < L_y_cut ) ) return false;
  if( strictTOF_cut == 1 && pi_hasTOFinfo == 0 ) return false; //TOF matched pions
  if( strictTOF_cut == 2 && (pi_hasTOFinfo == 0 || p_hasTOFinfo == 0) ) return false; //TOF matched pions and protons
  if(cos(L_theta) < L_cos_theta_cut) return false;
  if(L_decayL > L_decayL_cut) return false;

  return true;

}

//analyze invariant mass spectra
//arguments are vectors of bin numbers of the invariant mass peak
//the bins are determined via fit to the invariant mass spectra
//bool InvMass(TTree *L_tree, vector<int> &invMassBins_L, vector<int> &invMassBins_L0, vector<int> &invMassBins_L0bar, const int readMode)
bool InvMass(TChain *L_tree, double (&invMassRange_L)[2][nPtBins+1][nEtaBins+1], double (&invMassRange_L0)[2][nPtBins+1][nEtaBins+1], double (&invMassRange_L0bar)[2][nPtBins+1][nEtaBins+1], double (&L0_kappa)[nPtBins+1][nEtaBins+1], double (&L0bar_kappa)[nPtBins+1][nEtaBins+1], const int ReadMode, const int energy = 510, const int year = 2017)
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
  Float_t L_mass, L_pt, L_eta, L_decayL, L_theta, L_DCAdaughters;

  Float_t pi_pt, p_pt;
  Float_t pi_eta, p_eta;
  Float_t pi_phi, p_phi;
  Float_t pi_ch;
  Float_t pi_dca, p_dca;
  Int_t pi_hasTOFinfo, p_hasTOFinfo;

  Float_t thetaProdPlane;

  Int_t eventId;

  if(ReadMode == 0)
  {

    //new variable names
    //---------------SET BARANCH ADDRESSES------------------------
    L_tree->SetBranchAddress("pair_charge", &charge);
    L_tree->SetBranchAddress("pair_mass", &L_mass);
    L_tree->SetBranchAddress("pair_pt", &L_pt);
    L_tree->SetBranchAddress("pair_eta", &L_eta);
    L_tree->SetBranchAddress("pair_decayL", &L_decayL);
    L_tree->SetBranchAddress("pair_theta", &L_theta);
    L_tree->SetBranchAddress("pair_DCAdaughters", &L_DCAdaughters);

    //pion is particle 2 in the pair niside the TTree
    L_tree->SetBranchAddress("p2_pt", &pi_pt);
    L_tree->SetBranchAddress("p2_eta", &pi_eta);
    L_tree->SetBranchAddress("p2_phi", &pi_phi);
    L_tree->SetBranchAddress("p2_ch", &pi_ch);
    L_tree->SetBranchAddress("p2_dca", &pi_dca);
    L_tree->SetBranchAddress("p2_hasTOFinfo", &pi_hasTOFinfo);

    //proton is particle 1 in the pair niside the TTree
    L_tree->SetBranchAddress("p1_pt", &p_pt);
    L_tree->SetBranchAddress("p1_eta", &p_eta);
    L_tree->SetBranchAddress("p1_phi", &p_phi);
    L_tree->SetBranchAddress("p1_dca", &p_dca);
    L_tree->SetBranchAddress("p1_hasTOFinfo", &p_hasTOFinfo);

    L_tree->SetBranchAddress("thetaProdPlane", &thetaProdPlane);

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
    if( !cuts(pi_hasTOFinfo, p_hasTOFinfo, L_y, L_theta, L_decayL) ) continue;



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

      if( pi_ch == -1)
      {
        L0_inv_mass_US[pT_bin][eta_bin]->Fill(L_mass);
        L0_inv_mass_US[nPtBins][eta_bin]->Fill(L_mass);
        L0_inv_mass_US[pT_bin][nEtaBins]->Fill(L_mass);
        L0_inv_mass_US[nPtBins][nEtaBins]->Fill(L_mass);
      }

      if(pi_ch == 1)
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



        if( pi_ch == 1)
        {
          L0_inv_mass_LS[pT_bin][eta_bin]->Fill(L_mass);
          L0_inv_mass_LS[nPtBins][eta_bin]->Fill(L_mass);
          L0_inv_mass_LS[pT_bin][nEtaBins]->Fill(L_mass);
          L0_inv_mass_LS[nPtBins][nEtaBins]->Fill(L_mass);
        }

        if(pi_ch == -1)
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

        cout<< nEtaBins+1<<endl;

        InvMassCan_US_LS[pTbin][etaBin] = new TCanvas(Form("InvMassCan_US_LS_pT_%i_eta_%i", pTbin, etaBin), Form("InvMassCan_US_LS_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);
        InvMassCan[pTbin][etaBin] = new TCanvas(Form("InvMassCan_pT_%i_eta_%i", pTbin, etaBin), Form("InvMassCan_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);

        InvMassCanL0_US_LS[pTbin][etaBin] = new TCanvas(Form("L0_InvMassCan_US_LS_pT_%i_eta_%i", pTbin, etaBin), Form("L0_InvMassCan_US_LS_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);
        InvMassCanL0[pTbin][etaBin] = new TCanvas(Form("L0_InvMassCan_pT_%i_eta_%i", pTbin, etaBin), Form("L0_InvMassCan_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);

        InvMassCanL0bar_US_LS[pTbin][etaBin] = new TCanvas(Form("L0bar_InvMassCan_US_LS_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_InvMassCan_US_LS_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);
        InvMassCanL0bar[pTbin][etaBin] = new TCanvas(Form("L0bar_InvMassCan_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_InvMassCan_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);


        TString *pT_range = new TString();
        if(pTbin < nPtBins) pT_range->Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c}", pT_bins[pTbin], pT_bins[pTbin+1]);
        else pT_range->Form("p_{T} integrated");


        TString *eta_range = new TString();
        if(etaBin < nEtaBins) eta_range->Form("%0.1f < #eta < %0.1f", eta_bins[etaBin], eta_bins[etaBin+1]);
        else eta_range->Form("-1 < #eta < 1");



        TF1 *fitGauss_L = new TF1("fitGauss_L", "gaus(0)", 1.07, 2.);

        TF1 *fitGauss_L0 = new TF1("fitGaussL0", "gaus(0)", 1.07, 2.);

        TF1 *fitGauss_L0bar = new TF1("fitGauss_L0bar", "gaus(0)", 1.07, 2.);


        InvMassCan_US_LS[pTbin][etaBin]->cd();

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

        L_inv_mass_LS[pTbin][etaBin]->SetMarkerStyle(20);
        L_inv_mass_LS[pTbin][etaBin]->SetMarkerColor(kBlue);
        L_inv_mass_LS[pTbin][etaBin]->SetLineColor(kBlue);
        L_inv_mass_LS[pTbin][etaBin]->Draw("p e same");

        TLegend *SigAndBckgLeg = new TLegend(0.2, 0.65, 0.4, 0.89 );
        SigAndBckgLeg->AddEntry(L_inv_mass_US[pTbin][etaBin], "Unlike-sign");
        SigAndBckgLeg->AddEntry(L_inv_mass_LS[pTbin][etaBin], "Like-sign");
        SigAndBckgLeg->SetBorderSize(0);
        SigAndBckgLeg->Draw("same");

        TPaveText *cent_text = new TPaveText(0.6, 0.6, 0.89, 0.75, "NDC");
        //cent_text->AddText("STAR preliminary");
        //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
        cent_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
        cent_text->AddText("Minimum bias");
        cent_text->AddText("#Lambda^{0} and #bar{#Lambda^{0}}");
        cent_text->AddText(eta_range->Data());
        cent_text->AddText(pT_range->Data());
        cent_text->SetFillColorAlpha(0, 0.01);
        cent_text->Draw("same");

        InvMassCan_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_inv_mass_US_LS_pT_%i_eta_%i.png", pTbin, etaBin));
        //InvMassCan_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L_inv_mass_US_LS_pT_%i_eta_%i.pdf", pTbin, etaBin));


        InvMassCan[pTbin][etaBin]->cd();

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

          if(/*!fit_res_gaus_L ||*/ !fit_res_gaus_L->IsValid())
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

        TLegend *SigMinusBckgLeg = new TLegend(0.2, 0.65, 0.4, 0.89 );
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

        L0_inv_mass_US[pTbin][etaBin]->SetMarkerStyle(20);
        L0_inv_mass_US[pTbin][etaBin]->SetMarkerColor(kRed);
        L0_inv_mass_US[pTbin][etaBin]->SetLineColor(kRed);
        L0_inv_mass_US[pTbin][etaBin]->GetXaxis()->SetTitle("M_{inv}^{p#pi^{-}} (GeV/#it{c}^{2})");
        L0_inv_mass_US[pTbin][etaBin]->GetXaxis()->CenterTitle();
        L0_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
        L0_inv_mass_US[pTbin][etaBin]->GetYaxis()->CenterTitle();
        //L0_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 7000);
        //L0_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 3000);
        L0_inv_mass_US[pTbin][etaBin]->Draw("p e");

        L0_inv_mass_LS[pTbin][etaBin]->SetMarkerStyle(20);
        L0_inv_mass_LS[pTbin][etaBin]->SetMarkerColor(kBlue);
        L0_inv_mass_LS[pTbin][etaBin]->SetLineColor(kBlue);
        L0_inv_mass_LS[pTbin][etaBin]->Draw("p e same");


        SigAndBckgLeg->Draw("same");

        TPaveText *cent_text_2 = new TPaveText(0.6, 0.6, 0.89, 0.75, "NDC");
        //cent_text->AddText("STAR preliminary");
        //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
        cent_text_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
        cent_text_2->AddText("Minimum bias");
        cent_text_2->AddText("#Lambda^{0}");
        cent_text_2->AddText(eta_range->Data());
        cent_text_2->AddText(pT_range->Data());
        cent_text_2->SetFillColorAlpha(0, 0.01);
        cent_text_2->Draw("same");

        InvMassCanL0_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0_inv_mass_US_LS_pT_%i_eta_%i.png", pTbin, etaBin));
        //InvMassCanL0_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0_inv_mass_US_LS_pT_%i_eta_%i.pdf", pTbin, etaBin));


        InvMassCanL0[pTbin][etaBin]->cd();

        L0_inv_mass[pTbin][etaBin] = (TH1D*)L0_inv_mass_US[pTbin][etaBin]->Clone();
        L0_inv_mass[pTbin][etaBin]->Sumw2();
        L0_inv_mass[pTbin][etaBin]->Add(L0_inv_mass_LS[pTbin][etaBin], -1);
        L0_inv_mass[pTbin][etaBin]->GetXaxis()->SetTitle("M_{inv}^{p#pi^{-}} (GeV/#it{c}^{2})");
        L0_inv_mass[pTbin][etaBin]->GetXaxis()->CenterTitle();
        L0_inv_mass[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
        L0_inv_mass[pTbin][etaBin]->GetYaxis()->CenterTitle();
        L0_inv_mass[pTbin][etaBin]->SetMarkerStyle(20);
        L0_inv_mass[pTbin][etaBin]->SetMarkerColor(kRed);
        L0_inv_mass[pTbin][etaBin]->SetLineColor(kRed);
        L0_inv_mass[pTbin][etaBin]->Draw("p e");


        if( ( (energy == 510 || strictTOF_cut == 2) && pTbin > 0 ) || energy == 200 )//skip first pT bin for 510 GeV - no peak visisble
        {
          //if( pTbin == 1 ) fitGauss_L0->SetParameters(2000, 1.116, 0.002);
          if( ( ( energy == 510 || strictTOF_cut == 2 ) && pTbin == 1) ) fitGauss_L0->SetParameters(2000, 1.116, 0.002);
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

        L0bar_inv_mass_US[pTbin][etaBin]->SetMarkerStyle(20);
        L0bar_inv_mass_US[pTbin][etaBin]->SetMarkerColor(kRed);
        L0bar_inv_mass_US[pTbin][etaBin]->SetLineColor(kRed);
        L0bar_inv_mass_US[pTbin][etaBin]->GetXaxis()->SetTitle("M_{inv}^{#bar{p}#pi^{+}} (GeV/#it{c}^{2})");
        L0bar_inv_mass_US[pTbin][etaBin]->GetXaxis()->CenterTitle();
        L0bar_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
        L0bar_inv_mass_US[pTbin][etaBin]->GetYaxis()->CenterTitle();
        //L0bar_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 7000);
        //L0bar_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 3000);
        L0bar_inv_mass_US[pTbin][etaBin]->Draw("p e");

        L0bar_inv_mass_LS[pTbin][etaBin]->SetMarkerStyle(20);
        L0bar_inv_mass_LS[pTbin][etaBin]->SetMarkerColor(kBlue);
        L0bar_inv_mass_LS[pTbin][etaBin]->SetLineColor(kBlue);
        L0bar_inv_mass_LS[pTbin][etaBin]->Draw("p e same");


        SigAndBckgLeg->Draw("same");

        TPaveText *cent_text_3 = new TPaveText(0.6, 0.6, 0.89, 0.75, "NDC");
        //cent_text->AddText("STAR preliminary");
        //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
        cent_text_3->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
        cent_text_3->AddText("Minimum bias");
        cent_text_3->AddText("#bar{#Lambda^{0}}");
        cent_text_3->AddText(eta_range->Data());
        cent_text_3->AddText(pT_range->Data());
        cent_text_3->SetFillColorAlpha(0, 0.01);
        cent_text_3->Draw("same");

        InvMassCanL0bar_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0bar_inv_mass_US_LS_pT_%i_eta_%i.png", pTbin, etaBin));
        //InvMassCanL0bar_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0bar_inv_mass_US_LS_pT_%i_eta_%i.pdf", pTbin, etaBin));


        InvMassCanL0bar[pTbin][etaBin]->cd();

        L0bar_inv_mass[pTbin][etaBin] = (TH1D*)L0bar_inv_mass_US[pTbin][etaBin]->Clone();
        L0bar_inv_mass[pTbin][etaBin]->Sumw2();
        L0bar_inv_mass[pTbin][etaBin]->Add(L0bar_inv_mass_LS[pTbin][etaBin], -1);
        L0bar_inv_mass[pTbin][etaBin]->GetXaxis()->SetTitle("M_{inv}^{#bar{p}#pi^{+}} (GeV/#it{c}^{2})");
        L0bar_inv_mass[pTbin][etaBin]->GetXaxis()->CenterTitle();
        L0bar_inv_mass[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
        L0bar_inv_mass[pTbin][etaBin]->GetYaxis()->CenterTitle();
        L0bar_inv_mass[pTbin][etaBin]->SetMarkerStyle(20);
        L0bar_inv_mass[pTbin][etaBin]->SetMarkerColor(kRed);
        L0bar_inv_mass[pTbin][etaBin]->SetLineColor(kRed);
        L0bar_inv_mass[pTbin][etaBin]->Draw("p e");


        if( ( (energy == 510 || strictTOF_cut == 2) && pTbin > 0 ) || energy == 200 )//skip first pT bin for 510 GeV - no peak visisble
        {
          //if(   pTbin == 1 ) fitGauss_L0bar->SetParameters(2000, 1.116, 0.002);
          if( ( (energy == 510 || strictTOF_cut == 2) && pTbin == 1) ) fitGauss_L0bar->SetParameters(2000, 1.116, 0.002);
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
void LambdaPolarization(TChain *L_tree, double invMassRange_L0[2][nPtBins+1][nEtaBins+1], double invMassRange_L0bar[2][nPtBins+1][nEtaBins+1], double L0_kappa[nPtBins+1][nEtaBins+1], double L0bar_kappa[nPtBins+1][nEtaBins+1], const int ReadMode, const int energy = 510, const int year = 2017)
{
  const float L0_alpha = 0.732; //decay parameter of L0
  const float L0bar_alpha = -0.758; //decay paramteter of L0bar

  TFile *ProdPlaneOutFile; //output file to store production plane histograms

  if(ReadMode == 0) //create production plane file from nTuple - run in this mode first
  {
    ProdPlaneOutFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/ProdPlane_Lambda_work.root", "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT
  }
  else  //read file to store wrong and good sign M_inv distributions and daughter pT
  {

    ProdPlaneOutFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/ProdPlane_Lambda_work.root", "read"); //old, non-optimized cuts

    if( !(ProdPlaneOutFile->IsOpen()) )
    {
      cout<<"Unable to open file with production plane histograms!"<<endl;
      return;
    }
  }


  TFile *EffFile;

  if(energy == 510) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run17_new_200M_LL.root", "read");
  else if(energy == 200) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/L_cosThetaStar_eff_Run12_new_10M.root", "read");
  else
  {
    cout<<"Not a valid colliison energy! Abborting!"<<endl;
    return;
  }

  //_______________________________________________________________________________________________________________________________________________


  TH1D *L0_thetaProdPlane_US_hist[nPtBins+1][nEtaBins+1];
  TH1D *L0_thetaProdPlane_LS_hist[nPtBins+1][nEtaBins+1];


  TH1D *L0_cosThetaProdPlane_US_hist[nPtBins+1][nEtaBins+1];
  TH1D *L0_cosThetaProdPlane_LS_hist[nPtBins+1][nEtaBins+1];


  TH1D *L0bar_thetaProdPlane_US_hist[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_thetaProdPlane_LS_hist[nPtBins+1][nEtaBins+1];

  TH1D *L0bar_cosThetaProdPlane_US_hist[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_cosThetaProdPlane_LS_hist[nPtBins+1][nEtaBins+1];

  TH1D *L0_polarization[nEtaBins+1];
  TH1D *L0bar_polarization[nEtaBins+1];

  TH1D *L0_xF_US[nEtaBins+1];
  TH1D *L0_xF_LS[nEtaBins+1];

  TH1D *L0bar_xF_US[nEtaBins+1];
  TH1D *L0bar_xF_LS[nEtaBins+1];
  //________________________________________________________________

  //efficiency histograms
  TH1D *L0_cosThetaProdPlane_eff_hist[nPtBins+1][nEtaBins+1];
  TH1D *L0bar_cosThetaProdPlane_eff_hist[nPtBins+1][nEtaBins+1];

  //______________________________________________________________


  TH1D *nLambdaPairsInEvent_US_hist = new TH1D("nLambdaPairsInEvent_US_hist", "nLambdaPairsInEvent_US_hist", 10, 0, 10);
  TH1D *nLambdasInEvent_US_hist = new TH1D("nLambdasInEvent_US_hist", "nLambdasInEvent_US_hist", 10, 0, 10);
  TH1D *nLambdaBarsInEvent_US_hist = new TH1D("nLambdaBarsInEvent_US_hist", "nLambdaBarsInEvent_US_hist", 10, 0, 10);

  TH1D *nEventsWithLambdaPair_US_hist = new TH1D("nEventsWithLambdaPair_US_hist", "nEventsWithLambdaPair_US_hist", 2, 0, 2);

  TH1D *nLambdaPairsInEvent_LS_hist = new TH1D("nLambdaPairsInEvent_LS_hist", "nLambdaPairsInEvent_LS_hist", 10, 0, 10);
  TH1D *nLambdasInEvent_LS_hist = new TH1D("nLambdasInEvent_LS_hist", "nLambdasInEvent_LS_hist", 10, 0, 10);
  TH1D *nLambdaBarsInEvent_LS_hist = new TH1D("nLambdaBarsInEvent_LS_hist", "nLambdaBarsInEvent_LS_hist", 10, 0, 10);

  TH1D *nEventsWithLambdaPair_LS_hist = new TH1D("nEventsWithLambdaPair_LS_hist", "nEventsWithLambdaPair_LS_hist", 2, 0, 2);


  if(ReadMode == 0) //create histograms to be saved into file
  {
    for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
    {
      for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
      {
        if(pTbin == 0 )
        {
          L0_polarization[etaBin] = new TH1D(Form("L0_polarization_eta_%i", etaBin), Form("L0_polarization_eta_%i", etaBin), nPtBins, pT_bins);
          L0bar_polarization[etaBin] = new TH1D(Form("L0bar_polarization_eta_%i", etaBin), Form("L0bar_polarization_eta_%i", etaBin), nPtBins, pT_bins);

          L0_xF_US[etaBin] = new TH1D(Form("L0_xF_US_eta_%i", etaBin), Form("L0_xF_US_eta_%i", etaBin), 100, 0, 0.01);
          L0_xF_LS[etaBin] = new TH1D(Form("L0_xF_LS_eta_%i", etaBin), Form("L0_xF_LS_eta_%i", etaBin), 100, 0, 0.01);

          L0bar_xF_US[etaBin] = new TH1D(Form("L0bar_xF_US_eta_%i", etaBin), Form("L0bar_xF_US_eta_%i", etaBin), 100, 0, 0.01);
          L0bar_xF_LS[etaBin] = new TH1D(Form("L0bar_xF_LS_eta_%i", etaBin), Form("L0bar_xF_LS_eta_%i", etaBin), 100, 0, 0.01);
        }


        //theta star histograms
        L0_thetaProdPlane_US_hist[pTbin][etaBin] = new TH1D(Form("L0_thetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin), Form("L0_thetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());
        L0_thetaProdPlane_LS_hist[pTbin][etaBin] = new TH1D(Form("L0_thetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin), Form("L0_thetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());

        L0_cosThetaProdPlane_US_hist[pTbin][etaBin] = new TH1D(Form("L0_cosThetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin), Form("L0_cosThetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);
        L0_cosThetaProdPlane_LS_hist[pTbin][etaBin] = new TH1D(Form("L0_cosThetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin), Form("L0_cosThetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);

        L0bar_thetaProdPlane_US_hist[pTbin][etaBin] = new TH1D(Form("L0bar_thetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_thetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());
        L0bar_thetaProdPlane_LS_hist[pTbin][etaBin] = new TH1D(Form("L0bar_thetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_thetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());

        L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin] = new TH1D(Form("L0bar_cosThetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_cosThetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);
        L0bar_cosThetaProdPlane_LS_hist[pTbin][etaBin] = new TH1D(Form("L0bar_cosThetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_cosThetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);

        //efficiency histograms
        L0_cosThetaProdPlane_eff_hist[pTbin][etaBin] = (TH1D*)EffFile->Get(Form("L0_cosThetaProdPlane_eff_pT_%i_eta_%i", pTbin, etaBin));
        L0bar_cosThetaProdPlane_eff_hist[pTbin][etaBin] = (TH1D*)EffFile->Get(Form("L0bar_cosThetaProdPlane_eff_pT_%i_eta_%i", pTbin, etaBin));




      }

    }

  }
  else //load histograms from file
  {
    for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
    {
      for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
      {

        L0_thetaProdPlane_US_hist[pTbin][etaBin] = (TH1D*)ProdPlaneOutFile->Get(Form("L0_thetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin));
        L0_thetaProdPlane_LS_hist[pTbin][etaBin] = (TH1D*)ProdPlaneOutFile->Get(Form("L0_thetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin));

        L0_cosThetaProdPlane_US_hist[pTbin][etaBin] = (TH1D*)ProdPlaneOutFile->Get(Form("L0_cosThetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin));
        L0_cosThetaProdPlane_LS_hist[pTbin][etaBin] = (TH1D*)ProdPlaneOutFile->Get(Form("L0_cosThetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin));

        L0bar_thetaProdPlane_US_hist[pTbin][etaBin] = (TH1D*)ProdPlaneOutFile->Get(Form("L0bar_thetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin));
        L0bar_thetaProdPlane_LS_hist[pTbin][etaBin] = (TH1D*)ProdPlaneOutFile->Get(Form("L0bar_thetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin));

        L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin] = (TH1D*)ProdPlaneOutFile->Get(Form("L0bar_cosThetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin));
        L0bar_cosThetaProdPlane_LS_hist[pTbin][etaBin] = (TH1D*)ProdPlaneOutFile->Get(Form("L0bar_cosThetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin));


        //efficiency histograms
        L0_cosThetaProdPlane_eff_hist[pTbin][etaBin] = (TH1D*)EffFile->Get(Form("L0_cosThetaProdPlane_eff_pT_%i_eta_%i", pTbin, etaBin));
        L0bar_cosThetaProdPlane_eff_hist[pTbin][etaBin] = (TH1D*)EffFile->Get(Form("L0bar_cosThetaProdPlane_eff_pT_%i_eta_%i", pTbin, etaBin));

      }

    }

  }

  //________________________________________________________________________________________

  //Labmda and Lambda_bar stats

  int nLambda = -1; //default value
  int nLambdaBar = -1;
  //int nLambdaPairs = -1;

  int nLambda_background = -1; //default value
  int nLambdaBar_background = -1;

  int eventID_last = -1; //to store event ID from last L candidate
  int eventID_last_background = -1; //to store event ID from last L candidate


  Long64_t nEntries = 0; //total nEntries

  Int_t charge;
  Float_t L_mass, L_pt, L_eta, L_decayL, L_theta, L_DCAdaughters;

  Float_t pi_pt, p_pt;
  Float_t pi_eta, p_eta;
  Float_t pi_phi, p_phi;
  Float_t pi_ch;
  Float_t pi_dca, p_dca;
  Int_t pi_hasTOFinfo, p_hasTOFinfo;

  Float_t thetaProdPlane;

  Int_t eventId;

  if(ReadMode == 0)
  {
    //new variable names
    //---------------SET BARANCH ADDRESSES------------------------
    L_tree->SetBranchAddress("pair_charge", &charge);
    L_tree->SetBranchAddress("pair_mass", &L_mass);
    L_tree->SetBranchAddress("pair_pt", &L_pt);
    L_tree->SetBranchAddress("pair_eta", &L_eta);
    L_tree->SetBranchAddress("pair_decayL", &L_decayL);
    L_tree->SetBranchAddress("pair_theta", &L_theta);
    L_tree->SetBranchAddress("pair_DCAdaughters", &L_DCAdaughters);

    //pion is particle 2 in the pair niside the TTree
    L_tree->SetBranchAddress("p2_pt", &pi_pt);
    L_tree->SetBranchAddress("p2_eta", &pi_eta);
    L_tree->SetBranchAddress("p2_phi", &pi_phi);
    L_tree->SetBranchAddress("p2_ch", &pi_ch);
    L_tree->SetBranchAddress("p2_dca", &pi_dca);
    L_tree->SetBranchAddress("p2_hasTOFinfo", &pi_hasTOFinfo);

    //proton is particle 1 in the pair niside the TTree
    L_tree->SetBranchAddress("p1_pt", &p_pt);
    L_tree->SetBranchAddress("p1_eta", &p_eta);
    L_tree->SetBranchAddress("p1_phi", &p_phi);
    L_tree->SetBranchAddress("p1_dca", &p_dca);
    L_tree->SetBranchAddress("p1_hasTOFinfo", &p_hasTOFinfo);

    L_tree->SetBranchAddress("thetaProdPlane", &thetaProdPlane);

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



    double L_xF = fabs(pz(L_pt, L_eta))/energy*2.; //energy is in CMS, need energz of one proton

    //calculate Lambda rapidity y
    double L_y = rapidity(L_pt, L_eta, L_mass_PDG);




    if(nLambda == -1 || nLambdaBar == -1) //first iteration
    {
      eventID_last = eventId;
      eventID_last_background = eventId;

      nLambda = 0;
      nLambdaBar = 0;
      //nLambdaPairs = 0;

      nLambda_background = 0;
      nLambdaBar_background = 0;
    }

    //------------------------------------------------------------------------------------------------------------------

    //cuts
    if( !cuts(pi_hasTOFinfo, p_hasTOFinfo, L_y, L_theta, L_decayL) ) continue;

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

    if( pT_bin == -1 ) continue;


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



    if(charge == 0 ) //like-sign combinations
    {

      if(L_mass > invMassRange_L0[0][pT_bin][eta_bin] && L_mass < invMassRange_L0[1][pT_bin][eta_bin])
      {

        if( pi_ch == -1)
        {
          L0_thetaProdPlane_US_hist[pT_bin][eta_bin]->Fill(thetaProdPlane);
          L0_thetaProdPlane_US_hist[nPtBins][eta_bin]->Fill(thetaProdPlane);
          L0_thetaProdPlane_US_hist[pT_bin][nEtaBins]->Fill(thetaProdPlane);
          L0_thetaProdPlane_US_hist[nPtBins][nEtaBins]->Fill(thetaProdPlane);

          L0_cosThetaProdPlane_US_hist[pT_bin][eta_bin]->Fill(cos(thetaProdPlane));
          L0_cosThetaProdPlane_US_hist[nPtBins][eta_bin]->Fill(cos(thetaProdPlane));
          L0_cosThetaProdPlane_US_hist[pT_bin][nEtaBins]->Fill(cos(thetaProdPlane));
          L0_cosThetaProdPlane_US_hist[nPtBins][nEtaBins]->Fill(cos(thetaProdPlane));

          L0_xF_US[eta_bin]->Fill(L_xF);
          L0_xF_US[nEtaBins]->Fill(L_xF);

        }

        if( pi_ch == 1)
        {
          L0bar_thetaProdPlane_US_hist[pT_bin][eta_bin]->Fill(thetaProdPlane);
          L0bar_thetaProdPlane_US_hist[nPtBins][eta_bin]->Fill(thetaProdPlane);
          L0bar_thetaProdPlane_US_hist[pT_bin][nEtaBins]->Fill(thetaProdPlane);
          L0bar_thetaProdPlane_US_hist[nPtBins][nEtaBins]->Fill(thetaProdPlane);

          L0bar_cosThetaProdPlane_US_hist[pT_bin][eta_bin]->Fill(cos(thetaProdPlane));
          L0bar_cosThetaProdPlane_US_hist[nPtBins][eta_bin]->Fill(cos(thetaProdPlane));
          L0bar_cosThetaProdPlane_US_hist[pT_bin][nEtaBins]->Fill(cos(thetaProdPlane));
          L0bar_cosThetaProdPlane_US_hist[nPtBins][nEtaBins]->Fill(cos(thetaProdPlane));

          L0bar_xF_US[eta_bin]->Fill(L_xF);
          L0bar_xF_US[nEtaBins]->Fill(L_xF);
        }



        if(eventId == eventID_last) //same event as in previous iteration
        {

          if( pi_ch == -1)
          {
            nLambda++;
          }
          else if( pi_ch == 1)
          {
            nLambdaBar++;
          }

        }
        else if(eventId != eventID_last) //new event
        {
          //cout<<"new event"<<endl;

          eventID_last = eventId;

          if(nLambda > 0 && nLambdaBar > 0)
          {
            nEventsWithLambdaPair_US_hist->Fill(1.5);
          }
          else
          {
            nEventsWithLambdaPair_US_hist->Fill(0.5);
          }

          nLambdasInEvent_US_hist->Fill(nLambda);
          nLambdaBarsInEvent_US_hist->Fill(nLambdaBar);

          if( pi_ch == -1)
          {
            nLambda = 1;
            nLambdaBar = 0;
          }
          else if( pi_ch == 1)
          {
            nLambda = 0;
            nLambdaBar = 1;
          }

        }
      }
    }
    else //unlike-sign combinations
    {

      if(L_mass > invMassRange_L0[0][pT_bin][eta_bin] && L_mass < invMassRange_L0[1][pT_bin][eta_bin])
      {


        if( pi_ch == 1 )
        {
          L0_thetaProdPlane_LS_hist[pT_bin][eta_bin]->Fill(thetaProdPlane);
          L0_thetaProdPlane_LS_hist[nPtBins][eta_bin]->Fill(thetaProdPlane);
          L0_thetaProdPlane_LS_hist[pT_bin][nEtaBins]->Fill(thetaProdPlane);
          L0_thetaProdPlane_LS_hist[nPtBins][nEtaBins]->Fill(thetaProdPlane);

          L0_cosThetaProdPlane_LS_hist[pT_bin][eta_bin]->Fill(cos(thetaProdPlane));
          L0_cosThetaProdPlane_LS_hist[nPtBins][eta_bin]->Fill(cos(thetaProdPlane));
          L0_cosThetaProdPlane_LS_hist[pT_bin][nEtaBins]->Fill(cos(thetaProdPlane));
          L0_cosThetaProdPlane_LS_hist[nPtBins][nEtaBins]->Fill(cos(thetaProdPlane));

          L0_xF_LS[eta_bin]->Fill(L_xF);
          L0_xF_LS[nEtaBins]->Fill(L_xF);
        }

        if( pi_ch == -1 )
        {
          L0bar_thetaProdPlane_LS_hist[pT_bin][eta_bin]->Fill(thetaProdPlane);
          L0bar_thetaProdPlane_LS_hist[nPtBins][eta_bin]->Fill(thetaProdPlane);
          L0bar_thetaProdPlane_LS_hist[pT_bin][nEtaBins]->Fill(thetaProdPlane);
          L0bar_thetaProdPlane_LS_hist[nPtBins][nEtaBins]->Fill(thetaProdPlane);

          L0bar_cosThetaProdPlane_LS_hist[pT_bin][eta_bin]->Fill(cos(thetaProdPlane));
          L0bar_cosThetaProdPlane_LS_hist[nPtBins][eta_bin]->Fill(cos(thetaProdPlane));
          L0bar_cosThetaProdPlane_LS_hist[pT_bin][nEtaBins]->Fill(cos(thetaProdPlane));
          L0bar_cosThetaProdPlane_LS_hist[nPtBins][nEtaBins]->Fill(cos(thetaProdPlane));

          L0bar_xF_LS[eta_bin]->Fill(L_xF);
          L0bar_xF_LS[nEtaBins]->Fill(L_xF);
        }


        if(eventId == eventID_last_background) //same event as in previous iteration
        {
          if( pi_ch == -1)
          {
            nLambda_background++;
          }
          else if( pi_ch == 1)
          {
            nLambdaBar_background++;
          }
        }
        else if(eventId != eventID_last_background) //new event
        {
          eventID_last_background = eventId;

          if(nLambda_background > 0 && nLambdaBar_background > 0)
          {
            nEventsWithLambdaPair_LS_hist->Fill(1.5);
          }
          else
          {
            nEventsWithLambdaPair_LS_hist->Fill(0.5);
          }

          nLambdasInEvent_LS_hist->Fill(nLambda_background);
          nLambdaBarsInEvent_LS_hist->Fill(nLambdaBar_background);

          if( pi_ch == -1)
          {
            nLambda_background = 1;
            nLambdaBar_background = 0;
          }
          else if( pi_ch == 1)
          {
            nLambda_background = 0;
            nLambdaBar_background = 1;
          }

        }

      }//end if L_mass
    }//end else for flag_int check
  }//end loop over entries in NTuple

  //________________________________________________________________________________________________________


  //----------------------------------Histograms, pT binnig, main analyisis---------------------------------


  TCanvas *Lambda_stats_can = new TCanvas("Lambda_stats_can", "Lambda_stats_can", 1200, 1000);

  Lambda_stats_can->cd();
  nLambdasInEvent_US_hist->Sumw2();
  nLambdasInEvent_US_hist->Add(nLambdasInEvent_LS_hist, -1);
  nLambdasInEvent_US_hist->Draw();

  Lambda_stats_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/nLambdasInEvent.png");


  TCanvas *LambdaBar_stats_can = new TCanvas("LambdaBar_stats_can", "LambdaBar_stats_can", 1200, 1000);

  LambdaBar_stats_can->cd();
  nLambdaBarsInEvent_US_hist->Sumw2();
  nLambdaBarsInEvent_US_hist->Add(nLambdaBarsInEvent_LS_hist, -1);
  nLambdaBarsInEvent_US_hist->Draw();

  LambdaBar_stats_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/nLambdaBarsInEvent.png");


  TCanvas *nEventsWithLambdaPair_can = new TCanvas("nEventsWithLambdaPair_can", "nEventsWithLambdaPair_can", 1200, 1000);

  nEventsWithLambdaPair_can->cd();
  nEventsWithLambdaPair_US_hist->Sumw2();
  nEventsWithLambdaPair_US_hist->Add(nEventsWithLambdaPair_LS_hist, -1);
  nEventsWithLambdaPair_US_hist->Draw();

  nEventsWithLambdaPair_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/nEventsWithLambdaPair.png");


  //analysis hostograms in pT and eta bins
  TCanvas *L0_thetaProdPlane_US_LS_can[nPtBins+1][nEtaBins+1];
  TCanvas *L0_cosThetaProdPlane_US_LS_can[nPtBins+1][nEtaBins+1];

  TCanvas *L0bar_thetaProdPlane_US_LS_can[nPtBins+1][nEtaBins+1];
  TCanvas *L0bar_cosThetaProdPlane_US_LS_can[nPtBins+1][nEtaBins+1];

  TCanvas *L0_polarization_can[nEtaBins+1];
  TCanvas *L0bar_polarization_can[nEtaBins+1];

  TCanvas *L0_xF_can[nEtaBins+1];
  TCanvas *L0bar_xF_can[nEtaBins+1];


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
  {

    L0_polarization_can[etaBin] = new TCanvas(Form("L0_polarization_can_eta_%i", etaBin), Form("L0_polarization_can_eta_%i", etaBin), 1200, 1000 );
    L0bar_polarization_can[etaBin] = new TCanvas(Form("L0bar_polarization_can_eta_%i", etaBin), Form("L0bar_polarization_can_eta_%i", etaBin), 1200, 1000 );

    L0_xF_can[etaBin] = new TCanvas(Form("L0_xF_can_eta_%i", etaBin), Form("L0_xF_can_eta_%i", etaBin), 1200, 1000);
    L0bar_xF_can[etaBin] = new TCanvas(Form("L0bar_xF_can_eta_%i", etaBin), Form("L0bar_xF_can_eta_%i", etaBin), 1200, 1000);

    TString *eta_range = new TString();
    if(etaBin < nEtaBins) eta_range->Form("%0.1f < #eta < %0.1f", eta_bins[etaBin], eta_bins[etaBin+1]);
    else eta_range->Form("-1 < #eta < 1");


    for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
    {
      if( energy == 510  && pTbin == 0 ) continue;
      if( energy == 200 && strictTOF_cut == 2 && pTbin == 0 ) continue;

      L0_thetaProdPlane_US_LS_can[pTbin][etaBin] = new TCanvas(Form("L0_thetaProdPlane_US_LS_can_pt_%i_eta_%i", pTbin, etaBin), Form("L0_thetaProdPlane_US_LS_can_pt_%i_eta_%i", pTbin, etaBin), 1200, 1000);
      L0_cosThetaProdPlane_US_LS_can[pTbin][etaBin] = new TCanvas(Form("L0_cosThetaProdPlane_US_LS_can_pT_%i_eta_%i", pTbin, etaBin), Form("L0_cosThetaProdPlane_US_LS_can_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);

      L0bar_thetaProdPlane_US_LS_can[pTbin][etaBin] = new TCanvas(Form("L0bar_thetaProdPlane_US_LS_can_pt_%i_eta_%i", pTbin, etaBin), Form("L0bar_thetaProdPlane_US_LS_can_pt_%i_eta_%i", pTbin, etaBin), 1200, 1000);
      L0bar_cosThetaProdPlane_US_LS_can[pTbin][etaBin] = new TCanvas(Form("L0bar_cosThetaProdPlane_US_LS_can_pT_%i_eta_%i", pTbin, etaBin), Form("L0bar_cosThetaProdPlane_US_LS_can_pT_%i_eta_%i", pTbin, etaBin), 1200, 1000);


      TString *pT_range = new TString();
      if(pTbin < nPtBins) pT_range->Form("%0.1f < #it{p}_{T} < %0.1f GeV/#it{c}", pT_bins[pTbin], pT_bins[pTbin+1]);
      else pT_range->Form("p_{T} integrated");



      TLegend *SigAndBckgLeg_theta = new TLegend(0.2, 0.65, 0.4, 0.89 );
      SigAndBckgLeg_theta->AddEntry(L0_thetaProdPlane_US_hist[pTbin][etaBin], "Unlike-sign");
      SigAndBckgLeg_theta->AddEntry(L0_thetaProdPlane_LS_hist[pTbin][etaBin], "Like-sign");
      SigAndBckgLeg_theta->SetBorderSize(0);
      //SigAndBckgLeg->Draw("same");


      ProdPlaneOutFile->cd();

      L0_thetaProdPlane_US_LS_can[pTbin][etaBin]->cd();

      L0_thetaProdPlane_US_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0_thetaProdPlane_US_hist[pTbin][etaBin]->SetMarkerColor(kRed);
      L0_thetaProdPlane_US_hist[pTbin][etaBin]->SetLineColor(kRed);
      L0_thetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->SetTitle("#theta*");
      L0_thetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->CenterTitle();
      L0_thetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->SetTitle("d#{N}/d#theta*");
      L0_thetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0_thetaProdPlane_US_hist[pTbin][etaBin]->Sumw2();
      L0_thetaProdPlane_US_hist[pTbin][etaBin]->Scale(1./L0_thetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(1)); //scale by bin width to obtain d N/d theta*

      L0_thetaProdPlane_LS_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0_thetaProdPlane_LS_hist[pTbin][etaBin]->SetMarkerColor(kBlue);
      L0_thetaProdPlane_LS_hist[pTbin][etaBin]->SetLineColor(kBlue);
      L0_thetaProdPlane_LS_hist[pTbin][etaBin]->Sumw2();
      L0_thetaProdPlane_LS_hist[pTbin][etaBin]->Scale(1./L0_thetaProdPlane_LS_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(1));

      //L0_thetaProdPlane_LS_hist[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 7000);


      L0_thetaProdPlane_US_hist[pTbin][etaBin]->Draw("p e");
      L0_thetaProdPlane_LS_hist[pTbin][etaBin]->Draw("p e same");


      TPaveText *L0_text = new TPaveText(0.35, 0.75, 0.75, 0.89, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_text->AddText("Minimum bias");
      L0_text->AddText("#Lambda^{0}");
      L0_text->AddText(eta_range->Data());
      L0_text->AddText(pT_range->Data());
      L0_text->SetFillColorAlpha(0, 0.01);
      L0_text->Draw("same");


      L0_thetaProdPlane_US_LS_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0_thetaProdPlane_US_LS_can_pT_%i_eta_%i.png", pTbin, etaBin));
      //____________________________________________________________________________________________________________________________________________________________________________________________________________



      L0_cosThetaProdPlane_US_LS_can[pTbin][etaBin]->cd();

      L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->SetMarkerColor(kRed);
      L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->SetLineColor(kRed);
      L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->SetTitle("cos(#theta*)");
      L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->CenterTitle();
      L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->Sumw2();
      L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->Scale(1./L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(1)); //scale by vin width to obtain d N/d cos(theta*)
      L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->Divide(L0_cosThetaProdPlane_eff_hist[pTbin][etaBin]);
      L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->SetMinimum(0);

      L0_cosThetaProdPlane_LS_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0_cosThetaProdPlane_LS_hist[pTbin][etaBin]->SetMarkerColor(kBlue);
      L0_cosThetaProdPlane_LS_hist[pTbin][etaBin]->SetLineColor(kBlue);
      L0_cosThetaProdPlane_LS_hist[pTbin][etaBin]->Sumw2();
      L0_cosThetaProdPlane_LS_hist[pTbin][etaBin]->Scale(1./L0_cosThetaProdPlane_LS_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(1));
      //L0_cosThetaProdPlane_LS_hist[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 7000);

      L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->Draw("p e");
      L0_cosThetaProdPlane_LS_hist[pTbin][etaBin]->Draw("p e same");



      L0_text->Draw("same");


      TFitResultPtr fit_res_L0_US_theta_star;

      TF1 *fitL0_US_ThetaStar = new TF1("fitL0_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      fit_res_L0_US_theta_star = L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->Fit(fitL0_US_ThetaStar, "s i 0 r");

      fitL0_US_ThetaStar->SetLineColor(1);
      fitL0_US_ThetaStar->Draw("same");


      TFitResultPtr fit_res_L0_LS_theta_star;

      TF1 *fitL0_LS_ThetaStar = new TF1("fitL0_LS_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitL0_LS_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaLS_wrong_sign = L_inv_mass_LS->Fit(fitGaLSsBack, "s i 0", "", 1.07, 1.4);
      fit_res_L0_LS_theta_star = L0_cosThetaProdPlane_LS_hist[pTbin][etaBin]->Fit(fitL0_LS_ThetaStar, "s i 0 r");

      fitL0_LS_ThetaStar->SetLineColor(1);
      fitL0_LS_ThetaStar->Draw("same");


      L0_cosThetaProdPlane_US_LS_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0_cosThetaProdPlane_US_LS_can_pT_%i_eta_%i.png", pTbin, etaBin));

      //fill polarization histogram




      double L0_pol_tot = fit_res_L0_US_theta_star->Parameter(1)/L0_alpha;
      double L0_pol_back = fit_res_L0_LS_theta_star->Parameter(1)/L0_alpha;

      double L0_pol_sig = ( L0_pol_tot - ( 1 - L0_kappa[pTbin][etaBin] )*L0_pol_back )/L0_kappa[pTbin][etaBin];

      L0_polarization[etaBin]->SetBinContent( pTbin+1, L0_pol_sig );


      double L0_pol_tot_err = fit_res_L0_US_theta_star->Error(1)/L0_alpha/L0_kappa[pTbin][etaBin];
      double L0_pol_back_err = fit_res_L0_LS_theta_star->Error(1)/L0_alpha * (1-L0_kappa[pTbin][etaBin])/L0_kappa[pTbin][etaBin];

      double L0_pol_sig_err = sqrt( L0_pol_tot_err*L0_pol_tot_err + L0_pol_back_err*L0_pol_back_err);

      L0_polarization[etaBin]->SetBinError( pTbin+1, L0_pol_sig_err );







      //_________________________________________________________________________________________________________________________________________________



      L0bar_thetaProdPlane_US_LS_can[pTbin][etaBin]->cd();

      L0bar_thetaProdPlane_US_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0bar_thetaProdPlane_US_hist[pTbin][etaBin]->SetMarkerColor(kRed);
      L0bar_thetaProdPlane_US_hist[pTbin][etaBin]->SetLineColor(kRed);
      L0bar_thetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->SetTitle("#theta*");
      L0bar_thetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->CenterTitle();
      L0bar_thetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->SetTitle("d#{N}/d#theta*");
      L0bar_thetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0bar_thetaProdPlane_US_hist[pTbin][etaBin]->Sumw2();
      L0bar_thetaProdPlane_US_hist[pTbin][etaBin]->Scale(1./L0bar_thetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(1)); //scale by bin width to obtain d N/d theta*

      L0bar_thetaProdPlane_LS_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0bar_thetaProdPlane_LS_hist[pTbin][etaBin]->SetMarkerColor(kBlue);
      L0bar_thetaProdPlane_LS_hist[pTbin][etaBin]->SetLineColor(kBlue);
      L0bar_thetaProdPlane_LS_hist[pTbin][etaBin]->Sumw2(kBlue);
      L0bar_thetaProdPlane_LS_hist[pTbin][etaBin]->Scale(1./L0bar_thetaProdPlane_LS_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(1));
      //L0bar_thetaProdPlane_LS_hist[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 7000);

      L0bar_thetaProdPlane_US_hist[pTbin][etaBin]->Draw("p e");
      L0bar_thetaProdPlane_LS_hist[pTbin][etaBin]->Draw("p e same");



      TPaveText *L0bar_text = new TPaveText(0.35, 0.75, 0.75, 0.89, "NDC");
      //cent_text->AddText("STAR preliminary");
      //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_text->AddText("Minimum bias");
      L0bar_text->AddText("#bar{#Lambda^{0}}");
      L0bar_text->AddText(eta_range->Data());
      L0bar_text->AddText(pT_range->Data());
      L0bar_text->SetFillColorAlpha(0, 0.01);
      L0bar_text->Draw("same");

      L0bar_thetaProdPlane_US_LS_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0bar_thetaProdPlane_US_LS_can_pT_%i_eta_%i.png", pTbin, etaBin));




      L0bar_cosThetaProdPlane_US_LS_can[pTbin][etaBin]->cd();

      L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->SetMarkerColor(kRed);
      L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->SetLineColor(kRed);
      L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->SetTitle("cos(#theta*)");
      L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->CenterTitle();
      L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->Sumw2();
      L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->Scale(1./L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(1));
      L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->Divide(L0bar_cosThetaProdPlane_eff_hist[pTbin][etaBin]);
      L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->SetMinimum(0);
      //L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 7000);

      L0bar_cosThetaProdPlane_LS_hist[pTbin][etaBin]->SetMarkerStyle(20);
      L0bar_cosThetaProdPlane_LS_hist[pTbin][etaBin]->SetMarkerColor(kBlue);
      L0bar_cosThetaProdPlane_LS_hist[pTbin][etaBin]->SetLineColor(kBlue);
      L0bar_cosThetaProdPlane_LS_hist[pTbin][etaBin]->Sumw2();
      L0bar_cosThetaProdPlane_LS_hist[pTbin][etaBin]->Scale(1./L0bar_cosThetaProdPlane_LS_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(1));

      L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->Draw("p e");
      L0bar_cosThetaProdPlane_LS_hist[pTbin][etaBin]->Draw("p e same");



      L0bar_text->Draw("same");


      TFitResultPtr fit_res_L0bar_US_theta_star;

      TF1 *fitL0bar_US_ThetaStar = new TF1("fitL0bar_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_US_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      fit_res_L0bar_US_theta_star = L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->Fit(fitL0bar_US_ThetaStar, "s i 0 r");

      fitL0bar_US_ThetaStar->SetLineColor(1);
      fitL0bar_US_ThetaStar->Draw("same");


      TFitResultPtr fit_res_L0bar_LS_theta_star;

      TF1 *fitL0bar_LS_ThetaStar = new TF1("fitL0bar_LS_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_LS_ThetaStar->SetParameters(100, 0.5);

      //fit_res_gaLS_wrong_sign = L_inv_mass_LS->Fit(fitGaLSsBack, "s i 0", "", 1.07, 1.4);
      fit_res_L0bar_LS_theta_star = L0bar_cosThetaProdPlane_LS_hist[pTbin][etaBin]->Fit(fitL0bar_LS_ThetaStar, "s i 0 r");

      fitL0bar_LS_ThetaStar->SetLineColor(1);
      fitL0bar_LS_ThetaStar->Draw("same");


      L0bar_cosThetaProdPlane_US_LS_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0bar_cosThetaProdPlane_US_LS_can_pT_%i_eta_%i.png", pTbin, etaBin));


      //L0bar_polarization[etaBin]->SetBinContent( pTbin+1, fit_res_L0bar_US_theta_star->Parameter(1)/L0bar_alpha );
      //L0bar_polarization[etaBin]->SetBinError( pTbin+1, fit_res_L0bar_US_theta_star->Error(1)/L0bar_alpha );


      double L0bar_pol_tot = fit_res_L0bar_US_theta_star->Parameter(1)/L0bar_alpha;
      double L0bar_pol_back = fit_res_L0bar_LS_theta_star->Parameter(1)/L0bar_alpha;

      double L0bar_pol_sig = ( L0bar_pol_tot - ( 1 - L0bar_kappa[pTbin][etaBin] )*L0bar_pol_back )/L0bar_kappa[pTbin][etaBin];

      L0bar_polarization[etaBin]->SetBinContent( pTbin+1, L0bar_pol_sig );


      double L0bar_pol_tot_err = fit_res_L0bar_US_theta_star->Error(1)/L0bar_alpha/L0bar_kappa[pTbin][etaBin];
      double L0bar_pol_back_err = fit_res_L0bar_LS_theta_star->Error(1)/L0bar_alpha * (1-L0bar_kappa[pTbin][etaBin])/L0bar_kappa[pTbin][etaBin];

      double L0bar_pol_sig_err = sqrt( L0bar_pol_tot_err*L0bar_pol_tot_err + L0bar_pol_back_err*L0bar_pol_back_err);

      L0bar_polarization[etaBin]->SetBinError( pTbin+1, L0bar_pol_sig_err );





      //save histograms to multi-page PDF
      if(pTbin == 0 && etaBin == 0)
      {
        L0_cosThetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0_cosThetaProdPlane_US_LS.pdf(", "pdf");
        L0_thetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0_thetaProdPlane_US_LS.pdf(", "pdf");

        L0bar_cosThetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0bar_cosThetaProdPlane_US_LS.pdf(", "pdf");
        L0bar_thetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0bar_thetaProdPlane_US_LS.pdf(", "pdf");



      }
      else if( pTbin == nPtBins && etaBin == nEtaBins)
      {
        L0_cosThetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0_cosThetaProdPlane_US_LS.pdf)", "pdf");
        L0_thetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0_thetaProdPlane_US_LS.pdf)", "pdf");

        L0bar_cosThetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0bar_cosThetaProdPlane_US_LS.pdf)", "pdf");
        L0bar_thetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0bar_thetaProdPlane_US_LS.pdf)", "pdf");


      }
      else
      {
        L0_cosThetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0_cosThetaProdPlane_US_LS.pdf", "pdf");
        L0_thetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0_thetaProdPlane_US_LS.pdf", "pdf");

        L0bar_cosThetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0bar_cosThetaProdPlane_US_LS.pdf", "pdf");
        L0bar_thetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0bar_thetaProdPlane_US_LS.pdf", "pdf");
      }


      if(ReadMode == 0)
      {
        L0_thetaProdPlane_US_hist[pTbin][etaBin]->Write();
        L0_thetaProdPlane_LS_hist[pTbin][etaBin]->Write();

        L0_cosThetaProdPlane_US_hist[pTbin][etaBin]->Write();
        L0_cosThetaProdPlane_LS_hist[pTbin][etaBin]->Write();

        L0bar_thetaProdPlane_US_hist[pTbin][etaBin]->Write();
        L0bar_thetaProdPlane_LS_hist[pTbin][etaBin]->Write();

        L0bar_cosThetaProdPlane_US_hist[pTbin][etaBin]->Write();
        L0bar_cosThetaProdPlane_LS_hist[pTbin][etaBin]->Write();


      }

    }//end for pTbin

    L0_polarization_can[etaBin]->cd();

    L0_polarization[etaBin]->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
    L0_polarization[etaBin]->GetXaxis()->CenterTitle();
    L0_polarization[etaBin]->GetYaxis()->SetTitle("P_{#Lambda}");
    L0_polarization[etaBin]->GetYaxis()->CenterTitle();
    L0_polarization[etaBin]->SetMarkerStyle(20);
    L0_polarization[etaBin]->SetMarkerColor(kRed);
    L0_polarization[etaBin]->SetLineColor(kRed);
    L0_polarization[etaBin]->Draw("p e");

    TPaveText *L0_pol_text = new TPaveText(0.35, 0.75, 0.75, 0.89, "NDC");
    //cent_text->AddText("STAR preliminary");
    //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
    L0_pol_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0_pol_text->AddText("Minimum bias");
    L0_pol_text->AddText("#Lambda^{0}");
    L0_pol_text->AddText(eta_range->Data());
    L0_pol_text->SetFillColorAlpha(0, 0.01);
    L0_pol_text->Draw("same");

    L0_polarization_can[etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0_polarization_eta_%i.png", etaBin));
    L0_polarization_can[etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0_polarization_eta_%i.pdf", etaBin));



    L0bar_polarization_can[etaBin]->cd();

    L0bar_polarization[etaBin]->GetXaxis()->SetTitle("p_{T} (GeV/#it{c})");
    L0bar_polarization[etaBin]->GetXaxis()->CenterTitle();
    L0bar_polarization[etaBin]->GetYaxis()->SetTitle("P_{#bar{#Lambda}}");
    L0bar_polarization[etaBin]->GetYaxis()->CenterTitle();
    L0bar_polarization[etaBin]->SetMarkerStyle(20);
    L0bar_polarization[etaBin]->SetMarkerColor(kRed);
    L0bar_polarization[etaBin]->SetLineColor(kRed);
    L0bar_polarization[etaBin]->Draw("p e");


    TPaveText *L0bar_pol_text = new TPaveText(0.35, 0.75, 0.75, 0.89, "NDC");
    //cent_text->AddText("STAR preliminary");
    //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
    L0bar_pol_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0bar_pol_text->AddText("Minimum bias");
    L0bar_pol_text->AddText("#bar{#Lambda^{0}}");
    L0bar_pol_text->AddText(eta_range->Data());
    L0bar_pol_text->SetFillColorAlpha(0, 0.01);
    L0bar_pol_text->Draw("same");

    L0bar_polarization_can[etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0bar_polarization_eta_%i.png", etaBin));
    L0bar_polarization_can[etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/thetaProdPlane/L0bar_polarization_eta_%i.pdf", etaBin));

    //____________________________________________________________________________________

    L0_xF_can[etaBin]->cd();

    //gPad->SetLogy();

    L0_xF_US[etaBin]->Add(L0_xF_LS[etaBin]);
    L0_xF_US[etaBin]->GetXaxis()->SetTitle("x_{F}");
    L0_xF_US[etaBin]->GetXaxis()->CenterTitle();
    L0_xF_US[etaBin]->GetYaxis()->SetTitle("Counts");
    L0_xF_US[etaBin]->GetYaxis()->CenterTitle();
    L0_xF_US[etaBin]->SetMarkerStyle(20);
    L0_xF_US[etaBin]->SetMarkerColor(kRed);
    L0_xF_US[etaBin]->SetLineColor(kRed);
    L0_xF_US[etaBin]->Draw("p e");


    L0_pol_text->Draw("same");

    L0_xF_can[etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0_xF_eta_%i.png", etaBin));


    L0bar_xF_can[etaBin]->cd();

    //gPad->SetLogy();

    L0bar_xF_US[etaBin]->Add(L0bar_xF_LS[etaBin]);
    L0bar_xF_US[etaBin]->GetXaxis()->SetTitle("x_{F}");
    L0bar_xF_US[etaBin]->GetXaxis()->CenterTitle();
    L0bar_xF_US[etaBin]->GetYaxis()->SetTitle("Counts");
    L0bar_xF_US[etaBin]->GetYaxis()->CenterTitle();
    L0bar_xF_US[etaBin]->SetMarkerStyle(20);
    L0bar_xF_US[etaBin]->SetMarkerColor(kRed);
    L0bar_xF_US[etaBin]->SetLineColor(kRed);
    L0bar_xF_US[etaBin]->Draw("p e");

    L0bar_pol_text->Draw("same");


    L0bar_xF_can[etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/L0bar_xF_eta_%i.png", etaBin));



  }//end for etaBin

  cout<<"N events with L-Lbar pair: "<<nEventsWithLambdaPair_US_hist->GetBinContent(2)<<endl;

  cout<<"N events with 2+ L in event: "<<nLambdasInEvent_US_hist->Integral(3, 10)<<endl;
  cout<<"N events with 2+ Lbar in event: "<<nLambdaBarsInEvent_US_hist->Integral(3, 10)<<endl;

  ProdPlaneOutFile->Close();


  return;

}
//__________________________________________________________________________________________________________________________

//ReadMode = 0 - read TTree, ReadMode = 1 - read histograms - First run in ReadMode = 0 to save relevant histograms, then can run in ReadMode = 1 to read just histograms and save time
//energy - collision energy in GeV
void Ana001_Lambda(const int ReadMode = 0, const int energy = 510, const int year = 2017)
{
  ifstream fileList;

  if(energy == 510) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run17/fileList.list");
  else if(energy == 200) fileList.open("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Run12/fileList.list");
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


  bool invMassFinish = InvMass(myChain, invMassRange_L, invMassRange_L0, invMassRange_L0bar, L0_kappa, L0bar_kappa, ReadMode, energy , year);

  if(!invMassFinish)
  {
    cout<<"Analysis of invariant spectra ended abnormally. Abborting!"<<endl;

    return;
  }

  LambdaPolarization(myChain, invMassRange_L0, invMassRange_L0bar, L0_kappa, L0bar_kappa, ReadMode, energy , year);

  cout<<endl;
  cout<<"Nubmer of accepted events: "<<hEventStat1->GetBinContent(6)<<endl;


  return;
}
