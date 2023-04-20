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
#include"TLatex.h"


using namespace std;



const int nPtBins = 8;
float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};

const int nEtaBins = 3;
float const eta_bins[nEtaBins+1] = { -1, -0.2, 0.2, 1 };

const float K0s_mass_PDG = 0.497611; //PDG mass in GeV

const float y_cut = 1;
const int strictTOF_cut = 2; //0 - hybrid TOF for both daughters, 1 - strict TOF for pions, 2 - strict TOF for both pion and proton
const float cos_theta_cut = 0.996;
const float decayL_cut = 25;


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

bool cuts(int pi_hasTOFinfo, int p_hasTOFinfo, float y, float theta, float decayL)
{

  if( !( TMath::Abs(y) < y_cut ) ) return false;
  if( strictTOF_cut == 1 && pi_hasTOFinfo == 0 ) return false; //TOF matched pions
  if( strictTOF_cut == 2 && (pi_hasTOFinfo == 0 || p_hasTOFinfo == 0) ) return false; //TOF matched pions and protons
  if(cos(theta) < cos_theta_cut) return false;
  if(decayL > decayL_cut) return false;

  return true;

}


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
  Float_t K0s_mass, K0s_pt, K0s_eta, K0s_decayL, K0s_theta, K0s_DCAdaughters;


  Float_t pi1_pt, pi2_pt;
  Float_t pi1_eta, pi2_eta;
  Float_t pi1_phi, pi2_phi;
  //Float_t pi1_ch;
  Int_t pi1_hasTOFinfo, pi2_hasTOFinfo;


  Float_t thetaProdPlane;


  //---------------SET BARANCH ADDRESSES------------------------
  K0s_tree->SetBranchAddress("pair_charge", &charge);
  K0s_tree->SetBranchAddress("pair_mass", &K0s_mass);
  K0s_tree->SetBranchAddress("pair_pt", &K0s_pt);
  K0s_tree->SetBranchAddress("pair_eta", &K0s_eta);
  K0s_tree->SetBranchAddress("pair_decayL", &K0s_decayL);
  K0s_tree->SetBranchAddress("pair_theta", &K0s_theta);
  K0s_tree->SetBranchAddress("pair_DCAdaughters", &K0s_DCAdaughters);

  K0s_tree->SetBranchAddress("p1_pt", &pi1_pt);
  K0s_tree->SetBranchAddress("p1_eta", &pi1_eta);
  K0s_tree->SetBranchAddress("p1_phi", &pi1_phi);
  //K0s_tree->SetBranchAddress("p1_ch", &pi1_ch);
  K0s_tree->SetBranchAddress("p1_hasTOFinfo", &pi1_hasTOFinfo);

  K0s_tree->SetBranchAddress("p2_pt", &pi2_pt);
  K0s_tree->SetBranchAddress("p2_eta", &pi2_eta);
  K0s_tree->SetBranchAddress("p2_phi", &pi2_phi);
  K0s_tree->SetBranchAddress("p2_hasTOFinfo", &pi2_hasTOFinfo);

  K0s_tree->SetBranchAddress("thetaProdPlane", &thetaProdPlane);

  //--------------------------------------------------------------------------


  TH1D *K0s_inv_mass_US[nPtBins+1][nEtaBins+1];
  TH1D *K0s_inv_mass_LS[nPtBins+1][nEtaBins+1];
  TH1D *K0s_inv_mass[nPtBins+1][nEtaBins+1];



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
    if(i%100000 == 0)
    {
      cout<<i<<endl;
    }

    K0s_tree->GetEntry(i);



    double K0s_y = rapidity(K0s_pt, K0s_eta, K0s_mass_PDG);

    //------------------------------------------------------------------------------------------------------------------

    //cuts
    if( !cuts(pi1_hasTOFinfo, pi2_hasTOFinfo, K0s_y, K0s_theta, K0s_decayL) ) continue;

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

    //find eta bin of Lambda
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


  TCanvas *InvMassCan_US_LS[nPtBins+1][nEtaBins+1];
  TCanvas *InvMassCan[nPtBins+1][nEtaBins+1];

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  TFitResultPtr fit_res_gaus;

  for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
  {
    for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
    {

      InvMassCan_US_LS[pTbin][etaBin] = new TCanvas(Form("InvMassCan_US_LS_pT_%i_eta_%i_can", pTbin, etaBin), Form("InvMassCan_US_LS_pT_%i_eta_%i_can", pTbin, etaBin), 1200, 1000);
      InvMassCan[pTbin][etaBin] = new TCanvas(Form("InvMass_pT_%i_eta_%i_Can", pTbin, etaBin), Form("InvMass_pT_%i_eta_%i_Can", pTbin, etaBin), 1200, 1000);

      TString *pT_range = new TString();
      if(pTbin < nPtBins) pT_range->Form("%0.1f < #it{p}_{T} < %0.1f GeV/c", pT_bins[pTbin], pT_bins[pTbin+1]);
      else pT_range->Form("p_{T} integrated");

      TString *eta_range = new TString();
      if(etaBin < nEtaBins) eta_range->Form("%0.1f < #eta < %0.1f", eta_bins[etaBin], eta_bins[etaBin+1]);
      else eta_range->Form("-1 < #eta < 1");


      TPaveText *text = new TPaveText(0.35, 0.75, 0.75, 0.89, "NDC");
      //text[pTbin][etaBin]->AddText("STAR preliminary");
      //((TText*)text[pTbin][etaBin]->GetListOfLines()[pTbin][etaBin]->Last())[pTbin][etaBin]->SetTextColor(2);
      text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      text->AddText("Minimum bias");
      text->AddText("K_{s}^{0} and #bar{K_{s}^{0}}");
      text->AddText(eta_range->Data());
      text->AddText(pT_range->Data());
      text->SetFillColorAlpha(0, 0.01);





      InvMassCan_US_LS[pTbin][etaBin]->cd();

      K0s_inv_mass_US[pTbin][etaBin]->SetMarkerStyle(20);
      K0s_inv_mass_US[pTbin][etaBin]->SetMarkerColor(kRed);
      K0s_inv_mass_US[pTbin][etaBin]->SetLineColor(kRed);
      K0s_inv_mass_US[pTbin][etaBin]->GetXaxis()->SetTitle("M_{inv}^{p#pi} (GeV/#it{c}^{2})");
      K0s_inv_mass_US[pTbin][etaBin]->GetXaxis()->CenterTitle();
      K0s_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
      K0s_inv_mass_US[pTbin][etaBin]->GetYaxis()->CenterTitle();
      //K0s_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 7000);
      //K0s_inv_mass_US[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 6000);
      K0s_inv_mass_US[pTbin][etaBin]->SetMinimum(0);
      K0s_inv_mass_US[pTbin][etaBin]->Draw("p e");

      K0s_inv_mass_LS[pTbin][etaBin]->SetMarkerStyle(20);
      K0s_inv_mass_LS[pTbin][etaBin]->SetMarkerColor(kBlue);
      K0s_inv_mass_LS[pTbin][etaBin]->SetLineColor(kBlue);
      K0s_inv_mass_LS[pTbin][etaBin]->Draw("p e same");

      TLegend *SigAndBckgLeg = new TLegend(0.2, 0.65, 0.4, 0.89 );
      SigAndBckgLeg->AddEntry(K0s_inv_mass_US[pTbin][etaBin], "Unlike-sign");
      SigAndBckgLeg->AddEntry(K0s_inv_mass_LS[pTbin][etaBin], "Like-sign");
      SigAndBckgLeg->SetBorderSize(0);
      SigAndBckgLeg->Draw("same");

      TPaveText *cent_text = new TPaveText(0.6, 0.6, 0.89, 0.75, "NDC");
      //cent_text[pTbin][etaBin]->AddText("STAR preliminary");
      //((TText*)cent_text[pTbin][etaBin]->GetListOfLines()[pTbin][etaBin]->Last())[pTbin][etaBin]->SetTextColor(2);
      cent_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      cent_text->AddText("Minimum bias");
      cent_text->AddText("K_{s}^{0} and #bar{K_{s}^{0}}");
      cent_text->AddText(eta_range->Data());
      cent_text->AddText(pT_range->Data());
      cent_text->SetFillColorAlpha(0, 0.01);
      cent_text->Draw("same");

      InvMassCan_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/K0s_inv_mass_US_LS_pT_%i_eta_%i.png", pTbin, etaBin));
      //InvMassCan_US_LS[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/K0s_inv_mass_US_LS_pT_%i_eta_%i.pdf", pTbin, etaBin));




      InvMassCan[pTbin][etaBin]->cd();

      K0s_inv_mass[pTbin][etaBin] = (TH1D*)K0s_inv_mass_US[pTbin][etaBin]->Clone();
      K0s_inv_mass[pTbin][etaBin]->Sumw2();
      K0s_inv_mass_LS[pTbin][etaBin]->Sumw2();
      K0s_inv_mass[pTbin][etaBin]->Add(K0s_inv_mass_LS[pTbin][etaBin], -1);
      K0s_inv_mass[pTbin][etaBin]->GetXaxis()->SetTitle("M_{inv}^{p#pi} (GeV/#it{c}^{2})");
      K0s_inv_mass[pTbin][etaBin]->GetXaxis()->CenterTitle();
      K0s_inv_mass[pTbin][etaBin]->GetYaxis()->SetTitle("Counts");
      K0s_inv_mass[pTbin][etaBin]->GetYaxis()->CenterTitle();
      K0s_inv_mass[pTbin][etaBin]->SetMarkerStyle(20);
      K0s_inv_mass[pTbin][etaBin]->SetMarkerColor(kRed);
      K0s_inv_mass[pTbin][etaBin]->SetLineColor(kRed);
      K0s_inv_mass[pTbin][etaBin]->Draw("p e");

      TF1 *zeroLine = new TF1("zeroLine", "[0]", 1, 1.2);
      zeroLine->SetParameter(0, 0);
      zeroLine->SetLineColor(1);
      zeroLine->SetLineStyle(9);
      zeroLine->Draw("same");

      TLegend *SigMinusBckgLeg = new TLegend(0.2, 0.65, 0.4, 0.89 );
      SigMinusBckgLeg->AddEntry(K0s_inv_mass[pTbin][etaBin], "US-LS");
      SigMinusBckgLeg->SetBorderSize(0);
      SigMinusBckgLeg->Draw("same");

      cent_text->Draw("same");





      TF1 *fitGauss = new TF1("fitGauss", "gaus(0)", 0.45, 0.55);

      if( pTbin == 0 ) fitGauss->SetParameters(10000, 1.116, 0.002);
      else
      {
        //cout<<"Valid"<<endl;
        fitGauss->SetParameters(fit_res_gaus->Parameter(0), fit_res_gaus->Parameter(1), fit_res_gaus->Parameter(2));
      }


      //fit_res_gaus_wrong_sign = K0s_inv_mass_US[pTbin][etaBin]->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
      fit_res_gaus = K0s_inv_mass_US[pTbin][etaBin]->Fit(fitGauss, "s i l 0", "", 0.49, 0.51);

      fitGauss->SetLineColor(1);
      fitGauss->Draw("same");


      if(!fit_res_gaus->IsValid())
      {
        cout<<"Fit not valid for pT "<<pTbin<<" eta "<<etaBin<<endl;
        return false;
      }


      float peak_range_min = fit_res_gaus->Parameter(1) - 3*fit_res_gaus->Parameter(2);
      float peak_range_max = fit_res_gaus->Parameter(1) + 3*fit_res_gaus->Parameter(2);

      invMassRange[0][pTbin][etaBin] = peak_range_min;
      invMassRange[1][pTbin][etaBin] = peak_range_max;

      int peak_range_min_bin = K0s_inv_mass[pTbin][etaBin]->GetXaxis()->FindBin(peak_range_min);
      int peak_range_max_bin = K0s_inv_mass[pTbin][etaBin]->GetXaxis()->FindBin(peak_range_max);

      kappa[pTbin][etaBin] = K0s_inv_mass[pTbin][etaBin]->Integral(peak_range_min_bin, peak_range_max_bin)/K0s_inv_mass[pTbin][etaBin]->Integral(peak_range_min_bin, peak_range_max_bin);



      InvMassCan[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/K0s_inv_mass_signal_pT_%i_eta_%i.png", pTbin, etaBin));
      //InvMassCan[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/K0s_inv_mass_signal_pT_%i_eta_%i.pdf", pTbin, etaBin));



      //save plots
      if(pTbin == 0 && etaBin == 0)
      {

        InvMassCan_US_LS[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/K0s_inv_mass_US_LS.pdf(", "pdf");

        InvMassCan[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/K0s_inv_mass.pdf(", "pdf");

      }
      else if( pTbin == nPtBins && etaBin == nEtaBins)
      {

        InvMassCan_US_LS[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/K0s_inv_mass_US_LS.pdf)", "pdf");

        InvMassCan[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/K0s_inv_mass.pdf)", "pdf");

      }
      else
      {

        InvMassCan_US_LS[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/K0s_inv_mass_US_LS.pdf", "pdf");

        InvMassCan[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/K0s_inv_mass.pdf", "pdf");

      }

    }//end for eta

  }//end for pT

  return true;

}


//for analysis of Lambda polarization
void K0sPolarization(TChain *K0s_tree, double invMassRange[2][nPtBins+1][nEtaBins+1], double kappa[nPtBins+1][nEtaBins+1], const int ReadMode, const int energy = 510, const int year = 2017)
{

  TFile *outFile;

  if(ReadMode == 0) //create invariant mass file Form nTuple - run in this mode first
  {
    outFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/ProdPlane_K0s_work.root", "recreate"); //create file to store wrong and good sign M_inv distributions and daughter pT
  }
  else  //read file to store wrong and good sign M_inv distributions and daughter pT
  {

    outFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/ProdPlane_K0s_work.root", "read"); //old, non-optimized cuts

    if( !(outFile->IsOpen()) )
    {
      cout<<"Unable to open file with invariant mass histograms!"<<endl;
      return;
    }
  }



  TFile *EffFile;

  //if(energy == 510) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/K0s_cosThetaStar_eff_Run17_new.root", "read");
  //else if(energy == 200) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/K0s_cosThetaStar_eff_Run12.root", "read");
  if(energy == 510) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/K0s_cosThetaStar_eff_Run17_1B_tight_eta.root", "read");
  else if(energy == 200) EffFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/Eff/K0s_cosThetaStar_eff_Run12_1B_tight_eta.root", "read");
  else
  {
    cout<<"Not a valid colliison energy! Abborting!"<<endl;
    return;
  }
  //_______________________________________________________________________________________________________________________________________________

  //define cuts variables
  //Lambda

  Int_t charge;
  Float_t K0s_mass, K0s_pt, K0s_eta, K0s_decayL, K0s_theta, K0s_DCAdaughters;


  Float_t pi1_pt, pi2_pt;
  Float_t pi1_eta, pi2_eta;
  Float_t pi1_phi, pi2_phi;
  //Float_t pi1_ch;
  Int_t pi1_hasTOFinfo, pi2_hasTOFinfo;

  Float_t thetaProdPlane;

  //---------------SET BARANCH ADDRESSES------------------------
  K0s_tree->SetBranchAddress("pair_charge", &charge);
  K0s_tree->SetBranchAddress("pair_mass", &K0s_mass);
  K0s_tree->SetBranchAddress("pair_pt", &K0s_pt);
  K0s_tree->SetBranchAddress("pair_eta", &K0s_eta);
  K0s_tree->SetBranchAddress("pair_decayL", &K0s_decayL);
  K0s_tree->SetBranchAddress("pair_theta", &K0s_theta);
  K0s_tree->SetBranchAddress("pair_DCAdaughters", &K0s_DCAdaughters);

  K0s_tree->SetBranchAddress("p1_pt", &pi1_pt);
  K0s_tree->SetBranchAddress("p1_eta", &pi1_eta);
  K0s_tree->SetBranchAddress("p1_phi", &pi1_phi);
  //K0s_tree->SetBranchAddress("p1_ch", &pi1_ch);
  K0s_tree->SetBranchAddress("p1_hasTOFinfo", &pi1_hasTOFinfo);

  K0s_tree->SetBranchAddress("p2_pt", &pi2_pt);
  K0s_tree->SetBranchAddress("p2_eta", &pi2_eta);
  K0s_tree->SetBranchAddress("p2_phi", &pi2_phi);
  K0s_tree->SetBranchAddress("p2_hasTOFinfo", &pi2_hasTOFinfo);

  K0s_tree->SetBranchAddress("thetaProdPlane", &thetaProdPlane);

  //--------------------------------------------------------------------------


  TH1D *thetaProdPlane_US_hist[nPtBins+1][nEtaBins+1];
  TH1D *thetaProdPlane_LS_hist[nPtBins+1][nEtaBins+1];


  TH1D *cosThetaProdPlane_US_hist[nPtBins+1][nEtaBins+1];
  TH1D *cosThetaProdPlane_LS_hist[nPtBins+1][nEtaBins+1];

  TH1D *K0s_cosThetaProdPlane_eff_hist[nPtBins+1][nEtaBins+1];


  if(ReadMode == 0) //create histograms to be saved into file
  {

    for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
    {
      for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
      {
        thetaProdPlane_US_hist[pTbin][etaBin] = new TH1D(Form("thetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin), Form("thetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());
        thetaProdPlane_LS_hist[pTbin][etaBin] = new TH1D(Form("thetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin), Form("thetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin), 20, 0, TMath::Pi());

        cosThetaProdPlane_US_hist[pTbin][etaBin] = new TH1D(Form("cosThetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin), Form("cosThetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);
        cosThetaProdPlane_LS_hist[pTbin][etaBin] = new TH1D(Form("cosThetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin), Form("cosThetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin), 20, -1, 1);

        //efficiency histograms
        K0s_cosThetaProdPlane_eff_hist[pTbin][etaBin] = (TH1D*)EffFile->Get(Form("K0s_cosThetaProdPlane_eff_pT_%i_eta_%i", pTbin, etaBin));
      }
    }
  }
  else //load histograms Form file
  {

    for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
    {
      for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
      {
        thetaProdPlane_US_hist[pTbin][etaBin] = (TH1D*)outFile->Get(Form("thetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin));
        thetaProdPlane_LS_hist[pTbin][etaBin] = (TH1D*)outFile->Get(Form("thetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin));

        cosThetaProdPlane_US_hist[pTbin][etaBin] = (TH1D*)outFile->Get(Form("cosThetaProdPlane_US_hist_pT_%i_eta_%i", pTbin, etaBin));
        cosThetaProdPlane_LS_hist[pTbin][etaBin] = (TH1D*)outFile->Get(Form("cosThetaProdPlane_LS_hist_pT_%i_eta_%i", pTbin, etaBin));

        //efficiency histograms
        K0s_cosThetaProdPlane_eff_hist[pTbin][etaBin] = (TH1D*)EffFile->Get(Form("K0s_cosThetaProdPlane_eff_pT_%i_eta_%i", pTbin, etaBin));
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
    if(i%100000 == 0)
    {
      cout<<i<<endl;
    }

    K0s_tree->GetEntry(i);


    double K0s_y = rapidity(K0s_pt, K0s_eta, K0s_mass_PDG);

    //------------------------------------------------------------------------------------------------------------------

    //cuts
    if( !cuts(pi1_hasTOFinfo, pi2_hasTOFinfo, K0s_y, K0s_theta, K0s_decayL) ) continue;

    //------------------------------------------------------------------------------------------------------------------


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

    //find eta bin of Lambda
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
      if(K0s_mass > invMassRange[0][pT_bin][eta_bin] && K0s_mass < invMassRange[1][pT_bin][eta_bin])
      {
        thetaProdPlane_US_hist[pT_bin][eta_bin]->Fill(thetaProdPlane);
        thetaProdPlane_US_hist[nPtBins][eta_bin]->Fill(thetaProdPlane);
        thetaProdPlane_US_hist[pT_bin][nEtaBins]->Fill(thetaProdPlane);
        thetaProdPlane_US_hist[nPtBins][nEtaBins]->Fill(thetaProdPlane);

        cosThetaProdPlane_US_hist[pT_bin][eta_bin]->Fill(cos(thetaProdPlane));
        cosThetaProdPlane_US_hist[nPtBins][eta_bin]->Fill(cos(thetaProdPlane));
        cosThetaProdPlane_US_hist[pT_bin][nEtaBins]->Fill(cos(thetaProdPlane));
        cosThetaProdPlane_US_hist[nPtBins][nEtaBins]->Fill(cos(thetaProdPlane));
      }
    }
    else //like-sign combinations
    {
      if(K0s_mass > invMassRange[0][pT_bin][eta_bin] && K0s_mass < invMassRange[1][pT_bin][eta_bin])
      {
        thetaProdPlane_LS_hist[pT_bin][eta_bin]->Fill(thetaProdPlane);
        thetaProdPlane_LS_hist[nPtBins][eta_bin]->Fill(thetaProdPlane);
        thetaProdPlane_LS_hist[pT_bin][nEtaBins]->Fill(thetaProdPlane);
        thetaProdPlane_LS_hist[nPtBins][nEtaBins]->Fill(thetaProdPlane);

        cosThetaProdPlane_LS_hist[pT_bin][eta_bin]->Fill(cos(thetaProdPlane));
        cosThetaProdPlane_LS_hist[nPtBins][eta_bin]->Fill(cos(thetaProdPlane));
        cosThetaProdPlane_LS_hist[pT_bin][nEtaBins]->Fill(cos(thetaProdPlane));
        cosThetaProdPlane_LS_hist[nPtBins][nEtaBins]->Fill(cos(thetaProdPlane));
      }
    }//end else for flag_int check
  }//end loop over entries in NTuple

  //_________________________________________________________________________________________________________

  if(ReadMode == 0) //if in ReadMode = 0, save all histograms to output file
  {
    outFile->cd();

    //hEventStat1->Write();

    for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
    {
      for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
      {
        thetaProdPlane_US_hist[pTbin][etaBin]->Write();
        thetaProdPlane_LS_hist[pTbin][etaBin]->Write();

        cosThetaProdPlane_US_hist[pTbin][etaBin]->Write();
        cosThetaProdPlane_LS_hist[pTbin][etaBin]->Write();
      }
    }

  } //end if RreadMode = 0

  //InvMassFile->Close();




  TCanvas *thetaProdPlane_US_LS_can[nPtBins+1][nEtaBins+1];

  TCanvas *cosThetaProdPlane_US_LS_can[nPtBins+1][nEtaBins+1];


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);


  for(unsigned int pTbin = 0; pTbin < nPtBins+1; pTbin++)
  {
    for(unsigned int etaBin = 0; etaBin < nEtaBins+1; etaBin++)
    {
      thetaProdPlane_US_LS_can[pTbin][etaBin] = new TCanvas(Form("thetaProdPlane_US_LS_pT_%i_eta_%i_can", pTbin, etaBin), Form("thetaProdPlane_US_LS_pT_%i_eta_%i_can", pTbin, etaBin), 1200, 1000);

      cosThetaProdPlane_US_LS_can[pTbin][etaBin] = new TCanvas(Form("cosThetaProdPlane_US_LS_pT_%i_eta_%i_can", pTbin, etaBin), Form("cosThetaProdPlane_US_LS_pT_%i_eta_%i_can", pTbin, etaBin), 1200, 1000);


      TString *pT_range = new TString();
      if(pTbin < nPtBins) pT_range->Form("%0.1f < #it{p}_{T} < %0.1f GeV/c", pT_bins[pTbin], pT_bins[pTbin+1]);
      else pT_range->Form("p_{T} integrated");

      TString *eta_range = new TString();
      if(etaBin < nEtaBins) eta_range->Form("%0.1f < #eta < %0.1f", eta_bins[etaBin], eta_bins[etaBin+1]);
      else eta_range->Form("-1 < #eta < 1");


      TPaveText *text = new TPaveText(0.35, 0.75, 0.75, 0.89, "NDC");
      //text[pTbin][etaBin]->AddText("STAR preliminary");
      //((TText*)text[pTbin][etaBin]->GetListOfLines()[pTbin][etaBin]->Last())[pTbin][etaBin]->SetTextColor(2);
      text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      text->AddText("Minimum bias");
      text->AddText("K_{s}^{0} and #bar{K_{s}^{0}}");
      text->AddText(eta_range->Data());
      text->AddText(pT_range->Data());
      text->SetFillColorAlpha(0, 0.01);



      thetaProdPlane_US_LS_can[pTbin][etaBin]->cd();

      thetaProdPlane_US_hist[pTbin][etaBin]->SetMarkerStyle(20);
      thetaProdPlane_US_hist[pTbin][etaBin]->SetMarkerColor(kRed);
      thetaProdPlane_US_hist[pTbin][etaBin]->SetLineColor(kRed);
      thetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->SetTitle("#theta*");
      thetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->CenterTitle();
      thetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->SetTitle("d#it{N}/d#theta*");
      thetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      thetaProdPlane_US_hist[pTbin][etaBin]->Sumw2();
      thetaProdPlane_US_hist[pTbin][etaBin]->Scale(1./thetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(1)); //scale by bin width to obtain d N/d theta*
      //thetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 7000);
      thetaProdPlane_US_hist[pTbin][etaBin]->Draw("p e");


      thetaProdPlane_LS_hist[pTbin][etaBin]->SetMarkerStyle(20);
      thetaProdPlane_LS_hist[pTbin][etaBin]->SetMarkerColor(kBlue);
      thetaProdPlane_LS_hist[pTbin][etaBin]->SetLineColor(kBlue);
      thetaProdPlane_LS_hist[pTbin][etaBin]->Sumw2();
      thetaProdPlane_LS_hist[pTbin][etaBin]->Scale(1./thetaProdPlane_LS_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(1)); //scale by bin width to obtain d N/d theta*
      thetaProdPlane_LS_hist[pTbin][etaBin]->Draw("p e same");

      text->Draw("same");

      thetaProdPlane_US_LS_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/thetaProdPlane/thetaProdPlane_US_LS_pT_%i_eta_%i.png", pTbin, etaBin));


      cosThetaProdPlane_US_LS_can[pTbin][etaBin]->cd();

      cosThetaProdPlane_US_hist[pTbin][etaBin]->SetMarkerStyle(20);
      cosThetaProdPlane_US_hist[pTbin][etaBin]->SetMarkerColor(kRed);
      cosThetaProdPlane_US_hist[pTbin][etaBin]->SetLineColor(kRed);
      cosThetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->SetTitle("cos(#theta*)");
      cosThetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->CenterTitle();
      cosThetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
      cosThetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->CenterTitle();
      cosThetaProdPlane_US_hist[pTbin][etaBin]->Sumw2();
      cosThetaProdPlane_US_hist[pTbin][etaBin]->Scale(1./cosThetaProdPlane_US_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(1));
      cosThetaProdPlane_US_hist[pTbin][etaBin]->Divide(K0s_cosThetaProdPlane_eff_hist[pTbin][etaBin]);
      //cosThetaProdPlane_US_hist[pTbin][etaBin]->GetYaxis()->SetRangeUser(0, 7000);
      cosThetaProdPlane_US_hist[pTbin][etaBin]->SetMinimum(0);
      cosThetaProdPlane_US_hist[pTbin][etaBin]->Draw("p e");


      cosThetaProdPlane_LS_hist[pTbin][etaBin]->SetMarkerStyle(20);
      cosThetaProdPlane_LS_hist[pTbin][etaBin]->SetMarkerColor(kBlue);
      cosThetaProdPlane_LS_hist[pTbin][etaBin]->SetLineColor(kBlue);
      cosThetaProdPlane_LS_hist[pTbin][etaBin]->Sumw2();
      cosThetaProdPlane_LS_hist[pTbin][etaBin]->Scale(1./cosThetaProdPlane_LS_hist[pTbin][etaBin]->GetXaxis()->GetBinWidth(1));
      cosThetaProdPlane_LS_hist[pTbin][etaBin]->Draw("p e same");

      text->Draw("same");

      cosThetaProdPlane_US_LS_can[pTbin][etaBin]->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/thetaProdPlane/cosThetaProdPlane_US_LS_pT_%i_eta_%i.png", pTbin, etaBin));




      if(pTbin == 0 && etaBin == 0)
      {
        cosThetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/thetaProdPlane/cosThetaProdPlane_US_LS.pdf(", "pdf");
        thetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/thetaProdPlane/thetaProdPlane_US_LS.pdf(", "pdf");
      }
      else if( pTbin == nPtBins && etaBin == nEtaBins)
      {
        cosThetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/thetaProdPlane/cosThetaProdPlane_US_LS.pdf)", "pdf");
        thetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/thetaProdPlane/thetaProdPlane_US_LS.pdf)", "pdf");
      }
      else
      {
        cosThetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/thetaProdPlane/cosThetaProdPlane_US_LS.pdf", "pdf");
        thetaProdPlane_US_LS_can[pTbin][etaBin]->Print("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/K0s/thetaProdPlane/thetaProdPlane_US_LS.pdf", "pdf");
      }


    }
  }

  //------------------------------------------------------------------------------------------------------------


}


//ReadMode = 0 - read TTree, ReadMode = 1 - read histograms - First run in ReadMode = 0 to save relevant histograms, then can run in ReadMode = 1 to read just histograms and save time
//energy - collision energy in GeV
void Ana002_K0s_reference(const int ReadMode = 0, const int energy = 510, const int year = 2017)
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


  //Lambda signal to signal+background ratio
  //for bacground subtraction for polarization measurement
  double kappa[nPtBins+1][nEtaBins+1];



  bool invMassFinish = InvMass(myChain, invMassRange, kappa, ReadMode, energy , year);

  if(!invMassFinish)
  {
    cout<<"Analysis of invariant spectra ended abnormally. Abborting!"<<endl;

    return;
  }

  K0sPolarization(myChain, invMassRange, kappa, ReadMode, energy , year);

  cout<<endl;
  cout<<"Nubmer of accepted events: "<<hEventStat1->GetBinContent(6)<<endl;


  return;
}
