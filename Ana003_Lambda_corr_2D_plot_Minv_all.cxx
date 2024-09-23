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
#include"TLine.h"


using namespace std;


const int nPtBins_corr = 2;
float const pT_bins_corr[nPtBins_corr+1] = { 0.5, 1.5, 5.};

//const int nEtaBins = 3;
//float const eta_bins[nEtaBins+1] = { -1, -0.2, 0.2, 1 };
//float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };

void Ana003_Lambda_corr_2D_plot_Minv_all(const int cut_type = 0, const int energy = 510, const int year = 2017)
{


  TFile *InvMassFile; //output file to store invariant mass histograms

  if(cut_type == 0)
  {
    if(year == 2012)
    {
      InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts_work.root", year), "read");

      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/2024_06_new_ME/InvMass_Lambda_ana_cuts_new_ME.root", year), "read");
      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/2024_06_new_ME/InvMass_Lambda_ana_cuts_new_ME.root", year), "read");

      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts_012024.root", year), "read");

      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_orig_full_prod.root", year), "read");
      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_new_prod_open_Vz.root", year), "read");

    }
    else if(year == 2015)
    {
      InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts_work.root", year), "read");
      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts_new_no_TOF.root", year), "read");
    }
    else
    {
      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts.root", year), "read");
      InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts_work.root", year), "read");
    }
  }

  if(cut_type == 1) InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_tight_topo_cuts_work.root", year), "read");
  if(cut_type == 2) InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_tight_pT_cut_work.root", year), "read");

  if( !(InvMassFile->IsOpen()) )
  {
    cout<<"Unable to open file with invariant mass histograms!"<<endl;
    return;
  }


  //_______________________________________________________________________________________________________________________________________________


  //need to change binning
  // 2D - pT1 vs. pT2
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_all_US = (TH2F*)InvMassFile->Get("L0_inv_mass_vs_L0bar_inv_mass_all_US"); //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_all_US_LS = (TH2F*)InvMassFile->Get("L0_inv_mass_vs_L0bar_inv_mass_all_US_LS"); //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_all_LS = (TH2F*)InvMassFile->Get("L0_inv_mass_vs_L0bar_inv_mass_all_LS"); //for LS-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass;



  TH2F *L0_inv_mass_vs_L0_inv_mass_all_US = (TH2F*)InvMassFile->Get("L0_inv_mass_vs_L0_inv_mass_all_US"); //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_all_US_LS = (TH2F*)InvMassFile->Get("L0_inv_mass_vs_L0_inv_mass_all_US_LS"); //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_all_LS = (TH2F*)InvMassFile->Get("L0_inv_mass_vs_L0_inv_mass_all_LS"); //for LS-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass;


  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_all_US = (TH2F*)InvMassFile->Get("L0bar_inv_mass_vs_L0bar_inv_mass_all_US"); //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS = (TH2F*)InvMassFile->Get("L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS"); //for US-LS Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_all_LS = (TH2F*)InvMassFile->Get("L0bar_inv_mass_vs_L0bar_inv_mass_all_LS"); //for LS-LS Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass;


  //________________________________________________________________________________________________________

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(43);
  gStyle->SetLegendTextSize(30);



  float sideBandScale_L_Lbar;
  float sideBandScale_L_L;
  float sideBandScale_Lbar_Lbar;


  //find bins for side band and for projection to x
  //projection bins just for testing, later can do projections with parameters of 2D fit
  int binLow = L0_inv_mass_vs_L0bar_inv_mass_all_US->GetXaxis()->FindBin(1.11);
  int binHigh = L0_inv_mass_vs_L0bar_inv_mass_all_US->GetXaxis()->FindBin(1.12);

  int binLow_sideBand = L0_inv_mass_vs_L0bar_inv_mass_all_US->GetXaxis()->FindBin(1.14);
  int binHigh_sideBand = L0_inv_mass_vs_L0bar_inv_mass_all_US->GetXaxis()->FindBin(1.18);

  //L-Lbar
  //US and LS
  TCanvas *L0_inv_mass_vs_L0bar_inv_mass_all_US_can = new TCanvas("L0_inv_mass_vs_L0bar_inv_mass_all_US_can", "L0_inv_mass_vs_L0bar_inv_mass_all_US_can", 1200, 1000);
  L0_inv_mass_vs_L0bar_inv_mass_all_US_can->cd();

  L0_inv_mass_vs_L0bar_inv_mass_all_US->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
  L0_inv_mass_vs_L0bar_inv_mass_all_US->GetXaxis()->SetTitleOffset(2);
  L0_inv_mass_vs_L0bar_inv_mass_all_US->GetXaxis()->CenterTitle();
  L0_inv_mass_vs_L0bar_inv_mass_all_US->GetXaxis()->SetRangeUser(1.07, 1.2);
  L0_inv_mass_vs_L0bar_inv_mass_all_US->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
  L0_inv_mass_vs_L0bar_inv_mass_all_US->GetYaxis()->SetTitleOffset(2.1);
  L0_inv_mass_vs_L0bar_inv_mass_all_US->GetYaxis()->CenterTitle();
  L0_inv_mass_vs_L0bar_inv_mass_all_US->GetYaxis()->SetRangeUser(1.07,1.2);
  float L_Lbar_US_ratio_integral = L0_inv_mass_vs_L0bar_inv_mass_all_US->Integral(binLow, binHigh, binLow, binHigh);
  float nLLbar_US_tot = L0_inv_mass_vs_L0bar_inv_mass_all_US->Integral();
  L0_inv_mass_vs_L0bar_inv_mass_all_US->Draw("surf1");

  TPaveText *L0_L0bar_text_Minv_2D_US = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
  L0_L0bar_text_Minv_2D_US->SetTextFont(43);
  L0_L0bar_text_Minv_2D_US->SetTextSize(30);
  L0_L0bar_text_Minv_2D_US->AddText("STAR");
  //L0_L0bar_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  //L0_L0bar_text_Minv_2D_US->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0_L0bar_text_Minv_2D_US->AddText(Form("p+p #sqrt{s} = %i GeV", energy));
  L0_L0bar_text_Minv_2D_US->AddText("Minimum bias");
  L0_L0bar_text_Minv_2D_US->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_Minv_2D_US->Draw("same");


  TPaveText *L0_L0bar_text_Minv_2D_kine = new TPaveText(0.73, 0.83, 0.98, 0.98, "NDC");
  L0_L0bar_text_Minv_2D_kine->SetTextFont(43);
  L0_L0bar_text_Minv_2D_kine->SetTextSize(30);
  L0_L0bar_text_Minv_2D_kine->AddText("|#it{y}| < 1");
  L0_L0bar_text_Minv_2D_kine->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_Minv_2D_kine->Draw("same");

  TPaveText *L0_L0bar_text_STAR = new TPaveText(0.4, 0.9, 0.6, 0.98, "NDC");
  L0_L0bar_text_STAR->SetTextFont(43);
  L0_L0bar_text_STAR->SetTextSize(30);
  L0_L0bar_text_STAR->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
  L0_L0bar_text_STAR->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
  L0_L0bar_text_STAR->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_STAR->Draw("same");

  L0_inv_mass_vs_L0bar_inv_mass_all_US_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0bar_inv_mass_all_US.png");
  L0_inv_mass_vs_L0bar_inv_mass_all_US_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0bar_inv_mass_all_US.pdf");



  TCanvas *L0_inv_mass_vs_L0bar_inv_mass_all_LS_can = new TCanvas("L0_inv_mass_vs_L0bar_inv_mass_all_LS_can", "L0_inv_mass_vs_L0bar_inv_mass_all_LS_can", 1200, 1000);
  L0_inv_mass_vs_L0bar_inv_mass_all_LS_can->cd();

  L0_inv_mass_vs_L0bar_inv_mass_all_LS->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
  L0_inv_mass_vs_L0bar_inv_mass_all_LS->GetXaxis()->SetTitleOffset(2);
  L0_inv_mass_vs_L0bar_inv_mass_all_LS->GetXaxis()->CenterTitle();
  L0_inv_mass_vs_L0bar_inv_mass_all_LS->GetXaxis()->SetRangeUser(1.07, 1.2);
  L0_inv_mass_vs_L0bar_inv_mass_all_LS->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
  L0_inv_mass_vs_L0bar_inv_mass_all_LS->GetYaxis()->SetTitleOffset(2.1);
  L0_inv_mass_vs_L0bar_inv_mass_all_LS->GetYaxis()->CenterTitle();
  L0_inv_mass_vs_L0bar_inv_mass_all_LS->GetYaxis()->SetRangeUser(1.07,1.2);
  L0_inv_mass_vs_L0bar_inv_mass_all_LS->Draw("surf1");

  TPaveText *L0_L0bar_text_Minv_2D_LS = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
  L0_L0bar_text_Minv_2D_LS->SetTextFont(42);
  //L0_L0bar_text_no_corr->AddText("STAR Internal");
  //L0_L0bar_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0bar_text_Minv_2D_LS->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0_L0bar_text_Minv_2D_LS->AddText("Minimum bias");
  L0_L0bar_text_Minv_2D_LS->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
  L0_L0bar_text_Minv_2D_LS->AddText("LS(#pi^{+}p) vs. LS(#pi^{-}#bar{p})");
  L0_L0bar_text_Minv_2D_LS->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_Minv_2D_LS->Draw("same");

  L0_L0bar_text_Minv_2D_kine->Draw("same");

  L0_inv_mass_vs_L0bar_inv_mass_all_LS_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0bar_inv_mass_all_LS.png");
  L0_inv_mass_vs_L0bar_inv_mass_all_LS_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0bar_inv_mass_all_LS.pdf");



  TCanvas *L0_inv_mass_vs_L0bar_inv_mass_all_US_LS_can = new TCanvas("L0_inv_mass_vs_L0bar_inv_mass_all_US_LS_can", "L0_inv_mass_vs_L0bar_inv_mass_all_US_LS_can", 1200, 1000);
  L0_inv_mass_vs_L0bar_inv_mass_all_US_LS_can->cd();

  //background already scaled from Ana003_Lambda_corr_2D.cxx
  //here just potting
  L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->Sumw2();
  L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
  L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetXaxis()->SetTitleOffset(2);
  L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetXaxis()->CenterTitle();
  L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetXaxis()->SetRangeUser(1.07, 1.2);
  L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
  L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetYaxis()->SetTitleOffset(2);
  L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetYaxis()->CenterTitle();
  L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetYaxis()->SetRangeUser(1.07,1.2);
  float L_Lbar_US_LS_ratio_integral = L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->Integral(binLow, binHigh, binLow, binHigh);
  float nLLbar_US_LS_tot = L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->Integral(1,180,1,180);
  L0_inv_mass_vs_L0bar_inv_mass_all_US_LS->Draw("surf1");

  TPaveText *L0_L0bar_text_Minv_2D_US_LS = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
  L0_L0bar_text_Minv_2D_US_LS->SetTextFont(42);
  //L0_L0bar_text_no_corr->AddText("STAR Internal");
  //L0_L0bar_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0bar_text_Minv_2D_US_LS->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0_L0bar_text_Minv_2D_US_LS->AddText("Minimum bias");
  L0_L0bar_text_Minv_2D_US_LS->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
  L0_L0bar_text_Minv_2D_US_LS->AddText("US(#pi^{-}p) vs. LS(#pi^{-}#bar{p})");
  L0_L0bar_text_Minv_2D_US_LS->AddText("LS(#pi^{+}p) vs. US(#pi^{+}#bar{p})");
  L0_L0bar_text_Minv_2D_US_LS->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_Minv_2D_US_LS->Draw("same");

  L0_L0bar_text_Minv_2D_kine->Draw("same");

  L0_inv_mass_vs_L0bar_inv_mass_all_US_LS_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0bar_inv_mass_all_US_LS.png");
  L0_inv_mass_vs_L0bar_inv_mass_all_US_LS_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0bar_inv_mass_all_US_LS.pdf");



  //------------------------------------------------------------------------------------------------------------------------------------

  //US-LS
  TF2 *doubleGauss_L_Lbar = new TF2("doubleGauss_L_Lbar", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 1.07, 1.2, 1.07, 1.2);
  doubleGauss_L_Lbar->SetParameters(1200, 1.116, 0.002, 1.116, 0.002);


  TCanvas *L0_inv_mass_vs_L0bar_inv_mass_all_can = new TCanvas("L0_inv_mass_vs_L0bar_inv_mass_all_can", "L0_inv_mass_vs_L0bar_inv_mass_all_can", 1200, 1000);
  L0_inv_mass_vs_L0bar_inv_mass_all_can->cd();

  L0_inv_mass_vs_L0bar_inv_mass = (TH2F*)L0_inv_mass_vs_L0bar_inv_mass_all_US->Clone("L0_inv_mass_vs_L0bar_inv_mass");
  L0_inv_mass_vs_L0bar_inv_mass->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
  L0_inv_mass_vs_L0bar_inv_mass->GetXaxis()->SetTitleOffset(2);
  L0_inv_mass_vs_L0bar_inv_mass->GetXaxis()->CenterTitle();
  L0_inv_mass_vs_L0bar_inv_mass->GetXaxis()->SetRangeUser(1.07, 1.2);
  L0_inv_mass_vs_L0bar_inv_mass->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
  L0_inv_mass_vs_L0bar_inv_mass->GetYaxis()->SetTitleOffset(2.1);
  L0_inv_mass_vs_L0bar_inv_mass->GetYaxis()->CenterTitle();
  L0_inv_mass_vs_L0bar_inv_mass->GetYaxis()->SetRangeUser(1.07,1.2);
  //L0_inv_mass_vs_L0bar_inv_mass->Add(L0_inv_mass_vs_L0bar_inv_mass_all_LS, -1);
  L0_inv_mass_vs_L0bar_inv_mass->Add(L0_inv_mass_vs_L0bar_inv_mass_all_US_LS, -1);
  L0_inv_mass_vs_L0bar_inv_mass->Fit(doubleGauss_L_Lbar, "S 0", "", 1.11, 1.125);
  L0_inv_mass_vs_L0bar_inv_mass->Draw("surf1");


  TPaveText *L0_L0bar_text_Minv_2D_sig = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
  L0_L0bar_text_Minv_2D_sig->SetTextFont(42);
  //L0_L0bar_text_no_corr->AddText("STAR Internal");
  //L0_L0bar_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0bar_text_Minv_2D_sig->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0_L0bar_text_Minv_2D_sig->AddText("Minimum bias");
  L0_L0bar_text_Minv_2D_sig->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
  L0_L0bar_text_Minv_2D_sig->AddText("US - Bckg.");
  L0_L0bar_text_Minv_2D_sig->SetFillColorAlpha(0, 0.01);
  L0_L0bar_text_Minv_2D_sig->Draw("same");

  L0_L0bar_text_Minv_2D_kine->Draw("same");

  //doubleGauss_L_Lbar->Draw("surf1");

  L0_inv_mass_vs_L0bar_inv_mass_all_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0bar_inv_mass.png");
  L0_inv_mass_vs_L0bar_inv_mass_all_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0bar_inv_mass.pdf");

  //________________________________________________________________________________________________________

  //L-L
  //US
  TCanvas *L0_inv_mass_vs_L0_inv_mass_all_US_can = new TCanvas("L0_inv_mass_vs_L0_inv_mass_all_US_can", "L0_inv_mass_vs_L0_inv_mass_all_US_can", 1200, 1000);
  L0_inv_mass_vs_L0_inv_mass_all_US_can->cd();

  L0_inv_mass_vs_L0_inv_mass_all_US->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
  L0_inv_mass_vs_L0_inv_mass_all_US->GetXaxis()->SetTitleOffset(2);
  L0_inv_mass_vs_L0_inv_mass_all_US->GetXaxis()->CenterTitle();
  L0_inv_mass_vs_L0_inv_mass_all_US->GetXaxis()->SetRangeUser(1.07, 1.2);
  L0_inv_mass_vs_L0_inv_mass_all_US->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
  L0_inv_mass_vs_L0_inv_mass_all_US->GetYaxis()->SetTitleOffset(2.1);
  L0_inv_mass_vs_L0_inv_mass_all_US->GetYaxis()->CenterTitle();
  L0_inv_mass_vs_L0_inv_mass_all_US->GetYaxis()->SetRangeUser(1.07,1.2);
  float L_L_US_ratio_integral = L0_inv_mass_vs_L0_inv_mass_all_US->Integral(binLow, binHigh, binLow, binHigh);
  float nLL_US_tot = L0_inv_mass_vs_L0_inv_mass_all_US->Integral();
  L0_inv_mass_vs_L0_inv_mass_all_US->Draw("surf1");


  TPaveText *L0_L0_text_Minv_2D_US = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
  L0_L0_text_Minv_2D_US->SetTextFont(42);
  //L0_L0_text_no_corr->AddText("STAR Internal");
  //L0_L0_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0_text_Minv_2D_US->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0_L0_text_Minv_2D_US->AddText("Minimum bias");
  L0_L0_text_Minv_2D_US->AddText("#Lambda^{0}-#Lambda^{0}");
  L0_L0_text_Minv_2D_US->AddText("US(#pi^{-}p) vs. US(#pi^{-}p)");
  L0_L0_text_Minv_2D_US->SetFillColorAlpha(0, 0.01);
  L0_L0_text_Minv_2D_US->Draw("same");

  L0_L0bar_text_Minv_2D_kine->Draw("same");

  L0_inv_mass_vs_L0_inv_mass_all_US_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0_inv_mass_all_US.png");
  L0_inv_mass_vs_L0_inv_mass_all_US_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0_inv_mass_all_US.pdf");



  TCanvas *L0_inv_mass_vs_L0_inv_mass_all_US_LS_can = new TCanvas("L0_inv_mass_vs_L0_inv_mass_all_US_LS_can", "L0_inv_mass_vs_L0_inv_mass_all_US_LS_can", 1200, 1000);
  L0_inv_mass_vs_L0_inv_mass_all_US_LS_can->cd();

  //background already scaled from Ana003_Lambda_corr_2D.cxx
  //here just potting
  L0_inv_mass_vs_L0_inv_mass_all_US_LS->Sumw2();
  L0_inv_mass_vs_L0_inv_mass_all_US_LS->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
  L0_inv_mass_vs_L0_inv_mass_all_US_LS->GetXaxis()->SetTitleOffset(2);
  L0_inv_mass_vs_L0_inv_mass_all_US_LS->GetXaxis()->CenterTitle();
  L0_inv_mass_vs_L0_inv_mass_all_US_LS->GetXaxis()->SetRangeUser(1.07, 1.2);
  L0_inv_mass_vs_L0_inv_mass_all_US_LS->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
  L0_inv_mass_vs_L0_inv_mass_all_US_LS->GetYaxis()->SetTitleOffset(2.1);
  L0_inv_mass_vs_L0_inv_mass_all_US_LS->GetYaxis()->CenterTitle();
  L0_inv_mass_vs_L0_inv_mass_all_US_LS->GetYaxis()->SetRangeUser(1.07,1.2);
  float L_L_US_LS_ratio_integral = L0_inv_mass_vs_L0_inv_mass_all_US_LS->Integral(binLow, binHigh, binLow, binHigh);
  float nLL_US_LS_tot = L0_inv_mass_vs_L0_inv_mass_all_US_LS->Integral();
  L0_inv_mass_vs_L0_inv_mass_all_US_LS->Draw("surf1");

  TPaveText *L0_L0_text_Minv_2D_US_LS = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
  L0_L0_text_Minv_2D_US_LS->SetTextFont(42);
  //L0_L0_text_no_corr->AddText("STAR Internal");
  //L0_L0_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0_text_Minv_2D_US_LS->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0_L0_text_Minv_2D_US_LS->AddText("Minimum bias");
  L0_L0_text_Minv_2D_US_LS->AddText("#Lambda^{0}-#Lambda^{0}");
  L0_L0_text_Minv_2D_US_LS->AddText("US(#pi^{-}p) vs. LS(#pi^{+}p)");
  L0_L0_text_Minv_2D_US_LS->AddText("LS(#pi^{+}p) vs. US(#pi^{-}p)");
  L0_L0_text_Minv_2D_US_LS->SetFillColorAlpha(0, 0.01);
  L0_L0_text_Minv_2D_US_LS->Draw("same");

  L0_L0bar_text_Minv_2D_kine->Draw("same");

  L0_inv_mass_vs_L0_inv_mass_all_US_LS_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0_inv_mass_all_US_LS.png");
  L0_inv_mass_vs_L0_inv_mass_all_US_LS_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0_inv_mass_all_US_LS.pdf");



  //US-LS
  TF2 *doubleGauss_L_L = new TF2("doubleGauss_L_L", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 1.07, 1.2, 1.07, 1.2);
  //doubleGauss_L_L->SetParameters(1200, 1.116, 0.002, 1.116, 0.002);
  if( energy == 200 ) doubleGauss_L_L->SetParameters(1000, 1.116, 0.002, 1.116, 0.002);


  TCanvas *L0_inv_mass_vs_L0_inv_mass_all_can = new TCanvas("L0_inv_mass_vs_L0_inv_mass_all_can", "L0_inv_mass_vs_L0_inv_mass_all_can", 1200, 1000);
  L0_inv_mass_vs_L0_inv_mass_all_can->cd();

  L0_inv_mass_vs_L0_inv_mass = (TH2F*)L0_inv_mass_vs_L0_inv_mass_all_US->Clone("L0_inv_mass_vs_L0_inv_mass");
  L0_inv_mass_vs_L0_inv_mass->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
  L0_inv_mass_vs_L0_inv_mass->GetXaxis()->SetTitleOffset(2);
  L0_inv_mass_vs_L0_inv_mass->GetXaxis()->CenterTitle();
  L0_inv_mass_vs_L0_inv_mass->GetXaxis()->SetRangeUser(1.07,1.2);
  L0_inv_mass_vs_L0_inv_mass->GetYaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
  L0_inv_mass_vs_L0_inv_mass->GetYaxis()->SetTitleOffset(2.1);
  L0_inv_mass_vs_L0_inv_mass->GetYaxis()->CenterTitle();
  L0_inv_mass_vs_L0_inv_mass->GetYaxis()->SetRangeUser(1.07,1.2);
  //L0_inv_mass_vs_L0_inv_mass->Add(L0_inv_mass_vs_L0_inv_mass_all_LS, -1);
  L0_inv_mass_vs_L0_inv_mass->Add(L0_inv_mass_vs_L0_inv_mass_all_US_LS, -1);
  L0_inv_mass_vs_L0_inv_mass->Fit(doubleGauss_L_L, "s 0", "", 1.11, 1.125);
  L0_inv_mass_vs_L0_inv_mass->Draw("surf1");


  TPaveText *L0_L0_text_Minv_2D_sig = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
  L0_L0_text_Minv_2D_sig->SetTextFont(42);
  //L0_L0_text_no_corr->AddText("STAR Internal");
  //L0_L0_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0_L0_text_Minv_2D_sig->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0_L0_text_Minv_2D_sig->AddText("Minimum bias");
  L0_L0_text_Minv_2D_sig->AddText("#Lambda^{0}-#Lambda^{0}");
  L0_L0_text_Minv_2D_sig->AddText("US - Bckg.");
  L0_L0_text_Minv_2D_sig->SetFillColorAlpha(0, 0.01);
  L0_L0_text_Minv_2D_sig->Draw("same");

  L0_L0bar_text_Minv_2D_kine->Draw("same");

  L0_inv_mass_vs_L0_inv_mass_all_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0_inv_mass.png");
  L0_inv_mass_vs_L0_inv_mass_all_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0_inv_mass_vs_L0_inv_mass.pdf");

  //________________________________________________________________________________________________________

  //Lbar-Lbar

  //US
  TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_all_US_can = new TCanvas("L0bar_inv_mass_vs_L0bar_inv_mass_all_US_can", "L0bar_inv_mass_vs_L0bar_inv_mass_all_US_can", 1200, 1000);
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_can->cd();

  L0bar_inv_mass_vs_L0bar_inv_mass_all_US->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US->GetXaxis()->SetTitleOffset(2);
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US->GetXaxis()->CenterTitle();
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US->GetXaxis()->SetRangeUser(1.07, 1.2);
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US->GetYaxis()->SetTitleOffset(2.1);
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US->GetYaxis()->CenterTitle();
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US->GetYaxis()->SetRangeUser(1.07,1.2);
  float Lbar_Lbar_US_ratio_integral = L0bar_inv_mass_vs_L0bar_inv_mass_all_US->Integral(binLow, binHigh, binLow, binHigh);
  float nLbarLbar_US_tot = L0bar_inv_mass_vs_L0bar_inv_mass_all_US->Integral();
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US->Draw("surf1");

  TPaveText *L0bar_L0bar_text_Minv_2D_US = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
  L0bar_L0bar_text_Minv_2D_US->SetTextFont(42);
  //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
  //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0bar_L0bar_text_Minv_2D_US->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0bar_L0bar_text_Minv_2D_US->AddText("Minimum bias");
  L0bar_L0bar_text_Minv_2D_US->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
  L0bar_L0bar_text_Minv_2D_US->AddText("US(#pi^{+}#bar{p}) vs. US(#pi^{+}#bar{p})");
  L0bar_L0bar_text_Minv_2D_US->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_Minv_2D_US->Draw("same");

  L0_L0bar_text_Minv_2D_kine->Draw("same");

  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0bar_inv_mass_vs_L0bar_inv_mass_all_US.png");
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0bar_inv_mass_vs_L0bar_inv_mass_all_US.pdf");



  TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS_can = new TCanvas("L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS_can", "L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS_can", 1200, 1000);
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS_can->cd();

  //background already scaled from Ana003_Lambda_corr_2D.cxx
  //here just potting

  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->Sumw2();
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetXaxis()->SetTitleOffset(2);
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetXaxis()->CenterTitle();
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetXaxis()->SetRangeUser(1.07, 1.2);
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetYaxis()->SetTitleOffset(2.1);
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetYaxis()->CenterTitle();
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->GetYaxis()->SetRangeUser(1.07,1.2);
  float Lbar_Lbar_US_LS_ratio_integral = L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->Integral(binLow, binHigh, binLow, binHigh);
  float nLbarLbar_US_LS_tot = L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->Integral();
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS->Draw("surf1");

  TPaveText *L0bar_L0bar_text_Minv_2D_US_LS = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
  L0bar_L0bar_text_Minv_2D_US_LS->SetTextFont(42);
  //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
  //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0bar_L0bar_text_Minv_2D_US_LS->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0bar_L0bar_text_Minv_2D_US_LS->AddText("Minimum bias");
  L0bar_L0bar_text_Minv_2D_US_LS->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
  L0bar_L0bar_text_Minv_2D_US_LS->AddText("US(#pi^{+}#bar{p}) vs. LS(#pi^{+}p)");
  L0bar_L0bar_text_Minv_2D_US_LS->AddText("LS(#pi^{+}p) vs. US(#pi^{+}#bar{p})");
  L0bar_L0bar_text_Minv_2D_US_LS->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_Minv_2D_US_LS->Draw("same");

  L0_L0bar_text_Minv_2D_kine->Draw("same");

  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS.png");
  L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS.pdf");



  //US-LS
  TF2 *doubleGauss_Lbar_Lbar = new TF2("doubleGauss_Lbar_Lbar", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 1.07, 1.2, 1.07, 1.2);
  //doubleGauss_Lbar_Lbar->SetParameters(1200, 1.116, 0.002, 1.116, 0.002);
  doubleGauss_Lbar_Lbar->SetParameters(1000, 1.116, 0.002, 1.116, 0.002);


  TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_all_can = new TCanvas("L0bar_inv_mass_vs_L0bar_inv_mass_all_can", "L0bar_inv_mass_vs_L0bar_inv_mass_all_can", 1200, 1000);
  L0bar_inv_mass_vs_L0bar_inv_mass_all_can->cd();

  L0bar_inv_mass_vs_L0bar_inv_mass = (TH2F*)L0bar_inv_mass_vs_L0bar_inv_mass_all_US->Clone("L0bar_inv_mass_vs_L0bar_inv_mass");
  L0bar_inv_mass_vs_L0bar_inv_mass->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
  L0bar_inv_mass_vs_L0bar_inv_mass->GetXaxis()->SetTitleOffset(2);
  L0bar_inv_mass_vs_L0bar_inv_mass->GetXaxis()->CenterTitle();
  L0bar_inv_mass_vs_L0bar_inv_mass->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
  L0bar_inv_mass_vs_L0bar_inv_mass->GetYaxis()->SetTitleOffset(2.1);
  L0bar_inv_mass_vs_L0bar_inv_mass->GetYaxis()->CenterTitle();
  //L0bar_inv_mass_vs_L0bar_inv_mass->Add(L0bar_inv_mass_vs_L0bar_inv_mass_all_LS, -1);
  L0bar_inv_mass_vs_L0bar_inv_mass->Add(L0bar_inv_mass_vs_L0bar_inv_mass_all_US_LS, -1);
  //if( pTbin1 == 0 && pTbin2 == 1 ) fit_res_gaus_L0Bar_L0Bar = L0bar_inv_mass_vs_L0bar_inv_mass->Fit(doubleGauss_Lbar_Lbar, "s 0", "", 1.1, 1.12);
  //else if( pTbin1 == 1 && pTbin2 == 1 ) fit_res_gaus_L0Bar_L0Bar = L0bar_inv_mass_vs_L0bar_inv_mass->Fit(doubleGauss_Lbar_Lbar, "s 0", "", 1.1, 1.12);
  //else fit_res_gaus_L0Bar_L0Bar = L0bar_inv_mass_vs_L0bar_inv_mass->Fit(doubleGauss_Lbar_Lbar, "s 0", "", 1.11, 1.125);
  L0bar_inv_mass_vs_L0bar_inv_mass->Fit(doubleGauss_Lbar_Lbar, "s 0", "", 1.113, 1.125);
  L0bar_inv_mass_vs_L0bar_inv_mass->Draw("surf1");


  TPaveText *L0bar_L0bar_text_Minv_2D_sig = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
  L0bar_L0bar_text_Minv_2D_sig->SetTextFont(42);
  //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
  //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
  //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
  L0bar_L0bar_text_Minv_2D_sig->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  L0bar_L0bar_text_Minv_2D_sig->AddText("Minimum bias");
  L0bar_L0bar_text_Minv_2D_sig->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
  L0bar_L0bar_text_Minv_2D_sig->AddText("US - Bckg.");
  L0bar_L0bar_text_Minv_2D_sig->SetFillColorAlpha(0, 0.01);
  L0bar_L0bar_text_Minv_2D_sig->Draw("same");

  L0_L0bar_text_Minv_2D_kine->Draw("same");

  L0bar_inv_mass_vs_L0bar_inv_mass_all_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0bar_inv_mass_vs_L0bar_inv_mass.png");
  L0bar_inv_mass_vs_L0bar_inv_mass_all_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_all/L0bar_inv_mass_vs_L0bar_inv_mass.pdf");

  //________________________________________________________________________________________________________

  cout<<endl;
  cout<<"All"<<endl;
  cout<<"N L-Lbar pairs: "<<nLLbar_US_tot<<endl;
  cout<<"N L-Lbar background pairs: "<<nLLbar_US_LS_tot<<endl;
  cout<<endl;
  cout<<"N L-L pairs: "<<nLL_US_tot<<endl;
  cout<<"N L-L background pairs: "<<nLL_US_LS_tot<<endl;
  cout<<endl;
  cout<<"N Lbar-Lbar pairs: "<<nLbarLbar_US_tot<<endl;
  cout<<"N Lbar-Lbar background pairs: "<<nLbarLbar_US_LS_tot<<endl;

  cout<<endl;
  cout<<"+-3sigma"<<endl;
  cout<<"N L-Lbar pairs: "<<L_Lbar_US_ratio_integral<<endl;
  cout<<"N L-Lbar background pairs: "<<L_Lbar_US_LS_ratio_integral<<endl;
  cout<<endl;
  cout<<"N L-L pairs: "<<L_L_US_ratio_integral<<endl;
  cout<<"N L-L background pairs: "<<L_L_US_LS_ratio_integral<<endl;
  cout<<endl;
  cout<<"N Lbar-Lbar pairs: "<<Lbar_Lbar_US_ratio_integral<<endl;
  cout<<"N Lbar-Lbar background pairs: "<<Lbar_Lbar_US_LS_ratio_integral<<endl;

  InvMassFile->Close();

  return;
}
