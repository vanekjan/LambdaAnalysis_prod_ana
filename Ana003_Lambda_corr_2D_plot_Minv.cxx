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

void Ana003_Lambda_corr_2D_plot_Minv(const int cut_type = 0, const int energy = 510, const int year = 2017)
{


  TFile *InvMassFile; //output file to store invariant mass histograms

  if(cut_type == 0)
  {
    if(year == 2012)
    {
      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts_work.root", year), "read");

      InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/2024_06_new_ME/InvMass_Lambda_ana_cuts_new_ME.root", year), "read");
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
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_LS[nPtBins_corr][nPtBins_corr]; //for LS-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass[nPtBins_corr][nPtBins_corr];

  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_zoom[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_zoom[nPtBins_corr][nPtBins_corr];


  TH2F *L0_inv_mass_vs_L0_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_LS[nPtBins_corr][nPtBins_corr]; //for LS-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass[nPtBins_corr][nPtBins_corr];

  TH2F *L0_inv_mass_vs_L0_inv_mass_US_zoom[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_US_LS_zoom[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_zoom[nPtBins_corr][nPtBins_corr];


  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_LS[nPtBins_corr][nPtBins_corr]; //for LS-LS Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass[nPtBins_corr][nPtBins_corr];

  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_zoom[nPtBins_corr][nPtBins_corr];


  //----------------------------------------------------------------------------


  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2));

      L0_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_US_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2));


      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2));

      L0_inv_mass_vs_L0_inv_mass_US_zoom[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_US_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass_US_LS_zoom[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2));


      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_LS_pT1_%i_pT2_%i", pTbin1, pTbin2));

      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_zoom_pT1_%i_pT2_%i", pTbin1, pTbin2));

      //______________________________________________________________________________________________________________________________
    }
  }

  //________________________________________________________________________________________________________

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(43);
  gStyle->SetLegendTextSize(30);

  //L-Lbar
  TH1F *L0_inv_mass_LLbar_projection_US[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_inv_mass_LLbar_projection_US[nPtBins_corr][nPtBins_corr];

  TH1F *L0_inv_mass_LLbar_projection_LS[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_inv_mass_LLbar_projection_LS[nPtBins_corr][nPtBins_corr];

  TH1F *L0_inv_mass_LLbar_projection[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_inv_mass_LLbar_projection[nPtBins_corr][nPtBins_corr];

  TH1F *L0_inv_mass_LLbar_projection_full[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar_inv_mass_LLbar_projection_full[nPtBins_corr][nPtBins_corr];

  //L-L
  TH1F *L01_inv_mass_LL_projection_US[nPtBins_corr][nPtBins_corr];
  TH1F *L02_inv_mass_LL_projection_US[nPtBins_corr][nPtBins_corr];

  TH1F *L01_inv_mass_LL_projection_LS[nPtBins_corr][nPtBins_corr];
  TH1F *L02_inv_mass_LL_projection_LS[nPtBins_corr][nPtBins_corr];

  TH1F *L01_inv_mass_LL_projection[nPtBins_corr][nPtBins_corr];
  TH1F *L02_inv_mass_LL_projection[nPtBins_corr][nPtBins_corr];

  TH1F *L01_inv_mass_LL_projection_full[nPtBins_corr][nPtBins_corr];
  TH1F *L02_inv_mass_LL_projection_full[nPtBins_corr][nPtBins_corr];

  //Lbar-Lbar
  TH1F *L0bar1_inv_mass_LbarLbar_projection_US[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar2_inv_mass_LbarLbar_projection_US[nPtBins_corr][nPtBins_corr];

  TH1F *L0bar1_inv_mass_LbarLbar_projection_LS[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar2_inv_mass_LbarLbar_projection_LS[nPtBins_corr][nPtBins_corr];

  TH1F *L0bar1_inv_mass_LbarLbar_projection[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar2_inv_mass_LbarLbar_projection[nPtBins_corr][nPtBins_corr];

  TH1F *L0bar1_inv_mass_LbarLbar_projection_full[nPtBins_corr][nPtBins_corr];
  TH1F *L0bar2_inv_mass_LbarLbar_projection_full[nPtBins_corr][nPtBins_corr];


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
      int binHigh = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(1.12);

      int binLow_sideBand = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(1.14);
      int binHigh_sideBand = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->FindBin(1.18);

      //L-Lbar
      //US and LS
      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_US_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_US_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_US_can->cd();

      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      float L_Lbar_US_ratio_integral = L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow, binHigh);
      L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Draw("surf1");

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
      L0_L0bar_text_Minv_2D_kine->AddText(Form("%.2f < p_{T,1} < %.2f GeV/c", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0bar_text_Minv_2D_kine->AddText(Form("%.2f < p_{T,2} < %.2f GeV/c", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0bar_text_Minv_2D_kine->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_2D_kine->Draw("same");

      TPaveText *L0_L0bar_text_STAR = new TPaveText(0.4, 0.9, 0.6, 0.98, "NDC");
      L0_L0bar_text_STAR->SetTextFont(43);
      L0_L0bar_text_STAR->SetTextSize(30);
      L0_L0bar_text_STAR->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_text_STAR->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0_L0bar_text_STAR->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_STAR->Draw("same");

      L0_inv_mass_vs_L0bar_inv_mass_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_US_zoom_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_US_zoom_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_zoom_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_US_zoom_can->cd();

      L0_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv,1}");
      L0_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      //L0_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv,2}");
      L0_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      //L0_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0bar_text_Minv_2D_US_zoom = new TPaveText(0.02, 0.84, 0.32, 0.98, "NDC");
      L0_L0bar_text_Minv_2D_US_zoom->SetTextFont(43);
      L0_L0bar_text_Minv_2D_US_zoom->SetTextSize(30);
      L0_L0bar_text_Minv_2D_US_zoom->AddText("STAR");
      //L0_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      //L0_L0bar_text_Minv_2D_US_zoom->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_text_Minv_2D_US_zoom->AddText(Form("p+p #sqrt{s} = %i GeV", energy));
      L0_L0bar_text_Minv_2D_US_zoom->AddText("Minimum bias");
      //L0_L0bar_text_Minv_2D_US_zoom->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      //L0_L0bar_text_Minv_2D_US_zoom->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0_L0bar_text_Minv_2D_US_zoom->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_2D_US_zoom->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0_L0bar_text_STAR->Draw("same");

      L0_inv_mass_vs_L0bar_inv_mass_US_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_US_zoom_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass_US_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_US_zoom_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_zoom_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_zoom_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_zoom_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_zoom_can->cd();

      L0_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->Add(L0_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[pTbin1][pTbin2], -1);
      L0_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0bar_text_Minv_2D_zoom = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0bar_text_Minv_2D_zoom->SetTextFont(42);
      //L0_L0bar_text_Minv_2D_zoom->AddText("STAR");
      //L0_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0bar_text_Minv_2D_zoom->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_text_Minv_2D_zoom->AddText("Minimum bias");
      L0_L0bar_text_Minv_2D_zoom->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_text_Minv_2D_zoom->AddText("Signal");
      L0_L0bar_text_Minv_2D_zoom->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_2D_zoom->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0bar_inv_mass_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0_inv_mass_vs_L0bar_inv_mass_zoom_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0_inv_mass_vs_L0bar_inv_mass_zoom_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_LS_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_LS_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_LS_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_LS_can->cd();

      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2]->Draw("surf1");

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

      L0_inv_mass_vs_L0bar_inv_mass_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_LS_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_LS_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_US_LS_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_can->cd();


      //background already scaled from Ana003_Lambda_corr_2D.cxx
      //here just potting
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Sumw2();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      float L_Lbar_US_LS_ratio_integral = L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow, binHigh);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Draw("surf1");

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

      L0_inv_mass_vs_L0bar_inv_mass_US_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));



      TCanvas *L0_inv_mass_LLbar_projection_US_can = new TCanvas(Form("L0_inv_mass_LLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_LLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_LLbar_projection_US_can->cd();

      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      //L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binHigh, -1); //projection in side-band for testing
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->Draw("p e");


      L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->SetMarkerStyle(24);
      L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->SetLineColor(kBlue);
      L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->Draw("p e same");

      TPaveText *L0_L0bar_text_Minv_L = new TPaveText(0.6, 0.55, 0.85, 0.9, "NDC");
      L0_L0bar_text_Minv_L->SetTextFont(42);
      //L0_L0bar_text_no_corr->AddText("STAR Internal");
      //L0_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0bar_text_Minv_L->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_text_Minv_L->AddText("Minimum bias");
      L0_L0bar_text_Minv_L->AddText("Projection to #Lambda^{0}");
      //L0_L0bar_text_Minv_L->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0_L0bar_text_Minv_L->AddText("|#it{y}| < 1");
      L0_L0bar_text_Minv_L->AddText(Form("%.2f < p_{T}^{1} < %.2f GeV/c", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0bar_text_Minv_L->AddText(Form("%.2f < p_{T}^{2} < %.2f GeV/c", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0bar_text_Minv_L->AddText(Form("Sig./Bckg. = %.2f", (L_Lbar_US_ratio_integral- L_Lbar_US_LS_ratio_integral)/L_Lbar_US_LS_ratio_integral));
      L0_L0bar_text_Minv_L->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_L->Draw("same");

      TLegend *L0_L0bar_leg_Minv_L = new TLegend(0.2, 0.35, 0.45, 0.54);
      L0_L0bar_leg_Minv_L->AddEntry(L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2], "US(#pi^{-}p)");
      L0_L0bar_leg_Minv_L->AddEntry(L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2], "Conbinatorial bckg.");
      L0_L0bar_leg_Minv_L->SetBorderSize(0);
      L0_L0bar_leg_Minv_L->Draw("same");

      L0_inv_mass_LLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_LLbar_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_LLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0_inv_mass_LLbar_projection_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L0bar_inv_mass_LLbar_projection_US_can = new TCanvas(Form("L0bar_inv_mass_LLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_LLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_LLbar_projection_US_can->cd();

      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->ProjectionY(Form("L0bar_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
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

      TPaveText *L0_L0bar_text_Minv_Lbar = new TPaveText(0.6, 0.55, 0.85, 0.9, "NDC");
      L0_L0bar_text_Minv_Lbar->SetTextFont(42);
      //L0_L0bar_text_no_corr->AddText("STAR Internal");
      //L0_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0bar_text_Minv_Lbar->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_text_Minv_Lbar->AddText("Minimum bias");
      L0_L0bar_text_Minv_Lbar->AddText("Projection to #bar{#Lambda^{0}}");
      //L0_L0bar_text_Minv_Lbar->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0_L0bar_text_Minv_Lbar->AddText("|#it{y}| < 1");
      L0_L0bar_text_Minv_Lbar->AddText(Form("%.2f < p_{T}^{1} < %.2f GeV/c", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0bar_text_Minv_Lbar->AddText(Form("%.2f < p_{T}^{2} < %.2f GeV/c", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0bar_text_Minv_Lbar->AddText(Form("Sig./Bckg. = %.2f", (L_Lbar_US_ratio_integral- L_Lbar_US_LS_ratio_integral)/L_Lbar_US_LS_ratio_integral));
      L0_L0bar_text_Minv_Lbar->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_Lbar->Draw("same");

      L0_L0bar_leg_Minv_L->Draw("same");

      L0bar_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->Draw("p e same");

      L0bar_inv_mass_LLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0bar_inv_mass_LLbar_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar_inv_mass_LLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US/L0bar_inv_mass_LLbar_projection_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      //------------------------------------------------------------------------------------------------------------------------------------

      //US-LS
      TF2 *doubleGauss_L_Lbar = new TF2("doubleGauss_L_Lbar", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 1.07, 1.2, 1.07, 1.2);
      if( pTbin1 == 0 && pTbin2 == 0 )
      {
        if( energy == 200 ) doubleGauss_L_Lbar->SetParameters(1200, 1.116, 0.002, 1.116, 0.002);
        if( energy == 510 ) doubleGauss_L_Lbar->SetParameters(3000, 1.116, 0.002, 1.116, 0.002);
      }
      else if( pTbin1 > 0 && fit_res_gaus_L0_L0Bar[pTbin1-1][0]->IsValid() ) doubleGauss_L_Lbar->SetParameters(fit_res_gaus_L0_L0Bar[pTbin1-1][0]->Parameter(0)/10, fit_res_gaus_L0_L0Bar[pTbin1-1][0]->Parameter(1), fit_res_gaus_L0_L0Bar[pTbin1-1][0]->Parameter(2), fit_res_gaus_L0_L0Bar[pTbin1-1][0]->Parameter(3), fit_res_gaus_L0_L0Bar[pTbin1-1][0]->Parameter(4));
      else if( pTbin2 > 0 && fit_res_gaus_L0_L0Bar[0][pTbin2-1]->IsValid() ) doubleGauss_L_Lbar->SetParameters(fit_res_gaus_L0_L0Bar[0][pTbin2-1]->Parameter(0)/10, fit_res_gaus_L0_L0Bar[0][pTbin2-1]->Parameter(1), fit_res_gaus_L0_L0Bar[0][pTbin2-1]->Parameter(2), fit_res_gaus_L0_L0Bar[0][pTbin2-1]->Parameter(3), fit_res_gaus_L0_L0Bar[0][pTbin2-1]->Parameter(4));
      else
      {
        cout<<"Fit not valid for L-Lbar pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
        return;
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
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      //L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Add(L0_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2], -1);
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Add(L0_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2], -1);
      fit_res_gaus_L0_L0Bar[pTbin1][pTbin2] = L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_L_Lbar, "S 0", "", 1.11, 1.125);
      L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Draw("surf1");


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

      L0_inv_mass_vs_L0bar_inv_mass_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0_inv_mass_vs_L0bar_inv_mass_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0_inv_mass_vs_L0bar_inv_mass_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //---------------------------------------------------------------------------------------

      TCanvas *L0_inv_mass_LLbar_projection_can = new TCanvas(Form("L0_inv_mass_LLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_LLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_LLbar_projection_can->cd();

      //int LLbar_projectionBin_L = L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(1));
      int LLbar_projectionBin_L = L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(3));

      //L0_inv_mass_LLbar_projection[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LLbar_projectionBin_L, LLbar_projectionBin_L);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_inv_mass_LLbar_projection[pTbin1][pTbin2]->Draw("p e");


      //TF12 *LLbar_Gaus_L0 = new TF12("LLbar_Gaus_L0", doubleGauss_L_Lbar, fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(1), "x");
      TF12 *LLbar_Gaus_L0 = new TF12("LLbar_Gaus_L0", doubleGauss_L_Lbar, fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(3), "x");
      LLbar_Gaus_L0->SetLineColor(1);
      LLbar_Gaus_L0->Draw("same");


      L0_L0bar_text_Minv_L->Draw("same");

      TLegend *L0_L0bar_leg_Minv_L_sig = new TLegend(0.6, 0.45, 0.85, 0.54);
      L0_L0bar_leg_Minv_L_sig->AddEntry(L0_inv_mass_LLbar_projection[pTbin1][pTbin2], "US - Bckg.");
      L0_L0bar_leg_Minv_L_sig->AddEntry(LLbar_Gaus_L0, "Gaussian fit");
      L0_L0bar_leg_Minv_L_sig->SetBorderSize(0);
      L0_L0bar_leg_Minv_L_sig->Draw("same");


      L0_inv_mass_LLbar_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0_inv_mass_LLbar_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_LLbar_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0_inv_mass_LLbar_projection_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L0bar_inv_mass_LLbar_projection_can = new TCanvas(Form("L0bar_inv_mass_LLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_LLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_LLbar_projection_can->cd();

      //int LLbar_projectionBin_Lbar = L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(3));
      int LLbar_projectionBin_Lbar = L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(1));

      //L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0barbar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0bar_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L0bar_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LLbar_projectionBin_Lbar, LLbar_projectionBin_Lbar);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar_inv_mass_LLbar_projection[pTbin1][pTbin2]->Draw("p e");


      //TF12 *LLbar_Gaus_L0bar = new TF12("LLbar_Gaus_L0bar", doubleGauss_L_Lbar, fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(3), "y");
      TF12 *LLbar_Gaus_L0bar = new TF12("LLbar_Gaus_L0bar", doubleGauss_L_Lbar, fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(1), "y");
      LLbar_Gaus_L0bar->SetLineColor(1);
      LLbar_Gaus_L0bar->Draw("same");

      L0_L0bar_text_Minv_Lbar->Draw("same");

      L0_L0bar_leg_Minv_L_sig->Draw("same");

      L0bar_inv_mass_LLbar_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0bar_inv_mass_LLbar_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar_inv_mass_LLbar_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0bar_inv_mass_LLbar_projection_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //---------------------------------------------------------------------------------------

      TCanvas *L0_inv_mass_LLbar_projection_full_can = new TCanvas(Form("L0_inv_mass_LLbar_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_LLbar_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_LLbar_projection_full_can->cd();

      int LLbar_projection_fullBin_L = L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(1));

      L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      //L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LLbar_projection_fullBin_L, LLbar_projection_fullBin_L);
      L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->SetMarkerSize(2);
      L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->SetLineColor(kRed);
      L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->Draw("p e");

      L0_L0bar_text_Minv_L->Draw("same");

      TLegend *L0_L0bar_leg_Minv_L_sig_full = new TLegend(0.6, 0.45, 0.85, 0.54);
      L0_L0bar_leg_Minv_L_sig_full->AddEntry(L0_inv_mass_LLbar_projection_full[pTbin1][pTbin2], "US - Bckg.");
      L0_L0bar_leg_Minv_L_sig_full->SetBorderSize(0);
      L0_L0bar_leg_Minv_L_sig_full->Draw("same");


      L0_inv_mass_LLbar_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0_inv_mass_LLbar_projection_full_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_LLbar_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0_inv_mass_LLbar_projection_full_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L0bar_inv_mass_LLbar_projection_full_can = new TCanvas(Form("L0bar_inv_mass_LLbar_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_LLbar_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_LLbar_projection_full_can->cd();

      int LLbar_projection_fullBin_Lbar = L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus_L0_L0Bar[pTbin1][pTbin2]->Parameter(3));

      L0bar_inv_mass_LLbar_projection_full[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L0bar_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      //L0bar_inv_mass_LLbar_projection_full[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L0bar_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LLbar_projection_fullBin_Lbar, LLbar_projection_fullBin_Lbar);
      L0bar_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L0bar_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L0bar_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar_inv_mass_LLbar_projection_full[pTbin1][pTbin2]->Draw("p e");


      L0_L0bar_text_Minv_Lbar->Draw("same");

      L0_L0bar_leg_Minv_L_sig_full->Draw("same");

      L0bar_inv_mass_LLbar_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0bar_inv_mass_LLbar_projection_full_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar_inv_mass_LLbar_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_Lbar/US-LS/L0bar_inv_mass_LLbar_projection_full_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //________________________________________________________________________________________________________

      //L-L
      //US
      TCanvas *L0_inv_mass_vs_L0_inv_mass_US_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_US_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_US_can->cd();

      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      float L_L_US_ratio_integral = L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow, binHigh);
      L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->Draw("surf1");


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

      L0_inv_mass_vs_L0_inv_mass_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L0_inv_mass_vs_L0_inv_mass_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L0_inv_mass_vs_L0_inv_mass_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L0_inv_mass_vs_L0_inv_mass_US_zoom_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_US_zoom_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_zoom_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_US_zoom_can->cd();

      L0_inv_mass_vs_L0_inv_mass_US_zoom[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv,1}");
      L0_inv_mass_vs_L0_inv_mass_US_zoom[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0_inv_mass_US_zoom[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      //L0_inv_mass_vs_L0_inv_mass_US_zoom[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0_inv_mass_US_zoom[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv,2}");
      L0_inv_mass_vs_L0_inv_mass_US_zoom[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0_inv_mass_US_zoom[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      //L0_inv_mass_vs_L0_inv_mass_US_zoom[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0_inv_mass_US_zoom[pTbin1][pTbin2]->Draw("surf1");


      TPaveText *L0_L0_text_Minv_2D_US_zoom = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0_text_Minv_2D_US_zoom->SetTextFont(42);
      //L0_L0_text_Minv_2D_US_zoom->AddText("STAR");
      //L0_L0_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_text_Minv_2D_US_zoom->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_text_Minv_2D_US_zoom->AddText("Minimum bias");
      L0_L0_text_Minv_2D_US_zoom->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_text_Minv_2D_US_zoom->AddText("US(#pi^{-}p) vs. US(#pi^{-}p)");
      L0_L0_text_Minv_2D_US_zoom->SetFillColorAlpha(0, 0.01);
      L0_L0_text_Minv_2D_US_zoom->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0_inv_mass_US_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L0_inv_mass_vs_L0_inv_mass_US_zoom_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass_US_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L0_inv_mass_vs_L0_inv_mass_US_zoom_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L0_inv_mass_vs_L0_inv_mass_zoom_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_zoom_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_zoom_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_zoom_can->cd();

      L0_inv_mass_vs_L0_inv_mass_US_zoom[pTbin1][pTbin2]->Add(L0_inv_mass_vs_L0_inv_mass_US_LS_zoom[pTbin1][pTbin2], -1);
      L0_inv_mass_vs_L0_inv_mass_US_zoom[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0_text_Minv_2D_zoom = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0_text_Minv_2D_zoom->SetTextFont(42);
      //L0_L0_text_Minv_2D_zoom->AddText("STAR");
      //L0_L0_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_text_Minv_2D_zoom->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_text_Minv_2D_zoom->AddText("Minimum bias");
      L0_L0_text_Minv_2D_zoom->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_text_Minv_2D_zoom->AddText("Signal");
      L0_L0_text_Minv_2D_zoom->SetFillColorAlpha(0, 0.01);
      L0_L0_text_Minv_2D_zoom->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0_inv_mass_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L0_inv_mass_vs_L0_inv_mass_zoom_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L0_inv_mass_vs_L0_inv_mass_zoom_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L0_inv_mass_vs_L0_inv_mass_US_LS_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_US_LS_can->cd();

      //background already scaled from Ana003_Lambda_corr_2D.cxx
      //here just potting
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->Sumw2();
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      float L_L_US_LS_ratio_integral = L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow, binHigh);
      L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2]->Draw("surf1");

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

      L0_inv_mass_vs_L0_inv_mass_US_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L0_inv_mass_vs_L0_inv_mass_US_LS_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass_US_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L0_inv_mass_vs_L0_inv_mass_US_LS_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));



      TCanvas *L01_inv_mass_LL_projection_US_can = new TCanvas(Form("L01_inv_mass_LL_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L01_inv_mass_LL_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L01_inv_mass_LL_projection_US_can->cd();

      L01_inv_mass_LL_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->ProjectionX(Form("L01_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L01_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L01_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L01_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L01_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
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


      TPaveText *L0_L0_text_Minv_L1 = new TPaveText(0.6, 0.55, 0.85, 0.9, "NDC");
      L0_L0_text_Minv_L1->SetTextFont(42);
      //L0_L0_text_no_corr->AddText("STAR Internal");
      //L0_L0_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_text_Minv_L1->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_text_Minv_L1->AddText("Minimum bias");
      L0_L0_text_Minv_L1->AddText("Projection to #Lambda^{0}_{1}");
      //L0_L0_text_Minv_L1->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0_L0_text_Minv_L1->AddText("|#it{y}| < 1");
      L0_L0_text_Minv_L1->AddText(Form("%.2f < p_{T}^{1} < %.2f GeV/c", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0_text_Minv_L1->AddText(Form("%.2f < p_{T}^{2} < %.2f GeV/c", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0_text_Minv_L1->AddText(Form("Sig./Bckg. = %.2f", (L_L_US_ratio_integral- L_L_US_LS_ratio_integral)/L_L_US_LS_ratio_integral));
      L0_L0_text_Minv_L1->SetFillColorAlpha(0, 0.01);
      L0_L0_text_Minv_L1->Draw("same");

      TLegend *L0_L0_leg_Minv_L1 = new TLegend(0.2, 0.35, 0.45, 0.54);
      L0_L0_leg_Minv_L1->AddEntry(L01_inv_mass_LL_projection_US[pTbin1][pTbin2], "US(#pi^{-}p)");
      L0_L0_leg_Minv_L1->AddEntry(L01_inv_mass_LL_projection_LS[pTbin1][pTbin2], "Conbinatorial bckg.");
      L0_L0_leg_Minv_L1->SetBorderSize(0);
      L0_L0_leg_Minv_L1->Draw("same");

      L01_inv_mass_LL_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L01_inv_mass_LL_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L01_inv_mass_LL_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L01_inv_mass_LL_projection_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L02_inv_mass_LL_projection_US_can = new TCanvas(Form("L02_inv_mass_LL_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L02_inv_mass_LL_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L02_inv_mass_LL_projection_US_can->cd();

      L02_inv_mass_LL_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->ProjectionY(Form("L02_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L02_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L02_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L02_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L02_inv_mass_LL_projection_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
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

      TPaveText *L0_L0_text_Minv_L2 = new TPaveText(0.6, 0.55, 0.85, 0.9, "NDC");
      L0_L0_text_Minv_L2->SetTextFont(42);
      //L0_L0_text_no_corr->AddText("STAR Internal");
      //L0_L0_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_text_Minv_L2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_text_Minv_L2->AddText("Minimum bias");
      L0_L0_text_Minv_L2->AddText("Projection to #Lambda^{0}_{2}");
      //L0_L0_text_Minv_L2->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0_L0_text_Minv_L2->AddText("|#it{y}| < 1");
      L0_L0_text_Minv_L2->AddText(Form("%.2f < p_{T}^{1} < %.2f GeV/c", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0_text_Minv_L2->AddText(Form("%.2f < p_{T}^{2} < %.2f GeV/c", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0_text_Minv_L2->AddText(Form("Sig./Bckg. = %.2f", (L_L_US_ratio_integral- L_L_US_LS_ratio_integral)/L_L_US_LS_ratio_integral));
      L0_L0_text_Minv_L2->SetFillColorAlpha(0, 0.01);
      L0_L0_text_Minv_L2->Draw("same");

      L0_L0_leg_Minv_L1->Draw("same");

      L02_inv_mass_LL_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L02_inv_mass_LL_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L02_inv_mass_LL_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US/L02_inv_mass_LL_projection_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      //US-LS
      TF2 *doubleGauss_L_L = new TF2("doubleGauss_L_L", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 1.07, 1.2, 1.07, 1.2);
      //doubleGauss_L_L->SetParameters(1200, 1.116, 0.002, 1.116, 0.002);
      if( pTbin1 == 0 && pTbin2 == 0 )
      {
        if( energy == 200 ) doubleGauss_L_L->SetParameters(1000, 1.116, 0.002, 1.116, 0.002);
        if( energy == 510 ) doubleGauss_L_L->SetParameters(1000, 1.116, 0.002, 1.116, 0.002);
      }
      else if( pTbin1 > 0 && fit_res_gaus_L0_L0[pTbin1-1][0]->IsValid() ) doubleGauss_L_L->SetParameters(fit_res_gaus_L0_L0[pTbin1-1][0]->Parameter(0)/10, fit_res_gaus_L0_L0[pTbin1-1][0]->Parameter(1), fit_res_gaus_L0_L0[pTbin1-1][0]->Parameter(2), fit_res_gaus_L0_L0[pTbin1-1][0]->Parameter(3), fit_res_gaus_L0_L0[pTbin1-1][0]->Parameter(4));
      else if( pTbin2 > 0 && fit_res_gaus_L0_L0[0][pTbin2-1]->IsValid() ) doubleGauss_L_L->SetParameters(fit_res_gaus_L0_L0[0][pTbin2-1]->Parameter(0)/10, fit_res_gaus_L0_L0[0][pTbin2-1]->Parameter(1), fit_res_gaus_L0_L0[0][pTbin2-1]->Parameter(2), fit_res_gaus_L0_L0[0][pTbin2-1]->Parameter(3), fit_res_gaus_L0_L0[0][pTbin2-1]->Parameter(4));
      else
      {
        cout<<"Fit not valid for L-L pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
        return;
      }

      TCanvas *L0_inv_mass_vs_L0_inv_mass_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_can->cd();

      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2] = (TH2F*)L0_inv_mass_vs_L0_inv_mass_US[pTbin1][pTbin2]->Clone(Form("L0_inv_mass_vs_L0_inv_mass_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07,1.2);
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      //L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->Add(L0_inv_mass_vs_L0_inv_mass_LS[pTbin1][pTbin2], -1);
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->Add(L0_inv_mass_vs_L0_inv_mass_US_LS[pTbin1][pTbin2], -1);
      fit_res_gaus_L0_L0[pTbin1][pTbin2] = L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_L_L, "s 0", "", 1.11, 1.125);
      L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->Draw("surf1");

      if( !fit_res_gaus_L0_L0[pTbin1][pTbin2]->IsValid())
      {
        cout<<"Fit not valid for L-L pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
        return;
      }

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

      L0_inv_mass_vs_L0_inv_mass_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L0_inv_mass_vs_L0_inv_mass_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L0_inv_mass_vs_L0_inv_mass_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //----------------------------------------------------------------------------------------------------------------------

      TCanvas *L01_inv_mass_LL_projection_can = new TCanvas(Form("L01_inv_mass_LL_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L01_inv_mass_LL_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L01_inv_mass_LL_projection_can->cd();

      //int LL_projectionBin_L1 = L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(1));
      int LL_projectionBin_L1 = L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(3));

      //L0_inv_mass_LL_projection[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L01_inv_mass_LL_projection[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L01_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LL_projectionBin_L1, LL_projectionBin_L1);
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      L01_inv_mass_LL_projection[pTbin1][pTbin2]->Draw("p e");


      //TF12 *LL_Gaus_L0 = new TF12("LL_Gaus_L0", doubleGauss_L_L, fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(1), "x");
      TF12 *LL_Gaus_L0 = new TF12("LL_Gaus_L0", doubleGauss_L_L, fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(3), "x");
      LL_Gaus_L0->SetLineColor(1);
      LL_Gaus_L0->Draw("same");

      L0_L0_text_Minv_L1->Draw("same");

      TLegend *L0_L0_leg_Minv_L_sig = new TLegend(0.6, 0.45, 0.85, 0.54);
      L0_L0_leg_Minv_L_sig->AddEntry(L01_inv_mass_LL_projection[pTbin1][pTbin2], "US - Bckg.");
      L0_L0_leg_Minv_L_sig->AddEntry(LL_Gaus_L0, "Gaussian fit");
      L0_L0_leg_Minv_L_sig->SetBorderSize(0);
      L0_L0_leg_Minv_L_sig->Draw("same");

      L01_inv_mass_LL_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L01_inv_mass_LL_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L01_inv_mass_LL_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L01_inv_mass_LL_projection_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));



      TCanvas *L02_inv_mass_LL_projection_can = new TCanvas(Form("L02_inv_mass_LL_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L02_inv_mass_LL_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L02_inv_mass_LL_projection_can->cd();

      //int LL_projectionBin_L2 = L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(3));
      int LL_projectionBin_L2 = L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(1));

      //L02_inv_mass_LL_projection[pTbin1][pTbin2] = (TH1F*)L02_inv_mass_vs_L02bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L02_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L02_inv_mass_LL_projection[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L02_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LL_projectionBin_L2, LL_projectionBin_L2);
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      L02_inv_mass_LL_projection[pTbin1][pTbin2]->Draw("p e");


      //TF12 *LL_Gaus_L02 = new TF12("LL_Gaus_L02", doubleGauss_L_L, fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(3), "y");
      TF12 *LL_Gaus_L02 = new TF12("LL_Gaus_L02", doubleGauss_L_L, fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(1), "y");
      LL_Gaus_L02->SetLineColor(1);
      LL_Gaus_L02->Draw("same");

      L0_L0_text_Minv_L2->Draw("same");

      L0_L0_leg_Minv_L_sig->Draw("same");

      L02_inv_mass_LL_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L02_inv_mass_LL_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L02_inv_mass_LL_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L02_inv_mass_LL_projection_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      //----------------------------------------------------------------------------------------------------------------------

      TCanvas *L01_inv_mass_LL_projection_full_can = new TCanvas(Form("L01_inv_mass_LL_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L01_inv_mass_LL_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L01_inv_mass_LL_projection_full_can->cd();

      int LL_projection_fullBin_L1 = L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(1));

      L01_inv_mass_LL_projection_full[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      //L01_inv_mass_LL_projection_full[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L01_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LL_projection_fullBin_L1, LL_projection_fullBin_L1);
      L01_inv_mass_LL_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L01_inv_mass_LL_projection_full[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L01_inv_mass_LL_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L01_inv_mass_LL_projection_full[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L01_inv_mass_LL_projection_full[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L01_inv_mass_LL_projection_full[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L01_inv_mass_LL_projection_full[pTbin1][pTbin2]->SetMarkerStyle(20);
      L01_inv_mass_LL_projection_full[pTbin1][pTbin2]->SetMarkerSize(2);
      L01_inv_mass_LL_projection_full[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L01_inv_mass_LL_projection_full[pTbin1][pTbin2]->SetLineColor(kRed);
      L01_inv_mass_LL_projection_full[pTbin1][pTbin2]->Draw("p e");

      L0_L0_text_Minv_L1->Draw("same");


      TLegend *L0_L0_leg_Minv_L_sig_full = new TLegend(0.6, 0.45, 0.85, 0.54);
      L0_L0_leg_Minv_L_sig_full->AddEntry(L01_inv_mass_LL_projection_full[pTbin1][pTbin2], "US - Bckg.");
      L0_L0_leg_Minv_L_sig_full->SetBorderSize(0);
      L0_L0_leg_Minv_L_sig_full->Draw("same");

      L01_inv_mass_LL_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L01_inv_mass_LL_projection_full_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L01_inv_mass_LL_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L01_inv_mass_LL_projection_full_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));



      TCanvas *L02_inv_mass_LL_projection_full_can = new TCanvas(Form("L02_inv_mass_LL_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L02_inv_mass_LL_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L02_inv_mass_LL_projection_full_can->cd();

      int LL_projection_fullBin_L2 = L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus_L0_L0[pTbin1][pTbin2]->Parameter(3));

      L02_inv_mass_LL_projection_full[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L02_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      //L02_inv_mass_LL_projection_full[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L02_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LL_projection_fullBin_L2, LL_projection_fullBin_L2);
      L02_inv_mass_LL_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L02_inv_mass_LL_projection_full[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L02_inv_mass_LL_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L02_inv_mass_LL_projection_full[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L02_inv_mass_LL_projection_full[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L02_inv_mass_LL_projection_full[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L02_inv_mass_LL_projection_full[pTbin1][pTbin2]->SetMarkerStyle(20);
      L02_inv_mass_LL_projection_full[pTbin1][pTbin2]->SetMarkerSize(2);
      L02_inv_mass_LL_projection_full[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L02_inv_mass_LL_projection_full[pTbin1][pTbin2]->SetLineColor(kRed);
      L02_inv_mass_LL_projection_full[pTbin1][pTbin2]->Draw("p e");


      L0_L0_text_Minv_L2->Draw("same");

      L0_L0_leg_Minv_L_sig_full->Draw("same");

      L02_inv_mass_LL_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L02_inv_mass_LL_projection_full_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L02_inv_mass_LL_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/L_L/US-LS/L02_inv_mass_LL_projection_full_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //________________________________________________________________________________________________________

      //Lbar-Lbar

      //US
      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_US_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_can->cd();

      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      float Lbar_Lbar_US_ratio_integral = L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow, binHigh);
      L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Draw("surf1");

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

      L0bar_inv_mass_vs_L0bar_inv_mass_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom_can->cd();

      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv,1}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      //L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv,2}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      //L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0bar_L0bar_text_Minv_2D_US_zoom = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0bar_L0bar_text_Minv_2D_US_zoom->SetTextFont(42);
      //L0bar_L0bar_text_Minv_2D_US_zoom->AddText("STAR");
      //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_text_Minv_2D_US_zoom->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_text_Minv_2D_US_zoom->AddText("Minimum bias");
      L0bar_L0bar_text_Minv_2D_US_zoom->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_text_Minv_2D_US_zoom->AddText("US(#pi^{+}#bar{p}) vs. US(#pi^{+}#bar{p})");
      L0bar_L0bar_text_Minv_2D_US_zoom->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_Minv_2D_US_zoom->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_zoom_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_zoom_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_zoom_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_zoom_can->cd();

      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->Add(L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_zoom[pTbin1][pTbin2], -1);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_zoom[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0bar_L0bar_text_Minv_2D_zoom = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0bar_L0bar_text_Minv_2D_zoom->SetTextFont(42);
      //L0bar_L0bar_text_Minv_2D_zoom->AddText("STAR");
      //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_text_Minv_2D_zoom->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_text_Minv_2D_zoom->AddText("Minimum bias");
      L0bar_L0bar_text_Minv_2D_zoom->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_text_Minv_2D_zoom->AddText("Signal");
      L0bar_L0bar_text_Minv_2D_zoom->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_Minv_2D_zoom->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0bar_inv_mass_vs_L0bar_inv_mass_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar_inv_mass_vs_L0bar_inv_mass_zoom_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass_zoom_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar_inv_mass_vs_L0bar_inv_mass_zoom_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can->cd();

      //background already scaled from Ana003_Lambda_corr_2D.cxx
      //here just potting

      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Sumw2();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      float Lbar_Lbar_US_LS_ratio_integral = L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow, binHigh);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->Draw("surf1");

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

      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));



      TCanvas *L0bar1_inv_mass_LbarLbar_projection_US_can = new TCanvas(Form("L0bar1_inv_mass_LbarLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar1_inv_mass_LbarLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar1_inv_mass_LbarLbar_projection_US_can->cd();

      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->ProjectionX(Form("L0bar1_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
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

      TPaveText *L0bar_L0bar_text_Minv_L1 = new TPaveText(0.6, 0.55, 0.85, 0.9, "NDC");
      L0bar_L0bar_text_Minv_L1->SetTextFont(42);
      //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
      //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_text_Minv_L1->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_text_Minv_L1->AddText("Minimum bias");
      L0bar_L0bar_text_Minv_L1->AddText("Projection to #bar{#Lambda^{0}}_{1}");
      //L0bar_L0bar_text_Minv_L1->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0bar_L0bar_text_Minv_L1->AddText("|#it{y}| < 1");
      L0bar_L0bar_text_Minv_L1->AddText(Form("%.2f < p_{T}^{1} < %.2f GeV/c", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0bar_L0bar_text_Minv_L1->AddText(Form("%.2f < p_{T}^{2} < %.2f GeV/c", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0bar_L0bar_text_Minv_L1->AddText(Form("Sig./Bckg. = %.2f", (Lbar_Lbar_US_ratio_integral- Lbar_Lbar_US_LS_ratio_integral)/Lbar_Lbar_US_LS_ratio_integral));
      L0bar_L0bar_text_Minv_L1->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_Minv_L1->Draw("same");

      TLegend *L0bar_L0bar_leg_Minv_L1 = new TLegend(0.2, 0.35, 0.45, 0.54);
      L0bar_L0bar_leg_Minv_L1->AddEntry(L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2], "US(#pi^{+}#bar{p})");
      L0bar_L0bar_leg_Minv_L1->AddEntry(L0bar1_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2], "Conbinatorial bckg.");
      L0bar_L0bar_leg_Minv_L1->SetBorderSize(0);
      L0bar_L0bar_leg_Minv_L1->Draw("same");


      L0bar1_inv_mass_LbarLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar1_inv_mass_LbarLbar_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar1_inv_mass_LbarLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar1_inv_mass_LbarLbar_projection_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));



      TCanvas *L0bar2_inv_mass_LbarLbar_projection_US_can = new TCanvas(Form("L0bar2_inv_mass_LbarLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar2_inv_mass_LbarLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar2_inv_mass_LbarLbar_projection_US_can->cd();

      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->ProjectionY(Form("L0bar2_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->Draw("p e");

      //L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0bar2_inv_mass_vs_L0bar2_inv_mass_LS[pTbin1][pTbin2]->ProjectionX(Form("L0bar2_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2]->ProjectionY(Form("L0bar2_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->SetMarkerStyle(24);
      L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->SetMarkerColor(kBlue);
      L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->SetLineColor(kBlue);
      L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2]->Draw("p e same");

      TPaveText *L0bar_L0bar_text_Minv_L2 = new TPaveText(0.6, 0.55, 0.85, 0.9, "NDC");
      L0bar_L0bar_text_Minv_L2->SetTextFont(42);
      //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
      //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_text_Minv_L2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_text_Minv_L2->AddText("Minimum bias");
      L0bar_L0bar_text_Minv_L2->AddText("Projection to #bar{#Lambda^{0}}_{2}");
      //L0bar_L0bar_text_Minv_L2->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0bar_L0bar_text_Minv_L2->AddText("|#it{y}| < 1");
      L0bar_L0bar_text_Minv_L2->AddText(Form("%.2f < p_{T}^{1} < %.2f GeV/c", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0bar_L0bar_text_Minv_L2->AddText(Form("%.2f < p_{T}^{2} < %.2f GeV/c", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0bar_L0bar_text_Minv_L2->AddText(Form("Sig./Bckg. = %.2f", (Lbar_Lbar_US_ratio_integral- Lbar_Lbar_US_LS_ratio_integral)/Lbar_Lbar_US_LS_ratio_integral));
      L0bar_L0bar_text_Minv_L2->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_Minv_L2->Draw("same");

      L0bar_L0bar_leg_Minv_L1->Draw("same");

      L0bar2_inv_mass_LbarLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar2_inv_mass_LbarLbar_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar2_inv_mass_LbarLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US/L0bar2_inv_mass_LbarLbar_projection_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //US-LS
      TF2 *doubleGauss_Lbar_Lbar = new TF2("doubleGauss_Lbar_Lbar", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", 1.07, 1.2, 1.07, 1.2);
      //doubleGauss_Lbar_Lbar->SetParameters(1200, 1.116, 0.002, 1.116, 0.002);

      if( pTbin1 == 0 && pTbin2 == 0 )
      {
        if( energy == 200 && year == 2012 )
        {
          if(cut_type == 2) doubleGauss_Lbar_Lbar->SetParameters(200, 1.116, 0.002, 1.116, 0.002);
          else doubleGauss_Lbar_Lbar->SetParameters(1000, 1.116, 0.002, 1.116, 0.002);
        }

        if( energy == 200 && year == 2015 ) doubleGauss_Lbar_Lbar->SetParameters(200, 1.116, 0.002, 1.116, 0.002);
        if( energy == 200 && year == 2016 ) doubleGauss_Lbar_Lbar->SetParameters(600, 1.116, 0.002, 1.116, 0.002);
        if( energy == 510 ) doubleGauss_Lbar_Lbar->SetParameters(600, 1.116, 0.002, 1.116, 0.002);
      }
      //else if( pTbin1 == 1 && pTbin2 == 0 ) doubleGauss_Lbar_Lbar->SetParameters(1000, 1.116, 0.002, 1.116, 0.002);
      else if( pTbin1 > 0 && fit_res_gaus_L0Bar_L0Bar[pTbin1-1][0]->IsValid() ) doubleGauss_Lbar_Lbar->SetParameters(fit_res_gaus_L0Bar_L0Bar[pTbin1-1][0]->Parameter(0)/10, fit_res_gaus_L0Bar_L0Bar[pTbin1-1][0]->Parameter(1), fit_res_gaus_L0Bar_L0Bar[pTbin1-1][0]->Parameter(2), fit_res_gaus_L0Bar_L0Bar[pTbin1-1][0]->Parameter(3), fit_res_gaus_L0Bar_L0Bar[pTbin1-1][0]->Parameter(4));
      else if( pTbin2 > 0 && fit_res_gaus_L0Bar_L0Bar[0][pTbin2-1]->IsValid() ) doubleGauss_Lbar_Lbar->SetParameters(fit_res_gaus_L0Bar_L0Bar[0][pTbin2-1]->Parameter(0)/10, fit_res_gaus_L0Bar_L0Bar[0][pTbin2-1]->Parameter(1), fit_res_gaus_L0Bar_L0Bar[0][pTbin2-1]->Parameter(2), fit_res_gaus_L0Bar_L0Bar[0][pTbin2-1]->Parameter(3), fit_res_gaus_L0Bar_L0Bar[0][pTbin2-1]->Parameter(4));
      else
      {
        cout<<"Fit not valid for Lbar-Lbar pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
        return;
      }

      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_can->cd();


      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2] = (TH2F*)L0bar_inv_mass_vs_L0bar_inv_mass_US[pTbin1][pTbin2]->Clone(Form("L0bar_inv_mass_vs_L0bar_inv_mass_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      //L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Add(L0bar_inv_mass_vs_L0bar_inv_mass_LS[pTbin1][pTbin2], -1);
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Add(L0bar_inv_mass_vs_L0bar_inv_mass_US_LS[pTbin1][pTbin2], -1);
      //if( pTbin1 == 0 && pTbin2 == 1 ) fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2] = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_Lbar_Lbar, "s 0", "", 1.1, 1.12);
      //else if( pTbin1 == 1 && pTbin2 == 1 ) fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2] = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_Lbar_Lbar, "s 0", "", 1.1, 1.12);
      //else fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2] = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_Lbar_Lbar, "s 0", "", 1.11, 1.125);
      fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2] = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Fit(doubleGauss_Lbar_Lbar, "s 0", "", 1.113, 1.125);
      L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->Draw("surf1");

      if( !fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->IsValid())
      {
        cout<<"Fit not valid for Lbar-Lbar pT1 "<<pTbin1<<", pT2 "<<pTbin2<<endl;
        return;
      }

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

      L0bar_inv_mass_vs_L0bar_inv_mass_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar_inv_mass_vs_L0bar_inv_mass_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar_inv_mass_vs_L0bar_inv_mass_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //---------------------------------------------------------------------------------------------------------------------------------

      TCanvas *L0bar1_inv_mass_LbarLbar_projection_can = new TCanvas(Form("L0bar1_inv_mass_LbarLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar1_inv_mass_LbarLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar1_inv_mass_LbarLbar_projection_can->cd();

      //int LbarLbar_projectionBin_Lbar1 = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(1));
      int LbarLbar_projectionBin_Lbar1 = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(3));

      //L0bar_inv_mass_LbarLbar_projection[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0barbar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0bar_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0bar1_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LbarLbar_projectionBin_Lbar1, LbarLbar_projectionBin_Lbar1);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->Draw("p e");


      //TF12 *LbarLbar_Gaus_L0bar = new TF12("LbarLbar_Gaus_L0bar", doubleGauss_Lbar_Lbar, fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(1), "x");
      TF12 *LbarLbar_Gaus_L0bar = new TF12("LbarLbar_Gaus_L0bar", doubleGauss_Lbar_Lbar, fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(3), "x");
      LbarLbar_Gaus_L0bar->SetLineColor(1);
      LbarLbar_Gaus_L0bar->Draw("same");

      L0bar_L0bar_text_Minv_L1->Draw("same");

      TLegend *L0bar_L0bar_leg_Minv_L_sig = new TLegend(0.6, 0.45, 0.85, 0.54);
      L0bar_L0bar_leg_Minv_L_sig->AddEntry(L0bar1_inv_mass_LbarLbar_projection[pTbin1][pTbin2], "US - Bckg.");
      L0bar_L0bar_leg_Minv_L_sig->AddEntry(LbarLbar_Gaus_L0bar, "Gaussian fit");
      L0bar_L0bar_leg_Minv_L_sig->SetBorderSize(0);
      L0bar_L0bar_leg_Minv_L_sig->Draw("same");

      L0bar1_inv_mass_LbarLbar_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar1_inv_mass_LbarLbar_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar1_inv_mass_LbarLbar_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar1_inv_mass_LbarLbar_projection_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));



      TCanvas *L0bar2_inv_mass_LbarLbar_projection_can = new TCanvas(Form("L0bar2_inv_mass_LbarLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar2_inv_mass_LbarLbar_projection_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar2_inv_mass_LbarLbar_projection_can->cd();

      //int LbarLbar_projectionBin_Lbar2 = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(3));
      int LbarLbar_projectionBin_Lbar2 = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(1));

      //L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L0bar2_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LbarLbar_projectionBin_Lbar2, LbarLbar_projectionBin_Lbar2);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L0bar2_inv_mass_projection_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LbarLbar_projectionBin_Lbar1, LbarLbar_projectionBin_Lbar1);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar2_inv_mass_LbarLbar_projection[pTbin1][pTbin2]->Draw("p e");


      //TF12 *LbarLbar_Gaus_L0bar2 = new TF12("LbarLbar_Gaus_L0bar2", doubleGauss_Lbar_Lbar, fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(3), "y");
      TF12 *LbarLbar_Gaus_L0bar2 = new TF12("LbarLbar_Gaus_L0bar2", doubleGauss_Lbar_Lbar, fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(1), "y");
      LbarLbar_Gaus_L0bar2->SetLineColor(1);
      LbarLbar_Gaus_L0bar2->Draw("same");

      L0bar_L0bar_text_Minv_L2->Draw("same");

      L0bar_L0bar_leg_Minv_L_sig->Draw("same");

      L0bar2_inv_mass_LbarLbar_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar2_inv_mass_LbarLbar_projection_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar2_inv_mass_LbarLbar_projection_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar2_inv_mass_LbarLbar_projection_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //---------------------------------------------------------------------------------------------------------------------------------

      TCanvas *L0bar1_inv_mass_LbarLbar_projection_full_can = new TCanvas(Form("L0bar1_inv_mass_LbarLbar_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar1_inv_mass_LbarLbar_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar1_inv_mass_LbarLbar_projection_full_can->cd();

      int LbarLbar_projection_fullBin_Lbar1 = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetXaxis()->FindBin(fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(1));

      L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0bar_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      //L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionX(Form("L0bar1_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LbarLbar_projection_fullBin_Lbar1, LbarLbar_projection_fullBin_Lbar1);
      L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->Draw("p e");

      L0bar_L0bar_text_Minv_L1->Draw("same");


      TLegend *L0bar_L0bar_leg_Minv_L_sig_full = new TLegend(0.6, 0.45, 0.85, 0.54);
      L0bar_L0bar_leg_Minv_L_sig_full->AddEntry(L0bar1_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2], "US - Bckg.");
      L0bar_L0bar_leg_Minv_L_sig_full->SetBorderSize(0);
      L0bar_L0bar_leg_Minv_L_sig_full->Draw("same");

      L0bar1_inv_mass_LbarLbar_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar1_inv_mass_LbarLbar_projection_full_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar1_inv_mass_LbarLbar_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar1_inv_mass_LbarLbar_projection_full_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));



      TCanvas *L0bar2_inv_mass_LbarLbar_projection_full_can = new TCanvas(Form("L0bar2_inv_mass_LbarLbar_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar2_inv_mass_LbarLbar_projection_full_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar2_inv_mass_LbarLbar_projection_full_can->cd();

      int LbarLbar_projection_fullBin_Lbar2 = L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->GetYaxis()->FindBin(fit_res_gaus_L0Bar_L0Bar[pTbin1][pTbin2]->Parameter(3));

      L0bar2_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L0bar2_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LbarLbar_projection_fullBin_Lbar2, LbarLbar_projection_fullBin_Lbar2);
      //L0bar2_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass[pTbin1][pTbin2]->ProjectionY(Form("L0bar2_inv_mass_projection_full_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), LbarLbar_projection_fullBin_Lbar1, LbarLbar_projection_fullBin_Lbar1);
      L0bar2_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar2_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar2_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(1.1);
      L0bar2_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.1,1.13);
      L0bar2_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar2_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar2_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar2_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar2_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar2_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar2_inv_mass_LbarLbar_projection_full[pTbin1][pTbin2]->Draw("p e");


      L0bar_L0bar_text_Minv_L2->Draw("same");

      L0bar_L0bar_leg_Minv_L_sig_full->Draw("same");

      L0bar2_inv_mass_LbarLbar_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar2_inv_mass_LbarLbar_projection_full_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar2_inv_mass_LbarLbar_projection_full_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv/Lbar_Lbar/US-LS/L0bar2_inv_mass_LbarLbar_projection_full_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //________________________________________________________________________________________________________

    }//end for pT bin 2
  }//end for pT bin 1

  InvMassFile->Close();

  return;
}
