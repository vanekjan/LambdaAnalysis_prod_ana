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

void Ana003_Lambda_corr_2D_plot_Minv_ME(const int cut_type = 0, const int energy = 510, const int year = 2017)
{


  TFile *InvMassFile; //ME invariant mass histograms are in the correlation file

  if(cut_type == 0)
  {
    if(year == 2012)
    {
      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_work.root", year), "read");

      InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/2024_07_new_ME_new_reweight/ProdPlane_Lambda_ana_cuts_new_ME_reweight.root", year), "read"); //analysis cuts

      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/20240729_tests_5sigma_signal/ProdPlane_Lambda_ana_cuts_2_sigma_signal.root", year), "read"); //analysis cuts

      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_012024.root", year), "read");
    }
    else if(year == 2015)
    {
      InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_work.root", year), "read");
    }
    else
    {
      InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_work.root", year), "read");
    }
  }

  if(cut_type == 1) InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_work.root", year), "read");
  if(cut_type == 2) InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_work.root", year), "read");

  if( !(InvMassFile->IsOpen()) )
  {
    cout<<"Unable to open file with invariant mass histograms!"<<endl;
    return;
  }


  //_______________________________________________________________________________________________________________________________________________


  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_ME[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  TH2F *L0_inv_mass_vs_L0bar_inv_mass_ME[nPtBins_corr][nPtBins_corr];


  TH2F *L0_inv_mass_vs_L0_inv_mass_US_ME[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0_inv_mass_vs_L0_inv_mass_US_LS_ME[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  TH2F *L0_inv_mass_vs_L0_inv_mass_ME[nPtBins_corr][nPtBins_corr];


  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[nPtBins_corr][nPtBins_corr]; //for US-US Lambda pairs
  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[nPtBins_corr][nPtBins_corr]; //for US-LS Lambda pairs

  TH2F *L0bar_inv_mass_vs_L0bar_inv_mass_ME[nPtBins_corr][nPtBins_corr];


  //----------------------------------------------------------------------------


  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2] =  (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_US_ME_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME_pT1_%i_pT2_%i", pTbin1, pTbin2));


      L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_US_ME_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_ME_pT1_%i_pT2_%i", pTbin1, pTbin2));


      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_ME_pT1_%i_pT2_%i", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2] = (TH2F*)InvMassFile->Get(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME_pT1_%i_pT2_%i", pTbin1, pTbin2));
    }
  }

  //________________________________________________________________________________________________________

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(42);

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



  for(unsigned int pTbin1 = 0; pTbin1 < nPtBins_corr; pTbin1++)
  {
    for(unsigned int pTbin2 = 0; pTbin2 < nPtBins_corr; pTbin2++)
    {
      //find bins for side band and for projection to x
      //projection bins just for testing, later can do projections with parameters of 2D fit
      int binLow = L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->FindBin(1.11);
      int binHigh = L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->FindBin(1.12);

      int binLow_sideBand = L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->FindBin(1.125);
      int binHigh_sideBand = L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->FindBin(1.13);

      //L-Lbar
      //US and LS
      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_US_ME_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_US_ME_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_ME_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_US_ME_can->cd();

      L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      float L_Lbar_US_scale_integral = L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow_sideBand, binLow_sideBand);
      float L_Lbar_US_ratio_integral = L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow, binHigh);
      L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0bar_text_Minv_2D_US_ME = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0bar_text_Minv_2D_US_ME->SetTextFont(42);
      //L0_L0bar_text_no_corr->AddText("STAR Internal");
      //L0_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0bar_text_Minv_2D_US_ME->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_text_Minv_2D_US_ME->AddText("Mixed event");
      L0_L0bar_text_Minv_2D_US_ME->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_text_Minv_2D_US_ME->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0_L0bar_text_Minv_2D_US_ME->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_2D_US_ME->Draw("same");


      TPaveText *L0_L0bar_text_Minv_2D_kine = new TPaveText(0.73, 0.83, 0.98, 0.98, "NDC");
      L0_L0bar_text_Minv_2D_kine->SetTextFont(42);
      L0_L0bar_text_Minv_2D_kine->AddText("|#it{y}| < 1");
      L0_L0bar_text_Minv_2D_kine->AddText(Form("%.2f < p_{T,1} < %.2f GeV/c", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0bar_text_Minv_2D_kine->AddText(Form("%.2f < p_{T,2} < %.2f GeV/c", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0bar_text_Minv_2D_kine->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0bar_inv_mass_US_ME_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass_US_ME_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      //------------------------------------------------------------------

      TCanvas *L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME_can = new TCanvas(Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME_can->cd();


      //background already scaled from Ana003_Lambda_corr_2D.cxx
      //here just potting
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->Sumw2();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      float L_Lbar_US_LS_scale_integral = L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow_sideBand, binLow_sideBand);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->Scale(L_Lbar_US_scale_integral/L_Lbar_US_LS_scale_integral);
      float L_Lbar_US_LS_ratio_integral = L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow, binHigh);
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0bar_text_Minv_2D_US_LS_ME = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0bar_text_Minv_2D_US_LS_ME->SetTextFont(42);
      //L0_L0bar_text_no_corr->AddText("STAR Internal");
      //L0_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0bar_text_Minv_2D_US_LS_ME->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0bar_text_Minv_2D_US_LS_ME->AddText("Mixed event");
      L0_L0bar_text_Minv_2D_US_LS_ME->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
      L0_L0bar_text_Minv_2D_US_LS_ME->AddText("US(#pi^{-}p) vs. LS(#pi^{-}#bar{p})");
      L0_L0bar_text_Minv_2D_US_LS_ME->AddText("LS(#pi^{+}p) vs. US(#pi^{+}#bar{p})");
      L0_L0bar_text_Minv_2D_US_LS_ME->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_2D_US_LS_ME->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_Lbar/US/L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //-----------------------------------------------------------------------------

      TCanvas *L0_inv_mass_LLbar_projection_US_can = new TCanvas(Form("L0_inv_mass_LLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0_inv_mass_LLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_LLbar_projection_US_can->cd();

      L0_inv_mass_LLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
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


      L0_inv_mass_LLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->ProjectionX(Form("L0_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
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
      L0_L0bar_text_Minv_L->AddText("Mixed event");
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

      L0_inv_mass_LLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_Lbar/US/L0_inv_mass_LLbar_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_LLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_Lbar/US/L0_inv_mass_LLbar_projection_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //-----------------------------------------------------------

      TCanvas *L0bar_inv_mass_LLbar_projection_US_can = new TCanvas(Form("L0bar_inv_mass_LLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_LLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_LLbar_projection_US_can->cd();

      L0bar_inv_mass_LLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->ProjectionY(Form("L0bar_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
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
      L0bar_inv_mass_LLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->ProjectionY(Form("L0bar_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
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
      L0_L0bar_text_Minv_Lbar->AddText("Mixed event");
      L0_L0bar_text_Minv_Lbar->AddText("Projection to #bar{#Lambda^{0}}");
      //L0_L0bar_text_Minv_Lbar->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0_L0bar_text_Minv_Lbar->AddText("|#it{y}| < 1");
      L0_L0bar_text_Minv_Lbar->AddText(Form("%.2f < p_{T}^{1} < %.2f GeV/c", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0bar_text_Minv_Lbar->AddText(Form("%.2f < p_{T}^{2} < %.2f GeV/c", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0bar_text_Minv_Lbar->AddText(Form("Sig./Bckg. = %.2f", (L_Lbar_US_ratio_integral - L_Lbar_US_LS_ratio_integral)/L_Lbar_US_LS_ratio_integral));
      L0_L0bar_text_Minv_Lbar->SetFillColorAlpha(0, 0.01);
      L0_L0bar_text_Minv_Lbar->Draw("same");

      L0_L0bar_leg_Minv_L->Draw("same");

      L0bar_inv_mass_LLbar_projection_LS[pTbin1][pTbin2]->Draw("p e same");

      L0bar_inv_mass_LLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_Lbar/US/L0bar_inv_mass_LLbar_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar_inv_mass_LLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_Lbar/US/L0bar_inv_mass_LLbar_projection_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      //________________________________________________________________________________________________________

      //L-L
      //US
      TCanvas *L0_inv_mass_vs_L0_inv_mass_US_ME_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_US_ME_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_ME_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_US_ME_can->cd();

      L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      float L_L_US_scale_integral = L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow_sideBand, binLow_sideBand);
      float L_L_US_ratio_integral = L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow, binHigh);
      L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2]->Draw("surf1");


      TPaveText *L0_L0_text_Minv_2D_US_ME = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0_text_Minv_2D_US_ME->SetTextFont(42);
      //L0_L0_text_no_corr->AddText("STAR Internal");
      //L0_L0_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_text_Minv_2D_US_ME->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_text_Minv_2D_US_ME->AddText("Mixed event");
      L0_L0_text_Minv_2D_US_ME->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_text_Minv_2D_US_ME->AddText("US(#pi^{-}p) vs. US(#pi^{-}p)");
      L0_L0_text_Minv_2D_US_ME->SetFillColorAlpha(0, 0.01);
      L0_L0_text_Minv_2D_US_ME->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0_inv_mass_US_ME_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_L/US/L0_inv_mass_vs_L0_inv_mass_US_ME_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass_US_ME_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_L/US/L0_inv_mass_vs_L0_inv_mass_US_ME_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //---------------------------------------------------------

      TCanvas *L0_inv_mass_vs_L0_inv_mass_US_LS_ME_can = new TCanvas(Form("L0_inv_mass_vs_L0_inv_mass_US_LS_ME_can_%i_%i", pTbin1, pTbin2), Form("L0_inv_mass_vs_L0_inv_mass_US_LS_ME_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME_can->cd();


      //here just potting
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->Sumw2();
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      float L_L_US_LS_scale_integral = L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow_sideBand, binLow_sideBand);
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->Scale(L_L_US_scale_integral/L_L_US_LS_scale_integral);
      float L_L_US_LS_ratio_integral = L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow, binHigh);
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0_L0_text_Minv_2D_US_LS_ME = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0_L0_text_Minv_2D_US_LS_ME->SetTextFont(42);
      //L0_L0_text_no_corr->AddText("STAR Internal");
      //L0_L0_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0_L0_text_Minv_2D_US_LS_ME->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0_L0_text_Minv_2D_US_LS_ME->AddText("Mixed event");
      L0_L0_text_Minv_2D_US_LS_ME->AddText("#Lambda^{0}-#Lambda^{0}");
      L0_L0_text_Minv_2D_US_LS_ME->AddText("US(#pi^{-}p) vs. LS(#pi^{+}p)");
      L0_L0_text_Minv_2D_US_LS_ME->AddText("LS(#pi^{+}p) vs. US(#pi^{-}p)");
      L0_L0_text_Minv_2D_US_LS_ME->SetFillColorAlpha(0, 0.01);
      L0_L0_text_Minv_2D_US_LS_ME->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0_inv_mass_vs_L0_inv_mass_US_LS_ME_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_L/US/L0_inv_mass_vs_L0_inv_mass_US_LS_ME_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0_inv_mass_vs_L0_inv_mass_US_LS_ME_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_L/US/L0_inv_mass_vs_L0_inv_mass_US_LS_ME_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //---------------------------------------------------------

      TCanvas *L01_inv_mass_LL_projection_US_can = new TCanvas(Form("L01_inv_mass_LL_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L01_inv_mass_LL_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L01_inv_mass_LL_projection_US_can->cd();

      L01_inv_mass_LL_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2]->ProjectionX(Form("L01_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
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
      L01_inv_mass_LL_projection_LS[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->ProjectionX(Form("L01_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
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
      L0_L0_text_Minv_L1->AddText("Mixed event");
      L0_L0_text_Minv_L1->AddText("Projection to #Lambda^{0}_{1}");
      //L0_L0_text_Minv_L1->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0_L0_text_Minv_L1->AddText("|#it{y}| < 1");
      L0_L0_text_Minv_L1->AddText(Form("%.2f < p_{T}^{1} < %.2f GeV/c", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0_text_Minv_L1->AddText(Form("%.2f < p_{T}^{2} < %.2f GeV/c", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0_text_Minv_L1->AddText(Form("Sig./Bckg. = %.2f", (L_L_US_ratio_integral - L_L_US_LS_ratio_integral)/L_L_US_LS_ratio_integral));
      L0_L0_text_Minv_L1->SetFillColorAlpha(0, 0.01);
      L0_L0_text_Minv_L1->Draw("same");

      TLegend *L0_L0_leg_Minv_L1 = new TLegend(0.2, 0.35, 0.45, 0.54);
      L0_L0_leg_Minv_L1->AddEntry(L01_inv_mass_LL_projection_US[pTbin1][pTbin2], "US(#pi^{-}p)");
      L0_L0_leg_Minv_L1->AddEntry(L01_inv_mass_LL_projection_LS[pTbin1][pTbin2], "Conbinatorial bckg.");
      L0_L0_leg_Minv_L1->SetBorderSize(0);
      L0_L0_leg_Minv_L1->Draw("same");

      L01_inv_mass_LL_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_L/US/L01_inv_mass_LL_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L01_inv_mass_LL_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_L/US/L01_inv_mass_LL_projection_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      TCanvas *L02_inv_mass_LL_projection_US_can = new TCanvas(Form("L02_inv_mass_LL_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L02_inv_mass_LL_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L02_inv_mass_LL_projection_US_can->cd();

      L02_inv_mass_LL_projection_US[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass_US_ME[pTbin1][pTbin2]->ProjectionY(Form("L02_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
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
      L02_inv_mass_LL_projection_LS[pTbin1][pTbin2] = (TH1F*)L0_inv_mass_vs_L0_inv_mass_US_LS_ME[pTbin1][pTbin2]->ProjectionY(Form("L02_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
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
      L0_L0_text_Minv_L2->AddText("Mixed event");
      L0_L0_text_Minv_L2->AddText("Projection to #Lambda^{0}_{2}");
      //L0_L0_text_Minv_L2->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0_L0_text_Minv_L2->AddText("|#it{y}| < 1");
      L0_L0_text_Minv_L2->AddText(Form("%.2f < p_{T}^{1} < %.2f GeV/c", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0_L0_text_Minv_L2->AddText(Form("%.2f < p_{T}^{2} < %.2f GeV/c", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0_L0_text_Minv_L2->AddText(Form("Sig./Bckg. = %.2f", (L_L_US_ratio_integral - L_L_US_LS_ratio_integral)/L_L_US_LS_ratio_integral));
      L0_L0_text_Minv_L2->SetFillColorAlpha(0, 0.01);
      L0_L0_text_Minv_L2->Draw("same");

      L0_L0_leg_Minv_L1->Draw("same");

      L02_inv_mass_LL_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_L/US/L02_inv_mass_LL_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L02_inv_mass_LL_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/L_L/US/L02_inv_mass_LL_projection_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));


      //________________________________________________________________________________________________________

      //Lbar-Lbar

      //US
      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_US_ME_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_ME_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_ME_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME_can->cd();

      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      float Lbar_Lbar_US_scale_integral = L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow_sideBand, binLow_sideBand);
      float Lbar_Lbar_US_ratio_integral = L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow, binHigh);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0bar_L0bar_text_Minv_2D_US_ME = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0bar_L0bar_text_Minv_2D_US_ME->SetTextFont(42);
      //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
      //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_text_Minv_2D_US_ME->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_text_Minv_2D_US_ME->AddText("Mixed event");
      L0bar_L0bar_text_Minv_2D_US_ME->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_text_Minv_2D_US_ME->AddText("US(#pi^{+}#bar{p}) vs. US(#pi^{+}#bar{p})");
      L0bar_L0bar_text_Minv_2D_US_ME->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_Minv_2D_US_ME->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/Lbar_Lbar/US/L0bar_inv_mass_vs_L0bar_inv_mass_US_ME_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass_US_ME_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/Lbar_Lbar/US/L0bar_inv_mass_vs_L0bar_inv_mass_US_ME_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //------------------------------------------------------------------------------------------------

      TCanvas *L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME_can = new TCanvas(Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME_can_%i_%i", pTbin1, pTbin2), Form("L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME_can_%i_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME_can->cd();


      //here just potting

      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->Sumw2();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{-}p}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetXaxis()->SetTitleOffset(2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetXaxis()->SetRangeUser(1.07, 1.2);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetYaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetYaxis()->SetTitleOffset(2.1);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->GetYaxis()->SetRangeUser(1.07,1.2);
      float Lbar_Lbar_US_LS_scale_integral = L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow_sideBand, binLow_sideBand);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->Scale(Lbar_Lbar_US_scale_integral/Lbar_Lbar_US_LS_scale_integral);
      float Lbar_Lbar_US_LS_ratio_integral = L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->Integral(binLow, binHigh, binLow, binHigh);
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->Draw("surf1");

      TPaveText *L0bar_L0bar_text_Minv_2D_US_LS_ME = new TPaveText(0.02, 0.83, 0.32, 0.98, "NDC");
      L0bar_L0bar_text_Minv_2D_US_LS_ME->SetTextFont(42);
      //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
      //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
      //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
      L0bar_L0bar_text_Minv_2D_US_LS_ME->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
      L0bar_L0bar_text_Minv_2D_US_LS_ME->AddText("Mixed event");
      L0bar_L0bar_text_Minv_2D_US_LS_ME->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
      L0bar_L0bar_text_Minv_2D_US_LS_ME->AddText("US(#pi^{+}#bar{p}) vs. LS(#pi^{+}p)");
      L0bar_L0bar_text_Minv_2D_US_LS_ME->AddText("LS(#pi^{+}p) vs. US(#pi^{+}#bar{p})");
      L0bar_L0bar_text_Minv_2D_US_LS_ME->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_Minv_2D_US_LS_ME->Draw("same");

      L0_L0bar_text_Minv_2D_kine->Draw("same");

      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/Lbar_Lbar/US/L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/Lbar_Lbar/US/L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //-------------------------------------------------------------------------

      TCanvas *L0bar1_inv_mass_LbarLbar_projection_US_can = new TCanvas(Form("L0bar1_inv_mass_LbarLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar1_inv_mass_LbarLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar1_inv_mass_LbarLbar_projection_US_can->cd();

      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->ProjectionX(Form("L0bar1_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->SetTitle("M_{inv}^{#pi^{+}#bar{p}}");
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetXaxis()->CenterTitle();
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->SetTitle("Counts");
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->GetYaxis()->CenterTitle();
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetMarkerStyle(20);
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetMarkerSize(2);
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetMarkerColor(kRed);
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->SetLineColor(kRed);
      L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2]->Draw("p e");

      //L0bar1_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0bar1_inv_mass_vs_L0bar1_inv_mass_LS[pTbin1][pTbin2]->ProjectionX(Form("L0bar1_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
      L0bar1_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->ProjectionX(Form("L0bar1_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
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
      L0bar_L0bar_text_Minv_L1->AddText("Mixed event");
      L0bar_L0bar_text_Minv_L1->AddText("Projection to #bar{#Lambda^{0}}_{1}");
      //L0bar_L0bar_text_Minv_L1->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0bar_L0bar_text_Minv_L1->AddText("|#it{y}| < 1");
      L0bar_L0bar_text_Minv_L1->AddText(Form("%.2f < p_{T}^{1} < %.2f GeV/c", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0bar_L0bar_text_Minv_L1->AddText(Form("%.2f < p_{T}^{2} < %.2f GeV/c", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0bar_L0bar_text_Minv_L1->AddText(Form("Sig./Bckg. = %.2f", (Lbar_Lbar_US_ratio_integral - Lbar_Lbar_US_LS_ratio_integral)/Lbar_Lbar_US_LS_ratio_integral));
      L0bar_L0bar_text_Minv_L1->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_Minv_L1->Draw("same");

      TLegend *L0bar_L0bar_leg_Minv_L1 = new TLegend(0.2, 0.35, 0.45, 0.54);
      L0bar_L0bar_leg_Minv_L1->AddEntry(L0bar1_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2], "US(#pi^{+}#bar{p})");
      L0bar_L0bar_leg_Minv_L1->AddEntry(L0bar1_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2], "Conbinatorial bckg.");
      L0bar_L0bar_leg_Minv_L1->SetBorderSize(0);
      L0bar_L0bar_leg_Minv_L1->Draw("same");


      L0bar1_inv_mass_LbarLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/Lbar_Lbar/US/L0bar1_inv_mass_LbarLbar_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar1_inv_mass_LbarLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/Lbar_Lbar/US/L0bar1_inv_mass_LbarLbar_projection_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //---------------------------------------------------------------------

      TCanvas *L0bar2_inv_mass_LbarLbar_projection_US_can = new TCanvas(Form("L0bar2_inv_mass_LbarLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), Form("L0bar2_inv_mass_LbarLbar_projection_US_can_pT1_%i_pT2_%i", pTbin1, pTbin2), 1200, 1000);
      L0bar2_inv_mass_LbarLbar_projection_US_can->cd();

      L0bar2_inv_mass_LbarLbar_projection_US[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass_US_ME[pTbin1][pTbin2]->ProjectionY(Form("L0bar2_inv_mass_projection_US_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
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
      L0bar2_inv_mass_LbarLbar_projection_LS[pTbin1][pTbin2] = (TH1F*)L0bar_inv_mass_vs_L0bar_inv_mass_US_LS_ME[pTbin1][pTbin2]->ProjectionY(Form("L0bar2_inv_mass_projection_LS_can_pT1_%i_pT2_%i.png", pTbin1, pTbin2), binLow, binHigh);
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
      L0bar_L0bar_text_Minv_L2->AddText("Mixed event");
      L0bar_L0bar_text_Minv_L2->AddText("Projection to #bar{#Lambda^{0}}_{2}");
      //L0bar_L0bar_text_Minv_L2->AddText("US(#pi^{-}p) vs. US(#pi^{+}#bar{p})");
      L0bar_L0bar_text_Minv_L2->AddText("|#it{y}| < 1");
      L0bar_L0bar_text_Minv_L2->AddText(Form("%.2f < p_{T}^{1} < %.2f GeV/c", pT_bins_corr[pTbin1], pT_bins_corr[pTbin1+1]));
      L0bar_L0bar_text_Minv_L2->AddText(Form("%.2f < p_{T}^{2} < %.2f GeV/c", pT_bins_corr[pTbin2], pT_bins_corr[pTbin2+1]));
      L0bar_L0bar_text_Minv_L2->AddText(Form("Sig./Bckg. = %.2f", (Lbar_Lbar_US_ratio_integral - Lbar_Lbar_US_LS_ratio_integral)/Lbar_Lbar_US_LS_ratio_integral));
      L0bar_L0bar_text_Minv_L2->SetFillColorAlpha(0, 0.01);
      L0bar_L0bar_text_Minv_L2->Draw("same");

      L0bar_L0bar_leg_Minv_L1->Draw("same");

      L0bar2_inv_mass_LbarLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/Lbar_Lbar/US/L0bar2_inv_mass_LbarLbar_projection_US_pT1_%i_pT2_%i.png", pTbin1, pTbin2));
      L0bar2_inv_mass_LbarLbar_projection_US_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/2D_analysis/Minv_ME/Lbar_Lbar/US/L0bar2_inv_mass_LbarLbar_projection_US_pT1_%i_pT2_%i.pdf", pTbin1, pTbin2));

      //________________________________________________________________________________________________________

    }//end for pT bin 2
  }//end for pT bin 1

  InvMassFile->Close();

  return;
}
