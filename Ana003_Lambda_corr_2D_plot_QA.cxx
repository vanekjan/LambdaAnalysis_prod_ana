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

const int nPtBins = 8;
float const pT_bins[nPtBins+1] = { 0.5, 0.75, 1.,1.5, 2., 2.5, 3., 4., 5.};

//const int nEtaBins = 3;
//float const eta_bins[nEtaBins+1] = { -1, -0.2, 0.2, 1 };
//float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };

void Ana003_Lambda_corr_2D_plot_QA(const int cut_type = 0, const int energy = 510, const int year = 2017)
{


  TFile *InvMassFile; //input file with invariant mass histograms

  TFile *InFile; //input file with dN/dcos(theta*) and QA histoigrams

  if(cut_type == 0)
  {
    if(year == 2012)
    {
      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_orig_full_prod.root", year), "read");
      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_new_prod_open_Vz.root", year), "read");
      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts_012124.root", year), "read");

      InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts_work.root", year), "read");

      //InFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_work.root", year), "read");

      //InFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/cuts_tests/ProdPlane_Lambda_ana_cuts_ana_cuts_R_scan_hist.root", year), "read");

      InFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/20240729_tests_n_sigma_signal/ProdPlane_Lambda_ana_cuts_2_sigma_signal_alt_ME_Delta_pi_third_new_scan_hists_fixed_bckg.root", year), "read"); //analysis cuts
      //InFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/20240729_tests_n_sigma_signal/ProdPlane_Lambda_ana_cuts_3_sigma_signal_alt_ME_Delta_pi_third_new_scan_hists_fixed_bckg.root", year), "read"); //analysis cuts


      //InFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/Delta_R_tests/ProdPlane_Lambda_ana_cuts_Delta_R_0.93.root", year), "read");

      //InFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/tests_Delta_phi_06_2024/ProdPlane_Lambda_ana_cuts_1_SE_in_ME.root", year), "read");
      //InFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/tests_Delta_phi_06_2024/ProdPlane_Lambda_ana_cuts_10_SE_in_ME.root", year), "read");
      //InFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/tests_Delta_phi_06_2024/ProdPlane_Lambda_ana_cuts_100_SE_in_ME.root", year), "read");

    }
    else if(year == 2015)
    {
      InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts_work.root", year), "read");
      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts_new_no_TOF.root", year), "read");
    }
    else if(year == 2016)
    {
      InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts_work.root", year), "read");

      //InFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_work.root", year), "read"); //analysis cuts

      //InFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_2_sigma.root", year), "read"); //analysis cuts
      InFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts_3_sigma.root", year), "read"); //analysis cuts
    }
    else
    {
      //InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts.root", year), "read");
      InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_ana_cuts_work.root", year), "read");
    }
  }

  if(cut_type == 1) InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_tight_topo_cuts.root", year), "read");
  if(cut_type == 2) InvMassFile = new TFile(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/invariant_mass/%i/InvMass_Lambda_tight_pT_cut.root", year), "read");

  if( !(InvMassFile->IsOpen()) )
  {
    cout<<"Unable to open file with invariant mass histograms!"<<endl;
    return;
  }


  if( !(InFile->IsOpen()) )
  {
    cout<<"Unable to open file with QA histograms!"<<endl;
    return;
  }


  //_______________________________________________________________________________________________________________________________________________

  //L and Lbar kinematics
  //delta phi and delta eta distributions

  //same event
  TH1F *L0_L0bar_delta_phi_US_hist = (TH1F*)InFile->Get("L0_L0bar_delta_phi_US_hist");
  TH1F *L0_L0bar_delta_phi_US_LS_hist = (TH1F*)InFile->Get("L0_L0bar_delta_phi_US_LS_hist");

  TH1F *L0_L0bar_delta_eta_US_hist = (TH1F*)InFile->Get("L0_L0bar_delta_eta_US_hist");
  TH1F *L0_L0bar_delta_eta_US_LS_hist = (TH1F*)InFile->Get("L0_L0bar_delta_eta_US_LS_hist");

  TH1F *L0_L0bar_delta_R_US_hist = (TH1F*)InFile->Get("L0_L0bar_delta_R_US_hist");
  TH1F *L0_L0bar_delta_R_US_LS_hist = (TH1F*)InFile->Get("L0_L0bar_delta_R_US_LS_hist");

  TH2F *L0_L0bar_delta_phi_vs_delta_eta_US_hist = (TH2F*)InFile->Get("L0_L0bar_delta_phi_vs_delta_eta_US_hist");
  TH2F *L0_L0bar_delta_phi_vs_delta_eta_US_LS_hist = (TH2F*)InFile->Get("L0_L0bar_delta_phi_vs_delta_eta_US_LS_hist");


  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_hist = (TH2F*)InFile->Get("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist = (TH2F*)InFile->Get("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_hist");

  TH2F *L0_L0bar_delta_phi_vs_nFill_US_ME_hist = (TH2F*)InFile->Get("L0_L0bar_delta_phi_vs_nFill_US_ME_hist");

  TH3F *L0_L0bar_cos_theta_star_vs_delta_phi_vs_nFill_US_ME_hist = (TH3F*)InFile->Get("L0_L0bar_cos_theta_star_vs_delta_phi_vs_nFill_US_ME_hist");


  //------------------------------

  TH1F *L0_L0_delta_phi_US_hist = (TH1F*)InFile->Get("L0_L0_delta_phi_US_hist");
  TH1F *L0_L0_delta_phi_US_LS_hist = (TH1F*)InFile->Get("L0_L0_delta_phi_US_LS_hist");

  //------------------------------

  TH1F *L0bar_L0bar_delta_phi_US_hist = (TH1F*)InFile->Get("L0bar_L0bar_delta_phi_US_hist");
  TH1F *L0bar_L0bar_delta_phi_US_LS_hist = (TH1F*)InFile->Get("L0bar_L0bar_delta_phi_US_LS_hist");

  //mixed event
  TH1F *L0_L0bar_delta_phi_US_ME_hist = (TH1F*)InFile->Get("L0_L0bar_delta_phi_US_ME_hist");
  TH1F *L0_L0bar_delta_phi_US_LS_ME_hist = (TH1F*)InFile->Get("L0_L0bar_delta_phi_US_LS_ME_hist");

  TH1F *L0_L0bar_delta_eta_US_ME_hist = (TH1F*)InFile->Get("L0_L0bar_delta_eta_US_ME_hist");
  TH1F *L0_L0bar_delta_eta_US_LS_ME_hist = (TH1F*)InFile->Get("L0_L0bar_delta_eta_US_LS_ME_hist");


  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist = (TH2F*)InFile->Get("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist = (TH2F*)InFile->Get("L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_LS_ME_hist");

  //---------------------

  TH1F *L0_L0_delta_phi_US_ME_hist = (TH1F*)InFile->Get("L0_L0_delta_phi_US_ME_hist");
  TH1F *L0_L0_delta_phi_US_LS_ME_hist = (TH1F*)InFile->Get("L0_L0_delta_phi_US_LS_ME_hist");

  //----------------------

  TH1F *L0bar_L0bar_delta_phi_US_ME_hist = (TH1F*)InFile->Get("L0bar_L0bar_delta_phi_US_ME_hist");
  TH1F *L0bar_L0bar_delta_phi_US_LS_ME_hist = (TH1F*)InFile->Get("L0bar_L0bar_delta_phi_US_LS_ME_hist");



  //__________________________________________________________________________________________________________________________________________________________________________________________________________

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //L and Lbar kinematics

  //delta phi

  //L-Lbar
  TCanvas *L0_L0bar_delta_phi_ME_can = new TCanvas("L0_L0bar_delta_phi_ME_can", "L0_L0bar_delta_phi_ME_can", 1200, 1000);
  L0_L0bar_delta_phi_ME_can->cd();

  L0_L0bar_delta_phi_US_hist->SetLineColor(1);
  L0_L0bar_delta_phi_US_hist->GetXaxis()->SetTitle("#Delta#phi");
  L0_L0bar_delta_phi_US_hist->GetXaxis()->CenterTitle();
  //L0_L0bar_delta_phi_US_hist->Rebin(2);
  L0_L0bar_delta_phi_US_hist->SetMinimum(0);
  L0_L0bar_delta_phi_US_hist->Draw("hist");

  L0_L0bar_delta_phi_US_ME_hist->SetLineColor(kBlue);
  L0_L0bar_delta_phi_US_ME_hist->Scale(L0_L0bar_delta_phi_US_hist->Integral()/L0_L0bar_delta_phi_US_ME_hist->Integral());
  //L0_L0bar_delta_phi_US_LS_hist->Rebin(2);
  L0_L0bar_delta_phi_US_ME_hist->Draw("hist same");


  L0_L0bar_delta_phi_US_LS_hist->SetLineColor(2);
  //L0_L0bar_delta_phi_US_LS_hist->Rebin(2);
  L0_L0bar_delta_phi_US_LS_hist->Draw("hist same");


  L0_L0bar_delta_phi_US_LS_ME_hist->SetLineColor(6);
  L0_L0bar_delta_phi_US_LS_ME_hist->Scale(L0_L0bar_delta_phi_US_LS_hist->Integral()/L0_L0bar_delta_phi_US_LS_ME_hist->Integral());
  //L0_L0bar_delta_phi_US_LS_ME_hist->Rebin(2);
  L0_L0bar_delta_phi_US_LS_ME_hist->Draw("hist same");

  TPaveText *text_LLbar = new TPaveText(0.4, 0.6, 0.8, 0.8, "NDC");
  //cent_text->AddText("STAR preliminary");
  //((TText*)cent_text->GetListOfLines()->Last())->SetTextColor(2);
  text_LLbar->SetTextFont(43);
  text_LLbar->SetTextSize(33);
  text_LLbar->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  text_LLbar->AddText("Minimum bias");
  text_LLbar->AddText("#Lambda#bar{#Lambda}");
  text_LLbar->SetFillColorAlpha(0, 0.01);
  text_LLbar->Draw("same");

  L0_L0bar_delta_phi_ME_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0_L0bar_delta_phi_ME.png");



  TCanvas *L0_L0bar_delta_phi_separate_can = new TCanvas("L0_L0bar_delta_phi_separate_can", "L0_L0bar_delta_phi_separate_can", 1200, 1000);
  L0_L0bar_delta_phi_separate_can->cd();

  L0_L0bar_delta_phi_US_hist->SetLineColor(1);
  L0_L0bar_delta_phi_US_hist->GetXaxis()->SetTitle("#Delta#phi");
  L0_L0bar_delta_phi_US_hist->GetXaxis()->CenterTitle();
  //L0_L0bar_delta_phi_US_hist->Rebin(2);
  L0_L0bar_delta_phi_US_hist->SetMinimum(0);
  L0_L0bar_delta_phi_US_hist->Draw("hist");

  L0_L0bar_delta_phi_US_LS_hist->SetLineColor(2);
  //L0_L0bar_delta_phi_US_LS_hist->Rebin(2);
  L0_L0bar_delta_phi_US_LS_hist->Draw("hist same");

  text_LLbar->Draw("same");

  L0_L0bar_delta_phi_separate_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0_L0bar_delta_phi_separate.png");


  TCanvas *L0_L0bar_delta_phi_can = new TCanvas("L0_L0bar_delta_phi_can", "L0_L0bar_delta_phi_can", 1200, 1000);
  L0_L0bar_delta_phi_can->cd();

  L0_L0bar_delta_phi_US_hist->Add(L0_L0bar_delta_phi_US_LS_hist, -1);
  L0_L0bar_delta_phi_US_hist->SetMinimum(0);
  L0_L0bar_delta_phi_US_hist->Draw("hist");

  text_LLbar->Draw("same");

  L0_L0bar_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0_L0bar_delta_phi.png");


  TCanvas *L0_L0bar_delta_phi_ME_subtract_can = new TCanvas("L0_L0bar_delta_phi_ME_subtract_can", "L0_L0bar_delta_phi_ME_subtract_can", 1200, 1000);
  L0_L0bar_delta_phi_ME_subtract_can->cd();

  L0_L0bar_delta_phi_US_hist->Draw("hist");

  L0_L0bar_delta_phi_US_ME_hist->Add(L0_L0bar_delta_phi_US_LS_ME_hist, -1);
  L0_L0bar_delta_phi_US_ME_hist->SetMinimum(0);
  L0_L0bar_delta_phi_US_ME_hist->Draw("hist same");

  text_LLbar->Draw("same");

  L0_L0bar_delta_phi_ME_subtract_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0_L0bar_delta_phi_ME_subtract.png");

  //---------------------------------------------------------------

  TCanvas *L0_L0bar_delta_eta_ME_can = new TCanvas("L0_L0bar_delta_eta_ME_can", "L0_L0bar_delta_eta_ME_can", 1200, 1000);
  L0_L0bar_delta_eta_ME_can->cd();

  L0_L0bar_delta_eta_US_hist->SetLineColor(1);
  L0_L0bar_delta_eta_US_hist->GetXaxis()->SetTitle("#Deltay");
  L0_L0bar_delta_eta_US_hist->GetXaxis()->CenterTitle();
  //L0_L0bar_delta_eta_US_hist->Rebin(2);
  L0_L0bar_delta_eta_US_hist->SetMinimum(0);
  L0_L0bar_delta_eta_US_hist->Draw("hist");

  L0_L0bar_delta_eta_US_ME_hist->SetLineColor(kBlue);
  L0_L0bar_delta_eta_US_ME_hist->Scale(L0_L0bar_delta_eta_US_hist->Integral()/L0_L0bar_delta_eta_US_ME_hist->Integral());
  //L0_L0bar_delta_eta_US_LS_hist->Rebin(2);
  L0_L0bar_delta_eta_US_ME_hist->Draw("hist same");

  text_LLbar->Draw("same");

  L0_L0bar_delta_eta_ME_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0_L0bar_delta_eta_ME.png");



  TCanvas *L0_L0bar_delta_eta_separate_can = new TCanvas("L0_L0bar_delta_eta_separate_can", "L0_L0bar_delta_eta_separate_can", 1200, 1000);
  L0_L0bar_delta_eta_separate_can->cd();

  L0_L0bar_delta_eta_US_hist->SetLineColor(1);
  L0_L0bar_delta_eta_US_hist->GetXaxis()->SetTitle("#Deltay");
  L0_L0bar_delta_eta_US_hist->GetXaxis()->CenterTitle();
  //L0_L0bar_delta_eta_US_hist->Rebin(2);
  L0_L0bar_delta_eta_US_hist->SetMinimum(0);
  L0_L0bar_delta_eta_US_hist->Draw("hist");

  L0_L0bar_delta_eta_US_LS_hist->SetLineColor(2);
  //L0_L0bar_delta_eta_US_LS_hist->Rebin(2);
  L0_L0bar_delta_eta_US_LS_hist->Draw("hist same");

  text_LLbar->Draw("same");

  L0_L0bar_delta_eta_separate_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0_L0bar_delta_eta_separate.png");


  TCanvas *L0_L0bar_delta_eta_can = new TCanvas("L0_L0bar_delta_eta_can", "L0_L0bar_delta_eta_can", 1200, 1000);
  L0_L0bar_delta_eta_can->cd();

  L0_L0bar_delta_eta_US_hist->Add(L0_L0bar_delta_eta_US_LS_hist, -1);
  L0_L0bar_delta_eta_US_hist->SetMinimum(0);
  L0_L0bar_delta_eta_US_hist->Draw("hist");

  text_LLbar->Draw("same");

  L0_L0bar_delta_eta_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0_L0bar_delta_eta.png");

  //--------------------------------------------------------------

  TCanvas *L0_L0bar_delta_R_separate_can = new TCanvas("L0_L0bar_delta_R_separate_can", "L0_L0bar_delta_R_separate_can", 1200, 1000);
  L0_L0bar_delta_R_separate_can->cd();

  L0_L0bar_delta_R_US_hist->SetLineColor(1);
  L0_L0bar_delta_R_US_hist->GetXaxis()->SetTitle("#DeltaR");
  L0_L0bar_delta_R_US_hist->GetXaxis()->CenterTitle();
  //L0_L0bar_delta_R_US_hist->Rebin(2);
  L0_L0bar_delta_R_US_hist->SetMinimum(0);
  L0_L0bar_delta_R_US_hist->Draw("hist");

  L0_L0bar_delta_R_US_LS_hist->SetLineColor(2);
  //L0_L0bar_delta_R_US_LS_hist->Rebin(2);
  L0_L0bar_delta_R_US_LS_hist->Draw("hist same");

  L0_L0bar_delta_R_separate_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0_L0bar_delta_R_separate.png");


  TCanvas *L0_L0bar_delta_R_can = new TCanvas("L0_L0bar_delta_R_can", "L0_L0bar_delta_R_can", 1200, 1000);
  L0_L0bar_delta_R_can->cd();

  L0_L0bar_delta_R_US_hist->Add(L0_L0bar_delta_R_US_LS_hist, -1);
  L0_L0bar_delta_R_US_hist->SetMinimum(0);
  L0_L0bar_delta_R_US_hist->Draw("hist");

  L0_L0bar_delta_R_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0_L0bar_delta_R.png");

  //--------------------------------------------------------------

  TCanvas *L0_L0bar_delta_phi_vs_delta_eta_can = new TCanvas("L0_L0bar_delta_phi_vs_delta_eta_can", "L0_L0bar_delta_phi_vs_delta_eta_can", 1200, 1000);
  L0_L0bar_delta_phi_vs_delta_eta_can->cd();

  L0_L0bar_delta_phi_vs_delta_eta_US_hist->Add(L0_L0bar_delta_phi_vs_delta_eta_US_LS_hist, -1);
  L0_L0bar_delta_phi_vs_delta_eta_US_hist->GetXaxis()->SetTitle("#Delta#phi");
  L0_L0bar_delta_phi_vs_delta_eta_US_hist->GetXaxis()->CenterTitle();
  L0_L0bar_delta_phi_vs_delta_eta_US_hist->GetYaxis()->SetTitle("#Delta#eta");
  L0_L0bar_delta_phi_vs_delta_eta_US_hist->GetYaxis()->CenterTitle();
  L0_L0bar_delta_phi_vs_delta_eta_US_hist->Draw("colz");

  L0_L0bar_delta_phi_vs_delta_eta_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0_L0bar_delta_phi_vs_delta_eta.png");

  //-------------------------------------------------------------------------------------------------------------------------------------------

  //L-L
  TCanvas *L0_L0_delta_phi_separate_can = new TCanvas("L0_L0_delta_phi_separate_can", "L0_L0_delta_phi_separate_can", 1200, 1000);
  L0_L0_delta_phi_separate_can->cd();

  L0_L0_delta_phi_US_hist->SetLineColor(1);
  L0_L0_delta_phi_US_hist->GetXaxis()->SetTitle("#Delta#phi");
  L0_L0_delta_phi_US_hist->GetXaxis()->CenterTitle();
  //L0_L0_delta_phi_US_hist->Rebin(2);
  L0_L0_delta_phi_US_hist->SetMinimum(0);
  L0_L0_delta_phi_US_hist->Draw("hist");

  L0_L0_delta_phi_US_LS_hist->SetLineColor(2);
  //L0_L0_delta_phi_US_LS_hist->Rebin(2);
  L0_L0_delta_phi_US_LS_hist->Draw("hist same");

  L0_L0_delta_phi_separate_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0_L0_delta_phi_separate.png");


  TCanvas *L0_L0_delta_phi_can = new TCanvas("L0_L0_delta_phi_can", "L0_L0_delta_phi_can", 1200, 1000);
  L0_L0_delta_phi_can->cd();

  L0_L0_delta_phi_US_hist->Add(L0_L0_delta_phi_US_LS_hist, -1);
  L0_L0_delta_phi_US_hist->SetMinimum(0);
  L0_L0_delta_phi_US_hist->Draw("hist");

  L0_L0_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0_L0_delta_phi.png");

  //-------------------------------------------------------------------------------------------------------------------------------------------

  //Lbar-Lbar
  TCanvas *L0bar_L0bar_delta_phi_separate_can = new TCanvas("L0bar_L0bar_delta_phi_separate_can", "L0bar_L0bar_delta_phi_separate_can", 1200, 1000);
  L0bar_L0bar_delta_phi_separate_can->cd();

  L0bar_L0bar_delta_phi_US_hist->SetLineColor(1);
  L0bar_L0bar_delta_phi_US_hist->GetXaxis()->SetTitle("#Delta#phi");
  L0bar_L0bar_delta_phi_US_hist->GetXaxis()->CenterTitle();
  //L0bar_L0bar_delta_phi_US_hist->Rebin(2);
  L0bar_L0bar_delta_phi_US_hist->SetMinimum(0);
  L0bar_L0bar_delta_phi_US_hist->Draw("hist");

  L0bar_L0bar_delta_phi_US_LS_hist->SetLineColor(2);
  //L0bar_L0bar_delta_phi_US_LS_hist->Rebin(2);
  L0bar_L0bar_delta_phi_US_LS_hist->Draw("hist same");

  L0bar_L0bar_delta_phi_separate_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0bar_L0bar_delta_phi_separate.png");


  TCanvas *L0bar_L0bar_delta_phi_can = new TCanvas("L0bar_L0bar_delta_phi_can", "L0bar_L0bar_delta_phi_can", 1200, 1000);
  L0bar_L0bar_delta_phi_can->cd();

  L0bar_L0bar_delta_phi_US_hist->Add(L0bar_L0bar_delta_phi_US_LS_hist, -1);
  L0bar_L0bar_delta_phi_US_hist->SetMinimum(0);
  L0bar_L0bar_delta_phi_US_hist->Draw("hist");

  L0bar_L0bar_delta_phi_can->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/L_kine/L0bar_L0bar_delta_phi.png");


  //--------------------------------------------------------------------------

  //projections to Delta phi in individual dN/d cos(theta*) bins
  for(unsigned int bin = 0; bin < 10; bin++)
  {
    TH1D* Delta_phi_proj = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_hist->ProjectionY(Form("Proj_%i", bin), bin+1, bin+1);

    TH1D* Delta_phi_ME_proj = L0_L0bar_cos_theta_star_vs_delta_phi_for_corr_US_ME_hist->ProjectionY(Form("Proj_ME_%i", bin), bin+1, bin+1);

    TCanvas *Delta_phi_proj_can = new TCanvas(Form("Delta_phi_proj_can_%i", bin), Form("Delta_phi_proj_can_%i", bin), 1200, 1000);
    Delta_phi_proj_can->cd();

    Delta_phi_proj->GetXaxis()->SetTitle("#Delta#phi");
    Delta_phi_proj->GetXaxis()->CenterTitle();
    Delta_phi_proj->SetMinimum(0);
    Delta_phi_proj->SetLineColor(1);
    Delta_phi_proj->Draw("hist");

    Delta_phi_ME_proj->SetLineColor(kRed);
    Delta_phi_ME_proj->Scale(Delta_phi_proj->Integral()/Delta_phi_ME_proj->Integral());
    Delta_phi_ME_proj->Draw("hist same");

    TPaveText *text = new TPaveText(0.5, 0.7, 0.8, 0.9, "NDC");
    text->AddText(Form("%.1f < dN/dcos(#theta*) < %.1f", -1+bin*0.2, -0.8+bin*0.2));
    text->SetFillColorAlpha(0, 0.01);
    text->Draw("same");

    Delta_phi_proj_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/Delta_phi/L_Lbar/L0_L0bar_Delta_phi_%i.png", bin));


  }

  //--------------------------------------------------------------------------

  //Delta phi dstributions for different N_match ranges

  for(unsigned int Nfill_bin = 0; Nfill_bin < 5; Nfill_bin++)
  {
    TH1D *L0_L0bar_delta_phi_proj_hist = (TH1D*)L0_L0bar_delta_phi_vs_nFill_US_ME_hist->ProjectionX(Form("proj_nFill_%i", Nfill_bin), Nfill_bin*20+1, Nfill_bin*20+20 );


    TCanvas *Delta_phi_proj_can = new TCanvas(Form("Delta_phi_proj_nFill_can_%i", Nfill_bin), Form("Delta_phi_proj_nFill_can_%i", Nfill_bin), 1200, 1000);
    Delta_phi_proj_can->cd();


    L0_L0bar_delta_phi_proj_hist->SetLineColor(1);
    L0_L0bar_delta_phi_proj_hist->SetMinimum(0);
    L0_L0bar_delta_phi_proj_hist->Draw("hist");


    Delta_phi_proj_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/Delta_phi_nFill_proj/L_Lbar/L0_L0bar_Delta_phi_nFill_bin_%i.png", Nfill_bin));


    TH1D *L0_L0bar_cos_theta_star_proj_hist = (TH1D*)L0_L0bar_cos_theta_star_vs_delta_phi_vs_nFill_US_ME_hist->ProjectionX(Form("proj_cos_t_%i", Nfill_bin), 1, 5, Nfill_bin*2+1, Nfill_bin*2+2);

    TCanvas *cos_theta_star_can = new TCanvas(Form("cos_theta_star_can_nFill_can_%i", Nfill_bin), Form("cos_theta_star_can_nFill_can_%i", Nfill_bin), 1200, 1000);
    cos_theta_star_can->cd();

    L0_L0bar_cos_theta_star_proj_hist->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0bar_cos_theta_star_proj_hist->GetXaxis()->CenterTitle();
    //L0_L0bar_cos_theta_star_proj_hist->GetXaxis()->SetTextSizePixels(30);
    //L0_L0bar_cos_theta_star_proj_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0bar_cos_theta_star_proj_hist->GetYaxis()->SetTitle("[d#it{N}/d cos(#theta*)]_{meas}");
    L0_L0bar_cos_theta_star_proj_hist->GetYaxis()->CenterTitle();
    //L0_L0bar_cos_theta_star_proj_hist->GetYaxis()->SetTextSizePixels(30);
    L0_L0bar_cos_theta_star_proj_hist->GetYaxis()->SetMaxDigits(3);
    L0_L0bar_cos_theta_star_proj_hist->SetMarkerSize(1.5);
    L0_L0bar_cos_theta_star_proj_hist->SetMarkerStyle(20);
    L0_L0bar_cos_theta_star_proj_hist->SetMarkerColor(kRed);
    L0_L0bar_cos_theta_star_proj_hist->SetLineColor(kRed);
    L0_L0bar_cos_theta_star_proj_hist->SetMinimum(0);
    L0_L0bar_cos_theta_star_proj_hist->Draw("p e");

    cos_theta_star_can->SaveAs(Form("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/Lambda/QA/cos_theta_star_nFill_proj/L_Lbar/L0_L0bar_cos_theta_star_proj_bin_%i.png", Nfill_bin));

  }

  //--------------------------------------------------------------------------

  //dN/dcos(theta*) distributions for different N_match ranges




  //_______________________________________________________________________________________________________________________________________________


  InvMassFile->Close();

  return;
}
