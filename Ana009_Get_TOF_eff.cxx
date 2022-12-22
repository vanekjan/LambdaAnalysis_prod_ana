/*
  TOF matching efficiency
  Macro for production of tof_eff_Dmp_run16.root input for the data-driven fast simulator
*/

#include<iostream>
#include<fstream>
#include<vector>
#include"TH1.h"
#include"TH2.h"
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

using namespace std;

void Ana009_Get_TOF_eff()
{
  //load all files
  TFile *fileData = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/input/TOF_match/output.root", "READ");

  TFile *outFile = new TFile("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/output/TOF_eff/tof_eff_Run17_pp510GeV.root", "RECREATE");
  //TFile *outFile = new TFile("./TOF_eff_output/physics_stream/tof_eff_Dmp_run16_HFT_1sig_20_hist_const_fit.root", "RECREATE");

//------------------LOAD HISTOGRAMS FROM OUTPUT FILE-------------------------------------------
  TList *histoList = (TList*) fileData->Get("picoLambdaAnaMaker");

  TH2F *pi_TPC_TOF = static_cast<TH2F*>(histoList->FindObject("h_piTOF_20"));
  TH2F *p_TPC_TOF = static_cast<TH2F*>(histoList->FindObject("h_pTOF_20"));

  TH1F *piPtTPC = (TH1F*)pi_TPC_TOF->ProjectionX("pa", 1, 2); //pion TPC+TOF tracks
  TH1F *piPtTOF = (TH1F*)pi_TPC_TOF->ProjectionX("pb", 2, 2); //pion TOF tracks

  TH1F *pPtTPC = (TH1F*)p_TPC_TOF->ProjectionX("pc", 1, 2); //kaon TPC+TOF tracks
  TH1F *pPtTOF = (TH1F*)p_TPC_TOF->ProjectionX("pd", 2, 2); //kaon TOF tracksk



  TH2F *pi_TPC_TOF_BEMC_match = static_cast<TH2F*>(histoList->FindObject("h_piTOF_BEMC_match"));
  TH2F *p_TPC_TOF_BEMC_match = static_cast<TH2F*>(histoList->FindObject("h_pTOF_BEMC_match"));

  TH1F *piPtTPC_BEMC_match = (TH1F*)pi_TPC_TOF_BEMC_match->ProjectionX("pa_BEMC", 1, 2); //pion TPC+TOF tracks
  TH1F *piPtTOF_BEMC_match = (TH1F*)pi_TPC_TOF_BEMC_match->ProjectionX("pb_BEMC", 2, 2); //pion TOF tracks

  TH1F *pPtTPC_BEMC_match = (TH1F*)p_TPC_TOF_BEMC_match->ProjectionX("pc_BEMC", 1, 2); //kaon TPC+TOF tracks
  TH1F *pPtTOF_BEMC_match = (TH1F*)p_TPC_TOF_BEMC_match->ProjectionX("pd_BEMC", 2, 2); //kaon TOF tracksk



  //_______________________________________________________________________________________

//-------------------------TOF matching efficiency---------------------------------------------------
  TCanvas *TOF_match_canvas = new TCanvas("TOF_match_canvas", "TOF_match_canvas", 1000, 1000);
  TOF_match_canvas->Divide(1,2);

  piPtTOF->Sumw2();
  piPtTPC->Sumw2();
  pPtTOF->Sumw2();
  pPtTPC->Sumw2();

  piPtTOF->Divide(piPtTPC);
  pPtTOF->Divide(pPtTPC);

  piPtTOF->SetNameTitle("h_pi_tof_eff", "h_pi_tof_eff");
  pPtTOF->SetNameTitle("h_p_tof_eff", "h_p_tof_eff");

  TOF_match_canvas->cd(1);
  piPtTOF->Draw();
  piPtTOF->Write();

  TOF_match_canvas->cd(2);
  pPtTOF->Draw();
  pPtTOF->Write();

  TOF_match_canvas->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/TOF_eff/TOF_match_Run17_pp510GeV.png");
  TOF_match_canvas->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/TOF_eff/TOF_match_Run17_pp510GeV.pdf");
  //__________________________________________________________________________________________

  TCanvas *TOF_match_BEMC_canvas = new TCanvas("TOF_match_BEMC_canvas", "TOF_match_BEMC_canvas", 1000, 1000);
  TOF_match_BEMC_canvas->Divide(1,2);

  piPtTOF_BEMC_match->Sumw2();
  piPtTPC_BEMC_match->Sumw2();
  pPtTOF_BEMC_match->Sumw2();
  pPtTPC_BEMC_match->Sumw2();

  piPtTOF_BEMC_match->Divide(piPtTPC_BEMC_match);
  pPtTOF_BEMC_match->Divide(pPtTPC_BEMC_match);

  piPtTOF_BEMC_match->SetNameTitle("h_pi_tof_eff_BEMC_match", "h_pi_tof_eff_BEMC_match");
  pPtTOF_BEMC_match->SetNameTitle("h_p_tof_eff_BEMC_match", "h_p_tof_eff_BEMC_match");

  TOF_match_BEMC_canvas->cd(1);
  piPtTOF_BEMC_match->Draw();
  piPtTOF_BEMC_match->Write();

  TOF_match_BEMC_canvas->cd(2);
  pPtTOF_BEMC_match->Draw();
  pPtTOF_BEMC_match->Write();

  TOF_match_BEMC_canvas->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/TOF_eff/TOF_match_BEMC_Run17_pp510GeV.png");
  TOF_match_BEMC_canvas->SaveAs("/home/jvanek/C_drive_windows/Work/Analysis/STAR/Production/plots/TOF_eff/TOF_match_BEMC_Run17_pp510GeV.pdf");
  //__________________________________________________________________________________________







  //fitPiTOF->Write();
  //fitKTOF->Write();

  fileData->Close();
  outFile->Close();



}
