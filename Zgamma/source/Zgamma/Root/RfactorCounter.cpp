#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include "ConfigReader.h"

// #include "AtlasUtils.C"
// #include "AtlasLabels.C"
// #include "AtlasStyle.C"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "RfactorCounter.h"
#include <TH2F.h>
#include <TProfile.h>

using namespace std;

 // double RfactorCounter::RfactorMC(){
  double RfactorMC(const char* fnameMC[], int NFiles){
   //SetAtlasStyle();
   ConfigReader config = ConfigReader("../source/Zgamma/config.cfg");
   TCanvas *c0 = new TCanvas("c0", "c0",147,128,700,458);
   c0->SetRightMargin(0.15);
   //gStyle->SetPalette(kRainBow);

   TH2F *hist = new TH2F("hist", "Hist", 30, 0.00001, 0.3, 30, 100, 400);
   auto hprof  = new TProfile("hprof","Profile of pz versus px", 30, 0.00001, 0.3, 100, 400);
   //hist->GetXaxis()->SetTitle("p_{T}^{cone20}/p_{T}, [GeV]");
   hist->GetYaxis()->SetTitle("E_{T}^{miss}, [GeV]");

   bool WithNjets = true;
   //bool WithNjets = false;

   hist->GetXaxis()->SetLabelSize(0.035);
   hist->GetYaxis()->SetLabelSize(0.04);  /// set the size of the font ( <0.15 smaller, >0.15 bigger)

   hist->GetYaxis()->SetTitleOffset(0.9);
   hist->GetXaxis()->SetTitleOffset(1.01);

   // Электрослабый анализ (используются только предотброры для увеличеия статистики):
   // (fabs(ph_z_point)>=250 || fabs(weight)>=100) continue;
   // (mc_ph_type >= 13 && mc_ph_type <= 15) continue;
   // (metTST_pt <= 120) continue; // EWK analysis
   // (ph_pt <= 150) continue;
   // (n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;

   // Инклюзивный анализ:
   // preselections
   // (fabs(ph_z_point)>=250 || fabs(weight)>=100) continue;
   // (metTST_pt <= 130) continue; // EWK analysis
   // (ph_pt <= 150) continue;
   // (n_ph !=1 || n_mu !=0 || n_e_medium != 0) continue;
   // selections
   // (metTSTsignif <= 11) continue;
   // if(fabs(met.DeltaPhi(jet)) <= 0.4 ) continue;
   // if(fabs(met.DeltaPhi(jet2)) <= 0.3 ) continue;
   // if(fabs(met.DeltaPhi(ph)) <= 0.7) continue;

   // if(PhotonIsolationName.Contains("FixedCutTight_Tight") || PhotonIsolationName.Contains("FixedCutTightCaloOnly")){
   //   //hist->GetXaxis()->SetTitle("E_{T}^{cone40} - 0.022p_{T}, [GeV]");
   //   hist->GetXaxis()->SetTitle("p_{T}^{cone20}/p_{T}");
   //   //hist->GetXaxis()->SetTitle("#Delta#phi(p_{T}^{miss}, j_{2}), [GeV]");
   // } else //hist->GetXaxis()->SetTitle("#Delta#phi(p_{T}^{miss}, j_{2}), [GeV]");
   // hist->GetXaxis()->SetTitle("E_{T}^{cone20} - 0.065p_{T}, [GeV]");

   TString PhotonIsolationName = config.getString("Isolation");
   bool EWK = config.getBool("EWK");
   bool Inclusive = config.getBool("Inclusive");
   int CoefInversion = config.getInt("TrackInversion");

   int Loose = config.getInt("LoosePrime");
   switch(Loose){
       case 2: Loose = 0x27fc01; break;
       case 3: Loose = 0x25fc01; break;
       case 4: Loose = 0x5fc01; break;
       case 5: Loose = 0x1fc01; break;
   }

   double R_sum_1;

   double MinCut = config.getDouble("MinCut");
   double MediumCut = config.getDouble("MediumCut");
   double MaxCut = config.getDouble("MaxCut");

   double ArrayR[60]{};
   double ArrayRError[60]{};

   for(int j = 0; j < 1; j++){

     double sum_A = 0, sum_B = 0, sum_C = 0, R_sum = 0, sum_D = 0, del_R_sum = 0;
     double sum_err_A = 0, sum_err_B = 0, sum_err_C = 0, sum_err_D = 0;

   for(int i = 0; i< NFiles; i++){

     char ftempname[120]{};
     sprintf( ftempname, "%s", fnameMC[i] );
     TFile *file = new TFile(ftempname, "READ");
     cout<<ftempname<<endl;

   double sum_of_weights_bk_xAOD, sumw_MC16a = 0, weight, sum = 0, ph_pt, sum_koef = 0, koef;
   double ph_iso_pt, ph_iso_et40, ph_z_point, metTST_pt, ph_iso_et20, metTSTsignif;
   double ph_phi, jet_lead_phi, jet_sublead_phi, metTST_phi, ph_eta, softTerm;
   UInt_t ph_isem, n_ph, n_mu, n_e_medium, n_jet;
   Int_t mc_ph_type;
   double jet_lead_eta, jet_lead_pt, jet_lead_E,jet_sublead_pt, jet_sublead_eta, jet_sublead_E;
   TLorentzVector met, ph, jet, jet2;

   TTree *tree_MC_sw = (TTree*)file->Get("output_tree_sw");
   TTree *tree = (TTree*)file->Get("output_tree");
   TTree *tree_norm = (TTree*)file->Get("norm_tree");
   tree_MC_sw->SetBranchAddress("sum_of_weights_bk_xAOD",&sum_of_weights_bk_xAOD);
   tree->SetBranchAddress("weight",&weight);
   tree->SetBranchAddress("ph_pt",&ph_pt);
   tree->SetBranchAddress("ph_phi",&ph_phi);
   tree->SetBranchAddress("ph_eta",&ph_eta);

   tree->SetBranchAddress("jet_lead_pt", &jet_lead_pt);  //leading jet p_x
   tree->SetBranchAddress("jet_lead_eta", &jet_lead_eta);  //p_y
   tree->SetBranchAddress("jet_lead_phi", &jet_lead_phi);  //p_z
   tree->SetBranchAddress("jet_lead_E", &jet_lead_E);    //E

   tree->SetBranchAddress("jet_sublead_pt", &jet_sublead_pt);  //leading jet p_x
   tree->SetBranchAddress("jet_sublead_eta", &jet_sublead_eta);  //p_y
   tree->SetBranchAddress("jet_sublead_phi", &jet_sublead_phi);  //p_z
   tree->SetBranchAddress("jet_sublead_E", &jet_sublead_E);    //E

   tree->SetBranchAddress("metTST_pt", &metTST_pt);  //MET p_x
   tree->SetBranchAddress("metTST_phi", &metTST_phi);  //p_y

   tree->SetBranchAddress("ph_iso_et40", &ph_iso_et40);
   tree->SetBranchAddress("ph_iso_et20", &ph_iso_et20);
   tree->SetBranchAddress("ph_iso_pt", &ph_iso_pt);
   tree->SetBranchAddress("weight", &weight);
   tree->SetBranchAddress("n_ph", &n_ph);
   tree->SetBranchAddress("n_mu", &n_mu);
   tree->SetBranchAddress("n_jet", &n_jet);
   tree->SetBranchAddress("n_e_looseBL", &n_e_medium);
   tree->SetBranchAddress("ph_isem", &ph_isem);
   tree->SetBranchAddress("ph_z_point", &ph_z_point);
   tree->SetBranchAddress("mc_ph_type", &mc_ph_type);
   tree->SetBranchAddress("metTSTsignif", &metTSTsignif);
   tree->SetBranchAddress("soft_term_pt", &softTerm);
   tree_norm->SetBranchAddress("koef",&koef);

   int entry = (int)tree_MC_sw->GetEntries();
   int N = (int)tree->GetEntries();
   int N_koef = (int)tree_norm->GetEntries();
   for (int i=0; i<entry; i++) {
    tree_MC_sw->GetEntry(i);
    sumw_MC16a += sum_of_weights_bk_xAOD;
   }
   for(int i = 0; i < N; i++){
     tree->GetEntry(i);
     sum += weight;
   }
   for (int i=0; i<1; i++) {
    tree_norm->GetEntry(i);
    sum_koef += koef;
   }


   TH1F *hist_A = new TH1F ("hist_A", "hist_A", 100, 0, 50);
   TH1F *hist_B = new TH1F ("hist_B", "hist_B", 100, 0, 50);
   TH1F *hist_C = new TH1F ("hist_C", "hist_C", 100, 0, 50);
   TH1F *hist_D = new TH1F ("hist_D", "hist_D", 100, 0, 50);


   Double_t lumi_mc16a = 36214.96;
   Double_t lumi_mc16d = 44307.4;
   Double_t lumi_mc16e = 58450.1;

   double ph_z_point_var = config.getDouble("Zpointing");
   double weight_var = config.getDouble("Weight");
   int mc_ph_type_var_min = config.getInt("PhTypeMin");
   int mc_ph_type_var_max = config.getInt("PhTypeMax");
   double ph_pt_var = config.getDouble("PhPt");
   int n_jet_var = config.getInt("NJet");
   int n_ph_var = config.getInt("NPh");
   int n_mu_var = config.getInt("Nmu");
   int n_e_medium_var = config.getInt("Ne");
   int metTST_var_inc = config.getInt("metTSTinc");
   int metTST_var_ewk = config.getInt("metTSTewk");
   int metsingif_var_inc = config.getInt("metsingifinc");
   int metsingif_var_ewk = config.getInt("metsingifewk");
   double deltaPhiPh_var_inc = config.getDouble("deltaPhiPhinc");
   double deltaPhiPh_var_ewk = config.getDouble("deltaPhiPhewk");
   double deltaPhiJet_var_inc = config.getDouble("deltaPhiJetinc");
   double deltaPhiJet_var_ewk = config.getDouble("deltaPhiJetewk");
   double deltaPhiSubJet_var_inc = config.getDouble("deltaPhiSubJetinc");
   double deltaPhiSubJet_var_ewk = config.getDouble("deltaPhiSubJetewk");
   int SoftTerm_var = config.getInt("SoftTerm");


  for(int i = 0; i < N; i++){

    tree->GetEntry(i);
    jet.SetPtEtaPhiE(jet_lead_pt,jet_lead_eta,jet_lead_phi,jet_lead_E);
    jet2.SetPtEtaPhiE(jet_sublead_pt,jet_sublead_eta,jet_sublead_phi,jet_sublead_E);
    met.SetPtEtaPhiM(metTST_pt,0,metTST_phi,0);
    ph.SetPtEtaPhiE(ph_pt,ph_eta,ph_phi,ph_iso_et40);

    if(fabs(ph_z_point)>=ph_z_point_var || fabs(weight)>=weight_var) continue;
    if(ph_pt <= ph_pt_var) continue;
    //if(mc_ph_type >= mc_ph_type_var_min && mc_ph_type <= mc_ph_type_var_max) continue;
    //if(n_jet < n_jet_var) continue;
    if(n_ph != n_ph_var || n_mu != n_mu_var || n_e_medium != n_e_medium_var) continue;

    //Inclusive
    if(Inclusive){
    if(metTST_pt <= metTST_var_inc) continue;
    if(metTSTsignif <= metsingif_var_inc) continue;
    if(fabs(met.DeltaPhi(jet)) <= deltaPhiJet_var_inc ) continue;
    //if(fabs(met.DeltaPhi(jet2)) <= deltaPhiSubJet_var_inc ) continue;
    if(fabs(met.DeltaPhi(ph)) <= deltaPhiPh_var_inc) continue;
  }

    //EWK
    if(EWK){
    // if(metTST_pt <= metTST_var_ewk) continue;
    // if(metTSTsignif <= metsingif_var_ewk) continue;
    // if(fabs(met.DeltaPhi(jet)) <= deltaPhiJet_var_ewk ) continue;
    // if(fabs(met.DeltaPhi(jet2)) <= deltaPhiSubJet_var_ewk ) continue;
    // if(fabs(met.DeltaPhi(ph)) <= deltaPhiPh_var_ewk) continue;
    // if(softTerm >= SoftTerm_var) continue;
  }
    TString new_ftempname = TString(ftempname);
    if(PhotonIsolationName.Contains("FixedCutTight_Tight")){
     if(new_ftempname.Contains("MC16a")){
      if((ph_iso_et40 - 0.022*ph_pt ) < MinCut && ph_isem == 0 && ph_iso_pt/ph_pt < 0.05) hist_A->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_B->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && ph_iso_pt/ph_pt < 0.05) hist_C->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_D->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
    }else if(new_ftempname.Contains("MC16d")){
      if((ph_iso_et40 - 0.022*ph_pt ) < MinCut && ph_isem == 0 && ph_iso_pt/ph_pt < 0.05) hist_A->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_B->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && ph_iso_pt/ph_pt < 0.05) hist_C->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_D->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
    }else if(new_ftempname.Contains("MC16e")){
      if((ph_iso_et40 - 0.022*ph_pt ) < MinCut && ph_isem == 0 && ph_iso_pt/ph_pt < 0.05) hist_A->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_B->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && ph_iso_pt/ph_pt < 0.05) hist_C->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_D->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
    }
  }else if(PhotonIsolationName.Contains("FixedCutTightCaloOnly")){

    if(new_ftempname.Contains("MC16a")){
     if((ph_iso_et40 - 0.022*ph_pt ) < MinCut && ph_isem == 0) hist_A->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
     else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_B->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
     else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_C->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
     else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_D->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
   }else if(new_ftempname.Contains("MC16d")){
     if((ph_iso_et40 - 0.022*ph_pt ) < MinCut && ph_isem == 0) hist_A->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
     else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_B->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
     else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_C->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
     else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_D->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
   }else if(new_ftempname.Contains("MC16e")){
     if((ph_iso_et40 - 0.022*ph_pt ) < MinCut && ph_isem == 0) hist_A->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
     else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_B->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
     else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_C->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
     else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_D->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
   }
     }else if(PhotonIsolationName.Contains("FixedCutLoose")){
      if(new_ftempname.Contains("MC16a")){
       if((ph_iso_et20 - 0.065*ph_pt) < MinCut && ph_isem == 0 && ph_iso_pt/ph_pt < 0.05) hist_A->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_B->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et20 - 0.065*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && ph_iso_pt/ph_pt < 0.05) hist_C->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_D->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
    }else if(new_ftempname.Contains("MC16d")){
       if((ph_iso_et20 - 0.065*ph_pt) < MinCut && ph_isem == 0 && ph_iso_pt/ph_pt < 0.05) hist_A->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_B->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et20 - 0.065*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && ph_iso_pt/ph_pt < 0.05) hist_C->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
       else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_D->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
    }else if(new_ftempname.Contains("MC16e")){
      if((ph_iso_et20 - 0.065*ph_pt) < MinCut && ph_isem == 0 && ph_iso_pt/ph_pt < 0.05) hist_A->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_B->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et20 - 0.065*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && ph_iso_pt/ph_pt < 0.05) hist_C->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
      else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_D->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
    }
      }

  }


   Double_t errA, errB, errC, errD;

   double N_A = hist_A->IntegralAndError(1, hist_A->GetNbinsX(), errA, "");
   double N_B = hist_B->IntegralAndError(1, hist_B->GetNbinsX(), errB, "");
   double N_C = hist_C->IntegralAndError(1, hist_C->GetNbinsX(), errC, "");
   double N_D = hist_D->IntegralAndError(1, hist_D->GetNbinsX(), errD, "");



   double R;
   R = N_A*N_D/(N_C*N_B);

   cout<<"N_A = "<<N_A<<" +- "<<errA<<endl;
   cout<<"N_B = "<<N_B<<" +- "<<errB<<endl;
   cout<<"N_C = "<<N_C<<" +- "<<errC<<endl;
   cout<<"N_D = "<<N_D<<" +- "<<errD<<endl;
   double deltaR;
   deltaR = sqrt(pow(errA*N_D/(N_B*N_C) , 2) + pow(errD*N_A/(N_B*N_C), 2) + pow(errB*N_D*N_A/(N_B*N_C*N_B), 2) + pow(errC*N_D*N_A/(N_B*N_C*N_C) , 2));
   cout<<"R factor = "<<R<<" +- "<<deltaR<<endl;


   /// couting sum of events with weights
   sum_A += N_A;
   sum_B += N_B;
   sum_C += N_C;
   sum_D += N_D;

   sum_err_A += errA*errA;
   sum_err_B += errB*errB;
   sum_err_C += errC*errC;
   sum_err_D += errD*errD;

   cout<<"Photons isolation: "<<PhotonIsolationName<<endl;
   cout<<"Events in region A = "<<sum_A<<" +- "<<sqrt(sum_err_A)<<endl;
   cout<<"Events in region B = "<<sum_B<<" +- "<<sqrt(sum_err_B)<<endl;
   cout<<"Events in region C = "<<sum_C<<" +- "<<sqrt(sum_err_C)<<endl;
   cout<<"Events in region D = "<<sum_D<<" +- "<<sqrt(sum_err_D)<<endl;
   R_sum = sum_A*sum_D/(sum_C*sum_B);
   del_R_sum = sqrt(pow(sqrt(sum_err_A)*sum_D/(sum_B*sum_C) , 2) + pow(sqrt(sum_err_D)*sum_A/(sum_B*sum_C), 2) + pow(sqrt(sum_err_B)*sum_D*sum_A/(sum_B*sum_C*sum_B), 2) + pow(sqrt(sum_err_C)*sum_D*sum_A/(sum_B*sum_C*sum_C) , 2));

   cout<<"R factor = "<<R_sum<<" +- "<<del_R_sum<<endl;
   R_sum_1 = R_sum;
   // if(i == 68){
   // ArrayR[j] = R_sum;
   // } else continue;
   file->Close();
   if(i == 68){
   ArrayR[j] = R_sum;
   ArrayRError[j] = del_R_sum;
 }

 }
     //MaxCut+=2.0;
     if(j == 0) break;
  }
    hist->Draw("COLZ");
    hprof->Draw("SAME");

    // TLegend *leg1 = new TLegend(0.7042607,0.6578261,0.839599,0.8330435);
    //  leg1->SetShadowColor(10);
    //  leg1->SetBorderSize(0);
    //  leg1->SetTextSize(0.05217391);
    //  leg1->SetFillStyle(1002);
    //  leg1->SetFillColor(10);
    //  ATLASLabel(0.19,0.90,"Internal");
    //  leg1->Draw();

 // if(PhotonIsolationName.Contains("FixedCutTight_Tight")){
 //     TLatex *tex = new TLatex(-21.31286,0.399021,"FixedCutTight (MC)");
 //      tex->SetTextColor(kWhite);
 //      tex->SetTextSize(0.054);
 //      tex->SetTextAngle(0.0);
 //      //tex->Draw();
 //  } else if(PhotonIsolationName.Conttex1->SetTextAngle(0.0);
 //        //tex1->Draw();
 //     }
 //
 //     for(int j = 0; j< 60; j++){
 //       if(j == 0) MaxCut = 5.0;
 //       cout<<"MaxCut = "<<MaxCut<<" R factor = "<<ArrayR[j]<<" +- "<<ArrayRError[j]<<endl;
 //       MaxCut+=0.5;
 //     }
 //
 //     for(int j = 0; j< 60; j++){
 //       if(j == 0) MaxCut = 5.0;
 //       cout<<ArrayR[j]<<",";
 //       MaxCut+=0.5;
 //     }
 //
 //     for(int j = 0; j< 60; j++){
 //       if(j == 0) MaxCut = 5.0;
 //       cout<<ArrayRError[j]<<",";
 //       MaxCut+=0.5;
 //     }ains("FixedCutTightCaloOnly")){
 //    TLatex *tex = new TLatex(-21.31286,0.399021,"FixedCutCaloOnly (MC)");
 //     tex->SetTextColor(kWhite);
 //     tex->SetTextSize(0.054);
 //     tex->SetTextAngle(0.0);
 //     //tex->Draw();
 //  } else {
 //    TLatex *tex = new TLatex(-21.31286,0.399021,"FixedCutLoose (MC)");
 //     tex->SetTextColor(kWhite);
 //     tex->SetTextSize(0.054);
 //     tex->SetTextAngle(0.0);
 //     //tex->Draw();
 //  }
 //
 //   if(WithNjets){
 //      TLatex *tex1 = new TLatex(-21.31286,0.2067037,"with cut on N_{jets}");
 //       tex1->SetTextColor(kWhite);
 //       tex1->SetTextSize(0.045);
 //       tex1->SetTextAngle(0.0);
 //       //tex1->Draw();
 //     }else {
 //       TLatex *tex1 = new TLatex(-21.31286,0.2067037,"w/o cut on N_{jets}");
 //        tex1->SetTextColor(kWhite);
 //        tex1->SetTextSize(0.045);
 //        tex1->SetTextAngle(0.0);
 //        //tex1->Draw();
 //     }
 //
 //     for(int j = 0; j< 60; j++){
 //       if(j == 0) MaxCut = 5.0;
 //       cout<<"MaxCut = "<<MaxCut<<" R factor = "<<ArrayR[j]<<" +- "<<ArrayRError[j]<<endl;
 //       MaxCut+=0.5;
 //     }
 //
 //     for(int j = 0; j< 60; j++){
 //       if(j == 0) MaxCut = 5.0;
 //       cout<<ArrayR[j]<<",";
 //       MaxCut+=0.5;
 //     }
 //
 //     for(int j = 0; j< 60; j++){
 //       if(j == 0) MaxCut = 5.0;
 //       cout<<ArrayRError[j]<<",";
 //       MaxCut+=0.5;
 //     }
  return R_sum_1;
 }
