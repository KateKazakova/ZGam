#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <TLorentzVector.h>
#include "LeakageParameters.h"
#include "ConfigReader.h"

// #include "AtlasUtils.C"
// #include "AtlasLabels.C"
// #include "AtlasStyle.C"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include <TH2F.h>
#include <TProfile.h>

using namespace std;

   // for signal files we use seletion on mc_ph_type [13, 16]
   // for background we don't use this seletion

 double* Leakage(const char* fnameLeakage[], int NFiles){
   //SetAtlasStyle();
   ConfigReader config = ConfigReader("../source/Zgamma/config.cfg");
   double sum_A = 0, sum_B = 0, sum_C = 0, sum_D = 0;
   double sum_BE = 0, sum_E = 0, sum_DF = 0, sum_F = 0;
   double sum_err_A = 0, sum_err_B = 0, sum_err_C = 0, sum_err_D = 0;

   bool EWK = config.getBool("EWK");
   bool Inclusive = config.getBool("Inclusive");
   int CoefInversion = config.getInt("TrackInversion");

   double* LeakageMass = new double[6];
   double c_A, c_B, c_C, c_D;
   double err_c_A, err_c_B, err_c_C, err_c_D;
   double c_BE, c_E, c_DF, c_F;
   double Leakage;

   // double MinCut = 2.45, MediumCut = 2.45 + 2.0, MaxCut = 25.45;
   double MinCut = config.getDouble("MinCut");
   double MediumCut = config.getDouble("MediumCut");
   double MaxCut = config.getDouble("MaxCut");
   int Loose = config.getInt("LoosePrime");
   switch(Loose){
       case 2: Loose = 0x27fc01; break;
       case 3: Loose = 0x25fc01; break;
       case 4: Loose = 0x5fc01; break;
       case 5: Loose = 0x1fc01; break;
   }

   for(int i = 0; i< NFiles; i++){

     char ftempname[104]{};
     sprintf( ftempname, "%s", fnameLeakage[i] );
     TFile *file = new TFile(ftempname, "READ");
     cout<<ftempname<<endl;

   double sum_of_error_A = 0, sum_of_error_B = 0, sum_of_error_C = 0, sum_of_error_D = 0;
   double sum_of_weights_bk_xAOD, sumw_MC16 = 0, weight, sum = 0, ph_pt, sum_koef = 0, koef;
   double ph_iso_pt, ph_iso_et40, ph_z_point, metTST_pt, ph_iso_et20, ph_phi, ph_eta;
   UInt_t ph_isem, n_ph, n_mu, n_e_medium, n_jet;
   Int_t mc_ph_type;
   double jet_lead_phi, jet_sublead_phi, metTST_phi, metTSTsignif, soft_term_pt;
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
   tree->SetBranchAddress("n_jet", &n_jet);
   tree->SetBranchAddress("n_mu", &n_mu);
   tree->SetBranchAddress("metTST_pt", &metTST_pt);
   tree->SetBranchAddress("n_e_looseBL", &n_e_medium);
   tree->SetBranchAddress("ph_isem", &ph_isem);
   tree->SetBranchAddress("ph_z_point", &ph_z_point);
   tree->SetBranchAddress("mc_ph_type", &mc_ph_type);
   tree_norm->SetBranchAddress("koef",&koef);
   tree->SetBranchAddress("metTSTsignif", &metTSTsignif);
   tree->SetBranchAddress("soft_term_pt", &soft_term_pt);

   int entry = (int)tree_MC_sw->GetEntries();
   int N = (int)tree->GetEntries();
   int N_koef = (int)tree_norm->GetEntries();
   for (int i=0; i<entry; i++) {
    tree_MC_sw->GetEntry(i);
    sumw_MC16 += sum_of_weights_bk_xAOD;
   }

    for(int i = 0; i < 1; i++){
      tree_norm->GetEntry(i);
      sum_koef = koef;
    }

   TH1D *hist_A = new TH1D ("hist_A", "hist_A", 12, -3, 3);
   TH1D *hist_B = new TH1D ("hist_B", "hist_B", 12, -3, 3);
   TH1D *hist_C = new TH1D ("hist_C", "hist_C", 12, -3, 3);
   TH1D *hist_D = new TH1D ("hist_D", "hist_D", 12, -3, 3);
   TH1D *hist_BE = new TH1D ("hist_BE", "hist_BE", 12, -3, 3);
   TH1D *hist_DF = new TH1D ("hist_DF", "hist_DF", 12, -3, 3);
   TH1D *hist_E = new TH1D ("hist_E", "hist_E", 12, -3, 3);
   TH1D *hist_F = new TH1D ("hist_F", "hist_F", 12, -3, 3);

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
     //
     if(fabs(ph_z_point)>=ph_z_point_var || fabs(weight)>=weight_var) continue;
     if(mc_ph_type < mc_ph_type_var_min || mc_ph_type > mc_ph_type_var_max) continue;
     if(ph_pt <= ph_pt_var) continue;
     if(n_jet < n_jet_var) continue;
     if(n_ph !=n_ph_var || n_mu !=n_mu_var || n_e_medium != n_e_medium_var) continue;

     //Inclusive
     if(Inclusive){
     if(metTST_pt <= metTST_var_inc) continue;
     if(metTSTsignif <= metsingif_var_inc) continue;
     if(fabs(met.DeltaPhi(jet)) <= deltaPhiJet_var_inc ) continue;
     if(fabs(met.DeltaPhi(jet2)) <= deltaPhiSubJet_var_inc ) continue;
     if(fabs(met.DeltaPhi(ph)) <= deltaPhiPh_var_inc) continue;
   }

     //EWK
     if(EWK){
     if(metTST_pt <= metTST_var_ewk) continue;
     if(metTSTsignif <= metsingif_var_ewk) continue;
     if(fabs(met.DeltaPhi(jet)) <= deltaPhiJet_var_ewk ) continue;
     if(fabs(met.DeltaPhi(jet2)) <= deltaPhiSubJet_var_ewk ) continue;
     if(fabs(met.DeltaPhi(ph)) <= deltaPhiPh_var_ewk) continue;
     if(soft_term_pt >= SoftTerm_var) continue;
   }


     TString new_ftempname = TString(ftempname);
     if(new_ftempname.Contains("MC16a")){
             if ((ph_iso_et40 - 0.022*ph_pt) < MinCut && ph_isem == 0 && ph_iso_pt/ph_pt < 0.05) hist_A->Fill(ph_eta, lumi_mc16a*sum_koef*weight/(sumw_MC16));
             else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_B->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16));
             else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && ph_iso_pt/ph_pt < 0.05) hist_C->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16));
             else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_D->Fill(n_ph, lumi_mc16a*sum_koef*weight/(sumw_MC16));
     }
    else if(new_ftempname.Contains("MC16d")){
            if ((ph_iso_et40 - 0.022*ph_pt) < MinCut && ph_isem == 0 && ph_iso_pt/ph_pt < 0.05) hist_A->Fill(ph_eta, lumi_mc16d*sum_koef*weight/(sumw_MC16));
            else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_B->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16));
            else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && ph_iso_pt/ph_pt < 0.05) hist_C->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16));
            else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_D->Fill(n_ph, lumi_mc16d*sum_koef*weight/(sumw_MC16));
     }
    else if(new_ftempname.Contains("MC16e")){
            if ((ph_iso_et40 - 0.022*ph_pt) < MinCut && ph_isem == 0 && ph_iso_pt/ph_pt < 0.05) hist_A->Fill(ph_eta, lumi_mc16e*sum_koef*weight/(sumw_MC16));
            else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_B->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16));
            else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && ph_iso_pt/ph_pt < 0.05) hist_C->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16));
            else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_D->Fill(n_ph, lumi_mc16e*sum_koef*weight/(sumw_MC16));
    }

  }


     Double_t errA, errB, errC, errD;

     double N_A = hist_A->IntegralAndError(1, hist_A->GetNbinsX(), errA, "");
     double N_B = hist_B->IntegralAndError(1, hist_B->GetNbinsX(), errB, "");
     double N_C = hist_C->IntegralAndError(1, hist_C->GetNbinsX(), errC, "");
     double N_D = hist_D->IntegralAndError(1, hist_D->GetNbinsX(), errD, "");


     // cout<<"N_A = "<<N_A<<" +- "<<errA<<endl;
     // cout<<"N_B = "<<N_B<<" +- "<<errB<<endl;
     // cout<<"N_C = "<<N_C<<" +- "<<errC<<endl;
     // cout<<"N_D = "<<N_D<<" +- "<<errD<<endl;


     /// couting sum of events with weights
     sum_A += N_A;
     sum_B += N_B;
     sum_C += N_C;
     sum_D += N_D;

     sum_err_A += errA*errA;
     sum_err_B += errB*errB;
     sum_err_C += errC*errC;
     sum_err_D += errD*errD;

     // cout<<"loose'2:"<<endl;
     // cout<<"-----Result-----"<<endl;
     // cout<<"Sum in region A = "<<sum_A<<" +- "<<sqrt(sum_err_A)<<endl;
     // cout<<"Sum in region B = "<<sum_B<<" +- "<<sqrt(sum_err_B)<<endl;
     // cout<<"Sum in region C = "<<sum_C<<" +- "<<sqrt(sum_err_C)<<endl;
     // cout<<"Sum in region D = "<<sum_D<<" +- "<<sqrt(sum_err_D)<<endl;

     c_B = sum_B/sum_A;
     c_C = sum_C/sum_A;
     c_D = sum_D/sum_A;

     err_c_B = sqrt(pow(sqrt(sum_err_A)*sum_B/(sum_A*sum_A) , 2) + pow(sqrt(sum_err_B)/(sum_A), 2));
     err_c_C = sqrt(pow(sqrt(sum_err_A)*sum_C/(sum_A*sum_A) , 2) + pow(sqrt(sum_err_C)/(sum_A), 2));
     err_c_D = sqrt(pow(sqrt(sum_err_A)*sum_D/(sum_A*sum_A) , 2) + pow(sqrt(sum_err_D)/(sum_A), 2));


     // cout<<"Events: "<<endl;
     // cout<<"N_A = "<<sum_A<<endl;
     // cout<<"N_B = "<<sum_B<<endl;
     // cout<<"N_C = "<<sum_C<<endl;
     // cout<<"N_D = "<<sum_D<<endl;

     // cout<<"Leakage parameters:"<<endl;
     // cout<<"c_B = "<<c_B<<" +- "<<err_c_B<<endl;
     // cout<<"c_C = "<<c_C<<" +- "<<err_c_C<<endl;
     // cout<<"c_D = "<<c_D<<" +- "<<err_c_D<<endl;
     LeakageMass[0] = c_B, LeakageMass[1] = c_C, LeakageMass[2] = c_D;
     LeakageMass[3] = err_c_B, LeakageMass[4] = err_c_C, LeakageMass[5] = err_c_D;

   file->Close();
  }
  cout<<"Leakage parameters:"<<endl;
  cout<<"c_B = "<<c_B<<" +- "<<err_c_B<<endl;
  cout<<"c_C = "<<c_C<<" +- "<<err_c_C<<endl;
  cout<<"c_D = "<<c_D<<" +- "<<err_c_D<<endl;
  return LeakageMass;

 }

  double Sum(double LeakageMass[3], int longArr){
    double sum = 0;
    for(int i = 0; i < longArr; i++){
      sum+= LeakageMass[i];
    }
    return sum;
  }
