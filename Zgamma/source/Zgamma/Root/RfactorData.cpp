#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <TLorentzVector.h>
#include "RfactorData.h"
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

  ConfigReader config = ConfigReader("/home/katet/Programs/Znunugamma/Zgamma/source/Zgamma/config.cfg");

  double MinCut = config.getDouble("MinCut");
  double MediumCut = config.getDouble("MediumCut");
  double MaxCut = config.getDouble("MaxCut");
  TString PhotonIsolationName = config.getString("Isolation");
  bool EWK = config.getBool("EWK");
  bool Inclusive = config.getBool("Inclusive");
  int CoefInversion = config.getInt("TrackInversion");

  double* DataEventCounting(const char fnameData[]){

    double N_BE_data, N_DF_data, N_E_data, N_F_data;
    int Loose = config.getInt("LoosePrime");
    switch(Loose){
        case 2: Loose = 0x27fc01; break;
        case 3: Loose = 0x25fc01; break;
        case 4: Loose = 0x5fc01; break;
        case 5: Loose = 0x1fc01; break;
    }

    TFile *file_data = new TFile(fnameData, "READ");

    TTree *tree = (TTree*)file_data->Get("output_tree");

    TH1F *hist_BE_data = new TH1F ("hist_BE_data", "hist_BE_data", 50, -1, 49);
    TH1F *hist_DF_data = new TH1F ("hist_DF_data", "hist_DF_data", 50, -1, 49);
    TH1F *hist_E_data = new TH1F ("hist_E_data", "hist_E_data", 50, -1, 49);
    TH1F *hist_F_data = new TH1F ("hist_F_data", "hist_F_data", 50, -1, 49);

    double ph_pt, errBE_data, errDF_data, errE_data, errF_data;
    double ph_iso_pt, ph_iso_et40, ph_z_point, metTST_pt, ph_iso_et20, ph_phi, ph_eta, soft_term_pt;
    UInt_t ph_isem, n_ph, n_mu, n_e_medium, n_jets;
    double jet_lead_phi, jet_sublead_phi, metTST_phi, metTSTsignif, weight;
    double jet_lead_eta, jet_lead_pt, jet_lead_E,jet_sublead_pt, jet_sublead_eta, jet_sublead_E;
    TLorentzVector met, ph, jet, jet2;

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
    tree->SetBranchAddress("jet_sublead_E", &jet_sublead_E);    //Edouble MinCut = 4.45, MediumCut = 10.45, MaxCut = 25.45;

    tree->SetBranchAddress("metTST_pt", &metTST_pt);  //MET p_x
    tree->SetBranchAddress("metTST_phi", &metTST_phi);  //p_y

    tree->SetBranchAddress("ph_iso_et40", &ph_iso_et40);
    tree->SetBranchAddress("ph_iso_et20", &ph_iso_et20);
    tree->SetBranchAddress("ph_iso_pt", &ph_iso_pt);
    tree->SetBranchAddress("n_ph", &n_ph);
    tree->SetBranchAddress("n_mu", &n_mu);
    tree->SetBranchAddress("n_jet", &n_jets);
    tree->SetBranchAddress("n_e_looseBL", &n_e_medium);
    tree->SetBranchAddress("ph_isem", &ph_isem);
    tree->SetBranchAddress("ph_z_point", &ph_z_point);
    tree->SetBranchAddress("metTSTsignif", &metTSTsignif);
    tree->SetBranchAddress("soft_term_pt", &soft_term_pt);

    int N_data = (int)tree->GetEntries();

    double ph_z_point_var = config.getDouble("Zpointing");
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

    for(int i = 0; i < N_data; i++){

      tree->GetEntry(i);
      jet.SetPtEtaPhiE(jet_lead_pt,jet_lead_eta,jet_lead_phi,jet_lead_E);
      jet2.SetPtEtaPhiE(jet_sublead_pt,jet_sublead_eta,jet_sublead_phi,jet_sublead_E);
      met.SetPtEtaPhiM(metTST_pt,0,metTST_phi,0);
      ph.SetPtEtaPhiE(ph_pt,ph_eta,ph_phi,ph_iso_et40);

      if(fabs(ph_z_point)>=ph_z_point_var) continue;
      if(ph_pt <= ph_pt_var) continue;
      //if(n_jets < n_jet_var) continue;
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
        // if(metTSTsignif <= metsingif_var_ewk) continue;
        // if(fabs(met.DeltaPhi(jet)) <= deltaPhiJet_var_ewk ) continue;
        // if(fabs(met.DeltaPhi(jet2)) <= deltaPhiSubJet_var_ewk ) continue;
        // if(fabs(met.DeltaPhi(ph)) <= deltaPhiPh_var_ewk) continue;
        // if(soft_term_pt >= SoftTerm_var) continue;
    }

      if(PhotonIsolationName.Contains("FixedCutTight_Tight")){
         if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_BE_data->Fill(n_ph, 1.0);
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_E_data->Fill(n_ph, 1.0);
         else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_DF_data->Fill(n_ph, 1.0);
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_F_data->Fill(n_ph, 1.0);

      }else if(PhotonIsolationName.Contains("FixedCutTightCaloOnly")){

        if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0) hist_BE_data->Fill(n_ph, 1.0);
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_E_data->Fill(n_ph, 1.0);
        else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_DF_data->Fill(n_ph, 1.0);
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_F_data->Fill(n_ph, 1.0);

       }else if(PhotonIsolationName.Contains("FixedCutLoose")){

         if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_BE_data->Fill(n_ph, 1.0);
         else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_E_data->Fill(n_ph, 1.0);
         else if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_DF_data->Fill(n_ph, 1.0);
         else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_F_data->Fill(n_ph, 1.0);

        }
    }


    N_BE_data =  hist_BE_data->IntegralAndError(-1, 49, errBE_data, "");
    N_DF_data =  hist_DF_data->IntegralAndError(-1, 49, errDF_data, "");
    N_E_data =  hist_E_data->IntegralAndError(-1, 49, errE_data, "");
    N_F_data =  hist_F_data->IntegralAndError(-1, 49, errF_data, "");

    //cout<<"Photons isolation: "<<PhotonIsolationName<<endl;
    cout<<"N_B-E_data: "<<N_BE_data<<" +- "<<errBE_data<<endl;
    cout<<"N_D-F_data: "<<N_DF_data<<" +- "<<errDF_data<<endl;
    cout<<"N_E_data: "<<N_E_data<<" +- "<<errE_data<<endl;
    cout<<"N_F_data: "<<N_F_data<<" +- "<<errF_data<<endl;

    file_data->Close();

    double* DataEventArray = new double[8];
    DataEventArray[0] = N_BE_data;
    DataEventArray[1] = N_DF_data;
    DataEventArray[2] = N_E_data;
    DataEventArray[3] = N_F_data;
    DataEventArray[4] = errBE_data;
    DataEventArray[5] = errDF_data;
    DataEventArray[6] = errE_data;
    DataEventArray[7] = errF_data;
    return DataEventArray;
  }



  void RfactorDataCounting(const char *fnameBkg[], double DataEventArray[8], int longArr, int NFilesBkg){

    int Loose = config.getInt("LoosePrime");
    switch(Loose){
        case 2: Loose = 0x27fc01; break;
        case 3: Loose = 0x25fc01; break;
        case 4: Loose = 0x5fc01; break;
        case 5: Loose = 0x1fc01; break;
    }

    double N_BE_data, N_E_data, N_DF_data, N_F_data;
    double errBE_data, errE_data, errDF_data, errF_data;

    N_BE_data = DataEventArray[0], N_DF_data = DataEventArray[1], N_E_data = DataEventArray[2], N_F_data = DataEventArray[3];
    errBE_data = DataEventArray[4], errDF_data = DataEventArray[5], errE_data = DataEventArray[6], errF_data = DataEventArray[7];

    float sum_BE = 0, sum_DF = 0, sum_E = 0, sum_F = 0, R_sum = 0, del_R_sum = 0;
    float sum_err_BE = 0, sum_err_DF = 0, sum_err_E = 0, sum_err_F = 0;
    float sum_B = 0, sum_D = 0, sum_err_B = 0, sum_err_D = 0;

    for(int i = 0; i<NFilesBkg; i++){

       char ftempname[111]{};
       sprintf(ftempname, "%s", fnameBkg[i] );
       TFile *file = new TFile(ftempname, "READ");
       cout<<ftempname<<endl;

     double sum_of_weights_bk_xAOD, sumw_MC16a = 0, weight, sum = 0, ph_pt, sum_koef = 0, koef;
     double ph_iso_pt, ph_iso_et40, ph_z_point, metTST_pt, ph_iso_et20, ph_phi, soft_term_pt, ph_eta;
     UInt_t ph_isem, n_ph, n_mu, n_e_medium, n_jets;
     double jet_lead_phi, jet_sublead_phi, metTST_phi, metTSTsignif;
     Int_t mc_ph_type;
     double jet_lead_eta, jet_lead_pt, jet_lead_E,jet_sublead_pt, jet_sublead_eta, jet_sublead_E;
     TLorentzVector met, ph, jet, jet2;

     double EtoGam_BE = config.getDouble("EtoGamBE");
     double EtoGam_BE_err = config.getDouble("EtoGamBEerr");
     double EtoGam_E = config.getDouble("EtoGamE");
     double EtoGam_E_err = config.getDouble("EtoGamEerr");
     double EtoGam_DF = config.getDouble("EtoGamDF");
     double EtoGam_DF_err = config.getDouble("EtoGamDFerr");
     double EtoGam_F = config.getDouble("EtoGamF");
     double EtoGam_F_err = config.getDouble("EtoGamFerr");

     TTree *tree_MC_sw = (TTree*)file->Get("output_tree_sw");
     TTree *tree = (TTree*)file->Get("output_tree");
     TTree *tree_norm = (TTree*)file->Get("norm_tree");
     tree_MC_sw->SetBranchAddress("sum_of_weights_bk_xAOD", &sum_of_weights_bk_xAOD);
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

     tree->SetBranchAddress("metTSTsignif", &metTSTsignif);
     tree->SetBranchAddress("ph_iso_et40", &ph_iso_et40);
     tree->SetBranchAddress("ph_iso_et20", &ph_iso_et20);
     tree->SetBranchAddress("ph_iso_pt", &ph_iso_pt);
     tree->SetBranchAddress("weight", &weight);
     tree->SetBranchAddress("n_ph", &n_ph);
     tree->SetBranchAddress("n_mu", &n_mu);
     tree->SetBranchAddress("n_jet", &n_jets);
     tree->SetBranchAddress("metTST_pt", &metTST_pt);
     tree->SetBranchAddress("n_e_looseBL", &n_e_medium);
     tree->SetBranchAddress("ph_isem", &ph_isem);
     tree->SetBranchAddress("ph_z_point", &ph_z_point);
     tree->SetBranchAddress("mc_ph_type", &mc_ph_type);
     tree_norm->SetBranchAddress("koef",&koef);
     tree->SetBranchAddress("soft_term_pt", &soft_term_pt);

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

      for(int i = 0; i < 1; i++){
        tree_norm->GetEntry(i);
        sum_koef = koef;
      }


     TH1F *hist_BE = new TH1F ("hist_BE", "hist_BE", 50, -1, 49);
     TH1F *hist_DF = new TH1F ("hist_DF", "hist_DF", 50, -1, 49);
     TH1F *hist_E = new TH1F ("hist_E", "hist_E", 50, -1, 49);
     TH1F *hist_F = new TH1F ("hist_F", "hist_F", 50, -1, 49);

     Double_t lumi_mc16a = 36214.96;
     Double_t lumi_mc16d = 44307.4;
     Double_t lumi_mc16e = 58450.1;

   for(int i = 0; i < N; i++){

     tree->GetEntry(i);
     jet.SetPtEtaPhiE(jet_lead_pt,jet_lead_eta,jet_lead_phi,jet_lead_E);
     jet2.SetPtEtaPhiE(jet_sublead_pt,jet_sublead_eta,jet_sublead_phi,jet_sublead_E);
     met.SetPtEtaPhiM(metTST_pt,0,metTST_phi,0);
     ph.SetPtEtaPhiE(ph_pt,ph_eta,ph_phi,ph_iso_et40);

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
     double SoftTerm_var = config.getDouble("SoftTerm");

     if(fabs(ph_z_point)>=ph_z_point_var || fabs(weight)>=weight_var) continue;
     if(ph_pt <= ph_pt_var) continue;
     //if(n_jets < n_jet_var) continue;
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
       if(metTSTsignif <= metsingif_var_ewk) continue;
       // if(fabs(met.DeltaPhi(jet)) <= deltaPhiJet_var_ewk ) continue;
       // if(fabs(met.DeltaPhi(jet2)) <= deltaPhiSubJet_var_ewk ) continue;
       // if(fabs(met.DeltaPhi(ph)) <= deltaPhiPh_var_ewk) continue;
       // if(soft_term_pt >= SoftTerm_var) continue;
     }

     TString new_ftempname = TString(ftempname);
     if(PhotonIsolationName.Contains("FixedCutTight_Tight")){
       if(new_ftempname.Contains("MC16a")){
        if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0 && (CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion)) hist_BE->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0 && (CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion)) hist_E->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && (CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion)) hist_DF->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && (CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion)) hist_F->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
      }else if(new_ftempname.Contains("MC16d")){
        if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0 && (CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion)) hist_BE->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0 && (CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion)) hist_E->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && (CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion)) hist_DF->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && (CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion)) hist_F->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
      }else if(new_ftempname.Contains("MC16e")){
        if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0 && (CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion)) hist_BE->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0 && (CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion)) hist_E->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && (CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion)) hist_DF->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && (CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion)) hist_F->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
      }

     }else if(PhotonIsolationName.Contains("FixedCutTightCaloOnly")){

       if(new_ftempname.Contains("MC16a")){
        if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0) hist_BE->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_E->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_DF->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_F->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
      }else if(new_ftempname.Contains("MC16d")){
        if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0) hist_BE->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_E->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_DF->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_F->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
      }else if(new_ftempname.Contains("MC16e")){
        if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && ph_isem == 0) hist_BE->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_E->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) < MediumCut && (ph_iso_et40 - 0.022*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_DF->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_F->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
      }

      }else if(PhotonIsolationName.Contains("FixedCutLoose")){

        if(new_ftempname.Contains("MC16a")){
          if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_BE->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
          else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_E->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
          else if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_DF->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
          else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_F->Fill(ph_phi, lumi_mc16a*sum_koef*weight/(sumw_MC16a));
        }else if(new_ftempname.Contains("MC16d")){
          if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_BE->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
          else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_E->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
          else if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_DF->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
          else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_F->Fill(ph_phi, lumi_mc16d*sum_koef*weight/(sumw_MC16a));
        }else if(new_ftempname.Contains("MC16e")){
          if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_BE->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
          else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_E->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
          else if((ph_iso_et20 - 0.065*ph_pt) < MediumCut && (ph_iso_et20 - 0.065*ph_pt) > MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_DF->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));
          else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_F->Fill(ph_phi, lumi_mc16e*sum_koef*weight/(sumw_MC16a));

        }


       }

   }
       double errBE, errDF, errE, errF, errD;

       double N_BE = hist_BE->IntegralAndError(-1, 50, errBE, "");
       double N_DF = hist_DF->IntegralAndError(-1, 50, errDF, "");
       double N_E = hist_E->IntegralAndError(-1, 50, errE, "");
       double N_F = hist_F->IntegralAndError(-1, 50, errF, "");

       double R;
       R = N_BE*N_F/(N_E*N_DF);

       sum_BE += N_BE;
       sum_DF += N_DF;
       sum_E += N_E;
       sum_F += N_F;

       sum_err_BE += errBE*errBE;
       sum_err_DF += errDF*errDF;
       sum_err_E += errE*errE;
       sum_err_F += errF*errF;


       cout<<"Sum in region B-E = "<<fabs(sum_BE)<<" +- "<<sqrt(sum_err_BE)<<endl;
       cout<<"Sum in region E = "<<fabs(sum_E)<<" +- "<<sqrt(sum_err_E)<<endl;
       cout<<"Sum in region D-F = "<<fabs(sum_DF)<<" +- "<<sqrt(sum_err_DF)<<endl;
       cout<<"Sum in region F = "<<fabs(sum_F)<<" +- "<<sqrt(sum_err_F)<<endl;


       if(i == (NFilesBkg-1)){

         sum_BE = - sum_BE + N_BE_data - EtoGam_BE;
         sum_err_BE = sum_err_BE + errBE_data*errBE_data + EtoGam_BE_err*EtoGam_BE_err;
         sum_E = - sum_E + N_E_data - EtoGam_E;
         sum_err_E = sum_err_E + errE_data*errE_data + EtoGam_E_err*EtoGam_E_err;
         sum_DF = - sum_DF + N_DF_data - EtoGam_DF;
         sum_err_DF = sum_err_DF + errDF_data*errDF_data + EtoGam_DF_err*EtoGam_DF_err;
         sum_F = - sum_F + N_F_data - EtoGam_F;
         sum_err_F = sum_err_F + errF_data*errF_data + EtoGam_F_err*EtoGam_F_err;

       }

       file->Close();

  }

  R_sum = sum_BE*sum_F/(sum_E*sum_DF);
  del_R_sum = sqrt(pow(sqrt(sum_err_BE)*sum_F/(sum_DF*sum_E) , 2) + pow(sqrt(sum_err_F)*sum_BE/(sum_DF*sum_E), 2)
    + pow(sqrt(sum_err_DF)*sum_F*sum_BE/(sum_DF*sum_E*sum_DF), 2) + pow(sqrt(sum_err_E)*sum_F*sum_BE/(sum_DF*sum_E*sum_E) , 2));

    cout<<"----Results----"<<endl;
    cout<<"Photons isolation: "<<PhotonIsolationName<<endl;
    cout<<"Sum in region B-E = "<<fabs(sum_BE)<<" +- "<<sqrt(sum_err_BE)<<endl;
    cout<<"Sum in region E = "<<fabs(sum_E)<<" +- "<<sqrt(sum_err_E)<<endl;
    cout<<"Sum in region D-F = "<<fabs(sum_DF)<<" +- "<<sqrt(sum_err_DF)<<endl;
    cout<<"Sum in region F = "<<fabs(sum_F)<<" +- "<<sqrt(sum_err_F)<<endl;
    cout<<"Summing for R factor = "<<R_sum<<" +- "<<del_R_sum<<endl;

}
