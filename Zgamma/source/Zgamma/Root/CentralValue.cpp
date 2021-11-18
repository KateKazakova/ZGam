#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <TLorentzVector.h>
#include "CentralValue.h"
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

  double* DataEventCounterCenterV(const char fnameData[]){

    ConfigReader config = ConfigReader("/home/katet/Programs/Znunugamma/Zgamma/source/Zgamma/config.cfg");

    double MinCut = config.getDouble("MinCut");
    double MediumCut = config.getDouble("MediumCut");
    double MaxCut = config.getDouble("MaxCut");
    TString PhotonIsolationName = config.getString("Isolation");
    TString VetoJets = config.getString("VetoJets");
    bool EWK = config.getBool("EWK");
    bool Inclusive = config.getBool("Inclusive");
    int CoefInversion = config.getInt("TrackInversion");

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

    double EtoGam_A = config.getDouble("EtoGamA");
    double EtoGam_A_err = config.getDouble("EtoGamAerr");
    double EtoGam_B = config.getDouble("EtoGamB");
    double EtoGam_B_err = config.getDouble("EtoGamBerr");
    double EtoGam_C = config.getDouble("EtoGamC");
    double EtoGam_C_err = config.getDouble("EtoGamCerr");
    double EtoGam_D = config.getDouble("EtoGamD");
    double EtoGam_D_err = config.getDouble("EtoGamDerr");

    double errA_data, errB_data, errC_data, errD_data;
    double N_A_data, N_B_data, N_C_data, N_D_data;

    int Loose = config.getInt("LoosePrime");
    switch(Loose){
        case 2: Loose = 0x27fc01; break;
        case 3: Loose = 0x25fc01; break;
        case 4: Loose = 0x5fc01; break;
        case 5: Loose = 0x1fc01; break;
    }

    TFile *file_data = new TFile(fnameData, "READ");
    TTree *tree = (TTree*)file_data->Get("output_tree");

    TH1D *hist_A_data = new TH1D ("hist_A_data", "hist_A_data", 50, 1, 50);
    TH1D *hist_B_data = new TH1D ("hist_B_data", "hist_B_data", 50, 1, 50);
    TH1D *hist_C_data = new TH1D ("hist_C_data", "hist_C_data", 50, 1, 50);
    TH1D *hist_D_data = new TH1D ("hist_D_data", "hist_D_data", 50, 1, 50);

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

    for(int i = 0; i < N_data; i++){

      tree->GetEntry(i);
      jet.SetPtEtaPhiE(jet_lead_pt,jet_lead_eta,jet_lead_phi,jet_lead_E);
      jet2.SetPtEtaPhiE(jet_sublead_pt,jet_sublead_eta,jet_sublead_phi,jet_sublead_E);
      met.SetPtEtaPhiM(metTST_pt,0,metTST_phi,0);
      ph.SetPtEtaPhiE(ph_pt,ph_eta,ph_phi,ph_iso_et40);

      if(fabs(ph_z_point)>=ph_z_point_var) continue;
      if(ph_pt <= ph_pt_var) continue;
      if(VetoJets.Contains("YES")){ if(n_jets < n_jet_var) continue; }
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

      if(PhotonIsolationName.Contains("FixedCutTight_Tight")){
         if((ph_iso_et40 - 0.022*ph_pt) < MinCut && ph_isem == 0 && (ph_iso_pt/ph_pt) < 0.05) hist_A_data->Fill(n_ph, 1.0);
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_B_data->Fill(n_ph, 1.0);
         else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && (ph_iso_pt/ph_pt) < 0.05) hist_C_data->Fill(n_ph, 1.0);
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_D_data->Fill(n_ph, 1.0);

      }else if(PhotonIsolationName.Contains("FixedCutTightCaloOnly")){

        if((ph_iso_et40 - 0.022*ph_pt) < MinCut && ph_isem == 0) hist_A_data->Fill(n_ph, 1.0);
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0) hist_B_data->Fill(n_ph, 1.0);
        else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_C_data->Fill(n_ph, 1.0);
        else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 )) hist_D_data->Fill(n_ph, 1.0);

       }else if(PhotonIsolationName.Contains("FixedCutLoose")){

         if((ph_iso_et40 - 0.022*ph_pt) < MinCut && ph_isem == 0 && (ph_iso_pt/ph_pt) < 0.05) hist_A_data->Fill(n_ph, 1.0);
         else if((ph_iso_et40 - 0.022*ph_pt) > MediumCut && (ph_iso_et40 - 0.022*ph_pt) < MaxCut && ph_isem == 0 && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_B_data->Fill(n_ph, 1.0);
         else if((ph_iso_et40 - 0.022*ph_pt) < MinCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && (ph_iso_pt/ph_pt) < 0.05) hist_C_data->Fill(n_ph, 1.0);
         else if((ph_iso_et20 - 0.065*ph_pt) > MediumCut && (ph_iso_et20 - 0.065*ph_pt) < MaxCut && (ph_isem != 0 && (ph_isem & Loose) == 0 ) && CoefInversion*(ph_iso_pt/ph_pt) < 0.05*CoefInversion) hist_D_data->Fill(n_ph, 1.0);

        }
    }


    N_A_data =  hist_A_data->IntegralAndError(1, hist_A_data->GetNbinsX(), errA_data, "");
    N_B_data =  hist_B_data->IntegralAndError(1, hist_B_data->GetNbinsX(), errB_data, "");
    N_C_data =  hist_C_data->IntegralAndError(1, hist_C_data->GetNbinsX(), errC_data, "");
    N_D_data =  hist_D_data->IntegralAndError(1, hist_D_data->GetNbinsX(), errD_data, "");

    //cout<<"Photons isolation: "<<PhotonIsolationName<<endl;
    cout<<"N_A_data: "<<N_A_data<<" +- "<<errA_data<<endl;
    cout<<"N_B_data: "<<N_B_data<<" +- "<<errB_data<<endl;
    cout<<"N_C_data: "<<N_C_data<<" +- "<<errC_data<<endl;
    cout<<"N_D_data: "<<N_D_data<<" +- "<<errD_data<<endl;

    file_data->Close();

    double* DataEventArrayCentralV = new double[8];
    DataEventArrayCentralV[0] = N_A_data;
    DataEventArrayCentralV[1] = N_B_data;
    DataEventArrayCentralV[2] = N_C_data;
    DataEventArrayCentralV[3] = N_D_data;
    DataEventArrayCentralV[4] = errA_data;
    DataEventArrayCentralV[5] = errB_data;
    DataEventArrayCentralV[6] = errC_data;
    DataEventArrayCentralV[7] = errD_data;
    return DataEventArrayCentralV;
  }



  void CentralValueCounter(const char *fnameBkg[], double DataEventArrayCentralV[8], double LeakageMass[6], int longArr, int NFilesBkg){

   ConfigReader config = ConfigReader("/home/katet/Programs/Znunugamma/Zgamma/source/Zgamma/config.cfg");

    int Loose = config.getInt("LoosePrime");
        switch(Loose){
            case 2: Loose = 0x27fc01; break;
            case 3: Loose = 0x25fc01; break;
            case 4: Loose = 0x5fc01; break;
            case 5: Loose = 0x1fc01; break;
        }

   double errA_data, errB_data, errC_data, errD_data;
   double N_A_data, N_B_data, N_C_data, N_D_data;

   double c_B, c_C, c_D;
   double err_c_B, err_c_C, err_c_D;
   double sum_A = 0, sum_B = 0, sum_D = 0, sum_C = 0, R_sum = 0, del_R_sum = 0;
   double sum_err_A = 0, sum_err_B = 0, sum_err_C = 0, sum_err_D = 0;

   c_B = LeakageMass[0], c_C = LeakageMass[1], c_D = LeakageMass[2];
   err_c_B = LeakageMass[3], err_c_C = LeakageMass[4], err_c_D = LeakageMass[5];

   N_A_data = DataEventArrayCentralV[0], N_B_data = DataEventArrayCentralV[1], N_C_data = DataEventArrayCentralV[2], N_D_data = DataEventArrayCentralV[3];
   errA_data = DataEventArrayCentralV[4], errB_data = DataEventArrayCentralV[5], errC_data = DataEventArrayCentralV[6], errD_data = DataEventArrayCentralV[7];

   cout<<"FROM MAIN FUNCTION"<<endl;
   cout<<"N_A_data: "<<N_A_data<<" +- "<<errA_data<<endl;
   cout<<"N_B_data: "<<N_B_data<<" +- "<<errB_data<<endl;
   cout<<"N_C_data: "<<N_C_data<<" +- "<<errC_data<<endl;
   cout<<"N_D_data: "<<N_D_data<<" +- "<<errD_data<<endl;


   for(int i = 0; i<NFilesBkg; i++){

       char ftempname[140]{};
       sprintf(ftempname, "%s", fnameBkg[i] );
       TFile *file = new TFile(ftempname, "READ");
       cout<<ftempname<<endl;

       double sum_of_weights_bk_xAOD, sumw_MC16a = 0, weight, sum = 0, ph_pt, sum_koef = 0, koef;
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
       sumw_MC16a += sum_of_weights_bk_xAOD;
      }
       for(int i = 0; i < 1; i++){
         tree_norm->GetEntry(i);
         sum_koef = koef;
       }


       TH1D *hist_A = new TH1D ("hist_A", "hist_A", 50, 1, 50);
       TH1D *hist_B = new TH1D ("hist_B", "hist_B", 50, 1, 50);
       TH1D *hist_C = new TH1D ("hist_C", "hist_C", 50, 1, 50);
       TH1D *hist_D = new TH1D ("hist_D", "hist_D", 50, 1, 50);

       double MinCut = config.getDouble("MinCut");
       double MediumCut = config.getDouble("MediumCut");
       double MaxCut = config.getDouble("MaxCut");
       TString PhotonIsolationName = config.getString("Isolation");
       TString VetoJets = config.getString("VetoJets");
       bool EWK = config.getBool("EWK");
       bool Inclusive = config.getBool("Inclusive");
       int CoefInversion = config.getInt("TrackInversion");

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


        double EtoGam_A = config.getDouble("EtoGamA");
        double EtoGam_A_err = config.getDouble("EtoGamAerr");
        double EtoGam_B = config.getDouble("EtoGamB");
        double EtoGam_B_err = config.getDouble("EtoGamBerr");
        double EtoGam_C = config.getDouble("EtoGamC");
        double EtoGam_C_err = config.getDouble("EtoGamCerr");
        double EtoGam_D = config.getDouble("EtoGamD");
        double EtoGam_D_err = config.getDouble("EtoGamDerr");

        double sigma_A = config.getDouble("sigmaA");
        double sigma_B = config.getDouble("sigmaB");
        double sigma_C = config.getDouble("sigmaC");
        double sigma_D = config.getDouble("sigmaD");


       Double_t lumi_mc16a = 36214.96;
       Double_t lumi_mc16d = 44307.4;
       Double_t lumi_mc16e = 58450.1;


      for(int i = 0; i < N; i++){

        tree->GetEntry(i);
        jet.SetPtEtaPhiE(jet_lead_pt,jet_lead_eta,jet_lead_phi,jet_lead_E);
        jet2.SetPtEtaPhiE(jet_sublead_pt,jet_sublead_eta,jet_sublead_phi,jet_sublead_E);
        met.SetPtEtaPhiM(metTST_pt,0,metTST_phi,0);
        ph.SetPtEtaPhiE(ph_pt,ph_eta,ph_phi,ph_iso_et40);

        if(fabs(ph_z_point)>=ph_z_point_var || fabs(weight)>=weight_var) continue;
        if(ph_pt <= ph_pt_var) continue;
        if(VetoJets.Contains("YES")){ if(n_jet < n_jet_var) continue; }
        if(n_ph != n_ph_var || n_mu != n_mu_var || n_e_medium != n_e_medium_var) continue;

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


      cout<<"N_A = "<<N_A<<" +- "<<errA<<endl;
      cout<<"N_B = "<<N_B<<" +- "<<errB<<endl;
      cout<<"N_C = "<<N_C<<" +- "<<errC<<endl;
      cout<<"N_D = "<<N_D<<" +- "<<errD<<endl;

      /// couting sum of events with weights
      sum_A += N_A;
      sum_B += N_B;
      sum_C += N_C;
      sum_D += N_D;

      sum_err_A += errA*errA;
      sum_err_B += errB*errB;
      sum_err_C += errC*errC;
      sum_err_D += errD*errD;

      cout<<"loose'2:"<<endl;
      cout<<"Sum in region A = "<<sum_A<<" +- "<<sqrt(sum_err_A)<<endl;
      cout<<"Sum in region B = "<<sum_B<<" +- "<<sqrt(sum_err_B)<<endl;
      cout<<"Sum in region C = "<<sum_C<<" +- "<<sqrt(sum_err_C)<<endl;
      cout<<"Sum in region D = "<<sum_D<<" +- "<<sqrt(sum_err_D)<<endl;

      cout<<"Leakage parameters:"<<endl;
      cout<<"c_B = "<<c_B<<" +- "<<err_c_B<<endl;
      cout<<"c_C = "<<c_C<<" +- "<<err_c_C<<endl;
      cout<<"c_D = "<<c_D<<" +- "<<err_c_D<<endl;


      if(i == (NFilesBkg -1)){

        sum_A = -sum_A + N_A_data - EtoGam_A + sigma_A;
        sum_err_A = sum_err_A + errA_data*errA_data + EtoGam_A_err*EtoGam_A_err;
        sum_B = -sum_B + N_B_data - EtoGam_B + sigma_B;
        sum_err_B = sum_err_B + errB_data*errB_data + EtoGam_B_err*EtoGam_B_err;
        sum_C = -sum_C + N_C_data - EtoGam_C + sigma_C;
        sum_err_C = sum_err_C + errC_data*errC_data + EtoGam_C_err*EtoGam_C_err;
        sum_D = -sum_D + N_D_data - EtoGam_D + sigma_C;
        sum_err_D = sum_err_D + errD_data*errD_data + EtoGam_D_err*EtoGam_D_err;
     }
        file->Close();

    }

    double N_jet_to_gam, N_jet_to_gam_R, N_jet_to_gam_R_data;
    double b, c, a;
    double b_R, c_R, a_R;
    double b_R_data, c_R_data, a_R_data;
    double RfactorMC = config.getDouble("RfactorMC");
    double RfactorData = config.getDouble("RfactorData");
    a = c_D - c_B*c_C;
    b = sum_D + (c_D*sum_A) - (c_B*sum_C + c_C*sum_B);
    c = (sum_D*sum_A) - (sum_C*sum_B);

    a_R = -RfactorMC * c_B*c_C + c_D;
    b_R = sum_D + (c_D*sum_A) - RfactorMC*(c_B*sum_C + c_C*sum_B);
    c_R = (sum_D*sum_A) - RfactorMC*(sum_C*sum_B);

    a_R_data = -RfactorData * c_B*c_C + c_D;
    b_R_data = sum_D + (c_D*sum_A) - RfactorData*(c_B*sum_C + c_C*sum_B);
    c_R_data = (sum_D*sum_A) - RfactorData*(sum_C*sum_B);

    N_jet_to_gam = (sum_A) - ((b - sqrt(b*b - (4*a*c)))/(2*a));
    N_jet_to_gam_R = (sum_A) - ((b_R - sqrt(b_R*b_R - (4*a_R*c_R)))/(2*a_R));
    N_jet_to_gam_R_data = (sum_A) - ((b_R_data - sqrt(b_R_data*b_R_data - (4*a_R_data*c_R_data)))/(2*a_R_data));

    cout<<"Jet to Gam without R = "<<N_jet_to_gam<<endl;
    cout<<"a = "<<a<<endl;
    cout<<"b = "<<b<<endl;
    cout<<"c = "<<c<<endl;

    cout<<"---------------------------------"<<N_jet_to_gam_R<<endl;
    cout<<"Jet to Gam with R MK = "<<endl;
    cout<<"a = "<<a_R<<endl;
    cout<<"b = "<<b_R<<endl;
    cout<<"c = "<<c_R<<endl;

    cout<<"---------------------------------"<<N_jet_to_gam_R_data<<endl;
    cout<<"Jet to Gam with R data = "<<endl;
    cout<<"a = "<<a_R_data<<endl;
    cout<<"b = "<<b_R_data<<endl;
    cout<<"c = "<<c_R_data<<endl;


    cout<<"----Results----"<<endl;
    cout<<"Sum in region A = "<<(sum_A)<<" +- "<<sqrt(sum_err_A)<<endl;
    cout<<"Sum in region B = "<<(sum_B)<<" +- "<<sqrt(sum_err_B)<<endl;
    cout<<"Sum in region C = "<<(sum_C)<<" +- "<<sqrt(sum_err_C)<<endl;
    cout<<"Sum in region D = "<<(sum_D)<<" +- "<<sqrt(sum_err_D)<<endl;


  }
