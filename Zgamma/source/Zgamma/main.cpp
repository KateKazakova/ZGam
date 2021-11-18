#include <iostream>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include "RfactorCounter.h"
#include "LeakageParameters.h"
#include "RfactorData.h"
#include "ConfigReader.h"
#include "CentralValue.h"
#include <map>
// #include "AtlasUtils.C"
// #include "AtlasLabels.C"
// #include "AtlasStyle.C"
#include <TCanvas.h>
#include <TH2F.h>
#include <TProfile.h>
using namespace std;

const char *fnameData = "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Data.root";

const char* fnameMC[120] = {
  /*"/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet32.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet33.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet34.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet35.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet36.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet37.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet38.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet39.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet40.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet41.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Multijet42.root",*/
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj364222.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj364223.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366011.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366012.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366013.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366014.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366015.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366016.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366017.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366020.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366021.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366022.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366023.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366024.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366025.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366026.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366029.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366030.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366031.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366032.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366033.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366034.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16a_Zj366035.root",
      /*"/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet31.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet32.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet33.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet34.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet35.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet36.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet37.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet38.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet39.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet40.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet41.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Multijet42.root",*/
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj364222.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj364223.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366011.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366012.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366013.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366014.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366015.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366016.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366017.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366020.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366021.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366022.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366023.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366024.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366025.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366026.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366029.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366030.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366031.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366032.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366033.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366034.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16d_Zj366035.root",
      /*"/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet31.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet32.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet33.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet34.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet35.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet36.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet37.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet38.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet39.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet40.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet41.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Multijet42.root",*/
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj364222.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj364223.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366011.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366012.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366013.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366014.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366015.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366016.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366017.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366020.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366021.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366022.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366023.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366024.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366025.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366026.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366029.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366030.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366031.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366032.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366033.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366034.root",
      "/home/katet/Programs/Znunugamma/Samples/user.akurova.MxAOD_MC16e_Zj366035.root"};

      const char *fnameBkg[111] = {
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361045.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361046.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361047.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361048.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361049.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361050.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361051.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361052.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361053.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361054.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361055.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16a_361056.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361045.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361046.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361047.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361048.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361049.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361050.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361051.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361052.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361053.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361054.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361055.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16d_361056.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361045.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361046.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361047.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361048.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361049.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361050.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361051.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361052.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361053.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361054.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361055.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Gammajet_MC16e_361056.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_ttgamma_MC16a.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_ttgamma_MC16d.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_ttgamma_MC16e.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16a_361273.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16a_361274.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16a_361275.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16d_361273.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16d_361274.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16d_361275.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16e_361273.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16e_361274.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaEWK_MC16e_361275.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16a_364525.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16a_364530.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16a_364535.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16d_364525.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16d_364530.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16d_364535.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16e_364525.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16e_364530.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_WgammaQCD_MC16e_364535.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364184.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364185.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364186.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364187.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364188.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364189.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364190.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364191.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364192.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364193.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364194.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364195.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16a_364196.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364184.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364185.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364186.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364187.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364188.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364189.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364190.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364191.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364192.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364193.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364194.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364195.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16d_364196.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364184.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364185.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364186.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364187.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364188.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364189.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364190.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364191.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364192.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364193.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364194.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364195.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Wtaunu_MC16e_364196.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16a_364504.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16a_364509.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16a_364514.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16d_364504.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16d_364509.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16d_364514.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16e_364504.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16e_364509.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Zllgamma_MC16e_364514.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jQCD_MC16a.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jQCD_MC16d.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jQCD_MC16e.root",
        // "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jEWK_MC16a.root",
        // "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jEWK_MC16d.root",
        // "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jEWK_MC16e.root",
      };


  const char *fnameLeakage[104] = {
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jEWK_MC16a.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jEWK_MC16d.root",
        "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jEWK_MC16e.root",
        // "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jQCD_MC16a.root",
        // "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jQCD_MC16d.root",
        // "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jQCD_MC16e.root",
        // "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jQCD_Pythia_MC16a.root",
        // "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jQCD_Pythia_MC16d.root",
        // "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jQCD_Pythia_MC16e.root",
        // "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jEWK_Herwig_MC16a.root",
        // "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jEWK_Herwig_MC16d.root",
        // "/home/katet/Programs/Znunugamma/ForDataSamples/user.akurova.MxAOD_Znunugamma2jEWK_Herwig_MC16e.root",
      };


 int main(int argc, char **argv){

  unsigned int key;

  cout<<" RFactorMC: 1\n RFactorData: 2\n LeakageParameters: 3\n NJetToGamma:\n";
  cin>>key;
  switch(key){

          case 1: RfactorMC(fnameMC, 69); break;
          case 2: {
                    double *ptr = DataEventCounting(fnameData);
                    RfactorDataCounting(fnameBkg, ptr, 8, 111);
                    break;}
          case 3: {
            cout<<"Do you need the CentralValue?"<<endl;
            int choise;
            cin>>choise;
            double *ptr_leakage = Leakage(fnameLeakage, 3);
            if(choise){
                double *ptr_data = DataEventCounterCenterV(fnameData);
                CentralValueCounter(fnameBkg, ptr_data, ptr_leakage, 8, 108);
                delete[] ptr_data;
                delete[] ptr_leakage;} else break;
           } break;
          default: return 0;
  }

 return 0;


}
