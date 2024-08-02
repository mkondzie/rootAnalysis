//  calibrationSi.c
//  Created by Giulia Colucci on 28/02/24
//  Modified by Maria Kondzielska on 24/04/24

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TString.h"
#include <RtypesCore.h>
#include <TFile.h>
#include <TGraph.h>
#include <TPolyMarker.h>
#include <TStyle.h>
#include <TTree.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <stdio.h>
#include <vector>

template <typename T> void setPlottingOptions(T *plot) {

  plot->SetMarkerSize(0.5);
  plot->SetMarkerStyle(21);
  plot->SetMarkerColor(kBlack);
  plot->SetLineColor(kAzure + 2);
  plot->SetFillColor(kAzure + 2);
  plot->SetLineWidth(1);

  plot->GetYaxis()->SetNdivisions(505, 4);
  plot->GetXaxis()->SetNdivisions(505, 4);

  plot->GetXaxis()->CenterTitle(true);
  plot->GetYaxis()->CenterTitle(true);

  gStyle->SetCanvasColor(kWhite);
  gStyle->SetOptStat(0);
}

template <typename T, typename P> P functionFormula(T *x, P *p) {
  return p[0] + p[1] * (*x) + p[2] * TMath::Power(*x, 2) +
         p[3] * TMath::Power(*x, 2.5) + p[4] * TMath::Power(*x, 3);
}

void readParameters(std::map<int, TF1 *> fittingFunction, int lDetector,
                    bool usePresetParameters = true,
                    TString filename = "kal_bk_det_run_cal_si_02.dat") {
  std::ifstream file(filename.Data());

  if (!file.is_open() || usePresetParameters) {
    fittingFunction[lDetector]->SetParameter(0, -0.025);
    fittingFunction[lDetector]->SetParameter(1, 0.00018);
    fittingFunction[lDetector]->FixParameter(2, 0.0);
    fittingFunction[lDetector]->FixParameter(3, 0.0);
    fittingFunction[lDetector]->FixParameter(4, 0.0);
  }

  else {
    // find the end of the first line and move to the parameters
    std::string firstLine;
    std::getline(file, firstLine);
    size_t pos = firstLine.find('\n');

    if (pos != std::string::npos) {
      file.seekg(pos + 1, std::ios_base::beg);
    }

    TString line;
    TString format("%d %lf %lf %lf %lf %lf %lf");

    while (line.ReadLine(file)) {
      Int_t lInitial;
      Double_t initialA, initialB, initialC, initialD, initialE, initialK;

      if (sscanf(line.Data(), format.Data(), &lInitial, &initialA, &initialB,
                 &initialC, &initialD, &initialE, &initialK) != 7) {
        continue;
      }

      if (lInitial == lDetector) {
        fittingFunction[lDetector]->SetParameter(0, initialA);
        fittingFunction[lDetector]->SetParameter(1, initialB);
        fittingFunction[lDetector]->SetParameter(2, initialC);
        fittingFunction[lDetector]->SetParameter(3, initialD);
        fittingFunction[lDetector]->SetParameter(4, initialE);
      }
    }

    file.close();
  }
}

void calibrationSi() {

  std::vector<int> workingDetectors = {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15};
  const int nSpectra =
      *std::max_element(workingDetectors.begin(), workingDetectors.end()) + 1;

  //--------------------spectra--------------------//
  std::map<int, TH1F *> histSpectra;

  //---------------alpha variables----------------//
  std::map<int, TH1F *> histAlpha;

  std::map<int, TSpectrum *> spectrumAlpha;
  std::map<int, Int_t> nFoundAlpha;
  Double_t *peaksAlphaPositionX[nSpectra];
  TPolyMarker *polyMarkerAlpha;

  //--------------pulser variables---------------//
  std::map<int, TH1F *> histPulser;

  std::map<int, TSpectrum *> spectrumPulser;
  std::map<int, Int_t> nFoundPulser;
  Double_t *peaksPulserPositionX[nSpectra];
  TPolyMarker *polyMarkerPulser;

  std::map<int, TF1 *> fChVeryk;
  std::map<int, TGraph *> graph;

  Double_t voltage[] = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
                        0.4, 0.45, 0.5, 0.55, 0.6, 0.65};

  Double_t A[nSpectra], B[nSpectra], C[nSpectra], D[nSpectra], E[nSpectra], K[nSpectra];

  //////////////////////////////////////////////////
  UShort_t energyChannel;
  UChar_t detectorN;

  TFile *file = new TFile("cal_Si_02.root");
  TTree *cc21 = (TTree *)file->Get("cc21");

  cc21->SetBranchAddress("det", &detectorN);
  cc21->SetBranchAddress("e", &energyChannel);
  //////////////////////////////////////////////////

  UShort_t lowerBound = 250;
  UShort_t lowerBoundAlpha = 310;
  UShort_t upperBoundAlpha = 385;

  for (const auto &l : workingDetectors) {

    Int_t nEntries = (Int_t)cc21->GetEntries();

    //-----------------Create and fill spectra-----------------//
    histSpectra[l] = new TH1F(
        Form("h%i", l), Form("Energy spectrum Silicon %i", l), 4092, 0, 4092);
    histAlpha[l] =
        new TH1F(Form("halpha%i", l),
                 Form("Energy alpha spectrum Silicon %i", l), 4092, 0, 4092);
    histPulser[l] = new TH1F(
        Form("hp%i", l), Form("Pulser spectrum Silicon %i", l), 4092, 0, 4092);

    for (Long_t i = 0; i < nEntries; i++) {

      cc21->GetEntry(i);
      if (detectorN == l) {

        histSpectra[l]->Fill(energyChannel);

        if (energyChannel > lowerBound) {

          if (energyChannel > lowerBoundAlpha &&
              energyChannel < upperBoundAlpha) {

            histAlpha[l]->Fill(energyChannel);

          }

          else if (energyChannel < lowerBoundAlpha ||
                   energyChannel > upperBoundAlpha) {

            histPulser[l]->Fill(energyChannel);
          }
        }
      }
    }

    //---------------Find alpha----------------//
    // Use TSpectrum to find the peak candidates//

    spectrumAlpha[l] = new TSpectrum();

    nFoundAlpha[l] = spectrumAlpha[l]->Search(histAlpha[l], 2, "nodraw", 0.50);
    TList *functionsAlpha = histAlpha[l]->GetListOfFunctions();
    polyMarkerAlpha = (TPolyMarker *)functionsAlpha->FindObject("TPolyMarker");

    std::cout << "----Silicon detector : " << l << std::endl;
    printf("Found %d candidate alpha peaks to fit\n", nFoundAlpha[l]);
    peaksAlphaPositionX[l] = spectrumAlpha[l]->GetPositionX();

    //----------Initialize indices array for sorting----------//
    long long spectrumAlphaIndices[nFoundAlpha[l]];

    for (int k = 0; k < nFoundAlpha[l]; ++k) {
      spectrumAlphaIndices[k] = 0;
    }
    //----------------Sort alpha peaks positions---------------//
    TMath::Sort((long long)nFoundAlpha[l], peaksAlphaPositionX[l],
                spectrumAlphaIndices, kFALSE);
    for (int j = 0; j < nFoundAlpha[l]; ++j) {
      std::cout << peaksAlphaPositionX[l][spectrumAlphaIndices[j]] << std::endl;
    }

    //---------------Find pulser----------------//
    // Use TSpectrum to find the peak candidates//

    spectrumPulser[l] = new TSpectrum();

    nFoundPulser[l] =
        spectrumPulser[l]->Search(histPulser[l], 8, "nodraw", 0.10);
    TList *functionsPulser = histPulser[l]->GetListOfFunctions();
    polyMarkerPulser =
        (TPolyMarker *)functionsPulser->FindObject("TPolyMarker");
    printf("Found %d candidate PULSER peaks to fit\n", nFoundPulser[l]);

    peaksPulserPositionX[l] = spectrumPulser[l]->GetPositionX();
    Double_t pulserPositionsX[nFoundPulser[l]];

    //----------Initialize indices array for sorting----------//
    long long sortedIndices[nFoundPulser[l]];
    for (int k = 0; k < nFoundPulser[l]; ++k) {
      sortedIndices[k] = 0;
    }

    //---------------Sort pulser peaks positions--------------//
    TMath::Sort((long long)nFoundPulser[l], peaksPulserPositionX[l],
                sortedIndices, kFALSE);
    std::cout << "Sorted pulser positions:" << std::endl;
    for (int j = 0; j < nFoundPulser[l]; ++j) {
      pulserPositionsX[j] = peaksPulserPositionX[l][sortedIndices[j]];
      std::cout << pulserPositionsX[j] << std::endl;
    }

    std::cout << "---------------" << std::endl;

    //---------------Fit ~linear formula to the pulser peaks----------------//

    fChVeryk[l] = new TF1(Form("fChVeryk%i", l), functionFormula, 500, 4000, 5);
    readParameters(fChVeryk, l, false);

    graph[l] = new TGraph(nFoundPulser[l], pulserPositionsX, voltage);

    graph[l]->Fit(fChVeryk[l], "RMEN", "", 500, 4000);

    A[l] = fChVeryk[l]->GetParameter(0);
    B[l] = fChVeryk[l]->GetParameter(1);
    C[l] = fChVeryk[l]->GetParameter(2);
    D[l] = fChVeryk[l]->GetParameter(3);
    E[l] = fChVeryk[l]->GetParameter(4);

    double calibrationParameters[5] = {A[l], B[l], C[l], D[l], E[l]};

    Double_t xAlphaPeakChannel;

    xAlphaPeakChannel = peaksAlphaPositionX[l][0];
    double energyAlphaMeV = 5.4395;

    K[l] = energyAlphaMeV /
           functionFormula(&xAlphaPeakChannel, calibrationParameters);
  }
  //-----------Save calibration parameters-----------//
  TString calibrationParametersFilename = "kal_bk_det_run_cal_si_02.dat";

  FILE *calibrationParametersOutfile =
      fopen(std::string(calibrationParametersFilename).c_str(), "w");

  fprintf(calibrationParametersOutfile,
          "%-8s %-15s %-15s %-15s %-15s %-15s %-10s\n", "l", "A[l]", "B[l]",
          "C[l]", "D[l]", "E[l]", "K[l]");

  //----Check reference Ne scattered on Au target----//
  std::vector<double> NeAuChannel = {
      0.0,    0.0,    0.0,    0.0,    2170.0, 2022.0, 2027.0, 0.0,
      2072.0, 2017.0, 2075.0, 2010.0, 2046.0, 2170.0, 0.0,    2095.0};

  TString NeAuTargetFilename = "NeAuTargetcal_si_02.dat";

  FILE *NeAuTargetOutfile = fopen(std::string(NeAuTargetFilename).c_str(), "w");

  fprintf(NeAuTargetOutfile, "%-8s %-15s %-15s\n", "l", "NeAuChannel[l]",
          "NeAuEnergyMeV[l]");
  std::cout << "l\t\t"
            << "NeAuChannel[l]\t"
            << "NeAuEnergyMeV[l]" << std::endl;
  //////////////////////////////////////////////////

  for (const auto &l : workingDetectors) {

    fprintf(calibrationParametersOutfile,
            "%-8d %-15g %-15g %-15g %-15g %-15g %-10g\n", l, A[l], B[l], C[l],
            D[l], E[l], K[l]);

    double calibrationParameters[5] = {A[l], B[l], C[l], D[l], E[l]};

    if (l != 7) {
      Double_t NeAuEnergyMeV =
          K[l] * functionFormula(&NeAuChannel.at(l), calibrationParameters);

      std::cout << l << "\t\t" << NeAuChannel.at(l) << "\t\t" << NeAuEnergyMeV
                << std::endl;

      fprintf(NeAuTargetOutfile, "%-8d %-15g %-15g\n", l, NeAuChannel.at(l),
              NeAuEnergyMeV);
    }
  }

  fclose(calibrationParametersOutfile);
  fclose(NeAuTargetOutfile);

  //-----------Save calibrated energy-----------//
  TFile *outFile = new TFile("cc21CalibrationResult.root", "RECREATE");
  TTree *cc21Calibration = cc21->CloneTree(0);
  Double_t energyMeV;
  cc21Calibration->Branch("eMeV", &energyMeV, "eMeV/d");

  for (auto &l : workingDetectors) {
    double calibrationParameters[5] = {A[l], B[l], C[l], D[l], E[l]};
    for (Long64_t i = 0; i < cc21->GetEntries(); i++) {
      cc21->GetEntry(i);
      if (detectorN == l) {

        energyMeV = K[l] * functionFormula(
                               &energyChannel, calibrationParameters);
      
        cc21Calibration->Fill(); 
      }
    }
  }
  cc21Calibration->Write();
  outFile->Close();
  //-----------------Plot results-----------------//
  TCanvas *canvas = new TCanvas("canvas", "canvas", 10, 10, 1000, 900);
  canvas->Divide(5, 3);

  TCanvas *canvasSpectra =
      new TCanvas("canvasSpectra", "canvasSpectra", 10, 10, 1000, 900);
  canvasSpectra->Divide(5, 3);

  TCanvas *canvasAlpha =
      new TCanvas("canvasAlpha", "canvasAlpha", 10, 10, 1000, 900);
  canvasAlpha->Divide(5, 3);

  TCanvas *canvasPulser =
      new TCanvas("canvasPulser", "canvasPulser", 10, 10, 1000, 900);
  canvasPulser->Divide(5, 3);
  //////////////////////////////////////////////////

  for (const auto &l : workingDetectors) {

    canvas->cd(l);
    graph[l]->GetXaxis()->SetTitle("Channel (V)");
    graph[l]->GetYaxis()->SetTitle("Voltage (V)");
    graph[l]->SetTitle(Form("Voltage(Channel) - Fit Eryk -  Silicon %i", l));
    setPlottingOptions(graph[l]);
    graph[l]->GetXaxis()->SetRange(1, 4092);
    graph[l]->Draw("AP");
    setPlottingOptions(fChVeryk[l]);
    fChVeryk[l]->Draw("same");

    canvasSpectra->cd(l);
    histSpectra[l]->Draw("Hist");
    histSpectra[l]->GetXaxis()->SetTitle("Channels");
    histSpectra[l]->SetTitle(Form("Spectra -  Silicon %i", l));
    setPlottingOptions(histSpectra[l]);
    polyMarkerAlpha->Draw();
    polyMarkerPulser->Draw();

    canvasAlpha->cd(l);
    histAlpha[l]->Draw("Hist");
    histAlpha[l]->GetXaxis()->SetTitle("Channels");
    histAlpha[l]->SetTitle(Form("#alpha Spectra -  Silicon %i", l));
    histAlpha[l]->GetXaxis()->SetRange(310, 380);
    setPlottingOptions(histAlpha[l]);

    canvasPulser->cd(l);
    histPulser[l]->Draw("Hist");
    histPulser[l]->GetXaxis()->SetTitle("Channels");
    histPulser[l]->SetTitle(Form("Spectra Pulser -  Silicon %i", l));
    histPulser[l]->GetXaxis()->SetRange(1, 4092);
    setPlottingOptions(histPulser[l]);
    histPulser[l]->SetLineColor(kAlpine);
  }
 
  for (auto &l : workingDetectors) {

    delete spectrumAlpha[l];
    delete spectrumPulser[l];
  }
}
