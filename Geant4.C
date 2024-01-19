#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TCanvas.h>
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TColor.h"

void Geant4(){
TFile *rootFile = new TFile("merged28_2e10.root", "READ");
TTree* scattering = (TTree*)rootFile->Get("Scattering");
TH2D* hist2 = new TH2D("", "", 70, -0.55, 0.55, 70, -0.55, 0.55);
double_t x, y, E, z;
scattering->SetBranchAddress("fX", &x);
scattering->SetBranchAddress("fY", &y);

for (Long64_t i = 0; i < scattering->GetEntries(); i++) {
	scattering->GetEntry(i);
	hist2->Fill(x, y);

}

TCanvas *canvas = new TCanvas("canvas", "", 670, 600);
canvas->Print("Geant.pdf[");
gStyle->SetPalette(kCandy);
gStyle->SetCanvasColor(kWhite);
gStyle->SetOptStat(0);
hist2->GetXaxis()->SetNdivisions(409);
hist2->GetYaxis()->SetNdivisions(409);
hist2->GetXaxis()->SetTitle("x (mm)");
hist2->GetXaxis()->SetLabelSize(0.04);
hist2->GetYaxis()->SetLabelSize(0.04);
hist2->GetYaxis()->SetLabelOffset(0.005);
hist2->GetXaxis()->SetLabelOffset(0.005);
hist2->GetXaxis()->SetTitleOffset(1.2);
hist2->GetYaxis()->SetTitleOffset(1);
hist2->GetYaxis()->SetTitle("y (mm)");
hist2->GetZaxis()->SetRangeUser(0, 8);
hist2->GetXaxis()->CenterTitle(true);
hist2->GetYaxis()->CenterTitle(true);
hist2->Draw("COLZ");

canvas->Print("Geant.pdf");
//////////////////////////

TH1F* hist1 = new TH1F("", "", 100, 9940, 10060);
scattering->SetBranchAddress("fZ", &z);

for (Long64_t i = 0; i < scattering->GetEntries(); i++) {
	scattering->GetEntry(i);
	hist1->Fill(z);
}

//n of major divisions to 5 and n of minor divisions to 4
hist1->GetXaxis()->SetNdivisions(505, 4);
hist1->GetXaxis()->SetTitle("z (mm)"); 
hist1->GetYaxis()->SetNdivisions(409);
hist1->GetXaxis()->SetLabelSize(0.04);
hist1->GetYaxis()->SetLabelSize(0.04);
hist1->GetYaxis()->SetLabelOffset(0.005);
hist1->GetXaxis()->SetLabelOffset(0.005);
hist1->GetXaxis()->CenterTitle(true);
hist1->GetYaxis()->CenterTitle(true);
hist1->GetXaxis()->SetTitleOffset(1.2);


hist1->SetLineColor(kAzure + 2);
hist1->SetFillColor(kAzure + 2);
hist1->SetFillStyle(3002);
hist1->SetTitle(" ");

gStyle->SetCanvasColor(kWhite);
gStyle->SetOptStat(0);
hist1->Draw("");
canvas->Print("Geant.pdf");
///////////////////////////////////////////

TH1D * hist4 = new TH1D("", "", 100, 0.01, 0.55);

for (Long64_t i = 0; i < scattering->GetEntries(); i++) {
	scattering->GetEntry(i);
	hist4->Fill(std::sqrt(x * x + y * y));

}


hist4->GetXaxis()->SetNdivisions(505, 5);
hist4->GetYaxis()->SetNdivisions(505, 4);
hist4->GetXaxis()->SetLabelSize(0.04);
hist4->GetYaxis()->SetLabelSize(0.04);
hist4->GetYaxis()->SetLabelOffset(0.005);
hist4->GetXaxis()->SetLabelOffset(0.005);
hist4->GetXaxis()->SetTitle("r (mm)");
hist4->GetXaxis()->SetTitleOffset(1.2);
hist4->GetXaxis()->CenterTitle(true);
hist4->GetYaxis()->CenterTitle(true);
hist4->SetLineColor(kPink - 9);
hist4->SetFillColor(kPink - 9);
hist4->SetFillStyle(3002);
hist4->Draw("hist");

canvas->Print("Geant.pdf");

///////////////////////////////////////////
TTree *hits = (TTree *)rootFile->Get("Hits"); 

auto hist = new TProfile2D("", "", 25, -5.99, 5.9, 25, -5.99, 5.9, 0, 1);

hits->SetBranchAddress("fX", &x);
hits->SetBranchAddress("fY", &y);
hits->SetBranchAddress("fE", &E);

Long64_t nEntries = hits->GetEntries();
for (Long64_t i = 0; i < nEntries; i++) {
    hits->GetEntry(i);
	hist->Fill(x, y, (10 - E) * 1e6); 
   
}


gStyle->SetPalette(kViridis);
TColor::InvertPalette();
hist->GetZaxis()->SetRangeUser(0, 1.);
hist->GetXaxis()->SetNdivisions(409);
hist->GetYaxis()->SetNdivisions(409);
hist->GetXaxis()->SetLabelSize(0.04);
hist->GetYaxis()->SetLabelSize(0.04);
hist->GetYaxis()->SetLabelOffset(0.005);
hist->GetXaxis()->SetLabelOffset(0.005);
hist->GetXaxis()->SetTitle("x (mm)");
hist->GetYaxis()->SetTitle("y (mm)");
hist->SetTitle("#DeltaE (eV)");
gStyle->SetTitleAlign(33);
gStyle->SetTitleX(.99);
hist->GetXaxis()->SetTitleOffset(1.2);
hist->GetXaxis()->CenterTitle(true);
hist->GetYaxis()->CenterTitle(true);
hist->Draw("COLZ");

canvas->Print("Geant.pdf");

//////////////////////////

auto hist3 = new TProfile("", "", 28, 0.01, 6, 0, 1);

for (Long64_t i = 0; i < nEntries; i++) {
	hits->GetEntry(i);
	hist3->Fill(std::sqrt(x*x + y*y), (10 - E) * 1e6); 

}

hist3->GetXaxis()->SetNdivisions(505, 4);
hist3->GetYaxis()->SetNdivisions(505, 4);
hist3->GetXaxis()->SetLabelSize(0.04);
hist3->GetYaxis()->SetLabelSize(0.04);
hist3->GetYaxis()->SetLabelOffset(0.005);
hist3->GetXaxis()->SetLabelOffset(0.005);
hist3->GetXaxis()->SetTitle("r (mm)");
hist3->GetYaxis()->SetTitle("#DeltaE (eV)");
hist3->GetXaxis()->SetTitleOffset(1.2);
hist3->GetXaxis()->CenterTitle(true);
hist3->GetYaxis()->CenterTitle(true);
hist3->SetLineColor(kSpring - 7);
hist3->SetFillColor(kSpring - 7);
hist3->SetFillStyle(3002);
hist3->Draw("hist");

canvas->Print("Geant.pdf");


///////////////////////////
TH1D* hist5 = new TH1D("", "", 35, 0.01, 6);

for (Long64_t i = 0; i < hits->GetEntries(); i++) {
	hits->GetEntry(i);
	hist5->Fill(std::sqrt(x * x + y * y));

}

hist5->GetXaxis()->SetNdivisions(505, 5);
hist5->GetYaxis()->SetNdivisions(505, 4);
hist5->GetXaxis()->SetLabelSize(0.04);
hist5->GetYaxis()->SetLabelSize(0.04);
hist5->GetYaxis()->SetLabelOffset(0.005);
hist5->GetXaxis()->SetLabelOffset(0.005);
hist5->GetXaxis()->SetTitle("r (mm)");
hist5->GetXaxis()->SetTitleOffset(1.2);
hist5->GetXaxis()->CenterTitle(true);
hist5->GetYaxis()->CenterTitle(true);
hist5->SetLineColor(kAzure - 9);
hist5->SetFillColor(kAzure - 9);
hist5->SetFillStyle(3002);
hist5->Draw("hist");

canvas->Print("Geant.pdf");
////////////////////////

canvas->Print("Geant.pdf]");




}