#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TLegend.h>
#include <THStack.h>
#include <TGraph.h>
#include <TF1.h>
#include <TProfile.h>
#include <TLatex.h>

#include <iostream>
#include <fstream>
#include <string>
#include <list>

#include "common.cpp"

// TODO Error bars


// Hist Binning Parameters
const int bins_1d = 80;
const int bins_2d = 30;
const int bins_resolution = 15;
const int min_bin = 0;
const int e_max = 30;
const int ge_max = 80;
const int norm_min = -2;
const int norm_max = 2;
const int norm_resolution = 60;

const float slice_energy = 3; // GeV

// Cuts
const double r = 0.5;   // r^2 < dphi^2 + deta^2
const double r2 = r * r;

// Plotting
// const std::string secondPlot("energyScale");
const std::string secondPlot("normalizedEnergyHist");
// const std::string secondPlot("etaEnergyHist");


class jetEnergyData: public jetData {
    public:
        TH2F *truthEnergyHist;
        TH2F *normalizedEnergyHist;
        uint32_t fullBins;
        TProfile *profile;
        TH1D *projection;
        float *energy;
        float *scale;
        float *resolution;
        TGraph *jetScale;
        TGraph *jetResolution;
        TH1D *slice;
};


void plotJetEnergyScale(std::string centralFileList = "", std::string forwardFileList = "", std::string backwardFileList = "") {
    // Initialization, i.e. loading file list and creating histogram
    jetEnergyData jets[NUM_REGIONS];
    if (centralFileList != "") {
        jets[CENTRAL].loaded = true;
        jets[CENTRAL].descriptiveName = std::string("Central");
        std::cout << "loaded " << readFileList(centralFileList, jets[CENTRAL].files) << " files in central region" << std::endl;
        jets[CENTRAL].color = kRed;
        jets[CENTRAL].secondaryColor = kMagenta;
        jets[CENTRAL].marker = 21;
        jets[CENTRAL].markerSize = 1;
        jets[CENTRAL].truthEnergyHist = new TH2F(Form("energy_ratio, %s", jets[CENTRAL].descriptiveName.c_str()), "", bins_2d, min_bin, e_max, bins_2d, min_bin, e_max);
        jets[CENTRAL].normalizedEnergyHist = new TH2F(Form("reco-truth/truth, %s", jets[CENTRAL].descriptiveName.c_str()), "", bins_resolution, min_bin, e_max, norm_resolution, norm_min, norm_max);
    }
    if (forwardFileList != "") {
        jets[FORWARD].loaded = true;
        jets[FORWARD].descriptiveName = std::string("Forward");
        std::cout << "loaded " << readFileList(forwardFileList, jets[FORWARD].files) << " files in forward region" << std::endl;
        jets[FORWARD].color = kBlue;
        jets[FORWARD].secondaryColor = kCyan;
        jets[FORWARD].marker = 22;
        jets[FORWARD].markerSize = 1.5;
        jets[FORWARD].truthEnergyHist = new TH2F(Form("energy_ratio, %s", jets[FORWARD].descriptiveName.c_str()), "", bins_2d, min_bin, e_max, bins_2d, min_bin, e_max);
        jets[FORWARD].normalizedEnergyHist = new TH2F(Form("reco-truth/truth, %s", jets[FORWARD].descriptiveName.c_str()), "", bins_resolution, min_bin, e_max, norm_resolution, norm_min, norm_max);
    }
    if (backwardFileList != "") {
        jets[BACKWARD].loaded = true;
        jets[BACKWARD].descriptiveName = std::string("Backward");
        std::cout << "loaded " << readFileList(backwardFileList, jets[BACKWARD].files) << " files in backward region" << std::endl;
        jets[BACKWARD].color = kGreen;
        jets[BACKWARD].secondaryColor = kYellow + 1;
        jets[BACKWARD].marker = 21;
        jets[BACKWARD].markerSize = 1;
        jets[BACKWARD].truthEnergyHist = new TH2F(Form("energy_ratio, %s", jets[BACKWARD].descriptiveName.c_str()), "", bins_2d, min_bin, e_max, bins_2d, min_bin, e_max);
        jets[BACKWARD].normalizedEnergyHist = new TH2F(Form("reco-truth/truth, %s", jets[BACKWARD].descriptiveName.c_str()), "", bins_resolution, min_bin, e_max, norm_resolution, norm_min, norm_max);
    }

    // Loop over files
    for (uint8_t jetRegion = 0; jetRegion < NUM_REGIONS; jetRegion++) {
        if (!jets[jetRegion].loaded) {
            continue;
        }
        for (std::list<std::string>::iterator iter = jets[jetRegion].files.begin(); iter != jets[jetRegion].files.end(); ++iter) {
            TFile *inFile = TFile::Open((*iter).c_str());
            TTree *jetTree = (TTree*) inFile->Get("ntp_truthjet");
            if (jetTree == nullptr) {
                std::cout << "Could not find jet tree" << std::endl;
            }
            float truthE, recoE;
            float pos[4];
            jetTree->SetBranchAddress("ge", &truthE);
            jetTree->SetBranchAddress("e", &recoE);
            jetTree->SetBranchAddress("geta", &pos[0]);
            jetTree->SetBranchAddress("eta", &pos[1]);
            jetTree->SetBranchAddress("gphi", &pos[2]);
            jetTree->SetBranchAddress("phi", &pos[3]);


            // Create histograms
            for (uint32_t i = 0; i < jetTree->GetEntries(); i++) {
                jetTree->GetEntry(i);
                if (r2 < calculateDistance(pos)) {
                    continue;
                }
                // Filling Histograms
                if (!std::isnan(recoE) && !std::isnan(truthE)) {
                    jets[jetRegion].truthEnergyHist->Fill(truthE, recoE);
                    jets[jetRegion].normalizedEnergyHist->Fill(truthE, (recoE - truthE) / truthE);
                }
                // std::cout << truthE << "\t" << recoE << std::endl;
            
            }
            inFile->Close();
        }
    }

    
    // Calculate energy scale and resolution
    // TProfile *profile = truthEnergyHist->ProfileX();
    for (uint8_t jetRegion = 0; jetRegion < NUM_REGIONS; jetRegion++) {
        if (!jets[jetRegion].loaded) {
            continue;
        }
        jets[jetRegion].profile = jets[jetRegion].normalizedEnergyHist->ProfileX();
        jets[jetRegion].projection = jets[jetRegion].normalizedEnergyHist->ProjectionX();
        jets[jetRegion].profile->BuildOptions(0, 0, "s");
        jets[jetRegion].energy = (float*)malloc(bins_resolution * sizeof(float));
        jets[jetRegion].scale = (float*)malloc(bins_resolution * sizeof(float));
        jets[jetRegion].resolution = (float*)malloc(bins_resolution * sizeof(float));
        uint32_t fullBins = 0;
        for (uint32_t i = 1; i <= bins_resolution; i++) {
            if (jets[jetRegion].projection->GetBinContent(i) == 0) {
                continue;
            }
            jets[jetRegion].energy[fullBins] = jets[jetRegion].profile->GetBinCenter(i);
            // scale[fullBins] = (profile->GetBinContent(i) - profile->GetBinCenter(i)) / profile->GetBinCenter(i);
            jets[jetRegion].scale[fullBins] = jets[jetRegion].profile->GetBinContent(i);
            jets[jetRegion].resolution[fullBins] = jets[jetRegion].profile->GetBinError(i);
            fullBins++;
            // std::cout << projection->GetBinContent(i) << "\t" << energy[i] << "\t" << scale[i] << std::endl;
        }
        jets[jetRegion].fullBins = fullBins;
    }

    // gStyle->SetPadLeftMargin(0.15);


    // Plotting
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadTopMargin(0.12);
    gStyle->SetPadBottomMargin(0.12);
    TCanvas *jetEnergy = new TCanvas("jet_energy", "", 1000, 500);
    jetEnergy->Divide(2, 1);
    
    jetEnergy->cd(1);
    TMultiGraph *jetScale = new TMultiGraph("jet_energy_scale", "Jet Energy Scale");
    TLegend *jetScaleLegend = new TLegend(0.65, 0.8, 0.87, 0.87);
    for (uint8_t jetRegion = 0; jetRegion < NUM_REGIONS; jetRegion++) {
        if (!jets[jetRegion].loaded) {
            continue;
        }
        jets[jetRegion].jetScale = new TGraph(jets[jetRegion].fullBins, jets[jetRegion].energy, jets[jetRegion].scale);
        jets[jetRegion].jetScale->SetMarkerColor(jets[jetRegion].color);
        jets[jetRegion].jetScale->SetMarkerStyle(jets[jetRegion].marker);
        jets[jetRegion].jetScale->SetMarkerSize(jets[jetRegion].markerSize);
        jetScale->Add(jets[jetRegion].jetScale);
        jetScaleLegend->AddEntry(jets[jetRegion].jetScale, Form("%s Jets", jets[jetRegion].descriptiveName.c_str()));
    }
    jetScale->Draw("ALP");
    jetScale->GetXaxis()->SetTitle("Energy");
    jetScale->GetYaxis()->SetTitle("Scale (Mean((reco-truth)/truth))");
    jetScale->SetTitle("Jet Energy Scale");
    jetScaleLegend->Draw();


    
    jetEnergy->cd(2);
    TMultiGraph *jetResolution = new TMultiGraph("jet_energy_resolution", "Jet Energy Resolution");
    TLegend *jetResolutionLegend = new TLegend(0.65, 0.8, 0.87, 0.87);
    for (uint8_t jetRegion = 0; jetRegion < NUM_REGIONS; jetRegion++) {
        if (!jets[jetRegion].loaded) {
            continue;
        }
        jets[jetRegion].jetResolution = new TGraph(jets[jetRegion].fullBins, jets[jetRegion].energy, jets[jetRegion].resolution);
        jets[jetRegion].jetResolution->SetMarkerColor(jets[jetRegion].color);
        jets[jetRegion].jetResolution->SetMarkerStyle(jets[jetRegion].marker);
        jets[jetRegion].jetResolution->SetMarkerSize(jets[jetRegion].markerSize);
        jetResolution->Add(jets[jetRegion].jetResolution);
        jetResolutionLegend->AddEntry(jets[jetRegion].jetResolution, Form("%s Jets", jets[jetRegion].descriptiveName.c_str()));
    }
    jetResolution->Draw("ALP");
    jetResolution->GetXaxis()->SetTitle("Energy");
    jetResolution->GetYaxis()->SetTitle("Resolution (RMS((reco-truth)/truth))"); 
    jetResolution->SetTitle("Jet Energy Resolution");
    jetResolutionLegend->Draw();

    jetEnergy->SaveAs("canvas.png");
    jetEnergy->SaveAs("canvas.c");


    TCanvas *sliceCanvas = new TCanvas("slice", "", 500, 500);
    THStack *sliceStack = new THStack();
    TLegend *sliceLegend = new TLegend(0.65, 0.8, 0.87, 0.87);
    for (uint8_t jetRegion = 0; jetRegion < NUM_REGIONS; jetRegion++) {
        if (!jets[jetRegion].loaded) {
            continue;
        }
        jets[jetRegion].slice = jets[jetRegion].normalizedEnergyHist->ProjectionY(Form("%s slice", jets[jetRegion].descriptiveName.c_str()), 
                                                                                  slice_energy, slice_energy);
        jets[jetRegion].slice->SetLineColor(jets[jetRegion].color);
        sliceStack->Add(jets[jetRegion].slice);
        sliceLegend->AddEntry(jets[jetRegion].slice, Form("%s Jets", jets[jetRegion].descriptiveName.c_str()));
    }
    sliceStack->Draw("nostack");
    sliceStack->SetTitle(Form("%.0f GeV Profile", 2 * slice_energy)); // TODO Make this independent of binning
    sliceStack->GetXaxis()->SetTitle("(reco - truth) / truth");
    sliceStack->GetYaxis()->SetTitle("Counts");
    sliceLegend->Draw();
    sliceCanvas->SaveAs("canvas2.png");
    sliceCanvas->SaveAs("canvas2.c");



    // Some cleanup
    // pft, who needs cleanup
    // delete jetEnergy;
    // delete jetScale;
    // delete truthEnergyHist;
    // delete truthGeHist;
    // delete truthEHist;

    // free(energy);
    // free(scale);
}