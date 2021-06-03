#include <TROOT.h>
#include <TH1F.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <THStack.h>
#include <TLegend.h>

#include <list>
#include <string>
#include <iostream>
#include <fstream>

// Binning
const int num_bins = 50;
const int min_energy = 0;
const int max_energy = 50;

// Cuts
const float r = 0.5;

// Shouldn't need to touch these
const float r2 = r * r;

// Translate file list into list of file paths
// Really should make a personal library of functions...
int readFileList(std::string fileList, std::list<std::string> &list) {
    int numFiles = 0;
    std::ifstream files(fileList);
    std::string filePath;
    
    while (std::getline(files, filePath)) {
        list.push_back(filePath);
        numFiles++;
    }
    files.close();
    return numFiles;
}

// pos = [truthEta, truthPhi, recoEta, recoPhi]
// returns R2 = dEta * dEta + dPhi * dPhi
float calculateDistance(float *pos) {
    for (uint8_t i = 0; i < 4; i++) {
        if (std::isnan(pos[i])) {
            return 9999;
        }
    }
    float dEta, dPhi;
    dEta = pos[0] - pos[2];
    dPhi = pos[1] - pos[3];
    if (dPhi > TMath::Pi()) {
        dPhi -= TMath::TwoPi();
    }
    if (dPhi < -1 * TMath::Pi()) {
        dPhi += TMath::TwoPi();
    }
    return dEta *dEta + dPhi * dPhi;
}

void jetEfficiency(std::string jets = "") {
    // Load files
    std::list<std::string> filePaths;
    std::cout << "read " << readFileList(jets, filePaths) << " file paths" << std::endl;

    // 1D histograms to store number of truth jets and reco jets for each energy bin
    TH1F *truthEnergy = new TH1F("truth_energy", "", num_bins, min_energy, max_energy);
    TH1F *matchedEnergy  = new TH1F("reco_energy",  "", num_bins, min_energy, max_energy);

    // Loop over all the files
    for (std::list<std::string>::iterator iter = filePaths.begin(); iter != filePaths.end(); ++iter) {
        TFile *inFile = TFile::Open((*iter).c_str());       // open root file
        if (inFile == nullptr) {
            std::cerr << "Could not open file " << *iter << std::endl;
            continue;
        }
        TTree *jets = (TTree*) inFile->Get("ntp_truthjet"); // get truthjet tree
        if (jets == nullptr) {
            std::cerr << "Could not file jet tree" << std::endl;
            continue;
        }

        float truthE, recoE;
        float pos[4];

        jets->SetBranchAddress("ge", &truthE);
        jets->SetBranchAddress("e", &recoE);
        jets->SetBranchAddress("geta", &pos[0]);
        jets->SetBranchAddress("gphi", &pos[1]);
        jets->SetBranchAddress("eta", &pos[2]);
        jets->SetBranchAddress("phi", &pos[3]);

        for (uint32_t i = 0; i < jets->GetEntries(); i++) {
            jets->GetEntry(i);
            if (std::isnan(pos[0]) || abs(pos[0]) > 1.5) {
                continue;
            }
            if (std::isnan(truthE)) {
                continue;
            }
            // Do we filter on R for efficiency? Probably
            
            truthEnergy->Fill(truthE);
            if (r2 < calculateDistance(pos)) {
                continue;
            }
            if (std::isnan(truthE) || std::isnan(recoE)) {
                continue;
            }
            matchedEnergy->Fill(truthE);
            // std::cout << truthE << "\t" << recoE << std::endl;
            // std::cout << pos[0] << "\t" << pos[1] << std::endl;
        }
        inFile->Close();
    }

    // Calculate efficiencies
    // efficiency = (num matched) / (num truth)
    float *energy = (float*)malloc(num_bins * sizeof(float));
    float *efficiency = (float*)malloc(num_bins * sizeof(float));
    uint32_t fullBins = 0;
    for (uint32_t i = 1; i < num_bins; i++) {
        if (truthEnergy->GetBinContent(i) == 0 || matchedEnergy->GetBinContent(i) == 0) {
            continue;
        }
        energy[fullBins] = truthEnergy->GetBinCenter(i);
        efficiency[fullBins] = matchedEnergy->GetBinContent(i) / truthEnergy->GetBinContent(i);
        std::cout << matchedEnergy->GetBinContent(i) << "\t" << truthEnergy->GetBinContent(i) << std::endl;
        fullBins++;
    }
    std::cout << "filled " << fullBins << " bins" << std::endl;

    // plotting
    TCanvas *efficiencyCanvas = new TCanvas("jet_efficiency", "", 1000, 500);
    efficiencyCanvas->Divide(2, 1);

    efficiencyCanvas->cd(1);
    THStack *stack = new THStack("jet_energy", "");
    truthEnergy->SetLineColor(kRed);
    matchedEnergy->SetLineColor(kBlue);
    stack->Add(truthEnergy);
    stack->Add(matchedEnergy);
    stack->Draw("nostack");
    stack->SetTitle("Jet Energy");
    stack->GetXaxis()->SetTitle("Jet Energy");
    stack->GetYaxis()->SetTitle("Counts");
    TLegend *legend = new TLegend(0.50, 0.80, 0.9, 0.9);
    legend->AddEntry(truthEnergy, "Truth Jets");
    legend->AddEntry(matchedEnergy, "Matched Jets");
    legend->SetTextSize(0.035);
    legend->Draw();
    std::cout << "integrals: "  << truthEnergy->Integral() << "\t" << matchedEnergy->Integral() << std::endl;
    gPad->SetLogy();
    

    efficiencyCanvas->cd(2);
    // gPad->SetLeftMargin(0.1);
    TGraph *efficiencyGraph = new TGraph(fullBins, energy, efficiency);
    efficiencyGraph->SetTitle("Jet Efficiency");
    efficiencyGraph->GetXaxis()->SetTitle("Jet Energy");
    efficiencyGraph->GetYaxis()->SetTitle("Efficiency");
    efficiencyGraph->Draw("A*");
    // gPad->SetLogy();
    efficiencyCanvas->SaveAs("canvas.png");

    delete truthEnergy;
    delete matchedEnergy;
    delete efficiencyCanvas;
    delete efficiencyGraph;
    delete stack;
    delete legend;

    free(energy);
    free(efficiency);


}