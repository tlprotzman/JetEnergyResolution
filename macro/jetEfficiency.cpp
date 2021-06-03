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
#include <TMultiGraph.h>

#include <list>
#include <string>
#include <iostream>
#include <fstream>

#include "common.cpp"

// Binning
const int num_bins = 50;
const int min_energy = 0;
const int max_energy = 50;

// Cuts
const float r = 0.5;

// Shouldn't need to touch these
const float r2 = r * r;


// Jet Regions
const int CENTRAL = 0;
const int FORWARD = 1;
const int BACKWARD = 2;

const int NUM_REGIONS = 3;

class jetEfficiencyData: public jetData {
    public:
        uint32_t fullBins;
        TH1F *truthEnergy;
        TH1F *matchedEnergy;
        float *energy;
        float *efficiency;
        TGraph *efficiencyGraph;
};

void jetEfficiency(std::string centralFileList = "", std::string forwardFileList = "", std::string backwardFileList = "") {
    // Load files
    jetEfficiencyData jets[NUM_REGIONS];
    if (centralFileList != "") {
        jets[CENTRAL].loaded = true;
        jets[CENTRAL].descriptiveName = std::string("Central");
        std::cout << "loaded " << readFileList(centralFileList, jets[CENTRAL].files) << " files in central region" << std::endl;
        jets[CENTRAL].minEta = -1.5;
        jets[CENTRAL].maxEta = 1.5;
        jets[CENTRAL].color = kRed;
        jets[CENTRAL].secondaryColor = kMagenta;
        jets[CENTRAL].marker = 21;
        jets[CENTRAL].markerSize = 1;
    }
    if (forwardFileList != "") {
        jets[FORWARD].loaded = true;
        jets[FORWARD].descriptiveName = std::string("Forward");
        std::cout << "loaded " << readFileList(forwardFileList, jets[FORWARD].files) << " files in forward region" << std::endl;
        jets[FORWARD].minEta = 1.5;
        jets[FORWARD].maxEta = 3;
        jets[FORWARD].color = kBlue;
        jets[FORWARD].secondaryColor = kCyan;
        jets[FORWARD].marker = 22;
        jets[FORWARD].markerSize = 1.4;
    }
    if (backwardFileList != "") {
        jets[BACKWARD].loaded = true;
        jets[BACKWARD].descriptiveName = std::string("Backward");
        std::cout << "loaded " << readFileList(backwardFileList, jets[BACKWARD].files) << " files in backward region" << std::endl;
        jets[BACKWARD].minEta = 1.5;    // THESE VALUES ARE NOT SET
        jets[BACKWARD].maxEta = 3;
        jets[BACKWARD].color = kGreen + 2;
        jets[BACKWARD].secondaryColor = kYellow + 1;
        jets[BACKWARD].marker = 23;
        jets[BACKWARD].markerSize = 1.4;
    }

    // Loop over all the files
    for (uint8_t jetRegion = 0; jetRegion < NUM_REGIONS; jetRegion++) {
        if (!jets[jetRegion].loaded) {  // Skip over regions we aren't studying
            continue;
        }
            // 1D histograms to store number of truth jets and reco jets for each energy bin
        jets[jetRegion].truthEnergy = new TH1F(Form("truth_energy_%s", jets[jetRegion].descriptiveName.c_str()), "", num_bins, min_energy, max_energy);
        jets[jetRegion].matchedEnergy  = new TH1F(Form("reco_energy_%s", jets[jetRegion].descriptiveName.c_str()),  "", num_bins, min_energy, max_energy);
        for (std::list<std::string>::iterator iter = jets[jetRegion].files.begin(); iter != jets[jetRegion].files.end(); ++iter) {
            TFile *inFile = TFile::Open((*iter).c_str());       // open root file
            if (inFile == nullptr) {
                std::cerr << "Could not open file " << *iter << std::endl;
                continue;
            }
            TTree *jetTree = (TTree*) inFile->Get("ntp_truthjet"); // get truthjet tree
            if (jetTree == nullptr) {
                std::cerr << "Could not file jet tree" << std::endl;
                continue;
            }

            float truthE, recoE;
            float pos[4];

            jetTree->SetBranchAddress("ge", &truthE);
            jetTree->SetBranchAddress("e", &recoE);
            jetTree->SetBranchAddress("geta", &pos[0]);
            jetTree->SetBranchAddress("gphi", &pos[1]);
            jetTree->SetBranchAddress("eta", &pos[2]);
            jetTree->SetBranchAddress("phi", &pos[3]);

            for (uint32_t i = 0; i < jetTree->GetEntries(); i++) {
                jetTree->GetEntry(i);
                if (std::isnan(truthE)) {
                    continue;
                }
                if (pos[0] < jets[jetRegion].minEta || pos[0] > jets[jetRegion].maxEta) {
                    continue;
                }
                // Do we filter on R for efficiency? Probably
                
                jets[jetRegion].truthEnergy->Fill(truthE);
                if (r2 < calculateDistance(pos)) {
                    continue;
                }
                if (std::isnan(truthE) || std::isnan(recoE)) {
                    continue;
                }
                jets[jetRegion].matchedEnergy->Fill(truthE);
                // std::cout << truthE << "\t" << recoE << std::endl;
                // std::cout << pos[0] << "\t" << pos[1] << std::endl;
            }
            inFile->Close();
        }
    }

    // Calculate efficiencies
    // efficiency = (num matched) / (num truth)
    for (uint8_t jetRegion = 0; jetRegion < NUM_REGIONS; jetRegion++) {
        if (!jets[jetRegion].loaded) {
            continue;
        }
        jets[jetRegion].energy = (float*)malloc(num_bins * sizeof(float));
        jets[jetRegion].efficiency = (float*)malloc(num_bins * sizeof(float));
        uint32_t fullBins = 0;
        for (uint32_t i = 1; i < num_bins; i++) {
            if (jets[jetRegion].truthEnergy->GetBinContent(i) == 0 || jets[jetRegion].matchedEnergy->GetBinContent(i) == 0) {
                continue;
            }
            jets[jetRegion].energy[fullBins] = jets[jetRegion].truthEnergy->GetBinCenter(i);
            jets[jetRegion].efficiency[fullBins] = jets[jetRegion].matchedEnergy->GetBinContent(i) / jets[jetRegion].truthEnergy->GetBinContent(i);
            // std::cout << jets[jetRegion].matchedEnergy->GetBinContent(i) << "\t" << jets[jetRegion].truthEnergy->GetBinContent(i) << std::endl;
            fullBins++;
        }
        jets[jetRegion].fullBins = fullBins;
        std::cout << "filled " << fullBins << " bins in " << jets[jetRegion].descriptiveName << " region" << std::endl;
    }

    // plotting
    TCanvas *efficiencyCanvas = new TCanvas("jet_efficiency", "", 1000, 500);
    efficiencyCanvas->Divide(2, 1);

    // TODO make this handle different selections of regions better
    efficiencyCanvas->cd(1);
    THStack *stack = new THStack("jet_energy", "");
    TLegend *histLegend = new TLegend(0.45, 0.70, 0.9, 0.9);
    for (uint8_t jetRegion = 0; jetRegion < NUM_REGIONS; jetRegion++) {
        if (!jets[jetRegion].loaded) {
            continue;
        }
        jets[jetRegion].truthEnergy->SetLineColor(jets[jetRegion].color);
        jets[jetRegion].matchedEnergy->SetLineColor(jets[jetRegion].secondaryColor);
        stack->Add(jets[jetRegion].truthEnergy);
        stack->Add(jets[jetRegion].matchedEnergy);
        histLegend->AddEntry(jets[jetRegion].truthEnergy, Form("%s Truth Jets", jets[jetRegion].descriptiveName.c_str()));
        histLegend->AddEntry(jets[jetRegion].matchedEnergy, Form("%s Matched Jets", jets[jetRegion].descriptiveName.c_str()));
    }
    stack->Draw("nostack");
    stack->SetTitle("Jet Energy");
    stack->GetXaxis()->SetTitle("Jet Energy");
    stack->GetYaxis()->SetTitle("Counts");
    histLegend->SetTextSize(0.035);
    histLegend->Draw();
    gPad->SetLogy();
    

    efficiencyCanvas->cd(2);
    // gPad->SetLeftMargin(0.1);
    TMultiGraph *mGraph = new TMultiGraph();
    TLegend *graphLegend = new TLegend(0.45, 0.8, 0.9, 0.9);
    for (uint8_t jetRegion = 0; jetRegion < NUM_REGIONS; jetRegion++) {
        if (!jets[jetRegion].loaded) {
            continue;
        }
        jets[jetRegion].efficiencyGraph = new TGraph(jets[jetRegion].fullBins, jets[jetRegion].energy, jets[jetRegion].efficiency);
        jets[jetRegion].efficiencyGraph->SetLineColor(jets[jetRegion].color);
        jets[jetRegion].efficiencyGraph->SetMarkerColor(jets[jetRegion].color + 2);
        jets[jetRegion].efficiencyGraph->SetMarkerSize(jets[jetRegion].markerSize);
        jets[jetRegion].efficiencyGraph->SetMarkerStyle(jets[jetRegion].marker);
        mGraph->Add(jets[jetRegion].efficiencyGraph);
        graphLegend->AddEntry(jets[jetRegion].efficiencyGraph, Form("%s Jet Efficiency", jets[jetRegion].descriptiveName.c_str()));
    }
    mGraph->SetTitle("Jet Efficiency");
    mGraph->GetXaxis()->SetTitle("Jet Energy");
    mGraph->GetYaxis()->SetTitle("Efficiency");
    mGraph->Draw("ALP");
    graphLegend->SetTextSize(0.035);
    graphLegend->Draw();
    // gPad->SetLogy();
    efficiencyCanvas->SaveAs("canvas.png");
    efficiencyCanvas->SaveAs("canvas.c");

    for (uint8_t jetRegion; jetRegion < NUM_REGIONS; jetRegion++) {
        if (!jets[jetRegion].loaded) {
            continue;
        }
        delete jets[jetRegion].truthEnergy;
        delete jets[jetRegion].matchedEnergy;
        delete jets[jetRegion].efficiencyGraph;
        free(jets[jetRegion].energy);
        free(jets[jetRegion].efficiency);

    }
    delete efficiencyCanvas;
    delete stack;
    delete histLegend;
    delete mGraph;
    delete graphLegend;



}