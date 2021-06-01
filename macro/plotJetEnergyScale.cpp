#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
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


// Hist Binning Parameters
const int bins_1d = 80;
const int bins_2d = 100;
const int bins_resolution = 30;
const int min_bin = 0;
const int e_max = 30;
const int ge_max = 80;

// Cuts
const double r = 0.5;   // r^2 < dphi^2 + deta^2
const double r2 = r * r;

// Plotting
// const std::string secondPlot("energyScale");
const std::string secondPlot("normalizedEnergyHist");
// const std::string secondPlot("etaEnergyHist");


void plotJetEnergyScale(std::string centralJets = "", std::string forwardJets = "") {
    // Initialization, i.e. loading file list and creating histogram
    std::list<std::string> fileList;

    std::ifstream centralFiles(centralJets);        // Load central jets
    std::string filePath;
    while (std::getline(centralFiles, filePath)) {
        fileList.push_back(filePath);
    }
    centralFiles.close();
    
    std::ifstream forwardFiles(forwardJets);        // Load forward jets
    while (std::getline(forwardFiles, filePath)) {
        fileList.push_back(filePath);
    }
    forwardFiles.close();

    TH2D *truthEnergyHist = new TH2D("energy_ratio, truth->reco", "", bins_2d, min_bin, e_max, bins_2d, min_bin, e_max);
    TH2D *normalizedEnergyHist = new TH2D("reco-truth/truth, truth->reco", "", bins_resolution, min_bin, e_max, bins_resolution, -4, e_max);
    TH2D *etaEnergyHist = new TH2D("eta vs energy, truth->reco", "", bins_2d, -4.5, 4.5, bins_2d, min_bin, ge_max);
    
    TH1D *truthGeHist = new TH1D("truth_ge", "", bins_1d, min_bin, ge_max);
    TH1D *truthEHist = new TH1D("truth_e", "", bins_1d, min_bin, e_max);

    // Loop over files
    for (std::list<std::string>::iterator iter = fileList.begin(); iter != fileList.end(); ++iter) {
        TFile *inFile = TFile::Open((*iter).c_str());
        TTree *jets = (TTree*) inFile->Get("ntp_truthjet");
        if (jets == nullptr) {
            std::cout << "Could not find jet tree" << std::endl;
        }
        float truthE, recoE;
        float truthEta, recoEta;
        float truthPhi, recoPhi;
        jets->SetBranchAddress("ge", &truthE);
        jets->SetBranchAddress("e", &recoE);
        jets->SetBranchAddress("geta", &truthEta);
        jets->SetBranchAddress("eta", &recoEta);
        jets->SetBranchAddress("gphi", &truthPhi);
        jets->SetBranchAddress("phi", &recoPhi);


        // Create histograms
        for (uint32_t i = 0; i < jets->GetEntries(); i++) {
            jets->GetEntry(i);

            float dEta = truthEta - recoEta;
            float dPhi = truthPhi - recoPhi;
            if (dPhi > TMath::Pi()) {
                dPhi -= TMath::TwoPi();
            }
            if (dPhi < -1 * TMath::Pi()) {
                dPhi += TMath::TwoPi();
            }
            if (r2 < dPhi * dPhi + dEta * dEta) {
                continue;
            }
            // Filling Histograms
            truthGeHist->Fill(recoE);
            truthEHist->Fill(truthE);
            if (false || (!std::isnan(truthEta) && !std::isnan(recoE))) {
                etaEnergyHist->Fill(truthEta, recoE);
            }
            if (!std::isnan(recoE) && !std::isnan(truthE)) {
                truthEnergyHist->Fill(truthE, recoE);
                normalizedEnergyHist->Fill(truthE, (recoE - truthE) / truthE);
            }
            // std::cout << truthE << "\t" << recoE << std::endl;
        
        }
        inFile->Close();
    }

    
    // Calculate energy scale and resolution
    // TProfile *profile = truthEnergyHist->ProfileX();
    TProfile *profile = normalizedEnergyHist->ProfileX();
    TH1D *projection = normalizedEnergyHist->ProjectionX();
    profile->BuildOptions(0, 0, "s");
    double *energy = (double*)malloc(bins_resolution * sizeof(double));
    double *scale = (double*)malloc(bins_resolution * sizeof(double));
    double *resolution = (double*)malloc(bins_resolution * sizeof(double));
    int fullBins = 0;
    for (uint32_t i = 1; i <= bins_resolution; i++) {
        if (projection->GetBinContent(i) == 0) {
            continue;
        }
        energy[fullBins] = profile->GetBinCenter(i);
        // scale[fullBins] = (profile->GetBinContent(i) - profile->GetBinCenter(i)) / profile->GetBinCenter(i);
        scale[fullBins] = profile->GetBinContent(i);
        resolution[fullBins] = profile->GetBinError(i);
        fullBins++;
        std::cout << projection->GetBinContent(i) << "\t" << energy[i] << "\t" << scale[i] << std::endl;
    }

    // gStyle->SetPadLeftMargin(0.15);


    // Plotting
    TCanvas *jetEnergy = new TCanvas("jet_energy", "", 1000, 1000);
    jetEnergy->Divide(2, 2);
    jetEnergy->cd(1);
    TLatex *cutList = new TLatex(0.2, 0.8, Form("r<%.1f", r));
    cutList->SetNDC();
    cutList->SetTextFont(43);
    cutList->SetTextSize(20);
    truthEnergyHist->Draw("colz");
    cutList->Draw();
    truthEnergyHist->SetXTitle("ge");
    truthEnergyHist->SetYTitle("e");
    truthEnergyHist->SetTitle("ep, 10 GeV x 100 GeV, truth->reco, |eta| < 1.5");
    truthEnergyHist->SetStats(false);
    gPad->SetLogz();
    
    

    jetEnergy->cd(2);
    if (secondPlot == "normalizedEnergyHist") {
        normalizedEnergyHist->Draw("colz");
        normalizedEnergyHist->SetXTitle("ge");
        normalizedEnergyHist->SetYTitle("(reco - truth) / truth");
        normalizedEnergyHist->SetTitle("ep, 10x100 GeV, truth->reco, Jet Scale");
        normalizedEnergyHist->SetStats(false);
    }
    else if (secondPlot == "etaEnergyHist") {
        gStyle->SetStatX(.4);
        gStyle->SetStatY(0.92);
        etaEnergyHist->Draw("colz");
        etaEnergyHist->SetXTitle("eta");
        etaEnergyHist->SetYTitle("ge");
        etaEnergyHist->SetTitle("ep, 10x100 GeV, truth->reco, Truth Energy vs Eta");
        etaEnergyHist->SetStats(false);

    }
    else if (secondPlot == "energyScale") {
        profile->Draw("");
        profile->SetXTitle("ge");
        profile->SetYTitle("(reco-truth)/truth");
        profile->SetTitle("Profile on X");
    }
    gPad->SetLogz();
    

    jetEnergy->cd(3);
    TGraph *jetScale = new TGraph(fullBins, energy, scale);
    jetScale->GetXaxis()->SetTitle("Energy");
    jetScale->GetYaxis()->SetTitle("Scale (Mean((reco-truth)/truth))");
    jetScale->SetTitle("Jet Scale");
    jetScale->SetMarkerSize(2.5);
    // jetScale->GetXaxis()->SetRangeUser(0, 20);
    // jetScale->GetYaxis()->SetRangeUser(-0.1, 0.2);
    jetScale->Draw("A*");


    
    jetEnergy->cd(4);
    TGraph *jetResolution = new TGraph(fullBins, energy, resolution);
    jetResolution->GetXaxis()->SetTitle("Energy");
    jetResolution->GetYaxis()->SetTitle("Resolution (RMS((reco-truth)/truth))"); 
    jetResolution->SetTitle("Jet Resolution");
    jetResolution->SetMarkerSize(2.5);
    jetResolution->GetYaxis()->SetRangeUser(-0.1, 1);
    jetResolution->Draw("A*");

    jetEnergy->SaveAs("canvas.png");


    // Some cleanup
    delete jetEnergy;
    delete cutList;
    delete jetScale;
    delete truthEnergyHist;
    delete truthGeHist;
    delete truthEHist;

    free(energy);
    free(scale);
}