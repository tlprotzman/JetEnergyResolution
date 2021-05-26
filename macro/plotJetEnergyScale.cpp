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
#include <TF1.h>

#include <iostream>
#include <fstream>
#include <string>
#include <list>

const int bins_2d = 400;
const int bins_1d = 80;
const int min_bin = 0;
const int e_max = 20;
const int ge_max = 80;


void plotJetEnergyScale(std::string inFilePath = "smallfilelist.txt") {
    // Initialization, i.e. loading file list and creating histogram
    std::list<std::string> fileList;

    std::ifstream files(inFilePath);
    std::string filePath;
    while (std::getline(files, filePath)) {
        fileList.push_back(filePath);
    }

    TH2D *truthEnergyHist = new TH2D("energy_ratio, truth->reco", "", bins_2d, min_bin, e_max, bins_2d, min_bin, e_max);
    TH2D *recoEnergyHist = new TH2D("energy_ratio, reco->truth", "", bins_2d, min_bin, e_max, bins_2d, min_bin, e_max);
    
    TH1D *truthGeHist = new TH1D("truth_ge", "", bins_1d, min_bin, ge_max);
    TH1D *truthEHist = new TH1D("truth_e", "", bins_1d, min_bin, e_max);
    TH1D *recoGeHist = new TH1D("reco_ge", "", bins_1d, min_bin, ge_max);
    TH1D *recoEHist = new TH1D("reco_e", "", bins_1d, min_bin, e_max);

    // Loop over files
    for (std::list<std::string>::iterator iter = fileList.begin(); iter != fileList.end(); ++iter) {
        TFile *inFile = TFile::Open((*iter).c_str());
        TTree *truthJets = (TTree*) inFile->Get("ntp_truthjet");
        TTree *recoJets = (TTree*) inFile->Get("ntp_recojet");
        if (truthJets == nullptr || recoJets == nullptr) {
            std::cout << "Could not find jet tree" << std::endl;
        }
        float truthE, truthGe;
        float recoE, recoGe;
        truthJets->SetBranchAddress("e", &truthE);
        truthJets->SetBranchAddress("ge", &truthGe);
        recoJets->SetBranchAddress("e", &recoE);
        recoJets->SetBranchAddress("ge", &recoGe);


        // Create histograms
        for (uint32_t i = 0; i < truthJets->GetEntries(); i++) {
            truthJets->GetEntry(i);
            recoJets->GetEntry(i);
            if (!std::isnan(truthGe) && !std::isnan(truthE)) {truthEnergyHist->Fill(truthGe, truthE);}
            if (!std::isnan(recoGe) && !std::isnan(recoE))   {recoEnergyHist->Fill(recoGe, recoE);}
            truthGeHist->Fill(truthGe);
            truthEHist->Fill(truthE);
            recoGeHist->Fill(recoGe);
            recoEHist->Fill(recoE);
            // std::cout << truthE << "\t" << truthGe << std::endl;
            // std::cout << recoE << "\t" << recoGe << std::endl;
        
        }
        inFile->Close();
    }

    TCanvas *jetEnergy = new TCanvas("jet_energy", "", 1000, 1000);
    jetEnergy->Divide(2, 2);
    jetEnergy->cd(1);
    truthEnergyHist->Draw("colz");
    truthEnergyHist->SetXTitle("ge");
    truthEnergyHist->SetYTitle("e");
    truthEnergyHist->SetTitle("ep, 10 GeV x 100 GeV, truth->reco");
    gPad->SetLogz();
    

    jetEnergy->cd(2);
    recoEnergyHist->Draw("colz");
    recoEnergyHist->SetXTitle("ge");
    recoEnergyHist->SetYTitle("e");
    recoEnergyHist->SetTitle("ep, 10 GeV x 100 GeV, reco->truth");
    // line.Draw("lsame");
    gPad->SetLogz();
    

    jetEnergy->cd(3);
    truthGeHist->SetLineColor(2);
    recoGeHist->SetLineColor(1);

    THStack *geStack = new THStack("ge", "");
    geStack->Add(truthGeHist);
    geStack->Add(recoGeHist);
    geStack->Draw("nostack hist lp");
    gPad->SetLogy();
    geStack->GetXaxis()->SetTitle("ge");
    geStack->GetYaxis()->SetTitle("counts");
    geStack->SetTitle("ge distribution");

    TLegend *geLegend = new TLegend(0.65, 0.55, 0.85, 0.75);
    geLegend->AddEntry(truthGeHist, "truth->reco");
    geLegend->AddEntry(recoGeHist, "reco->truth");
    geLegend->Draw();
    

    jetEnergy->cd(4);
    truthEHist->SetLineColor(2);
    recoEHist->SetLineColor(1);

    THStack *eStack = new THStack("e", "");
    eStack->Add(truthEHist);
    eStack->Add(recoEHist);
    eStack->Draw("nostack hist lp");
    gPad->SetLogy();
    eStack->GetXaxis()->SetTitle("e");
    eStack->GetYaxis()->SetTitle("counts");
    eStack->SetTitle("e distribution");

    TLegend *eLegend = new TLegend(0.65, 0.55, 0.85, 0.75);
    eLegend->AddEntry(truthEHist, "truth->reco");
    eLegend->AddEntry(recoEHist, "reco->truth");
    eLegend->Draw();
    jetEnergy->Draw();
}