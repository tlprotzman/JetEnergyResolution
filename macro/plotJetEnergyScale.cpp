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


void plotJetEnergyScale() {
    // Initialization, i.e. loading file list and creating histogram
    std::list<std::string> fileList;

    // std::ifstream files("filelist.txt");
    std::ifstream files("smallfilelist.txt");
    std::string filePath;
    while (std::getline(files, filePath)) {
        fileList.push_back(filePath);
    }

    TH2D *truthHist = new TH2D("energy_ratio, truth->reco", "", bins_2d, min_bin, e_max, bins_2d, min_bin, e_max);
    TH2D *recoHist = new TH2D("energy_ratio, reco->truth", "", bins_2d, min_bin, e_max, bins_2d, min_bin, e_max);
    
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
            std::cout << "?!?!?!?" << std::endl;
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
            if (!std::isnan(truthGe) && !std::isnan(truthE)) {truthHist->Fill(truthGe, truthE);}
            if (!std::isnan(recoGe) && !std::isnan(recoE))   {recoHist->Fill(recoGe, recoE);}
            truthGeHist->Fill(truthGe);
            truthEHist->Fill(truthE);
            recoGeHist->Fill(recoGe);
            recoEHist->Fill(recoE);
            // std::cout << truthE << "\t" << truthGe << std::endl;
            // std::cout << recoE << "\t" << recoGe << std::endl;
        
        }
        inFile->Close();
    }

    TCanvas *canvas = new TCanvas("canvas", "", 1000, 1000);
    canvas->Divide(2, 2);
    canvas->cd(1);
    truthHist->Draw("colz");
    truthHist->SetXTitle("ge");
    truthHist->SetYTitle("e");
    truthHist->SetTitle("ep, 10 GeV x 100 GeV, truth->reco");
    line.Draw("lsame");
    gPad->SetLogz();
    

    canvas->cd(2);
    recoHist->Draw("colz");
    recoHist->SetXTitle("ge");
    recoHist->SetYTitle("e");
    recoHist->SetTitle("ep, 10 GeV x 100 GeV, reco->truth");
    // line.Draw("lsame");
    gPad->SetLogz();
    

    canvas->cd(3);
    truthGeHist->SetLineColor(2);
    recoGeHist->SetLineColor(1);

    THStack *stack = new THStack("ge", "");
    stack->Add(truthGeHist);
    stack->Add(recoGeHist);
    stack->Draw("nostack hist lp");
    gPad->SetLogy();
    stack->GetXaxis()->SetTitle("ge");
    stack->GetYaxis()->SetTitle("counts");
    stack->SetTitle("ge distribution");

    TLegend *legend = new TLegend(0.65, 0.55, 0.85, 0.75);
    legend->AddEntry(truthGeHist, "truth->reco");
    legend->AddEntry(recoGeHist, "reco->truth");
    legend->Draw();
    

    canvas->cd(4);
    truthEHist->SetLineColor(2);
    recoEHist->SetLineColor(1);

    THStack *stack2 = new THStack("e", "");
    stack2->Add(truthEHist);
    stack2->Add(recoEHist);
    stack2->Draw("nostack hist lp");
    gPad->SetLogy();
    stack2->GetXaxis()->SetTitle("e");
    stack2->GetYaxis()->SetTitle("counts");
    stack2->SetTitle("e distribution");

    TLegend *legend2 = new TLegend(0.65, 0.55, 0.85, 0.75);
    legend2->AddEntry(truthEHist, "truth->reco");
    legend2->AddEntry(recoEHist, "reco->truth");
    legend2->Draw();
    canvas->Draw();
}