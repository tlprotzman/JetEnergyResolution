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
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <string>
#include <list>

// Hist Binning Parameters
const int bins_1d = 150;
const int bins_2d = 400;
const double phiRange = TMath::Pi() + 0.1;
const double truthEtaRange = 4.5;
const double recoEtaRange = 1.7;


const double phiMin = -1 * phiRange;
const double phiMax = phiRange;
const double truthEtaMin = -1 * truthEtaRange;
const double truthEtaMax = truthEtaRange;
const double recoEtaMin = -1 * recoEtaRange;
const double recoEtaMax = recoEtaRange;


void plotJetAngularResolution(std::string inFilePath = "smallfilelist.txt") {
    // Initialization, i.e. loading file list and creating histogram
    std::list<std::string> fileList;

    std::ifstream files(inFilePath);
    std::string filePath;
    while (std::getline(files, filePath)) {
        fileList.push_back(filePath);
    }
    files.close();

    TH2D *phiHist2d = new TH2D("phi, truth->reco", "", bins_2d, phiMin, phiMax, bins_2d, phiMin, phiMax);
    TH2D *etaHist2d = new TH2D("eta, truth->reco", "", bins_2d, truthEtaMin, truthEtaMax, bins_2d, recoEtaMin, recoEtaMax);
    
    TH1D *gphiHist = new TH1D("gphi", "", bins_1d, phiMin, phiMax);
    TH1D *phiHist = new TH1D("phi", "", bins_1d, phiMin, phiMax);
    TH1D *getaHist = new TH1D("eta", "", bins_1d, truthEtaMin, truthEtaMax);
    TH1D *etaHist = new TH1D("geta", "", bins_1d, recoEtaMin, recoEtaMax);

    // Loop over files
    for (std::list<std::string>::iterator iter = fileList.begin(); iter != fileList.end(); ++iter) {
        TFile *inFile = TFile::Open((*iter).c_str());
        TTree *truthJets = (TTree*) inFile->Get("ntp_truthjet");
        TTree *recoJets = (TTree*) inFile->Get("ntp_recojet");
        if (truthJets == nullptr || recoJets == nullptr) {
            std::cout << "Could not find jet tree" << std::endl;
        }
        float truthPhi, truthEta;
        float recoPhi, recoEta;
        truthJets->SetBranchAddress("gphi", &truthPhi);
        truthJets->SetBranchAddress("geta", &truthEta);
        truthJets->SetBranchAddress("phi", &recoPhi);
        truthJets->SetBranchAddress("eta", &recoEta);


        // Create histograms
        for (uint32_t i = 0; i < truthJets->GetEntries(); i++) {
            truthJets->GetEntry(i);
            recoJets->GetEntry(i);
            if (!std::isnan(truthPhi) && !std::isnan(recoPhi)) {phiHist2d->Fill(truthPhi, recoPhi);}
            if (!std::isnan(truthEta) && !std::isnan(recoEta))   {etaHist2d->Fill(truthEta, recoEta);}
            gphiHist->Fill(truthPhi);
            phiHist->Fill(recoPhi);
            getaHist->Fill(truthEta);
            etaHist->Fill(recoEta);
            // std::cout << truthPhi << "\t" << truthEta << std::endl;
            // std::cout << recoPhi << "\t" << recoEta << std::endl;
        
        }
        inFile->Close();
    }


    // Drawing
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.42);
    TCanvas *jetEnergy = new TCanvas("jet_energy", "", 1000, 1000);
    jetEnergy->Divide(2, 2);
    jetEnergy->cd(1);
    phiHist2d->Draw("colz");
    phiHist2d->SetXTitle("truth phi");
    phiHist2d->SetYTitle("reco phi");
    phiHist2d->SetTitle("phi, 10 GeV x 100 GeV, truth->reco");
    gPad->SetLogz();
    

    jetEnergy->cd(2);
    etaHist2d->Draw("colz");
    etaHist2d->SetXTitle("truth eta");
    etaHist2d->SetYTitle("reco eta");
    etaHist2d->SetTitle("eta, 10 GeV x 100 GeV, truth->reco");
    // line.Draw("lsame");
    gPad->SetLogz();
    

    jetEnergy->cd(3);
    gphiHist->SetLineColor(2);
    phiHist->SetLineColor(1);

    THStack *geStack = new THStack("ge", "");
    geStack->Add(gphiHist);
    geStack->Add(phiHist);
    geStack->Draw("nostack hist lp");
    gPad->SetLogy();
    geStack->GetXaxis()->SetTitle("phi");
    geStack->GetYaxis()->SetTitle("counts");
    geStack->SetTitle("phi distribution");

    TLegend *geLegend = new TLegend(0.80, 0.80, 0.9, 0.9);
    geLegend->AddEntry(gphiHist, "gphi");
    geLegend->AddEntry(phiHist, "phi");
    geLegend->Draw();
    

    jetEnergy->cd(4);
    getaHist->SetLineColor(2);
    etaHist->SetLineColor(1);

    THStack *eStack = new THStack("e", "");
    eStack->Add(getaHist);
    eStack->Add(etaHist);
    eStack->Draw("nostack hist lp");
    gPad->SetLogy();
    eStack->GetXaxis()->SetTitle("eta");
    eStack->GetYaxis()->SetTitle("counts");
    eStack->SetTitle("eta distribution");

    TLegend *eLegend = new TLegend(0.80, 0.60, 0.9, 0.7);
    eLegend->AddEntry(getaHist, "geta");
    eLegend->AddEntry(etaHist, "eta");
    eLegend->Draw();
    jetEnergy->Draw();


    // Some cleanup
    delete phiHist2d;
    delete etaHist2d;
    delete gphiHist;
    delete phiHist;
    delete getaHist;
    delete etaHist;
}