#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TLegend.h>
#include <THStack.h>
#include <TF1.h>
#include <TMath.h>
#include <TProfile.h>
#include <TMultiGraph.h>

#include <iostream>
#include <fstream>
#include <string>
#include <list>

// Hist Binning Parameters
const int bins_1d = 150;
const int bins_2d = 100;
const int bin_resolution = 15;
const double phiRange = TMath::Pi() + 0.1;
const double truthEtaRange = 5;
const double recoEtaRange = 5;
const double r = 0.4;

const double phiMin = -1 * phiRange;
const double phiMax = phiRange;
const double truthEtaMin = -1.7;//-1 * truthEtaRange;
const double truthEtaMax = 4;//truthEtaRange;
const double recoEtaMin = -1.7;//-1 * recoEtaRange;
const double recoEtaMax = 4;//recoEtaRange;


float calculateAngularity(float eta, float pt, float e, float parameter=-2.0) {
    return 0;
}

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
    TH2D *normalizedEtaHist = new TH2D("eta, (reco-truth)/truth", "", bin_resolution, recoEtaMin, recoEtaMax, bin_resolution, recoEtaMin, recoEtaMax);
    TH2D *normalizedPhiHist = new TH2D("phi, (reco-truth)/truth", "", bin_resolution, phiMin, phiMax, bin_resolution, phiMin, phiMax);
    
    
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
            float dEta, dPhi;
            dEta = truthEta - recoEta;
            dPhi = truthPhi - recoPhi;
            if (dPhi > TMath::Pi()) {
                dPhi -= TMath::TwoPi();
            }
            if (dPhi < -1 * TMath::Pi()) {
                dPhi += TMath::TwoPi();
            }
            if (r * r < dEta * dEta + dPhi * dPhi) {
                continue;
            }
            // if (abs(truthEta) > 1.5) {
            //     continue;
            // }
            if (!std::isnan(truthPhi) && !std::isnan(recoPhi)) {
                phiHist2d->Fill(truthPhi, recoPhi);
                normalizedPhiHist->Fill(truthPhi, (recoPhi - truthPhi));
            }
            if (!std::isnan(truthEta) && !std::isnan(recoEta))   {
                etaHist2d->Fill(truthEta, recoEta);
                normalizedEtaHist->Fill(truthEta, (recoEta - truthEta));
            }
            gphiHist->Fill(truthPhi);
            phiHist->Fill(recoPhi);
            getaHist->Fill(truthEta);
            etaHist->Fill(recoEta);
            // std::cout << truthPhi << "\t" << truthEta << std::endl;
            // std::cout << recoPhi << "\t" << recoEta << std::endl;
        
        }
        inFile->Close();
    }


    // Calculate scale and resolution of the jet angularity measurement 
    TProfile *etaProfile = normalizedEtaHist->ProfileX();
    TH1D *etaProjection = normalizedEtaHist->ProjectionX();
    etaProfile->BuildOptions(0, 0, "s");
    double *eta = (double*)malloc(bin_resolution * sizeof(double));
    double *etaScale = (double*)malloc(bin_resolution * sizeof(double));
    double *etaResolution = (double*)malloc(bin_resolution * sizeof(double));
    int fullEtaBins = 0;
    for (uint32_t i = 1; i <= bin_resolution; i++) {
        if (etaProjection->GetBinContent(i) == 0) {
            continue;
        }
        eta[fullEtaBins] = etaProfile->GetBinCenter(i);
        etaScale[fullEtaBins] = etaProfile->GetBinContent(i);
        etaResolution[fullEtaBins] = etaProfile->GetBinError(i);
        fullEtaBins++;
    }
    
    TProfile *phiProfile = normalizedPhiHist->ProfileX();
    TH1D *phiProjection = normalizedPhiHist->ProjectionX();
    phiProfile->BuildOptions(0, 0, "s");
    double *phi = (double*)malloc(bin_resolution * sizeof(double));
    double *phiScale = (double*)malloc(bin_resolution * sizeof(double));
    double *phiResolution = (double*)malloc(bin_resolution * sizeof(double));
    int fullPhiBins = 0;
    for (uint32_t i = 1; i <= bin_resolution; i++) {
        if (phiProjection->GetBinContent(i) == 0) {
            continue;
        }
        phi[fullPhiBins] = phiProfile->GetBinCenter(i);
        phiScale[fullPhiBins] = phiProfile->GetBinContent(i);
        phiResolution[fullPhiBins] = phiProfile->GetBinError(i);
        fullPhiBins++;
    }


    // Drawing
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.42);
    gStyle->SetPadLeftMargin(0.15);
    TCanvas *jetEnergy = new TCanvas("jet_energy", "", 1000, 1000);
    jetEnergy->Divide(2, 2);
    jetEnergy->cd(1);
    // phiHist2d->Draw("colz");
    // phiHist2d->SetXTitle("truth phi");
    // phiHist2d->SetYTitle("reco phi");
    // phiHist2d->SetTitle("phi, 10 GeV x 100 GeV, truth->reco");
    etaHist2d->Draw("colz");
    etaHist2d->SetXTitle("truth eta");
    etaHist2d->SetYTitle("reco eta");
    etaHist2d->SetTitle("eta, 10x100 GeV");
    etaHist2d->SetStats(false);
    gPad->SetLogz();
    

    jetEnergy->cd(2);
    // normalizedEtaHist->Draw("colz");
    // normalizedEtaHist->SetXTitle("truth eta");
    // normalizedEtaHist->SetYTitle("(reco - truth) / truth eta");
    phiHist2d->Draw("colz");
    phiHist2d->SetXTitle("truth phi");
    phiHist2d->SetYTitle("reco phi");
    phiHist2d->SetTitle("phi, 10x100 GeV");
    phiHist2d->SetStats(false);
    gPad->SetLogz();
    


    if (false) {
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
    }
    if (true) {
        jetEnergy->cd(3);
        TGraph *etaScaleGraph = new TGraph(fullEtaBins, eta, etaScale);
        etaScaleGraph->GetXaxis()->SetTitle("Eta");
        etaScaleGraph->GetYaxis()->SetTitle("Scale (Mean((reco-truth)/truth))");
        etaScaleGraph->SetTitle("Eta Scale");
        etaScaleGraph->SetMarkerStyle(3);
        etaScaleGraph->SetMarkerSize(2.5);
        // etaScaleGraph->GetXaxis()->SetRangeUser(0, 20);
        // etaScaleGraph->GetYaxis()->SetRangeUser(-0.5, 0.5);
        // etaScaleGraph->Draw("A*");

        TGraph *etaResolutionGraph = new TGraph(fullEtaBins, eta, etaResolution);
        etaResolutionGraph->GetXaxis()->SetTitle("Eta");
        etaResolutionGraph->GetYaxis()->SetTitle("Resolution (RMS((reco-truth)/truth))"); 
        etaResolutionGraph->SetTitle("Eta Resolution");
        etaResolutionGraph->SetMarkerStyle(2);
        etaResolutionGraph->SetMarkerColor(kRed);
        etaResolutionGraph->SetMarkerSize(2.5);
        // etaResolutionGraph->GetYaxis()->SetRangeUser(-0.1, 1);
        // etaResolutionGraph->Draw("A* same");

        TMultiGraph *etaMGraph = new TMultiGraph();
        etaMGraph->Add(etaScaleGraph);
        etaMGraph->Add(etaResolutionGraph);
        etaMGraph->SetTitle("Eta Scale and Resolution");
        etaMGraph->GetXaxis()->SetTitle("Eta");
        etaMGraph->GetYaxis()->SetRangeUser(-0.15,0.2);
        etaMGraph->Draw("ap");

        TLegend *etaLegend = new TLegend();
        etaLegend->AddEntry(etaScaleGraph);
        etaLegend->AddEntry(etaResolutionGraph);
        etaLegend->Draw();

        jetEnergy->cd(4);
        TGraph *phiScaleGraph = new TGraph(fullPhiBins, phi, phiScale);
        phiScaleGraph->GetXaxis()->SetTitle("phi");
        phiScaleGraph->GetYaxis()->SetTitle("Scale (Mean((reco-truth)/truth))");
        phiScaleGraph->SetTitle("Phi Scale");
        phiScaleGraph->SetMarkerStyle(3);
        phiScaleGraph->SetMarkerSize(2.5);
        // phiScaleGraph->GetXaxis()->SetRangeUser(0, 20);
        // phiScaleGraph->GetYaxis()->SetRangeUser(-0.5, 0.5);
        // phiScaleGraph->Draw("A*");

        TGraph *phiResolutionGraph = new TGraph(fullPhiBins, phi, phiResolution);
        phiResolutionGraph->GetXaxis()->SetTitle("phi");
        phiResolutionGraph->GetYaxis()->SetTitle("Resolution (RMS((reco-truth)/truth))"); 
        phiResolutionGraph->SetTitle("Phi Resolution");
        phiResolutionGraph->SetMarkerStyle(2);
        phiResolutionGraph->SetMarkerColor(kRed);
        phiResolutionGraph->SetMarkerSize(2.5);
        phiResolutionGraph->GetYaxis()->SetRangeUser(-0.2,0.2);
        // phiResolutionGraph->Draw("A* same");

        TMultiGraph *phiMGraph = new TMultiGraph();
        phiMGraph->Add(phiScaleGraph);
        phiMGraph->Add(phiResolutionGraph);
        phiMGraph->SetTitle("Phi Scale and Resolution");
        phiMGraph->GetXaxis()->SetTitle("Phi");
        phiMGraph->GetYaxis()->SetRangeUser(-0.15,0.2);
        phiMGraph->Draw("ap");

        TLegend *phiLegend = new TLegend();
        phiLegend->AddEntry(phiScaleGraph);
        phiLegend->AddEntry(phiResolutionGraph);
        phiLegend->Draw();


    }
    // jetEnergy->Draw();
    jetEnergy->SaveAs("canvas.png");


    // Some cleanup
    delete phiHist2d;
    delete etaHist2d;
    delete gphiHist;
    delete phiHist;
    delete getaHist;
    delete etaHist;
}