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

#include "common.cpp"

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

class jetAngularData: public jetData {
    public:
        TH2F *phiHist;
        TH2F *etaHist;
        TH2F *normalizedPhiHist;
        TH2F *normalizedEtaHist;

        TH1D *etaProjection;
        TProfile *etaProfile;
        double *eta, *etaScale, *etaResolution;
        int fullEtaBins;
        TGraph *etaScaleGraph;
        TGraph *etaResolutionGraph;

        TH1D *phiProjection;
        TProfile *phiProfile;
        double *phi, *phiScale, *phiResolution;
        int fullPhiBins;
        TGraph *phiScaleGraph;
        TGraph *phiResolutionGraph;
};

void plotJetAngularResolution(std::string centralFileList = "", std::string forwardFileList = "", std::string backwardFileList = "") {
    // Initialization, i.e. loading file list and creating histogram
    jetAngularData jets[NUM_REGIONS];
    if (centralFileList != "") {
        jets[CENTRAL].loaded = true;
        jets[CENTRAL].descriptiveName = std::string("Central");
        std::cout << "loaded " << readFileList(centralFileList, jets[CENTRAL].files) << " files in central region" << std::endl;
        jets[CENTRAL].color = kRed;
        jets[CENTRAL].secondaryColor = kMagenta;
        jets[CENTRAL].marker = 21;
        jets[CENTRAL].markerSize = 1;
    }
    if (forwardFileList != "") {
        jets[FORWARD].loaded = true;
        jets[FORWARD].descriptiveName = std::string("Forward");
        std::cout << "loaded " << readFileList(forwardFileList, jets[FORWARD].files) << " files in forward region" << std::endl;
        jets[FORWARD].color = kBlue;
        jets[FORWARD].secondaryColor = kCyan;
        jets[FORWARD].marker = 21;
        jets[FORWARD].markerSize = 1;
    }
    if (backwardFileList != "") {
        jets[BACKWARD].loaded = true;
        jets[BACKWARD].descriptiveName = std::string("Backward");
        std::cout << "loaded " << readFileList(backwardFileList, jets[BACKWARD].files) << " files in backward region" << std::endl;
        jets[BACKWARD].color = kGreen + 2;
        jets[BACKWARD].secondaryColor = kYellow + 1;
        jets[BACKWARD].marker = 21;
        jets[BACKWARD].markerSize = 1;
    }


    for (uint8_t jetRegion = 0; jetRegion < NUM_REGIONS; jetRegion++) {
        if (!jets[jetRegion].loaded) {
            continue;
        }
        jets[jetRegion].phiHist = new TH2F(Form("%s phi", jets[jetRegion].descriptiveName.c_str()), "", bins_2d, phiMin, phiMax, bins_2d, phiMin, phiMax);
        jets[jetRegion].etaHist = new TH2F(Form("%s eta", jets[jetRegion].descriptiveName.c_str()), "", bins_2d, truthEtaMin, truthEtaMax, bins_2d, recoEtaMin, recoEtaMax);
        jets[jetRegion].normalizedEtaHist = new TH2F(Form("%s eta, (reco-truth)/truth", jets[jetRegion].descriptiveName.c_str()), "", bin_resolution, recoEtaMin, recoEtaMax, bin_resolution, recoEtaMin, recoEtaMax);
        jets[jetRegion].normalizedPhiHist = new TH2F(Form("%s phi, (reco-truth)/truth", jets[jetRegion].descriptiveName.c_str()), "", bin_resolution, phiMin, phiMax, bin_resolution, phiMin, phiMax);

        // Loop over files
        for (std::list<std::string>::iterator iter = jets[jetRegion].files.begin(); iter != jets[jetRegion].files.end(); ++iter) {
            TFile *inFile = TFile::Open((*iter).c_str());
            TTree *truthJets = (TTree*) inFile->Get("ntp_truthjet");
            if (truthJets == nullptr) {
                std::cout << "Could not find jet tree" << std::endl;
            }
            float pos[4];
            truthJets->SetBranchAddress("geta", &pos[0]);
            truthJets->SetBranchAddress("gphi", &pos[1]);
            truthJets->SetBranchAddress("eta", &pos[2]);
            truthJets->SetBranchAddress("phi", &pos[3]);


            // Create histograms
            for (uint32_t i = 0; i < truthJets->GetEntries(); i++) {
                truthJets->GetEntry(i);
                if (r * r < calculateDistance(pos)) {
                    continue;
                }
                // if (abs(truthEta) > 1.5) {
                //     continue;
                // }
                if (!std::isnan(pos[1]) && !std::isnan(pos[3])) {
                    jets[jetRegion].phiHist->Fill(pos[1], pos[3]);
                    jets[jetRegion].normalizedPhiHist->Fill(pos[1], (pos[3] - pos[1]));
                }
                if (!std::isnan(pos[0]) && !std::isnan(pos[2]))   {
                    jets[jetRegion].etaHist->Fill(pos[0], pos[2]);
                    jets[jetRegion].normalizedEtaHist->Fill(pos[0], (pos[2] - pos[0]));
                }
            
            }
            inFile->Close();
        }
    }


    // Calculate scale and resolution of the jet angularity measurement 
    for (uint8_t jetRegion = 0; jetRegion < NUM_REGIONS; jetRegion++) {
        if (!jets[jetRegion].loaded) {
            continue;
        }
        jets[jetRegion].etaProfile = jets[jetRegion].normalizedEtaHist->ProfileX();
        jets[jetRegion].etaProjection = jets[jetRegion].normalizedEtaHist->ProjectionX();
        jets[jetRegion].etaProfile->BuildOptions(0, 0, "s");
        jets[jetRegion].eta = (double*)malloc(bin_resolution * sizeof(double));
        jets[jetRegion].etaScale = (double*)malloc(bin_resolution * sizeof(double));
        jets[jetRegion].etaResolution = (double*)malloc(bin_resolution * sizeof(double));
        int fullEtaBins = 0;
        for (uint32_t i = 1; i <= bin_resolution; i++) {
            if (jets[jetRegion].etaProjection->GetBinContent(i) == 0) {
                continue;
            }
            jets[jetRegion].eta[fullEtaBins] = jets[jetRegion].etaProfile->GetBinCenter(i);
            jets[jetRegion].etaScale[fullEtaBins] = jets[jetRegion].etaProfile->GetBinContent(i);
            jets[jetRegion].etaResolution[fullEtaBins] = jets[jetRegion].etaProfile->GetBinError(i);
            fullEtaBins++;
        }
        jets[jetRegion].fullEtaBins = fullEtaBins;
        
        jets[jetRegion].phiProfile = jets[jetRegion].normalizedPhiHist->ProfileX();
        jets[jetRegion].phiProjection = jets[jetRegion].normalizedPhiHist->ProjectionX();
        jets[jetRegion].phiProfile->BuildOptions(0, 0, "s");
        jets[jetRegion].phi = (double*)malloc(bin_resolution * sizeof(double));
        jets[jetRegion].phiScale = (double*)malloc(bin_resolution * sizeof(double));
        jets[jetRegion].phiResolution = (double*)malloc(bin_resolution * sizeof(double));
        int fullPhiBins = 0;
        for (uint32_t i = 1; i <= bin_resolution; i++) {
            if (jets[jetRegion].phiProjection->GetBinContent(i) == 0) {
                continue;
            }
            jets[jetRegion].phi[fullPhiBins] = jets[jetRegion].phiProfile->GetBinCenter(i);
            jets[jetRegion].phiScale[fullPhiBins] = jets[jetRegion].phiProfile->GetBinContent(i);
            jets[jetRegion].phiResolution[fullPhiBins] = jets[jetRegion].phiProfile->GetBinError(i);
            fullPhiBins++;
        }
        jets[jetRegion].fullPhiBins = fullPhiBins;
    }


    // Drawing
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.42);
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadTopMargin(0.12);
    gStyle->SetPadBottomMargin(0.12);
    TCanvas *jetEnergy = new TCanvas("jet_energy", "", 1000, 500);
    jetEnergy->Divide(2, 1);
    // jetEnergy->cd(1);
    // // phiHist2d->Draw("colz");
    // // phiHist2d->SetXTitle("truth phi");
    // // phiHist2d->SetYTitle("reco phi");
    // // phiHist2d->SetTitle("phi, 10 GeV x 100 GeV, truth->reco");
    // etaHist2d->Draw("colz");
    // etaHist2d->SetXTitle("truth eta");
    // etaHist2d->SetYTitle("reco eta");
    // etaHist2d->SetTitle("eta, 10x100 GeV");
    // etaHist2d->SetStats(false);
    // gPad->SetLogz();
    

    // jetEnergy->cd(2);
    // // normalizedEtaHist->Draw("colz");
    // // normalizedEtaHist->SetXTitle("truth eta");
    // // normalizedEtaHist->SetYTitle("(reco - truth) / truth eta");
    // phiHist2d->Draw("colz");
    // phiHist2d->SetXTitle("truth phi");
    // phiHist2d->SetYTitle("reco phi");
    // phiHist2d->SetTitle("phi, 10x100 GeV");
    // phiHist2d->SetStats(false);
    // gPad->SetLogz();
    


    // if (false) {
    //     jetEnergy->cd(3);
    //     gphiHist->SetLineColor(2);
    //     phiHist->SetLineColor(1);
    //     THStack *geStack = new THStack("ge", "");
    //     geStack->Add(gphiHist);
    //     geStack->Add(phiHist);
    //     geStack->Draw("nostack hist lp");
    //     gPad->SetLogy();
    //     geStack->GetXaxis()->SetTitle("phi");
    //     geStack->GetYaxis()->SetTitle("counts");
    //     geStack->SetTitle("phi distribution");

    //     TLegend *geLegend = new TLegend(0.80, 0.80, 0.9, 0.9);
    //     geLegend->AddEntry(gphiHist, "gphi");
    //     geLegend->AddEntry(phiHist, "phi");
    //     geLegend->Draw();
        

    //     jetEnergy->cd(4);
    //     getaHist->SetLineColor(2);
    //     etaHist->SetLineColor(1);

    //     THStack *eStack = new THStack("e", "");
    //     eStack->Add(getaHist);
    //     eStack->Add(etaHist);
    //     eStack->Draw("nostack hist lp");
    //     gPad->SetLogy();
    //     eStack->GetXaxis()->SetTitle("eta");
    //     eStack->GetYaxis()->SetTitle("counts");
    //     eStack->SetTitle("eta distribution");

    //     TLegend *eLegend = new TLegend(0.80, 0.60, 0.9, 0.7);
    //     eLegend->AddEntry(getaHist, "geta");
    //     eLegend->AddEntry(etaHist, "eta");
    //     eLegend->Draw();
    // }
    if (true) {
        jetEnergy->cd(1);
        TMultiGraph *etaMGraph = new TMultiGraph();
        TLegend *etaLegend = new TLegend();
        for (uint8_t jetRegion = 0; jetRegion < NUM_REGIONS; jetRegion++) {
            if (!jets[jetRegion].loaded) {
                continue;
            }
            jets[jetRegion].etaScaleGraph = new TGraph(jets[jetRegion].fullEtaBins, jets[jetRegion].eta, jets[jetRegion].etaScale);
            jets[jetRegion].etaScaleGraph->SetMarkerStyle(3);
            jets[jetRegion].etaScaleGraph->SetMarkerSize(2.5);
            jets[jetRegion].etaScaleGraph->SetMarkerColor(jets[jetRegion].color);


            jets[jetRegion].etaResolutionGraph = new TGraph(jets[jetRegion].fullEtaBins, jets[jetRegion].eta, jets[jetRegion].etaResolution);
            jets[jetRegion].etaResolutionGraph->SetMarkerStyle(5);
            jets[jetRegion].etaResolutionGraph->SetMarkerColor(jets[jetRegion].secondaryColor);
            jets[jetRegion].etaResolutionGraph->SetMarkerSize(2.5);

            etaMGraph->Add(jets[jetRegion].etaScaleGraph);
            etaMGraph->Add(jets[jetRegion].etaResolutionGraph);
            
            etaLegend->AddEntry(jets[jetRegion].etaScaleGraph, Form("%s Eta Scale", jets[jetRegion].descriptiveName.c_str()));
            etaLegend->AddEntry(jets[jetRegion].etaResolutionGraph, Form("%s Eta Resolution", jets[jetRegion].descriptiveName.c_str()));
        }

        etaMGraph->SetTitle("Eta Scale and Resolution");
        etaMGraph->GetXaxis()->SetTitle("Eta");
        etaMGraph->GetYaxis()->SetRangeUser(-0.15,0.2);
        etaMGraph->Draw("ap");

        etaLegend->Draw();

        jetEnergy->cd(2);
        TMultiGraph *phiMGraph = new TMultiGraph();
        TLegend *phiLegend = new TLegend();
        for (uint8_t jetRegion = 0; jetRegion < NUM_REGIONS; jetRegion++) {
            if (!jets[jetRegion].loaded) {
                continue;
            }
            jets[jetRegion].phiScaleGraph = new TGraph(jets[jetRegion].fullPhiBins, jets[jetRegion].phi, jets[jetRegion].phiScale);
            jets[jetRegion].phiScaleGraph->SetMarkerStyle(3);
            jets[jetRegion].phiScaleGraph->SetMarkerSize(2.5);
            jets[jetRegion].phiScaleGraph->SetMarkerColor(jets[jetRegion].color);

            jets[jetRegion].phiResolutionGraph = new TGraph(jets[jetRegion].fullPhiBins, jets[jetRegion].phi, jets[jetRegion].phiResolution);
            jets[jetRegion].phiResolutionGraph->SetMarkerStyle(5);
            jets[jetRegion].phiResolutionGraph->SetMarkerSize(2.5);
            jets[jetRegion].phiResolutionGraph->SetMarkerColor(jets[jetRegion].secondaryColor);

            phiMGraph->Add(jets[jetRegion].phiScaleGraph);
            phiMGraph->Add(jets[jetRegion].phiResolutionGraph);

            phiLegend->AddEntry(jets[jetRegion].phiScaleGraph, Form("%s Phi Scale", jets[jetRegion].descriptiveName.c_str()));
            phiLegend->AddEntry(jets[jetRegion].phiResolutionGraph, Form("%s Phi Resolution", jets[jetRegion].descriptiveName.c_str()));
        }

        phiMGraph->SetTitle("Phi Scale and Resolution");
        phiMGraph->GetXaxis()->SetTitle("Phi");
        phiMGraph->GetYaxis()->SetRangeUser(-0.15,0.2);
        phiMGraph->Draw("ap");

        phiLegend->Draw();


    }
    // jetEnergy->Draw();
    jetEnergy->SaveAs("jetAngularScale.png");
    jetEnergy->SaveAs("jetAngularScale.c");


    // Some cleanup
    // delete phiHist2d;
    // delete etaHist2d;
    // delete gphiHist;
    // delete phiHist;
    // delete getaHist;
    // delete etaHist;
}