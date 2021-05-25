//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in JetEnergyResolution.h.
//
// JetEnergyResolution(const std::string &name = "JetEnergyResolution")
// everything is keyed to JetEnergyResolution, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// JetEnergyResolution::~JetEnergyResolution()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int JetEnergyResolution::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int JetEnergyResolution::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int JetEnergyResolution::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT; 
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT; 
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int JetEnergyResolution::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int JetEnergyResolution::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int JetEnergyResolution::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int JetEnergyResolution::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void JetEnergyResolution::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

#include "JetEnergyResolution.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>

#include <g4jets/Jet.h>
#include <g4jets/JetMap.h>

#include <g4eval/JetEvalStack.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

//____________________________________________________________________________..
JetEnergyResolution::JetEnergyResolution(const std::string &name):
 SubsysReco(name)
{
  std::cout << "JetEnergyResolution::JetEnergyResolution(const std::string &name) Calling ctor" << std::endl;
  outfile = new TFile();
  recoJetTree = new TTree("RecoJetTree", "A tree containing reconstructed jets");
  recoJetTree->Branch("recoPt", &recoPt, "recoPt/D");
  recoJetTree->Branch("recoEnergy", &recoEnergy, "recoEnergy/D");
  recoJetTree->Branch("truthPt", &truthPt, "truthPt/D");
  recoJetTree->Branch("truthEnergy", &truthEnergy, "truthEnergy/D");
  recoJetTree->Branch("dR", &dR, "dR/D");
}

//____________________________________________________________________________..
JetEnergyResolution::~JetEnergyResolution()
{
  std::cout << "JetEnergyResolution::~JetEnergyResolution() Calling dtor" << std::endl;
  delete recoJetTree;
}

//____________________________________________________________________________..
int JetEnergyResolution::Init(PHCompositeNode *topNode)
{
  std::cout << "JetEnergyResolution::Init(PHCompositeNode *topNode) Initializing" << std::endl;
  outfile = new TFile("out.root", "RECREATE"); // Create file for writing to
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
// int JetEnergyResolution::InitRun(PHCompositeNode *topNode)
// {
//   std::cout << "JetEnergyResolution::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
//   return Fun4AllReturnCodes::EVENT_OK;
// }

//____________________________________________________________________________..
int JetEnergyResolution::process_event(PHCompositeNode *topNode)
{
  std::cout << "JetEnergyResolution::process_event(PHCompositeNode *topNode) Processing Event" << std::endl;
  // Copied from AnaTutorial->getReconstructedJets
  JetMap *recoJets = findNode::getClass<JetMap>(topNode, "AntiKt_Tower_r04");
  JetMap *truthJets = findNode::getClass<JetMap>(topNode, "AntiKt_Truth_r04");
  if (!jetEvalStack) {
    jetEvalStack = new JetEvalStack(topNode, "AntiKt_Tower_r04", "AntiKt_Truth_r04");
  }
  jetEvalStack->next_event(topNode);
  JetRecoEval *recoEval = jetEvalStack->get_reco_eval();
  if (!recoJets) {
    std::cout << "No reconstructed jet node: " << PHWHERE << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }
  for (JetMap::Iter recoIter = recoJets->begin(); recoIter != recoJets->end(); ++recoIter) {
    Jet *recoJet = recoIter->second;
    recoPt = recoJet->get_pt();
    recoEnergy = recoJet->get_e();

    Jet *truthJet = recoEval->max_truth_jet_by_energy(recoJet);
    if (truthJet->get_pt() < 5) {
      continue;
    }
    if (truthJet) {
      truthPt = truthJet->get_pt();
      truthEnergy = truthJet->get_e();
      dR = -100;
    }
    else if (truthJets) {
      float closestJet = 9999;
      float recoEta = recoJet->get_eta();
      float recoPhi = recoJet->get_phi();

      for (JetMap::Iter truthIter = truthJets->begin(); truthIter != truthJets->end(); ++truthIter) {
        const Jet *truthJet = truthIter->second;
        truthPt = truthJet->get_pt();
        truthEnergy = truthJet->get_e();
        float truthEta = truthJet->get_eta();
        float truthPhi = truthJet->get_phi();

        float dPhi = recoPhi - truthPhi;
        if (dPhi > TMath::Pi()) {
          dPhi -= TMath::TwoPi();
        }
        if (dPhi < -1 * TMath::Pi()) {
          dPhi += TMath::TwoPi();
        }

        float dEta = recoEta - truthEta;
        // Calculate distance in eta phi space
        dR = sqrt(pow(dPhi, 2) + pow(dEta, 2));

        if (dR < recoJets->get_par() && dR < closestJet) {
          closestJet = dR;
        }
      }
      dR = closestJet;
    }
    std::cout << "filling tree" << std::endl;
    if (dR < 9998) {
      recoJetTree->Fill();
    }
  }
  std::cout << "about to return from here" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
// int JetEnergyResolution::ResetEvent(PHCompositeNode *topNode)
// {
//   std::cout << "JetEnergyResolution::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
//   return Fun4AllReturnCodes::EVENT_OK;
// }

//____________________________________________________________________________..
// int JetEnergyResolution::EndRun(const int runnumber)
// {
//   std::cout << "JetEnergyResolution::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
//   return Fun4AllReturnCodes::EVENT_OK;
// }

//____________________________________________________________________________..
int JetEnergyResolution::End(PHCompositeNode *topNode)
{
  std::cout << "JetEnergyResolution::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  outfile->cd();
  recoJetTree->Write();
  outfile->Write();
  outfile->Close();
  std::cout << "is this actually running??" << std::endl;
  delete outfile;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
// int JetEnergyResolution::Reset(PHCompositeNode *topNode)
// {
//  std::cout << "JetEnergyResolution::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
//   return Fun4AllReturnCodes::EVENT_OK;
// }

//____________________________________________________________________________..
// void JetEnergyResolution::Print(const std::string &what) const
// {
//   std::cout << "JetEnergyResolution::Print(const std::string &what) const Printing info for " << what << std::endl;
// }
