//R__LOAD_LIBRARY(Event_h.so)
#include "TCanvas.h"
#include "TH1F.h"
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TTreeReader.h>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Event.h"
#include "GlobalActorReader.hh"

using namespace std;

void addEntryToEvent(const GlobalActorReader &gar, Event *outEvent);
void clearEvent(Event *outEvent);
void transformToEventTree(const std::string &inFileName,
                          const std::string &outFileName);

void transformToEventTree(const std::string &inFileName,
                          const std::string &outFileName) {
  TFile fileOut(outFileName.c_str(), "RECREATE");
  TTree *tree = new TTree("Tree", "Tree");
  Event *event = nullptr;
  tree->Branch("Event", &event, 16000, 99);
  try {
    event = new Event;
    GlobalActorReader gar;
    if (gar.LoadFile(inFileName.c_str())) {
      bool isNewEvent = false;
      bool isFirstEvent = false;
      auto previousID = event->fEventID;
      auto currentID = previousID;
      while (gar.Read()) {
        currentID = gar.GetEventID();
        isFirstEvent = (previousID < 0) && (currentID > 0);
        isNewEvent = currentID != previousID;

        if (isFirstEvent) {
          addEntryToEvent(gar, event);
        } else {
          if (isNewEvent) {
            tree->Fill();
            clearEvent(event);
          }
          addEntryToEvent(gar, event);
        }
        previousID = currentID;
      }
      if (event->fEventID > 0) {
        tree->Fill();
        clearEvent(event);
      }
    } else {
      std::cerr << "Loading file failed." << std::endl;
    }
  } catch (const std::logic_error &e) {
    std::cerr << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Udefined exception" << std::endl;
  }
  fileOut.cd();
  assert(tree);
  fileOut.Write();
}

void simpleExampleHowToWorkWithEventTree(const std::string &inFile) {

  TFile testOut("testOutSimple.root", "RECREATE");
  TH1F h("h", "h", 500, 0, 1200);
  TFile file(inFile.c_str(), "READ");
  TTreeReader reader("Tree", &file);
  TTreeReaderValue<Event> event(reader, "Event");
  /// Let's assume I want to plot the deposited energy of second scattering from
  /// photons which are scattered
  /// at least twice above the given energy threshold.
  double cut = 200;
  while (reader.Next()) {
    for (const auto &track : event->fTracks) {
      auto &steps = track.fTrackInteractions;
      if (steps.size() > 1) {
        int counterAboveCut = 0;
        for (auto i = 0u; i < steps.size(); i++) {
          auto &step = steps[i];
          if (step.fEnergyDeposition > cut)
            counterAboveCut++;
          if (counterAboveCut == 2) {
            h.Fill(step.fEnergyDeposition);
            break;
          }
        }
      }
    }
  }
  testOut.cd();
  h.Write();
  testOut.Close();
}

bool isScatteringInPhantom(const TrackInteraction& step) {
  return step.fVolumeName == "detector1";
}

bool isScattering511(const TrackInteraction &hit) {
    return (hit.fEnergyBeforeProcess == 511);
};

void simpleExampleHowToWorkWithEventTree2(const std::string &inFile) {

  TFile testOut("testOutSimple.root", "RECREATE");
  TH1F h("h", "h", 500, 0, 1200);
  TFile file(inFile.c_str(), "READ");
  TTreeReader reader("Tree", &file);
  TTreeReaderValue<Event> event(reader, "Event");
  double cut = 200;
  while (reader.Next()) {
    for (const auto &track : event->fTracks) {

      bool wasInPhantom = false;
      bool isInPhantom = false;
      bool isInDetector = false;
      bool is511 = false;

      auto &steps = track.fTrackInteractions;
      int counterAboveCut = 0;
      for (auto i = 0u; i < steps.size(); i++) {
        auto &hit = steps[i];
        is511 =isScattering511(hit);
        isInPhantom = isScatteringInPhantom(hit);
        isInDetector = !isInPhantom;
        if(is511){
          std::cout << "It is 511! " << std::endl;
        }
        if(isInDetector) {
          if(wasInPhantom) {
          h.Fill(hit.fEnergyDeposition);
        } else {
          ///second histogram 
          }
        }

        if (isInPhantom) {
          wasInPhantom =  true;
        }
      }
    }
  }
  testOut.cd();
  h.Write();
  testOut.Close();
}

void hardcoreExampleHowToWorkWithEventTree(const std::string &inFile) {

  auto exactlyOneHit = [](const Track &track) -> bool {
    return (track.fTrackInteractions.size() == 1);
  };
  auto exactlyTwoHitsInEvent = [&exactlyOneHit](const Event &event) -> bool {
    return (event.fTracks.size() == 2) &&
           std::all_of(event.fTracks.begin(), event.fTracks.end(),
                       exactlyOneHit);
  };

  auto scattering511 = [](const TrackInteraction &hit) -> bool {
    return (hit.fEnergyBeforeProcess == 511);
  };
  auto onlyfirstScattering511 = [scattering511](const Track &track) -> bool {
    return (track.fTrackInteractions.size() == 1) &&
           scattering511(track.fTrackInteractions.front());
  };

  double energyCut = 200;
  auto aboveEnergy = [energyCut](const TrackInteraction &hit) -> bool {
    return (hit.fEnergyDeposition > energyCut);
  };
  auto atLeastOneAbove = [energyCut, &aboveEnergy](const Track &track) -> bool {
    return std::any_of(track.fTrackInteractions.begin(),
                       track.fTrackInteractions.end(), aboveEnergy);
  };
  auto allHitsAboveSomeEnergyCut =
      [energyCut, &atLeastOneAbove](const Event &event) -> bool {
    return std::all_of(event.fTracks.begin(), event.fTracks.end(),
                       atLeastOneAbove);
  };

  TFile testOut("testOut.root", "RECREATE");
  TH1F h("h", "h", 500, 0, 1200);
  TFile file(inFile.c_str(), "READ");
  TTreeReader reader("Tree", &file);
  TTreeReaderValue<Event> event(reader, "Event");
  while (reader.Next()) {
    if (exactlyTwoHitsInEvent(*event) && allHitsAboveSomeEnergyCut(*event)) {
      h.Fill(event->fTracks[0].fTrackInteractions[0].fEnergyDeposition);
      h.Fill(event->fTracks[1].fTrackInteractions[0].fEnergyDeposition);
    }
  }
  testOut.cd();
  h.Write();
  testOut.Close();
}

void addEntryToEvent(const GlobalActorReader &gar, Event *outEvent) {
  assert(outEvent);
  outEvent->fEventID = gar.GetEventID();

  TrackInteraction trkStep;
  trkStep.fHitPosition = gar.GetProcessPosition();
  trkStep.fEnergyDeposition = gar.GetEnergyLossDuringProcess();
  trkStep.fEnergyBeforeProcess = gar.GetEnergyBeforeProcess();
  trkStep.fVolumeName = gar.GetVolumeName();

  int currentTrackID = gar.GetTrackID();
  if (!outEvent->fTracks.empty()) {
    auto &lastTrack = outEvent->fTracks.back();
    if (lastTrack.fTrackID == currentTrackID) {
      lastTrack.fTrackInteractions.push_back(trkStep);
    } else {
      Track trk;
      trk.fEmissionEnergy = gar.GetEmissionEnergyFromSource();
      trk.fTrackID = currentTrackID;
      trk.fTrackInteractions.push_back(trkStep);
      outEvent->fTracks.push_back(trk);
    }
  } else {
    Track trk;
    trk.fEmissionEnergy = gar.GetEmissionEnergyFromSource();
    trk.fTrackID = currentTrackID;
    trk.fTrackInteractions.push_back(trkStep);
    outEvent->fTracks.push_back(trk);
  }
}

void clearEvent(Event *outEvent) {
  assert(outEvent);
  outEvent->fEventID = -1;
  outEvent->fTracks.clear();
}

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cerr << "Invalid number of variables." << std::endl;
    std::cerr << "usage: ./GAR inputFile.root outputFile.root" << std::endl;
  } else {
    std::string in_file_name(argv[1]);
    std::string out_file_name(argv[2]);
    transformToEventTree(in_file_name, out_file_name);
    //hardcoreExampleHowToWorkWithEventTree(out_file_name);
    simpleExampleHowToWorkWithEventTree(out_file_name);
    //simpleExampleHowToWorkWithEventTree2(out_file_name);
  }
  return 0;
}
