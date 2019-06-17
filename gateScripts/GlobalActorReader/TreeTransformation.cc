/**
 *  @copyright Copyright 2018 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  @file TreeTransformation.cc
 */

#include "TreeTransformation.h"
#include <cassert>

namespace tree_transformation
{

void transformToEventTree(const std::string& inFileName,
                          const std::string& outFileName,
                          int maxNumEvents)
{
  TFile fileOut(outFileName.c_str(), "RECREATE");
  TTree* tree = new TTree("Tree", "Tree");
  Event* event = nullptr;
  tree->Branch("Event", &event, 16000, 99);
  try {
    event = new Event;
    GlobalActorReader gar;
    if (gar.LoadFile(inFileName.c_str())) {
      bool isNewEvent = false;
      bool isFirstEvent = false;
      auto previousID = event->fEventID;
      auto currentID = previousID;
      if (maxNumEvents >= 0) {
        int currSeq = 0;
        while (gar.Read() && (currSeq < maxNumEvents)) {
          currSeq++;
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
      } else {
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
      }
      if (event->fEventID > 0) {
        tree->Fill();
        clearEvent(event);
      }
    } else {
      std::cerr << "Loading file failed." << std::endl;
    }
  } catch (const std::logic_error& e) {
    std::cerr << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Udefined exception" << std::endl;
  }
  fileOut.cd();
  assert(tree);
  fileOut.Write();
}

void addEntryToEvent(const GlobalActorReader& gar, Event* outEvent)
{
  assert(outEvent);
  outEvent->fEventID = gar.GetEventID();

  Hit trkStep;
  trkStep.fHitPosition = gar.GetProcessPosition();
  trkStep.fEnergyDeposition = gar.GetEnergyLossDuringProcess();
  trkStep.fEnergyBeforeProcess = gar.GetEnergyBeforeProcess();
  trkStep.fVolumeName = gar.GetVolumeName();

  int currentTrackID = gar.GetTrackID();
  if (!outEvent->fTracks.empty()) {
    auto& lastTrack = outEvent->fTracks.back();
    if (lastTrack.fTrackID == currentTrackID) {
      lastTrack.fHits.push_back(trkStep);
    } else {
      Track trk;
      trk.fEmissionEnergy = gar.GetEmissionEnergyFromSource();
      trk.fTrackID = currentTrackID;
      trk.fHits.push_back(trkStep);
      outEvent->fTracks.push_back(trk);
    }
  } else {
    Track trk;
    trk.fEmissionEnergy = gar.GetEmissionEnergyFromSource();
    trk.fTrackID = currentTrackID;
    trk.fHits.push_back(trkStep);
    outEvent->fTracks.push_back(trk);
  }
}

void clearEvent(Event* outEvent)
{
  assert(outEvent);
  outEvent->fEventID = -1;
  outEvent->fTracks.clear();
}

}
