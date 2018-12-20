
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

void addEntryToEvent(const GlobalActorReader& gar, Event* outEvent);
void clearEvent(Event* outEvent);
void transformToEventTree(const std::string& inFileName,
                          const std::string& outFileName);

void transformToEventTree(const std::string& inFileName,
                          const std::string& outFileName)
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
  } catch (const std::logic_error& e) {
    std::cerr << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Udefined exception" << std::endl;
  }
  fileOut.cd();
  assert(tree);
  fileOut.Write();
}

void simpleExampleHowToWorkWithEventTree(const std::string& inFile)
{

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
    for (const auto& track : event->fTracks) {
      auto& steps = track.fTrackInteractions;
      if (steps.size() > 1) {
        int counterAboveCut = 0;
        for (auto i = 0u; i < steps.size(); i++) {
          auto& step = steps[i];
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


random_device rd;
mt19937 gen(rd());

Double_t sigmaE(Double_t E, Double_t coeff = 0.0444)
{
  return coeff / TMath::Sqrt(E) * E;
}

double r_norm(double mean, double sigmaE)
{
  normal_distribution<double> d(mean, sigmaE);
  return d(gen);
}

bool isEqual(double x, double y, double epsilon = 10e-9)
{
  return std::abs(x - y) < epsilon;
}


double smearEnergy(double energy)
{
  return r_norm(energy, 1000. * sigmaE((energy) * 1. / 1000.));
}

double calculateDistance(double x1, double y1, double x2, double y2)
{
  double distance = 0;
  distance = abs(x2 * y1 - y2 * x1) / sqrt(pow((y2 - y1), 2) + pow((x2 - x1), 2));
  return distance;
}

auto exactlyOneHit = [](const Track& track) -> bool {
  return (track.fTrackInteractions.size() == 1);
};

auto exactlyTwoHitsInEvent = [](const Event& event) -> bool {
  return (event.fTracks.size() == 2) &&
  std::all_of(event.fTracks.begin(), event.fTracks.end(),
  exactlyOneHit);
};

auto exactlyThreeHitsInEvent = [](const Event& event) -> bool {
  return (event.fTracks.size() == 3) &&
  std::all_of(event.fTracks.begin(), event.fTracks.end(),
  exactlyOneHit);
};
TGraph* gamma1gamma2_wykres = 0;

void simpleExampleHowToWorkWithEventTree2(const std::string& inFile)
{

  TFile testOut("testOutSimple.root", "RECREATE");
  TH1F hgamma1X("hgamma1X", "h1X", 500, -1200, 1200);
  TH1F hgamma2X("hgamma2X", "hgamma2X", 500, -1200, 1200);
  TH1F hpromptX("hpromptX", "hpromptX", 500, -1200, 1200);
  TH1F hgamma1Y("hgamma1Y", "hgamma1Y", 500, -1200, 1200);
  TH1F hgamma1Z("hgamma1Z", "hgamma1Z", 500, -1200, 1200);
  TH1F hpromptY("hpromptY", "hpromptY", 500, -1200, 1200);
  TH1F hpromptZ("hpromptZ", "hpromptZ", 500, -1200, 1200);
  TH1F hgamma2Y("hgamma2Y", "hgamma2Y", 500, -1200, 1200);
  TH1F hgamma2Z("hgamma2Z", "hgamma2Z", 500, -1200, 1200);
  TH1F hgamma1prompt("hgamma1prompt", "hpromptZ", 500, -1200, 1200);
  TH1F hgamma2prompt("hgamma2prompt", "hgamma2Y", 500, -1200, 1200);
  TH1F hgamma1gamma2("hgamma1gamma2", "hgamma2Z", 500, -1200, 1200);


  TLorentzVector gammaPrompt;
  TLorentzVector gamma1;
  TLorentzVector gamma2;
  TFile file(inFile.c_str(), "READ");
  TTreeReader reader("Tree", &file);
  TTreeReaderValue<Event> event(reader, "Event");
  double cut = 200;

  std::vector <TLorentzVector> gammaPromptPos;
  std::vector <TLorentzVector> gamma511Pos1;
  std::vector <TLorentzVector> gamma511Pos2;
  assert(isEqual(calculateDistance(1, sqrt(3), 1, -sqrt(3)), 1));
  while (reader.Next()) {
    for (const auto& track : event->fTracks) {

      double gamma1X, gamma2X;

      auto& rozp = track.fTrackID;
      double emissionEnergy = track.fEmissionEnergy;
      auto& steps = track.fTrackInteractions;

      for (auto i = 0u; i < steps.size(); i++) {
        auto& hit = steps[i];

        if (exactlyThreeHitsInEvent(*event)) {
          if (rozp == 2 && isEqual(emissionEnergy, 511)) {

            gamma1 = TLorentzVector(hit.fHitPosition, hit.fEnergyDeposition);
            gamma1X = gamma1.Vect().X();
            gamma511Pos1.push_back(gamma1);

            hgamma1X.Fill(gamma1X);
          }
          if (rozp == 3 && isEqual(emissionEnergy, 511)) {
            gamma2 = TLorentzVector(hit.fHitPosition, hit.fEnergyDeposition);
            gamma2X = gamma2.Vect().X();
            gamma511Pos2.push_back(gamma2);
            hgamma2X.Fill(gamma2X);
          }

          if (rozp == 1 && isEqual(emissionEnergy, 1157)) {
            gammaPrompt = TLorentzVector(hit.fHitPosition, hit.fEnergyDeposition);
            gammaPromptPos.push_back(gammaPrompt);
            hpromptX.Fill(gammaPrompt.X());

          }
        }

      }
    }


  }

  int eventstep = gammaPromptPos.size();
  //assert(eventstep==110);
  //Double_t number[eventstep];
  //Double_t gamma1gamma2 [eventstep];
  for (Int_t i = 0; i < eventstep ; i++) {
    //gamma1gamma2[i] = calculateDistance(gamma511Pos1[i].X(), gamma511Pos1[i].Y(), gamma511Pos2[i].X(), gamma511Pos2[i].Y());
    //number[i] = i;
    hgamma1prompt.Fill(calculateDistance(gammaPromptPos[i].X(), gammaPromptPos[i].Y(), gamma511Pos2[i].X(), gamma511Pos2[i].Y()));
    hgamma2prompt.Fill(calculateDistance(gammaPromptPos[i].X(), gammaPromptPos[i].Y(), gamma511Pos1[i].X(), gamma511Pos1[i].Y()));
    hgamma1gamma2.Fill(calculateDistance(gamma511Pos1[i].X(), gamma511Pos1[i].Y(), gamma511Pos2[i].X(), gamma511Pos2[i].Y()));

  }
  //  gamma1gamma2_wykres = new TGraph(eventstep, gamma1gamma2, number);

  testOut.cd();
  hgamma1X.Write();
  hgamma2X.Write();
  hpromptX.Write();
  hgamma1gamma2.Write();
  hgamma2prompt.Write();
  hgamma1prompt.Write();
  //gamma1gamma2_wykres->Write();
  testOut.Close();
}

//void hardcoreExampleHowToWorkWithEventTree(const std::string& inFile)
//{

  //auto exactlyOneHit = [](const Track & track) -> bool {
    //return (track.fTrackInteractions.size() == 1);
  //};
  //auto exactlyTwoHitsInEvent = [&exactlyOneHit](const Event & event) -> bool {
    //return (event.fTracks.size() == 2) &&
    //std::all_of(event.fTracks.begin(), event.fTracks.end(),
    //exactlyOneHit);
  //};

  //auto scattering511 = [](const TrackInteraction & hit) -> bool {
    //return (hit.fEnergyBeforeProcess == 511);
  //};
  //auto onlyfirstScattering511 = [scattering511](const Track & track) -> bool {
    //return (track.fTrackInteractions.size() == 1) &&
    //scattering511(track.fTrackInteractions.front());
  //};

  //double energyCut = 200;
  //auto aboveEnergy = [energyCut](const TrackInteraction & hit) -> bool {
    //return (hit.fEnergyDeposition > energyCut);
  //};
  //auto atLeastOneAbove = [energyCut, &aboveEnergy](const Track & track) -> bool {
    //return std::any_of(track.fTrackInteractions.begin(),
    //track.fTrackInteractions.end(), aboveEnergy);
  //};
  //auto allHitsAboveSomeEnergyCut =
  //[energyCut, &atLeastOneAbove](const Event & event) -> bool {
    //return std::all_of(event.fTracks.begin(), event.fTracks.end(),
    //atLeastOneAbove);
  //};

  //TFile testOut("testOut.root", "RECREATE");
  //TH1F h("h", "h", 500, 0, 1200);
  //TFile file(inFile.c_str(), "READ");
  //TTreeReader reader("Tree", &file);
  //TTreeReaderValue<Event> event(reader, "Event");
  ////auto hitPosition = fHitPosition.Vect().X();
  //while (reader.Next()) {
    ////auto &steps = hit.fTrackInteractions;

    //if (exactlyTwoHitsInEvent(*event) && allHitsAboveSomeEnergyCut(*event)) {

      //h.Fill(event->fTracks[0].fTrackInteractions[0].fEnergyDeposition);
      //h.Fill(event->fTracks[1].fTrackInteractions[0].fEnergyDeposition);
    //}
  //}
  //testOut.cd();
  //h.Write();

  //TCanvas c8("c", "c", 2000, 1200);
  //h.SetTitle("Background 511: 1-n fantom && n detector");
  //h.GetXaxis()->SetTitle("Energy [keV]");
  //h.GetYaxis()->SetTitle("Events");
  //h.SetLineColor(kBlack);
  //h.Draw();
  //c8.SaveAs("b2_511.png");
  //testOut.Close();
//}

void addEntryToEvent(const GlobalActorReader& gar, Event* outEvent)
{
  assert(outEvent);
  outEvent->fEventID = gar.GetEventID();

  TrackInteraction trkStep;
  trkStep.fHitPosition = gar.GetProcessPosition();
  trkStep.fEnergyDeposition = gar.GetEnergyLossDuringProcess();
  trkStep.fEnergyBeforeProcess = gar.GetEnergyBeforeProcess();
  trkStep.fVolumeName = gar.GetVolumeName();

  int currentTrackID = gar.GetTrackID();
  if (!outEvent->fTracks.empty()) {
    auto& lastTrack = outEvent->fTracks.back();
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

void clearEvent(Event* outEvent)
{
  assert(outEvent);
  outEvent->fEventID = -1;
  outEvent->fTracks.clear();
}

int main(int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Invalid number of variables." << std::endl;
    std::cerr << "usage: ./GAR inputFile.root outputFile.root" << std::endl;
  } else {
    std::string in_file_name(argv[1]);
    std::string out_file_name(argv[2]);
    transformToEventTree(in_file_name, out_file_name);
    //hardcoreExampleHowToWorkWithEventTree(out_file_name);
    //simpleExampleHowToWorkWithEventTree(out_file_name);
    simpleExampleHowToWorkWithEventTree2(out_file_name);
  }
  return 0;
}
