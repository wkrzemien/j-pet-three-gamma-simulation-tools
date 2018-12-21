//R__LOAD_LIBRARY(Event_h.so)
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

#include <TH1F.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>

#include "Event.h"
#include "GlobalActorReader.hh"
#include "TreeTransformation.h"
#include "HelperFunctions.h"

using namespace std;
using namespace helper_functions;


double calculateDistance(double x1, double y1, double x2, double y2)
{
  if (x1 == x2 && y1 == y2) return 0;
  double distance = 0;
  distance = abs(x2 * y1 - y2 * x1) / sqrt(pow((y2 - y1), 2) + pow((x2 - x1), 2));
  return distance;
}

auto exactlyOneHit = [](const Track& track) -> bool {
  return (track.fHits.size() == 1);
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

void analyse(const std::string& inFile)
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

  std::vector <TLorentzVector> gammaPromptPos;
  std::vector <TLorentzVector> gamma511Pos1;
  std::vector <TLorentzVector> gamma511Pos2;
  while (reader.Next()) {
    for (const auto& track : event->fTracks) {

      double gamma1X, gamma2X;

      auto& rozp = track.fTrackID;
      double emissionEnergy = track.fEmissionEnergy;
      auto& steps = track.fHits;

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


void runTests()
{
  assert(isEqual(calculateDistance(1, sqrt(3), 1, -sqrt(3)), 1));
}

int main(int argc, char* argv[])
{
  runTests();
  using namespace tree_transformation;
  if (argc != 3) {
    std::cerr << "Invalid number of variables." << std::endl;
    std::cerr << "usage: ./GAR inputFile.root outputFile.root" << std::endl;
  } else {
    std::string in_file_name(argv[1]);
    std::string out_file_name(argv[2]);
    transformToEventTree(in_file_name, out_file_name);
    analyse(out_file_name);
  }
  return 0;
}
