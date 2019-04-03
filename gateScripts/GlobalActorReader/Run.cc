//R__LOAD_LIBRARY(Event_h.so)
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <initializer_list>
#include <stdlib.h>
#include <stdio.h>

#include <TH1F.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>

#include "Event.h"
#include "GlobalActorReader.hh"
#include "TreeTransformation.h"
#include "HelperFunctions.h"



struct Counter {
  int numericprompt = 0;
  int numericgamma = 0;
  int numericgammamax = 0;
  int numericgammamin = 0;
  int numericgammamid = 0;
  int g1g2 = 0;
  int pg1 = 0;
  int pg2 = 0;
};

using namespace std;
using namespace helper_functions;




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

void printResults(int numOfEvents, const Counter& counter, double percent)
{
  std::cout << "Results for back-to-back: " << std::endl;
  std::cout << "Liczba gamma dla min odleglosci  " << counter.numericgammamin << std::endl;
  std::cout << "Liczba gamma dla max odleglosci  " << counter.numericgammamax << std::endl;
  std::cout << "Liczba gamma dla mid odleglosci  " << counter.numericgammamid << std::endl;
  std::cout << "Liczba zaakceptowanych prompt:  " << counter.numericprompt << std::endl;
  std::cout << "Liczba zaakceptowanych gamma jako prompt:  " << counter.numericgamma << std::endl;
  std::cout << "Percent of false " << percent << "%" << std::endl;
  std::cout << "Event steps  " << numOfEvents << std::endl;
  std::cout << "gamma1prompt  " << counter.pg1 << std::endl;
  std::cout << "gamma2prompt  " << counter.pg2 << std::endl;
  std::cout << "gamma1gamma2  " << counter.g1g2 << std::endl;
}

void analyse(const std::string& inFile)
{

  TFile testOut("testOutSimple.root", "RECREATE");
  TH1F histMultiplicyt("histMultiplicyt", "histMultiplicyt", 500, -1200, 1200);
  TH1F hgamma1energy("hgamma1energy", "h1energy", 500, 0, 1200);
  TH1F hgamma2energy("hgamma2energy", "h2energy", 500, 0, 1200);
  TH1F hgammapromptenergy("hgammapromptenergy", "hpromptenergy", 500, 0, 1200);
  TH1F hgamma1X("hgamma1X", "h1X", 500, -1200, 1200);
  TH1F hgamma2X("hgamma2X", "hgamma2X", 500, -1200, 1200);
  TH1F hpromptX("hpromptX", "hpromptX", 500, -1200, 1200);
  TH1F hgamma1Y("hgamma1Y", "hgamma1Y", 500, -1200, 1200);
  TH1F hgamma1Z("hgamma1Z", "hgamma1Z", 500, -1200, 1200);
  TH1F hpromptY("hpromptY", "hpromptY", 500, -1200, 1200);
  TH1F hpromptZ("hpromptZ", "hpromptZ", 500, -1200, 1200);
  TH1F hgamma2Y("hgamma2Y", "hgamma2Y", 500, -1200, 1200);
  TH1F hgamma2Z("hgamma2Z", "hgamma2Z", 500, -1200, 1200);
  TH1F hgamma1prompt("hgamma1prompt", "hpromptZ", 500, 0, 600);
  TH1F hgamma2prompt("hgamma2prompt", "hgamma2Y", 500, 0, 600);
  TH1F hgamma1gamma2("hgamma1gamma2", "hgamma2Z", 500, 0, 600);

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
        histMultiplicyt.Fill(event->fTracks.size());
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

  int numOfEvents = gammaPromptPos.size();
  double distgamma1prompt[numOfEvents], distgamma2prompt[numOfEvents], distgamma1gamma2[numOfEvents];
  Counter counter;

  //bool gamma1 = false;
  //assert(numOfEvents==110);
  //Double_t number[numOfEvents];
  //Double_t gamma1gamma2 [numOfEvents];
  for (Int_t i = 0; i < numOfEvents ; i++) {

    distgamma1prompt[i] = calculateDistance2D(gammaPromptPos[i].X(), gammaPromptPos[i].Y(), gamma511Pos2[i].X(), gamma511Pos2[i].Y());
    distgamma2prompt[i] = calculateDistance2D(gammaPromptPos[i].X(), gammaPromptPos[i].Y(), gamma511Pos1[i].X(), gamma511Pos1[i].Y());
    distgamma1gamma2[i] = calculateDistance2D(gamma511Pos1[i].X(), gamma511Pos1[i].Y(), gamma511Pos2[i].X(), gamma511Pos2[i].Y());
    hgamma1prompt.Fill(distgamma1prompt[i]);
    hgamma2prompt.Fill(distgamma2prompt[i]);
    hgamma1gamma2.Fill(distgamma1gamma2[i]);

    auto maxVal = std::max({distgamma1gamma2[i], distgamma1prompt[i], distgamma2prompt[i]  });
    auto minVal = std::min({distgamma1gamma2[i], distgamma1prompt[i], distgamma2prompt[i]  });
    if (isEqual(maxVal, distgamma1gamma2[i])) {
      //gammaPromptPos[i].E();
      counter.numericgammamax++;
    } else {
      if (isEqual(minVal, distgamma1gamma2[i])) {
        counter.numericgammamin++;

      } else {
        counter.numericgammamid++;
      }
    }

    if (isEqual(maxVal, distgamma1prompt[i]) || isEqual(maxVal, distgamma2prompt[i])) {
      counter.numericprompt++;
    } else {
      counter.numericgamma++;
    }

    if (isEqual(maxVal, distgamma1prompt[i])) {
      counter.pg1++;
      hgamma1energy.Fill(gamma511Pos1[i].E());
      hgammapromptenergy.Fill(gammaPromptPos[i].E());
    }
    if (isEqual(maxVal, distgamma2prompt[i])) {
      counter.pg2++;
      hgamma1energy.Fill(gamma511Pos2[i].E());
      hgammapromptenergy.Fill(gammaPromptPos[i].E());
    }
    if (isEqual(maxVal, distgamma1gamma2[i])) {
      counter.g1g2++;
      // hgamma1energy.Fill(gamma511Pos1[i].E());
      // hgamma2energy.Fill(gamma511Pos2[i].E());

    }

    //koniec pętli for
  }
  double percent = counter.numericgamma * 100. / numOfEvents;

  printResults(numOfEvents, counter, percent);


  testOut.cd();
  histMultiplicyt.Write();
  hgamma1energy.Write();
  hgamma2energy.Write();
  hgammapromptenergy.Write();
  hgamma1X.Write();
  hgamma2X.Write();
  hpromptX.Write();
  hgamma1gamma2.Write();
  hgamma2prompt.Write();
  hgamma1prompt.Write();
  testOut.Close();
}

using LOR = std::pair<TLorentzVector, TLorentzVector> ;

double calculateDistance(const LOR& lor)
{
  return calculateDistance3D(lor.first.X(), lor.first.Y(), lor.first.Z(), lor.second.X(), lor.second.Y(), lor.second.Z());
}

bool isInEnergyRange(double E)
{
  double Emin = 200;
  double Ecut = 338;
  return ((E >= Emin) && (E <= Ecut));
}

std::vector<LOR> select(const TLorentzVector& gamma1, const TLorentzVector& gamma2, const TLorentzVector& gamma3)
{
  LOR lor1 = {gamma1, gamma2};
  LOR lor2 = {gamma1, gamma3};
  LOR lor3 = {gamma2, gamma3};
  std::vector<LOR> all = {lor1, lor2, lor3};
  std::sort(all.begin(), all.end(), [](LOR a, LOR b)->bool {return calculateDistance(a) < calculateDistance(b);});
  // now all= {lor_min, lor_med, lor_max}
  // we reject lor_max
  std::vector<LOR> selectedAfterGeomCut = {all[0], all[1]};
  std::vector<LOR> finalSelection;


  for (auto& lor : selectedAfterGeomCut ) {
    double E1 = lor.first.Energy();
    double E2 = lor.second.Energy();
    if (isInEnergyRange(E1) && isInEnergyRange(E2)) {
      finalSelection.push_back(lor);
    }
  }
  return finalSelection;
}



void runTests()
{
  assert(isEqual(calculateDistance2D(1, sqrt(3), 1, -sqrt(3)), 1));
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
