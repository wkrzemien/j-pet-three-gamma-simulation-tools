//R__LOAD_LIBRARY(Event_h.so)
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
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
  int distgamma1prompt[eventstep], distgamma2prompt[eventstep], distgamma1gamma2[eventstep];
  int numeric = 0;
  int numericprompt = 0;
  int numericgamma = 0;
  int numericgammamax = 0;
  int numericgammamin = 0;
  int numericgammamid = 0;
  int g1g2 = 0;
  int pg1 = 0;
  int pg2 = 0;
  double percent =0;

  //bool gamma1 = false;
  //assert(eventstep==110);
  //Double_t number[eventstep];
  //Double_t gamma1gamma2 [eventstep];
  for (Int_t i = 0; i < eventstep ; i++) {

    distgamma1prompt[i]=calculateDistance(gammaPromptPos[i].X(), gammaPromptPos[i].Y(), gamma511Pos2[i].X(), gamma511Pos2[i].Y());
    distgamma2prompt[i]=calculateDistance(gammaPromptPos[i].X(), gammaPromptPos[i].Y(), gamma511Pos1[i].X(), gamma511Pos1[i].Y());
    distgamma1gamma2[i]=calculateDistance(gamma511Pos1[i].X(), gamma511Pos1[i].Y(), gamma511Pos2[i].X(), gamma511Pos2[i].Y());
    hgamma1prompt.Fill(distgamma1prompt[i]);
    hgamma2prompt.Fill(distgamma2prompt[i]);
    hgamma1gamma2.Fill(distgamma1gamma2[i]);

    auto maxVal= std::max({distgamma1gamma2[i],distgamma1prompt[i],distgamma2prompt[i]  });
    auto minVal= std::min({distgamma1gamma2[i],distgamma1prompt[i],distgamma2prompt[i]  });
    if (isEqual(maxVal, distgamma1gamma2[i])) {
    	numericgammamax++;
    } else {
    	if (isEqual(minVal, distgamma1gamma2[i])) {
    		numericgammamin++;

    	} else{
    		numericgammamid++;
    	}
    }

    if (isEqual(maxVal, distgamma1prompt[i])||isEqual(maxVal, distgamma2prompt[i])) {
      numericprompt++;
    }
    else
    {
      numericgamma++;
    }

    /*if (isEqual(maxVal, distgamma1prompt[i])) {
      pg1++;
    }
    if (isEqual(maxVal, distgamma2prompt[i])) {
      pg2++;
    }
    if (isEqual(maxVal, distgamma1gamma2[i])) {
      g1g2++;
    }*/

    //koniec pętli for
  }
  percent=numericgamma*100 / eventstep;

    std::cout << "Results for back-to-back: " << std::endl;
    std::cout << "Liczba gamma dla min odleglosci  " << numericgammamin <<std::endl;
    std::cout << "Liczba gamma dla max odleglosci  " << numericgammamax <<std::endl;
    std::cout << "Liczba gamma dla mid odleglosci  " << numericgammamid <<std::endl;
    std::cout << "Liczba zaakceptowanych prompt:  " << numericprompt <<std::endl;
    std::cout << "Liczba zaakceptowanych gamma jako prompt:  " << numericgamma <<std::endl;
    std::cout << "Percent of false " << percent << "%" << std::endl;
    std::cout << "Event steps  " << eventstep <<std::endl;
    std::cout << "gamma1prompt  " << pg1 <<std::endl;
    std::cout << "gamma2prompt  " << pg2 <<std::endl;
    std::cout << "gamma1gamma2  " << g1g2 <<std::endl;


  testOut.cd();
  hgamma1X.Write();
  hgamma2X.Write();
  hpromptX.Write();
  hgamma1gamma2.Write();
  hgamma2prompt.Write();
  hgamma1prompt.Write();
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