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
#include <TH3F.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>

#include "TCanvas.h"
#include "TH1F.h"
#include <TLegend.h>
#include <TMath.h>
#include <TRandom3.h>

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

using LOR = std::pair<TLorentzVector, TLorentzVector> ;

std::vector<LOR> select2(const TLorentzVector& gamma1, const TLorentzVector& gamma2, const TLorentzVector& gamma3);
std::vector<LOR> select(const TLorentzVector& gamma1, const TLorentzVector& gamma2, const TLorentzVector& gamma3);

auto exactlyOneHit = [](const Track& track) -> bool {
  return (track.fHits.size() == 1);
};

auto oneOrMoreHits = [](const Track& track) -> bool {
  return (track.fHits.size() >= 1);
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

auto nOrMoreHitsInEvent = [](const Event& event, int N) -> bool {
  int hits = 0;
  for (auto track : event.fTracks)
  {
    hits = hits + track.fHits.size();
  }
  return (hits >= N);
};

int numOfHitsInEvent(const Event& event)
{
  int hits = 0;
  for (const auto& track : event.fTracks)  {
    hits = hits + track.fHits.size();
  }
  return hits;
}

bool isScatteringInDetector1(const Hit& step)
{
  return step.fVolumeName == "crystal1";
}

bool isScatteringInDetector2(const Hit& step)
{
  return step.fVolumeName == "crystal2";
}

bool isScatteringInDetector3(const Hit& step)
{
  return step.fVolumeName == "crystal3";
}

auto scattering511 = [](const Hit& hit) -> bool {
  return (isEqual(hit.fEnergyBeforeProcess, 511));
};

auto scatteringprompt = [](const Hit& hit) -> bool {
  return (isEqual(hit.fEnergyBeforeProcess, 1157));
};

TGraph* purity_prompt = 0;
TGraph* efficiency_prompt = 0;
TGraph* ROC_prompt = 0;
TGraph* purity_511 = 0;
TGraph* efficiency_511 = 0;
TGraph* ROC_511 = 0;
TGraph* F_511_g = 0;
TGraph* F_511_g2 = 0;
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
  std::cout << "finalSelection  " << counter.g1g2 << std::endl;
}

void analyse(const std::string& inFile)
{

  TFile testOut("testOutSimple.root", "RECREATE");
  TH1F histMultiplicity("histMultiplicity", "histMultiplicity", 20, 0, 20);
  TH1F histMultiplicity2("histMultiplicity2", "histMultiplicity2", 20, 0, 20);
  TH1F hgammaenergyMax("hgammaenergyMax", "hgammaenergyMax", 500, 0, 1200);
  TH1F hgammapromptenergyMax("hgammapromptenergyMax", "hgammapromptenergyMax", 500, 0, 1200);
  TH1F hgammaenergyMid("hgammaenergyMid", "hgammaenergyMid", 500, 0, 1200);
  TH1F hgammapromptenergyMid("hgammapromptenergyMid", "hgammapromptenergyMid", 500, 0, 1200);
  TH1F hgammaenergyMin("hgammaenergyMin", "hgammaenergyMin", 500, 0, 1200);
  TH1F hgammapromptenergyMin("hgammapromptenergyMin", "hgammapromptenergyMin", 500, 0, 1200);
  TH1F hgamma1X("hgamma1X", "h1X", 500, -1200, 1200);
  TH1F hgamma2X("hgamma2X", "hgamma2X", 500, -1200, 1200);
  TH1F hpromptX("hpromptX", "hpromptX", 500, -1200, 1200);
  TH1F hgamma1prompt("hgamma1prompt", "hpromptZ", 500, 0, 600);
  TH1F hgamma2prompt("hgamma2prompt", "hgamma2Y", 500, 0, 600);
  TH1F hgamma1gamma2("hgamma1gamma2", "hgamma2Z", 500, 0, 600);
  TH3F hXYZprompt1("hXYZprompt1", "pr1", 20, -4, 4, 20, -4, 4, 20, 0, 20);
  TH3F hXYZprompt2("hXYZprompt2", "pr2", 20, -4, 4, 20, -4, 4, 20, 0, 20);
  TH3F hXYZgamma12("hXYZgamma12", "12", 20, -4, 4, 20, -4, 4, 20, 0, 20);
  TH1F h511All("h511All", "h", 500, 0, 1200);
  TH1F hpromptAll("hpromptAll", "h", 500, 0, 1200);
  TH1F h511Sygnal("h511Sygnal", "h", 500, 0, 1200); // 0 fantom && 1 detector
  TH1F hpromptSygnal("hpromptSygnal", "h", 500, 0, 1200); //0 fantom && 1 detector
  TH1F h511B1B2("h511B1B2", "h", 500, 0, 1200); //0 fantom && n detector + n fantom && n detector
  TH1F hpromptB1B2("hpromptB1B2", "h", 500, 0, 1200); //0 fantom && n detector + n fantom && n detector
  TH1F h511B1("h511B1", "h", 500, 0, 1200); //0 fantom && n detector
  TH1F hpromptB1("hpromptB1", "h", 500, 0, 1200); //0 fantom && n detector
  TH1F h511B2("h511B2", "h", 500, 0, 1200); //n fantom && n detector
  TH1F hpromptB2("hpromptB2", "h", 500, 0, 1200); //n fantom && n detector
  TH1F h511SygnalB1B2("h511SygnalB1B2", "h detector", 500, 0, 1200); //all types fantom && all types detector
  TH1F hpromptSygnalB1B2("hpromptSygnalB1B2", "h", 500, 0, 1200); //all types fantom && all types detector
  TH1F hostalos("hostalos", "h", 500, 0, 1200);
  TH1F hWTF("hWTF", "h", 500, 0, 1200);

  TLorentzVector gammaPrompt;
  TLorentzVector gamma1;
  TLorentzVector gamma2;
  TFile file(inFile.c_str(), "READ");
  TTreeReader reader("Tree", &file);
  TTreeReaderValue<Event> event(reader, "Event");

  std::vector <TLorentzVector> gammaPromptPos;
  std::vector <TLorentzVector> gamma511Pos1;
  std::vector <TLorentzVector> gamma511Pos2;

  //struct Triples{
  //std::vector <TLorentzVector> trackPrompt;
  //std::vector <TLorentzVector> track511_1;
  //std::vector <TLorentzVector> track511_2;
  //};
  //std::vector<Triples> registeredEvents;
  //registeredEvents.reserve(10000);

  while (reader.Next()) {
    histMultiplicity.Fill(event->fTracks.size());
    histMultiplicity2.Fill(numOfHitsInEvent(*event));
    if (numOfHitsInEvent(*event) < 3) {
      continue;
    }

    //Triples triples;

    if (exactlyThreeHitsInEvent(*event)) {
      for (const auto& track : event->fTracks) {

        double emissionEnergy = track.fEmissionEnergy;
        auto& rozp = track.fTrackID;
        auto& steps = track.fHits;
        for (auto i = 0u; i < steps.size(); i++) {
          auto& hit = steps[i];
          if (rozp == 2 && isEqual(emissionEnergy, 511)) {
            gamma1 = TLorentzVector(hit.fHitPosition, hit.fEnergyDeposition);
            gamma511Pos1.push_back(gamma1);
            //triples.track511_1.push_back(gamma1);
          } else {
            if (rozp == 3 && isEqual(emissionEnergy, 511)) {
              gamma2 = TLorentzVector(hit.fHitPosition, hit.fEnergyDeposition);
              gamma511Pos2.push_back(gamma2);
              //triples.track511_2.push_back(gamma2);
            } else {
              if (rozp == 1 && isEqual(emissionEnergy, 1157)) {
                gammaPrompt = TLorentzVector(hit.fHitPosition, hit.fEnergyDeposition);
                gammaPromptPos.push_back(gammaPrompt);
                //triples.trackPrompt.push_back(gammaPrompt);

              } else {
                std::cerr << "This should never happen" << std::endl;
                assert(1 == 0);
              }
            }
          }
        }
      }
    }
    //registeredEvents.push_back(triples);
  }


  int numOfEvents = gammaPromptPos.size();
  double distgamma1prompt[numOfEvents], distgamma2prompt[numOfEvents], distgamma1gamma2[numOfEvents];
  Counter counter;

  std::vector<LOR> result;
//int A=0;
  for (Int_t i = 0; i < numOfEvents ; i++) {

    LOR trueLOR = {gamma511Pos1[i], gamma511Pos2[i]};
    result = select2(gamma511Pos1[i], gamma511Pos2[i], gammaPromptPos[i]);
    // if

    distgamma1prompt[i] = calculateDistance3D(gammaPromptPos[i].X(), gammaPromptPos[i].Y(), gamma511Pos2[i].X(), gamma511Pos2[i].Y(), gammaPromptPos[i].Z(), gamma511Pos2[i].Z());
    distgamma2prompt[i] = calculateDistance3D(gammaPromptPos[i].X(), gammaPromptPos[i].Y(), gamma511Pos1[i].X(), gamma511Pos1[i].Y(), gammaPromptPos[i].Z(), gamma511Pos1[i].Z());
    distgamma1gamma2[i] = calculateDistance3D(gamma511Pos1[i].X(), gamma511Pos1[i].Y(), gamma511Pos2[i].X(), gamma511Pos2[i].Y(), gamma511Pos1[i].Z(), gamma511Pos2[i].Z());
    hgamma1prompt.Fill(distgamma1prompt[i]);
    hgamma2prompt.Fill(distgamma2prompt[i]);
    hgamma1gamma2.Fill(distgamma1gamma2[i]);
    hXYZprompt1.Fill(gammaPromptPos[i].X(), gammaPromptPos[i].Y(), gammaPromptPos[i].Z());

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
      hgammaenergyMax.Fill(gamma511Pos1[i].E());
      hgammapromptenergyMax.Fill(gammaPromptPos[i].E());
    } else if (isEqual(minVal, distgamma1prompt[i])) {
      hgammaenergyMin.Fill(gamma511Pos1[i].E());
      hgammapromptenergyMin.Fill(gammaPromptPos[i].E());
    } else {
      hgammaenergyMid.Fill(gamma511Pos1[i].E());
      hgammapromptenergyMid.Fill(gammaPromptPos[i].E());
    }

    if (isEqual(maxVal, distgamma2prompt[i])) {
      counter.pg2++;
      hgammaenergyMax.Fill(gamma511Pos2[i].E());
      hgammapromptenergyMax.Fill(gammaPromptPos[i].E());
    } else if (isEqual(minVal, distgamma2prompt[i])) {
      hgammaenergyMin.Fill(gamma511Pos2[i].E());
      hgammapromptenergyMin.Fill(gammaPromptPos[i].E());
    } else {
      hgammaenergyMid.Fill(gamma511Pos2[i].E());
      hgammapromptenergyMid.Fill(gammaPromptPos[i].E());
    }

    if (isEqual(maxVal, distgamma1gamma2[i])) {
      counter.g1g2++;
      hgammaenergyMax.Fill(gamma511Pos1[i].E());
      hgammaenergyMax.Fill(gamma511Pos2[i].E());
    } else if (isEqual(minVal, distgamma1gamma2[i])) {
      hgammaenergyMin.Fill(gamma511Pos1[i].E());
      hgammaenergyMin.Fill(gamma511Pos2[i].E());
    } else {
      hgammaenergyMid.Fill(gamma511Pos1[i].E());
      hgammaenergyMid.Fill(gamma511Pos2[i].E());
    }

    //koniec pętli for
  }
  double percent = counter.numericgamma * 100. / numOfEvents;

  printResults(numOfEvents, counter, percent);

  testOut.cd();
  histMultiplicity.Write();
  histMultiplicity2.Write();
  hgammaenergyMax.Write();
  hgammapromptenergyMax.Write();
  hgammaenergyMid.Write();
  hgammapromptenergyMid.Write();
  hgammaenergyMin.Write();
  hgammapromptenergyMin.Write();
  hgamma1X.Write();
  hgamma2X.Write();
  hpromptX.Write();
  hXYZgamma12.Write();
  hXYZprompt2.Write();
  hXYZprompt1.Write();
  hgamma1gamma2.Write();
  hgamma2prompt.Write();
  hgamma1prompt.Write();
  testOut.Close();
}



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

bool isEqualLor(const LOR& lor1, const LOR& lor2)
{
  return ((lor1.first == lor2.first)  && (lor1.second == lor2.second)) || ((lor1.first == lor2.second)  && (lor1.second == lor2.first));
}

// 3 hity N trojek
// select
// A - wybral dobrego LORa(511 *2)
// B - wybrala zle
// C - odrzucila trzy LORY
// D - obie zaakceptowal 2 LORY
// N =A+B+C


// A/N; B/N;

/// }

std::vector<LOR> select(const TLorentzVector& gamma1, const TLorentzVector& gamma2, const TLorentzVector& gamma3)
{
  LOR lor1 = {gamma1, gamma2};
  LOR lor2 = {gamma1, gamma3};
  LOR lor3 = {gamma2, gamma3};
  std::vector<LOR> all = {lor1, lor2, lor3};
  std::sort(all.begin(), all.end(), [](LOR a, LOR b)->bool {return calculateDistance(a) < calculateDistance(b);});
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

std::vector<LOR> select2(const TLorentzVector& gamma1, const TLorentzVector& gamma2, const TLorentzVector& gamma3)
{
  std::vector<LOR> results;
  int A = 0;
  int B = 0;
  int C = 0;
  int D = 0;
  int A_res = 0;
  int B_res = 0;
  int C_res = 0;
  int D_res = 0;
  std::vector <TLorentzVector> gamma12;
  gamma12.push_back(gamma1);
  int numOfEvents = gamma12.size();
  for (Int_t i = 0; i < numOfEvents ; i++) {
    LOR lorTrue;
    auto result = select(gamma1, gamma2, gamma3);
    if (result.empty()) {
      C++;
    } else {
      if (result.size() == 2) {
        D++;
      } else {
        auto lor = result[0];
        if (isEqualLor(lorTrue, lor)) {
          A++;
        } else {
          B++;
        }
      }
    }
  }
  A_res = A / numOfEvents;
  B_res = B / numOfEvents;
  C_res = C / numOfEvents;
  D_res = D / numOfEvents;
  return results;
}

void runTests()
{
  assert(isEqual(calculateDistance2D(1, sqrt(3), 1, -sqrt(3)), 1));
}

int main(int argc, char* argv[])
{
  runTests();
  using namespace tree_transformation;
  if ((argc != 3) && (argc != 4)) {
    std::cerr << "Invalid number of variables." << std::endl;
    std::cerr << "usage 1: ./GAR inputFile.root outputFile.root  " << std::endl;
    std::cerr << "usage 2: ./GAR inputFile.root outputFile.root  10000" << std::endl;
  } else {
    std::string in_file_name(argv[1]);
    std::string out_file_name(argv[2]);
    int  maxEvents = -1;
    if (argc >= 4) {
      maxEvents = std::stoi((argv[3]));
    }

    transformToEventTree(in_file_name, out_file_name, maxEvents);
    analyse(out_file_name);
    //std::cout << "finalSelection  " << finalSelection << std::endl;
  }
  return 0;
}
