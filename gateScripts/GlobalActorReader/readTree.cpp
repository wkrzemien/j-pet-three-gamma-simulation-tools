#include <TROOT.h>
#include <TChain.h>
/// This is some completely awful hack to make RDF work
/// on CIS cluster
#ifdef R__HAS_VDT
#undef R__HAS_VDT
#endif
#include "ROOT/RDataFrame.hxx"
#ifdef R__HAS_VDT
#undef R__HAS_VDT
#endif

#include "Event.h"
#include <iostream>
#include <string>
#include <TCanvas.h>
#include "HelperFunctions.h"

using namespace std;
using namespace helper_functions;
using TH1DPtr = ROOT::RDF::RResultPtr<TH1D>;

std::string kFirstLayerName = "crystal1";
std::string kSecondLayerName = "crystal2";
std::string kThirdLayerName = "crystal3";
std::string kPhantomName = "NEMA_IQ";


void readTree()
{
  //const char* inFile = "NEMA_git.root"; 
  //const char* inFile = "NEMA_git2.root"; 
  const char* inFile = "NEMA_git_Na_1mln2.root"; 
  //const char* inFile = "NEMA_git_1mln_zle.root"; 
  const char* chainName = "GateGlobalActorTree";
  int nEvents = 0;
  TChain chain(chainName);
  chain.Add(inFile);
  ROOT::RDataFrame df(chain);


  auto ranged_df = df.Range(nEvents);

  std::unique_ptr<TCanvas> canv(new TCanvas("canv", "canv", 1920, 1080));
  auto df2 = ranged_df.Alias("Edep", "EnergyLossDuringProcess");
  auto hist = df2.Histo1D({"EnergyBeforeProcess", "EnergyBeforeProcess;E [keV];nevents",1000, 0 ,4000}, "EnergyBeforeProcess");
  auto hist2 = df2.Histo1D({"EmissionEnergyFromSource", "EmissionEnergy;E [keV];nevents",2000, 0 ,2000}, "EmissionEnergyFromSource");
  //auto hist3 = df2.Filter("Edep >200").Filter("EmissionEnergyFromSource ==1157").Filter("VolumeName == \"crystal\"").Histo1D({"Edep", "Edep;Edep [keV];nevents",500, 0 ,4000}, "Edep");
  auto hist3 = df2.Filter("Edep >200").Filter("EmissionEnergyFromSource >1273 && EmissionEnergyFromSource < 1275").Filter("VolumeName == \"crystal1\"").Histo1D({"Edep", "Edep;Edep [keV];nevents",500, 0 ,4000}, "Edep");
  hist->DrawClone();
  hist2->DrawClone();
  hist3->DrawClone();

}

/// General algorithm:
// 1. We take 3 (or more hits) that deposited more than 200 keV
// 2. We form 3 LORs (more?)
// 3. We apply geometrical criterium 
// 4. We apply energetical criterium
// 5. We should apply temporal criterium

using LOR = std::pair<TLorentzVector, TLorentzVector>;
double calculateDistance(const LOR& lor)
{
  return calculateDistance3D(lor.first.X(), lor.first.Y(), lor.first.Z(), lor.second.X(), lor.second.Y(), lor.second.Z());
}

bool isInEnergyRange(double E, double Emin, double Ecut) { return ((E >= Ecut) && (E <= Emin)); }

bool isEqualLor(const LOR& lor1, const LOR& lor2)
{
  return ((lor1.first == lor2.first) && (lor1.second == lor2.second)) || ((lor1.first == lor2.second) && (lor1.second == lor2.first));
}

//kryterium geometryczny
std::vector<LOR> geometry(const LOR& lor1, const LOR& lor2, const LOR& lor3)
{
  std::vector<LOR> all = {lor1, lor2, lor3};
  //sortujemy po długości
  std::sort(all.begin(), all.end(), [](LOR a, LOR b)->bool {
    return calculateDistance(a) < calculateDistance(b);
    });
  //wyrzucamy najdłuższy LOR
  std::vector<LOR> selectedAfterGeomCut = {all[0], all[1]};  
  return selectedAfterGeomCut;
}

// kryterium energetyczny
std::vector<LOR> select2(const TLorentzVector& gamma1, const TLorentzVector& gamma2, const TLorentzVector& gamma3, double Emin, double Ecut)
{
  // są 3 LORa
  LOR lor1 = {gamma1, gamma2};
  LOR lor2 = {gamma1, gamma3};
  LOR lor3 = {gamma2, gamma3};

  std::vector<LOR> selectedAfterGeomCut;

  // sortujemy po długości
  // cout << " " <<endl;
  selectedAfterGeomCut = geometry(lor1, lor2, lor3);
  int licznik = 0;
  int licznik2 = 0;
  std::vector<LOR> finalSelection;
  // kryteruim energetyczny
  for (auto& lor : selectedAfterGeomCut) {
    double E1 = (lor.first.Energy());
    double E2 = (lor.second.Energy());
    // cout << "E1: " <<E1<<endl;
    // cout << "E2: " <<E2<<endl;
    if (isInEnergyRange(E1, Emin, Ecut) && isInEnergyRange(E2, Emin, Ecut)) {
      finalSelection.push_back(lor);
      licznik++;
    }
    licznik2++;
  }
  // cout << "procent energetyczny: " <<licznik*100/licznik2<<endl;
  return finalSelection;
}

void analyzeTree(const std::string &inFile = "outNew2.root",
               const std::string &outFile = "out.root") {
  const char *chainName = "Tree";
  int nEvents = 0;
  TChain chain(chainName);
  chain.Add(inFile.c_str());
  ROOT::RDataFrame df(chain);
  auto getHitsPerEvent = [](const Event &event) -> int {
    int hits = 0;
    for (const auto &track : event.fTracks) {
      for (const auto &hit : track.fHits) {
        if ((!isScatteringInFantom(hit)))
          hits = hits + 1;
      }
    }
    return hits;
  };

  auto getHitsPerEventWithThreshold = [](const Event &event) -> int {
    int hits = 0;
    for (const auto &track : event.fTracks) {
      for (const auto &hit : track.fHits) {
        if (hit.fEnergyBeforeProcess > 200 && (!isScatteringInFantom(hit)))
          hits = hits + 1;
      }
    }
    return hits;
  };
  auto getTracksPerEvent = [](const Event &event) -> int {
    return event.fTracks.size();
  };
  /// Does it make any sense?

  /// This is probably wrong
  auto isTrueThreeHitsInEvent = [](const Event &event) -> bool {
    for (const auto &track : event.fTracks) {
      // auto &steps = track.fHits;
      if (isEqual(track.fEmissionEnergy, 511) ||
          ((track.fEmissionEnergy > 1273) && (track.fEmissionEnergy < 1275))) {
      } else {
        return false;
      }
    }
    return true;
  };

  // auto isEventScatteredInPhantom= [](const Event& event)->int {
  // bool scattered = false;
  // One should define properly scattered in phantom as a hit scattered in
  // phantom and then the track should  contain next hit with energy higher than
  // 200 keV and registered in detector
  // In addition one should probably define a minimum energy loss in scattering
  // in phantom to remove events with almost no scattering
  // for (const auto& track : event.fTracks) {
  // scattered = scattered ||
  // std::any(track.fHits(track.fHits.begin(),track.fHits(track.fHits.end(),
  // isScatteringInFantom);
  //}
  // return scattered;
  //}

  // auto a =  df.Filter([] (const Event& event) { return event.fTracks.size()
  // >0; }, {"Event"}).Count().GetValue();

  auto ranged_df = df.Range(nEvents);
  auto df2 = df.Define("TracksPerEvent", getTracksPerEvent, {"Event"})
                 .Define("HitsPerEvent", getHitsPerEvent, {"Event"})
                 .Define("HitsPerEventWithThreshold",
                         getHitsPerEventWithThreshold, {"Event"});
  auto hist = df2.Histo1D(
      {"HitsPerEvent", "HitsPerEvent;nevents;multiplicity", 30, 0, 30},
      "HitsPerEvent");
  auto hist2 =
      df2.Histo1D({"HitsPerEventWithThreshold",
                   "HitsPerEventWithThreshold;nevents;multiplicity", 30, 0, 30},
                  "HitsPerEventWithThreshold");
  std::unique_ptr<TCanvas> canv(new TCanvas("canv", "canv", 1920, 1080));
  canv->Divide(2, 1);
  canv->cd(1);
  hist->DrawClone();
  canv->cd(2);
  hist2->DrawClone();
  canv->SaveAs(outFile.c_str());

  df2.Filter("HitsPerEventWithThreshold >=3");

  // auto df2 = ranged_df.Alias("Edep", "EnergyLossDuringProcess");
  // auto hist = df2.Histo1D({"EnergyBeforeProcess", "EnergyBeforeProcess;E
  // [keV];nevents",1000, 0 ,4000}, "EnergyBeforeProcess");
  // auto hist2 = df2.Histo1D({"EmissionEnergyFromSource", "EmissionEnergy;E
  // [keV];nevents",500, 0 ,4000}, "EmissionEnergyFromSource");
  // auto hist3 = df2.Filter("Edep >200").Filter("EmissionEnergyFromSource
  // ==1157").Filter("VolumeName == \"crystal\"").Histo1D({"Edep", "Edep;Edep
  // [keV];nevents",500, 0 ,4000}, "Edep");
  // hist->DrawClone();
  // hist2->DrawClone();
}

int main(int argc, char **argv) {
    std::string inFileName;
    std::string outFileName;
  if (argc != 3) {
    std::cerr << "Invalid number of variables." << std::endl;
    std::cerr << "usage : ./start.ext inputFile.root outputFile.root  "
              << std::endl;
    return -1;
  } else {
    inFileName = argv[1];
    outFileName = argv[2];
  }
  analyzeTree(inFileName, outFileName);
  // readTree();
}
