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

//TODO: plot energies and energies after selections
/// the main question is why they are so few true events registered from 

using LOR = std::pair<TLorentzVector, TLorentzVector>;
double calculateDistance(const LOR& lor)
{
  return calculateDistance3D(lor.first.X(), lor.first.Y(), lor.first.Z(), lor.second.X(), lor.second.Y(), lor.second.Z());
}

bool isInEnergyRange(double E, double Emin, double Emax) { return ((E >= Emin) && (E <= Emax)); }

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

std::vector<LOR> select(const TLorentzVector& gamma1, const TLorentzVector& gamma2, const TLorentzVector& gamma3, double Ethreshold, double Ecut)
{
  LOR lor1 = {gamma1, gamma2};
  LOR lor2 = {gamma1, gamma3};
  LOR lor3 = {gamma2, gamma3};

  std::vector<LOR> selectedAfterGeomCut;

  selectedAfterGeomCut = geometry(lor1, lor2, lor3);
  std::vector<LOR> finalSelection;
  for (auto& lor : selectedAfterGeomCut) {
    double E1 = (lor.first.Energy());
    double E2 = (lor.second.Energy());
    if (isInEnergyRange(E1, Ethreshold, Ecut) && isInEnergyRange(E2, Ethreshold, Ecut)) {
      finalSelection.push_back(lor);
    }
  }
  return finalSelection;
}

const double kEthreshold = 200;
const double kEcut = 400;
std::vector<LOR> select3(const std::vector<TLorentzVector>& gammas)
{
  if (gammas.size() == 2) 
  {
    std::vector<LOR> lors;
    double E1 = gammas[0].Energy();
    double E2 = gammas[1].Energy();
    if (isInEnergyRange(E1, kEthreshold, kEcut) && isInEnergyRange(E2, kEthreshold, kEcut)) {
      lors= {{gammas[0], gammas[1]}};
    }
    return lors;
  }
  assert(gammas.size() ==3);
  return select(gammas[0], gammas[1], gammas[2], kEthreshold , kEcut);
}

/// Lambdas for selection
auto isEmissionEnergyCloseTo511or1274 = [](double emissionEnergy)->bool 
{
  return isEqual(emissionEnergy, 511, 0.01) ||isEqual(emissionEnergy, 1274.54, 0.01);
};


auto getHitsInDetectorPerEvent = [](const Event& event) -> int {
  int hits = 0;
  for (const auto& track : event.fTracks)
  {
      if (!isEmissionEnergyCloseTo511or1274(track.fEmissionEnergy)) continue;
      for (const auto &hit : track.fHits) {
      if ((!isScatteringInFantom(hit)))
        hits = hits + 1;
      }
  }
  return hits;
};

auto getHitsInDetectorPerEventWithThreshold = [&](const Event& event) -> int {
  int hits = 0;
  for (const auto& track : event.fTracks)
  {
      if (!isEmissionEnergyCloseTo511or1274(track.fEmissionEnergy)) continue;
      for (const auto &hit : track.fHits) {
      if ((hit.fEnergyDeposition > kEthreshold) && (!isScatteringInFantom(hit)) && (isScatteringInAnyDetector(hit)))
      {
        hits = hits + 1;
      }
      }
  }
  return hits;
};

auto getTracksPerEvent = [](const Event& event) -> int { return event.fTracks.size(); };

auto getHitsPositionsAndEnergyInDetector = [&](const Event& event) -> std::vector<TLorentzVector> {
  std::vector<TLorentzVector> photons;
  for (const auto& track : event.fTracks)
  {
    if (!isEmissionEnergyCloseTo511or1274(track.fEmissionEnergy)) continue;
    for (const auto& hit : track.fHits)
    {
      if ((hit.fEnergyDeposition > kEthreshold) && (isScatteringInAnyDetector(hit)))
      {
        TLorentzVector photon(hit.fHitPosition, hit.fEnergyDeposition);
        photons.push_back(photon);
      }
    }
  }
  return photons;
};

auto select511Lors = [](const std::vector<TLorentzVector>& hits) { return select3(hits); };

auto selectTrue511s = [](const Event& event) -> std::vector<TLorentzVector> {
  std::vector<TLorentzVector> photons;
  for (const auto& track : event.fTracks)
  {
    if (!isEmissionEnergyCloseTo511or1274(track.fEmissionEnergy)) continue;
    if (!isEqual(track.fEmissionEnergy, 511, 0.01))
    {
      continue;
    }

    for (const auto& hit : track.fHits)
    {
      if ((hit.fEnergyDeposition > kEthreshold) && (isScatteringInAnyDetector(hit)))
      {
        TLorentzVector photon(hit.fHitPosition, hit.fEnergyDeposition);
        photons.push_back(photon);
      }
    }
  }
  return photons;
};

auto isSelectionCorrect = [](const std::vector<LOR>& lors, std::vector<TLorentzVector>& true511s) -> int {
  if (true511s.size() < 2)
    return 1;
  if (lors.size() == 0)
    return 2;
  if (lors.size() > 1)
    return 3;
  LOR trueLOR = {true511s[0], true511s[1]};
  auto recoLOR = lors[0];
  if (!isEqualLor(trueLOR, recoLOR))
    return 4;
  return 0;
};

// auto isEventScatteredInPhantom= [](const Event& event)->int {
// bool scattered = false;
// One should define properly scattered in phantom as a hit scattered
// in
// phantom and then the track should  contain next hit with energy
// higher than
// 200 keV and registered in detector
// In addition one should probably define a minimum energy loss in
// scattering
// in phantom to remove events with almost no scattering
// for (const auto& track : event.fTracks) {
// scattered = scattered ||
// std::any(track.fHits(track.fHits.begin(),track.fHits(track.fHits.end(),
// isScatteringInFantom);
//}
// return scattered;
//}

void analyzeTree(const std::string &inFile = "outNew2.root",
               const std::string &outFile = "out.root") {

  const char *chainName = "Tree";
  int nEvents = 0;
  TChain chain(chainName);
  chain.Add(inFile.c_str());
  ROOT::RDataFrame df(chain);


  auto ranged_df = df.Range(nEvents);
  auto df2 = df.Define("TracksPerEvent", getTracksPerEvent, {"Event"})
                 .Define("HitsPerEvent", getHitsInDetectorPerEvent, {"Event"})
                 .Define("HitsPerEventWithThreshold",
                         getHitsInDetectorPerEventWithThreshold, {"Event"});

  auto oneHit = df2.Filter("HitsPerEventWithThreshold ==1").Count().GetValue();
  auto twoHit = df2.Filter("HitsPerEventWithThreshold ==2").Count().GetValue();
  auto threeHit = df2.Filter("HitsPerEventWithThreshold ==3").Count().GetValue();
  auto fourHit = df2.Filter("HitsPerEventWithThreshold ==4").Count().GetValue();
  
  std::cout <<  oneHit << std::endl;
  std::cout << twoHit << std::endl;
  std::cout << threeHit << std::endl;
  std::cout << fourHit << std::endl;

  auto filtered =
      df2.Filter("(HitsPerEventWithThreshold ==2) ||(HitsPerEventWithThreshold ==3) ")
          .Define("HitsVectors", getHitsPositionsAndEnergyInDetector, {"Event"})
          .Define("true511", selectTrue511s, {"Event"})
          .Define("selected511Lors", select511Lors, {"HitsVectors"})
          .Define("selectionResults", isSelectionCorrect,
                  {"selected511Lors", "true511"});

  auto filtered2H = filtered.Filter("HitsPerEventWithThreshold ==2");
  auto filtered3H = filtered.Filter("HitsPerEventWithThreshold ==3");
  auto histSelectionResults2 = filtered2H.Histo1D({"selectionResults2Hits", ";nevents;results flag", 5, -0.5, 4.5}, "selectionResults");
  auto histSelectionResults3 = filtered3H.Histo1D({"selectionResults3Hits", ";nevents;results flag", 5, -0.5, 4.5}, "selectionResults");
  auto histX = new TH1F("Edep", "Edep",650, 0, 1300); 
  auto fillEnergyHist = [&](const std::vector<TLorentzVector>& hits)->void { if(hits.size()!=2) return; for(const auto& h: hits) histX->Fill(h.Energy());
  };
  filtered2H.Foreach(fillEnergyHist,{"true511"});
  
  //auto hist = df2.Histo1D(
      //{"HitsPerEvent", "HitsPerEvent;nevents;multiplicity", 30, 0, 30},
      //"HitsPerEvent");
  //auto hist2 =
      //df2.Histo1D({"HitsPerEventWithThreshold",
                   //"HitsPerEventWithThreshold;nevents;multiplicity", 30, 0, 30},
                  //"HitsPerEventWithThreshold");

  std::unique_ptr<TCanvas> canv(new TCanvas("canv", "canv", 1920, 1080));
  canv->Divide(3, 1);
  canv->cd(1);
  histSelectionResults2->Draw(); 
  canv->cd(2);
  histSelectionResults3->Draw(); 
  canv->cd(3);
  histX->Draw();
  //hist->DrawClone();
  //hist2->DrawClone();
  canv->SaveAs(outFile.c_str());
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
