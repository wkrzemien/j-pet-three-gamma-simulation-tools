#include <iostream>
#include "GlobalActorReader.hh"
#include <stdexcept>
#include "TCanvas.h"
#include "TH1F.h"
#include <string>
#include <TRandom3.h>
#include <TMath.h>
#include <algorithm>
#include <TGraph.h>
#include <TLegend.h>
#include <TFile.h>
#include <cmath>
#include <TLorentzVector.h>
#include <cassert>

using namespace std;

struct FullEvent {
  FullEvent()
  {
    reset();
  }

  int fEventID = -1;
  TLorentzVector gammaPrompt;
  TLorentzVector gamma1;
  TLorentzVector gamma2;
  bool gammaPromptScatteredInWater = false;
  bool gamma1ScatteredInWater = false;
  bool gamma2ScatteredInWater = false;

  void reset()
  {
    gamma1.SetXYZT(0, 0, 0, -1);
    gamma2.SetXYZT(0, 0, 0, -1);
    gammaPrompt.SetXYZT(0, 0, 0, -1);
    gammaPromptScatteredInWater = false;
    gamma1ScatteredInWater = false;
    gamma2ScatteredInWater = false;
  }
};

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

/// return true if full event has been read.
/// Track1 ==> prompt
/// Track2 ==> 511 keV
/// Track3 ==> 511 keV
bool readEvent(const GlobalActorReader& gar, const std::string& water_volume_name, FullEvent& outEvent)
{
  int trackID = gar.GetTrackID();
  std::string volume_name = gar.GetVolumeName();
  double energyDeposition = gar.GetEnergyLossDuringProcess();
  double emissionEnergy = gar.GetEmissionEnergyFromSource();
  auto hitPosition = gar.GetProcessPosition(); /// I am not sure here

  int currentEventID = gar.GetEventID();
  int previousEventID = outEvent.fEventID;
  if ( previousEventID != currentEventID ) {
    if ( previousEventID > -1 ) {
      outEvent.fEventID = currentEventID; /// we save new eventID to outEvent
      return true;
    } else {
      outEvent.fEventID = currentEventID; /// this is the first event starting
      return false;
    }
  }

/// if there is multiple scattering we save only the last part
  if ( volume_name == water_volume_name ) {
    /// scatter in volume
    if (trackID == 1) {
      assert(!isEqual(emissionEnergy, 511));
      assert(isEqual(emissionEnergy, 1157));
      outEvent.gammaPromptScatteredInWater = true;
    }
    if (trackID == 2) {
      assert(isEqual(emissionEnergy, 511));
      assert(!isEqual(emissionEnergy, 1157));
      outEvent.gamma1ScatteredInWater = true;
    }
    if (trackID == 3) {
      assert(isEqual(emissionEnergy, 511));
      assert(!isEqual(emissionEnergy, 1157));
      outEvent.gamma2ScatteredInWater = true;
    }
  } else {
    /// scatter in detector
    if (trackID == 1) {
      assert(!isEqual(emissionEnergy, 511));
      assert(isEqual(emissionEnergy, 1157));
      outEvent.gammaPrompt = TLorentzVector(hitPosition, energyDeposition);
    }
    if (trackID == 2) {
      assert(isEqual(emissionEnergy, 511));
      assert(!isEqual(emissionEnergy, 1157));
      outEvent.gamma1 = TLorentzVector(hitPosition, energyDeposition);
    }
    if (trackID == 3) {
      assert(isEqual(emissionEnergy, 511));
      assert(!isEqual(emissionEnergy, 1157));
      outEvent.gamma2 = TLorentzVector(hitPosition, energyDeposition);
    }
  }
  return false;
}

TH1F* hAllPrompt = 0;
TH1F* hScatteredinPhantomPrompt  = 0;
TH1F* hNotScatteredinPhantomPrompt  = 0;
TH1F* hAll511 = 0;
TH1F* hScatteredinPhantom511  = 0;
TH1F* hNotScatteredinPhantom511  = 0;

void createHistograms()
{
  hAllPrompt = new TH1F( "hAllPrompt", "All ", 1200, 0, 1200);
  hScatteredinPhantomPrompt = new TH1F( "hScatteredinPhantomPrompt", "hScatteredinPhantomPrompt", 1200, 0, 1200);
  hNotScatteredinPhantomPrompt = new TH1F( "hNotScatteredinPhantomPrompt", "hNotScatteredinPhantomPrompt", 1200, 0, 1200);

  hAll511 = new TH1F( "hAll511", "All ", 1200, 0, 1200);
  hScatteredinPhantom511 = new TH1F( "hScatteredinPhantom511", "hScatteredinPhantom511", 1200, 0, 1200);
  hNotScatteredinPhantom511 = new TH1F( "hNotScatteredinPhantom511", "hNotScatteredinPhantom511", 1200, 0, 1200);
}

void fillHistograms(const FullEvent& event)
{
  bool isSmearingOn = false;  //change it to true to turn on the smearing
  /// if for gamma energy is less than 0, that means that this this gamma was not registered at all in the scanner
  double promptEnergy = event.gammaPrompt.Energy();
  double gamma1Energy = event.gamma1.Energy();
  double gamma2Energy = event.gamma2.Energy();

  if (isSmearingOn) {
    if (promptEnergy >= 0) {
      promptEnergy = smearEnergy(promptEnergy);
    }
    if (gamma1Energy  >= 0) {
      gamma1Energy = smearEnergy(gamma1Energy);
    }
    if (gamma2Energy  >= 0) {
      gamma2Energy = smearEnergy(gamma2Energy);
    }
  }

  if (promptEnergy > 0) {
    hAllPrompt->Fill(promptEnergy);
    if (event.gammaPromptScatteredInWater) {
      hScatteredinPhantomPrompt->Fill(promptEnergy);
    } else {
      hNotScatteredinPhantomPrompt->Fill(promptEnergy);
    }
  }

  if (gamma1Energy > 0) {
    hAll511->Fill(gamma1Energy);
    if (event.gamma1ScatteredInWater) {
      hScatteredinPhantom511->Fill(gamma1Energy);
    } else {
      hNotScatteredinPhantom511->Fill(gamma1Energy);
    }
  }
  if (gamma2Energy > 0) {
    hAll511->Fill(gamma2Energy);
    if (event.gamma2ScatteredInWater) {
      hScatteredinPhantom511->Fill(gamma2Energy);
    } else {
      hNotScatteredinPhantom511->Fill(gamma2Energy);
    }
  }


  ///To get hit position probably you need:
  /// event.gamma1.Vect().X()
  /// event.gamma1.Vect().Y()
  /// event.gamma1.Vect().Z()
  // Again check if energy>0 otherwise that means that gamma was not registered
}

void saveHistograms()
{

  //histograms prompt
  TCanvas c1("c", "c", 2000, 2000);
  hAllPrompt->SetLineColor(kBlack);
  hAllPrompt->Draw();
  hScatteredinPhantomPrompt->SetLineColor(kRed);
  hScatteredinPhantomPrompt->Draw("same");
  hNotScatteredinPhantomPrompt->SetLineColor(kBlue);
  hNotScatteredinPhantomPrompt->Draw("same");
  c1.SaveAs("prompt.png");

  //histograms 511
  TCanvas c2("c", "c", 2000, 2000);
  hAll511->SetLineColor(kBlack);
  hAll511->Draw();
  hScatteredinPhantom511->SetLineColor(kRed);
  hScatteredinPhantom511->Draw("same");
  hNotScatteredinPhantom511->SetLineColor(kBlue);
  hNotScatteredinPhantom511->Draw("same");
  c2.SaveAs("511.png");

  TFile f("histograms_prompt.root", "recreate");
  hAllPrompt->Write();
  hScatteredinPhantomPrompt->SetLineColor(kRed);
  hScatteredinPhantomPrompt->Write();
  hNotScatteredinPhantomPrompt->SetLineColor(kBlue);
  hNotScatteredinPhantomPrompt->Write();
  f.Close();
  TFile f2("histograms_511.root", "recreate");
  hAll511->Write();
  hScatteredinPhantom511->SetLineColor(kRed);
  hScatteredinPhantom511->Write();
  hNotScatteredinPhantom511->SetLineColor(kBlue);
  hNotScatteredinPhantom511->Write();
  f2.Close();
}

void clearHistograms()
{
  if (hAllPrompt) {
    delete hAllPrompt;
    hAllPrompt = 0;
  }
  if (hScatteredinPhantomPrompt ) {
    delete hScatteredinPhantomPrompt ;
    hScatteredinPhantomPrompt = 0;
  }
  if (hNotScatteredinPhantomPrompt ) {
    delete hNotScatteredinPhantomPrompt ;
    hNotScatteredinPhantomPrompt = 0;
  }
  if (hAll511) {
    delete hAll511;
    hAll511 = 0;
  }
  if (hScatteredinPhantom511 ) {
    delete hScatteredinPhantom511 ;
    hScatteredinPhantom511 = 0;
  }
  if (hNotScatteredinPhantom511 ) {
    delete hNotScatteredinPhantom511 ;
    hNotScatteredinPhantom511  = 0;
  }
}

int main(int argc, char* argv[])
{
  if (argc != 4) {
    std::cerr << "Invalid number of variables." << std::endl;
  } else {
    std::string file_name( argv[1] );
    std::string water_volume_name( argv[2] );
    std::string emissionEnergy(argv[3]);


    FullEvent event;
    bool isNewEvent = false;

    createHistograms();
    try {
      GlobalActorReader gar;
      if (gar.LoadFile(file_name)) {
        while (gar.Read()) {
          isNewEvent = readEvent(gar, water_volume_name, event);
          if (isNewEvent) {
            fillHistograms(event);
            event.reset();
          }
        }
      } else {
        std::cerr << "Loading file failed." << std::endl;
      }
      saveHistograms();

      //Int_t n = 1500, a = 0;
      //Double_t m1, m2, m3, m4;
      //Double_t TPR[n], PPV[n], FPR[n], w[n];
      //for (Int_t i = 0; i < n; i++) {
      //int wbini = hAllPrompt.GetXaxis()->FindBin(i);//zwraca numer binu dla ktrego energia jest i keV
      //int wbinn = hAllPrompt.GetXaxis()->FindBin(n);
      //int gini = hFantomPrompt.GetXaxis()->FindBin(i);
      //int ginn = hAllPrompt.GetXaxis()->FindBin(n);
      //int wbin0 = hFantomPrompt.GetXaxis()->FindBin(a);
      //int gin0 = hFantomPrompt.GetXaxis()->FindBin(a);
      //m4 = hAllPrompt.Integral(wbini, wbinn) - hFantomPrompt.Integral(gini, ginn); //hSignal
      //TPR[i] = (m4) / (hAllPrompt.Integral(wbin0, wbinn) - hFantomPrompt.Integral(gin0, ginn)); //effi
      //PPV[i] = (m4) / (hAllPrompt.Integral(wbini, wbinn)); //purity
      //FPR[i] = (hFantomPrompt.Integral(gin0, gini)) / (hFantomPrompt.Integral(wbin0, wbini));
      //w[i] = i;
      //}

      //TGraph	TFF(n, w, TPR); //efficiency
      //TGraph	TFF2(n, TPR, FPR); //purity
      //TGraph	TFF3(n, w, PPV); //false alarm


      //purity i efficiency
      //TCanvas c2("c2");
      //TFF.Draw();
      ////c2.Print("efficiency_20.png");
      ////TCanvas c4("c4");
      //TFF3.SetLineColor(kRed);
      //TFF3.Draw("same");
      //TLegend legend(0.1, 0.1, 0.1, 0.1);
      ////legend.SetHeader("The Legend Title","C"); // option "C" allows to center the header
      //legend.AddEntry("TFF", "efficiency", "l");
      //legend.AddEntry("TFF3", "purity (czerwony)", "same");
      //legend.Draw();
      //c2.Print("purity_and_efficiency_25.png");

      ////krzywa ROC
      //TCanvas c3("c3");
      //TFF2.Draw("AC+");
      //c3.Print("ROC_25.png");

      clearHistograms();
    } catch (const std::logic_error& e ) {
      std::cerr << e.what() << std::endl;
    } catch (...) {
      std::cerr << "Udefined exception" << std::endl;
    }
  }
  return 0;
}
