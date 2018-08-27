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

int gAllPrompt = 0;
int gAll511 = 0;

struct FullEvent {
  FullEvent()
  {
    reset();
  }

  int fEventID = -1;
  TLorentzVector gammaPrompt;
  TLorentzVector gamma1;
  TLorentzVector gamma2;
  //bool gammaPromptScatteredInWater = false;
  //bool gamma1ScatteredInWater = false;
  //bool gamma2ScatteredInWater = false;

  void reset()
  {
    gamma1.SetXYZT(0, 0, 0, -1);
    gamma2.SetXYZT(0, 0, 0, -1);
    gammaPrompt.SetXYZT(0, 0, 0, -1);
    //gammaPromptScatteredInWater = false;
    //gamma1ScatteredInWater = false;
    //gamma2ScatteredInWater = false;
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

/*double calculateDistance(double x1, double y1, double x2, double y2, double radius)
{
double distance =0;
distance= sqrt(pow((radius/2),2)-pow((x2-x1),2)-pow((y2-y1),2));
return distance;
}*/

/// returns pair<bool,FullEvent >  if true the return FullEvent can be saved
/// Track1 ==> prompt
/// Track2 ==> 511 keV
/// Track3 ==> 511 keV
std::pair<bool, FullEvent> readEvent(const GlobalActorReader& gar, FullEvent& outEvent)
{
  int trackID = gar.GetTrackID();
  //std::string volume_name = gar.GetVolumeName();
  double energyBeforeProcess = gar.GetEnergyBeforeProcess();
  double energyDeposition = gar.GetEnergyLossDuringProcess();
  double emissionEnergy = gar.GetEmissionEnergyFromSource();
  auto hitPosition = gar.GetProcessPosition(); /// I am not sure here

  int currentEventID = gar.GetEventID();
  int previousEventID = outEvent.fEventID;
  bool isNewEvent = false;
  FullEvent lastEvent = outEvent;
  if ( previousEventID == currentEventID ) {
    isNewEvent = false;
  } else {
    if ( previousEventID > -1 ) {
      isNewEvent = true;
      outEvent.reset();
      outEvent.fEventID = currentEventID; /// we save new eventID to outEvent
    } else {
      isNewEvent = false;
      outEvent.fEventID = currentEventID; /// this is the first event starting
    }
  }
/// we save only the first scattering!!!
  if (isEqual(emissionEnergy, 511) || isEqual(emissionEnergy, 1157)) {
    /// scatter in detector
    if (trackID == 1) {
      assert(!isEqual(emissionEnergy, 511));
      assert(isEqual(emissionEnergy, 1157));
      ///This means it is the first time in the detector
      /// This condition will not work with in-phantom scattering!!!
      if (isEqual(energyBeforeProcess, 1157)) {
        outEvent.gammaPrompt = TLorentzVector(hitPosition, energyDeposition);
      }
      /// else this means that prompt scatters for the second time in the detector
    }
    if (trackID == 2) {
      assert(isEqual(emissionEnergy, 511));
      assert(!isEqual(emissionEnergy, 1157));
      /// This condition will not work with in-phantom scattering!!!
      if (isEqual(energyBeforeProcess, 511)) {
        outEvent.gamma1 = TLorentzVector(hitPosition, energyDeposition);
      }
    }
    if (trackID == 3) {
      assert(isEqual(emissionEnergy, 511));
      assert(!isEqual(emissionEnergy, 1157));
      /// This condition will not work with in-phantom scattering!!!
      if (isEqual(energyBeforeProcess, 511)) {
        outEvent.gamma2 = TLorentzVector(hitPosition, energyDeposition);
      }
    }
  }
  return std::make_pair(isNewEvent, lastEvent);
}

TH1F* hAllPrompt = 0;
TH1F* hAll511 = 0;
TH1F* h3detPrompt = 0;
TH1F* h2det1Prompt = 0;
TH1F* h2det2Prompt = 0;
TH1F* hPromptMult = 0;
TH1F* h511Mult = 0;

TGraph* purity_prompt = 0;
TGraph* efficiency_prompt = 0;
TGraph* ROC_prompt = 0;
TGraph* purity_511 = 0;
TGraph* efficiency_511 = 0;
TGraph* ROC_511 = 0;

void createHistograms()
{
  hAllPrompt = new TH1F( "hAllPrompt", "All ", 500, 0, 1200);
  hAllPrompt->GetXaxis()->SetCanExtend(true);
  hAll511 = new TH1F( "hAll511", "All ", 500, 0, 1200);
  h3detPrompt = new TH1F( "h3detPrompt", "Prompt with 3 det", 500, 0, 1200);
  h2det1Prompt = new TH1F( "h2det1Prompt", "Prompt with det prompt and gamma1 ", 500, 0, 1200);
  h2det2Prompt = new TH1F( "h2det2Prompt", "Prompt with det prompt and gamma2  ", 500, 0, 1200);
  hPromptMult = new TH1F("hPromptMult","Prompt multiplicity", 2, -0.5, 1.5); 
  h511Mult = new TH1F("h511Mult","511 gamma multiplicity", 3, -0.5, 2.5);

}

void fillHistograms(const FullEvent& event)
{
  bool isSmearingOn = false;  //change it to true to turn on the smearing
  bool isLowEnergyCutOn = false; //chane it to true to turn on the low energy cut
  const double kLowEnergyThreshold = 100; /// can be also 50 keV
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

  if (isLowEnergyCutOn) {
    if (promptEnergy < kLowEnergyThreshold) {
      promptEnergy = -1;

    }
    if (gamma1Energy  < kLowEnergyThreshold ) {
      gamma1Energy = -1;

    }
    if (gamma2Energy  < kLowEnergyThreshold) {
      gamma2Energy = -1;

    }
  }


  if (promptEnergy >= 0) {
    hAllPrompt->Fill(promptEnergy);
    hPromptMult->Fill(1);
    gAllPrompt++;
  } else{
    hPromptMult->Fill(0);
  }

  if (gamma1Energy > 0) {
    hAll511->Fill(gamma1Energy);
    gAll511++;
  }
  if (gamma2Energy > 0) {
    hAll511->Fill(gamma2Energy);
    gAll511++;
  }

  if (promptEnergy > 0 && gamma1Energy > 0 && gamma2Energy > 0) {
    h3detPrompt->Fill(promptEnergy);
  }
  if (promptEnergy > 0 && gamma1Energy > 0) {
    h2det1Prompt->Fill(promptEnergy);
    h511Mult->Fill(2);  
  }
  if (promptEnergy > 0 && gamma2Energy > 0) {
    h2det2Prompt->Fill(promptEnergy);
  }
  
  //only 1 511 detected
  if((gamma1Energy >=0 && gamma2Energy < 0) ||(gamma1Energy < 0 && gamma2Energy >= 0)){
    h511Mult->Fill(1);  
  }
  //no 511 gammas detected
  if(gamma1Energy < 0 && gamma2Energy < 0){
    h511Mult->Fill(0);  
  }

  //new TH1F if(promptEnergy>100&&gamma1Energy>100&&gamma2Energy>100)

}

void saveHistograms()
{

  TCanvas c1("c", "c", 2000, 2000);
  hAll511->SetLineColor(kBlack);
  hAll511->Draw();
  hAllPrompt->SetLineColor(kRed);
  hAllPrompt->Draw("same");
  c1.SaveAs("all.png");

  TCanvas c12("c", "c", 2000, 2000);
  h3detPrompt->SetLineColor(kBlack);
  h3detPrompt->Draw();
  h3detPrompt->SetLineColor(kRed);
  c12.SaveAs("h3detPrompt.png");

  TCanvas c13("c", "c", 2000, 2000);
  h2det1Prompt->SetLineColor(kBlack);
  h2det1Prompt->Draw();
  c13.SaveAs("h2det1Prompt.png");

  TCanvas c14("c", "c", 2000, 2000);
  h2det2Prompt->SetLineColor(kBlack);
  h2det2Prompt->Draw();
  c14.SaveAs("h2det2Prompt.png");

  TFile f("histograms_all.root", "recreate");
  hAllPrompt->Write();
  hAll511->Write();
  h3detPrompt->Write();
  h2det1Prompt->Write();
  h2det2Prompt->Write();
  f.Close();
}
void drawPlot()
{
  const double kEnergyMin = 0;
  const double kEnergyMax = 1158;
  const int kNumberOfSteps = 1158;
  const double kEnergyStep = (kEnergyMax - kEnergyMin) / kNumberOfSteps;
  Double_t TPR_prompt[kNumberOfSteps];
  Double_t PPV_prompt[kNumberOfSteps];
  Double_t FPR_prompt[kNumberOfSteps];
  Double_t TPR_511[kNumberOfSteps];
  Double_t PPV_511[kNumberOfSteps];
  Double_t FPR_511[kNumberOfSteps];
  Double_t energies[kNumberOfSteps];


  int minBinPrompt = hAllPrompt->GetXaxis()->FindBin(kEnergyMin); /// min bin for Promtp
  int maxBinPrompt = hAllPrompt->GetXaxis()->FindBin(kEnergyMax); /// max bin for Prompt
  int minBin511 = hAll511->GetXaxis()->FindBin(kEnergyMin); // min bin for 511
  int maxBin511 = hAll511->GetXaxis()->FindBin(kEnergyMax); /// max bin for 511

  double allPrompt = hAllPrompt->Integral(minBinPrompt, maxBinPrompt);
  double all511 = hAll511->Integral(minBin511, maxBin511);

  for (Int_t i = 0; i < kNumberOfSteps ; i++) {
    double currentEnergy =  kEnergyStep * i;
    energies[i] = currentEnergy;

    int currentBin511 = hAll511->GetXaxis()->FindBin(currentEnergy);
    int currentBinPrompt = hAllPrompt->GetXaxis()->FindBin(currentEnergy);

    double allAcceptedAsPrompt = hAllPrompt->Integral(currentBinPrompt, maxBinPrompt) + hAll511->Integral(currentBin511, maxBin511);
    double promptAcceptedAsPrompt = hAllPrompt->Integral(currentBinPrompt, maxBinPrompt);
    double v511AcceptedAsPrompt = hAll511->Integral(currentBin511, maxBin511);

    if (allAcceptedAsPrompt !=0) {
      PPV_prompt[i] = promptAcceptedAsPrompt/allAcceptedAsPrompt; //purity
    } else {
      PPV_prompt[i] = 0;
    }
    TPR_prompt[i] = promptAcceptedAsPrompt/allPrompt; //effi
    /// This is the porcentage of 511 that is classifed falsely as prompt
    FPR_prompt[i] = v511AcceptedAsPrompt / all511;

    double allAcceptedAs511 = (hAllPrompt->Integral(minBinPrompt, currentBinPrompt) + hAll511->Integral(minBin511, currentBin511));
    double v511AcceptedAs511 = hAll511->Integral(minBin511, currentBin511);
    double promptAcceptedAs511 = hAllPrompt->Integral(minBinPrompt, currentBinPrompt);

    if (allAcceptedAs511 !=0) {
      PPV_511[i] = v511AcceptedAs511 /allAcceptedAs511; //purity
    } else {
      PPV_511[i] = 0;
    }
    TPR_511[i] = v511AcceptedAs511 / all511; //effi
    /// This is the porcentage of prompt that is classifed falsely as 511
    FPR_511[i] = promptAcceptedAs511  / allPrompt;
  }

  ROC_prompt = new TGraph(kNumberOfSteps, FPR_prompt, TPR_prompt);
  efficiency_prompt = new TGraph(kNumberOfSteps, energies, TPR_prompt);
  purity_prompt = new TGraph(kNumberOfSteps, energies, PPV_prompt);

  ROC_511 = new TGraph(kNumberOfSteps, FPR_511, TPR_511);
  efficiency_511 = new TGraph(kNumberOfSteps, energies, TPR_511);
  purity_511 = new TGraph(kNumberOfSteps, energies, PPV_511);
}

void savePlot()
{
  TCanvas c2("c", "c", 2000, 2000);
  purity_prompt->SetLineColor(kBlack);
  purity_prompt->SetTitle("purity prompt");
  purity_prompt->GetXaxis()->SetTitle("Energy cut [MeV]");
  purity_prompt->GetYaxis()->SetTitle("PPV");
  purity_prompt->Draw();
  c2.SaveAs("purity_prompt.png");

  TCanvas c3("c", "c", 2000, 2000);
  purity_511->SetLineColor(kBlack);
  purity_511->SetTitle("purity 511");
  purity_511->GetXaxis()->SetTitle("Energy cut [MeV]");
  purity_511->GetYaxis()->SetTitle("PPV");
  purity_511->Draw();
  c3.SaveAs("purity_511.png");

  TCanvas c4("c", "c", 2000, 2000);
  efficiency_511->SetLineColor(kRed);
  efficiency_511->SetTitle("efficiency 511");
  efficiency_511->GetXaxis()->SetTitle("Energy cut [MeV]");
  efficiency_511->GetYaxis()->SetTitle("TPR");
  efficiency_511->Draw();
  c4.SaveAs("efficiency_511.png");

  TCanvas c5("c", "c", 2000, 2000);
  efficiency_prompt->SetLineColor(kRed);
  efficiency_prompt->SetTitle("efficiency prompt");
  efficiency_prompt->GetXaxis()->SetTitle("Energy cut [MeV]");
  efficiency_prompt->GetYaxis()->SetTitle("TPR");
  efficiency_prompt->Draw();
  c5.SaveAs("efficiency_prompt.png");

  TCanvas c6("c", "c", 2000, 2000);
  ROC_prompt->SetLineColor(kBlack);
  ROC_prompt->SetTitle("ROC prompt");
  ROC_prompt->GetXaxis()->SetTitle("FPR");
  ROC_prompt->GetYaxis()->SetTitle("TPR");
  ROC_prompt->Draw();
  c6.SaveAs("ROc_prompt.png");

  TCanvas c7("c", "c", 2000, 2000);
  ROC_511->SetLineColor(kBlack);
  ROC_511->SetTitle("ROC 511");
  ROC_511->GetXaxis()->SetTitle("FPR");
  ROC_511->GetYaxis()->SetTitle("TPR");
  ROC_511->Draw();
  c7.SaveAs("ROc_511.png");

  TFile f("histograms_all.root", "recreate");
  efficiency_511->Write();
  purity_511->Write();
  ROC_511->Write();
  efficiency_prompt->Write();
  purity_prompt->Write();
  ROC_prompt->Write();
  hPromptMult->Write();
  h511Mult->Write();
  f.Close();
}

void clearHistograms()
{
  if (hAllPrompt) {
    delete hAllPrompt;
    hAllPrompt = 0;
  }

  if (hAll511) {
    delete hAll511;
    hAll511 = 0;
  }

  if (efficiency_prompt) {
    delete efficiency_prompt;
    efficiency_prompt = 0;
  }

  if (purity_prompt) {
    delete purity_prompt;
    purity_prompt = 0;
  }

  if (ROC_prompt) {
    delete ROC_prompt;
    ROC_prompt = 0;
  }

  if (efficiency_511) {
    delete efficiency_511;
    efficiency_511 = 0;
  }

  if (purity_511) {
    delete purity_511;
    purity_511 = 0;
  }

  if (ROC_511) {
    delete ROC_511;
    ROC_511 = 0;
  }

}

int main(int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Invalid number of variables." << std::endl;
  } else {
    std::string file_name( argv[1] );
    std::string emissionEnergy(argv[2]);

    FullEvent event;

    createHistograms();
    try {
      GlobalActorReader gar;
      if (gar.LoadFile(file_name)) {
        while (gar.Read()) {
          auto res = readEvent(gar, event);
          if (res.first) {
            fillHistograms(res.second);
          }
        }
      } else {
        std::cerr << "Loading file failed." << std::endl;
      }
      drawPlot();
      saveHistograms();
      savePlot();
      clearHistograms();
      
      //calculating efficiencies of registration 
      const double  kNumberOfGeneratedEvents = 5e6;
      double PromptEff = (double)gAllPrompt/kNumberOfGeneratedEvents;
      double Eff511 = (double)gAll511/(2*kNumberOfGeneratedEvents);
      std::cout << "Efficiencies: " << std::endl;
      std::cout << "GammaPrompt[%]  " << 100.*PromptEff <<std::endl;
      std::cout << "Gamma 511[%] " << 100.*Eff511 <<std::endl;
    } catch (const std::logic_error& e ) {
      std::cerr << e.what() << std::endl;
    } catch (...) {
      std::cerr << "Udefined exception" << std::endl;
    }
  }
  return 0;
}
