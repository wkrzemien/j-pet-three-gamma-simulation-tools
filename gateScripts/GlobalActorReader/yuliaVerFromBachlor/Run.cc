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

using namespace std;
using namespace helper_functions;

struct Counter {
  int numericprompt = 0;
  int numericgamma = 0;
  int numeric_dist_gamma1gamma2_max = 0;
  int numeric_dist_gamma1gamma2_min = 0;
  int numeric_dist_gamma1gamma2_mid = 0;
  int numeric_dist_gamma1prompt_max = 0;
  int numeric_dist_gamma2prompt_max = 0;
}; 

TGraph* purityPromptEcut0 = 0 ;
TGraph* efficiencyPromptEcut0 = 0;
TGraph* ROCpromptEcut0 = 0;
TGraph* F1PromptEcut0 = 0;
TGraph* purityPromptEcutN = 0 ;
TGraph* efficiencyPromptEcutN = 0;
TGraph* ROCpromptEcutN = 0;
TGraph* F1PromptEcutN = 0;

TGraph* purity511Ecut0 = 0;
TGraph* efficiency511Ecut0 = 0;
TGraph* ROC511Ecut0 = 0;
TGraph* F1_511Ecut0 = 0;
TGraph* purity511EcutN = 0;
TGraph* efficiency511EcutN = 0;
TGraph* ROC511EcutN = 0;
TGraph* F1_511EcutN = 0;

TGraph* F1_max_511 = 0;
TGraph* F1_max_prompt = 0;

TGraph* xy511_1 = 0;
TGraph* xy511_2 = 0;
TGraph* xz511_1 = 0;
TGraph* xz511_2 = 0;
TGraph* zy511_1 = 0;
TGraph* zy511_2 = 0;
TGraph* xzprompt = 0;
TGraph* xyprompt = 0;
TGraph* zyprompt = 0;

TGraph* klas0 = 0;
TGraph* klas1 = 0;
TGraph* klas2 = 0;
TGraph* klas3 = 0;

std::vector <TLorentzVector> gammaPromptPos;
std::vector <TLorentzVector> gamma511Pos1;
std::vector <TLorentzVector> gamma511Pos2;

TLorentzVector gammaPrompt;
TLorentzVector gamma1;
TLorentzVector gamma2;

/*----------------------------------------------------OPIS FUNKCJI--------------------------------------------------------------------
(double energy)----------------------------------------------------------------- rozmycie energii
isEqual(double x, double y, double epsilon = 10e-9)---------------------------------------- porównanie wartości
calculateDistance3D(double x1, double y1, double x2, double y2, double z1, double z2)------ distance dla LOR 3D
double calculateDistance2D(double x1, double y1, double x2, double y2)--------------------- distance dla LOR 2D
exactlyOneHit = [](const Track& track)----------------------------------------------------- 1 hit dla jednego eventa
oneOrMoreHits = [](const Track& track)----------------------------------------------------- 1 albo więcej hitów dla jednego eventa
exactlyTwoHitsInEvent = [](const Event& event)--------------------------------------------- 2 hity dla jednego eventa
exactlyThreeHitsInEvent = [](const Event& event)------------------------------------------- 3 hity dla jednego eventa
nOrMoreHitsInEvent = [](const Event& event, int N)----------------------------------------- dowolna liczba hitów dla jednego eventa
numOfHitsInEvent(const Event& event)------------------------------------------------------- liczba hitów w jednym evencie
isScatteringInDetector1(const Hit& step)--------------------------------------------------- hity zarejestrowane w 1 warstwie detektora
isScatteringInDetector2(const Hit& step)--------------------------------------------------- hity zarejestrowane w 2 warstwie detektora
isScatteringInDetector3(const Hit& step)--------------------------------------------------- hity zarejestrowane w 3 warstwie detektora
scattering511 = [](const Hit& hit)--------------------------------------------------------- emitowany foton ma energie 511
scatteringprompt = [](const Hit& hit)------------------------------------------------------ emitowany foton ma energie 1157
--------------------------------------------------------------------------------------------------------------------------------------*/
using LOR = std::pair<TLorentzVector, TLorentzVector> ;

std::vector<LOR> select2(const TLorentzVector& gamma1, const TLorentzVector& gamma2, const TLorentzVector& gamma3, double Emin, double Ecut);

bool isEqualLor(const LOR& lor1, const LOR& lor2);

void printResults(int numOfEvents, vector<int> resultNum, int Emin, int k)
{
  std::cout <<"------------------------------------------ "<<endl;
  std::cout <<"Wynik: "<<endl;
  std::cout <<""<<endl;
  std::cout << "dla E cut= " << k  <<endl;
  std::cout << "best E treshold= " << Emin  <<endl;
  std::cout <<" "<<endl;
  std::cout <<"Definicja wyników działania algorythmu selekcyjnego:"<<endl;
  std::cout <<"Klasa 0 - wybral dobrego LORa(511 *2)"<<endl;
  std::cout <<"klasa 1 - wybrala zle"<<endl;
  std::cout <<"klasa 2 - odrzucila trzy LORY"<<endl;
  std::cout <<"klasa 3 - obie zaakceptowal 2 LORY"<<endl;
  std::cout <<" "<<endl;
  int licznik=0;
  for (auto r:resultNum) {
  std::cout <<"klasa:"<<licznik<< ", ilosc:"<<r<< ", procent:"<<float(r)/numOfEvents * 100.<<endl;
  licznik++;
  }
    std::cout <<" "<<endl;
}

void printevents(int numOfEvents, const Counter& counter, double percent){
  std::cout <<" "<<endl;std::cout << "Results for back-to-back: " << std::endl;
  std::cout << "Liczba gamma dla min odleglosci  " << counter.numeric_dist_gamma1gamma2_min << std::endl;
  std::cout << "Liczba gamma dla max odleglosci  " << counter.numeric_dist_gamma1gamma2_max << std::endl;
  std::cout << "Liczba gamma dla mid odleglosci  " << counter.numeric_dist_gamma1gamma2_mid << std::endl;
  std::cout << "Liczba zaakceptowanych prompt:  " << counter.numericprompt << std::endl;
  std::cout << "Liczba zaakceptowanych gamma jako prompt:  " << counter.numericgamma << std::endl;
  std::cout << "Event steps  " << numOfEvents << std::endl;
  std::cout <<" "<<endl;
  std::cout << "maksymalna odległość dla gamma1prompt:  " << counter.numeric_dist_gamma1prompt_max << std::endl;
  std::cout << "maksymalna odległość dla gamma2prompt:  " << counter.numeric_dist_gamma2prompt_max << std::endl;
  std::cout << "maksymalna odległość dla gamma1gamma2: " << counter.numeric_dist_gamma1gamma2_max << std::endl;
  std::cout << "finalSelection dla 511-511 = max " << counter.numeric_dist_gamma1gamma2_max << std::endl;
  std::cout << "Percent of false " << percent << "%" << std::endl;
}

void analyse(const std::string& inFile)
{

  //inicjalizujemy histogramy
  TFile testOut("testOutSimple.root", "RECREATE");
  TH1F histMultiplicity("histMultiplicity", "histMultiplicity", 20, 0, 20);
  TH1F histMultiplicity2("histMultiplicity2", "histMultiplicity2", 20, 0, 20);
  TH1F hgammaenergyMax("hgammaenergyMax", "hgammaenergyMax", 250, 0, 1200);
  TH1F hgammapromptenergyMax("hgammapromptenergyMax", "hgammapromptenergyMax", 250, 0, 1200);
  TH1F hgammaenergyMid("hgammaenergyMid", "hgammaenergyMid", 250, 0, 1200);
  TH1F hgammapromptenergyMid("hgammapromptenergyMid", "hgammapromptenergyMid", 250, 0, 1200);
  TH1F hgammaenergyMin("hgammaenergyMin", "hgammaenergyMin", 250, 0, 1200);
  TH1F hgammapromptenergyMin("hgammapromptenergyMin", "hgammapromptenergyMin", 250, 0, 1200);
  TH1F hDISTANCEgamma1prompt("hDISTANCEgamma1prompt", "hpromptZ", 250, 0, 600);
  TH1F hDISTANCEgamma2prompt("hDISTANCEgamma2prompt", "hgamma2Y", 250, 0, 600);
  TH1F hDISTANCEgamma1gamma2("hDISTANCEgamma1gamma2", "hgamma2Z", 250, 0, 600);
  TH1F h511Sygnal("h511Sygnal", "h",100, 0, 1200); // 0 fantom && 1 detector
  TH1F hpromptSygnal("hpromptSygnal", "h", 100, 0, 1200); //0 fantom && 1 detector
  TH1F h511B1B2("h511B1B2", "h", 100, 0, 1200); //1 fantom && 0 detector + n fantom && 0 detector
  TH1F hpromptB1B2("hpromptB1B2", "h", 100, 0, 1200); //1 fantom && 0 detector + n fantom && 0 detector
  TH1F h511B1("h511B1", "h", 100, 0, 1200); //1 fantom && 0 detector
  TH1F hpromptB1("hpromptB1", "h", 100, 0, 1200); //1 fantom && 0 detector
  TH1F h511B2("h511B2", "h", 100, 0, 1200); //n fantom && 0 detector
  TH1F hpromptB2("hpromptB2", "h", 100, 0, 1200); //n fantom && 0 detector
  TH1F h511B3("h511B3", "h", 100, 0, 1200); //0 fantom && n detector
  TH1F h511B4("h511B4", "h", 100, 0, 1200); //0 fantom && n detector
  TH1F hpromptB3("hpromptB3", "h", 100, 0, 1200); //0 fantom && n detector
  TH1F h511SygnalB1B2("h511SygnalB1B2", "h detector", 100, 0, 1200); //all types fantom && all types detector
  TH1F hpromptSygnalB1B2("hpromptSygnalB1B2", "h", 100, 0, 1200); //all types fantom && all types detector
  TH1F hostalos("hostalos", "h", 100, 0, 1200);
  TH1F hDwa("hDwa", "h", 100, 0, 1200);
  TH1F hDwaFantom("hDwaFantom", "h", 100, 0, 1200);
  TH1F hDwaAll("hDwaAll", "h", 100, 0, 1200);

  //wczytujemy file
  TFile file(inFile.c_str(), "READ");
  TTreeReader reader("Tree", &file);
  TTreeReaderValue<Event> event(reader, "Event");
  int treevent = 0;
  int dwaevent = 0;
  int truedwaevent = 0;
  int eve1 = 0;
  int eve4 = 0;
  int eve5 = 0;
  int evemore=0;
  double smearEnergy511;
  double smearEnergy5112;
  double smearEnergyPrompt;

  Double_t x5111[16000];
  Double_t y5112[16000];
  Double_t z5112[16000];
  Double_t z5111[16000];
  Double_t x5112[16000];
  Double_t y5111[16000];
  Double_t xprompt[16000];
  Double_t yprompt[16000];
  Double_t zprompt[16000];

//----------------------------------------------------------SORTOWANIE FOTONÓW PO GRUPACH-----------------------------------------------------------------------
  while (reader.Next()) {
    histMultiplicity.Fill(event->fTracks.size()); //histogram z liczbą rozproszeń w evencie
    histMultiplicity2.Fill(numOfHitsInEvent(*event)); //histogram z 3 hitami w jednym evencie
    if (exactlyTwoHitsInEvent(*event)) {
      dwaevent++;
      if (trueTwoHitsInEvent(*event)){
        truedwaevent++;
      }
      for (const auto& track : event->fTracks) {
      double emissionEnergy = track.fEmissionEnergy;
      auto &steps = track.fHits;
      auto& rozp = track.fTrackID;
      for (auto i = 0u; i < steps.size(); i++) 
      {
        auto &hit = steps[i];
        hDwaAll.Fill(hit.fEnergyDeposition);
        if(isScatteringInDetector1(hit) || isScatteringInDetector2(hit) || isScatteringInDetector3(hit))
        {
          hDwa.Fill(hit.fEnergyDeposition);
        }
        if(isScatteringInFantom2(hit))
        {
          hDwaFantom.Fill(hit.fEnergyDeposition);
        }
      }
    }
    }
     if (exactly1HitsInEvent(*event)) {
eve1++;
    }
      if (exactlyMoreHitsInEvent(*event)) {
        evemore++;
      }

    if (numOfHitsInEvent(*event) < 3) { 
      continue;
    }
    if (exactlyThreeHitsInEvent(*event)) {
      treevent++;
      
      cout << " " << endl;
      cout << " " << endl;
      cout << "Exactly 3 hits per event:" << endl;
      cout << "Nr zdarzenia: " << treevent << endl;
      if (trueThreeHitsInEvent(*event)){
      for (const auto& track : event->fTracks) {
      bool isInPhantom = false;
      bool isInPhantom2 = false;
      bool isInDetector = false;
      bool isInDetector1 = false;
      bool isInDetector2 = false;
      bool isInDetector3 = false;
      bool is511 = false;
      bool isprompt = false;
      double emissionEnergy = track.fEmissionEnergy;
      auto &steps = track.fHits;
      auto& rozp = track.fTrackID;
      for (auto i = 0u; i < steps.size(); i++) 
      {
        auto &hit = steps[i];
        is511 =scattering511(hit);
        isprompt=scatteringprompt(hit);
        isInDetector1 = isScatteringInDetector1(hit);
        isInDetector2 = isScatteringInDetector2(hit);
        isInDetector3 = isScatteringInDetector3(hit);
        isInPhantom = isScatteringInFantom(hit); //fantom NEMA
        isInPhantom2 = isScatteringInFantom2(hit); //fantom wodny
        isInDetector = isInDetector1  || isInDetector2 || isInDetector3;
        isInPhantom2 = !isInDetector;
        isInPhantom = !isInDetector;

          if (isEqual(emissionEnergy, 511) && gamma511Pos1.size()< gamma511Pos2.size() ) {
            smearEnergy511 = smearEnergy(hit.fEnergyDeposition);
            gamma1 = TLorentzVector(hit.fHitPosition, smearEnergy511);
            gamma511Pos1.push_back(gamma1);
            x5111[treevent] = hit.fHitPosition.x();
            y5111[treevent] = hit.fHitPosition.y();
            z5111[treevent] = hit.fHitPosition.z();
            cout << " " << endl;
            cout << "First gamma-511" << endl; 
            cout << "energy deposition= " << hit.fEnergyDeposition << endl;
            cout << "emission energy = " << emissionEnergy << endl;
            cout << "number of Track = " << rozp << endl;
            cout << "energy after smearing = " << gamma511Pos1[gamma511Pos1.size()-1].Energy() << endl;
            cout << "number of event = " << gamma511Pos1.size() << endl;
            h511SygnalB1B2.Fill(smearEnergy511); //zapisujemy wszystkie 511 
            if(isInDetector){
              if(is511) //energy before process
              {
                h511Sygnal.Fill(smearEnergy511); // 1 rozproszenie w detectorze/rejestracja fotonu w detektorze
              }
              else
              {
                h511B1B2.Fill(smearEnergy511);
                h511B3.Fill(smearEnergy511); //n rozproszenie w detektorze
              }
           } else if (isInPhantom || isInPhantom2)
            {h511B4.Fill(hit.fEnergyBeforeProcess);
              if(is511)
              {
                h511B1B2.Fill(smearEnergy511);
                h511B1.Fill(smearEnergy511); //1 w fantomie n w det
                //h511B4.Fill(hit.fEnergyBeforeProcess);
              }
              else
              {
                h511B1B2.Fill(smearEnergy511);
                h511B2.Fill(smearEnergy511); //n w fantomie 0 w det
              }
            } else
            {
             hostalos.Fill((hit.fEnergyDeposition)); 
            }

          } else {
            if ( isEqual(emissionEnergy, 511 )) {
            smearEnergy5112 = smearEnergy(hit.fEnergyDeposition);
            gamma2 = TLorentzVector(hit.fHitPosition, smearEnergy5112);
              gamma511Pos2.push_back(gamma2);
              x5112[treevent] = hit.fHitPosition.x();
            y5112[treevent] = hit.fHitPosition.y();
            z5111[treevent] = hit.fHitPosition.z();
              cout << " " << endl;
              cout << "Second gamma-511" << endl; 
            cout << "energy deposition= " << hit.fEnergyDeposition << endl;
            cout << "emission energy = " << emissionEnergy << endl;
            cout << "number of hit = " << rozp << endl;
            cout << "energy = " << gamma511Pos2[gamma511Pos2.size()-1].Energy() << endl;
            cout << "liczba = " << gamma511Pos2.size() << endl;
            h511SygnalB1B2.Fill(smearEnergy5112); //zapisujemy wszystkie 511 
              if(isInDetector){
              if(is511)
              {
                h511Sygnal.Fill(smearEnergy5112); // 1 rozproszenie w detectorze
              }
              else
              {
                h511B1B2.Fill(smearEnergy5112);
                h511B3.Fill(smearEnergy5112); //n rozproszenie w detektorze
              }
            }else if (isInPhantom || isInPhantom2)
            {h511B4.Fill(hit.fEnergyBeforeProcess);
              if(is511)
              {
                h511B1B2.Fill(smearEnergy5112);
                h511B1.Fill(smearEnergy5112); //1 w fantomie n w det
                              }
              else
              {
                h511B1B2.Fill(smearEnergy5112);
                h511B2.Fill(smearEnergy5112); //n w fantomie 0 w det

              }
            } else
            {
             hostalos.Fill(smearEnergy5112); 
            }
            } else {
              if ( isEqual(emissionEnergy, 1157)) {
                smearEnergyPrompt = smearEnergy(hit.fEnergyDeposition);
                gammaPrompt = TLorentzVector(hit.fHitPosition, smearEnergyPrompt);
                gammaPromptPos.push_back(gammaPrompt);
                xprompt[treevent] = hit.fHitPosition.x();
            yprompt[treevent] = hit.fHitPosition.y();
            zprompt[treevent] = hit.fHitPosition.y();
                cout << " " << endl;
                cout << "Gamma-prompt" << endl; 
            cout << "energy deposition= " << hit.fEnergyDeposition << endl;
            cout << "emission energy = " << emissionEnergy << endl;
            cout << "number of hit = " << rozp << endl;
            cout << "energy = " << gammaPromptPos[gammaPromptPos.size()-1].Energy() << endl;
            cout << "liczba = " << gammaPromptPos.size() << endl;
            hpromptSygnalB1B2.Fill(smearEnergyPrompt); //zapisujemy wszystkie 511 
                if(isInDetector){

              if(isprompt)
              {
                hpromptSygnal.Fill(smearEnergyPrompt); // 1 rozproszenie w detectorze

              }
              else
              {
                hpromptB1B2.Fill(smearEnergyPrompt);
                hpromptB3.Fill(smearEnergyPrompt); //n rozproszenie w detektorze
                
              }
           } else if (isInPhantom || isInPhantom2)
            {
              
              if(isprompt)
              {
                hpromptB1B2.Fill(smearEnergyPrompt);
                hpromptB1.Fill(smearEnergyPrompt); //1 w fantomie n w det
                
              }
              else
              {
                hpromptB1B2.Fill(smearEnergyPrompt);
                hpromptB2.Fill(smearEnergyPrompt); //n w fantomie 0 w det
                
              }
             }else
            {
             hostalos.Fill((hit.fEnergyDeposition)); 
            }
              } else {
                std::cerr << "This should never happen !" << std::endl;
                assert(1 == 0);
              }
            }
          }
        }
      }}
    }
    
  }
  cout << "Double coincidence = " << dwaevent << endl;
  cout << "True double coincidence = " << truedwaevent << endl;
  cout << "Onie hit in event = " << eve1 << endl;
  cout << "More than three hits in event = " << evemore << endl;
  auto xy511_1 = new TGraph(treevent, x5111, y5111);
  auto xy511_2 = new TGraph(treevent, x5112, y5112);
  auto xyprompt = new TGraph(treevent, xprompt, yprompt);

  auto xz511_1 = new TGraph(treevent, x5111, z5111);
  auto xz511_2 = new TGraph(treevent, x5112, z5112);
  auto xzprompt = new TGraph(treevent, xprompt, zprompt);

  auto zy511_1 = new TGraph(treevent, z5111, y5111);
  auto zy511_2 = new TGraph(treevent, z5112, y5112);
  auto zyprompt = new TGraph(treevent, zprompt, yprompt);

  TCanvas c120("c", "Positions", 600, 600);
  xy511_1->SetLineColor(kBlack);
  xy511_1->SetTitle("Positions XY");
  xy511_1->GetXaxis()->SetTitle("X");
  xy511_1->GetYaxis()->SetTitle("Y");
  xy511_1->Draw("ap");
  xy511_2->Draw("ap");
  xyprompt->Draw("ap");
  c120.SaveAs("positions.png");

  TCanvas c121("c", "Positions", 600, 600);
  xz511_1->SetLineColor(kBlack);
  xz511_1->SetTitle("Positions XZ");
  xz511_1->GetXaxis()->SetTitle("X");
  xz511_1->GetYaxis()->SetTitle("Y");
  xz511_1->Draw("ap");
  xz511_2->Draw("ap");
  xzprompt->Draw("ap");
  c121.SaveAs("positionsXZ.png");

  TCanvas c122("c", "Positions", 600, 600);
  zy511_1->SetLineColor(kBlack);
  zy511_1->SetTitle("Positions YZ");
  zy511_1->GetXaxis()->SetTitle("X");
  zy511_1->GetYaxis()->SetTitle("Y");
  zy511_1->Draw("ap");
  zy511_2->Draw("ap");
  zyprompt->Draw("ap");
  c122.SaveAs("positionsYZ.png");

  //------------------------------------------------------------------------ENERGIA-------------------------------------------------------------------
  const double kEnergyMin = 0;
  const double kEnergyMax = 1157;
  const int kNumberOfSteps = 1157;
  const double kEnergyStep = (kEnergyMax - kEnergyMin) / kNumberOfSteps;
  double Emin;
  double Ecut;
  int k_save =0;
  int i_save=0;
  
  Double_t TPR_prompt[kNumberOfSteps]; //true positive rate (efficirncy/wydajność)
  Double_t PPV_prompt[kNumberOfSteps]; //positive predictive value (purity/czystość)
  Double_t FPR_prompt[kNumberOfSteps]; //false positive rate (false alarm)
  Double_t TPR_511[kNumberOfSteps];
  Double_t PPV_511[kNumberOfSteps];
  Double_t FPR_511[kNumberOfSteps];
  Double_t energies_treshold[kNumberOfSteps];
  Double_t energies_cut[kNumberOfSteps];
  Double_t treshold_511[kNumberOfSteps];
  Double_t treshold_prompt[kNumberOfSteps];
  Double_t F_511[kNumberOfSteps];
  Double_t F_511_manual[kNumberOfSteps];
  Double_t F_prompt[kNumberOfSteps];
  Double_t F_prompt_manual[kNumberOfSteps];

  std::vector<Double_t> Emin_save;
  std::vector<Double_t> F_511_save;
  std::vector<Double_t> F_prompt_save;
  std::vector<Double_t> F_511_save2;
  vector<Double_t>::iterator result_F_511;
  vector<Double_t>::iterator result_F_511_end;

  std::vector<LOR> result;
  //F_511_save.clear();

  //ograniczamy maksymalną ilość cięcia energetycznego i tresholda
  const int kMaxEnergyCutStep = 201;
  int resultnumsave0;
  int resultnumsave1;
  int licz=0;
  Double_t klas0save[kMaxEnergyCutStep];
  Double_t klas1save[kMaxEnergyCutStep];
  Double_t klas2save[kMaxEnergyCutStep];
  Double_t klas3save[kMaxEnergyCutStep];
  
  int minBinSygnalPrompt = hpromptSygnal.GetXaxis()->FindBin(kEnergyMin); /// min bin for Prompt
  int maxBinSygnalPrompt = hpromptSygnal.GetXaxis()->FindBin(kEnergyMax); /// max bin for Prompt
  int minBinSygnal511 = h511Sygnal.GetXaxis()->FindBin(kEnergyMin); // min bin for 511
  int maxBinSygnal511 = h511Sygnal.GetXaxis()->FindBin(kEnergyMax); /// max bin for 511
  int maxBinAllPrompt = hpromptSygnalB1B2.GetXaxis()->FindBin(kEnergyMax); /// max bin for Prompt
  int maxBinAll511 = h511SygnalB1B2.GetXaxis()->FindBin(kEnergyMax); /// max bin for 511
  int maxBinB1B2Prompt = hpromptB1B2.GetXaxis()->FindBin(kEnergyMax); /// max bin for Prompt
  
  //liczba wszystkich gamma-prompt i 511 z jednym zdarzeniem detektorze
  double allPrompt = hpromptSygnal.Integral(minBinSygnalPrompt, maxBinSygnalPrompt);
  double all511 = h511Sygnal.Integral(minBinSygnal511, maxBinSygnal511);
  
  //---------------------------------------ENERGY CUT-----------------------------------------------------------------------------------------
  
  for (Int_t k = 0; k < kMaxEnergyCutStep; k++){

    double currentEnergy_cut =  kEnergyStep * k;
    energies_cut[k] = currentEnergy_cut;
    //cout << "k: " << k <<endl;

    //przeliczamy aktualne biny zmienione po zmianie energy cut
    int currentBinSygnal511_cut = h511Sygnal.GetXaxis()->FindBin(currentEnergy_cut);
    //int currentBinSygnalPrompt_cut = hpromptSygnal.GetXaxis()->FindBin(currentEnergy_cut);
    int currentBinAll511_cut = h511SygnalB1B2.GetXaxis()->FindBin(currentEnergy_cut);
    int currentBinAllPrompt_cut = hpromptSygnalB1B2.GetXaxis()->FindBin(currentEnergy_cut);
    int currentBinB1B2511_cut = h511B1B2.GetXaxis()->FindBin(currentEnergy_cut);
    int currentBinB1B2Prompt_cut = hpromptB1B2.GetXaxis()->FindBin(currentEnergy_cut);

    for (int zero=0; zero<kEnergyStep; zero++)
    {
      TPR_511[zero]=0;
      TPR_prompt[zero]=0;
      PPV_prompt[zero]=0;
      FPR_prompt[zero]=0; 
      PPV_511[zero]=0;
      FPR_511[zero]=0;
      energies_treshold[zero]=0;
      //energies_cut[zero]=0;
      treshold_511[zero]=0;
      treshold_prompt[zero]=0;
      F_511[zero]=0;
      F_511_manual[zero]=0;
      F_prompt[zero]=0;
      F_prompt_manual[zero]=0;
    }
    F_511_save.clear();

    //wyliczamy treshold
    i_save=0;
    for (Int_t i = k; i < kNumberOfSteps ; i++) {
      double currentEnergy =  kEnergyStep * i;
      energies_treshold[i] = currentEnergy;

      //przeliczamy aktualne biny zmienione po zmianie treshold
      int currentBinSygnal511 = h511Sygnal.GetXaxis()->FindBin(currentEnergy);
      int currentBinSygnalPrompt = hpromptSygnal.GetXaxis()->FindBin(currentEnergy);
      int currentBinAll511 = h511SygnalB1B2.GetXaxis()->FindBin(currentEnergy);
      int currentBinAllPrompt = hpromptSygnalB1B2.GetXaxis()->FindBin(currentEnergy);
      int currentBinB1B2511 = h511B1B2.GetXaxis()->FindBin(currentEnergy);
      int currentBinB1B2Prompt = hpromptB1B2.GetXaxis()->FindBin(currentEnergy);

      //for prompt
      double allAcceptedAsPrompt = hpromptSygnalB1B2.Integral(currentBinAllPrompt, maxBinAllPrompt) + h511SygnalB1B2.Integral(currentBinAll511, maxBinAll511);
      double promptAcceptedAsPrompt = hpromptSygnal.Integral(currentBinSygnalPrompt, maxBinSygnalPrompt);
      double remainingAcceptedAsPrompt = hpromptB1B2.Integral(currentBinB1B2Prompt, maxBinB1B2Prompt)+h511SygnalB1B2.Integral(currentBinAll511, maxBinAll511);

      //for 511
      double allAcceptedAs511 = hpromptSygnalB1B2.Integral(currentBinAllPrompt_cut, currentBinAllPrompt) + h511SygnalB1B2.Integral(currentBinAll511_cut, currentBinAll511);
      double v511AcceptedAs511 = h511Sygnal.Integral(currentBinSygnal511_cut, currentBinSygnal511);
      double remainingAcceptedAs511 = hpromptSygnalB1B2.Integral(currentBinAllPrompt_cut, currentBinAllPrompt) + h511B1B2.Integral(currentBinB1B2511_cut, currentBinB1B2511);

      //purity prompt
      if (allAcceptedAsPrompt !=0) 
      {
        PPV_prompt[i] = promptAcceptedAsPrompt/allAcceptedAsPrompt; 
      } else {
        PPV_prompt[i] = 0;
      }
      
      //efficiency ptompt
      TPR_prompt[i] = promptAcceptedAsPrompt/allPrompt; 

      /// This is the porcentage of 511 that is classifed falsely as prompt
      FPR_prompt[i] = remainingAcceptedAsPrompt / all511;

      //purity 511
      
      if (allAcceptedAs511 !=0) {
        PPV_511[i] = v511AcceptedAs511 /allAcceptedAs511;
      } else {
        PPV_511[i] = 0;
      }
      //std::cout << "v511AcceptedAs511 " << v511AcceptedAs511 << " allAcceptedAs511 " << allAcceptedAs511  << " PPV_511[i] " << PPV_511[i] <<std::endl;

      //efficiency 511
      TPR_511[i] = v511AcceptedAs511 / all511;
      //std::cout << "TPR511 " << TPR_511[i] <<std::endl;
      /// This is the porcentage of prompt that is classifed falsely as 511

      FPR_511[i] = remainingAcceptedAs511  / allPrompt;
      //std::cout << "FPR_511 " << FPR_511[i] <<std::endl;

      if((PPV_511[i]+TPR_511[i]!=0)){
        F_511[i]=(2*PPV_511[i]*TPR_511[i])/(PPV_511[i]+TPR_511[i]);
      } else {
        F_511[i]=0;
      }
      if ((PPV_prompt[i]+TPR_prompt[i])!=0){
        F_prompt[i]=(2*PPV_prompt[i]*TPR_prompt[i])/(PPV_prompt[i]+TPR_prompt[i]);
      }else{
        F_prompt[i]=0;
      }

    //cout << " "  <<endl; 
    //cout << "F_511 " << F_511[i] <<endl; 
    //cout << "i " << i <<endl; 
          //--------------MANUAL------------------------
      if (F_511[i]>=F_511[i_save] || isEqual(i,0))
      {
        treshold_511[k]=energies_treshold[i];
        F_511_manual[k]=F_511[i];
        i_save=i;
        //cout << "Energy step dla F_511[i]>F_511[i-1]: " << F_511_test[k] <<endl;               
      }

      if (F_prompt[i]>F_prompt[i-1])
      {
        treshold_prompt[k]=energies_treshold[i];
        F_prompt_manual[k]=F_prompt[i];
        
      }

      //---------------------------------AUTO------------------
      F_511_save.push_back(F_511[i]);
      F_prompt_save.push_back(F_prompt[i]);
    } // END OF FOR------------------------------------------------------------------------------------------------------------
    Emin=treshold_511[k];

    if (isEqual(k,0)){
    
    efficiency511Ecut0 = new TGraph(kNumberOfSteps, energies_treshold, TPR_511);
    purity511Ecut0 = new TGraph(kNumberOfSteps, energies_treshold, PPV_511);
    ROC511Ecut0 = new TGraph(kNumberOfSteps, FPR_511, TPR_511);
    F1_511Ecut0 = new TGraph(kNumberOfSteps, energies_treshold, F_511);

    ROCpromptEcut0 = new TGraph(kNumberOfSteps, FPR_prompt, TPR_prompt);
    efficiencyPromptEcut0 = new TGraph(kNumberOfSteps, energies_treshold, TPR_prompt);
    purityPromptEcut0 = new TGraph(kNumberOfSteps, energies_treshold, PPV_prompt);
    F1PromptEcut0 = new TGraph(kNumberOfSteps, energies_treshold, F_prompt);

    TCanvas c12("c", "Purity and efficiency for gamma-prompt with energy cut = 0", 2000, 1200);
  efficiencyPromptEcut0->SetLineColor(kBlack);
  efficiencyPromptEcut0->SetTitle("Purity and efficiency for gamma-prompt with energy cut = 0");
  efficiencyPromptEcut0->GetXaxis()->SetTitle("Energy threshold value [keV]");
  efficiencyPromptEcut0->GetYaxis()->SetTitle("PPV/TPR");
  efficiencyPromptEcut0->GetYaxis()->SetRange(0,1.2);
  efficiencyPromptEcut0->Draw();
  purityPromptEcut0->SetLineColor(kRed);
  purityPromptEcut0->Draw("same");
  TLegend legend(0.7,0.75,0.9,0.9);
  legend.AddEntry(purityPromptEcut0,"purity","l");
  legend.AddEntry(efficiencyPromptEcut0,"efficiency","l");
  legend.Draw();
  c12.SaveAs("purity_efficiency_prompt_0.png");

  TCanvas c13("c", "c", 2000, 1200);
  efficiency511Ecut0->SetLineColor(kBlack);
  efficiency511Ecut0->SetTitle("Purity and efficiency for 511 keV with energy cut = 0");
  efficiency511Ecut0->GetXaxis()->SetTitle("Energy threshold value [keV]");
  efficiency511Ecut0->GetYaxis()->SetTitle("PPV/TPR");
  efficiency511Ecut0->GetYaxis()->SetRange(0,1.2);
  efficiency511Ecut0->Draw();
  purity511Ecut0->SetLineColor(kRed);
  purity511Ecut0->Draw("same");
  TLegend legend2(0.7,0.75,0.9,0.9);
  legend2.AddEntry(purity511Ecut0,"purity","l");
  legend2.AddEntry(efficiency511Ecut0,"efficiency","l");
  legend2.Draw();
  c13.SaveAs("purity_efficiency_511_0.png");

  TCanvas c16("c", "c", 1000, 1000);
  ROCpromptEcut0->SetLineColor(kBlack);
  ROCpromptEcut0->SetTitle("ROC for gamma-prompt and 511 keV with energy cut = 0");
  ROCpromptEcut0->GetXaxis()->SetTitle("FPR");
  ROCpromptEcut0->GetYaxis()->SetTitle("TPR");
  ROCpromptEcut0->GetXaxis()->SetRangeUser(0,1);
  //ROCpromptEcut0->GetYaxis()->SetRange(0,1);
  ROCpromptEcut0->Draw();
  ROC511Ecut0->SetLineColor(kRed);
  ROC511Ecut0->Draw("same");
  TLegend legend3(0.7,0.75,0.9,0.9);
  legend3.AddEntry(ROCpromptEcut0,"gamma-prompt","l");
  legend3.AddEntry(ROC511Ecut0,"2-511 KeV","l");
  legend3.Draw();
  c16.SaveAs("ROC_prompt_511_0.png");

  TCanvas c17("c", "c", 2000, 1200);
  F1_511Ecut0->SetLineColor(kBlack);
  F1_511Ecut0->SetTitle("F1 score for gamma-prompt and 511 keV with energy cut = 0");
  F1_511Ecut0->GetXaxis()->SetTitle("Energy treshold value [keV]");
  F1_511Ecut0->GetYaxis()->SetTitle("F1 [0-1]");
  //F1_511Ecut0->GetYaxis()->SetRange(0,1.2);
  F1_511Ecut0->Draw();
  F1PromptEcut0->SetLineColor(kRed);
  F1PromptEcut0->Draw("same");
  TLegend legend4(0.7,0.75,0.9,0.9);
  legend4.AddEntry(F1PromptEcut0,"gamma-prompt","l");
  legend4.AddEntry(F1_511Ecut0,"2-511 KeV","l");
  legend4.Draw();
  c17.SaveAs("F1_511Ecut0.png");
    }

    // ------------------------------------------------------------DISTANCE----------------------------------------------------------------------------------------
//inicjalizacja zmiennych
  int numOfEvents = gamma511Pos1.size();
  double distgamma1prompt[numOfEvents], distgamma2prompt[numOfEvents], distgamma1gamma2[numOfEvents];
  Counter counter;

  std::vector<int> resultNum = {0,0,0,0,};
  //resultNum.clear();
//główna pętla for
  for (Int_t i = 0; i < numOfEvents; i++) {
    
    //wyliczamy długość pomiędzy wskazanymi cząstkami (511-511 albo 511-prompt)
    distgamma1prompt[i] = calculateDistance3D(gammaPromptPos[i].X(), gammaPromptPos[i].Y(), gamma511Pos2[i].X(), gamma511Pos2[i].Y(), gammaPromptPos[i].Z(), gamma511Pos2[i].Z());
    distgamma2prompt[i] = calculateDistance3D(gammaPromptPos[i].X(), gammaPromptPos[i].Y(), gamma511Pos1[i].X(), gamma511Pos1[i].Y(), gammaPromptPos[i].Z(), gamma511Pos1[i].Z());
    distgamma1gamma2[i] = calculateDistance3D(gamma511Pos1[i].X(), gamma511Pos1[i].Y(), gamma511Pos2[i].X(), gamma511Pos2[i].Y(), gamma511Pos1[i].Z(), gamma511Pos2[i].Z());
    //wsadzamy długość do histogramu
    hDISTANCEgamma1prompt.Fill(distgamma1prompt[i]);
    hDISTANCEgamma2prompt.Fill(distgamma2prompt[i]);
    hDISTANCEgamma1gamma2.Fill(distgamma1gamma2[i]);

    //wyliczamy największy/najmniejszy distance
    auto maxVal = std::max({distgamma1gamma2[i], distgamma1prompt[i], distgamma2prompt[i]  });
    auto minVal = std::min({distgamma1gamma2[i], distgamma1prompt[i], distgamma2prompt[i]  });

    //wyluczamy ile jest 511-prompt i 511-511 z największym distancem
    if (isEqual(maxVal, distgamma1prompt[i]) || isEqual(maxVal, distgamma2prompt[i])) {
      counter.numericprompt++;
    } else {
      counter.numericgamma++;
    }

    //sprawdzamy do jakiej długości należy najwiecej par 511-511
      if (isEqual(maxVal, distgamma1gamma2[i])) {
      counter.numeric_dist_gamma1gamma2_max++;
      hgammaenergyMax.Fill(gamma511Pos1[i].E());
      hgammaenergyMax.Fill(gamma511Pos2[i].E());
    } else if (isEqual(minVal, distgamma1gamma2[i])) {
      counter.numeric_dist_gamma1gamma2_min++;
      hgammaenergyMin.Fill(gamma511Pos1[i].E());
      hgammaenergyMin.Fill(gamma511Pos2[i].E());
    } else {
      counter.numeric_dist_gamma1gamma2_mid++;
      hgammaenergyMid.Fill(gamma511Pos1[i].E());
      hgammaenergyMid.Fill(gamma511Pos2[i].E());
    }

    //sprawdzamy do jakiej długości należy najwiecej par 511-prompt
    if (isEqual(maxVal, distgamma1prompt[i])) {
      counter.numeric_dist_gamma1prompt_max++;
      hgammaenergyMax.Fill(gamma511Pos1[i].E());
      hgammapromptenergyMax.Fill(gammaPromptPos[i].E());
    } else if (isEqual(minVal, distgamma1prompt[i])) {
      hgammaenergyMin.Fill(gamma511Pos1[i].E());
      hgammapromptenergyMin.Fill(gammaPromptPos[i].E());
    } else {
      hgammaenergyMid.Fill(gamma511Pos1[i].E());
      hgammapromptenergyMid.Fill(gammaPromptPos[i].E());
    }

    //sprawdzamy do jakiej długości należy najwiecej par 511-prompt
    if (isEqual(maxVal, distgamma2prompt[i])) {
      counter.numeric_dist_gamma2prompt_max++;
      hgammaenergyMax.Fill(gamma511Pos2[i].E());
      hgammapromptenergyMax.Fill(gammaPromptPos[i].E());
    } else if (isEqual(minVal, distgamma2prompt[i])) {
      hgammaenergyMin.Fill(gamma511Pos2[i].E());
      hgammapromptenergyMin.Fill(gammaPromptPos[i].E());
    } else {
      hgammaenergyMid.Fill(gamma511Pos2[i].E());
      hgammapromptenergyMid.Fill(gammaPromptPos[i].E());
    }

//-------------------------------------------------------------------------ALGORYTHM------------------------------------------------------------------------------
    //oznaczamy prawidłowy LOR (511-511)
    LOR trueLOR = {gamma511Pos1[i], gamma511Pos2[i]};

    /*cout << "e1 " << gamma511Pos1[i].Energy() <<endl;
    cout << "e2 " << gamma511Pos2[i].Energy()<<endl;

    cout << "  " <<i<<endl;*/
    
    //włączamy algorythm selekcji
    result = select2(gamma511Pos1[i], gamma511Pos2[i], gammaPromptPos[i], treshold_511[k], k);
      
    //Statystyka wyników:
    //żadny LOR nie został zaakceptowany prawidłowo
    if (result.empty()) {
      resultNum[2]++;
    } else {
      //2 LORy zostały uznane za prawidłowy
      if (result.size() == 2) {
        resultNum[3]++;
      } else {
        auto lor = result[0];
        //algorythm znalazł prawidłowy LOR
        if (isEqualLor(trueLOR, lor)) {
          resultNum[0]++;
        } else {
          //pozostałe błędne przypadki
          resultNum[1]++;
        }
      }
    }

  }


cout << "licz " << licz <<endl;
  //procent zaakceptowanych 511-511 jako prompt
  double percent = counter.numericgamma * 100. / numOfEvents;
  
  if ((resultNum[0]>=resultnumsave0 && resultNum[1]<=resultnumsave1)|| isEqual(k,0))
      {
        resultnumsave0=resultNum[0];
        resultnumsave1=resultNum[1];
        cout << "res 0 " << resultnumsave0 <<endl;
        cout << "res 1 " << resultnumsave1 <<endl;
        k_save=k;              
    } 
        printResults(numOfEvents, resultNum, Emin, k);
        klas0save[k]=resultNum[0];
        klas1save[k]=resultNum[1];
        klas2save[k]=resultNum[2];
        klas3save[k]=resultNum[3];

  resultNum = {0,0,0,0,};
  cout << "K save " << k_save <<endl;
      
  }
  klas0 = new TGraph(kMaxEnergyCutStep, energies_cut, klas0save);
  klas1 = new TGraph(kMaxEnergyCutStep, energies_cut, klas1save);
  klas2 = new TGraph(kMaxEnergyCutStep, energies_cut, klas2save);
  klas3 = new TGraph(kMaxEnergyCutStep, energies_cut, klas3save);


  TCanvas c18("c", "c", 2000, 1200);
    //c18.SetLogy();
  klas2->SetLineColor(kBlack);
  klas2->SetTitle("Algorithm results depending on energy cut value");
  klas2->GetXaxis()->SetTitle("Energy cut value");
  klas2->GetYaxis()->SetTitle("Number of events");
  klas2->GetYaxis()->SetRangeUser(0,3500);
  klas2->Draw();
  klas1->SetLineColor(kRed);
  klas1->Draw("same");
  klas0->SetLineColor(kBlue);
  klas0->Draw("same");
  klas3->SetLineColor(kGray);
  klas3->Draw("same");
  TLegend legend0(0.7,0.75,0.9,0.9);
  legend0.AddEntry(klas0,"True LOR accepted","l");
  legend0.AddEntry(klas1,"False LOR accepted","l");
  legend0.AddEntry(klas2,"all LORs rejected","l");
  legend0.AddEntry(klas3,"all LORs accepted","l");
  legend0.Draw();
  c18.SaveAs("klasy.png");
    //konic for K
  //printevents(numOfEvents, counter);

  /*F1_max_511 = new TGraph(kMaxEnergyCutStep, energies_cut, gF_511);
  F1_max_prompt = new TGraph(kMaxEnergyCutStep, energies_cut, gFprompt);

  TCanvas c11("c", "c", 2000, 1200);
  F1_max_prompt->GetYaxis()->SetRangeUser(0,1);
  //gr->GetXaxis()->SetRangeUser(0.2,0.4);
  F1_max_prompt->SetLineColor(kBlack);
  F1_max_prompt->SetTitle("Max F1 score for gamma-prompt and 511 keV for different energy cut");
  F1_max_prompt->GetXaxis()->SetTitle("Energy cut value [keV]");
  //F1_max_prompt->SetMarkerStyle(21);
  //F1_max->GetXaxis()->SetRange(0,34.1);
  F1_max_prompt->GetYaxis()->SetTitle("F1 [0-1]");
  
  F1_max_prompt->Draw();
  F1_max_511->SetLineColor(kRed);
  F1_max_511->Draw("same");
  TLegend legend0(0.7,0.75,0.9,0.9);
  legend0.AddEntry(F1_max_prompt,"gamma-prompt","l");
  legend0.AddEntry(F1_max_511,"2-511 keV","l");
  legend0.Draw();
  c11.SaveAs("F1_max.png");*/



  //procent zaakceptowanych 511-511 jako prompt
  //double percent = counter.numericgamma * 100. / numOfEvents;

//------------------------------------------------------------------DRUKOWANIE WYNIKÓW-------------------------------------------------------------------------------
  //printResults(numOfEvents, counter, percent, resultNum, Emin, Ecut);

  //zapisujemy histogramy wynikowe w pliku
    testOut.cd();
  histMultiplicity.Write();
  histMultiplicity2.Write();
  hgammaenergyMax.Write();
  hgammapromptenergyMax.Write();
  hgammaenergyMid.Write();
  hgammapromptenergyMid.Write();
  hgammaenergyMin.Write();
  hgammapromptenergyMin.Write();
  hDISTANCEgamma1gamma2.Write();
  hDISTANCEgamma2prompt.Write();
  hDISTANCEgamma1prompt.Write();

  //--------------------------------------------------------------Wyniki w postaci odrazków--------------------------------------------------------------------
  TCanvas c1("c", "c", 2000, 1200);
  //hpromptSygnalB1B2.SetStats(kFALSE);
  h511SygnalB1B2.SetStats(kFALSE);
  h511SygnalB1B2.SetTitle("All types of 511 keV");
  h511SygnalB1B2.GetXaxis()->SetTitle("Energy [keV]");
  h511SygnalB1B2.GetYaxis()->SetTitle("Events");
  h511SygnalB1B2.SetLineColor(kGreen);
  h511SygnalB1B2.Draw();
  h511Sygnal.SetLineColor(kBlack);
  h511Sygnal.Draw("same");
  h511B1.SetLineColor(kRed);
  h511B1.Draw("same");
  h511B2.SetLineColor(kBlue);
  h511B2.Draw("same");
  h511B3.SetLineColor(kGray);
  h511B3.Draw("same");
  TLegend legend260(0.7,0.75,0.9,0.9);
  legend260.AddEntry("h511SygnalB1B2","All events","l");
  legend260.AddEntry("h511Sygnal","0 phantom && 1 detector","l");
  legend260.AddEntry("h511B1","1 phantom && 0 detector","l");
  legend260.AddEntry("h511B2","n phantom && 0 detector","l");
  legend260.AddEntry("h511B3","n phantom && n detector","l");
  legend260.Draw();
  c1.SaveAs("all511.png");

  TCanvas c2("c", "c", 2000, 1200);
  hpromptSygnalB1B2.SetStats(kFALSE);
  hpromptSygnalB1B2.SetTitle("All types of gamma-prompt");
  hpromptSygnalB1B2.GetXaxis()->SetTitle("Energy [keV]");
  hpromptSygnalB1B2.GetYaxis()->SetTitle("Events");
  hpromptSygnalB1B2.SetLineColor(kBlue);
  hpromptSygnalB1B2.Draw();
  hpromptSygnal.SetLineColor(kRed);
  hpromptSygnal.Draw("same");
  hpromptB1.SetLineColor(kBlack);
  hpromptB1.Draw("same");
  hpromptB2.SetLineColor(kGreen);
  hpromptB2.Draw("same");
   hpromptB3.SetLineColor(kGray);
  hpromptB3.Draw("same");
  TLegend legend_allprompt(0.7,0.75,0.9,0.9);
  legend_allprompt.AddEntry("hpromptSygnalB1B2","All events","l");
  legend_allprompt.AddEntry("hpromptSygnal","0 phantom && 1 detector","l");
  legend_allprompt.AddEntry("hpromptB1","1 phantom && 0 detector","l");
  legend_allprompt.AddEntry("hpromptB2","n phantom && 0 detector","l");
  legend_allprompt.AddEntry("hpromptB3","n phantom && n detector","l");
  legend_allprompt.Draw();
  c2.SaveAs("allprompt.png");

  TCanvas c9("c", "c", 2000, 1200);
  h511SygnalB1B2.SetStats(kFALSE);
  h511SygnalB1B2.SetTitle("Histograms for calculation purity and efficiency for gamma-prompt");
  h511SygnalB1B2.GetXaxis()->SetTitle("Energy [keV]");
  h511SygnalB1B2.GetYaxis()->SetTitle("Events");
  h511SygnalB1B2.SetLineColor(kRed);
  h511SygnalB1B2.Draw();
  hpromptSygnalB1B2.SetStats(kFALSE);
  hpromptSygnalB1B2.SetLineColor(kBlack);
  hpromptSygnalB1B2.Draw("same");
TLegend legend_calcprompt(0.7,0.75,0.9,0.9);
  legend_calcprompt.AddEntry("h511SygnalB1B2","All 2-gamma 511 keV","l");
  legend_calcprompt.AddEntry("hpromptSygnalB1B2","All gamma-prompt","l");
  legend_calcprompt.Draw();
  c9.SaveAs("calculation_purityprompt.png");

TCanvas c10("c", "c", 2000, 1200);
  h511B1.SetStats(kFALSE);
  h511B1.SetTitle("Histograms for calculation purity and efficiency for 511 keV");
  h511B1.GetXaxis()->SetTitle("Energy [keV]");
  h511B1.GetYaxis()->SetTitle("Events");
  c10.SetLogy();
  h511B1.SetLineColor(kRed);
  h511B1.Draw();
  h511Sygnal.SetStats(kFALSE);
  h511Sygnal.SetLineColor(kBlue);
  h511Sygnal.Draw("same");
  h511B2.SetLineColor(kGreen);
  h511B2.Draw("same");
  h511B3.SetLineColor(kGray);
  h511B3.Draw("same");
  hpromptSygnalB1B2.SetLineColor(kBlack);
  hpromptSygnalB1B2.Draw("same");
TLegend legend_calc511(0.7,0.75,0.9,0.9);
  legend_calc511.AddEntry("hpromptSygnalB1B2","all prompt","l");
  legend_calc511.AddEntry("h511Sygnal","0 phantom && 1 detector","l");
  legend_calc511.AddEntry("h511B1","1 phantom && 0 detector","l");
  legend_calc511.AddEntry("h511B2","n phantom && 0 detector","l");
  legend_calc511.AddEntry("h511B3","n phantom && n detector","l");
  legend_calc511.Draw();
  c10.SaveAs("calculation_purity511.png");

  TCanvas c3("c", "c", 2000, 1200);
  hpromptSygnal.SetTitle("Sygnal prompt: 0 fantom && 1 detector");
  hpromptSygnal.GetXaxis()->SetTitle("Energy [keV]");
  hpromptSygnal.GetYaxis()->SetTitle("Events");
  hpromptSygnal.SetLineColor(kBlack);
  hpromptSygnal.Draw();
  c3.SaveAs("sygnal_prompt.png");

  TCanvas c6("c", "c", 2000, 1200);
  h511Sygnal.SetTitle("Sygnal 511 keV: 0 fantom && 1 detector");
  h511Sygnal.GetXaxis()->SetTitle("Energy [keV]");
  h511Sygnal.GetYaxis()->SetTitle("Events");
  h511Sygnal.SetLineColor(kBlack);
  h511Sygnal.Draw();
  c6.SaveAs("sygnal_511.png");

  TCanvas c60("c", "c", 2000, 1200);
  h511B4.SetTitle("B4");
  h511B4.GetXaxis()->SetTitle("Energy [keV]");
  h511B4.GetYaxis()->SetTitle("Events");
  h511B4.SetLineColor(kBlack);
  h511B4.Draw();
  c60.SaveAs("B4.png");

  h511B1B2.Write();
  hpromptB1B2.Write();
  hpromptSygnalB1B2.Write();
  h511Sygnal.Write();
  hpromptSygnal.Write();
  h511SygnalB1B2.Write();
  h511B2.Write();
  h511B1.Write();
  h511B3.Write();
  hpromptB3.Write();
  hpromptB2.Write();
  hpromptB1.Write();
  hostalos.Write();
  hDwa.Write();
  hDwaFantom.Write();
  hDwaAll.Write();

  c1.Write();
  c2.Write();
  c3.Write();
  //c4.Write();
  //c5.Write();
  c6.Write();
  c60.Write();
  //c8.Write();
  c9.Write();
  c10.Write();
}

double calculateDistance(const LOR& lor)
{
  return calculateDistance3D(lor.first.X(), lor.first.Y(), lor.first.Z(), lor.second.X(), lor.second.Y(), lor.second.Z());
}

bool isInEnergyRange(double E, double Emin, double Ecut)
{

  return ((E >= Ecut) && (E <= Emin));
}

bool isEqualLor(const LOR& lor1, const LOR& lor2)
{
  return ((lor1.first == lor2.first)  && (lor1.second == lor2.second)) || ((lor1.first == lor2.second)  && (lor1.second == lor2.first));
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

//kryterium energetyczny
std::vector<LOR> select2(const TLorentzVector& gamma1, const TLorentzVector& gamma2, const TLorentzVector& gamma3, double Emin, double Ecut)
{
  //są 3 LORa
  LOR lor1 = {gamma1, gamma2};
  LOR lor2 = {gamma1, gamma3};
  LOR lor3 = {gamma2, gamma3};

  std::vector<LOR> selectedAfterGeomCut;
  
  //sortujemy po długości
  //cout << " " <<endl;
  selectedAfterGeomCut = geometry(lor1, lor2, lor3);
  int licznik =0;
  int licznik2 =0;
  std::vector<LOR> finalSelection;
  //kryteruim energetyczny
  for (auto& lor : selectedAfterGeomCut ) {
    double E1 = (lor.first.Energy());
    double E2 = (lor.second.Energy());
    //cout << "E1: " <<E1<<endl;
    //cout << "E2: " <<E2<<endl;
    if (isInEnergyRange(E1, Emin, Ecut) && isInEnergyRange(E2, Emin, Ecut)) {
      finalSelection.push_back(lor);
      licznik++;
    }
    licznik2++;
  }
  //cout << "procent energetyczny: " <<licznik*100/licznik2<<endl;
  return finalSelection;
}


void runTests()
{
  assert(isEqual(calculateDistance2D(1, sqrt(3), 1, -sqrt(3)), 1));
}

void clearHistograms()
{
 

  if (ROC511EcutN) {
    delete ROC511EcutN;
    ROC511EcutN = 0;
  }

  if (ROC511Ecut0) {
    delete ROC511Ecut0;
    ROC511Ecut0 = 0;
  }

  if (efficiency511EcutN) {
    delete efficiency511EcutN;
    efficiency511EcutN = 0;
  }

    if (efficiency511Ecut0) {
    delete efficiency511Ecut0;
    efficiency511Ecut0 = 0;
  }

  if (purity511EcutN) {
    delete purity511EcutN;
    purity511EcutN = 0;
  }

  if (purity511Ecut0) {
    delete purity511Ecut0;
    purity511Ecut0 = 0;
  }

    if (F1_511EcutN) {
    delete F1_511EcutN;
    F1_511EcutN = 0;
  }

  if (F1_511Ecut0) {
    delete F1_511Ecut0;
    F1_511Ecut0 = 0;
  }
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
    clearHistograms();
    //std::cout << "finalSelection  " << finalSelection << std::endl;
  }
  return 0;
}