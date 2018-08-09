#include <iostream>
#include "GlobalActorReader.hh"
#include <stdexcept>
#include "TCanvas.h"
#include "TH1F.h"
#include<string>
#include <TRandom3.h>
#include <TMath.h>
#include <algorithm>
#include <TGraph.h>
#include <TLegend.h>

using namespace std;


/***
CO ROBI TEN PRZYKŁADOWY KOD:
Rozróżni czy rozproszenie nastąpiło w objętości z wodą czy w geometrii.
Na podstawie tych informacji wykręsla histogramy:
1. depozycja energii w objętości z wodą - patrzymy tylko na pierwsze rozproszenie w wodzie - jak chcecie mieć mięcej rozproszeń zarejestrowanych to trzeba odpowiedni kod napisać
2. depozycja ogólnie w geometrii detektora podczas pierwszego rozproszenie Comptonowskiego w tym detektrorze - czyli może ale nie musi być tu gamma rozproszona w wodzie ( jeśli byście chcieli badać kolejne rozproszenia to konieczne jest sprawdzanie wektora położenia scyntylatora - bo mogą i zachodzą wielokrotne rozproszenia w scyntylatorach nawet dla tej geometrii)
3. depozycja energii w geometrii detektora tylko przez gammy które rozproszyły się wcześniej w objętości z wodą

Zauważcie, że tutaj nie interesuje nas w której konkretnej warstwie nastąpiło rozproszenie - dlatego nie trzeba sprawdzać nazw warstw detektora. Jeśli chcielibyście robić studium dla konkretnych warstw detektora to już trzeba dodać odpowiednie sprawdzania.

Ten kod ma taką a nie inną formę ponieważ wiemy na podstawie konfiguracji actora, że:
1. rejestrujemy tylko gammy
2. rejestrujemy tylko Comptona
***/
Double_t sigmaE(Double_t E, Double_t coeff = 0.0444)
{
  return coeff / TMath::Sqrt(E) * E;
}


struct AnalysisVariables {
  int fEventID;

  bool fGammaScattered511ForTrack1; //gamma rozproszone przed
  bool fGammaScatteredPromptForTrack1;

  bool fGammaScattered511ForTrack2;
//bool fGammaScatteredInDetectorForTrack2;

  double fEnergyDeposition511ForTrack1;
  double fEnergyDepositionPromptForTrack1;

  double fEnergyDeposition511ForTrack2;
//double fEnergyDepositionPromptForTrack2;

  void init();
  void reset( int eventID );
};

void AnalysisVariables::init()
{
  reset( -1 );
}

void AnalysisVariables::reset( int eventID )
{
  fEventID = eventID;
  fGammaScattered511ForTrack1 = fGammaScatteredPromptForTrack1 = fGammaScattered511ForTrack2 = false;
  fEnergyDeposition511ForTrack1 = fEnergyDepositionPromptForTrack1 = fEnergyDeposition511ForTrack2 = 0;
}

random_device rd;
mt19937 gen(rd());

double r_norm(double mean, double sigmaE)
{
  normal_distribution<double> d(mean, sigmaE);
  return d(gen);
}

void Check( GlobalActorReader& gar, AnalysisVariables& av,  TH1F& h511, TH1F& hPrompt, TH1F& h5112, TH1F& hPrompt2)
{
  int eventID = gar.GetEventID();
  int trackID = gar.GetTrackID();
  double energyDeposition = gar.GetEnergyLossDuringProcess();
  double emissionEnergy = gar.GetEmissionEnergyFromSource();

  if ( av.fEventID != eventID ) {
    //Skoro tutaj jesteśmy to znaczy, że zaczynamy nowy event
    if ( av.fEventID > -1 ) {
      //A jeśli tu jesteśmy to znaczy, że nie jest to pierwszy event, który jest w pliku, lecz kolejny i można próbować go zapisać

      if ( av.fGammaScattered511ForTrack1 )
        h5112.Fill(av.fEnergyDeposition511ForTrack1);
      h511.Fill(r_norm(av.fEnergyDeposition511ForTrack1, 1000 * sigmaE((av.fEnergyDeposition511ForTrack1) * 1 / 1000)));
      if ( av.fGammaScatteredPromptForTrack1 ) {
        hPrompt2.Fill(av.fEnergyDepositionPromptForTrack1);
        hPrompt.Fill(r_norm(av.fEnergyDepositionPromptForTrack1, 1000 * sigmaE((av.fEnergyDepositionPromptForTrack1) * 1 / 1000)));
      }

      if ( av.fGammaScattered511ForTrack2 ) {
        h5112.Fill( av.fEnergyDeposition511ForTrack2 );
        h511.Fill(r_norm(av.fEnergyDeposition511ForTrack2, 1000 * sigmaE((av.fEnergyDeposition511ForTrack2) * 1 / 1000)));
      }
      if ( av.fGammaScattered511ForTrack2 ) {
        h5112.Fill( av.fEnergyDeposition511ForTrack2 );
        h511.Fill(r_norm(av.fEnergyDeposition511ForTrack2, 1000 * sigmaE((av.fEnergyDeposition511ForTrack2) * 1 / 1000)));
      }


      //Zapisaliśmy dane do histogramu - o ile było co zapisywać
    }
    //Resetujemy dane dla nowego eventu aby poprawnie rozpoznawać zdarzenia
    av.reset( eventID );
  }

//Sprawdzamy z jakim rozproszeniem mamy doczynienie
  if ( emissionEnergy == 511 ) {
    //Czyli mamy do czynienia z 511

    if ( trackID == 1 && !av.fGammaScattered511ForTrack1 ) {
      //Czyli mamy pierwsze rozporszenie w wodzie dla track ID = 1
      av.fGammaScattered511ForTrack1 = true; //w ten sposób oznaczam, że już zarajestrowałem pierwszą gamma rozproszenie w wodzie dla trackID=1
      av.fEnergyDeposition511ForTrack1 = energyDeposition;
    } else if ( trackID == 2 && !av.fGammaScattered511ForTrack2 ) {
      //Czyli mamy pierwsze rozporszenie w wodzie dla track ID = 2
      av.fGammaScattered511ForTrack2 = true; //w ten sposób oznaczam, że już zarajestrowałem pierwszą gamma rozproszenie w wodzie dla trackID=1
      av.fEnergyDeposition511ForTrack2 = energyDeposition;
    }
  } else {
    //Skoro nie 511 to musi być to 1147

    if ( trackID == 1 && !av.fGammaScatteredPromptForTrack1 ) {
      av.fGammaScatteredPromptForTrack1 = true;
      av.fEnergyDepositionPromptForTrack1 = energyDeposition;
    }
  }

}


int main(int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Invalid number of variables." << std::endl;
  } else {
    std::string file_name( argv[1] );
    std::string emissionEnergy( argv[2] );

    AnalysisVariables av;
    av.init();

    //Histogramy bez rozmycia
    TH1F h511( "h511", "h511", 200, 0, 1400 ); //Wszystkie zarajestrowane
    TH1F hPrompt( "hPrompt", "hPrompt", 200, 0, 1400 );
    //TH1F hGeoFromWater( "hGeoFromWater", "hGeoFromWater", 1000, 0, 400 ); //Wszystkie w fantomie
    //Histogramy z rozmyciem
    TH1F h5112("h5112", "h5112", 200, 0, 400); //Wszystkie zarajestrowane
    TH1F hPrompt2("hPrompt2", "hPrompt2", 200, 0, 1400);
    //TH1F hGeoFromWater2("hGeoFromWater2", "hGeoFromWater2", 1000, 0, 400); //Wszystkie w fantomie


    try {
      GlobalActorReader gar;
      if (gar.LoadFile(file_name))
        while (gar.Read())
          Check(gar, av, h511, hPrompt, h5112, hPrompt2);
      else
        std::cerr << "Loading file failed." << std::endl;

      TFile file("histos.root", "recreate");
      hPrompt.Write();
      h511.Write();
      file.Close();

      Int_t n = 1500;
      Double_t m1, m2, m3, m4;
      Double_t TPR[n], PPV[n], FPR[n], w[n];
      for (Int_t i = 0; i < n; i++) {
        m4 = hPrompt.Integral(i, n); //hSignal
        TPR[i] = hPrompt.Integral(i, n) / hPrompt.Integral(0, n); //effi
        PPV[i] = (m4) / (hPrompt.Integral(i, n) + h511.Integral(i, n)); //purity
        FPR[i] = (h511.Integral(n, i)) / (h511.Integral(0, i));
        // purity = 511 not scattered in phantom/all registered
        // hSignal = hGeoAll2 -hWater2;  TP
        // purity = hSignal(i,n)/hGeoAll2(i,n)
        // hEff = hSignal(i,n)/hSignal(0,n)

        // One plot: hWater from smeared Energy

        // 511 keV -not scattered in phantonm
        // ~1100 keV not scattered in phantom
        // background:
        // 511 keV scattered in phantom
        // ~1100 keV scattered in phanton
        w[i] = i;
      }

      TGraph	TFF(n, w, TPR); //efficiency
      TGraph	TFF2(n, TPR, FPR); //purity
      TGraph	TFF3(n, w, PPV); //false alarm

      //histogramy bez rozmycia
      TCanvas c("c", "c", 1000, 1000);
      h511.Draw("hist l");
      //c.Print("hWater.png");
      hPrompt.Draw("same");
      //c.Print("hGeoAll.png");
      //hGeoFromWater.Draw("same");
      c.Print("bez_rozmycia_30.png");

      //histogramy z rozmyciem
      TCanvas c1("c1");
      h5112.Draw("hist l");
      //c1.Print("hWater2.png");
      //hGeoFromWater2.Draw("hist l");
      hPrompt2.Draw("same");
      //c1.Print("hGeoAll2.png");
      //hGeoFromWater2.Draw("same");
      c1.Print("z_rozmyciem_30.png");

      //purity i efficiency
      TCanvas c2("c2", "c2", 1000, 1000);
      TFF.Draw();
      //c2.Print("efficiency_20.png");
      //TCanvas c4("c4");
      //TFF.SetLineColor(kRed);
      //TFF.Draw("same");
      TLegend legend(0.1, 0.1, 0.1, 0.1);
      //legend.SetHeader("The Legend Title","C"); // option "C" allows to center the header
      legend.AddEntry("TFF", "efficiency", "l");
      legend.AddEntry("TFF3", "purity (czerwony)", "l");
      legend.Draw();
      c2.Print("purity_and_efficiency_30.png");

      //krzywa ROC
      TCanvas c3("c3");
      TFF2.Draw("AC+");
      c3.Print("ROC.png");


    }

    // hWater.Write();
    //myFile.Close();
    catch (const std::logic_error& e ) {
      std::cerr << e.what() << std::endl;
    } catch (...) {
      std::cerr << "Udefined exception" << std::endl;
    }
  }
  return 0;
}
