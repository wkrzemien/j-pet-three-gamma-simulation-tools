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
#include <TFile.h>

using namespace std;



Double_t sigmaE(Double_t E, Double_t coeff = 0.0444)
{
	return coeff / TMath::Sqrt(E) * E;
}


struct AnalysisVariables
{
	int fEventID;

 //fantom

 bool fGammaScatteredBeforeInWater511ForTrack1; //gamma rozproszone przed  
 bool fGammaScatteredInDetector511ForTrack1;

 bool fGammaScatteredBeforeInWater511ForTrack2;
 bool fGammaScatteredInDetector511ForTrack2;

 bool fGammaScatteredBeforeInWaterPromptForTrack1; //gamma rozproszone przed  
 bool fGammaScatteredInDetectorPromptForTrack1;

 double fEnergyDepositionInWater511ForTrack1;
 double fEnergyDepositionInWaterPromptForTrack1;
 double fEnergyDepositionInDetetcor511ForTrack1;
 double fEnergyDepositionInDetetcorPromptForTrack1;

 double fEnergyDepositionInWater511ForTrack2;
 double fEnergyDepositionInDetetcor511ForTrack2;
 
 //511 i prompt
 bool fGammaScattered511ForTrack1; //gamma rozproszone przed  
 bool fGammaScatteredPromptForTrack1;

 bool fGammaScattered511ForTrack2;

 double fEnergyDeposition511ForTrack1;
 double fEnergyDepositionPromptForTrack1;

 double fEnergyDeposition511ForTrack2;

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

double r_norm(double mean, double sigmaE){
	normal_distribution<double> d(mean, sigmaE);
	return d(gen);
}

void Check( GlobalActorReader& gar, AnalysisVariables& av, TH1F& h511, TH1F& hPrompt, TH1F& hWaterPrompt2, TH1F& hGeoAllPrompt2, TH1F& hGeoFromWaterPrompt2, TH1F& hWater5112, TH1F& hGeoAll5112, TH1F& hGeoFromWater5112, std::string& water_volume_name)
{
	int eventID = gar.GetEventID();
	int trackID = gar.GetTrackID();
	std::string volume_name = gar.GetVolumeName();
	double energyDeposition = gar.GetEnergyLossDuringProcess();
	double emissionEnergy = gar.GetEmissionEnergyFromSource();

	if ( av.fEventID != eventID )
	{
	  //Skoro tutaj jesteśmy to znaczy, że zaczynamy nowy event
		if ( av.fEventID > -1 )
		{
		   //A jeśli tu jesteśmy to znaczy, że nie jest to pierwszy event, który jest w pliku, lecz kolejny i można próbować go zapisać
		   if ( av.fGammaScattered511ForTrack1 ) //511 track 1
		   {
			     h511.Fill(r_norm(av.fEnergyDeposition511ForTrack1, 1000 * sigmaE((av.fEnergyDeposition511ForTrack1) * 1 / 1000)));
			    if ( av.fGammaScatteredBeforeInWater511ForTrack1 ) //w wodzie
			    	hWater5112.Fill(r_norm(av.fEnergyDepositionInWater511ForTrack1, 1000 * sigmaE((av.fEnergyDepositionInWater511ForTrack1) * 1 / 1000)));
			    if ( av.fGammaScatteredInDetector511ForTrack1 )
			    {
			    	hGeoAll5112.Fill(r_norm(av.fEnergyDepositionInDetetcor511ForTrack1, 1000 * sigmaE((av.fEnergyDepositionInDetetcor511ForTrack1) * 1 / 1000)));
			    	if( av.fGammaScatteredBeforeInWater511ForTrack1 )
			    	{
			    		hGeoFromWater5112.Fill(r_norm(av.fEnergyDepositionInDetetcor511ForTrack1, 1000 * sigmaE((av.fEnergyDepositionInDetetcor511ForTrack1) * 1 / 1000)));
			    	}
			    }
			}
		   if ( av.fGammaScatteredPromptForTrack1 ) //prompt track 1
		   {
			   	hPrompt.Fill(r_norm(av.fEnergyDepositionPromptForTrack1, 1000 * sigmaE((av.fEnergyDepositionPromptForTrack1) * 1 / 1000)));
		  		if ( av.fGammaScatteredBeforeInWaterPromptForTrack1 ) //w wodzie
		  			hWaterPrompt2.Fill(r_norm(av.fEnergyDepositionInWaterPromptForTrack1, 1000 * sigmaE((av.fEnergyDepositionInWaterPromptForTrack1) * 1 / 1000)));
		  		if ( av.fGammaScatteredInDetectorPromptForTrack1 )
		  		{
		  			hGeoAllPrompt2.Fill(r_norm(av.fEnergyDepositionInDetetcorPromptForTrack1, 1000 * sigmaE((av.fEnergyDepositionInDetetcorPromptForTrack1) * 1 / 1000)));
		  			if( av.fGammaScatteredBeforeInWaterPromptForTrack1 )
		  			{
		  				hGeoFromWaterPrompt2.Fill(r_norm(av.fEnergyDepositionInDetetcorPromptForTrack1, 1000 * sigmaE((av.fEnergyDepositionInDetetcorPromptForTrack1) * 1 / 1000)));
		  			}
		  		}
		  	}

		   if ( av.fGammaScattered511ForTrack2 ) //511 track 2
		   {
		       h511.Fill(r_norm(av.fEnergyDeposition511ForTrack2, 1000 * sigmaE((av.fEnergyDeposition511ForTrack2) * 1 / 1000)));
		   		if ( av.fGammaScatteredBeforeInWater511ForTrack2 )
		   		hWater5112.Fill(r_norm(av.fEnergyDepositionInWater511ForTrack2, 1000 * sigmaE((av.fEnergyDepositionInWater511ForTrack2) * 1 / 1000)));
		   		if ( av.fGammaScatteredInDetector511ForTrack2 )
			   	{
			   		hGeoAll5112.Fill(r_norm(av.fEnergyDepositionInDetetcor511ForTrack2, 1000 * sigmaE((av.fEnergyDepositionInDetetcor511ForTrack2) * 1 / 1000)));
			   		if( av.fGammaScatteredBeforeInWater511ForTrack2 )
			   		{
			   			hGeoFromWater5112.Fill(r_norm(av.fEnergyDepositionInDetetcor511ForTrack2, 1000 * sigmaE((av.fEnergyDepositionInDetetcor511ForTrack2) * 1 / 1000)));
			   		}
			   	}
	   		}	
	   		//Zapisaliśmy dane do histogramu - o ile było co zapisywać
		}
	  //Resetujemy dane dla nowego eventu aby poprawnie rozpoznawać zdarzenia
		av.reset( eventID );
	}



	if ( emissionEnergy == 511 )
	{
	  //Czyli mamy do czynienia z 511
		if ( trackID == 1 && !av.fGammaScattered511ForTrack1 )
		{
	   	//Czyli mamy pierwsze rozporszenie 511 track ID = 1
		//Sprawdzamy czy mamay do czynienia z rozproszeniem w wodzie
			if ( volume_name == water_volume_name )
			{
	  		//Czyli mamy do czynienia z 1 rozproszeniem w wodzie
				if ( trackID == 1 && !av.fGammaScatteredBeforeInWater511ForTrack1 )
				{
	   			//Czyli mamy pierwsze rozporszenie w wodzie dla track ID = 1
	   				av.fGammaScatteredBeforeInWater511ForTrack1 = true; //w ten sposób oznaczam, że już zarajestrowałem pierwszą gamma rozproszenie w wodzie dla trackID=1
	   				av.fEnergyDepositionInWater511ForTrack1 = energyDeposition;
	   			}
	   		}
	   	}
	   	else
	   	{
	  		//Skoro nie woda to musi być to detektor
	   		if ( trackID == 1 && !av.fGammaScatteredInDetector511ForTrack1 )
	   		{
	   			av.fGammaScatteredInDetector511ForTrack1 = true;
	   			av.fEnergyDepositionInDetetcor511ForTrack1 = energyDeposition;
	   		}
	   	}
	}
	else if ( trackID == 2 && !av.fGammaScattered511ForTrack2 )
	{
	   	//Czyli mamy 2 rozporszenie 511 dla track ID = 2
	 	//Sprawdzamy czy mamay do czynienia z rozproszeniem w wodzie
	   	if ( volume_name == water_volume_name )
	   	{
	  		//Czyli mamy do czynienia z rozproszeniem w wodzie
	   		if ( trackID == 2 && !av.fGammaScatteredBeforeInWater511ForTrack2 )
	   		{
	   			//Czyli mamy pierwsze rozporszenie w wodzie dla track ID = 2
	   			av.fGammaScatteredBeforeInWater511ForTrack2 = true; //w ten sposób oznaczam, że już zarajestrowałem pierwszą gamma rozproszenie w wodzie dla trackID=1
	   			av.fEnergyDepositionInWater511ForTrack2 = energyDeposition;
	   		}
	   	}
	   	else
	   	{
	  		//Skoro nie woda to musi być to detektor
	   		if ( trackID == 2 && !av.fGammaScatteredInDetector511ForTrack2 )
	   		{
	   			av.fGammaScatteredInDetector511ForTrack2 = true;
	   			av.fEnergyDepositionInDetetcor511ForTrack2 = energyDeposition;
	   		}
	   	}

	}

	else
	{
	  	//Skoro nie 511 to musi być to 1147
	   	if ( trackID == 1 && !av.fGammaScatteredPromptForTrack1 )
	   	{
	   		if ( volume_name == water_volume_name )
	   		{
	  			//Czyli mamy do czynienia z rozproszeniem w wodzie
	   			if ( trackID == 1 && !av.fGammaScatteredBeforeInWaterPromptForTrack1 )
	   			{
	   				//Czyli mamy pierwsze rozporszenie w wodzie dla track ID = 1
	   				av.fGammaScatteredBeforeInWaterPromptForTrack1 = true; //w ten sposób oznaczam, że już zarajestrowałem pierwszą gamma rozproszenie w wodzie dla trackID=1
	   				av.fEnergyDepositionInWaterPromptForTrack1 = energyDeposition;
	   			}
	   		}
	   		else
	   		{
	  			//Skoro nie woda to musi być to detektor

	   			if ( trackID == 1 && !av.fGammaScatteredInDetectorPromptForTrack1 )
	   			{
	   				av.fGammaScatteredInDetectorPromptForTrack1 = true;
	   				av.fEnergyDepositionInDetetcorPromptForTrack1 = energyDeposition;
	   			}
	   		}
	   	}
	}
}


int main(int argc, char *argv[])
{
   	if(argc != 3)
   	{
   		std::cerr<<"Invalid number of variables."<<std::endl;
   	}
   	else
   	{
   		std::string file_name( argv[1] ); 
   		std::string water_volume_name( argv[2] );
   		double emissionEnergy();
			

   		AnalysisVariables av;
   		av.init();

	//TFile myFile("bla.root");
  //Histogramy bez rozmycia
	TH1F hWater5112( "hWater", "hWater", 200, 0, 400 ); //Wszystkie zarajestrowane
	TH1F hGeoAll5112( "hGeoAll", "hGeoAll", 200, 0, 400 );
	TH1F hGeoFromWater5112( "hGeoFromWater", "hGeoFromWater", 200, 0, 400 ); //Wszystkie w fantomie
  //Histogramy z rozmyciem
	TH1F hWaterPrompt2("hWater2", "hWater2", 200, 0, 400); //Wszystkie zarajestrowane
	TH1F hGeoAllPrompt2("hGeoAll2", "hGeoAll2", 200, 0, 400);
	TH1F hGeoFromWaterPrompt2("hGeoFromWater2", "hGeoFromWater2", 200, 0, 400); //Wszystkie w fantomie
	TH1F h511("h511", "h511", 200, 0, 400);
	TH1F hPrompt("hPrompt", "hPrompt", 200, 0, 400);
	

	try
	{
		GlobalActorReader gar;
		if(gar.LoadFile(file_name))
			while(gar.Read())
				Check(gar, av, h511, hPrompt, hWater5112, hGeoAll5112, hGeoFromWater5112, hWaterPrompt2, hGeoAllPrompt2, hGeoFromWaterPrompt2, water_volume_name);
			else
				std::cerr<<"Loading file failed."<<std::endl;
			Int_t n= 1500;
			Double_t m1, m2, m3, m4;
			Double_t TPR[n], PPV[n], FPR[n], w[n];
			for (Int_t i=0; i < n; i++) 
			{                    
            /*int pbini = hPrompt.GetXaxis()->FindBin(i);//zwraca numer binu dla ktrego energia jest i keV
            int pbinn = hPrompt.GetXaxis()->FindBin(n);
            int bini = h511.GetXaxis()->FindBin(i);
            int binn = h511.GetXaxis()->FindBin(n);
            int pbin0 = hPrompt.GetXaxis()->FindBin(a);
            int bin0 = h511.GetXaxis()->FindBin(a); */
            //PPV[i]=hPrompt.Integral(pbini,pbinn)/(hPrompt.Integral(pbini,pbinn)+h511.Integral(bini,binn));//purity
            //TPR[i]=hPrompt.Integral(pbini,pbinn)/hPrompt.Integral(pbin0,pbinn); //effi
            //FPR[i]=(h511.Integral(binn,bini))/(h511.Integral(bin0,bini));
            
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
            w[i]=i;
        }

		//TGraph	TFF(n, w, TPR); //efficiency
		//TGraph	TFF2(n, TPR, FPR); //purity
		//TGraph	TFF3(n, w, PPV); //false alarm

		//histogramy bez rozmycia
        TCanvas c("c", "c", 1000, 1000);
        hWater5112.Draw("hist l");
		//c.Print("hWater.png");
        hGeoAll5112.Draw("same");
		//c.Print("hGeoAll.png");
        hGeoFromWater5112.Draw("same");
        c.Print("bez_rozmycia_25.png");

    //histogramy z rozmyciem
        TCanvas c1("c1");
        hWaterPrompt2.Draw("hist l");
		//c1.Print("hWater2.png");
		//hGeoFromWater2.Draw("hist l");
        hGeoAllPrompt2.Draw("same");
		//c1.Print("hGeoAll2.png");
        hGeoFromWaterPrompt2.Draw("same");
        c1.Print("z_rozmyciem_25.png");

        TCanvas c2("c2");
        h511.Draw("hist l");
		//c1.Print("hWater2.png");
		//hGeoFromWater2.Draw("hist l");
        hPrompt.Draw("same");
	
        c1.Print("z_rozmyciem_205.png");
    //zapisujemy
        TFile f("histogram.root", "update");
        hGeoFromWaterPrompt2.Write("same");
        f.Close();

    //purity i efficiency
    //TCanvas c2("c2");
    //TFF.Draw();
		//c2.Print("efficiency_20.png");
		//TCanvas c4("c4");
		//TFF3.SetLineColor(kRed);
		//TFF3.Draw("same");
		//TLegend legend(0.1,0.1,0.1,0.1);
   	//legend.SetHeader("The Legend Title","C"); // option "C" allows to center the header
    //legend.AddEntry("TFF","efficiency","l");
   	//legend.AddEntry("TFF3","purity (czerwony)","l");
   	//legend.Draw();
		//2.Print("purity_and_efficiency_25.png");

    //krzywa ROC
		//TCanvas c3("c3");
		//TFF2.Draw("AC+");
		//c3.Print("ROC_25.png");


    }

       // hWater.Write();
        //myFile.Close(); 

    catch(const std::logic_error& e )
    {
    	std::cerr<<e.what()<<std::endl;
    }
    catch(...)
    {
    	std::cerr<<"Udefined exception"<<std::endl;
    }
}
return 0;
}
