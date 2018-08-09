#include <TGraph.h>
#include <TCanvas.h>
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TVector3.h"
#include <iostream>
#include <TH2.h>
#include <TRandom3.h>
#include <TMath.h>
#include <algorithm>
using namespace std;
//#include "comptonscattering.h"

//unsigned objectID_= 1;
using namespace std;
Double_t sigmaE(Double_t E, Double_t coeff=0.0444) 
{
    return coeff/TMath::Sqrt(E) * E;
}

void hsimpleReader() 
{
   // Create a histogram for the values we read.
   auto myHist = new TH1D("h1","Rozklad energii",1000,0,500);
   auto myHist2 = new TH1D("h2","Rozklad energii smeared",1000,0,500);
      TH2F *h1 = new TH2F("h1","",40,-4,4,40,-20,20);
   h1->SetFillColor(kBlue);
   // Open the file containing the tree.
   auto myFile = TFile::Open("example7.root");
   if (!myFile || myFile->IsZombie()) 
   {
      cout <<"File not opened correctly" <<endl;
      return;
   }
   // Create a TTreeReader for the tree, for instance by passing the
   // TTree's name and the TDirectory / TFile it is in.
   TTreeReader myReader("GateGlobalActorTree", myFile);
   // The branch "px" contains floats; access them as myPx.
   TTreeReaderValue<Double_t> energyLoss(myReader, "EnergyLossDuringProcess");
  TTreeReaderValue<TVector3> vect(myReader, "EmissionPointFromSource");
   Int_t n= 500;
   Double_t q[n], w[n];
   //sigma z ~ 1cm
   //x,y ~0.5cm
  double energySmeared;
  double threshold = 200;
   while (myReader.Next()) 
   {
            energySmeared = gRandom->Gaus(*energyLoss, 1000* sigmaE((*energyLoss)*1/1000));        
           // if(ProcessName =="CompotonInWater") {
            //    hist->Fill();
            //}
            //first = gRandom->Gaus(x(), 0.5));
//thurd = gRandom->Gaus(y(), 0.5));
            //h1->Fill(vect->first, vect->thurd);
              h1->Fill(gRandom->Gaus(vect->x(), 0.5), gRandom->Gaus(vect->y(), 0.5));//, gRandom->Gaus(vect->z(), 1));
            if(energySmeared < threshold) {

            }
            //cout << "energy:" <<*energyLoss <<endl;
            //cout << "sigmaE:" <<sigmaE(*energyLoss) <<endl;
            //cout << "smeared:" << energySmeared <<endl;

            myHist->Fill(*energyLoss);
            myHist2->Fill(energySmeared);
   }
             for (Int_t i=0; i < n; i++) 
      {
            q[i]=myHist->Integral(0,i);
            w[i]=i;
      }

   TGraph *gr  = new TGraph(n,q,w);
   TCanvas* c1 = new TCanvas;
   c1->SetGrid();  //wydlad pola grafiku
   // draw the graph with axis, continuous line, and put
   // a * at each point
   c1->Divide(1,3);
   c1->cd(1);
   gr->SetTitle("Zaleznosc ilosci zaakceptowanych zdarzen od progu");
   gr->GetXaxis()->SetTitle("prog");
   gr->GetYaxis()->SetTitle("zaakceptowanych zdarzen");
   gr->Draw("AC*");
   c1->cd(2);
   myHist2->GetXaxis()->SetTitle("energia");
   myHist2->GetYaxis()->SetTitle("ilosc zarezestrowanych");
   myHist2->Draw();
   c1->cd(3);
   h1->Draw();
   c1->Update();
   cout << myHist->Integral(0,1) <<endl;
   //save file
   TFile f("histogram.root","recreate");
   c1->Write();
   gr->Write();
   f.Close();
}