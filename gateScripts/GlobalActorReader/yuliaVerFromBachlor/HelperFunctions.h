/**
 *  @copyright Copyright 2018 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  @file HelperFunctions.h
 */

#ifndef HELPERFUNCTIONS_H_H
#define HELPERFUNCTIONS_H_H
#include <random>
#include <TMath.h>
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

namespace helper_functions
{

std::random_device rd;
std::mt19937 gen(rd());

Double_t sigmaE(Double_t E, Double_t coeff = 0.0444)
{
  return coeff / TMath::Sqrt(E) * E;
}

double r_norm(double mean, double sigmaE)
{
  std::normal_distribution<double> d(mean, sigmaE);
  return d(gen);
}

double smearEnergy(double energy)
{
  return r_norm(energy, 1000. * sigmaE((energy) * 1. / 1000.));
}

bool isEqual(double x, double y, double epsilon = 10e-3)
{
  return std::abs(x - y) < epsilon;
}

double calculateDistance3D(double x1, double y1, double x2, double y2, double z1, double z2)
{
  if (x1 == x2 && y1 == y2 && z1 == z2) return 0;
  // sqrt(pow(y1*(z2-z1)-z1*(y2-y1),2)+pow(x1*(z2-z1)-z1*(x2-x1),2)+pow(x1*(y2-y1)-y1*(x2-x1),2))/ sqrt(pow((y2 - y1), 2) + pow((x2 - x1), 2) + pow((z2 - z1), 2));
  return sqrt(pow(y1 * z2 - y2 * z1, 2) + pow(x1 * z2 - z1 * x2, 2) + pow(x1 * y2 - y1 * x2, 2)) / sqrt(pow((y2 - y1), 2) + pow((x2 - x1), 2) + pow((z2 - z1), 2));
}

double calculateDistance2D(double x1, double y1, double x2, double y2)
{
  return calculateDistance3D(x1, y1, x2, y2, 0, 0);
}

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

auto exactly4HitsInEvent = [](const Event& event) -> bool {
  return (event.fTracks.size() == 4) && 
  std::all_of(event.fTracks.begin(), event.fTracks.end(),
  exactlyOneHit);
};

auto exactly1HitsInEvent = [](const Event& event) -> bool {
  return (event.fTracks.size() == 1);
};

auto exactly5HitsInEvent = [](const Event& event) -> bool {
  return (event.fTracks.size() == 5);
};

auto exactlyMoreHitsInEvent = [](const Event& event) -> bool {
  return (event.fTracks.size() >= 5);
};



bool trueThreeHitsInEvent(const Event& event){
  int licz =0;
  int licz2=0;
for (const auto& track: event.fTracks)
{
  if (isEqual(track.fEmissionEnergy, 511)  ){
    licz++;
  }
  if (isEqual(track.fEmissionEnergy, 1157))
  {
    licz2++;
  }
}

  if (licz == 2 && licz2==1)
  {
    return true;
  }
  else{
    return false;
  }
};

bool trueTwoHitsInEvent(const Event& event){
  int licz =0;
for (const auto& track: event.fTracks)
{
  if (isEqual(track.fEmissionEnergy, 511)  ){
    licz++;
  }
  
}

  if (licz == 2 )
  {
    return true;
  }
  else{
    return false;
  }
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
  return step.fVolumeName == "crystal" || step.fVolumeName == "crystal1";
}

bool isScatteringInDetector2(const Hit& step)
{
  return step.fVolumeName == "crystal2";
}

bool isScatteringInDetector3(const Hit& step)
{
  return step.fVolumeName == "crystal3";
}

bool isScatteringInFantom(const Hit& step)
{
  return step.fVolumeName == "NEMA_IQ";
}

bool isScatteringInFantom2(const Hit& step)
{
  return step.fVolumeName == "cylinder_with_water";
}



auto scattering511 = [](const Hit& hit) -> bool {
  return (isEqual(hit.fEnergyBeforeProcess, 511));
};

auto scatteringprompt = [](const Hit& hit) -> bool {
  return (isEqual(hit.fEnergyBeforeProcess,1157));
};

}

#endif /*  !HELPERFUNCTIONS_H_H */
