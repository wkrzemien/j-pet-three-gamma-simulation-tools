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

bool isEqual(double x, double y, double epsilon = 10e-9)
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

}

#endif /*  !HELPERFUNCTIONS_H_H */