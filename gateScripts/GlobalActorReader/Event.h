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
 *  @file Event.h
 */

#ifndef EVENT_H 
#define EVENT_H 
#include <TObject.h>
#include <string>
#include <vector>
#include <TVector3.h>


class Hit : public TObject {
public:
  Hit();
  double fEnergyBeforeProcess = -1;
  double fEnergyDeposition = -1;
  double fLocalTime = -1;
  double fGlobalTime = -1;
  TVector3 fHitPosition;
  TVector3 fVolumeCenter;
  std::string fVolumeName;
  ClassDef(Hit,1)
};

class Track : public TObject {
public:
  Track();
  int fTrackID = -1;
  double fEmissionEnergy = -1;
  std::vector<Hit> fHits;
  ClassDef(Track,1)
};


class Event: public TObject {
public:
  Event();
  int fEventID = -1;
  std::vector<Track> fTracks;
  ClassDef(Event,1)
};
#endif /*  !EVENT_H */