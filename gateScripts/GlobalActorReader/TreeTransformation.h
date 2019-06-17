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
 *  @file TreeTransformation.h
 */

#ifndef TREETRANFORMATION_H
#define TREETRANFORMATION_H
#include <string>
#include "Event.h"
#include "GlobalActorReader.hh"

namespace tree_transformation
{
void addEntryToEvent(const GlobalActorReader& gar, Event* outEvent);
void clearEvent(Event* outEvent);
void transformToEventTree(const std::string& inFileName,
                          const std::string& outFileName, int maxNumEvents = -1);

}

#endif /*  !TREETRANFORMATION_H */
