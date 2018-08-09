/**
 *  @copyright Copyright 2018 The J-PET Gate Authors. All rights reserved.
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  @file GlobalActorReader.hh
 */
#ifndef GLOBALACTORREADER_HH
#define GLOBALACTORREADER_HH
#include "Variable.hh"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

/** @About: This class is a reader of data from GateGlobalActor.
 * Thanks to simple implementation you just need to declar variable and add initialization of variable in InitVariable() to obtain access to data from file.
 * @Author: Mateusz Bała
 * @Email: bala.mateusz@gmail.com
*/
class GlobalActorReader
{
public:
    GlobalActorReader();
    ~GlobalActorReader();

    /** Function try open file and get acces to tree.
     * @param: file_name - ROOT file name
     * @return: false - when file does not exist, tree does not exist or tree is empty, otherwise true
    */
    bool LoadFile(const std::string& file_name);

    /** Function read data from tree. Each call of this function call next entity.
     * @return : false - when function reach maximal number of entries in tree, true - otherwise
    */
    bool Read();

    /** Resturn number number of entries in tree.
    */
    int GetEntriesNumber();

    /** Return number of initialized variables
    */
    int GetInitializedVariablesNumber();

    /** Reset counters and close file.
    */
    void Reset();

private:
    /** Try init variables. Use this function to read your variable.
     * You just need to add in this function : Notice(yourVariable.TryAttachToBranch(mTree, "yourBranchName");
    */
    void InitVariables();

    /** This function counter number of initialized variables.
     * @param: value - value from clas Variable function TryAttachToBranch
    */
    void Notice(const bool& value);

    /** Convert radians to degree
    */
    double deg(const double& angle_radians);

    //ROOT file
    TFile* pFile;
    //Tree from mFile
    TTree* pTree;
    //Number of initialized variables
    int mInitializedVariablesNumber;
    //Number of entries in tree
    int mEntriesNumber;
    //Current index of entry
    int mCurrentEntryIndex;

public:
/**GET-functions.
	This is only way how you return values;
    For more informations check variables descriptions.
*/
    std::string GetVolumeName() const;

    TVector3 GetScintilatorPosition() const;

    int GetEventID() const;

    int GetTrackID() const;

    double GetEnergyBeforeProcess() const;

    double GetEnergyAfterProcess() const;

    double GetEnergyLossDuringProcess() const;

    TVector3 GetMomentumDirectionBeforeProcess() const;

    TVector3 GetMomentumDirectionAfterProcess() const;

    TVector3 GetProcessPosition() const;

    TVector3 GetEmissionPointFromSource() const;

    TVector3 GetEmissionMomentumDirectionFromSource() const;

    double GetEmissionEnergyFromSource() const;

    std::string GetParticleName() const;

    int GetParticlePGDCoding() const;

    double GetProcessAngle() const;

    TVector3 GetPolarizationBeforeProcess() const;

    TVector3 GetPolarizationAfterProcess() const;

    std::string GetProcessName() const;

    int GetParentID() const;

    double GetInteractionTime() const;

    double GetLocalTime() const;

    double GetGlobalTime() const;

    double GetProperTime() const;

private:
/** VARIABLES
 * Declar here your variable by using template class Variable<yourType>
 * Please add description of your variable here
**/
    //Volume name (layer name) - this variable tell in which volume (layer) something happened e.g. gamma scattered in volume wit name "LayerA"
    Variable<std::string, true> VolumeName;
    //Scintilator translation vector (if detector does not rotate it is equal scintilator centrum position) - this is useful when you want to know where exacly process happened (from VolumeName you know layer and from ScintilatorPosition you know which scintilator)
    Variable<TVector3, true> ScintilatorPosition;
    //Event ID - this variable tell you which event is related to particle e.g. you have e+ and gamma - both created during event no 1 (both has the same EventID)
    Variable<int, false> EventID;
    //Track ID - each EventID has one or more tracks e.g. e+ and 2 gammas from e+ annihilation - both has EventID=1, e+ has trackID=1, gamma no 1 has trackID=2, gamma no 2 has trackID = 3
    Variable<int, false> TrackID;
    //Particle energy before process (e.g. before scattering) - the unit of energy is keV
    Variable<double, false> EnergyBeforeProcess;
    //Particle energy after process (e.g. after scattering) - the unit of energy is keV
    Variable<double, false> EnergyAfterProcess;
    //Particle energy loss during process (e.g. during scattering) - this is exactly this: EnergyLossDuringProcess = EnergyBeforeProcess - EnergyAfterProcess; the unit of energy is keV
    Variable<double, false> EnergyLossDuringProcess;
    //Particle momemntum direction before process (e.g. before scattering)
    Variable<TVector3, true> MomentumDirectionBeforeProcess;
    //Particle momemntum direction after process (e.g. after scattering)
    Variable<TVector3, true> MomentumDirectionAfterProcess;
    //Process position in lab coordinate system - by this variable you know where exacly process happened
    Variable<TVector3, true> ProcessPosition;
    //Emission point from source - by this variable you know position of source which emitted particle
    Variable<TVector3, true> EmissionPointFromSource;
    //Particle momentum direction just after emission from source
    Variable<TVector3, true> EmissionMomentumDirectionFromSource;
    //Particle energy just after emission from source - the unit of energy is keV
    Variable<double, false> EmissionEnergyFromSource;
    //Particle name e.g. gamma, e+, e-, etc. You can find particle name by checking Geant4 particle class and looking at parameter "aName" in G4ParticleDefinition constructor (it's 1st argument of constructor).  For some examples check this: https://agenda.infn.it/getFile.py/access?sessionId=26&resId=0&materialId=0&confId=12061
    Variable<std::string, true> ParticleName;
    //Particle PDG code - e.g. gamma has 0. You can find pdg code of particle by checking Geant4 particle class and looking at parameter "encoding" in G4ParticleDefinition constructor (it's 14th argument of constructor). For some examples check this: https://agenda.infn.it/getFile.py/access?sessionId=26&resId=0&materialId=0&confId=12061
    Variable<int, false> ParticlePGDCoding;
    //Angle betwean particle momentum before and after process
    Variable<double, false> ProcessAngle;
    //Particle polarization before process happened
    Variable<TVector3, true> PolarizationBeforeProcess;
    //Particle polarization after process happened
    Variable<TVector3, true> PolarizationAfterProcess;
    //Process name e.g. "Compton", "Ionization", etc. - if you use Gate physics lists, "compt", etc. - if you use Geant4 physics lists. For GATE check this page: http://wiki.opengatecollaboration.org/index.php/Users_Guide:Setting_up_the_physics ; for Geant4 find specific class and/or check Geant4 documentation
    Variable<std::string, true> ProcessName;
    //Parent ID - parent ID tell which track created current track e.g e+ has EventID = 1 and TrackID=1 and e+ decay to 2 gammas, gamma has EventID=1, TrackID=2, ParentID=1 etc.
    Variable<int, false> ParentID;
    //Interaction time - difference between local time of pre and post step of track
    Variable<double, false> InteractionTime;
    //Local time - time since the track was created
    Variable<double, false> LocalTime;
    //Global time -  time since the event was created
    Variable<double, false> GlobalTime;
    //Proper time -  time in its rest frame since the track was created
    Variable<double, false> ProperTime;

public:
/** EXTRA FUNCTIONS
 * If you need variable which is based on other variables lisetd in InitVariables().
 * Please add here your GetSthValue() e.g. GetComptonThetaValue();
 * Remember add decription of file
*/
    /** Return Compton scaterring angle (theta angle).
     * This is angle betwean prime and scattered gamma momentum direction.
     * Angle unit: degree
    */
    double GetComptonThetaValue();

    /** Return Compton phi angle.
     * This is angle betwean prime gamma polarization vector and scattered gamma momentum vector throw on
     * the plane perpendicular to prime gamma momentum direction.
     * Angle unit: degree
    */
   double GetComptonPhiValue();

};

#endif // GLOBALACTORREADER_HH
