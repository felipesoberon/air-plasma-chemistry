#ifndef GLOBALMODEL_H
#define GLOBALMODEL_H

#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

#include "reaction.h"
#include "species.h"
#include "commandline.h"


class globalmodel
{
 private:

  const static int NoSpecies = 53;
  const static int NoReactions = 673;

  species SPECIES[NoSpecies+1];
  reaction REACTION[NoReactions+1];

  constexpr static int listOfSpeciesContainingH[16] = \
                       {11,12,13,14,15,16,26,27,32,44,45,46,47,48,49,50};
  
  double peakElectronTemperatureKelvin = 298;
  double electronTemperatureKelvin;
  double gasTemperatureKelvin = 298;

  unsigned long stepCount = 0;
  double dt = 50E-12;
  bool adaptiveDtEnabled = true;
  double relativeChangeMinLimit = 0.01;
  double relativeChangeMaxLimit = 0.05;
  double lastMaxRelativeChange = 0.0;
  double simulationTime = 0;
  double lastSavedSimulationTime = -1;
  double plasmaTime = 1E-9;
  double totalTime = 1E-6;

  commandline commandLine;
  
 public:

  void setParametersFromCommandLineInput(int numberOfArguments, char* valueOfArgument[]);

  void setSpeciesFormula(void);
  
  void setDefaultSpeciesDensities(void);
  void setH2ODensity(double);
  void setNO2Density(double);
  void setO3Density(double);
  double returnSaturationPressurePaAtTemperatureKelvin(double);
  double returnH2ODensityFromRelativeHumidityPercent(double);
  
  void setElectronTemperatureeV(double);
  void setElectronTemperatureKelvin(double);
  void setGasTemperatureKelvin(double);

  void setReactionRates(void);
  void setReactionReactantAndProductSpecies(void);
  void setBalanceEquations(void);
  void processBalanceEquations(void);
  double processTimeStepSpeciesDensities(void);
  void updateDtFromRelativeChangeMetric(double);
  
  void processMainLoop(void);
  
  void printSpeciesFormulaAndDensity(void);
  void printListOfSourcesAndLossesReactions(void);
  void printListOfReactions(void);
  
  void readSpeciesDensityDataFile(void);
  
  double returnElectronTemperatureKelvinAtTime(void);
  
  unsigned long returnSaveIntervalStep();
  
};
#endif
