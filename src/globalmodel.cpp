#include "globalmodel.h"


void globalmodel::setParametersFromCommandLineInput(int numberOfArguments, char* valueOfArgument[])
{
  cout << "\nCOMMAND_LINE_INPUT_PARAMETERS\n\n";

  commandLine.setArgumentList(numberOfArguments, valueOfArgument);
  commandLine.printArgumentList();

  commandLine.setFlagName("-Te", "Electron temperature in eV");
  commandLine.setFlagName("-[H2O]", "Density of water in m-3");
  commandLine.setFlagName("-RH", "Relative humidity in %");
  commandLine.setFlagName("-totaltime", "Total simulation time is s");
  commandLine.setFlagName("-plasmatime", "Plasma pulse time in s");
  commandLine.setFlagName("-dt","Simulation time step in s");
  commandLine.setFlagName("-metricmin","Adaptive dt lower relative-change limit");
  commandLine.setFlagName("-metricmax","Adaptive dt upper relative-change limit");
  commandLine.printFlagNameList();

  commandLine.setFlagValues();
  commandLine.printFlagValues();

  if (commandLine.flagValueIsNumber(0))
    peakElectronTemperatureKelvin = 11605. * commandLine.returnFloatFlagValue(0);
  
  if (commandLine.flagValueIsNumber(1))
    setH2ODensity( commandLine.returnFloatFlagValue(1) );

  if (commandLine.flagValueIsNumber(2))
    setH2ODensity( returnH2ODensityFromRelativeHumidityPercent(commandLine.returnFloatFlagValue(2)) );
  
  if (commandLine.flagValueIsNumber(3))
    totalTime =  commandLine.returnFloatFlagValue(3);
  
  if (commandLine.flagValueIsNumber(4))
    plasmaTime =  commandLine.returnFloatFlagValue(4);

  if (commandLine.flagValueIsNumber(5))
    {
      dt =  commandLine.returnFloatFlagValue(5);
      adaptiveDtEnabled = false;
    }

  if (commandLine.flagValueIsNumber(6))
    relativeChangeMinLimit = commandLine.returnFloatFlagValue(6);

  if (commandLine.flagValueIsNumber(7))
    relativeChangeMaxLimit = commandLine.returnFloatFlagValue(7);

  if (relativeChangeMinLimit <= 0.0 || relativeChangeMaxLimit <= 0.0 || relativeChangeMinLimit >= relativeChangeMaxLimit)
    {
      cout << "WARNING: invalid adaptive metric limits; resetting to defaults (0.01, 0.05)." << endl;
      relativeChangeMinLimit = 0.01;
      relativeChangeMaxLimit = 0.05;
    }

  cout << endl;
 }






void globalmodel::setSpeciesFormula(void)
{
  int i;
  for (i=0; i<=NoSpecies; i++)
    {
      SPECIES[i].setFormula(i);
    }
}









void globalmodel::setDefaultSpeciesDensities(void)
{
  int i;
  for (i=0; i<=NoSpecies; i++)
    {
      SPECIES[i].setDensity(0.0);
    } 
  SPECIES[0].setDensity(  2.40E25 ); //M (e.g., O2 + N2)
  SPECIES[17].setDensity( 1.00E3  ); //e
  SPECIES[36].setDensity( 0.00E0  ); //O3
  SPECIES[39].setDensity( 0.00E0  ); //NO2
  SPECIES[51].setDensity( 1.92E25 ); //N2
  SPECIES[52].setDensity( 4.80E24 ); //O2
  SPECIES[53].setDensity( 1.20E24 ); //H2O  (5% of 2.4E25)
}





void globalmodel::setH2ODensity(double densityValue)
{ SPECIES[53].setDensity( densityValue ); }


void globalmodel::setNO2Density(double densityValue)
{ SPECIES[39].setDensity( densityValue ); }


void globalmodel::setO3Density(double densityValue)
{ SPECIES[36].setDensity( densityValue ); }

double globalmodel::returnSaturationPressurePaAtTemperatureKelvin(double temperatureValueKelvin)
{
  // Tetens approximation (water over liquid), valid near room conditions.
  double temperatureValueCelsius = temperatureValueKelvin - 273.15;
  return 610.94 * exp( (17.625 * temperatureValueCelsius) / (temperatureValueCelsius + 243.04) );
}

double globalmodel::returnH2ODensityFromRelativeHumidityPercent(double relativeHumidityPercent)
{
  constexpr double BoltzmannConstant = 1.380649E-23;
  double relativeHumidityFraction = relativeHumidityPercent / 100.0;
  if (relativeHumidityFraction < 0.0) relativeHumidityFraction = 0.0;
  if (relativeHumidityFraction > 1.0) relativeHumidityFraction = 1.0;

  double saturationPressurePa = returnSaturationPressurePaAtTemperatureKelvin(gasTemperatureKelvin);
  double h2oPartialPressurePa = relativeHumidityFraction * saturationPressurePa;
  return h2oPartialPressurePa / (BoltzmannConstant * gasTemperatureKelvin);
}


void globalmodel::setElectronTemperatureeV(double temperatureValueeV)
{ electronTemperatureKelvin = temperatureValueeV * 11605.; }

void globalmodel::setElectronTemperatureKelvin(double temperatureValue)
{ electronTemperatureKelvin = temperatureValue; }

void globalmodel::setGasTemperatureKelvin(double temperatureValue)
{ gasTemperatureKelvin = temperatureValue; }



void globalmodel::setReactionRates(void)
{
  if (!reaction::bolsigLoaded)
    reaction::loadBolsigTable("bolsigplus072024-win/02_output/rates_Peng.csv");

  int j;
  for (j=1;  j<=NoReactions; j++)
    {
      if (simulationTime > 1E-4 && j >= 625)
	REACTION[j].setReactionRate(0, gasTemperatureKelvin, electronTemperatureKelvin);
      else
	REACTION[j].setReactionRate(j, gasTemperatureKelvin, electronTemperatureKelvin);
    }
}



void globalmodel::setReactionReactantAndProductSpecies()
{
  int j;
  for (j=1;  j<=NoReactions; j++)
    {
      REACTION[j].setReactantAndProductSpecies(j);
    }
}








void globalmodel::setBalanceEquations(void)
{
  int i, j, k;
  int NoReactants, NoProducts;
  int itemIndex;
  int repeatsource, repeatloss; //number of times species in reaction 	
  
  for (i=1; i<=NoSpecies; i++)
    {		
      for (j=1; j<=NoReactions; j++)
	{
	  
	  /* LOSSES */
	  NoReactants = REACTION[j].returnNumberOfReactants();
	  repeatloss = 0; 
	  for (k=1; k<=NoReactants; k++)
	    {
	      if (REACTION[j].returnReactant(k) == i) { repeatloss++; } 
	    } 
	  if (repeatloss>0) 
            {
	      SPECIES[i].incrementNumberOfReactionNoLosses();
	      itemIndex = SPECIES[i].returnNumberOfReactionNoLosses();
	      SPECIES[i].setReactionNoLossesItem( itemIndex, j);
	      SPECIES[i].setReactionNoLossesMultiplierItem( itemIndex, repeatloss);
	    }
	  
	  /* SOURCES */  
	  NoProducts = REACTION[j].returnNumberOfProducts(); 
	  repeatsource = 0; 
	  for (k=1; k<=NoProducts; k++) 
	    {
	      if (REACTION[j].returnProduct(k) == i) { repeatsource++; }
	    }
	  if (repeatsource>0)
	    {
	      SPECIES[i].incrementNumberOfReactionNoSources();
	      itemIndex = SPECIES[i].returnNumberOfReactionNoSources();
	      SPECIES[i].setReactionNoSourcesItem( itemIndex, j);
	      SPECIES[i].setReactionNoSourcesMultiplierItem( itemIndex, repeatsource);
	    }
	  
	}//for j (run over all reactions)

      SPECIES[i].processReduceListsOfReactionNo();
      
    }//for i (run over species)	
}





void globalmodel::processBalanceEquations(void)
{
  int i; //all species index 
  int j; //reaction number
  int g; //index of species containing hydrogen
  int h; //reaction list index
  int ri;//reactants index
  
  int multiplier;
  int NoReactions;
  int NoReactants;
  int reactantIndex;	
  
  double auxDensity;
  double auxLoss;
  double auxSource;
  double reactionRate; 
  
  bool flagSpeciesContainingH;  
  
  /*Run over all species, except N2, O2, H2O*/
  for (i=1; i<=50; i++)
    {
      SPECIES[i].setLoss(0.0);
      SPECIES[i].setSource(0.0);
      
      /*[H20] == 0 ?*/
      if ( SPECIES[53].returnDensity() == 0 ) 
	{
	  flagSpeciesContainingH = false;
	  for (g=0; g<16; g++)
	    {
	      if ( i == listOfSpeciesContainingH[g] )
		{
		  flagSpeciesContainingH = true;
		  break;
		}
	    }
	  if ( flagSpeciesContainingH ) continue;
	}
      
      /* LOSSES */
      if (SPECIES[i].returnDensity() > 0)
	{
	  NoReactions = SPECIES[i].returnNumberOfReactionNoLosses();
	  for (h=1; h<=NoReactions; h++)
	    {
	      j = SPECIES[i].returnReactionNoLosses(h); 
	      NoReactants = REACTION[j].returnNumberOfReactants();				
	      auxDensity = 1.0;
	      for (ri=1; ri<=NoReactants; ri++) 
		{ 
		  reactantIndex = REACTION[j].returnReactant(ri);
		  auxDensity = auxDensity * SPECIES[reactantIndex].returnDensity(); 
		}
	      multiplier = SPECIES[i].returnReactionNoLossesMultiplier(h);
	      reactionRate = REACTION[j].returnReactionRate();
	      auxLoss = SPECIES[i].returnLoss();
	      
	      auxLoss = auxLoss + reactionRate * auxDensity * multiplier;
	      
	      SPECIES[i].setLoss(auxLoss);
	    }	  
	}
      
      
      /* SOURCES */  
      NoReactions = SPECIES[i].returnNumberOfReactionNoSources();
      for (h=1; h<=NoReactions; h++)
	{	  
	  j = SPECIES[i].returnReactionNoSources(h); 
	  NoReactants = REACTION[j].returnNumberOfReactants();
	  auxDensity = 1.0;
	  for (ri=1; ri<=NoReactants; ri++) 
	    { 
	      reactantIndex = REACTION[j].returnReactant(ri);
	      auxDensity = auxDensity * SPECIES[reactantIndex].returnDensity(); 
	    }
	  multiplier = SPECIES[i].returnReactionNoSourcesMultiplier(h);
	  reactionRate = REACTION[j].returnReactionRate();
	  auxSource = SPECIES[i].returnSource();
	  
	  auxSource = auxSource + reactionRate * auxDensity * multiplier;
	  
	  SPECIES[i].setSource(auxSource);			
	}
    }//for i (run over species)	
}







double globalmodel::processTimeStepSpeciesDensities(void)
{
  int i;
  double SourcesMinusLosses;
  double currentDensity;
  double deltaDensity;
  double relativeChange;
  double maxRelativeChange = 0.0;
  double auxDensity;
  
  /*Run over all species, except N2, O2, H2O*/
  for (i=1; i<=50; i++)
    {
      if ( SPECIES[i].returnLoss() >0 || SPECIES[i].returnSource() >0) 
	{ 
	  currentDensity = SPECIES[i].returnDensity();
	  SourcesMinusLosses = SPECIES[i].returnSource() - SPECIES[i].returnLoss();
	  deltaDensity = SourcesMinusLosses * dt;

	  if (currentDensity > 0.0)
	    {
	      relativeChange = fabs(deltaDensity/currentDensity);
	      if (relativeChange > maxRelativeChange) maxRelativeChange = relativeChange;
	    }
	  
	  if (SourcesMinusLosses>0 && deltaDensity>currentDensity && currentDensity>1E5)
	    {
	      cout << "WARNING: species [";
	      cout << SPECIES[i].returnFormula();
	      cout << "] > x2 density at time ";
	      cout << simulationTime;
	      cout << endl;
	    }
	  
	  auxDensity = currentDensity + deltaDensity;
	  
	  if (SourcesMinusLosses != 0.0) 
	    { SPECIES[i].setDensity( auxDensity ); }
	}
    }//for i (run over species)

  lastMaxRelativeChange = maxRelativeChange;
  return maxRelativeChange;
}







void globalmodel::updateDtFromRelativeChangeMetric(double maxRelativeChange)
{
  if (!adaptiveDtEnabled) return;
  if (maxRelativeChange > relativeChangeMaxLimit) dt = dt * 0.5;
  else if (maxRelativeChange < relativeChangeMinLimit) dt = dt * 2.0;
}




void globalmodel::processMainLoop(void)
{
  int i;
  double maxRelativeChange;
  ofstream dumpfile;             
  dumpfile.open("output.csv", ios_base::app);	
  
  if (simulationTime == 0)
    {		
      dumpfile << "#"; 
      for (i=1; i<=NoSpecies; i++) 
	{
	  dumpfile << SPECIES[i].returnFormula();
	  dumpfile << ",";
	}
      dumpfile << "Time(s),StepNo";
      dumpfile << endl;
    }
  
  /* Plasma Pulse*/
  cout << "\nPLASMA PULSE: (duration = " << plasmaTime << ")\n";
  while (true)
    {
      if (simulationTime >= totalTime || simulationTime >= 10*plasmaTime) break;
      
      electronTemperatureKelvin = returnElectronTemperatureKelvinAtTime();
      setReactionRates();		
      
      processBalanceEquations();
      maxRelativeChange = processTimeStepSpeciesDensities();
      updateDtFromRelativeChangeMetric(maxRelativeChange);
      
      if ( stepCount % returnSaveIntervalStep() == 0 && simulationTime > lastSavedSimulationTime ) 
	{
	  for (i=1; i<=NoSpecies; i++) dumpfile << SPECIES[i].returnDensity() << ",";
	  dumpfile << simulationTime << "," << stepCount << endl;
	  cout << simulationTime<< "\t";
	  cout << stepCount << "\t";
	  cout << electronTemperatureKelvin;
	  cout << "\t";
	  cout << dt;
	  cout << endl;
	}
      simulationTime = simulationTime + dt;
      stepCount++;
    }//plasma pulse
  
  electronTemperatureKelvin = returnElectronTemperatureKelvinAtTime();
  setReactionRates();		
  
  /* the MAIN LOOP: Afterglow*/
  cout << "\nAFTERGLOW:\n" << endl;
  while (true)
    {
      /*exit the simulation if time */
      if (simulationTime >= totalTime) break;
      
      processBalanceEquations();
      maxRelativeChange = processTimeStepSpeciesDensities();
      updateDtFromRelativeChangeMetric(maxRelativeChange);
      
      if ( stepCount % returnSaveIntervalStep() == 0 && simulationTime > lastSavedSimulationTime ) 
	{
	  for (i=1; i<=NoSpecies; i++) dumpfile << SPECIES[i].returnDensity() << ",";
	  dumpfile << simulationTime << "," << stepCount << endl;
	  cout << simulationTime<< "\t";
	  cout << stepCount << "\t";
	  cout << electronTemperatureKelvin;
	  cout << "\t";
	  cout << dt;
	  cout << endl;
	}		
      simulationTime = simulationTime + dt;
      stepCount++;
    }//afterglow
  
  dumpfile.close();	
}







void globalmodel::printSpeciesFormulaAndDensity(void)
{
  int i;
  cout << "\nSPECIES\n\n";
  for (i=0; i<=NoSpecies; i++)
    {
      cout << i << ",";
      cout << SPECIES[i].returnFormula();
      cout << ",";
      cout << SPECIES[i].returnDensity();
      cout << endl;
    }
}





void globalmodel::printListOfSourcesAndLossesReactions(void)
{
  int i;
  cout << "\nSPECIES_AND_THEIR_LIST_OF_SOURCE_AND_LOSSES_REACTIONS\n\n";
  for (i=0; i<=NoSpecies; i++)
    {
      SPECIES[i].printListOfSourcesAndLossesReactions();
    }
}






void globalmodel::printListOfReactions(void)
{
  int i, j;
  int NoReactants, NoProducts;
  int ReactantIndex, ProductIndex;
  double reactionRateValue;
	
  cout << "\nREACTION_LIST\n\n";
  cout << "No.,Rate,r1,r2,r3,r4,--->,p1,p2,p3,p4\n";
	
  for (i=1; i<=NoReactions; i++)
    {
      reactionRateValue = REACTION[i].returnReactionRate();
		
      NoReactants= REACTION[i].returnNumberOfReactants();
      NoProducts = REACTION[i].returnNumberOfProducts();
		
      cout << i << "," << reactionRateValue << ",";
		
      /*Reactants*/
      for (j=1; j<=4; j++) 
	{
	  if (j <= NoReactants) 
	    {
	      ReactantIndex = REACTION[i].returnReactant(j);
	      cout << SPECIES[ReactantIndex].returnFormula(); 
	      cout << ",";
	    }
	  else 
	    {
	      cout << ",";
	    }
	}
		
      cout << "--->,";
		
      /*Products*/
      for (j=1; j<=4; j++) 
	{
	  if (j <= NoProducts) 
	    {
	      ProductIndex = REACTION[i].returnProduct(j);
	      cout << SPECIES[ProductIndex].returnFormula();
	      cout << ",";
	    }
	  else 
	    {
	      cout << ",";
	    }
	}
		
      cout << endl;
    }
}





void globalmodel::readSpeciesDensityDataFile(void)
{
  int numberOfLines;
  string line;  
  string dataLine;  
	
  fstream file; 
  file.open("output.csv");
  if(!file.is_open())
    {
      cout << "WARNING: Problem opening <output.csv> file.\n"; 
      cout << "NOTICE: Program will default to initial species density values.\n";
    }
  else
    {
      numberOfLines = 0;
      while(!file.eof())
	{
	  getline(file, line);
	  if (line.length() > 50 ) 
	    {
	      dataLine = line;
	    }
	  numberOfLines++;	  
	}		
      file.close();
		
      stringstream dataLineString(dataLine);
      double val[NoSpecies+3]; //species:1-53, simulationTime:54, StepCount:55 
      int i=1;
      while ( dataLineString >> val[i] && i<=(NoSpecies+2) )
	{
	  if (dataLineString.peek() == ',') dataLineString.ignore();
	  i++;
	}
      cout << "\nNOTICE: <output.csv> file read... " << i-1 << endl;
		
      double densityOfSpeciesM = 0;
      for (i=1; i<=NoSpecies; i++) 
	{
	  densityOfSpeciesM = densityOfSpeciesM + val[i];
	  SPECIES[i].setDensity( val[i] );
	}
      SPECIES[0].setDensity( densityOfSpeciesM );		
      simulationTime = val[54];
      lastSavedSimulationTime = simulationTime;
      stepCount = (unsigned long)val[55];
    }
}





double globalmodel::returnElectronTemperatureKelvinAtTime(void)
{
  double temperatureValue = gasTemperatureKelvin;
  if ( simulationTime < 10*plasmaTime )
    {
      temperatureValue = gasTemperatureKelvin + ( peakElectronTemperatureKelvin \
						  - gasTemperatureKelvin ) * exp( -0.5 \
										  * pow( (simulationTime-5*plasmaTime)/plasmaTime,2) );
    }
  return temperatureValue;
}		





unsigned long globalmodel::returnSaveIntervalStep()
{
  unsigned long result =1;
  double ratioFactor; //ratio between dt and 50ps
  ratioFactor = dt/50E-12;
	
  /* step intervals for simualtions up to 1000 seconds */
  if (simulationTime <= 1E3 ) result = 1000000000000;
  if (simulationTime <= 1E2 ) result = 100000000000;
  if (simulationTime <= 1E1 ) result = 10000000000;
  if (simulationTime <= 1E-0) result = 1000000000;
  if (simulationTime <= 1E-1) result = 100000000;
  if (simulationTime <= 1E-2) result = 10000000;
  if (simulationTime <= 1E-3) result = 1000000;
  if (simulationTime <= 1E-4) result = 100000;
  if (simulationTime <= 1E-5) result = 10000;
  if (simulationTime <= 1E-6) result = 1000;
  if (simulationTime <= 1E-7) result = 100;
  if (simulationTime <= 1E-8) result = 10;
  if (simulationTime <= 1E-9) result = 1;
	
  ratioFactor = ( (double)result ) / ratioFactor;
  result  = (unsigned long)ratioFactor; 
  if (result < 1) 
    {
      result = 1;
    }
  return result;	
}

