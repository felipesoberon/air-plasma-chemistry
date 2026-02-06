#ifndef SPECIES_H
#define SPECIES_H

#include <iostream>
#include <string>
using namespace std;


class species
{
 private:

  const static int maxNoReactions = 624; 
  string formula;
  double density;
  constexpr static double minimumDensity = 1.0E-20;
  double loss;
  double source;

  int listOfReactionNoSources[maxNoReactions+1];
  int listOfReactionNoLosses[maxNoReactions+1]; 
  
  int listOfReactionNoSourcesMultiplier[maxNoReactions+1];
  int listOfReactionNoLossesMultiplier[maxNoReactions+1];

 public:

  species();

  void setFormula(int);
  void setDensity(double);
  void setLoss(double);
  void setSource(double);

  void incrementNumberOfReactionNoSources(void);
  void incrementNumberOfReactionNoLosses(void);
  void setReactionNoSourcesItem(int, int);
  void setReactionNoLossesItem(int, int);
  void setReactionNoSourcesMultiplierItem(int, int);
  void setReactionNoLossesMultiplierItem(int, int);

  int returnNumberOfReactionNoSources(void);
  int returnNumberOfReactionNoLosses(void);
  
  int returnReactionNoSources(int);
  int returnReactionNoLosses(int);

  int returnReactionNoSourcesMultiplier(int);
  int returnReactionNoLossesMultiplier(int);

  string returnFormula(void);
  double returnDensity(void);
  double returnLoss(void);
  double returnSource(void);

  void printListOfSourcesAndLossesReactions(void);

  void processReduceListsOfReactionNo(void);
  
};
#endif
