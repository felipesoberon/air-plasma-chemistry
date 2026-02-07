#ifndef REACTION_H
#define REACTION_H


#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
using namespace std;


class reaction
{
 private:
  
  const static int maxNoReactionSpecies = 4;

  double reactionRate;

  int reactantSpeciesList[maxNoReactionSpecies+1];
  int productSpeciesList[maxNoReactionSpecies+1];
  
  double reactionRateFunction(int j, double Tgas, double Telectron);
  double interpolateBolsigRate(int pengReaction, double TeEv);


 public:

  // Static BOLSIG+ table (shared across all instances)
  static const int bolsigMaxPoints = 50;
  static const int bolsigMaxReactions = 50;
  static bool bolsigLoaded;
  static int bolsigNPoints;
  static double bolsigTeEv[bolsigMaxPoints];
  static double bolsigRate[bolsigMaxPoints][bolsigMaxReactions];
  static int bolsigNCols;
  static int bolsigPengCol[674];  // maps Peng# -> CSV column index (-1 if absent)

  static void loadBolsigTable(const char* filename);
  
  reaction(void);

  void setReactionRate(int j, double Tgas, double Telectron);
  void setReactantAndProductSpecies(int j);
  
  int returnNumberOfReactants(void);
  int returnNumberOfProducts(void);
  
  int returnReactant(int);
  int returnProduct(int);
  
  double returnReactionRate();
  
};
#endif
