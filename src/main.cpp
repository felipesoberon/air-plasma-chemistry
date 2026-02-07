#include "globalmodel.h"

void printProgramDescription(void);

int main (int argc, char* argv[])
{
  printProgramDescription();
  
  globalmodel MicrodischargeModel;
    
  MicrodischargeModel.setSpeciesFormula();
  MicrodischargeModel.setDefaultSpeciesDensities();
  MicrodischargeModel.setParametersFromCommandLineInput(argc, argv);

  MicrodischargeModel.readSpeciesDensityDataFile(); //override defaults 
  
  MicrodischargeModel.setReactionRates();
  MicrodischargeModel.setReactionReactantAndProductSpecies();
  MicrodischargeModel.setBalanceEquations();
  
    
  //MicrodischargeModel.printSpeciesFormulaAndDensity();
  //MicrodischargeModel.printListOfSourcesAndLossesReactions();
  //MicrodischargeModel.printListOfReactions();
  

  MicrodischargeModel.processMainLoop();
  
  
}






void printProgramDescription(void)
{
  cout << "**********************************************" << endl;
  cout << "*                                            *" << endl;
  cout << "* Felipe Soberon (felipe.soberon@gmail.com)  *" << endl;
  cout << "* 2026                                       *" << endl;
  cout << "*                                            *" << endl;
  cout << "* Global_Model_2.1 of atmospheric pressure   *" << endl;
  cout << "* plasma discharge in (humid) air;           *" << endl;
  cout << "* using Sakiyama et al., 2012 (R1-R624)      *" << endl;
  cout << "* and Peng et al., 2022 (R625-R673)          *" << endl;
  cout << "* reaction data.                             *" << endl;
  cout << "*                                            *" << endl;
  cout << "**********************************************" << endl;
  
  cout << "Usage:" << endl;
  cout << "$ airGM2.1 <-flag> <flag value>";
  cout << endl;
}
