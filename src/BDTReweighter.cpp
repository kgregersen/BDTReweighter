//local includes
#include "DecisionTree.h"
#include "Variable.h"
#include "Variables.h"
#include "Event.h"
#include "HistDefs.h"
#include "Store.h"
#include "Log.h"

// stl includes
#include <vector>
#include <string>
#include <fstream>
#include <ctime>

// ROOT includes
#include "TFile.h"
#include "TTree.h"



int main(int argc, char * argv[]) {

  // get confiuration file
  std::string configpath = "config.txt";
  Store * config = Store::createStore(configpath.c_str());

  // initialize log
  Log log("BDTRweighter");
  std::string str_level;
  config->getif<std::string>("PrintLevel", str_level);
  if (str_level.length() > 0) {
    Log::LEVEL level = Log::StringToLEVEL(str_level);
    log.SetLevel(level);
  }

  // get input file
  const std::string & inputfilename = config->get<std::string>("InputFileName");
  TFile * f = new TFile(inputfilename.c_str(), "read");
  if ( ! f->IsOpen() ) {
    log << Log::ERROR << "Couldn't open file : " << inputfilename << Log::endl();
    return 0;
  }
  
  // get initial tree
  const std::string & treenameinitial = config->get<std::string>("InputTreeNameInitial");
  TTree * initial = static_cast<TTree *>(f->Get(treenameinitial.c_str()));
  if ( ! initial ) {
    log << Log::ERROR << "Couldn't get TTree : " << treenameinitial << Log::endl();
    return 0;
  }

  // get target tree
  const std::string & treenametarget = config->get<std::string>("InputTreeNameTarget");
  TTree * target = static_cast<TTree *>(f->Get(treenametarget.c_str()));
  if ( ! target ) {
    log << Log::ERROR << "Couldn't get TTree : " << treenametarget << Log::endl();
    return 0;
  }

  // create event object and connect TTrees
  Event & event = Event::Instance(config);
  event.ConnectAllVariables(initial);
  event.ConnectAllVariables(target);

  // initialize variables
  Variables::Initialize();
  
  // get histogram definitions
  HistDefs histDefs(config);
  histDefs.Initialize();
  histDefs.UpdateVariableRanges(target);
  histDefs.UpdateVariableRanges(initial);
  for (const HistDefs::Entry & entry : histDefs.GetEntries()) {
    log << Log::INFO << "Histogram name : " << entry.Name() << ", range = ( " << entry.Xmin() << " , " << entry.Xmax() << " )" << Log::endl();
  }
   
  // grow decision trees
  int ntree = config->get<int>("NumberOfTrees");
  std::vector<const DecisionTree *> decisionTrees;
  for (int itree = 0; itree < ntree; ++itree) {

    // create tree
    DecisionTree * dtree = new DecisionTree(initial, target, config, histDefs);
    dtree->GrowTree(decisionTrees);
    
    // add tree to forest
    decisionTrees.push_back(dtree);

  }

  // open ouput file
  std::string outfilename = "BDTweights.txt";
  config->getif<std::string>("OutputFileName", outfilename);
  std::ofstream outfile;
  outfile.open(outfilename.c_str());

  // print timestamp to file
  std::time_t now= std::time(0);
  std::tm* now_tm= std::gmtime(&now);
  char buf[200];
  std::strftime(buf, 200, "%a, %d %b %y %T %z", now_tm);
  outfile << "# Time stamp : " << buf << "\n";

  // print variables to file
  outfile << " # Variables  : ";
  const std::vector<const Variable *> & variables = Variables::Get();
  for (const Variable * var : variables) {
    outfile << var->Name() << ",";
  }
  outfile << "\n";
  
  // write decision trees to file
  for (const DecisionTree * dtree : decisionTrees) {
    dtree->Write( outfile );
  }

  // close file
  outfile.close();
  
  
  // and we're done!
  return 0;
  
}
