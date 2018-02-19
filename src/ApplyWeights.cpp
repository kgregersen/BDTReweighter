//local includes
#include "DecisionTree.h"
#include "Branch.h"
#include "Variable.h"
#include "Variables.h"
#include "Event.h"
#include "Store.h"
#include "Log.h"
#include "Method.h"

// stl includes
#include <vector>
#include <string>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <utility>

// ROOT includes
#include "TFile.h"
#include "TTree.h"



int main(int argc, char * argv[]) {

  // get confiuration file
  std::string configpath = "config_apply.txt";
  Store * config = Store::createStore(configpath.c_str());

  // initialize log
  Log log("ApplyWeights");
  std::string str_level;
  config->getif<std::string>("PrintLevel", str_level);
  if (str_level.length() > 0) {
    Log::LEVEL level = Log::StringToLEVEL(str_level);
    log.SetLevel(level);
  }

  // set method
  std::string str_method;
  config->getif<std::string>("Method", str_method);
  if (str_level.length() == 0) {
    log << Log::ERROR << "Method not specified! Syntax : 'string Method = <method-name>'. Available methods: BDT, RF." << Log::endl();
    return 0;    
  }
  Method::TYPE method = Method::Type(str_method);

  // initialize variables
  Variables::Initialize();

  // open weights file
  std::string weightsfilename = config->get<std::string>("WeightsFileName");
  std::ifstream weightsfile;
  log << Log::INFO << "Opening file " << weightsfilename << Log::endl();
  weightsfile.open(weightsfilename.c_str());
 
  // forest of decision trees
  std::vector<const DecisionTree *> trees;

  // single tree (collection of final node weights and corresponding cuts)
  std::vector<std::pair<float, std::vector<const Branch::Cut *> > > treeReadIn;
  
  // read lines
  log << Log::INFO << "Reading file " << weightsfilename << Log::endl();
  std::string line;
  int lineNumber = 0;
  while ( ! weightsfile.eof() ) {

    // increment line counter
    ++lineNumber;

    // get line 
    std::getline( weightsfile , line );
    log << Log::DEBUG << line << Log::endl();

    // check if we are at new tree
    if (line.size() >= 14 && line.substr(2,13) == "Decision Tree") {
      if (treeReadIn.size()) {
	trees.push_back( new DecisionTree(treeReadIn, config) );
      }
      treeReadIn.clear();
      continue;
    }
    
    // only consider lines that starts with 'weight='
    if ( ! (line.size() >= 7 && line.substr(0,7) == "weight=") ) continue;

    // retrieve weight
    size_t pos = line.find(':') + 1;
    std::string buffer = line.substr(7, pos - 7);
    std::istringstream iss(buffer);
    float weight;
    iss >> weight;

    // branches for this weight
    std::vector<const Branch::Cut *> cuts;

    // retrieve cuts from this line and convert them to branches
    while (pos != std::string::npos) {

      size_t next = line.find_first_of('|',pos);
      buffer = line.substr(pos, next - pos);
      
      size_t posLT = buffer.find('<');
      size_t posGT = buffer.find('>');
      float value;
      if (posLT != std::string::npos && posGT != std::string::npos) {
	log << Log::ERROR << "There is both a '<' and a '>' in line " << lineNumber << Log::endl();
	return 0;
      }
      else if (posLT != std::string::npos) {
	std::string name  = buffer.substr(0, posLT);
	std::string valueStr = buffer.substr(posLT + 1);
	std::istringstream issV(valueStr);
	issV >> value;
	cuts.push_back( new Branch::Smaller(Variables::Get(name), value) );
      }
      else if (posGT != std::string::npos) {
	std::string name  = buffer.substr(0, posGT);
	std::string valueStr = buffer.substr(posGT + 1);
	std::istringstream issV(valueStr);
	issV >> value;
	cuts.push_back( new Branch::Greater(Variables::Get(name), value) );
      }
      
      pos = (next == std::string::npos ? next : next + 1);

    }

    // reverse vector of branches so input branch is first
    std::reverse(cuts.begin(), cuts.end());

    // add final node to tree (weight, branches)
    treeReadIn.push_back( std::make_pair(weight, cuts) );
        
  }
  
  // remember to add last tree
  trees.push_back( new DecisionTree(treeReadIn, config) );

  log << Log::INFO << "Weights succesfully read from file!" << Log::endl();

  // open input file in 'update' mode
  const std::string & inputfilename = config->get<std::string>("InputFileName");
  TFile * f = new TFile(inputfilename.c_str(), "update");
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

  // create event object and connect TTrees
  Event & event = Event::Instance(config);
  event.ConnectAllVariables(initial, false);

  // BDT event weight
  float weight;
  std::string weightName = "BDTWeight";
  config->getif<std::string>("WeightName", weightName);
  TBranch * b_weight = initial->Branch(weightName.c_str(), &weight);

  // prepare for loop over tree entries
  long maxEvent = initial->GetEntries();
  long reportFrac = maxEvent/(maxEvent > 100000 ? 10 : 1) + 1;
  log << Log::INFO << "FillHistograms() : Looping over events (" << initial->GetName() << ") : "  << maxEvent << Log::endl();
  std::clock_t start = std::clock();

  // Loop over ree entries
  for (long ievent = 0; ievent < maxEvent; ++ievent) {

    // print progress
    if( ievent > 0 && ievent % reportFrac == 0 ) {
      double duration     = (std::clock() - start)/static_cast<double>(CLOCKS_PER_SEC);    
      double frequency    = static_cast<double>(ievent) / duration;
      double timeEstimate = static_cast<double>(maxEvent - ievent) / frequency;
      log << Log::INFO << "FillHistograms() : ---> processed : " << std::setw(8) << 100*ievent/maxEvent << "\%  ---  frequency : " << std::setw(7) << static_cast<int>(frequency) << " events/sec  ---  time : " << std::setw(4) << static_cast<int>(duration) << " sec  ---  remaining time : " << std::setw(4) << static_cast<int>(timeEstimate) << " sec"<< Log::endl(); 
    }
    
    // get event
    initial->GetEntry( ievent );

    // initialise weight
    if      (method == Method::BDT) weight = 1.;
    else if (method == Method::RF ) weight = 0.;

    // loop over trees
    for (const DecisionTree * t : trees) {
      
      float w = t->GetWeight();
      
      // use weight if event falls on this node, and break out of this tree
      if      (method == Method::BDT) weight *= w;
      else if (method == Method::RF ) weight += w;

    }

    // finalise weight
    if (method == Method::RF) {
      weight /= static_cast<float>( trees.size() );
    }
    
    // fill weight
    b_weight->Fill();
    
  }

  // print out
  double duration  = (std::clock() - start)/static_cast<double>(CLOCKS_PER_SEC);    
  double frequency = static_cast<double>(maxEvent) / duration;
  log << Log::INFO << "FillHistograms() : ---> processed :  100\%  ---  frequency : " << std::setw(7) << static_cast<int>(frequency) << " events/sec  ---  time : " << std::setw(4) << static_cast<int>(duration) << " sec  ---  remaining time :    0 sec"<< Log::endl(); 

  // write tree
  initial->Write();  
  
  // and we're done!
  return 0;
  
}
