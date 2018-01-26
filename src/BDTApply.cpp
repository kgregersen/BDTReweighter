//local includes
#include "Branch.h"
#include "Variable.h"
#include "Variables.h"
#include "Event.h"
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
  Log log("BDTApply");
  std::string str_level;
  config->getif<std::string>("PrintLevel", str_level);
  if (str_level.length() > 0) {
    Log::LEVEL level = Log::StringToLEVEL(str_level);
    log.SetLevel(level);
  }

  // initialize variables
  Variables::Initialize();

  // open BDT weights file
  std::string weightsfilename = "BDTweights.txt";
  config->getif<std::string>("WeightsFileName", weightsfilename);
  std::ifstream weightsfile;
  log << Log::INFO << "Opening file " << weightsfilename << Log::endl();
  weightsfile.open(weightsfilename.c_str());

  // vector containing weights and associated cuts
  std::vector<std::pair<float, std::vector<const Branch *> > > weights;

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

    // only consider lines that starts with 'weight='
    if (line.size() == 0 || line.substr(0,7) != "weight=") continue;

    // retrieve weight
    size_t pos = line.find(':') + 1;
    std::string buffer = line.substr(7, pos - 7);
    std::istringstream iss(buffer);
    float weight;
    iss >> weight;

    std::vector<const Branch *> branches;
    
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
	branches.push_back( new Branch(config, 0, name, value, false, 0, 0) );
      }
      else if (posGT != std::string::npos) {
	std::string name  = buffer.substr(0, posGT);
	std::string valueStr = buffer.substr(posGT + 1);
	std::istringstream issV(valueStr);
	issV >> value;
	branches.push_back( new Branch(config, 0, name, value, true, 0, 0) );
      }
      
      pos = next == std::string::npos ? next : next + 1;

    }

    weights.push_back( std::make_pair(weight, branches) );
        
  }

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
  event.ConnectAllVariables(initial);

  // BDT event weight
  float BDTEventWeight;
  std::string BDTEventWeightName = "BDTWeight";
  config->getif<std::string>("BDTEventWeightName", BDTEventWeightName);
  TBranch * b_BDTEventWeight = initial->Branch(BDTEventWeightName.c_str(), &BDTEventWeight);

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
      log << Log::INFO << "FillHistograms() : ---> processed : " << std::setw(4) << 100*ievent/maxEvent << "\%  ---  frequency : " << std::setw(7) << static_cast<int>(frequency) << " events/sec  ---  time : " << std::setw(4) << static_cast<int>(duration) << " sec  ---  remaining time : " << std::setw(4) << static_cast<int>(timeEstimate) << " sec"<< Log::endl(); 
    }
    
    // get event
    initial->GetEntry( ievent );

    // apply cuts
    BDTEventWeight = 1.;
    for (unsigned int i = 0; i < weights.size(); ++i) {

      // check if event falls on the i'th node
      bool pass = true;
      for (const Branch * b : weights.at(i).second) {
	if ( ! b->Pass() ) {
	  pass = false;
	  break;
	}
      }

      // use weight if event falls on this node
      if ( pass ) BDTEventWeight *= weights.at(i).first;     
      
    }

    // fill BDT event weight
    b_BDTEventWeight->Fill();
    
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
