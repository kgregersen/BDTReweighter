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


// 
class Forest {

public:
  
  // constructor
  Forest(const std::vector<const DecisionTree *> trees, Method::TYPE method) : m_trees(trees), m_method(method), m_weight(0), m_error(0) {}

  // destructor
  ~Forest() {}

  // get weight
  void GetWeight(float & weight, float & error) 
  {   
    CalculateWeight();
    weight = m_weight;
    error  = m_error;
  }


private:

  // calculate weight
  void CalculateWeight() {
    if (m_method == Method::BDT) {
      for (unsigned int i = 0; i < m_trees.size(); ++i) {
	float w = m_trees.at(i)->GetWeight();
	m_weight *= w;
      }
      m_error = 0.;
    }
    else if (m_method == Method::RF || m_method == Method::ET) {
      static std::vector<float> error_vec;
      error_vec.resize(m_trees.size());
      for (unsigned int i = 0; i < m_trees.size(); ++i) {
	float w = m_trees.at(i)->GetWeight();
	m_weight += w;
	error_vec.at(i) = w;
      }
      m_weight /= static_cast<float>( m_trees.size() );
      for (unsigned int i = 0; i < m_trees.size(); ++i) {
	m_error += pow(error_vec.at(i) - m_weight, 2);
      }
      m_error = sqrt( m_error/( m_trees.size() > 1 ? m_trees.size() - 1 : 1 ) );
    }
  }
  
  // trees
  const std::vector<const DecisionTree *> m_trees;
 
  // method
  const Method::TYPE m_method;

  // weight/error
  float m_weight;
  float m_error;

};




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
    log << Log::ERROR << "Method not specified! Syntax : 'string Method = <method-name>'. Available methods: BDT, RF, ET (see ./inc/Methods.h)." << Log::endl();
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
 
  // vector of forests
  std::vector<Forest *> forests;

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

    // check if we are at new forest
    if ( line.size() >= 10 && line.substr(0,10) == "Time stamp" ) {
      if (trees.size()) {
	forests.push_back( new Forest(trees, method) );
	trees.clear();
      }
    }

    // check if we are at new tree
    if ( (line.size() >= 14 && line.substr(2,13) == "Decision Tree") || (line.size() >= 5 && line.substr(2,3) == "End") ) {
      if (treeReadIn.size()) {
	trees.push_back( new DecisionTree(treeReadIn, config) );
      }
      for (unsigned int i = 0; i < treeReadIn.size(); ++i) {
	std::vector<const Branch::Cut *> cuts = treeReadIn.at(i).second;
	for (unsigned int j = 0; j < cuts.size(); ++j) {
	  delete cuts.at(j);
	}
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
   
  // remember to add last forest
  forests.push_back( new Forest(trees, method) );
  
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

  // BDT/RF/ET event weight
  float weight;
  std::string weightName = "Weight";
  config->getif<std::string>("WeightName", weightName);
  TBranch * b_weight = initial->Branch(weightName.c_str(), &weight);

  // RF/ET event weight error
  float weight_err;
  std::string weightErrName = weightName + "_err";
  TBranch * b_weight_err = 0;
  bool bagging = false;
  config->getif<bool>("Bagging", bagging);
  if (method == Method::RF || method == Method::ET || bagging) b_weight_err = initial->Branch(weightErrName.c_str(), &weight_err);
  std::vector<float> weight_err_vec(forests.size());
  
  // prepare for loop over tree entries
  long maxEvent = initial->GetEntries();
  long reportFrac = maxEvent/(maxEvent > 100000 ? 100 : 1) + 1;
  log << Log::INFO << "Looping over events (" << initial->GetName() << ") : "  << maxEvent << Log::endl();
  std::clock_t start = std::clock();

  // Loop over ree entries
  for (long ievent = 0; ievent < maxEvent; ++ievent) {

    // print progress
    if( ievent > 0 && ievent % reportFrac == 0 ) {
      double duration     = (std::clock() - start)/static_cast<double>(CLOCKS_PER_SEC);    
      double frequency    = static_cast<double>(ievent) / duration;
      double timeEstimate = static_cast<double>(maxEvent - ievent) / frequency;
      log << Log::INFO << "---> processed : " << std::setw(8) << 100*ievent/maxEvent << "\%  ---  frequency : " << std::setw(7) << static_cast<int>(frequency) << " events/sec  ---  time : " << std::setw(4) << static_cast<int>(duration) << " sec  ---  remaining time : " << std::setw(4) << static_cast<int>(timeEstimate) << " sec"<< Log::endl(); 
    }
    
    // get event
    initial->GetEntry( ievent );

    // get weight/error
    weight = 0;
    weight_err = 0;
    for (unsigned int f = 0; f < forests.size(); ++f) {
      Forest * forest = forests.at(f);
      float w = 0;
      float e = 0;
      forest->GetWeight(w, e);
      weight += w;
      weight_err_vec.at(f) = (e > 0 ? e : w);
    }
    
    // finalise weight
    weight /= static_cast<float>( forests.size() );
    for (unsigned int f = 0; f < forests.size(); ++f) {
      weight_err += pow(weight_err_vec.at(f) - weight, 2);
    }
    weight_err = sqrt( weight_err/( forests.size() > 1 ? forests.size() - 1 : 1 ) );
    
    // fill weight
    b_weight->Fill();
    if (b_weight_err) {
      b_weight_err->Fill();
    }
    
  }

  // print out
  double duration  = (std::clock() - start)/static_cast<double>(CLOCKS_PER_SEC);    
  double frequency = static_cast<double>(maxEvent) / duration;
  log << Log::INFO << "---> processed : " << std::setw(8) << 100 << "\%  ---  frequency : " << std::setw(7) << static_cast<int>(frequency) << " events/sec  ---  time : " << std::setw(4) << static_cast<int>(duration) << " sec  ---  remaining time :    0 sec"<< Log::endl(); 

  // write tree
  initial->Write();  
  
  // and we're done!
  return 0;
  
}
