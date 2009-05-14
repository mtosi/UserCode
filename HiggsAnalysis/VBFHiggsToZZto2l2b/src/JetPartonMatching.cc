#include "HiggsAnalysis/VBFHiggsToZZto2l2b//interface/JetPartonMatching.h"

#include "HiggsAnalysis/VBFHiggsToZZto2l2b/interface/VBFHZZllbbUtils.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <Math/VectorUtil.h>

JetPartonMatching::JetPartonMatching(const std::vector<const reco::Candidate*>& partonVec, const std::vector<reco::GenJet>& jetVec,
				     const int algorithm = totalMinDist, const bool useMaxDist = true, 
				     const bool useDeltaR = true, const double maxDist = 0.3)
  : partons(partonVec), algorithm_(algorithm), useMaxDist_(useMaxDist), useDeltaR_(useDeltaR), maxDist_(maxDist)
{
  std::vector<const reco::Candidate*> js;
  for(unsigned int index=0; index != jetVec.size(); index++)
    js.push_back( &(jetVec[index]) );
  jets = js;
  calculate();
}

JetPartonMatching::JetPartonMatching(const std::vector<const reco::Candidate*>& partonVec, const std::vector<reco::CaloJet>& jetVec,
				     const int algorithm = totalMinDist, const bool useMaxDist = true,
				     const bool useDeltaR = true, const double maxDist = 0.3)
  : partons(partonVec), algorithm_(algorithm), useMaxDist_(useMaxDist), useDeltaR_(useDeltaR), maxDist_(maxDist) {

  std::vector<const reco::Candidate*> js;
  for(unsigned int index = 0; index < jetVec.size(); ++index)
    js.push_back( &(jetVec[index]) );
  jets = js; 
  calculate(); 
}

JetPartonMatching::JetPartonMatching(const std::vector<const reco::Candidate*>& partonVec, const std::vector<const reco::Candidate*>& jetVec,
				     const int algorithm = totalMinDist, const bool useMaxDist = true,
				     const bool useDeltaR = true, const double maxDist = 0.3)
  : partons(partonVec), jets(jetVec), algorithm_(algorithm), useMaxDist_(useMaxDist), useDeltaR_(useDeltaR), maxDist_(maxDist) {

  calculate();
}

void 
JetPartonMatching::calculate() {
  // use maximal distance between objects in case of unambiguousOnly algorithmm
  if(algorithm_ == unambiguousOnly) useMaxDist_=true;

  // check if there are empty partons in the partonVector
  bool emptyParton = false;
  for ( unsigned int partonIndex = 0; partonIndex != partons.size(); partonIndex++ ){
    if ( partons[partonIndex]->pdgId() == 0 ) {
      emptyParton = true;
      break;
    }
  }

  // switch algorithm,
  // default is to match on the minimal sum of the distance 
  // (if jets or a parton is empty fill match with blanks)
  if( jets.empty() || emptyParton ) {
    MatchingCollection dummyMatch;
    for ( unsigned int partonIndex = 0; partonIndex != partons.size(); partonIndex++ )
      dummyMatch.push_back( std::make_pair(partonIndex, -1) );
    matching.push_back( dummyMatch );
  }
  else {
    switch( algorithm_ ) {
    case totalMinDist: 
      matchingTotalMinDist();    
      break;
      
    case minSumDist: 
      matchingMinSumDist();
      break;
      
    case ptOrderedMinDist: 
      matchingPtOrderedMinDist();
      break;
      
    case unambiguousOnly:
      matchingUnambiguousOnly();
      break;
      
    default:
      matchingMinSumDist();
    }
  }

  numberOfUnmatchedPartons.clear();
  sumDeltaE .clear();
  sumDeltaPt.clear();
  sumDeltaR .clear();
  for(unsigned int comb = 0; comb < matching.size(); comb++) {
    MatchingCollection match = matching[comb];
    std::sort(match.begin(), match.end());
    matching[comb] = match;
    
    int nUnmatchedPartons = partons.size();
    for(unsigned int partonIndex = 0; partonIndex != partons.size(); partonIndex++ )
      if(getMatchForParton(partonIndex,comb) >= 0) --nUnmatchedPartons;

    double sumDE  = -999.;
    double sumDPt = -999.;
    double sumDR  = -999.;
    if(nUnmatchedPartons == 0){
      sumDE  = 0;
      sumDPt = 0;
      sumDR  = 0;
      for(unsigned int matchIndex = 0; matchIndex != match.size(); matchIndex++ ){
	sumDE  += fabs(partons[match[matchIndex].first]->energy() - jets[match[matchIndex].second]->energy());
	sumDPt += fabs(partons[match[matchIndex].first]->pt()     - jets[match[matchIndex].second]->pt());
	sumDR  += vbfhzz2l2b::distance(partons[match[matchIndex].first]->p4(), jets[match[matchIndex].second]->p4(), useDeltaR_);
      }
    }

    numberOfUnmatchedPartons.push_back( nUnmatchedPartons );
    sumDeltaE .push_back( sumDE  );
    sumDeltaPt.push_back( sumDPt );
    sumDeltaR .push_back( sumDR  );
  }
}
//-----------------------------------------------------------------------------

// match parton to jet
// with shortest distance starting with the shortest distance available
// apply some outlier rejection if desired
void 
JetPartonMatching::matchingTotalMinDist() {

  // prepare vector of pairs
  // with distances between all partons to all jets
  // in the input vectors
  std::vector< std::pair<double, unsigned int> > distances;
  for(unsigned int partonIndex = 0; partonIndex != partons.size(); partonIndex++) {
    for(unsigned int jetIndex = 0; jetIndex != jets.size(); jetIndex++ ) { 
      double dist = vbfhzz2l2b::distance(jets[jetIndex]->p4(), partons[partonIndex]->p4(),useDeltaR_);
      distances.push_back(std::pair<double, unsigned int>(dist, partonIndex*jets.size()+jetIndex));
    }
  }
  std::sort(distances.begin(), distances.end());

  MatchingCollection match;

  while(match.size() < partons.size()){
    unsigned int partonIndex = distances[0].second/jets.size();
    int jetIndex = distances[0].second-jets.size()*partonIndex;
    
    // use primitive outlier rejection if desired
    if( useMaxDist_&& distances[0].first > maxDist_ ) jetIndex = -1;

    // prevent underflow in case of too few jets
    if( distances.empty() )
      match.push_back(std::make_pair(partonIndex, -1));
    else
      match.push_back(std::make_pair(partonIndex, jetIndex));
    
    // remove all values for the matched parton and the matched jet
    for(unsigned int distIndex = 0; distIndex < distances.size(); distIndex++) {
      unsigned int pIndex = distances[distIndex].second/jets.size();
      int jIndex = distances[distIndex].second-jets.size()*pIndex;
      if((pIndex == partonIndex) || (jIndex == jetIndex)) {
	distances.erase(distances.begin()+distIndex, distances.begin()+distIndex+1); 
	--distIndex;
      }
    }
  }
  matching.clear();
  matching.push_back( match );
  return;
}
//-----------------------------------------------------------------------------
void 
JetPartonMatching::minSumDist_recursion(const unsigned int partonIndex,
					std::vector<unsigned int> & jetIndices,
					std::vector<bool> & usedJets,
					std::vector<std::pair<double, MatchingCollection> > & distMatchVec) {

  // build up jet combinations recursively
  if(partonIndex != partons.size()){
    for(unsigned int jetIndex = 0; jetIndex != jets.size(); jetIndex++ ){
      if(usedJets[jetIndex]) continue;
      usedJets[jetIndex] = true;
      jetIndices[partonIndex] = jetIndex;
      minSumDist_recursion(partonIndex+1, jetIndices, usedJets, distMatchVec);
      usedJets[jetIndex] = false;
    }
    return;
  }

  // calculate sumDist for each completed combination
  double sumDist = 0;
  MatchingCollection match;
  for(unsigned int partonIndex = 0; partonIndex != partons.size(); partonIndex++){
    double dist  = vbfhzz2l2b::distance(partons[partonIndex]->p4(), jets[jetIndices[partonIndex]]->p4(),useDeltaR_);
    if(useMaxDist_ && dist > maxDist_) return; // outlier rejection
    sumDist += vbfhzz2l2b::distance(partons[partonIndex]->p4(), jets[jetIndices[partonIndex]]->p4(),useDeltaR_);
    match.push_back(std::make_pair(partonIndex, jetIndices[partonIndex]));
  }

  distMatchVec.push_back( std::make_pair(sumDist, match)  );
  return;
}
//-----------------------------------------------------------------------------

// match partons to jets
// with minimal sum of the distances between all partons and jets
void JetPartonMatching::matchingMinSumDist() {

  std::vector<std::pair<double, MatchingCollection> > distMatchVec;

  std::vector<bool> usedJets;
  for(unsigned int jetIndex = 0; jetIndex < jets.size(); jetIndex++ )
    usedJets.push_back(false);

  std::vector<unsigned int> jetIndices;
  jetIndices.reserve(partons.size());

  minSumDist_recursion(0, jetIndices, usedJets, distMatchVec);
  std::sort(distMatchVec.begin(), distMatchVec.end());

  matching.clear();
  for(unsigned int i=0; i<distMatchVec.size(); ++i)
    matching.push_back( distMatchVec[i].second );  
  return;
}
//-----------------------------------------------------------------------------

// match partons to jets
// with minimal sum of the distances between all partons and jets
// order partons in pt first
void 
JetPartonMatching::matchingPtOrderedMinDist() {
  std::vector<std::pair <double, unsigned int> > ptOrderedPartons;

  for(unsigned int partonIndex=0; partonIndex != partons.size(); partonIndex++)
    ptOrderedPartons.push_back(std::make_pair(partons[partonIndex]->pt(), partonIndex));

  std::sort(ptOrderedPartons.begin(), ptOrderedPartons.end());
  std::reverse(ptOrderedPartons.begin(), ptOrderedPartons.end());

  std::vector<unsigned int> jetIndices;
  for(unsigned int jetIndex = 0; jetIndex != jets.size(); jetIndex++)
    jetIndices.push_back(jetIndex);

  MatchingCollection match;

  for(unsigned int partonIndex = 0; partonIndex != ptOrderedPartons.size(); partonIndex++ ) {
    double minDist = 999.;
    int jetIndexMin = -1;

    for(unsigned int jetIndex = 0; jetIndex != jetIndices.size(); jetIndex++ ){
      double dist = vbfhzz2l2b::distance(partons[ptOrderedPartons[partonIndex].second]->p4(), jets[jetIndices[jetIndex]]->p4(),useDeltaR_);
      if(dist < minDist){
	if(!useMaxDist_ || dist <= maxDist_) {
	  minDist = dist;
	  jetIndexMin = jetIndex;
	}
      }
    }
    
    if(jetIndexMin >= 0){
      match.push_back( std::make_pair(ptOrderedPartons[partonIndex].second, jetIndices[jetIndexMin]) );
      jetIndices.erase(jetIndices.begin() + jetIndexMin, jetIndices.begin() + jetIndexMin + 1);
    }
    else
      match.push_back( std::make_pair(ptOrderedPartons[partonIndex].second, -1) );
  }

  matching.clear();
  matching.push_back( match );
  return;
}
//-----------------------------------------------------------------------------

// match partons to jets,
// only accept event 
// if there are no ambiguities
void 
JetPartonMatching::matchingUnambiguousOnly() {
  std::vector<bool> jetMatched;
  for(unsigned int jetIndex = 0; jetIndex != jets.size(); jetIndex++ )
    jetMatched.push_back(false);

  MatchingCollection match;
  
  for(unsigned int partonIndex = 0; partonIndex != partons.size(); partonIndex++){
    int iMatch = -1;
    for(unsigned int jetIndex = 0; jetIndex != jets.size(); jetIndex++ ){
      double dist = vbfhzz2l2b::distance(partons[partonIndex]->p4(), jets[jetIndex]->p4(),useDeltaR_);
      if(dist <= maxDist_){
	if(!jetMatched[jetIndex]){ // new match for jet
	  jetMatched[jetIndex] = true;
	  if(iMatch == -1) // new match for parton and jet
	    iMatch = jetIndex;
	  else // parton already matched: ambiguity!
	    iMatch = -2;
	}
	else // jet already matched: ambiguity!
	  iMatch = -2;
      }
    }
    match.push_back(std::make_pair(partonIndex, iMatch));
  }

  matching.clear();
  matching.push_back( match );
  return;
}
//-----------------------------------------------------------------------------

int
JetPartonMatching::getMatchForParton(const unsigned int part, const unsigned int comb) {
  if(comb >= matching.size()) return -9;
  if(part >= matching[comb].size()) return -9;
  return (matching[comb])[part].second;
}

// return a vector with the indices of the matched jets
// (ordered according to the vector of partons)
std::vector<int>
JetPartonMatching::getMatchesForPartons(const unsigned int comb) {
  std::vector<int> jetIndices;
  for(unsigned int partonIndex = 0; partonIndex != partons.size(); partonIndex ++)
    jetIndices.push_back( getMatchForParton(partonIndex, comb) );
  return jetIndices;
}

double 
JetPartonMatching::getDistanceForParton(const unsigned int part, const unsigned int comb) {
  // get the distance between parton and its best matched jet
  if(getMatchForParton(part, comb) < 0) return -999.;
  return vbfhzz2l2b::distance( jets[getMatchForParton(part,comb)]->p4(), partons[part]->p4(), useDeltaR_ );
}

double 	
JetPartonMatching::getSumDistances(const unsigned int comb) {
  // get sum of distances between partons and matched jets
  double sumDists = 0.;
  for(unsigned int partonIndex = 0; partonIndex != partons.size(); partonIndex++){
    double dist = getDistanceForParton(partonIndex, comb);
    if(dist < 0.) return -999.;
    sumDists += dist;
  }
  return sumDists;
}

void
JetPartonMatching::print() {
  //report using MessageLogger
  edm::LogInfo log("JetPartonMatching");
  log << "++++++++++++++++++++++++++++++++++++++++++++++ \n";
  log << " algorithm : ";
  switch(algorithm_) {
    case totalMinDist     : log << "totalMinDist    "; break;
    case minSumDist       : log << "minSumDist      "; break;
    case ptOrderedMinDist : log << "ptOrderedMinDist"; break;
    case unambiguousOnly  : log << "unambiguousOnly "; break;
    default               : log << "UNKNOWN         ";
  }
  log << "\n";
  log << " useDeltaR : ";
  switch(useDeltaR_) {
    case false : log << "false"; break;
    case true  : log << "true ";
  }
  log << "\n";
  log << " useMaxDist: ";
  switch(useMaxDist_) {
    case false : log << "false"; break;
    case true  : log << "true ";
  }
  log << "      maxDist: " << maxDist_ << "\n";
  log << " number of partons / jets: " << partons.size() << " / " << jets.size() << "\n";
  log << " number of available combinations: " << getNumberOfAvailableCombinations() << "\n";
  for(unsigned int comb = 0; comb < matching.size(); ++comb) {
    log << " -------------------------------------------- \n";
    log << " ind. of matched jets:";
    for(unsigned int partonIndex = 0; partonIndex != partons.size(); partonIndex++ )
      log << std::setw(4) << getMatchForParton(partonIndex, comb);
    log << "\n";
    log << " sumDeltaR             : " << getSumDeltaR(comb) << "\n";
    log << " sumDeltaPt / sumDeltaE: " << getSumDeltaPt(comb) << " / "  << getSumDeltaE(comb);
    log << "\n";
  }
  log << "++++++++++++++++++++++++++++++++++++++++++++++";
}
