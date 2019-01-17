/*
  Substructure Searching over Reactions
*/

#include <fstream>
#include <string>
#include <iostream>

#include <memory>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::ofstream;

#include "cmdline.h"
#include "misc.h"

#include "istream_and_type.h"
#include "rwsubstructure.h"
#include "misc2.h"
#include "molecule.h"
#include "aromatic.h"
#include "iwstandard.h"
#include "substructure.h"
#include "target.h"
#include "path.h"

typedef enum RxnQueryType
{
	Smiles,
	Smarts,
	QryFile
} RxnQueryType;
    
class Reaction_To_Search
{
  private:
    vector<Molecule *> reactants;
    vector<Molecule *> agents;
    vector<Molecule *> products;

    IWString _name;
    
    Chemical_Standardisation &chemicalStandardization;
		bool &reduceToLargestFragment;


//  private functions

    int buildRxnPart(const_IWSubstring s, vector<Molecule *> & m);

  public:
    Reaction_To_Search(Chemical_Standardisation &chemicalStandardizationInit
    								, bool reduceToLargestFragmentInit ):
    	    chemicalStandardization(chemicalStandardizationInit)
					,reduceToLargestFragment(reduceToLargestFragmentInit)

    {
    }
    ~Reaction_To_Search()
    {
    	for (vector<Molecule *>::iterator thisIter = reactants.begin() ;  thisIter != reactants.end() ; ++thisIter)
    		delete *thisIter;
    	for (vector<Molecule *>::iterator thisIter = agents.begin() ;  thisIter != agents.end() ; ++thisIter)
    		delete *thisIter;
    	for (vector<Molecule *>::iterator thisIter = products.begin() ;  thisIter != products.end() ; ++thisIter)
    		delete *thisIter;
    }

    int build(const const_IWSubstring & buffer);

    vector<Molecule *> & getReactants() { return reactants;}
    vector<Molecule *> & getAgents() { return agents;}
    vector<Molecule *> & getProducts() { return products;}
};

int
Reaction_To_Search::build(const const_IWSubstring & buffer)
{
  const_IWSubstring s, n;
  if (! buffer.split(s, ' ', n))
  {
  	s = buffer;
  	n = "";
  }
  
  if (0 == s.length())
  {
    cerr << "Reaction_To_Search::build:invalid input '" << buffer << "'\n";
    return 0;
  }

  if (2 != s.ccount('>'))
  {
    cerr << "Reaction_To_Search::build:reaction must have three components separated by >, " << buffer << " invalid\n";
    return 0;
  }

  _name = n;

  const_IWSubstring r, a, p;

  int i = 0;
  s.nextword_single_delimiter(r, i, '>');
  s.nextword_single_delimiter(a, i, '>');
  s.nextword_single_delimiter(p, i, '>');

  if (! buildRxnPart(r, reactants) || ! buildRxnPart(a, agents) || ! buildRxnPart(p, products))
  {
    //cerr << "Reaction_To_Search::build:invalid input '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

int
Reaction_To_Search::buildRxnPart(const_IWSubstring s, vector<Molecule *> & m)
{
	m.clear();
	
  if (0 == s.length())    // nothing to do
    return 1;

//  resizable_array<int> pos;
//
//  identify_plus_positions(s, pos);
//  pos.add(s.length() - 1);
//
//  const int n = pos.number_elements();
  //iwaray<const_IWSubstring>parts;
	resizable_array_p<const_IWSubstring> parts;
//	s.split(parts,'+');
//  const int n = parts.number_elements();
 	const int n=splitOnPlusses(s,parts);
 	
  for (int i = 0; i < n; ++i)
  {
//    int cstart;         
//    if (0 == i)
//      cstart = 0;
//    else
//      cstart = pos[i-1] + 1;
//
//    int cstop;
//
//    if (i == n - 1)
//      cstop = s.length() - 1;
//    else
//      cstop = pos[i] - 1;
//
//    const_IWSubstring smiles;
//
//    s.from_to(cstart, cstop, smiles);

    Molecule * newMol = new Molecule();

//    if (! newMol->build_from_smiles(smiles))
    if (! newMol->build_from_smiles(*(parts[i])))
    {
      cerr << "buildRxnPart:invalid smiles " << s << "\n";
      delete newMol;
      return 0;
    }
 		
 		if (reduceToLargestFragment)
    	newMol->reduce_to_largest_fragment();
  	if (chemicalStandardization.active())
    	chemicalStandardization.process(*newMol);
    m.push_back(newMol);
  }

  return n;
}

// The Reaction_Query_Component class is one of zero or more components of a reaction 
// query part (reactants, agents, products).  This component supports two types of input:
// 1) the original form form Ian Watson that had the >'s and +'s of the reaction specification on 
// 		the command line in the form -q "qry+qry+qry>qry+qry>qry+qry+qry"
//		each part can have any number of qry's, including none.  These has to be at least one qry
//		somewhere in the overall specification.
//		Each qry represents a list of Substructure_Query class instances.  A hit reaction to the overall 
//		query must have a match to at least one of the members of the list of Substructure_Query for each
//		qry specified on the command line.
//	  the qry itself can be represented as one of:
//					 a smarts string 
//					 a qry file
//					 a file of smarts strings
//					 a file of qry filenames
//		Use of multiple Substructure_Query in any one component can be easily specified, but care must 
//		be taken to avoid unexpected hits.  For example, if the user wanted to search for hits to two 
//		different reactions, like:
//			C-[NH]-CCC>>C-[N](=O)-CCC
//			CN=[N+]=[N-]>>CN
//
// 		and this was formulated as -q 'reactant1.smarts>>product1.smarts',
//		with reactant1.smarts as:
//				C-[NH]-CCC R1
//				CN=[N+]=[N-] R2
//		and product1.smarts as:
//				C-[N](=O)-CCC  P1
//				CN P2
//		then the desired hits would be found, but so would the unexpected reaction hits that contain
//		C-[NH]-CCC in the reactants and CN in the products. For cases like this, the user should use
//		the newer, more flexible entry method (vide infra).  This input type is maintained for backward
//		compatibility.
//
//	2) the revised entry method, in which the >'s and +'s of the reaction specification are part of
//		representation of the query (as a string or the lines of a file).  
//		This is formulated as -s qry2, where qry2 can be:
//			a Reaction SMARTS string of the form 'smarts+smarts+smarts>smarts+smarts>smarts+smarts'
//			a Reactions SMILES string of the form 's:smiles+smiles+smiles>smiles+smiles>smiles+smiles'
//			a file of Reaction SMARTS specified as 'S:filename
//			a file of Reaction SMILES specified as 'M:filename
//			a Reaction Qry filename, as 'Q:filename.  The Reaction Qry file contains lines with the  
//				format:
//				file.qry+file.qry+file.qry>file.qry+file.qry>file.qry+file.qry
//				where file.qry is a qry filename
//
//	the eaction_Query_Component class and it partent Reaction_Query class supports both entry methods.
//
//  In the original method, the Reaction_Query_Component 
//	object could contain multiple Substructure_Query items and there will be only one (1) 
//	master Reaction_Query object.  
//
//	In the newer method, the Reaction_Query_Component will contain only one (1) 
//	Substructure_Query item, but can contain multipe  master Reaction_Query items.
		

class Reaction_Query_Component
{
  private:
    vector<Substructure_Query *> queries;
												
  public:
    Reaction_Query_Component()
    {
    }
    ~Reaction_Query_Component()
    {
    	for (vector<Substructure_Query *>::iterator thisIter = queries.begin() ;  thisIter != queries.end() ; ++thisIter)
    		delete *thisIter;
    }
    
		int addQueryFromString(const_IWSubstring s, RxnQueryType typeFlag);
		int addQueryOldFormat(const const_IWSubstring & buffer);
	
  	int matches(Molecule &m);

};


int Reaction_Query_Component::matches(Molecule &m)
{	 			
 		Molecule_to_Match target(&m);
		Substructure_Results sresults;
		for(vector<Substructure_Query *>::iterator queryIter = queries.begin() 
				; queryIter != queries.end()
				; ++ queryIter)
		{
	  	if ((*queryIter)->substructure_search(target, sresults))
	  		return 1;
		}
  	
  	return 0;
}

int Reaction_Query_Component::addQueryFromString(const_IWSubstring buffer, RxnQueryType typeFlag)
{
	// the buffer might contain a name for the query - if so discard it
	
	const_IWSubstring s,n;
	
	if (!buffer.split(s, ' ', n))
    s = buffer;
  else  if (0 == s.length())
  {
    cerr << "RReaction_Query_Component::addQueryFromString:invalid input '" << buffer << "'\n";
    return 0;
  }
  
  // make the new query
  
  Substructure_Query * q = new Substructure_Query();

  if (typeFlag == RxnQueryType::Smiles)
  {
    if (! q->create_from_smiles(s))
    {
      cerr << "Reaction_Query_Component::addQueryFromString:invalid smiles '" << s << "'\n";
      delete q;
      return 0;    	
    }
	}
	else if(typeFlag == RxnQueryType::Smarts)
	{
		if (! q->create_from_smarts(s))
	  {
	    cerr << "Reaction_Query_Component::addQueryFromString:invalid smarts '" << s << "'\n";
	    delete q;
	    return 0;
	  }			
	}
	else if (typeFlag == RxnQueryType::QryFile)
	{
		if (! q->read(s))
	  {
	    cerr << "Reaction_Query_Component::addQueryFromString:invalid qry file '" << s << "'\n";
	    delete q;
	    return 0;
	  }		
	} 

  queries.push_back(q);

  return 1;
}

int Reaction_Query_Component::addQueryOldFormat(const const_IWSubstring & buffer)
{
  // the buffer might contain a name for the query - if so discard it
	
  const_IWSubstring s,n;
	
  if (!buffer.split(s, ' ', n))
    s = buffer;
  else  if (0 == s.length())
  {
    cerr << "Reaction_Query_Component::addQueryFromString:invalid input '" << buffer << "'\n";
    return 0;
  }
  
  // make the new query
  //create the query based on the flags at the front of the string
	
  RxnQueryType typeFlag;
  if (s.starts_with("Q:") || s.starts_with("q:"))  
    typeFlag = RxnQueryType::QryFile;
  else
    typeFlag = RxnQueryType::Smarts;
 	
  // if this is a file type, open the file and read all the lines
  // parse them according to the query type determined above

  if (s.starts_with("Q:") || s.starts_with("S:")) 
  {
    s.remove_leading_chars(2); 
    
    iwstring_data_source input(s);

    if (! input.good())
    {
      cerr << ": cannot open '" << s << "'\n";
      return 0;
    }

    const_IWSubstring buffer2;  // buffer is input, should not be redeclared here

    while (input.next_record(buffer2))
    {
      buffer2.remove_line_terminators();
      // OK - we have a Smarts or a full path
      if (typeFlag == RxnQueryType::Smarts || buffer2.starts_with('/'))
      {
        if (! addQueryFromString(buffer2, typeFlag)) {
    	  cerr << ": cannot add query from string '" << buffer2 << "' from file '"<< s << "'\n";
    	  return 0;
        } 
      }
      else  // only a file/sub path, use the parent query files path here
      {
        IWString buffer3(s);
        buffer3.truncate_at_last('/');
        buffer3 += '/';
        buffer3 += buffer2;
        if (! addQueryFromString(buffer3, typeFlag))
        {
    	  cerr << ": cannot add query from string '" << buffer3 << "' from file '"<< s << "'\n";
    	  return 0;
        }
      }
    }
    input.do_close();
  }
  else // either "q:" or no leading tag
  {
    if (s.starts_with("q:"))
      s.remove_leading_chars(2);
    if (! addQueryFromString(s, typeFlag))
    {
      cerr << ": cannot add query from string '" << s << "'\n";
      return 0;
    }
  }
 
  return 1;
}


class Reaction_Query
{
  private:
    vector<Reaction_Query_Component *> reactants;
    vector<Reaction_Query_Component *> agents;
    vector<Reaction_Query_Component *> products;

//  private functions

		int addQueryPartFromString(const_IWSubstring s,
                                            vector<Reaction_Query_Component *> & queryPart,
                                            RxnQueryType typeFlag);
                                            
		int matchesOnePart(vector<Molecule *> &m, vector<Reaction_Query_Component *> &queries);
		int matchesOneQuery(vector<Molecule *> &m
												, vector<int> &molDoneFlags
												, vector<Reaction_Query_Component *>&queries
												, int queryIndex);
		int addQueryPartOldFormat(const_IWSubstring &s, vector<Reaction_Query_Component *> &queryPart);
												
  public:
    Reaction_Query()
    {
    }
    ~Reaction_Query()
    {
    	for (vector<Reaction_Query_Component *>::iterator thisIter = reactants.begin() ;  thisIter != reactants.end() ; ++thisIter)
    		delete *thisIter;
    	for (vector<Reaction_Query_Component *>::iterator thisIter = agents.begin() ;  thisIter != agents.end() ; ++thisIter)
    		delete *thisIter;
    	for (vector<Reaction_Query_Component *>::iterator thisIter = products.begin() ;  thisIter != products.end() ; ++thisIter)
    		delete *thisIter;
    }
    
		int addQuery(const const_IWSubstring & buffer, RxnQueryType typeFlag);
		
		int addQueryOldFormat(const const_IWSubstring & buffer);
	
  	int matches(Reaction_To_Search &m);
  	 	
};


// recursive routine to check for a match to a query

int Reaction_Query::matchesOneQuery(vector<Molecule *> &m
																	, vector<int> &molDoneFlags
																	, vector<Reaction_Query_Component *>&queries
																	, int queryIndex)
{
	 	
	for (int i=0; i < m.size(); ++i)
	{
		if (molDoneFlags[i])
			continue;	// already used - not available
			 	
  	if (!(queries[queryIndex])->matches(*m[i]))
  		continue;  // this query does not match this target
  		
    // if we have matches all queries, it is a match!
    
    if (queryIndex + 1 == queries.size())
    	return 1;
    	
    // see if matches can be found to the subsquent queries
    
	  molDoneFlags[i] = 1;  // in use in this incarnation
	  int recursiveMatch = matchesOneQuery(m, molDoneFlags, queries, queryIndex+1); // try the next query
	  molDoneFlags[i] = 0;
	  
	  if (recursiveMatch)
	  	return 1;
  }

  // if we get here, no choice of molecules found a match to the query and allowed matches to
  //	all subsequent queries
  
  return 0;
}

int Reaction_Query::matchesOnePart(vector<Molecule *> &m, vector<Reaction_Query_Component *> &queries)
{
	// first  check that the query does not have more components that the molecule does.  If it does, 
	//  it cannot hit!
	
	if (m.size() < queries.size())
		return false;
		
		
  if (queries.size() == 0)
  	return true;  // there are no queries for this part, so it does not fail to match
		
	// each query in the part must match one of the component molecules in the Molecule part
	// each molecule component can only be used ONCE
	
	vector<int> molDoneFlags(m.size(),0);   // flags for which mol components are already used
	
	// check the first query - the rest are checked recursively

  return matchesOneQuery(m, molDoneFlags, queries, 0);
}

int Reaction_Query::addQuery(const const_IWSubstring & buffer, RxnQueryType typeFlag)
{
 	// the buffer might contain a name for the query - if so discard it
	
	const_IWSubstring s,n;
	
	if (!buffer.split(s, ' ', n))
    s = buffer;
  else if (0 == s.length())
  {
    cerr << "Reaction_Query::addQuery:invalid input '" << buffer << "'\n";
    return 0;
  }

  if (2 != s.ccount('>'))
  {
    cerr << "Reaction_Query::addQuery:reaction must have three components separated by >, " << buffer << " invalid\n";
    return 0;
  }

  const_IWSubstring r, a, p;

  int i = 0;
  s.nextword_single_delimiter(r, i, '>');
  s.nextword_single_delimiter(a, i, '>');
  s.nextword_single_delimiter(p, i, '>');

  if (! addQueryPartFromString(r, reactants, typeFlag) || 
  	  ! addQueryPartFromString(a, agents, typeFlag) || 
  	  ! addQueryPartFromString(p, products, typeFlag))
  {
    cerr << "Reaction_Query::addQuery:invalid input '" << buffer << "'\n";
    return 0;
  }

  return 1;
}

int Reaction_Query::addQueryOldFormat(const const_IWSubstring & s)
{
	// this is the format where the >s and +s are on the command line
	
	//make sure there are two > in the string

  if (2 != s.ccount('>'))
  {
    cerr << "Reaction_Query::addQueryOldFormat:parameter must have three components separated by >, "
     << s << " is invalid\n";
    return 0;
  }

  const_IWSubstring r, a, p;

  int i = 0;
  s.nextword_single_delimiter(r, i, '>');
  s.nextword_single_delimiter(a, i, '>');
  s.nextword_single_delimiter(p, i, '>');

  if (! addQueryPartOldFormat(r, reactants) || 
  	  ! addQueryPartOldFormat(a, agents) || 
  	  ! addQueryPartOldFormat(p, products))
  {
    cerr << "Reaction_Query::addQueryOldFormat:invalid input '" << s << "'\n";
    return 0;
  }

  return 1;
}


int Reaction_Query::addQueryPartFromString(const_IWSubstring s,
                                            vector<Reaction_Query_Component *> & queryPart,
                                            RxnQueryType typeFlag)
{
	queryPart.clear();
	
  if (0 == s.length())    // nothing to do
    return 1;

  resizable_array<int> pos;

  identify_plus_positions(s, pos);
  pos.add(s.length() - 1);

  const int n = pos.number_elements();

  for (int i = 0; i < n; ++i)
  {
    int cstart;          // bunch of code lifted from build_molecules in rxnfile2.cc
    if (0 == i)
      cstart = 0;
    else
      cstart = pos[i-1] + 1;

    int cstop;

    if (i == n - 1)
      cstop = s.length() - 1;
    else
      cstop = pos[i] - 1;

    const_IWSubstring componentStr;

    s.from_to(cstart, cstop, componentStr);

    Reaction_Query_Component * q = new Reaction_Query_Component();
 
    if (! q->addQueryFromString(componentStr,typeFlag))
    {
      cerr << "Reaction_Query::addQueryPartFromString:invalid query string  '" << s << "'\n";
      delete q;
      return 0;    	
    }
	
    queryPart.push_back(q);
  }

  return n;
}

int Reaction_Query::addQueryPartOldFormat(const_IWSubstring &s, vector<Reaction_Query_Component *> &queryPart)
{
	queryPart.clear();
	
  if (0 == s.length())    // nothing to do
    return 1;

  resizable_array<int> pos;

  identify_plus_positions(s, pos);
  pos.add(s.length() - 1);

  const int n = pos.number_elements();

  for (int i = 0; i < n; ++i)
  {
    int cstart;          // bunch of code lifted from build_molecules in rxnfile2.cc
    if (0 == i)
      cstart = 0;
    else
      cstart = pos[i-1] + 1;

    int cstop;

    if (i == n - 1)
      cstop = s.length() - 1;
    else
      cstop = pos[i] - 1;

    const_IWSubstring componentStr;

    s.from_to(cstart, cstop, componentStr);

    Reaction_Query_Component * q = new Reaction_Query_Component();
 
    if (! q->addQueryOldFormat(componentStr))
    {
      cerr << "Reaction_Query::addQueryPartOldFormat:invalid query string '" << s << "'\n";
      delete q;
      return 0;    	
    }
	
    queryPart.push_back(q);
  }

  return n;
}



class RXN_Substructure_Search
{
  private:
  	// these are the reaction queries, that can have one or more quries for reactants, agents or products
  	// these come from reaction smiles or smarts
  	// example: c1ccccc1Cl<<   this is a smiles query for a reaction that has a reactant that is a Cl-benzene derivitive
  	
    vector<Reaction_Query *>  queries;
    vector<IWString> filesToSearch;
    
		ofstream streamForMatches;
		ofstream streamForNonMatches;
		Chemical_Standardisation chemicalStandardization;
		bool reduceToLargestFragment;
		int reactionsRead;
		int reactionsMatching;
  	int verbose;
    
	  		   	
  public:

    RXN_Substructure_Search();
    ~RXN_Substructure_Search();

    int addQueryFromSmiles(const const_IWSubstring & buffer);
    int addQueryFromSmarts(const const_IWSubstring & buffer);
    int addQueryFromRxnQry(const const_IWSubstring & buffer);
    
    int addQueriesFromSmilesFile(const const_IWSubstring & fname);
    int addQueriesFromSmartsFile(const const_IWSubstring & fname);
    int addQueriesFromRxnQryFile(const const_IWSubstring & fname);
    int addQueriesOldFormat(const const_IWSubstring & buffer);
    int isQueryValid() const
    {
      return queries.size() > 0;
    }

    void addSmilesFile (const char *fname)
    {
      IWString fnameToAdd(fname);
      filesToSearch.push_back(fnameToAdd);
    }
		    
    int setStreamForMatches(IWString filename)
    {
      if (streamForMatches.is_open())
        streamForMatches.close();
  
      streamForMatches.open(filename);
      return streamForMatches.is_open();
    }
    
    int setStreamForNonMatches(IWString filename)
    {
      if (streamForNonMatches.is_open())
        streamForNonMatches.close();
    
      streamForNonMatches.open(filename);
      return streamForNonMatches.is_open();
    }
    
    int setChemicalStandardisation(Command_Line &cl, char thisOption)
    {
      return chemicalStandardization.construct_from_command_line(cl, verbose > 1,thisOption);
    }
    
    void setReduceToLargestFragment(bool valueToSet)
    {
    	reduceToLargestFragment = valueToSet;
    }
    
    void setVerbose(int valueToSet)
    {
    	verbose = valueToSet;
    }
    
    int  getReactionsRead()
    {
    	return reactionsRead;
    }
    
    int  getReactionsMatching()
    {
    	return reactionsMatching;
    } 
    
    int search();
      
};

RXN_Substructure_Search::RXN_Substructure_Search():
		reduceToLargestFragment(false)
		,reactionsRead(0)
		,reactionsMatching(0)
  	,verbose(0)
{
}

RXN_Substructure_Search::~RXN_Substructure_Search()
{
   	for (vector<Reaction_Query *>::iterator thisIter = queries.begin() ;  
   														thisIter != queries.end() ; ++thisIter)
  		delete *thisIter;
  		

		if (streamForMatches.is_open())
			streamForMatches.close();

  	
  	if (streamForNonMatches.is_open())
  		streamForNonMatches.close();

}

int
RXN_Substructure_Search::addQueryFromSmiles(const const_IWSubstring & buffer)
{

	Reaction_Query *newOne = new Reaction_Query();
	if (!newOne->addQuery(buffer,RxnQueryType::Smiles))
  {
  	delete newOne;
  	return 0;
  }
   
  queries.push_back(newOne);
  return 1;
}

int
RXN_Substructure_Search::addQueryFromRxnQry(const const_IWSubstring & buffer)
{

	Reaction_Query *newOne = new Reaction_Query();
	if (!newOne->addQuery(buffer,RxnQueryType::QryFile))
  {
  	delete newOne;
  	return 0;
  }
  
  queries.push_back(newOne);
  return 1;
}

int RXN_Substructure_Search::addQueryFromSmarts(const const_IWSubstring & buffer)
{

	Reaction_Query *newOne = new Reaction_Query();
	if (!newOne->addQuery(buffer,RxnQueryType::Smarts))
  {
  	delete newOne;
  	return 0;
  }
  
  queries.push_back(newOne);
  return 1;
}

int RXN_Substructure_Search::addQueriesFromSmilesFile(const const_IWSubstring & fname)
{

  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << ": cannot open '" << fname << "'\n";
    return 0;
  }
                          
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {  
  	 buffer.remove_line_terminators();

     if (!addQueryFromSmiles(buffer))
     {
  		input.do_close();
     	return 0;
     }  	
  }
  input.do_close();
  return 1;
}
    
int RXN_Substructure_Search::addQueriesFromSmartsFile(const const_IWSubstring & fname)
{

  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << ": cannot open '" << fname << "'\n";
    return 0;
  }
                          
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {   
    buffer.remove_line_terminators();

    if (!addQueryFromSmarts(buffer))
    {
      input.do_close();
      return 0;
    }
  }
  input.do_close();
  return 1;
}

int RXN_Substructure_Search::addQueriesFromRxnQryFile(const const_IWSubstring & fname)
{

  iwstring_data_source input(fname);

  if (! input.good())
  {
    cerr << ": cannot open '" << fname << "'\n";
    return 0;
  }
                          
  const_IWSubstring buffer;

  while (input.next_record(buffer))
  {
    buffer.remove_line_terminators();
    if (buffer.starts_with('/')) // OK - we have full path(s)
    {
      if (!addQueryFromRxnQry(buffer))
      {
        input.do_close();
        return 0;
      }
    } 
    else  // only a file/sub path, need to token string and put it back
          // together, for each string token prepend the parent query files
          // path here
    {
      IWString buffer2(fname);
      buffer2.truncate_at_last('/');
      buffer2 += '/';

      IWString query;
      std::vector<std::string> tokens;
      std::string sep(">+");
      buffer.split(tokens, sep);
      for (auto tok : tokens)
      {
        if (tok.front() == '>' || tok.front() == ('+'))
        {
          query += tok;
        }
        else
        {
          query += buffer2;
          query += tok;
        }
      }
 
      if (!addQueryFromRxnQry(query))
      {
        input.do_close();
        return 0;
      }
      buffer.make_empty();  // Seems necessary here
    }
  }
  input.do_close();
  return 1;
}
 
int RXN_Substructure_Search::addQueriesOldFormat(const const_IWSubstring & buffer)
{

	Reaction_Query *newOne = new Reaction_Query();
	if (! newOne->addQueryOldFormat(buffer))
		return 0;
   
  queries.push_back(newOne);
  return 1;
}

int RXN_Substructure_Search::search()
{
  
  // loop over all of the molecule files
  
  reactionsRead = 0;
  for (vector<IWString>::iterator filenameIter = filesToSearch.begin() ; 
  							filenameIter != filesToSearch.end(); ++filenameIter)
  {
  	iwstring_data_source input(*filenameIter);

	  if (! input.good())
	  {
	    cerr << "cannot open '" << *filenameIter << "'\n";
	    return 0;
	  }
	                          
	  const_IWSubstring reactionToSearchStr;

	  while (input.next_record(reactionToSearchStr))
	  { 
	  	reactionToSearchStr.remove_line_terminators();
 
	  	reactionsRead++;

	  	Reaction_To_Search reactionToSearch(chemicalStandardization,reduceToLargestFragment);
	  	if (! reactionToSearch.build(reactionToSearchStr))
	    {
	      cerr << "Cannot read reaction " << reactionToSearchStr << "\n";
	      return 0;
	    }
	    
	    bool thisReactionIsAHit = false;
		  for (vector<Reaction_Query *>::iterator queryIter = queries.begin() ; 
		  							queryIter != queries.end(); ++queryIter)
		  {
	  		if ((*queryIter)->matches(reactionToSearch))
  			{
  				reactionsMatching++;
  				thisReactionIsAHit = true;
				 
		      if (streamForMatches.is_open() )
		      {
		        streamForMatches << reactionToSearchStr << endl;
		      }

  				break;
  			}
		  }

			if (!thisReactionIsAHit)
			{
	      if (streamForNonMatches.is_open())
	      {
	        streamForNonMatches << reactionToSearchStr << endl;
				}
      }		
    }
  	input.do_close();
  }
  
  if (verbose || (!streamForMatches.is_open() && !streamForNonMatches.is_open()))
  {
    cerr << "Read " << reactionsRead << " reactions, " << reactionsMatching << " matched\n";
  }

 	return 1;
}


int Reaction_Query::matches(Reaction_To_Search &m)
{
  // for each part (reactants, agents, products) check the corresponding parts of the 
  //Reaction to search.  Each part in the query must match something for a match to occur
  
  if (!matchesOnePart(m.getReactants(), this->reactants) ||
			!matchesOnePart(m.getAgents(), this->agents) ||
			!matchesOnePart(m.getProducts(), this->products))
	  return 0;
	return 1;
}

void usage ()
{
  cerr << __FILE__ << " compiled " << __DATE__ << " " << __TIME__ << endl;
  cerr << "Substructure Searching over reactions\n";
  cerr << "  -s <q>        Reaction q where q can be\n";
  cerr << "                reaction SMARTS query of the format \"rs1+rs1+rs3>as1+as2>ps1+ps2+ps3\"\n";
  cerr << "                s:rxnSmiles reaction Smiles query of the format\n";
  cerr << "                                 \"rs1+rs1+rs3>as1+as2>ps1+ps2+ps3\"\n";
  cerr << "                f:RxnQry    Reaction Qry string of the format:\n";
  cerr << "                            r1.qry+r2.qry+r3.qry>a1.qry+a2.qry>p1.qry+P2.qry+P3.qry\n";
  cerr << "                            where the *.qry are reactant, agent and product query files\n";
  cerr << "                S:fname     reaction SMARTS file\n";
  cerr << "                M:fname     reaction SMILES file\n";
  cerr << "                F:fname     reaction QRY file\n";;
  cerr << "  -q <q>        queries in the form: \"qry+qry+qry>qry>qry+qry\"\n";
  cerr << "                where qry can be\n";
  cerr << "                smarts      a smarts query\n";
  cerr << "                S:fname     file of smarts\n";
  cerr << "                Q:fname     file of query file names\n";
  cerr << "                q:fname     single query file\n";
  cerr << "  -m <fname>    write reactions that        match to <fname>\n";
  cerr << "  -n <fname>    write reactions that do not match to <fname>\n";
  cerr << "  -l            reduce to largest fragment\n";
  cerr << "  -g ...        chemical standardisation options (enter -g help for details\n";
  cerr << "  -E ...        standard element specificationsenter -E help for details\n";
  cerr << "  -A ...        standard aromaticity specificationsenter -A help for details\n";
  cerr << "  -v            verbose output\n";
}

int
main (int argc, char ** argv)
{
  RXN_Substructure_Search rxnSubstructureSearch;
  
   Command_Line cl(argc, argv, "vA:E:g:lq:n:m:s:");

  if (cl.unrecognised_options_encountered())
  {
    cerr << "Unrecognised options encountered\n";
    usage();
    exit(1);
  }

  int verbose = cl.option_count('v');
  rxnSubstructureSearch.setVerbose(verbose);

  if (cl.option_present('A'))
  {
    if (! process_standard_aromaticity_options(cl, verbose, 'A'))
    {
      cerr << "Cannot initialise aromaticity specifications\n";
	    usage();
	    exit(1);
    }
  }
  else 
    set_global_aromaticity_type(Daylight);

  if (cl.option_present('E'))
  {
    if (! process_elements(cl, verbose, 'E'))
    {
      cerr << "Cannot process -E option\n";
 			usage();
    	exit(1);   
    }
  }
 
  if (cl.option_present('g'))
  {
    if (! rxnSubstructureSearch.setChemicalStandardisation(cl,  'g'))
    {
      cerr << "Cannot process chemical standardisation options (-g)\n";
 			usage();
    	exit(1);   
    }
  }

  if (cl.option_present('l'))
  {
    rxnSubstructureSearch.setReduceToLargestFragment(true);

    if (verbose)
      cerr << "Will reduce to largest fragment\n";
  }

  if (! cl.option_present('q') &&	! cl.option_present('s') )
  {
    	cerr << "Must specify query via the -q  or -s options\n";
 			usage();
    	exit(1);   
  }
  
  if (cl.option_present('q') &&	cl.option_present('s')  )
  {
    cerr << "The options -q and -s cannot both be specified\n";
 			usage();
    	exit(1);   
  } 


  if (cl.option_present('s'))
  {
    const_IWSubstring q = cl.string_value('s');

    if (q.starts_with("s:"))
    {
      q.remove_leading_chars(2); 
	    if (!rxnSubstructureSearch.addQueryFromSmiles(q))
	    {
		   	cerr << "-q s: must be followed by a valid SMILES.  Found: " << q << "\n";
	    	exit(1); 
	    }
    }   
    else if (q.starts_with("f:"))
    {
      q.remove_leading_chars(2); 
	    if (!rxnSubstructureSearch.addQueryFromRxnQry(q))
	    {
		   	cerr << "-q f: must be followed by a valid SMILES.  Found: " << q << "\n";
	    	exit(1); 
	    }
    }    
    else if (q.starts_with("M:"))
    {
      q.remove_leading_chars(2);
	    if (!rxnSubstructureSearch.addQueriesFromSmilesFile(q))
	    {
		   	cerr << "-q M: must be followed by a valid Reactions SMILES.  Found: " << q << "\n";
	    	exit(1); 
	    }
    }
    else if (q.starts_with("S:"))
    {
      q.remove_leading_chars(2);
	    if (!rxnSubstructureSearch.addQueriesFromSmartsFile(q))
	    {
		   	cerr << "-q M: must be followed by a valid Reactions Smarts.  Found: " << q << "\n";
	    	exit(1); 
	    }
    }
    else if (q.starts_with("F:"))
    {
    	q.remove_leading_chars(2);
	    if (!rxnSubstructureSearch.addQueriesFromRxnQryFile(q))
	    {
		   	cerr << "-q M: must be followed by a valid Reactions Qry file.  Found: " << q << "\n";
	    	exit(1); 
	    }
    }
    else
    {
	    if (! rxnSubstructureSearch.addQueryFromSmarts(q))
	    {
	      cerr << "Cannot build reaction query for the reaction SMARTS '" << q << "'\n";
	 			usage();
	    	exit(1);   
	    }
    }
  }
 

  if (cl.option_present('q'))
  {
    const_IWSubstring q = cl.string_value('q');

		// make one query
		
    if (!rxnSubstructureSearch.addQueriesOldFormat(q))
    {
	   	cerr << "-q : must be followed a valid query string  Found: " << q << "\n";
    	usage();
    	exit(1); 
    }		
  
  }	
  // see how many input files are specified
   
  if (0 == cl.number_elements())
  {
      cerr << "Insufficient arguments\n";
 			usage();
    	exit(1); 
  }


  if (cl.option_present('m'))
  {
    IWString tmp = cl.string_value('m');

    if (! tmp.ends_with(".rxnsmi"))
      tmp << ".rxnsmi";

    if (! rxnSubstructureSearch.setStreamForMatches(tmp))
    {
      cerr << "Cannot open stream for matches '" << tmp << "'\n";
 			usage();
    	exit(1);    
    }

    if (verbose)
      cerr << "Reactions matching written to '" << tmp << "'\n";

  }

  if (cl.option_present('n'))
  {
    IWString tmp = cl.string_value('n');

    if (! tmp.ends_with(".rxnsmi"))
      tmp << ".rxnsmi";

  	if (! rxnSubstructureSearch.setStreamForNonMatches(tmp))
    {
      cerr << "Cannot open stream for non matches '" << tmp << "'\n";
      exit(1);  
    }

    if (verbose)
      cerr << "Reactions not matching written to '" << tmp << "'\n";

  }
  
  // check to make sure there is at least one query is some part
  
  if (!rxnSubstructureSearch.isQueryValid())
  {
  	cerr << "the rxn query must contain some valid query query in the reactants, products or agents\n";
  	exit(1);  
  }
  	
  
  //Here is the actual searching part

  for (int i = 0; i < cl.number_elements(); i++)
  {
  	
  	rxnSubstructureSearch.addSmilesFile(cl[i]);
  }
  
  if (! rxnSubstructureSearch.search())
		exit(1);
  exit(0);
}
