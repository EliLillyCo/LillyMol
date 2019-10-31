// The format read by the Enrichment class is as follows:
// ----------------------------
// 8 1
// 15 2
// 128 3
// 231 4
// 347 5
// 1090 6
// 1456 7
// 2500 7
// ----------------------------
// 
// The first column corresponds to the rank of each active
// and the second column to the active count. The last row
// contains the total number of compounds ranked and the 
// total number of actives. In the situation where not all
// actives are kept by the ranking method (discarded before
// they are scored), it is possible to have (active 6 discarded):
// ----------------------------
// 8 1
// 15 2
// 128 3
// 231 4
// 347 5
// 2500 7
//-----------------------------
// In fact, the last line is not used as a rank, but to 
// determine n (number of actives) and N (the total number
// of compounds).
//

#include "Enrichment.h"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <iostream>
#include <set>
#include <time.h>
Enrichment::Enrichment( void )
{
  initialize();
}

Enrichment::~Enrichment()
{
}

void Enrichment::initialize( void )
{
  m_nbr_actives = 0;
  m_nbr_compounds = 0;
  m_rank.clear();
}

void Enrichment::setNbrActives( int p_nbr_actives )
{
  assert( p_nbr_actives > 0 );
  m_nbr_actives = p_nbr_actives;
}

void Enrichment::setNbrCompounds( int p_nbr_compounds )
{
  assert( p_nbr_compounds > 0 );
  m_nbr_compounds = p_nbr_compounds;
}

std::ostream& Enrichment::printRanks( std::ostream& p_out ) const
{
  rank_type::const_iterator iter;
  int count = 1;
  for( iter = m_rank.begin(); iter < m_rank.end(); iter++ )
  {
    p_out << *iter << " " << count << std::endl;
    count++;
  }
  p_out << getNbrCompounds() << " " << getNbrActives() << std::endl;
  return p_out;
}

std::ostream& Enrichment::printROCCoordinates( std::ostream & p_out ) const
{
  rank_type::const_iterator rank_iter = m_rank.begin();
  int count_actives = 0;
  int count_decoys = 0;
  int active_rank = m_rank[0];
  p_out << "0.0 0.0" << std::endl;
  for( int k = 1; k <= getNbrCompounds() - 1; k++ )
  {
    if( active_rank > k )
    {
      count_decoys++;
    }
    else if( rank_iter == ( m_rank.end() - 1 ) )
    {
      count_decoys++;
    }
    else
    {
      count_actives++;
      rank_iter++;
      active_rank = *rank_iter;
    }
    double x = static_cast<double>( count_decoys ) / static_cast<double>( getNbrInactives() );
    double y = static_cast<double>( count_actives ) / static_cast<double>( getNbrActives() );
    p_out << x << " " << y << std::endl;
  }
  p_out << "1.0 1.0" << std::endl;
  
  return p_out;
}


void Enrichment::buildFromRandom( const int p_nbr_actives, const int p_nbr_compounds )
{
  initialize();
  assert( p_nbr_actives > 0 );
  assert( p_nbr_compounds > p_nbr_actives );

  setNbrActives( p_nbr_actives );
  setNbrCompounds( p_nbr_compounds );

  std::set<int> bag;
  std::set<int>::iterator bag_iter;

  while( bag.size() < getNbrActives() )
  {
    int rank = static_cast<int>( static_cast<double>( rand() ) / static_cast<double>( RAND_MAX ) * static_cast<double>( p_nbr_compounds - 1 ) + 1.0 );
    bag.insert( rank );
  }
  
  m_rank.reserve( p_nbr_actives + 1 );
  for( bag_iter = bag.begin(); bag_iter != bag.end(); bag_iter++ )
  {
    m_rank.push_back( *bag_iter );
  }
}

void Enrichment::buildFromExponential( const int p_nbr_actives, const int p_nbr_compounds, const double p_exponential )
{
  initialize();
  assert( p_nbr_compounds >= p_nbr_actives );
  assert( p_exponential > 0.0 );
  setNbrActives( p_nbr_actives );
  setNbrCompounds( p_nbr_compounds );

  std::set<int> bag;
  std::set<int>::iterator bag_iter;

  const double alpha = p_exponential;
  const double minus_inv_alpha = -1.0 / alpha;
  const double one_minus_exp_alpha = 1.0 - exp( -alpha );
  while( bag.size() < getNbrActives() )
  {
    double u = ( static_cast<double>( rand() ) / static_cast<double>( RAND_MAX ) );
    double x = minus_inv_alpha * log( 1.0 - u * one_minus_exp_alpha );
    int rank = static_cast<int>( x * static_cast<double>( getNbrCompounds() - 1 ) ) + 1;
    bag.insert( rank );
  }
  
  for( bag_iter = bag.begin(); bag_iter != bag.end(); bag_iter++ )
  {
    m_rank.push_back( *bag_iter );
  }
}

void Enrichment::buildSubset( const double p_fraction, Enrichment& p_enrich ) const
{
  assert( p_fraction > 0.0 );
  int nbr_to_delete = static_cast<int>( p_fraction * static_cast<double>( getNbrCompounds() ) );
  
  // Generate nbr_to_delete different random positions
  std::set<int> bag;
  std::set<int>::iterator i_bag;
  srand( time(NULL) );

  while( bag.size() < nbr_to_delete )
  {
    int rank = ( rand() % getNbrCompounds() ) + 1;
    bag.insert( rank );
  }
  
  rank_type new_rank;

  // Update rank of actives found in m_ramk
  rank_type::const_iterator i_rank;
  i_bag = bag.begin();
  int previous_count = 0;

  // Build offset of rank to delete
  bool active_to_delete = false;
  rank_type result;
  int position = 0;
  for( i_rank = m_rank.begin(); i_rank < m_rank.end(); i_rank++ )
  {
    while( *i_bag <= *i_rank && i_bag != bag.end() )
    {
      if( *i_bag == *i_rank )
      {
        active_to_delete = true;
      }
      previous_count++;
      i_bag++;
    }
    if( ! active_to_delete )
    {
      result.push_back( *i_rank - previous_count );
    }
    position++;
    active_to_delete = false;
  }
  p_enrich.initialize();
  p_enrich.setNbrCompounds( getNbrCompounds() - nbr_to_delete );
  p_enrich.setNbrActives( static_cast<int>( result.size() ) );
  p_enrich.m_rank = result;
}

int Enrichment::buildFromSortedArray(float * activity, int n ,float cutoff)
{
  if(n <= 1) 
    return 0;
  initialize();
  setNbrCompounds( n );
  int actives =0;
  for(int i=0;i<n;i++)
  {
    if(activity[i] >= cutoff)
    {
      actives++;
      m_rank.push_back( i+1 );
    }
  }
  if(actives >= n || actives <=0 )
    return 0;
  setNbrActives( actives );
    return 1;
}

void Enrichment::readFromStream( std::istream& p_istream )
{
  initialize();
  int rank = -1;
  int previous_rank = -1;
  int count = -1;
  int previous_count = -1;
  while( p_istream >> rank >> count )
  {
    assert( rank > 0 );
    assert( count > 0 );
    assert( rank >= previous_rank );
    assert( count >= previous_count );
    setNbrCompounds( rank );
    setNbrActives( count );
    m_rank.push_back( rank );
    previous_rank = rank;
    previous_count = count;
  }
  m_rank.pop_back();
}

double Enrichment::calculateRandVarRIE( const int p_nbr_actives, const int p_nbr_compounds, const double p_exponential ) 
{
  assert( p_nbr_actives > 0 );
  assert( p_nbr_compounds >= p_nbr_actives );
  assert( p_exponential > 0.0 );
  double n = static_cast<double>( p_nbr_actives );
  double N = static_cast<double>( p_nbr_compounds );
  double a = p_exponential;
  double exp_aN = exp( a * N );
  double exp_a  = exp( a );
  double exp_2a = exp_a * exp_a;
  double exp_a_N = exp( a / N );
  double exp_2a_N = exp_a_N * exp_a_N;
  double factor = N/n*(1.0/exp_a_N-1.0)*(1.0/exp_a_N-1.0)/(1.0-exp_a)/(1.0-exp_a);
  double term1 = ( 1.0 - exp(-2.0*a ))/( exp(2.0*a/N)-1.0 );
  double term2 = 2.0*(n-1.0)/(N-1.0)*exp(-2.0*a)*(exp(a/N)-exp(a))*(1.0-exp(a))/(exp(a/N)-1.0)/(exp(a/N)-1.0)/(1.0+exp(a/N));
  double term3 = n/N*(1.0-exp(-a))/(exp(a/N)-1.0)*(1.0-exp(-a))/(exp(a/N)-1.0);
  return (term1 + term2)/term3 - 1.0;
}

double Enrichment::calculateRandVarAUAC( const int p_nbr_actives, const int p_nbr_compounds ) 
{
  assert( p_nbr_actives > 0 );
  assert( p_nbr_compounds >= p_nbr_actives );
  double n = static_cast<double>( p_nbr_actives );
  double N = static_cast<double>( p_nbr_compounds );
  return ( N - n )*( N + 1 ) / 12.0 / n / N / N;
}

double Enrichment::calculateRandVarROC( const int p_nbr_actives, const int p_nbr_compounds ) 
{
  assert( p_nbr_actives > 0 );
  assert( p_nbr_compounds >= p_nbr_actives );
  double ri = static_cast<double>( p_nbr_compounds - p_nbr_actives )/static_cast<double>( p_nbr_compounds );
  return calculateRandVarAUAC( p_nbr_actives, p_nbr_compounds ) / ri / ri;
}

double Enrichment::calculateRandVarEF( const int p_nbr_actives, const int p_nbr_compounds, const double p_fraction )
{
  assert( p_nbr_actives > 0 );
  assert( p_nbr_compounds >= p_nbr_actives );
  assert( p_fraction > 0.0 );
  double w = floor( p_fraction * static_cast<double>( p_nbr_compounds ) );
  double n = static_cast<double>( p_nbr_actives );
  double N = static_cast<double>( p_nbr_compounds );
  return w/n/N/p_fraction/p_fraction * ( 1.0 + (n-1.0)/(N-1.0)*(w-1.0)) - w*w/p_fraction/p_fraction/N/N;
}
