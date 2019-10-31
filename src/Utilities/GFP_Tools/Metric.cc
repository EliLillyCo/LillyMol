#include "Metric.h"
#include <cassert>
#include <cmath>
#include <iostream>

void calculateRIE::setAlpha( double p_alpha )
{
  assert( p_alpha > 0.0 );
  m_alpha = p_alpha;
}

void calculateBEDROC::setAlpha( double p_alpha )
{
  assert( p_alpha > 0.0 );
  m_alpha = p_alpha;
}

void calculateEF::setFraction( double p_fraction )
{
  assert( p_fraction > 0.0 );
  m_fraction = p_fraction;
}

double calculateROC::operator()( const Enrichment& p_enrich ) const
{
  int sum = 0;
  int count_active = 0;
  int count_inactive = 0;
  
  for( std::vector<int>::const_iterator i_rank = p_enrich.m_rank.begin(); i_rank < p_enrich.m_rank.end(); ++i_rank )
  {
    count_active++;
    if( i_rank == p_enrich.m_rank.begin() )
    {
      count_inactive = *i_rank - 1;
    }
    else
    {
      count_inactive = *i_rank - count_active;
    }
    sum += p_enrich.getNbrInactives() - count_inactive;
  }
  return static_cast<double>( sum ) / static_cast<double>( p_enrich.getNbrInactives() ) / static_cast<double>( p_enrich.getNbrActives() );
};

double calculateBEDROC::operator()( const Enrichment& p_enrich ) const
{
  double sum = 0.0;
  const double N = static_cast<double>( p_enrich.getNbrCompounds() );
  const double n = static_cast<double>( p_enrich.getNbrActives() );
  double a = getAlpha();
  for(std::vector<int>::const_iterator iter = p_enrich.m_rank.begin(); iter < p_enrich.m_rank.end(); ++iter )
  {
    sum += exp( -a * static_cast<double>( *iter ) / N );
  }
  double ra  = n/N;
  double ri  = (N-n)/N;
  double random_sum = ra*exp( -a/N )*(1.0 - exp( -a ) )/ ( 1.0 - exp( -a / N ) );
  double factor = ra*sinh(a/2.0)/(cosh(a/2.0)-cosh(a/2.0-a*ra));
  double cte = 1.0/(1.0-exp(a*ri));
  return sum/random_sum*factor + cte;
}

double calculateEF::operator()( const Enrichment& p_enrich ) const
{
  int count_active = 0;
  int rank_limit = static_cast<int>( floor( getFraction() * static_cast<double>( p_enrich.getNbrCompounds() ) ) );
  for( std::vector<int>::const_iterator i_rank = p_enrich.m_rank.begin(); i_rank < p_enrich.m_rank.end(); ++i_rank )
  {
    if( *i_rank > rank_limit )
    {
      break;
    }
    count_active++;
  }
  return static_cast<double>( count_active ) / static_cast<double>( p_enrich.getNbrActives() ) / getFraction(); 
}

double calculateRIE::operator()( const Enrichment& p_enrich ) const
{
  double sum = 0.0;
  const double N = static_cast<double>( p_enrich.getNbrCompounds() );
  const double n = static_cast<double>( p_enrich.getNbrActives() );
  double a = getAlpha();
  for (std::vector<int>::const_iterator iter = p_enrich.m_rank.begin(); iter < p_enrich.m_rank.end(); ++iter )
  {
    sum += exp( -a * static_cast<double>( *iter ) / N );
  }
  double ra  = n/N;
  double ri  = (N-n)/N;
  double random_sum = ra*exp( -a/N )*(1.0 - exp( -a ) )/ ( 1.0 - exp( -a / N ) );
  return sum/random_sum;
}

double calculateAUAC::operator()( const Enrichment& p_enrich ) const
{
  double sum = 0.0;
  double N = static_cast<double>( p_enrich.getNbrCompounds() );
  double n = static_cast<double>( p_enrich.getNbrActives() );
  for (std::vector<int>::const_iterator iter = p_enrich.m_rank.begin(); iter < p_enrich.m_rank.end(); ++iter)
  {
    sum += *iter;
  }
  return 1.0 - sum / n / N + 0.5/N;// N /n ;
};
