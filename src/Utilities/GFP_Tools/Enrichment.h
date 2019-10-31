#ifndef ENRICHMENT_H
#define ENRICHMENT_H

#include <vector>
#include <iostream>
//class ostream;

/**
 * This class contains the information concerning the rank ordered list, or equivalently
 * the accumulation curve information.
 */ 
class Enrichment
{
 public:

  /**
   * Default constructor.
   */
  Enrichment();

  /**
   * Default destructor.
   */
  ~Enrichment();

  /**
   * Return the number of actives.
   */
  int getNbrActives( void ) const{ return m_nbr_actives; };

  /**
   * Return the total number of compounds (actives + inactives).
   */
  int getNbrCompounds( void ) const { return m_nbr_compounds; };

  /**
   * Return the number of inactives.
   */
  int getNbrInactives( void ) const { return m_nbr_compounds - m_nbr_actives; };

  /**
   * Output the current state of the object following the format outlined hereafter.
   */
  std::ostream& printRanks( std::ostream & p_out ) const;

  /**
   * Prints the ROC curve coordinates for the current state of the object.
   */
  std::ostream& printROCCoordinates( std::ostream & p_out ) const;
  
  /**
   * Build the rank of each actives from a ranked array of the activities
   * @param activity is experimental values
   * @param n is the nbr of compounds in the array
   * @param cutoff is the cutoff of experimental values for active and inactive
   */
  int buildFromSortedArray(float * activity, int n,float cutoff);
  /**
   * Read the rank of each actives according to the format described hereafter
   */
  void readFromStream( std::istream& p_istream );

  /**
   * Generate a uniformly distributed set of p_nbr_actives considering that there are p_nbr_compounds overall. In most random number generator, the maximum generated number is 32000, so p_nbr_compounds should not be larger than 32000 for a correct behavior. The previous state of this object is lost.
   */
  void buildFromRandom( const int p_nbr_actives, const int p_nbr_compounds );

  /**
   * Generate p_nbr_actives positions according to an exponential distribution of parameter p_exponent: PDF(x) = p_exponent * exp( -p_exponent * x ) / ( 1 - exp(-p_exponent) ) where 'x' is the normalized position ( rank / p_nbr_compounds ). The previous state of this object is lost. 
   * @param p_exponent Must be greater than 0.
   */
  void buildFromExponential( const int p_nbr_actives, const int p_nbr_compounds, const double p_exponent );

  /** 
   * This takes the actual object and randomly remove  1 - p_fraction of the actives/inactives  (no distinction) from the actual enrichment list to generate a new subset of size p_fraction * getNbrCompounds() keeping the order actives/inactives not deleted.
   */
  void buildSubset( const double p_fraction, Enrichment& p_enrich ) const;

  
  // These methods calculate the analytical variance of the specified metrics given that the active are uniformly distributed in the list. Each pair of actives (p_nbr_actives) and total number of compound (p_nbr_compounds) lead to a different value. 
  static double calculateRandVarBEDROC( const int p_nbr_actives, const int p_nbr_compounds, const double p_exponential );
  static double calculateRandVarRIE( const int p_nbr_actives, const int p_nbr_compounds, const double p_exponential );
  static double calculateRandVarAUAC( const int p_nbr_actives, const int p_nbr_compounds ); 
  static double calculateRandVarROC( const int p_nbr_actives, const int p_nbr_compounds ); 
  static double calculateRandVarEF( const int p_nbr_actives, const int p_nbr_compounds, const double p_fraction );

  // Function object classes that calculate metrics on Enrichment objects
  friend class calculateROC;
  friend class calculateBEDROC;
  friend class calculateEF;
  friend class calculateRIE;
  friend class calculateAUAC;

  // Container type that contain the rank of each active in the list
  typedef std::vector< int > rank_type;

 protected:

  void setNbrActives( int p_nbr_actives );
  void setNbrCompounds( int p_nbr_compounds );
  
  // Only the number of actives and total number of compounds are kept
  // in memory, the number of inactives is easily calculated on the fly.
  int m_nbr_actives;
  int m_nbr_compounds;

  // Positions from 1 to N of the actives in the rank ordered list.
  rank_type m_rank;
  
  // Reset the state of the object: empty ranks, 0 actives and inactives.
  void initialize( void );
  
 private:

  // Copy constructor not defined and should not be used.
  Enrichment( const Enrichment& p_enrich ){};
  Enrichment& operator=( const Enrichment& p_enrich ){ return *this; };

};

#endif

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
