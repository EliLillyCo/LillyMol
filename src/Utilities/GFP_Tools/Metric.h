#ifndef METRIC_H
#define METRIC_H

#include "Enrichment.h"
#include <string>
//class ostream;

/**
 * This is an abstract class difining the interface for a derived metric class. The intention here is to provide a set of function object that can operate on the Enrichment class using polymorphism on the metric calculation method.
 */ 
class Metric
{
 public:

  /**
   * Return the scalar value calculated on the ranking list contained by p_enrich.
   */
  virtual double operator()( const Enrichment& p_enrich ) const = 0;
  
  /**
   * A way to return the name of the metric.
   */
  virtual std::string whichMetric( void ) const = 0;
};

/**
 * Calculate the area under the receiver operating characteristic curve (ROC).
 */
class calculateROC : public Metric
{
 public:
  double operator()( const Enrichment& p_enrich ) const;
  std::string whichMetric( void ) const {return "ROC";};
};

/**
 * Calculate the Boltzmann-averaged receiver operating characteristic (BEDROC) metric.
 */
class calculateBEDROC : public Metric
{
 public:
  calculateBEDROC( double p_alpha = 20.0 ){ setAlpha( p_alpha );};
  void setAlpha( double p_alpha );
  double getAlpha( void ) const { return m_alpha; };
  double operator()( const Enrichment& p_enrich ) const;
  std::string whichMetric( void ) const {return "BEDROC";};
  
 protected:
  double m_alpha;
};

/**
 * Calculate the Enrichment Factor (EF) metric.
 */
class calculateEF : public Metric
{
 public:
  calculateEF( double p_fraction ){ setFraction( p_fraction ); };
  void setFraction( double p_fraction );
  double getFraction( void ) const{ return m_fraction;};
  double operator()( const Enrichment& p_enrich ) const;
  std::string whichMetric( void ) const {return "EF";};

 protected:
  double m_fraction;
};

/**
 * Calculate the Robust Initial Enhancement (RIE) metric.
 */
class calculateRIE : public Metric
{
 public:
  calculateRIE( double p_alpha ){ setAlpha( p_alpha );};
  void setAlpha( double p_alpha );
  double getAlpha( void ) const { return m_alpha; };
  double operator()( const Enrichment& p_enrich ) const;
  std::string whichMetric( void ) const {return "RIE";};
  
 protected:
  double m_alpha;
};

/**
 * Calculate the area under the accumulation curve metric (AUAC).
 */
class calculateAUAC : public Metric
{
 public:
  double operator()( const Enrichment& p_enrich ) const;
  std::string whichMetric( void ) const {return "AUAC";};
};

#endif
