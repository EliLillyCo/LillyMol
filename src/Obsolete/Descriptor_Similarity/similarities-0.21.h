#include <math.h>
#include "nrutil.h"

/* Binary (presence/absence) data similarity measures */
double **diceS_row(int **data_matrix, int rows, int columns);                 /* Dice similarity */
double **jaccardS_row(int **data_matrix, int rows, int columns);              /* Jaccard similarity */
double **otsukaS_row(int **data_matrix, int rows, int columns);               /* Otsuka similarity */
double **matchingS_row(int **data_matrix, int rows, int columns);             /* Simple matching similarity */
double **simpsonS_row(int **data_matrix, int rows, int columns);              /* Simspon similarity */
double **soerensenS_row(int **data_matrix, int rows, int columns);            /* Soerensen similarity */
double **baroniUrbaniS_row(int **data_matrix, int rows, int columns);         /* Baroni-Urbani similarity */
double **jaccardS2_row(int **data_matrix, int rows, int columns);             /* Jaccard similarity */
double **soerensenS2_row(int **data_matrix, int rows, int columns);           /* Soerensen similarity */
double **simpleMatchingS_row(int **data_matrix, int rows, int columns);       /* Simple matching similarity */
double **simpleMatching2S_row(int **data_matrix, int rows, int columns);      /* Simple matching similarity */
double **dice2S_row(int **data_matrix, int rows, int columns);                /* Dice similarity */
double **yuleQS_row(int **data_matrix, int rows, int columns);                /* Yule's Q similarity */
double **binCosineS_row(int **data_matrix, int rows, int columns);            /* Cosine similarity */
double **binOverlapS_row(int **data_matrix, int rows, int columns);           /* Binary overlap similarity */      
double **binDispersionS_row(int **data_matrix, int rows, int columns);        /* Binary dispersion similarity */
double **equivalenceS_row(int **data_matrix, int rows, int columns);          /* Equivalence similarity */
double **faith83S_row(int **data_matrix, int rows, int columns);              /* Faith 83 similarity */
double **inclusionS_row(int **data_matrix, int rows, int columns);            /* Inclusion similarity */
double **kulczynski1S_row(int **data_matrix, int rows, int columns);          /* Kulczynski 1 similarity */
double **kulczynski2S_row(int **data_matrix, int rows, int columns);          /* Kulczynski 2 similarity */
double **binPearsonCorrS_row(int **data_matrix, int rows, int columns);       /* Binary Pearson correlation similarity */
double **binOchiai1S_row(int **data_matrix, int rows, int columns);           /* Ochiai 1 similarity */
double **binOchiai2S_row(int **data_matrix, int rows, int columns);           /* Ochiai 2 similarity */
double **binRogersTanimotoS_row(int **data_matrix, int rows, int columns);    /* Rogers & Tanimoto similarity */
double **binRusselRaoS_row(int **data_matrix, int rows, int columns);         /* Russel & Rao similarity */
double **binSokalSneath1S_row(int **data_matrix, int rows, int columns);      /* Sokal & Sneath 1 similarity */
double **binSokalSneath2S_row(int **data_matrix, int rows, int columns);      /* Sokal & Sneath 2 similarity */
double **binSokalSneath3S_row(int **data_matrix, int rows, int columns);      /* Sokal & Sneath 3 similarity */
double **binSokalSneath4S_row(int **data_matrix, int rows, int columns);      /* Sokal & Sneath 4 similarity */
double **binSokalSneath5S_row(int **data_matrix, int rows, int columns);      /* Sokal & Sneath 5 similarity */
double **spearmanRhoS_row(int **data_matrix, int rows, int columns);          /* Spearman rho rank correlation coefficient */
double **kendallTauS_row(int **data_matrix, int rows, int columns);           /* Kendall tau */
double **goodmanKruskalGammaS_row(int **data_matrix, int rows, int columns);  /* Goodman-Kruskal gamma */    
double **binCramerPhiS_row(int **data_matrix, int rows, int columns);         /* Cramer phi */
double **binCramerPhiSquaredS_row(int **data_matrix, int rows, int columns);  /* Cramer phi-squared */
double **binOGES_row(int **data_matrix, int rows, int columns);               /* OGE similarity */
double **binCAS_row(int **data_matrix, int rows, int columns);                /* CA similarity */

/* Binary (presence/absence) data dissimilarity measures */
double **binSqEuclidD_row(int **data_matrix, int rows, int columns);          /* binary squared Euclidean distance */
double **binEuclidD_row(int **data_matrix, int rows, int columns);            /* binary Euclidean distance */
double **binPatternDiffD_row(int **data_matrix, int rows, int columns);       /* binary pattern difference */
double **binVarianceD_row(int **data_matrix, int rows, int columns);          /* binary variance */
double **lanceWilliamsD_row(int **data_matrix, int rows, int columns);        /* Lance & Williams distance */
double **binDiceD_row(int **data_matrix, int rows, int columns);              /* Dice distance */
double **binShapeDiffD_row(int **data_matrix, int rows, int columns);         /* binary shape difference */
double **binSizeDiffD_row(int **data_matrix, int rows, int columns);          /* binary size difference */
double **binEuclid2D_row(int **data_matrix, int rows, int columns);           /* binary Euclidean distance */
double **binAverageSquaredD_row(int **data_matrix, int rows, int columns);    /* binary average squared distance */
double **binShannonD_row(int **data_matrix, int rows, int columns);           /* binary Shannon distance */

/* Interval data similarity measures */
double **circleProductS_row(float **data_matrix, int rows, int columns);      /* circle product similarity */
double **robinsonS_row(float **data_matrix, int rows, int columns);           /* Robinson similarity */
double **scaledTaxonomicS_row(float **data_matrix, int rows, int columns);    /* scaled taxonomic similarity */
double **simRatioS_row(float **data_matrix, int rows, int columns);           /* similarity ratio */
double **czekanowskiS_row(float **data_matrix, int rows, int columns);        /* Czekanowski similarity */
double **pearsonCorrS_row(float **data_matrix, int rows, int columns);        /* Pearson correlation coefficient */

/* Interval data dissimilarity measures */
double **euclidD_row(float **data_matrix, int rows, int columns);             /* Euclidean distance */
double **meanEuclidD_row(float **data_matrix, int rows, int columns);         /* mean Euclidean distance */
double **manhattanD_row(float **data_matrix, int rows, int columns);          /* Manhattan distance */
double **drennanD_row(float **data_matrix, int rows, int columns);            /* Drennan's distance */
double **taxonomicD_row(float **data_matrix, int rows, int columns);          /* taxonomic distance */
double **canberraD_row(float **data_matrix, int rows, int columns);           /* Canberra metric */
double **cosineD_row(float **data_matrix, int rows, int columns);             /* cosine theta distance */
double **minkowskiInfiniteD_row(float **data_matrix, int rows, int columns);  /* Minkowski infinite-power distance */
double **brayCurtisD_row(float **data_matrix, int rows, int columns);         /* Bray-Curtis distance */
double **meanCensoredEuclidD_row(float **data_matrix, int rows, int columns); /* mean censored Euclidean distance */
double **minkowskiD_row(int power,float **data_matrix, int rows, int columns);/* Minkowski distance */

/* Similarity matrix transformations */
double **pairwiseRankComparisonT_row(float **orig_sim_matrix, int rows, int columns);     /* pairwise rank comparison transformation */
double **circleProductT_row(float **orig_sim_matrix, int rows, int columns);              /* circle product transformation */


#else
/* Binary (presence/absence) data similarity measures */
extern double **diceS_row(int **data_matrix, int rows, int columns);                 /* Dice similarity */
extern double **jaccardS_row(int **data_matrix, int rows, int columns);              /* Jaccard similarity */
extern double **otsukaS_row(int **data_matrix, int rows, int columns);               /* Otsuka similarity */
extern double **matchingS_row(int **data_matrix, int rows, int columns);             /* Simple matching similarity */
extern double **simpsonS_row(int **data_matrix, int rows, int columns);              /* Simspon similarity */
extern double **soerensenS_row(int **data_matrix, int rows, int columns);            /* Soerensen similarity */
extern double **baroniUrbaniS_row(int **data_matrix, int rows, int columns);         /* Baroni-Urbani similarity */
extern double **jaccardS2_row(int **data_matrix, int rows, int columns);             /* Jaccard similarity */
extern double **soerensenS2_row(int **data_matrix, int rows, int columns);           /* Soerensen similarity */
extern double **simpleMatchingS_row(int **data_matrix, int rows, int columns);       /* Simple matching similarity */
extern double **simpleMatching2S_row(int **data_matrix, int rows, int columns);      /* Simple matching similarity */
extern double **dice2S_row(int **data_matrix, int rows, int columns);                /* Dice similarity */
extern double **yuleQS_row(int **data_matrix, int rows, int columns);                /* Yule's Q similarity */
extern double **binCosineS_row(int **data_matrix, int rows, int columns);            /* Cosine similarity */
extern double **binOverlapS_row(int **data_matrix, int rows, int columns);           /* Binary overlap similarity */      
extern double **binDispersionS_row(int **data_matrix, int rows, int columns);        /* Binary dispersion similarity */
extern double **equivalenceS_row(int **data_matrix, int rows, int columns);          /* Equivalence similarity */
extern double **faith83S_row(int **data_matrix, int rows, int columns);              /* Faith 83 similarity */
extern double **inclusionS_row(int **data_matrix, int rows, int columns);            /* Inclusion similarity */
extern double **kulczynski1S_row(int **data_matrix, int rows, int columns);          /* Kulczynski 1 similarity */
extern double **kulczynski2S_row(int **data_matrix, int rows, int columns);          /* Kulczynski 2 similarity */
extern double **binPearsonCorrS_row(int **data_matrix, int rows, int columns);       /* Binary Pearson correlation similarity */
extern double **binOchiai1S_row(int **data_matrix, int rows, int columns);           /* Ochiai 1 similarity */
extern double **binOchiai2S_row(int **data_matrix, int rows, int columns);           /* Ochiai 2 similarity */
extern double **binRogersTanimotoS_row(int **data_matrix, int rows, int columns);    /* Rogers & Tanimoto similarity */
extern double **binRusselRaoS_row(int **data_matrix, int rows, int columns);         /* Russel & Rao similarity */
extern double **binSokalSneath1S_row(int **data_matrix, int rows, int columns);      /* Sokal & Sneath 1 similarity */
extern double **binSokalSneath2S_row(int **data_matrix, int rows, int columns);      /* Sokal & Sneath 2 similarity */
extern double **binSokalSneath3S_row(int **data_matrix, int rows, int columns);      /* Sokal & Sneath 3 similarity */
extern double **binSokalSneath4S_row(int **data_matrix, int rows, int columns);      /* Sokal & Sneath 4 similarity */
extern double **binSokalSneath5S_row(int **data_matrix, int rows, int columns);      /* Sokal & Sneath 5 similarity */
extern double **spearmanRhoS_row(int **data_matrix, int rows, int columns);          /* Spearman rho rank correlation coefficient */
extern double **kendallTauS_row(int **data_matrix, int rows, int columns);           /* Kendall tau */
extern double **goodmanKruskalGammaS_row(int **data_matrix, int rows, int columns);  /* Goodman-Kruskal gamma */    
extern double **binCramerPhiS_row(int **data_matrix, int rows, int columns);         /* Cramer phi */
extern double **binCramerPhiSquaredS_row(int **data_matrix, int rows, int columns);  /* Cramer phi-squared */
extern double **binOGES_row(int **data_matrix, int rows, int columns);               /* OGE similarity */
extern double **binCAS_row(int **data_matrix, int rows, int columns);                /* CA similarity */

/* Binary (presence/absence) data dissimilarity measures */
extern double **binSqEuclidD_row(int **data_matrix, int rows, int columns);          /* binary squared Euclidean distance */
extern double **binEuclidD_row(int **data_matrix, int rows, int columns);            /* binary Euclidean distance */
extern double **binPatternDiffD_row(int **data_matrix, int rows, int columns);       /* binary pattern difference */
extern double **binVarianceD_row(int **data_matrix, int rows, int columns);          /* binary variance */
extern double **lanceWilliamsD_row(int **data_matrix, int rows, int columns);        /* Lance & Williams distance */
extern double **binDiceD_row(int **data_matrix, int rows, int columns);              /* Dice distance */
extern double **binShapeDiffD_row(int **data_matrix, int rows, int columns);         /* binary shape difference */
extern double **binSizeDiffD_row(int **data_matrix, int rows, int columns);          /* binary size difference */
extern double **binEuclid2D_row(int **data_matrix, int rows, int columns);           /* binary Euclidean distance */
extern double **binAverageSquaredD_row(int **data_matrix, int rows, int columns);    /* binary average squared distance */
extern double **binShannonD_row(int **data_matrix, int rows, int columns);           /* binary Shannon distance */

/* Interval data similarity measures */
extern double **circleProductS_row(float **data_matrix, int rows, int columns);      /* circle product similarity */
extern double **robinsonS_row(float **data_matrix, int rows, int columns);           /* Robinson similarity */
extern double **scaledTaxonomicS_row(float **data_matrix, int rows, int columns);    /* scaled taxonomic similarity */
extern double **simRatioS_row(float **data_matrix, int rows, int columns);           /* similarity ratio */
extern double **czekanowskiS_row(float **data_matrix, int rows, int columns);        /* Czekanowski similarity */
extern double **pearsonCorrS_row(float **data_matrix, int rows, int columns);        /* Pearson correlation coefficient */

/* Interval data dissimilarity measures */
extern double **euclidD_row(float **data_matrix, int rows, int columns);             /* Euclidean distance */
extern double **meanEuclidD_row(float **data_matrix, int rows, int columns);         /* mean Euclidean distance */
extern double **manhattanD_row(float **data_matrix, int rows, int columns);          /* Manhattan distance */
extern double **drennanD_row(float **data_matrix, int rows, int columns);            /* Drennan's distance */
extern double **taxonomicD_row(float **data_matrix, int rows, int columns);          /* taxonomic distance */
extern double **canberraD_row(float **data_matrix, int rows, int columns);           /* Canberra metric */
extern double **cosineD_row(float **data_matrix, int rows, int columns);             /* cosine theta distance */
extern double **minkowskiInfiniteD_row(float **data_matrix, int rows, int columns);  /* Minkowski infinite-power distance */
extern double **brayCurtisD_row(float **data_matrix, int rows, int columns);         /* Bray-Curtis distance */
extern double **meanCensoredEuclidD_row(float **data_matrix, int rows, int columns); /* mean censored Euclidean distance */
extern double **minkowskiD_row(int power,float **data_matrix, int rows, int columns);/* Minkowski distance */

/* Similarity matrix transformations */
extern double **pairwiseRankComparisonT_row(float **orig_sim_matrix, int rows, int columns); /* pairwise rank comparison transformation */
extern double **circleProductT_row(float **orig_sim_matrix, int rows, int columns);          /* circle product transformation */

