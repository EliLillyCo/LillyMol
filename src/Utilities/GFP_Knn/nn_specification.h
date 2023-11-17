#ifndef NN_SPECIFICATION_H
#define NN_SPECIFICATION_H

#include "Foundational/iwstring/iwstring.h"

/*
  Common code for parsing the specification of a model
*/

#include "Foundational/iwstring/iwstring.h"

#define MODEL_COMPONENT_SEPARATOR ';'

#define KNN_MODEL_PREFIX "knn_model:"

#define KNN_MODEL_TYPE "mtype:"

#define TR_GFP_DIRECTIVE "TRgfp:"

#define TR_GFP_UNKNOWN "%UNKTRgfp%"

#define KNN_CONTINUOUS "continuous"
#define KNN_CLASSIFICATION "classification"

#define KNN_MODEL_FILE_VERSION_1 "KNN MODEL 1"

#define EXPT_SOURCE "expt:"

extern int
parse_model_specifications (const const_IWSubstring & c,
                            const_IWSubstring & weight_function,
                            int & special_class,
                            double & special_class_probability,
                            double & distance,
                            int & number_neighbours);

#endif

