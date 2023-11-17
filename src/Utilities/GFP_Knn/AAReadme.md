# Knn with gfp

This is very old code and should be re-implemented.

It does work and is used in some production systems, but needs to be re-done
from the ground up.

nn_training builds tunes the K and/or distance parameters of a knn model to decide on an
optimal set, and instantiates that model. That model file can then be used by
nn_predictions to make predictions.

The code was never satisfactory, but does seem to work and often builds a
helpful baseline model.
