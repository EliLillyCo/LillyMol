`in/neighbours.tfrecord` is a TFDataRecord file of `nnbr::NearNeighbours` protos,
generated from the adjacent human-readable neighbour file:

```
nn2proto -T in/neighbours.tfrecord in/neighbours.nn
```

It covers multiple neighbours, empty lists, self-neighbours, distinct zero-distance
neighbours, and a missing query ID. Regenerate with the optimized LillyMol tool.
