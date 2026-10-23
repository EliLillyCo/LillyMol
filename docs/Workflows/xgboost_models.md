# xgboost models
Many people have observed the excellent training times for xgboost.
The xgboost authors have gone to heroic efforts on the efficiency of that
tool.

Perhaps unsurprisingly, inference is also very fast. This raises the
possibility of having models that might be able to process very
large numbers of molecules in reasonable times.

Note that xgboost does work in parallel, in highly optimised ways, so
all the times reported here are based on whatever parallel processing
xgboost has been able to perform.

While svmfp models have proven to be very accurate, scoring can be
slow. The speed of scoring depends on the number of support vectors
in the model, and the number and type of fingerprints used.

For example, evaluating an svmfp model that contains 13k support
vectors and with the fingerprints
```
-CATSP13 -ECC3 -ecbig -CLOGP
```
takes almost 4:17 minutes to score 50k molecules - this is on hardware
from 2014. From this we infer a rate of about 85 minutes per million
molecules.

As an aside, that fingerprint computation uses BioByte for the logP
fingerprint. Just generating the fingerprints on 50k molecules takes
24 seconds. If instead the -ALOGP fingerprint is used, that time drops to
under 6 seconds. But clearly the model evaluation would remain a significant
bottleneck.

A much larger model, this one with 88k support vectors, takes 1 hour and 57
minutes to score 50k molecules. Clearly such a model cannot be applied to
sets of molecules in the millions.

## Experiment

Can we build a fast-to-score model that is good enough that it will enable scoring
of much larger datasets. To do this, we will take advantage of the very fast
scoring available in xgboost.

In order to do that, we then need descriptors that are fast to compute. That
pretty much rules out all 3D descriptors and some of the 2D descriptors as well.
Start with iwdescr.sh, which can compute as many as 270 mostly interpretable
descriptors quickly.

If asked to compute all available descriptors, it takes 16 seconds to process 50k
molecules. If all optional features are disabled, `iwdescr.sh -O none rand50k.smi` that
takes just 3.3 seconds, generating 104 features.

Some of the less expensive features can be enabled at little extra cost. So running
```
iwdescr.sh -O none -O adjring -O bbr -O charge -O complex -O crowd -O hbond -O shbond -O alogp\
           -O ncon -O pbond -O psa -O ramey -O rcj -O rfuse -O rss -O rssr -O spch -O spcgrp
```
takes about 6 seconds.
Omitting the `charge` features drops this to 4.4 seconds, but those features are often
important. It all depends.

So if descriptors can be computed at the rate of approximately 88 seconds per
million, and scoring is essentially free, that would see a scoring rate much
better than the 85 minutes per million noted above for the same model with svmfp.

The question becomes whether it is possible to build an adequate model with
this fast-to-compute set of features.

## Model Performance
For these experiments we use a single 80/20 split of an ADME dataset. The whole
set contains 32k molecules, so in the split train.smi contains 25800 molecules
and test.smi contains 6460.

The production model for this target generates the following predictions for
the test set
```
 SVMFP:ACTIVITY RMS = 0.478801
 SVMFP:ACTIVITY R2 = 0.716938
 SVMFP:ACTIVITY Q2 = 0.715249
 SVMFP:ACTIVITY Bias = 0.0128124
 SVMFP:ACTIVITY AE50 = 0.267501
 SVMFP:ACTIVITY AE75 = 0.461546
 SVMFP:ACTIVITY AE95 = 0.952051
 SVMFP:ACTIVITY max_error = 3.50656
 SVMFP:ACTIVITY cMSR95_sd = 8.6994
 SVMFP:ACTIVITY cMSR95_mad = 5.37719
 SVMFP:ACTIVITY cMSR90_sd = 6.15692
 SVMFP:ACTIVITY cMSR90_mad = 4.10893
 SVMFP:ACTIVITY B2.input = 0.844993

```
generating those results in 36 seconds.

The xgboost model reports
```
 XGBD:adme RMS = 0.566686
 XGBD:adme R2 = 0.602107
 XGBD:adme Q2 = 0.601533
 XGBD:adme Bias = -0.00240129
 XGBD:adme AE50 = 0.295358
 XGBD:adme AE75 = 0.560729
 XGBD:adme AE95 = 1.19482
 XGBD:adme max_error = 3.68697
 XGBD:adme cMSR95_sd = 13.0512
 XGBD:adme cMSR95_mad = 6.72523
 XGBD:adme cMSR90_sd = 8.66478
 XGBD:adme cMSR90_mad = 4.9626
 XGBD:adme B2.input = 0.792273
```
producing that result in 0.8 seconds for descriptor computation and 2.6 seconds
for evaluation via python. See below for how this can be improved...

Not surprisingly we find that indeed the model built with a very limited
number of features not as good as the SVMFP model built with a carefully selected
combination of sophisticated fingerprints and properties. But it did that an order
of magnitude faster.

Note too that a more accurate SVMFP model, built with `-w 0.01` has many more
support vectors and takes almost 4 minutes to score the 6400 molecules in the
test split, but reporting enhanced performance
```
 SVMFP:adme RMS = 0.466653
 SVMFP:adme R2 = 0.73105
 SVMFP:adme Q2 = 0.729794
 SVMFP:adme Bias = -0.0315132
 SVMFP:adme AE50 = 0.217225
 SVMFP:adme AE75 = 0.425649
 SVMFP:adme AE95 = 0.967693
 SVMFP:adme max_error = 3.99069
 SVMFP:adme cMSR95_sd = 8.16318
 SVMFP:adme cMSR95_mad = 4.20616
 SVMFP:adme cMSR90_sd = 5.82863
 SVMFP:adme cMSR90_mad = 3.33818
 SVMFP:adme B2.input = 0.863363
``` 
Although surprisingly, AE95 is not better. We see that the median error, AE50,
has dropped from 0.267 to 0.217, AE75 drops from 0.461 to 0.425, while the
95'th percentile error is worse, 0.952 vs 0.967. Perhaps there is something
going on here to be better understood.

A slightly more expensive descriptor calculation, using all the -w features and adding the -abr features,
takes about 3 seconds to compute and results in
```
 XGBD_adme RMS = 0.527344
 XGBD_adme R2 = 0.65496
 XGBD_adme Q2 = 0.654939
 XGBD_adme Bias = 0.000660324
 XGBD_adme AE50 = 0.27875
 XGBD_adme AE75 = 0.522334
 XGBD_adme AE95 = 1.10442
 XGBD_adme max_error = 3.66587
 XGBD_adme cMSR95_sd = 10.9326
 XGBD_adme cMSR95_mad = 6.15278
 XGBD_adme cMSR90_sd = 7.46692
 XGBD_adme cMSR90_mad = 4.60556
 XGBD_adme B2.input = 0.819743
```
which does seem like a worthwhile improvement in performance at a moderate increase
in descriptor computation time - and still remaining an order of magnitude faster
than the production SVMFP model. During calibrate, the range of RMS values observed
was between 0.551 (worst) and 0.499 (best) so the descriptor model above (0.527) is quite
competitive, despite the small number and crudeness of the features.

## Summary
This highlights an important trade-off between accuracy and time. We can have
a more accurate model, but at much greater cost. We can save the computational
cost and get a model that is not quite as good.

Here we might think about two users of such models. The retail consumer is
a scientist with mouse in hand, checking the results for (say) 10 molecules.
Important decisions are going to be made based on those 10 results.
The wholesale consumer is the Cheminformatics person who is generating
100M virtual molecules and needs them evaluated. The retail customer is
best served by the most accurate model we know how to produce. But that
model might be infeasible for the wholesale user. They would be better
served by a faster, but less accurate, model.

## Pipelined Descriptor Computation and Model Evaluation
Some descriptor computation scripts have been modified to enable them to participate
in a pipelined descriptor computation - where descriptors are generated by one programme
and then passed to the next tool which adds its features to the output. This enables
reasonable parallel usage of multiple cores, while avoiding all use of temporary files.

In addition, the tool `xgboost_model_evaluate` evaluates xgboost models via the
C interface and is therefore an order of magnitude faster than doing model evaluation
in python. 

The following command
```
descriptor_pipeline.sh -w -abr rand50k.smi | xgboost_model_evaluate.sh -mdir MODEL -
```
scores 50k molecules in about 18 seconds, or a rate of 1 million in 6 minutes, about
10 million molecules per hour.

What would make sense is to enable multiple -mdir options to xgboost_model_evaluate
so that multiple models could be scored.

What is especially interesting about this is that it entails no temporary
files at all, and could therefore run indefinitely. Output could be piped
to something like `dfilefilter` which filters a descriptor file based on
cutoffs.

So something like
```
generator ... | descriptor_pipeline.sh -w -abr - | 
              xgboost_model_evaluate.sh -mdir MODEL1 -mdir MODEL2 -mdir MODEL3 - |
              dfilefilter -e 'm1<1.0 && m2>3.14 && m3<2.4' -
```
should soon be a reality, with speeds of millions per hour and no temporary
files. It is doable today with just one model...
