# Leave One Out Models

We recently encountered a dataset with fewer than 50 pieces of available
data. It is quite unclear whether or not a useful model could be built
given such paucity of data. But we did want to explore whether or not
there was any chance of success.

We decided to generate all possible leave one out models and then across
those models, see what the average performance looked like. If it was
poor, then this endeavour seems unfruitful. If there is some kind of
encouraging result, then this should be explored.

## Splits
This part is really the main reason for this article. Using `dopattern`
and standard Linux commands, we can easily generate the needed files.

For the smiles file, we need to generate, for each line in the file,
a test set file containing that line, and training set file containing
all the other records.
```
dopattern -o 40 'sed -n %p all.smi > test%.smi
```
Conversely, for each training set
```
dopattern -o 40 'sed %d all.smi > train%.smi
```
generates a training set file with all records except one.

## Building models
Theoretically some kind of fingerprint calibration process could have
been run here, but that is not reported here. In reality, just a couple of
fingerprint types were tested.

```
dopattern.sh -a 1 -o 40 -parallel 8 'svmfp_make -gfp -CATSP12 -RS -EC3:APT -MAP12:CPY1 -IW -w 20 -w -v -gfp -mdir MODEL% -A all.activity  train%.smi > train%.log 2>&1'
```
builds all models, 8 way parallel on the local machine. This is fast.

Next the models must be evaluated.
```
dopattern -parallel 8 -a 1 -o 40 'svmfp_evaluate -mdir MODEL% test%.smi > test%.pred
```
Generates a bunch of predicted value files, each with 1 record. Concatenate them.
```
descriptor_file_cat.sh test*pred > all.pred
```
and get summary stats
```
iwstats.sh -p 2 -E all.activity  all.pred
```
Building any kind of ML model with fewer than 50 data points must always be
met with considerable skepticism. There just isn't enough data to lear much.
But doing this simple leave-one-out exercise can hopefully give some idea of
whether or not further attempts at model building are likely to be fruitful
or not - or is the only solution to get more data, or completely different
features..
