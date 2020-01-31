# The Eli Lilly (LillyMol) Matched Molecular Pairs Toolkit

## Contact
[JamesALumley](https://github.com/JamesALumley)

## About
A toolkit of scripts and methods for Matched Molecular Pairs (MMP) analysis based on an in-memory implementation of the Hussain and Rea fragment indexing algorithm (JCIM 2010). 

An MMP is defined as a pair of molecules that differ in only a minor single point of chemical change. MMP Analysis is a method that compares the properties of these two molecules to derive a change or delta value associated with the chemical (fragment) change.

A Transform can be defined as a specified chemical change and its statistical significance across many pairs of molecules (many MMPs). An example transformation is the replacement of one functional group by another (e.g.: Cl => Br) or the addition of a phenyl ring (e.g.: H => Phenyl). A useful Molecular transformation in a specified context is termed a significant transformation if it is seen to systematically decrease or increase a desired property (i.e.: The specified fragment change or Transform is similar in magnitude and direction for many example MMPs. Different metrics can be used to define significant.  Transforms or design rules can be used as an aid to decision-making in the multiple parameter optimisation problem of small molecule drug design.

## Getting Started

### Prerequisites

This python code depends on multiple binary executables from the core LillyMol code base such as dicer for fragmentation.  Please follow the instructions to install and compile LillyMol first as per the root folder README file.  Once LillyMol is installed and working the file pymo.py can be edited to correctly reference all binaries, typically using the environment variable LILLYMOL_HOME to point to the top directory containing the _ _contrib_ _ subdir.  The pymo.py file acts as a file IO based python interface to the main LillyMol executables. The python path should then be set as below:

```
export PYTHONPATH=${LILLYMOL_HOME}/contrib/script/py
cd $LILLYMOL_HOME
cd bin
ln -s Linux-gcc-7.2.1 Linux
```

### Simple MMP / MMS Generation
The following scripts are included in this project:
* **getMMPfromSMI.py** - For an input SMILES file (Molecule SMILES and numeric ID) generate a CSV of Matched Molecular Pairs
* **getMMPStatsfromCSV.py** - For an input CSV file (Molecule SMILES, numeric ID and numerical data column) generate a CSV of Matched Molecular Pairs with associated delta for the given data column(s). Summarised (aggregated) transforms can additionally be requested. A second CSV file will be generated containing these transforms. Options are available for different statistical aggregation methods.   
  * MEAN_DIFF: The simplest form of aggregation would ensure the input bioactivity data is in log form, then calculate the mean of the difference (-A MEAN_DIFF). The [scipy students t-test](http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.t.html) is used to calculate statistical significance.  
  * FOLD_CHANGE: For data that is already in the log scale, such as solubility (LogS  =  log10(Solubility mg/ml)), this method can be used to return the reverse log of the mean (10^x). This represents the average solubility ratio of substituted over the parent molecules for a given matched-pair.  Lower values indicate the transformation is likely to reduce solubility, positive values vice-versa and values close to 1 indicate minimal change expected.
  * DIFFXX:  The diffxx metric is a scaled Index with upper and lower boundaries determined by the t.interval. It is for use with data like % Microsomal Metabolic Turnover (DIFF60) or Est IC50 for CYP Inhibition (DIFF50) and will be described elsewhere.
  * CATEGORICAL: For categorical data, the method assumes two class categorical data with values 0 and 1.  A textual value is output describing the number of pairs that have or have not seen a change in class.
  * Atom Context: The aggregation type switch (-g) allows you to include the atom context in the aggregation and will most likely produce a greater number of summarised pairs.  In some cases transforms will have improved predictive power using atom attachment points.
  * Property Data: can be added via the -p switch (default is no property differences). This is based on the output from the LillyMol iwdescr binary and calculates difference (FRAG_R – FRAG_L) regardless of -A.
  * Quality Metric: The quality estimate was designed for use with the commercial Biobyte clogP algorithm.  As an alternative, the opensource LillyMol Abrahams AlogP can be used. It is used to estimate the quality of a transfrom from the perspective of chemical diversity and uses logP as a surrogate.  Where L.clogp is the clogp of the parent (Left) compound and d.clogp is the delta clogp (Right – Left), we can defined the following:  A good transform has >=15 contributing MMPs, a stdev of L.clogp >= 0.5 and when |d.clogp|>0, stdev of d.clogp>0.  A medium quality transform has >=6 contributing MMPs, a stdev of L.clogp >= 0.3 and when |d.clogp>|0, stdev of d.clogp>0.  Any other transform is classed as poor quality.
* **getMMPbasedMCSSfromSMI.py** - For an input SMILES file (Molecule SMILES and numeric ID) generate a CSV of the MMP that has the smallest change between the two molecules (by atom count). This is roughly equivalent to the Maximum Common Substructure (MCSS) between two molecules but is based on SMILES and is not as expressive as more flexible SMARTS based definition.
* **getMMPEnumeratedNewMols.py** - For an input SMI file (Molecule SMILES and numeric ID) mutate the molecule to create new idea molecules based on a prebuilt or custom set of Transforms.  This code is useful in evolving molecules in an automated design scenario by applying a file of suggested/selected transforms.
* **getMMPSeriesfromCSV.py** - For an input CSV file (Molecule SMILES, numeric ID and numerical data column) generate a CSV of Matched Molecular Series. The file will print each series as a ordered series with one molecule per line with associated data. This is an automated approach to deriving new idea compounds with improved activity data. It is also an approach for SAR transfer.
* **getMMPSeriesSuggestfromCSV.py** - For an input SMI file and a second .series file (or a directory of .series files) search the input SMI for matched series that can be extended by the series in the second file(s). This is an automated approach to deriving new idea compounds with improved activity data and for SAR transfer.

### Implementation Details
The python MMP code is based on the LillyMol dicer executable, a c-code based executable for molecule fragmentation.  As molecule smiles are parsed by dicer, the resulting canonicalised fragments are processed and stored in a dictionary.  The larger fragment is termed the context smiles and the smaller part the fragment.  The dicer switch MAXFF=0.50001 is used and means that the ‘context’ or the bit of the matched pair that does not change, is defined as the larger part.  Larger is >50.001% of the whole molecule (fraction of atoms).  The resulting dictionary is keyed by the context smiles and an associated iterator can return all MMP's from the fully populated dictionary.  Atoms involved in bond split/fragmentation are tagged with an isotopic label e.g.: the smiles [1CH3] which represents a terminal methyl group that has been cut/fragmented from a molecule with the C atom labelled with an isotopic label 1.  Support for atom attachment points is included.  The concept of denormalisation is used on the stored fragment/context smiles strings in order to conserve memory.  A lookup dictionary is used to assign an arbitary id to all smiles strings and the id stored in the mmp dictionary to avoid repetitively storing larger smiles strings.  The result is a fast mmp implementation that is not database coupled and therefore lighterweight than some other opensource implementations. The obvious downside is greater memory usage.

### Examples

```
getMMPfromSMI.py --help

getMMPfromSMI.py –i input.smi –o output.pairs

getMMPStatsfromCSV.py --help

getMMPStatsfromCSV.py –i input.csv –o output.pairs -s SMI_COLUMN -n ID_COLUMN -a DATA_COLUMN -A AGGREGATION_METHOD

```


## Running the tests

Each class has it's own set of unit tests that can be triggered independently e.g.:
```
python mmp_objects.py
python mmp_stats_functions.py
```
Alternatively, integration tests that demonstrate the usage of the code can be run from the LillyMol test suite:
```
cd $LILLYMOL_HOME/tests/getMMPStatsfromCSV/
./run_test
```

## License

Please see the LICENSE file in the root of the repository
