# TooT-SC.v1.svm.11

This tool predicts the substrate class of a given transporter. The class can belong to one of the following eleven categories 
[1] nonselective inorganic molecule
[2] water
[3] inorganic cation
[4] inorganic anion
[5] organic anion
[6] organooxogyn
[7] amino acid and derivatives
[8] other Organonitrogen compound
[9] nucleotide
[10] organic heterocyclic
[11] miscellaneous.


Details on ChEBI mapping is found on Substrate_ChEBI.csv file.
 
Input: transporter proteins sequences in Fasta format
Output: the substrate class with highest probability, and the probabilities of the other classes.

## Installation

```bash
git clone https://github.com/bioinformatics-group/TooT-SC.git
cd TooT-SC
wget https://tootsuite.encs.concordia.ca/databases/TooT-SC/v1.svm.11_db.tar.gz
wget https://tootsuite.encs.concordia.ca/datasets/TooT-SC/v1.svm.11_datasets.tar.gz
tar -xzf v1.svm.11_db.tar.gz
tar -xzf v1.svm.11_datasets.tar.gz
rm v1*
```
Use `curl` if you don't have `wget`.

## FOLDERS
There are a number of folders that support the running of TranCEP and its outputs.

### dataset/transporter_substrate_class
Contains both the training and the independent testing dataset if fasta format, the identifiers of independent testset is found under substrate_classes_indep.csv.
The models were trained using the training dataset =(all dataset - independent testset)


### Working Directory
When run, the tool will create a `work` directory, common among all the TooT Suite. TooT-SC will create a directory named `TooT-SC` which will hold intermediate files created during the running of these tools.
The top-level working directory contains the homology details needed to extract the features. Details of the `Blast hits` for each sequence are found here, one directory per sequence.

The `TooT-SC` directory will also contain a `Compostions` folder which contains the extracted `MSA_TPAAC` features of the test set

The working directory can be quite large, so you may want to use `-work` to specify a better location for it if you're short on space, and be sure to clean it out from time to time.

### db
The database to be used when performing BLAST.

By default, if you `tar -xzf` the contents of [this](https://tootsuite.encs.concordia.ca/databases/SwissOct18.tar.gz) into a `db` folder adjacent to the source folder, it should be in the default location that this tool checks. You may also specify the location manually via `-db`.

### src
The scripts needed to use the tool.

## HOW TO USE
 - This tool requires that `BLAST` be pre-installed
 - Usage: `src/TooT_SCTool.R -query=<input> [-out=<outdir>] [-db=<path to db>] [-work=<Workdir>] [-TooTSC=<TooTSCdir>]`
  - `<input>` is your sequence input file in fasta format
  - `<out>` is the output directory where you want the predicted 	results, formatted as csv
  - `<db> is the directory where the database is stored`
  - `<Workdir>` is the directory where intermediate files will be stored after each run
  - `<TooTSCdir>` is the directory where the base TooT-SC files 	are located
 - `MSA_PAAC` features of each sequence in the test set are found under the working directory `Compositions/MSA_PAAC.csv`

