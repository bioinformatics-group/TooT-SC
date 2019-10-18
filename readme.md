# TooT-SC.v1.svm.11

This tool predicts the substrate class of a given transporter. The class can belong to one of the following eleven categories 
[1] Nonselective
[2] water
[3] inorganic cation
[4] inorganic anion
[5] organic anion
[6] organooxogyn
[7] amino acid and derivatives
[8] other Organonitrogen compound
[9] nucleotide
[10] Organic heterocyclic
[11] Miscellaneous

 
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



### intermediate_files
Contains the homology details needed to extract the features. Details of the  `Blast hits` for each sequence is found here.

intermediate_files/Compostions: Contains the extracted `MSA_TPAAC` features of test set

### db
Contains the database to be used when performing BLAST.


### src
The scripts needed to use the tool.

## HOW TO USE
 - This tool requires that `BLAST` be pre-installed
 - Usage: `Rscript src/TooT_SC_V1_11.R -query=<input> [-TooTSC=<TooTSCdir>] [-out=<outdir>]`
  - `<input>` is your sequence input file in fasta format
  - `<out>` is the output directory where you want the predicted 	results, formatted as csv
  - `<TooTSCdir>` is the directory where the base TooT-SC files 	are located
 - `MSA_PAAC` features of each sequence in the test set is  found under [intermediate_files/Compositions/MSA_PAAC.csv](intermediate_files/Compositions/MSA_PAAC.csv)

