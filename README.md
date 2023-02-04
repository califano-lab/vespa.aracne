# vespa.aracne - ARACNe module for the VESPA algorithm
Gene regulatory or signaling network reverse engineering through inference of mutual information by hybrid partitioning.

## Overview
vespa.aracne (Algorithm for the Reconstruction of Accurate Cellular Networks with Adaptive Partitioning) is a signaling network specific extension of the original ARANCe and ARACNE-AP algorithms.

Lachmann A, Giorgi FM, Lopez G, Califano A. *ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information.* **Bioinformatics.** 2016 Jul 15;32(14):2233-5. doi: [10.1093/bioinformatics/btw216](https://dx.doi.org/10.1093/bioinformatics/btw216). Epub 2016 Apr 23.

Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G, Dalla Favera R, Califano A. *ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context.* **BMC Bioinformatics.** 2006 Mar 20;7 Suppl 1:S7. doi: [10.1186/1471-2105-7-S1-S7](https://dx.doi.org/10.1186/1471-2105-7-S1-S7)

## Installation & Building
``vespa.aracne`` requires JDK > 1.8 and ANT. Use the following command in the repository root directory to build the ``jar`` and documentation:

```
ant main
```

The jar will be placed in ``dist/aracne.jar``. The documentation can be found in ``docs/index.html``.

## Usage
For further information on usage, refer to the [`vespa.tutorial`](https://github.com/califano-lab/vespa.tutorial) documentation.

### Input files needed to run vespa.aracne
See below for file format specification (or download the test files from our repository)
1.	Gene expression or phosphoproteomic matrix.
2.	List of regulators (e.g. Kinases or Transcription Factors)

### Steps required to run vespa.aracne
1.	Calculate a threshold for Mutual Information
2.	Run vespa.aracne on bootstraps of the input matrix
3.	Consolidate, i.e. combine the bootstraps into a final network file

### Optional ways to run vespa.aracne
1.	Removing signal transduction (stDPI) or standard DPI (Data Process Inequality) will preserve every edge that passes the Mutual Information threshold.
2.	vespa.aracne by default operates on a bootstrapped version of the input matrix. It is possible to turn this feature off.
3.	During the consolidation step, a Bonferroni correction is applied to the p-values obtained from the Poisson distribution used to determine the significance of each edge according to their appearance in different bootstraps. It is possible to turn this correction off, with the result of having a slightly bigger output network.

### Output of vespa.aracne
Apart from the individual bootstraps, the consolidation step of ARACNe-AP will produce a file called network.txt in the output folder provided by the user. This file shows every significant interaction in four columns
1.	The regulator.
2.	The target.
3.	The MI (Mutual Information) of the pair.
4.  Spearman's correlation of the pair.
5.  The prior probability of the pair.
6.	The pvalue of the pair (assessed during the consolidation step by integrating the bootstraps).

## Input file format
### Gene / phosphosite lists
A text file, containing one gene symbol per line, e.g.
```
g165
g196
g257
g367
g401
g1390
```

### Dataset
A text file, tab separated, with genes / phosphosites on rows and samples on columns
```
gene    Sample1   Sample2   Sample3
g1   1.8 5.2 4.1
g2   5.7 8.3 2.0
g3   6.2 3.1 9.2
g4   7.2 9.1 0.6
```

## Parameters
``-c`` Run ARACNe in consolidation mode

``-t`` Run ARACNe in MI threshold estimation mode

``-nd`` Run ARACNe without DPI

``-nb`` Run ARACNe without bootstrapping

``-e`` Expression Matrix (M x N); M=genes, N=samples; Designate missing values with NA

``-o`` Output directory

``-r`` Regulator identifier file (e.g. transcription factors or kinases & phosphatases)

``-a`` Activator identifier file (e.g. kinases)

``-tg`` Target identifier file (e.g. genes) [default: use all genes or proteins]

``-i`` Protein-protein interaction file (e.g. HSM/P physical interaction predictions)

``-f`` Threshold estimation mode: family-wise error-rate [default: 0.05]

``-mi`` Threshold estimation mode: maximum interactions to assess [default: 1000000]

``-ct`` Absolute correlation threshold to trust mode of interaction [default: 0.0]

``-s`` Optional seed for reproducible results [default: random]

``-j`` Number of threads to use [default: 1]

``-p`` P-value threshold for the Poisson test of edge significance [default: 0.05]

``-m`` Method for multiple-testing correction [BH, Bonferroni, none; default: BH]

## Examples
Note: the examples have been written based on the provided test sets: ``test/matrix.txt`` (the gene expression matrix) and ``test/tfs.txt`` (the list of regulators). Also, example 3 (the running of 100 bootstraps) is written using a “for loop” as a useful method to run 100 bootstraps with a controlled seed.

### Example 1: calculate threshold with a fixed seed
```
java -Xmx5G -jar aracne.jar -e test/matrix.txt  -o outputFolder -r test/tfs.txt -t -s 1
```

### Example 2: run vespa.aracne on a single bootstrap
```
java -Xmx5G -jar aracne.jar -e test/matrix.txt  -o outputFolder -r test/tfs.txt -s 1
```

### Example 3: run 100 reproducible bootstraps
#### UNIX loop
```
for i in {1..100}
do
java -Xmx5G -jar aracne.jar -e test/matrix.txt  -o outputFolder -r test/tfs.txt -s $i
done
```
#### Windows loop
```
for /l %i in (1, 1, 100) do java -Xmx5G -jar aracne.jar -e test/matrix.txt  -o outputFolder -r test/tfs.txt -s %i
```

### Example 4: consolidate bootstraps in the output folder
```
java -Xmx5G -jar aracne.jar -o outputFolder -c
```

### Example 5: run a single ARACNE with no bootstrap and no DPI
```
java -Xmx5G -jar aracne.jar -e test/matrix.txt  -o outputFolder -r test/tfs.txt -s 1 -nd -nb
```

### Example 6: consolidate bootstraps without Bonferroni correction
```
java -Xmx5G -jar aracne.jar -o outputFolder -c -m none
```
