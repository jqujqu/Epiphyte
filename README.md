# Epiphyte
Epigenetic models on phylogenetic tree


# Epiphyte

Epiphyte is a software package for estimating phylogenetic tree from epigenetic modification (DNA methylation primarily) profiles of multiple species.  It assumes

  - known phylogenetic relationship between species
  - epigenetic modification has binary state

## Programs in Epiphyte

### simulation

##### Input: 

The `simulation` requires two input files: `CpGs.bed` -- a `bed` format text file containing the location of CpG sites in the genome, for which methylation read proprotions are to be simulated; and `params.txt` -- parameters required for simulation. The format of `params.txt` is shown below. The first column contains the required key word, the rest of each line contains user defined parameters for simulation.

    NEWICKTREE  ((Human:0.2,Chimp:0.2)HC:0.5,(Mouse:0.4,Rat:0.4)MR:0.3,Dog:0.4)HCMRD;
    PI0	0.3
    G-MATRIX  0.9 0.95
    Q-MATRIX  0.4 0.6
    COV-PARAM 20  0.5
    FG-BETA-PARAMS  0.443 4.314
    BG-BETA-PARAMS  2.576 0.491

```sh
simulation -s -o Name -c CpGs.bed params.txt
```
##### Output 
The output files all share the same prefix as specified in the `-o` option. A `methcount` format text file is generated for each leaf species. The `methcount` format is defined in
    chr1	19464625	+	T	0.210526	19

To create a table containing the methylation information of the simulated methylomes of all species, you can use the `merge-methcounts` program from the [methpipe] package.
```sh
merge-methcounts -t test_Human test_Chimp test_Mouse   \
  test_Rat test_Dog | awk 'NR==1{OFS="\t"; 
  print "Human","Chimp", "Mouse", "Rat", "Dog"} 
  NR>1{print}' > test_meth_table
```




License
----



[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [methpipe]:<https://github.com/smithlabcode/methpipe>

