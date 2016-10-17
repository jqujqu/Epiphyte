# Epiphyte

Epiphyte is a software package for estimating phylogenetic tree from epigenetic modification (DNA methylation primarily) profiles of multiple species. It assumes

  - known phylogenetic relationship between species
  - epigenetic modification has binary state

## Programs in Epiphyte

### PROGRAM: *simulation*
simulate methylomes according to a phylogenetic tree

Usage:`simulation [OPTIONS] <parameter file>`

|Option| Long tag   | Type   | Default | Description |
| ---- | :--------- |:-------| :-------| :---------- |
| -c   | -CpG       | string | Null    | methcount file or bed file (required) |
| -s   | -singlesite|boolean | false   | simulate sites independently |
| -d   | -desert    |integer | 1000    | desert size|
| -o   | -output    |string  | stdout    | name of output file (default: stdout)|
| -v   | -verbose   |boolean | false   | print more run info |

To see the list of options, use "-?" or "-help".

Example:
```sh
simulation -o Name -c CpGs.bed params.txt
```
##### Required Input:
 - Option `-c` requires a `bed` format text file containing the location of CpG
 sites in the genome, for which methylation read proportions are to be simulated;
 - Parameters required for simulation. The format of `params.txt` is shown
 below. The first column contains the required key word, the rest of each line
 contains user defined parameters for simulation.


    NEWICKTREE  ((Human:0.2,Chimp:0.2)HC:0.5,(Mouse:0.4,Rat:0.4)MR:0.3,Dog:0.4)HCMRD;
    PI0	0.3
    G-MATRIX  0.9 0.95
    Q-MATRIX  0.4 0.6
    COV-PARAM 20  0.5
    FG-BETA-PARAMS  0.443 4.314
    BG-BETA-PARAMS  2.576 0.491

##### Option:

By default, methylomes are simulated assuming interdependent sites. When the
option `-s` is specified, sites in methylomes are simulated independently,
i.e. the G-MATRIX parameter is ignored.

##### Output
- Simulated methylation levels for individual extant species in `methcount`
format: ``<chromosome> <position> <strand> <simulated state (T/C)>
<methylation level> <coverage>``
```sh    
chr1	19464625	+	T	0.210526	19
```
To create a table containing the methylation information of the simulated
methylomes of all species, you can use the `merge-methcounts` program from the
[methpipe] package.
```sh
merge-methcounts -t Name_Human Name_Chimp Name_Mouse   \
  Name_Rat Name_Dog | awk 'NR==1{OFS="\t";
  print "Human","Chimp", "Mouse", "Rat", "Dog"}
  NR>1{print}' > Name_meth_table
```
- Simulated methylation states of the entire tree: ``<Name>_treestates``. The
header contains the newick tree used for simulation. Subsequent lines are
methylation states of the whole tree (pre-order traversal) for all CpG sites.
```sh
##((Human:0.2,Chimp:0.2)HC:0.5,(Mouse:0.4,Rat:0.4)MR:0.3,Dog:0.4)HCMRD:0;
chr1	19464625	+	TTTTTTTT
chr1	19464645	+	CTTTCCCC
chr1	19464648	+	CCCCCCCC
```

- Simulated methylation states presented in "hypomethylation probability" (
1 for fully unmethylated state, 0 for fully methylated state). The first line
contains Tab-separated extant species names defined in the input newick tree.
Subsequent lines contains site name (`<crhom>:<pos>`), and Tab-separated
hypomethylation probabilities.
```sh
Human	Chimp	Mouse	Rat	Dog
chr1:19464625	1.0	1.0	1.0	1.0	1.0
chr1:19464645	1.0	1.0	0.0	0.0	0.0
chr1:19464648	0.0	0.0	0.0	0.0	0.0
```


We can estimate hypomethylation probabilities from simulated methylation levels
for each extant species.
### PROGRAM: *indep-sites-hypo*
posteriors for hypo meth state assuming independent sites

Usage: ```indep-sites-hypo [OPTIONS] <meth-tab>```


|Option| Long tag   | Type    | Default | Description |
| ---- | :--------- |:--------| :-------| :---------- |
| -o   | -out       |str  | Null    | output file (default: stdout) |
|      | -params-in |str   | Null    | parameters file (no training) |
|      | -params-out|str   | Null    | write estimated parameters to file |
| -v   | -verbose   |bool  | false   | print more run info (default: false) |

To see the list of options, use "-?" or "-help".


If sites were simulated as independent from each other using the `-s` option,
we estimate the hypomethylation probabilities using a two-component
beta-binomial mixture model.
```sh
for i in `echo Chimp Human Mouse Rat Dog`; do
  echo $i
  indep-sites-hypo -o Name_${i}.hypoprob Name_${i} \
    -params-out Name_${i}.hypoprob.params
done
```

If sites were simulated as interdependent, we use the program `hmr_posterior`
from the [methpipe] package, which estimates posterior probabilities of
individual CpG sites being hypomethylated.
```sh
for i in `echo Chimp Human Mouse Rat Dog`; do
echo $i
~/git/methpipe/bin/hmr_posterior -v -o Name_${i}.hypoprob Name_${i} \
  -p Name_${i}.params
done
```

##### merge probabilities into a single file
```sh
paste Name_Human.hypoprob Name_Chimp.hypoprob Name_Mouse.hypoprob \
  Name_Rat.hypoprob Name_Dog.hypoprob |  \
  awk 'BEGIN{OFS="\t"; print "Human","Chimp","Mouse","Rat","Dog"}
       {print  $1":"$2,$5,$11,$17,$23,$29, $35}' > Name.hypoprob_table
```


### PROGRAM: *indep-epi-phylo*

Estimate tree shape and mutation rates assuming site-independency

Usage: ``indep-epi-phylo [OPTIONS] assuming independent sites<newick> <meth-tab>``

|Option| Long tag    | Type| Default | Description |
| ---- | :---------- |:----| :-------| :---------- |  
| -c   | -counts    | bool | false   | meth-table contains read counts |
| -p   | -params    | str  | Null    | use given parameters and skip optimization |
|  -i  | -iteration | int | 10      | maximum number of iteration  
|  -n  | -nodemap   | bool | false   | output MAP states of individual nodes|
|  -o  | -out       | str  | stdout  | output file |
|  -v  | -verbose   | bool | false   | print more run info |

To see the list of options, use "-?" or "-help".


```sh
indep-epi-phylo -c -v -o Name_counts.out init.nwk  Name_meth_table
```

### PROGRAM: *epiphyte*
Estimate phylogeny shape and methylation state transition rates for methylome
evolution

Usage: ``epiphyte [OPTIONS] <newick> <hypoprob-tab>``

|Option| Long tag    | Type| Default | Description |
| ---- | :---------- |:---- | :-------| :---------- |  
|  -m  | -minCpG     |int | 10 | minimum #CpGs per block|
|  -i  | -maxiter    |int  | 10 | maximum iteration|
|  -c  | -complete   |bool | false| input contains complete observations |
|  -v  | -verbose    |bool | false | print more run info |
|  -p  | -params     |str  | Null  | given parameters in file |
|  -P  | -outparams  |str  | Null  | output parameters to file |
|  -o  | -output     |str  | Null  | output file name |
|  -f  | -minfragCpG |int  | 5 | minimum #CpGs per fragment to output   |
|  -s  | -single     |bool | false | also output states by sites (use with -o)|

To see the list of options, use "-?" or "-help".

```sh
epiphyte init.nwk Name_hypoprob_table -i 30 -v -s   \
  -o Name_hme.out -P Name_hme.params 2> Name_RunInfo.txt
```
License
----



[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [methpipe]:<https://github.com/smithlabcode/methpipe>
