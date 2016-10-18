# Epiphyte

Epiphyte is a software package for estimating phylogenetic tree from epigenetic
modification (DNA methylation primarily) profiles of multiple species. It
assumes

  - known phylogenetic relationship between species
  - epigenetic modification has binary state

## Programs in Epiphyte

### PROGRAM: *simulation*
simulate methylomes according to a phylogenetic tree

Usage:`simulation [OPTIONS] <parameter file>`

|Option| Long tag   | Type   | Default | Description |
| ---- | :--------- |:-------| :-------| :---------- |
| -c   | -CpG       |str | Null    | methcount file or bed file (required) |
| -s   | -singlesite|bool | false   | simulate sites independently |
| -d   | -desert    |int | 1000    | desert size|
| -o   | -output    |str  | stdout    | name of output file (default: stdout)|
| -v   | -verbose   |bool | false   | print more run info |

To see the list of options, use "-?" or "-help".


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


### PROGRAM: *indep-epi-phylo*

Estimate tree shape and mutation rates assuming site-independency

Usage: ``indep-epi-phylo [OPTIONS] assuming independent sites<newick> <meth-tab>``

|Option| Long tag    | Type| Default | Description |
| ---- | :---------- |:----| :-------| :---------- |  
| -c   | -counts    | bool | false   | meth-table contains read counts |
| -p   | -params    | str  | Null    | use given parameters and skip optimization |
|  -i  | -iteration | int | 20      | maximum number of iteration  
|  -n  | -nodemap   | bool | false   | output MAP states of individual nodes|
|  -o  | -out       | str  | stdout  | output file |
|  -v  | -verbose   | bool | false   | print more run info |

To see the list of options, use "-?" or "-help".


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


## Example 1
In this example, we start by simulate a site-independent evolution process for  
a phylogenetic tree containing 5 extant species, then estimate the model
parameters from the simulated data of extant species.

##### Simulate methylome evolution (independent sites)

Two input files are required to start the simulation:

 - `CpGs.bed` contains the location of CpG sites in the genome, for example


     chr1	19464625	19464626	X	 0	 +
     chr1	19464645	19464646	X	 0	 +
     chr1	19464648	19464649	X	 0   +
     chr1	19464651	19464652	X	 0	 +


 - `params.txt` specifies the parameters for methylome evolution.
 The first column contains the required key word, the rest of each line
 contains user defined parameters for simulation.


    NEWICKTREE  ((Human:0.2,Chimp:0.2)HC:0.5,(Mouse:0.4,Rat:0.4)MR:0.3,Dog:0.4)HCMRD;
    PI0	0.3
    G-MATRIX  0.9 0.95
    Q-MATRIX  0.4 0.6
    COV-PARAM 20  0.5
    FG-BETA-PARAMS  0.443 4.314
    BG-BETA-PARAMS  2.576 0.491

```sh
$ simulation -s -o Name -c CpGs.bed params.txt
```

Output files from this command include:

- Simulated methylation levels for individual extant species in `methcount`
format: ``<chromosome> <position> <strand> <simulated state (T/C)>
<methylation level> <coverage>``
```sh    
chr1	19464625	+	T	0.210526	19
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

##### Make table of methylomes of extant species

Two types of tables are accepted by `indep-epi-phylo` and `epiphyte`

- Table of read counts.
To create a table containing the methylation information of the simulated
methylomes of all species, you can use the `merge-methcounts` program from the
[methpipe] package.
```sh
$ merge-methcounts -t Name_Human Name_Chimp Name_Mouse   \
  Name_Rat Name_Dog | awk 'NR==1{OFS="\t";
  print "Human","Chimp", "Mouse", "Rat", "Dog"}
  NR>1{print}' > Name_meth_table
```
The resulting file contains 11 columns (except for the first row), site name
followed by two columns for each species indicating the read-coverage and number
 of methylated reads. The top few lines look like below:
```sh
    Human	Chimp	Mouse	Rat	Dog
    chr1:19464625:+:T	22	0   27  0   23  3   27  0   13  4
    chr1:19464645:+:C	24	1   18  6   27  16  30  28  26  12
    chr1:19464648:+:C	13	13  37  37  31  31  28  27  24  24
```
- Table of hypomethylation probabilities.
```sh
$ for i in `echo Chimp Human Mouse Rat Dog`; do
    echo $i
    indep-sites-hypo -o Name_${i}.hypoprob Name_${i} \
    -params-out Name_${i}.hypoprob.params
  done
$ # merge probabilities into a single file
$ paste Name_Human.hypoprob Name_Chimp.hypoprob Name_Mouse.hypoprob \
  Name_Rat.hypoprob Name_Dog.hypoprob |  \
  awk 'BEGIN{OFS="\t"; print "Human","Chimp","Mouse","Rat","Dog"}
       {print  $1":"$2,$5,$11,$17,$23,$29, $35}' > Name.hypoprob_table
```

##### Estimate phylo-epigenetic tree
With a initial tree specificed by `init.nwk`
```sh
((Human:0.1,Chimp:0.1)HC:0.1,(Mouse:0.1,Rat:0.1)MR:0.1,Dog:0.1)HCMRD;
```
and the table of hypomethylation probabilities ,`Name.hypoprob_table`,
we estimate site-independent phylo-epigenetic model parameters with command:
```sh
$ indep-epi-phylo -v -o Name_hypoprob.out init.nwk Name.hypoprob_table
```
##### Examine the results
In the output file `Name_hypoprob.out`, the first line contains the estimated
branch lengths in the newick format, hypomethylation-to-methylation transition
rate, root node hypomethylation probability. Subsequent lines contains site name
(chromosome:position) and methylation pattern on the whole tree in a pre-order
traversal. `0` indicates hypomethylation and `1` indicates methylation state.
```sh
#((Human:0.223032,Chimp:0.233348)HC:0.504895,(Mouse:0.441373,Rat:0.41714)MR:0.300187,Dog:0.435095)HCMRD:0;	Rate = 0.389063;Root unmeth prob = 0.298923	log-likelihood = -212593
chr1:19464625	00000000
chr1:19464645	10001111
chr1:19464648	11111111
chr1:19464651	00010000
```
Check the rate of correct reconstructions:
```sh
$ paste Name_treestates Name_hypoprob.out | grep -v "#" |cut -f 4,6 | \
  tr '0' 'T' < /dev/stdin | tr '1' 'C' < /dev/stdin |  \
  awk '{if($1==$2) i+=1} END{print i/NR}'
# 0.662636
```

Check the accuracy of methylation state estimates at individual nodes:
```sh  
$ for i in `seq 1 8`; do
    rate=`paste Name_treestates Name_hypoprob.out |  \
      grep -v "#" |cut -f 4,6| tr '0' 'T' < /dev/stdin |  \
      tr '1' 'C' < /dev/stdin |  \
      awk -v node=${i} '{if(substr($1,node,1)==substr($2,node,1)) i+=1}
                        END{print i/NR}' `
  echo node $i accuracy $rate
  done
# node 1 accuracy 0.872479
# node 2 accuracy 0.921771
# node 3 accuracy 0.961141
# node 4 accuracy 0.96066
# node 5 accuracy 0.882314
# node 6 accuracy 0.956536
# node 7 accuracy 0.958562
# node 8 accuracy 0.956041
```

## Example 2
In this example, we simulate a evolution process with interdependent-sites and
use `epiphyte` to estimate the model parameters from the simulated data of
extant species.

##### Simulate methylome evolution (interdependent-sites)
```sh
$ simulation -o Name -c CpGs.bed params.txt
```

As sites were simulated as interdependent, we use the program `hmr_posterior`
from the [methpipe] package, which estimates posterior probabilities of
individual CpG sites being hypomethylated using a hidden Markov model.
```sh
$ for i in `echo Chimp Human Mouse Rat Dog`; do
    echo $i
    ~/git/methpipe/bin/hmr_posterior -v -o Name_${i}.hypoprob Name_${i} \
    -p Name_${i}.hypoprob.params
  done
```

##### Make table of methylomes of extant species
```sh
$ paste Name_Human.hypoprob Name_Chimp.hypoprob Name_Mouse.hypoprob \
  Name_Rat.hypoprob Name_Dog.hypoprob |  \
  awk 'BEGIN{OFS="\t"; print "Human","Chimp","Mouse","Rat","Dog"}
       {print  $1":"$2,$5,$11,$17,$23,$29, $35}' > Name.hypoprob_table
```

##### Estimate phylo-epigenetic tree
```sh
$ epiphyte init.nwk Name.hypoprob_table -i 50 -v -s   \
  -o Name_hme.out -P Name_hme.params 2> Name_RunInfo.txt
```

##### Examine the results
```sh
$ head Name_hme.out
#((Human:0.233754,Chimp:0.254006)HC:0.525329,(Mouse:0.488647,Rat:0.457306)MR:0.316842,Dog:0.454417)HCMRD:0;	pi0=0.262087	Rate=0.394385	g0=0.885621	g1=0.94534
chr1	19465257	19465388	01110000	6	+
chr1	19465595	19465652	11101111	9	+
chr1	19465820	19466136	00000000	22	+
chr1	19466442	19466517	11111011	8	+
chr1	19466522	19467614	11111111	40	+
chr1	19468403	19468464	00000010	6	+
chr1	19470816	19470899	00000000	6	+
chr1	19470973	19471158	00100000	7	+
chr1	19471176	19471369	00000000	8	+
```

Check the accuracy of methylation state estimates at individual nodes:
```sh
$ for i in `seq 1 8`; do   
  rate=`grep -v "#" Name_treestates | paste /dev/stdin Name_hme.out_bysite | \
    cut -f 4,8| tr '0' 'T' < /dev/stdin |  \
    tr '1' 'C' < /dev/stdin |  \
    awk -v node=${i} '{ if(substr($1,node,1)=="C")
      {j+=1; if(substr($1,node,1)==substr($2,node,1)) i+=1}} END{print i/j}' `;
  echo node $i accuracy $rate;
  done
# node 1 accuracy 0.948621
# node 2 accuracy 0.973123
# node 3 accuracy 0.993424
# node 4 accuracy 0.992879
# node 5 accuracy 0.956937
# node 6 accuracy 0.992716
# node 7 accuracy 0.991833
# node 8 accuracy 0.992711
```


License
----



[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [methpipe]:<https://github.com/smithlabcode/methpipe>
