# Male Mutation Bias in Primates
## To use the scripts:

	—> Sample input file: “hg18.farFromGenes.20kb.nosep.schweinfurthii.bothsexes.reheader”
	—> The input file has data for all the chromosomes (chr1-22 & chrX)

## Format of the input file:
	chrom   start   stop    cm2genes        	pi      		NumBases        Dist2Root
	chr1    32450   52450   3.98951152150029e-08    0       		1240    11
	chr1    62254   82254   1.46132807076949e-08    7.44238095238095        1785    52.2083333333333
	chr1    82254   102254  1.37304716842559e-07    0       		0       0
	chr1    102254  122254  1.46071458287778e-08    0       		0       0
	chr1    132091  152091  1.01900332726624e-08    0       		0       0
	chr1    152091  172091  1.32887604286096e-07    0       		0       0
	chr1    172091  192091  2.55585175299423e-07    0       		0       0
	chr1    192091  212091  3.78282746313204e-07    0       		0       0
	chr1    212091  232091  5.00980317326777e-07    0       		0       0
	chr1    232091  252091  3.78295016070769e-07    0       		0       0
	chr1    252091  272091  2.55597445058046e-07    0       		0       0
	chr1    272091  292091  1.32899874045323e-07    0       		0       0
	chr1    292091  312091  1.02023030326007e-08    0       		0       0

**Below are the fields given in the input data file:**

**chrom**:     	name of the chromosome	
**start**:     	start position of the chromosome		
**stop**:      	stop position of the chromosome		
**cm2genes**:  	distance of the defined region of the chromosome in centiMorgans from the nearest gene		
**pi**:        	diversity value of the defined region		
**NumBases**:  	number of callable bases in the defined region	
**Dist2Root**: 	divergence value of the defined region	

***************************************************************************************************************************************

##Mean and Confidence Intervals are calculated in two steps:

1. Generate bootstrap data (the program generates 1000 bootstrap files): **generate_bootstrap.R**	
	This script creates 1000 bootstrap files for autosomes and 1000 bootstrap files for X chromosome	

2. Calculate confidence intervals: **CalculateMean_CI.R**	 
 	This script calculates the following:	
	(i)   Mean and 95% confidence intervals for X and autosomal Divergence for each bin 	
	(ii)  Mean and 95% confidence intervals for X to autosomal divergence ratios, both uncorrected and corrected for ancestral polymorphism, for each bin 	
	(iii) Mean and 95% confidence intervals for Alpha, both uncorrected and corrected for ancestral polymorphism, for each bin	

Defined Bins:
	Increasing genetic distances from the nearest genes (in centiMorgans): ([0-0.05], [0.05-0.1], [0.1-0.2], [0.2-0.4], [0.4-0.8], [0.8-2.0]). 
