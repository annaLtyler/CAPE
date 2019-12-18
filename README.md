# CAPE with example

This repository contains code to run the Combined Analysis of Epistasis and 
Pleiotropy, or `cape`. This method, described in [1,2] combines information 
across multiple quantitative traits to infer directed genetic interactions. By 
combining information across multiple traits, CAPE not only increases power 
to detect genetic interactions, but also interprets these interactions across 
traits to identify a single interaction that is consistent across all observed data. 
CAPE can be applied to a variety of genetic variants, such as single nucleotide 
polymorphisms (SNPs), copy number variations (CNVs) or structural variations 
(SVs). Included here is also a small example data set with code to infer a predictive 
network between quantitative trait loci (QTL) in a BXD mouse population assayed
for three immune phenotypes. We used phenotypes 16320, 10062, and 13011 from 
GeneNetwork (www.genenetwork.org).
In this example, CAPE generates a genetic interaction network that describes how 
variants interact with each other to influence this group of related traits.

The traits used here are as follows:
 
* **16320** Immune system: ELISA-3x, IgG class antibody binding to TSHR A-subunit 
protein in ELISA 4 weeks after 3 immunizations with TSHR A-subunit adenovirus [OD490 nm]

* **10062** Immune system,  pulmonary system: Tumor necrosis factor (TNF)-alpha level 
in lung after aerosolized lipopolysaccharide (LPS) exposure [pg/ml]

* **13011** Infectious disease, immune function: H1N1 (PR8) influenza A virus (2x10E3 FFU), 
median body weight loss day 7 after infection in 9-14 week-old females [\%]

## Parameters
The parameters for CAPE are supplied with a parameter file. An example is supplied
to show the basic format of the file. The parameters specified through the parameter 
file are as follows:

* **Traits:** Names of traits to be analyzed
* **scan.what** Whether to analyze eigentraits or original traits
* **traits.scaled** Wether to scale traits 
* **traits.normalized**  Whether to normalize trait
* **eig.which** If eigentraits are being analyzed, which of the eigentraits to use
* **pval.correction** What type of p value correction method to use
* **use.kinship** Wether a kinship correction should be implemented
* **kinship.type** What type of kinship correction to implement. Either "overall" or "LTCO"
* **pop**: What type of population is being analyzed:  2PP = two-parent, MPP - multi-parent, RIL = recombinant inbred lines
* **ref.allele** Which allele to use as the reference allele. Typically this is A, but may be different in a multi-parent population.
* **singlescan.perm** How many permutations to run for the single-locus scan. This is only to test single-locus
significance for the user. CAPE does not consider significance when selecting markers for the pair scan.
* **marker.selection.method** The marker selection method to use. Typically top.effects. But from.file can also be used to
specify specific markers to test.
* **peak.density** When *marker.selection.method* is top.effects, this parameter indicates how densely to sample 
markers under effect size peaks. A value of 0.5 indicates that half the markers under an effect size peak will be selected 
for pairwise testing.
* **tolerance** How many markers away from the specified number (*num.alelles.in.pairscan*) can be tolerated.
* **num.alleles.in.pairscan** How many markers should be selected for pairwise testing.
* **max.pair.cor** The maximum correlation between markers for pairwise testing. Testing markers that are highly
 correlated may lead to false positives.
* **pairscan.null.size**  The size of the null distribution to generate for significance testing in the pairwise scan.
We recommend at least 500k for significance testing.

 
 ## References

1. Carter, G. W., Hays, M., Sherman, A. & Galitski, T. Use of pleiotropy to model genetic 
interactions in a population. PLoS Genet. 8, e1003010 (2012).

2. Tyler, A. L., Lu, W., Hendrick, J. J., Philip, V. M. & Carter, G. W. CAPE: an R package 
for combined analysis of pleiotropy and epistasis. PLoS Comput. Biol. 9, e1003270 (2013).