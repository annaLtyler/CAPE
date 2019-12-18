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

These traits encompass both molecular traits and higher level physiological traits, which is 
idea


1. Carter, G. W., Hays, M., Sherman, A. & Galitski, T. Use of pleiotropy to model genetic 
interactions in a population. PLoS Genet. 8, e1003010 (2012).

2. Tyler, A. L., Lu, W., Hendrick, J. J., Philip, V. M. & Carter, G. W. CAPE: an R package 
for combined analysis of pleiotropy and epistasis. PLoS Comput. Biol. 9, e1003270 (2013).