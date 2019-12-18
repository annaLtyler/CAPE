# CAPE with example

This repository contains code to run the Combined Analysis of Epistasis and 
Pleiotropy, or `cape`. This code implements a method, originally described
in [@Carter:2012fd]. The method infers directed interaction networks between 
genetic variants for predicting the influence of genetic perturbations on complex 
traits. This method takes advantage of complementary information in partially 
pleiotropic genetic variants to resolve directional influences between variants
that interact epistatically. `cape` can be applied to a variety of genetic 
variants, such as single nucleotide polymorphisms (SNPs), copy number variations 
(CNVs) or structural variations (SVs). Here we demonstrate the functionality of 
`cape` by inferring a predictive network between quantitative trait loci (QTL) 
in a BXD mouse population assayed for three immune phenotypes. We used
phenotypes 16320, 10062, and 13011 from GeneNetwork (www.genenetwork.org).

The trait descriptions are as follows:
 
* **16320** Immune system: ELISA-3x, IgG class antibody binding to TSHR A-subunit 
protein in ELISA 4 weeks after 3 immunizations with TSHR A-subunit adenovirus [OD490 nm]

* **10062** Immune system,  pulmonary system: Tumor necrosis factor (TNF)-alpha level 
in lung after aerosolized lipopolysaccharide (LPS) exposure [pg/ml]

* **13011** Infectious disease, immune function: H1N1 (PR8) influenza A virus (2x10E3 FFU), 
median body weight loss day 7 after infection in 9-14 week-old females [\%]

These traits encompass both molecular traits and higher level physiological traits and
are moderately correlated.
