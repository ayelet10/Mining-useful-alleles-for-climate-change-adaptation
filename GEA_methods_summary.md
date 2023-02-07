# Genotype-environment association (GEA) Methods #

In this document, the following methods for GEA are reviewed:
categorical tests, logistic regressions, matrix correlations, general
linear models, mixed effects models, machine learning algorithms.

# Introduction to General GEA approaches (also called environmental association analysis (EAA) or landscape genomics)

GEA methods can identify adaptive loci based on correlations between
genetic and environmental data. The environmental association analysis
can detect signals of adaptive patterns that might not be detected by
the traditional population genetic differentiation methods. They may be
split into two main approaches, Univariate and Multivariate methods. The
univariate methods are very common even though the data is high
dimensional. The multivariate approaches analyze many loci together and
they should be advantageous in estimating how sets of markers covary in
response to environment.

## Sampling design 

In Rellstab et al. 2015 they talk about the importance of sampling
design that is tricky in case you want to replicate the experiment. If
you sample along a gradient, it could be hard to track the exact
conditions in a replication experiment. It is suggested by Manal et al.
2012 that a modeled-based approach should be considered in designing the
sampling strategy, instead of random sampling that is usually done.
Another sampling method that showed effective is to "**sample scattered
and random pairs of closely situated populations that exhibit
substantial differences in environmental conditions while being within
geneflow distance"** it was introduced by Lotterhos & Whitlock 2015 that
showed by simulations that it increased power in detecting true
positives compared to random or transect designs.

# Methods description

## Logistic regressions

### **SAMᵦADA** 

This landscape genomic software designed to run *univariate or
multivariate* logistic regression between the presence of a genotype and
one or several environmental variables.

It is an extended version of SAM (Spatial Analysis method), developed to
overcome some of the limitations of SAM. The software includes the
possibility of multivariate analyses testing, enabling the introduction
of neutral genetic structure as an additional factor. SAMᵦADA can
further quantify the level of spatial autocorrelation of genotypes.

Advantages: It is claimed to be substantially faster than BAYENV2 and
LFMM with the univariate model (i.e. not including neutral genetic
structure) and faster than BAYENV2 with a bivariate model. SAMᵦADA comes
with a module that can split and remerge large data files. Hence,
analyses can be run on different processors in parallel.

Can be downloaded as R package 'R.SamBada':
<https://cran.r-project.org/web/packages/R.SamBada/index.html>

Reference: <https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12629>

Manual and examples can be found here:
<https://github.com/Sylvie/sambada> ;
<https://cran.r-project.org/web/packages/R.SamBada/R.SamBada.pdf>

## Matrix Correlations 

### **Partial Mantel tests** 

In matrix correlations, one aims to test for correlation between
matrices that express distances or dissimilarities between sampling
units. A simple Mantel test estimates the strength of correlation
(linear or rank linear) between two distance matrices (Mantel 1967) and
computes a P-value for the correlation coefficient in a permutation
procedure. As an extension, the partial Mantel test checks if there is a
correlation between two distance matrices given a third matrix (Smouse
et al. 1986). In EAA, partial Mantel tests can be used with individual
or population data. The first matrix includes pairwise genetic distances
or differentiation among individuals or populations at particular loci,
the second matrix consists of environmental distances between sampling
locations, and the third matrix can be used to control for genetic
structure with neutral pairwise genetic distances.

Advantages: it deals with distances and does not rely on any parametric
assumptions.

Disadvantages: if there is spatial autocorrelation in the two matrices,
Mantel tests result in P-values that are not well calibrated, because
the permutation procedure fails to produce a valid null hypothesis.
(This can be overcome by ignoring P-values and concentrating on effect
sizes instead (i.e. the correlation coefficient r) when identifying top
associations between loci and environmental factors.)

R package: 'vegan'

Tutorial:
<https://rfunctions.blogspot.com/2016/10/simple-and-partial-mantel-tests.html>

## General Linear Models 

General linear models are statistical models in which a response
variable is modelled as a linear function of some set of explanatory
variables.

Advantages: These models can account for neutral genetic structure and
include statistical methods largely familiar to biologists.

## Multiple linear regressions and univariate general linear models

### Standard linear modelling in R. 

### **TASSEL** 

TASSEL is a software package used to evaluate traits associations,
evolutionary patterns, and linkage disequilibrium.

Advantages:

1\. The opportunity for several new and powerful statistical approaches
to association mapping such as a General Linear Model (GLM) and Mixed
Linear Model (MLM).

2\. An ability to handle a wide range of indels (insertion & deletions).
Most software ignore this type of polymorphism; however, in some species
(like maize), this is the most common type of polymorphism.

Reference: <https://pubmed.ncbi.nlm.nih.gov/17586829/>

Downloads and manual - <https://www.maizegenetics.net/tassel>

rTASSEL - R package that connects a variety of highly used TASSEL
methods and analytical tools. See this link for demo, installation, and
user instructions: <https://github.com/maize-genetics/rTASSEL>

## Canonical correlations and multivariate linear regressions

The general linear model framework can be extended to models with
*multivariate* response variables to [account for the polygenic
architecture of adaptive traits]{.underline}.

### **CCA** 

The most popular method is canonical correlation analysis (CCA), which
finds the linear combinations of two sets of variables -- multiple loci
and multiple environmental factors -- that are maximally correlated
(Legendre & Legendre 2012). The results are orthogonal sets of canonical
variables that can be tested for significance.

Disadvantages:

-   Strong patterns of multicollinearity could skew the results.
    Moreover, as CCA does not allow missing data, global deletion of
    samples or imputation of missing values is often required.

Advantages:

-   Mosca et al. (2012) were able to show with only several hundred
    SNPs, how geographic factors shape the population genetic structure
    of four subalpine conifer tree species in the European Alps.

### **Constrained ordinations (****RDA)** 

RDA is a *multivariate* ordination technique that can be used to analyze
many loci and environmental predictors simultaneously. RDA determines
how groups of loci covary in response to the multivariate environment,
and can detect processes that result in weak, multilocus molecular
signatures (Rellstab et al., 2015; Forester et al., 2018).

RDA is a two-step analysis in which genetic and environmental data are
analyzed using **multivariate linear regression**, producing a matrix of
fitted values. Then PCA of the fitted values is used to produce
canonical axes, which are linear combinations of the predictors
(Legendre & Legendre, 2012). RDA can be used to analyze genomic data
derived from both individual and population-based sampling designs.

The RDA method was recommended as a useful approach to test hypotheses
about **specific environmental factors** is redundancy analysis
(Legendre & Legendre 2012, Forester et al. 2018).

Advantage: It allows for building and testing models of varying
complexity, including those that condition results based on neutral
genetic structure or spatial effects, referred to as partial RDA (pRDA).
Significance of the model, each synthetic orthogonal axis and each
explanatory variable can be tested using a permutation-based analysis of
variance (Legendre & Legendre 2012).

In a simulation study, RDA showed a superior combination of low false
positive and high true positive rates across weak, moderate, and strong
multilocus selection. These results were robust across the levels of
population structure, demographic histories, sampling designs, and
sample sizes tested (Forester et al., 2018).

Disadvantage:

Lotterhos (2022 preprint) analyzed simulated data and concluded that RDA
analysis is not suited for detecting QTNs since it had a high false
positive rate and low power, regardless of whether structure was
accounted for or not.

R package for constrained ordinations (CCA, RDA) is 'vegan'

Tutorials:
<https://www.mooreecology.com/uploads/2/4/2/1/24213970/constrained_ordination.pdf>

<https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html>

<https://link.springer.com/content/pdf/10.1007%2F978-1-4419-7976-6_6>

<https://popgen.nescent.org/2018-03-27_RDA_GEA.html#data-packages>

References:

Borcard D, Gillet F, Legendre P (2011) [*Numerical Ecology with
R*](http://www.springer.com/us/book/9781441979759). Springer, New York.

Legendre P, Legendre L (2012) [*Numerical Ecology*, 3rd
edition](https://www.elsevier.com/books/numerical-ecology/legendre/978-0-444-53868-0).
Elsevier, Amsterdam.

## Mixed effects models 

In general**: allele frequencies** of individuals or populations are
treated as response variables, environmental factors are used as fixed
factors, whereas neutral genetic structure is incorporated as a random
factor.

The difference between the approaches is in how significance is tested,
how neutral genetic structure is incorporated, and which type of
genotype-- environment association (linear/rank-linear/logistic) is
assumed.

General advantage: they provide a unified statistical framework for
controlling for the effects of neutral genetic structure.

## Mixed effects models -- specific approaches

### **BAYENV** and **BAYENV2** 

A Bayesian method that estimates the empirical pattern of covariance in
allele frequencies between populations from a set of markers, and then
uses this as a null model for a test at individual SNPs.

For a given genetic variant, BAYENV tests whether a model that includes
an environmental factor has an improved fit to the data compared to a
null model that includes only neutral genetic structure, which is
represented by a covariance matrix of estimated allele frequencies.
BAYENV delivers Bayes factors for each locus--variable combination

Advantages:

> -BAYENV allows for the incorporation of uncertainty of allele
> frequencies that arises from differences in sample sizes.
>
> -Low rate of false positives (De Mita et al. 2013)
>
> \- Perform best under scenarios with weak hierarchical genetic
> structure (de Villemereuil et al. 2014)
>
> \- Only BAYENV2 (Gunther & Coop, 2013) yet accounts for the variance
> introduced by variation in sequencing coverage in Pool-Seq. (Rellstab
> et al. 2015)

Disadvantages:

> \- BAYENV is slow because it is computationally very intensive.
>
> \- It is not applicable to individual and scattered sampling designs.
>
> \- The output factors may not be directly compared across
> environmental variables because of variable-specific value ranges.
>
> \- Run-to-run variation of BAYENV (version 1) can be large (Blair et
> al 2014), so it is advised to average Bayes factors among multiple
> runs to produce more stable and reliable results.

Reference:

Günther, Torsten, and Graham Coop. "Robust identification of local
adaptation from allele frequencies." Genetics vol. 195,1 (2013): 205-20.
doi:10.1534/genetics.113.152462

Download: <https://bitbucket.org/tguenther/bayenv2_public/src/master/>

### **GINLAND** 

gINLAnd is a spatial generalized mixed model (SGLMM) method, that
accounts for spatial autocorrelation.

Software: gINLAnd (Guillot et al. 2014), a spatial generalized mixed
model (SGLMM) which uses a Markov chain Monte Carlo (MCMC)-free approach
with shorter computing time (compared to BAYENV).

The program essentially quantifies the magnitude of a statistical
dependence between a loci displaying an outstanding statistical
dependence with a certain environmental variable and computes a
statistical measure of how likely it is to observe this dependence by
chance. Computations are based on the Integrated Nested Laplace
Approximation (INLA) method proposed by Rue et al. (2009) and the
connection between Stochastic Partial Differential Equations (SPDEs) and
Gaussian Markov Random Field (GMRF) introduced by Lindgren et al (2011).

Download and documentation:
<https://i-pri.org/special/Biostatistics/Software/gINLAnd/gINLAnd_doc.html>

Advantages:

-   Shorter computing time

-   Considers pure spatial autocorrelation based on a geographical
    distance matrix.

-   The main tasks, namely inference and Bayes factors computation, can
    be carried out via the graphical user interface that does not
    require any knowledge about R.

Reference:

Guillot, G., Vitalis, R., le Rouzic, A., & Gautier, M. (2014). Detecting
correlation between allele frequencies and environmental variables as a
signature of selection. A fast computational approach for genome-wide
studies. Spatial Statistics, 8, 145-155.
<https://www.sciencedirect.com/science/article/pii/S2211675313000468>

## Latent factor mixed models (LFMMs)

In **LFMMs**, neutral genetic structure is introduced as a random factor
with the so-called latent factors, which are similar to principal
components and calculated from all available markers.

It can be used for fitting latent factor mixed models and evaluating
association between a response matrix (SNP genotype or methylation
levels) and a variable of interest (phenotype or exposure levels) in
genome-wide (GW), genome-environment (GE), epigenome-wide (EW)
association studies.

Before starting the final analysis, the number of latent factors (K) has
to be chosen, either by an analysis of histograms of test P-values for
different K-values (i.e. it should look similar to a uniform
distribution), by performing a Tracy--Widom test on the eigenvalues of a
PCA on the genetic data, or using programs such as STRUCTURE (Pritchard
et al. 2000) to determine plausible values for K.

Advantages:

-   The advantage of this linear approach is that the effects of
    environmental factors and neutral genetic structure on allele
    frequencies are simultaneously estimated

-   Computing time is reasonably fast, making LFMM attractive for EAA
    with whole genomes or subsets of large random batches of SNPs in
    parallel.

-   LFMM is suited for both population based and scattered,
    individual-based sampling designs (although the graphical
    user-interface is not available for handling population allele
    frequencies so that should be done by using the command-line
    version)

-   Low rates of false positives and negatives and that it performs
    slightly better than BAYENV in detecting weak selection (Frichot et
    al., 2013)

-   LFMM provides the best compromise between detection power and error
    rates in situations with complex hierarchical neutral genetic
    structure and polygenic selection. (de Villemereuil et al.,2014)

-   LFMM is quite robust to a variety of sampling designs and underlying
    demographic models. (Lotterhos & Whitlock ,2015).

> Disadvantages:

-   As the stochastic algorithm of LFMM (MCMC) does not provide exact
    results, Frichot et al. (2013) recommend performing **multiple
    runs**, use the median of the resulting Z-scores and adjust their
    P-values as described in the software manual.

-   Power was scenario-dependent and was reduced (relative to
    ordinations) for loci under weak selection (Forester et al., 2018,
    de Villemereuil et al.,2014).

Input formats and command line version tutorial:
<http://membres-timc.imag.fr/Olivier.Francois/lfmm/files/note.pdf>

Download R package 'lfmm': <https://bcm-uga.github.io/lfmm/>

Tutorial:
<https://cran.r-project.org/web/packages/lfmm/vignettes/lfmm.html>

Another function, 'lfmm2' is available in R package 'LEA'
<https://rdrr.io/bioc/LEA/man/lfmm2.html>

Authors Caye K, Jumentier B, Lepeule J, Francois O. (2019). LFMM 2: fast
and accurate inference of gene-environment associations in genome-wide
studies. Molecular biology and evolution, 36(4), 852-860.

References: Frichot et al. 2013; Lotterhos & Whitlock ,2015; de
Villemereuil et al.,2014

## GWAS mixed models 

Methods based on GWAS, replacing the response variable phenotype by
environment.

### **EMMA** 

Kang et al. (2008) developed an efficient mixed-model association (EMMA)
method that includes a simple identity-by-state allele sharing kinship
matrix to control for neutral genetic background

Disadvantages**:**

-   **EMMA** is optimized to test associations of only one allele with
    climate

-   Allowing heterozygous genotypes of outbred individuals is possible,
    but complex and computationally intensive (Kang et al. 2008).

-   The use of a kinship matrix to describe neutral genetic structure of
    populations may be inappropriate.

-   Similarly, a linear mixed-model method is implemented in the
    software TASSEL (Bradbury et al. 2007).

-   GWAS mixed models are designed for individual rather than population
    sampling, making them best suited for analyses with samples
    continuously distributed across a study region.

Downloads: <http://mouse.cs.ucla.edu/emma_jemdoc/>

Documentation is found here:
<https://cran.r-project.org/web/packages/EMMAgeo/EMMAgeo.pdf>

Reference:

Kang, H. M., Zaitlen, N. A., Wade, C. M., Kirby, A., Heckerman, D.,
Daly, M. J., & Eskin, E. (2008). Efficient control of population
structure in model organism association mapping. Genetics, 178(3),
1709--1723. <https://doi.org/10.1534/genetics.107.080101>

## *F~ST~*-based approach 

### **BayeScEnv** 

*F~ST~* -based genome-scan method, which incorporates environmental
information in the form of 'environmental differentiation'. It is based
on the *F* model, but, as opposed to existing approaches, it considers
two locus-specific effects: one due to divergent selection and the other
due to various other processes different from local adaptation (e.g.
range expansions, differences in mutation rates across loci or
background selection).

Advantages:

-   Lower false positive rate than an existing *F~ST~*-based method,
    BayeScan, under a wide range of demographic scenarios.

-   Was shown in human data that it can be used successfully to study
    local adaptation.

> Code is available here, in C++ :
> <http://github.com/devillemereuil/bayescenv>
>
> Reference:
>
> de Villemereuil, P. and Gaggiotti, O.E. (2015), A new *F~ST~*-based
> method to uncover local adaptation using environmental variables.
> Methods Ecol Evol, 6:
> 1248-1258. <https://doi.org/10.1111/2041-210X.12418>

## Machine-learning algorithm

### Random Forest

Random Forest (RF) is a machine-learning algorithm that is designed to
identify structure in complex data and generate accurate predictive
models. It is based on classification and regression trees (CART), which
recursively partition data into response groups based on splits in
predictors variables. CART models can capture interactions,
contingencies and nonlinear relationships among variables,
differentiating them from linear models (De'ath & Fabricius, 2000). RF
reduces some of the problems associated with CART models (e.g.,
overfitting and instability) by building a "forest" of classification or
regression trees with two layers of stochasticity: random bootstrap
sampling of the data and random subsetting of predictors at each node
(Breiman, 2001). This provides a built-in assessment of predictive
accuracy (based on data left out of the bootstrap sample) and variable
importance (based on the change in accuracy when covariates are
permuted). For GEA, variable importance is the focal statistic, where
the predictor variables used at each split in the tree are molecular
markers, and the goal is to sort individuals into groups based on an
environmental category (classification) or to predict home environmental
conditions (regression). Markers with high variable importance are best
able to sort individuals or predict environments.

Advantages: RF has been used in several GEA and GWAS studies (e.g.,
Brieuc, Ono, Drinan, & Naish, 2015; Holliday, Wang, & Aitken, 2012;
Laporte et al., 2016; Pavey et al., 2015).

Disadvantages: Random Forest had lower detection rates overall and
performed poorly compared to univariate GEAs and constrained ordinations
(RDA). It was previously noted that RF is not suited to identifying weak
multilocus selection or interaction effects in these large data sets.
The reason is that most SNPs are not under selection, which results in
detection that is dominated by strongly selected loci (Forester et al.,
2018).

R package: 'randomForest'

Tutorial: <https://www.tutorialspoint.com/r/r_random_forest.htm>

Reference: Brieuc, MSO, Waters, CD, Drinan, DP, Naish, KA. A practical
introduction to Random Forest for genetic association studies in ecology
and evolution. Mol Ecol
Resour. 2018; 18: 755-- 766. <https://doi.org/10.1111/1755-0998.12773>

Book chapter on random forest in genomic predictions:

Montesinos López, O.A., Montesinos López, A., Crossa, J. (2022). Random
Forest for Genomic Prediction. In: Multivariate Statistical Machine
Learning Methods for Genomic Prediction. Springer, Cham.
<https://doi.org/10.1007/978-3-030-89010-0_15>

### Gradient Forest

Gradient Forest (GF) is a machine learning algorithm designed to analyze
spatial patterns of biodiversity as a function of environmental
gradients.

This prediction-based approach can be used to predict how specific
genotypes perform in different environments. In Láruson et al. 2021,
they use simulations to test GF use as an offset measure that predicts
the loss of environmentally adapted alleles under rapid environmental
change. This GF Offset is commonly used to measure the offset between
the GF-predicted environmental association of adapted alleles and a new
environment.

Disadvantages:

\- Not yet been tested for increased demographic complexities.

\- Sampling scheme might bias environmental predictor importance values.

Download "gradientForest" R package (Ellis et al., 2012)
<https://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf>

Examples:
<https://gradientforest.r-forge.r-project.org/biodiversity-survey.pdf>

Reference: Láruson et al. 2021
<https://onlinelibrary.wiley.com/doi/full/10.1111/eva.13354>


<table>
<colgroup>
<col style="width: 14%" />
<col style="width: 8%" />
<col style="width: 8%" />
<col style="width: 8%" />
<col style="width: 10%" />
<col style="width: 11%" />
<col style="width: 8%" />
<col style="width: 7%" />
<col style="width: 8%" />
<col style="width: 11%" />
</colgroup>
<thead>
<tr class="header">
<th><strong>Method</strong></th>
<th><strong>Reference</strong></th>
<th><strong>Association type</strong></th>
<th><strong>Sampling design</strong></th>
<th><strong>Incorporation of neutral genetic structure</strong></th>
<th><strong>Incorporation of spatial autocorrelation</strong></th>
<th><strong>Individual / population data</strong></th>
<th><strong>Mode for pooled data</strong></th>
<th><strong>Correction for sample size</strong></th>
<th><strong>Software/ R package</strong></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Categories</td>
<td></td>
<td>Categorical</td>
<td>Categorical</td>
<td>Possible</td>
<td>Possible</td>
<td>Both</td>
<td>Possible</td>
<td>Possible</td>
<td>Various<br />
statistical<br />
methods</td>
</tr>
<tr class="even">
<td>Spatial analysis method<br />
(SAM)</td>
<td>Joost et al. (2007)</td>
<td>Logistic</td>
<td>Gradient/<br />
scattered</td>
<td>Possible (in<br />
SAMᵦADA)</td>
<td>Possible (in<br />
SAMᵦADA)</td>
<td>Individual</td>
<td>No</td>
<td>No</td>
<td>SAM (Joost et al.<br />
2008), SAMᵦADA<br />
(Stucki et al.<br />
submitted)</td>
</tr>
<tr class="odd">
<td>Multiple logistic regression</td>
<td></td>
<td>Logistic</td>
<td>Gradient/<br />
scattered</td>
<td>Possible</td>
<td>Possible</td>
<td>Individual</td>
<td>No</td>
<td>No</td>
<td>R (R<br />
Development<br />
Core Team<br />
2011)</td>
</tr>
<tr class="even">
<td><p>Generalized estimating equations</p>
<p>(GEEs)</p></td>
<td>Carl &amp; Kuhn<br />
(2007),<br />
Poncet et al.<br />
(2010)</td>
<td>Logistic</td>
<td>Gradient/<br />
scattered</td>
<td>No</td>
<td>Yes</td>
<td>Individual</td>
<td>No</td>
<td>No</td>
<td>GEEPACK (Yan &amp;<br />
Fine 2004)</td>
</tr>
<tr class="odd">
<td>Partial Mantel test</td>
<td>Smouse<br />
et al. (1986)</td>
<td>Linear/<br />
rank-linear</td>
<td>Gradient/<br />
scattered</td>
<td>Yes</td>
<td>Possible</td>
<td>Both</td>
<td>No</td>
<td>No</td>
<td>ECODIST (Goslee<br />
&amp; Urban 2007),<br />
VEGAN<br />
(Oksanen et al.<br />
2013)</td>
</tr>
<tr class="even">
<td>Multiple linear regression/General linear models</td>
<td></td>
<td>Linear</td>
<td>Gradient/<br />
scattered</td>
<td>Possible</td>
<td>Possible</td>
<td>Both</td>
<td>No</td>
<td>No</td>
<td>R (R<br />
Development<br />
Core Team<br />
2011), TASSEL<br />
(Bradbury et al.<br />
2007)</td>
</tr>
<tr class="odd">
<td><p>Canonical correlation analysis</p>
<p>(CCA)</p></td>
<td>Legendre &amp;<br />
Legendre<br />
(2012)</td>
<td>Linear</td>
<td>Gradient/<br />
scattered</td>
<td>Possible</td>
<td>Possible</td>
<td>Both</td>
<td>No</td>
<td>No</td>
<td>VEGAN (Oksanen<br />
et al. 2013)</td>
</tr>
<tr class="even">
<td>(Partial) redundancy analysis<br />
(RDA)</td>
<td>Legendre &amp;<br />
Legendre<br />
(2012)</td>
<td>Linear</td>
<td>Gradient/<br />
scattered</td>
<td>Possible</td>
<td>Possible</td>
<td>Both</td>
<td>No</td>
<td>No</td>
<td>VEGAN (Oksanen<br />
et al. 2013)</td>
</tr>
<tr class="odd">
<td>BAYENV</td>
<td>Coop et al.<br />
(2010)</td>
<td>Linear/<br />
rank-linear</td>
<td>Gradient/<br />
scattered</td>
<td>Yes</td>
<td>No</td>
<td>Population</td>
<td>Yes (in<br />
BAYENV2)</td>
<td>Yes</td>
<td>BAYENV (Coop<br />
et al. 2010),<br />
BAYENV2<br />
(Gunther &amp; €<br />
Coop 2013)</td>
</tr>
<tr class="even">
<td>Spatial generalized linear mixed model<br />
(SGLMM)</td>
<td>Guillot et al.<br />
(2014)</td>
<td>Linear</td>
<td>Gradient/<br />
scattered</td>
<td>Yes</td>
<td>Yes</td>
<td>Both</td>
<td>No</td>
<td>Yes</td>
<td>GINLAND (Guillot<br />
et al. 2014)</td>
</tr>
<tr class="odd">
<td>Latent factor mixed models<br />
(LFMMs)</td>
<td>Frichot et al.<br />
(2013)</td>
<td>Linear</td>
<td>Gradient/<br />
scattered</td>
<td>Yes</td>
<td>No</td>
<td>Both</td>
<td>No</td>
<td>No</td>
<td><p>LFMM (Frichot<br />
et al. 2013),</p>
<p>LEA<br />
(Frichot &amp;<br />
Francois 2015)</p></td>
</tr>
<tr class="even">
<td>GWAS mixed models</td>
<td></td>
<td>Linear</td>
<td>Gradient/<br />
scattered</td>
<td>Yes</td>
<td>No</td>
<td>Individual</td>
<td>No</td>
<td>No</td>
<td><p>EMMA</p>
<p>(Kang et al. 2008),<br />
TASSEL<br />
(Bradbury et al.<br />
2007),</p>
<p>LME4 (Bates et al. 2014)</p></td>
</tr>
<tr class="odd">
<td><em>F<sub>ST</sub></em>-based methods</td>
<td>Frichot et al.<br />
(2013)</td>
<td>Differentiation-based</td>
<td>Gradient/<br />
scattered</td>
<td>Yes</td>
<td>No</td>
<td>Both</td>
<td>No</td>
<td>Yes</td>
<td><p>BAYESCENV</p>
<p>(de Villemereuil &amp;<br />
Gaggiotti in<br />
press)</p></td>
</tr>
<tr class="even">
<td>Random forest Machine-learning algorithms</td>
<td>Brieuc et al. 2018</td>
<td>Classification and regression</td>
<td>Gradient/ scattered</td>
<td>No</td>
<td>No</td>
<td>Population</td>
<td>No</td>
<td>No</td>
<td><p>randomForest</p>
<p>Breiman (2001)</p></td>
</tr>
<tr class="odd">
<td>Gradient Forest algorithm</td>
<td>Láruson et al. 2021</td>
<td>Classification and regression</td>
<td>Gradient/ scattered</td>
<td>No</td>
<td>No</td>
<td>Population</td>
<td>No</td>
<td>No</td>
<td><p>gradientForest</p>
<p>(Ellis et al., 2012)</p></td>
</tr>
</tbody>
</table>

Table reproduced from Rellstab et al. 2015

# References: 

Rellstab, C, Gugerli, F, Eckert, AJ, Hancock, AM, Holderegge,r R. A
practical guide to environmental association analysis in landscape
genomics. Mol Ecol. 2015 Sep;24(17):4348-70. doi:10.1111/mec.13322.
PMID:
26184487.<https://www.umt.edu/ces/conferences/congen/congen2019/images/all%20together/TGraves.pdf>

Forester, BR, Lasky, JR, Wagner, HH, Urban, DL. Comparing methods for
detecting multilocus adaptation with multivariate genotype--environment
associations. Mol
Ecol. 2018; 27: 2215-- 2233. <https://doi.org/10.1111/mec.14584>

Láruson, ÁJ, Fitzpatrick, MC, Keller, SR, Haller, BC, Lotterhos, KE.
Seeing the forest for the trees: Assessing genetic offset predictions
from gradient forest. Evolutionary
Applications, 2022; 15, 403-- 416. <https://doi.org/10.1111/eva.13354>
