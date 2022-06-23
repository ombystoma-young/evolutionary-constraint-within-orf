![[logo-bi-18-7| width=10px]](https://user-images.githubusercontent.com/90496643/169656572-a93ad3c6-2e70-481a-b749-470e02f84e7e.svg#gh-dark-mode-only)
![logo-bi-18-3](https://user-images.githubusercontent.com/90496643/169656574-08b10a55-abe4-401b-bdd2-c9518c4c4f38.svg#gh-light-mode-only)

</br>

# Analysis of variable evolutionary constraint within a single ORF



**Author**  
- Oksana Kotovskaya

**Supervisors**  
- Yury Barbitoff (Bioinformatics Institute)
- Mikhail Skoblov (Research Center of Medical Genetics)


  
### Introduction

Genetic variants leading to loss of function are not found in all genes. If a gene is found under selection pressure, protein truncation variants (PTV) are much less common in them (Cassa C., 2017). Most often, such genes have important functions, and a such catastrophic change in the protein leads to various diseases or death (Samocha K., 2014). In this work, we are interested in the case when the division of genes into conservative (that is, under selection) and non-conservative (that is, free from selection) becomes less unambiguous, namely, cases when non-conservative genes are found in relatively conservative genes. This work is devoted to implementation of algorithm to the search for such sequences.


**Goal**: to estimate the evolutionary conservativeness of individual regions within single ORF.  

  
The key task to achieve this goal was to implement an algorithm based on the hidden Markov model (HMM), which allows to determine the conservativeness of individual regions of the protein-coding sequence (CDS).
  
### Data
In this work, the GRCh37/hg19 assembly of the human genome was used (Frankish A. 2021). Sequence `.fa` and annotation `.gff3` were obtained from this [link](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/).  
Data on existing human PTVs was provided by Yury Barbitoff (Skitchenko R.K., 2020) and was taken from the Genome Aggregation Database (GnomAD) v2.1.1 (Karczewski K.J. et al, 2020).  
Also mutation rates in the trinucleotide context were taken from the work of Karczewski K.J. _et al_ (2020).




### Workflow

First, we needed to obtain data on the coding DNA sequence (CDS) for the gene of interest, namely the nucleotide sequence itself, the CDS coordinates, as well as information about the observed PTV for this gene. The data obtaining scheme is presented below.

![data extraction_dark_2](https://user-images.githubusercontent.com/90496643/169664239-98854a08-111e-4a45-816c-85a10cba0ecc.svg#gh-dark-mode-only)
![data extraction_light_2](https://user-images.githubusercontent.com/90496643/169664240-1bbfa53d-faf8-49c7-ba5b-751b4dd90c42.svg#gh-light-mode-only)

  
#### Stage 1. Sequence extraction  
  We started by getting the CDS coordinates from `.gff3` file. To do this, we used a combination of `grep`, `agrep` and `awk` (presented in the file `get_cds.sh`). In order to calculate the mutation rate in the trinucleotide context (described below), we also needed nucleotides before and after each region of coding DNA sequence (CDS).  
After that, we indexed the file .fa using `samtools faidx` (Danecek P., 2021):  
  
```sh  
! samtools faidx ./data/GRCh37.p13.genome.fa  
``` 
The indexed file was intersected with the previously obtained coordinates to get sequence with `bedtools getfasta` (Quinlan A.R., 2010):  
  
```sh
! bedtools getfasta -fi ./data/GRCh37.p13.genome.fa -bed {gene_name}.gff -name > {gene}.fa
```  

In this way, we obtained the DNA sequences necessary to calculate the expected probability of PTV.
  

#### Stage 2. Counting of PTV mutation rate at the codon in a trinucleotide context  
  
The next step is to calculate the mutation rate in the trinucleotide context for each codon. According to the work of Samocha K., 2014, the best context for determining the variability of a single nucleotide is the inclusion of both 5’ and 3’ flanking nucleotides. 
We used mutation frequencies for each possible variant (G, T, C for A), specified in Supplement Materials to the work of Karczewski K.J., 2020.  
Thus, for each nucleotide, we considered three possible substitutions in a trinucleotide context. At the same time, we consider this base as part of the codon and determine what such a variant leads to: synonymous, missense or nonsense mutation ([frequencies_count](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/5803ca531e6e99400d1802899676b257c62a029a/freqcounter.py#L151) function in file [`freqcounter.py`](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/main/freqcounter.py)).
  
#### Stage 3. Computation of PTV mutation rate per window
  
Next, we divided the gene into regions of fixed length, scanning window ([frequencies_per_window](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/5803ca531e6e99400d1802899676b257c62a029a/freqcounter.py#L278) function in file [`freqcounter.py`](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/main/freqcounter.py)). We did not fix _window size_ it for all genes, because the length of genes varies very much, and in the case of too small values of the window size, it becomes problematic to resolve long non-conservative regions (if there are any). For each window, the sum of the expected frequencies was calculated (example for gene ARFGEF1 presented below).
  
![test_rasterization_dark](https://user-images.githubusercontent.com/90496643/169657497-cd9d0be2-a419-49ee-a576-bae184c2005d.svg#gh-dark-mode-only)
![test_rasterization_light](https://user-images.githubusercontent.com/90496643/169657498-dfb865f8-253e-4b86-8a87-2b2800ae27a6.svg#gh-light-mode-only)



 The resulting distribution of mutation rates corresponds to the expected: the mutation level for missense mutations is higher than the mutation level for synonymous variants, while loss-of-function (LoF) mutations are lower than all. 
  For further analysis, we needed LoF mutation rates $(U)$.
  

#### Stage 4. Counting observed PTV 
Next, we calculated the sum of the observed variants (allele count, $n$) per each window, as well as the average value of allele number ($N$) using function [observed_ptv](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/5803ca531e6e99400d1802899676b257c62a029a/ploffinder.py#L5) in [`ploffinder.py`](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/main/ploffinder.py). `.tsv` with high-confidence pLOF variants were provided by Yury Barbitoff (file [`variants_coords_ac_an_name.tsv`](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/main/data/variants_coords_ac_an_name.tsv)).  
For regions in which no protein truncation variants were detected, we took the $N$ value averaged over the entire gene (since we considered only genes containing at least one PTV).
  
#### Stage 5. Hidden Markov model
  Based on the data obtained, we built a hidden Markov model (presented below).  
  For the sake of simplicity, we have built a model in which only two states are possible: _conservative_ (Cons) and _non-conservative_ (Not) (in the future, the number of states can be increased).

As observations, we choose the allele count per window $n$ (since it is finite, we can consider this quantity discrete).

To assess conservativeness, we used an estimation of the selective effect against heterozygous PTVs ($s_{het}$) that takes into account for each region the observed value of protein truncation variants (PTV) $n$, the allele numbers $N$ and the expected mutation rate.
 

![Evolutionary conservativeness of ORF regions_dark](https://user-images.githubusercontent.com/90496643/169657643-81083110-c36f-4316-9505-a6ef2bf2486f.svg#gh-dark-mode-only)
![Evolutionary conservativeness of ORF regions_light](https://user-images.githubusercontent.com/90496643/169657645-9908d55b-fdb0-49f0-a839-f7b8f93d52f3.svg#gh-light-mode-only)


#### Stage 5.1. Search of HMM parameters. 
As mentioned above, we have chosen $s_{het}$ as an estimate of conservativeness. Similarly to the work of Cassa C. 2017, we assumed the observed distribution of PTV counts across $i$-th region:
  
  $$P(n_i | \alpha, \beta; \nu_i) = \int P(n_i|s_{het};\nu_i) P(s_{het};\alpha, \beta)ds_{het}, $$  
 where $\nu_i = N_iU_i$ - expected genic PTV counts. 
  
  Further, $P(n_i|s_{het};\nu_i) = Pois(n_i, \lambda_i),\ \mbox{where}\ \lambda_i=\nu_i / s_{het}$.  
  $P(s_{het};\alpha, \beta) = IG(s_{het};\alpha, \beta)$ is inverse Gaussian distribution with mean $\alpha$ and shape $\beta$ parameters (calculated for the gene as a whole). 
  Thus,
  
   $$P(n_i | \alpha, \beta; \nu_i) = \int  Pois(n_i, \nu_i / s_{het}) IG(s_{het};\alpha, \beta) ds_{het} $$  
  
  We have chosen as the emission probabilities ($e_k(n_i)$ for $n_i$ in $k$-th state) the normalized probabilities obtained by taking a certain integral ([emission_probabilities](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/5803ca531e6e99400d1802899676b257c62a029a/emissionprobgetter.py#L87) function in [`emissionprobgetter.py`](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/main/emissionprobgetter.py) file):
  
  $$ e_k(n_i) = \frac{P(n_i | \alpha, \beta; \nu_i; a, b)}{P(n_i | \alpha, \beta; \nu_i)} =  \frac{\int\limits_{a}^{b}  Pois(n_i, \nu_i / s_{het}) IG(s_{het};\alpha, \beta) ds_{het}}{P(n_i | \alpha, \beta; \nu_i)},$$
  
  where $a, b = [0, 0.01]$ for $k=Not$ and  $a, b = [0.01, 1]$ for $k=Cons$. The choice of such values is also due to the results obtained in the work  Cassa C. 2017.
  
  Since we had no assumptions about the transition probabilities, we used the Baum–Welch algorithm to find transition probabilities corresponding to the maximum likelihood of the model (function [baum_welch_algo](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/5803ca531e6e99400d1802899676b257c62a029a/hmmalg.py#L41) in file [`hmmalg.py`](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/main/hmmalg.py)).
  
  #### Stage 5.2. Decoding sequence
  
  The decoding of the path of states was carried out using the Viterbi algorithm (function [viterbi_algo](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/5803ca531e6e99400d1802899676b257c62a029a/hmmalg.py#L73) in file [`hmmalg.py`](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/main/hmmalg.py)), it consists in finding a path that also meets the maximum likelihood of the model.
   
  
### Results and future plans
  
 An algorithm based on the hidden Markov model for the estimation of conservativeness of individual regions of the protein-coding sequence was invented and implemented ([conservativeness_estimator](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/5803ca531e6e99400d1802899676b257c62a029a/getconservativeness.py#L7) function in file [`getconservativeness.py`](https://github.com/ombystoma-young/evolutionary-constraint-within-orf/blob/main/getconservativeness.py) gathers all the previous ones and gives the result). Also, a complete pipeline for finding conservative regions for one gene was implemented with `Snakemake` v6.10.0 (Mölder F. 2021). To use it, it's needed to download and index the genome, as described above, and, specifying the desired gene name, window size and initial approximations of the transition probabilities at top of `Snakefile`, run `Snakemake`.
  
The algorithm was tested on set of genes: mostly conservative, mostly non-conservative, presumably possessing regions of conservative, and it solves the first two cases and not so much sensitive for the third.
  
There is a graph below showing the dependence of the allele counts per window for the regions of the gene ARFGEF1. In the case of a large number of alleles, the algorithm with the specified parameters copes with finding the most conservative region.

![ARFGEF1r_tst2_dark](https://user-images.githubusercontent.com/90496643/169666714-04c972ab-888f-4323-9f2c-168913935aaf.svg#gh-dark-mode-only)
![ARFGEF1r_tst2_light](https://user-images.githubusercontent.com/90496643/169666735-c2461c30-aac0-4926-966b-a6b0498091de.svg#gh-light-mode-only)


  
In the future, it is planned to establish whether it is possible to fix transition probabilities based on the theory of population evolution. 
Also change the approach to calculating the allele number for regions where no PTV is observed (the current estimate is rather rough, as well as $a$, $b$ for emission probabilities). 



  

 **Requirements for the pipeline**: 
 - agrep tool;
 - Python v3.9.10 (the requirements for computing packages are specified in `requirements.txt`);
 - samtools v1.14;
 - bedtools v2.30.0;
 - Snakemake v6.10.0.

  

  
  
### Literature  
- Cassa, C., Weghorn, D., Balick, D. _et al_. Estimating the selective effects of heterozygous protein-truncating variants from human exome data. _Nat Genet_ 49, 806–810 (2017). 
- Danecek, P., Bonfield, J.K., Liddle, J., Li, H. _et al_. Twelve years of SAMtools and BCFtools. _GigaScience_, 10(2), (2021).  
- Davis, R.I. and Lovell, B.C. Comparing and evaluating HMM ensemble training algorithms using train and test and condition number criteria. _Pattern Anal. Appl._ 6, 4, pp. 327–336, (2003).
- Frankish, A., Diekhans, M., Jungreis, I., _et al_. GENCODE 2021. _Nucleic Acids Research_. 49(D1):D916-D923 (2021). 
- Karczewski, K.J., Francioli, L.C., Tiao, G. _et al_. The mutational constraint spectrum quantified from variation in 141,456 humans. _Nature_, 581, pp. 434–443 (2020).
- Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J.. Sustainable data analysis with Snakemake. _F1000Res_ 10, 33. (2021)
- Quinlan, A.R., Hall, I.M., BEDTools: a flexible suite of utilities for comparing genomic features, _Bioinformatics_, 26(6), pp. 841–842 (2010).
- Samocha, K., Robinson, E., Sanders, S. _et al_. A framework for the interpretation of _de novo_ mutation in human disease. _Nat Genet_ 46, 944–950 (2014).
- Skitchenko, R.K., Kornienko, J.S., Maksiutenko, E.M., Glotov, A.S., Predeus, A.V., Barbitoff, Y.A. Harnessing population-specific protein truncating variants to improve the annotation of loss-of-function alleles. _bioRxiv_ 2020.08.17.254904, (2020).
