# Introduction

In Mammals, life starts when a male sperm cell fertilizes a female oocyte to form a zygote. This unique totipotent cell can give rise to a complete multicellular organism. This formation relies on a series of cellular differentiation and spatial organization events. Much insight into Mammalian development has come from the mouse, which benefits from a rapid development (21 days) and provides unrivalled accessibility to embryos. Regarding human development however, differences have emerged, and many questions remain unsolved ([Mole *et al*, 2020](#mole_2020)). Thanks to the advent of single-cell multi-omics techniques, we have now access to the very intimate molecular processes occuring in the early human embryo. However, exploration of these data has been optimized, and notably, a careful integrated analysis of the different omics levels has not been conducted so far.

In this study, we propose to apply integrated bioinformatics analyses to transcriptomic and methylomic data of the peri-implantation human embryo, in order to gain insight into the interconnection of these different molecular levels. For this purpose, we use single-cell RNA-Seq and PBAT data from ([Zhou *et al*, 2019](#zhou_2019)), extracted from *in vitro* cultured human embryos from day 6 to 12 of development. These data are particularly suited for our analysis, because transcriptomic and DNA methylation data were obtained from the same cell.

First, we normalize the data to make them suitable for secondary analysis algorithms. Then, we start investigating each dataset separately. We determine inter-individual similarities by dimension reduction approaches such as Principal Component Analysis (PCA). We cross-validate the results with unsupervised hierarchical clustering.

Next, we perform integrated analysis of the two datasets with different approaches:

- Multi-Omics Factor Analysis (MOFA)
- Kernel computation

We follow the analysis by investigating similarities between gene expression and DNA methylation patterns through the sample this time. For this purpose, we use Weighted Gene Correlation Network Analysis (WGCNA) in order to identify correlated gene modules.

Finally, we consult selected databases for functional enrichment analysis, and draw conclusions and perspectives

# Data preprocessing

## Data exploration

These two datasets consist of 130 cells from the three earliest cell lineages:

- the epiblast (EPI), which is pluripotent and give rise to the fetus
- the primitive endoderm (PE), which produces the yolk sac
- the trophectoderm (TE), which yields the placenta.

The table below summarizes the distribution of cells into lineages and embryonic days.
## Normalisation of data distribution

Input data have been pre-processed to some extent: methylation levels range from 0 to 1 while RNA-Seq counts are given in FPKM. However, distribution of feature values do not follow a Gaussian law. We need to normalise the data prior to further processing, in order to increase the robustness of secondary analyses.

By looking at total measures per individual for each dataset, we observe that library sizes for DNA methyaltion differ within sample, likely due to technical issues related to sequencing depth. So we decide to standardize this size by scaling units to 1e6. In contrast, each library from RNA-Seq dataset has the same size as it corresponds to FPKM.
We can also filter out null and low-variance features, as these might not contain relevant information for differential analysis we plan to conduct.

We start creating tables that summarize quantitative measures per individual and per feature.
When we plot input data for the two datasets, we observe a 0 inflation, which can have two origins:

- biological: 0 expression or methylation, unlikely
- technical: no detection of expression or methylation, likely

Moreover, we observe that distribution is also skewed by high values.

So for better visualization of the latent distribution, we apply log transformation of the data in order to remove 0 and reduce dynamic range of variables.

By doing this, we reveal that biologically related distribution of gene expression and DNA methylation levels at gene promoters, after filtering out low-variance features, both approximate the normal/Gaussian law.

It is interesting to mention that low-variance levels correspond to low methylated regions, and if not filtered out, yield a bimodal distribution.
We choose to calculate Pearson correlation distance between DNA methylated gene promoters, because this method produces the best Rand Adjusted Index ([Jaskowiak *et al*, 2014](#jaskowiak_2014)). Then we apply hierarchical clustering with Ward method ("ward.D2"), as this most robustly identifies groups of individuals with similar quantitative traits. Indeed, this groups individuals according to variance instead of mean, contrary to other methods such as "average" ([Lawlor *et al*, 2016](#lawlor_2016)).

For methylation levels of gene promoters, we observe that the two main clusters segregate day 6 from later developmental stages. However, within the second main cluster, we see that subclusters well separate cell lineages. For subsequent analyses, we will consider analyzing day6 either together or separately from later stages, depending on the question to be addressed.

For gene expression, we observe that clusters nicely segregate cell lineages.

# Dimensionality reduction based analysis of inter-individual vicinity

In order to investigate individual vicinity within sample, we perform a Principal Component Analysis (PCA), which reduces dimensions while keeping maximum of sample variance summarised in "eigenvectors" ([Lever *et al*, 2017](#lever_2017)).

We also perform another dimensionality reduction analysis, the Uniform Manifold Approximation and Projection (UMAP), which is complementary to PCA, as it is nonlinear, and outperforms t-SNE ([Kobak *et al*, 2019](#kobak_2019)). A key parameter for the UMAP is the number of neighbours for each individual. In order to determine the most appropriate number of neighbours according to the latent structure of the data, we use the result of the hierarchical clustering of cells, cut at three clusters recapitulating developmental days or cell types relative to the dataset, and apply the mean number of cells per cluster.
We observe that applying the "correlation" metric based on Pearson distance to umap calculation provides much better results than euclidean distance.

PCA and UMAP yield groups of individuals in line with hierarchical clustering. For both analyses, we observe that gene promoter methylation mainly recapitulates developmental progression (embryonic days) while gene expression recapitulates cell lineages(EPI, PE, TE). This observation is in line with previously published results from ([Zhou *et al*, 2019](#zhou_2019)), as shown in the figure below:

"Principal component analysis 
of DNA methylome data showed that these 130 cells formed 4 major 
clusters (Extended Data Fig. 8g, h), with a combination of the EPI, PE 
and TE at the blastocyst stage (day 6) as a single cluster, and the EPI, 
PE and TE beyond the blastocyst stage as another 3 separate clusters, 
suggesting that all of the 3 lineages showed considerable changes in 
DNA methylation soon after implantation."

In terms of biology, this could mean that at day 6, global DNA methylation state results from epigenetic waves at early embryonic stages, prior to lineage specification. Later stages show rewriting of epigenetic marks specific to cell lineages, which have distinct transcriptomic signatures.

This observation is interesting as a combination of DNA methylation and gene expression levels might allow to assign developmental stage and lineage to single cells of the human embryo.

# WGCNA

After considering inter-individual vicinity, we now investigate correlation within sample at gene promoter methylation and gene expression levels. For this purpose, we perform Weighted Gene Network Correlation Analysis (WGCNA). This algorithm searches for a latent structure within features, i.e modules of correlated genes at promoter methylation or transcriptional levels ([Langfelder *et al*, 2008](#langfelder_2008)).

## Parameter settings & WGCNA computation

First, we need to determine the most suitable WGCNA parameters to apply for each dataset. We notably need to choose the soft-power (soft thresholding power) which is used to power the adjacency matrix of genes, in order to reduce signal noise.

According to [WGCNA package guidelines](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html), it is not recommended filtering data prior to WGCNA. Indeed, this is designed to be an unsupervised analysis method that clusters genes based on their expression profiles (or methylation patterns, by extension). For instance, filtering genes by differential expression would lead to a set of correlated genes that will essentially form a few highly correlated modules. It would also completely invalidate the scale-free topology assumption, so choosing soft thresholding power by scale-free topology fit would fail. Therefore, we provide input data to the WGCNA algorithm.

We choose a soft-power of 12 for both datasets, as this yields a > 0.9 scale-free topology fitting index (R<sup>2</sup>). Intuitively, we would have rather chosen a soft-power of 4 for the scPBAT and 7 for the scRNA-Seq datasets, as these correspond to the inflection points of the curves that also yield a > 0.9 R<sup>2</sup>. However, it is reported in [WGCNA package guidelines](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html) that a minimum soft-power of 12 is advised for signed networks. This soft-power thresholding allows the correlation maximization with scale-free topology producing low mean connectivity.

One pitfall of WGCNA is overfitting the model. This could happen by setting soft-power too high.

Given the high number of genes, we fix a threshold of 100 genes per module, to avoid over-resolution yielding inflated small modules.

WGCNA on gene promoter methylation produces 50 modules, while it yields 39 modules of correlated gene expression. The grey groups correspond to genes that have not been assigned to any module.

## WGCNA calculation
## Visualization of WGCNA networks

Instead of Cytoscape, we use plotly to see it in 3D. What could be the third dimension ? interaction too (?)

## Gene module membership & connectivity
We calculate membership and connectivity of module genes and rename WGCNA modules with eigengenes defined by the highest membership to module.

## Correlation heatmap of WGCNA modules & day of treatment 

We observe that modules of gene promoters and transcripts correlate nicely with development stage or cell lineage, showing almost exclusive patterns. 

# Integration of correlated gene and DNA methylation modules

We now perform an integrative approach on correlation networks between the two datasets, to investigate potential interconnections between genes and DNA methylation  modules.
We select the top most correlated or anti-correlated gene promoter and transcript modules for further investigation. We foucs on 3 pairs of top correlated transcript-pomoter modules, that characterize some lineage-stage traits:

- ME_RPS23-ME_SNORD116.8 (Epi-day6)
- ME_VWCE-ME_SPACA5B (PE-day8)
- ME_UBE2E2-ME_MXI1 (TE-day8)

co-expression similarity and "co-methylation similarity"

## Visualization of WGCNA module network with Cytoscape

## Functional enrichment on WGCNA selected modules

We set a p-value threshold of 0.1 for significant process enrichment, as modules ~ 40 genes fails to produce significant results under a 0.05 threshold.

# Network visualization

As they are the main drivers of cell fate progression, we focus on transcription factor interactions. We query the [TF2DNA database](http://fiserlab.org/tf2dna_db//index.html) on computed interactions in human based on DNA binding motives. Transcription factors are preceded with the "TF" symbol and yield a list of regulated genes they interact with, while other classes of genes do not.
## Mix-Kernel analysis


We want now to go further and use a kernel based approach, which computes a multi-dimension space in which embryo cells will place given their vicinity across gene promoter methylation and gene expression. This exploratory non-supervised analysis traces a non-linear frontier of decision that allocates individuals into groups ([Mariette and Vialaneix, 2018](#mariette_2018)).

We needed a space within cells segregate by developmental stage and lineage progression -> kernel almost succeeded.
Indeed, contrary to single PCA or UMAP, kernel recapitulates both developmental stage progression and lineage specification.

We can use this kernel space to visualize methylation and expression of genes, notably WGCNA module genes, throughout early development of the human embryo.
Clearly, the kernel PCA (KPCA) shows a space in which embryonic cell vicinity recapitulates both lineages and developmental day. This is of utmost interest as we could not combine these two traits with sufficient resolution.

Therefore, the kernel based approach integrates gene promoter methylation with gene expression and computes a space recapitulating developmental progression and lineage specification of the early human embryo.


# PERSPECTIVE : Machine learning -> algorithm that can assign robustly cell lineages and developmental time
We can use this model to build an algorithm that assign a given embryo cell to lineage and developmental stage, provided DNA methylation and gene expression measures.

Machine learning model:

- create
- train
- test (on other datasets too, should work for single dataset and both dataset input)
- use

A predictive model could help to resolve inter-embryo developmental delays or misleaded annotation, and thus accurate developmental stage and lineage maturation identification. This also could help to unify the diverse datasets on a common annotation, thus helping to compare these data.

For this purpose, we could subset the current dataset and use the kernel space and the WGCNA modules, to train the model on one set of cells, and test it on the other set, both of which representative of developmental day and lineage distribution.

However, it would be better that once new datasets are produced, we test our model on these.

# Difficulties

- no access to fastq files
- difficulties to normalise data properly
- results strongly depend on input data distribution
- hard to select appropriate soft-power threshold for WGCNA
- functional enrichment databases seem inaccurate for the human embryo

# Session info

# References

###### Argelaguet, R., Velten, B., Arnol, D., Dietrich, S., Zenz, T., Marioni, J. C., Buettner, F., Huber, W., & Stegle, O. (2018). Multi-Omics Factor Analysis-a framework for unsupervised integration of multi-omics data sets. Molecular systems biology, 14(6), e8124. [https://doi.org/10.15252/msb.20178124](https://doi.org/10.15252/msb.20178124) {#argelaguet_2018}

###### Jaskowiak, P. A., Campello, R. J., & Costa, I. G. (2014). On the selection of appropriate distances for gene expression data clustering. BMC bioinformatics, 15 Suppl 2(Suppl 2), S2. [https://doi.org/10.1186/1471-2105-15-S2-S2](https://doi.org/10.1186/1471-2105-15-S2-S2) {#jaskowiak_2014}


###### Kobak, D., Linderman, G.C. Initialization is critical for preserving global data structure in both t-SNE and UMAP. Nat Biotechnol 39, 156–157 (2021). [https://doi-org.proxy.insermbiblio.inist.fr/10.1038/s41587-020-00809-z](https://doi-org.proxy.insermbiblio.inist.fr/10.1038/s41587-020-00809-z) {#kobak_2021}

###### Langfelder, P., Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008). [https://doi.org/10.1186/1471-2105-9-559](https://doi.org/10.1186/1471-2105-9-559) {#langfelder_2008}

###### Lawlor, N., Fabbri, A., Guan, P., George, J., & Karuturi, R. K. (2016). multiClust: An R-package for Identifying Biologically Relevant Clusters in Cancer Transcriptome Profiles. Cancer informatics, 15, 103–114. [https://doi.org/10.4137/CIN.S38000](https://doi.org/10.4137/CIN.S38000) {#lawlor_2016}

###### Lever, J., Krzywinski, M. & Altman, N. Principal component analysis. Nat Methods 14, 641–642 (2017). [https://doi.org/10.1038/nmeth.4346](https://doi.org/10.1038/nmeth.4346) {#lever_2017}

###### Mariette, J., & Villa-Vialaneix, N. (2018). Unsupervised multiple kernel learning for heterogeneous data integration. Bioinformatics (Oxford, England), 34(6), 1009–1015.[https://doi.org/10.1093/bioinformatics/btx682](https://doi.org/10.1093/bioinformatics/btx682) {#mariette_2018}

###### Molè, M. A., Weberling, A., & Zernicka-Goetz, M. (2020). Comparative analysis of human and mouse development: From zygote to pre-gastrulation. Current topics in developmental biology, 136, 113–138. [https://doi.org/10.1016/bs.ctdb.2019.10.002](https://doi.org/10.1016/bs.ctdb.2019.10.002) {#mole_2020}

###### Zhou, F., Wang, R., Yuan, P., Ren, Y., Mao, Y., Li, R., Lian, Y., Li, J., Wen, L., Yan, L., Qiao, J., & Tang, F. (2019). Reconstituting the transcriptome and DNA methylome landscapes of human implantation. Nature, 572(7771), 660–664. [https://doi.org/10.1038/s41586-019-1500-0](https://doi.org/10.1038/s41586-019-1500-0) {#zhou_2019}

###### Zoppi, J., Guillaume, JF., Neunlist, M. et al. MiBiOmics: an interactive web application for multi-omics data exploration and integration. BMC Bioinformatics 22, 6 (2021). [https://doi.org/10.1186/s12859-020-03921-8](https://doi.org/10.1186/s12859-020-03921-8) {#zoppi_2021}