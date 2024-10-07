# **Introduction**

Adult diffuse gliomas are the most frequent malignant tumors affecting the central nervous system. Survival rates differ significantly based on subtype. Low-grade gliomas have a high 5-year survival rate, whereas high-grade gliomas are associated with much lower 5-year survival rates. _(Molinaro AM et al. 2019)_

Although histopathologic classification is well-established it is prone to significant intra- and inter-observer variability, especially in grade II-III tumors. _(Whitfield BT et al.,2022)_ 

One key identifier is The isocitrate dehydrogenase family, which comprises three isoforms located in the cytoplasm and peroxisomes (IDH1), and mitochondria (IDH2 and IDH3). These enzymes are involved in cellular processes including mitochondrial oxidative phosphorylation, glutamine metabolism, and more. _(Cairns R et al.,2013)_

Research suggested classifying gliomas into three groups: 

1. IDH wild-type

2. IDH mutant-codel

3. IDH mutant-non-codel

We sought to determine pan-glioma expression subtypes through unsupervised clustering analysis of  RNA-seq profiles by clustering Gliomas Identifies into four clusters labeled LGr1–4. _(Ceccarelli et al.,2016)_

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcd4r5ukJjy-WqjyeXdn_iws4mSBVhG9tQE-WJOF_2n3v7M7MONO_M2wQdzAGRYRv5DFMsMNpW8TAj87Oc78DG5JpsalGsMezikPhKmhR6zQe0Syf9n48jeB9ZHbAgnHCxkhqpKhjLfQ_IUahu23jVQMtK4?key=AgyQUxUvRaYsce19bm41Iw)****

Figure 1: An Overview of the Packages Used for Biomarker Discovery Pipeline Created in BioRender.com

## **Data Preprocessing and Analysis**

TCGA-LGG data is uploaded onto RStudio, then filtered according to IDH status and type. For this analysis we decided to compare between LGr1 & LGr2, and LGr1 & LGr4. Normalization and filtering were performed with a cutoff of 0.25.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeu7Ot6HEVSKTVcFBgfk4uU3Loag7xyqQ6QmxZaBVRiSQNR6Rv0hRkDjoXjIaR8cVgAc2clWsP7fZ1-B3G7PMHA-9ErYFjFVTGmr3DRLU_sVqDtLEnFWITLgM8fLYgveRZd410OjQILhsqgAgt-EfxkRUo?key=AgyQUxUvRaYsce19bm41Iw)

Figure 2: Data subset for this Analysis

## **Differential Expression Analysis**

We performed DEA using edgeR to identify DEGs in glioma samples with significance criteria set at an adjusted log2 fold change (log FC) > 2 and a false discovery rate (FDR) < 0.01.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfKlK8AYdpmw1iZDPqmwCz2kBTOMDALFTeufIgJzyQ8B3l5zGKpZgSoI3DGMY3M2ce9_rUUg4f7R3UVsEuzTzE1VUUrfWMUyNKtcI2LfzKR3V9oo22cBwDrB-bK3zw_xpSPSoAjmbhy7W5_BbyLa-QaA56o?key=AgyQUxUvRaYsce19bm41Iw)

Figure 4: Heatmap showcasing the DEA results of WT glioma in LGr 1 and 4

## **Functional Enrichment Analysis**

We extracted the DEGs for each subset then used biomaRt to convert EnsembleIDs to GeneIDS to perform FEA, and more specifically GO analysis to pinpoint the biological processes, molecular functions, and cellular components linked to these genes.


# **Comparison of results**

Our analysis showed that the up and down regulated pathways align with the pathways described in the paper. For example, some of the down regulated pathways in wildtype gliomas were related to cell adhesion, collagen metabolism, and extracellular matrix remodeling, which is consistent with the poor prognosis and invasive behavior described in the paper. 
Additionally, immune-related pathways are significantly upregulated in mutant gliomas. Both LGr1 and LGr4 analyses show enhanced lymphocyte activation, T-cell activation, and antigen presentation via MHC class II, suggesting active immune responses. While the paper primarily focuses on the less aggressive nature of mutant gliomas, these immunoupregulations hint at complex interactions that affect immune response or tumor suppression.
Also, from our analysis, in corroboration with the findings of the paper by Ceccarelli et al. (2016), we discovered that the comparison between the LGr1 and LGr4 showed that the was the IDH wild type samples contained more LGr4 pan-glioma expression sub-types, while the IDH mutant samples contained more LGr1


# **Conclusion**

In summary, IDH-wildtype gliomas exhibit down regulated pathways linked to invasion and

metastasis, reflecting their poor prognosis, while mutant gliomas show downregulation

of proliferative pathways and upregulation of immune-related processes, consistent with slower growth and better clinical outcomes.

**Citations:**

1. Molinaro AM, Taylor JW, Wiencke JK, Wrensch MR. Genetic and molecular epidemiology of adult diffuse glioma. Nat Rev Neurol. 2019 Jul;15(7):405-417. doi: 10.1038/s41582-019-0220-2. Epub 2019 Jun 21. PMID: 31227792; PMCID:PMC7286557.

2. Whitfield BT, Huse JT. Classification of adult-type diffuse gliomas: Impact of the World Health Organization 2021 update. Brain Pathol. 2022 Jul;32(4):e13062. doi: 10.1111/bpa.13062. Epub 2022 Mar 14. PMID: 35289001;PMCID: PMC9245936.

3. Cairns R, Mak T. Oncogenic isocitrate dehydrogenase mutations: mechanisms, models, and clinical opportunities. Cancer Discovery 2013; 3: 730-41.

4. Ceccarelli et. al, (2016). Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma. Cell, 164(3), 550–563. <https://doi.org/10.1016/j.cell.2015.12.028>


<!--EndFragment-->
