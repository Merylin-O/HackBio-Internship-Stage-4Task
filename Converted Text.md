Introduction\
Adult diffuse gliomas are the most frequent malignant tumors affecting the central\
nervous system. Survival rates differ significantly based on the glioma sub-type.\
Low-grade gliomas have a relatively high 5-year survival rate of up to 80%,\
whereas high-grade gliomas are associated with much lower 5-year survival rates (Molinaro AM et al.,2019).\
Although histopathologic classification is well-established and forms the basis of\
the World Health Organization (WHO) classification of CNS tumors it is prone to\
significant intra- and inter-observer variability, especially in grade II-III tumors (Whitfield BT et al.,2022).\
Mutations in the isocitrate dehydrogenase genes 1 and 2 (IDH1/IDH2) identify a\
distinct subset of glioblastoma. The isocitrate dehydrogenase (IDH) family of\
enzymes comprises three isoforms located in the cytoplasm and peroxysomes\
(IDH1), and mitochondria (IDH2 and IDH3). IDH enzymes are involved in a\
number of cellular processes including mitochondrial oxidative phosphorylation,\
glutamine metabolism, lipogenesis, glucose sensing, and regulation of cellular\
redox status (Cairns R et al.,2013).\
Dataset comprising all TCGA newly diagnosed diffuse glioma consisting of 1,122\
patients .Research has suggested classifying gliomas into three groups:\
·1. IDH wild-type\
·2. IDH mutant with 1p/19q codeletion (IDH mutant-codel)\
·3. IDH mutant without 1p/19q codeletion (IDH mutant-non-codel)\
they sought to determine pan-glioma expression sub-types through unsupervised\
clustering analysis of RNA-seq profiles by Unsupervised Clustering of Gliomas\
Identifies four main clusters labeled LGr1–4 ((Ceccarelli et al.,2016)

Data Preprocessing and Analysis\
This preprocessing step aimed to prepare the dataset for downstream analysis,\
ensuring data quality and reliability. We retrieved a dataset of LGG gliomas\
samples from The Cancer Genome Atlas (TCGA) website

Methodology\
Our analysis employed the TCGAbiolinks pipeline, leveraging various R packages:\
1\. We uploaded TCGA-LGG project data.\
2\. filter the data according to IDH satus : Mutant and Wild Type Comparison\
between two expression clusters in IDH status was done : (LGr1 and LGr4)\
3\. Normalization and filtering were performed using TCGAbiolinks libraries: -\
TCGAanalyze\_Normalization, TCGAanalyze\_Filtering Quantile method\
with a cutoff of 0.25.

Differential Expression Analysis\
Following data preprocessing, we performed differential expression analysis using\
the edgeR package to identify differentially expressed genes in glioma samples.\
We compared results between two IDH status groups—wild-type and mutant—\
across two clusters (LGr1 vs LGr4). The analysis was conducted using the\
TCGAanalyze\_DEA function, with significance criteria set at an adjusted log2 fold\
change (log FC) > 2 and a false discovery rate (FDR) < 0.01.

Functional Enrichment Analysis\
We then used the biomaRt library to map Ensembl IDs to gene IDs through the\
getBM function, which enabled us to carry out functional enrichment analysis.\
Employing the TCGAanalyze\_EAcomplete function, we performed Gene\
Ontology (GO) analysis on the differentially expressed genes. This helped us\
pinpoint the biological processes, molecular functions, and cellular components\
linked to these genes.

Comparison of results\
From our analysis, in corroboration with the findings of the paper by Ceccarelli _et al. _(2016), we discovered that the comparison between the LGr1 and LGr4 showed that the was the IDH wild type samples contained more LGr4 pan-glioma expression sub-types, while the IDH mutant samples contained more LGr1.

In IDH-wildtype gliomas, downregulation of pathways related to cell adhesion, collagen metabolism, and extracellular matrix (ECM) remodeling is observed. This is consistent with the aggressive and metastatic nature of these tumors, as the disruption of cell adhesion and\
ECM regulation promotes invasion. Additionally, the downregulation of colorectal cancer metastasis signaling and leukocyte extravasation signaling suggests a role in tumor migration, which aligns with the poor prognosis and invasive behavior of IDH-wildtype gliomas described in the paper. This molecular profile of aggressive tumor progression mirrors glioblastoma, further supporting the poor prognosis associated with IDH-wildtype gliomas.\
Conversely, IDH-mutant gliomas show a downregulation of cell cycle, mitotic, and nuclear division pathways, pointing to slower tumor proliferation. This downregulation is in line with the less aggressive nature of these tumors, which are associated with better clinical outcomes due to slower growth and epigenetic silencing of proliferative pathways via G-CIMP hypermethylation. The IDH-mutant gliomas also demonstrate a more stable genomic profile, which contributes to their better prognosis. These findings corroborate the paper’s portrayal of IDH-mutant gliomas as slower-growing tumors with reduced invasiveness.\
Interestingly, immune-related pathways are significantly upregulated in IDH-mutant gliomas. Both LGr1 and LGr4 analyses show enhanced lymphocyte activation, T-cell activation, and antigen presentation via MHC class II, suggesting active immune responses. Pathways such as dendritic cell maturation and T-helper cell differentiation are also enriched, indicating a role in immune surveillance and activation. While the paper primarily focuses on the less aggressive nature of IDH-mutant gliomas, these immune-related upregulations hint at a complex interaction between the tumor and its microenvironment, potentially affecting immune response or tumor suppression.

Conclusion

In summary, IDH-wildtype gliomas exhibit downregulated pathways linked to invasion and metastasis, reflecting their poor prognosis. In contrast, IDH-mutant gliomas show downregulation of proliferative pathways and upregulation of immune-related processes, consistent with slower growth and better clinical outcomes.\
There is potential for more clusters based on newer datasets. The upregulated results show significant activation of immune-related pathways (e.g., T-cell activation, antigen presentation, dendritic cell maturation). These findings suggest that immune activity could be a distinguishing feature in IDH-mutant gliomas and may warrant further clustering based on the immune microenvironment or tumor-infiltrating lymphocytes.

References\
Whitfield BT, Huse JT. Classification of adult-type diffuse gliomas: Impact of the World Health Organization 2021 update. Brain Pathol. 2022 Jul;32(4):e13062. doi:\
10.1111/bpa.13062. Epub 2022 Mar 14. PMID: 35289001; PMCID: PMC9245936.\
Molinaro AM, Taylor JW, Wiencke JK, Wrensch MR. Genetic and molecular epidemiology of adult diffuse glioma. Nat RevNeurol. 2019 Jul;15(7):405-417. doi: 10.1038/s41582-019-0220-2. Epub 2019 Jun 21. PMID: 31227792; PMCID: PMC7286557.\
Cairns R, Mak T. Oncogenic isocitrate dehydrogenase mutations: mechanisms, models, and clinical opportunities. Cancer Discovery 2013; 3: 730-41.\
Ceccarelli, M., Barthel, F. P., Malta, T. M., Sabedot, T. S., Salama, S. R., Murray, B. A., Morozova, O., Newton, Y., Radenbaugh, A., Pagnotta, S. M., Anjum, S., Wang, J., Manyam, G., Zoppoli, P., Ling, S., Rao, A. A., Grifford, M., Cherniack, A. D., Zhang, H., . . . Zmuda, E. (2016). Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma. Cell, 164(3), 550–563. https\://doi.org/10.1016/j.cell.2015.12.028\
Louis, D. N., Ohgaki, H., Wiestler, O. D., Cavenee, W. K., Burger, P. C., Jouvet, A., Scheithauer, B. W., & Kleihues, P. (2007). The 2007 WHO Classification of Tumours of the Central Nervous System. Acta Neuropathologica, 114(2),\
97–109. https\://doi.org/10.1007/s00401-007-0243-4

 

 
