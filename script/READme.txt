[Project Name] Byst_IPF_proj
[Paper Name] Type 2 innate-like pathogenic function of PD-1high CD4 T cells aggravate pulmonary fibrosis
[Created] 2025-04-17
[Owner] Jaehyun Choi

[Data]
- raw_data/: GSE295241
- processed: filtered and normalized with Scanpy and R,  see 01_preprocessing.ipynb

[Main Steps]
1. Preprocessing: 01_preprocessing.ipynb [conda env test02]
2. Clustering: 02_clustering.ipynb [conda env test02]
3. DEG Analysis: 03_DEG_analysis.md [R] 
GO and KEGG analyses were performed using differentially expressed genes (log₂ fold change > 0.58, adjusted P-value < 0.05) through the DAVID functional annotation tool
4. Velocity analysis: 04_deepvelo.ipynb [conda env deepvelo and scvelo_backup]
5. psuedotime and Fate probabilities analysis: 05_fate_anlaysis.ipynb [conda env cellrank2]
6. SCENIC_TF interaction analysis : 06_pySCENIC.md and 07_SCENIC_analysis.ipynb [conda env test02]


[Note]
