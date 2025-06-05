
# PySCENIC analysis

#pre process adata to analyze with SCENIC
import numpy as np

adata.X = adata.layers["matrix"]

adata.layers = {}
adata.obsm = {}
adata.varm = {}
adata.obsp = {}

adata.obs = adata.obs[["n_counts"]] 
adata.var = adata.var[[]]

adata.var_names_make_unique()  
adata.var["Gene"] = adata.var_names
adata.obs["CellID"] = adata.obs.index
adata.write_loom("/home/jaehyunchoi/pySCENIC/data/16clustername_241019.loom", write_obsm_varm=False)

import loompy

with loompy.connect("/home/jaehyunchoi/pySCENIC/data/16clustername_241019.loom") as ds:
    print("Row attributes:", ds.ra.keys())  # 'Gene'이 포함되어야 함
    print("Column attributes:", ds.ca.keys())  # 'CellID'가 포함되어야 함


# pySCENIC analysis with Docker
# 1. SCENIC : Run GRN inference

```bash
docker run --rm -v /home/jaehyunchoi/pySCENIC/data:/data -w /data aertslab/pyscenic:0.12.1
pyscenic grn /data/16clustername_241019.loom
/data/allTFs_mm.txt
-o /data/Byst_SCENIC.adjacencies.tsv
--num_workers 20
```

# 2. Run regulon prediction
```bash
docker run --rm -v /home/jaehyunchoi/pySCENIC/data:/data -w /data aertslab/pyscenic:0.12.1
pyscenic ctx /data/Byst_SCENIC.adjacencies.tsv
/data/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
/data/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
--annotations_fname /data/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
--expression_mtx_fname /data/16clustername_241019.loom
--output /data/Byst_SCENIC.motifs.csv
--num_workers 20
```
# 3. Run cellular enrichment (AU cell)
```bash
docker run --rm -v /home/jaehyunchoi/pySCENIC/data:/data -w /data aertslab/pyscenic:0.12.1
pyscenic aucell /data/16clustername_241019.loom
/data/Byst_SCENIC.motifs.csv
--output /data/Byst_SCENIC.auc.csv
--num_workers 20
```
# -------------------------------------------------------------------------- #

adj_matrix = pd.read_csv('/home/jaehyunchoi/pySCENIC/data/Byst_SCENIC.adjacencies.tsv', index_col=0)
motif_matrix = pd.read_csv('/home/jaehyunchoi/pySCENIC/data/Byst_SCENIC.motifs.csv', index_col=0)
auc_matrix = pd.read_csv('/home/jaehyunchoi/pySCENIC/data/Byst_SCENIC.auc.csv', index_col=0)
