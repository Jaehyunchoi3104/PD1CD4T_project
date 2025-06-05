
# PySCENIC 분석에 필요한 데이터 유지

#필요한 레이어 설정 (raw count 데이터가 'matrix' 레이어에 있다고 가정)
adata.X = adata.layers["matrix"]

#불필요한 레이어, obsm, varm, obsp 제거
adata.layers = {}
adata.obsm = {}
adata.varm = {}
adata.obsp = {}

#불필요한 obs와 var 컬럼 제거
adata.obs = adata.obs[["n_counts"]]  # 필요한 경우 "cluster" 같은 컬럼 추가
adata.var = adata.var[[]]

#PySCENIC 분석에 필요한 데이터 유지
adata

import numpy as np

#유전자 이름을 row_attrs에 추가
adata.var_names_make_unique()  # 유전자 이름이 중복되지 않도록 처리
adata.var["Gene"] = adata.var_names

#세포 ID를 col_attrs에 추가
adata.obs["CellID"] = adata.obs.index

#Loom 파일로 저장
adata.write_loom("/home/jaehyunchoi/pySCENIC/data/16clustername_241019.loom", write_obsm_varm=False)

import loompy

with loompy.connect("/home/jaehyunchoi/pySCENIC/data/16clustername_241019.loom") as ds:
    print("Row attributes:", ds.ra.keys())  # 'Gene'이 포함되어야 함
    print("Column attributes:", ds.ca.keys())  # 'CellID'가 포함되어야 함


# Docker를 통한 pyscenic 진행
# 1. SCENIC : Run GRN inference
Docker 로 진행하기

docker run --rm -v /home/jaehyunchoi/pySCENIC/data:/data -w /data aertslab/pyscenic:0.12.1
pyscenic grn /data/16clustername_241019.loom
/data/allTFs_mm.txt
-o /data/Byst_SCENIC.adjacencies.tsv
--num_workers 20

# 2. Run regulon prediction
docker run --rm -v /home/jaehyunchoi/pySCENIC/data:/data -w /data aertslab/pyscenic:0.12.1
pyscenic ctx /data/Byst_SCENIC.adjacencies.tsv
/data/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
/data/mm10_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
--annotations_fname /data/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl
--expression_mtx_fname /data/16clustername_241019.loom
--output /data/Byst_SCENIC.motifs.csv
--num_workers 20

# 3. Run cellular enrichment (AU cell)
docker run --rm -v /home/jaehyunchoi/pySCENIC/data:/data -w /data aertslab/pyscenic:0.12.1
pyscenic aucell /data/16clustername_241019.loom
/data/Byst_SCENIC.motifs.csv
--output /data/Byst_SCENIC.auc.csv
--num_workers 20

# -------------------------------------------------------------------------- #

adj_matrix = pd.read_csv('/home/jaehyunchoi/pySCENIC/data/Byst_SCENIC.adjacencies.tsv', index_col=0)
motif_matrix = pd.read_csv('/home/jaehyunchoi/pySCENIC/data/Byst_SCENIC.motifs.csv', index_col=0)
auc_matrix = pd.read_csv('/home/jaehyunchoi/pySCENIC/data/Byst_SCENIC.auc.csv', index_col=0)