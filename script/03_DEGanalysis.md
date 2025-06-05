# Load R package
```{r}
library(Seurat)
library(MAST)
library(dittoSeq)
library(caret)
library(ranger)
library(feseR)
library(devtools)
library(pROC)
library(reticulate)
library(SingleR)
library(slingshot)
library(celldex)
library(scCATCH)
library(dplyr)
```
# Cytokine differentiation DEG analysis
```{r}
Cord_diff = readRDS("/home/jaehyunchoi/bystander_collaboration/Figure_script_IPF/Data/Cord_diff.rds")
table(Cord_diff@meta.data$celltype)
  #Resting  Th1-like  Th2-like Treg-like Undefined 
  #11183      2421      4051       350       885 

 Idents(Cord_diff)<-"celltype"
 Diff_Th2_markers <- FindAllMarkers(Cord_diff,min.pct = 0.1,logfc.threshold = 0.25, only.pos = TRUE)
 Diff_Th2_DEG <- FindMarkers(Cord_diff, ident.1 = "Th2-like", only.pos = FALSE)
 write.csv(Diff_Th2_markers, file = "/home/jaehyunchoibystander_collaboration/Figure_script_IPF/Results/Diff_Th2_markers_fam.csv", row.names = TRUE)
 write.csv(Diff_Th2_DEG, file = "/home/jaehyunchoi/bystander_collaboration/Figure_script_IPF/Results/Diff_Th2_markers_fm.csv", row.names = TRUE)

EnhancedVolcano(Diff_Th2_DEG,
                lab = rownames(Diff_Th2_DEG),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = c("Areg","Il5","Il13","Batf","Junb","Gata3","Lgals1","Epas1","Il1rl1"),  
                max.overlaps = 10,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                labSize = 5.0,
                labFace = 'bold',
                drawConnectors = TRUE,
                xlim = c(-1.5, 2.5)
)
```
```{r}
# MPEC CD4 T cells include DEG analysis
Cord_byst = readRDS("/home/jaehyunchoi/bystander_collaboration/Figure_script_IPF/Data/Cord_byst.rds")
table(Cord_byst@meta.data$celltype)

 #CCR6_high IFN-I_high  PD-1_high    Resting   Th1-like   Th2-like  Treg-like  Undefined 
 #294        123        487      14281       3253       4051        755        885 
 
Idents(Cord_byst)<-"celltype"
Byst_Th2_markers <- FindAllMarkers(Cord_byst,min.pct = 0.1,logfc.threshold = 0.25, only.pos = TRUE)
Byst_Th2_DEG <- FindMarkers(Cord_byst, ident.1 = "Th2-like", only.pos = FALSE)
 write.csv(Byst_Th2_markers, file = "/home/jaehyunchoi/bystander_collaboration/Figure_script_IPF/Results/Byst_Th2_markers_fam.csv", row.names = TRUE)
 write.csv(Byst_Th2_DEG, file = "/home/jaehyunchoi/bystander_collaboration/Figure_script_IPF/Results/Byst_Th2_markers_fm.csv", row.names = TRUE)

EnhancedVolcano(Byst_Th2_DEG,
                lab = rownames(Byst_Th2_DEG),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                selectLab = c("Areg","Il5","Il13","Batf","Junb","Gata3","Lgals1","Epas1","Il1rl1"),  
                max.overlaps = 10,
                pCutoff = 0.05,
                FCcutoff = 0.58,
                labFace = 'bold',
                labSize = 5.0,
                drawConnectors = TRUE,
                xlim = c(-1.5, 2.5)
)

Cord_byst@meta.data$Th2vsother <- ifelse(Idents(Cord_byst) == "Th2-like", "Th2-like", "Other")
Idents(Cord_byst) <- "Th2vsother"

p <- DotPlot(Cord_byst, features = Th2_gene, group.by = "Th2vsother")
p$data$id <- factor(p$data$id, levels = c("Th2-like", "Other"))

p + 
  geom_point(
    data = p$data,
    aes(x = features.plot, y = id, size = pct.exp, fill = avg.exp.scaled),
    shape = 21,         # 내부 채움과 테두리 색 구분 가능
    color = "black",    # 테두리 검정색
    stroke = 0.5        # 테두리 두께
  ) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    #limits = c(-0.4, 0.4),
    oob = scales::squish,
    name = "Scaled\nExpression"
  ) +
  scale_size(range = c(1, 8), name = "% Expressed") +
  coord_flip() +
  theme_minimal() +
  theme(text = element_text(size = 12))

dotdata <- DotPlot(Cord_byst, features = Th2_gene, group.by = "Th2vsother")$data

dotdata$id <- factor(dotdata$id, levels = c("Th2-like", "Other"))

ggplot(dotdata, aes(x = features.plot, y = id)) +
    geom_point(
        aes(size = pct.exp, fill = avg.exp.scaled),
        shape = 21,
        color = "black",
        stroke = 0.5
    ) +
    scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0,
        oob = scales::squish,
        name = "Scaled Exp"
    ) +
    scale_size(range = c(1, 8), name = "% Cells") +
    coord_flip() +
    theme_minimal() +
    theme(
        text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)
    )
```
# Marker gene comparison between Tibbit et al., Khan et al., pTh2 cell groups
```{r}
library(fgsea)
 #1. Tibbit et al.,
gene_string <- "Igfbp7 Il13 Il1rl1 Plac8 Bhlhe40 Gata3 Nfkb1 Rbpj Fgl2 Ltb4r1 Epas1 Il5 Zcchc10 Hlf Mns1 Pparg Cxcr6 Tagln2 Vdr Gclc Pdcd1 Rnf128 Il6 Paxbp1 Id2 4932438H23Rik Mboat7 Il4 Litaf Il17rb Csf1 Adck5 Ccr8 Rgs1 Serpinb6a Lgals1 Cd200r1 Dusp6 Sh3bgrl3 Ldha Sorl1 Rgcc Ccdc50 Tmsb4x Cd44 Samsn1 Ctla2a Cotl1 Lgals3 Tnfsf11 Tusc3 S100a11 Icosl Ramp1 Sdf2 Dync1li1 Rgs2 Aff4 Aprt Cd3g Prnp Bcl2a1d Lig1 Plp2 Impa2 Degs1 L1cam Cmpk1 Acvr2a Neurl3 Vipr2 Pdia3 Plk3 Lgals3bp Ehbp1l1 Rb1 Fam213a Atpif1 Arl1 Cib2 Adam8 4930404I05Rik Baiap2 Dennd5a Uhrf2 Ech1 Pcmt1 Lztfl1 Tmem64 Abhd14a Junb Padi2 Lman1 Slc25a4 Btg1 Il9r Lmna Ell2 Etfb Hist1h1c Pld3 D630039A03Rik Gsn Arg1 1190002N15Rik Mrpl20 Nrp1 Ccnd2 Endod1 Zyx Gm17745 Capg Gm5127 Slmo2 Psmc2 Serinc3 Ptms Ttc3 Ero1l Cyfip1 Cyp51 Galnt1 Tpi1 Cysltr1 2700089E24Rik Il18rap Rtca Cd52 Cish Ube2l6 Osbpl5 Inpp5k Gna15 Orc6 AA467197 Abi1 Msmo1 Snx12 Kdelr1 Cd40lg Crip1 Fdx1 Rexo2 Ndfip1 Jak2 Psma1 Crip2 Mrpl33 Rnf135 4930503L19Rik P4hb Plin2 Egr1 Gmds Lyn Cd8b1 Stx12 Tspan31 F2r Ccdc109b Cnih1 Ly6g5b Prelid2 Suox Pmpcb Mvk Sqle Gnb4 Aamdc Pros1 Cox19 Blcap Slc25a17 Nagk 1700021F05Rik Psat1 Pradc1 Drg2 Snrnp25 Spcs3 Crcp Tm2d2 Acp6 Nt5c3 Eef1e1 Cept1 Myadm Farsa Arl4a Zcchc9 Gdap2 Igsf8 Tmem238 Rnf11 Far1 Rbbp7 Crem Rbfa Klf10 Chst12 Cass4 Pdcl Coq2 Lxn Nceh1 Raly Tspan13 Fos Glipr1 Slc52a2 Bad Hmgcs1 Ttc39c Gpr171 Guk1 Fam69a Ptcd3 Trp53i13 Anxa5 Rtcb Acaa2 Fam179b Mrpl36 Cd82 Ndufv1 Cers5 Fundc2 Cetn2 Ash2l 9530068E07Rik Elovl5 Pop5 Scamp2 Rtn3 Tsc22d3 Tsg101 Cd69 Skap1 Sept7 Tmem128 Hmgcr Vapa Hspa4 Tmem30a Tnfaip3 Hnrnpab Zap70 Hspa5"

genes <- strsplit(gene_string, "\\s+")[[1]]
Tib_th2 <- genes
gene_list <- Byst_Th2_DEG$avg_log2FC
names(gene_list)<-rownames(Byst_Th2_DEG)
gene_list <- sort(gene_list, decreasing = TRUE)
pathways_list <- list(peTh2_signature = Tib_th2)

fgseaRes <- fgsea(pathways = pathways_list,
                  stats = gene_list)

plotEnrichment(pathways_list$peTh2_signature, gene_list) +
    ggtitle("GSEA: peTh2 Signature")
#NES와 FDR 값을 추출
my_nes <- round(fgseaRes$NES[fgseaRes$pathway == "peTh2_signature"], 2)
my_fdr <- signif(fgseaRes$padj[fgseaRes$pathway == "peTh2_signature"], 2)
#유전자 포함 여부 계산
n_total <- length(geneset_genes)
n_found <- sum(geneset_genes %in% names(gene_list))


geneset_name <- "peTh2_signature"
geneset_genes <- pathways_list[[geneset_name]]
annotation_text <- paste0("NES = ", my_nes,
                          "\nP Value = ", my_fdr,
                          "\n", n_found, "/", n_total, " genes found")
#enrichment plot + NES/FDR annotation
plotEnrichment(pathways_list$peTh2_signature, gene_list) +
    ggtitle("GSEA: peTh2 Signature") +
    annotate("text",
             x = Inf, y = Inf,
             label = paste0("NES = ", my_nes, "\nP Value = ", my_fdr),
             hjust = 1.1, vjust = 1.5, size = 4)

plotEnrichment(geneset_genes, gene_list) +
    ggtitle(paste("Airway pTh2 gene set(Tibbitt et al.)")) +
    annotate("text",
             x = Inf, y = Inf,
             label = annotation_text,
             hjust = 1.1, vjust = 1.5,
             size = 4
             )


 #2. Khan et al.,
Khan_Th2 <- read.csv("/home/jaehyunchoi/bystander_collaboration/Figure_script_IPF/Results/Khanetal_Th2.csv")
Khan_Th2 <- Khan_Th2$gene
pathways_list <- list(peTh2_signature = Khan_Th2)
gene_list <- Byst_Th2_DEG$avg_log2FC
names(gene_list)<-rownames(Byst_Th2_DEG)
gene_list <- sort(gene_list, decreasing = TRUE)
fgseaRes <- fgsea(pathways = pathways_list,
                  stats = gene_list)
plotEnrichment(pathways_list$peTh2_signature, gene_list) +
    ggtitle("GSEA: peTh2 Signature")

#NES와 FDR 값을 추출
my_nes <- round(fgseaRes$NES[fgseaRes$pathway == "peTh2_signature"], 2)
my_fdr <- signif(fgseaRes$padj[fgseaRes$pathway == "peTh2_signature"], 2)
#필요한 값 추출
geneset_name <- "peTh2_signature"
geneset_genes <- pathways_list[[geneset_name]]
annotation_text <- paste0("NES = ", my_nes,
                          "\nP Value = ", my_fdr,
                          "\n", n_found, "/", n_total, " genes found")
#유전자 포함 여부 계산
n_total <- length(geneset_genes)
n_found <- sum(geneset_genes %in% names(gene_list))
#enrichment plot + NES/FDR annotation
plotEnrichment(pathways_list$peTh2_signature, gene_list) +
    ggtitle("GSEA: peTh2 Signature") +
    annotate("text",
             x = Inf, y = Inf,
             label = paste0("NES = ", my_nes, "\nFDR = ", my_fdr),
             hjust = 1.1, vjust = 1.5, size = 4)

plotEnrichment(geneset_genes, gene_list) +
    ggtitle(paste("Airway pTh2 gene set(Khan et al.)")) +
    annotate("text",
             x = Inf, y = Inf,
             label = annotation_text,
             hjust = 1.1, vjust = 1.5,
             size = 4
             )
```
