
## non-kinetics

### data

- Constructs: GFP_RBCC, GFP_VHH, HuR_RBCC, HuR_VHH, RBCC_GFP, RBCC_HuR
- Plates: R1, R2, R3

Comparisons to make:
- RBCC_HuR vs HuR_VHH
- RBCC_HuR vs RBCC_GFP
- RBCC_HuR vs GFP_VHH
- HuR_VHH vs RBCC_GFP
- RBCC_GFP vs HuR_VHH
- HuR_VHH vs GFP_VHH
- RBCC_GFP vs GFP_VHH

Name keys for manuscript:
- RBCC-HuR = TRIM21-VHHHuR
- HuR-RBCC = VHHHuR-TRIM21
- RBCC-GFP = TRIM21-VHHGFP
- GFP-RBCC = VHHGFP-TRIM21
- HuR VHH = VHHHuR
- GFP VHH = VHHGFP



### limma

```{r}
library(devtools)
load_all("/projects/qbio/bifo/R/3.6.0/pxanalytics")
library(tidyverse)
library(yaml)
library(data.table)
library(limma)
library(ggplot2)
library(ggrepel)


# Enlarge the view width when printing tables
options(width = 300)


# Load data
data <- fread("20210414_214732_DOA078_Protac_IVEB_Trim21_041421_ReportProteins.txt")
data <- unique(data, by = "genes")
data <- data.frame(data)
rownames(data) <- data$genes


# Assay data
assay <- na.omit(data[,2:37])
nrow(data) # 6455


# Annotations
anno <- data.frame("ID" = rownames(assay), stringsAsFactors = FALSE)


# Metadata
metadata <- data.frame("ID" = colnames(assay), "Sample" = sapply(colnames(assay), function (x) paste(unlist(strsplit(x, "_"))[1:2], collapse="_")), "Plate" = sapply(colnames(assay), function (x) unlist(strsplit(x, "_"))[3]), "Replicate" = sapply(colnames(assay), function (x) unlist(strsplit(x, "_"))[4]), stringsAsFactors = FALSE)


# pxdata object initialisation
px_raw <- pxinit(assay, annotations = anno, metadata = metadata)


# Normalisation
px_norm <- px_raw %>% pxnormalise(method = "total", batch_var = "Plate")


# model.matrix, lmFit, makeContrasts, contrasts.fit, eBayes
des <- model.matrix(~ 0 + metadata$Sample)
colnames(des) <- levels(factor(metadata$Sample))

fit <- lmFit(pxdata_logged(px_norm), des)

contrast.matrix <- makeContrasts(
  "RBCC_HuR_vs_HuR_VHH" = RBCC_HuR-HuR_VHH,
  "RBCC_HuR_vs_RBCC_GFP" = RBCC_HuR-RBCC_GFP,
  "RBCC_HuR_vs_GFP_VHH" = RBCC_HuR-GFP_VHH,
  "HuR_VHH_vs_RBCC_GFP" = HuR_VHH-RBCC_GFP,
  "RBCC_GFP_vs_HuR_VHH" = RBCC_GFP-HuR_VHH,
  "HuR_VHH_vs_GFP_VHH" = HuR_VHH-GFP_VHH,
  "RBCC_GFP_vs_GFP_VHH" = RBCC_GFP-GFP_VHH,
  levels=des)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)


# biogrid data
biogrid <- fread("BioGRID_List_of_Proteins.csv")
biogrid_vector <- setdiff(unique(c(biogrid[,8][[1]], biogrid[,9][[1]])), c("ELAVL1"))
biogrid_vector <- append(biogrid_vector, c("SGO1", "BRF1", "NR2F2", "GOT1", "FH", "SAT2", ";TUBB3"))


#######################
# RBCC_HuR_vs_HuR_VHH #
#######################
detable <- data.table(topTable(fit2, coef="RBCC_HuR_vs_HuR_VHH", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn    logFC  AveExpr         t      P.Value    adj.P.Val        B
#1: ELAVL1 -2.81906 9.491571 -42.17192 1.325216e-29 8.554269e-26 54.10872

detable[rn == "ELAVL2"]
#       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
#1: ELAVL2 -2.639226 10.96142 -38.65639 2.034756e-28 6.567175e-25 52.04169

nrow(detable[adj.P.Val < 0.05]) # 2353
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 111
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 9

detable[logFC < -1 & adj.P.Val < 0.05]
#          rn     logFC   AveExpr          t      P.Value    adj.P.Val         B
#  1:  ELAVL1 -2.819060  9.491571 -42.171918 1.325216e-29 8.554269e-26 54.108716
#  2:  ELAVL2 -2.639226 10.961423 -38.656385 2.034756e-28 6.567175e-25 52.041693
#  3:   TFAP4 -2.483983  7.111919 -34.210852 9.251127e-27 1.990534e-23 49.004644
#  4: IGF2BP3 -1.227526 10.038029 -29.658186 7.770201e-25 1.253916e-21 45.288042
#  5:   YLPM1 -1.992173  8.225543 -24.960952 1.543015e-22 1.992033e-19 40.624987
# ---                                                                           
#107:  MINDY4 -1.211854  4.227627  -2.733921 1.011545e-02 3.148275e-02 -3.645650
#108:    NANP -1.370740  7.769224  -2.704756 1.086325e-02 3.337566e-02 -3.711514
#109:     SF1 -1.184486  6.878176  -2.560893 1.536279e-02 4.407414e-02 -4.029808
#110: B4GALT5 -1.047678  6.906644  -2.495547 1.792822e-02 4.945584e-02 -4.170623
#111: SLC35B3 -1.589751  7.209701  -2.494121 1.798835e-02 4.950969e-02 -4.173668

detable[logFC > 1 & adj.P.Val < 0.05]
#        rn    logFC   AveExpr         t      P.Value    adj.P.Val          B
#1: CCDC171 5.048946  8.859662 20.559459 5.349100e-20 5.754740e-17 35.2463530
#2:    NME3 1.149213  8.651534 12.989454 2.672913e-14 1.568514e-11 22.6135618
#3:   PCBP4 1.094157 11.197335  8.195913 2.328513e-09 1.330137e-07 11.3387095
#4:  FYTTD1 1.077030  5.745152  4.666012 5.239546e-05 4.551987e-04  1.4026828
#5:  CALML5 1.301974  6.134296  3.875820 4.962189e-04 2.841087e-03 -0.7860154
#6:    KRT9 1.154140 10.239632  3.282377 2.494200e-03 1.021578e-02 -2.3329412
#7:  OXNAD1 1.166505  6.943503  2.985357 5.394426e-03 1.914295e-02 -3.0604369
#8: TSPAN13 1.207254  5.689078  2.974545 5.544950e-03 1.956952e-02 -3.0862140
#9:  SYNGR3 2.174762  6.738369  2.926179 6.268142e-03 2.155613e-02 -3.2008696


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_HuR_vs_HuR_VHH.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_HuR_vs_HuR_VHH.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_HuR_vs_HuR_VHH.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_HuR_vs_HuR_VHH.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR vs VHHHuR") +
theme_bw() +
coord_cartesian(xlim = c(-5, 5), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_RBCC_HuR_vs_HuR_VHH.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR vs VHHHuR") +
theme_bw() +
coord_cartesian(xlim = c(-5, 5), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5))

ggsave('figures/volcano_RBCC_HuR_vs_HuR_VHH_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)


########################
# RBCC_HuR_vs_RBCC_GFP #
########################
detable <- data.table(topTable(fit2, coef="RBCC_HuR_vs_RBCC_GFP", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn     logFC  AveExpr         t      P.Value    adj.P.Val       B
#1: ELAVL1 -2.895363 9.491571 -43.31337 5.722468e-30 3.693853e-26 51.9705

detable[rn == "ELAVL2"]
#       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
#1: ELAVL2 -2.848272 10.96142 -41.71825 1.861626e-29 6.008398e-26 51.24771

nrow(detable[adj.P.Val < 0.05]) # 836
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 13
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 6

detable[logFC < -1 & adj.P.Val < 0.05]
#            rn     logFC   AveExpr          t      P.Value    adj.P.Val         B
# 1:     ELAVL1 -2.895363  9.491571 -43.313372 5.722468e-30 3.693853e-26 51.970499
# 2:     ELAVL2 -2.848272 10.961423 -41.718248 1.861626e-29 6.008398e-26 51.247713
# 3:    IGF2BP3 -1.165338 10.038029 -28.155670 3.860022e-24 8.305481e-21 42.563760
# 4:      TFAP4 -1.801600  7.111919 -24.812675 1.849997e-22 2.985433e-19 39.429406
# 5:        DSP -1.115766 10.065480  -8.696967 6.178842e-10 2.658961e-07 12.740353
# 6:      PEG10 -1.109217  6.005135  -6.609307 1.889217e-07 2.540603e-05  7.120344
# 7:       SDC4 -1.517687  6.439895  -5.715496 2.483256e-06 1.908264e-04  4.591440
# 8:     RNF115 -1.215109  5.974537  -5.165758 1.229266e-05 6.399122e-04  3.026155
# 9:     MRFAP1 -1.056535  6.494097  -4.732994 4.317329e-05 1.671465e-03  1.801501
#10:    CCDC171 -1.015862  8.859662  -4.136622 2.384651e-04 5.786814e-03  0.145861
#11: ;MCM8;MCM8 -1.020628  6.331991  -4.010836 3.400418e-04 7.329533e-03 -0.195861
#12:     MAGEA2 -2.252860  8.581502  -3.648775 9.294773e-04 1.383187e-02 -1.159117
#13:    B4GALT5 -1.230316  6.906644  -2.930587 6.198727e-03 4.855920e-02 -2.947455

detable[logFC > 1 & adj.P.Val < 0.05]
#        rn    logFC   AveExpr        t      P.Value   adj.P.Val          B
#1:  CEP350 1.395232  6.617072 5.975319 1.168920e-06 0.000106273  5.3304385
#2:    CTSL 1.066260  7.183053 4.792707 3.632086e-05 0.001482888  1.9696693
#3:   CREB1 1.597339  6.832917 4.314020 1.440135e-04 0.004095187  0.6328717
#4:    KRT9 1.243430 10.239632 3.536318 1.263020e-03 0.016706552 -1.4510801
#5: CNEP1R1 1.107842  5.592712 3.396084 1.843109e-03 0.021245127 -1.8095806
#6:   ACTA1 1.325237 12.416964 2.920820 6.353531e-03 0.049360775 -2.9703556


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_HuR_vs_RBCC_GFP.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_HuR_vs_RBCC_GFP.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_HuR_vs_RBCC_GFP.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_HuR_vs_RBCC_GFP.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR vs TRIM21-VHHGFP") +
theme_bw() +
coord_cartesian(xlim = c(-5, 5), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_RBCC_HuR_vs_RBCC_GFP.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR vs TRIM21-VHHGFP") +
theme_bw() +
coord_cartesian(xlim = c(-5, 5), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5))

ggsave('figures/volcano_RBCC_HuR_vs_RBCC_GFP_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)


#######################
# RBCC_HuR_vs_GFP_VHH #
#######################
detable <- data.table(topTable(fit2, coef="RBCC_HuR_vs_GFP_VHH", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
#1: ELAVL1 -2.846842 9.491571 -42.58753 9.736491e-30 6.284905e-26 54.36229

detable[rn == "ELAVL2"]
#       rn     logFC  AveExpr        t      P.Value    adj.P.Val        B
#1: ELAVL2 -2.435586 10.96142 -35.6737 2.506799e-27 8.090694e-24 50.08132

nrow(detable[adj.P.Val < 0.05]) # 2615
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 113
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 11

detable[logFC < -1 & adj.P.Val < 0.05]
#          rn     logFC   AveExpr          t      P.Value    adj.P.Val         B
#  1:  ELAVL1 -2.846842  9.491571 -42.587532 9.736491e-30 6.284905e-26 54.362286
#  2:  ELAVL2 -2.435586 10.961423 -35.673697 2.506799e-27 8.090694e-24 50.081321
#  3:   TFAP4 -2.552946  7.111919 -35.160641 3.939557e-27 8.476614e-24 49.716845
#  4: IGF2BP3 -1.196532 10.038029 -28.909339 1.710780e-24 2.760771e-21 44.619662
#  5:   YLPM1 -1.937996  8.225543 -24.282147 3.569808e-22 4.608622e-19 39.875260
# ---                                                                           
#109:   KLF16 -1.064159  4.355561  -2.753708 9.635734e-03 2.774249e-02 -3.603721
#110:     SF1 -1.222862  6.878176  -2.643863 1.259317e-02 3.441528e-02 -3.850614
#111:    BRF1 -1.210674  6.116357  -2.579196 1.470742e-02 3.890744e-02 -3.992962
#112:   CDK17 -1.891126  7.386061  -2.557365 1.549219e-02 4.053590e-02 -4.040501
#113: TGFB1I1 -2.197078  6.490737  -2.537615 1.623517e-02 4.205377e-02 -4.083277

detable[logFC > 1 & adj.P.Val < 0.05]
#         rn    logFC   AveExpr         t      P.Value    adj.P.Val          B
# 1: CCDC171 5.323255  8.859662 21.676452 1.099794e-20 1.183195e-17 36.7273238
# 2:   AIF1L 1.111773  6.328864  4.419099 1.066286e-04 7.207197e-04  0.7049697
# 3:    KRT9 1.379875 10.239632  3.924369 4.333262e-04 2.248489e-03 -0.6578683
# 4:    MT1X 1.291303  8.019070  3.783278 6.416940e-04 3.112047e-03 -1.0369098
# 5:  SYNGR3 2.729961  6.738369  3.673208 8.692332e-04 3.943008e-03 -1.3289810
# 6:    MT1E 1.345569  9.107433  3.644832 9.395727e-04 4.220558e-03 -1.4037237
# 7: COX7A2L 1.184289  5.902787  2.865264 7.305707e-03 2.224450e-02 -3.3467227
# 8: METTL2B 1.826493  6.555492  2.744935 9.845730e-03 2.822122e-02 -3.6236674
# 9:    DAD1 1.730326  5.984755  2.733194 1.013349e-02 2.885385e-02 -3.6503040
#10:    HRNR 1.112490  6.765955  2.641228 1.267350e-02 3.456165e-02 -3.8564578
#11:  ZDHHC3 1.028406  9.229493  2.551282 1.571762e-02 4.094785e-02 -4.0536996


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_HuR_vs_GFP_VHH.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_HuR_vs_GFP_VHH.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_HuR_vs_GFP_VHH.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_HuR_vs_GFP_VHH.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR vs VHHGFP") +
theme_bw() +
coord_cartesian(xlim = c(-5, 5), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_RBCC_HuR_vs_GFP_VHH.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR vs VHHGFP") +
theme_bw() +
coord_cartesian(xlim = c(-5, 5), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5))

ggsave('figures/volcano_RBCC_HuR_vs_GFP_VHH_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)


#######################
# HuR_VHH_vs_RBCC_GFP #
#######################
detable <- data.table(topTable(fit2, coef="HuR_VHH_vs_RBCC_GFP", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn       logFC  AveExpr         t   P.Value adj.P.Val         B
#1: ELAVL1 -0.07630259 9.491571 -1.141454 0.2621571 0.4800636 -6.194031

detable[rn == "ELAVL2"]
#       rn      logFC  AveExpr         t     P.Value  adj.P.Val         B
#1: ELAVL2 -0.2090457 10.96142 -3.061863 0.004434693 0.02792775 -2.659425

nrow(detable[adj.P.Val < 0.05]) # 1271
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 9
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 37

detable[logFC < -1 & adj.P.Val < 0.05]
#        rn     logFC   AveExpr          t      P.Value    adj.P.Val         B
#1: CCDC171 -6.064809  8.859662 -24.696080 2.135180e-22 1.378259e-18 39.441877
#2:    NME3 -1.046601  8.651534 -11.829631 3.224776e-13 2.312881e-10 20.106009
#3:   PTOV1 -1.207874  7.253424 -11.412624 8.207586e-13 4.414997e-10 19.202223
#4:   PCBP4 -1.046355 11.197335  -7.837842 6.126086e-09 7.908776e-07 10.481990
#5:     UBC -1.144932  9.372804  -7.229558 3.279973e-08 2.550871e-06  8.830629
#6:  FYTTD1 -1.250840  5.745152  -5.419005 5.881924e-06 1.586371e-04  3.729170
#7:  SNAPC1 -1.112054  4.203166  -4.911477 2.574008e-05 5.143955e-04  2.285735
#8:    IGKC -1.799652  3.744130  -3.227670 2.881088e-03 2.036913e-02 -2.254715
#9:  OXNAD1 -1.115737  6.943503  -2.855432 7.487590e-03 4.041170e-02 -3.146850

detable[logFC > 1 & adj.P.Val < 0.05]


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/HuR_VHH_vs_RBCC_GFP.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/HuR_VHH_vs_RBCC_GFP.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/HuR_VHH_vs_RBCC_GFP.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/HuR_VHH_vs_RBCC_GFP.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("VHHHuR vs TRIM21-VHHGFP") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_HuR_VHH_vs_RBCC_GFP.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("VHHHuR vs TRIM21-VHHGFP") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5))

ggsave('figures/volcano_HuR_VHH_vs_RBCC_GFP_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)


#######################
# RBCC_GFP_vs_HuR_VHH #
#######################
detable <- data.table(topTable(fit2, coef="RBCC_GFP_vs_HuR_VHH", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn       logFC  AveExpr         t   P.Value adj.P.Val         B
#1: ELAVL1 -0.07630259 9.491571 -1.141454 0.2621571 0.4800636 -6.194031

detable[rn == "ELAVL2"]
#       rn      logFC  AveExpr         t     P.Value  adj.P.Val         B
#1: ELAVL2 -0.2090457 10.96142 -3.061863 0.004434693 0.02792775 -2.659425

nrow(detable[adj.P.Val < 0.05]) # 1271
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 37
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 9

detable[logFC < -1 & adj.P.Val < 0.05]

detable[logFC > 1 & adj.P.Val < 0.05]
#        rn    logFC   AveExpr         t      P.Value    adj.P.Val         B
#1: CCDC171 6.064809  8.859662 24.696080 2.135180e-22 1.378259e-18 39.441877
#2:    NME3 1.046601  8.651534 11.829631 3.224776e-13 2.312881e-10 20.106009
#3:   PTOV1 1.207874  7.253424 11.412624 8.207586e-13 4.414997e-10 19.202223
#4:   PCBP4 1.046355 11.197335  7.837842 6.126086e-09 7.908776e-07 10.481990
#5:     UBC 1.144932  9.372804  7.229558 3.279973e-08 2.550871e-06  8.830629
#6:  FYTTD1 1.250840  5.745152  5.419005 5.881924e-06 1.586371e-04  3.729170
#7:  SNAPC1 1.112054  4.203166  4.911477 2.574008e-05 5.143955e-04  2.285735
#8:    IGKC 1.799652  3.744130  3.227670 2.881088e-03 2.036913e-02 -2.254715
#9:  OXNAD1 1.115737  6.943503  2.855432 7.487590e-03 4.041170e-02 -3.146850


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_GFP_vs_HuR_VHH.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_GFP_vs_HuR_VHH.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_GFP_vs_HuR_VHH.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_GFP_vs_HuR_VHH.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHGFP vs VHHHuR") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_RBCC_GFP_vs_HuR_VHH.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHGFP vs VHHHuR") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5))

ggsave('figures/volcano_RBCC_GFP_vs_HuR_VHH_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)


######################
# HuR_VHH_vs_GFP_VHH #
######################
detable <- data.table(topTable(fit2, coef="HuR_VHH_vs_GFP_VHH", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn       logFC  AveExpr          t   P.Value adj.P.Val         B
#1: ELAVL1 -0.02778247 9.491571 -0.4156137 0.6804685 0.9997142 -6.038788

detable[rn == "ELAVL2"]
#       rn     logFC  AveExpr        t     P.Value adj.P.Val         B
#1: ELAVL2 0.2036401 10.96142 2.982688 0.005431222 0.5393622 -2.290516

nrow(detable[adj.P.Val < 0.05]) # 26
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 4
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 2

detable[logFC < -1 & adj.P.Val < 0.05]
#      rn     logFC   AveExpr          t      P.Value    adj.P.Val         B
#1:  ISCU -1.198841  8.145434 -10.302813 1.095924e-11 2.358063e-08 15.609096
#2: ADAT2 -1.147498  6.931987  -9.979075 2.404461e-11 3.880198e-08 14.944379
#3: PCBP4 -1.032542 11.197335  -7.734379 8.124987e-09 6.555849e-06  9.874131
#4:  SPEG -2.872693  9.396652  -4.212071 1.925289e-04 4.779901e-02  0.729110

detable[logFC > 1 & adj.P.Val < 0.05]
#       rn    logFC  AveExpr        t      P.Value  adj.P.Val        B
#1:  GON4L 1.239836 6.168955 4.782518 3.740855e-05 0.01270906 2.227674
#2: TECPR2 1.054601 6.268936 4.577066 6.772036e-05 0.02185675 1.684422


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/HuR_VHH_vs_GFP_VHH.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/HuR_VHH_vs_GFP_VHH.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/HuR_VHH_vs_GFP_VHH.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/HuR_VHH_vs_GFP_VHH.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("VHHHuR vs VHHGFP") +
theme_bw() +
coord_cartesian(xlim = c(-5, 5), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_HuR_VHH_vs_GFP_VHH.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("VHHHuR vs VHHGFP") +
theme_bw() +
coord_cartesian(xlim = c(-5, 5), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5))

ggsave('figures/volcano_HuR_VHH_vs_GFP_VHH_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)


#######################
# RBCC_GFP_vs_GFP_VHH #
#######################
detable <- data.table(topTable(fit2, coef="RBCC_GFP_vs_GFP_VHH", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn      logFC  AveExpr         t   P.Value adj.P.Val         B
#1: ELAVL1 0.04852013 9.491571 0.7258401 0.4732161 0.6637403 -6.578296

detable[rn == "ELAVL2"]
#       rn     logFC  AveExpr        t      P.Value    adj.P.Val        B
#1: ELAVL2 0.4126857 10.96142 6.044552 9.567837e-07 3.529165e-05 5.511707

nrow(detable[adj.P.Val < 0.05]) # 1396
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 43
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 5

detable[logFC < -1 & adj.P.Val < 0.05]

detable[logFC > 1 & adj.P.Val < 0.05]
#        rn    logFC  AveExpr         t      P.Value    adj.P.Val         B
#1: CCDC171 6.339117 8.859662 25.813073 5.540906e-23 3.576655e-19 40.564604
#2:   PTOV1 1.086350 7.253424 10.264404 1.202144e-11 4.849898e-09 16.592761
#3:     UBC 1.125356 9.372804  7.105942 4.635692e-08 3.324821e-06  8.490208
#4:  SNAPC1 1.033420 4.203166  4.564187 7.027982e-05 1.103786e-03  1.307954
#5:  SYNGR3 2.406225 6.738369  3.237616 2.806721e-03 1.795578e-02 -2.229779


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_GFP_vs_GFP_VHH.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_GFP_vs_GFP_VHH.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_GFP_vs_GFP_VHH.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/RBCC_GFP_vs_GFP_VHH.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHGFP vs VHHGFP") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_RBCC_GFP_vs_GFP_VHH.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHGFP vs VHHGFP") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 30), clip = "off") +
annotate("text", x = -4, y = 30, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 30, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=14), axis.text = element_text(size=14, color = "black"), plot.title = element_text(size=16, hjust = 0.5))

ggsave('figures/volcano_RBCC_GFP_vs_GFP_VHH_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)
```



### venn diagrams and intersecting lists

triple venn diagram 1:
- RBCC_HuR vs HuR_VHH
- RBCC_HuR vs RBCC_GFP
- RBCC_HuR vs GFP_VHH

double venn diagram 1:
- RBCC_HuR vs HuR_VHH
- RBCC_HuR vs GFP_VHH

double venn diagram 2:
- RBCC_GFP vs HuR_VHH
- RBCC_GFP vs GFP_VHH

double venn diagram 3:
- intersection double venn diagram 1
- intersection double venn diagram 2

double venn diagram 4:
- RBCC_HuR vs HuR_VHH
- RBCC_GFP vs GFP_VHH

triple venn diagram 2:
- RBCC_HuR vs HuR_VHH
- HuR_VHH vs RBCC_GFP
- HuR_VHH vs GFP_VHH


Name keys:
- RBCC-HuR = TRIM21-VHHHuR
- HuR-RBCC = VHHHuR-TRIM21
- RBCC-GFP = TRIM21-VHHGFP
- GFP-RBCC = VHHGFP-TRIM21
- HuR VHH = VHHHuR
- GFP VHH = VHHGFP


```r
library(data.table)
library(VennDiagram)


# Enlarge the view width when printing tables
options(width = 300)


#triple venn diagram 1:
#- RBCC_HuR vs HuR_VHH
#- RBCC_HuR vs RBCC_GFP
#- RBCC_HuR vs GFP_VHH

RBCC_HuR_vs_HuR_VHH <- fread("RBCC_HuR_vs_HuR_VHH.down.txt")
RBCC_HuR_vs_RBCC_GFP <- fread("RBCC_HuR_vs_RBCC_GFP.down.txt")
RBCC_HuR_vs_GFP_VHH <- fread("RBCC_HuR_vs_GFP_VHH.down.txt")

venn.plot <- draw.triple.venn(
  area1 = length(RBCC_HuR_vs_HuR_VHH$rn),
  area2 = length(RBCC_HuR_vs_RBCC_GFP$rn),
  area3 = length(RBCC_HuR_vs_GFP_VHH$rn),
  n12 = length(intersect(RBCC_HuR_vs_HuR_VHH$rn,RBCC_HuR_vs_RBCC_GFP$rn)),
  n23 = length(intersect(RBCC_HuR_vs_RBCC_GFP$rn,RBCC_HuR_vs_GFP_VHH$rn)),
  n13 = length(intersect(RBCC_HuR_vs_HuR_VHH$rn,RBCC_HuR_vs_GFP_VHH$rn)),
  n123 = length(intersect(intersect(RBCC_HuR_vs_HuR_VHH$rn,RBCC_HuR_vs_RBCC_GFP$rn),RBCC_HuR_vs_GFP_VHH$rn)),
  category = c(sprintf("TRIM21-VHHHuR vs VHHHuR\n (%s)", length(RBCC_HuR_vs_HuR_VHH$rn)), sprintf("TRIM21-VHHHuR vs TRIM21-VHHGFP\n (%s)", length(RBCC_HuR_vs_RBCC_GFP$rn)), sprintf("TRIM21-VHHHuR vs VHHGFP\n (%s)", length(RBCC_HuR_vs_GFP_VHH$rn))),
  fill = c("#830051", "#C4D600", "#003865"),
  cex = 1.5,
  fontfamily = "sans",
  cat.dist = c(0.15, 0.15, 0.15),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  print.mode = "raw",
  margin = 0.25)

pdf("../figures/venn_RBCC_HuR_vs_HuR_VHH_RBCC_HuR_vs_RBCC_GFP_RBCC_HuR_vs_GFP_VHH_downregulated.pdf", width = 24/2.54, height = 24/2.54, useDingbats = FALSE)
g <- grid.draw(venn.plot)
dev.off()

# lists
write.table(setdiff(intersect(RBCC_HuR_vs_HuR_VHH$rn,RBCC_HuR_vs_RBCC_GFP$rn),RBCC_HuR_vs_GFP_VHH$rn), "venn_RBCC_HuR_vs_HuR_VHH_RBCC_HuR_vs_RBCC_GFP_notRBCC_HuR_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(intersect(RBCC_HuR_vs_RBCC_GFP$rn,RBCC_HuR_vs_GFP_VHH$rn),RBCC_HuR_vs_HuR_VHH$rn), "venn_notRBCC_HuR_vs_HuR_VHH_RBCC_HuR_vs_RBCC_GFP_RBCC_HuR_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(intersect(RBCC_HuR_vs_HuR_VHH$rn,RBCC_HuR_vs_GFP_VHH$rn),RBCC_HuR_vs_RBCC_GFP$rn), "venn_RBCC_HuR_vs_HuR_VHH_notRBCC_HuR_vs_RBCC_GFP_RBCC_HuR_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(intersect(intersect(RBCC_HuR_vs_HuR_VHH$rn,RBCC_HuR_vs_RBCC_GFP$rn),RBCC_HuR_vs_GFP_VHH$rn), "venn_RBCC_HuR_vs_HuR_VHH_RBCC_HuR_vs_RBCC_GFP_RBCC_HuR_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)


#double venn diagram 1:
#- RBCC_HuR vs HuR_VHH
#- RBCC_HuR vs GFP_VHH

RBCC_HuR_vs_HuR_VHH <- fread("RBCC_HuR_vs_HuR_VHH.down.txt")
RBCC_HuR_vs_GFP_VHH <- fread("RBCC_HuR_vs_GFP_VHH.down.txt")

venn.plot <- draw.pairwise.venn(
  area1 = length(RBCC_HuR_vs_HuR_VHH$rn),
  area2 = length(RBCC_HuR_vs_GFP_VHH$rn),
  cross.area = length(intersect(RBCC_HuR_vs_HuR_VHH$rn, RBCC_HuR_vs_GFP_VHH$rn)),
  category = c(sprintf("TRIM21-VHHHuR vs VHHHuR\n (%s)", length(RBCC_HuR_vs_HuR_VHH$rn)), sprintf("TRIM21-VHHHuR vs VHHGFP\n (%s)", length(RBCC_HuR_vs_GFP_VHH$rn))),
  euler.d = FALSE,
  scaled = FALSE,
  inverted = FALSE,
  fill = c("#830051", "#003865"),
  cex = 1.5,
  fontfamily = "sans",
  cat.pos = c(45, -45),
  cat.dist = c(0.15, 0.15),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  print.mode = "raw",
  margin = 0.25)

pdf("../figures/venn_RBCC_HuR_vs_HuR_VHH_RBCC_HuR_vs_GFP_VHH_downregulated.pdf", width = 24/2.54, height = 24/2.54, useDingbats = FALSE)
g <- grid.draw(venn.plot)
dev.off()

# lists
write.table(intersect(RBCC_HuR_vs_HuR_VHH$rn,RBCC_HuR_vs_GFP_VHH$rn), "venn_RBCC_HuR_vs_HuR_VHH_RBCC_HuR_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(RBCC_HuR_vs_HuR_VHH$rn,RBCC_HuR_vs_GFP_VHH$rn), "venn_RBCC_HuR_vs_HuR_VHH_notRBCC_HuR_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(RBCC_HuR_vs_GFP_VHH$rn,RBCC_HuR_vs_HuR_VHH$rn), "venn_notRBCC_HuR_vs_HuR_VHH_RBCC_HuR_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)


#double venn diagram 2:
#- RBCC_GFP vs HuR_VHH
#- RBCC_GFP vs GFP_VHH

RBCC_GFP_vs_HuR_VHH <- fread("RBCC_GFP_vs_HuR_VHH.down.txt")
RBCC_GFP_vs_GFP_VHH <- fread("RBCC_GFP_vs_GFP_VHH.down.txt")

venn.plot <- draw.pairwise.venn(
  area1 = length(RBCC_GFP_vs_HuR_VHH$rn),
  area2 = length(RBCC_GFP_vs_GFP_VHH$rn),
  cross.area = length(intersect(RBCC_GFP_vs_HuR_VHH$rn, RBCC_GFP_vs_GFP_VHH$rn)),
  category = c(sprintf("TRIM21-VHHGFP vs VHHHuR\n (%s)", length(RBCC_GFP_vs_HuR_VHH$rn)), sprintf("TRIM21-VHHGFP vs VHHGFP\n (%s)", length(RBCC_GFP_vs_GFP_VHH$rn))),
  euler.d = FALSE,
  scaled = FALSE,
  inverted = FALSE,
  fill = c("#830051", "#003865"),
  cex = 1.5,
  fontfamily = "sans",
  cat.pos = c(45, -45),
  cat.dist = c(0.15, 0.15),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  print.mode = "raw",
  margin = 0.25)

pdf("../figures/venn_RBCC_GFP_vs_HuR_VHH_RBCC_GFP_vs_GFP_VHH_downregulated.pdf", width = 24/2.54, height = 24/2.54, useDingbats = FALSE)
g <- grid.draw(venn.plot)
dev.off()

# lists
write.table(intersect(RBCC_GFP_vs_HuR_VHH$rn,RBCC_GFP_vs_GFP_VHH$rn), "venn_RBCC_GFP_vs_HuR_VHH_RBCC_GFP_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(RBCC_GFP_vs_HuR_VHH$rn,RBCC_GFP_vs_GFP_VHH$rn), "venn_RBCC_GFP_vs_HuR_VHH_notRBCC_GFP_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(RBCC_GFP_vs_GFP_VHH$rn,RBCC_GFP_vs_HuR_VHH$rn), "venn_notRBCC_GFP_vs_HuR_GFP_RBCC_HuR_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)


#double venn diagram 3:
#- intersection double venn diagram 1
#- intersection double venn diagram 2

venn1 <- fread("venn_RBCC_HuR_vs_HuR_VHH_RBCC_HuR_vs_GFP_VHH_down.txt", header = FALSE)
venn2 <- fread("venn_RBCC_GFP_vs_HuR_VHH_RBCC_GFP_vs_GFP_VHH_down.txt", header = FALSE)

venn.plot <- draw.pairwise.venn(
  area1 = length(venn1$V1),
  area2 = length(venn2$V1),
  cross.area = length(intersect(venn1$V1, venn2$V1)),
  category = c(sprintf("Venn 1\n (%s)", length(venn1$V1)), sprintf("Venn 2\n (%s)", length(venn2$V1))),
  euler.d = FALSE,
  scaled = FALSE,
  inverted = FALSE,
  fill = c("#830051", "#003865"),
  cex = 1.5,
  fontfamily = "sans",
  cat.pos = c(-45, 45),
  cat.dist = c(0.15, 0.15),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  print.mode = "raw",
  margin = 0.25)

pdf("../figures/venn_venn1_venn2_downregulated.pdf", width = 24/2.54, height = 24/2.54, useDingbats = FALSE)
g <- grid.draw(venn.plot)
dev.off()

# lists
write.table(intersect(venn1$V1,venn2$V1), "venn_venn1_venn2_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(venn1$V1,venn2$V1), "venn_venn1_notvenn2_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(venn2$V1,venn1$V1), "venn_notvenn1_venn2_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)


#double venn diagram 4:
#- RBCC_HuR vs HuR_VHH
#- RBCC_GFP vs GFP_VHH

RBCC_HuR_vs_HuR_VHH <- fread("RBCC_HuR_vs_HuR_VHH.down.txt")
RBCC_GFP_vs_GFP_VHH <- fread("RBCC_GFP_vs_GFP_VHH.down.txt")

venn.plot <- draw.pairwise.venn(
  area1 = length(RBCC_HuR_vs_HuR_VHH$rn),
  area2 = length(RBCC_GFP_vs_GFP_VHH$rn),
  cross.area = length(intersect(RBCC_HuR_vs_HuR_VHH$rn, RBCC_GFP_vs_GFP_VHH$rn)),
  category = c(sprintf("TRIM21-VHHHuR vs VHHHuR\n (%s)", length(RBCC_HuR_vs_HuR_VHH$rn)), sprintf("TRIM21-VHHGFP vs VHHGFP\n (%s)", length(RBCC_GFP_vs_GFP_VHH$rn))),
  euler.d = FALSE,
  scaled = FALSE,
  inverted = FALSE,
  fill = c("#830051", "#003865"),
  cex = 1.5,
  fontfamily = "sans",
  cat.pos = c(-45, 45),
  cat.dist = c(0.15, 0.15),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  print.mode = "raw",
  margin = 0.25)

pdf("../figures/venn_RBCC_HuR_vs_HuR_VHH_RBCC_GFP_vs_GFP_VHH_downregulated.pdf", width = 24/2.54, height = 24/2.54, useDingbats = FALSE)
g <- grid.draw(venn.plot)
dev.off()

# lists
write.table(intersect(RBCC_HuR_vs_HuR_VHH$rn,RBCC_GFP_vs_GFP_VHH$rn), "venn_RBCC_HuR_vs_HuR_VHH_RBCC_GFP_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(RBCC_HuR_vs_HuR_VHH$rn,RBCC_GFP_vs_GFP_VHH$rn), "venn_RBCC_HuR_vs_HuR_VHH_notRBCC_GFP_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(RBCC_GFP_vs_GFP_VHH$rn,RBCC_HuR_vs_HuR_VHH$rn), "venn_notRBCC_GFP_vs_HuR_VHH_RBCC_HuR_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)


#triple venn diagram 2:
#- RBCC_HuR vs HuR_VHH
#- HuR_VHH vs RBCC_GFP
#- HuR_VHH vs GFP_VHH

RBCC_HuR_vs_HuR_VHH <- fread("RBCC_HuR_vs_HuR_VHH.down.txt")
HuR_VHH_vs_RBCC_GFP <- fread("HuR_VHH_vs_RBCC_GFP.down.txt")
HuR_VHH_vs_GFP_VHH <- fread("HuR_VHH_vs_GFP_VHH.down.txt")

venn.plot <- draw.triple.venn(
  area1 = length(RBCC_HuR_vs_HuR_VHH$rn),
  area2 = length(HuR_VHH_vs_RBCC_GFP$rn),
  area3 = length(HuR_VHH_vs_GFP_VHH$rn),
  n12 = length(intersect(RBCC_HuR_vs_HuR_VHH$rn,HuR_VHH_vs_RBCC_GFP$rn)),
  n23 = length(intersect(HuR_VHH_vs_RBCC_GFP$rn,HuR_VHH_vs_GFP_VHH$rn)),
  n13 = length(intersect(RBCC_HuR_vs_HuR_VHH$rn,HuR_VHH_vs_GFP_VHH$rn)),
  n123 = length(intersect(intersect(RBCC_HuR_vs_HuR_VHH$rn,HuR_VHH_vs_RBCC_GFP$rn),HuR_VHH_vs_GFP_VHH$rn)),
  category = c(sprintf("TRIM21-VHHHuR vs VHHHuR\n (%s)", length(RBCC_HuR_vs_HuR_VHH$rn)), sprintf("VHHHuR vs TRIM21-VHHGFP\n (%s)", length(HuR_VHH_vs_RBCC_GFP$rn)), sprintf("VHHHuR vs VHHGFP\n (%s)", length(HuR_VHH_vs_GFP_VHH$rn))),
  fill = c("#830051", "#C4D600", "#003865"),
  cex = 1.5,
  fontfamily = "sans",
  cat.dist = c(0.15, 0.15, 0.15),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  print.mode = "raw",
  margin = 0.25)

pdf("../figures/venn_RBCC_HuR_vs_HuR_VHH_HuR_VHH_vs_RBCC_GFP_HuR_VHH_vs_GFP_VHH_downregulated.pdf", width = 24/2.54, height = 24/2.54, useDingbats = FALSE)
g <- grid.draw(venn.plot)
dev.off()

# lists
write.table(setdiff(intersect(RBCC_HuR_vs_HuR_VHH$rn,HuR_VHH_vs_RBCC_GFP$rn),HuR_VHH_vs_GFP_VHH$rn), "venn_RRBCC_HuR_vs_HuR_VHH_HuR_VHH_vs_RBCC_GFP_notHuR_VHH_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(intersect(HuR_VHH_vs_RBCC_GFP$rn,HuR_VHH_vs_GFP_VHH$rn),RBCC_HuR_vs_HuR_VHH$rn), "venn_notRRBCC_HuR_vs_HuR_VHH_HuR_VHH_vs_RBCC_GFP_HuR_VHH_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(intersect(RBCC_HuR_vs_HuR_VHH$rn,HuR_VHH_vs_GFP_VHH$rn),HuR_VHH_vs_RBCC_GFP$rn), "venn_RRBCC_HuR_vs_HuR_VHH_notHuR_VHH_vs_RBCC_GFP_HuR_VHH_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(intersect(intersect(RBCC_HuR_vs_HuR_VHH$rn,HuR_VHH_vs_RBCC_GFP$rn),HuR_VHH_vs_GFP_VHH$rn), "venn_RRBCC_HuR_vs_HuR_VHH_HuR_VHH_vs_RBCC_GFP_HuR_VHH_vs_GFP_VHH_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
```





## kinetics

### data

- Time points: 0, 24, 48, 72h
- Treatments: +Dox, -Dox
- Plate: A, B

- Comparisons to be made within plate A:
- 0h +DOX vs 0h -DOX
- 24h +DOX vs 0h +DOX
- 48h +DOX vs 0h +DOX
- 72h +DOX vs 0h +DOX
- 48h +DOX vs 24h +DOX
- 72h +DOX vs 24h +DOX
- 72h +DOX vs 48h +DOX



### limma

```{r}
library(devtools)
load_all("/projects/qbio/bifo/R/3.6.0/pxanalytics")
library(tidyverse)
library(yaml)
library(data.table)
library(limma)
library(ggplot2)
library(ggrepel)


# Enlarge the view width when printing tables
options(width = 300)


# Load data
data <- fread("20210415_094022_DOA078_041021_Protac_IVEB_Trim21_kinetic_HybridLib0411_041421_ReportProteins.txt")
data <- unique(data, by = "genes")
data <- data.frame(data)
rownames(data) <- data$genes


# Assay data
assay <- na.omit(data[, grep("_A_", colnames(data))]) # selecting plate A only
nrow(data) # 6494


# Annotations
anno <- data.frame("ID" = rownames(assay), stringsAsFactors = FALSE)


# Metadata
metadata <- data.frame("ID" = colnames(assay), "Sample" = sapply(colnames(assay), function (x) paste(unlist(strsplit(x, "_"))[1:3], collapse="_")), "Time" = sapply(colnames(assay), function (x) unlist(strsplit(x, "_"))[1]), "Dox" = sapply(colnames(assay), function (x) unlist(strsplit(x, "_"))[2]), "Plate" = sapply(colnames(assay), function (x) unlist(strsplit(x, "_"))[3]), "Replicate" = sapply(colnames(assay), function (x) unlist(strsplit(x, "_"))[4]), stringsAsFactors = FALSE)


# pxdata object initialisation
px_raw <- pxinit(assay, annotations = anno, metadata = metadata)


# Normalisation
## Which DIA software / normalisation they used if any?
px_norm <- px_raw %>% pxnormalise(method = "total")


# model.matrix, lmFit, makeContrasts, contrasts.fit, eBayes
des <- model.matrix(~ 0 + metadata$Sample)
colnames(des) <- levels(factor(metadata$Sample))

fit <- lmFit(pxdata_logged(px_norm), des)

contrast.matrix <- makeContrasts(
  "t0h_pDox_A_vs_t0h_mDox_A" = t0h_pDox_A-t0h_mDox_A,
  "t24h_pDox_A_vs_t0h_pDox_A" = t24h_pDox_A-t0h_pDox_A,
  "t48h_pDox_A_vs_t0h_pDox_A" = t48h_pDox_A-t0h_pDox_A,
  "t72h_pDox_A_vs_t0h_pDox_A" = t72h_pDox_A-t0h_pDox_A,
  "t48h_pDox_A_vs_t24h_pDox_A" = t48h_pDox_A-t24h_pDox_A,
  "t72h_pDox_A_vs_t24h_pDox_A" = t72h_pDox_A-t24h_pDox_A,
  "t72h_pDox_A_vs_t48h_pDox_A" = t72h_pDox_A-t48h_pDox_A,
  levels=des)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)


# biogrid data
biogrid <- fread("../materials/BioGRID_List_of_Proteins.csv")
biogrid_vector <- setdiff(unique(c(biogrid[,8][[1]], biogrid[,9][[1]])), c("ELAVL1"))
biogrid_vector <- append(biogrid_vector, c("SGO1", "BRF1", "NR2F2", "GOT1", "FH", "SAT2", ";TUBB3"))


############################
# t0h_pDox_A_vs_t0h_mDox_A #
############################
detable <- data.table(topTable(fit2, coef="t0h_pDox_A_vs_t0h_mDox_A", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn      logFC  AveExpr         t   P.Value adj.P.Val         B
#1: ELAVL1 -0.2659035 9.208881 -2.180214 0.0627846 0.4863435 -4.254519

detable[rn == "ELAVL2"]
#       rn      logFC  AveExpr         t     P.Value adj.P.Val         B
#1: ELAVL2 -0.2693578 9.913181 -3.681973 0.006846156 0.3393812 -2.122298

nrow(detable[adj.P.Val < 0.05]) # 11
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 2
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 5

detable[logFC < -1 & adj.P.Val < 0.05]
#        rn     logFC  AveExpr         t      P.Value  adj.P.Val        B
#1: CFAP298 -3.383978 9.878656 -9.594861 1.693959e-05 0.01805490 3.201224
#2:    FLG2 -1.074976 8.376421 -7.686019 7.872091e-05 0.04647396 1.984695

detable[logFC > 1 & adj.P.Val < 0.05]
#         rn    logFC   AveExpr         t      P.Value   adj.P.Val        B
#1:   RAB3IP 3.813450  6.973100 16.992738 2.694250e-07 0.001749646 5.689137
#2: C18orf21 2.102045  6.906558 14.998903 6.754821e-07 0.002193290 5.244289
#3:  NDUFAB1 1.055649  8.704001  9.551119 1.749362e-05 0.018054903 3.177195
#4:    RPLP1 1.249345 11.112244  8.674736 3.428823e-05 0.027833468 2.659962
#5:    PYGO2 1.122941  7.127258  8.475401 4.027969e-05 0.029064035 2.532146


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t0h_pDox_A_vs_t0h_mDox_A.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t0h_pDox_A_vs_t0h_mDox_A.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t0h_pDox_A_vs_t0h_mDox_A.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t0h_pDox_A_vs_t0h_mDox_A.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 0h +DOX vs TRIM21-VHHHuR 0h -DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_t0h_pDox_A_vs_t0h_mDox_A.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 0h +DOX vs TRIM21-VHHHuR 0h -DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5))

ggsave('figures/volcano_t0h_pDox_A_vs_t0h_mDox_A_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)


#############################
# t24h_pDox_A_vs_t0h_pDox_A #
#############################
detable <- data.table(topTable(fit2, coef="t24h_pDox_A_vs_t0h_pDox_A", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn      logFC  AveExpr        t      P.Value   adj.P.Val        B
#1: ELAVL1 -0.9856565 9.208881 -8.08166 5.589076e-05 0.009485456 2.450862

detable[rn == "ELAVL2"]
#       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
#1: ELAVL2 -1.177401 9.913181 -16.09443 4.022168e-07 0.0003264995 7.068696

nrow(detable[adj.P.Val < 0.05]) # 149
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 11
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 22

detable[logFC < -1 & adj.P.Val < 0.05]
#          rn     logFC  AveExpr          t      P.Value    adj.P.Val         B
# 1: C18orf21 -2.277612 6.906558 -16.251640 3.744110e-07 0.0003264995 7.1276115
# 2:   ELAVL2 -1.177401 9.913181 -16.094431 4.022168e-07 0.0003264995 7.0686963
# 3: TNFRSF6B -1.720296 7.556737 -14.327475 9.449291e-07 0.0005578518 6.3441307
# 4:    TRPM4 -1.322897 5.705686 -12.138233 3.158016e-06 0.0014648682 5.2578888
# 5:    YLPM1 -1.035213 8.197082 -11.087571 6.060059e-06 0.0026236014 4.6449426
# 6:  ZDHHC20 -1.999544 6.300255 -10.259598 1.054496e-05 0.0032609030 4.1116814
# 7:  DPY19L3 -1.681464 5.188724  -8.036597 5.807391e-05 0.0094854559 2.4119007
# 8:    RUNX1 -1.486513 8.223450  -7.799943 7.122575e-05 0.0105122736 2.2038299
# 9: C11orf68 -1.085339 8.973387  -6.794023 1.799095e-04 0.0192323832 1.2500216
#10:   TCF7L2 -1.042845 6.052560  -6.024942 3.928042e-04 0.0277268510 0.4369422
#11:   PHLDA1 -1.398805 7.844815  -6.004334 4.014932e-04 0.0280354516 0.4140681

detable[logFC > 1 & adj.P.Val < 0.05]


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t24h_pDox_A_vs_t0h_pDox_A.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t24h_pDox_A_vs_t0h_pDox_A.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t24h_pDox_A_vs_t0h_pDox_A.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t24h_pDox_A_vs_t0h_pDox_A.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 24h +DOX vs TRIM21-VHHHuR 0h +DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_t24h_pDox_A_vs_t0h_pDox_A.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 24h +DOX vs TRIM21-VHHHuR 0h +DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5))

ggsave('figures/volcano_t24h_pDox_A_vs_t0h_pDox_A_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)


#############################
# t48h_pDox_A_vs_t0h_pDox_A #
#############################
detable <- data.table(topTable(fit2, coef="t48h_pDox_A_vs_t0h_pDox_A", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
#1: ELAVL1 -1.853956 9.208881 -15.20108 6.122353e-07 0.0001757083 6.925642

detable[rn == "ELAVL2"]
#       rn     logFC  AveExpr        t      P.Value   adj.P.Val        B
#1: ELAVL2 -1.861626 9.913181 -25.4474 1.335105e-08 1.24943e-05 10.38878

nrow(detable[adj.P.Val < 0.05]) # 1042
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 95
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 76

detable[logFC < -1 & adj.P.Val < 0.05]

detable[logFC > 1 & adj.P.Val < 0.05]


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t48h_pDox_A_vs_t0h_pDox_A.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t48h_pDox_A_vs_t0h_pDox_A.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t48h_pDox_A_vs_t0h_pDox_A.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t48h_pDox_A_vs_t0h_pDox_A.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 48h +DOX vs TRIM21-VHHHuR 0h +DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_t48h_pDox_A_vs_t0h_pDox_A.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 48h +DOX vs TRIM21-VHHHuR 0h +DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5))

ggsave('figures/volcano_t48h_pDox_A_vs_t0h_pDox_A_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)


#############################
# t72h_pDox_A_vs_t0h_pDox_A #
#############################
detable <- data.table(topTable(fit2, coef="t72h_pDox_A_vs_t0h_pDox_A", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn     logFC  AveExpr         t      P.Value    adj.P.Val        B
#1: ELAVL1 -1.903302 9.208881 -15.60568 5.047299e-07 0.0001869938 7.107672

detable[rn == "ELAVL2"]
#       rn    logFC  AveExpr         t      P.Value    adj.P.Val        B
#1: ELAVL2 -1.74659 9.913181 -23.87492 2.150949e-08 2.328044e-05 9.934146

nrow(detable[adj.P.Val < 0.05]) # 1221
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 95
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 77

detable[logFC < -1 & adj.P.Val < 0.05]

detable[logFC > 1 & adj.P.Val < 0.05]


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t72h_pDox_A_vs_t0h_pDox_A.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t72h_pDox_A_vs_t0h_pDox_A.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t72h_pDox_A_vs_t0h_pDox_A.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t72h_pDox_A_vs_t0h_pDox_A.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 72h +DOX vs TRIM21-VHHHuR 0h +DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_t72h_pDox_A_vs_t0h_pDox_A.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 72h +DOX vs TRIM21-VHHHuR 0h +DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5))

ggsave('figures/volcano_t72h_pDox_A_vs_t0h_pDox_A_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)


##############################
# t48h_pDox_A_vs_t24h_pDox_A #
##############################
detable <- data.table(topTable(fit2, coef="t48h_pDox_A_vs_t24h_pDox_A", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn      logFC  AveExpr         t      P.Value   adj.P.Val        B
#1: ELAVL1 -0.8682991 9.208881 -7.119415 0.0001318563 0.008477971 1.571003

detable[rn == "ELAVL2"]
#       rn      logFC  AveExpr         t      P.Value   adj.P.Val        B
#1: ELAVL2 -0.6842245 9.913181 -9.352973 2.027185e-05 0.003330626 3.469669

nrow(detable[adj.P.Val < 0.05]) # 367
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 32
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 21

detable[logFC < -1 & adj.P.Val < 0.05]

detable[logFC > 1 & adj.P.Val < 0.05]


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t48h_pDox_A_vs_t24h_pDox_A.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t48h_pDox_A_vs_t24h_pDox_A.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t48h_pDox_A_vs_t24h_pDox_A.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t48h_pDox_A_vs_t24h_pDox_A.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 48h +DOX vs TRIM21-VHHHuR 24h +DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_t48h_pDox_A_vs_t24h_pDox_A.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 48h +DOX vs TRIM21-VHHHuR 24h +DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5))

ggsave('figures/volcano_t48h_pDox_A_vs_t24h_pDox_A_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)


##############################
# t72h_pDox_A_vs_t24h_pDox_A #
##############################
detable <- data.table(topTable(fit2, coef="t72h_pDox_A_vs_t24h_pDox_A", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn      logFC  AveExpr         t      P.Value   adj.P.Val        B
#1: ELAVL1 -0.9176453 9.208881 -7.524018 9.094409e-05 0.006359914 1.961061

detable[rn == "ELAVL2"]
#       rn      logFC  AveExpr         t      P.Value   adj.P.Val        B
#1: ELAVL2 -0.5691881 9.913181 -7.780489 7.244709e-05 0.006097685 2.193054

nrow(detable[adj.P.Val < 0.05]) # 655
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 37
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 24

detable[logFC < -1 & adj.P.Val < 0.05]

detable[logFC > 1 & adj.P.Val < 0.05]


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t72h_pDox_A_vs_t24h_pDox_A.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t72h_pDox_A_vs_t24h_pDox_A.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t72h_pDox_A_vs_t24h_pDox_A.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t72h_pDox_A_vs_t24h_pDox_A.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 72h +DOX vs TRIM21-VHHHuR 24h +DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_t72h_pDox_A_vs_t24h_pDox_A.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 72h +DOX vs TRIM21-VHHHuR 24h +DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5))

ggsave('figures/volcano_t72h_pDox_A_vs_t24h_pDox_A_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)


##############################
# t72h_pDox_A_vs_t48h_pDox_A #
##############################
detable <- data.table(topTable(fit2, coef="t72h_pDox_A_vs_t48h_pDox_A", adjust="BH", number = Inf, sort.by = 'P'), keep.rownames = TRUE)

detable[rn == "ELAVL1"]
#       rn       logFC  AveExpr          t   P.Value adj.P.Val         B
#1: ELAVL1 -0.04934619 9.208881 -0.4046026 0.6969743 0.9398155 -6.594986

detable[rn == "ELAVL2"]
#       rn     logFC  AveExpr        t   P.Value adj.P.Val         B
#1: ELAVL2 0.1150364 9.913181 1.572485 0.1566692 0.5490593 -5.495464

nrow(detable[adj.P.Val < 0.05]) # 118
nrow(detable[logFC < -1 & adj.P.Val < 0.05]) # 6
nrow(detable[logFC > 1 & adj.P.Val < 0.05]) # 20

detable[logFC < -1 & adj.P.Val < 0.05]

detable[logFC > 1 & adj.P.Val < 0.05]


# write tables
## down
write.table(detable[logFC < -1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t72h_pDox_A_vs_t48h_pDox_A.down.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t72h_pDox_A_vs_t48h_pDox_A.down.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
## up
write.table(detable[logFC > 1 & adj.P.Val < 0.05][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t72h_pDox_A_vs_t48h_pDox_A.up.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
write.table(detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector][,.(rn, logFC, P.Value, adj.P.Val)], "tables/t72h_pDox_A_vs_t48h_pDox_A.up.biogrid.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


# volcano
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "maroon") +
geom_vline(xintercept=-1, linetype="longdash", color = "maroon") +
geom_vline(xintercept=1, linetype="longdash", color = "maroon") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 72h +DOX vs TRIM21-VHHHuR 48h +DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5)) +
geom_text_repel(data = head(detable[logFC < -1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3) +
geom_text_repel(data = head(detable[logFC > 1 & adj.P.Val < 0.05], 10), aes(label = rn), force = 10, size = 3)

ggsave('figures/volcano_t72h_pDox_A_vs_t48h_pDox_A.pdf', width = 6, height = 6, useDingbats = FALSE)


# volcano biogrid
## FDR
gg <- ggplot(detable, aes(x = logFC, y = -log10(adj.P.Val))) +
geom_point(size = 1, col = "lightgray") +
geom_hline(yintercept=-log10(0.05), linetype="longdash", color = "black") +
geom_vline(xintercept=-1, linetype="longdash", color = "black") +
geom_vline(xintercept=1, linetype="longdash", color = "black") +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05], aes(x = logFC, y = -log10(adj.P.Val)), color='black') +
geom_point(data = detable[logFC < -1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[logFC > 1 & adj.P.Val < 0.05 & rn %in% biogrid_vector], aes(x = logFC, y = -log10(adj.P.Val)), color='forestgreen') +
geom_point(data = detable[rn=="ELAVL1"], aes(x = logFC, y = -log10(adj.P.Val)), color='maroon', size = 3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
ggtitle("TRIM21-VHHHuR 72h +DOX vs TRIM21-VHHHuR 48h +DOX") +
theme_bw() +
coord_cartesian(xlim = c(-6, 6), ylim = c(0, 10), clip = "off") +
annotate("text", x = -4, y = 10, label = sprintf("n = %s", nrow(detable[logFC < -1 & adj.P.Val < 0.05])), size = 5) +
annotate("text", x = 4, y = 10, label = sprintf("n = %s", nrow(detable[logFC > 1 & adj.P.Val < 0.05])), size = 5) +
theme(axis.title = element_text(size=12), axis.text = element_text(size=12, color = "black"), plot.title = element_text(size=12, hjust = 0.5))

ggsave('figures/volcano_t72h_pDox_A_vs_t48h_pDox_A_biogrid.pdf', width = 6, height = 6, useDingbats = FALSE)
```


### venn diagrams and intersecting lists

triple venn diagram down:
- 24h +DOX vs 0h +DOX
- 48h +DOX vs 0h +DOX
- 72h +DOX vs 0h +DOX

triple venn diagram up:
- 24h +DOX vs 0h +DOX
- 48h +DOX vs 0h +DOX
- 72h +DOX vs 0h +DOX

```r
library(data.table)
library(VennDiagram)


# Enlarge the view width when printing tables
options(width = 300)


#triple venn diagram down:
#- 24h +DOX vs 0h +DOX
#- 48h +DOX vs 0h +DOX
#- 72h +DOX vs 0h +DOX

t24h_pDox_A_vs_t0h_pDox_A <- fread("t24h_pDox_A_vs_t0h_pDox_A.down.txt")
t48h_pDox_A_vs_t0h_pDox_A <- fread("t48h_pDox_A_vs_t0h_pDox_A.down.txt")
t72h_pDox_A_vs_t0h_pDox_A <- fread("t72h_pDox_A_vs_t0h_pDox_A.down.txt")

venn.plot <- draw.triple.venn(
  area1 = length(t24h_pDox_A_vs_t0h_pDox_A$rn),
  area2 = length(t48h_pDox_A_vs_t0h_pDox_A$rn),
  area3 = length(t72h_pDox_A_vs_t0h_pDox_A$rn),
  n12 = length(intersect(t24h_pDox_A_vs_t0h_pDox_A$rn,t48h_pDox_A_vs_t0h_pDox_A$rn)),
  n23 = length(intersect(t48h_pDox_A_vs_t0h_pDox_A$rn,t72h_pDox_A_vs_t0h_pDox_A$rn)),
  n13 = length(intersect(t24h_pDox_A_vs_t0h_pDox_A$rn,t72h_pDox_A_vs_t0h_pDox_A$rn)),
  n123 = length(intersect(intersect(t24h_pDox_A_vs_t0h_pDox_A$rn,t48h_pDox_A_vs_t0h_pDox_A$rn),t72h_pDox_A_vs_t0h_pDox_A$rn)),
  category = c(sprintf("24h +DOX vs 0h +DOX\n (%s)", length(t24h_pDox_A_vs_t0h_pDox_A$rn)), sprintf("48h +DOX vs 0h +DOX\n (%s)", length(t48h_pDox_A_vs_t0h_pDox_A$rn)), sprintf("72h +DOX vs 0h +DOX\n (%s)", length(t72h_pDox_A_vs_t0h_pDox_A$rn))),
  fill = c("#830051", "#C4D600", "#003865"),
  cex = 1.5,
  fontfamily = "sans",
  cat.dist = c(0.15, 0.15, 0.15),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  print.mode = "raw",
  margin = 0.25)

pdf("../figures/venn_t24h_pDox_A_vs_t0h_pDox_A_t48h_pDox_A_vs_t0h_pDox_A_t72h_pDox_A_vs_t0h_pDox_A_downregulated.pdf", width = 24/2.54, height = 24/2.54, useDingbats = FALSE)
g <- grid.draw(venn.plot)
dev.off()

# lists
write.table(setdiff(intersect(t24h_pDox_A_vs_t0h_pDox_A$rn,t48h_pDox_A_vs_t0h_pDox_A$rn),t72h_pDox_A_vs_t0h_pDox_A$rn), "venn_t24h_pDox_A_vs_t0h_pDox_A_t48h_pDox_A_vs_t0h_pDox_A_not72h_pDox_A_vs_t0h_pDox_A_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

write.table(setdiff(intersect(t48h_pDox_A_vs_t0h_pDox_A$rn,t72h_pDox_A_vs_t0h_pDox_A$rn),t24h_pDox_A_vs_t0h_pDox_A$rn), "venn_not24h_pDox_A_vs_t0h_pDox_A_t48h_pDox_A_vs_t0h_pDox_A_t72h_pDox_A_vs_t0h_pDox_A_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

write.table(setdiff(intersect(t24h_pDox_A_vs_t0h_pDox_A$rn,t72h_pDox_A_vs_t0h_pDox_A$rn),t48h_pDox_A_vs_t0h_pDox_A$rn), "venn_t24h_pDox_A_vs_t0h_pDox_A_not48h_pDox_A_vs_t0h_pDox_A_t72h_pDox_A_vs_t0h_pDox_A_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

write.table(intersect(intersect(t24h_pDox_A_vs_t0h_pDox_A$rn,t48h_pDox_A_vs_t0h_pDox_A$rn),t72h_pDox_A_vs_t0h_pDox_A$rn), "venn_t24h_pDox_A_vs_t0h_pDox_A_t48h_pDox_A_vs_t0h_pDox_A_t72h_pDox_A_vs_t0h_pDox_A_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)



#triple venn diagram up:
#- 24h +DOX vs 0h +DOX
#- 48h +DOX vs 0h +DOX
#- 72h +DOX vs 0h +DOX

t24h_pDox_A_vs_t0h_pDox_A <- fread("t24h_pDox_A_vs_t0h_pDox_A.up.txt")
t48h_pDox_A_vs_t0h_pDox_A <- fread("t48h_pDox_A_vs_t0h_pDox_A.up.txt")
t72h_pDox_A_vs_t0h_pDox_A <- fread("t72h_pDox_A_vs_t0h_pDox_A.up.txt")

venn.plot <- draw.triple.venn(
  area1 = length(t24h_pDox_A_vs_t0h_pDox_A$rn),
  area2 = length(t48h_pDox_A_vs_t0h_pDox_A$rn),
  area3 = length(t72h_pDox_A_vs_t0h_pDox_A$rn),
  n12 = length(intersect(t24h_pDox_A_vs_t0h_pDox_A$rn,t48h_pDox_A_vs_t0h_pDox_A$rn)),
  n23 = length(intersect(t48h_pDox_A_vs_t0h_pDox_A$rn,t72h_pDox_A_vs_t0h_pDox_A$rn)),
  n13 = length(intersect(t24h_pDox_A_vs_t0h_pDox_A$rn,t72h_pDox_A_vs_t0h_pDox_A$rn)),
  n123 = length(intersect(intersect(t24h_pDox_A_vs_t0h_pDox_A$rn,t48h_pDox_A_vs_t0h_pDox_A$rn),t72h_pDox_A_vs_t0h_pDox_A$rn)),
  category = c(sprintf("24h +DOX vs 0h +DOX\n (%s)", length(t24h_pDox_A_vs_t0h_pDox_A$rn)), sprintf("48h +DOX vs 0h +DOX\n (%s)", length(t48h_pDox_A_vs_t0h_pDox_A$rn)), sprintf("72h +DOX vs 0h +DOX\n (%s)", length(t72h_pDox_A_vs_t0h_pDox_A$rn))),
  fill = c("#830051", "#C4D600", "#003865"),
  cex = 1.5,
  fontfamily = "sans",
  cat.dist = c(0.15, 0.15, 0.15),
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  print.mode = "raw",
  margin = 0.25)

pdf("../figures/venn_t24h_pDox_A_vs_t0h_pDox_A_t48h_pDox_A_vs_t0h_pDox_A_t72h_pDox_A_vs_t0h_pDox_A_upregulated.pdf", width = 24/2.54, height = 24/2.54, useDingbats = FALSE)
g <- grid.draw(venn.plot)
dev.off()

# lists
write.table(setdiff(intersect(t24h_pDox_A_vs_t0h_pDox_A$rn,t48h_pDox_A_vs_t0h_pDox_A$rn),t72h_pDox_A_vs_t0h_pDox_A$rn), "venn_t24h_pDox_A_vs_t0h_pDox_A_t48h_pDox_A_vs_t0h_pDox_A_not72h_pDox_A_vs_t0h_pDox_A_up.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

write.table(setdiff(intersect(t48h_pDox_A_vs_t0h_pDox_A$rn,t72h_pDox_A_vs_t0h_pDox_A$rn),t24h_pDox_A_vs_t0h_pDox_A$rn), "venn_not24h_pDox_A_vs_t0h_pDox_A_t48h_pDox_A_vs_t0h_pDox_A_t72h_pDox_A_vs_t0h_pDox_A_up.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

write.table(setdiff(intersect(t24h_pDox_A_vs_t0h_pDox_A$rn,t72h_pDox_A_vs_t0h_pDox_A$rn),t48h_pDox_A_vs_t0h_pDox_A$rn), "venn_t24h_pDox_A_vs_t0h_pDox_A_not48h_pDox_A_vs_t0h_pDox_A_t72h_pDox_A_vs_t0h_pDox_A_up.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)

write.table(intersect(intersect(t24h_pDox_A_vs_t0h_pDox_A$rn,t48h_pDox_A_vs_t0h_pDox_A$rn),t72h_pDox_A_vs_t0h_pDox_A$rn), "venn_t24h_pDox_A_vs_t0h_pDox_A_t48h_pDox_A_vs_t0h_pDox_A_t72h_pDox_A_vs_t0h_pDox_A_up.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
```



## uniprot keyword analyses

Four separate lists were created from `211117 List of Proteins for Keywords Analysis.xlsx` and then do uniprot keyword analysis similar following to what performed on [this manuscript](https://github.com/sblab-bioinformatics/cmpp)


### combine columns and upload to uniprot to retrieve info

```r
library(data.table)


# Enlarge the view width when printing tables
options(width = 300)

data <- fread("211117_List_of_Proteins_for_Keywords_Analysis.csv", skip=1)

write.table(setdiff(unique(data$V1), c("")), "211117_List_of_Proteins_for_Keywords_Analysis_nonkinetic.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(unique(data$V3), c("")), "211117_List_of_Proteins_for_Keywords_Analysis_kinetic_down.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(unique(data$V6), c("")), "211117_List_of_Proteins_for_Keywords_Analysis_kinetic_up.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(unique(data$V8), c("")), "211117_List_of_Proteins_for_Keywords_Analysis_kinetic_down_up.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
write.table(setdiff(unique(data$V10), c("")), "211117_List_of_Proteins_for_Keywords_Analysis_kinetic_down_up_72honly.txt", row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
```

I used [uniprot](https://www.uniprot.org/) ID mapping to obtain the uniprot keywords and go terms to result in:

- `211117_List_of_Proteins_for_Keywords_Analysis_nonkinetic_fixed_uniprot.txt`
- `211117_List_of_Proteins_for_Keywords_Analysis_kinetic_down_fixed_uniprot.txt`
- `211117_List_of_Proteins_for_Keywords_Analysis_kinetic_up_fixed_uniprot.txt`
- `211117_List_of_Proteins_for_Keywords_Analysis_kinetic_down_up_fixed_uniprot.txt` (used in downstream analyses)
- `211117_List_of_Proteins_for_Keywords_Analysis_kinetic_down_up_72honly_fixed_uniprot.txt`


### barplots biological function

prepare files:

```python

### 211111_List_of_Proteins_for_Keywords_Analysis_kinetic_down_up_fixed_uniprot.txt
ifile = open("211111_List_of_Proteins_for_Keywords_Analysis_kinetic_down_up_fixed_uniprot.txt", "r")
ilines = ifile.readlines()
ifile.close()

uniprot_keyword = []
uniprot_go = []

for l in ilines[1:]:
  fields = l.replace("\n", "").split("\t")
  uniprot = fields[0]
  gene_name = fields[1]
  ks = fields[2].split(";")
  gos = fields[4].split("; ")
  for k in ks:
    uniprot_keyword.append((uniprot, k))
  for go in gos:
    go_m = go.split(" [")[0]
    uniprot_go.append((uniprot, go_m))

ofile_k = open("211111_List_of_Proteins_for_Keywords_Analysis_kinetic_down_up_fixed_uniprot_keyword.txt", "w")
ofile_k.write("\n".join(["\t".join(e) for e in uniprot_keyword]))
ofile_k.close()

ofile_g = open("211111_List_of_Proteins_for_Keywords_Analysis_kinetic_down_up_fixed_uniprot_go.txt", "w")
ofile_g.write("\n".join(["\t".join(e) for e in uniprot_go]))
ofile_g.close()
```

barplots uniprot keyword:

```bash
library(data.table)
library(ggplot2)

### 211111_List_of_Proteins_for_Keywords_Analysis_kinetic_down_up_fixed_uniprot_keyword.txt
data_keyword <- fread("211111_List_of_Proteins_for_Keywords_Analysis_kinetic_down_up_fixed_uniprot_keyword.txt", header = FALSE)
setnames(data_keyword, c("uniprot", "keyword"))

data.table(sort(table(data_keyword$keyword), decreasing = T)[1:50])
# 5:      Alternative splicing 31
#18:                 Transport 15
#20:                    Signal 12
#29:             Transcription  7
#36:           Differentiation  5
#39:             Ion transport  5
#44:        Mental retardation  4
#47:        Biological rhythms  3
#49:             Cell adhesion  3
#50:           Cell projection  3

data_keyword_table <- data.table(sort(table(data_keyword$keyword), decreasing = T)[1:50])[V1 %in% c("Alternative splicing", "Transport", "Signal", "Transcription", "Differentiation", "Ion transport", "Mental retardation", "Biological rhythms", "Cell adhesion", "Cell projection")]

gg <- ggplot(data = data_keyword_table, aes(x = N, y = reorder(V1, N))) +
geom_bar(stat="identity") +
xlab("Number of proteins") +
ylab("") +
ggtitle("UniprotKB keyword") +
theme_bw() +
theme(axis.title = element_text(size=16), axis.text = element_text(size=12, color = "black"), axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1), plot.title = element_text(face="bold", size=16, hjust = 0.5))

ggsave("211111_List_of_Proteins_for_Keywords_Analysis_kinetic_down_up_fixed_uniprot_keyword.pdf", height = 4, width = 6)
```

For interpretation of uniprot keywords check this: https://www.uniprot.org/docs/keywlist
