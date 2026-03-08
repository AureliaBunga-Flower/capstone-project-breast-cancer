#library() digunakan agar fungsi di dalam package bisa digunakan 
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(hgu133a.db)
library(AnnotationDbi)
library(umap)
gset <- getGEO("GSE15852", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
#ExpressionSet berisi:
# - exprs() : matriks ekspresi gen
# - pData() : metadata sampel
# - fData() : metadata fitur (probe / gen)

#PART D. PRE-PROCESSING DATA EKSPRESI 

# exprs(): mengambil matriks ekspresi gen
# Baris  = probe/gen
# Kolom  = sampel
ex <- exprs(gset)
#quantile(): menghitung nilai kuantil (persentil)
#as.numeric(): mengubah hasil quantile (yang berupa named vector)
#menjadi vektor numerik biasa agar mudah dibandingkan
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

#LogTransform adalah variabel logika (TRUE / FALSE)
#Operator logika:
#>  : lebih besar dari
#| | : OR (atau)
#&& : AND (dan)
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

#IF statement:
#Jika LogTransform = TRUE, maka lakukan log2
if (LogTransform) {
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}


#PART E. DEFINISI KELOMPOK SAMPEL 

#pData(): metadata sampel
#source_name_ch1 berisi informasi kondisi biologis sampel
group_info <- pData(gset)[["source_name_ch1"]]

#make.names(): mengubah teks menjadi format valid untuk R
groups <- make.names(group_info)

#factor():
#Mengubah data kategorik menjadi faktor
#Faktor sangat penting untuk analisis statistik di R
gset$group <- factor(groups)

#levels(): melihat kategori unik dalam faktor
nama_grup <- levels(gset$group)
print(nama_grup)

#PART F. DESIGN MATRIX (KERANGKA STATISTIK) 

#model.matrix():
#Membuat matriks desain untuk model linear
#~0 berarti TANPA intercept (best practice limma)
design <- model.matrix(~0 + gset$group)

#colnames(): memberi nama kolom agar mudah dibaca
colnames(design) <- levels(gset$group)

#Menentukan perbandingan biologis
grup_kanker <- nama_grup[1]
grup_normal <- nama_grup[2]

contrast_formula <- paste(grup_kanker, "-", grup_normal)
print(paste("Kontras yang dianalisis:", contrast_formula))

#PART G. ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)

#lmFit():
#Membangun model linear untuk setiap gen
fit <- lmFit(ex, design)

#makeContrasts(): mendefinisikan perbandingan antar grup
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

#contrasts.fit(): menerapkan kontras ke model
fit2 <- contrasts.fit(fit, contrast_matrix)

#eBayes():
#Empirical Bayes untuk menstabilkan estimasi varians
fit2 <- eBayes(fit2)

#topTable():
#Mengambil hasil akhir DEG
#adjust = "fdr" -> koreksi multiple testing
#p.value = 0.01  -> gen sangat signifikan
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.01
)

head(topTableResults)



#PART H. ANOTASI NAMA GEN 

#Mengambil ID probe dari hasil DEG
probe_ids <- rownames(topTableResults)

#Mapping probe -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(
  hgu133a.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

#Gabungkan dengan hasil limma
topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

#Cek hasil anotasi
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])

#PART I.1 BOXPLOT DISTRIBUSI NILAI EKSPRESI 

#Boxplot digunakan untuk:
#- Mengecek distribusi nilai ekspresi antar sampel
#- Melihat apakah ada batch effect
#- Mengevaluasi apakah normalisasi/log-transform sudah wajar

#Set warna berdasarkan grup
group_colors <- as.numeric(gset$group)

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)


#PART I.2 DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT) 

#Density plot menunjukkan sebaran global nilai ekspresi gen
#Digunakan untuk:
#- Mengecek efek log-transform
#- Membandingkan distribusi antar grup

#Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )


#PART I.3 UMAP (VISUALISASI DIMENSI RENDAH)

#Transpose matriks ekspresi:
#UMAP bekerja pada OBSERVATION = sampel
umap_input <- t(ex)

#Jalankan UMAP
umap_result <- umap(umap_input)

#Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )


#PART J.1 VISUALISASI VOLCANO PLOT 

#Volcano plot menggabungkan:
#- Log fold change (efek biologis)
#- Signifikansi statistik

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

#Klasifikasi status gen
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.01] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.01] <- "DOWN"

#Visualisasi
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG Kanker Payudara")


#PART J.2 VISUALISASI HEATMAP 

#Heatmap digunakan untuk melihat pola ekspresi gen
#antar sampel berdasarkan gen-gen paling signifikan

#Pilih 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)

#Ambil matriks ekspresi untuk gen terpilih
mat_heatmap <- ex[top50$PROBEID, ]

#Gunakan Gene Symbol (fallback ke Probe ID)
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,      # jika SYMBOL kosong → probe ID
  top50$SYMBOL        # jika ada → gene symbol
)

rownames(mat_heatmap) <- gene_label

#Pembersihan data (WAJIB agar tidak error hclust)
#Hapus baris dengan NA
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

#Hapus gen dengan varians nol
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

#Anotasi kolom (kelompok sampel)
annotation_col <- data.frame(
  Group = gset$group
)

rownames(annotation_col) <- colnames(mat_heatmap)

#Visualisasi heatmap 
pheatmap(
  mat_heatmap,
  scale = "row",                 # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE,         # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)


#PART K. MENYIMPAN HASIL 

# write.csv(): menyimpan hasil analisis ke file CSV
write.csv(topTableResults, "Hasil_GSE15852_DEG.csv")

message("Analisis selesai. File hasil telah disimpan.")

# Ambil gen yang signifikan (misal logFC > 1 atau < -1)
genes_to_test <- tT$Gene.symbol[tT$adj.P.Val < 0.05]
genes_to_test <- topTableResults$Gene.symbol[topTableResults$adj.P.Val < 0.05]

# Install paket jika belum ada
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

library(clusterProfiler)
library(org.Hs.eg.db) # Database untuk manusia (Homo sapiens)
library(enrichplot)

# Konversi Symbol ke Entrez ID
entrez_ids <- bitr(genes_to_test, fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

go_enrich <- enrichGO(gene         = entrez_ids$ENTREZID,
                      OrgDb        = org.Hs.eg.db,
                      ont          = "BP", # Bisa diganti "MF" atau "CC"
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      readable      = TRUE)

# Visualisasi sederhana
dotplot(go_enrich, showCategory=20) + ggtitle("GO Enrichment - Biological Process")
kegg_enrich <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                          organism     = 'hsa', # 'hsa' untuk Homo sapiens
                          pvalueCutoff  = 0.05)

# Visualisasi
dotplot(kegg_enrich, showCategory=20) + ggtitle("KEGG Pathway Analysis")