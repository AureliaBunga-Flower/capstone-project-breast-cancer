# 1. Load Library
library(GEOquery)
library(limma)

# 2. Ambil Data (Ganti GSE ID sesuai kebutuhan)
gse_id <- "GSE15852" 
gse <- getGEO(gse_id, destdir = ".", getGPL = TRUE)
if (length(gse) > 1) idx <- grep("GPL", attr(gse, "names")) else idx <- 1
gse <- gse[[idx]]

# 3. Ekstrak Matriks Ekspresi dan Data Klinis (Metadata)
exprs_mat <- exprs(gse)
pheno <- pData(gse)

# --- BAGIAN KRUSIAL: Identifikasi Grup ---
# Kita cari kolom mana yang berisi informasi "tumor" atau "normal"
# Biasanya ada di kolom 'title' atau 'characteristics_ch1'
target_col <- "title" # Sesuaikan jika kolomnya berbeda

pheno$Group <- "Unknown"
# Gunakan grepl (L kecil), ignore.case = TRUE agar tidak sensitif huruf besar/kecil
pheno$Group[grepl("normal", pheno[[target_col]], ignore.case = TRUE)] <- "Normal"
pheno$Group[grepl("tumor|cancer|carcinoma", pheno[[target_col]], ignore.case = TRUE)] <- "Tumor"

# Filter hanya yang kita kenali (Normal & Tumor)
keep <- pheno$Group != "Unknown"
pheno <- pheno[keep, ]
exprs_mat <- exprs_mat[, keep]

# Ubah menjadi Factor (Wajib agar tidak error contrasts)
pheno$Group <- factor(pheno$Group, levels = c("Normal", "Tumor"))

# Cek jumlah sampel (Pastikan muncul di konsol)
cat("Statistik Sampel:\n")
print(table(pheno$Group))

# 4. Normalisasi Log2 (Jika data belum di-log)
qx <- as.numeric(quantile(exprs_mat, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  exprs_mat[exprs_mat <= 0] <- NaN
  exprs_mat <- log2(exprs_mat) 
  cat("Data telah di-log2 transformasi.\n")
}

# 5. Differential Expression Analysis (LIMMA)
design <- model.matrix(~0 + pheno$Group)
colnames(design) <- c("Normal", "Tumor")

fit <- lmFit(exprs_mat, design)
# Membuat kontras: Tumor vs Normal
cont.matrix <- makeContrasts(Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 6. Ambil Hasil
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
head(tT)

# Visualisasi Sederhana: Volcano Plot
volcanoplot(fit2, main="Volcano Plot: Tumor vs Normal")