# capstone-project-breast-cancer
Analisis Transcriptomics GSE15852: Breast Cancer vs Normal Tissue
Repository ini berisi proyek akhir (Capstone Project) untuk analisis data ekspresi gen pada kanker payudara menggunakan dataset publik GSE15852. Proyek ini mencakup alur kerja bioinformatika lengkap, mulai dari pra-pemrosesan data hingga interpretasi jalur biologis.

 Deskripsi Proyek
Tujuan utama dari proyek ini adalah mengidentifikasi Differentially Expressed Genes (DEGs) yang secara signifikan membedakan jaringan tumor payudara dengan jaringan payudara normal. Melalui pendekatan statistik dan fungsional, proyek ini memetakan perubahan molekuler yang terjadi pada sel kanker.

 Metodologi & Alat
Analisis dilakukan menggunakan bahasa pemrograman R (v4.5.2) dengan library utama sebagai berikut:

limma: Digunakan untuk analisis ekspresi diferensial menggunakan model linear.

clusterProfiler: Digunakan untuk analisis pengayaan fungsional (GO & KEGG).

ggplot2 & enrichplot: Digunakan untuk visualisasi data berkualitas tinggi.

UMAP: Digunakan untuk reduksi dimensi dan visualisasi pengelompokan sampel.

Hasil Utama
1. Pengelompokan Sampel (Clustering)
Berdasarkan plot UMAP, terdapat pemisahan yang sangat jelas antara kelompok Normal dan Cancer. Hal ini menunjukkan bahwa dataset memiliki profil transkriptomik yang konsisten dan kontras di antara kedua kondisi.

2. Identifikasi DEGs
Ditemukan sebanyak 1.775 gen yang terekspresi secara diferensial dengan kriteria adj.P.Val < 0.05.

Up-regulated: Termasuk gen penanda kanker seperti ERBB2 dan MMP9.

Down-regulated: Termasuk gen pengatur metabolisme seperti ADIPOQ dan LEP.

3. Analisis Jalur (Pathway Analysis)
Hasil analisis KEGG menunjukkan gangguan signifikan pada PPAR Signaling Pathway. Jalur ini berperan dalam metabolisme lipid dan energi, yang mengindikasikan adanya pemrograman ulang metabolik pada sel tumor payudara.
