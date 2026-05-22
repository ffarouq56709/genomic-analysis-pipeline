# Load the required libraries
suppressPackageStartupMessages(library(karyoploteR))
suppressPackageStartupMessages(library(GenomicRanges))

patients <- c("PAPKHU", "PASJJH", "PASNSF", "PATNBV")

# Map T2T accessions ONLY for autosomes (1-22)
# NC_060925.1 (Chr1) to NC_060946.1 (Chr22)
chrom_map <- c(
  "NC_060925.1" = "chr1", "NC_060926.1" = "chr2", "NC_060927.1" = "chr3", "NC_060928.1" = "chr4",
  "NC_060929.1" = "chr5", "NC_060930.1" = "chr6", "NC_060931.1" = "chr7", "NC_060932.1" = "chr8",
  "NC_060933.1" = "chr9", "NC_060934.1" = "chr10", "NC_060935.1" = "chr11", "NC_060936.1" = "chr12",
  "NC_060937.1" = "chr13", "NC_060938.1" = "chr14", "NC_060939.1" = "chr15", "NC_060940.1" = "chr16",
  "NC_060941.1" = "chr17", "NC_060942.1" = "chr18", "NC_060943.1" = "chr19", "NC_060944.1" = "chr20",
  "NC_060945.1" = "chr21", "NC_060946.1" = "chr22"
)

# Define the autosome set for plotting
autosomes_to_plot <- paste0("chr", 1:22)

for (patient in patients) {
  filename <- paste0(patient, "_100kb_depth.regions.bed.gz")
  
  if (!file.exists(filename)) {
    message(paste("Skipping", patient, "- file not found"))
    next
  }
  
  message(paste("Processing", patient, "..."))
  
  # 1. Load the Data
  df <- read.table(gzfile(filename), sep="\t", col.names=c("chrom", "start", "end", "depth"))
  
  # 2. Map Chromosomes (Filtering out X and Y here)
  df$chrom <- chrom_map[df$chrom]
  df <- df[!is.na(df$chrom), ]
  
  # 3. Calculate Global Log2 Baseline Math
  med_depth <- median(df$depth[df$depth > 0], na.rm=TRUE)
  df$log2_ratio <- log2((df$depth + 1) / med_depth)
  
  # Filter noise outliers
  df <- df[df$log2_ratio > -2.5 & df$log2_ratio < 2.5, ]
  
  # 4. Convert to GenomicRanges
  gr <- GRanges(seqnames = df$chrom, 
                ranges = IRanges(start = df$start, end = df$end), 
                log2_ratio = df$log2_ratio)
  
  # Assign Colors
  dot_cols <- ifelse(gr$log2_ratio > 0.58, "#e74c3c", 
              ifelse(gr$log2_ratio < -1.0, "#3498db", "#95a5a6"))
  
  # ==========================================
  # PLOT 1: Whole Genome (Autosomes Only)
  # ==========================================
  wg_out <- paste0(patient, "_WG_Autosomes_karyo.png")
  png(wg_out, width=2000, height=700, res=150)
  
  # chromosomes=autosomes_to_plot removes X and Y
  # main= adds the title above the graph
  kp <- plotKaryotype(genome="hg38", 
                      plot.type=4, 
                      chromosomes=autosomes_to_plot, 
                      labels.plotter=NULL, 
                      main=paste("Whole Genome CNV Profile (Autosomes) -", patient))
  
  kpAddChromosomeNames(kp, srt=45, cex=0.9)
  kpAxis(kp, ymin=-2, ymax=2, numticks=5, cex=0.8)
  
  # Add baselines and plot dots
  kpAbline(kp, h=0, col="black", ymin=-2, ymax=2)
  kpAbline(kp, h=c(0.58, -1), col=c("#e74c3c", "#3498db"), lty=2, ymin=-2, ymax=2)
  kpPoints(kp, data=gr, y=gr$log2_ratio, ymin=-2, ymax=2, col=dot_cols, cex=0.4)
  
  dev.off()
  
  # ==========================================
  # PLOT 2: Chromosome 21 Detail
  # ==========================================
  chr21_out <- paste0(patient, "_Chr21_Detail_karyo.png")
  png(chr21_out, width=1200, height=550, res=150)
  
  kp2 <- plotKaryotype(genome="hg38", 
                       chromosomes="chr21", 
                       main=paste("Chromosome 21 iAMP21 Landscape -", patient))
  
  kpDataBackground(kp2, data.panel=1, color="#f8f9fa") 
  kpAxis(kp2, ymin=-2, ymax=2, numticks=5)
  
  # Add baseline and plot dots
  kpAbline(kp2, h=0, col="black", ymin=-2, ymax=2)
  kpPoints(kp2, data=gr, y=gr$log2_ratio, ymin=-2, ymax=2, col=dot_cols, cex=0.8)
  
  dev.off()
}

message("Done! Check your folder for the polished graphics.")
