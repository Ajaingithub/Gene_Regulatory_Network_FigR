# Directory where data will be downloaded to
# workingDir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/FigR/"
# shareseqZip <- "https://s3.us-east-1.amazonaws.com/vkartha/FigR/FigR_SHAREseq.zip"
# download.file(url = shareseqZip, destfile = paste0(workingDir,basename(shareseqZip)))
# unzip(paste0(workingDir,basename(shareseqZip)),exdir = workingDir,overwrite = FALSE)
# setwd(paste0(workingDir))

# We require scATAC-seq peak counts as a SummarizedExperiment object, where the rowRanges are fixed-width accessibility peak ranges (typically 300 bp
# or so), and the assay object is a sparseMatrix of Tn5-instertion counts per peak per cell. Additionally, the scRNA-seq matrix can be provided as a 
# sparseMatrix of the log-normalized gene expression levels for the same cells (e.g. processed using Seurat).

# In this example, we show these two inputs pertaining to a SHARE-seq experiment performed on mouse skin cells. For simplicity, we only run this on a
# subset of cells (quality filtered, and sampled to n=10,000 cells for faster run-time). As a result, note that the results will be similar to what 
# was presented in Kartha et al. , but not identical.

FigR_regulome <- function(obj, savedir, objname, peak_cutoff = 5){
  library(ggplot2)
  library(dplyr)
  library(FigR)
  library(SummarizedExperiment)
  library(Seurat)
  library(BuenColors)
  rest_act <- obj
  RNAmat <- rest_act@assays$RNA@data
 
 message("Making a summarized experiment")
 counts <- rest_act@assays$ATAC@counts
 rowRanges <- rest_act@assays$ATAC@ranges
 colData <- rest_act@meta.data
 ATAC_se <- SummarizedExperiment(assays=list(counts=counts),rowRanges=rowRanges, colData=colData)
 wnn_umap <- as.data.frame(rest_act@reductions$wnn.umap@cell.embeddings)
 all(rownames(colData(ATAC_se)) == rownames(wnn_umap))
 colData(ATAC_se)$wnnUMAP_1 <- wnn_umap$wnnUMAP_1
 colData(ATAC_se)$wnnUMAP_2 <- wnn_umap$wnnUMAP_2
 message("ATAC dimension")
 print(dim(ATAC_se)) # Peaks x Cells
 
 message("Remove genes with zero expression across all cells")
 RNAmat <- RNAmat[Matrix::rowSums(RNAmat)!=0,]
 
 dir.create(paste(savedir,"saveRDS",sep = ""), showWarnings = FALSE)
 saveRDS(RNAmat,paste(savedir,"saveRDS/RNAmat.RDS",sep = ""))
 
 message("RNA dimension")
 print(dim(RNAmat)) # Genes x Cells
 
 # Additionally, we load an matrix object of cisTopic probability scores that was generated as part of the dimensionality reduction process implemented
 # on the scATAC-seq count data using the lsi performed in Seurat. This will enable us to derive a cell kNN matrix (k-nearest neighbors graph) that is used
 # for smoothing the data. Note, if you used an alternative method for dimensionality reduction, you can derive this cell kNN matrix using alternative 
 # inputs (e.g. using the top N LSI / PCA components)
 lsi_ATAC <- rest_act@reductions$lsi@cell.embeddings
 stopifnot(all(rownames(lsi_ATAC)==rownames(colData(ATAC_se)))) # if TRUE move forward
 
 message("Deriving cell kNN")
 set.seed(123)
 cellkNN <- FNN::get.knn(lsi_ATAC,k = 30)$nn.index
 dim(cellkNN)
 
 stopifnot(all(rownames(colData(ATAC_se))==rownames(rest_act@meta.data)))
 colData(ATAC_se)$cluster <- rest_act@meta.data$seurat_clusters
 
 saveRDS(ATAC_se,paste(savedir,"saveRDS/ATAC_se.RDS",sep = ""))
 
 dir.create(paste(savedir,"Figures",sep = ""),showWarnings = FALSE)
 pdf(paste(savedir,"Figures/FigR_UMAP_ATAC.pdf",sep = ""),width = 7, height = 6)
 print(colData(ATAC_se) %>% as.data.frame() %>% ggplot(aes(wnnUMAP_1,wnnUMAP_2,color=cluster)) + 
   geom_point(size=1) +
   theme_classic())
 dev.off()
 
 
 message("Running Peak Gene Correlation. This will take some time")
 # Peak-gene association testing
 # Next, we take a default (10 kb) window around each gene's TSS, and compute the Spearman correlation across all cells between their peak accessibility
 # counts (mean-centered) and the normalized RNA expression. For each peak-gene pair correlation, we use (default n=100) background peaks matched for
 # GC content and accessibility to correlate to the same gene, so that we can test for significance (permutation p-value).
 # Don't run interactively
 cisCorr <- FigR::runGenePeakcorr(ATAC.se = ATAC_se,
                                  RNAmat = RNAmat,
                                  genome = "hg38", # One of hg19, mm10 or hg38 
                                  nCores = 8,
                                  p.cut = NULL, # Set this to NULL and we can filter later
                                  n_bg = 100)
 print(head(cisCorr))
 
 dir.create(paste(savedir,"saveRDS",sep = ""),showWarnings = FALSE)
 saveRDS(cisCorr,paste(savedir,"saveRDS/cisCorr.RDS",sep = ""))
 
 message("Filtering the cis correlated with pvalue < 0.05")
 cisCorr.filt <- cisCorr %>% dplyr::filter(pvalZ <= 0.05)
 
 dir.create(paste(savedir,"Figures",sep = ""),showWarnings = FALSE)
 pdf(paste(savedir,"Figures/DORC_correlated.pdf",sep = ""), width = 6, height = 6)
 print(dorcGenes <- dorcJPlot(dorcTab = cisCorr.filt,
                        cutoff = peak_cutoff, # No. sig peaks needed to be called a DORC
                        labelTop = 30,
                        returnGeneList = TRUE, # Set this to FALSE for just the plot
                        force=5))
 dev.off()
 
 dir.create(paste(savedir,"saveRDS",sep = ""),showWarnings = FALSE)
 saveRDS(dorcGenes,paste(savedir,"saveRDS/dorcGenes.RDS",sep = ""))
 
 
 # We can also just get the full ranked table of genes and number of significantly associated peaks, instead of only getting the gene names (so you can threshold yourself)
 # Unfiltered
 numDorcs <- cisCorr.filt %>% group_by(Gene) %>% tally() %>% arrange(desc(n))
 numDorcs_df <- as.data.frame(numDorcs)
 
 dir.create(paste(savedir,"saveRDS",sep = ""),showWarnings = FALSE)
 saveRDS(numDorcs_df,paste(savedir,"saveRDS/numDorcs_df.RDS",sep = ""))
 
 dorcMat <- getDORCScores(ATAC.se = ATAC_se, # Has to be same SE as used in previous step
                          dorcTab = cisCorr.filt,
                          geneList = dorcGenes,
                          nCores = 8)
 dir.create(paste(savedir,"saveRDS",sep = ""),showWarnings = FALSE)
 saveRDS(dorcMat,paste(savedir,"saveRDS/dorcMat.RDS",sep = ""))
 
 ## It will get the DorcMat only for the dorcGene which were significant
 # We can then smooth these (sparse) DORC counts, and visualize side-by-side with the expression, on the same UMAP. In this particular case, we use the 
 # original Latent Semantic Indexing (LSI) cell embeddings that were pre-determined using all cells (using which the 2D UMAP was determined), subset to 
 # the cells being used, and then derive cell k-nearest neighbors, which we will use to smooth both the DORC and paired RNA expression, for the same UMAP:
 # Smooth dorc scores using cell KNNs (k=30)
 rownames(cellkNN) <- rownames(lsi_ATAC)
 dir.create(paste(savedir,"saveRDS",sep = ""),showWarnings = FALSE)
 saveRDS(cellkNN,paste(savedir,"saveRDS/cellkNN.RDS",sep = ""))
 
 dorcMat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = dorcMat,nCores = 8)
 # Smooth RNA using cell KNNs
 # This takes longer since it's all genes
 RNAmat.s <- smoothScoresNN(NNmat = cellkNN[,1:30],mat = RNAmat,nCores = 8)
 
 saveRDS(dorcMat.s,paste(savedir,"saveRDS/dorcMat_smooth.RDS",sep = ""))
 saveRDS(RNAmat.s,paste(savedir,"saveRDS/RNAmat_smooth.RDS",sep = ""))
 
 # # Visualize on pre-computed UMAP
 # umap.d <- as.data.frame(colData(ATAC_se)[,c("wnnUMAP_1","wnnUMAP_2")])
 # 
 # # DORC score for Dlx3
 # dorcg <- plotMarker2D(umap.d,dorcMat.s,markers = c("PREX1","IFNG"),maxCutoff = "q0.99",colorPalette = "brewer_heat")
 # 
 # #> Plotting  Dlx3
 # 
 # # RNA for Dlx3
 # rnag <- plotMarker2D(umap.d,RNAmat.s,markers = c("PREX1","IFNG"),maxCutoff = "q0.99",colorPalette = "brewer_purple")
 # 
 # #> Plotting  Dlx3
 # 
 # library(patchwork)
 # pdf(paste(savedir,"PREX1_IFNG_dorc_rna.pdf",sep = ""), width = 10, height = 5)
 # dorcg + rnag
 # dev.off()
 
 message("Running TF-gene associations")
 # 
 # The core component of FigR is to determine TFs that are putative regulators (acivators or repressors) of DORCs. By specifying a built-in reference
 # motif database (see Methods of associated manuscript), and providing smoothed DORC accessibility and RNA count matrices as input, we determine 
 # which DORCs are enriched for different TF binding motifs, in addition to testing their correlations to TF RNA expression. This can be done using a
 # single command in FigR
 # For 708 TFs against 586 DORCs (testing with n=50 background iterations), this took ~ 1.5 hrs using 4 cores (for the same 10k cells).
 
 figR.d <- runFigRGRN(ATAC.se = ATAC_se, 
                      dorcTab = cisCorr.filt, 
                      genome = "hg38",
                      dorcMat = dorcMat.s,
                      rnaMat = RNAmat.s, 
                      nCores = 8)
 saveRDS(figR.d,paste(savedir,"saveRDS/FigR.d.RDS",sep = ""))
 
 message("Making Global regulation profile")
 # First, we can simply plot all TF-DORC associations (all tested TFs by all tested DORCs) as a scatter plot, where each point is a TF-DORC pair. The 
 # x-axis provides the significance of the correlation of each TF's expression to the DORC accessibility (Z-test relative to background), and the
 # y-axis shows the relative enrichment of the TF's motif among the DORC-associated peaks. Points are colored by their determined FigR regulation 
 # score. Putative activating associations are highlighted in red, while repressive associations are in blue

 # figR.d_filt <- figR.d[figR.d$Corr.log10P < 10,]
 
 pdf(paste(savedir,"Figures/global_regulation_profile.pdf",sep = ""))
 print(figR.d %>% 
   ggplot(aes(Corr.log10P,Enrichment.log10P,color=Score)) + 
   ggrastr::geom_point_rast(size=0.01,shape=16) + 
   theme_classic() + 
   scale_color_gradientn(colours = jdb_palette("solar_extra"),limits=c(-3,3),oob = scales::squish,breaks=scales::breaks_pretty(n=3)))
 dev.off()
 
 # # We can rank TFs by the overall (mean) regulation score across all DORCs. This helps give a snapshot of the major activators and repressors across all cells in the data given the DORCs
 pdf(paste(savedir,"Figures/Rank_TF_based_on_meanscore.pdf",sep = ""), width = 10, height = 6)
 print(rankDrivers(figR.d, rankBy = "meanScore", interactive = FALSE))
 dev.off()

 # Specify the myLabels parameter to label the plots above for select TFs (and not just those in the top and bottom percentiles by default). Useful to see where your TFs of interest fall among all TFs
 pdf(paste(savedir,"Figures/Rank_TF_based_on_ntargets.pdf",sep = ""), width = 10, height = 6)
 print(rankDrivers(figR.d, score.cut = 1, rankBy = "nTargets", interactive = FALSE))
 dev.off()
 
 # plotDrivers(figR.d,score.cut = 2,marker = "PREX1")
 # plotDrivers(figR.d,score.cut = 2,marker = "CD3Z")
 
 library(ComplexHeatmap)
 pdf(paste(savedir,"Figures/heatmap_TF_and_DORC.pdf",sep = ""), width = 8, height = 6)
 print(plotfigRHeatmap(figR.d = figR.d,
                 score.cut = 1.5,
                 column_names_gp = grid::gpar(fontsize=8), # from ComplexHeatmap
                 show_row_dend = FALSE # from ComplexHeatmap
 ))
 dev.off()
 
 # library(session)
 # save.session(paste(savedir,objname,".Rda",sep = ""))
 # 
 # library(networkD3)
 # FigR::plotfigRNetwork(figR.d,
 #                       score.cut = 1.5,
 #                       TFs = c("TBX21"),
 #                       weight.edges = TRUE)
 
 # figR.d_filt <- figR.d[figR.d$Score < -1.5 | figR.d$Score > 1.5,]
 # figR.d_filt_selected <- dplyr::select(figR.d_filt, c("DORC","Motif","Score"))
 # figR.d_filt_selected_recast <- reshape::cast(figR.d_filt_selected, DORC ~ Motif)
 # figR.d_filt_selected_recast[is.na(figR.d_filt_selected_recast)] <- 0
 # 
 # library(circlize)
 # col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
 # col_fun(seq(-3, 3))
 # Heatmap(figR.d_filt_selected_recast, name = "mat", col = col_fun)
 
}



