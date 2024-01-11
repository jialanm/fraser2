IMPLEMENTATION = "PCA"


def count_split_reads_single_sample_r(num_of_cpu, sample_id, bam_path):
    return f"""
    library(data.table)
    library(FRASER)
    library(BiocParallel)

    if(.Platform$OS.type == "unix") {{
        register(MulticoreParam(workers=min({num_of_cpu}, multicoreWorkers())))
    }} else {{
        register(SnowParam(workers=min({num_of_cpu}, multicoreWorkers())))
    }}

    # create sample dataset for FRASER with columns of sampleID and bamFile
    annotation_dat <- data.table(sampleID="{sample_id}",
                                 bamFile="{bam_path}",
                                 pairedEnd=TRUE)
    print(annotation_dat)

    # fds <- FraserDataSet(colData=annotation_dat, workingDir=".")
    # 
    # # count split reads
    # split_counts <- getSplitReadCountsForAllSamples(fds, BPPARAM=bpparam())

    settings <- FraserDataSet(colData=annotation_dat, workingDir=".")
    # count the split and non-split reads
    fds <- countRNAData(settings, BPPARAM=bpparam())
    """


def count_reads_all_samples_r(num_of_cpu):
    return f"""
    library(data.table)
    library(FRASER)
    library(BiocParallel)

    args <- commandArgs(trailingOnly = TRUE)
    sample_ids <- unlist(strsplit(args[1], split=":"))
    bam_paths <- unlist(strsplit(args[2], split=":"))

    # create sample dataset for FRASER with columns of sampleID and bamFile
    annotation_dat <- data.table(sampleID=sample_ids,
                                 bamFile=bam_paths,
                                 pairedEnd=TRUE)
    print(annotation_dat)

    settings <- FraserDataSet(colData=annotation_dat, workingDir=".")

    if(.Platform$OS.type == "unix") {{
        register(MulticoreParam(workers=min({num_of_cpu}, multicoreWorkers())))
    }} else {{
        register(SnowParam(workers=min({num_of_cpu}, multicoreWorkers())))
    }}

    # count the split and non-split reads
    fds <- countRNAData(settings, BPPARAM = bpparam())

    # calculate Jaccard intron index and other metrics
    fds <- calculatePSIValues(fds, BPPARAM=bpparam())

    saveFraserDataSet(fds, dir=".")
    """


def run_fraser_r(psitype, num_of_cpu, result_table_filename, heatmap_before_ae,
                 heatmap_after_ae, enc_dim_auc,
                 enc_dim_loss, aberrant, delta_psi_threshold, padj_threshold, min_reads,
                 gene_models_gff_path):
    return f"""
    library(data.table)
    library(FRASER)
    library(ggplot2)
    library(BiocParallel)
    library(GenomicFeatures)
    library(org.Hs.eg.db)

# library(TxDb.Hsapiens.UCSC.hg38.knownGene)

    args <- commandArgs(trailingOnly = TRUE)
    sample_ids <- unlist(strsplit(args[1], split=":"))

    if(.Platform$OS.type == "unix") {{
        register(MulticoreParam(workers=min({num_of_cpu}, multicoreWorkers())))
    }} else {{
        register(SnowParam(workers=min({num_of_cpu}, multicoreWorkers())))
    }}


    fds = loadFraserDataSet(".")
    num_of_iter = 15

    print(fds)

    # change splice metrics
    fitMetrics(fds) <- "{psitype}"  # not available in FRASER1

    # plot color heatmap for samples before autoencoder correction for sample covariance
    # before_ae <- plotCountCorHeatmap(fds, type="{psitype}", logit=TRUE, plotType="sampleCorrelation")
    # ggsave(filename = "{heatmap_before_ae}", plot = before_ae, device = "png")

    # filter junctions with low expressions
    fds <- filterExpressionAndVariability(fds,  minExpressionInOneSample={min_reads}, minDeltaPsi={delta_psi_threshold}, filter=TRUE, BPPARAM=bpparam())

    # annotate with genes
    txdb_obj <- makeTxDbFromGFF("./{gene_models_gff_path}")
    fds <- annotateRangesWithTxDb(fds, txdb=txdb_obj, feature="ENSEMBL", keytype="ENSEMBL")
    # fds <- annotateRangesWithTxDb(fds, txdb=txdb_obj, keytype="ENSEMBL")
    print(fds)

    # get the optimal dimension of the latent space
    fds <- optimHyperParams(fds, type="{psitype}", implementation="{IMPLEMENTATION}", BPPARAM=bpparam())
    best_q = bestQ(fds, type="{psitype}")
    # best_q = 9
    print("Best Q is: ")
    print(best_q)

    # plot the encoding dimension search auc and loss
    # enc_auc = plotEncDimSearch(fds, type="{psitype}", plotType="auc") 
    # ggsave(filename = "{enc_dim_auc}", plot = enc_auc, device="png")
    # enc_loss = plotEncDimSearch(fds, type="{psitype}", plotType="loss") 
    # ggsave(filename = "{enc_dim_loss}", plot = enc_loss, device="png")

    # run FRASER pipeline
    print("Run FRASER pipeline... ")
    print(fds)

    fds <- FRASER(fds, q=best_q, type="{psitype}", implementation="{IMPLEMENTATION}", iterations=num_of_iter, BPPARAM=bpparam())
    fds <- calculateZscore(fds, type="{psitype}")
    print(fds)

    # plot heatmap after confounder correction
    # after_ae <- plotCountCorHeatmap(fds, type="{psitype}", logit=TRUE, normalized=TRUE, plotType="sampleCorrelation")
    # ggsave(filename = "{heatmap_after_ae}", plot = after_ae, device = "png")

    # aberrant <- plotAberrantPerSample(fds, type="{psitype}", sample_id, deltaPsiCutoff={delta_psi_threshold}, padjCutoff={padj_threshold})
    # ggsave(filename = "{aberrant}", plot = aberrant, device = "png")

    res_filtered <- as.data.table(results(fds, padjCutoff={padj_threshold}, 
    deltaPsiCutoff={delta_psi_threshold}))
    res <- as.data.table(results(fds, padjCutoff=1, deltaPsiCutoff=0))
    # res <- as.data.table(results(fds, padjCutoff={padj_threshold}, deltaPsiCutoff=0))
    # res <- as.data.table(results(fds))
    print(res)

    write.table(res_filtered, file="filtered_{result_table_filename}", quote=FALSE, row.names=FALSE)
    write.table(res, file="{result_table_filename}", quote=FALSE, row.names=FALSE)

    # result visualization
    for(sample_id in sample_ids) {{
        cur_plot_name <- paste("volcano_", sample_id, ".png", sep="")
        p <- plotVolcano(fds, sampleID=sample_id, type="{psitype}", deltaPsiCutoff = {delta_psi_threshold}, padjCutoff = {padj_threshold})
        ggsave(filename = cur_plot_name, plot = p, device = "png")
    }}
    """