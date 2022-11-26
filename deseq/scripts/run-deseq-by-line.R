# This script loads the input counts data into RDS format for downstream analysis
# Alastair Droop, 2022-08-08

# Load the necessary libraries:
library(tidyverse)
library(DESeq2)
library(GenomicRanges)
library(BiocParallel)

# Set the run IO:
io <- list(
    "counts" = file.path("..", "..", "transcript-counts", "summarised", "2581-feature-count-v103.txt"),
    "metacounts" = file.path("..", "..", "transcript-counts", "summarised", "2581-metafeature-count-v103.txt"),
    "metadata" = file.path("..", "..", "metadata", "2581-metadata.txt"),
    "gene_metadata" = file.path("..", "..", "metadata", "genes", "gene-length-GRCh38_103.txt"),
    "chr_map" = file.path("..", "..", "metadata", "genes", "chr-map.txt"),
    "plot_dir" = file.path("..", "data", "by-line", "qc"),
    "rdata_dir" = file.path("..", "data", "by-line", "rdata"),
    "results_dir" = file.path("..", "data", "by-line", "results"),
    "gsea_dir" = file.path("..", "data", "by-line", "gsea")
)

# Set the run parameters:
params <- list(
    "n_cpu" = 8, # The number of CPUs to use
    "dataset_prefix" = "2581-by-line", # This is the dataset prefix to use when saving these data to RDS format
    "count_threshold" = 10, # This is the threshold count for a gene inclusion summed across all samples
    "invalid_clones" = c("C18", "C21"), # This is the list of clones to remove from the data *before* analysis 
    "test_genes" = c("SIN3B", "MYC", "BMI1"), # The gene names to plot
    "coefficient" = "condition_knockout_vs_control"
)

# Register the CPUs:
register(MulticoreParam(params$n_cpu))

# Load & tidy the metadata:
metadata <- read.table(io$metadata, sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE)
metadata$line <- factor(metadata$line)
metadata$clone <- factor(metadata$clone, levels=sort(unique(metadata$clone)))
metadata$replicate <- factor(metadata$replicate, levels=sort(unique(metadata$replicate)))
metadata$condition <- factor(c('CONTROL'='control', 'SIN3B_KO'='knockout')[metadata$condition], levels=c('control', 'knockout'))

# Pre-load the counts data to test for counts:
counts <- read.table(io$counts, sep="\t", header=TRUE, row.names=1, nrows=1)

# Check the valid samples:
valid_samples <- sort(intersect(rownames(metadata), colnames(counts)))
stopifnot(
    "count & metadata sample mismatch" = length(setdiff(valid_samples, colnames(counts))) == 0,
    "count & metadata sample mismatch" = length(setdiff(valid_samples, rownames(metadata))) == 0
)

# Load the complete count & metacount data as a matrix of integers:
counts <- as.matrix(read.table(io$counts, sep="\t", header=TRUE, row.names=1, colClasses=c("character", rep("integer", length(valid_samples)))))
metacounts <- as.matrix(read.table(io$metacounts, sep="\t", header=TRUE, row.names=1, colClasses=c("character", rep("integer", length(valid_samples)))))
stopifnot("count & metacount sample mismatch" = all(colnames(counts) == colnames(metacounts)))

# Load the gene metadata:
valid_chrs <- read.table(io$chr_map, sep="\t", header=FALSE)[,2]
genes <- read.table(io$gene_metadata, sep="\t", header=TRUE, row.names=1, na.strings="-")
genes$chr <- factor(genes$chr, levels=valid_chrs)
genes$biotype <- factor(genes$biotype)

# Remove the gene set inclusions:
genes$MSK <- NULL
genes$ACE <- NULL

# Make sure the genes match:
stopifnot(
    "count & metadata gene mismatch" = length(setdiff(rownames(genes), rownames(counts))) == 0,
    "count & metadata gene mismatch" = length(setdiff(rownames(counts), rownames(genes))) == 0
)

# Remove genes without valid positions:
valid_genes <- sort(intersect(rownames(counts), rownames(genes[!is.na(genes$chr) & !is.na(genes$start) & !is.na(genes$end), ])))

# Reorder metadata & counts to ensure matching rows & columns:
metadata <- metadata[valid_samples,]
counts <- counts[valid_genes, valid_samples]
metacounts <- metacounts[, valid_samples]
genes <- genes[valid_genes, ]

# Remove any invalid clones:
deseq_samples <- rownames(metadata[metadata$clone %in% params$invalid_clones == FALSE,])

# Iterate through the lines, performing individual analyses on each one:
for(line in levels(metadata$line)) {
    message(sprintf("processing line %s...", line))
    # Split the metadata:
    line_metadata <- metadata[metadata$line == line, ]
    line_metadata$line  <- droplevels(line_metadata$line)
    # Get the valid samples:
    line_samples <- intersect(valid_samples, rownames(line_metadata))
    # Get the valid counts:
    line_counts <- counts[valid_genes, line_samples]
    # Get the valid metacounts:
    line_metacounts <- metacounts[, line_samples]
    # Build the DESeq2 object:
    deseq_d <- DESeqDataSetFromMatrix(
        countData = line_counts[, line_samples],
        colData = line_metadata[line_samples, ],
        rowData = makeGRangesFromDataFrame(genes, keep.extra.columns=TRUE),
        design = ~condition # This is the model we're fitting
    )
    # Clip to only genes with sufficient data:
    deseq_d <- deseq_d[rowSums(counts(deseq_d)) >= params$count_threshold,]
    # Perform DESeq processing:
    deseq_d <- DESeq(deseq_d, parallel=TRUE)
    # Make some pre-analysis plots:
    for(gene in params$test_genes) {
        gene_id <- rownames(genes[genes$gene == gene,])
        stopifnot("non-unique gene name" = length(gene_id) == 1)
        # Plot the raw counts:
        d <- data.frame("count"=line_counts[gene_id,])
        d$sample <- rownames(d)
        d$line <- line_metadata[rownames(d), "line"]
        d$condition <- line_metadata[rownames(d), "condition"]
        g <- ggplot(d, aes(x=sample, y=count, fill=condition, group=line))
        g <- g + geom_col()
        g <- g + theme(axis.ticks.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=4)) 
        g <- g + theme(panel.background=element_blank(), panel.border=element_rect(fill=NA))
        g <- g + labs(y="Raw Count", x="")
        g <- g + ggtitle(sprintf("%s Raw Counts (%s)", gene, gene_id))
        ggsave(g, file=file.path(io$plot_dir, sprintf("%s-%s-raw_count.pdf", line, gene)), width=12, height=4)
        # Plot the DESeq normalised data:
        d <- plotCounts(deseq_d, gene=gene_id, intgroup="condition", returnData=TRUE)
        g <- ggplot(d, aes(x=condition, y=count, colour=condition))
        g <- g + geom_point(position=position_jitter(w=0.1,h=0), show.legend=FALSE)
        g <- g + theme(axis.ticks.x=element_blank()) 
        g <- g + theme(panel.background=element_blank(), panel.border=element_rect(fill=NA))
        g <- g + labs(y="Normalised Count", x="")
        g <- g + scale_y_log10()
        g <- g + ggtitle(sprintf("%s Normalised Counts", gene))
        ggsave(g, file=file.path(io$plot_dir, sprintf("%s-%s-norm_count.pdf", line, gene)), width=4, height=6)
    }
    # Load the metacounts into a ggplot-able form:
    d_metacounts <- as.data.frame(line_metacounts)
    d_metacounts$metacount <- factor(rownames(d_metacounts))
    d_metacounts <- pivot_longer(data=d_metacounts, cols=-metacount, names_to="sample", values_to="count")
    d_metacounts$line <- droplevels(line_metadata[d_metacounts$sample, "line"])
    d_metacounts$condition <- line_metadata[d_metacounts$sample, "condition"]
    # Plot the metacounts:
    g <- ggplot(d_metacounts, aes(x=sample, y=count, fill=condition))
    g <- g + geom_col()
    g <- g + theme(axis.ticks.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, size=4)) 
    g <- g + theme(panel.background=element_blank(), panel.border=element_rect(fill=NA))
    g <- g + facet_grid(rows=vars(metacount))
    g <- g + labs(y="Raw Count", x="")
    ggsave(g, file=file.path(io$plot_dir, sprintf("%s-%s-metacounts.pdf", params$dataset_prefix, line)), width=12, height=8)
    # Save the data to RDS format:
    saveRDS(line_metadata, file=file.path(io$rdata_dir, sprintf("%s-%s-metadata.data", params$dataset_prefix, line)))
    saveRDS(genes, file=file.path(io$rdata_dir, sprintf("%s-%s-genes.data", params$dataset_prefix, line)))
    saveRDS(line_counts, file=file.path(io$rdata_dir, sprintf("%s-%s-counts.data", params$dataset_prefix, line)))
    saveRDS(line_metacounts, file=file.path(io$rdata_dir, sprintf("%s-%s-metacounts.data", params$dataset_prefix, line)))
    saveRDS(deseq_d, file=file.path(io$rdata_dir, sprintf("%s-%s-deseq.data", params$dataset_prefix, line)))
    # Save the DESeq2-normalised and scaled data to GSEA format:
    ## 1: Make the GCT file:
    deseq_counts <- counts(deseq_d, normalized=TRUE, replaced=TRUE)
    gct_file = file(file.path(io$gsea_dir, sprintf("%s-%s-normalized-counts.gct", params$dataset_prefix, line)), open="wt")
    writeLines("#1.2", con=gct_file)
    writeLines(sprintf("%d\t%d", nrow(deseq_counts), ncol(deseq_counts)), con=gct_file)
    writeLines(paste(c("Name", "Description", colnames(deseq_counts)), collapse="\t"), con=gct_file)
    write.table(cbind.data.frame(rownames(deseq_counts), genes[rownames(deseq_counts), 'gene'], deseq_counts), file=gct_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, na="")
    close(gct_file)
    ## 2: Make the CLS file:
    cls_file = file(file.path(io$gsea_dir, sprintf("%s-%s-normalized-counts.cls", params$dataset_prefix, line)), open="wt")
    writeLines(sprintf("%d\t%d\t1", nrow(colData(deseq_d)), nlevels(colData(deseq_d)$condition)), con=cls_file)
    writeLines(sprintf("#%s", paste(levels(colData(deseq_d)$condition), collapse="\t")), con=cls_file)
    writeLines(paste(as.character(colData(deseq_d)$condition), collapse="\t"), con=cls_file)
    close(cls_file)
    ## 3: Make and return the annotated results data:
    deseq_results <- lfcShrink(deseq_d, coef=params$coefficient, saveCols=c("gene", "biotype"), type="apeglm", parallel=TRUE)
    deseq_results <- deseq_results[order(deseq_results$pvalue, decreasing=FALSE, na.last=TRUE), ]
    write.csv(deseq_results, file=file.path(io$results_dir, sprintf("%s-%s-deseq-results.csv", params$dataset_prefix, line)))
}
