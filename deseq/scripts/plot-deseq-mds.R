# This script loads the DESeq-2 processed & raw counts and plots them using MDS

# Set the run IO:
io <- list(
    "input" = file.path("..", "data", "all-lines", "rdata", "2581-all-lines-counts.data"),
    "output_dir" = file.path("..", "data", "all-lines", "qc")
)

# Define colours for the lines:
colours <- c(
    "A2058" = "#66C2A5", 
    "A375" = "#FC8D62", 
    "SKMEL28" = "#8DA0CB"
)

# Load the raw counts data:
d <- readRDS(io$input)

# Calculate the MDS data for the raw counts:
sample_dist <- dist(t(d))
sample_mds <- as.data.frame(cmdscale(sample_dist))
sample_mds$col <- colours[gsub('^([^_]+)_.*', '\\1', rownames(sample_mds), perl=TRUE)]
pdf(file.path(io$output_dir, "2581-all-lines-raw-count-mds.pdf"), width=12, height=12)
plot(sample_mds[,1], sample_mds[,2], axes=FALSE, ann=FALSE, pch=16, col=sample_mds$col)
box()
text(sample_mds[,1], sample_mds[,2], labels=rownames(sample_mds), cex=0.5)
dev.off()

# Now, look at the data using ranked correlations. This focusses on the information in the signal, not the features of the signal itself.
dr <- cor(d, method='spearman')
# We need to convert a correlation into something that looks like a distance.  See https://stats.stackexchange.com/questions/165194/using-correlation-as-distance-metric-for-hierarchical-clustering
drc <- sqrt(2 * (1 - dr))
# We now need to make that thing look like a distance to R:
drc <- as.dist(drc)
# We can now plot the MDS just as before:
sample_mds <- as.data.frame(cmdscale(drc))
sample_mds$col <- colours[gsub('^([^_]+)_.*', '\\1', rownames(sample_mds), perl=TRUE)]
pdf(file.path(io$output_dir, "2581-all-lines-raw-count-mds-spearman.pdf"), width=12, height=12)
plot(sample_mds[,1], sample_mds[,2], axes=FALSE, ann=FALSE, pch=16, col=sample_mds$col)
box()
dev.off()
