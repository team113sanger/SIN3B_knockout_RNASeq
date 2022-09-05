us <- read.table('/nfs/users/nfs_a/ad33/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression/transcript-counts/summarised/2581-feature-count-v103.txt', sep='\t', header=TRUE, row.names=1)
casm <- read.table('/nfs/users/nfs_a/ad33/projects/6406_SIN3B_role_in_melanoma_resistance_and_progression/transcript-counts-casm/summarised/2581-feature-count-v103.txt', sep='\t', header=TRUE, row.names=1)

stopifnot(
    "rowname mismatch"=all(rownames(us) == rownames(casm)),
    "colname mismatch"=all(colnames(us) == colnames(casm))
)

samples <- colnames(us)


correlations <- do.call(rbind.data.frame, lapply(sample(rownames(us), 25), function(g){
    x <- as.numeric(us[g, samples])
    y <- as.numeric(casm[g, samples])
    return(data.frame(
        "gene" = g,
        "cor" = cor(x, y)
    ))
}))
