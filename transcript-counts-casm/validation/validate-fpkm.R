# This script loads the FPKM & TPM counts from the DERMATLAS pipeline and from canapps to compare them.
# Alastair Droop, 202-07-12

# Set the run arguments:
args <- list(
    'ctypes' = c('count', 'tpm', 'fpkm'),
    'samples' = file.path('..', '..', 'metadata', '2581-samples.txt'),
    'canapps.dir' = file.path('.', 'canapps-fpkm'),
    'dermatlas.input.base' = file.path('..', 'summarised', '2581-feature-%s-v103.txt'),
    'results' = file.path('.', 'validation-results.txt')
)

# Load the samples:
samples <- readLines(args$samples)

# Load the canapps input filenames:
c.filenames <- list.files(args$canapps.dir, pattern='*-converted-counts.txt')
names(c.filenames) <- gsub('-converted-counts.txt', '', c.filenames, fixed=TRUE)

# Open the output file:
res.file <- file(args$results, open='wt')

# Iterate through the count types:
for(ctype in args$ctypes){
    message(sprintf('testing %s...', ctype))
    # Load the DERMATLAS data:
    d.filename <- sprintf(args$dermatlas.input.base, ctype)
    d.counts <- read.table(d.filename, sep='\t', header=TRUE, row.names=1)

    # Get the list of valid samples:
    valid.samples <- intersect(colnames(d.counts), names(c.filenames))
    valid.samples <- sort(intersect(samples, valid.samples))
    message(sprintf('  %d/%d samples are valid', length(valid.samples), length(samples)))

    # Go through each sample in turn:
    max.delta <- 0
    for(sample in valid.samples){
        # Open the appropriate canapps data file:
        c.counts <- read.table(file.path(args$canapps.dir, c.filenames[sample]), sep='\t', header=FALSE, row.names=1)
        colnames(c.counts) <- c('gene', 'biotype', 'chr', 'longest_isoform', 'count', 'unfiltered_count', 'fpkm', 'fpkm_uq', 'tpm')

        # Pull out a list of matching genes, and extract matching data from DERMATLAS & canapps:
        valid.genes <- intersect(rownames(d.counts), rownames(c.counts))
        d.data <- d.counts[valid.genes, sample]
        c.data <- c.counts[valid.genes, ctype]

        # Calculate the maximum absolute difference between the two:
        data.diff <- abs(d.data - c.data)
        i <- which.max(data.diff)
        max.delta <- max(max.delta, max(data.diff))
        message(sprintf('  %s %s: delta_max = %f (%f : %f)', ctype, sample, max(data.diff), d.data[i], c.data[i]))
    }
    res.str <- sprintf('overall maximum delta: %f', max.delta)
    message(sprintf('  %s %s', ctype, res.str))
    writeLines(sprintf('%s\t%s', ctype, res.str), con=res.file)
}

# Close the output file:
close(res.file)
