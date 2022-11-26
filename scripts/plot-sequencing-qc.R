#! /usr/bin/env Rscript-base
# This script tests a set of samples against the thresholds defined in a JSON file

# Set the version:
.version <- '1.0.0'

# A function to write a given message to stderr:
log.message <- function(..., verbose=NA){
    if(is.na(verbose)){
        verb <- get0('.verbose', ifnotfound=TRUE)
    } else {
        verb <- verbose
    }
    if(identical(verb, TRUE)){
        message(sprintf(...))
        flush(stderr())
    }
}

# A function to exit with a given error:
error <- function(..., exit.code=1, verbose=TRUE){
    log.message(sprintf('ERROR: %s', sprintf(...)), verbose=verbose)
    if(!identical(interactive(), TRUE)){
        quit(save='no', status=exit.code)
    }
}

# A function to quietly load a vector libraries from character strings:
loadLibrary <- function(x, verbose=NA){
    for(l in x){
        log.message('loading library "%s"', l, verbose=verbose)
        res <- suppressWarnings(suppressPackageStartupMessages(require(l, character.only=TRUE, quietly=TRUE)))
        if(!identical(res, TRUE)) error('failed to load package "%s"', l)
    }
}

# Create pretty scientific labels:
# https://stackoverflow.com/questions/10762287/how-can-i-format-axis-labels-with-exponents-with-ggplot2-and-scales
scientific <- function(x){
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scales:::scientific_format()(x)))))
}

# A function to extract a series of single column values for each sample into ggplot format:
collateSingle <- function(x, sel){
    if(is.null(names(sel))){
        if(length(sel) > 1){
            names(sel) <- sprintf('value_%d', 1:length(sel))
        } else {
            names(sel) <- 'value'
        }
    }
    res <- data.frame(
        'sample' = factor(rownames(x), levels=rownames(x))
    )
    for(i in names(sel)){
        res[[i]] <- eval(parse(text=sel[i]), envir=x)
    }
    return(res)
}

# A function to return a pass/fail factor from a verctor of values and a threshold:
checkThreshold <- function(x, threshold){
    res <- factor(rep('pass', length(x)), levels=c('pass', 'fail'))
    if(!is.na(threshold[1])) res[x < threshold[1]] <- 'fail'
    if(!is.na(threshold[2])) res[x > threshold[2]] <- 'fail'
    return(res)
}

# A function to plot a single total with threshold:
plotSingleValueThreshold <- function(filename, x, col, threshold, ylab, title, width=50, axis=waiver(), colours=cols, ret.full=FALSE) {
    if(is.null(threshold)) return(invisible(NULL))
    # Build the ggplot input data:
    gd <- x[, c('sample', col)]
    gd$status <- checkThreshold(gd[[col]], threshold)
    #  Define the status colours:
    status.cols <- c('pass'=colours$pass, 'fail'=colours$fail)
    # Define the range:
    ylim <- range(gd[[col]])
    # Build the plot:
    g <- ggplot(gd, aes_string(x='sample', y=col, fill='status'))
    g <- g + geom_col(show.legend=FALSE)
    g <- g + theme(axis.ticks=element_blank(), panel.grid=element_blank(), panel.background=element_blank(), panel.border=element_blank())
    g <- g + theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, vjust=0.5, size=4, family='mono'), plot.subtitle=element_text(colour='grey25'), legend.title=element_blank())
    g <- g + scale_y_continuous(labels=axis)
    g <- g + ylab(ylab)
    g <- g + ggtitle(label=title)
    if(!is.na(threshold[1])){
        g <- g + geom_hline(yintercept=threshold[1], colour=colours$threshold)
        ylim <- c(min(ylim[1], threshold[1]), ylim[2])
    }
    if(!is.na(threshold[2])){
        g <- g + geom_hline(yintercept=threshold[2], colour=colours$threshold)
        ylim <- c(ylim[1], max(ylim[2], threshold[2]))
    }
    g <- g + coord_cartesian(ylim=ylim, expand=TRUE)
    g <- g + scale_fill_discrete(type=status.cols)
    ggsave(filename, plot=g, width=width, height=12, units='cm')
    if(!identical(ret.full, TRUE)){
        gd <- gd$status
    }
    return(invisible(gd))
}

# Build & process the CLI:
loadLibrary('argparse', verbose=FALSE)
# parser <- ArgumentParser(description='Generate QC plots from a bamstats data file')
# parser$add_argument('-V', '--version', dest='version', default=FALSE, action='store_true', help='print version information')
# parser$add_argument('-v', '--verbose', dest='verbose', default=FALSE, action='store_true', help='provide verbose output')
# parser$add_argument('-p', '--output-prefix', dest='output_prefix', metavar='prefix', default='', help='output plot prefix')
# parser$add_argument(dest='stats_file', metavar='file', help='sample statistics file')
# parser$add_argument(dest='threshold_file', metavar='file', help='threshold JSON file')
# parser$add_argument(dest='output_dir', metavar='dir', default='.', help='output plot folder')
# args <- parser$parse_args()
args <- list(
    'verbose' = TRUE,
    'version' = FALSE,
    'output_prefix' = 'res-', 
    'stats_file' = '/nfs/users/nfs_a/ad33/locations/scratch119/6406_SIN3B_role_in_melanoma_resistance_and_progression/qc/sequencing/summarised/2581-bamstats.txt',
    'threshold_file' = '/nfs/users/nfs_a/ad33/locations/scratch119/6406_SIN3B_role_in_melanoma_resistance_and_progression/qc/sequencing/2581-qc-thresholds.json',
    'output_dir' = '/nfs/users/nfs_a/ad33/locations/scratch119/6406_SIN3B_role_in_melanoma_resistance_and_progression/qc/sequencing/results'
)
.verbose <- args$verbose

# Process the version data if required:
if(identical(args$version, TRUE)){
    log.message('plot-bamstats %s', .version, verbose=TRUE)
    if(identical(.verbose, TRUE)){
        log.message('\nsession info:', verbose=TRUE)
        sink(stderr())
        print(sessionInfo())
        sink()
    }
    quit(save='no', status=0)
}

# Load the necessary libraries:
loadLibrary(c('jsonlite', 'ggplot2'))

# Read the threshold data:
# NB: The data are not checked here for consistency!
log.message('reading QC threshold data from "%s"', args$threshold_file)
thresholds <- jsonlite::fromJSON(args$threshold_file)

# Load the input data:
log.message('reading bamstats data from "%s"', args$stats_file)
d <- read.table(
    args$stats_file, 
    sep='\t', 
    header=TRUE,  
    comment.char='', 
    check.names=FALSE, 
    colClasses=c(rep('character', 1), rep('numeric', 16))
)
# colnames(d) <- gsub('#', 'n', colnames(d), fixed=TRUE) # Replace pound signs with n (for automatic expression parsing)
# rownames(d) <- sprintf('%s_%s', d$sample, d$platform_unit) # Make sure we have unique names

# Set the plot colours:
cols <- list(
    'pass' = '#66BD63',
    'fail' = '#F46D43',
    'threshold' = '#D73027'
)

# Convert ratios to percentages:
for(i in c('mapped_pairs', 'duplication')){
    d[[i]] <- d[[i]] * 100
}

# Check the thresholds present in the JSON file:
plotSingleValueThreshold(file.path(args$output_dir, sprintf('%sreads-total-count.pdf', args$output_prefix)), d, col='reads', threshold=thresholds$reads$total, ylab='Total Reads', title='Total Sample Reads', axis=scientific)
plotSingleValueThreshold(file.path(args$output_dir, sprintf('%sreads-mapped-percent.pdf', args$output_prefix)), d, col='mapped_pairs', threshold=thresholds$reads$mapped.percent, ylab='Percent of Total Reads', title='Sample Read Mapping Rate')
plotSingleValueThreshold(file.path(args$output_dir, sprintf('%sreads-duplicated-percent.pdf', args$output_prefix)), d, col='duplication', threshold=thresholds$reads$duplicate.percent, ylab='Percent of Total Reads', title='Read Duplication Rate')
plotSingleValueThreshold(file.path(args$output_dir, sprintf('%sreads-mean-insertsize.pdf', args$output_prefix)), d, col='insert_mean', threshold=thresholds$reads$insert.size, ylab='Mean Insert Size', title='Sample Mean Insert Size')
