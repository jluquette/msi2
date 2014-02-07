# Command line arguments:
#   output_file input_file1,sampleid1 [ input_file2,sampleid2 ... ]

if (length(commandArgs(TRUE)) < 2) {
    stop("usage: join_candidates.R output_file input_file1,sampleid1 [ input_file2,sampleid2 ... ]")
}

out.file = commandArgs(TRUE)[1]

# Columns that must exist in the file, used for joining
join_cols = c("chr", "start", "end", "ref_len", "unit", "region", "flank1", "flank2", "sequence")
data_cols = c("raw_alleles", "call", "genotype", "pval", "allele_summaries")

reader <- function(file) {
    cat(sprintf("Reading %s...\n", file))
    read.table(file, sep="\t", header=TRUE)
}

# Make an empty data.frame with the expected columns
df <- do.call(data.frame, as.list(c(join_cols, data_cols)))  # Dummy row
colnames(df) <- c(join_cols, data_cols)
df <- df[FALSE,]  # Remove the row of dummy values


all <- TRUE
#for (sample in samples) {
for (arg in commandArgs(TRUE)[-1]) {
    x <- strsplit(arg, ',')[[1]]
    if (length(x) != 2) {
        stop("input arguments must be of the form file,sampleid, with no spaces near the comma")
    }
    in.file <- x[1]
    sample <- x[2]
    t <- reader(in.file)
    df <- merge(df, t, sort=FALSE, by=join_cols, all=all,
                suffixes=c('', paste('.', sample, sep='')))
    all <- FALSE
}

# All of the data columns without a .XXX suffix are the original dummy
# data columns, which are all NA.
for (col in data_cols) {
    df[col] <- NULL
}

cat(sprintf("Writing combined output to %s..", out.file))
write.table(df, file=out.file, quote=FALSE, sep="\t", row.names=FALSE)
