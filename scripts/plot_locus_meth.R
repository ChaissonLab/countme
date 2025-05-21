
library(ggplot2)


library(readr)

# Load your file
file_path <- "region.cons.txt"  # Replace with actual path
df <- read_tsv(file_path, col_names = FALSE)

region <- read_file("region.txt")



df <- read.table(file_path, sep = "\t", stringsAsFactors = FALSE, header = FALSE)
blue_to_red <- function(val) {
  stopifnot(val >= 0 && val <= 1)  # enforce bounds

  # Define gradient from blue (low) to red (high)
  color_fn <- colorRamp(c("blue", "red"))
  
  rgb_vals <- color_fn(val)  # returns RGB in [0, 255]
  rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], maxColorValue = 255)
}

parse_triplets <- function(triplet_str, line_id) {
  triplets <- strsplit(triplet_str, "/")[[1]]
  if (length(triplets) <= 1) return(NULL)

  triplets <- triplets[-length(triplets)]  # remove the last triplet

  dat <- lapply(triplets, function(t) {
    vals <- as.numeric(strsplit(t, ",")[[1]])
    if (length(vals) != 3 || vals[2] < 3) return(NULL)  # exclude if count < 3
    data.frame(
      line_id = line_id,
      x = vals[1],
      count = vals[2],
      score = vals[3]
    )
  })

  dat <- do.call(rbind, dat)
  if (is.null(dat)) return(NULL)
  return(dat)
}

# Parse all rows
parsed_list <- lapply(1:nrow(df), function(i) parse_triplets(df[i, 5], i))
all_points <- do.call(rbind, parsed_list)


max_x_values <- sapply(parsed_list, function(df) if (!is.null(df)) max(df$x) else NA)
max_x <- max(max_x_values[!is.na(max_x_values)])
xOrder <- order(max_x_values)
pdf("tr_meth_plot.pdf")
plot(c(), xlim=c(0,max_x*1.1), ylim=c(0, length(xOrder)), ylab="Sample", xlab="Locus position", main=region)

for (j in seq(1, length(xOrder), by=10)) {
    i <- xOrder[j]
    cols <- sapply(parsed_list[i][[1]]$score, function(i) blue_to_red(i))
    points(parsed_list[i][[1]]$x, y= rep(j, length(parsed_list[i][[1]]$x)), pch=21,bg=cols)
}
dev.off()


