# Runs every plugin of a release's native library through ReVAMP on the 91 pilot recordings of the Acoustic Index
# User's Guide, at its defaults, and compares each value with the R function's value in the 0.3.0 benchmark
# (reference.csv). A value passes within 1e-4 of the R value, relative to it, plus 1e-6 absolute, as soundecology
# rounds some indices to six decimal places.
#   Rscript check_release.R <recordings folder> <results csv> [label]
# VAMP_PATH must name the folder holding the library. Exits 1 if a plugin is not found, a recording is not the
# benchmark's (md5) or a value does not pass.
suppressMessages(library(ReVAMP))
a <- commandArgs(TRUE)
dir <- a[1]; out <- a[2]; label <- if (length(a) >= 3) a[3] else R.version$platform
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
ref <- read.csv(file.path(here, "reference.csv"), colClasses = c(md5 = "character"), stringsAsFactors = FALSE)
REL <- 1e-4; ABS <- 1e-6
KEY <- function(p) paste0("vamp-ecoacoustics:", p)

cat(label, "\nReVAMP", format(packageVersion("ReVAMP")), "on", R.version.string, "\nVamp path:", vampPaths(), "\n")
installed <- vampPlugins()
missing <- setdiff(KEY(unique(ref$plugin)), installed$id)
if (length(missing)) stop("not found by ReVAMP: ", paste(missing, collapse = ", "))

wavs <- list.files(dir, pattern = "[.]wav$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
sums <- unname(tools::md5sum(wavs))
path_of <- function(f, md5) {
  hit <- wavs[basename(wavs) == f & sums == md5]
  if (!length(hit)) stop(f, " (md5 ", md5, ") is not in ", dir)
  hit[1]
}
first_value <- function(path, plugin, output) {
  v <- tryCatch(as.numeric(unlist(runPlugin(path, KEY(plugin))[[output]]$value))[1], error = function(e) NA_real_)
  if (length(v)) v else NA_real_
}

d <- do.call(rbind, lapply(split(ref, ref$file), function(r) {
  p <- path_of(r$file[1], r$md5[1])
  data.frame(host = label, file = r$file, plugin = r$plugin, output = r$output,
             value = mapply(first_value, p, r$plugin, r$output, USE.NAMES = FALSE),
             r_value = r$r_value, benchmark_value = r$benchmark_value)
}))
d$diff_r <- d$value - d$r_value
d$diff_benchmark <- d$value - d$benchmark_value
d$ok <- !is.na(d$value) & abs(d$diff_r) <= REL * abs(d$r_value) + ABS
write.csv(d, out, row.names = FALSE)

s <- do.call(rbind, lapply(unique(ref$plugin), function(p) {
  x <- d[d$plugin == p, ]
  rel <- abs(x$diff_r) / abs(x$r_value)
  data.frame(plugin = p, recordings = nrow(x), failed = sum(!x$ok), max_abs_diff_r = max(abs(x$diff_r)),
             max_rel_diff_r = max(rel[x$r_value != 0]), max_abs_diff_benchmark = max(abs(x$diff_benchmark)))
}))
print(s, row.names = FALSE, digits = 3)
md <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(md))
  cat(c(paste("###", label), "", paste("ReVAMP", format(packageVersion("ReVAMP")), "on", R.version.string), "",
        "| plugin | recordings | failed | largest difference from R | largest relative difference from R | largest difference from the benchmark build |",
        "|---|---|---|---|---|---|",
        sprintf("| %s | %d | %d | %.3g | %.3g | %.3g |", s$plugin, s$recordings, s$failed, s$max_abs_diff_r,
                s$max_rel_diff_r, s$max_abs_diff_benchmark), ""),
      file = md, sep = "\n", append = TRUE)
if (any(!d$ok)) {
  cat(sum(!d$ok), "of", nrow(d), "values outside the tolerance\n")
  quit(status = 1)
}
cat("all", nrow(d), "values within the tolerance\n")
