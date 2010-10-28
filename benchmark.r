#!/usr/bin/r -t

svnver <- system("svnversion", intern=TRUE)
cat("# At", format(Sys.time()), "\n# SVN ", svnver, "\n", sep="")
source("demo/SmallBenchmark.R")
