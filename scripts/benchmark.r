#!/usr/bin/r -t

svnver <- system("svnversion", intern=TRUE)
cat("# small benchmark at SVN", svnver, "\n")
source("demo/SmallBenchmark.R")
