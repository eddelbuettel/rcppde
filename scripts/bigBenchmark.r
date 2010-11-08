#!/usr/bin/r -t

svnver <- system("svnversion", intern=TRUE)
cat("# big benchmark at SVN", svnver, "\n")
source("demo/LargeBenchmark.R")
