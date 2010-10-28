#!/usr/bin/r -t

svnver <- system("svnversion", intern=TRUE)
cat("# SVN", svnver, "\n")
source("demo/LargeBenchmark.R")
