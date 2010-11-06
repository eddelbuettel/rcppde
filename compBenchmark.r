#!/usr/bin/r -t

svnver <- system("svnversion", intern=TRUE)
cat("# compiled benchmark at SVN", svnver, "\n")
source("demo/CompiledBenchmark.R")
