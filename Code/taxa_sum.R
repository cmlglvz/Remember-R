taxa.sum <- function(otu.table, taxa.table, taxa.lvl ){
  z <- NULL
  y <- NULL
  for (i in 1:length(unique(taxa.table[colnames(otu.table),taxa.lvl]))) {
    if (length(otu.table[,which(taxa.table[colnames(otu.table),taxa.lvl]==unique(taxa.table[colnames(otu.table),taxa.lvl])[i])])!=length(rownames(otu.table))) {
      z <- which(taxa.table[colnames(otu.table),taxa.lvl]==unique(taxa.table[colnames(otu.table),taxa.lvl])[i])
      y <- cbind(y,apply(otu.table[,which(taxa.table[colnames(otu.table),taxa.lvl]==unique(taxa.table[colnames(otu.table),taxa.lvl])[i])],1,function(x) sum(x)))
    } else {
      y <- cbind(y,otu.table[,which(taxa.table[colnames(otu.table),taxa.lvl]==unique(taxa.table[colnames(otu.table),taxa.lvl])[i])])
    }
  }
  colnames(y) <- unique(taxa.table[colnames(otu.table),taxa.lvl])
  invisible((y))
}