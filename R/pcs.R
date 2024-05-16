#' Two axes of genetic variation.
#'
#' @description
#'  The data frame contains two PCs. In this simulation example,we considered the 7 populations
#' from African ancestry which were created from Phase 3 of 1000 Genomes, incorporating two populations
#' from America continent and five populations from African continent. To derive the axes of genetic
#' variation, in each reference panel, we used a subset of genetic variants with minor allele frequency (MAF) >5\%
#' in all populations, which were separated by 1Mb, to derive the matrix of pairwise Euclidean distances
#' between the populations. Then the two axes of genetic variation would be derived from the multi-dimensional
#' scaling of the distance matrix in the corresponding reference panel.
#' @source 7 populations of African ancestry were collected from Phase 3 of 1000 Genomes \url{https://ctg.cncr.nl/software/MAGMA/ref_data/}.
"pcs"
