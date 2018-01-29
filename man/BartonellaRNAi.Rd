\name{BartonellaRNAi}
\alias{D}
\alias{KOmap}
\alias{siRNA_gene_Expts}
\docType{data}
\title{Perturbation with off targets}
\description{
  The dataset consists of binary data with 288 features and 35 experiments, a perturbation map for 35 experiments targetting 8 genes, and a dataframe of gene names corresponding to the siRNA IDS. 
}
\usage{data(BartonellaRNAi)}
\format{
  D: binary data matrix: 288 x 35\\
  KOmap: perturbation matrix:  35 x 8\\
  siRNA_gene_Expts:		data frame: 35 x 3\\
}
\details{
  The datawas derived from microscopy image-based infection assays where ATCC HeLa cells were transfected with a genome-wide single-siRNA
  library from Qiagen followed by infection with B. henselae. The cells were then fixed, stained and imaged. The images were corrected for 
  illumination distortion using CIDRE (Smith et al., 2015). Subsequently, the cell features (phenotypic effects) were extracted from the grid of 9
  images per knockdown experiment using screeningBee CellProfiler an in-house image analysis solution based on CellProfiler (Carpenter et al.,2006). 
  Features were grouped based on their source segmented objects (parts of the cell), which include: Cells (cell body), Nuclei (cell nuclei), and Perinuclei (perinuclear space).
  In addition, features derived from Voronoi segmentation of the images were included. The Qiagen siRNA library typically consists of four different siRNAs
  per gene, with the exception of talin1 and Cdc42, which had three and eight siRNAs, respectively.We used TargetScan’s (Lewis et al., 2005) predicted off-targets to define the 
  perturbation map for the eight genes across the 35 experiments.
  
  In order to convert single-cell data to gene-level binary data, we first
  applied B-score normalization to correct for row, column, and plate effects.
  Then we further normalised the data using MARS (multivariate adaptive
                                                  regression splines) and z-scoring. This entire process was performed using
  the R package singleCellFeatures (Bennet, 2015). We defined a
  control distribution because the biological controls were subject to strong
  edge effects from the experimental setup. The first and last two columns of
  wells constituted the control wells.We performedWilcoxon tests between
  all pairs of control populations and generated a distribution of p-values
  for each feature, choosing the lower 5th percentile as the critical p-value.
  Thiswas done to capture the differences across control populations. For the
  gene-level data, the knockdown populations were compared to six random
  control populations ( 10% of control wells) using a Wilcoxon test. The
  resulting p-values were combined using Fischer’s method, and this meta
  p-value was compared to the critical p-value. The feature was significant
  and took a value of 1 if the meta p-value was less than the critical p-value
  and 0 otherwise. The resulting gene-level binary data set consisted of 288
  features across 35 knockdown experiments.
}
\references{
  Bennet, N. (2015). Analysis of High Content RNA Interference Screens at Single Cell Level. Master’s thesis, ETH Zurich. \\
  Carpenter, A. E., Jones, T. R., Lamprecht, M. R., Clarke, C., Kang, I. H., Friman, O., Guertin, D. A., Chang, J. H., Lindquist, R. A., Moffat, J., et al. (2006). Cellprofiler: image analysis software for identifying and quantifying cell phenotypes. Genome biology, 7(10), R100.\\
  Lewis, B. P., Burge, C. B., and Bartel, D. P. (2005). Conserved seed pairing, often flanked by adenosines, indicates that thousands of human genes are microrna targets. cell, 120(1), 15–20.\\
  Smith, K., Li, Y., Piccinini, F., Csucs, G., Balazs, C., Bevilacqua, A., and Horvath, P. (2015). Cidre: an illumination-correction method for optical microscopy. Nature methods, 12(5), 404–406.\\
}
\seealso{\code{\link{BartonellaRNAi}}}
\examples{
  data("BartonellaRNAi")
  dim(D)          
}
\keyword{datasets}
