#' GBS data from Shirasawa et al (2017)
#'
#' Contains counts of reference alleles and total read counts from the GBS data of Shirasawa et al (2017) for
#' the three SNP's used as examples in Gerard, Ferrao, and Stephens (2017).
#'
#' @format A \code{\link[tibble]{tibble}} with 419 rows and 4 columns:
#' \describe{
#'     \item{id}{The identification label of the individuals.}
#'     \item{snp}{The SNP label.}
#'     \item{counts}{The number of read-counts that support the reference allele.}
#'     \item{size}{The total number of read-counts at a given SNP.}
#' }
#'
#' @source \url{http://sweetpotato-garden.kazusa.or.jp/}
#'
#' @references Shirasawa, Kenta, Masaru Tanaka, Yasuhiro Takahata, Daifu Ma, Qinghe Cao, Qingchang Liu, Hong Zhai et al. "A high-density SNP genetic map consisting of a complete set of homologous groups in autohexaploid sweetpotato (Ipomoea batatas)." Scientific Reports 7 (2017). DOI: 10.1038/srep44207
#'
#'   Gerard, David, Luis Felipe Ventorim Ferrão, and Matthew Stephens. 2017. "Harnessing Empirical Bayes and Mendelian Segregation for Genotyping Autopolyploids with Messy Sequencing Data." Overleaf Preprint.
#'
"snpdat"
