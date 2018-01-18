\name{SahinRNAi2008}
\alias{SahinRNAi2008}
\alias{dat.normalized}
\alias{dat.unnormalized}
\alias{map.int2node}
\docType{data}
\title{Combinatorial Protein Knockdowns in the ERBB Signaling Pathway}
\description{
Sixteen RNAi knockdowns (including 3 double knockdowns) of proteins in the ERBB signaling pathway of trastuzumab resistant breast cancer cells were conducted. Reverse Phase Protein Array (RPPA) measurments for 10 signaling intermediates are available before and after EGF stimulation with 4 technical and 3 biological replicates.
}
\usage{data(SahinRNAi2008)}
\format{
  dat.unnormalized: 408 x 17 matrix (rows = RPPA measurements for (16 KOs + MOCK) x 4 technical x 3 biological replicates, columns = 10 antibodies + 6 proteins without measurements + time)
  dat.normalized: 408 x 17 matrix (measurements from technical and biological replicates are quantile normalized for each RNAi experiment)
  map.int2node: list with names being interventions (=names of dat.normalized) and entries being node names (=column names of dat.normalized)
}
\details{
The cells were lysed on ice by scraping the cells in M-PER lysis buffer
(Pierce, Rockford, IL) containing protease inhibitor Complete Mini (Roche,
Basel), anti-phosphatase PhosSTOP (Roche, Basel), 10 mM NaF and 1mM
Na4VO3. Protein concentrations were determined with a BCA Protein Assay
Reagent Kit (Pierce, Rockford, IL). Lysates were mixed 1:2 with 2 times Protein Arraying Buffer (Whatman, Brentfort, UK) to obtain a final protein
concentration of 1.5 mug/muL. Briefly, these lysates were printed onto nitrocellulose
coated ONCYTE-slides (Grace Bio Labs, Bend, USA) using a
non-contact piezo spotter, sciFlexxarrayer S5 (Scienion, Berlin, Germany).
After primary and near-infrared (NIR)-dye labeled secondary antibodies applied,
spots were analysed using an Odyssey scanner (LI-COR, Lincoln,
USA) and signal intensities were quantified using Odyssey 2.0 software (For
detailed information and an antibody list, see Sahin et al., 2008). Since no
antibody against MEK1 was available, we measured protein expression of
pERK1/2, which is downstream of MEK1.
}
\references{
  Oezguer Sahin, Holger Froehlich, Christian Loebke, Ulrike Korf, Sara Burmester, Meher Majety, Jens Mattern, Ingo Schupp, Claudine Chaouiya, Denis Thieffry, Annemarie Poustka, Stefan Wiemann, Tim Beissbarth, Dorit Arlt, Modeling ERBB receptor-regulated G1/S transition to find novel targets for de novo trastuzumab resistance, BMC Systems Biology, 2008
}
\seealso{\code{\link{BoutrosRNAi2002}}}
\examples{
    data("SahinRNAi2008")
    dim(dat.normalized)
    dim(dat.unnormalized)      
}
\keyword{datasets}
