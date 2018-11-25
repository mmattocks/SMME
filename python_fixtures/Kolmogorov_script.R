#load acss libraries after, if necessary, installing
require(acss)
library(acss)

#Read in sequence files
empiricalData <- read.table("~/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/EOsequences", header=TRUE, colClasses=c(rep("numeric",1),"character"))
heData <- read.table("~/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/He", header=TRUE, colClasses=c(rep("numeric",2),"character"))
heRefitData <- read.table("~/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/HeRefit", header=TRUE, colClasses=c(rep("numeric",2),"character"))
deterministicData <- read.table("~/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/Deterministic", header=TRUE, colClasses=c(rep("numeric",2),"character"))
gomesData <- read.table("~/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/Gomes", header=TRUE, colClasses=c(rep("numeric",2),"character"))
boijeData <- read.table("~/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/Boije", header=TRUE, colClasses=c(rep("numeric",2),"character"))

empiricalACSS <- likelihood_d(empiricalData[,'Sequence'], alphabet=4)
heACSS <- likelihood_d(heData[,'Sequence'], alphabet=4)
heRefitACSS <- likelihood_d(heRefitData[,'Sequence'], alphabet=4)
deterministicACSS <- likelihood_d(deterministicData[,'Sequence'], alphabet=4)
gomesACSS <- likelihood_d(gomesData[,'Sequence'], alphabet=4)
boijeACSS <- likelihood_d(boijeData[,'Sequence'], alphabet=4)

write.table(empiricalACSS, "~/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/EO_acss", sep="\t", row.names=TRUE, col.names=TRUE)
write.table(heACSS, "~/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/he_acss", sep="\t", row.names=TRUE, col.names=TRUE)
write.table(heRefitACSS, "~/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/heRefit_acss", sep="\t", row.names=TRUE, col.names=TRUE)
write.table(deterministicACSS, "~/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/deterministic_acss", sep="\t", row.names=TRUE, col.names=TRUE)
write.table(gomesACSS, "~/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/gomes_acss", sep="\t", row.names=TRUE, col.names=TRUE)
write.table(boijeACSS, "~/git/chaste/projects/ISP/python_fixtures/testoutput/KolmogorovSequences/boije_acss", sep="\t", row.names=TRUE, col.names=TRUE)
