# R script
# Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>
# Thu Sep 14 17:04:24 EDT 2017

# script to label nodes for simmap_analysis

set.seed(420)

library(maps)
library(ape)
library(phytools)

mydir         <- "/Users/jfryan/Dropbox/00-JOE/06-Papers/2017/14-SASSON_ANCESTRAL_SEX/03-R_for2nd_revision/04-FINAL_R/02-PRINT_NODENUMBERS"
mynsim        <- 1
ctenosistree  <- "../ctenosis.nex"
spongesistree <- "../spongesis.nex"
ctenooutpdf   <- "ctenosis_nodenumbers.pdf"
spongeoutpdf  <- "spongesis_nodenumbers.pdf"
charmat_vals  <- "../charmat_vals.dat"
charmat_names <- "../charmat_names.dat"

setwd(mydir)

mysimmap <- function(mytree,mypdfout) {
    pdf(file=mypdfout,width=25.5,height=33)
    anitree <- read.nexus(mytree)

    anitree <- drop.tip(anitree, "Monosiga_brevicollis")
    anitree <- drop.tip(anitree, "Proterospongia")

    # we use anichrono to generate a pseudo-ultrametic tree
    anichrono <- chronopl(anitree, lambda = 0, age.min = 1)
    plot(anichrono,cex=0.8)
    nodelabels(bg="white")

}

mysimmap(ctenosistree,ctenooutpdf)
mysimmap(spongesistree,spongeoutpdf)
sessionInfo()

