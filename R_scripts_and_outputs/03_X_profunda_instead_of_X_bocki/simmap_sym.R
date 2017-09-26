# R script
# Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>
# Thu Sep 14 17:04:24 EDT 2017

# this version is very similar to ../simmap_sym.R but uses data files in this
# directory which have swapped out Xenoturbella bocki (hermaphrodite) for Xenoturbella profunda (separate sexes) to see the effect of this taxa

set.seed(420)

library(maps)
library(ape)
library(phytools)

mydir         <- .
mynsim        <- 1000
ctenosistree  <- "ctenosis.nex"
spongesistree <- "spongesis.nex"
ctenooutpdf   <- "simmap_ctenosis.pdf"
spongeoutpdf  <- "simmap_spongesis.pdf"
charmat_vals  <- "charmat_vals.dat"
charmat_names <- "charmat_names.dat"

setwd(mydir)

mysimmap <- function(mytree,mypdfout) {
    pdf(file=mypdfout,width=25.5,height=33)
    anitree <- read.nexus(mytree)
    anidata <- read.table(charmat_vals,sep=",")
    rn <- read.table(charmat_names)

    anitree <- drop.tip(anitree, "Monosiga_brevicollis")
    anitree <- drop.tip(anitree, "Proterospongia")

    animatrix <- as.matrix(anidata)
    rownames(animatrix) <- rn$V1
    colnames(animatrix) <- c("separatesexes","hermaphrodite", "asexual")

    SYM.simmap_trees <- make.simmap(anitree,animatrix,nsim=mynsim,model="SYM")
    SYM.simmap_trees$loglike

    res_simmap <- describe.simmap(SYM.simmap_trees)

    print(res_simmap)

    anichrono <- chronopl(anitree, lambda = 0, age.min = 1)
    plot(anichrono,label.offset=.01, cex=0.8)
    nodelabels(pie=res_simmap$ace,piecol=c("yellow","red","blue"),cex=0.2)
    tiplabels(pie=res_simmap$tips,piecol=c("yellow","red","blue"),cex=0.2)

    print(res_simmap$ace)
    cat("\n------------------------------------------------------\n")
}

mysimmap(ctenosistree,ctenooutpdf)
mysimmap(spongesistree,spongeoutpdf)
sessionInfo()

