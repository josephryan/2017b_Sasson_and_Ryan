# R script
# Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>
# Thu Sep 14 17:04:24 EDT 2017

# based on http://blog.phytools.org/2013/03/estimating-ancestral-states-when-tips.html
# script runs two stochastic character mappings, one with the ctenophore
#     sister tree and one with the sponge sister tree.

set.seed(420)

library(maps)
library(ape)
library(phytools)

#mydir         <- "/Users/jfryan/Dropbox/00-JOE/06-Papers/2017/14-SASSON_ANCESTRAL_SEX/03-R_for2nd_revision/04-FINAL_R"
mydir         <- "."
mynsim        <- 1000
ctenosistree  <- "ctenosis.nex"
spongesistree <- "spongesis.nex"
spongeoutpdf   <- "simmap_spongesis.pdf"
ctenooutpdf   <- "simmap_ctenosis.pdf"
ctenoplotsoutpdf   <- "simmap_ctenosis_plots.pdf"
spongeplotsoutpdf  <- "simmap_spongesis_plots.pdf"
charmat_vals  <- "charmat_vals.dat"
charmat_names <- "charmat_names.dat"

setwd(mydir)

mysimmap <- function(mytree,mypdfout,mypdfplotout) {
    pdf(file=mypdfout,width=8.5,height=11)
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

    # we use anichrono to generate a pseudo-ultrametic tree
    anichrono <- chronopl(anitree, lambda = 0, age.min = 1)
    plot(anichrono,label.offset=.01, cex=0.4)
    nodelabels(pie=res_simmap$ace,piecol=c("yellow","red","blue"),cex=0.2)
    tiplabels(pie=res_simmap$tips,piecol=c("yellow","red","blue"),cex=0.2)

    pdf(file=mypdfplotout,width=8.5,height=11)
    plot(anitree,label.offset=.01, cex=0.4)
    par(mfrow=c(5,5))
    colors <- setNames(c(c="yellow","red","blue"),c("asexual","hermaphrodite","separatesexes")) 
    plot(SYM.simmap_trees,lwd=1,ftype="off",colors=colors)

    print(res_simmap$ace)
    cat("\n------------------------------------------------------\n")
}

mysimmap(ctenosistree,ctenooutpdf,ctenoplotsoutpdf)
mysimmap(spongesistree,spongeoutpdf,spongeplotsoutpdf)
sessionInfo()

