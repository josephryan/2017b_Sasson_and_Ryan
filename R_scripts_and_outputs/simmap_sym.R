# R script
# Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>
# Sep 2017

# based on http://blog.phytools.org/2013/03/estimating-ancestral-states-when-tips.html
# and http://blog.phytools.org/2017/09/densitymap-like-plot-for-3-state.html
# script runs two stochastic character mappings, one with the ctenophore
#     sister tree and one with the sponge sister tree.

# set seed so script gives same answer each run
set.seed(420)

# if these packages are not installed 
# install.packages("maps")
# install.packages("ape")
# install.packages("phytools")

library(maps)
library(ape)
library(phytools)

# set variables; these can be adjusted if running on different data
mydir         <- "."
mynsim        <- 1000
ctenosistree  <- "ctenosis.nex"
spongesistree <- "spongesis.nex"
spongeoutpdf   <- "simmap_spongesis.pdf"
ctenooutpdf   <- "simmap_ctenosis.pdf"
ctenodensityoutpdf   <- "simmap_ctenosis_densityTree.pdf"
spongedensityoutpdf  <- "simmap_spongesis_densityTree.pdf"
charmat_vals  <- "charmat_vals.dat"
charmat_names <- "charmat_names.dat"

setwd(mydir)

# function runs make.simmap and densityTree
mysimmap <- function(mytree,mypdfout,mypdfdensityout) {

    # set output pdf file
    pdf(file=mypdfout,width=8.5,height=11)

    # read in tree and matrix; matrix is broken up into data and rownames
    # outgroups (Monosiga_brevicollis and Proterospongia) have been removed
    anitree <- read.nexus(mytree)
    anidata <- read.table(charmat_vals,sep=",")
    rn <- read.table(charmat_names)
    animatrix <- as.matrix(anidata)
    rownames(animatrix) <- rn$V1
    colnames(animatrix) <- c("separatesexes","hermaphrodite", "asexual")

    # outgroups are dropped because difficult to compare sexuality in 
    #    single-cellular organism vs. multi-cellular organism.
    #    If they are coded as unknown make.simmap will assign them a 
    #    sexual mode (separatesexes, hermaphrodite, or asexual).
    anitree <- drop.tip(anitree, "Monosiga_brevicollis")
    anitree <- drop.tip(anitree, "Proterospongia")

    # run make.simmap and print out information about the run
    SYM.simmap_trees <- make.simmap(anitree,animatrix,nsim=mynsim,model="SYM")
    SYM.simmap_trees$loglike
    res_simmap <- describe.simmap(SYM.simmap_trees)
    print(res_simmap)
    print (res_simmap$ace)

    # make psuedo chrono tree for plotting
    anichrono <- chronopl(anitree, lambda = 0, age.min = 1)

    # plot tree with labels slightly offset so pies can go at tips
    plot(anichrono,label.offset=.01, cex=0.4)
    nodelabels(pie=res_simmap$ace,piecol=c("yellow","red","blue"),cex=0.2)
    tiplabels(pie=res_simmap$tips,piecol=c("yellow","red","blue"),cex=0.2)

    # reset pdf file for densityTree (getting blank 1st page so onfile=FALSE)
    dev.off()
    pdf(file=mypdfdensityout,onefile=FALSE,width=8.5,height=11)

    # set colors; adjust margins and fontsize; 
    colors <- setNames(c(c="yellow","red","blue"),
                       c("asexual","hermaphrodite","separatesexes")) 
    par(mar=c(1,0,0,1),cex=0.6)

    # run densityTree
    densityTree(SYM.simmap_trees,method="plotSimmap",lwd=3,colors=colors);

    cat("\n------------------------------------------------------\n")
}

# call our simmap function twice with different tree topology
mysimmap(ctenosistree,ctenooutpdf,ctenodensityoutpdf)
mysimmap(spongesistree,spongeoutpdf,spongedensityoutpdf)

# print versions of all loaded modules so analysis can be reproduced
sessionInfo()

