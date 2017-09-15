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

mydir         <- "/Users/jfryan/Dropbox/00-JOE/06-Papers/2017/14-SASSON_ANCESTRAL_SEX/03-R_for2nd_revision/04-FINAL_R"
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

    # we use anichrono to generate a pseudo-ultrametic tree
    anichrono <- chronopl(anitree, lambda = 0, age.min = 1)
    plot(anichrono,label.offset=.01, cex=0.8)
    nodelabels(pie=res_simmap$ace,piecol=c("yellow","red","blue"),cex=0.2)
    tiplabels(pie=res_simmap$tips,piecol=c("yellow","red","blue"),cex=0.2)

    # FOR THE NEXT COMMAND THESE COMMANDS AND SUBSEQUENT VALUES WILL BE HELPFUL
    # use nodelabels(bg="white") to label nodes with corresponding numbers
    # key nodes:
    #                                           asexual       sep           herm
    #   animal LCA = 166                # 166   
    #   ctenophore LCA = 167            # 167   
    #   sponge+parahoxozoa LCA = 177    # 177   
    #   sponge LCA = 178                # 178   
    #   parahoxozoa LCA = 239           # 239   
    #   cnidarian-bilaterian LCA = 240  # 240   
    #   bilateria LCA = 271             # 271   
    #   cnidaria LCA = 241              # 241   

    # The following prints all node values
    # NOTE: HEADINGS ARE WRONG. SEE CTENOPHORE LCA 167 in ctenosis SHOULD BE: 
    #     ASEXUAL SEPARATESEXES HERMAPHRODITE
    #     asexual hermaphrodite separatesexes
    #166   0.000         0.542         0.458

    print(res_simmap$ace)
    cat("\n------------------------------------------------------\n")
}

mysimmap(ctenosistree,ctenooutpdf)
mysimmap(spongesistree,spongeoutpdf)
sessionInfo()

