# R script
# Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>
# Thu Sep 14 17:03:56 EDT 2017

# based on  https://informatics.nescent.org/wiki/R_Hackathon_1/Ancestral_State_Reconstruction#Reconstructing_Ancestral_States_for_Discrete_Variables
# script runs two likelihood ratio tests, one with the ctenophore sister tree 
#   and one with the sponge sister tree.

library(ape)

mydir         <- "/Users/jfryan/Dropbox/00-JOE/06-Papers/2017/14-SASSON_ANCESTRAL_SEX/03-R_for2nd_revision/04-FINAL_R"

ctenosistree  <- "ctenosis.nex"
spongesistree <- "spongesis.nex"
charmat       <- "charmat_nomissindata.dat"

# we remove missing data because they cause problems with ace
unknowns <- c("Labidiaster_annulatus","Membranipora_membranacae","Aegina_citrea","Botrynema_brucei","Keratoisidinae","Coeloplana_bannwarthii","Trichoplax","Acanthella_acuta","Agelas_clathrodes","Aplysina_fulva","Calyx_podatypa","Chalinula_molitba","Dictyonella_incisa","Dragmacidon_lunaecharta","Eunapius_fragilis","Haliclona_sarai","Heterochone_calyx","Hexadella_pruvoti","Lissodendoryx_colombiensis","Lophocalyx_profundum","Placospongia_intermedia","Nodastrella_asconemaoida","Tedania_ignis","Tethya1","Tethya2","Monosiga_brevicollis","Proterospongia")

setwd(mydir)

mylrt <- function(mytree) {
    anitree <- read.nexus(mytree)
    anidata <- read.table(charmat)

    for (unk in unknowns) {
        anitree <- drop.tip(anitree, unk)
    }

    ERreconstruction <- ace(anidata$V2, anitree, type="discrete", model="ER")
    SYMreconstruction <- ace(anidata$V2, anitree, type="discrete", model="SYM")
    ARDreconstruction <- ace(anidata$V2, anitree, type="discrete", model="ARD")

    cat("RESULTS FOR ",mytree,"\n")

    erlnl <- ERreconstruction$loglik
    cat("ER log likelihood: ",erlnl,"\n")
    
    symlnl <- SYMreconstruction$loglik
    cat("SYM log likelihood: ",symlnl,"\n")
    
    ardlnl <- ARDreconstruction$loglik
    cat("ARD log likelihood: ",ardlnl,"\n")
    
    #  For a three-state character, 
    #     ER is a one parameter model, 
    #     SYM a three parameter model,
    #     ARD a six parameter model.

    # df = 5 (i.e. ARD(6) - ER(1))
    erard <- 1-pchisq(2*abs(ERreconstruction$loglik - ARDreconstruction$loglik), 5)
    cat("ER vs. ARD: ", erard, "\n")
    cat("    if <= 0.05 than ARD is significantly better than ER\n")

    # df = 3 (i.e. ARD(6) - SYM(3))
    symard <- 1-pchisq(2*abs(SYMreconstruction$loglik - ARDreconstruction$loglik), 3)
    cat("SYM vs. ARD: ",symard,"\n")
    cat("    if <= 0.05 than ARD is significantly better than SYM\n")

    # df = 2 (i.e. SYM(3) - ER(1))
    ersym <- 1-pchisq(2*abs(ERreconstruction$loglik - SYMreconstruction$loglik), 2)
    cat("ER vs. SYM: ",ersym,"\n")
    cat("    if <= 0.05 than SYM is significantly better than ER\n")
    cat("\n---------------------------------------------------------------\n\n")
}    
    
mylrt(ctenosistree)
mylrt(spongesistree)
sessionInfo()

