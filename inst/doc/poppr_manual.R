### R code from vignette source 'poppr_manual.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: install (eval = FALSE)
###################################################
## install.packages("poppr", dependencies=TRUE)


###################################################
### code chunk number 2: install_depend (eval = FALSE)
###################################################
## install.packages(c("adegenet", "pegas", "vegan", "ggplot2", "phangorn", "ape", "igraph"))


###################################################
### code chunk number 3: install_source (eval = FALSE)
###################################################
## install.packages("/path/to/poppr.tar.gz", type="source", repos=NULL)


###################################################
### code chunk number 4: install_devtools (eval = FALSE)
###################################################
## install.packages("devtools")


###################################################
### code chunk number 5: install_github (eval = FALSE)
###################################################
## library(devtools)
## install_github(repo = "grunwaldlab/poppr")


###################################################
### code chunk number 6: install_devel (eval = FALSE)
###################################################
## library(devtools)
## install_github(repo = "grunwaldlab/poppr", ref = "devel")


###################################################
### code chunk number 7: poppr_manual.Rnw:132-134
###################################################
library(poppr)
x <- list(files="/path/to/R/poppr/files/rootrot.csv", path="/path/to/R/poppr/files")


###################################################
### code chunk number 8: getfilefunk (eval = FALSE)
###################################################
## library(poppr)


###################################################
### code chunk number 9: getfilefunk2 (eval = FALSE)
###################################################
## x <- getfile()


###################################################
### code chunk number 10: getfilex
###################################################
x


###################################################
### code chunk number 11: firstpoppr (eval = FALSE)
###################################################
## popdata <- poppr(x$files)


###################################################
### code chunk number 12: aflp
###################################################
options(width=90)
popprsoutput <- poppr(system.file("files/rootrot.csv", package="poppr"))


###################################################
### code chunk number 13: firstpoppr2 (eval = FALSE)
###################################################
## popdata


###################################################
### code chunk number 14: aflp2
###################################################
popprsoutput


###################################################
### code chunk number 15: multi_getfile (eval = FALSE)
###################################################
## x <- getfile(multi=TRUE)


###################################################
### code chunk number 16: show_multi_getfile1 (eval = FALSE)
###################################################
## x


###################################################
### code chunk number 17: show_multi_getfile
###################################################
x$files <- list.files(dirname(system.file("files/rootrot.csv", package="poppr")))
x$files <- paste(x$path, x$files, sep="/")
x


###################################################
### code chunk number 18: poppr.all (eval = FALSE)
###################################################
## poppr.all(x$files)


###################################################
### code chunk number 19: poppr.all_eval
###################################################
poppr.all(c(system.file("files/rootrot.csv", package="poppr"), system.file("files/rootrot2.csv", package="poppr"), system.file("files/simulated.dat", package="poppr")))


###################################################
### code chunk number 20: getfile_adegenet (eval = FALSE)
###################################################
## getfile(multi=TRUE)


###################################################
### code chunk number 21: list_adegenet_files
###################################################
nancylist <- dir(dirname(system.file("files/nancycats.gtx",package="adegenet")))
list(
  files = paste("/path/to/R/adegenet/files", nancylist, sep="/"),
  path = "/path/to/R/adegenet/files"
)


###################################################
### code chunk number 22: getfile_nancy (eval = FALSE)
###################################################
## getfile(multi=TRUE, pattern="nancy")


###################################################
### code chunk number 23: show_getfile_nancy
###################################################
nancylist <- list.files(dirname(system.file("files/nancycats.gtx",package="adegenet")), pattern="nancy")
list(
  files = paste("/path/to/R/adegenet/files", nancylist, sep="/"),
  path = "/path/to/R/adegenet/files"
)


###################################################
### code chunk number 24: getfile_gtx (eval = FALSE)
###################################################
## getfile(multi=TRUE, pattern="gtx")


###################################################
### code chunk number 25: show_getfile_gtx
###################################################
nancylist <- list.files(dirname(system.file("files/nancycats.gtx",package="adegenet")), pattern="gtx")
list(
  files = paste("/path/to/R/adegenet/files", nancylist, sep="/"),
  path = "/path/to/R/adegenet/files"
)


###################################################
### code chunk number 26: getfile_dat (eval = FALSE)
###################################################
## getfile(multi=TRUE, pattern="dat")


###################################################
### code chunk number 27: show_getfile_dat
###################################################
nancylist <- list.files(dirname(system.file("files/nancycats.gtx",package="adegenet")), pattern="dat")
list(
  files = paste("/path/to/R/adegenet/files", nancylist, sep="/"),
  path = "/path/to/R/adegenet/files"
)


###################################################
### code chunk number 28: getfile_datend (eval = FALSE)
###################################################
## getfile(multi=TRUE, pattern="dat$")


###################################################
### code chunk number 29: show_getfile_datend
###################################################
nancylist <- list.files(dirname(system.file("files/nancycats.gtx",package="adegenet")), pattern="dat$")
list(
  files = paste("/path/to/R/adegenet/files", nancylist, sep="/"),
  path = "/path/to/R/adegenet/files"
)


###################################################
### code chunk number 30: system_file_genalex (eval = FALSE)
###################################################
## system.file("files/rootrot.csv", package="poppr")


###################################################
### code chunk number 31: system_file_echo
###################################################
paste("/path/to/R/library/poppr/files/rootrot.csv")


###################################################
### code chunk number 32: read.genalex_ex
###################################################
rootrot <- read.genalex(system.file("files/rootrot.csv", package="poppr"))


###################################################
### code chunk number 33: read.genalex_ex2
###################################################
rootrot


###################################################
### code chunk number 34: microsave (eval = FALSE)
###################################################
## library(poppr)
## data(microbov)
## microbov@other$population_hierarchy <- data.frame(list(Country = microbov@other$coun, 
##   Species = microbov@other$spe, Breed = microbov@other$breed))
## microbov <- splitcombine(microbov, method=2, hier=c("Country", "Species", "Breed"))
## genind2genalex(microbov, file="~/Desktop/microbov.csv")


###################################################
### code chunk number 35: microsave_message
###################################################
cat("Extracting the table ... Writing the table to ~/Desktop/microbov.csv ... Done.")


###################################################
### code chunk number 36: genind2genalex (eval = FALSE)
###################################################
## genind2genalex(rootrot, "~/Desktop/rootrot.csv")


###################################################
### code chunk number 37: genind2genalex_cat
###################################################
cat("Extracting the table ... Writing the table to ~/Desktop/rootrot.csv ... Done.\n")


###################################################
### code chunk number 38: nancyxy
###################################################
data(nancycats)
nancycats@other$xy


###################################################
### code chunk number 39: genind2genalex_nancy (eval = FALSE)
###################################################
## genind2genalex(nancycats, "~/Desktop/nancycats_pop_xy.csv")


###################################################
### code chunk number 40: genind2genalex_cat2
###################################################
cat("Extracting the table ... Writing the table to ~/Desktop/nancycats_pop_xy.csv ... Done.\n")


###################################################
### code chunk number 41: nancy_grow_xy
###################################################
nan2 <- nancycats
nan2@other$xy <- nan2@other$xy[rep(1:17, table(pop(nan2))), ]
head(nan2@other$xy)


###################################################
### code chunk number 42: genind2genalex_nancy_grow (eval = FALSE)
###################################################
## genind2genalex(nan2, "~/Desktop/nancycats_inds_xy.csv")


###################################################
### code chunk number 43: genind2genalex_cat3
###################################################
cat("Extracting the table ... Writing the table to ~/Desktop/nancycats_inds_xy.csv ... Done.\n")


###################################################
### code chunk number 44: ex_data
###################################################
library(adegenet)
df <- data.frame(list(locus1=c("101/101", "102/103", "102/102"), 
                      locus2=c("201/201","202/203","203/204"), 
                      locus3=c("301/302", "301/303", "304/305")))
df


###################################################
### code chunk number 45: ex_genind
###################################################
df2genind(df, sep="/")$tab


###################################################
### code chunk number 46: example_data_frame
###################################################
df <- data.frame(list(locus1=c("101/101", "102/103", "102/102"), 
                      locus2=c("201/201", "202/203", "203/204"), 
                      locus3=c("301/302", "301/303", "304/305")
                      )
                 )
dfg <- df2genind(df, sep="/")


###################################################
### code chunk number 47: display_example
###################################################
dfg


###################################################
### code chunk number 48: other_slot
###################################################
# First off, how big is the object?
print(object.size(dfg), units="auto")
dfg$other$dfg <- dfg
dfg # we can now see that the @other slot is now filled.
dfg$other$dfg
print(object.size(dfg), units="auto") # How big is it now?


###################################################
### code chunk number 49: H3N2_data1
###################################################
data(H3N2)
H3N2
pop(H3N2)
H3N2$pop.names


###################################################
### code chunk number 50: H3N2_other_slot
###################################################
head(H3N2$other$x)
nrow(H3N2$other$x)


###################################################
### code chunk number 51: setting_H3N2_population
###################################################
pop(H3N2) <- H3N2$other$x$country
head(pop(H3N2))
H3N2$pop.names


###################################################
### code chunk number 52: initializing_poppr
###################################################
library(poppr)


###################################################
### code chunk number 53: initializing_nancycats
###################################################
data(nancycats)


###################################################
### code chunk number 54: nancy_summary
###################################################
summary(nancycats)


###################################################
### code chunk number 55: nancy_indiv
###################################################
nancycats$loc.names # Names of the loci
nancycats$tab[1:5, 8:13]


###################################################
### code chunk number 56: missingno_replace
###################################################
nanzero <- missingno(nancycats, type = "zero")
nanmean <- missingno(nancycats, type = "mean")
nanzero$tab[1:5, 8:13]
nanmean$tab[1:5, 8:13]


###################################################
### code chunk number 57: missingno_exclude
###################################################
nanloci <-  missingno(nancycats, "loci")
nangeno <-  missingno(nancycats, "geno")
nanloci$tab[1:5, 8:13]


###################################################
### code chunk number 58: missingno_loci
###################################################
length(nanloci$ind.names) # Individuals
nanloci$loc.names # Names of the loci


###################################################
### code chunk number 59: missingno_geno
###################################################
nangeno$tab[1:5, 8:13]
length(nangeno$ind.names) # Individuals
nangeno$loc.names # Names of the loci


###################################################
### code chunk number 60: Aeut
###################################################
Aeut <- read.genalex(system.file("files/rootrot.csv", package="poppr"))
summary(Aeut)


###################################################
### code chunk number 61: splitcombine_view
###################################################
head(Aeut$other$population_hierarchy)


###################################################
### code chunk number 62: splitcombine_split
###################################################
Aeut.pop <- splitcombine(Aeut, method=1, dfname="population_hierarchy", hier=c("Pop", "Subpop"), setpopulation=TRUE)
head(Aeut.pop$other$population_hierarchy)
summary(Aeut.pop)


###################################################
### code chunk number 63: splitcombine_poppr
###################################################
poppr(Aeut.pop, quiet=TRUE)


###################################################
### code chunk number 64: read_rootrot2
###################################################
Aeut2 <- read.genalex(system.file("files/rootrot2.csv", package="poppr"), region=TRUE)
head(Aeut2@other$population_hierarchy)
summary(Aeut2)


###################################################
### code chunk number 65: splitcombine_combine
###################################################
Aeut2.combine <- splitcombine(Aeut2, method=2, hier=2:1)
head(Aeut2.combine@other$population_hierarchy)
summary(Aeut2.combine)


###################################################
### code chunk number 66: popsub_sublist
###################################################
data(H3N2)
pop(H3N2) <- H3N2$other$x$country
H3N2$pop.names # Only two countries from North America.
H.na <- popsub(H3N2, sublist=c("USA", "Canada"))
H.na$pop.names


###################################################
### code chunk number 67: popsub_sizes
###################################################
nInd(H.na)
nInd(H3N2)


###################################################
### code chunk number 68: popsub_blacklist
###################################################
H.minus.na <- popsub(H3N2, blacklist=c("USA", "Canada"))
H.minus.na$pop.names


###################################################
### code chunk number 69: length_test
###################################################
(nInd(H.minus.na) + nInd(H.na)) == nInd(H3N2)


###################################################
### code chunk number 70: popsub_combine
###################################################
Hsort <- sort(H3N2$pop.names)[1:10]
Hsort
H.alph <- popsub(H3N2, sublist=Hsort, blacklist=c("USA", "Canada"))
H.alph$pop.names


###################################################
### code chunk number 71: clonecorrect
###################################################
data(Aeut)
A.cc <- clonecorrect(Aeut,  hier=c("Pop", "Subpop"), keep=1)
poppr(A.cc, quiet=TRUE)


###################################################
### code chunk number 72: clonecorrect_comparison
###################################################
poppr(Aeut, quiet=TRUE)


###################################################
### code chunk number 73: shuffle_mat
###################################################
exmat <- matrix(c(4,4,
         4,1,
         4,3,
         2,2,
         3,3), 5, byrow=TRUE)
exmat


###################################################
### code chunk number 74: multilocus_shuffle
###################################################
set.seed(1001)
exmat[sample(1:5), ]


###################################################
### code chunk number 75: permutation
###################################################
set.seed(1001)
matrix(sample(exmat), 5, byrow=T)


###################################################
### code chunk number 76: param_boot
###################################################
set.seed(1001)
cat("First Sample")
matrix(sample(1:4, 10, prob=c(0.1,0.2,0.3,0.4), replace=TRUE), 5, byrow=T)
cat("Second Sample")
matrix(sample(1:4, 10, prob=c(0.1,0.2,0.3,0.4), replace=TRUE), 5, byrow=T)


###################################################
### code chunk number 77: boot
###################################################
set.seed(1001)
matrix(sample(1:4, 10, prob=rep(1, 4), replace=TRUE), 5, byrow=T)


###################################################
### code chunk number 78: shuffle_ia
###################################################
data(nancycats)
nan1 <- popsub(nancycats, 1)
ia(nan1)
replicate(10, ia(shufflepop(nan1, method = 3), quiet=TRUE))


###################################################
### code chunk number 79: inform.H3N2.1 (eval = FALSE)
###################################################
## data(H3N2)
## H.five <- informloci(H3N2, cutoff = 0.05)


###################################################
### code chunk number 80: inform.H3N2.2
###################################################
res <- c(157,177,233,243,262,267,280,303,313,327,357,382,384,399,412,418,424,425,429,433,451,470,529,546,555,557,564,576,592,595,597,602,612,627,642,647,648,654,658,663,667,681,717,806,824,837,882)
cat("cutoff value: 5 percent ( 95 individuals ).\n","47 uninfomative loci found:", res, fill = 80)


###################################################
### code chunk number 81: inform.nancy
###################################################
data(nancycats)
naninform <- informloci(nancycats, cutoff = 0.05)


###################################################
### code chunk number 82: mlg
###################################################
data(H3N2)
mlg(H3N2, quiet=FALSE)


###################################################
### code chunk number 83: crosspop
###################################################
pop(H3N2) <- H3N2$other$x$country
H.dup <- mlg.crosspop(H3N2, quiet=TRUE)


###################################################
### code chunk number 84: crosspopout
###################################################
H.inds <- mlg.crosspop(H3N2, indexreturn=TRUE)
Hadoo <- mlg.crosspop(H3N2, mlgsub=H.inds[1:10])


###################################################
### code chunk number 85: crosspop2
###################################################
head(H.dup)
H.num <- sapply(H.dup, length) # count the number of populations each MLG crosses.
H.num


###################################################
### code chunk number 86: mlgbar (eval = FALSE)
###################################################
## H.tab <- mlg.table(H3N2, quiet=TRUE, bar=TRUE)
## H.tab[1:10, 1:10] # Showing the first 10 columns and rows of the table.


###################################################
### code chunk number 87: poppr_manual.Rnw:1013-1015
###################################################
H.tab <- mlg.table(H3N2, quiet=TRUE, bar=FALSE)
H.tab[1:10, 1:10]


###################################################
### code chunk number 88: mlgbarplot
###################################################
mlg.table(H3N2, sublist="Norway", quiet=TRUE, bar=TRUE)


###################################################
### code chunk number 89: mlgrare1
###################################################
H.year <- H3N2
pop(H.year) <- H.year$other$x$year
summary(H.year) # Check the data to make sure it's correct.


###################################################
### code chunk number 90: mlgrare2 (eval = FALSE)
###################################################
## library(vegan)
## H.year <- mlg.table(H.year, bar=FALSE)
## rarecurve(H.year, ylab="Multilocus genotypes expected", sample=min(rowSums(H.year)))


###################################################
### code chunk number 91: mlgrareplot
###################################################
library(vegan)
H.year <- mlg.table(H.year, bar=FALSE)
rarecurve(H.year, ylab="Multilocus Genotypes Expected", sample=min(rowSums(H.year)))


###################################################
### code chunk number 92: subcross
###################################################
UGNN.list <- c("United Kingdom", "Germany", "Netherlands", "Norway")
UGNN <- mlg.crosspop(H3N2, sublist=UGNN.list, indexreturn=TRUE)


###################################################
### code chunk number 93: subtable
###################################################
UGNN # Note that we have three numbers here. This will index the columns for us.
UGNN.list # And let's not forget that we have the population names.
H.tab[UGNN.list, UGNN]


###################################################
### code chunk number 94: subdata
###################################################
H.vec <- mlg.vector(H3N2)
H.sub <- H3N2[H.vec %in% UGNN, ]
mlg.table(H.sub, bar=FALSE)


###################################################
### code chunk number 95: mlgsub_flag (eval = FALSE)
###################################################
## mlg.table(H3N2, mlgsub=UGNN, bar=TRUE)


###################################################
### code chunk number 96: mlgsub_flagshow
###################################################
mlg.table(H3N2, mlgsub=UGNN, bar=FALSE)


###################################################
### code chunk number 97: subnor
###################################################
mlg.table(H3N2, sublist="Norway", mlgsub=UGNN)


###################################################
### code chunk number 98: poppr_manual.Rnw:1107-1109
###################################################
H.vec[1:22]
mlg.vector(H.sub)


###################################################
### code chunk number 99: poppr_manual.Rnw:1113-1116
###################################################
H3N2@other$MLG.vector <- H.vec
H.sub <- H3N2[H.vec %in% UGNN, ]
H.sub@other$MLG.vector


###################################################
### code chunk number 100: df hashing_echo (eval = FALSE)
###################################################
## H.df <- genind2df(H3N2)
## H.df[H.vec %in% UGNN, 1:15] # Showing only 15 columns becaus it is a large dataset.


###################################################
### code chunk number 101: df hashing_eval
###################################################
H.df <- genind2df(H3N2[, loc=names(H3N2@loc.names)[1:15]])
H.df[H.vec %in% UGNN, 1:15] # Showing only 15 columns becaus it is a large dataset.


###################################################
### code chunk number 102: subsetting
###################################################
UGNN
H.vec[H.vec %in% UGNN]


###################################################
### code chunk number 103: mlg_chart
###################################################
df <- mlg.crosspop(H3N2, df=TRUE, quiet=TRUE)
names(df)


###################################################
### code chunk number 104: poppr_manual.Rnw:1147-1150
###################################################
H.max <- names(sort(H.num, decreasing=TRUE)[1:10])
# Showing the data frame by the largest MLG complex.
df[df$MLG %in% H.max[1], ]


###################################################
### code chunk number 105: ggplotchart
###################################################
df2 <- df[df$MLG %in% H.max, ]
library(ggplot2)
qplot(y=MLG, x=Population, data=df2, color=Count, size=Count) +
  theme(axis.text.x = element_text(size = 10, angle = -45, hjust = 0))


###################################################
### code chunk number 106: ia_demo
###################################################
ia(nancycats)


###################################################
### code chunk number 107: ia_demo_sampling_dummy (eval = FALSE)
###################################################
## set.seed(1001)
## ia(popsub(nancycats, 5), sample=999)


###################################################
### code chunk number 108: ia_demo_sampling
###################################################
set.seed(1001)
simplenan <- ia(popsub(nancycats, 5))
invisible(lapply(1:999, function(x) {cat(ifelse(x %% 50 == 0, ".\n", "."))}))
c(simplenan[1], p.Ia = 0.589, simplenan[2], p.rD = 0.589)


###################################################
### code chunk number 109: poppr_manual.Rnw:1218-1220
###################################################
set.seed(1001)
nan5 <- ia(popsub(nancycats, 5), sample=999, quiet=TRUE)


###################################################
### code chunk number 110: ch_pval_one (eval = FALSE)
###################################################
## set.seed(1001)
## ia(popsub(nancycats, 5), sample=999, method=3, quiet=TRUE, hist=FALSE)


###################################################
### code chunk number 111: ch_pval_two
###################################################
set.seed(1001)
simplenan2 <- ia(popsub(nancycats, 5))
c(simplenan2[1], p.Ia = 0.589, simplenan2[2], p.rD = 0.596)


###################################################
### code chunk number 112: ia_Aeut_ex (eval = FALSE)
###################################################
## data(Aeut)
## set.seed(1001)
## ia(popsub(Aeut, 1), sample=999, method=3, quiet=TRUE, hist=FALSE)


###################################################
### code chunk number 113: ia_Aeut_ex_real
###################################################
data(Aeut)
set.seed(1001)
A.dum <- ia(popsub(Aeut, 1))
c(A.dum[1], p.Ia = 0.001, A.dum[2], p.rD = 0.001)


###################################################
### code chunk number 114: diss_dist (eval = FALSE)
###################################################
## data(Aeut)
## A.dist <- diss.dist(Aeut)
## heatmap(as.matrix(A.dist), symm=TRUE)


###################################################
### code chunk number 115: poppr_manual.Rnw:1276-1279
###################################################
data(Aeut)
A.dist <- diss.dist(Aeut)
heatmap(as.matrix(A.dist), symm=TRUE)


###################################################
### code chunk number 116: bruvo_matrix
###################################################
dist9 <- bruvo.dist(popsub(nancycats, 9), replen=rep(1,9))
dist9


###################################################
### code chunk number 117: bruvo_heat
###################################################
heatmap(as.matrix(dist9), symm=TRUE)


###################################################
### code chunk number 118: poppr_manual.Rnw:1313-1314
###################################################
heatmap(as.matrix(dist9), symm=TRUE)


###################################################
### code chunk number 119: attrLabels
###################################################
attr(dist9, "Labels")


###################################################
### code chunk number 120: attrLabel_replace
###################################################
dist9.attr <- attr(dist9, "Labels")
attr(dist9, "Labels") <- paste(rep("P09", 9), dist9.attr)
dist9


###################################################
### code chunk number 121: popcompare_bruvo1 (eval = FALSE)
###################################################
## dist9to8 <- bruvo.dist(popsub(nancycats, 8:9), replen=rep(1,9))
## dist9to8.attr <- attr(dist9to8, "Labels")
## nan9to8pop <- nancycats@pop[nancycats@pop %in% c("P08", "P09")]
## attr(dist9to8, "Labels") <- paste(nan9to8pop, dist9to8.attr)
## heatmap(as.matrix(dist9to8), symm=TRUE)


###################################################
### code chunk number 122: popcompare_bruvo2
###################################################
dist9to8 <- bruvo.dist(popsub(nancycats, 8:9), replen=rep(1,9))
dist9to8.attr <- attr(dist9to8, "Labels")
nan9to8pop <- nancycats@pop[nancycats@pop %in% c("P08", "P09")]
attr(dist9to8, "Labels") <- paste(nan9to8pop, dist9to8.attr)
heatmap(as.matrix(dist9to8), symm=TRUE)


###################################################
### code chunk number 123: bruvo_boot (eval = FALSE)
###################################################
## set.seed(1001)
## nan9tree <- bruvo.boot(popsub(nancycats, 8:9), replen=rep(1,9), sample=1000, cutoff=50)


###################################################
### code chunk number 124: progbar
###################################################
cat("|================================================================================| 100%\n")
cat("\nBootstrapping... (note: calculation of node labels can take a while even after the progress bar is full)\n\n")


###################################################
### code chunk number 125: bruvo_tree
###################################################
 set.seed(1001)
 nan9tree <- phangorn::upgma(bruvo.dist(popsub(nancycats, 8:9), replen = rep(1, 9)))
 nan9tree$node.labels <- c(100, NA, NA, NA, NA, 66, NA, NA, NA, NA, NA, 78, NA, 61, NA, NA, NA, NA)
 ape::plot.phylo(nan9tree, show.node.label=TRUE)
 ape::axisPhylo(3)


###################################################
### code chunk number 126: greycurve_normal
###################################################
greycurve()


###################################################
### code chunk number 127: greywidth_inverse
###################################################
greycurve(gweight = 2)


###################################################
### code chunk number 128: greycurve_small_heavy
###################################################
greycurve(glim = c(0.2, 0.9), gadj=15)


###################################################
### code chunk number 129: greywidth_large_heavy
###################################################
greycurve(glim = c(0.2, 0.9), gadj=15, gweight=2)


###################################################
### code chunk number 130: bruvo_msn_code (eval = FALSE)
###################################################
## data(partial_clone)
## set.seed(9005)
## pc.msn <- bruvo.msn(partial_clone, replen=rep(1, 10), vertex.label.cex=0.7, 
##           vertex.label.dist=-0.5, palette=colorRampPalette(c("blue", "yellow")))


###################################################
### code chunk number 131: bruvo_msn
###################################################
data(partial_clone)
set.seed(9005)
pc.msn <- bruvo.msn(partial_clone, replen=rep(1, 10), vertex.label.cex=0.7, 
          vertex.label.dist=-0.5, palette=colorRampPalette(c("blue", "yellow")))


###################################################
### code chunk number 132: bruvo_msn_igraph
###################################################
library(igraph)
pc.msn


###################################################
### code chunk number 133: bruvo_msn_reconstruct (eval = FALSE)
###################################################
## set.seed(9005)
## library(igraph)
## plot(pc.msn$graph, vertex.size = V(pc.msn$graph)$size * 3, vertex.label.cex=0.7, 
##      vertex.label.dist=-0.5,)
## legend(-1.55, 1, bty = "n", cex = 0.75, legend = pc.msn$populations, 
##        title = "Populations", fill = pc.msn$colors, border = NULL)


###################################################
### code chunk number 134: poppr_msn (eval = FALSE)
###################################################
## data(Aeut)
## A.dist <- diss.dist(Aeut)
## set.seed(9005)
## A.msn <- poppr.msn(Aeut, A.dist, vertex.label=NA, palette=rainbow, gadj=15)


###################################################
### code chunk number 135: poppr_msn_fig
###################################################
data(Aeut)
A.dist <- diss.dist(Aeut)
set.seed(9005)
A.msn <- poppr.msn(Aeut, A.dist, vertex.label=NA, palette=rainbow, gadj=15)


###################################################
### code chunk number 136: Aeut_rehash
###################################################
data(Aeut)
poppr(Aeut)


###################################################
### code chunk number 137: Aeut_sample (eval = FALSE)
###################################################
## poppr(Aeut, sample=999, hist=FALSE, quiet=TRUE)


###################################################
### code chunk number 138: Aeut_sample_show
###################################################
none1 <- poppr(Aeut, hist=FALSE, quiet=TRUE)
cbind(none1[1:10],list(p.Ia = rep(0.001, 3)), none1[11], list(p.rD = rep(0.001, 3)), none1[12])


###################################################
### code chunk number 139: Aeut_cc_sub_show (eval = FALSE)
###################################################
## poppr(Aeut, sample=999, clonecorrect=TRUE, hier=c("Pop","Subpop"), 
##       dfname="population_hierarchy", quiet=TRUE, hist=FALSE)


###################################################
### code chunk number 140: Aeut_cc_sub
###################################################
sub1 <- poppr(Aeut, clonecorrect=TRUE, hier=c("Pop","Subpop"), dfname="population_hierarchy", quiet=TRUE, hist=FALSE)
cbind(sub1[1:10],list(p.Ia = rep(0.001, 3)), sub1[11], list(p.rD = rep(0.001, 3)), sub1[12])


###################################################
### code chunk number 141: Aeut_cc_pop_show (eval = FALSE)
###################################################
## poppr(Aeut, sample=999, clonecorrect=TRUE, hier="Pop", 
##       dfname="population_hierarchy", quiet=TRUE, hist=FALSE)


###################################################
### code chunk number 142: Aeut_cc_pop
###################################################
pop1 <- poppr(Aeut, sample=0, clonecorrect=TRUE, hier="Pop", dfname="population_hierarchy", quiet=TRUE, hist=FALSE)
cbind(pop1[1:10],list(p.Ia = rep(0.001, 3)), pop1[11], list(p.rD = rep(0.001, 3)), pop1[12])


###################################################
### code chunk number 143: nancy_example_show (eval = FALSE)
###################################################
## set.seed(2001)
## poppr(nancycats, sublist=5:6, total=FALSE, sample=999, method=3, quiet=TRUE, hist=FALSE)


###################################################
### code chunk number 144: nancy_example_eval
###################################################
set.seed(2001)
nan_ex <- poppr(nancycats, sublist=5:6, total=FALSE, sample=0, method=3, quiet=TRUE, hist=FALSE)
cbind(nan_ex[1:10],list(p.Ia = c(0.599, 0.064), nan_ex[11], list(p.rD = c(0.599, 0.065), nan_ex[12])))


###################################################
### code chunk number 145: simulated_dist (eval = FALSE)
###################################################
## set.seed(2004)
## poppr(system.file("files/simulated.dat", package="poppr"), sample=999, method=1, quiet=TRUE)


###################################################
### code chunk number 146: simulated_dist_out
###################################################
set.seed(2004)
sim_ex <- poppr(system.file("files/simulated.dat", package="poppr"), sample=0, method=1, quiet=TRUE, hist=FALSE)
cbind(sim_ex[1:10],list(p.Ia = c(0.09), sim_ex[11], list(p.rD = c(0.09), sim_ex[12])))


###################################################
### code chunk number 147: poppr_manual.Rnw:1714-1718
###################################################
df <- data.frame(list(locus1=c("101/101", "102/103", "102/102"), 
                      locus2=c("201/201","202/203","203/204"), 
                      locus3=c("301/302", "301/303", "304/305")))
df


###################################################
### code chunk number 148: poppr_manual.Rnw:1721-1722
###################################################
dfg@tab[, 1:3]


###################################################
### code chunk number 149: poppr_manual.Rnw:1731-1734
###################################################
abs(dfg@tab[1, 1:3] - dfg@tab[2, 1:3])
abs(dfg@tab[1, 1:3] - dfg@tab[3, 1:3])
abs(dfg@tab[2, 1:3] - dfg@tab[3, 1:3])


###################################################
### code chunk number 150: ggsave1 (eval = FALSE)
###################################################
## data(nancycats) # Load the data set.
## poppr(nancycats, sublist=5, sample=999) # Produce a single plot.
## ggsave("nancy5.pdf")


###################################################
### code chunk number 151: png_save (eval = FALSE)
###################################################
## data(H3N2)
## pop(H3N2) <- H3N2$other$x$country
## ####
## png("H3N2_barchart%02d.png", width = 14, height = 14, units = "in", res = 300)
## H.tab <- mlg.table(H3N2)
## dev.off()
## ####


###################################################
### code chunk number 152: pdf_save (eval = FALSE)
###################################################
## pdf("H3N2_barcharts.png", width = 14, height = 14, compress = FALSE)
## H.tab <- mlg.table(H3N2)
## dev.off()


