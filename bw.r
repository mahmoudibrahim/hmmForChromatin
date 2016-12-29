rm(list = ls())
suppressPackageStartupMessages(library("mvtnorm"))
suppressPackageStartupMessages(library("parallel"))
options(scipen = 1000)


###define variables###
samplingSeed = 666 #seed used for randomized steps (used in data transformation, see Duttke et al. Mol. Cell, 2015.)
out = "" #path to output folder
inF = "" #path to input folder
marks = c("h3k4me1","h3k4me2","h3k4me3","h3k27ac") #which histone modifications to learn the HMM
clustNum = 2:10 #number of chromatin states to attempt
resolution = 10 #the resolution the wig files are produced with in basepairs
minimumSeqLen = 500 #minimum length in basepairs of a genome segment to be used in learning the HMM
chromo = c("chr1","chr2") #which chromosomes to use for the HMM learning
maxIter = 200 #maximum number of iterations to learn the HMM
cornum = 10 #number of cores to use during learning the HMM
color = colorRampPalette(c("white","blue"))(n=1000) #color palette to produce the chromatin state heatmap

##note: histone modification wig files should be arranged as follows (see line 145 of this script):
##inF "/" marks[f] "-" chromo[c] ".wig"
#########



###############
#normalize 0-to-1
score = function(x){
	x  = ((x-min(x))/(max(x)-min(x)))*1000
	return(x)
}
seqle = function(x) { 
	incr = 1
	n = length(x)
	y = x[-1L] != x[-n] + incr 
	i = c(which(y|is.na(y)),n)
	list(lengths = diff(c(0L,i)), values = x[head(c(0L,i)+1L,-1L)]) 
}
#read-in files
prasein = function(bedfile, chromo) {
	dat = vector('numeric')
	for (r in 1:length(chromo)) {
		dat = c(dat, read.table(bedfile[[r]])[[4]])
	}
	dat[dat < 0] = 0
	return(dat)
}
#data transformation
makeguass = function(mark) {
	l = length(mark[mark == 0])
	set.seed(samplingSeed)
	minis = sample(mini, l, replace = TRUE)
	mark[mark == 0] = minis
	mark = score(mark)
	mark = log(mark + 0.001)
	return(mark)
}
#get multivariate density
getDensityTemp = function(index, means, sds, dat, models) {
	x = index
	den = (dmvnorm(dat, means[,x], sigma = sds[[x]], log = TRUE))
}
#resturture a list to a matrix
makeMatrix = function(someList, numTracks) {
	return(matrix(someList, ncol = numTracks))
}
#forward algorithm
getfor = function(whichiswhich, datTemp, probsAll, seqseq, models, wei, trans) {
	
	dd = datTemp[[whichiswhich]]
	
	starting = (cumsum(seqseq[1:(whichiswhich+1)]))[whichiswhich] + 1
	ending = (cumsum(seqseq[1:(whichiswhich+1)]))[(whichiswhich+1)]
	normaF = rep(0, nrow(datTemp[[whichiswhich]]))
	probs = probsAll[(starting : ending), ]
	probsF = probs
	
	a = 1
	for (b in 1:models) {
		probsF[a,b] = probs[a,b] + ( log( sum( exp(wei) ) ) )
	}
	normaF[a] = log( sum( exp(probsF[a,]) ) )
	probsF[a,] = probsF[a,] - normaF[a]
	for (a in 2:(nrow(probs))) {
		for (b in 1:models) {
			probsF[a,b] = probs[a,b] + ( log( sum( exp(probsF[a-1,] + trans[b,]) ) ) )
		}
		normaF[a] = log( sum( exp(probsF[a,]) ) )
		probsF[a,] = probsF[a,] - normaF[a]
	}
	
	writethis = cbind(probsF, normaF)
	return(writethis)
}
#forward-backward algorithm
getforback = function(whichiswhich, datTemp, probsAll, seqseq, models, wei, trans) {
	
	dd = datTemp[[whichiswhich]]
	
	start = cumsum(seqseq[1:(whichiswhich+1)])[whichiswhich] + 1
	end = cumsum(seqseq[1:(whichiswhich+1)])[(whichiswhich+1)]
	
	normaF = rep(0, nrow(datTemp[[whichiswhich]]))
	probs = probsAll[(start : end), ]
	probsF = probs
	probsB = probs
	
	a = 1
	for (b in 1:models) {
		probsF[a,b] = probs[a,b] + ( log( sum( exp(wei) ) ) )
	}
	normaF[a] = log( sum( exp(probsF[a,]) ) )
	probsF[a,] = probsF[a,] - normaF[a]
	for (a in 2:(nrow(probs))) {
		for (b in 1:models) {
			probsF[a,b] = probs[a,b] + ( log( sum( exp(probsF[a-1,] + trans[b,]) ) ) )
		}
		normaF[a] = log( sum( exp(probsF[a,]) ) )
		probsF[a,] = probsF[a,] - normaF[a]
	}
	
	a = nrow(probs)
	probsB[a,] = log(sum(exp(probsB[a,] + probs[a,])))
	for (a in ((nrow(probs))-1):1) {
		for (b in 1:models) {
			probsB[a,b] = log(sum(exp(probsB[a+1,] + probs[a+1,] + trans[b,])))
		}
		probsB[a,] = probsB[a,] - normaF[a]
	}
	probs = probsF + probsB

	writethis = cbind(probs, normaF)
	return(writethis)
}
###############



###read in files###
message("Importing Data...")
files = list();
for (f in 1:length(marks)) {
	anathema = list()
	for (c in 1:length(chromo)) {
		anathema[[c]] = paste0(inF, "/", marks[f], "-", chromo[c], ".wig")
	}
	files[[f]] = anathema
}
dat = mclapply(files, prasein, chromo, mc.cores = cornum, mc.preschedule = FALSE)
datanr = length(files)
minimumSeq = minimumSeqLen / resolution
########

###Process Data###
message("Processing Data...")
dat = matrix(unlist(dat), ncol = datanr, byrow = FALSE)
zeros = rowSums(dat)
dat = cbind(zeros, dat)
dat = dat[dat[,1] != 0,, drop = FALSE]
dat = dat[,-1]
zeros = which(zeros != 0)
#########

###make the data look gaussian###
mmV = 0.1
mmM = 0.1
sigma = log(1+((mmV) / ((mmM)^2)))
mu = (log(mmM)) - (0.5 * (sigma))
set.seed(samplingSeed)
mini = rlnorm(5000000, mu, sqrt(sigma))
dat = apply(dat, 2, makeguass)
#########


for (m in 1:(length(clustNum))) {

models = clustNum[m]

##construct dataframe with training sequences
saltlines = seqle(zeros)
violins = which(saltlines$lengths > minimumSeq)
datTemp = split(dat, rep(saltlines$values, saltlines$lengths))
datTemp = datTemp[violins]
datTemp = mclapply(datTemp, makeMatrix, numTracks = length(marks), mc.cores = cornum, mc.preschedule = TRUE)
datTempMat = do.call(rbind, datTemp)
seqseq = saltlines$lengths[violins]
seqseq = c(0,seqseq)
message(paste0("Number of States: ", models))
message(paste0("Number of Training Sequences: ", length(datTemp), ", with median length of ", median(seqseq)))
##########





##initiatlize with kmeans
message("Initializing with K-means...")
means = matrix(0, nrow = length(marks), ncol = models)
sds = list(matrix(0, nrow = length(marks), ncol = length(marks)))
sds = rep(sds, models)
wei = rep(0, models)
trans = diag(rep(0.99, models))
trans[trans == 0] = (0.01 / (models - 1))
trans = log(trans)
datInd = 1:nrow(datTempMat)
datInd = dat[datInd,]
set.seed(samplingSeed)
k = (kmeans(datInd, models, iter.max = 10, nstart = 10))$cluster
for (x in 1:models) {
	means[,x] = apply(datInd[(which(k == x)),], 2, mean)
	sds[[x]] = cov(datInd[(which(k == x)),], method = "pearson")
	wei[x] = (length(which(k == x))) / (length(k))
}
pdf(paste0(out, "/means-", models, "-NT.pdf"))
a = t(apply(means, 1, score))
b = heatmap(a, Rowv = NA, Colv = NA, scale = "none", col = color, labRow = marks, cexRow = 2, cexCol = 2, margins = c(7,10))
b = dev.off()


probsAll = mclapply(1:models, getDensityTemp, means = means, sds = sds, dat = datTempMat, models = models, mc.cores = cornum)
probsAll = matrix(unlist(probsAll), ncol = models, byrow = FALSE)
probsAll = mclapply(1:(length(datTemp)), getfor, datTemp, probsAll, seqseq, models, wei, trans, mc.preschedule = TRUE, mc.cores = cornum)
probsAll = do.call(rbind, probsAll)
normaAll = probsAll[,(models + 1)]
probsAll = probsAll[,1:models]
probsAll = exp((probsAll) - (log(rowSums(exp(probsAll)))))
likelihood = sum(normaAll)
message(paste0("Initialization, Log-Likelihood: ", likelihood))
#############




message("Started EM...")
llPrev = likelihood

for (i in 1:maxIter) {


##
for (x in 1:models) {
	whoiswho = cov.wt(datTempMat, wt = probsAll[,x], cor = FALSE, center = TRUE, method = "ML")
	means[,x] = whoiswho$center
	sds[[x]] = whoiswho$cov
	wei[x] = (sum(probsAll[,x])) / (length(datTempMat[,1]))
}
########



#########
probsAll = mclapply(1:models, getDensityTemp, means = means, sds = sds, dat = datTempMat, models = models, mc.cores = cornum)
probsAll = matrix(unlist(probsAll), ncol = models, byrow = FALSE)
probsAll = mclapply(1:(length(datTemp)), getforback, datTemp, probsAll, seqseq, models, wei, trans, mc.preschedule = TRUE, mc.cores = cornum)
probsAll = do.call(rbind, probsAll)
normaAll = probsAll[,(models + 1)]
probsAll = probsAll[, 1:models]
probsAll = exp((probsAll) - (log(rowSums(exp(probsAll)))))
likelihood = sum(normaAll)
message(paste0("Iteration ", i, ", Log-Likelihood: ", likelihood))
###############


pdf(paste0(out, "/means-", models, "-NT.pdf"))
a = t(apply(means, 1, score))
heatmap(a, Rowv = NA, Colv = NA, scale = "none", col = color, labRow = marks, cexRow = 2, cexCol = 2, margins = c(7,10))
dev.off()


llDiff = ((llPrev - likelihood) / (likelihood))
if (llDiff < 0.000001) {
	break;
}

llPrev = likelihood
	
}


message("Saving Model to Disk...")
model = list()
model$means = means
model$covariance = sds
model$likelihood = likelihood
model$wei = wei
save(model, file = paste0(out, "/hmm-", models, ".RData"))


message("Done!")
}
