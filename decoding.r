rm(list = ls())
suppressPackageStartupMessages(library("mvtnorm"))
suppressPackageStartupMessages(library("parallel"))
options(scipen = 1000)


###define variables###
samplingSeed = 666 #seed used for randomized steps (used in data transformation, see Duttke et al. Mol. Cell, 2015.)
out = "" #path to output folder
inF = "" #path to input folder
resol = 10 #resolution used to produce the wig files in basepairs
marks = c("h3k4me1","h3k4me2","h3k4me3","h3k27ac") #which histone modifications were used to learn the HMM
inmodel = "" #path to RData file with the model (output of bw.r)
chromo = c("chr1","chr2","chr3") #which chromosomes to segment
cornum = 10 #how many cores to use

##note: histone modification bedGraph files should be arranged as follows (see line 180 of this script):
##inF "/" marks[f] "-" chromo[c] ".bedGraph"
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
	dat = read.table(bedfile)[[4]]
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
		probsF[a,b] = probs[a,b] + ( log( sum( exp(wei + trans[b,]) ) ) )
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
	starting = cumsum(seqseq[1:(whichiswhich+1)])[whichiswhich] + 1
	ending = cumsum(seqseq[1:(whichiswhich+1)])[(whichiswhich+1)]
	normaF = rep(0, nrow(datTemp[[whichiswhich]]))
	probs = probsAll[(starting : ending), ]
	probsF = probs
	probsB = probs
	
	
	if (nrow(dd) > 1) {
	

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
		
	} else {
		
		probs = t(as.matrix(probs))
		a = 1
		probsF = probs

		for (b in 1:models) {
			probsF[a,b] = probs[a,b] + ( log( sum( exp(wei) ) ) )
		}
		normaF[a] = log( sum( exp(probsF[a,]) ) )
		probsF[a,] = probsF[a,] - normaF[a]
		writethis = cbind(probs, normaF)
		return(writethis)
	}
}
whichmax = function(go, ppp) {
	return(which.max(ppp[go,]))
}
###############


#######
message(paste0("Loading model..."))
load(inmodel)
means = model$means
sds = model$covariance
likelihood = model$likelihood
wei = model$wei
models = length(wei)
trans = diag(rep(0.9, models))
trans[trans == 0] = (0.1 / (models - 1))
trans = log(trans)
#######







##############
for (xyz in 1:length(chromo)) {

chromName = paste0(chromo[xyz])


###read in files###
message(paste0("Importing Data for Chromosome ", chromName, "..."))
files = list();
for (f in 1:length(marks)) {
	anathema = paste0(inF, "/", marks[f], "-", chromName, ".bedGraph")
	files[[f]] = anathema
}
dat = mclapply(files, prasein, marks, mc.cores = cornum, mc.preschedule = FALSE)
datanr = length(files)
########

###Process Data###
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



##construct dataframe with training sequences
saltlines = seqle(zeros)
violins = which(saltlines$lengths > 0)
datTemp = split(dat, rep(saltlines$values, saltlines$lengths))
datTemp = datTemp[violins]
datTemp = mclapply(datTemp, makeMatrix, numTracks = length(marks), mc.cores = cornum, mc.preschedule = TRUE)
datTempMat = do.call(rbind, datTemp)
seqseq = saltlines$lengths[violins]
seqseq = c(0,seqseq)
message(paste0("Number of Sequences: ", length(datTemp), ", with median length of ", median(seqseq)))
##########




####
message("Started Decoding...")
probsAll = mclapply(1:models, getDensityTemp, means = means, sds = sds, dat = datTempMat, models = models, mc.cores = cornum)
probsAll = matrix(unlist(probsAll), ncol = models, byrow = FALSE)
probsAll = mclapply(1:length(datTemp), getforback, datTemp, probsAll, seqseq, models, wei, trans, mc.preschedule = TRUE, mc.cores = cornum)
probsAll = do.call(rbind, probsAll)
normaAll = probsAll[,(models + 1)]
probsAll = probsAll[, 1:models]
probsAll = exp((probsAll) - (log(rowSums(exp(probsAll)))))
###############






message("Saving Information to Disk...")
why = unlist(mclapply(1:nrow(probsAll), whichmax, probsAll, mc.preschedule = TRUE, mc.cores = cornum))
kc = cbind(why, probsAll, normaAll)
aa = rep(chromName, length(zeros))
bb = ((zeros*resol) - (resol - 1) - 1)
cc = (zeros*resol)
write.table(cbind(aa,bb,cc,why), sep = "\t", file = paste0(out, "/model-", chromName, "-", models, ".probs"), append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE)




}
message("Done!")
