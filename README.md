# hmmForChromatin

Those are R scripts implementing the Baum-Welch algorithm and posterior decoding for chromatin state annotation. The scripts take as an input wig files of histone modifications split by chromosome. It returns an R file with the learned model and the segmentation of the states.

The pipeline to create wig files suitable for those scripts is described in Duttke et al. Mol. Cell, 2015.
