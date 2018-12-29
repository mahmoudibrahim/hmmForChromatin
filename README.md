# hmmForChromatin

Those are R scripts implementing the Baum-Welch algorithm (initialized with K-means) and posterior decoding for chromatin state annotation. The scripts take as an input bedGraph files of histone modifications split by chromosome. It returns an R file with the learned model and the segmentation of the states.

If you use those scripts in a publication, please cite:
[Ibrahim et al., Nature Communications 2018](https://doi.org/10.1038/s41467-018-06962-z) (DOI: 10.1038/s41467-018-06962-z).

The main idea behind this approach to chromatin state annotation is to obtain signal for each ChIP-Seq data set only where there are peaks, and to clamp everywhere else at zeros. This allows for high-resolution chromatin state annotation using HMMs (similar to [ChromHMM](http://compbio.mit.edu/ChromHMM/) and [Segway](https://segway.hoffmanlab.org/)), taking into account differences in signal levels and the covariance between different histone modifications, but without requiring large memory or compute power resources. You can read more about this idea in [the supplemental methods section of Duttke et al. Molecular Cell 2015](https://ars.els-cdn.com/content/image/1-s2.0-S1097276514010077-mmc1.pdf), where this approach for chromatin state annotation was first introduced.

This approach was developed at the Ohler Lab, at [BIMSB, MDC in Berlin](https://www.mdc-berlin.de/bimsb).
