# README

This is a slight modification of bitbucket.org/afmagee/tree_convergence_code/

This repository is best interacted with via the Rproject.
The code in this repository depends on a number of R packages:

- treess (from [bitbucket.org/afmagee/treess](https://bitbucket.org/afmagee/treess/src/master/))
- coda
- ape
- phangorn
- viridis (for colors in several plots, not strictly necessary)

## simulation_study
Code that tests a variety of ESS metrics on a variety of realistic fake datasets, and that makes summary figures.

### data
To save space, does not contain the original data files, even tarred and gzipped these are too big for uploading.
Instead, contains files with only 95% HPD or first 4096 trees, whichever is smaller.
Trees are product of MrBayes Golden Runs of 1 billion total iterations, see original paper ([Whidden et al.,2020](https://academic.oup.com/sysbio/article/69/2/280/5555780)) for details.

### output
Currently contains the adjacency graphs for each dataset, named `<dataset>_adjacency_graph.Rdata`, required to run the "fake" MCMC.
Provided for those who wish to play with the examples without installing igraph.

Running the MCMC step produces files named in the format `<dataset>_MCMC_<run_length>_<sampling_frequency>.Rdata`.
All runs contain 100 independent chains, and are formatted as objects of class simulatedPosterior.

Computing the Monte Carlo error produces files named in the format `<dataset>_error_<run_length>_<sampling_frequency>.Rdata`.

Some summarizing and plotting requires computing the estimated tree/split probabilities, which produces files named in the format `<dataset>_probs_<run_length>_<sampling_frequency>.Rdata`.

### src
Scripts for running analyses and post-processing.
All scripts are designed to be callable from the top level of the repo with Rscript.

The more computationally intensive steps are running the "fake" MCMC and computing the raw error measures needed for Monte Carlo error calculations.
These steps come as a pair of scripts, one of which has `_slurm` in the name.
The `_slurm` scripts automatically loop over the entire pipeline, and thus contain all parameters, settings, and seeds used in the analyses in the paper.
They are designed to be called with Rscript inside a slurm environment (such as from a slurm scriptfile).

Running the full pipeline, all 45 dataset by run length combinations used in the paper, 100 replicate MCMC chains each, takes a long time and is not advised without using a cluster.
For debugging/testing/experimentation purposes, running the MCMC for only 1000 iterations, reducing the number of chains (but no less than 2), and possibly only running DS3, are recommended to speed up the process.

#### Simulating
The full simulation pipeline requires calling the scripts in the following order.
The product is the Monte Carlo error files, which can be used to determine how well the ESS measures worked, and which include ESS estimates.

- [condense_trees.R](simulation_study/src/condense_trees.R) Takes posterior tree samples and drops down to the 95% HPD (or first 4096 tree topologies, whichever is smaller). Then stores in a simple file without any branch lengths.
- [generate_adjacency_graphs.R](simulation_study/src/generate_adjacency_graphs.R) Takes in the best trees files and constructs the adjacency graph required for running fake MCMC.
- [run_pseudo_MCMC.R](simulation_study/src/run_pseudo_MCMC.R) Uses the adjacency graphs to simulate phylogenetic MCMC on the DS for a variety of run lengths. The `run_pseudo_MCMC_slurm.R` script is convenient for cluster job submission and contains all parameter settings used.
- [compute_errors.R](simulation_study/src/compute_errors.R) Computes the brute-force Monte Carlo error of the MCMC runs. These results are minimally processed. Specifically, for each ESS measure for each dataset and each run length, this script computes the ESS of each chain, draws iid samples, and computes error metrics on a per-chain and per-split (or per-tree) basis for both the MCMC and iid samples. The `compute_errors_slurm.R` script is convenient for cluster job submission and contains all parameter settings used.

#### Post-processing
- [compute_probs.R](simulation_study/src/compute_probs.R) The error computing step doesn't record probabilities of trees, so any later use of them (filtering, coloring points by, etc.) requires running this script to obtain them.
- [prepare_for_plotting.R](simulation_study/src/prepare_for_plotting.R) Aggregates the per-split/tree errors across each chain, creating per-split/tree Monte Carlo errors. Then puts all these into a list for plotting.
- [plot_main_results.R](simulation_study/src/plot_main_results.R) Makes the individual panels in all main and supplemental Monte Carlo error figures.


#### Misc
- [utils.R](simulation_study/src/utils.R) Contains miscellaneous utility functions.

## Normal_random_walk
Code for using Normal distributions in 1 or more dimensions.
These examples are used to understand features of the tree ESS measures and the simulation study setup.
The ESS measures used in the study are designed for trees but almost all are applicable to any object.
By recapitulating the simulation study approach for the mean of a single Normal random variable, we can establish a reference for how well standard approaches work to contextualize tree ESS performance.
And by running MCMC targeting a multivariate Normal and using the distance-based ESS measures, we can see how they scale with the sample size.

### src
- [compute_reference_performance.R](Normal_random_walk/src/compute_reference_performance.R) Recapitulates the simulation study approach for a 1-D Normal distribution and simple univariate ESS.
- [compute_mvn_ess.R](Normal_random_walk/src/compute_mvn_ess.R) Runs MCMC targeting a diagonal 10-D Normal, computes distance-based tree ESS measures. Iterates over a number of sample sizes.
- [compute_univariate_ess.R](Normal_random_walk/src/compute_univariate_ess.R) Runs MCMC targeting a diagonal 10-D Normal, computes classical univariate ESS measure separately on each dimension. Iterates over a number of sample sizes.
- [plot_reference_performance.R](Normal_random_walk/src/plot_reference_performance.R) Plots results of `compute_reference_performance.R`.
- [plot_mvn_ess.R](Normal_random_walk/src/plot_mvn_ess.R) Plots results of `compute_mvn_ess.R`.
- [plot_univariate_ess.R](Normal_random_walk/src/plot_univariate_ess.R) Plots results of `compute_univariate_ess.R`.
- [simpleMVRW.R](Normal_random_walk/src/simpleMVRW.R) MCMC on diagonal multivariate normals.

## MCMC_validation
Code for some preliminary testing of the fake phylogenetic MCMC algorithm.
Test cases include both fully synthetic and empirical-based analyses.

The fully synthetic datasets use the set of all (15) rooted 4-tip trees, one analysis assuming all trees have equal probability 1/15, and one with dirichlet-distributed probabilities.

The empirical test cases are the datasets from Scantlebury et al., run for 100,000 generations.
This assumes that the adjacency graphs have already been made (see section "Simulation study").
The empirical target distributions are natural posterior distributions with potential multimodalities and trees that may only be barely connected ([Whidden and Matsen, 2015](https://academic.oup.com/sysbio/article/64/3/472/1632660), [Whidden et al.,2020](https://academic.oup.com/sysbio/article/69/2/280/5555780)), so some noise is to expected in the estimates of the tree probabilities with only 100,000 samples taken.

In all cases, we run 100 chains and examine the mean tree probability for each tree, as well as the [5%,95%] quantile range.
Inspection of performance is done by making figures plotting true versus MCMC estimated tree probabilities.
Plots have means as "x"s, and ranges as bars (colored by whether they include the true value).

### src
- [4_taxon.R](MCMC_validation/src/5_taxon.R) Runs MCMC and summarizes/plots results for analyses targeting the 5-tip trees.
- [real_datasets.R](MCMC_validation/src/real_datasets.R) Runs MCMC on real datasets and plots/summarizes results.
- [general_plots.R](MCMC_validation/src/general_plots.R) Functions for plotting results.
