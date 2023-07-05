-----------------------------------------
# Distribution of Mutational Effects Code
-----------------------------------------

This repository contains code for calculating the likelihood of observations of colony growth rates in [Spontaneous single-nucleotide substitutions and microsatellite mutations have distinct distributions of fitness effects](https://www.biorxiv.org/content/10.1101/2023.07.04.547687v1).

The following abbreviations refer to models described in the text:

1. **poissonmut**: 'SNMs only' model
1. **lambdafit**: SNMs + unidentified mutations, ‘single DFE’ model
1. **mixedgauss**: SNMs + unidentified mutations, ‘Gaussian’ model
1. **mixed_digamma**: SNMs + unidentified mutations, ‘two-gamma’ model

+ `pre_MLE...` and `post_MLE...` functions perform functions that occur before and after optimization, respectively, that don't need to happen in each round of optimization (e.g. loading in data and organizing inputs in the case of `pre_MLE...` functions, organizing multiple outputs in the case of `post_MLE...` functions). These are essentially used for 'housekeeping' and are not necessary for understanding the implementation of the math described in the paper
+ `LL_...` functions calculate likelihoods of data; these are named directly in the setup files for each model. Note that some likelihood calculation functions perform optimization calling functions from [MutationEffectEstimation](https://github.com/plavskin/MutationEffectEstimation).
	- [`LL_mixef_calc`](LL_mixef_calc.m) performs the likelihood calculation for **Equation 5**
	- [`LL_calculator_within_pair_different_sigmas`](LL_calculator_within_pair_different_sigmas.m) calculates the likelihood of observed growth rate differences described in **Equation 8**; however, this code is run indirectly in order to allow for likelihood maximization across general parameters and multiple strains
		+ for estimation of mutational effect parameters relative to a single in-well control, [`LL_calculator_strain_looper_pairwise`](LL_calculator_strain_looper_pairwise.m) is used
		+ for estimation of mutation effect parameters relative to a non-marked strain present in a subset of wells, [`LL_calculator_strain_looper_pairwise_multiref_with_mixef`](LL_calculator_strain_looper_pairwise_multiref_with_mixef.m) is used
	- [`LL_calculator_me_estimate_single_dfe`](LL_from_mut_effect_estimate/LL_calculator_me_estimate_single_dfe.m) performs the likelihood calculation whose fourier transform is described in **Equation 16**
	- [`LL_calculator_observed_diffs_given_me_dist`](LL_calculator_observed_diffs_given_me_dist.m) performs the likelihood calculation in **Equation 17**; however, this code is run indirectly in order to allow for likelihood maximization across general parameters and multiple strains, through[`LL_calculator_single_DFE`](LL_calculator_single_DFE.m)
+ `sim...` functions (e.g. [`sim_mut_effect_poissonmut`](DFE_poissonmut/sim_mut_effect_poissonmut.m)) simulate mutational effect or growth rate data from a set of parameters (dimensionality of data simulated matches provided input data)
+ `nonlin_constraint...` functions (e.g. [`nonlin_constraint_poissonmut`](DFE_poissonmut/nonlin_constraint_poissonmut.m)) calculate non-linear constraints on parameter values to be used during optimization
+ `fourier_dfe_...` functions (e.g. [`fourier_dfe_poissonmut`](DFE_poissonmut/fourier_dfe_poissonmut.m)) calculate the characteristic function for a distribution of mutational effects

There are also a number of accessory functions used for e.g. safely converting between fourier space and mutational effect space, etc.