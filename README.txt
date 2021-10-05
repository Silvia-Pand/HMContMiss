This repository contains the R functions used for the paper
"Maximum likelihood estimation of hidden Markov models for
continuous longitudinal data with missing responses and dropout" by
- S.Pandolfi (University of Perugia, IT)  
- F.Bartolucci (University of Perugia, IT)
- F. Pennoni (University of Milano-Bicocca, IT)

lmbasic.cont.MISS.R 	--->	estimate the basic HM model for continuous outcomes with intermittent missingness and dropout using an extended EM algorithm 


lmcovlatent.cont.MISS.R --->	estimate the HM model for continuous outcomes with intermittent missingness and dropout including covariates in the distribution of the latent process
 
bootstrap.MISS.R --> perform non-parametric bootstrap procedure in order to compute standard errors of model parameters for the HM model with covariates
 	
lk_comp_cont_MISS.R ---> compute complete log-likelihood of the basic HM model for continuous outcomes (internal use)

lk_comp_latent_cont_MISS.R ---> compute the complete log-likelihood of the HM model for continuous outcomes with covariates in the distribution of the latent process (internal use)

lk_obs_latent_cont_MISS.R --->  compute the observable log-likelihood of the HM model with covariates in the latent model (internal use)

prob_multilogit.R ---> compute multinomial probabilities (internal function)

est_multilogit.R ---> perform maximum likelihood estimation of the multilogit model (internal function)

prob_post_cov_cont.R ---> use backward recursion to compute posterior probabilities (internal funtion)

example.R ---> example file that loads the workspace file "example_data.RData" and fits the basic HM model and the HM with covariates. It also perform non-parametric bootstrap using function bootstrap.MISS.R  				


example_data.RData  --> workspace file containing a random subsample from PBC data with 10 time occasion and 3 covariates 
