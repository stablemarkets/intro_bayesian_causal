/*
Author: Arman Oganisian
Purpose: Implement partial pooling example in SAS.
*/


/*  			Load data simulated in ..SAS_example/simulate_data.R  				*/
proc import datafile="C:\Users\aoganisi\Box Sync\Research\Analyses\IntroBayesCausal\partial_pool\SAS_example\simulated_data.csv"
	out=d
	dbms=CSV;
run;
proc freq data=d; 
	table V;
run;

/* create design matrix 	*/
proc transreg data=d design;
	model class( V / zero='1');  /* categorical, with A=0 as reference */
	id  Y W A; 					/* include out come in this dataset */
	output out=mcmc_data(drop=_:);
run;

/*  			Run MCMC 								 				*/

proc mcmc data=mcmc_data 
	/*run 10,000 iterations after 5,000 iteration burn-in. results stored in 
	  "posterior draws" dataset. 
	  Only monitor Psi - the causal effects of interest
      Output 2.5, 97.5 percentiles of posterior for credible interval. */
	outpost=posterior_draws seed=1 nmc=1000 nbi=1000 monitor=(mu) STATS(PERCENTAGE=(2.5 97.5))=SUMMARY ;

	/* Declare parameters */
	parms b0 bw bv2-bv5 t1-t5 mu;
	
	/* Set priors for nuisance pars */
	prior b0 ~ normal(0, sd=1);
	prior bw ~ normal(0, sd=1);
	prior bv: ~ normal(0, sd=1);
	
	/* partial pooling priors on V effects  */
	/* overall mean */
	prior mu ~ normal(0, sd=1);
	/* V-specific effects around mu */
	prior t1 ~ normal(mu, sd=.5);
	prior t2-t5 ~ normal(mu-t1, sd=.5);
	
	/* Specify Likelihood */
	eta = logistic( b0 + bw*W + bv2*V2 + bv3*V3 + bv4*V4 + bv5*V5 + (t1 + t2*V2 + t3*V3 + t4*V4 + t5*V5)*A );
	model Y ~ binomial(1, eta );

run;

/* Compute Bayesian Bootstrap Estimate of causal OR for each stratum of V */
proc iml;
	/* set seed */
	call randseed(1);

	/* Read in matrix of posterior draws from SAS to IML */
	use posterior_draws; read all var _ALL_ into pm; close posterior_draws;

	/* Read in model matrix from SAS to IML */
	use mcmc_data; read all var _ALL_ into X; close mcmc_data;

	n = nrow(X);
	n_iter = nrow(pm);

	/* shell to store posterior draws of Causal OR for each of the 5 strata */
	OR_mat = j(n_iter, 5, 0); 

	/* loop over posterior draws */
	do i=1 to n_iter; 
	
		/* loop over strata of V */
		do v = 1 to 5;
			nv = sum(X[,6] = v); /* find how many subjects in stratum v */
			idx = loc(X[,6] = v) ; /* find which obs are in stratum v */
			/* Draw from Dirichlet(1,...,1) distribution  
			   to do bayesian bootstrap estimate of P_v(W) */
			alpha= J(nv , 1 , 1);
			bb_w = RandDirichlet(1, alpha);
			bb_w = bb_w || 1-sum(bb_w);

			/* for each strata, compute logit of event 
			   under treatment 1 and 0: lp1, lp0*/

			/* compute reference group v=1 separately */
			if v=1 then do;
				lp1 = pm[i,2] + pm[i, 3]*X[ idx , 8] + pm[i, 8] ;
				lp0 = pm[i,2] + pm[i, 3]*X[ idx , 8] ;
			end;
			
			if v>1 then do;
				lp1 = pm[i,2] + pm[i, 3]*X[ idx , 8] + pm[i,2+v]  + (pm[i, 8] + pm[i,7+v] ) ;
				lp0 = pm[i,2] + pm[i, 3]*X[ idx , 8] + pm[i,2+v] ;
			end;
			
			/* inverse logit transform to convert to probability */
			p1 = exp(lp1)/(1+exp(lp1)); 
			p0 = exp(lp0)/(1+exp(lp0));
			
			/* bayesian bootstrap average of probability*/
			/* dot-product: bb_w is 1-X-n vector and p1, p0 are n-X-1  */ 
			mu1 = bb_w*p1; 
			mu0 = bb_w*p0;
			
			/* compute Odds Ratio for stratum v */
			OR_mat[i,v] = ( mu1/(1-mu1) ) / ( mu0/(1-mu0 ) )  ;
		end;

	end;
	
	/* Compute posterior means and 95% intervals for each of the 5 ORs */
	OR_means=mean(OR_mat);
	call qntl(OR_perc, OR_mat,{.025, .975});
	OR_out = t(OR_means // OR_perc) ;

	/* output the results from IML back to SAS */
	create iml_res from OR_out; append from OR_out ; close iml_res;
quit;

data OR_plot; set iml_res; 
	rename COL1 = OR_mean;
	rename COL2 = OR_lwr;
	rename COL3 = OR_upr;
	V= _n_;
run;

/* Estiamtes in OR_plot match Stan's estimates pretty closely. 
Point estimates are very close. Interval estimates are pretty close - 
but we likely need more samples before they start converging.
Small diffferences expected due to randomness. */ 

/* plot posterior mean curve along with 95% interval */
proc sgplot data=OR_plot;
	scatter x=V y=OR_mean / yerrorlower=OR_lwr yerrorupper=OR_upr ; 
run;
