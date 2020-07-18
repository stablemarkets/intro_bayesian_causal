/*
Author: Arman Oganisian 
Purpose: implement dose effect analysis in SAS. 
*/


/*  			Load data simulated in ..SAS_example/simulate_data.R  				*/
proc import datafile="C:\Users\aoganisi\Box Sync\Research\Analyses\IntroBayesCausal\dose_response\sas_example\simulatd_dose_data.csv"
	out=d
	dbms=CSV;
run;

/* create design matrix 	*/
proc transreg data=d design;
	model class( A / zero='0');  /* categorical, with A=0 as reference */
	id  Y L; 					/* include out come in this dataset */
	output out=mcmc_data(drop=_:);
run;

proc mcmc data=mcmc_data 
	/*run 10,000 iterations after 5,000 iteration burn-in. results stored in 
	  "posterior draws" dataset. 
	  Only monitor Psi - the causal effects of interest
      Output 2.5, 97.5 percentiles of posterior for credible interval. */
	outpost=posterior_draws seed=1 nmc=10000 nbi=5000 monitor=(Psi) STATS(PERCENTAGE=(2.5 97.5))=SUMMARY ;

	/* shell to store Psi draws in */
	array Psi[9];

	/* Declare parameters */
	parms phi t1-t9 bL b0;
	
	/* Set priors */
	prior bL ~ normal(0, sd=10);
	prior b0 ~ normal(0, sd=10);
	prior phi ~ cauchy( 0, 10, lower=0);
	
	/* AR-1 prior on dose effects  */
	prior t1 ~ normal(0, sd=10);
	prior t2 ~ normal(2*t1, sd=1);
	prior t3 ~ normal(2*t2 - t1, sd=1);
	prior t4 ~ normal(2*t3 - t2, sd=1);
	prior t5 ~ normal(2*t4 - t3, sd=1);
	prior t6 ~ normal(2*t5 - t4, sd=1);
	prior t7 ~ normal(2*t6 - t5, sd=1);
	prior t8 ~ normal(2*t7 - t6, sd=1);
	prior t9 ~ normal(2*t8 - t7, sd=1);
	
	muA = t1*A1 + t2*A2 + t3*A3 + t4*A4 + t5*A5 + t6*A6 + t7*A7 + t8*A8 + t9*A9;

	/* Specify Likelihood */
	model Y ~ normal( b0 + muA + bL*L  , sd=phi);

	/* compute points on dose-effect curve, Psi(k) */
	Psi[1] = t1 ;
	Psi[2] = t2 - t1;
	Psi[3] = t3 - t2;
	Psi[4] = t4 - t3;
	Psi[5] = t5 - t4;
	Psi[6] = t6 - t5;
	Psi[7] = t7 - t6;
	Psi[8] = t8 - t7;
	Psi[9] = t9 - t8;
	
	/* output posterior summaries */
	ods output PostSummaries=psi_sum;
run;


/*Add dose #, k, for plotting and plot */
data psi_sum; set psi_sum; dose = _n_; run;

/* plot posterior mean curve along with 95% interval */
proc sgplot data=psi_sum;
	band x=dose upper=p2_5 lower=P97_5;
	series X=dose Y=mean;
	yaxis grid values=(-4 to 4 by 1) label='Dose Effect - Psi(k)';
	xaxis label="Dose Level";
run;
