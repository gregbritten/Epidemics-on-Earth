data {
  real<lower=0> POP;          	
  int<lower=0>  N_obs;         
  int<lower=0>  y[N_obs]; 
  //real          x0[4];
}
parameters {
  real<lower=0.0> beta[N_obs];
  //real<lower=0.0> gamma;
  //real<lower=0.0> sigma;
  //real<lower=0.0,upper=1.0> I0;
  real<lower=0.0,upper=1.0> E0;
  real<lower=1E-15,upper=0.1> sigma_beta;
}
transformed parameters {
  real x[N_obs,4];
  //real R0[N_obs];
  x[1,1] = 0.0;
  x[1,2] = E0;
  x[1,3] = 1-E0;
  x[1,4] = 0.0;
  for(j in 2:N_obs){
	  x[j,1] = fmax(x[j-1,1] - beta[j]*x[j-1,1]*x[j-1,3],1E-15);
	  x[j,2] = fmax(x[j-1,2] + beta[j]*x[j-1,1]*x[j-1,3] - (1.0/7.0)*x[j-1,2],1E-15);
	  x[j,3] = fmax(x[j-1,3] + (1.0/7.0)*x[j-1,2]  - (1.0/7.0)*x[j-1,3],1E-15);
	  x[j,4] = fmax(x[j-1,4] + (1.0/7.0)*x[j-1,3],1E-15);
	  // x[j,1] = fmax(x[j-1,1] - beta[j]*x[j-1,1]*x[j-1,3],1E-15);
	  // x[j,2] = fmax(x[j-1,2] + beta[j]*x[j-1,1]*x[j-1,3] - (1.0/5.0)*x[j-1,2],1E-15);
	  // x[j,3] = fmax(x[j-1,3] + (1.0/5.0)*x[j-1,2]  - (1.0/7.0)*x[j-1,3],1E-15);
	  // x[j,4] = fmax(x[j-1,4] + (1.0/7.0)*x[j-1,3],1E-15);
  }  
//  for(i in 1:N_obs){
//	  R0[i] = beta[i]/gamma;
//  }
}
model {
  //gamma ~ normal(1.0/7.0,0.005);
  //sigma ~ normal(1.0/5.0,0.005);
  beta[1] ~ uniform(0.0,2.0);
  for(i in 2:N_obs){
	  beta[i] ~ normal(0.99*beta[i-1],sigma_beta);
  }
  for(j in 2:N_obs){
	y[j] ~ poisson(fmax((x[j,3]-x[j-1,3])*POP,1E-15));
  }
}
//generated quantities {
  // real x_pred[N_obs,4]; 
  // real y_pred[N_obs,4]; 
  // x_pred = integrate_ode_rk45(SIER, x0, 0, t_obs, theta, rep_array(0.0,0), rep_array(0,0),1e-8,1e-8,1e6) ;
  // for(i in 1:4){
	  // for(j in 1:N_obs){
		  // y_pred[j,i] = poisson_rng(fmax(x[j,i]*POP,0.0));
	  // }
  // }
// }


