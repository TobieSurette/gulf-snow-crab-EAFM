#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data:
   DATA_VECTOR(x_imm);          // Immature size measurements (ni).
   DATA_VECTOR(f_imm);          // Immature frequency observations (ni).
   // DATA_VECTOR(x_mat);          // Mature size measurements (nm).
   // DATA_VECTOR(f_mat);          // Mature frequency observations (nm).
   
   // Instar growth parameters:                       
   PARAMETER(mu0);                          // First instar mean size.
   PARAMETER(log_sigma0);                   // Log-scale standard error for first instar.
   PARAMETER_VECTOR(log_hiatt_slope);       // Hiatt slope parameters.
   PARAMETER_VECTOR(log_hiatt_intercept);   // Hiatt intercept parameters.
   PARAMETER_VECTOR(log_growth_error);      // Growth increment error inflation parameters.
   
   // Density parameters:
   PARAMETER(log_lambda_alpha);             // Log-scale global mean density.
   PARAMETER_VECTOR(log_lambda_instar);     // Log-scale instar mean density.
   PARAMETER(log_sigma_lambda_instar);      // Log-scale instar mean density error parameter.

   // Vector sizes:      
   int n_imm = x_imm.size();                // Number of immature observations.
   // int n_mat = x_mat.size();                // Number of mature observations.
   int n_instar = log_lambda_instar.size(); // Number instars in model.
   
   // Calculate instar means and errors:
   vector<Type> mu(n_instar);
   vector<Type> log_sigma(n_instar);
   mu[0] = mu0;
   log_sigma[0] = log_sigma0; 
   for (int k = 1; k < 5; k++){
      mu[k] = exp(log_hiatt_intercept[0]) + mu[k-1] + exp(log_hiatt_slope[0]) * mu[k-1];
      log_sigma[k] = log(1 + exp(log_hiatt_slope[0]) + exp(log_growth_error[0])) + log_sigma[k-1];
   }
   for (int k = 5; k < (n_instar); k++){
      mu[k] = exp(log_hiatt_intercept[1]) + mu[k-1] + exp(log_hiatt_slope[1]) * mu[k-1];
      log_sigma[k] = log(1 + exp(log_hiatt_slope[1]) + exp(log_growth_error[1])) + log_sigma[k-1];
   }
   vector<Type> sigma = exp(log_sigma);
   
   // Initialize log-likelihood:
   Type v = 0;
   
   // Random effects:
   v -= sum(dnorm(log_lambda_instar, 0, exp(log_sigma_lambda_instar), true));  // Instar-level random effect.
   
   // Calculate mean instar abundances:
   vector<Type> log_lambda = log_lambda_alpha + log_lambda_instar;

   // Likelihood evaluation for immatures:
   for (int i = 0; i < n_imm; i++){ 
      Type eta = 0;
      for (int k = 0; k < n_instar; k++){
         eta += exp(log_lambda[k]) * (pnorm(x_imm[i]+0.05, mu[k], sigma[k]) - pnorm(x_imm[i]-0.05, mu[k], sigma[k]));
      }
      v -= dpois(f_imm[i], eta, true);
   }
   
   REPORT(mu);
   REPORT(sigma);
   REPORT(log_lambda);
   REPORT(log_lambda_instar);
   
   return v;
}

