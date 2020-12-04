#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data:
   DATA_VECTOR(x);         // Size measurements.
   DATA_VECTOR(f);         // Frequency observations.
   DATA_INTEGER(n_instar); // Number of instars.
   
   // Parameters:                       
   PARAMETER(mu0);                        // First instar mean size.
   PARAMETER(log_sigma0);                 // Log-scale standard error for first instar.
   PARAMETER_VECTOR(log_hiatt_slope);     // Hiatt slope parameters.
   PARAMETER_VECTOR(log_hiatt_intercept); // Hiatt intercept parameters.
   PARAMETER_VECTOR(log_growth_error);    // Growth increment error inflation parameters.
   PARAMETER_VECTOR(logit_p_instar);      // Multi-logit-scale parameters for instar proportions.
   
   // Vector sizes:      
   int n = x.size();
   
   // Calculate instar means and errors:
   vector<Type> mu(n_instar);
   vector<Type> log_sigma(n_instar);
   mu[0] = mu0;
   log_sigma[0] = log_sigma0; 
   for (int i = 1; i < 5; i++){
      mu[i] = exp(log_hiatt_intercept[0]) + mu[i-1] + exp(log_hiatt_slope[0]) * mu[i-1];
      log_sigma[i] = log(1 + exp(log_hiatt_slope[0]) + exp(log_growth_error[0])) + log_sigma[i-1];
   }
   for (int i = 5; i < n_instar; i++){
      mu[i] = exp(log_hiatt_intercept[1]) + + mu[i-1] + exp(log_hiatt_slope[1]) * mu[i-1];
      log_sigma[i] = log(1 + exp(log_hiatt_slope[1]) + exp(log_growth_error[1])) + log_sigma[i-1];
   }
   vector<Type> sigma = exp(log_sigma);
      
   // Initialize log-likelihood:
   Type v = 0;
   
   // Instar proportions:
   vector<Type> p(n_instar);
   p[0] = 1 / (1 + sum(exp(logit_p_instar)));
   for (int j = 1; j < n_instar; j++){
      p[j] = exp(logit_p_instar[j-1]) / (1 + sum(exp(logit_p_instar)));
   }
   
   // Mixture likelihood:
   for (int i = 0; i < n; i++){ 
      Type d = 0;
      for (int j = 0; j < n_instar; j++){
         d += p[j] * dnorm(x[i], mu[j], sigma[j], false); 
      }
      v -= f[i] * log(d);
   }
   
   // Export instar stats:
   REPORT(mu);
   REPORT(sigma);
   REPORT(p);
   
   return v;
}
