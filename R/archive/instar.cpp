#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data:
   DATA_VECTOR(x); // Size measurements.
   DATA_VECTOR(f); // Frequency observations.
   
   // Instar means:
   PARAMETER(mu_instar_0);                // Size of the first instar.
   PARAMETER_VECTOR(log_increment);       // Vector of log-scale growth increments.
   PARAMETER(log_mu_increment);           // Log-scale mean parameter associated with instar growth increments.
   PARAMETER(log_sigma_increment);        // Log-scale error parameter associated with instar growth increments.

   // Instar errors:
   PARAMETER(mu_log_sigma_instar);        // Instar errors log-scale mean. 
   PARAMETER(log_sigma_log_sigma_instar); // Instar errors log-scale error. 
   PARAMETER_VECTOR(log_sigma_instar);    // Log-scale instar error parameters.
   
   // Instar proportions:
   PARAMETER_VECTOR(logit_p_instar);      // Multi-logit-scale parameters for instar proportions.
   
   // Vector sizes:      
   int n = x.size();
   int n_instar = log_sigma_instar.size();
   
   // Initialize log-likelihood:
   Type v = 0;
   
   // Instar mean sizes and growth increments:
   vector<Type> mu_instar(n_instar);
   mu_instar[0] = mu_instar_0;
   Type mu_increment = exp(log_mu_increment);
   Type sigma_increment = exp(log_sigma_increment);
   for (int j = 1; j < n_instar; j++){
      mu_instar[j] = mu_instar[j-1] + exp(log_increment[j-1]); // Instar means.
      v += -dnorm(mu_instar[j] - mu_instar[j-1], mu_increment, sigma_increment, true); // Growth increments.
   }

   // Instar standard errors:
   v += -sum(dnorm(log_sigma_instar, mu_log_sigma_instar, exp(log_sigma_log_sigma_instar), true));
   vector<Type> sigma_instar = exp(log_sigma_instar);
   
   // Instar proportions:
   vector<Type> p_instar(n_instar);
   p_instar[0] = 1 / (1 + sum(exp(logit_p_instar)));
   for (int j = 1; j < n_instar; j++){
      p_instar[j] = exp(logit_p_instar[j-1]) / (1 + sum(exp(logit_p_instar)));
   }

   // Mixture likelihood:
   for (int i = 0; i < n; i++){ 
      Type d = 0;
      for (int j = 0; j < n_instar; j++){
         d += p_instar[j] * dnorm(x[i], mu_instar[j], sigma_instar[j], false); 
      }
      v += -f[i] * log(d);
   }
   
   // Export parameters:
   REPORT(mu_instar);
   REPORT(sigma_instar);
   REPORT(p_instar);
   
   return v;
}
