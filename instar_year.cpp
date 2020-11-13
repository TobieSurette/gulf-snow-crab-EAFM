#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data:
   DATA_VECTOR(x);     // Size measurements (n).
   DATA_VECTOR(f);     // Frequency observations (n).
   DATA_IVECTOR(year); // Year identifier (n). 

   // Instar means:
   PARAMETER(mu_instar_0);                // Mean size of the first instar.
   PARAMETER_VECTOR(log_increment);       // Vector of log-scale growth increments.
   PARAMETER(log_mu_increment);           // Log-scale mean parameter associated with instar growth increments.
   PARAMETER(log_sigma_increment);        // Log-scale error parameter associated with instar growth increments.
   PARAMETER_VECTOR(mu_year);             // Annual deviations.
   PARAMETER(log_sigma_mu_year);          // Log-scale error for annual instar mean deviations.
      
   // Instar errors:
   PARAMETER(mu_log_sigma_instar);        // Instar errors log-scale mean. 
   PARAMETER(log_sigma_log_sigma_instar); // Instar errors log-scale error.
   PARAMETER_VECTOR(log_sigma_instar);    // Log-scale instar error parameters.
   
   // Instar proportions:
   PARAMETER(mu_logit_p);                 // Instar proportions log-scale mean parameter.
   PARAMETER(log_sigma_logit_p);          // Instar proportions log-scale error parameter.
   PARAMETER_MATRIX(logit_p_instar_year); // Multi-logit-scale parameters for instar proportions by year.

   // Vector sizes:      
   int n = x.size();
   int n_instar = log_sigma_instar.size();
   int n_year = mu_year.size();
   
   // Initialize log-likelihood:
   Type v = 0;
   
   // Instar mean sizes:
   vector<Type> mu_instar(n_instar);
   mu_instar[0] = mu_instar_0;
   Type mu_increment = exp(log_mu_increment);
   Type sigma_increment = exp(log_sigma_increment);
   for (int j = 1; j < n_instar; j++){
      mu_instar[j] = mu_instar[j-1] + exp(log_increment[j-1]); // Instar means.
      v += -dnorm(mu_instar[j] - mu_instar[j-1], mu_increment, sigma_increment, true); // Random effect on instar increments.
   }
   // Year effect:
   Type sigma_mu_year = exp(log_sigma_mu_year);
   v += -sum(dnorm(mu_year, 0, sigma_mu_year, true)); // Instar mean year effect.
   // Instar mean matrix:
   matrix<Type> mu_instar_year(n_instar,n_year);
   for (int i = 0; i < n_year; i++){
      for (int j = 0; j < n_instar; j++){
         mu_instar_year(j,i) =  mu_instar[j] + mu_year[i];
      }
   }   

   // Instar standard errors:
   v += -sum(dnorm(log_sigma_instar, mu_log_sigma_instar, exp(log_sigma_log_sigma_instar), true));
   vector<Type> sigma_instar = exp(log_sigma_instar);
   
   // Instar proportions:
   Type sigma_logit_p = exp(log_sigma_logit_p);
   vector<Type> sum_logit_p(n_year);
   for (int i = 0; i < n_year; i++){
      sum_logit_p[i] = 0;
      for (int j = 1; j < n_instar; j++){
         v += -dnorm(logit_p_instar_year(j-1,i), mu_logit_p, sigma_logit_p, true); // Proportions random effect.
         sum_logit_p[i] += exp(logit_p_instar_year(j-1,i));
      }
   }
   matrix<Type> p_instar_year(n_instar,n_year);
   for (int i = 0; i < n_year; i++){
      p_instar_year(0,i) = 1 / (1 + sum_logit_p[i]);
      for (int j = 1; j < n_instar; j++){
         p_instar_year(j,i) = exp(logit_p_instar_year(j-1,i)) / (1 + sum_logit_p[i]);
      }
   }
   
   // Mixture likelihood:
   for (int i = 0; i < n; i++){ 
      Type d = 0;
      for (int j = 0; j < n_instar; j++){
         d += p_instar_year(j,year[i]) * dnorm(x[i], mu_instar_year(j,year[i]), sigma_instar[j], false); 
      }
      v += -f[i] * log(d);
   }
   
   // Export parameters:
   //REPORT(mu_instar);
   //REPORT(mu_year);
   //REPORT(mu_instar_year);
   //REPORT(sigma_instar);
   //REPORT(p_instar_year);
   
   return v;
}
