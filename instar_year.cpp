#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data:
   DATA_VECTOR(x);
   DATA_VECTOR(n);
   DATA_VECTOR(year); 
   
   // Parameters:
   PARAMETER(log_mu_instar_0);  
   PARAMETER(log_sigma_mu);
   PARAMETER(log_sigma_sigma);
   PARAMETER_VECTOR(log_sigma);
   PARAMETER(log_mu_inc);
   PARAMETER(log_sigma_inc);
   PARAMETER_VECTOR(log_inc);
   PARAMETER_VECTOR(log_p);
   PARAMETER_MATRIX(mu_year);
      
   // Vector sizes:      
   int nx = x.size();
   int k = log_p.size() + 1;
   
   // Initialize log-likelihood:
   Type v = 0;
   
   // Random effect for increments:
   Type sigma_instar_inc = exp(log_sigma_instar_inc);
   v += -sum(dnorm(log_instar_inc, log_mu_instar_inc, sigma_instar_inc, true));
   
   // Mixture component proportions:
   matrix<Type> p_year_instar(n_year, n_instar);
   for (int i = 0; i < n_year; i++){
      p_year_instar(i,0)  = 1 / (1 + sum(exp(log_p_year_instar.row(i))));
      for (int j = 1; j < k; j++){
         p_year_instar(i,j) = exp(log_p_year_instar(i,j-1)) / (1 + sum(exp(log_p_year_instar.row(i))));
      }
   }
   
   // Mixture component means:
   vector<Type> mu_instar(k);
   mu_instar[0] = exp(log_mu_instar_0);
   for (int j = 1; j < n_instar; j++){
      mu_instar[j] = mu_instar[j-1] + exp(log_inc[j-1]);
   }
   for (int i = 0; i < n_year; i++){
      mu_year[i] ~ dnorm(0, sigma_year) 
   }
   Type sigma_instar_year = exp(log_sigma_instar_year);
   for (int i = 0; i < n_year; i++){
      for (int j = 1; j < n_instar; j++){
         mu_year_instar(i,j) ~ dnorm(0, sigma_instar_year) // Interaction.
      }
   }
   
   
   // Define mixture components:
   for (int i = 0; i < n_year; i++){
      for (int j = 0; j < n_instar; j++){
         mu(i,j) = mu_year[i] + mu_instar[j] + mu_year_instar(i,j);
         log_sigma(i,j) = log_sigma_year[i] + log_sigma_instar[j] + log_sigma_year_instar(i,j);
      }
      for (int j = 1; j < n_instar; j++){
         logit_p(i,j-1) = logit_p_year[i] + logit_p_instar[j-1] + logit_p_year_instar(i,j-1);
      }
   }

   
   Type sigma_year = exp(log_sigma_year);
   for (int i = 0; i < n_year; i++){
      for (int j = 0; j < n_instar; j++){
         v += -dnorm(mu_year_instar(i,j), mu_instar[j], sigma_year, true)
      }
   }
   
   // Mixture component standard errors:
   Type sigma_sigma_instar = exp(log_sigma_instar);
   v += -sum(dnorm(log_sigma_instar, log_sigma_mu_instar, sigma_sigma_instar, true));
   vector<Type> sigma_instar = exp(log_sigma_instar);
   
   // Mixture likelihood:
   for (int i = 0; i < n; i++){ 
      Type d = 0;
      for (int j = 0; j < n_instar; j++){
         d += p_year_instar(year[i],j) * dnorm(x[i], mu_year_instar(year[i],j), sigma_instar[j], false); 
      }
      v += -n[i] * log(d);
   }
   
   // Export parameters:
   REPORT(p);
   REPORT(mu);
   REPORT(sigma);
   
   return v;
}
