#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data:
   DATA_VECTOR(x);
   DATA_VECTOR(w);
   
   // Parameters:
   PARAMETER(log_mu_0);  
   PARAMETER(log_sigma_mu);
   PARAMETER(log_sigma_sigma);
   PARAMETER_VECTOR(log_sigma);
   PARAMETER(log_mu_inc);
   PARAMETER(log_sigma_inc);
   PARAMETER_VECTOR(log_inc);
   PARAMETER_VECTOR(log_p);
   
   // Vector sizes:      
   int n = x.size();
   int n_instar = log_p.size() + 1;
   
   // Initialize log-likelihood:
   Type v = 0;
   
   // Random effect for increments:
   Type sigma_inc = exp(log_sigma_inc);
   v += -sum(dnorm(log_inc, log_mu_inc, sigma_inc, true));
   
   // Mixture component proportions:
   vector<Type> p(n_instar);
   p[0] = 1 / (1 + sum(exp(log_p)));
   for (int j = 1; j < n_instar; j++){
      p[j] = exp(log_p[j-1]) / (1 + sum(exp(log_p)));
   }

   // Mixture component means:
   vector<Type> mu(n_instar);
   mu[0] = exp(log_mu_0);
   for (int j = 1; j < n_instar; j++){
      mu[j] = mu[j-1] + exp(log_inc[j-1]);
   }
   
   // Mixture component standard errors:
   Type sigma_sigma = exp(log_sigma_sigma);
   v += -sum(dnorm(log_sigma, log_sigma_mu, sigma_sigma, true));
   vector<Type> sigma = exp(log_sigma);
   
   // Mixture likelihood:
   for (int i = 0; i < n; i++){ 
      Type d = 0;
      for (int j = 0; j < n_instar; j++){
         d += p[j] * dnorm(x[i], mu[j], sigma[j], false); 
      }
      v += -w[i] * log(d);
   }
   
   // Export parameters:
   REPORT(p);
   REPORT(mu);
   REPORT(sigma);
   
   return v;
}
