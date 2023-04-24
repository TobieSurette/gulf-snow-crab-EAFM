#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data:
   DATA_VECTOR(x);                                // Log-size measurements (n).
   DATA_VECTOR(f);                                // Frequency observations (n).
   DATA_VECTOR(maturity);                         // Crab maturity (n).
   DATA_IVECTOR(group);                           // Group identifier (n). 
   
   // Instar growth parameters:
   PARAMETER(mu_0);                               // Mean size of the first instar.
   PARAMETER(a);                                  // Immature growth slope parameter.
   PARAMETER(b);                                  // Immature growth intercept parameter.
   PARAMETER(a_pubescent);                        // Pubescent growth slope parameter.
   PARAMETER(b_pubescent);                        // Pubescent growth intercept parameter.
   PARAMETER(a_mature);                           // Mature growth slope parameter.
   PARAMETER(b_mature);                           // Mature growth intercept parameter.
  
   // Random effect error parameters:
   PARAMETER(log_sigma_mu);                       // Instar mean random effect error.
   PARAMETER(log_sigma_p);                        // Instar proportion random effect error.
  
   // Instar mean group-level effect:
   PARAMETER_VECTOR(delta_mu_immature_group);     // Immature instar mean size deviations.
   PARAMETER_VECTOR(delta_mu_pubescent_group);    // Pubescent instar mean size deviations.
   PARAMETER_VECTOR(delta_mu_mature_group);       // Mature instar mean size deviations.
   
   // Instar proportion effects:
   PARAMETER_VECTOR(logit_p_immature);            // Immature instar global proportions.
   PARAMETER_VECTOR(logit_p_pubescent);           // Pubescent instar global proportions.
   PARAMETER_VECTOR(logit_p_mature);              // Mature instar global proportions.
   PARAMETER_VECTOR(delta_p_immature_group);      // Immature instar group-level proportion deviations.
   PARAMETER_VECTOR(delta_p_pubescent_group);     // Pubescent instar group-level proportion deviations.
   PARAMETER_VECTOR(delta_p_mature_group);        // Mature instar group-level proportion deviations.
   
   // Vector sizes:      
   int n = x.size();
   int n_instar = logit_p_immature.size() + 1;
   int n_group = delta_p_immature_group.size() / (n_instar-1);
    
   // Initialize log-likelihood:
   Type v = 0;
  
   // Global instar mean sizes:
   vector<Type> mu_immature(n_instar);
   vector<Type> mu_pubescent(n_instar);
   vector<Type> mu_mature(n_instar);
   mu_immature[0] = mu_0;
   for (int j = 1; j < n_instar; j++){
      mu_immature[j]  = a * mu_immature[j-1] + b;                     // Immature growth.
      mu_pubescent[j] = a_pubescent * mu_immature[j-1] + a_pubescent; // Pubescent growth.
      mu_mature[j]    = a_mature * mu_pubescent[j-1] + a_mature;      // Mature growth.
   }
   
   // Group-level instar mean sizes:
   matrix<type> mu_immature_group(n_instar,n_group);
   matrix<type> mu_pubescent_group(n_instar,n_group);
   matrix<type> mu_mature_group(n_instar,n_group);
   v -= sum(dnorm(delta_mu_immature_group, 0, exp(log_sigma_mu), true)); 
   v -= sum(dnorm(delta_mu_pubescent_group, 0, exp(log_sigma_mu), true)); 
   v -= sum(dnorm(delta_mu_mature_group, 0, exp(log_sigma_mu), true)); 
   for (int i = 0; i < n_group; i++){
      for (int j = 0; j < n_instar; j++){
         mu_immature_group(j,i)  =  mu_immature[j]  + delta_mu_immature_group[j * n_group + i];
         mu_pubescent_group(j,i) =  mu_pubescent[j] + delta_mu_pubescent_group[j * n_group + i];
         mu_mature_group(j,i)    =  mu_mature[j]    + delta_mu_mature_group[j * n_group + i];
     }
   }

   // Group-level instar proportions effects:
   matrix<type> logit_p_immature_group(n_instar-1,n_group);
   matrix<type> logit_p_pubescent_group(n_instar-1,n_group);
   matrix<type> logit_p_mature_group(n_instar-1,n_group);
   v -= sum(dnorm(delta_p_immature_group, 0, exp(log_sigma_p), true)); 
   v -= sum(dnorm(delta_p_pubescent_group, 0, exp(log_sigma_p), true)); 
   v -= sum(dnorm(delta_p_mature_group, 0, exp(log_sigma_p), true)); 
   for (int i = 0; i < n_group; i++){
     for (int j = 0; j < (n_instar-1); j++){
       logit_p_immature_group(j,i)  =  logit_p_immature[j]  + delta_p_immature_group[j * n_group + i];
       logit_p_pubescent_group(j,i) =  logit_p_pubescent[j] + delta_p_pubescent_group[j * n_group + i];
       logit_p_mature_group(j,i)    =  logit_p_mature[j]    + delta_p_mature_group[j * n_group + i];
     }
   }
   
   // Instar proportion effects:
   matrix<Type> p_immature(n_instar,n_group);
   matrix<Type> p_pubescent(n_instar,n_group);
   matrix<Type> p_mature(n_instar,n_group);
   vector<Type> sum_logit_p_immature(n_group);   // Immature proportion accumulator.
   vector<Type> sum_logit_p_pubescent(n_group);  // Pubescent proportion accumulator.
   vector<Type> sum_logit_p_mature(n_group);     // Mature proportion accumulator.
   for (int i = 0; i < n_group; i++){
      sum_logit_p_immature[i] = 0;
      sum_logit_p_pubescent[i] = 0;
      sum_logit_p_mature[i] = 0;
      for (int j = 1; j < n_instar; j++){
         sum_logit_p_immature[i]  += exp(logit_p_immature_group(j-1,i));  // Immature. 
         sum_logit_p_pubescent[i] += exp(logit_p_pubescent_group(j-1,i)); // Pubescent.
         sum_logit_p_mature[i]    += exp(logit_p_mature_group(j-1,i));    // Mature.
      }
   }
   
   // Calculate instar proportions by group level:
   for (int i = 0; i < n_group; i++){
      p_immature(0,i)  = 1 / (1 + sum_logit_p_immature[i]);
      p_pubescent(0,i) = 1 / (1 + sum_logit_p_pubescent[i]);
      p_mature(0,i)    = 1 / (1 + sum_logit_p_mature[i]);
      for (int j = 1; j < n_instar; j++){
         p_immature(j,i)  = exp(logit_p_immature_group(j-1,i))  / (1 + sum_logit_p_immature[i]);   // Immature group proportions.
         p_pubescent(j,i) = exp(logit_p_pubescent_group(j-1,i)) / (1 + sum_logit_p_pubescent[i]);  // Pubescent group proportions.
         p_mature(j,i)    = exp(logit_p_mature_group(j-1,i))    / (1 + sum_logit_p_mature[i]);     // Mature group proportions.
      }
   }
  
   // Mixture likelihood:
   for (int i = 0; i < n; i++){ 
      Type d = 0;
      for (int j = 0; j < n_instar; j++){
         if (maturity[i] == 1) d += p_immature(j,group[i])  * dnorm(x[i], mu_immature_group(j,group[i]),  exp(log_sigma))
         if (maturity[i] == 2) d += p_pubescent(j,group[i]) * dnorm(x[i], mu_pubescent_group(j,group[i]), exp(log_sigma))
         if (maturity[i] == 3) d += p_mature(j,group[i])    * dnorm(x[i], mu_mature_group(j,group[i]),    exp(log_sigma))
      }
      v -= f[i] * log(d);
   }
  
   return v;
}