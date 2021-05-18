#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data:
   DATA_VECTOR(x);           // Log-size measurements (n).
   DATA_VECTOR(f);           // Frequency observations (n).
   DATA_VECTOR(maturity);    // Crab maturity(n).
   DATA_VECTOR(precision);   // Precision of size measurements (n).
   DATA_IVECTOR(group);      // group identifier (n). 
   DATA_INTEGER(n_instar);   // Number of instars.
   DATA_INTEGER(n_group);    // Number of groups.

   // Instar mean parameters:
   PARAMETER(mu_instar_0);                        // Mean size of the first instar.
   PARAMETER(log_increment);                      // Vector of log-scale growth increments.
   PARAMETER(log_increment_delta);                // Vector of log-scale growth increments modifiers.
   PARAMETER_VECTOR(mu_instar_group);             // Instar-group means random effect.
   PARAMETER(log_sigma_mu_instar_group);          // Log-scale error for instar-group means random effect.
   
   // Instar error parameters:
   PARAMETER(log_sigma);                          // Log-scale instar standard error by instar.
   PARAMETER_VECTOR(log_sigma_instar_group);
   PARAMETER(log_sigma_sigma_instar_group);
   
   // Instar proportion parameters:
   PARAMETER_VECTOR(mu_logit_p);                  // Instar proportions log-scale global mean parameter (n_instar-1).
   PARAMETER(log_sigma_logit_p_instar_group);     // Instar proportions log-scale global error parameter.
   PARAMETER_VECTOR(logit_p_instar_group);        // Proportion deviations by instar and group.

   // Mature crab parameters:
   PARAMETER_VECTOR(mu_instar_group_mature);      // Instar-group means random effect for mature crab.
   PARAMETER_VECTOR(mu_logit_p_mature);           // Instar proportions log-scale global mean parameter (n_instar-1).
   PARAMETER_VECTOR(logit_p_instar_group_mature); // Proportion deviations by instar and group.

   // Vector sizes:      
   int n = x.size();
   
   // Initialize log-likelihood:
   Type v = 0;
   
   // Global instar mean sizes:
   vector<Type> mu_instar(n_instar);
   mu_instar[0] = mu_instar_0;
   for (int j = 1; j < n_instar; j++){
      mu_instar[j] = mu_instar_0 + j * exp(log_increment) - 0.5 * (j-1) * j * exp(log_increment_delta); 
      //mu_instar[j] = mu_instar[j-1] + exp(log_increment) - j * exp(log_increment_delta); // Instar means.
      //DATA_IVECTOR(ado);
      //mu_instar[j] = mu_instar[j-1] + exp(log_increment) - j * exp(log_increment_delta + ado[j] * log_increment_delta_ado); // Instar means.
   }
  
   // Group-level instar sizes and standard errors:
   v -= sum(dnorm(mu_instar_group,        0, exp(log_sigma_mu_instar_group), true)); 
   v -= sum(dnorm(mu_instar_group_mature, 0, exp(log_sigma_mu_instar_group), true)); 
   v -= sum(dnorm(log_sigma_instar_group, 0, exp(log_sigma_sigma_instar_group), true)); 
   matrix<Type> mu(n_instar,n_group);
   matrix<Type> mu_mature(n_instar,n_group);
   matrix<Type> sigma(n_instar,n_group);
   for (int i = 0; i < n_group; i++){
      for (int j = 0; j < n_instar; j++){
         mu(j,i)        = mu_instar[j] + mu_instar_group[j * n_group + i];          // Immature mean sizes.
         mu_mature(j,i) = mu_instar[j] + mu_instar_group_mature[j * n_group + i];   // Mature mean sizes.
         sigma(j,i)     = exp(log_sigma + log_sigma_instar_group[j * n_group + i]); // Instar standard errors.
      }
   }   

   // Instar proportions:
   v -= sum(dnorm(logit_p_instar_group,        0, exp(log_sigma_logit_p_instar_group), true));
   v -= sum(dnorm(logit_p_instar_group_mature, 0, exp(log_sigma_logit_p_instar_group), true));
   matrix<Type> p(n_instar,n_group);
   matrix<Type> p_mature(n_instar,n_group);
   vector<Type> sum_logit_p(n_group);
   vector<Type> sum_logit_p_mature(n_group);
   matrix<Type> logit_p(n_instar-1,n_group);
   matrix<Type> logit_p_mature(n_instar-1,n_group);
   for (int i = 0; i < n_group; i++){
      sum_logit_p[i] = 0;
      sum_logit_p_mature[i] = 0;
      for (int j = 1; j < n_instar; j++){
         // Normalizing constants:
         sum_logit_p[i]        += exp(logit_p(j-1,i)); 
         sum_logit_p_mature[i] += exp(logit_p_mature(j-1,i)); 
         
         // Logit-scale effects:
         logit_p(j-1,i)        = mu_logit_p[j-1]        + logit_p_instar_group[(j-1) * n_group + i]; 
         logit_p_mature(j-1,i) = mu_logit_p_mature[j-1] + logit_p_instar_group_mature[(j-1) * n_group + i]; 
      }
   }
   for (int i = 0; i < n_group; i++){
      p(0,i)        = 1 / (1 + sum_logit_p[i]);
      p_mature(0,i) = 1 / (1 + sum_logit_p_mature[i]);
      for (int j = 1; j < n_instar; j++){
         p(j,i)        = exp(logit_p(j-1,i))        / (1 + sum_logit_p[i]); 
         p_mature(j,i) = exp(logit_p_mature(j-1,i)) / (1 + sum_logit_p_mature[i]);
      }
   }
   
   // Mixture likelihood:
   for (int i = 0; i < n; i++){ 
      Type d = 0;
      Type xlower = log(exp(x[i]) - precision[i] / 2);
      Type xupper = log(exp(x[i]) + precision[i] / 2);
      for (int j = 0; j < n_instar; j++){
         Type pp = (1-maturity[i]) * p(j,group[i])  + maturity[i] * p_mature(j,group[i]);
         Type mm = (1-maturity[i]) * mu(j,group[i]) + maturity[i] * mu_mature(j,group[i]);
         d +=  pp * (pnorm(xupper, mm, sigma(j,group[i])) - pnorm(xlower, mm, sigma(j,group[i]))); 
      }
      v -= f[i] * log(d);
   }
   
   // Export variables:
   REPORT(mu_instar);
   REPORT(mu);
   REPORT(mu_mature);
   REPORT(sigma);   
   REPORT(p);
   REPORT(p_mature);

   return v;
}
