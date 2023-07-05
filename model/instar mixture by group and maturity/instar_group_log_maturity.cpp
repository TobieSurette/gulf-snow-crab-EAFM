#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data:
   DATA_VECTOR(x);                                // Log-size measurements (n).
   DATA_VECTOR(f);                                // Frequency observations (n).
   DATA_VECTOR(maturity);                         // Crab maturity(n).
   DATA_VECTOR(precision);                        // Precision of size measurements (n).
   DATA_IVECTOR(group);                           // group identifier (n). 
   DATA_INTEGER(n_instar);                        // Number of instars.
   DATA_INTEGER(n_group);                         // Number of groups.
   DATA_INTEGER(first_instar);                    // Number corresponding to first instar.
   DATA_INTEGER(first_mature_instar);             // Number corresponding to first mature instar.

   // Instar growth parameters:
   PARAMETER(mu_instar_0);                        // Mean size of the first instar.
   PARAMETER(log_increment);                      // Vector of log-scale growth increments.
   PARAMETER(log_increment_delta);                // Vector of log-scale growth increments modifiers.
   PARAMETER(log_increment_mature);               // Instar size moult-to-maturity growth scaling. 
   PARAMETER(log_delta_maturation);               // Instar size selectivity maturation bias.
    
   // Instar size random effects:
   PARAMETER_VECTOR(mu_instar_group);             // Instar-group means random effect.
   PARAMETER(log_sigma_mu_instar_group);          // Log-scale error for instar-group means random effect.

   // Instar error parameters:
   PARAMETER(log_sigma_0);                        // Instar IV error. 
   PARAMETER(log_sigma_increment);                // Instar error incremental effect. 
   PARAMETER(log_sigma_maturation);               // Instar error effect associated with maturation. 
      
   // Instar error random effects:
   PARAMETER_VECTOR(log_sigma_instar_group);      // Instar-group error random effect error parameter.
   PARAMETER(log_sigma_sigma_instar_group);       // Instar-group error random effect.

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
   
   // Global instar size growth equations:
   vector<Type> mu_instar(n_instar);
   vector<Type> mu_instar_mature(n_instar);
   vector<Type> log_sigma_instar(n_instar);  
   vector<Type> log_sigma_instar_mature(n_instar);   
   mu_instar[0] = mu_instar_0;
   mu_instar_mature[0] = mu_instar_0;
   log_sigma_instar[0] = log_sigma_0;
   for (int j = 1; j < (first_mature_instar - first_instar); j++){
      // Instar means:
      mu_instar[j] = mu_instar[j-1] + exp(log_increment) - j * exp(log_increment_delta); // Immature instar means.
      mu_instar_mature[j] = mu_instar[j];                                                // Mature instar means.
      
      // Instar errors:
      log_sigma_instar[j] = log_sigma_0 + j * log_sigma_increment;                       // Immature instar errors.
      log_sigma_instar_mature[j] = log_sigma_instar[j];                                  // Mature instar errors.
   }
   mu_instar[first_mature_instar - first_instar] = mu_instar[first_mature_instar - first_instar-1] + 
                                                   exp(log_increment) - 
                                                   (first_mature_instar - first_instar) * exp(log_increment_delta) - 
                                                   log_delta_maturation;
   mu_instar_mature[first_mature_instar - first_instar] = mu_instar[first_mature_instar - first_instar-1] + 
                                                          exp(log_increment_mature) + 
                                                          log_delta_maturation;

   log_sigma_instar[first_mature_instar - first_instar]        = log_sigma_instar[first_mature_instar - first_instar - 1] + log_sigma_maturation;
   log_sigma_instar_mature[first_mature_instar - first_instar] = log_sigma_instar[first_mature_instar - first_instar] + log_sigma_maturation;      
   for (int j = (first_mature_instar - first_instar + 1); j < n_instar; j++){
      mu_instar[j] = mu_instar[j-1] + exp(log_increment) - j * exp(log_increment_delta);
     
      mu_instar_mature[j] = mu_instar[j-1] + exp(log_increment_mature) + log_delta_maturation;

      log_sigma_instar[j]        = log_sigma_instar[first_mature_instar - first_instar];
      log_sigma_instar_mature[j] = log_sigma_instar[first_mature_instar - first_instar];
   }
   
   // Instar size and error by group:
   v -= sum(dnorm(mu_instar_group, 0, exp(log_sigma_mu_instar_group), true)); 
   v -= sum(dnorm(mu_instar_group_mature, 0, exp(log_sigma_mu_instar_group), true)); 
   v -= sum(dnorm(log_sigma_instar_group, 0, exp(log_sigma_sigma_instar_group), true));
   matrix<Type> mu(n_instar,n_group);
   matrix<Type> mu_mature(n_instar,n_group);
   matrix<Type> sigma(n_instar,n_group);
   matrix<Type> sigma_mature(n_instar,n_group);
   for (int i = 0; i < n_group; i++){
      for (int j = 0; j < n_instar; j++){
         // Instar mean matrix:
         mu(j,i)           =  mu_instar[j] + mu_instar_group[j * n_group + i]; 
         mu_mature(j,i)    =  mu_instar_mature[j] + mu_instar_group_mature[j * n_group + i]; 
         
         // Instar error matrix:
         sigma(j,i)        = exp(log_sigma_instar[j] + log_sigma_instar_group[j * n_group + i]); 
         sigma_mature(j,i) = exp(log_sigma_instar_mature[j] + log_sigma_instar_group[j * n_group + i]); 
      }
   }  

   // Instar proportion effects:
   v += -sum(dnorm(logit_p_instar_group, 0, exp(log_sigma_logit_p_instar_group), true));
   v += -sum(dnorm(logit_p_instar_group_mature, 0, exp(log_sigma_logit_p_instar_group), true));
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
         // Immatures:
         logit_p(j-1,i) = mu_logit_p[j-1] + logit_p_instar_group[(j-1) * n_group + i]; 
         sum_logit_p[i] += exp(logit_p(j-1,i)); 
         
         // Matures:
         logit_p_mature(j-1,i) = mu_logit_p_mature[j-1] + logit_p_instar_group_mature[(j-1) * n_group + i];
         sum_logit_p_mature[i] += exp(logit_p_mature(j-1,i));         
      }
   }
   // Calculate instar proportions:
   for (int i = 0; i < n_group; i++){
      p(0,i)        = 1 / (1 + sum_logit_p[i]);
      p_mature(0,i) = 1 / (1 + sum_logit_p_mature[i]);
      for (int j = 1; j < n_instar; j++){
         p(j,i)        = exp(logit_p(j-1,i)) / (1 + sum_logit_p[i]); 
         p_mature(j,i) = exp(logit_p_mature(j-1,i)) / (1 + sum_logit_p_mature[i]); // Proportions by instar and group.
      }
   }
   
   // Mixture likelihood:
   for (int i = 0; i < n; i++){ 
      Type d = 0;
      Type xlower = log(exp(x[i]) - precision[i] / 2);
      Type xupper = log(exp(x[i]) + precision[i] / 2);
      for (int j = 0; j < n_instar; j++){
         Type pp = (1-maturity[i]) * p(j,group[i])     + maturity[i] * p_mature(j,group[i]);
         Type mm = (1-maturity[i]) * mu(j,group[i])    + maturity[i] * mu_mature(j,group[i]);
         Type ss = (1-maturity[i]) * sigma(j,group[i]) + maturity[i] * sigma_mature(j,group[i]);
         d +=  pp * (pnorm(xupper, mm, ss) - pnorm(xlower, mm, ss)); 
      }
      v -= f[i] * log(d);
   }
   
   // Export values:
   REPORT(mu);
   REPORT(mu_mature);
   REPORT(p);
   REPORT(p_mature);
   REPORT(mu_instar);
   REPORT(mu_instar_mature);
   REPORT(log_sigma_instar);
   REPORT(log_sigma_instar_mature);
   REPORT(mu_instar_group);
   REPORT(sigma);
   REPORT(sigma_mature);
   REPORT(log_sigma_instar_group);
   
   return v;
}
