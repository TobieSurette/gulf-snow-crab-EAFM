#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data:
   DATA_VECTOR(x_imm);          // Immature size measurements (n_imm).
   DATA_VECTOR(f_imm);          // Immature frequency observations (n_imm).
   DATA_IVECTOR(year_imm);      // Immature survey year (n_imm).
   DATA_VECTOR(x_mat);          // Mature size measurements (n_mat).
   DATA_VECTOR(f_mat);          // Mature frequency observations (n_mat).
   DATA_IVECTOR(year_mat);      // Immature survey year (n_mat).
   DATA_SCALAR(dx);             // Size-bin width.   
  
   // Instar growth parameters:                       
   PARAMETER(mu0);                          // First instar mean size.
   PARAMETER(log_sigma0);                   // Log-scale standard error for first instar.
   PARAMETER_VECTOR(log_hiatt_slope);       // Hiatt slope parameters.
   PARAMETER_VECTOR(log_hiatt_intercept);   // Hiatt intercept parameters.
   PARAMETER_VECTOR(log_growth_error);      // Growth increment error inflation parameters.
   
   PARAMETER_VECTOR(log_mu_year);           // Log-scale instar mean year interaction.
   PARAMETER(log_sigma_mu_year);            // Instar mean year interaction error term.
   
   PARAMETER_VECTOR(log_n_imm_year_0);      // First year immature abundances (n_instar).
   PARAMETER_VECTOR(log_n_skip_instar_0);   // First year skip abundances (n_instar - 5).
   PARAMETER_VECTOR(log_n_imm_instar_0);    // First instar recruitment (n_year).
   PARAMETER_VECTOR(log_n_rec_instar_0);    // First year mature recruit abundances (n_instar - 5).
   PARAMETER_VECTOR(log_n_res_instar_0);    // First year mature residual abundances (n_instar - 5).
   
   PARAMETER(selectivity_x50);              // Logit-scale trawl selectivity at 50%.
   PARAMETER(logit_selectivity_slope);      // Logit-scale trawl selectivity slope.

   PARAMETER_VECTOR(logit_p_skip);          // Logit-scale skip-moulting probabilities (n_instar-1).
   PARAMETER_VECTOR(logit_p_mat);           // Logit-scale moult-to-maturity probabilities (n_instar-1).

   
   PARAMETER(logit_M_imm);                  // Logit-scale immature mortality.
   PARAMETER(logit_M_mat);                  // Logit-scale mature mortality.     
   
   // Vector sizes:      
   int ni = x_imm.size();                   // Number of immature observations.
   int nm  = x_mat.size();                  // Number of mature observations.
   int n_instar = log_n_imm_year_0.size();  // Number instars.
   int n_year = log_n_imm_instar_0.size();  // Number of years.
  
   // Initialize log-likelihood:
   Type v = 0;
  
   // Instar global mean and error:
   vector<Type> mu(n_instar+1);
   vector<Type> log_sigma(n_instar+1);
   mu[0] = mu0;
   log_sigma[0] = log_sigma0; 
   for (int k = 1; k < 5; k++){
      mu[k] = exp(log_hiatt_intercept[0]) + mu[k-1] + exp(log_hiatt_slope[0]) * mu[k-1];
      log_sigma[k] = log(1 + exp(log_hiatt_slope[0]) + exp(log_growth_error[0])) + log_sigma[k-1];
   }
   for (int k = 5; k < (n_instar+1); k++){
      mu[k] = exp(log_hiatt_intercept[1]) + mu[k-1] + exp(log_hiatt_slope[1]) * mu[k-1];
      log_sigma[k] = log(1 + exp(log_hiatt_slope[1]) + exp(log_growth_error[1])) + log_sigma[k-1];
   }
   vector<Type> sigma = exp(log_sigma);
  
   // Annual instar sizes:
   v -= sum(dnorm(log_mu_year, 0, exp(log_sigma_mu_year), true));
   matrix<Type> mu_imm(n_instar,n_year);
   matrix<Type> mu_mat(n_instar,n_year);
   for (int k = 0; k < (n_instar); k++){
      for (int y = 0; y < (n_year); y++){
         mu_imm(k,y) = exp(log(mu[k]) + log_mu_year[y * n_instar + k]);
         mu_mat(k,y) = mu_imm(k,y);
      }
   }
  
   // Population abundance variables:
   matrix<Type> n_imm(n_instar,n_year);  // Immatures.
   matrix<Type> n_skip(n_instar,n_year); // Immature skip-moulters.
   matrix<Type> n_rec(n_instar,n_year);  // Mature recruits.
   matrix<Type> n_res(n_instar,n_year);  // Mature residuals (old-shelled).
   matrix<Type> n_mat(n_instar,n_year);  // Total matures.
   
   // Initialize first year immature abundances:
   for (int k = 0; k < (n_instar); k++){ 
      n_imm(k,0) = exp(log_n_imm_year_0[k]);  
   }
   // Initialize first year skip abundances:
   for (int k = 0; k < 5; k++){
      n_skip(k,0) = 0; 
   } 
   for (int k = 5; k < n_instar; k++){
      n_skip(k,0) = exp(log_n_skip_instar_0[k-5]);  
   }
   // Initialize first instar recruitment:
   for (int y = 1; y < (n_year); y++){ 
     n_imm(0,y) = exp(log_n_imm_instar_0[y-1]); 
   } 
   // Initialize first-year mature recruits:
   for (int k = 0; k < 5; k++){   
      n_rec(k,0) = 0; 
   } 
   for (int k = 5; k < n_instar; k++){
      n_rec(k,0) = exp(log_n_rec_instar_0[k-5]);  
   }
   // Initialize first-year mature residuals:
   for (int k = 0; k < 5; k++){   
      n_res(k,0) = 0; 
   } 
   for (int k = 5; k < n_instar; k++){
      n_res(k,0) = exp(log_n_res_instar_0[k-5]);  
   }
   
   // Moulting probabilities:
   vector<Type> p_skip = Type(1) / (Type(1) + exp(-logit_p_skip)); // Skip-moulting probabilities.
   vector<Type> p_mat = Type(1) / (Type(1) + exp(-logit_p_mat));   // Moult-to-maturity probabilities.
   
   // Mortality probabilities:
   Type M_imm = Type(1) / (Type(1) + exp(-logit_M_imm));  // Immature mortality.
   Type M_mat = Type(1) / (Type(1) + exp(-logit_M_mat));  // Mature mortality.
  
   // Population dynamics equations:
   for (int k = 1; k < (n_instar); k++){
      for (int y = 1; y < (n_year); y++){
         n_imm(k,y)  = (1-p_mat[k-1]) * (1-p_skip[k-1]) * (1-M_imm) * n_imm(k-1,y-1); 
         n_skip(k,y) = (1-p_mat[k-1]) * p_skip[k-1] * (1-M_imm) * n_imm(k-1,y-1);    
         n_rec(k,y)  = (1-M_mat) * ((1-p_skip[k-1]) * p_mat[k-1] * n_imm(k-1,y-1) + n_skip(k-1,y-1)); 
         n_res(k,y)  = (1-M_mat) * n_rec(k,y-1);    
         n_mat(k,y)  = n_rec(k,y) + n_res(k,y); 
      }
   }
      
   // Likelihood evaluation for immatures:
   for (int i = 0; i < ni; i++){ 
      Type eta = 0;
      Type selectivity = Type(1) / (Type(1) + exp(-logit_selectivity_slope * (x_imm[i] - selectivity_x50)));  // Trawl fishing selectivity.
      for (int k = 0; k < n_instar; k++){
         eta += selectivity * (n_imm(k,year_imm[i]) + n_skip(k,year_imm[i])) * 
                (pnorm(x_imm[i] + dx / 2, mu_imm(k,year_imm[i]), sigma[k]) - 
                 pnorm(x_imm[i] - dx / 2, mu_imm(k,year_imm[i]), sigma[k]));
      }
      v -= dpois(f_imm[i], eta, true);
   }
   // Likelihood evaluation for matures:
   for (int i = 0; i < nm; i++){ 
      Type eta = 0;
      Type selectivity = Type(1) / (Type(1) + exp(-logit_selectivity_slope * (x_mat[i] - selectivity_x50)));  // Trawl fishing selectivity.
      for (int k = 0; k < n_instar; k++){
         eta += selectivity * 
                n_mat(k,year_mat[i]) * 
                (pnorm(x_mat[i] + dx / 2, mu_mat(k,year_mat[i]), sigma[k]) - 
                 pnorm(x_mat[i] - dx / 2, mu_mat(k,year_mat[i]), sigma[k]));
      }
      v -= dpois(f_mat[i], eta, true);
   }
   
   return v;
}
