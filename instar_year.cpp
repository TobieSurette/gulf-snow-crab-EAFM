#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data:
   DATA_VECTOR(x_imm);          // Immature size measurements (n_imm).
   DATA_VECTOR(f_imm);          // Immature frequency observations (n_imm).
   DATA_IVECTOR(year_imm);      // Immature survey year (n_imm).
   DATA_VECTOR(x_mat);          // Mature size measurements (n_mat).
   DATA_VECTOR(f_mat);          // Mature frequency observations (n_mat).
   DATA_IVECTOR(year_mat);      // Mature survey year (n_mat).
   DATA_SCALAR(delta_x);        // Size-measurement bin width.   
  
   // Instar growth parameters:                       
   PARAMETER(mu0);                           // First instar mean size.
   PARAMETER(log_sigma0);                    // Log-scale standard error for first instar.
   PARAMETER_VECTOR(log_hiatt_slope);        // Hiatt slope parameters.
   PARAMETER_VECTOR(log_hiatt_intercept);    // Hiatt intercept parameters.
   PARAMETER_VECTOR(log_growth_error);       // Growth increment error inflation parameters.
   PARAMETER_VECTOR(log_mu_year);            // Log-scale instar mean year interaction (n_instar x n_year).
   PARAMETER(log_sigma_mu_year);             // Instar mean year interaction error term.

   // Abundance parameters:
   PARAMETER_VECTOR(log_n_imm_year_0);       // First year immature instar abundances (n_instar-1).
   PARAMETER_VECTOR(log_n_imm_instar_0);     // First instar recruitment for all years (n_year).
   PARAMETER(log_sigma_n_imm_instar_0);      // Log-scale first instar annual recruitment error parameter.
   PARAMETER_VECTOR(log_n_skp_instar_0);     // First year skip abundances (n_instar - 5).
   PARAMETER_VECTOR(log_n_rec_instar_0);     // First year mature recruit abundances (n_instar - 5).
   PARAMETER_VECTOR(log_n_res_instar_0);     // First year mature residual abundances (n_instar - 5).
   
   // Selectivity parameters:
   PARAMETER_VECTOR(selectivity_x50);        // Size-at-50% trawl selectivity.
   PARAMETER_VECTOR(log_selectivity_slope);  // Log-scale trawl selectivity slope.
   PARAMETER(logit_selectivity_proportion);  
       
   // Moulting probability parameters:
   PARAMETER_VECTOR(logit_p_skp);            // Logit-scale skip-moulting probabilities (n_instar).
   PARAMETER_VECTOR(logit_p_mat);            // Logit-scale moult-to-maturity probabilities (n_instar).
   PARAMETER_VECTOR(logit_p_mat_year);       // Logit-scale mout-to-maturity instar x year interaction (n_instar x n_year).
   PARAMETER(log_sigma_p_mat_year);          // Moult-to-maturity instar x year interaction error term.
                   
   // Mortality parameters:
   PARAMETER(logit_M_imm);                   // Logit-scale immature mortality.
   PARAMETER(logit_M_mat);                   // Logit-scale mature mortality.     
   
   // Initialize log-likelihood:
   Type v = 0;
   
   // Vector sizes:      
   int ni = x_imm.size();                        // Number of immature observations.
   int nm = x_mat.size();                        // Number of mature observations.
   int n_instar = log_n_imm_year_0.size() + 1;   // Number instars.
   int n_year   = log_mu_year.size() / n_instar; // Number of years.
  
   // Instar global mean and error:
   vector<Type> mu(n_instar);
   vector<Type> log_sigma(n_instar);
   mu[0] = mu0;
   log_sigma[0] = log_sigma0; 
   for (int k = 1; k < 5; k++){
      mu[k] = exp(log_hiatt_intercept[0]) + mu[k-1] + exp(log_hiatt_slope[0]) * mu[k-1];
      log_sigma[k] = log(1 + exp(log_hiatt_slope[0]) + exp(log_growth_error[0])) + log_sigma[k-1];
   }
   for (int k = 5; k < n_instar; k++){
      mu[k] = exp(log_hiatt_intercept[1]) + mu[k-1] + exp(log_hiatt_slope[1]) * mu[k-1];
      log_sigma[k] = log(1 + exp(log_hiatt_slope[1]) + exp(log_growth_error[1])) + log_sigma[k-1];
   }
   vector<Type> sigma = exp(log_sigma);
  
   // Annual instar sizes:
   v -= sum(dnorm(log_mu_year, 0, exp(log_sigma_mu_year), true));
   matrix<Type> mu_imm(n_instar,n_year);
   matrix<Type> mu_mat(n_instar,n_year);
   for (int k = 0; k < n_instar; k++){
      for (int y = 0; y < n_year; y++){
         mu_imm(k,y) = exp(log(mu[k]) + log_mu_year[y * n_instar + k]);
         mu_mat(k,y) = mu_imm(k,y); 
      }
   }
  
   // Population abundance variables:
   matrix<Type> n_imm(n_instar,n_year); // Immatures.
   matrix<Type> n_skp(n_instar,n_year); // Immature skip-moulters.
   matrix<Type> n_rec(n_instar,n_year); // Mature recruits.
   matrix<Type> n_res(n_instar,n_year); // Mature residuals (old-shelled).
   matrix<Type> n_mat(n_instar,n_year); // Total matures.
   
   // Immature abundances for first year and first instar: 
   v -= sum(dnorm(log_n_imm_instar_0, 0, exp(log_sigma_n_imm_instar_0), true)); 
   for (int k = 1; k < n_instar; k++){ 
      n_imm(k,0) = exp(log_n_imm_year_0[k-1]);   // First-year immatures.
   }
   for (int y = 0; y < n_year; y++){ 
     n_imm(0,y) = exp(log_n_imm_instar_0[y]);  // First instar recruits.
   } 
   
   // Skip abundances for first year and first instar: 
   for (int k = 0; k < 5; k++){
      n_skp(k,0) = 0; 
   } 
   for (int k = 5; k < n_instar; k++){
      n_skp(k,0) = exp(log_n_skp_instar_0[k-5]);  
   }
   for (int y = 1; y < n_year; y++){ 
      n_skp(0,y) = 0;  
   }
   
   // Mature recruits for first year and first instar: 
   for (int k = 0; k < 5; k++){   
      n_rec(k,0) = 0; 
   } 
   for (int k = 5; k < n_instar; k++){
      n_rec(k,0) = exp(log_n_rec_instar_0[k-5]);  
   }
   for (int y = 1; y < n_year; y++){ 
      n_rec(0,y) = 0;  
   }
   
   // Mature residuals for first year and first instar:
   for (int k = 0; k < 5; k++){   
      n_res(k,0) = 0; 
   } 
   for (int k = 5; k < n_instar; k++){
      n_res(k,0) = exp(log_n_res_instar_0[k-5]);  
   }
   for (int y = 1; y < n_year; y++){ 
      n_res(0,y) = 0;  
   }   
   
   // Moulting probabilities:
   vector<Type> p_skp = Type(1) / (Type(1) + exp(-logit_p_skp)); // Skip-moulting probabilities.
   v -= sum(dnorm(logit_p_mat_year, 0, exp(log_sigma_p_mat_year), true)); 
   matrix<Type> p_mat(n_instar,n_year);
   for (int k = 0; k < n_instar; k++){
      for (int y = 0; y < n_year; y++){
         p_mat(k,y) = Type(1) / (Type(1) + exp(-logit_p_mat[k] - logit_p_mat_year[y * n_instar + k]));
      }
   }
   
   // Mortality probabilities:
   Type M_imm = Type(1) / (Type(1) + exp(-logit_M_imm));         // Immature annual mortality.
   Type M_mat = Type(1) / (Type(1) + exp(-logit_M_mat));         // Mature annual mortality.
  
   // Population dynamics equations:
   for (int k = 1; k < n_instar; k++){
      for (int y = 1; y < n_year; y++){
         n_imm(k,y) = (1-p_mat(k-1,y-1)) * (1-p_skp[k-1]) * (1-M_imm) * n_imm(k-1,y-1); 
         n_skp(k,y) = (1-p_mat(k-1,y-1)) * p_skp[k-1] * (1-M_imm) * n_imm(k,y-1);    
         n_rec(k,y) = (1-M_mat) * ((1-p_skp[k-1]) * p_mat(k-1,y-1) * n_imm(k-1,y-1) + n_skp(k-1,y-1)); 
         n_res(k,y) = (1-M_mat) * (n_rec(k,y-1) + n_res(k,y-1));    
      }
   }
   for (int k = 0; k < n_instar; k++){
      for (int y = 0; y < n_year; y++){
         n_mat(k,y) = n_rec(k,y) + n_res(k,y); 
      }
   }
      
   // Likelihood evaluation for immatures:
   vector<Type> eta_imm(ni);
   eta_imm.fill(0);
   for (int i = 0; i < ni; i++){ 
      // Define selectivity curve:
      Type p0 = Type(1) / (Type(1) + exp(-exp(log_selectivity_slope[0]) * (x_imm[i] - selectivity_x50[0]))); 
      Type p1 = Type(1) / (Type(1) + exp(-exp(log_selectivity_slope[1]) * (x_imm[i] - selectivity_x50[1]))); 
      Type w = Type(1) / (Type(1) + exp(-logit_selectivity_proportion));
      Type selectivity = w * p0 + (1-w) * p1;
      for (int k = 0; k < n_instar; k++){
         eta_imm[i] += selectivity * (n_imm(k,year_imm[i]) + n_skp(k,year_imm[i])) * 
                       (pnorm(x_imm[i] + delta_x / 2, mu_imm(k,year_imm[i]), sigma[k]) - 
                        pnorm(x_imm[i] - delta_x / 2, mu_imm(k,year_imm[i]), sigma[k]));
      }
      if (eta_imm[i] > 0){ 
         v -= dpois(f_imm[i], eta_imm[i], true);
      }
   }
   
   // Likelihood evaluation for matures:
   vector<Type> eta_mat(nm);
   eta_mat.fill(0);
   for (int i = 0; i < nm; i++){ 
      // Define selectivity curve:
      Type p0 = Type(1) / (Type(1) + exp(-exp(log_selectivity_slope[0]) * (x_mat[i] - selectivity_x50[0]))); 
      Type p1 = Type(1) / (Type(1) + exp(-exp(log_selectivity_slope[1]) * (x_mat[i] - selectivity_x50[1]))); 
      Type w = Type(1) / (Type(1) + exp(-logit_selectivity_proportion));
      Type selectivity = w * p0 + (1-w) * p1;
      for (int k = 0; k < n_instar; k++){
         eta_mat[i] += selectivity * 
                       n_mat(k,year_mat[i]) * 
                      (pnorm(x_mat[i] + delta_x / 2, mu_mat(k,year_mat[i]), sigma[k]) - 
                       pnorm(x_mat[i] - delta_x / 2, mu_mat(k,year_mat[i]), sigma[k]));
      }
      if (eta_mat[i] > 0){ 
         v -= dpois(f_mat[i], eta_mat[i], true);
      }
   }
   
   // Export results:
   REPORT(mu_imm);
   REPORT(sigma);
   REPORT(p_skp);
   REPORT(p_mat);
   REPORT(n_imm);
   REPORT(n_skp);
   REPORT(n_rec);
   REPORT(n_res);
   REPORT(n_mat);
   REPORT(eta_imm);
   REPORT(eta_mat);
   REPORT(M_imm);
   REPORT(M_mat);
   
   return v;
}
