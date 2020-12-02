#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
   // Data:
   DATA_VECTOR(x);          // Size measurements (n).
   DATA_VECTOR(f);          // Frequency observations (n).
   DATA_IVECTOR(tow);       // Tow identifiers (n).
   DATA_VECTOR(swept_area); // Tow swept area (n_tow).
   
   // Instar growth parameters:                       
   PARAMETER(mu0);                          // First instar mean size.
   PARAMETER(log_sigma0);                   // Log-scale standard error for first instar.
   PARAMETER_VECTOR(log_hiatt_slope);       // Hiatt slope parameters.
   PARAMETER_VECTOR(log_hiatt_intercept);   // Hiatt intercept parameters.
   PARAMETER_VECTOR(log_growth_error);      // Growth increment error inflation parameters.
   
   // Density parameters:
   PARAMETER(log_lambda_alpha);             // Log-scale global mean density.
   PARAMETER_VECTOR(log_lambda_instar);     // Log-scale instar mean density.
   PARAMETER_VECTOR(log_lambda_tow);        // Log-scale tow mean density.
   PARAMETER_VECTOR(log_lambda_instar_tow); // Log-scale instar x tow mean density.
   PARAMETER(log_sigma_lambda_instar);      // Log-scale instar mean density error parameter.
   PARAMETER(log_sigma_lambda_tow);         // Log-scale tow mean density error parameter.
   PARAMETER(log_sigma_lambda_instar_tow);  // Log-scale instar x tow mean density error parameter.
   
   // Vector sizes:      
   int n = x.size();                        // Number of data elements.
   int n_instar = log_lambda_instar.size(); // Number instars in model.
   int n_tow = log_lambda_tow.size();       // Number of tows.
   
   // Calculate instar means and errors:
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
   
   // Initialize log-likelihood:
   Type v = 0;
   
   // Random effects:
   v -= sum(dnorm(log_lambda_instar, 0, exp(log_sigma_lambda_instar), true));          // Tow-level random effect.
   v -= sum(dnorm(log_lambda_tow, 0, exp(log_sigma_lambda_tow), true));                // Tow-level random effect.
   v -= sum(dnorm(log_lambda_instar_tow, 0, exp(log_sigma_lambda_instar_tow), true));  // Instar x tow effect.
   
   // Calculate mean instar abundances:
   matrix<Type> log_lambda(n_tow,n_instar);
   for (int t = 0; t < n_tow; t++){
      for (int k = 0; k < n_instar; k++){
         log_lambda(t,k) = log_lambda_alpha + 
                           log_lambda_instar[k] + 
                           log_lambda_tow[t] + 
                           log_lambda_instar_tow[t * n_instar + k];
      }
   }

   // Likelihood evaluation:
   for (int i = 0; i < n; i++){ 
      Type eta = 0;
      for (int k = 0; k < n_instar; k++){
         eta += exp(log_lambda(tow[i],k) + dnorm(x[i], mu[k], sigma[k], true) + log(swept_area[tow[i]] / 1000));
      }
      v -= dpois(f[i], eta, true);
   }
   
   // Likelihood evaluation:
   vector<Type> lambda_size(140); 
   for (int i = 0; i < 140; i++){ 
      lambda_size[i] = 0;
      for (int k = 0; k < n_instar; k++){
         lambda_size[i] += exp(log_lambda_alpha + log_lambda_instar[k]) * dnorm(Type(i+1), mu[k], sigma[k], false);
      }
   }
   
   // Export instar stats:
   REPORT(mu);
   REPORT(sigma);
   REPORT(log_lambda);
   REPORT(log_lambda_alpha);
   REPORT(log_lambda_instar);
   REPORT(log_lambda_tow); 
   REPORT(log_lambda_instar_tow);
   REPORT(lambda_size);
   
   return v;
}
