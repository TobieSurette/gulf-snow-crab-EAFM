#include <TMB.hpp>
template<class Type> Type objective_function<Type>::operator()(){
  // Data:
  DATA_VECTOR(x_imm);          // Immature size measurements (n_imm).
  DATA_VECTOR(f_imm);          // Immature frequency observations (n_imm).
  DATA_VECTOR(year_imm);       // Immature survey year (n_imm).
  DATA_VECTOR(x_mat);          // Mature size measurements (n_mat).
  DATA_VECTOR(f_mat);          // Mature frequency observations (n_mat).
  DATA_VECTOR(year_mat);       // Immature survey year (n_mat).
  DATA_SCALAR(dx);             // Size-bin width.   
  
  // Instar growth parameters:                       
  PARAMETER(mu0);                          // First instar mean size.
  PARAMETER(log_sigma0);                   // Log-scale standard error for first instar.
  PARAMETER_VECTOR(log_hiatt_slope);       // Hiatt slope parameters.
  PARAMETER_VECTOR(log_hiatt_intercept);   // Hiatt intercept parameters.
  PARAMETER_VECTOR(log_growth_error);      // Growth increment error inflation parameters.
  
  PARAMETER_VECTOR(log_mu_year);           // Log-scale instar mean year interaction.
  PARAMETER(log_sigma_mu_year);            // Instar mean year interaction error term.
  
  // Density parameters:
  PARAMETER(log_lambda_alpha);             // Log-scale global mean density.
  PARAMETER_VECTOR(log_lambda);            // Log-scale instar mean density.
  PARAMETER(log_sigma_lambda);             // Log-scale instar mean density error parameter.
  PARAMETER_VECTOR(logit_scale);           // Logit-scale factors.
  PARAMETER_VECTOR(log_lambda_year);
  PARAMETER(log_sigma_lambda_year);
  
  // Vector sizes:      
  int n_imm = x_imm.size();                // Number of immature observations.
  int n_mat = x_mat.size();                // Number of mature observations.
  int n_instar = log_lambda_instar.size(); // Number instars in model.
  
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
        mu_mat(k,y) = mu_mat;
     }
  }
  
  // Annual instar abundances:
  v -= sum(dnorm(log_lambda, 0, exp(log_sigma_lambda), true));  
  v -= sum(dnorm(log_lambda_year, 0, exp(log_sigma_lambda_year), true));  
  matrix<Type> lambda_imm(n_instar,n_year);
  for (int k = 0; k < (n_instar); k++){
    for (int y = 0; y < (n_year); y++){
       lambda_imm(k,y) = exp(log_lambda_alpha + log_lambda[k] + log_lambda_year[y * n_instar + k]);
    }
  }  
    
  // Likelihood evaluation for immatures:
  for (int i = 0; i < n_imm; i++){ 
    Type eta_imm = 0;
    for (int k = 0; k < n_instar; k++){
        eta_imm += lambda_imm[k] * (pnorm(x_imm[i] + dx / 2, mu_imm[k], sigma[k]) - 
                                    pnorm(x_imm[i] - dx / 2, mu_imm[k], sigma[k]));
    }
    v -= dpois(f_imm[i], eta, true);
  }
  
  // Likelihood evaluation for matures:
  vector<Type> lambda_mat(n_instar + 1);
  
  vector<Type> scale = Type(1) / (Type(1) + exp(-logit_scale));
  
  lambda_mat[0] = 0;
  for (int k = 1; k < (n_instar+1); k++){
    lambda_mat[k] = (Type(1) / (Type(1) + exp(-logit_scale[k-1]))) * lambda_imm[k-1];
  }
  for (int i = 0; i < n_mat; i++){ 
    Type eta = 0;
    for (int k = 0; k < (n_instar+1); k++){
      eta += lambda_mat[k] * (pnorm(x_mat[i] + dx / 2, mu[k], sigma[k]) - 
        pnorm(x_mat[i] - dx / 2, mu[k], sigma[k]));
    }
    v -= dpois(f_mat[i], eta, true);
  }
  
  // Export results:
  REPORT(mu);
  REPORT(sigma);
  REPORT(log_lambda_instar);
  REPORT(logit_scale);
  REPORT(scale);
  REPORT(lambda_imm);
  REPORT(lambda_mat);
  
  return v;
}
