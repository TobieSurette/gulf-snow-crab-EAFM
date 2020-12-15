// Population abundance variables:
matrix<Type> n_imm(n_instar,n_year);   // Immatures.
matrix<Type> n_skp(n_instar,n_year);   // Immature skip-moulters.
array<Type>  n_mat(n_instar,n_year,6); // Matures.

// Mature recruits for first year and first instar: 
for (int k = 0; k < 5; k++){   
   n_mat(k,0,0) = 0; 
} 
for (int k = 5; k < n_instar; k++){
   n_mat(k,0,0) = exp(log_n_rec_instar_0[k-5]);  
}
for (int y = 1; y < n_year; y++){ 
   n_mat(0,y,0) = 0;  
}

// First year:
for (int k = 0; k < 5; k++){   
   n_mat(k,0,1) = 0; 
   n_mat(k,0,2) = 0;
   n_mat(k,0,3) = 0;
   n_mat(k,0,4) = 0;
   n_mat(k,0,5) = 0;
} 
// Matures first year:
for (int k = 5; k < n_instar; k++){
   n_mat(k,0,1) = exp(log_n_res_instar_0[k-5]);  
   n_mat(k,0,2) = exp(log_n_res_instar_0[k-5]);
   n_mat(k,0,3) = exp(log_n_res_instar_0[k-5]);
   n_mat(k,0,4) = exp(log_n_res_instar_0[k-5]);
   n_mat(k,0,5) = exp(log_n_res_instar_0[k-5]);
}
// First instar:
for (int y = 1; y < n_year; y++){ 
   n_mat(0,y,1) = 0;
   n_mat(0,y,2) = 0;
   n_mat(0,y,3) = 0;
   n_mat(0,y,4) = 0;
   n_mat(0,y,5) = 0;
} 

// Population dynamics equations:
for (int k = 1; k < n_instar; k++){
   for (int y = 1; y < n_year; y++){
      n_imm(k,y) = (Type(1)-p_mat(k-1,y-1)) * (Type(1)-p_skp[k-1]) * (1-M_imm) * n_imm(k-1,y-1); 
      n_skp(k,y) = (Type(1)-p_mat(k-1,y-1)) * p_skp[k-1] * (1-M_imm) * n_imm(k,y-1);    
      n_mat(k,y,0) = (Type(1)-M_mat[0]) * ((Type(1)-p_skp[k-1]) * p_mat(k-1,y-1) * n_imm(k-1,y-1) + n_skp(k-1,y-1)); // Recruits.
      n_mat(k,y,1) = (Type(1)-M_mat[1]) * n_mat(k,y-1,0); // First-year residual.
      n_mat(k,y,2) = (Type(1)-M_mat[1]) * n_mat(k,y-1,1); // Second-year residual.
      n_mat(k,y,3) = (Type(1)-M_mat[1]) * n_mat(k,y-1,2); // Third-year residual.
      n_mat(k,y,4) = (Type(1)-M_mat[1]) * n_mat(k,y-1,3); // Fourth-year residual.
      n_mat(k,y,5) = (Type(1)-M_mat[1]) * n_mat(k,y-1,4); // Fifth-year residual.
   }
}
