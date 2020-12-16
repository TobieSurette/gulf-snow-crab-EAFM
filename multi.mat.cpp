// Population abundance variables:
matrix<Type> n_imm(n_instar,n_year);   // Immatures.
matrix<Type> n_skp(n_instar,n_year);   // Immature skip-moulters.
array<Type>  n_mat(n_instar,n_year,6); // Matures.

n_ma
// Mature recruits for first year and first instar: 
for (int k = 0; k < 5; k++){   
   for (int a = 0; a < 6; a++){ 
      n_mat(k,0,a) = 0; // First year small mature recruits set to zero.
   }
} 
for (int y = 1; y < n_year; y++){ 
   n_mat(0,y,0) = 0; // First instar mature recruits set to zero.
}
for (int k = 5; k < n_instar; k++){
   for (int a = 0; a < 6; a++){ 
      
      n_mat(k,0,a) = exp(log_n_rec_instar_0();  // Mature abundances for first year. 
   }
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
