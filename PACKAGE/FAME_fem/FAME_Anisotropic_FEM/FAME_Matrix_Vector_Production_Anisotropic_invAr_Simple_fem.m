function vec_y = FAME_Matrix_Vector_Production_Anisotropic_invAr_Simple_fem( vec_x, B, Nx, Ny, Nz, N, Sigma_r, Pi_Qr, Pi_Qrs, D_k, D_ks, LS_tol )
    global inner_iter inner_count inner_cpu_time 
    inner_count = inner_count + 1;
    
    vec_y = Sigma_r.\vec_x;
    time_start_pcg = tic;
    [ vec_y, ~, ~, inner_iter(inner_count) ] = pcg(  @(x) FAME_Matrix_Vector_Production_Anisotropic_Ar_Simple_fem(x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks),...
              vec_y, LS_tol,1000 );
    inner_cpu_time(inner_count) = toc(time_start_pcg);
    vec_y = Sigma_r.\vec_y;
end