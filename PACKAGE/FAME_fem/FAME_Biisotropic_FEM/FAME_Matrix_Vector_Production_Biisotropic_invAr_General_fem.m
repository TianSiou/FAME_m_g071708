function vec_y = FAME_Matrix_Vector_Production_Biisotropic_invAr_General_fem( vec_x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Pr, Pi_Qrs, Pi_Prs, D_kx, D_ky, D_kz, LS_tol)
    global inner_iter inner_count inner_cpu_time 
    inner_count = inner_count + 1;
    time_start_bicgstabl = tic;
    [ vec_y, ~, ~, inner_iter(inner_count) ] = pcg(  @(x) FAME_Matrix_Vector_Production_Biisotropic_Ar_General_fem(x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Pr, Pi_Qrs, Pi_Prs, D_kx, D_ky, D_kz),...
              vec_x, LS_tol,1000 );
    inner_cpu_time(inner_count) = toc(time_start_bicgstabl);
end