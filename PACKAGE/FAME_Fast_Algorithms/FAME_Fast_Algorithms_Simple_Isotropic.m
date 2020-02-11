function [ freq, Ele_field, cpu_time, LS_iter, LS_cpu_time ] = FAME_Fast_Algorithms_Simple_Isotropic( grid_num, wave_vec, B, Lambdas, eigen_wanted)
    global inner_iter inner_count inner_cpu_time
    inner_count = 0;
    shift = 0 * pi^2;
    
    Nx = grid_num(1); Ny = grid_num(2); Nz = grid_num(3);
    
    N = Nx*Ny*Nz;
    if norm(wave_vec) == 0
        Nd = 2*N - 2;
    else 
        Nd = 2*N;
    end
    
%     Lambdas.D_k  = kron(Lambdas.D_kz, kron(Lambdas.D_ky, Lambdas.D_kx));
%     Lambdas.D_ks = conj(Lambdas.D_k);
%     invB_eps     = 1./B_eps;
    
    opt.isreal = 0;
    opt.issym  = 1;
    
    LS_tol  = 1e-12;
    no_orth = 5; 
    ev = [];
    ew = [];
    cpu_time_start = tic;
%     invAr_fun = @(x) FAME_Matrix_Vector_Production_Isotropic_invAr_Simple( x, B, Nx, Ny, Nz, N, Lambdas.Sigma_r, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, LS_tol);
    if shift == 0
        
        for i = 1 : (eigen_wanted / no_orth)
            if i == 1
                [pre_ev,pre_ew] = eigs( @(x) FAME_Matrix_Vector_Production_Isotropic_invAr_Simple( x, B, Nx, Ny, Nz, N,...
                    Lambdas.Sigma_r, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, LS_tol),...
                    Nd, no_orth, 'lm', opt );
            else 
                fun_mtx_uu = @(vec_x) ones(Nd,1) .* vec_x - (pre_ev * (pre_ev' * vec_x) );
                fun_mtx_A = @(x) FAME_Matrix_Vector_Production_Isotropic_invAr_Simple( x, B, Nx, Ny, Nz, N,...
                    Lambdas.Sigma_r, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, LS_tol);
                fun_mtx_uAu = @(vec_x) fun_mtx_uu( fun_mtx_A( fun_mtx_uu(vec_x) ) );
                
                [pre_ev,pre_ew] = eigs( @(x)  fun_mtx_uAu(x), Nd, no_orth, 'lm', opt );

            end
                ev = [ev, pre_ev];
                ew = [ew; diag(pre_ew)];
        end
        ew = diag(ew);
    else
    [ev,ew] = eigs( @(x) FAME_Matrix_Vector_Production_Isotropic_shift_invAr_Simple( x, B, Nx, Ny, Nz, N,...
                         Lambdas.Sigma_r, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, LS_tol, shift),...
                         Nd, eigen_wanted, 'lm', opt );        
    end
    cpu_time = toc(cpu_time_start);
    
    LS_iter     = inner_iter;
    LS_cpu_time = inner_cpu_time;
    
    clear('inner_iter','inner_cpu_time');
    
%     ew = 1./diag(ew);

    ew = diag(ew);
    ew = 1./ew + shift;
    [ Ele_field, freq ] = Eigen_Restoration_Simple_Isotropic( ev, ew, Nx, Ny, Nz, N, wave_vec, B, Lambdas, shift );
end

function [ Output_eigvec_mat, Output_eigval ] = Eigen_Restoration_Simple_Isotropic( Input_eigvec_mat, Input_eigval, Nx, Ny, Nz, N, wave_vec, B, Lambdas, shift )
%% Start to restore the eigenvectors
    for column_idx = 1 : size(Input_eigvec_mat,2)
        vec_y = Lambdas.Sigma_r.*Input_eigvec_mat(:,column_idx);
        vec_y = FAME_Matrix_Vector_Production_Qr_Simple(vec_y, Nx, Ny, Nz, N, Lambdas.Pi_Qr, Lambdas.Pi_Qrs, Lambdas.D_k, Lambdas.D_ks, 'normal');

%         Output_eigvec = invB_eps.*vec_y;     
        Output_eigvec = FAME_Matrix_Vector_Production_invB_Isotropic(vec_y, B);
% Normalizing the eigenvector     
        Output_eigvec = Output_eigvec / norm(Output_eigvec);           
        Output_eigvec_mat(:,column_idx) = Output_eigvec;        
    end
    Output_eigval = Input_eigval;
    if shift == 0
        if  norm(wave_vec) == 0  
            Output_eigvec_mat = [ ones(3*N,2), Output_eigvec_mat ];
            Output_eigval     = [0;0;Output_eigval];
            Output_eigvec_mat = Output_eigvec_mat(:,1:end-2);
            Output_eigval     = Output_eigval(1:end-2);
        end
    end
    Output_eigval = sqrt(Output_eigval);
    [Output_eigval,idx] = sort(real(Output_eigval));
    Output_eigvec_mat   = Output_eigvec_mat(:,idx);
end
function vec_y = FAME_Matrix_Vector_Production_Isotropic_preCM_Simple( vec_x, grid_num, B, Nx, Ny, Nz, N, C_sc, Cs_sc, shift)
    global inner_iter inner_count preCM_cpu_time 
    inner_count = inner_count + 1;
    k     = grid_num(3) / grid_num(1);
    P     = [];
 
    for i = 1:k
        ei    = sparse(k,1);
        ei(i) = 1;
        P     = [P,kron(speye(3),  kron(ei, speye(N) ) )];
    end
    vec_x = P' * vec_x;
    B_p   = P' * B * P;
    tic
    for i = 1 : k
        
    end
    preCM_cpu_time(inner_count) = toc
    
    [ vec_y, ~, ~, inner_iter(inner_count) ] = pcg(  @(x) FAME_Matrix_Vector_Production_Isotropic_Ar_Simple(x, B, Nx, Ny, Nz, N, Pi_Qr, Pi_Qrs, D_k, D_ks),...
              vec_y, LS_tol,1000 );
  
    
end

function vec_y = FAME_mxt_product_Iuu(vec_x, A, u, Nd)
         vec_y = ones(Nd,1) .* vec_x - (u .* (u' * vec_x) );
         vec_y = A * vec_y;
         vec_y = ones(Nd,1) .* vec_y - (u .* (u' * vec_x) );
end