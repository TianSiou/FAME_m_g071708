% clear
% clc
addpath('C:\Users\User\Desktop\FAME_m_g071708-test_for_git\SIRA_GEP');
addpath('C:\Users\User\Desktop\FAME_m_g071708-test_for_git\SIRA_SEP');
% x_end = 3*pi/2;
edge_len = [1,1,1];
grid_num = [16,16,16];
mesh_len = edge_len./grid_num;
%% Construct FDM matrices
ones_m = @(m) ones(m,1);
K_m    = @(m) spdiags(ones_m(m-1),0,m,m-1) - spdiags(ones_m(m-1),-1,m,m-1); % sparse matrix form
I_m    = @(m) speye(m);

n1 = grid_num(1); n2 = grid_num(2); n3 = grid_num(3);
Nr = (n1-1)* n2   * n3    +  n1   *(n2-1)* n3    +  n1   * n2   *(n3-1);
Nc =  n1   *(n2-1)*(n3-1) + (n1-1)* n2   *(n3-1) + (n1-1)*(n2-1)* n3   ;
N0_hat = (n1-1)*(n2-1)*(n3-1);


% [V_n1,Sigma_n1,U_n1] = svd(full(K_m(n1)));
% V_n1 = V_n1(:,1:end-1);
% Sigma_n1 = Sigma_n1(1:end-1,:) / mesh_len(1);
% [V_n2,Sigma_n2,U_n2] = svd(full(K_m(n2)));
% V_n2 = V_n2(:,1:end-1);
% Sigma_n2 = Sigma_n2(1:end-1,:) / mesh_len(2);
% [V_n3,Sigma_n3,U_n3] = svd(full(K_m(n3)));
% V_n3 = V_n3(:,1:end-1);
% Sigma_n3 = Sigma_n3(1:end-1,:) / mesh_len(3);
% 
% V_1_hat = kron( I_m(n3-1), kron( I_m(n2-1), V_n1        ) );
% V_2_hat = kron( I_m(n3-1), kron( V_n2     , I_m(n1-1)   ) );
% V_3_hat = kron( V_n3     , kron( I_m(n2-1), I_m(n1-1)   ) );
% U_1_hat = kron( I_m(n3-1), kron( I_m(n2-1), U_n1        ) );
% U_2_hat = kron( I_m(n3-1), kron( U_n2     , I_m(n1-1)   ) );
% U_3_hat = kron( U_n3     , kron( I_m(n2-1), I_m(n1-1)   ) );
% V_hat   = kron( V_n3     , kron( V_n2     , V_n1        ) );
% U_hat   = kron( U_n3     , kron( U_n2     , U_n1        ) );
% 
% 
% Sigma_1_hat = kron( I_m(n3-1), kron( I_m(n2-1), Sigma_n1        ) );
% Sigma_2_hat = kron( I_m(n3-1), kron( Sigma_n2 , I_m(n1-1)   ) );
% Sigma_3_hat = kron( Sigma_n3 , kron( I_m(n2-1), I_m(n1-1) ) );
% Sigma_1_til = kron( I_m(n3)  , kron( I_m(n2)  , Sigma_n1  ) );
% Sigma_2_til = kron( I_m(n3)  , kron( Sigma_n2 , I_m(n1)   ) );
% Sigma_3_til = kron( Sigma_n3 , kron( I_m(n2)  , I_m(n1)   ) );
% inv_Sigma_1_hat = spdiags(1./diag(Sigma_1_hat),0,size(Sigma_1_hat,1),size(Sigma_1_hat,1));
% inv_Sigma_2_hat = spdiags(1./diag(Sigma_2_hat),0,size(Sigma_2_hat,1),size(Sigma_2_hat,1));
% inv_Sigma_3_hat = spdiags(1./diag(Sigma_3_hat),0,size(Sigma_3_hat,1),size(Sigma_3_hat,1));
% inv_Sigma_1_til = spdiags(1./diag(Sigma_1_til),0,size(Sigma_1_hat,1),size(Sigma_1_hat,1));
% inv_Sigma_2_til = spdiags(1./diag(Sigma_2_til),0,size(Sigma_2_hat,1),size(Sigma_2_hat,1));
% inv_Sigma_3_til = spdiags(1./diag(Sigma_3_til),0,size(Sigma_3_hat,1),size(Sigma_3_hat,1));
% Sigma_q_hat     = sqrt(Sigma_1_hat.^2 + Sigma_2_hat.^2 + Sigma_3_hat.^2);
% % Sigma_q_til     = sqrt(Sigma_1_til.^2 + Sigma_2_til.^2 + Sigma_3_hat.^2);
% inv_Sigma_q_hat = spdiags(1./diag(Sigma_q_hat),0,size(Sigma_q_hat,1),size(Sigma_q_hat,1));
% % inv_Sigma_q_til = spdiags(1./diag(Sigma_q_til),0,size(Sigma_q_til,1),size(Sigma_q_til,1));

C_1_hat = kron( I_m(n3-1), kron( I_m(n2-1), K_m(n1)   ) ) / mesh_len(1);
C_1_til = kron( I_m(n3)  , kron( I_m(n2)  , K_m(n1)   ) ) / mesh_len(1);
C_1y    = kron( I_m(n3)  , kron( K_m(n2)  , I_m(n1-1) ) ) / mesh_len(2);
C_1z    = kron( K_m(n3)  , kron( I_m(n2)  , I_m(n1-1) ) ) / mesh_len(3);
C_2x    = kron( I_m(n3)  , kron( I_m(n2-1), K_m(n1)   ) ) / mesh_len(1);
C_2_hat = kron( I_m(n3-1), kron( K_m(n2)  , I_m(n1-1) ) ) / mesh_len(2);
C_2_til = kron( I_m(n3)  , kron( K_m(n2)  , I_m(n1)   ) ) / mesh_len(2);
C_2z    = kron( K_m(n3)  , kron( I_m(n2-1), I_m(n1)   ) ) / mesh_len(3);
C_3x    = kron( I_m(n3-1), kron( I_m(n2)  , K_m(n1)   ) ) / mesh_len(1);
C_3y    = kron( I_m(n3-1), kron( K_m(n2)  , I_m(n1)   ) ) / mesh_len(2);
C_3_hat = kron( K_m(n3)  , kron( I_m(n2-1), I_m(n1-1) ) ) / mesh_len(3);
C_3_til = kron( K_m(n3)  , kron( I_m(n2)  , I_m(n1)   ) ) / mesh_len(3);

O_1 = sparse(size(C_1z,1),size(C_2z,2));
O_2 = sparse(size(C_2x,1),size(C_3x,2));
O_3 = sparse(size(C_3y,1),size(C_1y,2));

C = [  O_1, -C_1z,  C_1y ;
       C_2z,  O_2, -C_2x ;
      -C_3y,  C_3x,  O_3  ];


% Q0 = [C_1_hat;C_2_hat;C_3_hat]*(U_hat*inv_Sigma_q_hat*U_hat'); 
% Q1_tilde = [V_1_hat*inv_Sigma_1_hat*U_1_hat';V_2_hat*inv_Sigma_2_hat*U_2_hat';V_3_hat*inv_Sigma_3_hat*U_3_hat'];
% Q1 = Q1_tilde - Q0*(Q0'*Q1_tilde);
% P0 = [C_1_til';C_2_til';C_3_til'];
% test_Q0 = norm(C *Q0,'fro')
% test_P0 = norm(C'*P0,'fro')

% Q1 = [V_1*sqrt(Sigma_1);V_2*sqrt(Sigma_2);V_3*sqrt(Sigma_3)] - Q0*(Q0'*[V_1*sqrt(Sigma_1);V_2*sqrt(Sigma_2);V_3*sqrt(Sigma_3)]);
% Q1_tilde = [V_1_hat*inv_Sigma_1_hat*U_1_hat';V_2_hat*inv_Sigma_2_hat*U_2_hat';V_3_hat*inv_Sigma_3_hat*U_3_hat'];
% Q1 = Q1_tilde - Q0*(Q0'*Q1_tilde);
% Q2 = C*Q1_tilde;
% % Q1 = [V_1;V_2;V_3] - Q0*(Q0'*[V_1;V_2;V_3]);
% % Q1 = [V_1;V_2;V_3];
% Q1 = orth(full(Q1));
% Q2 = orth(full(Q2));
% test_orth_Q0_Q1 = norm(Q0'*Q1)
% 
% tmp = sqrt(eig(C'*C));
% tmp1 = sqrt(eig(Q1'*(C'*C)*Q1));
% tmp2 = sqrt(eig(Q2'*(C'*C)*Q2));


% sqrt(eig(Q1'*(C'*C)*Q1))
% Q1 = [V_1 - V_1*Sigma_1*U_1'*U*(inv_Sigma_q.^2)*U'*(U_1*Sigma_1 + U_2*Sigma_2 + U_3*Sigma_3);
%       V_2 - V_2*Sigma_2*U_2'*U*(inv_Sigma_q.^2)*U'*(U_1*Sigma_1 + U_2*Sigma_2 + U_3*Sigma_3);
%       V_3 - V_3*Sigma_3*U_3'*U*(inv_Sigma_q.^2)*U'*(U_1*Sigma_1 + U_2*Sigma_2 + U_3*Sigma_3)];
% norm(Q1 - ([V_1;V_2;V_3]-Q0*(Q0'*[V_1*Sigma_1;V_2*Sigma_2;V_3*Sigma_3])),'fro')
%   
%   Q1_test = orth(full(Q1));
%%%%%%%%%%%%%%%%%%%%%%%%% test by siu
eigen_wanted = 500;
% B.B_eps     = ones(Nc, 1);
% B.invB_eps  = 1./ B.B_eps;
B.B_eps     = spdiags(B.B_eps, 0, Nc, Nc);
sigma = 0.1;
global count inner_count inner_iter
count = 0;
inner_count = 0;



%%%%%%%%%%%%%%% test by Lanczos
% LANCZOS_time_start = tic;
% % [ev, ew] = eigs( @(vec_x) inv_A(vec_x, C' * C, B, sigma, 1e-10, 1000), Nc,  inv_B , eigen_wanted, 'lm');
% [ev, ew] = eigs( @(vec_x) mxt_product_BAB(vec_x, C' * C, B.B_eps, sigma, 1e-10, 1000), Nc, 1, 'lm');%%,'SubspaceDimension',40);
% ew = diag(ew);
% ew = (1 ./ ew) + sigma;
% ev = (B.invB_eps .^ (1/2)) .* ev;
% cpu_time_LANCZOS = toc(LANCZOS_time_start);
%%%%%
%%%%%%%%%%%%%%%% test by sira_GEP
SIRA_time_start = tic;
        no_restart     = max(100, 2*eigen_wanted);
        stop_tolerance = 1.0e-10; %* sqrt(Nd); %1.0e-12;
        initial_vec  = ones(Nc,1);
%         initial_vec  = ev;
        initial_vec  = initial_vec / norm(initial_vec);
%         initial_vec  = ev;
        target         = sigma;
        target_type    = 'CL';%CL%%RGTC
        CorEq          = 'SIRA';
        LSinfo.precond = 'no';
        LSinfo.solver  = 'minres'; %'minres'; %'pcg';
        fid            = fopen('runsh_parameter.txt','r');
        [ew, ev] = GEP_AB_Herm_JDSIRA_Driver (C' * C, B.B_eps, Nc, no_restart, eigen_wanted, ...
                stop_tolerance, initial_vec, target, fid, target_type, CorEq, ...
                LSinfo);
%%%%%%
cpu_time_SIRA = toc(SIRA_time_start); 


for i = 1 : eigen_wanted
    Err(i) = norm( (C' * C * ev(:,i)  ) - ew(i) * (B.B_eps * ev(:,i) ) ); 
end
% function vec_y = inv_A(vec_x, A, B, sigma, LS_tol, maxit)
%     vec_y = minres(A - sigma * B, vec_x, LS_tol, maxit );
% end
% fclose(fid);
function vec_y = mxt_product_BAB(vec_x, A, B, sigma, LS_tol, maxit)
         global inner_count inner_iter
         inner_count = inner_count + 1;
         vec_y = (B .^(1/2) ) * vec_x;
         [vec_y, flag, res, inner_iter(inner_count)] = bicgstabl(A - sigma * B, vec_y, LS_tol, maxit);
         vec_y = (B .^(1/2) ) * vec_y;
         res
         flag
end
% function vec_y = mxt_product_BCCB(vec_x, C, B, sigma, LS_tol, maxit, alpha, Nc)
%          global inner_count inner_iter
%          inner_count = inner_count + 1;
%          vec_y = (B .^(1/2) ) * vec_x;
%          [vec_y, flag, res, inner_iter(inner_count)] = minres( @(x) mtx_product_shiftCCB(x, C, B, sigma), vec_y, LS_tol, maxit, C' * C + alpha*speye(Nc) );
%          vec_y = (B .^(1/2) ) * vec_y;
%          res
%          flag
% end
% 
% function vec_y = mtx_product_shiftCCB(vec_x, C, B, sigma)
%          vec_y = C' * (C * vec_x) - sigma * B *vec_x;
% end