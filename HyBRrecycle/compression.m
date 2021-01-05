function [kk,W,Y,R] = compression(A,B,V,trunc_options,Vy,vector,y,W)
%
% compression is a subroutine function for basis compression. Here we
% provide four compression approaches: SVD, solution, sparse and RBD.
%
% Inputs:
%  A:             given matrix
%  B:             previous bidiagonal matrix
%  V:             [W V_l], where V_l is the augmented bases and W is the bases we kept after previous compression
%  trunc_options: compression strategies: SVD, solution, sparse and RBD; the storage limit
%  Vy:            the initial solution
%  vector:        right side of the projected problem
%  y:             the solution of the projected problem (required by solution-oriented compression)
%  W:             the bases we kept after previous compression

% Outputs:
%  kk:             number of bases after compression
%  W:              compressed basis vectors
%  Y:              orthogonal matrix from QR factorization: AW = YR
%  R:              upper triangular matrix from QR factorization: AW = YR
%
% J. Chung, E. de Sturler, and J. Jiang, 2020

max_mm = trunc_options.max_mm; % maximum number of bases we can keep after compression
tol = 10^(-6); % tolerance for compression

switch trunc_options.compress
  case 'SVD' % TSVD compression
    [Ub, Sb, Vb] = svd(B,0);
    trun_index = sum(diag(Sb) > tol); % number of bases we can keep based on the tolerance
    kk = min(max_mm,trun_index);
    
    W = V*Vb(:,1:kk);
    Y = zeros(size(A,1),size(W,2));
    for j = 1:kk
      Y(:,j) = A*W(:,j);
    end
    
    [Y,R] = qr(Y,0);
    x_out = Vy;
    
    eta = W'*x_out;
    xhat = x_out - W*eta; etahat = norm(xhat);
    xhat = xhat / etahat;
    W = [W xhat];
    Axhat = A*xhat;
    zeta = Y'*Axhat;
    yhat = Axhat - Y*zeta;
    zetahat = norm(yhat);
    yhat  = yhat/zetahat;
    Y = [Y yhat];
    R = [R , zeta ; zeros(1,kk) , zetahat];
    kk = kk+1;
    
  case 'solution' % solution oriented compression
    idx = find(y~=0);
    kk = min(max_mm,length(idx));
    sort_abs_y = sort(abs(y),'descend');
    
    sort_abs_y = sort_abs_y(sort_abs_y >= tol);
    kk = min(kk,length(sort_abs_y));
    lowb = sort_abs_y(kk);
    idx = abs(y)>=lowb;
    
    W = V(:,idx);
    Y = zeros(size(A,1),size(W,2));
    for j = 1:kk
      Y(:,j) = A*W(:,j);
    end
    [Y,R] = qr(Y,0);
    
    x_out = Vy;
    eta = W'*x_out; xhat = x_out - W*eta; etahat = norm(xhat);
    xhat = xhat / etahat;
    W = [W xhat];
    
    Axhat = A*xhat;
    zeta = Y'*Axhat;
    yhat = Axhat - Y*zeta;
    zetahat = norm(yhat);
    yhat  = yhat/zetahat;
    Y = [Y yhat];
    R = [R , zeta ; zeros(1,kk) , zetahat];
    kk = kk+1;
    
  case 'sparse' % sparse solution is used
    lambda = .0005; % Regularization parameter for sparsity
    ysparse= SpaRSA(vector,B,lambda,'Verbose',0);
    idx = find(ysparse~=0);
    
    kk = min(max_mm,length(idx));
    sort_abs_ysparse = sort(abs(ysparse),'descend');
    lowb = sort_abs_ysparse(kk);
    idx = abs(ysparse)>=lowb;
    
    W = V(:,idx);
    Y = zeros(size(A,1),size(W,2));
    for j = 1:kk
      Y(:,j) = A*W(:,j);
    end
    [Y,R] = qr(Y,0);
    
    x_out = Vy;
    eta = W'*x_out; xhat = x_out - W*eta; etahat = norm(xhat);
    xhat = xhat / etahat;
    W = [W xhat];
    
    Axhat = A*xhat;
    zeta = Y'*Axhat;
    yhat = Axhat - Y*zeta;
    zetahat = norm(yhat);
    yhat  = yhat/zetahat;
    Y = [Y yhat];
    R = [R , zeta ; zeros(1,kk) , zetahat];
    kk = kk+1;
    
  case 'RBD' % Reduced Basis Decomposition compression
    [Bases, ~,~] = RBD(B', tol, max_mm);
    kk = size(Bases,2);
   
    W = V*Bases;
    Y = zeros(size(A,1),size(W,2));
    for j = 1:kk
      Y(:,j) = A*W(:,j);
    end
    [Y,R] = qr(Y,0);
    
    x_out = Vy;
    eta = W'*x_out; xhat = x_out - W*eta; etahat = norm(xhat);
    xhat = xhat / etahat;
    W = [W xhat];
    
    Axhat = A*xhat;
    zeta = Y'*Axhat;
    yhat = Axhat - Y*zeta;
    zetahat = norm(yhat);
    yhat  = yhat/zetahat;
    Y = [Y yhat];
    R = [R , zeta ; zeros(1,kk) , zetahat];
    kk = kk+1;
end

end