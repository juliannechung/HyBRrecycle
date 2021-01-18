function [x, trun_index, trun_tol] = TSVDsolver(U, S, V, b, options)
%
%      [x, trun_index, trun_tol] = TSVDsolver(U, S, V, b, options)
%
%  This function computes a TSVD regularized LS solution to the
%  PROJECTED problem, using the identity matrix as the regularization operator.
%
%  Input: 
%    U, S, V - SVD of A, where S is a column vector containing
%              singular (or spectral) values of A.
%          b - noisy image
%    options - structure with the following necessary fields:
%         RegPar: [value | GCV | {WGCV} | optimal] where
%                     (a) a value for the regularization parameter
%                     (b) a method for choosing a reg. parameter
%                                (i) 'GCV' - standard GCV
%                                (ii) 'WGCV' - weighted GCV (default)
%                     (c) 'optimal' - finds the optimal reg. parameter
%                                (requires x_true)
%         Omega: if RegPar is 'WGCV', then omega must be 
%                           [value | {adapt}]
%                     (a) a value for omega
%                     (b) 'adapt' - use the adaptive method (default)
%            Vx: V_k * x_true (where V_k comes from Lanczos bidiagonlization)
%                 This vector is needed to compute the optimal parameter.
%
%  Output: 
%           x - TSVD solution
%  trun_index - index of singular value where truncation occurs
%    trun_tol - truncation tolerance
%
% J.Chung and J. Nagy 4/2007

bhat = U'*b;
bhat = bhat(:);

n1 = length(bhat);
n2 = length(S);

% Get values from the options structure.
tol = HyBRget(options,'RegPar',[],'fast');

switch tol
  case 'gcv'
    tol = TSVDGCV(bhat, S);
  case 'wgcv'
    omega = HyBRget(options,'Omega',[],'fast');
    tol = TSVDGCV(bhat, S, omega);
  case 'optimal'
    [e, trun_index, tol] = OptTSVDTol2(U, S, V, b, options);
  otherwise
    if ~isa(alpha, 'double')
      error('Incorrect input argument for alpha.')
    end
end
mask = (abs(S) >= tol);
trun_index = sum(sum(mask));
mask2 = ~mask;
S = (mask .* (1./(S + mask2)) );

x = V * (S .* bhat(1:n2));
trun_tol = tol;

if isreal(b)
  x = real(x);
end
end
%% -----------------------SUBFUNCTIONS----------------------------------

function tol = TSVDGCV(beta, s, omega)
%
%     tol = TSVDGCV(beta, s, omega)
%
% This function uses the GCV function to choose a default truncation
% tolerance for TSVD regularization.
%
% Input:  
%      beta - vector U'*b
%         s - vector containing the singular values
%     omega -  used to implement the "weighted" GCV function
%              Default is omega = 1, which gives usual GCV.
%
% Output: tol - truncation tolerance
%

if nargin == 2
  omega = [];
end
if isempty(omega)
  omega = 1;
end

n = length(s);
m = length(beta);

%B = sum(beta(1:n).^2);
%for k = 1:n-1
%    B = B - beta(k)^2;
%    G(k) = n*B/(m - omega*k)^2;
%end

%[G_min, i] = min(G);
%tol = s(i);


[s, idx] = sort(abs(s));
s = flipud(s);

beta1 = beta(1:n);
beta2 = beta(n+1:m);
beta1 = abs( flipud( beta1(idx) ) );
beta = [beta1; abs(beta2)];

%
%  Here we use the GCV function to pick a tolerance.
%
rho2 = zeros(m-1,1);
rho2(m-1) = beta(m)^2;
for k=m-2:-1:1
  rho2(k) = rho2(k+1) + beta(k+1)^2;
end
G = zeros(m-1,1);
for k=1:m-1
  G(k) = n*rho2(k)/(m - omega*k)^2;
end

[minG,reg_min] = min(G);
tol = s(reg_min);
end

function [e, trun_index, trun_tol] = OptTSVDTol(U, S, V, b, x)
%
%   [e, trun_index, trun_tol] = OptTolTol(U, S, V, b, x)
%
%   This function computes the optimal regularization parameter
%       for TSVD regularization by minimizing:
%           || filtered solution - true solution ||_2
%   
%   Assume the problem has the form b = Ax + noise, with [U, S, V] = svd(A).
%
%   Input:
%       [U, S, V] - SVD of A, where S is a column vector containing
%                       singular or spectral values of A
%               b - noisy data
%               x - true data
%
%   Output:
%               e - relative error
%      trun_index - index of singular value corresponding to truncation
%        trun_tol - optimal truncation tolerance

xhat = V'*x;
bhat = U'*b;
y = bhat(1:length(S)) ./ S - xhat;

y2 = abs(y).^2;
xhat2 = abs(xhat).^2;

e = zeros(length(S),1);
e(1) = y2(1) + sum(xhat2(2:end));
for i = 2:length(e)
  e(i) = e(i-1) - xhat2(i) + y2(i);
end
e = sqrt(e) / norm(x(:));
[mm, trun_index] = min(e);
trun_tol = S(trun_index);
end

function [e, trun_index, trun_tol] = OptTSVDTol2(U, S, V, b, options)
%
%   [e, trun_index, trun_tol] = OptTolTol(U, S, V, b, x)
%
%   This function computes the optimal regularization parameter
%       for TSVD regularization by minimizing:
%           || filtered solution - true solution ||_2
%   
%   Assume the problem has the form b = Ax + noise, with [U, S, V] = svd(A).
%
%   Input:
%       [U, S, V] - SVD of A, where S is a column vector containing
%                       singular or spectral values of A
%               b - noisy data
%               x - true data
%
%   Output:
%               e - relative error
%      trun_index - index of singular value corresponding to truncation
%        trun_tol - optimal truncation tolerance

V_GK = HyBRget(options, 'Vx');
x_true = HyBRget(options,'x_true',[],'fast');
V = V_GK*V;
bhat = U'*b;

y = bhat(1:length(S)) ./ S;
evec = y(1)*V(:,1) -x_true(:);
e = zeros(length(S),1);
e(1) = norm(evec);
for i = 2:length(e)
  evec = evec + y(i)*V(:,i);
  e(i) = norm(evec);
end
e = e / norm(x_true(:));
[mm, trun_index] = min(e);
trun_tol = S(trun_index);
end