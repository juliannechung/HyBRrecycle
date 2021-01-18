function [x, alpha] = Tikhonovsolver(U, S, V, b, options, B, t, QV, m)
%
%         [x, alpha] = Tikhonovsolver(U, S, V, b, options)
%
%  This function computes a Tikhonov regularized LS solution to the
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
%            QV: Lanczos 
%             t: number of iteration
%
%  Output:
%           x - Tikhonov solution
%       alpha - regularization parameter
%
% J.Chung and J. Nagy 3/2007

bhat = U'*b;
bhat = bhat(:);

% Get values from the options structure.
alpha = HyBRget(options,'RegPar',[],'fast');
omega = HyBRget(options,'Omega',[],'fast');


switch alpha
    case 'dp'
        nLevel = HyBRget(options,'nLevel',[],'fast');
        alpha = fminbnd('TikDP', 0, S(1), [], V, bhat, S, nLevel, B, b, m);
    case 'gcv'
        alpha = fminbnd('TikGCV', 0, S(1), [], bhat, S);
    case 'wgcv'
        alpha = fminbnd('TikGCV', 0, S(1), [], bhat, S, omega);
    case 'nwgcv'
        new_omega = (t+1)/m;
        alpha = fminbnd('TikGCV', 0, S(1), [], bhat, S, new_omega);
    case 'upre'
        nLevel = HyBRget(options,'nLevel',[],'fast');
        alpha = fminbnd('TikUPRE', 0, S(1), [], bhat, S, nLevel);
        
    case 'optimal'

    x_true = HyBRget(options,'x_true',[],'fast');
    if isfield(options,'mu')
      mu = options.mu;
      errhan = @(lambda)TikOPT(lambda, V, bhat, S, QV, x_true(:)-mu);
    else
      errhan = @(lambda)TikOPT(lambda, V, bhat, S, QV, x_true(:));
    end
  
    options.TolX = eps;
    alpha = fminbnd(errhan, 0, S(1), options);
    
    
  otherwise
    if ~isa(alpha, 'double')
      error('Incorrect input argument for alpha.')
    end
end

D = abs(S).^2 + alpha^2;
bhat = conj(S) .* bhat(1:length(S));
xhat = bhat ./ D;
x = V * xhat;

end
%% -----------------------SUBFUNCTIONS----------------------------------
function [e, alpha] = OptTikTol(U, S, V, b, x)
%
%   [e, alpha] = OptTolTol(U, S, V, b, x)
%
%   This function computes the optimal regularization parameter
%       for Tikhonov regularization by minimizing:
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
%           alpha - optimal regularization parameter

bhat = U'*b;
k = length(S);
alpha = fminbnd('TikRelErrors',0,1,[],bhat(1:k),V'*x,S);
e = TikRelErrors(alpha, bhat(1:k), V'*x, S);
end

function [e, alpha] = OptTikTol2(U, S, V, b, options)
%
%   [e, alpha] = OptTolTol2(U, S, V, b, x)
% This computes the error in the image space.
%
%   This function computes the optimal regularization parameter
%       for Tikhonov regularization by minimizing:
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
%           alpha - optimal regularization parameter

V_GK = HyBRget(options, 'Vx');
x_true = HyBRget(options,'x_true',[],'fast');
bhat = U'*b;
k = length(S);
% alpha = fminbnd('TikRelErrors2',0,1,[],bhat(1:k),V_GK*V, x_true,S);
alpha = fminbnd('TikRelErrors2',0,S(1),[],bhat(1:k),V_GK*V, x_true,S);
e = TikRelErrors2(alpha, bhat(1:k), V_GK*V, x_true, S);
end