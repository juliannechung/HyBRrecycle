function U = TikUPRE(alpha, bhat, s, sigma)
%
%   U = TikUPRE(alpha, bhat, s, sigma)
%
%   This function evaluates the UPRE function for standard form Tikhonov.
%
%   Input:
%       alpha - regularization parameter
%        bhat - U'*b
%           s - singular values
%       sigma - noise estimate
%
%   Output:
%           U - value of the function U(alpha)
%
% J. Chung 6/18/14

% bhat = bhat(:).^2;
bhat = bhat(:);
s = s(:);
s2 = abs(s).^2; 
alpha2 = alpha^2;
sigma2 = sigma^2;

m = length(bhat);
n = length(s);

t0 = sum(abs(bhat(n+1:m)).^2);

tt = 1 ./ (s2 + alpha2);
t2 = abs(alpha2*(bhat(1:n) .* tt)) .^2;
t3 = s2 .* tt;

% Compute UPRE function
U = (sum(t2) + t0 + 2*sigma2*sum(t3))-n*sigma2;
