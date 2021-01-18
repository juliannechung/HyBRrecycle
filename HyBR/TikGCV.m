function G = TikGCV(alpha, bhat, s, omega)
%
%    G = TikGCV(alpha, bhat, s, omega)
%
%  This function evaluates the GCV function for Tikhonov
%  regularization.  
%
%  Input:  alpha -  regularization parameter
%           bhat -  vector U'*b, where U = left singular vectors
%              s -  vector containing the singular values
%          omega -  used to implement the weighted GCV function
%                   Default is omega = 1, which gives the standard GCV.
%
%  Output:     G -  the scalar G(alpha).
%
%  J.Chung and J. Nagy 3/2007

if nargin == 3
  omega = [];
end
if isempty(omega)
  omega = 1;
end

m = length(bhat);
n = length(s);
t0 = sum(abs(bhat(n+1:m)).^2);

s2 = abs(s) .^ 2;
alpha2 = alpha^2;

tt = 1 ./ (s2 + alpha2);
t1 = alpha2 .* tt;
t2 = abs(bhat(1:n) .* t1) .^2;
t3 = t1 + (1 - omega)*s2 .* tt;

G = n*(sum(t2) + t0) / (sum(t3) + m - n)^2;