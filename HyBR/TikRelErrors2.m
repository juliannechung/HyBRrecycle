function E = TikRelErrors2(alpha, bhat, V, xtrue, S)
%
%   E = TikRelErrors(alpha, bhat, xhat, S)
%
%   This function computes the relative error for Tikhonov Regularization.
%      
%   Assume the problem has the form b = Ax + noise, with [U, S, V] = svd(A).
%
%   Input:
%       alpha - regularization parameter
%        bhat - U'*b, where b is the noisy data
%        xhat - V'*x_true
%           S - column vector containing singular or spectral values of A
%   Output:
%       E - relative error at the given alpha
%
%   J.Chung and J. Nagy 3/2007

bhat = bhat(:);
xtrue = xtrue(:);
E = zeros(length(alpha), 1);
for i = 1:length(E)
  y = V*(( (abs(S).^2) .* bhat ./ S ) ./ ( abs(S).^2 + alpha(i)^2 )) - xtrue;
  E(i) = sum(abs(y).^2);
end
E = sqrt(E)/norm(xtrue(:));