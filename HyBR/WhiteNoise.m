function [N, sigma] = WhiteNoise(b, level, seed)
%
%      [N,sigma] = WhiteNoise(b, level, seed);
%
%  This function generates Gaussian white noise for the 
%  data b. 
%
%  Input:  b - array containing data
%      level - scalar in [0, 1] specifiying level (percentage) of
%              noise.  For example, level = 0.01 implies
%                  norm(N)/norm(b) = 0.01, or 1% noise
%              Default is level = 0.01.
%       seed - Used to set the random number generator.
%              Default is seed = 0.
%
%  Output: N - array same dimension as b, containing pseudo-random
%              values drawn from a normal distribution with mean zero
%              and standard deviation one, and scaled as described above.
%          sigma - standard deviation of the noise
%

% Check inputs and set default values.
if nargin == 1, level = [];, seed = [];, end
if nargin == 2, seed = [];, end
if isempty(level), level = 0.01;, end
if isempty(seed), seed = 0;, end

% Generate noise.
% randn('seed', seed);
rng(seed)
N = randn(size(b));
nN = norm(N(:));
N = N / nN;
N = level*norm(b(:))*N;
sigma = level*norm(b(:))/nN;
