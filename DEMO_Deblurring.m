% DEMO_Deblurring.m
%
% This script sets up and computes reconstructions for a 2D deblurring 
%   problem using hybrid projection methods with recycling.
%
%      Chung, de Sturler, and Jiang. "Hybrid Projection Methods with 
%           Recycling for Inverse Problems". SISC, 2020.
%
% Note: This code requires IRTools: https://github.com/jnagy1/IRtools
%
% J. Chung, E. de Sturler, J. Jiang 12/2020

%% Get problem
load Grain

gamma = [4,7];
[PSF, center] = Gauss( gamma , size(x_true,1));
PSF(center(1)+1:end,center(2)+1:end) = 0;
PSF = PSF/sum(PSF(:));
    
A = psfMatrix(PSF, 'reflexive', center);
b = A*x_true;

nLevel = 0.002; % noise level
[N, sigma] = WhiteNoise(b, nLevel, 0); % add noise 
bn = b+N;

figure(1); 
subplot(1,3,1), imshow(x_true,[]), title('True Image')
subplot(1,3,2), imshow(bn,[]), title('Noisy Blurred Image')
subplot(1,3,3), imshow(PSF,[]), title('PSF')

%% Standard LSQR (using HyBR code to get subproblem)
maxit = 250;
input = HyBRset('InSolv', 'none', 'x_true', x_true, 'Iter', maxit, 'ReOrth','on');
[x_lsqr, output_lsqr] = HyBR(A,bn,[],input);

figure(2), plot(output_lsqr.Enrm,'ko','LineWidth',2), hold on
figure(3), subplot(2,3,1), imshow(x_lsqr,[]), title('lsqr')

%% Standard HyBR with Tikhonov and dp
input = HyBRset('InSolv', 'Tikhonov', 'x_true', x_true,'Iter', 50,'RegPar','dp', 'ReOrth','on','nLevel', sigma);
[xhybr_tik, output_hybr_tik] = HyBR(A,bn,[],input);

figure(2), plot(output_hybr_tik.Enrm,'rs','LineWidth',2), hold on
figure(3), subplot(2,3,2), imshow(xhybr_tik,[]), title('HyBR')

%% HyBR recycle using TSVD for compression
% (1) Set trunction options and matrices
trunc_options.nOuter = 11; % number of outer iterations
trunc_options.nInner = 50; % maximum storage of solution vector space
trunc_options.max_mm = 30; % maminimum number of vectors to save at compression
trunc_options.compress = 'SVD'; 
trunc_mats = [];

% (2) Set options for hybrid method
input = HyBRset('InSolv', 'Tikhonov', 'x_true', x_true,'Iter', 50,'RegPar', 'dp','nLevel', sigma); 

% (3) Run HyBR recycle
[xhybr_recycle, toutput_recycle, trunc_mats] = HyBRrecycle(A,bn,[],input, trunc_options, trunc_mats);

figure(2), plot(toutput_recycle.Enrm,'b','LineWidth',2)
legend('LSQR','HyBR','HyBR-recycle-dp-tsvd')
figure(3), subplot(2,3,3), imshow(reshape(xhybr_recycle,size(x_true)),[]), title('HyBR-recycle-dp-tsvd')

X = [x_lsqr(:), xhybr_tik(:), xhybr_recycle(:)];
[E,cmin,cmax,ns,nt] = geterrorimages(X, x_true, gcf);

%% Function to construct a Gaussian blur
function [PSF, center] = Gauss( gamma , N )
%
%     [PSF, center] = Gauss( gamma , N )
%
% This function constructs the Gaussian blur PSF, and the center of the PSF
%
%   Input:
%        gamma - parameters defining the blur
%            N - size fo the PSF
%
%    Output:
%          PSF - Gaussian point spread function
%      center  - center of the PSF
%
%  J. Chung 1/24/09
%

m = N;
n = N;

%
% Set up grid points to evaluate the Gaussian function.
%
x = -fix(n/2):ceil(n/2)-1;
y = -fix(m/2):ceil(m/2)-1;
[X,Y] = meshgrid(x,y);

%
% Compute the Gaussian PSF
%
if length(gamma) == 1
  s1 = gamma(1);
  PSF = exp( -(x.^2)/(2*s1^2));
  PSFsum = sum(PSF(:));
elseif length(gamma) == 2
  s1 = gamma(1); s2 = gamma(2);
  PSF = exp( -(X.^2)/(2*s1^2) - (Y.^2)/(2*s2^2) );
  PSFsum = sum(PSF(:));
  
elseif length(gamma)==3
  s1 = gamma(1); s2 = gamma(2); s3 = gamma(3);
  num = -((X.^2)*(s1^2) + (Y.^2)*(s2^2) - 2*(s3^2)*(X.*Y));
  den = 2*(s1^2 * s2^2 - s3^4);
  PSF = exp( num / den );
  PSFsum = sum(PSF(:));
  
end
  
PSF = PSF/PSFsum;
center = [m/2+1,n/2+1];
end
