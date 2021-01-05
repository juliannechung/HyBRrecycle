% DEMO_Tomo.m
%
% This script sets up and computes reconstructions for a streaming 
%       tomography reconstruction problem with two datasets. We use 
%       hybrid projection methods with recycling as described in:
%
%      Chung, de Sturler, and Jiang. "Hybrid Projection Methods with 
%           Recycling for Inverse Problems". SISC, 2020.
%
% Note: This code requires IRTools: https://github.com/jnagy1/IRtools
%
% J. Chung, E. de Sturler, J. Jiang 12/2020

%% Problem Setup
% n = 1024;   % dimension of the image
n = 256;   % dimension of the image
nangles = 180; anglemax = 180;
theta = linspace(0,anglemax,nangles+1);
theta = theta(1:end-1); % all projection angles

%% First problem with projections from 0-89 degrees
theta1 = theta(1:90); % First set of angles
options.angles = theta1;
[A, b, x_true, ProbInfo] = PRtomo(n, options);
bsize = ProbInfo.bSize;

% Add noise to the data
nlevel= 0.02;
[N,sigma] = WhiteNoise(b(:),nlevel, 0); % seed is 0
bn = b(:) + N;
x_true = reshape(x_true,n,n);

%% Compute and save vectors for recycling
nsave = 10;
nW = 51;
trunc_options.nOuter = 1; % number of outer iterations
trunc_options.nInner = nW; % maximum storage of solution vector space
trunc_options.max_mm = nsave; % maminimum number of vectors to save at compression
trunc_options.compress = 'SVD'; 
trunc_mats = [];

input = HyBRset('InSolv', 'Tikhonov', 'x_true', x_true(:),'Iter', trunc_options.nInner,'Reorth','on','RegPar', 'dp','nLevel', sigma);

[x_0, ~, trunc_mats] = HyBRrecycle(A,bn,[],input, trunc_options, trunc_mats);
x_0 = reshape(x_0,n,n);
W = trunc_mats.W;

%% Second problem with projections from 90-179 degrees
theta2 = theta(91:end);
options.angles = theta2;
[A2, b2, ~, ProbInfo] = PRtomo(n, options);

% Add noise to the data
[N,sigma] = WhiteNoise(b2(:), nlevel, 100); % seed is 100
bn2 = b2(:) + N;

%% Display problem
figure, 
subplot(1,3,1), imshow(x_true,[]), title('Desired image')
subplot(1,3,2), imshow(reshape(bn,bsize),[]), title('Sinogram for problem 1')
subplot(1,3,3), imshow(reshape(bn2,bsize),[]), title('Sinogram for problem 2')

%% HyBR with recycling
% (1) Set trunction options
trunc_options.nOuter = 3; % number of outer iterations
trunc_options.nInner = 50; % maximum storage of solution vector space
trunc_options.max_mm = 10; % minimum number of vectors to save at compression
trunc_options.compress = 'SVD'; 
trunc_mats.Y = [];
trunc_mats.R = [];
trunc_mats.x = [];
trunc_mats.W = W;

% (2) Set hybrid options
input = HyBRset('InSolv', 'Tikhonov', 'x_true', x_true(:),'Iter', trunc_options.nInner,'RegPar', 'dp','nLevel', sigma);

tic
[xthybr_opt, toutput_opt, trunc_mats] = HyBRrecycle(A2, bn2(:),[],input, trunc_options, trunc_mats);
t1 = toc;
fprintf('Time for HyBR-recycle = %.4f\n', t1)

figure(100), set(gcf, 'Position',  [200, 100, 500, 300])
err_all = toutput_opt.Enrm;
plot(err_all,'-rx'), hold on
fontSize = 12;
maxit = length(err_all);
xlim([1 maxit]);
xlabel('iteration','fontsize',20)
ylabel('relative error','fontsize',20)

xthybr_opt = reshape(xthybr_opt,n,n);
figure(119), set(gcf, 'Position',  [200, 100, 600, 200])
subplot(2,4,1), imshow(xthybr_opt,[]), title('HyBR-recycle-dp-svd')

%% (2) Compare to starting from scratch on the second dataset

tic
input = HyBRset('InSolv', 'Tikhonov', 'x_true', x_true(:),'Iter', maxit, 'RegPar', 'dp', 'nLevel', sigma);
[x_2, output_2] = HyBR(A2, bn2(:), [],input); % From scratch
t2 = toc;
fprintf('Time for HyBR on 2nd dataset = %.4f\n', t2)

x_2 = reshape(x_2,n,n);
Enrm_2 = output_2.Enrm;
figure(100), hold on, plot(Enrm_2,'k--', 'LineWidth',2)
figure(119), subplot(2,4,2), imshow(x_2,[]), title('HyBR on 2nd dataset')


%% (3) Compare to HyBR on ALL data
input = HyBRset('InSolv', 'Tikhonov', 'x_true', x_true(:),'Iter', maxit, 'RegPar', 'dp', 'nLevel', sigma);
tic
[x_3, output_3] = HyBR([A; A2], [bn(:); bn2(:)], [], input);
t3 = toc;
fprintf('Time for HyBR on ALL data = %.4f\n', t3)

x_3 = reshape(x_3,n,n);
Enrm_3 = output_3.Enrm;
figure(100), hold on, plot(Enrm_3,'m', 'LineWidth',2), 
figure(119), subplot(2,4,3), imshow(x_3,[]), title('HyBR for all data')


%% (4) Add 2 solutions
x_4 = .5*( x_0 + x_2 );
Enrm_4 = norm(x_4(:)-x_true(:))/norm(x_true(:))*ones(maxit,1);
figure(100), hold on, plot(Enrm_4,'b:','LineWidth',2), 
legend('HyBR-recycle-dp-svd', 'HyBR on 2nd dataset','HyBR all data','Average of Images')

figure(119), subplot(2,4,4), imshow(x_4,[]), title('Average of images')

%% Display error images
figure(119),
[E,cmin,cmax,ns,nt] = geterrorimages([xthybr_opt(:), x_2(:), x_3(:), x_4(:)], x_true);
