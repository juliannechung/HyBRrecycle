%% Contents.m
%   This folder contains MATLAB code to accompany the paper:
%
%      Chung, de Sturler, and Jiang. "Hybrid Projection Methods with 
%           Recycling for Inverse Problems". SISC, 2021.
%
%   The DEMO codes require the following packages:
%       RestoreTools by James Nagy
%             http://www.mathcs.emory.edu/~nagy/RestoreTools/
%       IRTools by James Nagy, Silvia Gazzola, and Per Christian Hansen
%             https://github.com/jnagy1/IRtools
%       SpaRSA by Stephen Wright, Robert Nowak, and Mario Figueiredo
%             https://www.lx.it.pt/~mtf/SpaRSA/
%
%   First run startup_Recycle.m for setting paths
%
% Chung, de Sturler, and Jiang (2020)

%% DEMOs on how to use HyBRrecycle
%
%   DEMO_Deblurring.m     Sets up and runs a 2D image deblurring problem
%                           corresponding to the results in Section 4.1 
%                           of the paper
%
%   DEMO_Tomo.m           Sets up and runs a 2D streaming tomgraphy problem
%                           with two sets of data.  
%                           This example with n = 1024 corresponds to 
%                           Case 1 of Section 4.2.1 in the paper.
%
%% Supporting MATLAB functions
%
% HyBRrecycle.m           Main code to run HyBR with recycling
%                         Syntax: [x_out, output, trunc_mats] =
%                         HyBRrecycle(A, b, P, options, trunc_options, trunc_mats
% 
% recyclingGKB.m          Code to perform one step of the recycling GKB
%                         process
%
% compression.m           Code to perform basis compression. Here we
%                         provide four compression approaches: 
%                         SVD, solution, sparse and RBD. 
%
%
% To obtain full functionality of these codes, it is recommended to also
% install 
%   (3) SpaRSA by Stephen Wright, Robert Nowak, and Mario Figueiredo
%              https://www.lx.it.pt/~mtf/SpaRSA/
%   (4) RBD by Yanlai Chen 
%              http://yanlaichen.reawritingmath.com