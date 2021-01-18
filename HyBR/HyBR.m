function [x_out, output] = HyBR(A, b, P, options)
%
% [x_out, output] = HyBR(A, b, P, options)
%
% HyBR is a Hybrid Bidiagonalization Regularization method used for
% solving large-scale, ill-posed inverse problems of the form:
%               b = A*x + noise
% The method combines an iterative Lanczos Bidiagonalization (LBD) Method
% with a SVD-based regularization method to stabilize the semiconvergence
% behavior that is characteristic of many ill-posed problems.
%
% Inputs:
%                 A : either (a) a full or sparse matrix
%                            (b) a matrix object that performs matrix*vector
%                                 and matrix'*vector operations
%                b : rhs vector
%                P : left preconditioner, P_left, OR
%                  : cell containing left and right preconditioner (optional)
%                          {P_left, P_right}
% options : structure with the following fields (optional)
%         InSolv - solver for the inner problem: [none | TSVD | {Tikhonov}]
%         RegPar - a value or method to find the regularization parameter:
%                       [non-negative scalar | GCV | {WGCV} | {NWGCV}  |optimal]
%                   Note: 'optimal' requires x_true
%          Omega - if RegPar is 'WGCV', then omega must be
%                       [non-negative scalar | {adapt}]
%           Iter - maximum number of Lanczos iterations:
%                       [ positive integer | {min(m,n,100)} ]
%         Reorth - reorthogonalize Lanczos subspaces: [on | {off}]
%         x_true - True solution : [ array | {off} ]
%                Returns error norms with respect to x_true at each iteration
%                and is used to compute 'optimal' regularization parameters
%         BegReg - Begin regularization after this iteration:
%                   [ positive integer | {2} ]
%             Vx - extra space needed for finding optimal reg. parameters
%        FlatTol - Tolerance for detecting flatness in the GCV curve as a
%                    stopping criteria
%                   [ non-negative scalar | {10^-6}]
%         MinTol - Window of iterations for detecting a minimum of the GCV curve
%                    as a stopping criteria
%                   [ positive integer | {3}]
%         ResTol - Residual tolerance for stopping the LBD iterations,
%                    similar to the stopping criteria from [1]: [atol, btol]
%                   [non-negative scalar  | {[10^-6, 10^-6]}]
%
%       Note: options is a structure created using the function 'HyBRset'
%               (see 'HyBRset' for more details)
%
% Outputs:
%      x_out : computed solution
%     output : structure with the following fields:
%      iterations - stopping iteration (options.Iter | GCV-determined)
%         GCVstop - GCV curve used to find stopping iteration
%            Enrm - relative error norms (requires x_true)
%            Rnrm - relative residual norms
%            Xnrm - relative solution norms
%             U,V - Lanczos basis vectors
%               B - bidiagonal matrix from LBD
%            flag - a flag that describes the output/stopping condition:
%                       1 - flat GCV curve
%                       2 - min of GCV curve (within window of MinTol its)
%                       3 - performed max number of iterations
%                       4 - achieved residual tolerance
%           alpha - regularization parameter at (output.iterations) its
%           Alpha - vector of all regularization parameters computed
%
% References:
%   [1] Paige and Saunders, "LSQR an algorithm for sparse linear
%       equations an sparse least squares", ACM Trans. Math Software,
%       8 (1982), pp. 43-71.
%   [2] Bjorck, Grimme and Van Dooren, "An implicit shift bidiagonalization
%       algorithm for ill-posed systems", BIT 34 (11994), pp. 520-534.
%   [3] Chung, Nagy and O'Leary, "A Weighted-GCV Method for Lanczos-Hybrid
%       Regularization", ETNA 28 (2008), pp. 149-167.
%   [4] R. A. Renaut, S. Vatankhah and V. E. Ardestani, (2016), 
%       Hybrid and iteratively reweighted regularization by unbiased predictive risk and weighted GCV for projected systems. 
%       SIAM J. Sci. Comput. 39-2 (2017), pp. B221-B243. 
%
% J.Chung and J. Nagy 3/2007
% Updated by Jiahua and Julianne  2020
%

%% Initialization
defaultopt = struct('InSolv','tikhonov','RegPar','wgcv','Omega',...
    'adapt', 'Iter', [] , 'Reorth', 'off', 'x_true', 'off', 'BegReg', 2,...
    'Vx' , [], 'FlatTol', 10^-6, 'MinTol', 4, 'ResTol', [10^-6, 10^-6],'nLevel', 'est');

% If input is 'defaults,' return the default options in x_out
if nargin==1 && nargout <= 1 && isequal(A,'defaults')
    x_out = defaultopt;
    return;
end

% Check for acceptable number of input arguments
if nargin < 2
    error('HyBR: Not Enough Inputs')
elseif nargin < 3
    P = {[],[]}; options = [];
elseif nargin < 4
    options = [];
end
if isempty(options)
    options = defaultopt;
end

if iscell(P)
    P_le = P{1};
    P_ri = P{2};
else
    P_le = P;
    P_ri = [];
end

% Get options:
[m,n] = size(A);
defaultopt.Iter = min([m, n, 100]);
options = HyBRset(defaultopt, options);

solver = HyBRget(options,'InSolv',[],'fast');
regpar = HyBRget(options,'RegPar',[],'fast');
omega = HyBRget(options,'Omega',[],'fast');
maxiter = HyBRget(options,'Iter',[],'fast');
x_true = HyBRget(options,'x_true',[],'fast');
regstart = HyBRget(options,'BegReg',[],'fast');
degflat = HyBRget(options,'FlatTol',[],'fast');
mintol = HyBRget(options,'MinTol',[],'fast');
restol = HyBRget(options,'ResTol',[],'fast');
nLevel = HyBRget(options,'nLevel',[],'fast');

adaptWGCV = strcmp(regpar, {'wgcv'}) && strcmp(omega, {'adapt'});
notrue = strcmp(x_true,{'off'});

if (strcmpi(regpar,'dp') || strcmpi(regpar,'upre')) && strcmpi(nLevel,'est')
  % Estimate the noise level
  if size(b,2) == 1
    [cA, cD] = dwt(b,'db1');
    nLevel = median(abs(cD(:)))/.67;
    options = HyBR_lsmrset(options, 'nLevel', nLevel);
  else
    [cA2,cH2,cV2,cD2] = dwt2(b,'db1');
    nLevel = median(abs(cD2(:)))/.67;
    options = HyBR_lsmrset(options, 'nLevel', nLevel);
  end
end

%--------------------------------------------
%  The following is needed for RestoreTools:
%
if isa(A, 'psfMatrix')
    bSize = size(b);
    b = b(:);
    A.imsize = bSize;
    if ~notrue
        xTrueSize = size(x_true);
        x_true = x_true(:);
    end
end
%
%  End of new stuff needed for RestoreTools
%--------------------------------------------

if ~notrue
    nrmtrue = norm(x_true(:));
end

% Set-up output parameters:
outputparams = nargout>1;
if outputparams
    output.iterations = maxiter;
    output.GCVstop = [];
    output.Enrm = ones(maxiter,1);
    output.Rnrm = ones(maxiter,1);
    output.Xnrm = ones(maxiter,1);
    output.U = [];
    output.V = [];
    output.B = [];
    output.flag = 3;
    output.alpha = 0;
end

% Test for a left preconditioner and define solver:
if isempty(P_le)
    beta = norm(b);
    U = (1 / beta)*b;
    handle = @LBD;
    if ~isempty(P_ri)
        handle = @PLBD;
    end
else
    U = P_le\b;
    beta = norm(U); U = U / beta;
    handle = @PLBD;
end

switch solver
    case 'tsvd'
        solverhandle = @TSVDsolver;
    case 'tikhonov'
        solverhandle = @Tikhonovsolver;
end
%% Main Code Begins Here
B = []; V = []; GCV = []; Omega= []; x_out = [];
insolve = 'none'; terminate = 1; warning = 0; norma = 0; normr = beta;
h = waitbar(0, 'Beginning iterations: please wait ...');

for i = 1:maxiter+1 %Iteration (i=1) is just an initialization
    [U, B, V] = feval(handle, A, U, B, V, P_le, P_ri, options,1);
    vector = (beta*eye(size(B,2)+1,1));
    if ~notrue
%         options.Vx = V'*x_true;
      options.Vx = V;
    end
    
    if i >= 2 %Begin Lanczos iterations
        if i >= regstart %Begin to regularize projected problem
            insolve = solver;
        end
        switch insolve
            case {'tsvd', 'tikhonov'}
                
                [Ub, Sb, Vb] = svd(B);
                
                if adaptWGCV %Use the adaptive, weighted GCV method
                    Omega(i-1) = min(1, findomega(Ub'*vector, diag(Sb), insolve));
                    options.Omega = mean(Omega);
                end
                
                % Solve the projected problem with Tikhonov or TSVD
%                 [f, alpha] = feval(solverhandle, Ub, diag(Sb), Vb, vector, options);
              [f, alpha] = feval(solverhandle, Ub, diag(Sb), Vb, vector, options, B, i, V, m);
                Alpha(i-1) = alpha;
                
                % Compute the GCV value used to find the stopping criteria
                GCV(i-1) = GCVstopfun(alpha, Ub(1,:)', diag(Sb), beta, m, n, insolve);
                
                % Determine if GCV wants us to stop
                if i > 2 && terminate
                    %%-------- If GCV curve is flat, we stop -----------------------
                    if abs((GCV(i-1)-GCV(i-2)))/GCV(regstart-1) < degflat
                        x_out = V*f; % Return the solution at (i-1)st iteration
                        % Test for a right preconditioner:
                        if ~isempty(P_ri)
                            x_out = P_ri \ x_out;
                        end
                        if notrue %Set all the output parameters and return
                            if outputparams
                                output.U = U;
                                output.V = V;
                                output.B = B;
                                output.GCVstop = GCV(:);
                                output.iterations = i-1;
                                output.flag = 1;
                                output.alpha = alpha; % Reg Parameter at the (i-1)st iteration
                                output.Alpha = Alpha(1:i-1); % Reg Parameters
                            end
                            close(h)
                            %--------------------------------------------
                            %  The following is needed for RestoreTools:
                            %
                            if isa(A, 'psfMatrix')
                                x_out = reshape(x_out, bSize);
                            end
                            %
                            %  End of new stuff needed for RestoreTools
                            %--------------------------------------------
                            return;
                        else % Flat GCV curve means stop, but continue since have x_true
                            if outputparams
                                output.iterations = i-1; % GCV says stop at (i-1)st iteration
                                output.flag = 1;
                                output.alpha = alpha; % Reg Parameter at the (i-1)st iteration
                            end
                        end
                        terminate = 0; % Solution is already found!
                        
                        %%--- Have warning : Avoid bumps in the GCV curve by using a
                        %    window of (mintol+1) iterations --------------------
                    elseif warning && length(GCV) > iterations_save + mintol %Passed window
                        if GCV(iterations_save) < GCV(iterations_save+1:end)
                            % We should have stopped at iterations_save.
                            x_out = x_save;
                            % Test for a right preconditioner:
                            if ~isempty(P_ri)
                                x_out = P_ri \ x_out;
                            end
                            if notrue %Set all the output parameters and return
                                if outputparams
                                    output.U = U;
                                    output.V = V;
                                    output.B = B;
                                    output.GCVstop = GCV(:);
                                    output.iterations = iterations_save;
                                    output.flag = 2;
                                    output.alpha = alpha_save;
                                    output.Alpha = Alpha(1:iterations_save); % Reg Parameters
                                end
                                close(h)
                                %--------------------------------------------
                                %  The following is needed for RestoreTools:
                                %
                                if isa(A, 'psfMatrix')
                                    x_out = reshape(x_out, bSize);
                                end
                                %
                                %  End of new stuff needed for RestoreTools
                                %--------------------------------------------
                                return;
                            else % GCV says stop at iterations_save, but continue since have x_true
                                if outputparams
                                    output.iterations = iterations_save;
                                    output.flag = 2;
                                    output.alpha = alpha_save;
                                end
                            end
                            terminate = 0; % Solution is already found!
                            
                        else % It was just a bump... keep going
                            warning = 0;
                            x_out = [];
                            iterations_save = maxiter;
                            alpha_save = 0;
                        end
                        
                        %% ----- No warning yet: Check GCV function---------------------
                    elseif ~warning
                        if GCV(i-2) < GCV(i-1) %Potential minimum reached.
                            warning = 1;
                            % Save data just in case.
                            x_save = V*f;
                            iterations_save = i-1;
                            alpha_save = alpha;
                        end
                    end
                end
                
            case 'none'
                f = B \ vector;
                alpha = 0;
            case 'sparse'
                alpha = .001; % Regularization parameter for sparsity
                [f]= SpaRSA(vector,B,alpha);
            otherwise
                error('HyBR error: No inner solver!')
        end
        x = V*f;
        % Test for a right preconditioner:
        if ~isempty(P_ri)
            x = P_ri \ x;
        end
        r = b(:) - A*x(:);
        normr = norm(r(:));
        if outputparams
            if ~notrue
                output.Enrm(i-1,1) = norm(x(:)-x_true(:))/nrmtrue;
            end
            output.Rnrm(i-1,1) = normr;
            output.Xnrm(i-1,1) = norm(x(:));
        end
        
        norma = norm([norma B(i,i) B(i+1,i)]);
        normar = norm(A'*r);
        normx = norm(x(:));
        
        if normr <= restol(1)*beta+restol(2)*norma*normx || normar/(norma*normr) <= restol(2) && terminate
            if notrue %Set all the output parameters and return
                if outputparams
                    output.U = U;
                    output.V = V;
                    output.B = B;
                    output.GCVstop = GCV(:);
                    output.iterations = i-1;
                    output.flag = 4;
                    output.alpha = alpha; % Reg Parameter at the (i-1)st iteration
                    output.Alpha = Alpha(1:i-1); % Reg Parameters
                end
                close(h)
                %--------------------------------------------
                %  The following is needed for RestoreTools:
                %
                if isa(A, 'psfMatrix')
                    x_out = reshape(x, bSize);
                else
                    x_out = x;
                end
                %
                %  End of new stuff needed for RestoreTools
                %--------------------------------------------
                return
            else % Residual says stop, but continue since have x_true
                if outputparams
                    output.iterations = i-1;
                    output.flag = 4;
                    output.alpha = alpha;
                end
            end
            terminate = 0; % Solution is already found!
        end
    end
    waitbar(i/(maxiter+1), h)
end
close(h)

if isempty(x_out) % GCV did not stop the process, so we reached max. iterations
    x_out = x;
    % Test for a right preconditioner:
    if ~isempty(P_ri)
        x_out = P_ri \ x_out;
    end
end

%--------------------------------------------
%  The following is needed for RestoreTools:
%
if isa(A, 'psfMatrix')
    x_out = reshape(x, bSize);
end
%
%  End of new stuff needed for RestoreTools
%--------------------------------------------

if outputparams
    output.U = U;
    output.V = V;
    output.B = B;
    output.GCVstop = GCV(:);
    if output.alpha == 0
        output.alpha = alpha;
    end
    %  output.Alpha = Alpha;
end

%% -----------------------SUBFUNCTION---------------------------------------
function omega = findomega(bhat, s, insolv)
%
%   omega = findomega(bhat, s, insolv)
%
%  This function computes a value for the omega parameter.
%
%  The method: Assume the 'optimal' regularization parameter to be the
%  smallest singular value.  Then we take the derivative of the GCV
%  function with respect to alpha, evaluate it at alpha_opt, set the
%  derivative equal to zero and then solve for omega.
%
%  Input:   bhat -  vector U'*b, where U = left singular vectors
%              s -  vector containing the singular values
%         insolv -  inner solver method for HyBR
%
%  Output:     omega - computed value for the omega parameter.

%
%   First assume the 'optimal' regularization parameter to be the smallest
%   singular value.
%

%
% Compute the needed elements for the function.
%
m = length(bhat);
n = length(s);
switch insolv
    case 'tsvd'
        k_opt = n;
        omega = (m*bhat(k_opt)^2) / (k_opt*bhat(k_opt)^2 + 2*bhat(k_opt+1)^2);
        %  omega = ((m/2)*(bhat(k_opt)^2 + bhat(k_opt+1)^2)) / ((k_opt/2)*(bhat(k_opt)^2 +bhat(k_opt+1)^2) + 2*bhat(k_opt+1)^2)
    case 'tikhonov'
        t0 = sum(abs(bhat(n+1:m)).^2);
        alpha = s(end);
        s2 = abs(s) .^ 2;
        alpha2 = alpha^2;
        
        tt = 1 ./ (s2 + alpha2);
        
        t1 = sum(s2 .* tt);
        t2 = abs(bhat(1:n).*alpha.*s) .^2;
        t3 = sum(t2 .* abs((tt.^3)));
        
        t4 = sum((s.*tt) .^2);
        t5 = sum((abs(alpha2*bhat(1:n).*tt)).^2);
        
        v1 = abs(bhat(1:n).*s).^2;
        v2 = sum(v1.* abs((tt.^3)));
        
        %
        % Now compute omega.
        %
        omega = (m*alpha2*v2)/(t1*t3 + t4*(t5 + t0));
        
    otherwise
        error('Unknown solver');
end

%% ---------------SUBFUNCTION ---------------------------------------
function G = GCVstopfun(alpha, u, s, beta, m, n, insolv)
%
%  G = GCVstopfun(alpha, u, s, beta, n, insolv)
%  This function evaluates the GCV function G(i, alpha), that will be used
%     to determine a stopping iteration.
%
% Input:
%   alpha - regularization parameter at the kth iteration of HyBR
%       u - P_k^T e_1 where P_k contains the left singular vectors of B_k
%       s - singular values of bidiagonal matrix B_k
%    beta - norm of rhs b
%     m,n - size of the ORIGINAL problem (matrix A)
%  insolv - solver for the projected problem
%

k = length(s);
beta2 = beta^2;

switch insolv
    case 'tsvd'
        t2 = (abs(u(alpha+1:k+1))).^2;
        G = n*beta2*(sum(t2))/((m - alpha)^2);
    case 'tikhonov'
        s2 = abs(s) .^ 2;
        alpha2 = alpha^2;
        
        t1 = 1 ./ (s2 + alpha2);
        t2 = abs(alpha2*u(1:k) .* t1) .^2;
        t3 = s2 .* t1;
        
        num = beta2*(sum(t2) + abs(u(k+1))^2)/n;
        den = ( (m - sum(t3))/n )^2;
        G = num / den;
        
    otherwise
        error('Unknown solver');
end