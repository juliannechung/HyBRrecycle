function [U, B, V] = PLBD(A, U, B, V, P_le, P_ri, options)
%
%   [U, B, V] = PLBD(A, U, B, V, P_le, P_ri, options)
%
%  Perform one step of Lanczos bidiagonalization with or without
%  reorthogonalization, WITH preconditioner here.
%
% Input:
%        A - matrix
%     U, V - "orthogonal" matrix for bidiagonal factorization update
%        B - bidiagonal matrix
%     P_le - left preconditioner
%     P_ri - right preconditioner
%  options - structure from HyBR (see HyBRset)
%
% Output:
%     U, V - updated "orthogonal" matrix
%        B - updated bidiagonal matrix
%
%  Refs: 
%   [1] Paige and Saunders, "LSQR an algorithm for sparse linear
%       equations an sparse least squares", ACM Trans. Math Software,
%       8 (1982), pp. 43-71.
%   [2] Bjorck, Grimme and Van Dooren, "An implicit shift bidiagonalization
%       algorithm for ill-posed systems", BIT 34 (11994), pp. 520-534.
%
%   J.Chung and J. Nagy 3/2007

% Determine if we need to do reorthogonalization or not.
reorth = strcmp(HyBRget(options,'Reorth'), {'on'});

precleft = 0; precright = 0;
if ~isempty(P_le)
  precleft = 1;
end
if ~isempty(P_ri)
  precright = 1;
end

m = size(U,1);
k = size(B,2)+1;

if reorth % Need reorthogonalization
  v = U(:,k);
  if precleft
    v = (P_le'\v);
  end
  v = A'*v;
  if precright
    v = P_ri' \ v;
  end

  if k > 1
    v = v - B(k, k-1)*V(:,k-1); 
    for j = 1:k-1
      v = v - (V(:,j)'*v)*V(:,j);
    end
  end
  alpha = norm(v);
  v = v / alpha;
  
  if precright
    v = P_ri \ v;
  end  
  v = A*v;
  if precleft
    v = P_le \ v;
  end  
  
  u = v - alpha*U(:,k);
 
  for j = 1:k
    u = u - (U(:,j)'*u)*U(:,j);
  end

  beta = norm(u);
  u = u / beta;
  U = [U, u];
else
  v = U(:);
  if precleft
    v = (P_le' \ v);
  end
  v = A'*v;
  if precright
    v = P_ri' \ v;
  end
  if k > 1
    v = v - B(k, k-1)*V(:,k-1);
  end
  alpha = norm(v);
  v = v / alpha;
  u = v;
  if precright
    u = P_ri \ u;
  end  
  u = A*u;
  if precleft
    u = P_le \ u;
  end  
  
  u = u - alpha*U(:);

  beta = norm(u);
  u = u / beta;
  U = u(:);  
end

V = [V, v];
B = [B, [zeros(k-1,1); alpha]; [zeros(1,k-1), beta]];

