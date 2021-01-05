function [U, B, V, H] = recyclingGKB(A, U, B, V, H, P_le, P_ri, W, Y, options)
%
%     [U, B, V, H] = recyclingGKB(A, U, B, V, H, P_le, P_ri, options)
%
%  Perform one step of recycling GKB without reorthogonalization and 
%   without preconditioner.
%
% Input:
%          A - matrix
%       U, V - accumulation of vectors
%          H - ADD definitions
%          B - bidiagonal matrix
% P_le, P_ri - inputs ignored
%       W, Y - Additional vectors for recycling
%    options - structure from HyBR (see HyBRset), ignored here
%
% Output:
%       U, V - updated "orthogonal" matrix
%          B - updated bidiagonal matrix
%          H - updated upper triangular portion of the matrix
%
%  Reference:
%   Chung, de Sturler, and Jiang. "Hybrid Projection Methods with
%           Recycling for Inverse Problems". SISC, 2020.
%
% J. Chung, E. de Sturler, and J. Jiang, 2020

k = size(B,2)+1;

if k == 1
  Au = A'*U(:,k);
  WtAu = W'*Au;
  v = Au - W*WtAu;
else
  Au = A'*U(:,k);
  WtAu = W'*Au; 
  v = (Au - W*WtAu) - B(k, k-1)*V(:,k-1);
end
alpha = norm(v);
v = v / alpha;
Av = A*v; YtAv = Y'*Av; H(:,k) = YtAv;
u = Av - Y*YtAv - alpha*U(:,k);
beta = norm(u);
u = u / beta;
U = [U, u];

V = [V, v];
B = [B, [zeros(k-1,1); alpha]; [zeros(1,k-1), beta]];
