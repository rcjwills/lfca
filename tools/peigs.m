function [V, d, r] = peigs(A, rmax)
%PEIGS   Finds positive eigenvalues and corresponding eigenvectors. 
%  
%    [V, d, r] = PEIGS(A, rmax) determines the positive eigenvalues d
%    and corresponding eigenvectors V of a matrix A. The input
%    parameter rmax is an upper bound on the number of positive
%    eigenvalues of A, and the output r is an estimate of the actual
%    number of positive eigenvalues of A. The eigenvalues are returned
%    as the vector d(1:r). The eigenvectors are returned as the
%    columns of the matrix V(:, 1:r).
%
%    If rmax < min(size(A))/10, PEIGS calls the function EIGS to compute 
%    the first rmax eigenpairs of A by Arnoldi iterations. Eigenpairs 
%    corresponding eigenvalues that are nearly zero or less than zero are
%    subsequently discarded.
%
%    If rmax >= min(size(A))/10, PEIGS calls the function EIG to compute 
%    the full eigendecomposition of A (which is faster than computing 
%    only the first rmax eigenpairs unless A is strongly rank-deficient).
%
%    See also: EIGS, EIG.

  narginchk(2,2)                    % check number of input arguments 
  [m, n]  = size(A);
 
  if rmax > min(m,n)
    rmax  = min(m,n);                             % rank cannot exceed size of A 
  end
     
  if rmax < min(m, n)/10
    % get first rmax eigenvectors of A by Arnoldi iterations
    warning off
    [V, d]     = eigs(A, rmax, 'lm'); 
    warning on
  else
    % compute full eigendecomposition
    [V, d]     = eig(A);
  end
  
  % output format of eig/eigs differs (and it differs between Matlab
  % versions): fix so that d is always a vector
  if prod(size(d)) > max(size(d))
    d          = diag(d);                         % ensure d is vector
  end

  % ensure that eigenvalues are monotonically decreasing
  [d, I]       = sort(d, 'descend');
  V            = V(:, I);
  
  % estimate number of positive eigenvalues of A
  d_min        = max(d) * max(m,n) * eps; 
  r            = sum(d > d_min);
				 
  % discard eigenpairs with eigenvalues that are close to or less than zero
  d            = d(1:r);
  V            = V(:, 1:r);
  d            = d(:);				  % ensure d is column vector
