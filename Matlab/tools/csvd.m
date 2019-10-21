function [U, s, V] = csvd(A, flag)
%CSVD    Compact singular value decomposition.
%  
% CSVD(A) returns the singular values of the matrix A. If A is of
% size [m,n], min(m, n) singular values are computed.
%
% [U, s, V] = CSVD(A) returns the compact form of the singular
% value decomposition A = U*diag(s)*V'  of A. The matrix of left
% singular vectors U is of size [m, min(m,n)], the matrix of right
% singular vectors V is of size [n, min(m, n)], and the vector of
% singular values s is of length min(m,n).
%
% If a second input argument is present, as in CSVD(A, 'full'), the
% full matrices U and V are returned.

% Adapted from CSVD in Per Christian Hansen's Regularization Tools Toolbox.

  if nargin==1                    % economoy size SVD
    if nargout == 1           
      U = svd(full(A));           % return only singular values
    else
      [m,n] = size(A);            % return complete SVD
      if m >= n
	[U, s, V] = svd(full(A),0); 
      else
	[V, s, U] = svd(full(A)',0); 
      end
      s = diag(s);
    end
  else                            % full size SVD
    if nargout == 1
      U = svd(full(A));           % return only singular values
    else
      [U, s, V] = svd(full(A));   % return complete SVD
      s = diag(s);   
    end
  end

