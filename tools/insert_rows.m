function [x_aug] = insert_rows(x, irow_not_NaN, irow_NaN)
%INSERT_ROWS    Insert rows of NaNs into a data matrix.
%
%    X_AUG = INSERT_ROWS(X, IROW_NOT_NaN, IROW_NaN) inserts rows
%    of NaNs into the data matrix X, such that in the augmented
%    data matrix X_AUG the rows IROW_NOT_NaN are taken from the
%    original data matrix X and the columns IROW_NaN consist of
%    NaNs. 
%
%    The user must make sure that IROW_NOT_NaN and IROW_NaN have
%    no elements in common.

  narginchk(3,3)           % check number of input arguments
  
  % initialize augmented data matrix
  nins	  = length(irow_NaN);
  ncols   = size(x, 2);
  x_aug   = NaN .* ones(size(x, 1) + nins, ncols);
  
  x_aug(irow_not_NaN, :) = x;


