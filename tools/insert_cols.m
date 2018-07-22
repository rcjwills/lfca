function [x_aug] = insert_cols(x, icol_not_NaN, icol_NaN)
%INSERT_COLS    Insert columns of NaNs into a data matrix.
%
%    X_AUG = INSERT_COLS(X, ICOL_NOT_NaN, ICOL_NaN) inserts columns
%    of NaNs into the data matrix X, such that in the augmented
%    data matrix X_AUG the columns ICOL_NOT_NaN are taken from the
%    original data matrix X and the columns ICOL_NaN consist of
%    NaNs. 
%
%    The user must make sure that ICOL_NOT_NaN and ICOL_NaN have
%    no elements in common.

  narginchk(3,3)            % check number of input arguments
  
  % initialize augmented data matrix
  nins	  = length(icol_NaN);
  nrows   = size(x, 1);
  x_aug   = NaN .* ones(nrows, size(x, 2) + nins);
  
  x_aug(:, icol_not_NaN) = x;


