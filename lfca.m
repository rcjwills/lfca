function [LFCs, LFPs, weights, r, pvar, pcs, EOF, N, pvar_slow, pvar_lfc, r_eofs, pvar_slow_eofs] = lfca(X, cutoff, truncation, scale, Covtot)

%% LFCA  Truncated Low-Frequency Component Analysis
%     [LFCs,LFPs,WEIGHTS,R,PVAR,PCS,EOF,N,PVAR_SLOW,PVAR_LFC,R_EOFS,PVAR_SLOW_EOFS] = lfca(X,CUTOFF,TRUNCATION,SCALE,COVTOT)
%     performs low-frequency component analysis (LFCA) on the data in 
%     matrix X based on a the ratio of low-pass filtered to unfiltered 
%     variance, with a low-pass filter defined by CUTOFF (in # timesteps).
%     Another type of filter can be substituted for the Lanczos lowpass
%     filter in the code below.

%% INPUT
%     X is a 2D data matrix with time variantions along the first dimension
%     and spatial variations along the second dimension
%
%     CUTOFF is the lowpass cutoff for the LFCA in number of timesteps
%
%     TRUNCATION is the number of principal components / EOFs to include in
%     the LFCA
%
%     SCALE (optional) a scale vector, which for geospatial data should be 
%     equal to the square root of grid cell area. The default value is one  
%     for all grid points.
%
%     COVTOT (optional) the covariance matrix associated with the data in
%     X. If not specified, COVTOT will be computed from X.

%% OUTPUT

%     WEIGHTS is a matrix containing the canonical weight vectors as
%     columns. LFPs is a matrix containing the dual vectors of the
%     canonical weight vectors as rows. These are the so-called low-
%     frequency patterns (LFPs). R is a vector measuring the ratio of 
%     low-pass filtered to total variance for each low-frequency component
%
%     PVAR is the percentage of total sample variation accounted for
%     by each of the EOFs. PCS is a matrix containing the principal
%     component time series as columns. EOF is a matrix containing the
%     EOFs, the principal component patterns, as rows. The scalar N
%     is the rank at which the PCA was truncated.
%
%     R is a vector containing the ratio of low-frequency to
%     total variance of each LFC.
%     
%     PVAR_SLOW is a vector of the low-frequency variance associated with 
%     each LFC as a fraction of the total low-frequency variance. Note that
%     the LFPs are not orthogonal, so these values need not add to the
%     total low-frequency variance in the first N principal components.
%
%     PVAR_LFC is a vector of the variance associated with 
%     each LFC as a fraction of the total variance. Note that
%     the LFPs are not orthogonal, so these values need not add to the
%     total variance in the first N principal components.
%
%     R_EOFS, PVAR_SLOW_EOFS, and PVAR are equivalent to R, PVAR_SLOW,
%     and PVAR_LFC respectively, but for the original EOFs.

  narginchk(3,5)          % check number of input arguments 
  if ndims(X) ~= 2,  error('Data matrix must be 2-D.'); end

  disp(sprintf('\nLFCA:'))

  [n,p]         = size(X);
  
  % center data 
  if any(any(isnan(X)))               % there are missing values in x
    Xm  = nanmean(X);
  else                                % no missing values
    Xm  = mean(X);
  end
  X    = X - repmat(Xm, n, 1);  
      
  %% compute covariance matrix
  if nargin < 5
    % compute sample covariance if covariance is not specified
    Covtot               = cov(X);
  end
  if any(size(Covtot) ~= [p, p])
    error('Covariance matrix must have same dimension as data.')
  end
  
  %% scale vector (e.g. square root of normalized grid-cell area)
  if nargin > 3
    scale       = scale(:)';
    if length(scale) ~= p
      error('Scale vector must have same dimension as data.')
    end
    Xs           = X .* repmat(scale,n,1);
  else
    scale       = ones(1,p);
    Xs          = X;
  end
  clear X
  
  %% eigendecomposition of covariance matrix
  Covtot      = repmat(scale',1,p) .* Covtot .* repmat(scale,p,1);
  [pcvec,evl,rest] = peigs(Covtot, min(n-1, p));
  trCovtot    = trace(Covtot);
  
  % percent of total sample variation accounted for by each EOF
  pvar          = evl./trCovtot .* 100;
  % principal component time series
  pcs           = Xs*pcvec;
  % return EOFs in original scaling as patterns (row vectors)
  EOF           = pcvec' ./ repmat(scale,rest,1);
  
  %% truncation of EOFs
  if truncation < 1 
      % using basic % variance criterion, where truncation gives the
      % fraction of variance to be included in the EOF truncation
      truncation = truncation*100;
      cum_pvar = cumsum(pvar);
      N = closest(cum_pvar,truncation);
      disp(sprintf('\tChosen truncation level: %3i', N))
  else
      if (truncation-round(truncation))~=0
          error('Truncation must be fraction of total variance included in EOF truncation or integer number of EOFs.')
      end
      % using specified truncation level
      N = truncation;
      disp(sprintf('\tUsing specified truncation level: %3i', N))
  end
  
  % this section can be modified to use a specific EOF truncation
  % criterion, right now the truncation number is specified as input
  % TRUNCATION
  
  %% Whitening transformation
  % multiplication factor for principal components in whitening
  % transformation (such that they have unit variance)
  f             = sqrt(evl(1:N));
  
  % get transformation matrices that transform original variables to
  % whitened variables and back
  S		= pcvec(:, 1:N) * diag(1./f);
  Sadj	        = diag(f) * pcvec(:, 1:N)';

  %% filter data matrix
  % Lanczos lowpass filter after detrending (CUTOFF is in number of timesteps, reflected boundary conditions)
  t = 1:n; t = t';
  X_f = Xs;
  for i = 1:size(Xs,2)
      p = polyfit(t,Xs(:,i),1);
      tmp = Xs(t,i)-p(1)*t-p(2);
      tmp = lanczos([flipud(tmp); tmp; flipud(tmp)],1,cutoff);
      X_f(:,i) = tmp((end/3+1):2*end/3)+p(1)*t+p(2);
  end
  
  % Simple Lanczos lowpass filter (periodic boundary conditions)
  %X_f = lanczos(Xs,1,cutoff);
  
  % OR INSERT ALTERNATE FILTER HERE
  
  % Added options for highpass filter
  %X_f = lanczos(X_lp,1,cutoff,[],'high');
  
  % Bandpass filter
  %X_f = lanczos(X,1,cutoff1);
  %X_f = lanczos(X_f,1,cutoff2,[],'high');
  
  % A simple running mean filter with zeros on the edges
  %X_f = X;
  %for i = 1:size(X,2)
  %    X_f(cutoff/2:end-cutoff/2,i) = rmean(X(:,i),cutoff);
  %end
  
  %% whiten variables [such that cov(Y) = I * n/(n-1)]
  Y		= X_f * S;

  % slow covariance matrix of whitened variables
  % (i.e. covariance matrix of filtered and whitened principal components)
  Gamma = cov(Y);

  %% SVD of slow covariance matrix (such that r are eigenvalues and V are eigenvectors)
  [~, r, V]	= csvd(Gamma);

  %% weight vectors (canonical vectors) and patterns (LFPs) in original scaling
  % Note: canonical vectors are called u_k in Wills et al. (2018)
  weights	= repmat(scale', 1, N) .* (S * V);       % weights are columns
  LFPs	= (V' * Sadj) ./ repmat(scale, N, 1);    % patterns are rows 
  
  % choose signs of patterns, weights, eofs, and pcs such that the
  % scalar product of the vectors and the scale vector is positive
  for j=1:size(LFPs, 1)
    if LFPs(j, :)*scale' < 0
      LFPs(j, :) = -LFPs(j, :);
      weights(:, j)  = -weights(:, j);
    end
  end
  
  for j=1:size(EOF, 1)
    if EOF(j, :)*scale' < 0
      EOF(j, :)  = -EOF(j, :);
      pcs(:, j)  = -pcs(:, j);
    end
  end
  
%% Low-frequency components (LFCs)
  
if nargin > 3
    Xs = Xs./repmat(scale,n,1);
end

LFCs = Xs * weights;

%% slow covariance of untruncated state space 
Cov_slow = cov(X_f); % X_f has already been scaled
trCovslow = trace(Cov_slow);

w = weights./repmat(scale', 1, N);
p = LFPs.*repmat(scale, N, 1);

slow_var = diag(p*Cov_slow*w)./diag(p*w);
tot_var = diag(p*Covtot*w)./diag(p*w);

slow_var_eofs = diag(pcvec'*Cov_slow*pcvec)./diag(pcvec'*pcvec);
tot_var_eofs = diag(pcvec'*Covtot*pcvec)./diag(pcvec'*pcvec);

% slow variance and total variance in each LFC, as a fraction of total slow
% variance and total variance
pvar_slow = slow_var./trCovslow*100;
pvar_lfc = tot_var./trCovtot*100;

% low-frequency to total variance ratio of principal components and
% fraction of slow variance in each principal component
r_eofs = slow_var_eofs./tot_var_eofs;
pvar_slow_eofs = slow_var_eofs./trCovslow*100;
  