function [tk, FPs, fingerprints, s, pvar, pcs, EOF, N, pvar_FPs, s_eofs] = forced_pattern_analysis(X, Xe, truncation, scale, Covtot)

%% Truncated Forced Analysis Analysis
%     [TK,FPS,FINGERPRINTS,S,PVAR,PCS,EOF,N,PVAR_FPS,S_EOFS] = forced_component_analysis(X,XE,TRUNCATION,SCALE,COVTOT)
%     performs forced pattern analysis / signal-to-noise maximizing EOF analysis (see Wills et al. 2020) 
%     on the data in matrix X based on a the ratio of ensemble-mean variance
%     to total variance (Xe gives the mean over ensemble members, which are concatenated in time).

%% INPUT
%     X is a 2D data matrix with time variantions along the first dimension
%     and spatial variations along the second dimension. Ensemble members
%     are concatenated in time.
%
%     XE is a 2D ensemble-mean data matrix with time variantions along the 
%     first dimension and spatial variations along the second dimension
%
%     TRUNCATION is the number of principal components / EOFs to include in
%     the analysis, or (if less than 1) the fraction of variance to retain
%
%     SCALE (optional) a scale vector, which for geospatial data should be 
%     equal to the square root of grid cell area. The default value is one  
%     for all grid points.
%
%     COVTOT (optional) the covariance matrix associated with the data in
%     X. If not specified, COVTOT will be computed from X.

%% OUTPUT

%     FINGERPRINTS is a matrix containing the canonical weight vectors as
%     columns. FPs is a matrix containing the dual vectors of the
%     canonical weight vectors as rows. These are the so-called forced 
%     patterns (FPs). S is a vector measuring the ratio of ensemble-mean
%     signal to total variance for each forced pattern
%
%     PVAR is the percentage of total sample variation accounted for
%     by each of the EOFs. PCS is a matrix containing the principal
%     component time series as columns. EOF is a matrix containing the
%     EOFs, the principal component patterns, as rows. The scalar N
%     is the rank at which the PCA was truncated.
%
%     S is a vector containing the ratio of ensemble-mean signal to total
%     variance for each forced pattern
%
%     PVAR_FPS is a vector of the variance associated with 
%     each forced pattern as a fraction of the total variance. Note that
%     the FPs are not orthogonal, so these values need not add to the
%     total variance in the first N principal components.
%
%     S_EOFS and PVAR are equivalent to S and PVAR_FPS respectively, but 
%     for the original EOFs.

  narginchk(3,5)          % check number of input arguments 
  if ndims(X) ~= 2,  error('Data matrix must be 2-D.'); end

  disp(sprintf('\nForced Pattern Analysis:'))

  [n,p]         = size(Xe);
  ne            = size(X,1)/n;
  
  % center data 
  if any(any(isnan(X)))               % there are missing values in x
    Xm  = nanmean(X);
    Xem  = nanmean(Xe);
  else                                % no missing values
    Xm  = mean(X);
    Xem  = mean(Xe);
  end
  X    = X - repmat(Xm, n*ne, 1);  
  Xe    = Xe - repmat(Xem, n, 1);  
      
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
    Xs           = X .* repmat(scale,n*ne,1);
    Xes           = Xe .* repmat(scale,n,1);
  else
    scale       = ones(1,p);
    Xs          = X;
    Xes          = Xe;
  end
  clear X Xe
  
  %% eigendecomposition of covariance matrix
  Covtot      = repmat(scale',1,p) .* Covtot .* repmat(scale,p,1);
  [pcvec,evl,rest] = peigs(Covtot, min(n-1, p));
  trCovtot    = trace(Covtot);
  
  % percent of total sample variation accounted for by each EOF
  pvar          = evl./trCovtot .* 100;
  % principal component time series
  pcs           = Xs*pcvec;
  pces          = Xes*pcvec;
  s_eofs        = var(pces,0,1)./var(pcs,0,1);
  % return EOFs in original scaling as patterns (row vectors)
  EOF           = pcvec' ./ repmat(scale,rest,1);
  
  %% truncation of EOFs
  if truncation < 1 
      % using basic % variance criterion, where truncation gives the
      % fraction of variance to be included in the EOF truncation
      truncation = truncation*100;
      cum_pvar = cumsum(pvar);
      N = find(abs(cum_pvar-truncation) == min(abs(cum_pvar-truncation)),1,'first');
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
  
  %% whiten variables [such that cov(Y) = I * n/(n-1)]
  Y		= Xes * S;

  % slow covariance matrix of whitened variables
  % (i.e. covariance matrix of filtered and whitened principal components)
  Gamma = cov(Y);

  %% SVD of slow covariance matrix (such that r are eigenvalues and V are eigenvectors)
  [~, s, V]	= csvd(Gamma);

  %% fingerprint patterns (canonical vectors) and forced patterns (FPs) in original scaling
  % Note: canonical vectors are called u_k in Wills et al. (2020)
  fingerprints	= repmat(scale', 1, N) .* (S * V);       % weights are columns
  FPs	= (V' * Sadj) ./ repmat(scale, N, 1);    % patterns are rows 
  
  % choose signs of patterns, weights, eofs, and pcs such that the
  % scalar product of the vectors and the scale vector is positive
  for j=1:size(FPs, 1)
    if FPs(j, :)*scale' < 0
      FPs(j, :) = -FPs(j, :);
      fingerprints(:, j)  = -fingerprints(:, j);
    end
  end
  
  for j=1:size(EOF, 1)
    if EOF(j, :)*scale' < 0
      EOF(j, :)  = -EOF(j, :);
      pcs(:, j)  = -pcs(:, j);
    end
  end
  
%% timeseries (tk)
  
if nargin > 3
    Xs = Xs./repmat(scale,n*ne,1);
end

tk = Xs * fingerprints;

%% fraction of variance in forced patterns

w = fingerprints./repmat(scale', 1, N);
p = FPs.*repmat(scale, N, 1);

tot_var = diag(p*Covtot*w)./diag(p*w);

pvar_FPs = tot_var./trCovtot*100;
  