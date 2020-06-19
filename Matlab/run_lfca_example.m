load('ERSST_1900_2016.mat','LON_AXIS','LAT_AXIS','SST')
time = 1900:1/12:2016.99;

%% Parameters
cutoff = 120; % number of timesteps (months in this case)
truncation = 30; % number of EOFs

% also try truncation = 3

%% Preprocessing
% compute anomaly from annual cycle
[SST_anomalies,Mt] = monthly_anomalies(SST);

% reshape data for LFCA
s = size(SST_anomalies);
[Y,X] = meshgrid(LAT_AXIS,LON_AXIS);
area = cos(Y*pi/180);
area(isnan(mean(SST_anomalies,3))) = 0;

% Pacific domain
domain = ones(size(area));
domain(X<100) = 0;
domain(X<103 & Y<5) = 0;
domain(X<105 & Y<2) = 0;
domain(X<111 & Y<-6) = 0;
domain(X<114 & Y<-7) = 0;
domain(X<127 & Y<-8) = 0;
domain(X<147 & Y<-18) = 0;
domain(Y>70) = 0;
domain(Y>65 & (X<175 | X>200)) = 0;
domain(Y<-45) = 0;
domain(X>260 & Y>17) = 0;
domain(X>270 & Y<=17 & Y>14) = 0;
domain(X>276 & Y<=14 & Y>9) = 0;
domain(X>290 & Y<=9) = 0;

% Note: X is changing uses from this point forward. Was longitude array, now is data array.
X = reshape(SST_anomalies,s(1)*s(2),s(3))';
AREA_WEIGHTS = reshape(area,s(1)*s(2),1)';
domain = reshape(domain,s(1)*s(2),1)';

% icol_ret and icol_disc help reconstruct the data onto the original grid
icol_ret = find(AREA_WEIGHTS~=0 & domain);
icol_disc = find(AREA_WEIGHTS==0 | ~domain);
X = X(:,icol_ret);
AREA_WEIGHTS = AREA_WEIGHTS(icol_ret);

% scale by square root of grid cell area such that covariance is area
% weighted
normvec          = AREA_WEIGHTS' ./ sum(AREA_WEIGHTS);
scale    = sqrt(normvec);

%% Low-frequency component analysis (LFCA)
[LFCs, LFPs, weights, r, pvar, PCs, EOFs, N, pvar_slow, pvar_lfc, r_eofs, pvar_slow_eofs] = lfca(X, cutoff, truncation, scale);
LFPs       = insert_cols(LFPs, icol_ret, icol_disc);
EOFs       = insert_cols(EOFs, icol_ret, icol_disc);
weightsf        = insert_rows(weights, icol_ret, icol_disc);

%% plot results
if truncation < 5
    ctrs = linspace(-1,1,21);
else
    ctrs = linspace(-0.6,0.6,25);
end
plot_patterns(LFCs,LFPs,1,time,LON_AXIS,LAT_AXIS,ctrs,'Pacific');
plot_patterns(LFCs,LFPs,2,time,LON_AXIS,LAT_AXIS,ctrs,'Pacific');
plot_patterns(LFCs,LFPs,3,time,LON_AXIS,LAT_AXIS,ctrs,'Pacific');
if truncation > 3
    plot_patterns(LFCs,LFPs,4,time,LON_AXIS,LAT_AXIS,ctrs,'Pacific');
end
