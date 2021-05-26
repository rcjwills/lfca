load('ERSST_1900_2016.mat','LON_AXIS','LAT_AXIS','SST')
SST = SST(:,:,end-288+1:end);
load('AVISO_1993_2016.mat','lon','lat','ssh')
time = 1993:1/12:2016.99;

%% Parameters
cutoff = 120; % number of timesteps (months in this case)
truncation = 30; % number of EOFs

% also try truncation = 3

%% Preprocessing
% compute anomaly from annual cycle
[SST_anomalies,Mt] = monthly_anomalies(SST);
[SSH_anomalies,Mt2] = monthly_anomalies(ssh);

%% reshape data1 for LFCA
s = size(SST_anomalies);
[Y,X] = meshgrid(LAT_AXIS,LON_AXIS);
area = cos(Y*pi/180);
area(isnan(mean(SST_anomalies,3))) = 0;

% North_Atlantic domain
domain = ones(size(area));
domain(Y>70) = 0;
domain(Y<0) = 0;
domain(X<260 & X>43 & Y>17) = 0;
domain(X<260 & X>30 & Y<31) = 0;
domain(X<270 & X>36 & Y<=17 & Y>14) = 0;
domain(X<276 & X>36 & Y<=14  & Y>=9) = 0;
domain(X<290 & X>24 & Y<9) = 0;

xtf1 = reshape(SST_anomalies,s(1)*s(2),s(3))';
AREA_WEIGHTS_1 = reshape(area,s(1)*s(2),1)';
domain = reshape(domain,s(1)*s(2),1)';

% icol_ret and icol_disc help reconstruct the data onto the original grid
icol_ret_1 = find(AREA_WEIGHTS_1~=0 & domain);
icol_disc_1 = find(AREA_WEIGHTS_1==0 | ~domain);
xtf1 = xtf1(:,icol_ret_1);
AREA_WEIGHTS_1 = AREA_WEIGHTS_1(icol_ret_1);

% scale by square root of grid cell area such that covariance is area
% weighted
normvec          = AREA_WEIGHTS_1' ./ sum(AREA_WEIGHTS_1);
scale1    = sqrt(normvec);

%% reshape data2 for LFCA
s = size(SSH_anomalies);
[Y,X] = meshgrid(lat,lon);
area = cos(Y*pi/180);
area(isnan(mean(SSH_anomalies,3))) = 0;

% North_Atlantic domain
domain = ones(size(area));
domain(Y>70) = 0;
domain(Y<0) = 0;
domain(X<260 & X>43 & Y>17) = 0;
domain(X<260 & X>30 & Y<31) = 0;
domain(X<270 & X>36 & Y<=17 & Y>14) = 0;
domain(X<276 & X>36 & Y<=14  & Y>=9) = 0;
domain(X<290 & X>24 & Y<9) = 0;

xtf2 = reshape(SSH_anomalies,s(1)*s(2),s(3))';
AREA_WEIGHTS_2 = reshape(area,s(1)*s(2),1)';
domain = reshape(domain,s(1)*s(2),1)';

% icol_ret and icol_disc help reconstruct the data onto the original grid
icol_ret_2 = find(AREA_WEIGHTS_2~=0 & domain);
icol_disc_2 = find(AREA_WEIGHTS_2==0 | ~domain);
xtf2 = xtf2(:,icol_ret_2);
AREA_WEIGHTS_2 = AREA_WEIGHTS_2(icol_ret_2);

% scale by square root of grid cell area such that covariance is area
% weighted
normvec          = AREA_WEIGHTS_2' ./ sum(AREA_WEIGHTS_2);
scale2    = sqrt(normvec);

%% Combine variables for analysis

% scale each variable by sqrt of its total variance
[n,p]         = size(xtf1);
Cov1 = cov(xtf1);
Cov1      = repmat(scale1,1,p) .* Cov1 .* repmat(scale1',p,1);
scale_factor_1 = sqrt(sum(diag(Cov1)));

[n,p]         = size(xtf2);
Cov2 = cov(xtf2);
Cov2      = repmat(scale2,1,p) .* Cov2 .* repmat(scale2',p,1);
scale_factor_2 = sqrt(sum(diag(Cov2)));

xtf = [xtf1./scale_factor_1 xtf2./scale_factor_2]; 
scale = [scale1; scale2];

is1 = 1:length(scale1);
is2 = length(scale1)+1:size(xtf,2);

%% Low-frequency component analysis (LFCA)
[LFCs, LFPs, weights, r, pvar, PCs, EOFs, N, pvar_slow, pvar_lfc, r_eofs, pvar_slow_eofs] = lfca(xtf, cutoff, truncation, scale);

LFPs_1       = insert_cols(LFPs(:,is1), icol_ret_1, icol_disc_1);
EOFs_1       = insert_cols(EOFs(:,is1), icol_ret_1, icol_disc_1);
weightsf_1        = insert_rows(weights(is1,:), icol_ret_1, icol_disc_1);
LFPs_2       = insert_cols(LFPs(:,is2), icol_ret_2, icol_disc_2);
EOFs_2       = insert_cols(EOFs(:,is2), icol_ret_2, icol_disc_2);
weightsf_2        = insert_rows(weights(is2,:), icol_ret_2, icol_disc_2);

LFPs_1 = LFPs_1.*scale_factor_1;
EOFS_1 = EOFs_1.*scale_factor_1;
LFPs_2 = LFPs_2.*scale_factor_2;
EOFS_2 = EOFs_2.*scale_factor_2;

%% plot results

ctrs1 = linspace(-0.8,0.8,25);
ctrs2 = linspace(-6,6,25);
% SSH is multiplied by 100 to plot in cm

plot_patterns(LFCs,LFPs_1,1,time,LON_AXIS,LAT_AXIS,ctrs1,'Atlantic');
plot_patterns(LFCs,100*LFPs_2,1,time,lon,lat,ctrs2,'Atlantic');

plot_patterns(LFCs,LFPs_1,2,time,LON_AXIS,LAT_AXIS,ctrs1,'Atlantic');
plot_patterns(LFCs,100*LFPs_2,2,time,lon,lat,ctrs2,'Atlantic');

plot_patterns(LFCs,LFPs_1,3,time,LON_AXIS,LAT_AXIS,ctrs1,'Atlantic');
plot_patterns(LFCs,100*LFPs_2,3,time,lon,lat,ctrs2,'Atlantic');

plot_patterns(LFCs,LFPs_1,4,time,LON_AXIS,LAT_AXIS,ctrs1,'Atlantic');
plot_patterns(LFCs,1000*LFPs_2,4,time,lon,lat,ctrs2,'Atlantic');
