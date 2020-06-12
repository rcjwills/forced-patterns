load('CESM-LE_hist-rcp85_1920-2019_3monthly_ts.mat')
nlon = size(ts_CESM_all,1);
nlat = size(ts_CESM_all,2);
n = size(ts_CESM_all,3); 
ne = size(ts_CESM_all,4);
time = 1920:0.25:2019.75;

% remove the ensemble-mean seasonal cycle 
% -- Note that Wills et al. (2020) removes the mean seasonal cycle from each
% ensemble member, but this does not make a large difference
ts_CESM_anomalies_all = ts_CESM_all;
ts_CESM_clim = zeros([nlon nlat 4]);
for m = 1:4
    ts_CESM_clim(:,:,m) = squeeze(mean(mean(ts_CESM_all(:,:,m:4:end,:),4),3));
    ts_CESM_anomalies_all(:,:,m:4:end,:) = ts_CESM_anomalies_all(:,:,m:4:end,:) - repmat(ts_CESM_clim(:,:,m),[1 1 100 40]);
end

%% Parameters
truncation = 150; % number of EOFs
M = 8; % number of Forced Patterns to retain in FP filtering (higher than in paper due to differing resolution and differing eigenvalue spectrum)

%% Preprocessing

% reshape data
ts_CESM_anomalies_emean = mean(ts_CESM_anomalies_all,4);
ts_CESM_anomalies_1_20_emean = mean(ts_CESM_anomalies_all(:,:,:,1:20),4);
ts_CESM_anomalies_21_40_emean = mean(ts_CESM_anomalies_all(:,:,:,21:40),4);
ts_CESM_anomalies_all = reshape(ts_CESM_anomalies_all,[nlon nlat n*ne]); % concatenate ensemble members in time
s = size(ts_CESM_anomalies_all);

[Y,~] = meshgrid(LAT_AXIS,LON_AXIS);
area = cos(Y*pi/180);
%area(isnan(mean(ts_CESM_all,3))) = 0; % CESM-LE doesn't have NaNs,
%but enabling this line will ignore grid points that have NaNs

domain = ones(size(area));

% Can specify sub-domains to analyze using the following form:

% % Pacific domain
% domain(X<100) = 0;
% domain(X<103 & Y<5) = 0;
% domain(X<105 & Y<2) = 0;
% domain(X<111 & Y<-6) = 0;
% domain(X<114 & Y<-7) = 0;
% domain(X<127 & Y<-8) = 0;
% domain(X<147 & Y<-18) = 0;
% domain(Y>70) = 0;
% domain(Y>65 & (X<175 | X>200)) = 0;
% domain(Y<-45) = 0;
% domain(X>260 & Y>17) = 0;
% domain(X>270 & Y<=17 & Y>14) = 0;
% domain(X>276 & Y<=14 & Y>9) = 0;
% domain(X>290 & Y<=9) = 0;

X = reshape(ts_CESM_anomalies_all,s(1)*s(2),s(3))';
Xe_1_20 = reshape(ts_CESM_anomalies_1_20_emean,s(1)*s(2),n)';
Xe_21_40 = reshape(ts_CESM_anomalies_21_40_emean,s(1)*s(2),n)';
Xe = reshape(ts_CESM_anomalies_emean,s(1)*s(2),n)';
AREA_WEIGHTS = reshape(area,s(1)*s(2),1)';
domain = reshape(domain,s(1)*s(2),1)';

% icol_ret and icol_disc help reconstruct the data onto the original grid
icol_ret = find(AREA_WEIGHTS~=0 & domain);
icol_disc = find(AREA_WEIGHTS==0 | ~domain);
X = X(:,icol_ret);
Xe = Xe(:,icol_ret);
AREA_WEIGHTS = AREA_WEIGHTS(icol_ret);

% scale by square root of grid cell area such that covariance is area
% weighted
normvec          = AREA_WEIGHTS' ./ sum(AREA_WEIGHTS);
scale    = sqrt(normvec);

%% Calculate forced patterns / Signal-to-noise maximizing EOF analysis (sub-ensembles)
[tk, FPs, fingerprints, s, pvar, pcs, EOFs, N, pvar_FPs, s_eofs] = forced_pattern_analysis(X(1:end/2,:), Xe_1_20, truncation, scale);
FPs       = insert_cols(FPs, icol_ret, icol_disc);
Xe_1_20_filtered = tk(:,1:M)*FPs(1:M,:);
Xe_1_20_filtered = ensemble_average_timeseries(Xe_1_20_filtered,ne/2);

[tk, FPs, fingerprints, s, pvar, pcs, EOFs, N, pvar_FPs, s_eofs] = forced_pattern_analysis(X(end/2+1:end,:), Xe_21_40, truncation, scale);
FPs       = insert_cols(FPs, icol_ret, icol_disc);
Xe_21_40_filtered = tk(:,1:M)*FPs(1:M,:);
Xe_21_40_filtered = ensemble_average_timeseries(Xe_21_40_filtered,ne/2);

%% Calculate forced patterns / Signal-to-noise maximizing EOF analysis
[tk, FPs, fingerprints, s, pvar, pcs, EOFs, N, pvar_FPs, s_eofs] = forced_pattern_analysis(X, Xe, truncation, scale);
FPs       = insert_cols(FPs, icol_ret, icol_disc);
EOFs       = insert_cols(EOFs, icol_ret, icol_disc);
fingerprintsf        = insert_rows(fingerprints, icol_ret, icol_disc);

%% plot Forced Patterns

ndisc = 1; plot_ensemble_patterns(tk,FPs,ndisc,400,1920,LON_AXIS,LAT_AXIS,linspace(-0.8,0.8,25));
ndisc = 2; plot_ensemble_patterns(tk,FPs,ndisc,400,1920,LON_AXIS,LAT_AXIS,linspace(-0.8,0.8,25));

%% FP filtering (remove internal variability by truncating higher-order FPs)

load('landmasks.mat','ocean','land')

Xe_filtered = tk(:,1:M)*FPs(1:M,:);
Xe_filtered = ensemble_average_timeseries(Xe_filtered,ne);

field = reshape(Xe_filtered',[nlon,nlat n]);
Xe_EFCA_global = global_mean(LON_AXIS,LAT_AXIS,field);
Xe_EFCA_east_west = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[-6 6],[210 270]) - mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[-6 6],[120 180]);
Xe_EFCA_NA = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[40 60],[300 360]);
Xe_EFCA_USA = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*land,[30 45],[190 310]);
Xe_EFCA_Arctic = mean_in_a_box(LON_AXIS,LAT_AXIS,field,[65 90],[0 360]);
Xe_EFCA_Seattle = squeeze(field(closest(LON_AXIS,360-122.3),closest(LAT_AXIS,47.6),:));

field = reshape(Xe_1_20_filtered',[nlon,nlat n]);
Xe_EFCA_1_20_global = global_mean(LON_AXIS,LAT_AXIS,field);
Xe_EFCA_1_20_east_west = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[-6 6],[210 270]) - mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[-6 6],[120 180]);
Xe_EFCA_1_20_NA = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[40 60],[300 360]);
Xe_EFCA_1_20_USA = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*land,[30 45],[190 310]);
Xe_EFCA_1_20_Arctic = mean_in_a_box(LON_AXIS,LAT_AXIS,field,[65 90],[0 360]);
Xe_EFCA_1_20_Seattle = squeeze(field(closest(LON_AXIS,360-122.3),closest(LAT_AXIS,47.6),:));

field = reshape(Xe_21_40_filtered',[nlon,nlat n]);
Xe_EFCA_21_40_global = global_mean(LON_AXIS,LAT_AXIS,field);
Xe_EFCA_21_40_east_west = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[-6 6],[210 270]) - mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[-6 6],[120 180]);
Xe_EFCA_21_40_NA = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[40 60],[300 360]);
Xe_EFCA_21_40_USA = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*land,[30 45],[190 310]);
Xe_EFCA_21_40_Arctic = mean_in_a_box(LON_AXIS,LAT_AXIS,field,[65 90],[0 360]);
Xe_EFCA_21_40_Seattle = squeeze(field(closest(LON_AXIS,360-122.3),closest(LAT_AXIS,47.6),:));

field = ts_CESM_anomalies_emean;
Xe_global = global_mean(LON_AXIS,LAT_AXIS,field);
Xe_east_west = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[-6 6],[210 270]) - mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[-6 6],[120 180]);
Xe_NA = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[40 60],[300 360]);
Xe_USA = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*land,[30 45],[190 310]);
Xe_Arctic = mean_in_a_box(LON_AXIS,LAT_AXIS,field,[65 90],[0 360]);
Xe_Seattle = squeeze(field(closest(LON_AXIS,360-122.3),closest(LAT_AXIS,47.6),:));

field = ts_CESM_anomalies_1_20_emean;
Xe_1_20_global = global_mean(LON_AXIS,LAT_AXIS,field);
Xe_1_20_east_west = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[-6 6],[210 270]) - mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[-6 6],[120 180]);
Xe_1_20_NA = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[40 60],[300 360]);
Xe_1_20_USA = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*land,[30 45],[190 310]);
Xe_1_20_Arctic = mean_in_a_box(LON_AXIS,LAT_AXIS,field,[65 90],[0 360]);
Xe_1_20_Seattle = squeeze(field(closest(LON_AXIS,360-122.3),closest(LAT_AXIS,47.6),:));

field = ts_CESM_anomalies_21_40_emean;
Xe_21_40_global = global_mean(LON_AXIS,LAT_AXIS,field);
Xe_21_40_east_west = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[-6 6],[210 270]) - mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[-6 6],[120 180]);
Xe_21_40_NA = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*ocean,[40 60],[300 360]);
Xe_21_40_USA = mean_in_a_box(LON_AXIS,LAT_AXIS,field.*land,[30 45],[190 310]);
Xe_21_40_Arctic = mean_in_a_box(LON_AXIS,LAT_AXIS,field,[65 90],[0 360]);
Xe_21_40_Seattle = squeeze(field(closest(LON_AXIS,360-122.3),closest(LAT_AXIS,47.6),:));

%% plot pattern filtering results

t = 1920:0.25:2019.75;

figure; plot(t,Xe_1_20_global,'linewidth',0.7); 
hold on; plot(t,Xe_21_40_global,'linewidth',0.7); 
plot(t,Xe_global,'k','linewidth',1.5); 
pretty_figure(550,275,'Year','Global-Mean Temp. (°C)','none','none',16);
set(gca,'xlim',[1920 2020])
set(gca,'ygrid','on')
set(gca,'ylim',[-0.4 1]); set(gca,'ytick',-0.4:0.2:1)
%corr(Xe_1_20_global',Xe_21_40_global').^2

figure; plot(t,Xe_EFCA_1_20_global,'linewidth',0.7); 
hold on; plot(t,Xe_EFCA_21_40_global,'linewidth',0.7); 
plot(t,Xe_EFCA_global,'k','linewidth',1.5); 
pretty_figure(550,275,'Year','Global-Mean Temp. (°C)','none','none',16);
set(gca,'xlim',[1920 2020])
set(gca,'ygrid','on')
set(gca,'ylim',[-0.4 1]); set(gca,'ytick',-0.4:0.2:1)
%corr(Xe_EFCA_1_20_global',Xe_EFCA_21_40_global').^2

% %

figure; plot(t,Xe_1_20_NA,'linewidth',0.7); 
hold on; plot(t,Xe_21_40_NA,'linewidth',0.7); 
plot(t,Xe_NA,'k','linewidth',1.5); 
pretty_figure(550,275,'Year','NA SST 40-60°N (°C)','none','none',16);
set(gca,'xlim',[1920 2020])
set(gca,'ygrid','on')
set(gca,'ylim',[-0.6 0.6]); set(gca,'ytick',-0.6:0.3:0.6)
%corr(Xe_1_20_NA',Xe_21_40_NA').^2

figure; plot(t,Xe_EFCA_1_20_NA,'linewidth',0.7); 
hold on; plot(t,Xe_EFCA_21_40_NA,'linewidth',0.7); 
plot(t,Xe_EFCA_NA,'k','linewidth',1.5); 
pretty_figure(550,275,'Year','NA SST 40-60°N (°C)','none','none',16);
set(gca,'xlim',[1920 2020])
set(gca,'ygrid','on')
set(gca,'ylim',[-0.6 0.6]); set(gca,'ytick',-0.6:0.3:0.6)
%corr(Xe_EFCA_1_20_NA',Xe_EFCA_21_40_NA').^2

% %

figure; plot(t,Xe_1_20_east_west,'linewidth',0.7); 
hold on; plot(t,Xe_21_40_east_west,'linewidth',0.7); 
plot(t,Xe_east_west,'k','linewidth',1.5); 
pretty_figure(550,275,'Year','Pac. East-West SST Diff. (°C)','none','none',16);
set(gca,'xlim',[1920 2020])
set(gca,'ygrid','on')
set(gca,'ylim',[-0.8 1]); set(gca,'ytick',-0.8:0.4:1.2)
%corr(Xe_1_20_east_west',Xe_21_40_east_west').^2

figure; plot(t,Xe_EFCA_1_20_east_west,'linewidth',0.7); 
hold on; plot(t,Xe_EFCA_21_40_east_west,'linewidth',0.7); 
plot(t,Xe_EFCA_east_west,'k','linewidth',1.5); 
pretty_figure(550,275,'Year','Pac. East-West SST Diff. (°C)','none','none',16);
set(gca,'xlim',[1920 2020])
set(gca,'ygrid','on')
set(gca,'ylim',[-0.8 1]); set(gca,'ytick',-0.8:0.4:1.2)
%corr(Xe_EFCA_1_20_east_west',Xe_EFCA_21_40_east_west').^2

% %

figure; plot(t,Xe_1_20_USA,'linewidth',0.7); 
hold on; plot(t,Xe_21_40_USA,'linewidth',0.7); 
plot(t,Xe_USA,'k','linewidth',1.5); 
pretty_figure(550,275,'Year','US Temp. 30-45°N (°C)','none','none',16);
set(gca,'xlim',[1920 2020])
set(gca,'ygrid','on')
set(gca,'ylim',[-1.25 1.75]); set(gca,'ytick',-1:0.5:1.5)
%corr(Xe_1_20_USA',Xe_21_40_USA').^2

figure; plot(t,Xe_EFCA_1_20_USA,'linewidth',0.7); 
hold on; plot(t,Xe_EFCA_21_40_USA,'linewidth',0.7); 
plot(t,Xe_EFCA_USA,'k','linewidth',1.5); 
pretty_figure(550,275,'Year','US Temp. 30-45°N (°C)','none','none',16);
set(gca,'xlim',[1920 2020])
set(gca,'ygrid','on')
set(gca,'ylim',[-1.25 1.75]); set(gca,'ytick',-1:0.5:1.5)
%corr(Xe_EFCA_1_20_USA',Xe_EFCA_21_40_USA').^2

if(1) % extra plots
    figure; plot(t,Xe_1_20_Arctic,'linewidth',0.7);
    hold on; plot(t,Xe_21_40_Arctic,'linewidth',0.7);
    plot(t,Xe_Arctic,'k','linewidth',1.5);
    pretty_figure(550,275,'Year','Arctic Temp. 65-90°N (°C)','none','none',16);
    set(gca,'xlim',[1920 2020])
    set(gca,'ygrid','on')
    set(gca,'ylim',[-2 4]); set(gca,'ytick',-2:2:4)
    %corr(Xe_1_20_Arctic',Xe_21_40_Arctic').^2
    
    figure; plot(t,Xe_EFCA_1_20_Arctic,'linewidth',0.7);
    hold on; plot(t,Xe_EFCA_21_40_Arctic,'linewidth',0.7);
    plot(t,Xe_EFCA_Arctic,'k','linewidth',1.5);
    pretty_figure(550,275,'Year','Arctic Temp. 65-90°N (°C)','none','none',16);
    set(gca,'xlim',[1920 2020])
    set(gca,'ygrid','on')
    set(gca,'ylim',[-2 4]); set(gca,'ytick',-2:2:4)
    %corr(Xe_1_20_EFCA_Arctic',Xe_EFCA_21_40_Arctic').^2
    
    figure; plot(t,Xe_1_20_Seattle,'linewidth',0.7);
    hold on; plot(t,Xe_21_40_Seattle,'linewidth',0.7);
    plot(t,Xe_Seattle,'k','linewidth',1.5);
    pretty_figure(550,275,'Year','Seattle Temp. (°C)','none','none',16);
    set(gca,'xlim',[1920 2020])
    set(gca,'ygrid','on')
    set(gca,'ylim',[-0.8 1.2]); set(gca,'ytick',-0.8:0.4:1.2)
    %corr(Xe_1_20_Seattle',Xe_21_40_Seattle').^2
    
    figure; plot(t,Xe_EFCA_1_20_Seattle,'linewidth',0.7);
    hold on; plot(t,Xe_EFCA_21_40_Seattle,'linewidth',0.7);
    plot(t,Xe_EFCA_Seattle,'k','linewidth',1.5);
    pretty_figure(550,275,'Year','Seattle Temp. (°C)','none','none',16);
    set(gca,'xlim',[1920 2020])
    set(gca,'ygrid','on')
    set(gca,'ylim',[-0.8 1.2]); set(gca,'ytick',-0.8:0.4:1.2)
    %corr(Xe_EFCA_1_20_Seattle',Xe_EFCA_21_40_Seattle').^2
end
