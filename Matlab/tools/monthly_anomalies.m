function [anomalies,mean_seasonal_cycle] = monthly_anomalies(data,dim)

if nargin < 2
   dim = 3;
end

anomalies = data;

data(data>1e20) = nan;

switch dim
    case 4
        mean_seasonal_cycle = zeros(size(data(:,:,:,1:12)));
        for i = 1:12
            mean_seasonal_cycle(:,:,:,i) = nanmean(data(:,:,:,i:12:end),4);
            anomalies(:,:,:,i:12:end) = anomalies(:,:,:,i:12:end) - mean_seasonal_cycle(:,:,:,i);
        end
    case 3
        mean_seasonal_cycle = zeros(size(data(:,:,1:12)));
        for i = 1:12
            mean_seasonal_cycle(:,:,i) = nanmean(data(:,:,i:12:end),3);
            anomalies(:,:,i:12:end) = anomalies(:,:,i:12:end) - mean_seasonal_cycle(:,:,i);
        end
    case 2
        mean_seasonal_cycle = zeros(size(data(:,1:12)));
        for i = 1:12
            mean_seasonal_cycle(:,i) = nanmean(data(:,i:12:end),2);
            anomalies(:,i:12:end) = anomalies(:,i:12:end) - mean_seasonal_cycle(:,i);
        end
    case 1
        mean_seasonal_cycle = zeros(size(data(1:12,:)));
        for i = 1:12
            mean_seasonal_cycle(i,:) = nanmean(data(i:12:end,:),1);
            anomalies(i:12:end,:) = anomalies(i:12:end,:) - mean_seasonal_cycle(i,:);
        end
end

