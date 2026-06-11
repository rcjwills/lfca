function [detrended,delta_data_trend] = detrend(data,dim,order,t)

% note that this removes the trend, but does not remove constant offsets

if nargin<2
    dim = 1;
    order = 1;
end

if order > 4
    disp('Error: higher order than 4 not yet supported')
end

detrended = data;

switch dim
    case 3
        n = size(data,3);
        if nargin < 4
            t = 1:n; t = t';
        end
        if length(size(data))==4
            for k = 1:size(data,4)
                for i = 1:size(data,1)
                    for j = 1:size(data,2)
                        p = polyfit(t,squeeze(data(i,j,:,k)),order);
                        switch order
                            case 1
                                detrended(i,j,:,k) = squeeze(data(i,j,:,k)) - p(1).*t - p(2);
                            case 2
                                detrended(i,j,:,k) = squeeze(data(i,j,:,k)) - p(1).*t.^2 - p(2).*t - p(3);
                            case 3
                                detrended(i,j,:,k) = squeeze(data(i,j,:,k)) - p(1).*t.^3 - p(2).*t.^2 - p(3).*t - p(4);
                            case 4
                                detrended(i,j,:,k) = squeeze(data(i,j,:,k)) - p(1).*t.^4 - p(2).*t.^3 - p(3).*t.^2 - p(4).*t - p(5);
                        end
                    end
                end
            end
            trend = data-detrended;
            delta_data_trend = squeeze(trend(:,:,end,:)-trend(:,:,1,:));
        else
            for i = 1:size(data,1)
                for j = 1:size(data,2)
                    p = polyfit(t,squeeze(data(i,j,:)),order);
                    switch order
                        case 1
                            detrended(i,j,:) = squeeze(data(i,j,:)) - p(1).*t - p(2);
                        case 2
                            detrended(i,j,:) = squeeze(data(i,j,:)) - p(1).*t.^2 - p(2).*t - p(3);
                        case 3
                            detrended(i,j,:) = squeeze(data(i,j,:)) - p(1).*t.^3 - p(2).*t.^2 - p(3).*t - p(4);
                        case 4
                            detrended(i,j,:) = squeeze(data(i,j,:)) - p(1).*t.^4 - p(2).*t.^3 - p(3).*t.^2 - p(4).*t - p(5);
                    end
                end
            end
            trend = data-detrended;
            delta_data_trend = trend(:,:,end)-trend(:,:,1);
        end
        
    case 2
        n = size(data,2);
        if nargin < 4
            t = 1:n; t = t';
        end
        for i = 1:size(data,1)
                p = polyfit(t,squeeze(data(i,:)),order);
                switch order
                    case 1
                        detrended(i,:) = squeeze(data(i,:)) - p(1).*t;
                    case 2
                        detrended(i,:) = squeeze(data(i,:)) - p(1).*t.^2 - p(2).*t;
                    case 3
                        detrended(i,:) = squeeze(data(i,:)) - p(1).*t.^3 - p(2).*t.^2 - p(3).*t;
                    case 4
                        detrended(i,:) = squeeze(data(i,:)) - p(1).*t.^4 - p(2).*t.^3 - p(3).*t.^2 - p(4).*t;
                end
        end
        trend = data-detrended;
        delta_data_trend = trend(:,end)-trend(:,1);
    case 1
        if size(data,2) > 1
            if length(size(data)) > 2
                n = size(data,1);
                if nargin < 4
                    t = 1:n; t = t';
                end
                for i = 1:size(data,2)
                    for j = 1:size(data,3)
                        p = polyfit(t,squeeze(data(:,i,j)),order);
                        switch order
                            case 1
                                detrended(:,i,j) = squeeze(data(:,i,j)) - p(1).*t;
                            case 2
                                detrended(:,i,j) = squeeze(data(:,i,j)) - p(1).*t.^2 - p(2).*t;
                        end
                    end
                end
                trend = data-detrended;
                delta_data_trend = trend(end,:,:)-trend(1,:,:);
            else
                n = size(data,1);
                if nargin < 4
                    t = 1:n; t = t';
                end
                for i = 1:size(data,2)
                    p = polyfit(t,data(:,i),order);
                    switch order
                        case 1
                            detrended(:,i) = data(:,i) - p(1).*t;
                        case 2
                            detrended(:,i) = data(:,i) - p(1).*t.^2 - p(2).*t;
                    end
                end
                trend = data-detrended;
                delta_data_trend = trend(end,:)-trend(1,:);
            end
        else
            n = size(data,1);
            if nargin < 4
                t = 1:n; t = t';
            end
            p = polyfit(t,data,order);
            switch order
                case 1
                    detrended = data - p(1).*t;
                case 2
                    detrended = data - p(1).*t.^2 - p(2).*t;
                case 3
                    detrended = data - p(1).*t.^3 - p(2).*t.^2 - p(3).*t;
                case 4
                    detrended = data - p(1).*t.^4 - p(2).*t.^3 - p(3).*t.^2 - p(4).*t;
            end
            trend = data-detrended;
            delta_data_trend = trend(end)-trend(1);
        end
end