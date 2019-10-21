function [field_filt] = lanczos_lowpass(field,dim,cutoff,optional1,optional2)

% enter cutoff in number of timesteps

s = size(field);
field_filt = zeros(s);

if dim == 1
    for i = 1:s(2)
        if length(s) > 2
            for j = 1:s(3)
                if nargin > 3
                    field_filt(:,i,j) = lanczosfilter(squeeze(field(:,i,j)),1,1/cutoff,optional1,optional2);
                else
                    field_filt(:,i,j) = lanczosfilter(squeeze(field(:,i,j)),1,1/cutoff);
                end
            end
        else
            if nargin > 3
                field_filt(:,i) = lanczosfilter(squeeze(field(:,i)),1,1/cutoff,optional1,optional2);
            else
                field_filt(:,i) = lanczosfilter(squeeze(field(:,i)),1,1/cutoff);
            end
        end
%         if mod(i,25) == 0
%             perct = i/s(2)*100;
%             disp([num2str(perct),'% complete']);
%         end
    end
elseif dim == 2
    for i = 1:s(1)
        if nargin > 3
            field_filt(i,:) = lanczosfilter(squeeze(field(i,:)),1,1/cutoff,optional1,optional2);
        else
            field_filt(i,:) = lanczosfilter(squeeze(field(i,:)),1,1/cutoff);
        end
        if mod(i,25) == 0
            perct = i/s(2)*100;
            disp([num2str(perct),'% complete']);
        end
    end
elseif dim == 3
    for i = 1:s(1)
        for j = 1:s(2)
            if nargin > 3
                field_filt(i,j,:) = lanczosfilter(squeeze(field(i,j,:)),1,1/cutoff,optional1,optional2);
            else
                field_filt(i,j,:) = lanczosfilter(squeeze(field(i,j,:)),1,1/cutoff);
            end
        end
        if mod(i,50) == 0
            perct = i/s(2)*100;
            disp([num2str(perct),'% complete']);
        end
    end
end
                