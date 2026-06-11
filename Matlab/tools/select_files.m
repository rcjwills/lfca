function [files] = select_files(files,included_string,excluded_string)

is = [];

if nargin > 2
    for i = 1:length(files)
        if contains(files{i},included_string) && ~contains(files{i},excluded_string)
            is = [is, i];
        end
    end
else
    for i = 1:length(files)
        if contains(files{i},included_string)
            is = [is, i];
        end
    end
end

files = files(is);