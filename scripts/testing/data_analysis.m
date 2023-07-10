clear all
close all force hidden

% This is a not ready for production yet processing step that groups ALL of
% the produced csv's in the >output folder into num_conditions based off
% kmeans grouping of the actual strings 
% this is when you have multiple replicates of the the same condition

num_conditions = 8;

ovr_dir = fullfile(pwd,'output');

if ispc
    [~,message,~] = fileattrib([ovr_dir,'\*']);
else
    [~,message,~] = fileattrib([ovr_dir,'/*']);
end

fprintf('\nThere are %i total files & folders in the overarching folder.\n',numel(message));

allExts = cellfun(@(s) s(end-2:end), {message.Name},'uni',0); % Get exts

CSVidx = ismember(allExts,'csv');    % Search ext for "CSV" at the end
CSV_filepaths = {message(CSVidx).Name}';  % Use CSVidx to list all paths.

CSV_filepaths = natsort(CSV_filepaths);
CSV_filepaths = CSV_filepaths(~contains(CSV_filepaths,'Example_output'));
fprintf('There are %i files with *.CSV exts.\n',numel(CSV_filepaths));

kmeans_idx = kmeans(char(string(CSV_filepaths)),num_conditions);
kmeant_idx_unique_not_sorted = unique(kmeans_idx,'stable');

hold on
for i = 1:num_conditions
    
    these_csvs = CSV_filepaths(kmeans_idx==kmeant_idx_unique_not_sorted(i));
    
    [~,this_condition,~] = fileparts(these_csvs{1});
    condition{i} = replace(this_condition(2:end-7),'_','-'); %probably have to change this for specific use cases
    
    for j = 1:length(these_csvs)
        temp_table = readtable(these_csvs{j},'VariableNamingRule','preserve');
        if j == 1
            norm_data = temp_table.("Norm open");
        else
            norm_data = [norm_data,temp_table.("Norm open")];
        end
    end
    
    avg_norm_data = medfilt1(mean(norm_data,2),5);
    stds_data = std(norm_data,0,2);
    
    sums(i) = sum(avg_norm_data);
    
    errorbar(1:length(avg_norm_data),avg_norm_data,stds_data);
    
end
legend(condition)

[val,idx] = sort(sums);

figure;
bar(val);
xticklabels(condition(idx))





