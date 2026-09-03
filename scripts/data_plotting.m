clear all
close all force hidden
clc

% this is a simple plotting script 

% this finds all the produced CSV's from 'Entropy_segment_cells_Batch'
% and then groups them by the number of conditions (user given)
% finds the replicates and plots the normalized area per time point for
% each condition with standard error of the sample mean (SEM) error bars

% this is for conditions
% not number of repeats (repeats handled automatically)
num_conditions = 4;

experiment_name = 'DIC_analysis_output2';
ovr_dir = fullfile(pwd,'output',experiment_name);
path_to_export = fullfile(ovr_dir,[experiment_name '.png']);

if ispc
    [~,message,~] = fileattrib([char(ovr_dir),'\*']);
else
    [~,message,~] = fileattrib([char(ovr_dir),'/*']);
end

fprintf('representing %s\n', ovr_dir)
fprintf('There are %i total files & folders in the overarching folder.\n',numel(message));

allExts = cellfun(@(s) s(end-2:end), {message.Name},'uni',0); % Get exts

CSVidx = ismember(allExts,'csv');    % Search ext for "CSV" at the end
CSV_filepaths = {message(CSVidx).Name}';  % Use CSVidx to list all paths.

CSV_filepaths = natsort(CSV_filepaths);
CSV_filepaths = CSV_filepaths(~contains(CSV_filepaths,'Example_output'));
fprintf('There are %i files with *.CSV exts.\n\n',numel(CSV_filepaths));

kmeans_idx = kmeans(char(string(CSV_filepaths)),num_conditions);
kmeant_idx_unique_not_sorted = unique(kmeans_idx,'stable');

condition = cell(1,num_conditions);
f = figure('Units','normalized','Position',[0,0.4,0.5,0.5]);
hold on
for i = 1:num_conditions
    
    these_csvs = CSV_filepaths(kmeans_idx==kmeant_idx_unique_not_sorted(i));
    
    [~,this_condition,~] = fileparts(these_csvs{1});
    condition{i} = replace(this_condition(2:end-7),'_','-'); %probably have to change this for specific use cases
    
    disp(['Found coniditon:' condition{i}])
    for j = 1:length(these_csvs)
        temp_table = readtable(these_csvs{j},'VariableNamingRule','preserve');
        if j == 1
            norm_data = temp_table.("Norm open");
        else
            norm_data = [norm_data,temp_table.("Norm open")];
        end
    end
    
    avg_norm_data = medfilt1(mean(norm_data,2),5);
    stds_data = std(norm_data,0,2)/sqrt(size(norm_data,2));
        
    errorbar(1:length(avg_norm_data),avg_norm_data,stds_data, ...
        'Linewidth',0.5);
    
end
legend(condition)

axis('square');
ax = gca;                 % Get current axes
ax.Toolbar.Visible = 'off'; % Hide the toolbar completely
saveas(f,path_to_export)

fprintf('\nExport path:')
disp(path_to_export)

disp('end of script')
close all

% [val,idx] = sort(sums);
% 
% figure;
% bar(val);
% xticklabels(condition(idx))





