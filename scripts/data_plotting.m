clear all
close all force hidden
clc

% this is a simple plotting script 

% this finds all the produced CSV's from 'Entropy_segment_cells_Batch'
% and then groups them by the number of conditions (user given or automatic)
% finds the replicates and plots the normalized area per time point for
% each condition with standard error of the sample mean (SEM) error bars

% this is for conditions
% not number of repeats (repeats handled automatically)
% if all experiments are named exactly the same and just have numbers
% appended to the end XXXexp_1,XXXexp_2,XXXexp_3,...,XXXexp_N
% then you can leave this as 0 and it will automatically group them and
% then find the number of conditions
num_conditions = 0;

experiment_name = 'DIC_analysis_output';
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

% get the CSV file names and remove the *_data.csv suffix
[~,CSV_filenames,~] = fileparts(CSV_filepaths);
CSV_filenames = cellfun(@(s) s(1:end-5),CSV_filenames,'UniformOutput',false);

% determine the number of conditions
if num_conditions == 0
    % check to see if the last character is a repeat number
    condition_counter = cell(length(CSV_filenames),1);
    result_counter = 0;
    for i = 1:length(CSV_filenames)
        [result,matches] = is_last_part_numeric(CSV_filenames{i});
        if result
            condition_counter{i} = CSV_filenames{i}(1:end-length(matches));
            result_counter = result_counter + 1;
        end
    end
    if result_counter == length(CSV_filenames)
        found_conditions = unique(natsort(condition_counter));
        disp('Auto generated conditions')
        for i = 1:length(found_conditions)
            disp(found_conditions{i})
        end
        num_conditions = length(found_conditions);
        disp(['Number of generated conditions: ' num2str(num_conditions)])
        disp('')
    else
        error('%s\n%s', 'Incorrect naming scheme for auto generation of conditions/replicates', ...
                'please specify the number of conditions')
    end
end

kmeans_idx = kmeans(char(string(CSV_filepaths)),num_conditions);
kmeant_idx_unique_not_sorted = unique(kmeans_idx,'stable');

condition = cell(1,num_conditions);
f = figure('Units','normalized','Position',[0,0.4,0.5,0.5]);
title(experiment_name,'Interpreter','none')
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

function [result,matches] = is_last_part_numeric(inputStr)
% IS_LAST_PART_NUMERIC Determines if the final part of a string is numeric.
%   The function checks if the string ends with a sequence of one or more 
%   digits. This handles various delimiters (like '_' or '-') and cases 
%   where the number is directly appended to text, as well as leading zeros.
%
%   Example usage:
%       is_last_part_numeric('_RENCA_WT_10FBS_DMSO_400')  -> true
%       is_last_part_numeric('_RENCA_WT_10FBS_DMSO_004')  -> true
%       is_last_part_numeric('_RENCA_WT_10FBS_DMSO-1234') -> true
%       is_last_part_numeric('_RENCA_WT_10FBS_DMSO400')   -> true
%       is_last_part_numeric('_RENCA_WT_10FBS_DMSO_test')  -> false

    % 1. Validate input: ensure we are working with a character vector or string
    if ~ischar(inputStr) && ~isstring(inputStr)
        error('Input must be a string or a character vector.');
    end

    % 2. Use regexp to search for one or more digits (\d+) at the end of the string ($).
    % 'match' returns the actual substring that matches the pattern.
    matches = regexp(inputStr, '\d+$', 'match');

    % 3. If the match array is not empty, it means a numeric sequence was found at the end.
    result = ~isempty(matches);
end



