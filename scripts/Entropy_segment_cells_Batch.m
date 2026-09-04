clear all
close all force hidden
clc
warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% settings to change for specific use cases
containing_folder = "C:\Users\LabPC2\Desktop\DIC_image_data_test\Raul data\2022-02-16\ImageJ processed";
experiment_name = "DIC_analysis_output";

% image format setting
image_type_format = '*.tif'; % ****must**** be in format '*.xxx'
% this is the type of images that the system will use
% currently only tested on tif and png files

% verbosity settings user interactivity
show_progress_bar = 1; % makes a figure that gives a progress bar
show_progress_in_terminal = 1; % displays current step to the terminal
show_segmentation_plot = 0; % shows the segmentation plot (briefly) once calculated

% image/gif export settings
% determines the settings for the data visualization/export
export_gif = 1; % export the data into a gif format - binary yes(1) no(0)
export_frames = 0; % export the frames of the gif aswell - binary yes(1) no(0)
full_resolution_gif_output = 1; % exports gif (large) or not (recommended)

% mask settings 
% this assumes the largest inital mask (in the middle) is
% the initial start for the scratch and ignores all the other "voids" of
% cells if there are any.
use_initial_largest_mask = 1; % Use the largest intial mask for baseof the segmentation - binary yes(1) no(0)

% small object filter settings
% this is a small object filter that if toggled removes small 'blobs' that
% usually crop up as time goes on. It is applied after masking from the
% initial large mask and helps reduse noise as the scratch is closing
filter_tiny_regions = 1;

%% Data processing section
ovr_dir = dir(containing_folder);
ovr_dir(ismember( {ovr_dir.name}, {'.', '..'})) = [];  %remove . and ..
dir_flags = [ovr_dir.isdir];
ovr_dir = ovr_dir(dir_flags);
exp_names = string(natsort({ovr_dir.name}))';

output_path = fullfile(pwd,'output',experiment_name);

disp(['Starting analysis for: ' char(experiment_name)])

mkdir(output_path);
for i = 1:length(exp_names)

    % locate and sort the images in the given folder
    this_exp = char(exp_names(i));
    this_exp_display = replace(this_exp,'_','-');
    imgs_dir_path = fullfile(containing_folder,this_exp);
    img_paths_struct = dir(fullfile(imgs_dir_path,image_type_format));
    sorted_img_names = natsort({img_paths_struct.name});
    num_imgs = length(sorted_img_names);

    if show_progress_in_terminal
        disp(['Experiment: ' num2str(i) ' of ' num2str(length(exp_names))])
        disp(['Experiment name:' char(this_exp)])
    end

    % preallocate for images and subsequent data processes
    imgs = cell(1,num_imgs);
    E_imgs = cell(1,num_imgs);
    mean_vals = zeros(1,num_imgs);

    if show_progress_in_terminal
        disp('--- Loading data')
    end

    % read in the images and normalize them between 0->1
    progress_bar = 0;
    for j = 1:num_imgs
        if show_progress_bar
            progress_bar = progressbar_function(j,num_imgs,progress_bar,{'Loading data',char(this_exp_display)});
        end
        temp_img = double(imread(fullfile(imgs_dir_path,sorted_img_names{j})));
        imgs{j} = temp_img/max(temp_img(:));
        mean_vals(j) = mean2(imgs{j});
        if isequal(j,num_imgs)
            if show_progress_bar
                close_progressbar(progress_bar)
            end
        end
    end

    mean_of_stack = mean(mean_vals);

    if filter_tiny_regions
        max_img_size = max(unique(cell2mat(cellfun(@(x) size(x), imgs, 'UniformOutput', false)'),'rows'),1);
        % max img size is given as rows, columns
        max_img_pixels = prod(max_img_size);
        filter_size = ceil(max_img_pixels/5000);
    end

    % normalize the total brightness of each image
    % to the average brightness of all the images
    for j = 1:num_imgs
        imgs{j} = imgs{j}*mean_of_stack/mean_vals(j);
    end

    % entropy filter the normalized images
    % before processing the images are self normalized per experiment to
    % ensure exposure is at least similar

    if show_progress_in_terminal
        disp('---- Processing data')
    end

    progress_bar = 0;
    for j = 1:num_imgs
        if show_progress_bar
            progress_bar = progressbar_function(j,num_imgs,progress_bar,{'Processing data',char(this_exp_display)});
        end
        % entropy filter the image for amount of disturbance
        E_imgs{j} = entropyfilt(imgs{j},true(25));
        if isequal(j,num_imgs)
            % combine all the entropy images
            large_E_array = imtile(E_imgs);
            % get a histogram of the all the counts of the entropy
            [N,edges] = histcounts(nonzeros(large_E_array));
            if show_progress_bar
                close_progressbar(progress_bar)
            end
        end
    end

    % this get the segmentation threshold (inflection point)
    % smooth the histogram
    smooth_N1 = smooth(N,10,'rloess');
    smooth_N2 = smooth(N,10);

    if max(abs(diff(smooth_N1))) > max(abs(diff(smooth_N2)))
        smooth_N = smooth_N2;
    else
        smooth_N = smooth_N1;
    end

    % find local minimus and prominance (derivative to a point)
    [TF,P] = islocalmin(smooth_N);
    % find the largest inflextion point
    inflection_point = find(P==max(P));
    % get inflection point as a number
    E_sep_point = edges(inflection_point);

    % processing/segmentation
    E_mask_1 = (E_imgs{1} < E_sep_point);
    inital_largest_mask = imfill(imgaussfilt(double(bwareafilt(E_mask_1,1)),3)>0,'holes');

    % this is all exporting and very little data processing
    % the amount_open is the segmentation data
    mkdir(fullfile(output_path,this_exp));
    if export_gif
        E_max = max(cellfun(@(x) max(x(:)), E_imgs));
    end

    if show_progress_in_terminal
        disp('----- Writing data')
    end

    amount_open = zeros(1,num_imgs);
    progress_bar = 0;
    for j = 1:num_imgs
        if show_progress_bar
            progress_bar = progressbar_function(j,num_imgs,progress_bar,{'Writing data',char(this_exp_display)});
        end
        % get the mask from the entropy seperation
        E_mask = (E_imgs{j}<E_sep_point);

        if filter_tiny_regions
            E_mask = bwareaopen(E_mask,filter_size,4);
        end

        % if using the initial largest mask
        % then mask of the resulting masks
        if use_initial_largest_mask
            E_mask = inital_largest_mask.*E_mask;
        end
        amount_open(j) = sum(sum(E_mask));
        if export_gif

            % create the tiled array for the GIF
            img_norm = imgs{j};
            img_norm(img_norm>1) = 1;

            uint8_img = uint8(img_norm*255);
            uint8_mask = uint8(E_mask*255);
            uint8_E_img = uint8((E_imgs{j}/E_max)*255);

            rgb_img = cat(3,uint8_img,uint8_img,uint8_img);
            rgb_mask = cat(3,uint8_mask,uint8_mask,uint8_mask);
            rgb_E = cat(3,uint8_E_img,uint8_E_img,uint8_E_img);

            out_label = labeloverlay(uint8_img,E_mask,'Colormap','autumn','Transparency',0.75);
            if length(size(out_label)) < 3
                out_label = cat(3,out_label,out_label,out_label);
            end
            out_img = [rgb_img,rgb_mask;
                rgb_E,out_label];

            if ~full_resolution_gif_output
                out_img = imresize(out_img,1000/min(max_img_size*2));
            end

            [A,map] = rgb2ind(out_img,256);
            if j == 1
                imwrite(A,map,fullfile(output_path,this_exp,[this_exp '.gif']), ...
                    'gif','LoopCount',Inf,'DelayTime',1/24);
            else
                imwrite(A,map,fullfile(output_path,this_exp,[this_exp '.gif']), ...
                    'gif','WriteMode','append','DelayTime',1/24);
            end

            if export_frames
                imwrite(out_img,fullfile(output_path,this_exp,[num2str(j) '.jpg']))
            end
        end
        if isequal(j,num_imgs)
            if show_progress_bar
                close_progressbar(progress_bar)
            end
        end
    end

    time_points = 1:length(img_paths_struct);

    grad_amount_open = gradient(amount_open)';

    header = ["time point","Amount open (pixels)","Norm open","Speed of closing","Norm close speed"];

    out = [time_points',amount_open',amount_open'/max(amount_open),grad_amount_open,grad_amount_open/max(abs(grad_amount_open))];

    out = cell2table(num2cell(out),'VariableNames',header);

    plot_data_simple(this_exp,edges,smooth_N,inflection_point,amount_open,fullfile(output_path,this_exp,['_' this_exp '.png']),show_segmentation_plot)
    writetable(out,fullfile(output_path,this_exp,['_' this_exp '_data.csv']))

    clear amount_open imgs E_imgs

end

close all
disp('Finished processing data')

%% functions block

function plot_data_simple(this_exp,edges,smooth_N,inflection_point,amount_open,path_to_export,show_segmentation_plot)

e2 = edges(1:end-1);
if show_segmentation_plot
    f = figure('Units','normalized','Position',[0,0.4,0.5,0.5]);
else
    f = figure('Units','normalized','Position',[0,0.4,0.5,0.5],'Visible','off');
end
sgtitle(this_exp,'Interpreter','None');
subplot(1,2,1);
plot(e2,smooth_N,'b-',e2(inflection_point),smooth_N(inflection_point),'r*');
text(e2(inflection_point),smooth_N(inflection_point),{'inflex point', num2str(round(smooth_N(inflection_point)))});
title('Dual peak infereance');
xlabel('total entropy');
ylabel('measured entropy');

subplot(1,2,2);
plot(amount_open/max(amount_open));
title('Normalized amount of open space');
xlabel('Hours');

% axis;
axis('square');

ax = gca;                 % Get current axes
ax.Toolbar.Visible = 'off'; % Hide the toolbar completely

saveas(f,path_to_export)

close(f)

end


function progress_bar = progressbar_function(i,num_samples,progress_bar,title)

progress_ratio = i/num_samples;

if isequal(progress_bar,0)
    progress_bar = waitbar(progress_ratio,title);
else
    progress_bar = waitbar(progress_ratio,progress_bar,title);
end

end

function close_progressbar(progress_bar)

close(progress_bar)

end
