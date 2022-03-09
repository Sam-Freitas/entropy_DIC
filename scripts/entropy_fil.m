clear all
close all force hidden
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% imgs_dir = "Y:\Users\Raul Castro\Microscopes\Olympus Spining Disk\2022-02-16\ImageJ processed\RENCA_A28_5FBS_3HAA_1";
% containing_folder = "/Volumes/Sutphin server/Users/Raul Castro/Microscopes/Olympus Spining Disk/2022-02-16/ImageJ processed";
containing_folder = "Y:\Users\Raul Castro\Microscopes\Olympus Spining Disk\2022-02-16\ImageJ processed";
ovr_dir = dir(containing_folder);
ovr_dir(ismember( {ovr_dir.name}, {'.', '..'})) = [];  %remove . and ..

export_gif = 1;
export_frames = 0;

use_inital_largest_mask = 1;

mkdir('output');
for i = 1:length(ovr_dir)
        
    this_exp = ovr_dir(i).name;
    this_exp_display = replace(this_exp,'_','-');
    
    imgs_dir = fullfile(ovr_dir(1).folder,this_exp);
    
    img_paths = dir(fullfile(imgs_dir,'*.tif'));
    sorted_img_paths = natsort({img_paths.name});
    
    num_imgs = length(sorted_img_paths);
    
    imgs = cell(1,num_imgs);
    E_imgs = cell(1,num_imgs);
    
    progress_bar = 0;
    for j = 1:num_imgs
        progress_bar = progressbar_function(j,num_imgs,progress_bar,{'Loading data',char(this_exp_display)});
        imgs{j} = double(imread(fullfile(imgs_dir,sorted_img_paths{j})));
        imgs{j} = imgs{j}/max(imgs{j}(:));
        mean_vals(j) = mean2(imgs{j});
        if isequal(j,num_imgs)
            close_progressbar(progress_bar)
        end
    end
    
    mean_of_stack = mean(mean_vals);
    
    progress_bar = 0;
    for j = 1:num_imgs
        progress_bar = progressbar_function(j,num_imgs,progress_bar,{'Processing data',char(this_exp_display)});       
        % entropy filter the image for amount of disturbance
        imgs{j} = imgs{j}*mean_of_stack/mean_vals(j);
        E_imgs{j} = entropyfilt(imgs{j},true(25));
        if isequal(j,num_imgs)
            % combine all the entropy images
            large_E_array = imtile(E_imgs);
            % get a histogram of the all the counts of the entropy
            [N,edges] = histcounts(nonzeros(large_E_array));
            close_progressbar(progress_bar)
        end
    end
    
    % this get the segmentation threshold  
    % smooth the histogram 
    smooth_N = smooth(N,10,'rloess');
    % find local minimus and prominance (derivative to a point)
    [TF,P] = islocalmin(smooth_N);
    % find the largest inflextion point
    inflection_point = find(P==max(P));
    % get inflection point as a number
    E_sep_point = edges(inflection_point);
    
    % processing 
    E_mask_1 = (E_imgs{1} < E_sep_point);
    
    inital_largest_mask = imfill(imgaussfilt(bwareafilt(E_mask_1,1)*2.5,3)>0,'holes');
    
    mkdir(fullfile(pwd,'output',this_exp));
    progress_bar = 0;
    for j = 1:num_imgs
        progress_bar = progressbar_function(j,num_imgs,progress_bar,{'Writing data',char(this_exp_display)});
        
        amount_open(j) = sum(sum(E_imgs{j}<E_sep_point));
        if export_gif
            E_mask = (E_imgs{j} < E_sep_point);
            
            if use_inital_largest_mask
                E_mask = inital_largest_mask.*E_mask;
            end

            img_norm = imgs{j};
            img_norm(img_norm>1) = 1;

            uint8_img = uint8(img_norm*255);
            uint8_mask = uint8(E_mask*255);
            uint8_E_img = uint8(E_imgs{j}/max(nonzeros(E_imgs{j}))*255);

            rgb_img = cat(3,uint8_img,uint8_img,uint8_img);
            rgb_mask = cat(3,uint8_mask,uint8_mask,uint8_mask);
            rgb_E = cat(3,uint8_E_img,uint8_E_img,uint8_E_img);

            out_label = labeloverlay(uint8(img_norm*255),E_mask,'Colormap','autumn','Transparency',0.75);
            if length(size(out_label)) < 3
                out_label = cat(3,out_label,out_label,out_label);
            end
            out_img = [rgb_img,rgb_mask;rgb_E,out_label];
            
            if isequal(j,1)
                h = imshow(out_img);
                gif(fullfile(pwd,'output',this_exp,[this_exp '.gif']),'DelayTime',1/24,'overwrite',true)
            else
                set(h,'CData',out_img)
                gif
            end
            if export_frames
                imwrite(out_img,fullfile(pwd,'output',this_exp,[num2str(j) '.jpg']))
            end
        end
        if isequal(j,num_imgs)
            close_progressbar(progress_bar)
        end
    end
    
    time_points = 1:length(img_paths);
    
    grad_amount_open = gradient(amount_open)';
    
    header = ["time point","Amount open (pixels)","Norm open","Speed of closing","Norm close speed"];
    
    out = [time_points',amount_open',amount_open'/max(amount_open),grad_amount_open,grad_amount_open/max(abs(grad_amount_open))];
    
    out = cell2table(num2cell(out),'VariableNames',header);
    
    plot_data_simple(this_exp,edges,smooth_N,inflection_point,amount_open,fullfile(pwd,'output',this_exp,['_' this_exp '.png']))
    writetable(out,fullfile(pwd,'output',this_exp,['_' this_exp '_data.csv']))
        
    clear amount_open imgs E_imgs
    
end


function plot_data_simple(this_exp,edges,smooth_N,inflection_point,amount_open,path_to_export)

    e2 = edges(1:end-1);
        
    f = figure('Units','normalized','Position',[0,0.4,0.5,0.5]);
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
    
    axis;
    axis('square');
    
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
