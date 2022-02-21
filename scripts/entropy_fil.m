
% imgs_dir = "Y:\Users\Raul Castro\Microscopes\Olympus Spining Disk\2022-02-16\ImageJ processed\RENCA_A28_5FBS_3HAA_1";
containing_folder = "/Volumes/Sutphin server/Users/Raul Castro/Microscopes/Olympus Spining Disk/2022-01-25/ImageJ processed";

imgs_dir = "/Volumes/Sutphin server/Users/Raul Castro/Microscopes/Olympus Spining Disk/2022-01-25/ImageJ processed/RENCA_A28_DFO_27";
img_paths = dir(fullfile(imgs_dir,'*.tif'));


for i = 1:length(img_paths)
    
    imgs{i} = imread(fullfile(imgs_dir,img_paths(i).name));
    
    E_imgs{i} = entropyfilt(imgs{i},true(15));
    
end

large_E_array = imtile(E_imgs);
[N,edges] = histcounts(nonzeros(large_E_array));

smooth_N = smooth(N,10,'rloess');

[TF,P] = islocalmin(smooth_N);

inflection_point = find(P==max(P));

e2 = edges(1:end-1);
N_1000 = N>1000;

figure;
plot(e2,smooth_N,'b-',e2(inflection_point),smooth_N(inflection_point),'r*')

E_sep_point = e2(inflection_point);

for i = 1:length(img_paths)
    
    amount_open(i) = sum(sum(E_imgs{i}<E_sep_point));
    
end

figure;
plot(amount_open/max(amount_open))
