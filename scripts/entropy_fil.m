
imgs_dir = "Y:\Users\Raul Castro\Microscopes\Olympus Spining Disk\2022-02-16\ImageJ processed\RENCA_A28_5FBS_3HAA_1";

img_paths = dir(fullfile(imgs_dir,'*.tif'));

for i = 1:length(img_paths)
    
    imgs{i} = imread(fullfile(imgs_dir,img_paths(i).name));
    
    E_imgs{i} = entropyfilt(imgs{i});
    
end

large_E_array = imtile(E_imgs);
[N,edges] = histcounts(nonzeros(large_E_array));

smooth_N = smooth(N);

e2 = edges(1:end-1);
N_1000 = N>1000;