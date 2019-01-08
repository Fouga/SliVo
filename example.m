% example of use to register a 2D slice to a 3D volume
clear all
close all
% load 3D CT volume 
volume_dir = './Data/3dCT/';
imds =imageDatastore(fullfile(volume_dir), 'FileExtensions', {'.tif'});
[Xsize Ysize] = size(readimage(imds,1));
Zsize = size(imds.Files,1);
CTvolume = zeros(Xsize, Ysize, Zsize, 'uint8');
parfor i = 1:Zsize
    CTvolume(:,:,i) = readimage(imds,i);
end

% load histology
histologyName = './Data/Resized_histology.png';
Histology = imread(histologyName);

registerSliceToVolume(CTvolume,Histology,'lower_limit',40, 'upper_limit', 400,...
    'calculate_features', 1, 'optimization',1,'angle',pi/18,'rotation_invariance',1,...
    'method2dregistration','NCC_NMI', 'radius',2.2);


% load histology for NOT invariant to rotation detection
histologyName = './Data/Rotated_histology.png';
Histology = imread(histologyName);

% localize a histological slide in a 3D data
registerSliceToVolume(CTvolume,Histology,'lower_limit',40, 'upper_limit', 400,...
    'calculate_features', 1, 'optimization',1,'angle',pi/18,'rotation_invariance',0,...
     'radius',2.2);
%% MRI volume
clear all
close all
% load 3D MRI volume 
volume_dir = './Data/3dMRI/';
imds =imageDatastore(fullfile(volume_dir), 'FileExtensions', {'.tif'});
[Xsize Ysize] = size(readimage(imds,1));
Zsize = size(imds.Files,1);
MRIvolume = zeros(Xsize, Ysize, Zsize, 'uint8');
parfor i = 1:Zsize
    MRIvolume(:,:,i) = readimage(imds,i);
end

% load histology
histologyName = './Data/Resized_histology.png';
Histology = imread(histologyName);
tic;
registerSliceToVolume(MRIvolume,Histology,'lower_limit',40, 'upper_limit', 400,...
    'calculate_features', 1, 'optimization',1,'angle',pi/18,'rotation_invariance',1,...
    'method2dregistration','NCC', 'radius',2.2,'feature_detector','SL1');
t = toc/60

% load histology for NOT invariant to rotation detection
histologyName = './Data/Rotated_histology.png';
Histology = imread(histologyName);

% localize a histological slide in a 3D data
tic;
registerSliceToVolume(MRIvolume,Histology,'lower_limit',40, 'upper_limit', 400,...
    'calculate_features', 1, 'optimization',1,'angle',pi/18,'rotation_invariance',0,...
    'method2dregistration','NCC_NMI', 'radius',2.2,'feature_detector','SL1');
t = toc/60
% 155 minutes
