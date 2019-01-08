function options = find_features(Data3D, Slice2D, options)
% Calculates feature points and matches the descriptor vectors of these points
% for a input 2D slice (e.g. histology) and each image in the 3D volume (e.g. CT, MRI)
% 
% Usage:  OPTIONS =  find_features(Data3D, Slice2D, options)
%
% Input: Data3D - is a XxYxZ 3D uint8 array. The image volume data should be
%                 passed as a stack of 2D images 
%       Slice2D - is a 2D uint8 image which you want to localize in the
%                 Data3D, for example, histological slice.
%       options - is a structure containing various parameter values needed by
%                slice to volume registration 
%
% Output: OPTIONS a structure containing various parameter values needed by
%         slice to volume registration
%         The output feauture vectors and matched correspondences are saved
%         in the otput folder: options.folder_destination
%
% See also: match_features, save_keys, feature_load
%
% From the project, (https://github.com/Fouga/).
% Author: 2019 Natalia Chicherova.
%
% If you use this code please cite my paper:
% Chicherova, N., Hieber, S.E., Khimchenko, A., Bikis, C., Müuller, B.,
% Cattin, P.C., 2018. Automatic deformable registration of histological
% slides to μCT volume data. Journal of Microscopy 271, 49–61.

[X_size Y_size Z_size] = size(Data3D);
options.size =  [X_size Y_size Z_size];

if isempty(options.folder_destination) || options.calculate_features == false
    if ~exist('Results_of_2D3D_Registration')
        mkdir('Results_of_2D3D_Registration')
    end
    options.folder_destination = fullfile('Results_of_2D3D_Registration');
end


if ~exist(fullfile(options.folder_destination, 'Features','2Dslice'))
    mkdir(fullfile(options.folder_destination, 'Features','2Dslice'))
end
options.folder_features = fullfile(options.folder_destination, 'Features'); 


if isempty(options.lower_limit)
    options.lower_limit = 1;
end

if isempty(options.upper_limit)
    options.upper_limit = size(Data3D,3);
end


if options.calculate_features == true

    output_Slice2D=fullfile(options.folder_features, '2Dslice', 'Feature_2Dslice.txt');

    % calculate SURF feature points and save them in a txt file
    save_keys(output_Slice2D,Slice2D,options);

    % calculate SURF/SL1 features for every image in a 3D stack
    if ~exist(fullfile(options.folder_features, '3Dstack'))
        mkdir(fullfile(options.folder_features, '3Dstack'))
    end
    output_3Ddata=fullfile(options.folder_features, '3Dstack');

    for j=options.lower_limit:options.upper_limit 
        filename=sprintf('Feature_3Dslice%i.txt',j);
        output_file=fullfile(output_3Ddata, filename);
        save_keys(output_file,Data3D(:,:,j),options);
        fprintf('Feaures for a slice %i\n',j);
    end
end
    

%% Match the features
disp('Matching features...')
[des2D, locs2D]=feature_load(fullfile(options.folder_features, '2Dslice','Feature_2Dslice.txt'), options.feature_detector);

FeatureCoordinates_3D=cell(1,length(options.lower_limit:options.upper_limit));

for j=options.lower_limit:options.upper_limit 
    filename_3Dslice=sprintf('Feature_3Dslice%i.txt',j);
    [des3D,locs3D]=feature_load(fullfile(options.folder_features, '3Dstack', filename_3Dslice),options.feature_detector);
    [number,coor2Dslice,coor3D]=match_features(des2D,locs2D,des3D,locs3D,options); 
    FeatureCoordinates_3D{j}=coor3D(:,1:2);
    fprintf('***************Number of matched points for a slice %i is %i*************\n',j,number);
end

if ~exist(fullfile(options.folder_destination, 'Matches'))
    mkdir(fullfile(options.folder_destination, 'Matches'))
end
options.folder_matches =fullfile(options.folder_destination, 'Matches');        

save(fullfile(options.folder_matches, 'FeatureCoordinates_3D.mat'), 'FeatureCoordinates_3D');

