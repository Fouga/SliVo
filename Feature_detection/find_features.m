function options = find_features(Data3D, Slice2D, options)

[X_size Y_size Z_size] = size(Data3D);
options.size =  [X_size Y_size Z_size];

if isempty(options.folder_destination) || options.calculate_features == false
    if ~exist('./Results_of_2D3D_Registration/')
        mkdir('./Results_of_2D3D_Registration/')
    end
    options.folder_destination = './Results_of_2D3D_Registration/';
end


if ~exist([options.folder_destination 'Features/2Dslice/'])
    mkdir([options.folder_destination 'Features/2Dslice/'])
end
options.folder_features = [options.folder_destination 'Features/']; 

if isempty(options.feature_detector)
    options.feature_detector = 'SURF';
end

if isempty(options.lower_limit)
    options.lower_limit = 1;
end

if isempty(options.upper_limit)
    options.upper_limit = size(Data3D,3);
end


if options.calculate_features == true

    output_surf=[options.folder_features '2Dslice/Feature_2Dslice.txt'];

    % calculate SURF feature points and save them in a txt file
    thresh = 0.0002;
    save_keys(output_surf,Slice2D,options.feature_detector,thresh);

    % calculate SURF features for every image in a 3D stack
    if ~exist([options.folder_features '/3Dstack/'])
        mkdir([options.folder_features '/3Dstack/'])
    end
    output_3Ddata=[options.folder_features '/3Dstack/'];

    for j=options.lower_limit:options.upper_limit 
        filename=sprintf('Feature_3Dslice%i.txt',j);
        output_file=[output_3Ddata filename];
        save_keys(output_file,Data3D(:,:,j),options.feature_detector,thresh);
        fprintf('Feaures for a slice %i\n',j);
    end
end
    

%% Match the features
disp('Matching features...')
[des2D, locs2D]=feature_load([options.folder_features '2Dslice/Feature_2Dslice.txt'], options.feature_detector);

FeatureCoordinates_3D=cell(1,length(options.lower_limit:options.upper_limit));

for j=options.lower_limit:options.upper_limit 
    filename_3Dslice=sprintf('Feature_3Dslice%i.txt',j);
    [des3D,locs3D]=feature_load([options.folder_features '/3Dstack/' filename_3Dslice],options.feature_detector);
    [number,coor2Dslice,coor3D]=match_features(des2D,locs2D,des3D,locs3D); 
    FeatureCoordinates_3D{j}=coor3D(:,1:2);
%     fprintf('Number of matched point for a slice %i is %i\n',j,number);
end

if ~exist([options.folder_destination 'Matches/'])
    mkdir([options.folder_destination 'Matches/'])
end
options.folder_matches = [options.folder_destination 'Matches/'];        

save([options.folder_matches 'FeatureCoordinates_3D.mat'], 'FeatureCoordinates_3D');

