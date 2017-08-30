function Slice_ransac = fit_plane(Data3D,Slice2D,options)
 
% load coordinates
load([options.folder_matches 'FeatureCoordinates_3D.mat']);
[X_size Y_size Z_size] = size(Data3D);

% calculate normal vector
[B,~]=matching_plane(FeatureCoordinates_3D, options);
save([options.folder_matches 'Normal_vec_toPlane.mat'],'B','-mat')

% show cut image
[x_mesh, y_mesh]=meshgrid(1:X_size,1:Y_size);
Z=-(y_mesh.*B(1) + x_mesh.*B(2) + B(4))/B(3);

Slice_ransac = interp3(single(Data3D),x_mesh,y_mesh,Z); %cut a slice from 3D volume
imwrite(uint8(Slice_ransac), [options.folder_matches 'FoundMatch_in_3D.tif'],'Compression','None' )

figure,
subplot(1,2,1), imshow(Slice_ransac,[]), title('Found slice in 3D volume')
subplot(1,2,2), imshow(Slice2D,[]), title('Given 2D slice')

