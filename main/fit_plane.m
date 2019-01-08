function Slice_ransac = fit_plane(Data3D,Slice2D,options)
% Finds a plane in the 3D point cloud with the highest number of inliers
% or, in other words, fits a plane to the highest density of the points.
% The 3D could is a cloud of matching feature points between an input image
% (Slice2D) and each image in the 3D volume (Data3D). The details of the
% Ransac algorithm can be found here: 
% Fischler, M.A. & Bolles, R.C. (1981) Random sample consensus: a
% paradigm for model fitting with applications to image analysis and
% automated cartography. Commun. ACM 24, 381–395.
%
% 
% Usage:  Slice_ransac = fit_plane(Data3D,Slice2D,options)
%
% Input: Data3D - is a XxYxZ 3D uint8 array. The image volume data should be
%                 passed as a stack of 2D images 
%       Slice2D - is a 2D uint8 image which you want to localize in the
%                 Data3D, for example, histological slice.
%       options - is a structure containing various parameter values needed by
%                slice to volume registration 
%
% Output: Slice_ransac - is a 2D image from the Data3D found as a fit to the highest density of
%                        3D correspondences cloud.  
%         The output results are saved in the otput folder: options.folder_destination
%
% See also: matching_plane, MatlabFns
%
% From the project, (https://github.com/Fouga/).
% Author: 2019 Natalia Chicherova.
%
% If you use this code please cite my paper:
% Chicherova, N., Hieber, S.E., Khimchenko, A., Bikis, C., Mü ller, B.,
% Cattin, P.C., 2018. Automatic deformable registration of histological
% slides to μCT volume data. Journal of Microscopy 271, 49–61.

% load coordinates
load(fullfile(options.folder_matches, 'FeatureCoordinates_3D.mat'),'FeatureCoordinates_3D');
[X_size Y_size Z_size] = size(Data3D);

% calculate normal vector
[B,~]=matching_plane(FeatureCoordinates_3D, options);
save(fullfile(options.folder_matches, 'Normal_vec_toPlane.mat'),'B','-mat')

% show cut image
[x_mesh, y_mesh]=meshgrid(1:Y_size,1:X_size);
Z=-(y_mesh.*B(1) + x_mesh.*B(2) + B(4))/B(3);

Slice_ransac = interp3(single(Data3D),x_mesh,y_mesh,Z); %cut a slice from 3D volume
imwrite(uint8(Slice_ransac), fullfile(options.folder_matches, 'FoundMatch_in_3D.tif'),'Compression','None' )

figure,
subplot(1,2,1), imshow(Slice_ransac,[]), title('Found slice in 3D volume')
subplot(1,2,2), imshow(Slice2D,[]), title('Given 2D slice')

