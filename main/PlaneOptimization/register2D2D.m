function Slice2D_registered2D = register2D2D(Slice2D,Slice_ransac,options)
% Calculates feature points and matches the descriptor vectors of these points
% for a input 2D slice (e.g. histology) and each image in the 3D volume
% (e.g. CT, MRI).
% For invariant to rotation feature detector there are 3 methods can be
% used to aling the Slice2D with Slice_ransac: NCC, NMI, NCC_NMI.
%
% NCC - can be applied to images with the same SCALE and clear background 
%       because it uses Otsu thersholding to segment the ROI. It compemsates
%       for rotation difference by miximizing NCC. 
% NMI - can be applied to images with the same SCALE. It compensates for
%       rotation and translation difference between images.
% NCC_NMI - can be applied to images with unequal scale. It, first, finds
%           rotation and then compensates for slight similarity difference.
% 
% For NOT invariant to rotation feature detector the code compensates only
% for transaltion, the scale should be adjusted before. 
% 
% Usage:  OPTIONS =  find_features(Data3D, Slice2D, options)
%
% Input: Slice2D - is a 2D uint8 image which you want to localize in the
%                 Data3D, for example, histological slice.
%        Slice_ransac - is a 2D image that matches best the Slice2D in the 
%                       volume data (Data3D). It is found by using fitting
%                       a plane into the 3D feature point cloud.
%       options - is a structure containing various parameter values needed by
%                slice to volume registration 
%
% Output: OPTIONS a structure containing various parameter values needed by
%         slice to volume registration
%         The output feauture vectors and matched correspondences are saved
%         in the otput folder: options.folder_destination
%
% See also: fit_plane, find_features
%
% Note: the 2D-2D registration based on Ransac Homography matrix and correspondences is not
% implemented.
%
% From the project, (https://github.com/Fouga/).
% Author: 2019 Natalia Chicherova.
%
% If you use this code please cite my paper:
% Chicherova, N., Hieber, S.E., Khimchenko, A., Bikis, C., Mueler, B.,
% Cattin, P.C., 2018. Automatic deformable registration of histological
% slides to μCT volume data. Journal of Microscopy 271, 49–61.
   
% get rid of NaNs
Slice2D = double(Slice2D);
[X_size Y_size] = size(Slice_ransac);
if (isnan(Slice_ransac(:))==0)
else
   for i = 1:X_size
        for j = 1:Y_size
            if (isnan(Slice_ransac(i,j)) ==1)
                Slice_ransac(i,j) = 1;
            end
        end
    end
end
Slice_ransac = double(Slice_ransac);    
if options.rotation_invariance ==1
    switch options.method2dregistration
        case 'NCC'
            % use NCC to find rotation
            level = graythresh(Slice2D);
            slice2D_BW = im2bw(Slice2D,level);
            level = graythresh(Slice_ransac);
            Slice_ransac_BW = im2bw(Slice_ransac,level);

             m = 1; CC1 = zeros(1,360);

            for ang = 0:359
                im_rotate = imrotate(slice2D_BW,ang,'crop');
                c = normxcorr2_general(im_rotate,Slice_ransac_BW);
                CC1(m) = max(c(:));
                m = m+1;
            end
            % figure,plot(CC1)
            [max_value ind_val] = max(CC1);
            Slice2D_registered2D = imrotate(Slice2D,ind_val,'crop');
        case 'NMI'
            registration_type = 'rigid';
            % metric = registration.metric.MeanSquares;
            metric = registration.metric.MattesMutualInformation;
            optimizer = registration.optimizer.OnePlusOneEvolutionary;
            %optimizer = registration.optimizer.RegularStepGradientDescent;
            Slice2D_registered2D = imregister(Slice2D,Slice_ransac,registration_type,optimizer,metric);
        case 'NCC_NMI'
            % use NCC to find rotation
            level = graythresh(Slice2D);
            slice2D_BW = im2bw(Slice2D,level);
            level = graythresh(Slice_ransac);
            Slice_ransac_BW = im2bw(Slice_ransac,level);

             m = 1; CC1 = zeros(1,360);

            for ang = 0:359
                im_rotate = imrotate(slice2D_BW,ang,'crop');
                c = normxcorr2_general(im_rotate,Slice_ransac_BW);
                CC1(m) = max(c(:));
                m = m+1;
            end
            % figure,plot(CC1)
            [max_value ind_val] = max(CC1);
            Slice2D_rotated = imrotate(Slice2D,ind_val,'crop');
     
            registration_type = 'similarity';
            % metric = registration.metric.MeanSquares;
            metric = registration.metric.MattesMutualInformation;
            optimizer = registration.optimizer.OnePlusOneEvolutionary;
            %optimizer = registration.optimizer.RegularStepGradientDescent;
            Slice2D_registered2D = imregister(Slice2D_rotated,Slice_ransac,registration_type,optimizer,metric);
    end
            
    figure, subplot(1,2,1), imshow(uint8(Slice_ransac),[]), title('Found slice in 3D volume')
            subplot(1,2,2), imshow(uint8(Slice2D_registered2D),[]), title('Registered in 2D given 2D slice')

else
    % when the slices are already pre-alinged in case of not invariant to rotation feature detector
    % we need to compensate only for slight shift differences
    slice2D_rotated = Slice2D;
        
    % do only translation, because the slices already rotated
    registration_type = 'translation';%'similarity'
    % metric = registration.metric.MeanSquares;
    metric = registration.metric.MattesMutualInformation;
    optimizer = registration.optimizer.OnePlusOneEvolutionary;
    %optimizer = registration.optimizer.RegularStepGradientDescent;
    slice2D_translate = imregister(slice2D_rotated,Slice_ransac,registration_type,optimizer,metric);
    Slice_ransac = uint8(Slice_ransac);
    slice2D_translate = uint8(slice2D_translate);
    figure, subplot(1,3,1), imshow(slice2D_translate,[]), title('Translated given 2D slice');
                subplot(1,3,2), imshow(single(slice2D_translate)-single(Slice_ransac),[]), title('Image difference');
                    subplot(1,3,3), imshow(Slice_ransac,[]), title('Found slice in 3D volume');
    Slice2D_registered2D = slice2D_translate;    
    
end
    