clear
clc; close all
X_size=600;
Y_size=700;
%read histology
input_hist='/media/a31cf190-6597-45ee-8fc4-223ef970f9d5/natalia_HardDrive/Dental_project/data_anya/13/';
Hist=imread([input_hist 'DS Ilgenstein_SPIEE_Mac_13.10-1_00007.tif']);
Hist=imcomplement(imresize(rgb2gray(Hist),[X_size Y_size]));

%read first CT slice
n1=100; n2=890;
X_size=764;
Y_size= 764;
Z_size = 969;
read_dir= '/media/a31cf190-6597-45ee-8fc4-223ef970f9d5/natalia_HardDrive/Dental_project/data_anya/13/';
output_dir = '/media/a31cf190-6597-45ee-8fc4-223ef970f9d5/natalia_HardDrive/Project_multi_reg_Results/Dense_match_figures/13/';
filename=[read_dir 'bmc02_rawbin2_2b_0x1527_0y1527_0z1938.dat'];
  CT_file=fopen(filename, 'r');
    header=length(fgets(CT_file));
    fseek(CT_file,header,'bof');
    data_CT = fread(CT_file,X_size*Y_size*Z_size,'uint8=>uint8');
  fclose(CT_file);
image_global = reshape(data_CT, X_size, Y_size, Z_size); 
clear data_CT

% Load coordinates
surf_read_dir='/media/a31cf190-6597-45ee-8fc4-223ef970f9d5/natalia_HardDrive/Project_multi_reg_Results/MATCH_2_nearest_neighbor/13/SURF_results/';
im=1
    filename_coord=sprintf('Coord_3D_SURF-%i_0.8',im);
    load([surf_read_dir filename_coord '.mat']);
    [B,inliers]=matching_plane(surf_3D, n1,n2,X_size,Y_size);

    % show cut image
    [x_mesh, y_mesh]=meshgrid(1:X_size,1:Y_size);
    Z=-(y_mesh.*B(1) + x_mesh.*B(2) + B(4))/B(3); %set of coordinates for tilting and rotating
    %%%%% must change x and y
    Z_ = Z - floor(min(min(Z)));
    CT = interp3(single(image_global(:,:,floor(min(min(Z))):ceil(max(max(Z))))),x_mesh,y_mesh,Z_); %cut a slice from 3D volume
    CT=flipdim(CT,2);
    figure, subplot(3,2,1),
    imshow(uint8(CT));




% input_CT='/media/a31cf190-6597-45ee-8fc4-223ef970f9d5/natalia_HardDrive/Project_multi_reg_Results/Dense_match_figures/13/';
% CT=imread([input_CT 'CT13_histology_cut_01.tif']);
% CT=imresize(rgb2gray(CT),[X_size Y_size]);

%SURF keypoints
    Options.tresh=0.0001; %Hessian response treshold (default 0.0002)
    Ipts1=OpenSurf(CT,Options); %key points of the image
    Ipts2=OpenSurf(Hist,Options);
    des1_surf = reshape([Ipts1.descriptor],64,[]); 
    des1=des1_surf'; % matrix K by 64
    des2_surf = reshape([Ipts2.descriptor],64,[]);
    des2=des2_surf';
%     loc1=[reshape([Ipts1.y],1,[]);reshape([Ipts1.x],1,[])]; %locs of SIFT [column, row]
%     loc2=[reshape([Ipts2.y],1,[]);reshape([Ipts2.x],1,[])];

    loc1=[reshape([Ipts1.x],1,[]); reshape([Ipts1.y],1,[])]; %locs of SIFT [column, row]
    loc2=[reshape([Ipts2.x],1,[]); reshape([Ipts2.y],1,[])];
    


        distRatio = 0.9;   
k=1;
% For each descriptor in the first image, select its match to second image.
    des2t = des2';                          % Precompute matrix transpose
    dotprod_all=des1*des2t;
    for i = 1 : size(des1,1)
       %dotprods = des1(i,:) * des2t;        % Computes vector of dot products
       dotprods=dotprod_all(i,:);
       [vals,indx] = sort(acos(dotprods));  % Take inverse cosine and sort results

       % Check if nearest neighbor has angle less than distRatio times 2nd.
       if (vals(1) < distRatio * vals(2))
          match2(i) = indx(1);
          coor_hist(:,k)=loc2(:,indx(1));%Hist
          k=k+1;
          
       else
          match2(i) = 0;
       end
     
    end
    coor_CT=loc1(:,match2>0);
    number = sum(match2 > 0);
    
  
    t=0.001;
    [H, inliers] = ransacfithomography(coor_hist, coor_CT, t);
    %coor_CT=H*coor_hist, transform Hist so that we get CT
    tform=maketform('projective',H');
    imT = imtransform(Hist,tform,'bicubic');
    subplot(3,2,3), imshow(Hist,[])
    subplot(3,2,2), imshow(imT,[])
    
    
    [M N]=size(imT);
    imTpadd=padarray(imT, [ceil((X_size-M)) ceil((Y_size-N))],'post');
    CC=normxcorr2( imT, CT);
    figure, surf(CC), shading flat
    
    RC=corr(im2double(imTpadd),im2double(CT));
    figure, surf(RC), shading flat

    %%%%%%%%%%    
   
%    Find peak in cross-correlation.

% [ypeak, xpeak] = find(RC==max(RC(:)));
% 
% %Account for the padding that normxcorr2 adds.
% 
% yoffSet = ypeak-size(onion,1);
% xoffSet = xpeak-size(onion,2);
% 
% %Display matched area.
% 
% hFig = figure;
% hAx  = axes;
% imshow(peppers,'Parent', hAx);
% imrect(hAx, [xoffSet, yoffSet, size(onion,2), size(onion,1)]);
