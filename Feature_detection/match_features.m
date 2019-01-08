function [number,coor2Dslice,coor3D]=match_features(des2D,locs2D,des3D,locs3D,options)
% Matches feature points detector using second nearest neighbour in case of
% SURF detector and L1-norm for SL1 detector. More information about the
% truncated L1 norm outlier rejection can be found here: 
% Ask, E., Enqvist, O., Svam, L., Kahl, F., Lippolis, G., 2014. Tractable
% and reliable registration of 2D point sets, in: Computer Vision–ECCV
% 2014. Springer, pp. 393–406.
%
% I am deeply grateful to Benedikt Bitterli for writting the speed up
% versions of loss_function_fast and scripts for SS feature detection. 
%
% See also: find_features, feature_load, loss_function_fast
%
%
% From the project, (https://github.com/Fouga/).
% Author: 2019 Natalia Chicherova.
%
% If you use this code please cite my paper:
% Chicherova, N., Hieber, S.E., Khimchenko, A., Bikis, C., Mueler, B.,
% Cattin, P.C., 2018. Automatic deformable registration of histological
% slides to μCT volume data. Journal of Microscopy 271, 49–61.

method = options.feature_detector;
k = 1;
if (strcmp(method,'SURF'))

    distRatio = 0.8;   
    
    % For each descriptor in the first image, select its match to second image.
    des2t = des2D';                          % Precompute matrix transpose
    dotprod_all=des3D*des2t;
    for i = 1 : size(des3D,1)
       %dotprods = des1(i,:) * des2t;        % Computes vector of dot products
       dotprods=dotprod_all(i,:);
       [vals,indx] = sort(acos(dotprods));  % Take inverse cosine and sort results

       % Check if nearest neighbor has angle less than distRatio times 2nd.
       if (vals(1) < distRatio * vals(2))
          match(i) = indx(1);
          coor2Dslice(k,:)=locs2D(indx(1),:);
          k = k+1;
       else
          match(i) = 0;
       end
     
    end
    coor3D=locs3D(match>0,:);%coordinates of matched features
    
    number = sum(match > 0);
    if number == 0
        coor3D = double.empty(0,2);
        coor2Dslice = double.empty(0,2);
    end
    
elseif (strcmp(method,'SL1'))
    
    % first outlier rejection
    distRatio = 0.7;   
    
    % For each descriptor in the first image, select its match to second image.
    des2t = des2D';                          % Precompute matrix transpose
    dotprod_all=des3D*des2t;
    for i = 1 : size(des3D,1)
        dotprods=dotprod_all(i,:);
        % Take inverse cosine and sort results
        [vals,indx] = sort(acos(dotprods));  
        % Check if nearest neighbor has angle less than distRatio times 2nd.
        if (vals(1) < distRatio * vals(2))
           match(i) = indx(1);
           coor2Dslice(k,:)=locs2D(indx(1),:);
           k = k+1;
       else
          match(i) = 0;
       end
     
    end
    
    %coordinates of matched features
    coor3D=locs3D(match>0,:);  % reference CT
%     coorCT_1 = coor3D;
%     coor_hist_1 = coor2Dslice;
    if isempty(coor3D)
        coor3D = double.empty(0,2);
        coor2Dslice = double.empty(0,2);
    else
        % second outlier rejection
        X_size = options.size(1);
        Y_size = options.size(2);
        coor1 = coor2Dslice';% float
        coor2 = coor3D';% reference

        [R,alpha,t] = loss_function_fast(coor1-repmat([Y_size; X_size]/2,1,size(coor1,2)),coor2-repmat([Y_size; X_size]/2,1,size(coor1,2)),10);
        if isnan(R)==0
            Tr = repmat([Y_size/2 X_size/2 0]',1,size(coor1,2));

            coor1_flip = [coor1; ones(1,size(coor1,2))];
            coor_T = R*(coor1_flip-Tr)+Tr+repmat([t;0],1,size(coor1,2));

            %% Outlier rejection
            epsi = 10;
            dist = sum(abs(coor_T(1:2,:)-coor2));
            coor1_match = [];
            coor2_match = [];
            k = 1;
            for n = 1:size(coor1,2)
                if dist(n) < epsi
                    coor1_match(:,k) = coor1(:,n);
                    coor2_match(:,k) = coor2(:,n);
                    k = k+1;
                end
            end
            number = size(coor1_match,2);
        else number = 0
        end

        if number == 0
            coor3D = double.empty(0,2);
            coor2Dslice = double.empty(0,2);
        else
            
            % third outliers rejection

            % if they too close - no 
            N = size(coor1_match, 2);
            remove_elements = []; 

            for i = 1:N
                if min(max(abs(coor1_match(:,setdiff(1:N, i))-repmat(coor1_match(:,i), 1, N-1)))) <= 2
                    remove_elements = [remove_elements, i];
                    coor1_match(:,i) = [NaN,NaN];
                end
            end
            coor1_match(:, remove_elements) = [];
            coor2_match(:, remove_elements) = [];


            N = size(coor2_match, 2);
            remove_elements = []; 
            for i = 1:N
                if min(max(abs(coor2_match(:,setdiff(1:N, i))-repmat(coor2_match(:,i), 1, N-1)))) <= 2
                    remove_elements = [remove_elements, i];
                    coor2_match(:,i) = [NaN,NaN];
                end
            end
            coor1_match(:, remove_elements) = [];
            coor2_match(:, remove_elements) = [];


            coor3D =   coor2_match';
            coor2Dslice = coor1_match';
            number = size(coor3D,1);
        end
    end 
end


    
end
