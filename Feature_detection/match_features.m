function [number,coor2Dslice,coor3D]=match_features(des2D,locs2D,des3D,locs3D)

k = 1;


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
    
end
