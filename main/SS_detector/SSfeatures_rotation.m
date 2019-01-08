function [coords]= SSfeatures_rotation(image, rotations)%, method)
    % params for calculating the degree with higher angle interval
    % coords = [row, colum]
    parms = struct;
    parms.patch_size         = 5;
    parms.desc_rad           = 40;%40;
    parms.nrad               = 6;
    parms.nang               = 8;
    parms.var_noise          = 300000;
    parms.saliency_thresh    = 1;
    parms.homogeneity_thresh = 0.7;
    parms.snn_thresh         = 1;
    
    % for smaller angle interval
    % for final descriptor
    parms2 = struct;
    parms2.patch_size         = 5;
    parms2.desc_rad           = 40;%40;    
    parms2.nrad               = 4;   
    parms2.nang               = 10;
    parms2.var_noise          = 300000;
    parms2.saliency_thresh    = 1;                        
    parms2.homogeneity_thresh = 0.7;
    parms2.snn_thresh         = 1;
    % sample rate of the image for descriptor 
    sample = 5;

    
    [n, m]=size(image); % n = row, m-column
    doubleImage = double(image);

    margin = parms.desc_rad + (parms.patch_size-1)/2;
    x_range = roundInterval([margin + 1, m - margin], 3);
    y_range = roundInterval([margin + 1, n - margin], 3);

    coords = [];
    
    evalX = [1; 2; 3];
    fitMat = inv([evalX .* evalX, evalX, [1; 1; 1]]);
    fx = 1:0.01:3;
    fxSq = fx.^2;


        [respFull, respCoords, ~, ~, ~] = mexCalcSsdescs(doubleImage, parms,[x_range(1), y_range(1), x_range(2), y_range(2), sample]);
    if rotations == 1
    % 2 in the end means how much to skip
    for idx = 1:size(respCoords, 2)
        j = respCoords(1, idx);
        i = respCoords(2, idx);
        resp = respFull(:, idx);

        if image(i,j) > 2 % calculate only for non zero pixels
            resp = resp/norm(resp);
            histogram1 = sum(reshape(resp, [parms.nrad parms.nang]), 1);
            circ = [histogram1(end), histogram1, histogram1(1)];
            [~, index] = find(histogram1 > 0.8*max(histogram1) & ...
                              circ(2:parms.nang+1) > circ(1:parms.nang) & ...
                              circ(2:parms.nang+1) > circ(3:parms.nang + 2));
            values = [circ(index); circ(index + 1); circ(index + 2)];
            fs = fitMat*values;
            funs = fs(1, :)'*fxSq + fs(2,:)'*fx;
            [~, degs] = max(funs, [], 2);
            rots = (index' + (degs*0.01 - 1))*(360/parms.nang);
            coords = [coords; rots, repmat([i, j], length(index), 1)];
        end
    end

    coords = sortrows(coords);

    curAngle = 0;
    rothistology = doubleImage;
    des = [];
    for idx = 1:size(coords, 1)
        rot = coords(idx, 1);
        i = coords(idx, 2);
        j = coords(idx, 3);

        if curAngle ~= rot
            rothistology = double(imrotate(image,-rot,'crop'));
            curAngle = rot;
        end

        [resp, ~, ~, ~, ~] = mexCalcSsdescs(rothistology, parms2, [j, i]);
        if ~isempty(resp) && ~max(isnan(resp))
            des = [des, [i;j;rot;resp]];
        end
    end
    if ~isempty(des); 
        des = sortrows(des', [2, 1])';

        % save features as txt file
        header=[size(des,2), size(des,1)];
        dlmwrite('key.txt',header,'delimiter',' ');
        dlmwrite('key.txt',des','-append','delimiter',' ');
    end
    else 

        des = [];
        for idx = 1:size(respCoords, 2)
%             rot = coords2(idx, 1);
        j = respCoords(1, idx);
        i = respCoords(2, idx);
        resp = respFull(:, idx);
%             [resp, ~, ~, ~, ~] = mexCalcSsdescs(doubleImage, parms2, [j, i]);

            if ~isempty(resp) && ~max(isnan(resp))
                des = [des, [i;j;0;resp]];
            end
        
        end
        
        if ~isempty(des); 
        des = sortrows(des', [2, 1])';

        % save features as txt file
        header=[size(des,2), size(des,1)];
        dlmwrite('key.txt',header,'delimiter',' ');
        dlmwrite('key.txt',des','-append','delimiter',' ');
        end
    end
end

function rounded = roundInterval(interval, slot)
    rounded = [ceil(interval(1)/slot)*slot + 1, floor((interval(2) - 2)/slot)*slot + 1];
end