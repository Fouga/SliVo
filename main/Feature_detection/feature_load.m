function [descriptors,locs]=feature_load(filename,method)


% Open tmp.key and check its header
g = fopen(filename, 'r');
if g == -1
    error('Could not open file Feature.key or .txt.');
end
[header, count] = fscanf(g, '%d %d', [1 2]);
if count ~= 2
    error('Invalid keypoint file beginning.');
end

if (strcmp(method,'SURF'))
    num = header(1);
    len = header(2);


    % Creates the two output matrices (use known size for efficiency)
    locs = double(zeros(num, 4));
    descriptors = double(zeros(num, len));

    % Parse tmp.key
    for i = 1:num
        [vector, count] = fscanf(g, '%f %f %f %f', [1 4]); %row col scale ori
        if count ~= 4
            error('Invalid keypoint file format');
        end
        locs(i, :) = vector(1, :);

        [descrip, count] = fscanf(g, '%f', [1 len]);
        if (count ~= len)
            error('Invalid keypoint file value!');
        end
        % Normalize each input vector to unit length
        descrip = descrip / sqrt(sum(descrip.^2));
        descriptors(i, :) = descrip(1, :);
    end


elseif (strcmp(method,'SL1'))
    num = header(1);
    len = header(2)-3;

    % Creates the two output matrices (use known size for efficiency)
    locs = double(zeros(num, 2));
    orient = double(zeros(num,1));
    descriptors = double(zeros(num, len));

    % Parse tmp.key
    for i = 1:num
        [vector, count] = fscanf(g, '%f %f %f', [1 3]); %row col scale ori
        if count ~= 3
            error('Invalid keypoint file format');
        end
        locs(i, :) = vector(1, 1:2);
        orient(i) = vector(1,3);
        [descrip, count] = fscanf(g, '%f', [1 len]);
        if (count ~= len)
            error('Invalid keypoint file value!');
        end
        % Normalize each input vector to unit length
        descrip = descrip / sqrt(sum(descrip.^2));
        descriptors(i, :) = descrip(1, :);
    end
end

fclose(g);
