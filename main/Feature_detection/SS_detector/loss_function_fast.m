function [R, al,t] = loss_function_fast(coor1, coor2, epsilon)
    R = [1, 0, 0;
         0, 1, 0;
         0, 0, 1];
    al = 0;

    x1 = coor1(1,:);
    y1 = coor1(2,:);
    x2 = coor2(1,:);
    y2 = coor2(2,:);

    L_star = Inf;
    number = length(x1);
    
    for j = 1:number
        for k = 1:number
            % Step 1: Compute break points
          
            % Compute coefficients
            as = x1 - x1(j);
            bs = y1 - y1(j);
            cs = x2 - x2(j);
            ds = x1 - x1(k);
            es = y1 - y1(k);
            fs = y2 - y2(k);
            
            indexColumn = 1:number;
            
            % Solve for zero crossings of X term
            det = as.*as + bs.*bs - cs.*cs;
            indices = and(det >= 0, as + cs ~= 0);
            % Have two possible solutions
            a1 = 2*atan((-bs(indices) + sqrt(det(indices)))./(as(indices) + cs(indices)));
            a2 = 2*atan((-bs(indices) - sqrt(det(indices)))./(as(indices) + cs(indices)));
            alphas = [a1', indexColumn(indices)'; a2', indexColumn(indices)'];
            
            % Solve for zero crossings of Y term
            det = ds.*ds + es.*es - fs.*fs;
            indices = and(det >= 0, es + fs);
            % Have two possible solutions
            a1 = 2*atan((ds(indices) - sqrt(det(indices)))./(es(indices) + fs(indices)));
            a2 = 2*atan((ds(indices) + sqrt(det(indices)))./(es(indices) + fs(indices)));
            alphas = [alphas; a1', indexColumn(indices)'; a2', indexColumn(indices)'];
            
            % Solve for saturation points (i.e. sum of absolutes = epsilon)
            alphas = [alphas; solveTruncationVec( as,  bs,  cs,  ds,  es,  fs, epsilon, indexColumn)];
            alphas = [alphas; solveTruncationVec(-as, -bs, -cs,  ds,  es,  fs, epsilon, indexColumn)];
            alphas = [alphas; solveTruncationVec( as,  bs,  cs, -ds, -es, -fs, epsilon, indexColumn)];
            alphas = [alphas; solveTruncationVec(-as, -bs, -cs, -ds, -es, -fs, epsilon, indexColumn)];
            
            
            % Step 2: Compute L1 norm at each break- and stationary point
            
            if ~isempty(alphas)
                % Find uniques and sort the break points
                alpha_sort = unique(alphas, 'rows');
                
                % Make sure to include the interval that wraps from 2 pi -> 0
                alpha_circ = [alpha_sort; alpha_sort(1, :) + [2*pi, 0]];
                
                % Find midpoints of all intervals bordering a break point
                % on the left and right. This function skips empty intervals
                [left_midpoints, right_midpoints] = findAdjacentMidpoints(alpha_circ(:, 1)');
                
                % At each break point, there is only one i for which
                % w1,w2,w3 changes. Keep track of the w1,w2,w3 right before
                % and after the break point, then compute the difference to
                % know the delta of the total w1,w2,w3
                indices = alpha_circ(2:end - 1, 2);
                 preWs = normCoefficientsVec( left_midpoints(2:end - 1), as(indices), bs(indices), cs(indices), ds(indices), es(indices), fs(indices), epsilon);
                postWs = normCoefficientsVec(right_midpoints(2:end - 1), as(indices), bs(indices), cs(indices), ds(indices), es(indices), fs(indices), epsilon);
                
                % Compute initial w1,w2,w3 inside first interval
                ws = sum(normCoefficientsVec(ones(1, number)*right_midpoints(1), as, bs, cs, ds, es, fs, epsilon));
                % Add delta at each breakpoint -> now know w1,w2,w3 in all intervals
                wFull = cumsum([ws; postWs - preWs]);
                
                % Compute norm at each breakpoint
                break_L_local = wFull(:, 1).*cos(alpha_sort(:, 1)) + wFull(:, 2).*sin(alpha_sort(:, 1)) + wFull(:, 3);
                
                % Find minimum and track it if smaller than L*
                [L_local, idx_local] = min(break_L_local);
                if L_local < L_star
                    % Sanity check against (slow) ground truth L1 norm
                    ground_truth_L1 = truncatedL1(coor1, coor2, j, k, alpha_sort(idx_local, 1), epsilon);
                    % If L_local is much smaller than the ground truth L1,
                    % we've almost certainly hit a numerical edge case and
                    % we're better off not using it.
                    if ground_truth_L1 - L_local < 1e-10
                        L_star = L_local;
                        al = alpha_sort(idx_local, 1);
                        R = [cos(al), -sin(al); sin(al), cos(al)];
                        t = [(coor2(1, j) - R(1, :)*coor1(:, j));
                             (coor2(2, k) - R(2, :)*coor1(:, k))];
                        fprintf('L*: %f angle: %f t: [%f %f] (break point)\n', L_star, al*180/pi, t);
                    end
                end
                
                % Now: Solve for stationary points
                alpha_left  = alpha_circ(1:end - 1, 1);
                alpha_right = alpha_circ(2:end, 1);
                
                % Only use intervals where denominator is nonzero
                non_zeroes = wFull(:, 1).^2 + wFull(:, 2).^2 ~= 0;
                wFull = wFull(non_zeroes, :);
                alpha_left = alpha_left(non_zeroes);
                alpha_right = alpha_right(non_zeroes);
                
                % Solve with eq (7)
                sin_alphas = wFull(:, 2)./sqrt(wFull(:, 1).^2 + wFull(:, 2).^2);
                alpha_locals = asin(sin_alphas);
                
                % Only keep alphas that are inside bounds
                if isempty(alpha_locals); continue; end% !!! aded
                 inside_bounds = and(alpha_locals >= alpha_left, alpha_locals <= alpha_right);
                alpha_locals = alpha_locals(inside_bounds);
                wFull = wFull(inside_bounds, :);
                
                if isempty(wFull); continue; end% why not ~isempty
    
                % Now, compute L_local at each stationary point
                stationary_L_local = wFull(:, 1).*cos(alpha_locals) + wFull(:, 2).*sin(alpha_locals) + wFull(:, 3);
                % Find minimum and track it if smaller than L*
                [L_local, idx_local] = min(stationary_L_local);
                if L_local < L_star
                    % Sanity check against (slow) ground truth L1 norm
                    ground_truth_L1 = truncatedL1(coor1, coor2, j, k, alpha_locals(idx_local), epsilon);
                    % If L_local is much smaller than the ground truth L1,
                    % we've almost certainly hit a numerical edge case and
                    % we're better off not using it.
                    if ground_truth_L1 - L_local < 1e-10
                        L_star = L_local;
                        al = alpha_locals(idx_local);
                        R = [cos(al), -sin(al); sin(al), cos(al)];
                        t = [(coor2(1, j) - R(1, :)*coor1(:, j));
                             (coor2(2, k) - R(2, :)*coor1(:, k))];

                        fprintf('L*: %f angle: %f t: [%f %f] (stationary point)\n', L_star, al*180/pi, t);
                    end
                end
            end
        end
        %fprintf('Finished for a point j=%d\n', j);
    end
    
    % Wrap negative angles around to be in [0, 2*pi]
    if al < 0; al = al + 2*pi; end
    
    R = [cos(al), -sin(al), 0;
         sin(al),  cos(al), 0;
               0,        0, 1];
    if L_star == Inf
        R = NaN;
        t = NaN;
    end
end

function [left_mid, right_mid] = findAdjacentMidpoints(alphas)
    nonempty_intervals = alphas(2:end) - alphas(1:end - 1) > 1e-12;
    filtered_midpoints_l = (alphas([false, nonempty_intervals]) + alphas([nonempty_intervals, false]))/2;
    filtered_midpoints_r = circshift(filtered_midpoints_l, -1, 2);
    M = length(filtered_midpoints_l);
    
    indices = mod(cumsum(nonempty_intervals) + M - 1, M) + 1;
    indices = [indices(end), indices];
   
    left_mid = filtered_midpoints_l(indices);
    right_mid = filtered_midpoints_r(indices);
end

function error = truncatedL1(points1, points2, j, k, alpha, epsilon)
    error = 0;
    
    R = [cos(alpha), -sin(alpha);
         sin(alpha),  cos(alpha)];
    t = [(points2(1, j) - R(1, :)*points1(:, j));
         (points2(2, k) - R(2, :)*points1(:, k))];
    
    for i = 1:size(points1, 2)
        error = error + min(epsilon, sum(abs(R*points1(:, i) + t - points2(:, i))));
    end
end

function alphas = solveTruncationVec(a, b, c, d, e, f, epsilon, indexColumn)
    A = a + c + f + e + epsilon;
    B = 2*(b - d);
    C = -a + c + f - e + epsilon;

    det = B.^2 - 4*A.*C;
    
    indices = and(det >= 0, A ~= 0);
    det = det(indices);
    A = A(indices);
    B = B(indices);
    indexColumn = indexColumn(indices);

    tanX1 = (-B - sqrt(det))./(2*A);
    tanX2 = (-B + sqrt(det))./(2*A);
    
    alpha1 = 2*atan(tanX1);
    alpha2 = 2*atan(tanX2);
        
    term11 = a(indices).*cos(alpha1) - b(indices).*sin(alpha1) - c(indices);
    term12 = d(indices).*sin(alpha1) + e(indices).*cos(alpha1) - f(indices);
    term21 = a(indices).*cos(alpha2) - b(indices).*sin(alpha2) - c(indices);
    term22 = d(indices).*sin(alpha2) + e(indices).*cos(alpha2) - f(indices);
    
    index1 = and(term11 >= 0, term12 >= 0);
    index2 = and(term21 >= 0, term22 >= 0);
    
    alphas = [alpha1(index1)', indexColumn(index1)'; alpha2(index2)', indexColumn(index2)'];
end

function ws = normCoefficientsVec(alpha_half, a, b, c, d, e, f, epsilon)
    term1 = a.*cos(alpha_half) - b.*sin(alpha_half) - c;
    term2 = d.*sin(alpha_half) + e.*cos(alpha_half) - f;
    
    ws = repmat(sign(term1'), 1, 3).*[a', -b', -c'] + repmat(sign(term2'), 1, 3).*[e', d', -f'];
    saturation_indices = abs(term1) + abs(term2) >= epsilon;
    ws(saturation_indices, 1) = 0;
    ws(saturation_indices, 2) = 0;
    ws(saturation_indices, 3) = epsilon;
end
