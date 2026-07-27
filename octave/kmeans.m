function indProbes = kmeans(data, k, varargin)
    % A fast, dependency-free k-means implementation for GNU Octave
    % Bypasses the missing statistics package in Conda environments
    
    maxiter = 1000;
    [n_probes, n_samples] = size(data);
    if n_probes == 0
        error('kmeans: data has 0 rows. Gene intersection was empty — check Python-level guard.');
    end
    if n_probes < k
        error('kmeans: n_probes=%d < k=%d. Reduce k or increase gene count.', n_probes, k);
    end

    % 1. Randomly initialize k centroids from the data
    % (Using rand and sort to avoid missing randperm.m)
    [~, r_sort] = sort(rand(n_probes, 1));
    rand_idx = r_sort(1:k);
    centroids = data(rand_idx, :);
    
    indProbes = zeros(n_probes, 1);
    
    for iter = 1:maxiter
        % 2. Calculate squared Euclidean distance to each centroid
        distances = zeros(n_probes, k);
        for c = 1:k
            diff = data - centroids(c, :);
            distances(:, c) = sum(diff.^2, 2);
        end
        
        % 3. Assign each probe to the closest centroid
        prev_indProbes = indProbes;
        [~, indProbes] = min(distances, [], 2);
        
        % 4. Check for convergence using core math operators (bypasses isequal)
        if sum(indProbes ~= prev_indProbes) == 0
            break;
        end
        
        % 5. Update centroid positions
        for c = 1:k
            idx = (indProbes == c);
            count = sum(idx);
            if count > 0
                % Using sum/count instead of 'mean'
                centroids(c, :) = sum(data(idx, :), 1) / count;
            else
                % Handle empty cluster by picking a new random point
                % (Using rand and ceil to avoid missing randi.m)
                centroids(c, :) = data(ceil(rand() * n_probes), :);
            end
        end
    end
end