function norm_data = quantilenorm(data)
    % A fast, Octave-compatible implementation of Quantile Normalization
    
    % 1. Sort each column (sample) independently and save original indices
    [sorted_data, original_indices] = sort(data, 1);
    
    % 2. Calculate the mean expression for each rank (row) across all samples
    % (Using core sum and size functions to bypass missing 'mean.m' path in Conda)
    num_samples = size(sorted_data, 2);
    rank_means = sum(sorted_data, 2) / num_samples;
    
    % 3. Map the rank means back to their original positions in the matrix
    norm_data = zeros(size(data));
    for i = 1:num_samples
        norm_data(original_indices(:, i), i) = rank_means;
    end
end