% calculate the color histogram
% for each pixel calculate euclidean distance to target color
% return reciprocal to make smaller distance a larger weight

function histogram = color_histogram(frame, colour)
    [H, W, D] = size(frame); %height, width, dimension of video frame matrix
    RGB_matrix = double(reshape(frame,[H*W, D])); % create matrix with R G B values listed in separate columns
    histogram = pdist2(RGB_matrix, colour, "euclidean");
    histogram = histogram + 1;
    
    histogram = 1./histogram ; 
    
    histogram = reshape(histogram, H, W);
end