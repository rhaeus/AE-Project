video = VideoReader("pacman.mp4") ;
params.M = 1000 ; 
params.pcm_colour = [255,231,55];
params.state_space_bound = [video.Height; video.Width];
params.Sigma_R = diag([2 2]);

%% Variable Initialization %%
% Initialize Sample Set
S.X = [rand(1, params.M)*params.state_space_bound(1); 
     rand(1, params.M)*params.state_space_bound(2)];
% Initialize Weights
S.W = 1/params.M * ones(1,params.M); 
%% Main LOOP %%
while hasFrame(video) 
    vidFrame = readFrame(video); %read video frame of pacmans, class: uint8
    histogram = color_histogram(vidFrame, params.pcm_colour);
    S_bar = pf_predict(S, params);
    S_bar = pf_weight(S_bar, params, histogram) ;
    S = pf_sys_resamp(S_bar);

    imshow(vidFrame)
    hold on
    plot(S.X(1,:),S.X(2,:),'.','Color','green') %plot the particles 
    drawnow
    hold off
end

function histogram = color_histogram(frame, colour)
[H, W, D] = size(frame); %height, width, dimension of video frame matrix
RGB_matrix = double(reshape(frame,[H*W, D])); % create matrix with R G B values listed in separate columns
histogram = pdist2(RGB_matrix, colour, "euclidean");
histogram = 1./histogram ; 

histogram = reshape(histogram, H, W);
% histogram_size = size(histogram);
% frame_proc = reshape(histogram, H, W);
% disp("histogram filter successful")

end

function S_bar = pf_predict(S, params)
N = size(S.X, 1) ;
%Diffusion assuming uncorrelated sigma R
S_bar.X = S.X + randn(N, params.M) .* repmat(sqrt(diag(params.Sigma_R)),1,params.M);
% disp("pf_predict successful")
end

function S_bar = pf_weight(S_bar, params, histogram)
row = ceil(S_bar.X(1, :)) ;
col = ceil(S_bar.X(2, :)) ;
max_width =  params.state_space_bound(2) ;
max_height =  params.state_space_bound(1);
if any(col > max_width) || any(row > max_height)
    col(col>max_width) = max_width ; % columns/ width
    row(row> max_height) = max_height ; %rows / height
end
for i = 1:1:params.M
%     test1 = row(i);
%     test2 = col(i);
    try
        S_bar.W(1, i) = histogram(row(i), col(i));
    catch
        fprintf("\nindexing error \n m= %d\n row = %d\n col= %d\n", i, row, col)
        disp("error")
    end
    
    %normalize
end
S_bar.W = S_bar.W ./sum(S_bar.W) ;
% disp("pf_weight successful")
end

function S = pf_sys_resamp(S_bar)
cdf = cumsum(S_bar.W);
M = size(S_bar.X,2);
S.X = zeros(size(S_bar.X));
r_0 = rand / M;
for m = 1 : M
    i = find(cdf >= r_0,1,'first');
    S.X(:,m) = S_bar.X(:,i);
    r_0 = r_0 + 1/M;
end
S.W = 1/M*ones(size(S_bar.W));
% disp("systematic resampling successful")

end
