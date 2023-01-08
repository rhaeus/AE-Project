video = VideoReader("video/pacman.mp4") ;
video.CurrentTime = 28;
% video.FrameRate
params.M = 1000 ; 
% params.pcm_colour = [255,231,55];
params.pcm_colour = [255,255,0];

% params.state_space_bound = [video.Width; video.Height-100]; %1920 1080 
params.bounds = [100, 1000;450, 1450]; % height bounds; width bounds
params.Sigma_R = diag([400 400]);

params.cutoff_dist = 25;

%% Variable Initialization %%
% Initialize Sample Set

S.X = [params.bounds(2,1) + (params.bounds(2,2) - params.bounds(2,1)) * rand(1, params.M);% colums
    params.bounds(1,1) + (params.bounds(1,2) - params.bounds(1,1)) * rand(1, params.M) % row
    ]; 

% Initialize Weights
S.W = 1/params.M * ones(1,params.M); 

% imshow(vidFrame,Parent=gca)
% hold on
% plot(S.X(1,:),S.X(2,:),'.','Color','green') %plot the particles 
% drawnow
% hold off
% 
% pause(10)


%% Main LOOP %%
while hasFrame(video) 
    vidFrame = readFrame(video); %read video frame of pacmans, class: uint8

    histogram = calc_color_histogram(vidFrame, params.pcm_colour);
    S_bar = pf_predict(S, params);
    S_bar = pf_weight(S_bar, params, histogram) ;

    S = pf_sys_resamp(S_bar);
%     S = multinomial_resample(S_bar, params);

%     subplot(2,2,1);
    image(vidFrame,Parent=gca);
% 
%     subplot(2,2,2);
%     imshow(mat2gray(histogram),Parent=gca);

%     subplot(2,2,3);
    imshow(vidFrame,Parent=gca)
    plot_particles(S);
    % FIXME: ugly hack to make pacman center position work
    % look into how to make this work with original bounds in the find
    % funtion
    bounds = params.bounds;
    params.bounds = [1, video.Height - 100; 1, video.Width]; % height bounds; width bounds
    plot_pacman_center(vidFrame, params);
    params.bounds = bounds;


%     hold on
%     plot(params.bounds(2,1), params.bounds(1,1),'x','Color','yellow') %plot the particles avg 
%     drawnow
%     hold off
    pause(0.5)
end

function histogram = calc_color_histogram(frame, colour)
[H, W, D] = size(frame); %height, width, dimension of video frame matrix
RGB_matrix = double(reshape(frame,[H*W, D])); % create matrix with R G B values listed in separate columns
histogram = pdist2(RGB_matrix, colour, "euclidean");
histogram = 1./histogram ; 
% histogram = histogram / sum(sum(histogram)) ;

histogram = reshape(histogram, H, W);
% histogram_size = size(histogram);
% frame_proc = reshape(histogram, H, W);
% disp("histogram filter successful")

end

function S_bar = pf_predict(S, params)
N = size(S.X, 1) ;
%Diffusion assuming uncorrelated sigma R
S_bar.X = S.X + randn(N, params.M) .* repmat(sqrt(diag(params.Sigma_R)),1,params.M);

%wrap around bounds
S_bar.X(1,:) = wrapToInterval(S_bar.X(1,:), params.bounds(2,:));
S_bar.X(2,:) = wrapToInterval(S_bar.X(2,:), params.bounds(1,:));

S_bar.W = S.W;
% disp("pf_predict successful")
end

% histogram 1080x1920
% S.W 1xM
% S.X 2x1000
function S_bar = pf_weight(S, params, histogram)
    S_bar = S;
    x_coords = ceil(S_bar.X(1, :)) ; % 1x1000
    y_coords = ceil(S_bar.X(2, :)) ; % 1x1000


    for m=1:params.M
        x = x_coords(1,m);
        y = y_coords(1,m);

        if (x > params.bounds(2,2) || x <= params.bounds(2,1) || y > params.bounds(1,2) || y <= params.bounds(1,1))
            S_bar.W(m) = 0;
        else
            S_bar.W(m) = histogram(y, x);
        end
    end

    % normalize weights
    S_bar.W = S_bar.W ./ sum(S_bar.W) ;
end

% function S_bar = pf_weight(S_bar, params, histogram)
% row = ceil(S_bar.X(1, :)) ;
% col = ceil(S_bar.X(2, :)) ;
% max_width =  params.state_space_bound(2) ;
% max_height =  params.state_space_bound(1);
% if any(col > max_width) || any(row > max_height)
%     col(col>max_width) = max_width ; % columns/ width
%     row(row> max_height) = max_height ; %rows / height
% end
% for i = 1:1:params.M
% %     test1 = row(i);
% %     test2 = col(i);
%     try
%         S_bar.W(1, i) = histogram(row(i), col(i));
%     catch
%         fprintf("\nindexing error \n m= %d\n row = %d\n col= %d\n", i, row, col)
%         disp("error")
%     end
%     
%     %normalize
% end
% S_bar.W = S_bar.W ./sum(S_bar.W) ;
% % disp("pf_weight successful")
% end

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

function S = multinomial_resample(S_bar, params)

    cdf = cumsum(S_bar.W);
    S.X = zeros(size(S_bar.X));

    for m = 1 : params.M
        rm = rand;
        i = find(cdf >= rm,1,'first');
        S.X(:,m) = S_bar.X(:,i);
    end
    S.W = 1/params.M*ones(1,params.M);
end
