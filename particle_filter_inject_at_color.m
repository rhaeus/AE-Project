video = VideoReader("video/pacman.mp4") ;
video.CurrentTime = 27;
% video.FrameRate
params.M = 1000 ; 
% params.pcm_colour = [255,231,55];
params.pcm_colour = [255,255,0];

% params.state_space_bound = [video.Width; video.Height-100]; %1920 1080 
params.bounds = [1, video.Height - 100; 1, video.Width]; % height bounds; width bounds
params.Sigma_R = diag([400 400]);

params.cutoff_dist = 25;

params.random_particles = 10;

%% Variable Initialization %%
S = init(params);

warmup_it = 5;
it = 0;

avgs = [];


%% Main LOOP %%
while hasFrame(video) 
    it = it + 1;
    vidFrame = readFrame(video); %read video frame of pacmans, class: uint8
%     video.CurrentTime

    histogram = calc_color_histogram(vidFrame, params.pcm_colour);
    S_bar = pf_predict(S, params);
    [S_bar, weight_avg] = pf_weight(S_bar, params, histogram);

%     weight_avg
%     avgs = [avgs weight_avg];

    frame = vidFrame(params.bounds(1,1):params.bounds(1,2), params.bounds(2,1):params.bounds(2,2),:);
    dist_image = color_dist(frame, params.pcm_colour); %980x1920
    bin_image = dist_image < params.cutoff_dist;
    pos = find(bin_image, params.random_particles, 'first'); % get positions of particles at target color
    n = size(pos,1);
%     pos = ind2sub([params.bounds(1,2), params.bounds(2,2)], pos);
%     pos = ind2sub(size(bin_image), pos);

    p = zeros([2 n]);
    for i = 1 : n
        [c,r] = ind2sub(size(bin_image), pos(i));
        p(1,i) = r;
        p(2,i) = c;
    end
%     p

%     if it > warmup_it && weight_avg < 0.003
%         display('we think we are lost, re-init')
%         S = init(params);
%         it = 0;
%     else
        S = pf_sys_resamp(S_bar, params, p);
        % S = multinomial_resample(S_bar, params);
%     end



%     subplot(1,2,1);
    image(vidFrame,Parent=gca);
% 
%     subplot(2,2,2);
%     imshow(mat2gray(histogram),Parent=gca);

%     subplot(2,2,3);
    imshow(vidFrame,Parent=gca)
    hold on
    plot(S.X(1,:),S.X(2,:),'.','Color','green') %plot the particles 
    drawnow
    hold off

    hold on
    plot(mean(S.X(1,:)),mean(S.X(2,:)),'x','Color','red') %plot the particles avg 
    drawnow
    hold off

%     subplot(1,2,2);
%     plot(avgs);

%     hold on
%     plot(params.bounds(2,1), params.bounds(1,1),'x','Color','yellow') %plot the particles avg 
%     drawnow
%     hold off
    pause(0.2)
end

function dist_mat = color_dist(frame, colour)
    [H, W, D] = size(frame); %height, width, dimension of video frame matrix
    RGB_matrix = double(reshape(frame,[H*W, D])); % create matrix with R G B values listed in separate columns
    dist_mat = pdist2(RGB_matrix, colour, "euclidean");
    dist_mat = reshape(dist_mat, H, W);
end

function pos = getRandomPos(params, count)
    pos = [params.bounds(2,1) + (params.bounds(2,2) - params.bounds(2,1)) * rand(1, count);% colums
        params.bounds(1,1) + (params.bounds(1,2) - params.bounds(1,1)) * rand(1, count) % row
        ];
end

function S = init(params)
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
% S_bar.X(1,:) = wrapToInterval(S_bar.X(1,:), params.bounds(2,:));
% S_bar.X(2,:) = wrapToInterval(S_bar.X(2,:), params.bounds(1,:));

S_bar.W = S.W;
% disp("pf_predict successful")
end

% histogram 1080x1920
% S.W 1xM
% S.X 2x1000
function [S_bar, w_avg] = pf_weight(S, params, histogram)
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

    w_avg = mean(S_bar.W);

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


function S = pf_sys_resamp(S_bar, params, random_particles)
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
    
    % inject random particles
    for i = 1 : size(random_particles, 2)
        index = randi([1, params.M]);
    %     index = i;
        S.X(:, index) = random_particles(:,i);
    end

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

    % inject random particles
    for i = 1 : params.random_particles
        index = randi([1, params.M]);
        pos = getRandomPos(params);
        S.X(:, index) = pos;
    end
end
