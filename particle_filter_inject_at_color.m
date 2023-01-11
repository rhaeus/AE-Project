% Particle filter that detects kidnap event through weight average and then
% injects particles at position with color similar to pacman color

clc
clear all
close all

start = 50.5;
video = VideoReader("video/pacman_full.mp4") ;
params.Sigma_R = diag([50 50]);

% start = 5;
% video = VideoReader("video/pacman.mp4") ;
% params.Sigma_R = diag([200 200]);

startTime = start + 20;
stopTime = start + 40;
video.CurrentTime = startTime;


video_writer = VideoWriter("pf_inject_at_color.avi"); %create the video object
video_writer.FrameRate = video.FrameRate;
open(video_writer); %open the file for writing

bottom_cut = video.Height / 15;
% bottom_cut = 0;

% video.FrameRate
params.M = 1000 ; 
% params.pcm_colour = [255,231,55];
params.pcm_colour = [255,255,0];

% params.state_space_bound = [video.Width; video.Height-100]; %1920 1080 
params.bounds = [1, video.Height - bottom_cut; 1, video.Width]; % height bounds; width bounds


params.cutoff_dist = 25;

params.random_particles = 10;


%% Variable Initialization %%
S = init(params);

warmup_it = 5;
it = 0;

weight_avgs = [];
pos_estimates = [];
pos_groundtruths = [];
pos_errs = [];
particle_stddevs = [];
time = [];

mins = [];
maxs = [];
i = 1;

%% Main LOOP %%
while hasFrame(video) && video.CurrentTime < stopTime
    it = it + 1;
    i = i + 1;

    vidFrame = readFrame(video); %read video frame of pacmans, class: uint8
%     video.CurrentTime
%     if mod(i, 2) == 0
%         continue; % lower frame rate
%     end
    

    histogram = color_histogram(vidFrame, params.pcm_colour);
    mins = [mins min(min(histogram))];
    maxs = [maxs max(max(histogram))];
    S_bar = pf_predict(S, params);
    [S_bar, weight_avg] = pf_weight(S_bar, params, histogram);


    if it > warmup_it && weight_avg <= 0.003
        display('we think we are lost, inject at pacman color position')

        p = find_pacman_positions(vidFrame, params);
        S = pf_sys_resamp(S_bar, params, p);
        it = 0;
    else
        S = pf_sys_resamp(S_bar, params, []);
        % S = multinomial_resample(S_bar, params);
    end



%     subplot(1,2,1);
%     image(vidFrame,Parent=gca);
% 
%     subplot(2,2,2);
%     imshow(mat2gray(histogram),Parent=gca);

%     subplot(2,2,3);
    imshow(vidFrame,Parent=gca)
    pos_estimate = plot_particles(S);
    pos_groundtruth = get_pacman_center(vidFrame, params);
%     plot_pacman_center(vidFrame, params);

    time = [time video.CurrentTime];
    weight_avgs = [weight_avgs weight_avg];
    pos_estimates = [pos_estimates pos_estimate];
    pos_groundtruths = [pos_groundtruths pos_groundtruth];
    pos_err = norm(pos_estimate - pos_groundtruth);
    pos_errs = [pos_errs pos_err];
    stddev = [std(S.X(1,:)); std(S.X(2,:))];
    particle_stddevs = [particle_stddevs stddev];

    F = getframe(gcf);
    writeVideo(video_writer,F); % Write the image to file.

%     subplot(1,2,2);
%     frame = vidFrame(params.bounds(1,1):params.bounds(1,2),params.bounds(2,1):params.bounds(2,2),:);
%     dist_image = color_dist(frame, params.pcm_colour);
%     bin_image = dist_image > params.cutoff_dist;
%     imshow(mat2gray(histogram),Parent=gca);
%     plot(avgs);

%     pause(0.5)
end

pf_plot_stats(time,weight_avgs, pos_estimates, pos_groundtruths, pos_errs, particle_stddevs);

fprintf('min is %0.4f\n', min(mins));
fprintf('max is %0.4f\n', max(maxs));

close(video_writer);

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


function S = pf_sys_resamp(S_bar, params, random_particles)
    cdf = cumsum(S_bar.W);
    M = size(S_bar.X,2);
    S.X = zeros(size(S_bar.X));
    r_0 = rand / M;
    
    for m = 1 : M
        i = find(cdf >= r_0,1,'first');

        if isempty(i) % this is weird, sometimes find fails
            i = m;
        end
    
        S.X(:,m) = S_bar.X(:,i);
        r_0 = r_0 + 1/M;
    end
    S.W = 1/M*ones(size(S_bar.W));
    % disp("systematic resampling successful")
    
    % inject random particles
    if ~isempty(random_particles)
        for i = 1 : size(random_particles, 2)
            index = randi([1, params.M]);
        %     index = i;
            S.X(:, index) = random_particles(:,i);
        end
    end

end

% function S = multinomial_resample(S_bar, params)
% 
%     cdf = cumsum(S_bar.W);
%     S.X = zeros(size(S_bar.X));
% 
%     for m = 1 : params.M
%         rm = rand;
%         i = find(cdf >= rm,1,'first');
%         S.X(:,m) = S_bar.X(:,i);
%     end
%     S.W = 1/params.M*ones(1,params.M);
% end
