% Particle filter with extended motion model to wrap particles around state
% space if a bound is exceeded

clc
clear all
close all

start = 50.5;
video = VideoReader("video/pacman_full.mp4") ;
params.Sigma_R = diag([50 50]);

% start = 5;
% video = VideoReader("video/pacman.mp4") ;
% params.Sigma_R = diag([200 200]);

startTime = start + 25;
stopTime = start + 30;
video.CurrentTime = startTime;

bottom_cut = video.Height / 15;

video_writer = VideoWriter("pf_inject_at_color.avi"); %create the video object
video_writer.FrameRate = video.FrameRate;
open(video_writer); %open the file for writing

% video.FrameRate
params.M = 1000 ; 
% params.pcm_colour = [255,231,55];
params.pcm_colour = [255,255,0];

% params.state_space_bound = [video.Width; video.Height-100]; %1920 1080 
% params.bounds = [100, 1000;450, 1450]; % height bounds; width bounds
params.bounds = [1, video.Height - bottom_cut;1, video.Width]; % height bounds; width bounds
% params.Sigma_R = diag([400 400]);
params.Sigma_R = diag([20 20]);

params.cutoff_dist = 25;

%% Variable Initialization %%
% Initialize Sample Set

S.X = [params.bounds(2,1) + (params.bounds(2,2) - params.bounds(2,1)) * rand(1, params.M);% colums
    params.bounds(1,1) + (params.bounds(1,2) - params.bounds(1,1)) * rand(1, params.M) % row
    ]; 

% Initialize Weights
S.W = 1/params.M * ones(1,params.M); 


weight_avgs = [];
pos_estimates = [];
pos_groundtruths = [];
pos_errs = [];
particle_stddevs = [];
time = [];


%% Main LOOP %%
while hasFrame(video) && video.CurrentTime < stopTime
    vidFrame = readFrame(video); %read video frame of pacmans, class: uint8

    histogram = color_histogram(vidFrame, params.pcm_colour);
    S_bar = pf_predict(S, params);
    [S_bar, weight_avg] = pf_weight(S_bar, params, histogram);

    S = pf_sys_resamp(S_bar);
%     S = multinomial_resample(S_bar, params);

%     subplot(2,2,1);
%     image(vidFrame,Parent=gca);
% 
%     subplot(2,2,2);
%     imshow(mat2gray(histogram),Parent=gca);

%     subplot(2,2,3);
    imshow(vidFrame,Parent=gca)
    pos_estimate = plot_particles(S);
    % FIXME: ugly hack to make pacman center position work
    % look into how to make this work with original bounds in the find
    % funtion
    bounds = params.bounds;
    params.bounds = [1, video.Height - bottom_cut; 1, video.Width]; % height bounds; width bounds
    pos_groundtruth = get_pacman_center(vidFrame, params);
%     plot_pacman_center(vidFrame, params);
    params.bounds = bounds; 

    F = getframe(gcf);
    writeVideo(video_writer,F); % Write the image to file.

    time = [time video.CurrentTime];
    weight_avgs = [weight_avgs weight_avg];
    pos_estimates = [pos_estimates pos_estimate];
    pos_groundtruths = [pos_groundtruths pos_groundtruth];
    pos_err = norm(pos_estimate - pos_groundtruth);
    pos_errs = [pos_errs pos_err];
    stddev = [std(S.X(1,:)); std(S.X(2,:))];
    particle_stddevs = [particle_stddevs stddev];


%     hold on
%     plot(params.bounds(2,1), params.bounds(1,1),'x','Color','yellow') %plot the particles avg 
%     drawnow
%     hold off
%     pause(0.1)
end

close(video_writer);

pf_plot_stats(time,weight_avgs, pos_estimates, pos_groundtruths, pos_errs, particle_stddevs);


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
