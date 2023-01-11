% Particle filter

clc
clear all
close all

start = 50.5;
video = VideoReader("video/pacman_full.mp4") ;
params.Sigma_R = diag([20 20]);

% start = 5;
% video = VideoReader("video/pacman.mp4") ;
% params.Sigma_R = diag([200 200]);

startTime = start;
stopTime = start + 10;
video.CurrentTime = startTime;

bottom_cut = video.Height / 15;

% video.FrameRate
params.M = 1000 ; 
% params.pcm_colour = [255,231,55];
params.pcm_colour = [255,255,0];

params.state_space_bound = [video.Width; video.Height-bottom_cut]; %1920 1080 
params.bounds = [1, video.Height - bottom_cut; 1, video.Width]; % height bounds; width bounds


params.cutoff_dist = 25;

%% Variable Initialization %%
% Initialize Sample Set

S.X = [rand(1, params.M)*params.state_space_bound(1); % x coord of each particle in the set
    rand(1, params.M)*params.state_space_bound(2)]; % y coord of each particle in the set

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
    image(vidFrame,Parent=gca);
% 
%     subplot(2,2,2);
%     imshow(mat2gray(histogram),Parent=gca);

%     subplot(2,2,3);
    imshow(vidFrame,Parent=gca)
    pos_estimate = plot_particles(S);
    pos_groundtruth = plot_pacman_center(vidFrame, params);

    time = [time video.CurrentTime];
    weight_avgs = [weight_avgs weight_avg];
    pos_estimates = [pos_estimates pos_estimate];
    pos_groundtruths = [pos_groundtruths pos_groundtruth];
    pos_err = norm(pos_estimate - pos_groundtruth);
    pos_errs = [pos_errs pos_err];
    stddev = [std(S.X(1,:)); std(S.X(2,:))];
    particle_stddevs = [particle_stddevs stddev];

    pause(0.1)
end

pf_plot_stats(time,weight_avgs, pos_estimates, pos_groundtruths, pos_errs, particle_stddevs);

function S_bar = pf_predict(S, params)
N = size(S.X, 1) ;
%Diffusion assuming uncorrelated sigma R
S_bar.X = S.X + randn(N, params.M) .* repmat(sqrt(diag(params.Sigma_R)),1,params.M);
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
