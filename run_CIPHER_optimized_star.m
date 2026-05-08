%% run_CIPHER_optimized.m
% Single-file runner + CIPHER solvers (analytic-t OR fast full-grid over n1,n2,t)
% Optimized Version (if looking for deterministic version, check other script)
% Units:
%   lambda_um in micrometers
%   t in micrometers
%   n1 dimensionless
%   n2 in um^2 (Cauchy: n(lambda)=n1 + n2/lambda^2)
%   Authors: Ana Doblas, Carlos Trujillo

clear; close all; clc;

%% -------------------- INPUT DATA --------------------
% Define phi_map (radians). For example:
% load('star_input_data.mat','unwrap_star_background_zero');
% phi_map = unwrap_star_background_zero;
% If you already have it in workspace, just uncomment:

%phi_map = unwrap_star_background_zero(200:800, 200:800); %(For fast tests)
%phi_map = unwrap_star_background_zero(100:2000, 900:920); %(For fast tests)
phi_map = unwrap_star_background_zero;

lambda = 0.532;              % Wavelength in microns
RI = 1.5148255; %https://www.benchmarktech.com/quantitativephasemicroscop/
d_true = phi_map.*(lambda / (2*pi*(RI-1)));
figure; 
imagesc(d_true); axis image; colormap(jet); colorbar; title('True Thickness Map [\mum]');


%% --------------SYSTEM PHYSICS --------------------
lambda_um = 0.532;     % illumination wavelength [um]
nm        = 1.000;     % surrounding medium RI (air)

%% -------------------- SEARCH BOUNDS (SEEDS) --------------------
% Your star step height range:
seeds = struct();
seeds.t_min  = 0.0;    % [um]
seeds.t_max  = 0.8;    % [um]

% Broad bounds (adjust as needed)
seeds.n1_min = 1.30;
seeds.n1_max = 1.80;

% n2 is in um^2. Keep it tight if you see boundary-hugging.
seeds.n2_min = 0.001;
seeds.n2_max = 0.005;

%% -------------------- OPTIONS --------------------
opts = struct();

% Choose mode:
%   'analytic_t' : solves t* analytically, clamps to bounds (UNSTABLE)
%   'full_grid'  : Fast version of deterministic CIPHER (Choose this one!)
opts.mode = 'full_grid';

% Multi-level coarse-to-fine refinement (For the analytic_t option, ignore if running full_grid)
opts.levels        = 3;
opts.window_factor = 2;
opts.shrink        = 0.25;

% Initial grid steps (coarse). Full-grid cost grows fast: start coarse, then refine.
opts.n1_step0 = 0.002;      % try 0.001–0.005 (haven't tested others)
opts.n2_step0 = 0.0002;     % try 1e-4–5e-4 (um^2)
opts.t_step0  = 0.005;      % try 0.002–0.01 (um)

% Speed options
opts.use_gpu    = true;     % requires Parallel Computing Toolbox (gpuArray)
opts.use_parfor = true;     % parfor (often disable if using GPU)
if opts.use_gpu
    % Using GPU inside parfor is usually not what you want (it collapses)
    opts.use_parfor = false;
end

% Robustness: ignore pixels with tiny phase (background drives solutions to bounds. Important for analytic_t)
opts.phi_thresh = 0;     % radians (tune; try 0.05–0.3). If "0", then no threshold applied.
opts.den_thresh = 1e-9;     % avoid division by near-zero denom (mostly irrelevant here)

%% -------------------- RUN --------------------
tic
[n1_map, n2_map, t_map, res2_map] = cipher_solve_map(phi_map, lambda_um, nm, seeds, opts);
toc

n_map = n1_map + n2_map./(lambda_um.^2);    % effective RI at lambda

%% -------------------- DISPLAY --------------------
figure('Name','CIPHER outputs');
subplot(2,3,1); imagesc(n_map); axis image; colormap(jet); colorbar; clim([1.3 1.8]); title('n(\lambda) map');

%Once 't_map' is computed, then choose how to filter (smooth):

% Edge-preserving smoothing (keeps step edges better than Gaussian)
%degree = 15;  % smoothing strength (try 5–30)
%t_smooth = imbilatfilt(t_map, degree);

% Median filter (great for granular artifacts) I like this better
w = 7;  % window size: 3,5,7...
t_smooth = medfilt2(t_map, [w w], 'symmetric');

% Smooth thickness map with Gaussian filter
%sigma = 1.0;  % try 0.5–2.0
%t_smooth = imgaussfilt(t_map, sigma);

% Adaptive Wiener filter
%w = 5;
%t_smooth = wiener2(t_map, [w w]);

subplot(2,3,2); imagesc(t_map); axis image; colormap(jet); colorbar; clim([0 1]); title('t map [\mum]');
subplot(2,3,3); imagesc(d_true); axis image; colormap(jet); colorbar; clim([0 1]); title('t map (true) [\mum]');
subplot(2,3,4); imagesc(abs(t_smooth-d_true)); axis image; colormap(jet); colorbar; clim([0 1]); title('abs(t_smooth-d_true) map');
subplot(2,3,5); imagesc(t_smooth); axis image; colormap(jet); colorbar; clim([0 1]); title('t_map smoothed [\mum]');
subplot(2,3,6); imagesc(abs(t_map-d_true)); axis image; colormap(jet); colorbar; clim([0 1]); title('abs(t_map-d_true) map [\mum]');


figure('Name','Residual (phase-fit) map');
imagesc(res2_map); axis image; colormap(jet); colorbar;
title('min J = (phi - phi_model)^2');

figure('Name','Error (n_map-RI) map');
imagesc(abs(n_map-RI)); axis image; colormap(jet); colorbar;
title('Error (n_map-RI) map');


% compute the Otsu threshold of t_smooth (returns normalized threshold in [0,1])
% convert to uint8 image for graythresh compatibility if necessary
if ~isfloat(t_smooth) || max(t_smooth(:)) > 1 || min(t_smooth(:)) < 0
    ts_norm = (t_smooth - min(t_smooth(:)));
    if max(ts_norm(:)) > 0
        ts_norm = ts_norm / max(ts_norm(:));
    end
else
    ts_norm = t_smooth;
end
level = graythresh(ts_norm);   % normalized threshold [0,1]
otsu_threshold = level * (max(t_smooth(:)) - min(t_smooth(:))) + min(t_smooth(:));  % back to original scale

max_t = max(max(t_smooth));
min_t = min(min(t_smooth));
real_thres = (max_t-min_t)*otsu_threshold+min_t;

% Create mask from t_smooth > real_thres (otzu threshold)
mask = t_smooth > real_thres;

% Display the mask as an overlay and standalone image
figure('Name','Mask');
subplot(1,2,1);
imagesc(ts_norm); axis image; colormap(gray); colorbar; title('Normalized t\_smooth (for thresholding)');
hold on;
% Overlay mask contour in red
contour(mask, [0.5 0.5], 'r', 'LineWidth', 1.5);
hold off;

subplot(1,2,2);
imagesc(mask); axis image; colormap(gray); colorbar; title('Binary mask (t\_smooth > Otsu threshold)');

% Using this mask, extract corresponding values from n_map into n_masked (as a vector)
n_masked = n_map(mask);

% apply median-filter to n_masked
w = 7;
n_smooth = medfilt2(n_masked, [w w], 'symmetric');

%% -------------------- OPTIONAL: LINE PROFILES --------------------

% x-coordinate(s) for the vertical profile lines
x_coords = [1, 11, 20];

% Loop through each x-coordinate to plot the profile lines
for i = 1:length(x_coords)
    x = x_coords(i);
    
    % Create a new figure for each profile line
    figure; 
    hold on; % Hold on to plot multiple lines
    
    % Extract the vertical profiles from t_map and d_true_resized
    %t_profile = t_map(:, x);
    t_profile = t_smooth(:, x);
    d_true_profile = d_true(:, x);
    
    % Plot the profiles
    plot(t_profile, 'DisplayName', ['t\_map Profile at x = ' num2str(x)], 'LineWidth', 1.5);
    plot(d_true_profile, 'DisplayName', ['d\_true Profile at x = ' num2str(x)], 'LineWidth', 1.5);
    
    % Customize the plot
    xlabel('Position (pixels)');
    ylabel('Thickness (um)');
    title(['Vertical Profile Lines Comparison at x = ' num2str(x)]);
    legend show; % Show legend for the profiles
    grid on; % Add grid for better visibility
    hold off; % Release the hold on the current plot
end

% --- 90° Circle perimeter profile (TOP -> RIGHT), actual pixel values ---
R  = 1500;      % radius [px]
r0 = 1752;      % center row
c0 = 442;       % center col

[ny, nx] = size(t_map);

% Arc length = (pi/2)*R  -> ~1 sample per pixel
nSamples = round((pi/2)*R);

% Start at top point: theta = -pi/2
% Go 90 degrees to the right: theta -> 0
theta = linspace(-pi/2, 0, nSamples);

% Continuous circle coordinates
rr = r0 + R*sin(theta);
cc = c0 + R*cos(theta);

% Nearest pixel coordinates (actual pixel values)
ri = round(rr);
ci = round(cc);

% Keep only points inside the image
inside = (ri >= 1 & ri <= ny & ci >= 1 & ci <= nx);
ri = ri(inside);
ci = ci(inside);
theta = theta(inside);

% Remove duplicates caused by rounding (keep order)
pts = [ri(:) ci(:)];
[ptsUnique, ia] = unique(pts, 'rows', 'stable');
ri = ptsUnique(:,1);
ci = ptsUnique(:,2);
theta = theta(ia);

% Extract pixel values
idx = sub2ind([ny nx], ri, ci);
profile90 = t_smooth(idx);
profile90_d_true = d_true(idx); % Extract d_true values for the same points

% Plot vs angle from the top (0°..90°)
angleFromTop_deg = (theta + pi/2) * 180/pi;

figure;
subplot(1, 2, 1); % Create a subplot for the image
imagesc(t_smooth); axis image; colormap(jet); colorbar; title('t\_smooth Image');
hold on; % Hold on to plot the profile line
plot(ci, ri, 'r', 'LineWidth', 1.5); % Plot the profile line on the image
hold off; % Release hold after plotting

subplot(1, 2, 2); % Create a subplot for the profile
plot(angleFromTop_deg, profile90, 'LineWidth', 1.2, 'DisplayName', 't\_smooth');
hold on; % Hold on to plot both profiles
plot(angleFromTop_deg, profile90_d_true, 'LineWidth', 1.2, 'DisplayName', 'd\_true');
hold off; % Release hold after plotting

xlabel('Angle from top (deg)');
ylabel('Profile value on circle');
grid on;
title('Perimeter Profile (Top → Right, 90°)');
legend show; % Show legend for both profiles


%% Statistical Analysis 

mean_t_smooth = mean(t_smooth(:));
std_t_smooth = std(t_smooth(:));

mean_d_true = mean(d_true(:));
std_d_true = std(d_true(:));

mean_n_map = mean(n_map(:));
std_n_map = std(n_map(:));

n_map_filtered = n_map(abs(n_map) > 1.4); % Filter out values close to zero
mean_n_map_filtered = mean(n_map_filtered(:)); % Compute mean
std_n_map_filtered = std(n_map_filtered(:)); % Compute standard deviation
 
mean_n_smooth = mean(n_smooth(:));
std_n_smooth = std(n_smooth(:));

mean_n_masked = mean(n_masked(:));
std_n_masked = std(n_masked(:));

% Compute SSIM and correlation compared to d_true
ssim_value = ssim(t_smooth, d_true);
correlation_value = corr2(t_smooth, d_true);

disp(['SSIM between t_smooth and d_true: ', num2str(ssim_value)]);
disp(['Correlation between t_smooth and d_true: ', num2str(correlation_value)]);

disp(['Mean of t_smooth: ', num2str(mean_t_smooth)]);
disp(['Standard Deviation of t_smooth: ', num2str(std_t_smooth)]);
disp(['Mean of d_true: ', num2str(mean_d_true)]);
disp(['Standard Deviation of d_true: ', num2str(std_d_true)]);
disp(['Mean of n_map: ', num2str(mean_n_map)]);
disp(['Standard Deviation of n_map: ', num2str(std_n_map)]);
disp(['Mean of n_smooth: ', num2str(mean_n_smooth)]);
disp(['Standard Deviation of n_smooth: ', num2str(std_n_smooth)]);
disp(['Mean of n_masked: ', num2str(mean_n_masked)]);
disp(['Standard Deviation of n_masked: ', num2str(std_n_masked)]);
disp(['Mean of n_map filtered: ', num2str(mean_n_map_filtered)]);
disp(['Standard Deviation of n_map filtered: ', num2str(std_n_map_filtered)]);

% Analysis for t_map
mean_t_map = mean(t_map(:));
std_t_map = std(t_map(:));

ssim_value_t_map = ssim(t_map, d_true);
correlation_value_t_map = corr2(t_map, d_true);

disp(['Mean of t_map: ', num2str(mean_t_map)]);
disp(['Standard Deviation of t_map: ', num2str(std_t_map)]);
disp(['SSIM between t_map and d_true: ', num2str(ssim_value_t_map)]);
disp(['Correlation between t_map and d_true: ', num2str(correlation_value_t_map)]);

error_smooth = t_smooth - d_true;
figure; 
imagesc(error_smooth); 
axis image; 
colormap(jet); 
colorbar; 
title('Difference between t\_smooth and d\_true');

mean_error_smooth = mean(error_smooth(:));
std_error_smooth = std(error_smooth(:));

disp(['Mean of error (t_smooth): ', num2str(mean_error_smooth)]);
disp(['Standard Deviation of error (t_smooth): ', num2str(std_error_smooth)]);

error_map = t_map - d_true;
figure; 
imagesc(error_map); 
axis image; 
colormap(jet); 
colorbar; 
title('Difference between t\_map and d\_true');

mean_error_map = mean(error_map(:));
std_error_map = std(error_map(:));

disp(['Mean of error (t_map): ', num2str(mean_error_map)]);
disp(['Standard Deviation of error (t_map): ', num2str(std_error_map)]);

%% ---- Forward model: couple n_map and t_map back into phase ----

k = 2*pi/lambda_um;                     % [rad/um]
phi_pred = k .* (n_map - nm) .* t_map;  % predicted phase [rad]

% Residual (difference)
phi_residual = phi_map - phi_pred;

MAE_phi  = mean(abs(phi_residual(:)));
RMSE_phi = sqrt(mean(phi_residual(:).^2));
MAX_phi  = max(abs(phi_residual(:)));

fprintf('Phase MAE  = %.4f rad\n', MAE_phi);
fprintf('Phase RMSE = %.4f rad\n', RMSE_phi);
fprintf('Phase MAX  = %.4f rad\n', MAX_phi);

figure('Name','Phase coupling check');

subplot(1,3,1);
imagesc(phi_map); axis image; colormap(jet); colorbar;
title('\phi input (measured)');

subplot(1,3,2);
imagesc(phi_pred); axis image; colormap(jet); colorbar;
title('\phi reconstructed from (n,t)');

subplot(1,3,3);
imagesc(phi_residual); axis image; colormap(jet); colorbar;
title('\phi residual = \phi_{in} - \phi_{pred}');

x = 1700;  % choose any x index within your image width

figure('Name','Phase profiles');
plot(phi_map(:,x), 'LineWidth', 1.5); hold on;
plot(phi_pred(:,x), '--', 'LineWidth', 1.5);
plot(phi_residual(:,x), ':', 'LineWidth', 1.5);
grid on;
xlabel('y (pixels)');
ylabel('Phase (rad)');
legend('\phi input','\phi reconstructed','residual');
title(sprintf('Phase coupling check at x = %d', x));


%% ==================== LOCAL FUNCTIONS ====================

function [n1_map, n2_map, t_map, res2_map] = cipher_solve_map(phi_map, lambda_um, nm, seeds, opts)

phi = double(phi_map);
sz  = size(phi);
N   = numel(phi);

n1_v  = NaN(N,1);
n2_v  = NaN(N,1);
t_v   = NaN(N,1);
res2_v= NaN(N,1);

solve1 = @(p) cipher_solve_pixel(p, lambda_um, nm, seeds, opts);

if isfield(opts,'use_parfor') && opts.use_parfor
    parfor ii = 1:N
        out = solve1(phi(ii));
        n1_v(ii)  = out.n1;
        n2_v(ii)  = out.n2;
        t_v(ii)   = out.t_um;
        res2_v(ii)= out.res2;
    end
else
    for ii = 1:N
        out = solve1(phi(ii));
        n1_v(ii)  = out.n1;
        n2_v(ii)  = out.n2;
        t_v(ii)   = out.t_um;
        res2_v(ii)= out.res2;
    end
end

n1_map  = reshape(n1_v,  sz);
n2_map  = reshape(n2_v,  sz);
t_map   = reshape(t_v,   sz);
res2_map= reshape(res2_v,sz);

end

function out = cipher_solve_pixel(phi, lambda_um, nm, seeds, opts)

if nargin < 5, opts = struct; end
opts = setdefaults(opts, struct( ...
    'mode','full_grid', ...
    'levels',3,'n1_step0',0.01,'n2_step0',0.001,'t_step0',0.01, ...
    'shrink',0.25,'window_factor',2, ...
    'use_gpu',false,'phi_thresh',0,'den_thresh',1e-9));

% quick reject: background pixels bias solutions toward bounds (Important for Analytic_t)
if ~isfinite(phi) || abs(phi) < opts.phi_thresh
    out = struct('n1',NaN,'n2',NaN,'t_um',NaN,'res2',NaN,'history',{{}});
    return;
end

% bounds
n1_lo = seeds.n1_min; n1_hi = seeds.n1_max;
n2_lo = seeds.n2_min; n2_hi = seeds.n2_max;
t_lo  = seeds.t_min;  t_hi  = seeds.t_max;

k = 2*pi/lambda_um;

best = struct('n1',NaN,'n2',NaN,'t',NaN,'res2',Inf);
histL = cell(0,1);

% initialize coarse grids
n1_c = n1_lo:opts.n1_step0:n1_hi; if numel(n1_c)<3, n1_c=linspace(n1_lo,n1_hi,3); end
n2_c = n2_lo:opts.n2_step0:n2_hi; if numel(n2_c)<3, n2_c=linspace(n2_lo,n2_hi,3); end
t_c  = t_lo :opts.t_step0 :t_hi;  if numel(t_c) <3, t_c =linspace(t_lo,t_hi,3);  end

for lev = 1:opts.levels

    % 2D grid over n1,n2
    [N2,N1] = ndgrid(n2_c, n1_c);
    A     = (N1 + N2/(lambda_um^2) - nm);   % effective index-diff term
    denom = k .* A;

    if opts.use_gpu
        denom = gpuArray(denom);
    end

    if strcmpi(opts.mode,'analytic_t')
        % --- analytic t* (previous shortcut, unstable now) ---
        t_star = phi ./ denom;
        mask_bad = abs(denom) < opts.den_thresh;
        t_star(mask_bad) = NaN;

        t_hat = min(max(t_star, t_lo), t_hi);
        phi_model = denom .* t_hat;
        res2 = (phi - phi_model).^2;
        res2(mask_bad) = Inf;

        % --- tie-break regularization (prevents bound-hugging) ---
        t_mid = 0.6;           % um (midpoint of your known range, not sure of this)
        alpha = 1e-6;          % small; tune 1e-8 to 1e-4
        res2  = res2 + alpha*(t_hat - t_mid).^2;

        [rmin, idx] = min(gather(res2(:)));
        if isfinite(rmin) && rmin < best.res2
            [i2,i1] = ind2sub(size(res2), idx);
            best.n1 = n1_c(i1);
            best.n2 = n2_c(i2);
            best.t  = gather(t_hat(i2,i1));
            best.res2 = rmin;
        end

        else
        % --- FAST FULL GRID (exact for uniform t grid, which is our case!!) ---
        % Find nearest t grid point to t_star = phi/denom (clipped)
        mask_bad = abs(denom) < opts.den_thresh;

        % uniform t step for this level
        if numel(t_c) >= 2
            t_step = t_c(2) - t_c(1);
        else
            t_step = (t_hi - t_lo);
        end

        t_star = phi ./ denom;
        t_star(mask_bad) = NaN;

        % clamp to bounds
        t_star = min(max(t_star, t_lo), t_hi);

        % quantize to nearest grid value (replicates brute min over t_c)
        % This is the key point
        t_hat = t_lo + round((t_star - t_lo) ./ t_step) .* t_step;
        t_hat = min(max(t_hat, t_lo), t_hi);
        t_hat(mask_bad) = NaN;

        % compute residual
        phi_model = denom .* t_hat;
        res2 = (phi - phi_model).^2;
        res2(mask_bad) = Inf;

        [rmin, idx] = min(gather(res2(:)));
        if isfinite(rmin) && rmin < best.res2
            [i2,i1] = ind2sub(size(res2), idx);
            best.n1 = n1_c(i1);
            best.n2 = n2_c(i2);
            best.t  = gather(t_hat(i2,i1));
            best.res2 = rmin;
        end
    end

    histL{end+1} = struct('lev',lev,'n1_grid',[n1_c(1) n1_c(end)], ...
                          'n2_grid',[n2_c(1) n2_c(end)],'t_grid',[t_c(1) t_c(end)], ...
                          'n1_len',numel(n1_c),'n2_len',numel(n2_c),'t_len',numel(t_c), ...
                          'res2_min',best.res2);

    % refine window around current best
    if lev < opts.levels && isfinite(best.res2)

        n1_step = max( mean(diff(n1_c)), (n1_hi-n1_lo)/max(2,numel(n1_c)-1) );
        n2_step = max( mean(diff(n2_c)), (n2_hi-n2_lo)/max(2,numel(n2_c)-1) );
        t_step  = max( mean(diff(t_c )), (t_hi -t_lo )/max(2,numel(t_c )-1) );

        n1_step = n1_step * opts.shrink;
        n2_step = n2_step * opts.shrink;
        t_step  = t_step  * opts.shrink;

        n1_lo_r = max(n1_lo, best.n1 - opts.window_factor*n1_step);
        n1_hi_r = min(n1_hi, best.n1 + opts.window_factor*n1_step);
        n2_lo_r = max(n2_lo, best.n2 - opts.window_factor*n2_step);
        n2_hi_r = min(n2_hi, best.n2 + opts.window_factor*n2_step);
        t_lo_r  = max(t_lo,  best.t  - opts.window_factor*t_step );
        t_hi_r  = min(t_hi,  best.t  + opts.window_factor*t_step );

        n1_c = n1_lo_r:n1_step:n1_hi_r; if numel(n1_c)<3, n1_c=linspace(n1_lo_r,n1_hi_r,3); end
        n2_c = n2_lo_r:n2_step:n2_hi_r; if numel(n2_c)<3, n2_c=linspace(n2_lo_r,n2_hi_r,3); end
        t_c  = t_lo_r :t_step :t_hi_r;  if numel(t_c) <3, t_c =linspace(t_lo_r, t_hi_r, 3); end
    end
end

out = struct('n1',best.n1,'n2',best.n2,'t_um',best.t,'res2',best.res2,'history',{histL});

end

function s = setdefaults(s, d)
fn = fieldnames(d);
for i=1:numel(fn)
    if ~isfield(s,fn{i}) || isempty(s.(fn{i}))
        s.(fn{i}) = d.(fn{i});
    end
end
end
