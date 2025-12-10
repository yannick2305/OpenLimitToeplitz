%{
    ----------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2025]
    Description:  [Homotopy and Symbol collapse]
    ----------------------------------------------------------------------
%}

clc
clear;
close all;

% ==== Parameters ====
    m = 9;              % Truncation size for a_k
    p = 3.5;            % Decay rate upwards
    q = 4.8;            % Decay rate downwards
    DimT = 100;         % Dimension of finite Toeplitz matrix to simulate open limit
    num_lambda = 100;   % Number of plotting points (50-300)
    fs = 18;            % Fontsize for annotation
    

% ==== Generate m-bsned Dummy Toeplitz Matrix ====
    col = zeros(m,1);
    row = zeros(1,m+1);
    
    % --- Add noise to the coeffiients ---
    ai = 1;
    bi = 1;

    col(1) = 1; 
    row(1) = 1;

    % --- Populate above and below diagonals ---
    for k = 1:m
        r = ai + (bi-ai)*rand;
        col(k) = r / ((k+1)^q);
        r = ai + (bi-ai)*rand;
        row(k+1) = r / ((k+1)^p);
    end
   
    % --- Coefficients of the symbol function ---
    a = [ col(end:-1:1)', row];
       
    % --- Generate finite Toeplitz matrix ---
    T = fourier_to_toeplitz(a, DimT);
    eigT = sort(eig(T));

% ==== Approximate the open limit ====
    % --- Get initial guess inside the open limit ---
    lambda_start = ( min(real(eigT)) + max(real(eigT)) ) / 2;

    % --- Adaptive Computation of the Open Limit ---
     [lambda_interval(1), lambda_interval(2)]  = adaptive_equal_mod_interval(a, lambda_start);

% ==== Approximate conjugate root set Λ(f) ====
    % --- Get sampling which is close to the DOS (Uniform spacing on Λ(f) ---
    t = linspace(0, 1, num_lambda);
    lambda_vals = lambda_interval(1) + (lambda_interval(2) - lambda_interval(1)) * (0.5 * (1 - cos(pi*t)));
    
    % --- Plot the Density of states of the sampling ---
    %{
    histogram(lambda_vals, floor(num_lambda/4), 'Normalization', 'pdf');
    xlabel('Eigenvalue');
    ylabel('Density of States');
    %}

    % --- Preallocate for all polynomial coefficients [P(z) = Q(z) - lambda*z^m] ---
    P_coeffs_matrix = repmat(a(:), 1, num_lambda);
    
    % The coefficient of z^m is at position (m+1) in the array and substract lambda
    P_coeffs_matrix(m + 1, :) = P_coeffs_matrix(m + 1, :) - lambda_vals;
    
    % --- Compute all roots at once ---
    all_roots = zeros(2*m, num_lambda);
    for k = 1:num_lambda
        all_roots(:, k) = roots(P_coeffs_matrix(:, k));
    end
    
    % --- Sort roots by magnitude for each lambda ---
    [~, sort_idx] = sort(abs(all_roots), 1);
    linear_idx = sort_idx + (0:num_lambda-1) * (2*m);
    all_roots_sorted = all_roots(linear_idx);
    
    % --- Extract m-th and (m+1)-th roots ---
    candidate_roots = all_roots_sorted(m:m+1, :);

    % --- Filter to keep only pairs with similar modulus ---
    tolerance = 1e-8;  % keep in mind rootsolver has its limits
    mod_n  = abs(candidate_roots(1, :));
    mod_n1 = abs(candidate_roots(2, :));
    
    % Find pairs where moduli are approximately equal
    similar_modulus = abs(mod_n - mod_n1) ./ max(mod_n, mod_n1) < tolerance;
    
    % Keep only the filtered roots
    candidate_roots = candidate_roots(:, similar_modulus);

    % --- Check if they have approximately the same modulus ---
    openLimit = [candidate_roots(1, :), candidate_roots(2, :)];

    phase = angle(openLimit(:));
    [~, sortIdx] = sort(phase);
    openLimit_sorted = openLimit(sortIdx);

    % --- Remove duplicated as they mess up the weights in integral ---
    openLimit_sorted = merge_close_points(openLimit_sorted, 1e-4);
    phase_sorted = angle(openLimit_sorted(:));

    theta = [phase_sorted; phase_sorted(1)];
    fp    = [openLimit_sorted, openLimit_sorted(1)];


    % --- Homotopy to Λ(f) ---
    create_homotopy_movie(theta, fp, a, T, m, 'Filename', 'CollapsedSymbolHomotopy');


    %% --- Homotopy to scaled unit torus ---
    theta_unit_circle    = linspace(0, 2*pi, 100);
    r_scaled_unit_circle = 1.8*exp(1i*theta_unit_circle);

    create_homotopy_movie(theta_unit_circle, r_scaled_unit_circle, a, T, m, 'Filename', 'ScaledUnitCircle');
   

%%
function homotopy_movie = create_homotopy_movie(theta, fp, a, T, m, varargin)

p = inputParser;
addRequired(p, 'theta');
addRequired(p, 'fp');
addRequired(p, 'a');
addRequired(p, 'm');
addParameter(p, 'Filename', 'CollapsedSymbolHomotopy');
parse(p, theta, fp, a, m, varargin{:});

num_frames = 200;
filename = p.Results.Filename;
frame_rate = 10;

% Ensure column vectors
theta = theta(:);
fp = fp(:);
a = a(:);

% Validate Laurent polynomial coefficients
if length(a) ~= 2*m + 1
    error('Coefficient vector a must have length 2*m+1 = %d', 2*m+1);
end

% Sort by angle for proper curve plotting
[theta_sorted, sort_idx] = sort(theta);
fp_sorted = fp(sort_idx);

% Wraparound
theta_dense = theta_sorted; %[pi; theta_sorted]; 
fp_interp   = fp_sorted; %[fp_sorted(1); fp_sorted]; 

% Target function: e^(i*theta)
target = exp(1i * theta_dense);

fig = figure('Position', [100, 100, 1600, 500]);
set(fig, 'Color', 'w');

% Initialize movie
if ~isempty(filename)
    v = VideoWriter(filename, 'MPEG-4');
    v.FrameRate = frame_rate;
    open(v);
end
homotopy_movie(num_frames) = struct('cdata', [], 'colormap', []);

% Compute functional values over homotopy
t_values = linspace(1, 0, num_frames);
Q_values = zeros(num_frames, 1);
Q_curves = cell(num_frames, 1);

fprintf('Computing homotopy and Laurent polynomial...\n');
for k = 1:num_frames
    t = t_values(k);
    % Linear homotopy: H(theta, t) = (1-t)*fp(theta) + t*e^(i*theta)
    h_t = (1 - t) * fp_interp + t * target;
    
    % Evaluate Laurent polynomial Q(h_t)
    Q_h_t = evaluate_laurent_poly(h_t, a, m);
    Q_curves{k} = Q_h_t;
    
    % Compute functional: integral of Q over the curve
    Q_values(k) = trapz(theta_dense, Q_h_t) / (2*pi);
end

% Create animation
for k = 1:num_frames
    t = t_values(k);
    h_t = (1 - t) * fp_interp + t * target;
    Q_h_t = Q_curves{k};
    
% ==== First Subplot ====
    subplot(1, 2, 1);
    % Subplot 1: Curve in complex plane
    hold off;
    
    % Plot target circle
    plot(real(target), imag(target), 'k--', 'LineWidth', 1.5);
    hold on;
    
    % Plot current homotopy curve
    plot(real(h_t), imag(h_t), 'b-', 'LineWidth', 2.5);
    
    % Plot original sampled points
    %scatter(real(fp), imag(fp), 50, theta, 'filled', 'MarkerEdgeColor', 'k');
    
    % Plot current interpolated points (colored by angle)
    %scatter(real(h_t), imag(h_t), 20, theta_dense, 'filled', 'MarkerFaceAlpha', 0.3);
    
    axis equal;
    grid on;
    max_radius = max(abs([fp; target]));
    xlim([-1.5, 1.5] * max_radius);
    ylim([-1.5, 1.5] * max_radius);
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 18); 
    legend({'$e^{i\theta}$', '$p(e^{i\theta})$'}, 'Interpreter', 'latex', 'FontSize', 28);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 24);
   
% ==== Second Subplot ====
    subplot(1, 2, 2);

    hold off;
        
    % Plot Q(h_t) curve
    plot(real(Q_h_t), imag(Q_h_t), 'k-', 'LineWidth', 2.5);
    hold on;
        
    grid on;
    ylim([-0.07; 0.07]);
    xlim([min(real(eig(T)))-0.05, max(real(eig(T)))+0.05])
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 18); 
    legend({'$f\big(p(e^{i\theta})\big)$'}, 'Interpreter', 'latex', 'FontSize', 28);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 24);

    drawnow;
    
    % Capture frame
    homotopy_movie(k) = getframe(fig);
    if ~isempty(filename)
        writeVideo(v, homotopy_movie(k));
    end
end

if ~isempty(filename)
    close(v);
    fprintf('Movie saved to: %s\n ', filename);
end

gif_filename = 'homotopyCircle.gif';
frame_delay = 1 / frame_rate;

for k = 1:num_frames
    [img, cmap] = rgb2ind(frame2im(homotopy_movie(k)), 256);

    if k == 1
        imwrite(img, cmap, gif_filename, ...
                'gif', 'LoopCount', Inf, ...
                'DelayTime', frame_delay);
    else
        imwrite(img, cmap, gif_filename, ...
                'gif', 'WriteMode', 'append', ...
                'DelayTime', frame_delay);
    end
end


end

% ==== Other helper Functions ====

function Q_z = evaluate_laurent_poly(z, a, m)

    Q_z = zeros(size(z));
    
    for j = -m:m
        coeff = a(j + m + 1);
        Q_z = Q_z + coeff * z.^(-j);
    end
end


function T = fourier_to_toeplitz(a, dimT)
    K = (length(a) - 1) / 2;
    a_0 = a(K+1);          
    
    col = zeros(dimT, 1);
    col(1) = a_0;
    row = zeros(1, dimT);
    row(1) = a_0;

    for k = 1:min(K, dimT-1)
        col(k+1) = a(K+1-k);
        row(k+1) = a(K+1+k);
    end
   
    T = toeplitz(col, row);

end


function [lambdaL, lambdaR] = adaptive_equal_mod_interval(a, lambda0)

    % tolerance for detecting equality
    tol = 1e-12;

    % initial step size
    step0 = 0.1;
    min_step = 1e-12;

    % determine m
    m = (length(a)-1)/2;

    % build polynomial template
    coeff_Q = zeros(1,2*m+1);
    for j = -m:m
        coeff_Q(m-j+1) = a(j+m+1);
    end

    % helper: m-th and (m+1)-th modulus difference

    %
    function F = diff_mod(lambda)
        c = coeff_Q;
        c(m+1) = c(m+1) - lambda;
        r = roots(c);
        r = sort(abs(r));
        F = r(m+1) - r(m);
    end
    

    % -------------------------------
    % Find RIGHT endpoint
    % -------------------------------
    lambda = lambda0;
    step = step0;

    % ensure we start inside interval
    if abs(diff_mod(lambda)) > tol
        error('Initial guess lambda0 is not inside equal-modulus interval.');
    end

    while true
        % tentative step
        lambda_try = lambda + step;
        F_try = diff_mod(lambda_try);

        if abs(F_try) <= tol
            % still inside, check if stepping further breaks equality
            lambda_test = lambda_try + tol;
            if abs(diff_mod(lambda_test)) > tol
                lambdaR = lambda_try;
                break;
            else
                lambda = lambda_try;   % move forward
            end

        else
            % outside ⇒ shrink step
            step = step / 2;
            if step < min_step
                lambdaR = lambda_try;
                break;
            end
        end
    end

    % -------------------------------
    % Find LEFT endpoint
    % -------------------------------
    lambda = lambda0;
    step = step0;

    while true
        lambda_try = lambda - step;
        F_try = diff_mod(lambda_try);

        if abs(F_try) <= tol
            lambda_test = lambda_try - tol;
            if abs(diff_mod(lambda_test)) > tol
                lambdaL = lambda_try;
                break;
            else
                lambda = lambda_try;
            end

        else
            % outside ⇒ shrink step
            step = step / 2;
            if step < min_step
                lambdaL = lambda_try;
                break;
            end
        end
    end
end


function merged = merge_close_points(openLimit, tol)
    
    n = length(openLimit);
    if n == 0
        merged = [];
        return;
    end
    
    merged = openLimit;
    
    % Check if first and last elements should merge
    if abs(openLimit(end) - openLimit(1)) <= tol
        merged(1) = mean([openLimit(1), openLimit(end)]);
        merged(end) = [];  % Remove last element
    end
    
    % Check if two middle elements should merge
    if n >= 2
        mid1 = floor(n / 2);
        mid2 = mid1 + 1;
        
        if abs(merged(mid2) - merged(mid1)) <= tol
            merged(mid1) = mean([merged(mid1), merged(mid2)]);
            merged(mid2) = [];  % Remove second middle element
        end
    end
end

%{
function g_val = g(lambda, a)
    % Compute g(lambda) = |a(end)| * prod_{k=m+1}^{2m} |z_k(lambda)|
    %
    % Inputs:
    %   lambda - complex value
    %   a - coefficient array (must have odd length 2m+1)
    
    m = (length(a)-1)/2;
    
    % Build polynomial template
    coeff_Q = zeros(1, 2*m+1);
    for j = -m:m
        coeff_Q(m-j+1) = a(j+m+1);
    end
    
    % Modify polynomial: Q(z) - lambda
    c = coeff_Q;
    c(m+1) = c(m+1) - lambda;
    
    % Find roots and sort by magnitude
    r = roots(c);
    r_abs = abs(r);
    [~, idx] = sort(r_abs);
    r_sorted = r(idx);
    
    % Compute g(lambda) = |a(end)| * product of larger m roots
    larger_roots = r_sorted(m+1:2*m);
    g_val = abs(a(end)) * prod(abs(larger_roots));
end
%}


