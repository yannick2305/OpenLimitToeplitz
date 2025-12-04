%{
    ----------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2025]
    Description:  [Asymptotic Similarity Transform for Toeplitz matrices]
    ----------------------------------------------------------------------
%}

clc
clear all;
close all;

% --- Parameters ---
    n = 9;              % Truncation size for a_k
    p = 4.5;            % Decay rate upwards
    q = 4.5;            % Decay rate downwards
    DimT = 100;         % Dimension of finite Toeplitz matrix to simulate open limit
    num_lambda = 51;    % Number of plotting points (keep it large)
    fs = 18;            % Fontsize for annotation
    
% --- Generate sequence a_k ---
    col = zeros(n,1);
    row = zeros(1,n+1);
    
    % --- Add noise to the coeffiients ---
    ai = 1;
    bi = 2;

    col(1) = 1; % Diagonal element
    row(1) = 1;

    % Downward decay (sub-diagonals)
    for k = 1:n
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

    % --- Get initial Guess in the open limit ---
    lambda_range = [min(real(eigT)), max(real(eigT))];
    lambda_start = ( lambda_range(end) + lambda_range(1) ) / 2;

    % --- Open Limit ---
    [lambda_interval, diagnostics] = find_conjugate_interval(a, lambda_start);

% --- Approximate conjugate root set ---
    % --- Generate frequency range ---
    lambda_vals = linspace(lambda_interval(1), lambda_interval(2), num_lambda);

    % --- Preallocate for all polynomial coefficients ---
    % P(z) = Q(z) - lambda*z^m
    % We need to subtract lambda from the coefficient of z^m
    P_coeffs_matrix = repmat(a(:), 1, num_lambda);
    
    % The coefficient of z^m is at position (m+1) in the array
    % Subtract lambda values from that position
    P_coeffs_matrix(n + 1, :) = P_coeffs_matrix(n + 1, :) - lambda_vals;
    
    % --- Compute all roots at once ---
    all_roots = zeros(2*n, num_lambda);
    for k = 1:num_lambda
        all_roots(:, k) = roots(P_coeffs_matrix(:, k));
    end
    
    % --- Sort roots by magnitude for each lambda ---
    [~, sort_idx] = sort(abs(all_roots), 1);
    % Convert linear indices for sorted access
    linear_idx = sort_idx + (0:num_lambda-1) * (2*n);
    all_roots_sorted = all_roots(linear_idx);
    
    % --- Extract m-th and (m+1)-th roots ---
    candidate_roots = all_roots_sorted(n:n+1, :);
    
    % --- Check if they have approximately the same modulus ---
    openLimit = [candidate_roots(1, :), candidate_roots(2, :)];


% --- Plot the set Λ(f) ---
    figure;
    plot(real(openLimit), imag(openLimit), 'kx', 'LineWidth', 2.5)
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    grid on;
    axis equal;
    box on;
    hold off;


% --- Compute Fourier coefficients of f(p(z)) numerically  ---

    % --- Retrieve the phase ---
    phase = zeros(length(openLimit),1);
    for i = 1:length(openLimit)
        phase(i) = angle(openLimit(i));
    end
    
    % --- Evaluate f(p(z)) on the torus ---
    N = length(openLimit);
    fp_values = zeros(N, 1);
    
    for j = 1:N
        for k = -n:n
            k_idx = k + n + 1; 
            fp_values(j) = fp_values(j) + a(k_idx) * openLimit(j)^(-k);
        end
    end

    fp_values = real(fp_values);

    [phase_sorted, sort_idx] = sort(phase);
    fp_sorted = fp_values(sort_idx);
    
    [phase_unique, unique_idx] = unique(phase_sorted); 
    fp_unique = fp_sorted(unique_idx);  
    
    phase_unique = [-pi; phase_unique];   
    fp_unique    = [fp_unique(end); fp_unique];       

 
    figure;
    plot(phase_unique, fp_unique);
    hold off;

    % --- Compute the Fourier Transform of f(p(z)) ---

    F_range = n;
    FourierFP = fourier_coefficients_spectral(phase_unique, fp_unique, F_range);

    % --- Plot the decay of the Fouier Coefficients of f(p(z)) ---
    figure;
    loglog(1:F_range, abs(FourierFP(F_range+1: 2*F_range)), 'k.', 'LineWidth', 1.5)
    hold on
    plot(row)
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    grid on;
    hold off;

    % --- Toeplitz matrix for deformed path ---
    T_b = fourier_to_toeplitz(FourierFP, DimT);
    eigT_b = sort(eig(T_b));
    similar = norm(eigT - eigT_b);
    fprintf('Similarity Coefficient: %f\n', similar);

% --- Plot the eigenvalues before and after asymptotic Similarity transform ---
    figure;
    plot(real(eigT), imag(eigT), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold on;
    plot(real(eigT_b), imag(eigT_b), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    grid on;
    box on;
    xlabel('$\mathrm{Re}(\sigma(\mathbf{T}_N))$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}(\sigma(\mathbf{T}_N))$', 'Interpreter', 'latex', 'FontSize', 14); 
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 250]); 
    axis equal;
    hold off;



%% 

    % Diagonalize T: T = P * D * P^(-1)
    [P, D] = eig(T);
    [D_sorted, sort_idx] = sort(diag(D), 'descend'); % or 'ascend'
    D = diag(D_sorted);
    P = P(:, sort_idx);
    
    % Diagonalize T_b: T_b = Q * D_b * Q^(-1)
    [Q, D_b] = eig(T_b);
    [D_b_sorted, sort_idx_b] = sort(diag(D_b), 'descend'); % or 'ascend'
    D_b = diag(D_b_sorted);
    Q = Q(:, sort_idx_b);
    
    % Check that eigenvalues are the same
    disp('Are eigenvalues equal?');
    disp(['D diagonal  : ', num2str(real(sort(diag(D)))')]);
    disp(['D_b diagonal: ', num2str(real(sort(diag(D_b)))')]);
    
    fprintf('\n');
    invP = P \ eye(DimT);
    % Compute the similarity transformation matrix
    S = invP * Q;
    invS = S \ eye(DimT);

    % Verify: T = S * T_b * S^(-1)
    T_from_Tb = S * T_b * invS;
    MatrixDiff = T - T_from_Tb;


%% --- Defining functions ---

function ck = fourier_coefficients_spectral(phase, fp_values, K)
    % First, interpolate using barycentric trigonometric interpolation
    % to get uniformly-spaced values, then use FFT
    
    N = length(phase);
    M = max(2*K+1, 2*N);  % target uniform grid size
    
    % Create uniform grid
    theta_uniform = (0:M-1)' * 2*pi/M;
    
    % Barycentric trigonometric interpolation
    fp_uniform = trig_barycentric_interp(phase, fp_values, theta_uniform);
    
    % Now use standard FFT (spectral accuracy!)
    fft_result = fft(fp_uniform) / M;
    
    % Extract coefficients
    ck = zeros(2*K+1, 1);
    ck(K+1) = fft_result(1);
    ck(K+2:end) = fft_result(2:K+1);
    ck(1:K) = fft_result(end-K+1:end);
end

function f_interp = trig_barycentric_interp(x, f, x_target)
    % Barycentric formula for trigonometric interpolation
    N = length(x);
    M = length(x_target);
    f_interp = zeros(M, 1);
    
    for j = 1:M
        if any(abs(x_target(j) - x) < 1e-14)
            % Exact node
            idx = find(abs(x_target(j) - x) < 1e-14, 1);
            f_interp(j) = f(idx);
        else
            % Barycentric interpolation
            weights = (-1).^(0:N-1)' ./ sin((x_target(j) - x)/2);
            f_interp(j) = sum(weights .* f) / sum(weights);
        end
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


function [lambda_interval, diagnostics] = find_conjugate_interval(a, lambda_init, options)
    % Find the interval where m-th and (m+1)-th largest roots transition from phase 0 to ±π
    % Assumes roots are conjugates throughout and phase is monotone in lambda
    %
    % Inputs:
    %   a            - Coefficient vector of length 2*m+1
    %   lambda_init  - Initial guess for lambda (should be inside the interval)
    %   options      - (optional) struct with fields:
    %                  .step_size (default 0.1) - initial step size
    %                  .tol (default 1e-10) - tolerance for convergence
    %                  .check_radius (default 0.01) - radius for checking if phase is stuck
    %
    % Outputs:
    %   lambda_interval - [lambda_start, lambda_end] 
    %   diagnostics     - struct with additional information
    
    % Set default options
    if nargin < 3
        options = struct();
    end
    if ~isfield(options, 'step_size')
        options.step_size = 0.1;
    end
    if ~isfield(options, 'tol')
        options.tol = 1e-10;
    end
    if ~isfield(options, 'check_radius')
        options.check_radius = 1e-8;
    end
    
    % Determine m from the length of a
    m = (length(a) - 1) / 2;
    if mod(length(a), 2) == 0 || m ~= floor(m)
        error('Vector a must have odd length (2*m+1)');
    end
    
    % Build coeff_Q from a
    coeff_Q = zeros(1, 2*m + 1);
    for j = -m:m
        coeff_Q(m - j + 1) = a(j + m + 1);
    end
    
    idx_lambda_m = m + 1;  % Index where lambda appears
    
    % Helper function: compute phase of m-th largest root
    function [phase_val, root_m] = compute_root_phase(lambda)
        % Build polynomial coefficients for this lambda
        coeffs = coeff_Q;
        coeffs(idx_lambda_m) = coeffs(idx_lambda_m) - lambda;
        
        % Find roots
        r = roots(coeffs);
        
        % Sort by magnitude (ascending)
        [~, sort_idx] = sort(abs(r));
        r_sorted = r(sort_idx);
        
        % Pick m-th largest (end is largest)
        root_m = r_sorted(end - m + 1);
        
        % Compute phase
        phase_val = angle(root_m);
    end
    
    % Helper to compute distance to target phase
    function dist = phase_distance(phase_val, target)
        if abs(target) < 0.1
            % Target is 0
            dist = abs(phase_val);
        else
            % Target is ±π
            dist = min(abs(phase_val - pi), abs(phase_val + pi));
        end
    end
    
    % Helper to check if phase is stuck at boundary
    function is_stuck = check_if_stuck(lambda, target_phase)
        [phase_center, ~] = compute_root_phase(lambda);
        dist_center = phase_distance(phase_center, target_phase);
        
        % If not close to target, definitely not stuck
        if dist_center > options.tol * 10
            is_stuck = false;
            return;
        end
        
        % Check in a small neighborhood
        [phase_plus, ~] = compute_root_phase(lambda + options.check_radius);
        [phase_minus, ~] = compute_root_phase(lambda - options.check_radius);
        
        dist_plus = phase_distance(phase_plus, target_phase);
        dist_minus = phase_distance(phase_minus, target_phase);
        
        % If phase doesn't change in neighborhood, we're stuck (past boundary)
        if abs(dist_plus - dist_center) < options.tol && abs(dist_minus - dist_center) < options.tol
            is_stuck = true;
        else
            is_stuck = false;
        end
    end
    
    % Step 1: Evaluate phase at initial guess
    fprintf('Initial guess lambda = %.6f\n', lambda_init);
    [phase_init, ~] = compute_root_phase(lambda_init);
    fprintf('Phase at initial guess: %.6f\n', phase_init);
    
    % Determine which direction leads to which target
    delta = options.check_radius;
    [phase_plus, ~] = compute_root_phase(lambda_init + delta);
    [phase_minus, ~] = compute_root_phase(lambda_init - delta);
    
    dist_plus_to_zero = phase_distance(phase_plus, 0);
    dist_plus_to_pi = phase_distance(phase_plus, pi);
    dist_minus_to_zero = phase_distance(phase_minus, 0);
    dist_minus_to_pi = phase_distance(phase_minus, pi);
    
    % Determine directions
    if dist_plus_to_pi < dist_minus_to_pi
        % Increasing lambda goes toward ±π
        dir_to_pi = +1;
        dir_to_zero = -1;
    else
        % Decreasing lambda goes toward ±π
        dir_to_pi = -1;
        dir_to_zero = +1;
    end
    
    fprintf('Direction toward ±π: %+d\n', dir_to_pi);
    fprintf('Direction toward 0: %+d\n', dir_to_zero);
    
    % Step 2: Find boundary for phase ±π
    fprintf('\nSearching boundary where phase ≈ ±π...\n');
    lambda_pi = find_boundary_stepping(lambda_init, dir_to_pi, pi);
    
    % Step 3: Find boundary for phase 0
    fprintf('\nSearching boundary where phase ≈ 0...\n');
    lambda_zero = find_boundary_stepping(lambda_init, dir_to_zero, 0);
    
    % Order them correctly
    lambda_left = min(lambda_pi, lambda_zero);
    lambda_right = max(lambda_pi, lambda_zero);
    lambda_interval = [lambda_left, lambda_right];
    
    fprintf('\n=== Results ===\n');
    fprintf('Interval: [%.10f, %.10f]\n', lambda_left, lambda_right);
    
    % Verify
    [phase_left, ~] = compute_root_phase(lambda_left);
    [phase_right, ~] = compute_root_phase(lambda_right);
    fprintf('Verification - Left phase:  %.8f\n', phase_left);
    fprintf('Verification - Right phase: %.8f\n', phase_right);
    
    % Store diagnostics
    diagnostics.lambda_init = lambda_init;
    diagnostics.phase_init = phase_init;
    diagnostics.phase_left = phase_left;
    diagnostics.phase_right = phase_right;
    
    % Nested function: step in direction until phase reaches target
    function lambda_boundary = find_boundary_stepping(lambda_start, direction, target_phase)
        % direction: +1 to increase lambda, -1 to decrease lambda
        % target_phase: 0 or pi
        
        lambda_current = lambda_start;
        step = direction * options.step_size;
        
        % Step until we get close to the target
        max_steps = 10000;
        for i = 1:max_steps
            [phase_current, ~] = compute_root_phase(lambda_current);
            dist_current = phase_distance(phase_current, target_phase);
            
            % Check if we're stuck at the boundary
            if dist_current < options.tol
                if check_if_stuck(lambda_current, target_phase)
                    fprintf('  Warning: Phase stuck at boundary at lambda = %.10f\n', lambda_current);
                    fprintf('  Stepping back to find actual boundary...\n');
                    % Step back until phase starts changing
                    lambda_back = lambda_current - step;
                    for j = 1:100
                        if ~check_if_stuck(lambda_back, target_phase)
                            fprintf('  Found non-stuck point at lambda = %.10f\n', lambda_back);
                            % Bisect between this point and the stuck point
                            lambda_boundary = bisection_refine(lambda_back, lambda_current, target_phase);
                            return;
                        end
                        lambda_back = lambda_back - step;
                    end
                    warning('Could not find non-stuck point');
                    lambda_boundary = lambda_current;
                    return;
                else
                    fprintf('  Converged at lambda = %.10f (phase = %.8f) after %d steps\n', ...
                            lambda_current, phase_current, i);
                    lambda_boundary = lambda_current;
                    return;
                end
            end
            
            % Check if next step will overshoot
            lambda_next = lambda_current + step;
            [phase_next, ~] = compute_root_phase(lambda_next);
            dist_next = phase_distance(phase_next, target_phase);
            
            if dist_next < dist_current && dist_current < 0.1
                % Getting close, switch to smaller steps
                step = step * 0.5;
            elseif dist_next > dist_current && dist_current < 0.1
                % Overshot, refine with bisection
                fprintf('  Switching to bisection refinement...\n');
                lambda_boundary = bisection_refine(lambda_current, lambda_next, target_phase);
                return;
            end
            
            lambda_current = lambda_next;
            
            if mod(i, 100) == 0
                fprintf('  Step %d: lambda = %.6f, phase = %.6f, dist = %.6f\n', ...
                        i, lambda_current, phase_current, dist_current);
            end
        end
        
        warning('Maximum steps reached without convergence');
        lambda_boundary = lambda_current;
    end
    
    % Nested function: bisection refinement
    function lambda_boundary = bisection_refine(lambda_a, lambda_b, target_phase)
        max_iter = 100;
        for iter = 1:max_iter
            lambda_mid = (lambda_a + lambda_b) / 2;
            [phase_mid, ~] = compute_root_phase(lambda_mid);
            dist_mid = phase_distance(phase_mid, target_phase);
            
            % Check if stuck
            if dist_mid < options.tol && check_if_stuck(lambda_mid, target_phase)
                % Mid is stuck, boundary is between a and mid
                lambda_b = lambda_mid;
            elseif dist_mid < options.tol || abs(lambda_b - lambda_a) < 1e-12
                lambda_boundary = lambda_mid;
                fprintf('  Bisection converged: lambda = %.10f (phase = %.8f)\n', ...
                        lambda_mid, phase_mid);
                return;
            else
                % Standard bisection logic
                [phase_a, ~] = compute_root_phase(lambda_a);
                dist_a = phase_distance(phase_a, target_phase);
                
                if dist_mid < dist_a
                    lambda_a = lambda_mid;
                else
                    lambda_b = lambda_mid;
                end
            end
        end
        
        lambda_boundary = (lambda_a + lambda_b) / 2;
    end
end
