%{
    ----------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2025]
    Description:  [Asymptotic Similarity Transform for Toeplitz matrices]
    ----------------------------------------------------------------------
%}

clc
clear;
close all;

% ==== Parameters ====
    n = 23;             % Truncation size for a_k
    p = 6.5;            % Decay rate upwards
    q = 7.8;            % Decay rate downwards
    DimT = 80;          % Dimension of finite Toeplitz matrix to simulate open limit
    num_lambda = 101;   % Number of plotting points (keep it large)
    fs = 18;            % Fontsize for annotation
    

% ==== Generate Dummy Toeplitz Matrix ====
    col = zeros(n,1);
    row = zeros(1,n+1);
    
    % --- Add noise to the coeffiients ---
    ai = 1;
    bi = 1;

    col(1) = 1; 
    row(1) = 1;

    % --- Populate above and below diagonals ---
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


% ==== Approximate the open limit ====
    % --- Get initial guess inside the open limit ---
    lambda_start = ( min(real(eigT)) + max(real(eigT)) ) / 2;

    % --- Adaptive Computation of the Open Limit ---
     [lambda_interval(1), lambda_interval(2)]  = adaptive_equal_mod_interval(a, lambda_start);


% ==== Approximate conjugate root set Λ(f) ====
    % --- Get sampling which is close to the DOS (Uniform spacing on Λ(f) ---
    t = linspace(0, 1, num_lambda);
    lambda_vals = lambda_interval(1) + (lambda_interval(2) - lambda_interval(1)) * (0.5 * (1 - cos(2*pi*t)));

    %{
    histogram(lambda_vals, num_lambda, 'Normalization', 'pdf');
    xlabel('Eigenvalue');
    ylabel('Density of States');
    %}

    % --- Preallocate for all polynomial coefficients [P(z) = Q(z) - lambda*z^m] ---
    P_coeffs_matrix = repmat(a(:), 1, num_lambda);
    
    % The coefficient of z^m is at position (m+1) in the array and substract lambda
    P_coeffs_matrix(n + 1, :) = P_coeffs_matrix(n + 1, :) - lambda_vals;
    
    % --- Compute all roots at once ---
    all_roots = zeros(2*n, num_lambda);
    for k = 1:num_lambda
        all_roots(:, k) = roots(P_coeffs_matrix(:, k));
    end
    
    % --- Sort roots by magnitude for each lambda ---
    [~, sort_idx] = sort(abs(all_roots), 1);
    linear_idx = sort_idx + (0:num_lambda-1) * (2*n);
    all_roots_sorted = all_roots(linear_idx);
    
    % --- Extract m-th and (m+1)-th roots ---
    candidate_roots = all_roots_sorted(n:n+1, :);

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

    openLimit = merge_close_points(openLimit, 1e-2);

    phase = angle(openLimit(:));
    [phase_sorted, sortIdx] = sort(phase);
    openLimit_sorted = openLimit(sortIdx);

    % --- Plot the set Λ(f) ---
    %{
        wraparound_OpenLimit = [openLimit_sorted, openLimit_sorted(1)];
        figure;
        plot(real(wraparound_OpenLimit), imag(wraparound_OpenLimit), 'k-', 'LineWidth', 2.5)
        hold on;
        plot(real(openLimit_sorted), imag(openLimit_sorted), 'bx', 'LineWidth', 2.5)
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
        grid on;
        axis equal;
        box on;
        hold off;
    %}


% ==== Compute Fourier coefficients of f(p(z)) numerically  ====
    % --- Evaluate f(p(z)) on the torus ---
    k_values = -n:n;
    powers_matrix = openLimit_sorted(:).^(-k_values);  % N x (2n+1) matrix
    
    % --- Vectorized sum ---
    fp_values = powers_matrix * a(:);

    % --- Clean up data ---
    fp_values = real(fp_values);

    % --- Wrap around ---
    phase_sorted = [-pi; phase_sorted; pi];
    fp_values    = [ fp_values(end); fp_values; fp_values(end)];

    % --- Plot the function f(p(torus))
    %{
    figure;
    plot(phase_sorted, fp_values);
    hold off;
    %}

    % --- Compute the Fourier Transform of f(p(z)) ---
    F_range = n+8;
    FourierFP = fourier_coefficients_spectral(phase_sorted, fp_values, F_range);


% ==== Quasi Similarity Transformed Toeplitz matrix ====
    % --- Toeplitz matrix for deformed path ---
    T_b = fourier_to_toeplitz(FourierFP, DimT);
    eigT_b = sort(eig(T_b));

    % --- Plot the eigenvalues before and after asymptotic Similarity transform ---
    %
    figure;
    plot(real(eigT), imag(eigT), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold on;
    plot(real(eigT_b), imag(eigT_b), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    grid on;
    box on;
    xlabel('$\mathrm{Re}(\sigma(\mathbf{T}_N))$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}(\sigma(\mathbf{T}_N))$', 'Interpreter', 'latex', 'FontSize', 14); 
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 
    axis equal;
    hold off;
    %}

    %{
    histogram(eigT_b, floor(DimT /10), 'Normalization', 'pdf');
    xlabel('Eigenvalue');
    ylabel('Density of States');
    %}


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

    N = length(x);
    %M = length(x_target);

    % Precompute barycentric weights (independent of x_target)
    w = (-1).^(0:N-1)';   % column vector

    % Build matrix of differences: M×N
    % each row j is x_target(j) - x(:)'
    D = x_target(:) - x(:).';

    % Detect exact hits (machine precision)
    hit = abs(D) < 1e-14;

    % Avoid division by zero in denominator
    D(hit) = Inf;

    % Barycentric kernel: w / sin((x_target - x)/2)
    K = w.' ./ sin(D/2);     % w.' ensures broadcasting over rows

    % Compute numerator and denominator
    num = K * f;             % M×N times N×1 → M×1
    den = sum(K, 2);         % sum rows → M×1

    f_interp = num ./ den;

    % Insert exact-match values
    if any(hit, "all")
        [j_idx, x_idx] = find(hit);
        f_interp(j_idx) = f(x_idx);
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
    tol = 1e-9;

    % initial step size
    step0 = 0.1;
    min_step = 1e-10;

    % determine m
    m = (length(a)-1)/2;

    % build polynomial template
    coeff_Q = zeros(1,2*m+1);
    for j = -m:m
        coeff_Q(m-j+1) = a(j+m+1);
    end

    % helper: m-th and (m+1)-th modulus difference
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
    % Sort the array to simplify merging
    openLimit = sort(openLimit);
    
    merged = [];
    i = 1;
    while i <= length(openLimit)
        % Start a group with the current point
        group = openLimit(i);
        j = i + 1;
        
        % Collect all points close to the current point
        while j <= length(openLimit) && abs(openLimit(j) - openLimit(i)) <= tol
            group(end+1) = openLimit(j);
            j = j + 1;
        end
        
        % Store the average of the group
        merged(end+1) = mean(group);
        
        % Move to the next unprocessed point
        i = j;
    end
end

