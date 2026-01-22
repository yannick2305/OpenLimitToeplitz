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

itera = 1000;

ero = zeros(itera,1);

W1 = zeros(itera,1);

FourCo= zeros(itera,1);


for num_lambda = 5:itera
% ==== Parameters ====
    m = 5;              % Truncation size for a_k
    p = 3.5;            % Decay rate upwards
    q = 3.8;            % Decay rate downwards
    %DimT = 70;         % Dimension of finite Toeplitz matrix to simulate open limit
    %num_lambda = 1500;   % Number of plotting points (50-300)
    fs = 18;            % Fontsize for annotation

    DimT = 50;

    

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

    %a = [1 ,0, 1];
       
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

    FourierP = fourier_coefficients_spectral(phase_sorted, openLimit_sorted.', 10);

    theta = phase_sorted(:);       % ensure column vector
    p_vals = openLimit_sorted(:);  % ensure column vector
    N = length(theta);
    
    % Sort by theta (optional, helps trapezoid)
    [theta_sorted, idx] = sort(theta);
    p_sorted = p_vals(idx);

%{
% Fourier modes
M = 10;
n = -M:M;
c = zeros(length(n),1);

% Compute approximate Fourier coefficients using trapezoidal rule
for j = 1:length(n)
    % differences between consecutive theta for integration weight
    dtheta = diff([theta_sorted; theta_sorted(1)+2*pi]); 
    c(j) = sum( p_sorted .* exp(-1i*n(j)*theta_sorted) .* dtheta ) / (2*pi);
end


k = M;
m = -k:k;

E = exp(-1i * theta_sorted(:) * m);   % size N × (2k+1)
ct = (p_sorted(:).' * E).' / length(p_sorted);


    figure;
    semilogy(1:length(c), (abs(c)));
    hold on;
    semilogy(1:length(ct), (abs(ct)));
    hold off;

%%
%}

    % --- Plot the set Λ(f) ---
    %{
        wraparound_OpenLimit = [openLimit_sorted, openLimit_sorted(1)];
        figure;
        plot(real(wraparound_OpenLimit), imag(wraparound_OpenLimit), 'k-', 'LineWidth', 2.5)
        hold on;
        % ---- unit circle ----
        theta = linspace(0, 2*pi, 300);
        plot(cos(theta), sin(theta), 'r-', 'LineWidth', 2);
        plot(real(openLimit_sorted), imag(openLimit_sorted), 'bx', 'LineWidth', 2.5)
        %legend({'$\Lambda(f)$', 'Unit torus', ''}, 'Interpreter', 'latex', 'FontSize', 18);

        set(gcf, 'Position', [100, 100, 300, 300]);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
        xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 18);
        ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 18); 
        grid on;
        axis equal;
        box on;
        hold off;
    %}


% ==== Compute Fourier coefficients of f(p(z)) numerically  ====
    % --- Evaluate f(p(z)) on the torus ---
    k_values = -m:m;
    powers_matrix = openLimit_sorted(:).^(-k_values);  % N x (2n+1) matrix
    
    % --- Vectorized sum ---
    fp_values = powers_matrix * a(:);

    % --- Clean up data ---
    fp_values = real(fp_values);

    % --- Wrap around ---
    phase_sorted = [phase_sorted; pi];
    fp_values    = [fp_values; fp_values(1)];

    % --- Plot the function f(p(torus))
    %{
    figure;
    plot(phase_sorted, fp_values);
    hold off;
    %}

    % --- Compute the Fourier Transform of f(p(z)) ---
    F_range = DimT;
    FourierFP = fourier_coefficients_spectral(phase_sorted, fp_values, F_range);

    % --- Plot the decay in the Fourier Coefficients ---
    %{
    figure;
    semilogy(1:length(FourierFP), (abs(FourierFP)));
    hold off;
    %}

    
% ==== Quasi Similarity Transformed Toeplitz matrix ====
    % --- Toeplitz matrix for deformed path ---
    ba = FourierFP;

    %ba = [1,0, 0, 0, 1];
    T_b = fourier_to_toeplitz(ba, DimT);
    eigT_b = sort(eig(T_b));


    ero(num_lambda) =  hausdorffDistance(eigT, eigT_b);

    FourCo(num_lambda) = ba(1);

    % --- Wasserstein distance ---
    lam  = sort(real(eigT(:)));
    lamb = sort(real(eigT_b(:)));

    W1(num_lambda) = sum(abs(lam - lamb));

% ==== Hirschman density of states ====
    
    lambda_interval = linspace(lambda_interval(1), lambda_interval(2), 500);
    deriv_sums = zeros(size(lambda_interval));
    
    for i = 1:length(lambda_interval)
        deriv_sums(i) = 1/(2*pi) * 1/g(lambda_interval(i), a) * compute_g_derivative_sum(lambda_interval(i), a);
        deriv_sums(i) = min(20, deriv_sums(i));
    end


end


%% --- Plot the Hausdorff distance ---

close all;
nvals = 1:itera;
rate = 2;
trend = 1./nvals.^rate;
% Plot
    figure;
    loglog(1:itera, ero, 'k-', 'LineWidth', 1.5);
    hold on;
    %loglog(1:itera, W1, 'r', 'LineWidth', 1.5);
    %loglog(1:itera, abs(FourCo-FourCo(end)), 'g-', 'LineWidth', 1.5);

    loglog(1:itera, trend, 'r--', 'LineWidth', 1.5);
%loglog(kt, ero_fit, 'r--', 'LineWidth', 2);  % trendline
    xlabel('Sampling size $N$', 'Interpreter', 'latex', 'FontSize', 14);
    %ylabel('$\varepsilon_N$', 'Interpreter', 'latex', 'FontSize', 14); 
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 250]); 
    legend('$d_{\sigma}(\mathbf{T}_{50}(f), \mathbf{T}_{50}(f\circ p))$', '$N^{-2}$', 'Interpreter', 'latex', 'FontSize', 18);
    %ylim([1e-8, 1e-3]);
    
    grid on;
    hold off;

    
%%

% --- Plot Wasserstein Distance ---


    % --- Plot the eigenvalues before and after asymptotic Similarity transform ---
    %
    figure;
    plot(real(eigT), imag(eigT), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold on;
    plot(real(eigT_b), imag(eigT_b), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    grid on;
    box on;
    xlim([0.96*min(real(eigT)), 1.03*max(real(eigT))])
    ylim([-0.005, 0.005])
    legend({'$\sigma(\mathbf{T}_n(f))$', '$\sigma(\mathbf{T}_n(f\circ p))$'}, 'Interpreter', 'latex', 'FontSize', 18);
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14); 
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 
    axis equal;
    hold off;
    %}


        % --- Plot the DOS against the empirical measure ---
    %
    figure;
    histogram(eigT_b, floor(DimT/2), 'Normalization', 'pdf');
    hold on;
    histogram(eigT, floor(DimT/2), 'FaceColor', 'g', 'Normalization', 'pdf');
    %plot(lambda_interval, deriv_sums, 'LineWidth', 2);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 600, 250]); 
    xlabel('Eigenvalue', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Density of States', 'Interpreter', 'latex', 'FontSize', 14); 
    box on;
    hold off;
    %}
%% --- Defining functions ---


function ck = fourier_coefficients_spectral(phase, fp_values, K)
    
    N = length(phase);
    M = max(2*K+1, 2*N); 
    theta_uniform = (0:M-1)' * 2*pi/M;
    
    % --- Get uniform sampling ---
    fp_uniform = trig_barycentric_interp(phase, fp_values, theta_uniform);
    
    % --- Use standard (uniform) FFT ---
    fft_result = fft(fp_uniform) / M;
    
    % --- Extract coefficients ---
    ck = zeros(2*K+1, 1);
    ck(K+1) = fft_result(1);
    ck(K+2:end) = fft_result(2:K+1);
    ck(1:K) = fft_result(end-K+1:end);
end


function f_interp = trig_barycentric_interp(x, f, x_target)

    N = length(x);

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


function deriv_sum = compute_g_derivative_sum(lambda_real, a, h)
    % Compute ∂g/∂n₁ + ∂g/∂n₂ for Hirschman DOS
    
    if nargin < 3
        h = 1e-5;
    end
    
    lambda_complex = lambda_real + 0i;
    
    g_plus = g(lambda_complex + 1i*h, a);
    g_center = g(lambda_complex, a);
    g_minus = g(lambda_complex - 1i*h, a);
    
    deriv_sum = (g_plus - 2*g_center + g_minus) / h;
end


function dH = hausdorffDistance(A, B)
%HAUSDORFFDISTANCE Compute Hausdorff distance between two sets of complex numbers
%
%   dH = hausdorffDistance(A, B)
%
%   A, B : vectors of real or complex numbers
%   dH   : Hausdorff distance

    % Ensure column vectors
    A = A(:);
    B = B(:);

    % Pairwise distances |a - b|
    D = abs(A - B.');

    % Directed distances
    dAB = max(min(D, [], 2));  % from A to B
    dBA = max(min(D, [], 1));  % from B to A

    % Hausdorff distance
    dH = max(dAB, dBA);
end
