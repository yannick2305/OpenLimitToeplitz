%{
    ----------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2025]
    Description:  [Gamma curve crossing]
    ----------------------------------------------------------------------
%}

clc
clear;
close all;

% ==== Parameters ====
    m = 2;              % Truncation size for a_k
    p = 3.5;            % Decay rate upwards
    q = 4.8;            % Decay rate downwards
    DimT = 80;         % Dimension of finite Toeplitz matrix to simulate open limit
    num_lambda = 1000;   % Number of plotting points (50-300)
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
    %a = [ col(end:-1:1)', row];
    a = [0.6, 2.5, 3.2, 1.6, 0.0];
       

    %a = [0.6, 2.5, 4, 1.6, 0.0];
       
    % --- Generate finite Toeplitz matrix ---
    T = fourier_to_toeplitz(a, DimT);
    eigT = sort(eig(T));


    % Compute eigenvalues and eigenvectors
[V, D] = eig(T);

% --- Plot the eigenvectors ---
%{
    figure;
    hold on;
    grid on;
    
    % Small positive floor to avoid log(0)
    tol = eps;
    
    for k = 1:DimT
        vk = abs(V(:,k));
        vk(vk < tol) = tol;   % enforce positivity for log scale
        semilogy(1:length(vk), vk, 'k', 'LineWidth', 1.5);
    end
    
    set(gcf, 'Position', [100, 100, 500, 300]);
    
    set(gca, 'TickLabelInterpreter', 'latex', ...
             'FontSize', 18, ...
             'YScale', 'log');   % force log scale (important)
    
    xlabel('Index', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('$|\mathbf{v}_i|$', 'Interpreter', 'latex', 'FontSize', 18);
    
    box on;
    hold off;

 
%}

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


    % --- Plot the set Λ(f) ---
    %
% Suppose these are your points
A = 0.6827      + 0.126i;
B = 0.1211  + 0.7682i;
C = -0.749946 - 0.661467i;
%D = -1.3324   - 0.0159i;
D = -1.2922 + 0.1246i;

% Plot the curve and unit circle
wraparound_OpenLimit = [openLimit_sorted, openLimit_sorted(1)];
figure;
plot(real(wraparound_OpenLimit), imag(wraparound_OpenLimit), 'k-', 'LineWidth', 2.5)
hold on;

theta = linspace(0, 2*pi, 300);
plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1);

% ---- Add points A, B, C, D ----
plot(real(A), imag(A), 'rx', 'MarkerSize', 10, 'LineWidth', 2); % red circle
plot(real(B), imag(B), 'bs', 'MarkerSize', 10, 'LineWidth', 2); % blue square
plot(real(C), imag(C), 'g^', 'MarkerSize', 10, 'LineWidth', 2); % green triangle
plot(real(D), imag(D), 'md', 'MarkerSize', 10, 'LineWidth', 2); % magenta diamond

% ---- Add LaTeX labels next to points ----
text(real(A)+0.02, imag(A),      '$A$', 'FontSize', 18, 'Color', 'r', 'Interpreter', 'latex');
text(real(B)+0.02, imag(B)-0.18, '$B$', 'FontSize', 18, 'Color', 'b', 'Interpreter', 'latex');
text(real(C)+0.02, imag(C)+0.07, '$C$', 'FontSize', 18, 'Color', 'g', 'Interpreter', 'latex');
text(real(D)+0.05, imag(D),      '$D$', 'FontSize', 18, 'Color', 'm', 'Interpreter', 'latex');

% Formatting the figure
set(gcf, 'Position', [100, 100, 300, 300]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 18); 
grid on;
axis equal;
box on;
hold off;


    %}

%
% ==== Compute Fourier coefficients of f(p(z)) numerically  ====
    % --- Evaluate f(p(z)) on the torus ---
    k_values = -m:m;
    powers_matrix = openLimit_sorted(:).^(-k_values);  % N x (2n+1) matrix
    
    % --- Vectorized sum ---
    fp_values = powers_matrix * a(:);

    % --- Clean up data ---
    fp_values = real(fp_values);

    % --- Wrap around ---
    phase_sorted = [-pi; phase_sorted; pi];
    fp_values    = [ fp_values(end); fp_values; fp_values(end)];


points = [A, B, C, D];      % array of points
labels = {'A','B','C','D'}; % corresponding labels

% --- Compute phase (angle) of each point ---
point_phases = angle(points);  % returns values in [-pi, pi]

% --- phase_sorted and fp_values are given ---
% phase_sorted: array of phases where fp_values are evaluated
% fp_values: corresponding fp values

% --- Initialize array to store fp values corresponding to points ---
fp_points = zeros(size(points));

% --- Loop through each point, find closest phase in phase_sorted ---
for k = 1:length(points)
    [~, idx] = min(abs(phase_sorted - point_phases(k))); % closest index
    fp_points(k) = fp_values(idx);                       % pick fp value
end

% --- Optional: display results ---
for k = 1:length(points)
    fprintf('%s: phase = %.3f rad, closest fp = %.3f\n', labels{k}, point_phases(k), fp_points(k));
end

% Eigen-decomposition
[V, D] = eig(T);
eigvals = diag(D);      % eigenvalues
Neig = length(eigvals);

% Storage for matched eigenvectors
eigvecs_points = zeros(size(V,1), length(fp_points));
eigval_matches = zeros(size(fp_points));

for k = 1:length(fp_points)
    [~, idx] = min(abs(eigvals - fp_points(k)));  % closest eigenvalue
    eigvecs_points(:,k) = V(:,idx);               % eigenvector
    eigval_matches(k) = eigvals(idx);             % store eigenvalue
end

colors = {'r','b','g','m'};
labels = {'A','B','C','D'};

markers = {'x','s','^','d'};   % same order as A,B,C,D

figure;
hold on;

for k = 1:length(fp_points)
    v = eigvecs_points(:,k);

    % Normalize for visual comparison
    v = v / norm(v);

    semilogy(1:length(v), abs(v) + eps, ...
        markers{k}, ...                 % <-- marker style
        'Color', colors{k}, ...
        'LineWidth', 1.5, ...
        'MarkerSize', 6);
end

set(gca, 'TickLabelInterpreter', 'latex', ...
         'FontSize', 18, ...
         'YScale', 'log');

xlabel('Index', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$|\mathbf{v}_i|$', 'Interpreter', 'latex', 'FontSize', 18);

set(gcf, 'Position', [100, 100, 500, 300]);
grid on;
box on;
hold off;



%% --- Defining functions ---


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

