%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2025]
    Description:  [Reality of open limit]
    --------------------------------------------------------------
%}

clear all;
close all;

% --- Parameters ---
    n = 5;              % Truncation size for a_k
    p = 3.5;            % Decay rate upwards
    q = 4.8;            % Decay rate downwards
    DimT = 80;          % Dimension of finite Toeplitz matrix to simulate open limit
    num_lambda = 2000;  % Number of plotting points (keep it large)
    threshold = 0.8;    % Maximum distance to nearest point
    fs = 18;            % Fontsize for annotation
    
% --- Generate sequence a_k ---
    col = zeros(n,1);
    row = zeros(1,n+1);
    
    % --- Add noise to the coeffiients ---
    ai = 1;
    bi = 1;

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
    col = [1, col'];
    col = col';

    col(n+1 : DimT) = 0;
    col_band = col(1 : n);

    row(n+1 : DimT) = 0;
    row_band = row(1 : n);
    
    T = toeplitz(col, row);
    eigT = sort(eig(T));

    % --- Define interval for the CBS ---
    %lambda_range = [0.992*min(real(eigT)), 1.008*max(real(eigT)) ];
    lambda_range = [0.8, 1.3];

    % --- Plot the open spectrum ---
    figure;
    plot(real(eigT), imag(eigT), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold on;
    grid on; 
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 250]); 
    ylim([-max(1.2*max(imag(eigT)), 0.1), max(1.2*max(imag(eigT)), 0.1)]);
    hold off;


%% --- Generate the Complex band structure and Λ(f) path ---

    m = (length(a) - 1) / 2;

    % --- Generate the symbol function ---
    syms z
    Qz = sym(0);
    for j = -m:m
        Qz = Qz + a(j + m + 1) * z^(m - j);
    end

    % --- Generate frequency range ---
    lambda_vals = linspace(lambda_range(1), lambda_range(2), num_lambda);

    % --- Initiallise the CBS ---
    alpha_all  = [];
    beta_all   = [];
    lambda_all = [];

    openLimitRoots = [];

    % --- Loop through lambda values and collect roots ---
    for k = 1:length(lambda_vals)
        lambda = lambda_vals(k);
        Pz = Qz - lambda * z^m;
        roots_z = roots(sym2poly(Pz));
        
        % --- Sort roots by magnitude ---
        [~, sort_idx] = sort(abs(roots_z));
        roots_z_sorted = roots_z(sort_idx);
        
        % --- Keep the m-th and (m+1)-th roots ---
        candidate_roots = roots_z_sorted(m:m+1);
        
        % --- Check |z_m| = |z_{m+1}| ---
        tol = 1e-3;
        mag1 = abs(candidate_roots(1));
        mag2 = abs(candidate_roots(2));

        if abs(mag1 - mag2) < tol
            selected_roots = candidate_roots;
            openLimitRoots = [openLimitRoots, candidate_roots(1), candidate_roots(2)];
        end

        roots_z = roots_z(~(abs(roots_z) < 1e-12));
        
        % --- Compute the complex band structure ---
        alpha = angle(roots_z);
        beta = -log(abs(roots_z));
        
        alpha_all  = [alpha_all; alpha(:)];
        beta_all   = [beta_all; beta(:)];
        lambda_all = [lambda_all; lambda * ones(length(alpha), 1)];
    end
       

    % --- Check for confluent roots ---
        % --- Compute derivative of symbol function for confluent roots ---
        Qz_derivative = diff(Qz*z^(-m), z);
            
        % --- Confluent roots z satisfy f'(z) = 0 ---
        roots_Qz_derivative = double(solve(Qz_derivative, z)); 
    
        tol = 0.5; 
        confluentRoots = [];
            
        for i = 1:length(openLimitRoots)
            for j = 1:length(roots_Qz_derivative)
                if abs(openLimitRoots(i) - roots_Qz_derivative(j)) < tol
                    confluentRoots = [confluentRoots; roots_Qz_derivative(j)];
                    break;  % Move to next root in openLimitRoots
                end
            end
        end

    % --- Plot the set Λ(f) by Simple nearest-neighbor path ---
    x = real(openLimitRoots);
    y = imag(openLimitRoots);

    Limy  = max(y);
    Limxm = min(x);
    LimxM = max(x);

    figure;
    hold on;
    
    h_curves = [];  % Store handle for legend

    unvisited = true(size(x));
    while any(unvisited)
        % --- Start a new curve from an unvisited point --- 
        idx = find(unvisited, 1);
        currentCurve = idx;
        unvisited(idx) = false;
        
        % --- Build curve by adding nearest neighbors ---
        while true
            % --- Find nearest unvisited point to the last point in curve ---
            lastPt = currentCurve(end);
            distances = sqrt((x - x(lastPt)).^2 + (y - y(lastPt)).^2);
            distances(~unvisited) = inf;
            
            [minDist, nearestIdx] = min(distances);
            
            % --- Stop if no more nearby points ---
            if isinf(minDist)
                break;
            end

            if minDist > threshold
                break;
            end
            
            % --- Save the path ---
            currentCurve = [currentCurve; nearestIdx];
            unvisited(nearestIdx) = false;
        end
        
        % --- Plot the current branch ---
        h = plot(x(currentCurve), y(currentCurve), 'k-', 'LineWidth', 2.5);

        if isempty(h_curves)
            h_curves = h;  % handle for legend
        end
    end
    
    h_roots = plot(real(confluentRoots), imag(confluentRoots), 'rx', 'LineWidth', 3.0);

    xlabel('$\Re$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\Im$', 'Interpreter', 'latex', 'FontSize', 14);
    legend([h_curves, h_roots], '$\Lambda(f_m)$', 'Confluent roots', 'Interpreter', 'latex', 'FontSize', 16);
    set(gca, 'TickLabelInterpreter', 'latex');
    set(gca, 'FontSize', 18);
    ylim([-1.5 * Limy,  2 * Limy]);
    xlim([ 1.3 * Limxm, 1.3 * LimxM]);
    set(gcf, 'Position', [100, 100, 500, 300]);
    box on;
    hold off;



%%  --- Plot the corresponding complex band structure (Scatter) ---

    figure; 
    hold on;
    scatter(alpha_all, lambda_all, 10, 'k', 'filled'); 
    scatter(beta_all,  lambda_all, 10, 'r', 'filled'); 

    % Plot formatting
    xlabel('$\alpha$ and $\beta$ respectively', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', fs);
    
    set(gcf, 'Position', [100, 100, 500, 300]); 
    xlim([-pi - 0.01, pi + 0.01]);
    ylim([lambda_range(1), lambda_range(2)]);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs+2);
    box on;
    grid on;
