%{
    ----------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2025]
    Description:  [Polluted spectra]
    ----------------------------------------------------------------------
%}


clc
clear;
close all;

% ==== Parameters ====
    m = 9;              % Truncation size for a_k
    p = 3.5;            % Decay rate upwards
    q = 4.8;            % Decay rate downwards
    DimT = 2000;        % Dimension of finite Toeplitz matrix to simulate open limit
    fs = 18;            % Fontsize for annotation
    rad = 1.0;          % Radius of torus
    

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

    % --- Evaluate symbol on scaled unit torus ---
    theta = linspace(0, 2*pi, 300);
    r_torus = rad * exp(1i*theta);
    Laurent =  evaluate_laurent_poly(r_torus, a, m);

    % --- Plot the eigenvalues and the image of symbol function ---
    figure;
    plot(real(eigT), imag(eigT), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold on;
    plot(real(Laurent), imag(Laurent), 'r-', 'MarkerSize', 8, 'LineWidth', 3.5);
    grid on;
    box on;
    xlim([0.96*min(real(eigT)), 1.03*max(real(eigT))])
    ylim([-0.065, 0.1])
    xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 14); 
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    legend({'$\sigma(T(f))$', '$f(re^{i\theta})$'}, 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'northwest');
    set(gcf, 'Position', [100, 100, 500, 250]); 
    hold off;


%% --- Defining functions ---


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
