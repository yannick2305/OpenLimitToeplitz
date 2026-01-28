%{
    --------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [November 2025]
    Description:  [Illustrate Open spectrum and symbol collapse]
    --------------------------------------------------------------
%}

clear all;
close all;
clc;

% --- Parameters ---
    n = 6;          % Truncation size for a_k
    p = 1;          % Upwards decay rate
    q = 6;          % Downwards decay rate
    fs = 30;        % Fontsize for plot annotation
    DimT = 80;      % Dimension of finite Toeplitz matrix to simulate open limit
    gif_filename = 'Open_limit_complex.gif';
    delay_time = 0.05;  % Delay between frames in seconds 

% --- Generate the Toeplitz operator ---
    col = zeros(n,1);
    row = zeros(1,n);
    
    col(1) = 1;
    row(1) = 1;
    
    for k = 2:n
        col(k) = 1 / ((k)^q); % sub-diagonals
        row(k) = 1 / ((k)^p); % super-diagonals
    end

    col(n+1 : DimT) = 0;
    col_band = col(1 : n);
    
    row(n+1 : DimT) = 0;
    row_band = row(1 : n);
    
    a = [ col_band(end:-1:1)', row_band(2:end) ];
    
    % --- Generate finite Toeplitz matrix ---
    T = toeplitz(col, row);
    eigT = sort(eig(T));

 % --- Generate a movie for symbol collapse ---
    % --- Generate symbol function ---
    n = n-1;
    syms z
    Qz = sym(0);
    for j = -n:n
        Qz = Qz + a(j + n + 1) * z^(n - j);
    end
    
    Qcoeffs = sym2poly(Qz);

    % --- Define the unit circle ---
    N = 1000;
    theta = linspace(0, 2*pi, N);
    z_torus = exp(1i * theta);

    % --- Range of r values ---
    r_values = linspace(3.5, 5, 500);

    % --- Prepare figure ---
    fig = figure('Color','w', 'Position', [100, 100, 800, 600]);
    ylim([-max(1.3*max(imag(eigT)), 0.1), max(1.3*max(imag(eigT)), 0.1)]);
    grid on
    xlabel('$\Re$', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\Im$', 'Interpreter', 'latex', 'FontSize', fs);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs+2);

    hold on;

    % --- Generate approximate open limit ---

    x = real(eigT);
    y = imag(eigT);
    
    % Sort points: real axis first (by real part), then complex points (by distance from real axis, then real part)
    realAxisMask = abs(y) < 1e-6;  
    realAxisIdx = find(realAxisMask);
    complexIdx = find(~realAxisMask);
    
    % Sort real axis points by real part
    [~, sortReal] = sort(x(realAxisIdx));
    realAxisIdx = realAxisIdx(sortReal);
    
    % Sort complex points by real part (to group branches)
    [~, sortComplex] = sort(x(complexIdx));
    complexIdx = complexIdx(sortComplex);
    
    % Plot real axis segment
    if ~isempty(realAxisIdx)
        eig_plot = plot(x(realAxisIdx), y(realAxisIdx), 'b-', 'LineWidth', 1.5);
    else
        eig_plot = [];
    end
    
    % Plot branches from real axis into complex plane
    threshold = 0.15;  % Adjust based on spacing
    first_branch = isempty(eig_plot);
    
    for i = 1:length(complexIdx)
        idx = complexIdx(i);
        
        % Find nearest point on real axis or already plotted complex point
        distances = sqrt((x - x(idx)).^2 + (y - y(idx)).^2);
        distances(idx) = inf;
        
        % Prefer connecting to real axis or points with smaller imaginary part
        distances(~realAxisMask & abs(y) > abs(y(idx))) = inf;
        
        [minDist, nearestIdx] = min(distances);
        
        if minDist < threshold
            if first_branch
                eig_plot = plot([x(nearestIdx), x(idx)], [y(nearestIdx), y(idx)], 'b-', 'LineWidth', 2.5);
                first_branch = false;
            else
                plot([x(nearestIdx), x(idx)], [y(nearestIdx), y(idx)], 'b-', 'LineWidth', 2.5, 'HandleVisibility', 'off');
            end
        end
    end
   
    % --- Plot the symbol function on the scaled torus ---
    curve_plot = plot(nan, nan, 'k', 'LineWidth', 2);
    xlim([0.8, 1.3])
    ylim([-0.1, 0.1])
    box on;
    title(sprintf('r = %.3f', r_values(1)), 'FontSize', fs+7, 'Interpreter', 'latex');

    legend('$\lim_{n\to\infty} \sigma(T_n(f))$', '$f(r\mathbf{T})$', 'Interpreter', ...
        'latex', 'FontSize', 25., 'Location','northeast');

    % --- Animation loop ---
    for idx = 1:length(r_values)
        r = r_values(idx);
        z_torus = r * exp(1i * theta);
        Qz_vals = zeros(1, N);
    
        for j = 1:N
            z = z_torus(j);
            Qz = 0;
            for k = -n:n
                Qz = Qz + a(k + n + 1) * z^(n - k);
            end
            Qz_vals(j) = Qz / z^n;
        end
    
        % --- Update plot ---
        set(curve_plot, 'XData', real(Qz_vals), 'YData', imag(Qz_vals));
        title(sprintf('r = %.3f', r), 'FontSize', fs, 'Interpreter', 'latex');
        drawnow;
    
        % --- Capture frame and write to GIF ---
        frame = getframe(fig);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        
        if idx == 1
            imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
        else
            imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        end
    end
    
    close(fig);
    fprintf('GIF saved as: %s\n', gif_filename);   
