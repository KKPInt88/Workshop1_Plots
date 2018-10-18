%% Matlab Series Workshop
% Week 1: Data Visualization (Plots)
% Kan Kanjanapas (Ph.D.)

clc;
close all;
clear all;


% 1. Plot

fs = 100;      % Sampling Frequency [Hz]
Ts = 1/fs;      % Sampling Time [s]

t_vec = [0:Ts:5]';  % Vector of time stamps

% sinusoidal waveform: A*sin(omega*t + phase_shift)

% Assume the reading from sensor A is given by:
f_1 = 1;    

x_1 = 1*sin(2*pi*f_1*t_vec + 0);
x_2 = 2*sin(2*pi*f_1*t_vec + 0);
x_3 = 3*sin(2*pi*f_1*t_vec + 0);
x_4 = 4*sin(2*pi*f_1*t_vec + 0);
x_5 = 5*sin(2*pi*f_1*t_vec + 0);

x_ceil = [];
for ii = 1:5
    x_ceil{ii} = ii*sin(2*pi*f_1*t_vec + 0);
end



%% 1) Plot

% Your first plot



% Refine the first plot 1.0



% Refine the first plot: version 1.1 --------------------------------------------------------------------------------------------------




% Refine the first plot: version 1.2 --------------------------------------------------------------------------------------------------




% Refine the first plot: version 1.3 --------------------------------------------------------------------------------------------------





% Refine the first plot: version 1.4 --------------------------------------------------------------------------------------------------




%% 2) Subplot

% Subplot Version 2.1 ----------------------------------------------------------------------------------------------------------




% Subplot Version 2.2 ----------------------------------------------------------------------------------------------------------




% Subplot Version 2.3 ----------------------------------------------------------------------------------------------------------
% % % figure;
% % % set(gcf, 'Position', [0 0 2560 1280]/2);
% % % 
% % % 
% % % subplot(2,2,1);
% % % plot(t_vec, x_1, 'LineStyle', '-', 'LineWidth', 2, 'Color', [0.0  0.0  1.0], ...
% % %      'Marker', 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'y');
% % % xlabel('Time [s]');
% % % ylabel('Position [m]');
% % % title(sprintf('Position Measurement from Sensor A with sampling frequency of %.1f [Hz]', f_1));
% % % h_legend = legend('Amplitude = 1 [m]');
% % % set(h_legend, 'Location', 'NorthEast', 'Color', [1.0  1.0  0.9]);
% % % grid on;
% % % set(gca, 'FontSize', 14); 
% % % 
% % % 
% % % subplot(2,2,3);
% % % plot(t_vec, x_2, 'LineStyle', '-', 'LineWidth', 2, 'Color', [0.0  0.5  0.0], ...
% % %      'Marker', 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g');
% % % xlabel('Time [s]');
% % % ylabel('Position [m]');
% % % h_legend = legend('Amplitude = 2 [m]');
% % % set(h_legend, 'Location', 'NorthEast', 'Color', [1.0  1.0  0.9]);
% % % %title(sprintf('Position Measurement from Sensor A with sampling frequency of %.1f [Hz]', f_1));
% % % grid on;
% % % set(gca, 'FontSize', 14); 
% % % 
% % % subplot(2,2,[2 4]);
% % % plot(t_vec, x_3, 'LineStyle', '-', 'LineWidth', 2, 'Color', Color_Matrix(3,:), ...
% % %      'Marker', 'o', 'MarkerEdgeColor', MarkerEdgeColor_Ceil{3}, 'MarkerFaceColor', 'y');
% % % xlabel('Time [s]');
% % % ylabel('Position [m]');
% % % h_legend = legend('Amplitude = 3 [m]');
% % % set(h_legend, 'Location', 'NorthEast', 'Color', [1.0  1.0  0.9]);
% % % grid on;
% % % set(gca, 'FontSize', 14); 




%% 3) Plot3


theta  = [0:pi/100:20*pi]';
t_vec2 = [0:1:length(theta)-1]'*Ts;
lambda = 0.2;

r = exp(-lambda*t_vec2);
x = r.*cos(theta);
y = r.*sin(theta);
z = theta/4;

% 3.1: Plot3 ----------------------------------------------------------------------------------------------------------



% 3.2: Plot3 Multiple View ----------------------------------------------------------------------------------------------------------



%% 4) Surface

[X,Y] = meshgrid((-5:0.1:5)*pi, (-1:0.1:1)*pi);  % check size(X), size(Y)
Z = X.^2 + Y.^4;  % check size(Z)


% 4.1) Surf ----------------------------------------------------------------------------------------------------------



% 4.2) Surf (Cont) ----------------------------------------------------------------------------------------------------------




%% 5) Contour (contour3, countourc, contourf, quiver)

% [X,Y] = meshgrid((0:0.1:5)*pi, (0:0.1:1)*pi);  % check size(X), size(Y)
% Z = X.^2 + Y.^4;  % check size(Z)



%% 6) Semilogx, Semilogy, loglog

% % % X = logspace(1,3,100);
% % % Y = X.^2;
% % % 
% % % figure;
% % % set(gcf, 'Position', [0 0 2560 1280]/2);
% % % 
% % % subplot(1,3,1);
% % % semilogx(X,Y, 'LineWidth', 2, 'Color', 'b');
% % % xlabel('X [-]');
% % % ylabel('Y [-]');
% % % title('Semilogx');
% % % set(gca, 'FontSize', 14);
% % % grid on;
% % % 
% % % subplot(1,3,2);
% % % semilogy(X,Y, 'LineWidth', 2, 'Color', 'b');
% % % xlabel('X [-]');
% % % ylabel('Y [-]');
% % % title('Semilogy');
% % % set(gca, 'FontSize', 14);
% % % grid on;
% % % 
% % % subplot(1,3,3);
% % % loglog(X,Y, 'LineWidth', 2, 'Color', 'b');
% % % xlabel('X [-]');
% % % ylabel('Y [-]');
% % % title('LogLog');
% % % set(gca, 'FontSize', 14);
% % % grid on;


