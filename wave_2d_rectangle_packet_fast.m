% wave2d_rectangle_packet_fast.m
%
% 2D wave equation on a rectangle with fixed walls:
%     u_tt = c^2 (u_xx + u_yy)
%     u = 0 on the boundary
%
% Solved by double sine-series expansion.
% Initial condition: traveling Gaussian-modulated carrier wave packet.
%
% Uses exact modal propagator and faster graphics updates.

clear;
clc;
close all;

%% ========================================================================
%  OPTIONS
%  ========================================================================

opts.show_diagnostics_during_animation = false;   % true = slower
opts.show_final_amplitude_figure       = false;
opts.show_final_energy_figures         = false;
opts.print_summary                     = true;

opts.colormap_name = 'turbo';
opts.figure_name   = '2D Wave Packet on a Rectangle';

% Performance options
opts.visual_frame_stride = 1;    % 1 = draw every frame, 2 = every other frame, etc.
opts.use_drawnow_limitrate = true;
opts.reuse_graphics_handles = true;

%% ========================================================================
%  PHYSICAL / DOMAIN PARAMETERS
%  ========================================================================

Lx = 2.0;
Ly = 4.0;
c  = 2.0;

%% ========================================================================
%  MODE TRUNCATION
%  ========================================================================

M = 60;
N = 60;

%% ========================================================================
%  PROJECTION GRID
%  ========================================================================

Nxq = 200;
Nyq = 200;

xq = linspace(0, Lx, Nxq);
yq = linspace(0, Ly, Nyq);

dx = xq(2) - xq(1);
dy = yq(2) - yq(1);

[Xq, Yq] = meshgrid(xq, yq);

%% ========================================================================
%  INITIAL WAVE PACKET
%  ========================================================================

A0 = 1.0;
x0 = 0.28 * Lx;
y0 = 0.35 * Ly;

sigx = 0.08 * Lx;
sigy = 0.08 * Ly;

theta   = pi/6;
lambda0 = 0.12;
k0      = 2*pi/lambda0;

kx0 = k0 * cos(theta);
ky0 = k0 * sin(theta);

dirx = cos(theta);
diry = sin(theta);

Xc = Xq - x0;
Yc = Yq - y0;

env   = A0 * exp(-(Xc.^2)/(2*sigx^2) - (Yc.^2)/(2*sigy^2));
phase = kx0 * Xc + ky0 * Yc;

u0 = env .* cos(phase);

env_x = -(Xc / sigx^2) .* env;
env_y = -(Yc / sigy^2) .* env;

u0_x = env_x .* cos(phase) - env .* sin(phase) * kx0;
u0_y = env_y .* cos(phase) - env .* sin(phase) * ky0;

v0 = -c * (dirx * u0_x + diry * u0_y);

%% ========================================================================
%  SINE BASIS ON PROJECTION GRID
%  ========================================================================

mx = 1:M;
ny = 1:N;

Sxq = sin(xq(:) * (mx * pi / Lx));   % Nxq x M
Syq = sin(yq(:) * (ny * pi / Ly));   % Nyq x N

%% ========================================================================
%  TRAPEZOIDAL WEIGHTS
%  ========================================================================

wx = ones(Nxq, 1);
wy = ones(Nyq, 1);

wx(1)   = 0.5;
wx(end) = 0.5;
wy(1)   = 0.5;
wy(end) = 0.5;

W = wy * wx.';

U0w = u0 .* W;
V0w = v0 .* W;

%% ========================================================================
%  PROJECT TO MODAL SPACE
%  ========================================================================

Q = (4/(Lx*Ly)) * dx * dy * (Sxq.' * U0w.' * Syq);
P = (4/(Lx*Ly)) * dx * dy * (Sxq.' * V0w.' * Syq);

%% ========================================================================
%  MODAL FREQUENCIES
%  ========================================================================

[MM, NN] = ndgrid(mx, ny);
omega = c * pi * sqrt((MM/Lx).^2 + (NN/Ly).^2);

%% ========================================================================
%  DISPLAY GRID
%  ========================================================================

Nx = 300;   % lower than 400 for faster animation
Ny = 600;   % lower than 800 for faster animation

x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);

Sx = sin(x(:) * (mx * pi / Lx));   % Nx x M
Sy = sin(y(:) * (ny * pi / Ly));   % Ny x N

%% ========================================================================
%  TIME SETTINGS
%  ========================================================================

Tfinal  = 8.0;
nFrames = 120 * Tfinal;

dt    = Tfinal / (nFrames - 1);
tlist = (0:nFrames-1) * dt;

%% ========================================================================
%  EXACT ONE-STEP PROPAGATOR
%  ========================================================================

Cd        = cos(omega * dt);
Sd        = sin(omega * dt);
Sd_over_w = Sd ./ omega;
wSd       = omega .* Sd;

%% ========================================================================
%  FIXED COLOR LIMITS
%  ========================================================================

umax0 = max(abs(u0(:)));
if umax0 == 0
    umax0 = 1;
end
clim = 1.1 * umax0;

%% ========================================================================
%  DIAGNOSTICS STORAGE
%  ========================================================================

energy    = zeros(nFrames, 1);
rel_drift = zeros(nFrames, 1);
umax_phys = zeros(nFrames, 1);

Ekin = zeros(nFrames, 1);
Epot = zeros(nFrames, 1);

E0 = 0.5 * sum(sum(P.^2 + (omega.^2) .* (Q.^2)));
if E0 == 0
    E0 = 1;
end

%% ========================================================================
%  FIGURE / GRAPHICS SETUP
%  ========================================================================

fig = figure('Name', opts.figure_name);

try
    set(fig, 'Color', 'w');
catch
end

colormap(opts.colormap_name);

if opts.show_diagnostics_during_animation
    ax_field = subplot(2,2,[1 3]);
else
    ax_field = axes('Parent', fig);
end

% Create image once
U_init = Sx * Q * Sy.';
Uplot_init = U_init.';

hImg = imagesc(ax_field, x, y, Uplot_init);
axis(ax_field, 'xy');
axis(ax_field, 'equal');
axis(ax_field, 'tight');
caxis(ax_field, [-clim, clim]);
colorbar(ax_field);
xlabel(ax_field, 'x');
ylabel(ax_field, 'y');
hTitle = title(ax_field, sprintf('2D wave packet,  t = %.4f', 0));

if opts.show_diagnostics_during_animation
    ax_energy = subplot(2,2,2);
    hEnergy = plot(ax_energy, NaN, NaN, 'LineWidth', 1.5);
    grid(ax_energy, 'on');
    xlim(ax_energy, [tlist(1), tlist(end)]);
    xlabel(ax_energy, 't');
    ylabel(ax_energy, 'E(t)');
    title(ax_energy, 'Total modal energy');

    ax_drift = subplot(2,2,4);
    hDrift = plot(ax_drift, NaN, NaN, 'LineWidth', 1.5);
    grid(ax_drift, 'on');
    xlim(ax_drift, [tlist(1), tlist(end)]);
    xlabel(ax_drift, 't');
    ylabel(ax_drift, '(E(t)-E_0)/E_0');
    title(ax_drift, 'Relative energy drift');
end

%% ========================================================================
%  ANIMATION LOOP
%  ========================================================================

for it = 1:nFrames
    t = tlist(it);

    % Reconstruct field
    U = Sx * Q * Sy.';
    Uplot = U.';

    % Diagnostics
    Ekin(it)      = 0.5 * sum(sum(P.^2));
    Epot(it)      = 0.5 * sum(sum((omega.^2) .* (Q.^2)));
    energy(it)    = Ekin(it) + Epot(it);
    rel_drift(it) = (energy(it) - E0) / E0;
    umax_phys(it) = max(abs(Uplot(:)));

    % Update graphics only on selected visual frames
    if mod(it-1, opts.visual_frame_stride) == 0 || it == nFrames
        set(hImg, 'CData', Uplot);
        set(hTitle, 'String', sprintf('2D wave packet,  t = %.4f', t));

        if opts.show_diagnostics_during_animation
            set(hEnergy, 'XData', tlist(1:it), 'YData', energy(1:it));
            set(hDrift,  'XData', tlist(1:it), 'YData', rel_drift(1:it));
        end

        if opts.use_drawnow_limitrate
            try
                drawnow limitrate;
            catch
                drawnow;
            end
        else
            drawnow;
        end
    end

    % Advance one exact propagator step
    if it < nFrames
        Qnew = Cd .* Q + Sd_over_w .* P;
        Pnew = -wSd .* Q + Cd .* P;

        Q = Qnew;
        P = Pnew;
    end
end

%% ========================================================================
%  TEXT SUMMARY
%  ========================================================================

if opts.print_summary
    fprintf('Initial energy           : %.16e\n', energy(1));
    fprintf('Final energy             : %.16e\n', energy(end));
    fprintf('Max abs relative drift   : %.16e\n', max(abs(rel_drift)));
    fprintf('Max physical amplitude   : %.16e\n', max(umax_phys));
end

%% ========================================================================
%  OPTIONAL FINAL FIGURES
%  ========================================================================

if opts.show_final_energy_figures
    figure('Name', 'Energy diagnostics');
    plot(tlist, Ekin, 'LineWidth', 1.5);
    hold on;
    plot(tlist, Epot, 'LineWidth', 1.5);
    plot(tlist, energy, '--', 'LineWidth', 1.5);
    hold off;
    grid on;
    xlabel('t');
    ylabel('Energy');
    legend('Kinetic', 'Potential', 'Total');
    title('Energy partition');
end

if opts.show_final_amplitude_figure
    figure('Name', 'Amplitude diagnostic');
    plot(tlist, umax_phys, 'LineWidth', 1.5);
    grid on;
    xlabel('t');
    ylabel('max |u(x,y,t)|');
    title('Maximum physical-space amplitude');
end
