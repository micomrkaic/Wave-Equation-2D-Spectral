% wave2d_rectangle_packet_propagator_diag.m
%
% 2D wave equation on a rectangle with fixed walls:
%     u_tt = c^2 (u_xx + u_yy)
%     u = 0 on the boundary
%
% Solved by double sine-series expansion.
% Initial condition: traveling Gaussian-modulated carrier wave packet.
%
% This version advances the modal state with a constant exact propagator
% and records diagnostics:
%   1. total modal energy
%   2. relative energy drift
%   3. max |u(x,y,t)| on the display grid
%
% Written to run in both MATLAB and GNU Octave.

clear;
clc;
close all;

% =========================
% Physical/domain settings
% =========================
Lx = 2.0;
Ly = 4.0;
c  = 2.0;

% =========================
% Mode truncation
% =========================
M = 60;
N = 60;

% =========================
% Quadrature grid for projection
% =========================
Nxq = 200;
Nyq = 200;

xq = linspace(0, Lx, Nxq);
yq = linspace(0, Ly, Nyq);

dx = xq(2) - xq(1);
dy = yq(2) - yq(1);

[Xq, Yq] = meshgrid(xq, yq);

% =========================
% Packet parameters
% =========================
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

% =========================
% Initial displacement u(x,y,0)
% =========================
Xc = Xq - x0;
Yc = Yq - y0;

env   = A0 * exp(-(Xc.^2)/(2*sigx^2) - (Yc.^2)/(2*sigy^2));
phase = kx0 * Xc + ky0 * Yc;

u0 = env .* cos(phase);

% =========================
% Initial velocity ut(x,y,0)
% =========================
env_x = -(Xc / sigx^2) .* env;
env_y = -(Yc / sigy^2) .* env;

u0_x = env_x .* cos(phase) - env .* sin(phase) * kx0;
u0_y = env_y .* cos(phase) - env .* sin(phase) * ky0;

v0 = -c * (dirx * u0_x + diry * u0_y);

% =========================
% Build sine basis on quadrature grid
% =========================
mx = 1:M;
ny = 1:N;

Sxq = sin(xq(:) * (mx * pi / Lx));
Syq = sin(yq(:) * (ny * pi / Ly));

% =========================
% Trapezoidal weights
% =========================
wx = ones(Nxq, 1);
wy = ones(Nyq, 1);

wx(1)   = 0.5;
wx(end) = 0.5;
wy(1)   = 0.5;
wy(end) = 0.5;

W = wy * wx.';

U0w = u0 .* W;
V0w = v0 .* W;

% =========================
% Project initial data onto modes
% =========================
Q = (4/(Lx*Ly)) * dx * dy * (Sxq.' * U0w.' * Syq);
P = (4/(Lx*Ly)) * dx * dy * (Sxq.' * V0w.' * Syq);

% =========================
% Modal frequencies
% =========================
[MM, NN] = ndgrid(mx, ny);
omega = c * pi * sqrt((MM/Lx).^2 + (NN/Ly).^2);

% =========================
% Display grid
% =========================
Nx = 400;
Ny = 800;

x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);

Sx = sin(x(:) * (mx * pi / Lx));
Sy = sin(y(:) * (ny * pi / Ly));

% =========================
% Time settings
% =========================
Tfinal  = 2.0;
nFrames = 240;
dt      = Tfinal / (nFrames - 1);
tlist   = (0:nFrames-1) * dt;

% =========================
% Exact one-step propagator
% =========================
Cd        = cos(omega * dt);
Sd        = sin(omega * dt);
Sd_over_w = Sd ./ omega;
wSd       = omega .* Sd;

% =========================
% Fixed color scale
% =========================
umax0 = max(abs(u0(:)));
if umax0 == 0
    umax0 = 1;
end
clim = 1.1 * umax0;

% =========================
% Diagnostics storage
% =========================
energy      = zeros(nFrames, 1);
rel_drift   = zeros(nFrames, 1);
umax_phys   = zeros(nFrames, 1);

% Initial modal energy:
% E = 0.5 * sum_{m,n} ( P_mn^2 + omega_mn^2 * Q_mn^2 )
E0 = 0.5 * sum(sum(P.^2 + (omega.^2) .* (Q.^2)));
if E0 == 0
    E0 = 1;
end

% =========================
% Figure setup
% =========================
fig = figure('Name', '2D Wave Packet on a Rectangle');

try
    set(fig, 'Color', 'w');
catch
end

colormap(turbo);

% =========================
% Animation loop
% =========================
for it = 1:nFrames
    t = tlist(it);

    % Reconstruct physical-space field
    U = Sx * Q * Sy.';
    Uplot = U.';

    % Diagnostics
    energy(it)    = 0.5 * sum(sum(P.^2 + (omega.^2) .* (Q.^2)));
    rel_drift(it) = (energy(it) - E0) / E0;
    umax_phys(it) = max(abs(Uplot(:)));

    % Plot field
    subplot(2,2,[1 3]);
    imagesc(x, y, Uplot);
    axis xy;
    axis equal tight;
    caxis([-clim, clim]);
    colorbar;
    xlabel('x');
    ylabel('y');
    title(sprintf('2D wave packet,  t = %.4f', t));

    % Plot total energy
    subplot(2,2,2);
    plot(tlist(1:it), energy(1:it), 'LineWidth', 1.5);
    grid on;
    xlim([tlist(1), tlist(end)]);
    xlabel('t');
    ylabel('E(t)');
    title('Total modal energy');

    % Plot relative energy drift
    subplot(2,2,4);
    plot(tlist(1:it), rel_drift(1:it), 'LineWidth', 1.5);
    grid on;
    xlim([tlist(1), tlist(end)]);
    xlabel('t');
    ylabel('(E(t)-E_0)/E_0');
    title('Relative energy drift');

    drawnow;

    % Advance one exact propagator step
    if it < nFrames
        Qnew = Cd .* Q + Sd_over_w .* P;
        Pnew = -wSd .* Q + Cd .* P;
        Q = Qnew;
        P = Pnew;
    end
end

% =========================
% Final textual summary
% =========================
fprintf('Initial energy           : %.16e\n', energy(1));
fprintf('Final energy             : %.16e\n', energy(end));
fprintf('Max abs relative drift   : %.16e\n', max(abs(rel_drift)));
fprintf('Max physical amplitude   : %.16e\n', max(umax_phys));

% =========================
% Optional separate figure for amplitude history
% =========================
figure('Name', 'Amplitude diagnostic');
plot(tlist, umax_phys, 'LineWidth', 1.5);
grid on;
xlabel('t');
ylabel('max |u(x,y,t)|');
title('Maximum physical-space amplitude');

% =========================
% End of script
% =========================
