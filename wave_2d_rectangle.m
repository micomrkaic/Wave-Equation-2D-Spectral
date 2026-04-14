% wave2d_rectangle_packet.m
%
% 2D wave equation on a rectangle with fixed walls:
%     u_tt = c^2 (u_xx + u_yy)
%     u = 0 on the boundary
%
% Solved by double sine-series expansion.
% Initial condition: traveling Gaussian-modulated carrier wave packet.
%
% This script is written to run in both MATLAB and GNU Octave.
%
% -------------------------------------------------------------------------
clear;
clc;
close all;
% =========================
% Physical/domain settings
% =========================
Lx = 2.0;          % domain length in x
Ly = 3.8;          % domain length in y
c  = 2.0;          % wave speed
% =========================
% Mode truncation
% =========================
% Increase M,N for narrower packets or higher carrier wavenumber.
M = 40;            % number of sine modes in x
N = 40;            % number of sine modes in y
% =========================
% Quadrature grid for projection
% =========================
% This grid is used to project the initial conditions onto the sine basis.
Nxq = 220;
Nyq = 180;
xq = linspace(0, Lx, Nxq);
yq = linspace(0, Ly, Nyq);
dx = xq(2) - xq(1);
dy = yq(2) - yq(1);
[Xq, Yq] = meshgrid(xq, yq);   % size Nyq x Nxq
% =========================
% Packet parameters
% =========================
A0 = 1.0;          % amplitude
x0 = 0.28 * Lx;    % initial center x
y0 = 0.35 * Ly;    % initial center y
sigx = 0.08 * Lx;  % Gaussian width in x
sigy = 0.08 * Ly;  % Gaussian width in y
theta = pi/6;      % propagation direction angle
lambda0 = 0.12;    % carrier wavelength
k0 = 2*pi/lambda0; % carrier magnitude
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
% To make the packet start traveling in direction (dirx, diry),
% use ut = -c * (dir · grad u0).
%
% For u0 = env * cos(phase):
%   ux = env_x * cos(phase) - env * sin(phase) * kx0
%   uy = env_y * cos(phase) - env * sin(phase) * ky0
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
% Sxq: Nxq x M, Syq: Nyq x N
Sxq = sin((xq(:) * (mx * pi / Lx)));
Syq = sin((yq(:) * (ny * pi / Ly)));
% Trapezoidal-rule weights
wx = ones(Nxq, 1);
wy = ones(Nyq, 1);
wx(1)   = 0.5;
wx(end) = 0.5;
wy(1)   = 0.5;
wy(end) = 0.5;
% Put weights into matrices matching u0, v0 size (Nyq x Nxq)
% U-weighting is easiest if we keep x in columns and y in rows.
W = wy * wx.';   % size Nyq x Nxq
U0w = u0 .* W;
V0w = v0 .* W;
% =========================
% Project initial data onto modes
% =========================
% We want:
%   A_mn = (4/Lx/Ly) \int\int u0 sin(m pi x/Lx) sin(n pi y/Ly) dx dy
%   G_mn = (4/Lx/Ly) \int\int v0 sin(m pi x/Lx) sin(n pi y/Ly) dx dy
%   B_mn = G_mn / omega_mn
%
% Matrix form:
%   A = (4/Lx/Ly) * dx*dy * Sxq' * U0w' * Syq
%
% Dimensions:
%   Sxq'   : M x Nxq
%   U0w'   : Nxq x Nyq
%   Syq    : Nyq x N
% Result   : M x N
A = (4/(Lx*Ly)) * dx * dy * (Sxq.' * U0w.' * Syq);
G = (4/(Lx*Ly)) * dx * dy * (Sxq.' * V0w.' * Syq);
% =========================
% Modal frequencies
% =========================
[MM, NN] = ndgrid(mx, ny);
omega = c * pi * sqrt((MM/Lx).^2 + (NN/Ly).^2);
B = G ./ omega;
% =========================
% Display grid
% =========================
Nx = 180;
Ny = 150;
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
Sx = sin((x(:) * (mx * pi / Lx)));   % Nx x M
Sy = sin((y(:) * (ny * pi / Ly)));   % Ny x N
% =========================
% Time settings
% =========================
Tfinal = 2.0;
nFrames = 240;
tlist = linspace(0, Tfinal, nFrames);
% =========================
% Optional: estimate a color scale
% =========================
umax = max(abs(u0(:)));
if umax == 0
    umax = 1;
end
clim = 1.1 * umax;
% =========================
% Figure setup
% =========================
fig = figure('Name', '2D Wave Packet on a Rectangle');
% Some Octave versions ignore this harmlessly
try
    set(fig, 'Color', 'w');
catch
end
% =========================
% Animation loop
% =========================
for it = 1:nFrames
    t = tlist(it);
    % Modal time factor
    C = A .* cos(omega * t) + B .* sin(omega * t);
    % Reconstruct solution:
    % U = sum_{m,n} C_mn sin(m pi x/Lx) sin(n pi y/Ly)
    % Matrix form:
    % U = Sx * C * Sy'
    U = Sx * C * Sy.';     % size Nx x Ny
    Uplot = U.';           % transpose for plotting as y-rows, x-columns
    % Plot
    imagesc(x, y, Uplot);
    axis xy;
    axis equal tight;
    caxis([-clim, clim]);
    colorbar;
    xlabel('x');
    ylabel('y');
    title(sprintf('2D wave packet on a rectangle,  t = %.4f', t));
    drawnow;
end
% =========================
% End of script
% =========================
