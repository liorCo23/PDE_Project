%% Project 1: Solving Laplace and Poisson Equations in 2D using SOR
% ID: 314818345, so we set:
% n1 = 1, n2 = 1, n3 = 3, n4 = 3, n5 = 4.

clear; clc; close all;

%% Parameters from ID and Domain Setup
n1 = 1; n2 = 1; n3 = 3; n4 = 3; n5 = 4;

% Domain: x in [0,2.0], y in [0,1.0]
Lx = 2.0; Ly = 1.0;
Nx = 100;    % number of subintervals in x (=> Nx+1 grid points)
Ny = 50;     % number of subintervals in y (=> Ny+1 grid points)
dx = Lx / Nx; dy = Ly / Ny;  % here dx = dy = 0.02
h = dx;  % since dx = dy

% Create grid vectors
x = linspace(0, Lx, Nx+1);
y = linspace(0, Ly, Ny+1);

% SOR parameters
omega = 1.8;      % relaxation parameter (1 < omega < 2)
tol = 1e-6;       % tolerance for convergence
max_iter = 5000;  % maximum iterations

%% --- Solve the Laplace Equation (u_xx + u_yy = 0) ---
% Initialize solution U (rows correspond to x, columns to y)
U = zeros(Nx+1, Ny+1);

% Set boundary conditions:
% Bottom (y=0) and Top (y=1)
U(:,1) = 0;         % y = 0
U(:,end) = 0;       % y = 1

% Left boundary (x=0): u(0,y) = y^n1*(1-y)^n3 = y*(1-y)^3
U(1,:) = (y).^n1 .* (1 - y).^n3;
% Right boundary (x=2): u(2,y) = - y^n4*(1-y)^n2 = - y^3*(1-y)
U(end,:) = - (y).^n4 .* (1 - y).^n2;

% SOR iteration for Laplace (f=0)
iter = 0;
max_diff = 1;
while (max_diff > tol) && (iter < max_iter)
    max_diff = 0;
    iter = iter + 1;
    for i = 2:Nx        % interior in x
        for j = 2:Ny    % interior in y
            U_old = U(i,j);
            % Standard 5-point average (f=0 here)
            U_new = ( U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) ) / 4;
            % SOR update:
            U(i,j) = (1 - omega)*U_old + omega * U_new;
            max_diff = max(max_diff, abs(U(i,j) - U_old));
        end
    end
    if mod(iter,100)==0
        fprintf('Laplace iteration %d, max diff = %.2e\n', iter, max_diff);
    end
end
fprintf('Laplace converged in %d iterations with max diff = %.2e\n', iter, max_diff);

% Plot the Laplace solution
[X, Y] = meshgrid(y, x); % Note: for surf, rows correspond to x and columns to y
figure;
surf(X, Y, U);
xlabel('y'); ylabel('x'); zlabel('u(x,y)');
title('Laplace Equation Solution');
shading interp; colorbar;

%% Analysis: Identify x-values where u(x,y) (as a function of y)
% Preallocate enough space (the worst case is that every x in 2:Nx qualifies):
x_candidates = zeros(1, Nx-1);  
count = 0;  % how many x-values we actually store

for i = 2:Nx
    u_profile = U(i,:);
    if numel(u_profile) < 3
        continue;
    end
    du = diff(u_profile);

    has_max = any(du(1:end-1) > 0 & du(2:end) < 0);
    has_min = any(du(1:end-1) < 0 & du(2:end) > 0);

    if has_max && has_min
        count = count + 1;
        x_candidates(count) = x(i);  % store x(i)
    end
end

% Trim x_candidates to the actual number of stored elements
x_candidates = x_candidates(1:count);


if ~isempty(x_candidates)
    fprintf('For the Laplace problem, u(x,y) has both a local max and local min at x values:\n');
    disp(x_candidates);
else
    fprintf('No x values found with both a local max and local min in u(x,y).\n');
end

%% --- Solve the Poisson Equation (u_xx + u_yy = 1+10^n5) ---
% Here f = 1 + n5/10 = 1 + 4/10 = 1.4
f_val = 1 + (n5)/10;
fprintf('Solving Poisson with f = %d\n', f_val);

% Reinitialize solution U_p with the same boundary conditions
U_p = zeros(Nx+1, Ny+1);
U_p(:,1) = 0;
U_p(:,end) = 0;
U_p(1,:) = (y).^n1 .* (1 - y).^n3;
U_p(end,:) = - (y).^n4 .* (1 - y).^n2;

iter = 0;
max_diff = 1;
while (max_diff > tol) && (iter < max_iter)
    max_diff = 0;
    iter = iter + 1;
    for i = 2:Nx
        for j = 2:Ny
            U_old = U_p(i,j);
            % Update formula including the source term:
            % U_new = 1/4 * [ U(i+1,j) + U(i-1,j) + U(i,j+1) + U(i,j-1) - h^2 * f_val ]
            U_new = ( U_p(i+1,j) + U_p(i-1,j) + U_p(i,j+1) + U_p(i,j-1) - h^2 * f_val ) / 4;
            U_p(i,j) = (1 - omega)*U_old + omega * U_new;
            max_diff = max(max_diff, abs(U_p(i,j) - U_old));
        end
    end
    if mod(iter,100)==0
        fprintf('Poisson iteration %d, max diff = %.2e\n', iter, max_diff);
    end
end
fprintf('Poisson converged in %d iterations with max diff = %.2e\n', iter, max_diff);

% Find the minimum value of U_p and its location
[min_val, ind] = min(U_p(:));
[row_min, col_min] = ind2sub(size(U_p), ind);
x_min = x(row_min);
y_min = y(col_min);
fprintf('Poisson: Minimum u = %.4f at x = %.4f, y = %.4f\n', min_val, x_min, y_min);

% Plot the Poisson solution
figure;
surf(X, Y, U_p);
xlabel('y'); ylabel('x'); zlabel('u(x,y)');
title('Poisson Equation Solution');
shading interp; colorbar;

fprintf('Grid spacing: dx = %.4f, dy = %.4f (Nx = %d, Ny = %d)\n', dx, dy, Nx, Ny);
