%% Project 2: 1D Heat Equation (Homogeneous and Non-Homogeneous)
% ID: 314818345, so:
% n1 = 1, n2 = 1, n3 = 3, n4 = 3, n5 = 4.

clear; clc; close all;

%% Parameters and Grid Setup
n1 = 1; n2 = 1; n3 = 3; n4 = 3; n5 = 4;

% Spatial domain: x in [0,1]
L = 1;
Nx = 100; 
dx = L / Nx;
x = linspace(0, L, Nx+1);

% Stability for explicit Euler for the heat equation requires dt/dx^2 < 0.5.
% We choose r = dt/dx^2 = 0.4.
r = 0.4;
dt = r * dx^2;
fprintf('Using dt = %e with r = %f\n', dt, r);

%% --- Solve the Homogeneous Heat Equation: u_t = u_xx ---
% Initial condition: u(x,0) = 5*x^(n1)*(1-x)^(n3) = 5*x*((1-x)^3)*(1/2 - x)
u = 5 * (x.^n1) .* ((1-x).^n3) .*((1/2)-x);
u_new = u;  

% Apply boundary conditions at t=0 (though they will be updated in time)
% u(0,t) = t^(1+n4)/(1+t); at t = 0, we set u(0,0) = 0.
u(1) = 0;  % x=0
u(end) = n2;  % x=1, u(1,t) = 1

% Prepare to store solutions for plotting
store_interval = 50;   % store every 50 time steps
sol_matrix = [];
time_store = [];
time = 0;

T_monotonic = [];  % will record the first time when u(x) becomes monotonic increasing
found_T = false;

max_steps = 50000;  % maximum time steps

for step = 1:max_steps
    time = time + dt;
    
    % Update interior points using explicit Euler:
    % u_i^{n+1} = u_i^n + r*(u_{i+1} - 2*u_i + u_{i-1})
    for i = 2:Nx
        u_new(i) = u(i) + r * ( u(i+1) - 2*u(i) + u(i-1) );
    end
    
    % Update boundary conditions:
    % Left boundary: u(0,t) = t^(1+n4)/(1+t) = t^4/(1+t)
    u_new(1) = (time)/(1+n4);
    % Right boundary: u(1,t) = n2 = 1.
    u_new(end) = time/n2;
    
    % Update solution:
    u = u_new;
    
    % Store solution every "store_interval" steps
    if mod(step, store_interval)==0
        sol_matrix = [sol_matrix; u];
        time_store = [time_store; time];
    end
    
    % Check if u(x) is monotonic increasing (within a tolerance)
    if ~found_T
        if all(diff(u) >= -1e-8)
            T_monotonic = time;
            found_T = true;
            fprintf('Homogeneous: Monotonic increasing profile first detected at t = %.4f\n', T_monotonic);
            % Continue simulation until t >= 2*T_monotonic.
        end
    end
    
    if found_T && (time >= 2*T_monotonic)
        break;
    end
end

if ~found_T
    fprintf('Homogeneous: Monotonic increasing profile not detected within simulation time.\n');
else
    fprintf('Homogeneous: Final simulation time t = %.4f (T_monotonic = %.4f)\n', time, T_monotonic);
end

% Plot the homogeneous heat equation solution
figure;
[X_mesh, T_mesh] = meshgrid(x, time_store);
surf(X_mesh, T_mesh, sol_matrix);
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
title('1D Heat Equation Solution (Homogeneous)');
shading interp; colorbar;

%% --- Solve the Non-Homogeneous Heat Equation: u_t = u_xx + (n5*x)/5 ---
% Here, (n5*x)/5 = (4*x)/5.
u2 = 5 * (x.^n1) .* (((1-x).^n3)) .* ((1/2)-x) ;  % same initial condition
u2_new = u2;
% Set initial boundary values
u2(1) = 0;
u2(end) = n2;

sol_matrix2 = [];
time_store2 = [];
time = 0;
found_T2 = false;
T_monotonic2 = [];

for step = 1:max_steps
    time = time + dt;
    
    % Update interior points with source term:
    % u_i^{n+1} = u_i^n + r*(u_{i+1} - 2*u_i + u_{i-1}) + dt*(n5*x(i))/5
    for i = 2:Nx
        u2_new(i) = u2(i) + r*( u2(i+1) - 2*u2(i) + u2(i-1) ) + dt*( n5 * x(i) / 5 );
    end
    
    % Update boundaries:
    u2_new(1) = (time)/(1+n4);
    u2_new(end) = time/n2;
    
    % Update solution:
    u2 = u2_new;
    
    if mod(step, store_interval)==0
        sol_matrix2 = [sol_matrix2; u2];
        time_store2 = [time_store2; time];
    end
    
    % Check for monotonicity
    if ~found_T2
        if all(diff(u2) >= -1e-8)
            T_monotonic2 = time;
            found_T2 = true;
            fprintf('Non-homogeneous: Monotonic increasing profile first detected at t = %.4f\n', T_monotonic2);
        end
    end
    
    if (time >= 2*T_monotonic2)
        break;
    end
end

if ~found_T2
    fprintf('Non-homogeneous: Monotonic increasing profile not detected within simulation time.\n');
else
    fprintf('Non-homogeneous: Final simulation time t = %.4f (T_monotonic = %.4f)\n', time, T_monotonic2);
end

% Plot the non-homogeneous heat equation solution
figure;
[X_mesh2, T_mesh2] = meshgrid(x, time_store2);
surf(X_mesh2, T_mesh2, sol_matrix2);
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
title('1D Heat Equation Solution (Non-homogeneous)');
shading interp; colorbar;

fprintf('Spatial step dx = %.4f and time step dt = %e\n', dx, dt);
if found_T
    fprintf('Homogeneous: Monotonic increasing profile first detected at t = %.4f\n', T_monotonic);
end
