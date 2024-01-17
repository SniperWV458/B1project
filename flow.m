clear;
clc;

%initializing known values
N = 64;
U0 = 1; %m/s
L0 = 1; %m
rou = 1; %kg/m**3
re = 5;
mu = rou*U0*L0/re;
delta_x = L0/(N+1); %m
psi_grid = zeros(N+2, N+2);

%usage
%{
max_it = 100000;
extrapolation_factor = 1.919;
omega_grid = get_omega_grid(N);
psi_grid = force_iterate(omega_grid, psi_grid, N, delta_x, max_it, extrapolation_factor);
partial_ux = back_diff(psi_grid, N, delta_x);
force1 = get_force(partial_ux, N, mu);

N = 32;
delta_x = L0/(N+1);
psi_grid = zeros(N+2, N+2);
omega_grid = get_omega_grid(N);
psi_grid = force_iterate(omega_grid, psi_grid, N, delta_x, max_it, extrapolation_factor);
partial_ux = back_diff(psi_grid, N, delta_x);
force2 = get_force(partial_ux, N, mu);
ref_psi_grid = generate_ref_grid(N, delta_x);

N = 16;
max_it = 140000;
delta_x = L0/(N+1);
psi_grid = zeros(N+2, N+2);
omega_grid = get_omega_grid(N);
psi_grid = force_iterate(omega_grid, psi_grid, N, delta_x, max_it, extrapolation_factor);
partial_ux = back_diff(psi_grid, N, delta_x);
force3 = get_force(partial_ux, N, mu);

resolution = [1/65,1/33,1/17];
forces = [force1, force2, force3];

p = polyfit(resolution,forces,2);
x_dense = linspace(0, max(resolution), 100);
y_dense = polyval(p, x_dense);
plot(resolution,forces, 'ro', 'MarkerSize', 10, 'LineWidth', 2); hold on;
plot(x_dense, y_dense, 'b-', 'LineWidth', 2);
title('Best Fit Curve through Three Points')
xlabel('resolution')
ylabel('force')
%}

%{
test_omega_grid = generate_test_grid(N, delta_x);
[psi_grid, times] = iterate(test_omega_grid, psi_grid, N, delta_x, 10000, 1.9);
max_it = 60000;
test_omega_grid = generate_test_grid(16, delta_x);
low_psi_grid = force_iterate(test_omega_grid, zeros(18,18), 16, delta_x, max_it, 1.9);
contour_level = (0.1:0.1:0.9);
subplot(2,2,1)
plot_contour(psi_grid,N,'calculated',contour_level)
subplot(2,2,2)
plot_contour(low_psi_grid,16,'max it')
subplot(2,2,3)
ref_psi_grid = generate_ref_grid(N, delta_x);
plot_contour(ref_psi_grid,N,'theory',contour_level)
subplot(2,2,4)
plot_surf(ref_psi_grid-psi_grid,N,'error')
%}

%{
split = 100;
max_it = 60000;
test_omega_grid = generate_test_grid(N, delta_x);
[factor, times] = test_extrapolation(test_omega_grid,psi_grid,N,delta_x,max_it,split);
plot(factor,times)
xlabel('extrapolation factor')
ylabel('iterations')
%}


max_it = 120000;
extrapolation_factor = 1.919;
[N_size, abs_error] = test_grid_size(max_it, extrapolation_factor);
subplot(2,5,10)
loglog(N_size,abs_error)
xlabel('N size')
ylabel('avg abs error')


%jacobi function
function next_psi_grid = jacobi_sweep(omega_grid, psi_grid, N, delta_x)
new_psi_grid = ones(N+2,N+2);
    for i = 2:N+1
        for j = 2:N+1
            new_psi_grid(i,j) = 0.25*(psi_grid(i+1,j)+psi_grid(i,j+1)+psi_grid(i-1,j)+psi_grid(i,j-1)+omega_grid(i-1,j-1)*delta_x^2);
        end
    end
next_psi_grid = new_psi_grid;
end

%Gauss-Seidel function
function next_psi_grid = gs_sweep(omega_grid, psi_grid, N, delta_x)
    for i = 2:N+1
        for j = 2:N+1
            psi_grid(i,j) = 0.25*(psi_grid(i+1,j)+psi_grid(i,j+1)+psi_grid(i-1,j)+psi_grid(i,j-1)+omega_grid(i-1,j-1)*delta_x^2);
        end
    end
    next_psi_grid = psi_grid;
end

%Succesive over-relaxation function
function next_psi_grid = sor_sweep(omega_grid, psi_grid, N, delta_x, extrapolation_factor)
    for i = 2:N+1
        for j = 2:N+1
            gs_value = 0.25*(psi_grid(i+1,j)+psi_grid(i,j+1)+psi_grid(i-1,j)+psi_grid(i,j-1)+omega_grid(i-1,j-1)*delta_x^2);
            psi_grid(i,j) = extrapolation_factor*gs_value + (1-extrapolation_factor)*psi_grid(i,j);
        end
    end
    next_psi_grid = psi_grid;
end

%check diff
function [diff, convergence] = check_diff(omega_grid, psi_grid, N, delta_x, extrapolation_factor)
convergence_count = 0;
current_grid = psi_grid;
next_grid = sor_sweep(omega_grid, current_grid, N, delta_x, extrapolation_factor);
diff = zeros(N,N);
    for i = 2:N+1
        for j = 2:N+1
            diff(i-1,j-1) = next_grid(i,j) - current_grid(i,j);
            %set covergence diff
            if diff(i-1,j-1)<10^-6
                convergence_count = convergence_count + 1;
            end
        end
    end
    if convergence_count>N*N*2/3
        convergence = true;
    else
        convergence = false;
    end
end

%plotter
function plot_surf(psi_grid,N, graph_name)
x = linspace(0, 1, N);
y = linspace(0, 1, N);
[X, Y] = meshgrid(x, y);
surf(X, Y, psi_grid(2:N+1,2:N+1))
xlabel('x')
ylabel('y')
zlabel('ψ(x, y)')
title(graph_name)
end

function plot_contour(psi_grid, N, graph_name, contour_level)
x = linspace(0, 1, N+2);
y = linspace(0, 1, N+2);
[X, Y] = meshgrid(x, y);
if nargin == 4
    contour(X, Y, psi_grid,contour_level)
else
    contour(X, Y, psi_grid)
end
xlabel('x')
ylabel('y')
zlabel('ψ(x, y)')
title(graph_name)
end

%run till conversion
function [psi_grid, times] = iterate(omega_grid, psi_grid, N, delta_x, max_it,extrapolation_factor)
for times = 1:max_it
    psi_grid = sor_sweep(omega_grid, psi_grid, N, delta_x, extrapolation_factor);
    if mod(times, 50) == 0
        [diff, convergence] = check_diff(omega_grid, psi_grid,N, delta_x, extrapolation_factor);
        if convergence == true
            assignin("base", "diff", diff)
            break
        end
    end
end
end

%run set iterations
function psi_grid = force_iterate(omega_grid, psi_grid, N, delta_x, max_it,extrapolation_factor)
for times = 1:max_it
    psi_grid = sor_sweep(omega_grid, psi_grid, N, delta_x, extrapolation_factor);
end
end

%vary extrapolation
function [factor_list, times_list] = test_extrapolation(omega_grid, psi_grid, N, delta_x, max_it, split)
factor_list = linspace(0,2,split);
times_list = zeros(split,1);
count = 1;
for factor = factor_list
    [~, times] = iterate(omega_grid, psi_grid, N, delta_x, max_it, factor);
    times_list(count) = times;
    count = count+1;
end
end

%vary grid size
function [N_size, abs_error] = test_grid_size(max_it, extrapolation_factor)
N_size = linspace(4,128,32);
abs_error = zeros(32, 1);
count = 1;
plot_count = 1;
for N = N_size
    delta_x = 1/(N+1);
    omega_grid = generate_test_grid(N, delta_x);
    psi_grid = zeros(N+2,N+2);
    ref_psi_grid = generate_ref_grid(N, delta_x);

    [psi_grid ,~] = iterate(omega_grid, psi_grid, N, delta_x, max_it, extrapolation_factor);
    
    if mod(count,4) == 0
        subplot(2,5,plot_count)
        plot_surf(psi_grid,N,N)
        plot_count = plot_count + 1;
    end
    error_sum = 0;
    for i = 1:N+2
        for j = 1:N+2
            error_sum = error_sum +abs(psi_grid(i,j)-ref_psi_grid(i,j));
        end
    end
    abs_error(count) = error_sum/(N*N);
    count = count + 1;
    
end
subplot(2,5,9)
ref_psi_grid = generate_ref_grid(128, 1/129);

plot_surf(ref_psi_grid,128,'reference')
end

%setup test omega grid
function test_omega_grid = generate_test_grid(N, delta_x)
test_omega_grid = zeros(N, N);
for x = 1:N
    for y = 1:N
        x_value = x*delta_x;
        y_value = x*delta_x;
        test_omega_grid(x, y) = (2*pi^2)*sin(pi*x_value)*sin(pi*y_value);
    end
end
end

%setup referencs psi grid
function psi_grid = generate_ref_grid(N, delta_x)
%initialize reference
psi_grid = zeros(N+2, N+2);
for x = 2:N+1
    for y = 2:N+1
        x_value = x*delta_x;
        y_value = y*delta_x;
       psi_grid(x, y) = sin(pi*x_value)*sin(pi*y_value);
    end
end
end

%import function
function omega_grid = get_omega_grid(N)
file_name = "omegaN"+ string(N) +".dat";
fileid = fopen(file_name);
dat = fread (fileid, "single");
omega_grid = reshape(dat,[N+1,N+1]);
end

%second-order accurate backward difference
function partial_ux = back_diff(psi_grid, N, delta_x)
partial_ux = ones(N+2, 1);
for i = 2:N+1
    partial_ux(i) = (2*psi_grid(i,N+1)-5*psi_grid(i,N)+4*psi_grid(i,N-1)-psi_grid(i,N-2))/delta_x^2;
end
end

%force function
function force = get_force(partial_ux, N, mu)
x_values = linspace(0, 1, N+2);
force = mu * trapz(x_values, partial_ux);
end