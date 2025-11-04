close all; clear all;


%% Painleve II equation
t0 = 6;
tn = -8;
dydt = @(t, y) [y(2); t*y(1) + 2*y(1)^3; y(4); y(1)^2; -y(1)];
I0 = integral(@(x) airy(x).^2 .* (x - t0), t0, 10*t0);
J0 = integral(@(x) airy(x), t0, 10*t0);
y0 = [airy(t0); airy(1, t0); I0; airy(t0)^2; J0];
opts=odeset('reltol',1e-13,'abstol',1e-14);
[t, y] = ode45(dydt, [t0 tn], y0, opts);

%% beta = 2
F2 = exp(-y(:, 3));
f2 = -y(:, 4) .* F2;

num_trials = 100;
n = 1000;
tic
eig_vec= zeros(num_trials, 1);
for i=1:num_trials
    A = randn(n, n) + 1i * randn(n, n);
    A_sym = (A + A')/2;
    eig_vec(i) = n^(1/6) * (max(eigs(A_sym)) - 2*sqrt(n));
end
toc
figure
histogram(eig_vec, 'Normalization', 'pdf')
hold on
plot(t, f2)

%% beta = 1

F1 = sqrt(F2 .* exp(-y(:, 5)));
f1 = (f2 + y(:, 1) .* F2) .* exp(-y(:, 5)) ./ (2*F1);
eig_vec_2 = zeros(num_trials, 1);
tic
for i=1:num_trials
    A = randn(n, n);
    A_sym = (A + A')/(2);
    eig_vec_2(i) = n^(1/6) * (max(eigs(A_sym))*sqrt(2) - 2*sqrt(n));
end
toc
figure
histogram(eig_vec_2, 'Normalization', 'pdf')
hold on
plot(t, f1)

%% Laguerre beta = 1
m = n;
num_trials = 10000;
mu_mn = (sqrt(m-1)+sqrt(n))^2;
sigma_mn = (sqrt(m-1)+sqrt(n)) * ((1/sqrt(m-1)) + (1/n))^(1/3);
eig_vec_3 = zeros(num_trials, 1);
tic
for i=1:num_trials
    A = randn(m, n);
    A_sym = A*A';
    eig_vec_3(i) = (1/sigma_mn) * (max(eigs(A_sym)) - mu_mn);
end
toc
figure
histogram(eig_vec_3, 'Normalization', 'pdf')
hold on
plot(t, f1)

%% Laguerre beta = 2
m = n;
num_trials = 100;
mu_mn = (sqrt(m)+sqrt(n))^2;
sigma_mn = (sqrt(m)+sqrt(n)) * ((1/sqrt(m)) + (1/n))^(1/3);
eig_vec_4 = zeros(num_trials, 1);
tic
for i=1:num_trials
    %A = rand(m, n) + 1i * rand(m, n);
    A = gen_lag_bidiag(n, beta);
    A_sym = A*A';
    eig_vec_4(i) = (1/sigma_mn) * (max(eigs(A_sym)) - mu_mn);
end
toc
figure
histogram(eig_vec_4, 'Normalization', 'pdf')
hold on
plot(t, f1)





%% Painleve II sigma solution
p = y(:, 2) + y(:, 1).^2 + t;
H_II = -1/2 * (2*y(:, 1).^2 - p + t) - (y(:, 1))/2;
sigma_II = -2^(1/3) * H_II; % Have to plot with -2t^(1/3)
figure
plot(t/(-2^(1/3)), sigma_II)

%% Painleve V sigma form

dydt = @(t, y) [y(2); ...
    (2/t).*sqrt((y(1) - t.*y(2)).*(t.*y(2) - y(1) - (y(2)).^2)); ...
    y(1) ./ t];
t0 = 10^(-12);
tn = 16;
y0 = [(t0/pi) + (t0/pi)^2; (1/pi) + (2*t0/pi^2); (t0/pi) + (t0^2/2*pi^2)];
opts=odeset('reltol',1e-13,'abstol',1e-14);
[t, sigma] = ode45(dydt, [t0, tn], y0, opts);
plot(t, sigma(:, 1))
E = exp(-sigma(:, 3));
p = -(pi./t).^2 .* (t .* sigma(:, 2) - sigma(:, 1) - (sigma(:, 1)).^2) .* E;
plot(t/pi, p)
hold on
plot(t/pi, E)

n=100;
nrep=10000;
beta=2;
ds=zeros(1,nrep*n/2);
for ii=1:nrep
    chi_vals = sqrt(chi2rnd((n-1:-1:1)'*beta)/2);
    A = diag(randn(n,1)) + diag(chi_vals, -1) + diag(chi_vals, 1);
    l = eig(A);
    d=diff(l(n/4:3*n/4))/beta/pi.*sqrt(2*beta*n-l(n/4:3*n/4-1).^2);
    ds((ii-1)*n/2+1:ii*n/2)=d;
end
figure
histogram(ds, 'Normalization', 'pdf')
hold on
plot(t/pi, p)


%% Painleve III
N = 20;
a = 1;
dydt = @(t, y) [y(2); ...
    -(1/t) .* sqrt((a*y(2))^2 - 4*y(2)*(1+y(2))*(t*y(2)-y(1))); y(1) / t + a/(2*N)*y(2)];
t0 = 10^(-8);
tn = 10;
I0 = -(1+(a+1)*a/(2*N)) / (a+1) * t0^(a+1);
y0 = [-t0^(a+1); -t0^a * (a+1); I0] / gamma(a+1) / gamma(a+2);
opts=odeset('reltol',1e-3,'abstol',1e-4);
[t, y] = ode45(dydt, [t0, tn], y0, opts);

figure
plot(t, y(:, 3))
E_III = exp(y(:, 3));
figure
plot(t, E_III)

%% Painleve III sigma
a = 0.5;
dydt = @(t, y) [y(2); ...
    (1/t) .* sqrt((a*y(2))^2 - y(2)*(y(1) - t*y(2))*(4*y(2)-1)); ...
    y(1) / t];
t0 = 10^(-7);
tn = 10;
asymp_func = @(x, z) z / (gamma(a+1)*gamma(a+2)) * (x/4).^(a+1);
asymp_func_der = @(x, z) z / (gamma(a+1)^2) * x.^a / 4^(a+1);
I0 = 1/(gamma(a+2))^2 * (t0/4)^(a+1);
y0 = [asymp_func(t0, 1); asymp_func_der(t0, 1); I0];
opts=odeset('reltol',1e-3,'abstol',1e-4);
[t, sigma] = ode45(dydt, [t0, tn], y0, opts);

figure
plot(t, sigma(:, 1))

figure
plot(t, sigma(:, 3))
D_III = exp(sigma(:, 3));
figure
plot(t, D_III)


%% Painleve V Laguerre

N = 20;
alpha = 0;
dydt = @(t, y) [y(2); ...
    -(1/t).*sqrt((y(1) - y(2).*(t + 2*y(2) - 2*N - alpha)).^2 - ...
                4*((y(2)).^2).*(y(2) - N).*(y(2) - N - alpha)); ...
    y(1) ./ t];
t0 = 10^(-12);
tn = 5;
y0 = [N*t0/2; N/2; N*t0/2];
%{
y0 = [gamma(N+alpha+1)/(gamma(N)*gamma(alpha+1)*gamma(alpha+2)) * t0^(alpha+1); ...
    gamma(N+alpha+1)/(gamma(N)*gamma(alpha+1)*gamma(alpha+2))*(alpha+1) * t0^alpha; ...
    -(t0/pi) - (t0^2/2*pi^2)];
%}
opts=odeset('reltol',1e-3,'abstol',1e-4);
[t, sigma] = ode45(dydt, [t0, tn], y0, opts);
figure
plot(t, sigma(:, 1))
figure
plot(t, exp(-sigma(:, 3)))


%% Painleve IV sigma form

N = 10;
dydt = @(t, y) [y(2); -2*sqrt((t.*y(2)-y(1)).^2 - 4*(y(2)).^2 .* (y(2) + 2*N)); -y(1); -y(2)];
t0 = 1.8*sqrt(N-1);
tn = 0;
asymp_func = @(x) 2^(N-1) .* (x.^(2*N-2)) .* exp(-x.^2) / sqrt(pi) / gamma(N);
I0 = integral(asymp_func, t0, 10*t0);
asymp_func_der = @(x) 2^(N-1) * ((2*N-2) .* (x.^(2*N-3)) - 2*(x.^(2*N-1))) .* exp(-x.^2) / sqrt(pi) / gamma(N);
y0 = [asymp_func(t0); asymp_func_der(t0); I0; -asymp_func(t0)];
opts=odeset('reltol',1e-6,'abstol',1e-8);
[t, sigma] = ode45(dydt, [t0, tn], y0, opts);

plot(t, real(sigma(:, 1)))
p_test = sigma(:, 4) .* exp(-real(sigma(:, 3)));
plot(t, p_test)
%{
daeform = @(t, y, yp) [yp(1) - y(2);
    yp(2)^2 - 4*(y(1) - t*y(2)) + 4*(y(2))^2 * (y(2) + 2*N);
    yp(3) + y(1)];

asymp_func_der_der = @(x) 2^(N-1) * (((2*N-2) .* (x.^(2*N-3)) - 2*x.* (x.^(2*N-2))) .* (-2*x) + ...
                                     ((2*N-2) * (2*N-3) .* (x.^(2*N-4)) - 2*(2*N-1)*(x.^(2*N-2)))) .* exp(-x.^2) / sqrt(pi) / gamma(N);
yp0 = [asymp_func_der(t0); asymp_func_der_der(t0); -asymp_func(t0)];
%[y0_new, yp0_new] = decic(daeform, t0, y0, [1 1 1], yp0, [0 0 0]);

tspan = linspace(t0, tn, 1000);
%[t, sigma] = ode15i(daeform, tspan, y0, yp0, opts);

% Solving for Hamiltonian to get Painleve IV sigma



dhdt = @(t, y) [y(2); 2*sqrt((t.*y(2)-y(1)).^2 - 4*(y(2) - 2*N/3).^2 .* (y(2) + 4*N/3))];
h0 = [y0(1) - 2*N*t0/3; y0(2) - 2*N/3];

opts=odeset('reltol',1e-3,'abstol',1e-4);
[t, h_sol] = ode23(dhdt, [t0, tn], h0, opts);

hp0 = [h0(2); yp0(2)];
daeform_h = @(t, h, hp) [hp(1) - h(2);
    hp(2)^2 - 4*(t*h(2) - h(1))^2 + 4*(h(2) + 4*N/3)*(h(2)-2*N/3)^2];
%[h0_new, hp0_new] = decic(daeform_h, t0, h0, [1 0], hp0, [0 0]);
[t, h_IV] = ode15i(daeform_h, [t0, tn], h0, hp0, opts);

figure
plot(t, h_sol(:, 1))
plot(t, h_sol(:, 1) + 2*N*t/3)
%}

%% Generate Laguerre B matrix

function [B_n] = gen_lag_bidiag(n, beta)
    a = ceil(beta/2 * (n+1));
    diag_dofs = 2*a - beta*(1:n);
    off_diag_dofs = beta*((n-1):-1:1);
    main_diag = sqrt(chi2rnd(diag_dofs));
    off_diag = sqrt(chi2rnd(off_diag_dofs));
    B_n = diag(main_diag) + diag(off_diag, -1);
end
