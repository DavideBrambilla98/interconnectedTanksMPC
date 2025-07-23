clear;
clc;
close all;

set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked')
set(0,'defaulttextInterpreter','latex')
rng('default');

Ts = 5; % Tempo di campionamento

% Modello del sistema
addpath('funzioni');
modello;

% Matrici costo
Q = diag([1000, 1000, 1000, 1000]);
R = 1 * eye(2);

% Orizzonte
N = 15;

% Vincoli su stati e ingressi

% Costruzione MPC
mpc = MPC_uguaglianza(sys_d.A, sys_d.B, Hx, hx, Hu, hu, x_ref, u_ref, Q, R, N);

T_sim = 60;
x_log = zeros(4, T_sim+1);
u_log = zeros(2, T_sim);
flags = zeros(1, T_sim);

x_log(:, 1) = x_start - x_ref;  % stato centrato

for tt = 1:T_sim
    x_centrato = x_log(:, tt);

    f = real(mpc.f * x_centrato);
    b_ineq = real(mpc.b_ineq - ...
              [mpc.Hx_tilde * mpc.A_cal; zeros(size(mpc.Hu_tilde,1), size(sys_d.A,2))] * x_centrato);
    b_eq = real(-mpc.A_cal_n * x_centrato);

    [delta_u_seq, ~, exitflag] = quadprog(mpc.F, f, ...
        mpc.A_ineq, b_ineq, ...
        mpc.A_eq, b_eq);

    flags(tt) = exitflag;

    if isempty(delta_u_seq) || exitflag <= 0
        warning("Infeasible QP al passo %d", tt);
        delta_u_seq_first = zeros(2,1);
    else
        delta_u_seq_first = delta_u_seq(1:2);
    end

    u_real = u_ref + delta_u_seq_first;
    u_log(:, tt) = u_real;

    x_real = x_centrato + x_ref;
    dxdt = @(t,x) livSerbatoi(t, x, A, a, k, gamma, g, u_real);
    [~, xx] = ode45(dxdt, [0 Ts], x_real);

    x_log(:, tt+1) = xx(end, :)' - x_ref;
end

% Plot risultati
x1 = x_log(1,:) + x_ref(1);  x3 = x_log(3,:) + x_ref(3);
x2 = x_log(2,:) + x_ref(2);  x4 = x_log(4,:) + x_ref(4);

figure;
subplot(1,2,1)
plot(x1, x3, 'b'); hold on; scatter(x1, x3, 'c');
xlabel('$h_1$ [cm]'); ylabel('$h_3$ [cm]');
title('Traiettoria $h_1$ vs $h_3$');

subplot(1,2,2)
plot(x2, x4, 'r'); hold on; scatter(x2, x4, 'm');
xlabel('$h_2$ [cm]'); ylabel('$h_4$ [cm]');
title('Traiettoria $h_2$ vs $h_4$');

%% Andamento degli stati e degli ingressi

% Ricostruzione degli stati nel dominio reale
x_log_real = x_log + x_ref;  % [n x (T_sim+1)]

% Andamento degli stati
figure;
subplot(2,1,1);
plot((0:T_sim)*Ts/60, x_log_real');
title('Andamento degli stati');
xlabel('Tempo [min]');
ylabel('Stati [{$h_i$}]');
legend({'$h_1$', '$h_2$', '$h_3$', '$h_4$'}, 'Interpreter','latex');
grid on;

% Andamento degli ingressi
subplot(2,1,2);
plot((1:T_sim)*Ts/60, u_log');
title('Ingressi di controllo');
xlabel('Tempo [min]');
ylabel('Ingressi [$u_i$]');
legend({'$v_1$', '$v_2$'}, 'Interpreter','latex');
grid on;