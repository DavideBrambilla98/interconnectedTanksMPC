clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
set(0,'defaulttextInterpreter','latex')
rng('default');

%% Invariant set per il sistema controllato

% Tempo di campionamento
Ts = 60; %[s]

% Richiamo del modello del sistema dei serbatoio interconnessi
addpath('funzioni');
modello;

% Matrici del costo quadratico
Q = 1e1*eye(4); % Penalizza lo stato (quanto gli stati h1,h2,h3,h4 devono essere vicino al riferimento)
R = 1e3*eye(2); % Penalizza l'ingresso (quanto limitare l'uso degli ingressi v1 e v2)

% Control invariant set CIS_H*x <= CIS_h
[CIS_H, CIS_h] = cis(sys_d.A, sys_d.B, zeros(4,1), zeros(2,1), Hx, hx, Hu, hu, Q, R); % si passano zeri come riferimento poichè il sistema è traslato sul riferimento

%% Plot del CIS

CIS_G = Polyhedron(CIS_H, CIS_h);
CIS_G = CIS_G.minHRep();
CIS_G_H13 = projection(CIS_G, [1,3]);
CIS_G_H24 = projection(CIS_G, [2,4]);

dim = CIS_G.Dim;
disp(['La dimensione del poliedro è: ', num2str(dim)]);

% Plot delle proiezioni
figure

subplot(1 , 2 , 1)
CIS_G_H13.plot();
title("Proiezione del CIS dei livelli dei serbatoi 1 e 3")
grid on
axis equal
xlabel("$h_1$" , Interpreter="latex")
ylabel("$h_3$" , Interpreter="latex")

subplot(1 , 2 , 2)
CIS_G_H24.plot();
title("Proiezione del CIS dei livelli dei serbatoi 2 e 4")
grid on
axis equal
xlabel("$h_2$" , Interpreter="latex")
ylabel("$h_4$" , Interpreter="latex")

%% N-step controllable set

% Orizzonte di predizione fissato
N = 10;
fprintf('\n--- Calcolo del %d-step controllable set ---\n', N);
[Np_steps_H, Np_steps_h] = controllable_set(Hx, hx, Hu, hu, CIS_H, CIS_h, sys_d.A, sys_d.B, N);
fprintf('Vincoli nel %d-step set: %d\n', N, size(Np_steps_H,1));

% Costruzione poliedro
Np_steps_set = Polyhedron(Np_steps_H, Np_steps_h);
Np_steps_set = Np_steps_set.minHRep();

% Plot N-step-controllable set
figure()
subplot(1 , 2 , 1)
h_cis_13 = CIS_G_H13.plot();
h_nsteps_13 = projection(Np_steps_set, [1,3]);
hold on
h_nsteps_13 = h_nsteps_13.plot('Alpha', 0, 'LineWidth', 2);
title(sprintf('CIS e %d-step controllable set del sistema linearizzato serbatoi 1 e 3: %d\n', N))
xlabel("$h_1$" , Interpreter="latex")
ylabel("$h_3$" , Interpreter="latex")
legend([h_cis_13, h_nsteps_13], ...
    {'CIS', sprintf('%d-step set', N)},...
    'Interpreter','latex')
subplot(1 , 2 , 2)
h_cis_24 = CIS_G_H24.plot();
h_nsteps_24 = projection(Np_steps_set, [2,4]);
hold on
h_nsteps_24 = h_nsteps_24.plot('Alpha', 0, 'LineWidth', 2);
title(sprintf('CIS e %d-step controllable set del sistema linearizzato serbatoi 2 e 4: %d\n', N))
xlabel("$h_2$" , Interpreter="latex")
ylabel("$h_4$" , Interpreter="latex")
legend([h_cis_24, h_nsteps_24], ...
    {'CIS', sprintf('%d-step set', N)},...
    'Interpreter','latex')
%% mpc con vincolo di uguaglianza

T_sim = 60;
mpc_uguaglianza = MPC_uguaglianza(sys_d.A,sys_d.B,Hx,hx,Hu,hu,CIS_H,CIS_h,x_ref,u_ref,Q,R,N);

Hx_tilde = mpc_uguaglianza.Hx_tilde;
Hu_tilde = mpc_uguaglianza.Hu_tilde;
A_cal = mpc_uguaglianza.A_cal;
A_cal_n = mpc_uguaglianza.A_cal_n;

%% Verifica dei vincoli MPC con uguaglianza terminale

disp('--- Verifica vincoli MPC (vincolo terminale di uguaglianza) ---');
disp('Max valore b_ineq:');
disp(max(mpc_uguaglianza.b_ineq));
disp('Min valore b_ineq:');
disp(min(mpc_uguaglianza.b_ineq));

disp('Dimensioni A_ineq:');
disp(size(mpc_uguaglianza.A_ineq));
disp('Dimensioni b_ineq:');
disp(size(mpc_uguaglianza.b_ineq));

disp('Dimensioni A_eq:');
disp(size(mpc_uguaglianza.A_eq));
disp('Dimensioni b_eq:');
disp(size(mpc_uguaglianza.b_eq));

disp('hx centrati:');
disp(hx - Hx*x_ref);

disp('hu centrati:');
disp(hu - Hu*u_ref);


%% Simulazione MPC con vincolo terminale di uguaglianza

x_log = zeros(4, T_sim+1); % Stati centrati
u_log = zeros(2, T_sim);   % Ingressi reali
flags = zeros(1, T_sim);   % Exitflag di quadprog

x_log(:, 1) = x_start - x_ref; % Stato iniziale centrato

for tt = 1:T_sim
    % Stato attuale centrato
    x_centrato = x_log(:, tt);

    % Costruzione f e b_ineq
    f = real(mpc_uguaglianza.f * x_centrato);
    b_ineq = real(mpc_uguaglianza.b_ineq - ...
              [Hx_tilde * A_cal; zeros(size(Hu_tilde,1), size(A,2))] * x_centrato);

    % Costruzione b_eq centrato
   b_eq = real(-mpc_uguaglianza.A_cal_n * x_centrato);


    % Risoluzione QP
    [delta_u_seq, ~, exitflag] = quadprog(mpc_uguaglianza.F, f, ...
    mpc_uguaglianza.A_ineq, b_ineq, ...
    mpc_uguaglianza.A_eq, b_eq);

    flags(tt) = exitflag;

    if isempty(delta_u_seq) || exitflag <= 0
        warning("Quadprog non ha trovato una soluzione ammissibile al passo %d (exitflag = %d)", tt, exitflag);
        delta_u_seq_first = zeros(2,1); % fallback
    else
        delta_u_seq_first = delta_u_seq(1:2);
    end

    % Calcolo ingresso reale
    u_real = u_ref + delta_u_seq_first;
    u_log(:, tt) = u_real;

    % Simulazione del sistema non lineare nel dominio reale
    x_real = x_centrato + x_ref;

    dxdt = @(t,x) livSerbatoi(t, x, A, a, k, gamma, g, u_real);
    [~, xx] = ode45(dxdt, [0 Ts], x_real);

    % Aggiorna stato centrato
    x_log(:, tt+1) = xx(end, :)' - x_ref;
end

%% Plot risultati MPC nel dominio reale
% Traiettorie nel piano (h1, h3)
CIS_13_shifted = CIS_G_H13 + x_ref([1,3]);
Np_step_13_shifted = projection(Np_steps_set, [1,3]) + x_ref([1,3]);

x1_real = x_log(1,:) + x_ref(1);
x3_real = x_log(3,:) + x_ref(3);

figure()
subplot(1,2,1)
plot(Np_step_13_shifted, 'Alpha', 0, 'LineWidth', 2); hold on
plot(CIS_13_shifted);
plot(x1_real, x3_real, 'Color', [0 0 0.5]);
scatter(x1_real, x3_real, 'cyan');

xlabel('$h_1$ [cm]', 'Interpreter', 'latex');
ylabel('$h_3$ [cm]', 'Interpreter', 'latex');
title('Traiettoria con vincolo di uguaglianza: $(h_1, h_3)$', 'Interpreter', 'latex');
legend('N-step set', 'CIS', 'Traiettoria', 'Sample', 'Interpreter', 'latex');

% Traiettorie nel piano (h2, h4)
CIS_24_shifted = CIS_G_H24 + x_ref([2,4]);
Np_step_24_shifted = projection(Np_steps_set, [2,4]) + x_ref([2,4]);

x2_real = x_log(2,:) + x_ref(2);
x4_real = x_log(4,:) + x_ref(4);

subplot(1,2,2)
plot(Np_step_24_shifted, 'Alpha', 0, 'LineWidth', 2); hold on
plot(CIS_24_shifted);
plot(x2_real, x4_real, 'Color', [0 0 0.5]);
scatter(x2_real, x4_real, 'cyan');

xlabel('$h_2$ [cm]', 'Interpreter', 'latex');
ylabel('$h_4$ [cm]', 'Interpreter', 'latex');
title('Traiettoria con vincolo di uguaglianza: $(h_2, h_4)$', 'Interpreter', 'latex');
legend('N-step set', 'CIS', 'Traiettoria', 'Sample', 'Interpreter', 'latex');