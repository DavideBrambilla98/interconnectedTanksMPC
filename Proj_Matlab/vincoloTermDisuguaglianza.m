clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
set(0,'defaulttextInterpreter','latex')
rng('default');

%% Invariant set per il sistema controllato

% Richiamo del modello del sistema dei serbatoio interconnessi
addpath('funzioni');
modello;

% Matrici del costo quadratico
Q = diag([1, 1, 100, 100]);  % Penalizza lo stato (quanto gli stati h1,h2,h3,h4 devono essere vicino al riferimento)
R = 1*eye(2);                 % Penalizza l'ingresso (quanto limitare l'uso degli ingressi v1 e v2)
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
N = 3;
fprintf('\n--- Calcolo del %d-step controllable set ---\n', N);
[Np_steps_H, Np_steps_h] = controllable_set(Hx, hx, Hu, hu, CIS_H, CIS_h, sys_d.A, sys_d.B, N);
%fprintf('Vincoli nel %d-step set: %d\n', N, size(Np_steps_H,1));

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

%% MPC
T_sim = 120;
mpc = MPC(sys_d.A,sys_d.B,Hx,hx,Hu,hu,CIS_H,CIS_h,x_ref,u_ref,Q,R,N);

%% Simulazione MPC con sistema centrato e dinamica non lineare reale

x_log = zeros(4, T_sim+1); % Stati centrati
u_log = zeros(2, T_sim);   % Ingressi reali
flags = zeros(1, T_sim);   % Exitflag di quadprog

x_log(:, 1) = x_start - x_ref; % Stato iniziale centrato

for tt = 1:T_sim
    % Stato attuale centrato
    x_centrato = x_log(:, tt);

    % Calcolo del termine lineare del costo
    f = real(mpc.f_base * x_centrato);
    b_ineq = real(mpc.b_ineq_base - mpc.b_ineq_x0_factor * x_centrato);

    % Risoluzione del problema QP
    [delta_u_seq, ~, exitflag] = quadprog(mpc.F, f, mpc.A_ineq, b_ineq);
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

% Traslazione del CIS e dell'N-step set nelle coordinate originali
CIS_13_shifted = CIS_G_H13 + x_ref([1,3]);
h_nsteps_13 = projection(Np_steps_set, [1,3]);
Np_step_set_13_shifted = h_nsteps_13 + x_ref([1,3]);

% Traiettoria reale nel piano (h1, h3)
x1_real = x_log(1,:) + x_ref(1);
x3_real = x_log(3,:) + x_ref(3);

figure()
subplot(1 , 2 , 1)
h_npstep_13_shifted = Np_step_set_13_shifted.plot('Alpha', 0, 'LineWidth', 2);
hold on
h_cis_13_shifted = CIS_13_shifted.plot();
h_traj_13 = plot(x1_real, x3_real, 'Color',[0 0 0.5]);
h_traj_dots_13 = scatter(x1_real, x3_real, 'cyan');

title('Traiettoria del sistema nel piano $(h_1, h_3)$', 'Interpreter', 'latex')
xlabel("$h_1$ [cm]", 'Interpreter','latex')
ylabel("$h_3$ [cm]", 'Interpreter','latex')

legend([h_cis_13_shifted, h_npstep_13_shifted, h_traj_13, h_traj_dots_13], ...
    {'CIS', sprintf('%d-step set', N), 'Traiettoria', 'Sample'}, ...
    'Interpreter','latex')

CIS_24_shifted = CIS_G_H24 + x_ref([2,4]);
h_nsteps_24 = projection(Np_steps_set, [2,4]);
Np_step_set_24_shifted = h_nsteps_24 + x_ref([2,4]);

% Traiettoria reale nel piano (h2, h4)
x2_real = x_log(2,:) + x_ref(2);
x4_real = x_log(4,:) + x_ref(4);

subplot(1 , 2 , 2)
h_npstep_24_shifted = Np_step_set_24_shifted.plot('Alpha', 0, 'LineWidth', 2);
hold on
h_cis_24_shifted = CIS_24_shifted.plot();
h_traj_24 = plot(x2_real, x4_real, 'Color',[0 0 0.5]);
h_traj_dots_24 = scatter(x2_real, x4_real, 'cyan');

title('Traiettoria del sistema nel piano $(h_2, h_4)$', 'Interpreter', 'latex')
xlabel("$h_2$ [cm]", 'Interpreter','latex')
ylabel("$h_4$ [cm]", 'Interpreter','latex')

legend([h_cis_24_shifted, h_npstep_24_shifted, h_traj_24, h_traj_dots_24], ...
    {'CIS', sprintf('%d-step set', N), 'Traiettoria', 'Sample'}, ...
    'Interpreter','latex')

%% Andamento degli stati e degli ingressi

% Andamento degli stati
figure;
subplot(2,1,1);
hold on;

% Plot degli stati con riferimento sommato
plot((0:T_sim)*Ts/60, x_log' + x_ref');

% Aggiunta delle linee di riferimento
xline = (0:T_sim)*Ts/60;
for i = 1:4
    plot(xline, ones(size(xline)) * x_ref(i), '--', 'DisplayName', sprintf('x_%d ref', i));
end

title('Andamento degli stati');
xlabel('Tempo [min]');
ylabel('Stato');
legend('x_1','x_2','x_3','x_4','x_1 ref','x_2 ref','x_3 ref','x_4 ref');
grid on;


% Andamento degli ingressi
subplot(2,1,2);
plot((1:T_sim)*Ts/60, u_log');
title('Ingressi di controllo');
xlabel('Tempo [min]');
ylabel('Ingresso');
legend('u_1','u_2');
grid on
