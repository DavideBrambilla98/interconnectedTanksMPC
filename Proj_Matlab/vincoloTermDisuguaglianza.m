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
Q = 100*eye(4); % Penalizza lo stato (quanto gli stati h1,h2,h3,h4 devono essere vicino al riferimento)
R = eye(2); % Penalizza l'ingresso (quanto limitare l'uso degli ingressi v1 e v2)
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

%% MPC
T_sim = 60;
mpc = MPC(sys_d.A,sys_d.B,Hx,hx,Hu,hu,CIS_H,CIS_h,x_ref,u_ref,Q,R,N);

    
%% Verifica dei vincoli

disp('--- Verifica vincoli MPC ---');
disp('Max valore b_ineq_base:');
disp(max(mpc.b_ineq_base));
disp('Min valore b_ineq_base:');
disp(min(mpc.b_ineq_base));

disp('Dimensioni A_ineq:');
disp(size(mpc.A_ineq));
disp('Dimensioni b_ineq_base:');
disp(size(mpc.b_ineq_base));


%vincoli centrati
disp('hx originali:');
disp(hx);
disp('hx centrati:');
disp(hx - Hx*x_ref);

disp('hu originali:');
disp(hu);
disp('hu centrati:');
disp(hu - Hu*u_ref);

%%
% Log stati e ingressi del sistema
x_log = zeros(4, T_sim+1); % Matrice per salvare gli stati del sistema
u_log = zeros(2, T_sim); % Matrice per salvare gli ingressi di controllo
flags = zeros(1, T_sim); % Per salvare lo stato di uscita dell'ottimizzatore(quadprog)

x_log(:, 1) = x0_centrato; % Stato iniziale del sistema

for tt = 1:T_sim
    % Stato del sistema linearizzato
    x_current = x_log(:, tt);
    x_lin = x_current;
    x_lin_shifted = x_lin;

    % Impostare i vincoli MPC in base alla condizione iniziale
    f = mpc.f_base * x_lin_shifted;
    b_ineq = mpc.b_ineq_base - mpc.b_ineq_x0_factor*x_lin_shifted;

    % Risoluzione del problema di ottimizzazione
    [delta_u_seq, ~, exitflag] = quadprog(mpc.F, f, mpc.A_ineq, b_ineq);
    flags(tt) = exitflag;
    
    if isempty(delta_u_seq) || exitflag <= 0
        warning("Quadprog non ha trovato una soluzione ammissibile al passo %d (exitflag = %d)", tt, exitflag);
        delta_u_seq_first = zeros(2,1); % oppure: break; % per fermare la simulazione
    else
        delta_u_seq_first = delta_u_seq(1:2); % usa solo se esiste una soluzione
    end

    % Avanzamento/Risposta del sistema
    u_log(:,tt) = u_ref + delta_u_seq_first;

    dxdt = @(t,x) livSerbatoi(t, x, A, a, k, gamma, g, u_log(:,tt)); % u_log(:,tt) per passare entrambi gli ingressi
    [~, xx] = ode45(dxdt, [0 Ts], x_current);

    x_log(:, tt+1) = xx(end, :)';
end

%% Plot risultati MPC

%   Traslazione del CIS e dell'N-step set nelle coordinate originali
CIS_13_shifted = CIS_G_H13 + x_ref;
Np_step_set_13_shifted = h_nsteps_13 + x_ref;

figure(3)
h_npstep_13_shifted = Np_step_set_13_shifted.plot('Alpha', 0, 'LineWidth', 2);
hold on
h_cis_13_shifted = CIS_13_shifted.plot();
h_traj_13 = plot(x_log(1, :), x_log(2, :), 'Color',[0 0 0.5]);
h_traj_dots_13 = scatter(x_log(1,:),x_log(2,:),'cyan');
title('Traiettoria del sistema')
xlabel("$h_1$" , Interpreter="latex")
ylabel("$h_3$" , Interpreter="latex")

legend([h_cis_13_shifted, h_npstep_13_shifted, h_traj_13, h_traj_dots_13], ...
    {'CIS', sprintf('%d-step set', N), 'Traiettoria', 'Sample'},...
    'Interpreter','latex')
%% Andamento degli stati e degli ingressi

% Andamento degli stati
figure;
subplot(2,1,1);
plot(0:T_sim, x_log');
title('Andamento degli stati');
xlabel('Tempo [step]');
ylabel('Stato');
legend('x_1','x_2','x_3','x_4');

% Andamento degli ingressi
subplot(2,1,2);
plot(1:T_sim, u_log');
title('Ingressi di controllo');
xlabel('Tempo [step]');
ylabel('Ingresso');
legend('u_1','u_2');

