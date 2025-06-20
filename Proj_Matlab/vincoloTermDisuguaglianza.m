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
N = 5;
fprintf('\n--- Calcolo del %d-step controllable set ---\n', N);
[Np_steps_H, Np_steps_h] = controllable_set(Hx, hx, Hu, hu, CIS_H, CIS_h, sys_d.A, sys_d.B, N);
fprintf('Vincoli nel %d-step set: %d\n', N, size(Np_steps_H,1));

% Costruzione poliedro
Np_steps_set = Polyhedron(Np_steps_H, Np_steps_h);
Np_steps_set = Np_steps_set.minHRep();

% Plot N-step-controllable set
figure()
subplot(1 , 2 , 1)
h_cis_24 = CIS_G_H13.plot();
h_nsteps_24 = projection(Np_steps_set, [1,3]);
hold on
h_nsteps_24 = h_nsteps_24.plot('Alpha', 0, 'LineWidth', 2);
title(sprintf('CIS e %d-step controllable set del sistema linearizzato H1-H3: %d\n', N))
xlabel("$h_1$" , Interpreter="latex")
ylabel("$h_3$" , Interpreter="latex")
legend([h_cis_24, h_nsteps_24], ...
    {'CIS', sprintf('%d-step set', N)},...
    'Interpreter','latex')
subplot(1 , 2 , 2)
h_cis_24 = CIS_G_H24.plot();
h_nsteps_24 = projection(Np_steps_set, [2,4]);
hold on
h_nsteps_24 = h_nsteps_24.plot('Alpha', 0, 'LineWidth', 2);
title(sprintf('CIS e %d-step controllable set del sistema linearizzato H2-H4: %d\n', N))
xlabel("$h_2$" , Interpreter="latex")
ylabel("$h_4$" , Interpreter="latex")
legend([h_cis_24, h_nsteps_24], ...
    {'CIS', sprintf('%d-step set', N)},...
    'Interpreter','latex')
%% MPC

