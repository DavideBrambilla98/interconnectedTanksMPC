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
CIS = Polyhedron(CIS_H, CIS_h);
%
%% Plot del CIS
CIS_G = Polyhedron(CIS_H, CIS_h);
CIS_G = CIS_G.minHRep();
CIS_G_H12 = projection(CIS_G, 1:2);
CIS_G_H34 = projection(CIS_G, 3:4);

dim = CIS_G.Dim;
disp(['La dimensione del poliedro è: ', num2str(dim)]);

% Plot delle proiezioni
figure

subplot(1 , 2 , 1)
CIS_G_H12.plot();
title("Proiezione del CIS dei livelli dei serbatoi 1 e 2")
xlim([-20 20])
ylim([-20 20])
grid on
axis equal
xlabel("$h_1$" , Interpreter="latex")
ylabel("$h_2$" , Interpreter="latex")

subplot(1 , 2 , 2)
CIS_G_H34.plot();
title("Proiezione del CIS dei livelli dei serbatoi 3 e 4")
xlim([-20 20])
ylim([-20 20])
grid on
axis equal
xlabel("$h_3$" , Interpreter="latex")
ylabel("$h_4$" , Interpreter="latex")

%% N-step controllable set

% Orizzonte di predizione fissato
%N = 1;
%fprintf('\n--- Calcolo del %d-step controllable set ---\n', N);
%[Np_steps_H, Np_steps_h] = controllable_set(Hx, hx, Hu, hu, CIS_H, CIS_h, sys_d.A, sys_d.B, N);
%fprintf('Vincoli nel %d-step set: %d\n', N, size(Np_steps_H,1));

[Np_steps_H, Np_steps_h , Np] = controllable_set(Hx, hx, Hu, hu, CIS_H, CIS_h, sys_d.A, sys_d.B, x0_centrato);

% Costruzione poliedro
Np_steps_set = Polyhedron(Np_steps_H, Np_steps_h);
Np_steps_set = Np_steps_set.minHRep();

%% plot


%% MPC e simulazione
