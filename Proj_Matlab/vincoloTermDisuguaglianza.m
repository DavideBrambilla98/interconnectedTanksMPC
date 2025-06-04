clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
% set(0,'DefaultFigureWindowStyle', 'normal')
set(0,'defaulttextInterpreter','latex')
rng('default');

%% Invariant set per il sistema controllato

% Tempo di campionamento
Ts = 1; %[s]

% Richiamo del modello del sistema dei serbatoio interconnessi
modello

% Matrici del costo quadratico
Q = eye(4); % Penalizza lo stato (4 stati: h1,h2,h3,h4)
R = eye(2); % Penalizza l'ingresso (2 ingressi: v1,v2)

% Control invariant set CIS_H*x <= CIS_h
[CIS_H, CIS_h] = cis(sys_discretizzato.A, sys_discretizzato.B, x_ref, u_ref, Hx, hx, Hu, hu, Q, R);
CIS = Polyhedron(CIS_H, CIS_h);


% Costruzione del CIS
CIS_G = Polyhedron(CIS_H, CIS_h);
CIS_G.minHRep(); % Riduzione alla rappresentazione minimale

% Proiezione su h1 e h2
CIS_G_h12 = CIS_G.projection([1 2]);

% Proiezione su h3 e h4
CIS_G_h34 = CIS_G.projection([3 4]);

% Plot delle proiezioni
figure;

subplot(1, 2, 1);
plot(CIS_G_h12);
title("Proiezione del CIS su h1 e h2");
xlabel("h1 [cm]", 'Interpreter', 'latex');
ylabel("h2 [cm]", 'Interpreter', 'latex');
xlim([-10 10])
ylim([-10 10])
grid on;

subplot(1, 2, 2);
plot(CIS_G_h34);
title("Proiezione del CIS su h3 e h4");
xlabel("h3 [cm]", 'Interpreter', 'latex');
ylabel("h4 [cm]", 'Interpreter', 'latex');
xlim([-10 10])
ylim([-10 10])
grid on;


V = CIS_G_h12.V;
disp('Vertici della proiezione su h1 e h2:');
disp(V);





