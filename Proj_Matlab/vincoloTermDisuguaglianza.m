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
[CIS_H, CIS_h] = cis(sys_d.A, sys_d.B, zeros(4,1), zeros(2,1), Hx, hx, Hu, hu, Q, R);
CIS = Polyhedron(CIS_H, CIS_h);

%% Plot del CIS
CIS_G = Polyhedron(CIS_H, CIS_h);
CIS_G = minHRep(CIS_G);

CIS_G_H12 = projection(CIS_G, 1:2);
CIS_G_H34 = projection(CIS_G, 3:4);

figure
subplot(1 , 2 , 1)
CIS_G_H12.plot();
title("Proiezione del CIS dei livelli dei serbatoi 1 e 2")
xlim([-20 20])
ylim([-20 20])
xlabel("h1" , Interpreter="latex")
ylabel("h2" , Interpreter="latex")
