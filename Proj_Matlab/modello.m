%%
clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
set(0,'defaulttextInterpreter','latex')
rng('default');

%% Parametri del sistema

g = 9.81; % [m/s^2] accelerazione di gravità

% Sezione del serbatoio
A = [28, 32, 28, 32]; % [cm^2]

% Sezione del foro del serbatoio
a = [0.071, 0.057, 0.071, 0.057]; % [cm^2]

% Portata della pompa su Volt (flusso d'acqua generato dalla pompa per ogni volt)
k = [2.7, 3.2]; % [cm^3/(s*V)]

% Partizione flusso nei serbatoi
gamma = [0.3, 0.4]; 

% Tensione applicata alle pompe
v = [0, 0];

% Definizione obiettivi di controllo
x_start = [1.3767, 2.2772, 0.8386, 0.5604]'; % cond. iniziale
x_ref = [7.8253, 18.7323, 3.3545, 7.8801]';
u_ref = [0, 0]';

%% Equazioni dinamiche del sistema
% Linearizzazione
% Discretizzazione

%% Invariant Set

%% Progettazione MPC


