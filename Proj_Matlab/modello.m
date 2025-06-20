%% Setup iniziale
clear; clc; close all;
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
set(0,'defaulttextInterpreter','latex')
rng('default');

%richiamo funzioni
addpath('funzioni');

%% Parametri fisici del sistema
g = 981; % [cm/s^2] gravità
A = [28, 32, 28, 32]; % [cm^2] sezione dei serbatoi
a = [0.071, 0.057, 0.071, 0.057]; % [cm^2] sezione fori
k = [2.7, 3.2]; % [cm^3/(s·V)] costanti delle pompe
gamma = [0.3, 0.4]; % suddivisione flussi

% Tensione nominale (esempio, non usata direttamente)
u_start = [0,0]';          

%% Stato iniziale e di riferimento
x_start = [1.3767, 2.2772, 0.8386, 0.5604]';         % stato iniziale
x_ref   = [7.8253, 18.7323, 3.3545, 7.8801]';        % stato di equilibrio

% Calcolo delle tensioni in equilibrio
u2_ref = (a(3) * sqrt(2 * g * x_ref(3))) / ((1 - gamma(2)) * k(2));
u1_ref = (a(4) * sqrt(2 * g * x_ref(4))) / ((1 - gamma(1)) * k(1));
u_ref = [u1_ref; u2_ref];

disp('Tensioni di equilibrio u_ref:');
disp(u_ref);

%% Simulazione del sistema non lineare

sim_time = 60*2; % durata simulazione [s]

% ODE del sistema
dxdt = @(t, x) livSerbatoi(t, x, A, a, k, gamma, g, u_start);

% Integrazione numerica
[tt, xx] = ode45(dxdt, linspace(0, sim_time, sim_time+1), x_start);

% Conversione tempo in minuti
tt = tt / 60;

% Plot livelli nei serbatoi
figure;
sgtitle("Simulazione del Quadruple Tank Process") 

subplot(2,1,1)
plot(tt, xx(:, 1:4), 'LineWidth', 1.5); hold on;
yline(0.5, '--r', 'Min', 'LabelVerticalAlignment','bottom');
yline(20, '--r', 'Max', 'LabelVerticalAlignment','top');
ylim([0 30])
title("Livelli d'acqua nei serbatoi")
ylabel("Livello [cm]")
xlabel("Tempo [min]")
legend(["h1", "h2", "h3", "h4"], 'Location', 'best')
grid on

%% Linearizzazione simbolica

syms h1 h2 h3 h4 u1 u2 real

% Equazioni dinamiche
h1_dot = (-a(1)*sqrt(2*g*h1) + a(3)*sqrt(2*g*h3) + gamma(1)*k(1)*u1)/A(1);
h2_dot = (-a(2)*sqrt(2*g*h2) + a(4)*sqrt(2*g*h4) + gamma(2)*k(2)*u2)/A(2);
h3_dot = (-a(3)*sqrt(2*g*h3) + (1-gamma(2))*k(2)*u2)/A(3);
h4_dot = (-a(4)*sqrt(2*g*h4) + (1-gamma(1))*k(1)*u1)/A(4);

F = [h1_dot; h2_dot; h3_dot; h4_dot];
stati = [h1, h2, h3, h4];
ingressi = [u1, u2];

% Jacobiane
A_sym = jacobian(F, stati);
B_sym = jacobian(F, ingressi);

% Valutazione in x_ref e u_ref
A_lin = double(subs(A_sym, [h1 h2 h3 h4], x_ref.'));
B_lin = double(subs(B_sym, [u1 u2], u_ref.'));

C_lin = eye(4);
D_lin = zeros(4,2);
sys_lineare = ss(A_lin, B_lin, C_lin, D_lin);

x0_centrato = x_start - x_ref;

%% Stabilità

disp('--- Stabilità del sistema linearizzato ---');
disp('Autovalori A (continuo):');
disp(eig(A_lin));

%% Discretizzazione

Ts = 20; % [s]
sys_d = c2d(sys_lineare, Ts);

disp('--- Stabilità del sistema discretizzato ---');
disp('Moduli autovalori A_d:');
disp(abs(eig(sys_d.A)));

%% Raggiungibilità

disp('--- Raggiungibilità ---');
Mr_c = ctrb(sys_lineare);
Mr_d = ctrb(sys_d);

fprintf('Continuo: rango %d, dimensione %dx%d\n', rank(Mr_c), size(Mr_c));
fprintf('Discreto: rango %d, dimensione %dx%d\n', rank(Mr_d), size(Mr_d));

%% Vincoli di stato e ingresso

% Limiti assoluti
X_bounds = [0.5, 20];    % livelli serbatoi [cm]
U_bounds = [0.0, 4.5];   % tensioni pompe [V]

% Costruzione vettori max/min
x_max = X_bounds(2) * ones(4,1);
x_min = X_bounds(1) * ones(4,1);
u_max = U_bounds(2) * ones(2,1);
u_min = U_bounds(1) * ones(2,1);


h12_bounds = X_bounds;  % [min, max] per h1 e h2
h34_bounds = X_bounds;  % [min, max] per h3 e h4

% Creazione vincoli: max (prime righe), min (ultime righe)
X_vinc = [
    h12_bounds(2)*ones(2,1);  % h1 e h2 max
    h34_bounds(2)*ones(2,1);  % h3 e h4 max
    h12_bounds(1)*ones(2,1);  % h1 e h2 min
    h34_bounds(1)*ones(2,1)   % h3 e h4 min
];

% Vincoli sugli ingressi (u1, u2)
U_vinc = [
    U_bounds(2)*ones(2,1);  % u1 e u2 max
    U_bounds(1)*ones(2,1)   % u1 e u2 min
];

% Calcoliamo i vincoli centrati nel punto di equilibrio
X_vinc_lin = X_vinc - repmat(x_ref, 2, 1);  % duplica x_ref
U_vinc_lin = U_vinc - repmat(u_ref, 2, 1);  % duplica u_ref

% Definizione delle matrici H per i vincoli Hx*x <= hx
Hx = [eye(4); -eye(4)];  % 8 vincoli (4 max + 4 min)
hx = [ones(4,1); -ones(4,1)] .* X_vinc_lin;

Hu = [eye(2); -eye(2)];  % 4 vincoli (2 max + 2 min)
hu = [ones(2,1); -ones(2,1)] .* U_vinc_lin;


