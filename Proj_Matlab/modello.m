%% Setup iniziale
clear; clc; close all;
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 14)
set(0,'DefaultFigureWindowStyle', 'docked') 
set(0,'defaulttextInterpreter','latex')
rng('default');

%% Parametri fisici del sistema
g = 981;                         % [cm/s^2] gravità
A = [28, 32, 28, 32];            % [cm^2] sezione dei serbatoi
a = [0.071, 0.057, 0.071, 0.057];% [cm^2] sezione fori
k = [2.7, 3.2];                  % [cm^3/(s·V)] costanti delle pompe
gamma = [0.3, 0.4];             % suddivisione flussi

% Tensione nominale (esempio, non usata direttamente)
u_nom = [0.375, 0.3]';          

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
dxdt = @(t, x) livSerbatoi(t, x, A, a, k, gamma, g, u_nom); % NB: usa u_nom fissa

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

Ts = 1; % [s]
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

% Vincoli centrati
x_max_lin = x_max - x_ref;
x_min_lin = x_min - x_ref;
u_max_lin = u_max - u_ref;
u_min_lin = u_min - u_ref;

% Matrici vincoli: Hx*x <= hx
Hx = [eye(4); -eye(4)];
hx = [x_max_lin; -x_min_lin];

Hu = [eye(2); -eye(2)];
hu = [u_max_lin; -u_min_lin];
