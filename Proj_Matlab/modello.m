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
u = [0, 0];

% Definizione obiettivi di controllo
x_start = [1.3767, 2.2772, 0.8386, 0.5604]';    %cond. iniziale
x_ref = [7.8253, 18.7323, 3.3545, 7.8801]';     %equilibrio

% Calcolo delle tensioni di equilibrio
u2_ref = (a(3) * sqrt(2 * g * x_ref(3))) / ((1 - gamma(2)) * k(2));
u1_ref = (a(4) * sqrt(2 * g * x_ref(4))) / ((1 - gamma(1)) * k(1));

u_ref = [u1_ref, u2_ref]';

disp('Tensioni di equilibrio u_ref:');
disp(u_ref);

%% Simulazione del sistema

simulazione = 60*20; % durata simulazione in secondi

% Definizione dell'ODE
dxdt = @(t, x) livSerbatoi(t, x, A, a, k, gamma, g, u);

% Simulazione del comportamento del sistema
[tt, xx] = ode45(dxdt, linspace(0, simulazione, simulazione+1), x_start);

% Conversione tempo in minuti
tt = tt / 60;

% Plot simulazione (traiettorie)
figure
sgtitle("Simulazione del Quadruple Tank Process con controllo tensione alle pompe") 

% Livelli nei serbatoi
subplot(2,1,1)
plot(tt, xx(:, 1:4), 'LineWidth', 1.5);
hold on
yline(0.5, '--r', 'Min', 'LabelVerticalAlignment','bottom');
yline(5, '--r', 'Max', 'LabelVerticalAlignment','top');
title("Livelli d'acqua nei serbatoi")
ylabel("Livello $[cm]$", 'Interpreter', 'latex')
xlabel("Tempo $[min]$", 'Interpreter', 'latex')
legend(["h1", "h2", "h3", "h4"], 'Location', 'best')
grid on

%% Linearizzazione

% Variabili di stato e ingresso simboliche
syms h1 h2 h3 h4 u1 u2 real

h1_dot = (-a(1)*sqrt(2*g*h1) + a(3)*sqrt(2*g*h3) + gamma(1)*k(1)*u1)/A(1);
h2_dot = (-a(2)*sqrt(2*g*h2) + a(4)*sqrt(2*g*h4) + gamma(2)*k(2)*u2)/A(2);
h3_dot = (-a(3)*sqrt(2*g*h3) + (1-gamma(2))*k(2)*u2)/A(3);
h4_dot = (-a(4)*sqrt(2*g*h4) + (1-gamma(1))*k(1)*u1)/A(4);

F = [h1_dot; h2_dot; h3_dot; h4_dot];
stati = [h1 h2 h3 h4];
ingressi = [u1 u2];

% Jacobiane
A_sim(h1,h2,h3,h4) = jacobian(F, stati);
B_sim(u1,u2) = jacobian(F, ingressi);

% Valutazione nell'equilibrio
A_lin = double(A_sim(x_ref(1) , x_ref(2) , x_ref(3) ,x_ref(4)));
B_lin = double(B_sim(u_ref(1) , u_ref(2)));

C_lin = eye(4);
D_lin = zeros(4,2);

sys_lineare = ss(A_lin, B_lin, C_lin, D_lin);

x0_centrato = x_start - x_ref;

%% Verifica della Stabilità del sistema lineare

fprintf('\n--- Stabilità del sistema linearizzato ---\n');
fprintf('Autovalori della matrice A (continua):\n');
disp(eig(A_lin));

%% Discretizzazione

Ts = 1; % tempo di campionamento [s]
sys_discretizzato = c2d(sys_lineare, Ts);

fprintf('\n--- Stabilità del sistema discretizzato ---\n');
fprintf('Moduli degli autovalori della matrice A (discretizzata):\n');
disp(abs(eig(sys_discretizzato.A)));

%% Analisi della Raggiungibilità

fprintf('\n--- Analisi della Raggiungibilità ---\n');

Mr_lineare = ctrb(sys_lineare);
Mr_discretizzato = ctrb(sys_discretizzato);

fprintf('\nSistema Linearizzato (Continuo):\n');
fprintf('  - Rango: %d\n', rank(Mr_lineare));
fprintf('  - Dimensioni matrice: %d righe x %d colonne\n', size(Mr_lineare, 1), size(Mr_lineare, 2));

fprintf('\nSistema Discretizzato:\n');
fprintf('  - Rango: %d\n', rank(Mr_discretizzato));
fprintf('  - Dimensioni matrice: %d righe x %d colonne\n', size(Mr_discretizzato, 1), size(Mr_discretizzato, 2));


%% Vincoli sullo stato e sugli ingressi
U_vinc = [0, 4.5];   %vincolo sulle tensioni delle due pompe [V]
H_vinc = [0.5, 20];  %vincoli su livelli dell'acqua nei serbatoi [cm]

% vettori dei vincoli concatenati max/min

X_vinc = [ H_vinc(2) * ones(4,1); H_vinc(1) * ones(4,1)];
U_vinc = [U_vinc(2) * ones(2,1); U_vinc(1) * ones(2,1)]; 

% vincoli centrati nel punto di equilibrio

X_vinc_lin = X_vinc - [x_ref; x_ref];
U_vinc_lin = U_vinc - [u_ref; u_ref];

% Matrici dei vincoli: Hx * x <= hx e Hu * u <= hu

Hx = [eye(4); -eye(4)];
hx = X_vinc_lin;
Hu = [eye(2); -eye(2)];
hu = U_vinc_lin;
