function mpc = MPC(A,B,Hx,hx,Hu,hu,CIS_H,CIS_h,x_ref,u_ref,Q,R,Np)

%dimensioni
n = size(A,2);      %numero stati
m = size(B,2);      %numero ingressi
n_ter = length(CIS_h);   %numero di righe del vincolo terminale

%matrice costo terminale
[~, P] = dlqr(A, B, Q, R);

% Traslare i vincoli rispetto al riferimento
% Vincoli sullo stato traslato
Hx_shifted = Hx;
hx_shifted = hx - Hx*x_ref;

% Vincoli per l'ingresso
Hu_shifted = Hu;
hu_shifted = hu - Hu*u_ref;

%richiamo a funzione function [A_cal , A_cal_n ,B_cal , B_cal_n , Q_cal , R_cal] = Calligrafica(A , B , Q , R , riccati , N)

[A_cal , A_cal_n ,B_cal , B_cal_n , Q_cal , R_cal] = Calligrafica(A ,B , Q , R , P , Np)

%   Matrice hessiana costo quadratico
F = B_cal'*Q_cal*B_cal + R_cal;

%   Componente lineare costo quadratico
f = B_cal' * Q_cal * A_cal;

% Vincoli
Hx_tilde = kron(eye(Np+1), Hx_shifted);
hx_tilde = repmat(hx_shifted, [Np+1, 1]);

Hx_tilde = [Hx_tilde; zeros(n_ter, Np*n), CIS_H];
hx_tilde = [hx_tilde; CIS_h];

Hu_tilde = kron(eye(Np), Hu_shifted);
hu_tilde = repmat(hu, [Np, 1]);

% Set ammissibili di ingresso (inequalities)
A_ineq = [Hx_tilde*B_cal; Hu_tilde];
b_ineq = [hx_tilde; hu_tilde];

%   Creazione struttura mpc
mpc.F = F;
mpc.f_base = f;
mpc.A_ineq = A_ineq;
mpc.b_ineq_base = b_ineq;
mpc.Np = Np;
mpc.b_ineq_x0_factor = [Hx_tilde*A_cal; zeros(2*m*Np, n)];
end