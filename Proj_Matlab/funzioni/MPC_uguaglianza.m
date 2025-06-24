function mpc = MPC_uguaglianza(A, B, Hx, hx, Hu, hu, x_ref, u_ref, Q, R, Np)

% Dimensioni
n = size(A, 2);
m = size(B, 2);

[~, P] = dlqr(A, B, Q, R);

% Matrici estese 
[A_cal , A_cal_n ,B_cal , B_cal_n , Q_cal , R_cal] = Calligrafica(A ,B , Q , R , P , Np);

% Matrice hessiana del costo
F = (B_cal' * Q_cal * B_cal + R_cal);
F = (F + F') / 2;

% Componente lineare del costo
f = B_cal' * Q_cal * A_cal;

% Vincoli di disuguaglianza
Hx_shifted = Hx; hx_shifted = hx;
Hu_shifted = Hu; hu_shifted = hu;

Hx_tilde = kron(eye(Np+1), Hx_shifted);
hx_tilde = repmat(hx_shifted, [Np+1, 1]);

Hu_tilde = kron(eye(Np), Hu_shifted);
hu_tilde = repmat(hu_shifted, [Np, 1]);

A_ineq = [Hx_tilde * B_cal; Hu_tilde];
b_ineq = [hx_tilde; hu_tilde];

% Vincolo di uguaglianza terminale: x_N = 0
A_eq = B_cal_n;
b_eq = -A_cal_n;

% Output struct
mpc.F = F;
mpc.f = f;
mpc.A_ineq = A_ineq;
mpc.b_ineq = b_ineq;
mpc.A_eq = A_eq;
mpc.b_eq = b_eq;
mpc.Np = Np;

% Extra per runtime
mpc.Hx_tilde = Hx_tilde;
mpc.Hu_tilde = Hu_tilde;
mpc.A_cal = A_cal;
mpc.A_cal_n = A_cal_n;
mpc.B_cal = B_cal;
mpc.B_cal_n = B_cal_n;

end

