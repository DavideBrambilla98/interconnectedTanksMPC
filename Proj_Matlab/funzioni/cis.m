function [G,g] = cis(A,B,x_ref,u_ref,Fx,fx,Fu,fu,Q,R)
%  CIS Calcolo del control invariant set (CIS) di un sistema lineare
%   Questo metodo assume che un controllore LQR venga utilizzato
%   all'interno del CIS
%   Input:
%       - A, B: matrici del sistema
%       - x_ref: equilibrio attorno al quale costruire il CIS
%       - Fx*x <= fx: vincoli sullo stato
%       - Fu*u <= fu: vincoli sull'ingresso
%       - Q , R: matrici per LQR

%   Controllore LQR
K = -dlqr(A,B,Q,R);

%   Matrice A del sistema controllato con LQR
A_lqr = A+B*K;

%   Vincoli sullo stato e sull'ingresso
F = [Fx; Fu*K];
f = [fx; fu + Fu*(K*x_ref - u_ref)];

%   Calcolo del CIS (G*x <= g)
CIS_poly_prev = Polyhedron();
CIS_poly_curr = Polyhedron(F,f);

i = 0;
while CIS_poly_prev.isEmptySet || CIS_poly_prev ~= CIS_poly_curr
    i = i+1;
    %   Memorizza vecchio candidato
    CIS_poly_prev = CIS_poly_curr;
    
    %   Calcola nuovo candidato
    G_hat = [CIS_poly_curr.A * A_lqr; F];
    g_hat = [CIS_poly_prev.b + CIS_poly_curr.A*B*(K*x_ref - u_ref); f];
    CIS_poly_curr = Polyhedron(G_hat,g_hat);
fprintf('Iterazione %d\n', i);

end

%   Disequazioni che descrivono il CIS
G = CIS_poly_curr.A;
g = CIS_poly_curr.b;