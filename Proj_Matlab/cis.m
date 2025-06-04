function [G, g] = cis(A, B, x_ref, u_ref, Hx, hx, Hu, hu, Q, R)
%CIS Calcolo del control invariant set (CIS) di un sistema lineare
%   Input:
%       - A, B: matrici del sistema
%       - x_ref: equilibrio attorno al quale costruire il CIS
%       - u_ref: ingresso di equilibrio
%       - Hx*x <= hx: vincoli sullo stato
%       - Hu*u<= hu: vincoli sull'ingresso
%       - Q,R: matrici per LQR

%   Controllore LQR
K = -dlqr(A, B, Q, R);

% Matrice A del sistema controllato con LQR
A_lqr = A + B * K;

% Vincoli sullo stato e sull'ingresso (H*x <= h)
% x >> Hx <= hx
% u := Hu(K(x - x_ref) + u_ref) <= hu
H = [Hx; Hu * K];
h = [hx;hu + Hu*(K*x_ref - u_ref)]; 

%   Calcolo del CIS (G*x <= g)
CIS_poly_curr = Polyhedron(H, h);
primaIterazione = 1;
i = 0;
tic;
while  primaIterazione || CIS_poly_prev ~= CIS_poly_curr 
    i = i+1;

    primaIterazione = 0;

    %   Memorizza vecchio candidato
    CIS_poly_prev = Polyhedron(CIS_poly_curr.A , CIS_poly_curr.b);
    
    %   Calcola nuovo candidato (G_hat * x <= g_hat)
    % x_next := Ax + Bu = Ax + B[K(x - x_ref) + u_ref ]
    % >> Alqr * x - B(K * x_ref - u_ref)
    G_hat = [CIS_poly_curr.A * A_lqr;H];
    g_hat = [CIS_poly_curr.b + CIS_poly_curr.A * B * (K*x_ref - u_ref) ; h];
    CIS_poly_curr= Polyhedron(G_hat, g_hat);
    CIS_poly_curr= CIS_poly_curr.minHRep();

    if mod(i , 20) == 0
        disp("Iterazione numero: " + i + " tempo trascorso: " + toc);
    end
end
disp("Iterazioni:" + i)
disp("Tempo impiegato: " + toc);


%   Disequazioni che descrivono il CIS
G = CIS_poly_curr.A;
g = CIS_poly_curr.b;

end