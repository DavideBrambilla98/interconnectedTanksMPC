function [H_nsteps, h_nsteps] = controllable_set(Hx, hx, Hu, hu, H_target, h_target, A, B, N)
% Calcola il controllable set in N passi verso un insieme target
% Input:
%   Hu*u <= hu  → vincoli sull'ingresso
%   Hx*x <= hx  → vincoli sullo stato
%   H_target*x <= h_target  → regione target
%   A, B → dinamica del sistema

n = size(A,2);

% target iniziale
H_ii_steps = H_target;
h_ii_steps = h_target;

for ii=1:N
    %set ad un passo rispetto a quello preceente
    temp = Polyhedron('A',[H_ii_steps*A H_ii_steps*B; ...
        zeros(size(Hu,1),n) Hu],'b',[h_ii_steps; hu]);
    %   Proiezione in R^n
    temp = projection(temp,1:n);
    temp.minHRep();
    %   Intersezione con X := {x | Hx*x <= hx}
    H_ii_steps = [temp.A; Hx];
    h_ii_steps = [temp.b; hx];
    
end

H_nsteps = H_ii_steps;
h_nsteps = h_ii_steps;

end


