function [H_nsteps,h_nsteps] = controllable_set(Hx,hx,Hu,hu,H_target,h_target,A,B,N)
%CONTROLLABLE_SET Calcolo dell'N-step-controllable set ad un set target
%descritto da H_target * x <= h_target

n = size(A,2);
m = size(B,2);

%   Candidato iniziale
H_ii_steps = H_target;
h_ii_steps = h_target;

for ii=1:N
    
fprintf('Iterazione %d\n', ii);
    %   Computazione in R^(n+m)
    temp = Polyhedron('A',[H_ii_steps*A H_ii_steps*B; ...
        zeros(size(Hu,1),n) Hu],'b',[h_ii_steps; hu]);
    %   Proiezione in R^n
    temp = projection(temp,1:n);
    temp.minHRep();
    %   Intersezione con X := {x | Hx*x <= hx}
    H_ii_steps = [temp.A; Hx];
    h_ii_steps = [temp.b; hx];

    fprintf('Vincoli nel %d-step set: %d\n', N, size(H_ii_steps,1));
end

H_nsteps = H_ii_steps;
h_nsteps = h_ii_steps;

end
