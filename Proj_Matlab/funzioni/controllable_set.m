% Calcola il controllable set in N passi verso un insieme target N>1 si blocca
% function [H_nsteps, h_nsteps] = controllable_set(Hx, hx, Hu, hu, H_target, h_target, A, B, N)

function [H_nsteps, h_nsteps, Np] = controllable_set(Hx, hx, Hu, hu, H_target, h_target, A, B, start)
% Calcola l'N-step controllable set verso un insieme target.
% Restituisce la rappresentazione H-rep del set e l'orizzonte minimo Np tale che lo stato iniziale è controllabile.

    n = size(A, 2); % Dimensione dello stato
    m = size(B, 2); % Dimensione dell'ingresso

    % Parametri di configurazione
    N_max = 50;             % Limite massimo di iterazioni
    tol = 1e-6;             % Tolleranza numerica

    % Inizializzazione con il set target
    H_ii_steps = H_target;
    h_ii_steps = h_target;

    Np = NaN; % Inizializza Np come non trovato

    for ii = 1:N_max
        % Costruzione del poliedro vincolato nel passo corrente
        P = Polyhedron('A', H_ii_steps, 'b', h_ii_steps);

        % Applica la dinamica inversa
        temp = Polyhedron('A', [H_target*A, H_target*B; zeros(size(Hu,1),n), Hu], ...
                          'b', [h_ii_steps; hu]);

        % Proiezione sullo spazio degli stati x ∈ R^n
        temp = projection(temp, 1:n);
        temp = temp.minHRep(); 

        % Interseca con X := {x | Hx*x <= hx}
        H_ii_steps = [temp.A; Hx];
        h_ii_steps = [temp.b; hx];

        % Verifica se lo stato iniziale è contenuto nel set corrente
        current_set = Polyhedron('A', H_ii_steps, 'b', h_ii_steps);
        if current_set.contains(start, tol)
            Np = ii;
            break;
        end
    end

    % Assegna comunque un valore alle uscite
    if isnan(Np)
        H_nsteps = [];
        h_nsteps = [];
    else
        final_set = Polyhedron(H_ii_steps, h_ii_steps);
        final_set = final_set.minHRep();
        H_nsteps = final_set.A;
        h_nsteps = final_set.b;
    end

end