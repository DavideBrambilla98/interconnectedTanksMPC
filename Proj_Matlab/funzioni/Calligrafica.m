function [A_cal , A_cal_n ,B_cal , B_cal_n , Q_cal , R_cal] = Calligrafica(A , B , Q , R , S , N)

    %inizzializzazione matrici
    A_cal = [];
    B_cal = [];
    
    for i = 0:N
        mat = A^i;
        % Costruzione di A_cal da A_cal_0 a A_cal_n (stato dopo n passi)
        A_cal = [A_cal;mat];
        if i == N
            A_cal_n = mat;
        end
    end
    % Calcolo di B_cal
    % come input influenza stato a ciascun passo
    % B_cal_n = come tutti gli input influenzano lo stato al passo n
    for i = 0:N
        riga_B_cal = [];
        
        for j = 1:N
            if (i-j) < 0
                riga_B_cal = [riga_B_cal , zeros(height(B) , width(B))];
            else
                riga_B_cal= [riga_B_cal , A^(i-j) * B];
            end
        end
        B_cal = [B_cal ; riga_B_cal];
    
        if i == N
            B_cal_n = riga_B_cal;
        end

    end

    %Costruzione R_cal
    R_cal = kron(eye(N), R); %kron ripete R lungo la diagonale N volte  

    %Costruzione Q_cal
    Q_cal = diagonale(Q, S, N+1);
end


%creo una matrice composta da p sulla diagonale eccetto N-esimo = s
%i restanti elementi esterni alla diagonale = 0
function mat = diagonale(p, s, N)
    matrici = cell(1, N);
  
    % Inserisce la matrice p nelle prime N-1 posizioni
    for k = 1 : N-1
        matrici{k} = p;
    end
    
    % Inserisce la matrice s nell'ultima posizione(N)
    matrici{N} = s;
    
    % Costruisce la matrice diagonale a blocchi
    mat = blkdiag(matrici{:});
end
