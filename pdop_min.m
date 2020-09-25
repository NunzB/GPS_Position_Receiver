% FUNCAO - Calcula a constelacao de satelites que tem a menor PDOP

% RECEBE - Coordenadas satelites visiveis, numero de satelites visiveis, coordenadas do receptor e a dimensao da constelacao pretendida

% RETORNA - O valor do PDOP minimo e os ID's dos satelites

function[min_pdop, id_sat_min_pdop] = pdop_min(sat_vis, nr_sat, coordinates, dim_sub_const)


    % Inicializacao
    r = zeros(1,nr_sat);
    
    aux = nchoosek(1:nr_sat,dim_sub_const); % Corresponde a todas as combinacoes possiveis de satelites
    max = nchoosek(nr_sat,dim_sub_const); % Corresponde ao numero maximo de combinacoes
    
    G = zeros(dim_sub_const, 4);
    pdop = zeros(1,max);
    id_sat_min_pdop = zeros(1,dim_sub_const);
    
    for j=1:nr_sat
        % Calcula as distancias receptor - satelite cada combinacao
        r(j) = sqrt((coordinates(1)-sat_vis(j,1))^2 + (coordinates(2)-sat_vis(j,2))^2 + (coordinates(3)-sat_vis(j,3))^2);
    end
    
    for i=1:max
        for n=1:dim_sub_const
            % Matriz G para cada combinacao satelites
            G(n,1) = (coordinates(1)-sat_vis(aux(i,n),1))/r(aux(i,1));
            G(n,2) = (coordinates(2)-sat_vis(aux(i,n),2))/r(aux(i,2));
            G(n,3) = (coordinates(3)-sat_vis(aux(i,n),3))/r(aux(i,3));
            G(n,4) = 1;
        end
        
        % Calcula o PDOP da combinacao
        H = inv(G.'*G);
        pdop(i) = sqrt(H(1,1) + H(2,2) + H(3,3));
    end  
    
    % Encontra o PDOP minimo e a posicao na matriz
    [min_pdop, pos] = min(pdop(:));
    
    
    for k=1:dim_sub_const
        % Determina os IDs atraves da posicao na matriz onde estava o PDOP
        % minimo, e vai buscar a matriz das coordenadas dos satelites (que tambem continha os IDs)
        id_sat_min_pdop(k)= sat_vis(aux(pos,k),4);
    end
    
end