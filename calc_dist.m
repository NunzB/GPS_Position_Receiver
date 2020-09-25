% Calcula distancia entre duas coordenadas

function[Delta] = calc_dist(sub_const_matrix,coordinates,sat_dim_pret)

for k=1:sat_dim_pret
    dist(k) = sqrt((coordinates(1)-sub_const_matrix(k,1))^2 + (coordinates(2)-sub_const_matrix(k,2))^2 + (coordinates(3)-sub_const_matrix(k,3))^2);
end
Delta = dist;

end