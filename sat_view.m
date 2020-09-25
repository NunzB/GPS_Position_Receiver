% FUNCAO - determina quais os satelites que estao visiveis, ou seja, os que tem maior elevacao que o angulo de mascara

% RECEBE - Coordenadas ECEF, elevacao e azimuth dos satelites e o angulo mascara

% RETORNA - Numero satelites visiveis, elevacao e azimuth dos satelites
% visiveis e coordenadas ECEF,

function [nr_sat_mask,sat_mask_elev_azi,sat_mask_ecef] = sat_view(sat_matrix, sat_elev_azi, mask)

i=1;
for k=1:size(sat_elev_azi,1)
    if sat_elev_azi(k,1)>= mask
        sat_mask_elev_azi(i,:) = sat_elev_azi(k,:);
        sat_mask_ecef(i,:) = sat_matrix(k,:);
        i=i+1;
    end
end

nr_sat_mask=size(sat_mask_elev_azi,1);

end
