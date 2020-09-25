 % SATELLITE VIEW %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 % FUNCAO determina quais os satelites que estao visiveis, ou se ja , os que tem maior elevacao que oangulo de mascara

 % RECEBE Coordenadas ECEF , e l e v a c a o e azimuth dos s a t e l i t e s e o angulo mascara


 % RETORNA Numero s a t e l i t e s v i s i v e i s , coordenada s ECEF , e l e v a c a o e azimuth dos s a t e l i t e s v i s i v e i s

 function [nr_sat_mask , sat_mask_elev_azi , sat_mask_ecef] = satellite_view(sat_matrix , sat_elev_azi,mask_angle)

 i =1;

 for k=1: size(sat_elev_azi ,1)
 if sat_elev_azi(k,1)>= mask_angle
 sat_mask_elev_azi(i,:) = sat_elev_azi(k,:);
 sat_mask_ecef(i,:) = sat_matrix(k,:);
 i=i+1;
 end
 end

 nr_sat_mask=size(sat_mask_elev_azi , 1);

 end