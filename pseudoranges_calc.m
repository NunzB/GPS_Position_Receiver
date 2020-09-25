function[pseu_meas] = pseudoranges_calc(sub_const_matrix,coord_recep,sat_dim_wanted,iono,var_ruido)
pseu_meas(sat_dim_wanted)=0;
    for k=1:sat_dim_wanted
        pseu_meas(k) = sqrt((coord_recep(1)-sub_const_matrix(k,1))^2 + (coord_recep(2)-sub_const_matrix(k,2))^2 + (coord_recep(3)-sub_const_matrix(k,3))^2) + coord_recep(4);

        LLA = ecef2lla ([coord_recep(1) coord_recep(2) coord_recep(3)]);
        
        [sat_E,sat_N,sat_U] = ecef2enu(sub_const_matrix(k,1),sub_const_matrix(k,2),sub_const_matrix(k,3),LLA(1),LLA(2),LLA(3),referenceEllipsoid('wgs84'));
        [sat_elev_azi] = elev_azi(sat_E,sat_N,sat_U);

        elev = sat_elev_azi(1,1);
        noise=var_ruido*randn();     %Ruído
        pseu_meas(k) = pseu_meas(k) + noise;
        if iono == 1
            pseu_meas(k) = pseu_meas(k) + 10/sind(elev);
        end     
    end
end
