function coord_finais = least_sq_alg(sub_const_matrix,est_ini,coord_recep,sat_dim_pret,var_ruido)

  
%  for i=1:2
      
      pseu_meas = pseudoranges_calc(sub_const_matrix,coord_recep,sat_dim_pret,1,var_ruido);
      pseu_est = pseudoranges_calc(sub_const_matrix,est_ini,sat_dim_pret,0,var_ruido);
      dist = calc_dist(sub_const_matrix,est_ini,sat_dim_pret);

%Cálculo da matriz G
        for k=1:sat_dim_pret
            MG(k,1) = -(sub_const_matrix(k,1)-coord_recep(1))/dist(k);
            MG(k,2) = -(sub_const_matrix(k,2)-coord_recep(2))/dist(k);
            MG(k,3) = -(sub_const_matrix(k,3)-coord_recep(3))/dist(k);
            MG(k,4) = 1;
        end
      G = MG;
      
      Ginv = inv(G'*G);
      delta_pseu = pseu_meas - pseu_est;
      delta_pos = (Ginv * G') * delta_pseu';
      delta_pos_i = delta_pos';  
      est_ini = est_ini + delta_pos_i;
%  end
  coord_finais= est_ini;
  
end
