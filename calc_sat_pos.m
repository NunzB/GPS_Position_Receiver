%Função que dá a posição actual dos satélites
function[sat_matrix] = calc_sat_pos(square_A,mean_anom,e0,arg_perigeu,RAAN,rate_of_right_ascend,alfa,t_sim,t_st1,r,id)

    t_st = t_st1 + t_sim;      %Tempo de transmissão do sinal
    t_oe = 24*60*60*datenum(2000, 1, 23, 0, 0, 0); %Tempo de referencia da ephemeris
    cons_grav = 3.986005e14; %ctg de gravidade da terra

    delta_t = t_st-t_oe; 
    %Correção do delta t
    if delta_t > 302000 
        delta_t = delta_t - 604800;
    end
    
    %Matriz onde vai ser inserida a posição dos satélites
    sat_matrix = zeros(r,3);   

    %cálculo da constelação de satélites - coordenadas ECEF
    for k=1:r
 
        
        n(k) = sqrt(cons_grav/square_A(k)^6);      %Mean motion  
        mean_anom_aux = mean_anom(k) + n(k)*delta_t;  %Anomalia Média   
        eccent_anom = Calc_EA(mean_anom_aux,e0(k));  %Excentricidade
    %Anomalia verdadeira
        true_anom = atan2(sqrt(1-e0(k)^2)*sin(eccent_anom),cos(eccent_anom)-e0(k));   
        teta = true_anom + arg_perigeu(k); 
    %Longitude do nó ascendente
        omega = RAAN(k) + rate_of_right_ascend(k)*delta_t - 2*pi/86164*t_st;  
    %Raio da órbita
        orbit_rad = square_A(k)^2*(1-e0(k)*cos(eccent_anom));    

        sat_matrix(k,1) = orbit_rad*cos(teta)*cos(omega)-orbit_rad*sin(teta)*sin(omega)*cos(alfa(k));
        sat_matrix(k,2) = orbit_rad*cos(teta)*sin(omega)+orbit_rad*sin(teta)*cos(omega)*cos(alfa(k));
        sat_matrix(k,3) = orbit_rad*sin(teta)*sin(alfa(k));
        sat_matrix(k,4) = id(1,k);

    end

end