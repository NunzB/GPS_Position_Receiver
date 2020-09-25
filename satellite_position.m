%% FUNCAO Ca l c u l a a p o s i c a o dos s a t e l i t e s em coordenada s ECEF
 % RECEBE Dados dos almanaques e tempos de t r a n smi s s a o i n i c i a l do s i n a l e tempo de s imu l a c a o

 % RETORNA Coordenadas dos S a t e l i t e s ECEF

 function[sat_matrix] = satellite_position(square_A, mean_anom , e0 ,w,RAAN, rate_of_right_ascend,alfa, t_sim , t_st1 , r , id)

 % Tempo de transmissao sinal
 t_st = t_st1 + t_sim;
 % Tempo de referencia de Ephemeris
 t_oe = 24*60*60* datenum(2000 , 1 , 23 , 0 , 0 , 0);
 % Constante gravidade Terra
 gravity_cte=3.986005e14 ;

 delta_t = t_st - t_oe ;

 % Correccao ao delta_t
 if delta_t > 302000
 delta_t = delta_t  604800;
 end

 % Inicializacao da matriz posicao satelites
 sat_matrix = zeros(r,3) ;
 for k=1:r
 % Caso particular quando a excentricidade da orbita e nula
 if e0==0
 teta0 = mean_anom(k)+w(k) ;
 teta = teta0 + n(k) * d e l t a t ;
 else
 n(k) = sqrt(gravity_cte / square_A(k)^6) ; % Mean Motion n
 mean_anom_aux = mean_anom(k) + n(k)*delta_t; % Mean anomaly M
 eccent_anom = kepler(mean_anom_aux, e0(k)) ; % E c c e n t r i c anomaly E
 true_anom = atan2(sqrt(1 e0(k)^2)*sin(eccent_anom) , cos(eccent_anom) e0(k)); % True anomaly v
 teta = true_anom + w(k); % Argument o f l a t i t u d e t h e t t a
 end

 %Longitude of ascending node at time t_st
 omega = RAAN( k ) + rate_of_right_ascend(k)*delta_t  2*pi/86164*t_st;

 % Orbit Radius
 orbit_rad = square_A(k)^2*(1 e0(k)*cos(eccent_anom));

 sat_matrix(k,1) = orbit_rad*cos(teta)*cos(omega) orbit_rad*sin(teta)*sin(omega)*cos(alfa(k));
 sat_matrix(k,2) = orbit_rad*cos(teta)*sin(omega)+orbit_rad*sin(teta)*cos(omega)*cos(alfa(k));
 sat_matrix(k,3) = orbit_rad*sin(teta)*sin(alfa(k));
 sat_matrix(k,4) = id(1,k) ;
 end
 end
