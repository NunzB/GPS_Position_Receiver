%% read_almanaque.m
function [id,square_A,eccent,toa,alfa,rate_of_right_ascend,RAAN,arg_perigee,mean_anom,af0,af1] = read_almanaque (A)
    
    %Variáveis auxiliares
    n = 0;
    k = 1;
    i = 1;
    verif = 1;

    [n,m] = size(A.data);

    r = (n+1)/14;

    %Inicialização dos parâmetros contidos no almanaque
    id = (1/r);
    eccent = zeros(1,r);
    toa = zeros(1,r);
    alfa = zeros(1,r);
    rate_of_right_ascend = zeros(1,r);
    square_A = zeros(1,r);
    RAAN = zeros(1,r);
    arg_perigee = zeros(1,r);
    mean_anom = zeros(1,r);
    af0 = zeros(1,r);
    af1 = zeros(1,r);

    while(i<n)
            id(k) = A.data(i);
            i = i+2;
            eccent(k) = A.data(i);
            i=i+1;
            toa(k) = A.data(i);
            i=i+1;
            alfa(k) = A.data(i);
            i=i+1;
            rate_of_right_ascend(k) = A.data(i);
            i=i+1;
            square_A(k) = A.data(i);
            i=i+1;
            RAAN(k) = A.data(i);
            i=i+1;
            arg_perigee(k) = A.data(i);
            i=i+1;
            mean_anom(k) = A.data(i);
            i=i+1;
            af0(k) = A.data(i);
            i=i+1;
            af1(k) = A.data(i);
            i=i+3;

            k = k+1; 
    end
    
end


%% Calc_sat_pos.m

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


%% Sat_view.m

% Determina quais os satelites que estao visiveis, ou seja, os que tem maior 
%elevacao que o angulo de mascara

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


%% Sat_sub.m

function[sub_const_matrix] = sat_sub(dim_sub_const,sat_matrix,id_min,r)

        sub_const_matrix = zeros(dim_sub_const,3);
    
        for i=1:dim_sub_const
        for k=1:r
        if sat_matrix(k,4) == id_min(i)
           sub_const_matrix(i,1) = sat_matrix(k,1) ;
           sub_const_matrix(i,2) = sat_matrix(k,2);
           sub_const_matrix(i,3) = sat_matrix(k,3);    
        end
        end
        end
         
end

%% least_sq_alg.m

function coord_finais = least_sq_alg(sub_const_matrix,est_ini,coord_recep,sat_dim_pret,var_ruido)
      
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

  coord_finais= est_ini;
end


%% pdop_min.m

% Calcula a constelacao de satelites que tem a menor PDOP
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

%% pseudoranges_calc.m

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


%% path.m

function [P]=path()

    f=1; % Frequência a que as posições são actualizadas (1Hz)
    x_0=0;
    y_0=0;
    z_0=0;
    v_0 = 50; %Velocidade inicial (50m/s)
    m=10;     %Angulo de máscara
    P=[x_0 y_0 z_0 v_0 0 m]; % Inicialização do percurso
    
    %Percurso AB
    t=100; %duração do percurso em segundos
    a=1; %aceleração constante

    for i = 1:t
        x=x_0;
        y=y_0+(a*i^2)/2;
        z=z_0;
        v=v_0+a*i;
        P=[P;x y z v a m];
    end
    
   % Percuso BC
    t=50; % duração do percurso em segundos
    a=0;  % Velocidade constante a partir de B
    w=(3*pi/2)/t;   % velocidade angular
    r=v/w;          % raio do percuso
    c = [x, y+r];
    for i = 1:(t/f)
        x = c(1)-r*(1-cos(i*w));
        y = c(2)+r*sin(i*w);
        P=[P;x y z v a m];
    end
    % Percurso CD
    t=50; %duração do percurso em segundos
    x_c=x;
        for i = 1:t
        x=x_c+v*i;
        P=[P;x y z v a m];
    end
    % Percurso DE
    t=50; %duração do percurso em segundos
    x_d=x;
    m=28;
    for i = 1:t
        x=x_d+v*i;
        P=[P;x y z v a m];
    end
    
end


%% elev_azi.m

% CONVERTE DE ENU PARA AZIMUTE E ELEVAÇÃO
function sat_elev_azi=elev_azi(e,n,u)

for i =1:size(e(:,1))

    azim(i,1)=rad2deg(atan2(e(i,1),n(i,1)));
    elev(i,1)=rad2deg(asin(u(i,1)/sqrt(e(i,1)^2+n(i,1)^2+u(i,1)^2)));
    
end
%el=rad2deg(atan2(u,sqrt(n^2+e^2)));
sat_elev_azi = [azim, elev];
end


%% DMStoD.m

% Converte graus, minutos e segundos para graus.
function [D] = DMStoD (DMS)

    D = zeros(size(DMS,1),1);

    for i=1:size(DMS,1)
        
        for j=1:size(DMS,2)

            D(i) = D(i) + DMS(i,j) / 60^(j-1);
    
        end
    end

end


%% Calc_EA.m

function[E] = Calc_EA(M,e)

    tol = 10^-8;
    Etemp = M;
    ratio = 1;

    while abs(ratio) > tol
        f_E = Etemp - e*sin(Etemp) - M;
        f_Eprime = 1 - e*cos(Etemp);
        ratio = f_E/f_Eprime;
        if abs(ratio) > tol
            Etemp = Etemp - ratio;
        else
            E = Etemp;
        end
    end
    
end


%% calc_dist.m

% Calcula distancia entre duas coordenadas
function[Delta] = calc_dist(sub_const_matrix,coordinates,sat_dim_pret)

for k=1:sat_dim_pret
    dist(k) = sqrt((coordinates(1)-sub_const_matrix(k,1))^2 + (coordinates(2)-sub_const_matrix(k,2))^2 + (coordinates(3)-sub_const_matrix(k,3))^2);
end
Delta = dist;

end


%% isValidhour.m

function res = isValidhour( hour )

if (length(hour)~=8)
    res=0;
    return;
end

HHstr = hour(1:2);
HH=str2double(HHstr);

if (HH<0 || HH>23 || isnan(HH)==1)
    res=0;
    return;
end

MMstr = hour(4:5);
MM=str2double(MMstr);
if (isempty(MM)==1 || (60<MM)&& (MM<0) || isnan(MM)==1)
   res=0;
   return;
end

SSstr = hour(7:8);
SS=str2double(SSstr);
if (isempty(SS)==1 || (60<SS)&& (SS<0) || isnan(SS)==1)
   res=0;
   return;
end

if (hour(3)~=':'||hour(6)~=':')
    res=0;
    return;
end
res=1;
end


%% isValidDate.m

function res=isValidDate(d)
% d = Enter date in 'DD/MM/YYYY' format
if (length(d)~=10)
    res=0;
    return;
end
yearstr=d(7:10);
year=str2double(yearstr);
if (year<1993 || isnan(year)==1)
    res=0;
    return;
end
monthstr=d(4:5);
month=str2double(monthstr);
if (isempty(month)==1 || isnan(month)==1)
   res=0;
   return;
elseif (month<1||month>12)
      res=0;
      return;
end  
if (d(3)~='/'||d(6)~='/')
    res=0;
    return;
end
daystr=d(1:2);
day=str2double(daystr);
if (isempty(day)==1 || isnan(day)==1)
    res=0;
    return;
end
if (month==1||month==3||month==5||month==7||month==8||month==10||month==12)
    if (day<1||day>31)
        res=0;
        return;
    end
end
if (month==4||month==6||month==9||month==11)
    if (day<1||day>30)
        res=0;
        return;
    end
end

if (rem(year,400)==0)
    if (month==2)
        if (day<1||day>29)
            res=0;
            return;
        end
    end
else if (rem(year,100)==0)
        if (month==2)
            if (day<1||day>28)
                res=0;
                return;
            end
        end
else if (rem(year,4)==0)
        if (month==2)
            if (day<1||day>29)
                res=0;
                return;
            end
        end
    else
         if (month==2)
            if (day<1||day>28)
                res=0;
                return;
            end
         end
    end
    end
end

if day < 1 && month<7 && year < 1993
    res=0;
    return
end

res=1 ;
end


%% isValidcoordinates.m

function [test,coord_graus] = isValidcoordinates(coordinates)

if (length(coordinates)~=10)
    test=0;
    coord_graus=0;
    return;
end

if (coordinates(3)~='-'||coordinates(6)~='-'||coordinates(9)~='-')
    test=0;
    coord_graus=0;
    return;
end

dd = str2double(coordinates(1:2));
mm = str2double(coordinates(4:5));
ss = str2double(coordinates(7:8));
N_S= coordinates(10);

if N_S == 'N' || N_S =='E'
    coord_graus = DMStoD ([dd mm ss]);
    test=1;
    return
else if N_S == 'S' || N_S =='W'
    coord_graus = -DMStoD ([dd mm ss]);   
    test=1;
    return
    end
end

if (abs(coord_graus)>180)
if (abs(coord_graus)>90)
    test=0;
    coord_graus=0;
    return;
end
end
end
