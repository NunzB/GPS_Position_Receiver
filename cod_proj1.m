clc;

%Pr� aloca��o de vari�veis
erro_pos=zeros(251,1);
erro_quadratico=zeros(100,251);
est_enu=zeros(251,3);
vel_est=zeros(251,1);
erro_vel=zeros(252,1);
erro_quad_vel=zeros(100,251);

mse_AB=zeros(1,100);
mse_BC=zeros(1,100);
mse_CD=zeros(1,100);
mse_DE=zeros(1,100);
mse_vel=zeros(1,100);

%% Inicializa��o de vari�veis
    c = 3e8; 
    ano = 0;
    mes = 0;
    dia = 0;
    hora = 0;
    min = 0;
    sec = 0;
    mask = 10;
%% Abre o ficheiro que cont�m o almanaque
    afile = fopen('almanaque.txt','r');

    if afile == -1
        disp('Ficheiro n�o encontrado')
        disp(' ')
        return
    end

    A = importdata('almanaque.txt',':');
    [n,m] = size(A.data);
    numero_sat = (n+1)/14; %N�mero de sat�lites

    [id,square_A,e0,toa,alfa,rate_of_right_ascend,RAAN,arg_perigeu,mean_anom,af0,af1] = read_almanaque (A);
    
%% Data e hora de transmiss�o 
fprintf('\n-------------------------\n\n');
fprintf('Insira a data exacta de GPS superior a 6 de Janeiro de 2008:\n\n');
test=0;

while (test ~= 1)

    data_s = input('Data em formato (DD/MM/YYYY): ', 's');
    test = isValidDate(data_s); %d = Enter date in 'DD/MM/YYYY' format
    
    if test == 0
        disp('Introduza novamente a Data no formato referido.')
    end                              
       
end
test=0;
while (test ~= 1)
    disp('');
    disp('00 <Hour< 23  00<Minutes<60 0<Seconds<60');
    hour_s = input('Hora em formato (HH:MM:SS): ', 's');
    test = isValidhour(hour_s); 
    if test == 0
        disp('Introduza novamente a hora no formato referido')
    end    
end
test=0;
while (test ~= 1)
    disp('');
    disp('introduza um angulo de mascara minimo de referencia');
    mask_s = input('Angulo de mascara em graus:', 's');
    mask = str2double(mask_s);
    if isnan(mask) == 1 
        disp('Introduza novamente uma valor valido')
    else
        test = 1;
    end
end

ano= str2num(data_s(7:end));            %#ok<ST2NM>
mes= str2num(data_s(4:5));              %#ok<ST2NM>
dia= str2num (data_s (1:2));            %#ok<ST2NM>
hora= str2num (hour_s(1:2));            %#ok<ST2NM>
min= str2num (hour_s(4:5));             %#ok<ST2NM>
sec= str2num (hour_s(7:8));             %#ok<ST2NM>

t_st1 = 24*60*60*datenum(ano, mes, dia, hora, min, sec); %Tempo inicial de transmiss�o
t_sim = 0;


%% C�lculo da posi��o inicial dos sat�lites em ECEF
sat_matrix = calc_sat_pos(square_A,mean_anom,e0,arg_perigeu,RAAN,rate_of_right_ascend,alfa,t_sim,t_st1, numero_sat,id);
 
%% Coordenadas exatas do recetor em ECEF

% Coordenadas do recetor
% latitude = 40;
% longitude = -9;
% altitude = 2000;

fprintf('\n-------------------------');
fprintf('\n\nInsira a posi��o exacta do receptor:\n\n');
test=0;
while (test ~= 1)
        lat_recp_str = input('Latitude formato: DD-MM-SS-(N or S): ', 's');
        [test, latitude] = isValidcoordinates(lat_recp_str);
end
test=0;
while (test ~= 1)
        long_recp_str = input('Longitude formato: DD-MM-SS-(W or E): ', 's');
        [test, longitude] = isValidcoordinates(long_recp_str);
end
test=1;
while (test ~= 0)
        altitude_str = input('Altitude do receptor: ', 's');
        altitude = str2double(altitude_str);
        test = isnan(altitude) + imag(altitude);
end
test=1;

fprintf('Latitude (Graus):          %0.2f�\n', latitude);
fprintf('Longitude (Graus):         %0.2f�\n', longitude);
fprintf('Altitude (Graus):          %0.1f m\n', altitude);
fprintf('\nPrima uma tecla para continuar\n');
pause;

coordinates_receptor_e = lla2ecef([latitude longitude altitude]); 
coord_recep = (coordinates_receptor_e);


%% C�lculo da eleva��o e do azimute dos sat�lites em rela��o � posi��o inicial do recetor

 [sat_E , sat_N , sat_U] = ecef2enu(sat_matrix(:,1),sat_matrix(:,2),sat_matrix(:,3),coord_recep(1), coord_recep(2) , coord_recep(3) , referenceEllipsoid('wgs84'));
 sat_elev_azii = elev_azi(sat_E , sat_N , sat_U);

%% calculo dos satelites visiveis a partir da posi��o do receptor

[nr_sat_mask, sat_mask_elev_azii, sat_mask_ecef] = sat_view(sat_matrix, sat_elev_azii, mask);

% C�lculo dos ID's dos sat�lites que minimizam o PDOP
% sat_dim_pret=4; % Dimensao pretendida da sub constela��o

    verif=1;
    while(verif == 1)
        verif = 0;
        disp('');
        fprintf('Existem %i sat�lites vis�veis.', nr_sat_mask);
        disp('Quantos sat�lites pretende utilizar? (M�nimo 4)');

        sat_dim_pret = input('','s');
        sat_dim_pret = str2double(sat_dim_pret);
        disp('');

        if isnan(sat_dim_pret) == 1 || (sat_dim_pret > nr_sat_mask && sat_dim_pret < 4)
            disp ('Valor inv�lido. Insira um n� de sat�lites v�lido.');
            verif = 1;
        end
    end
    
[min_pdop , id_sat_min_pdop] = pdop_min(sat_mask_ecef, nr_sat_mask, coord_recep, sat_dim_pret);

%Constela��o dos sat�lites que minimizam o PDOP de dimens�o 4
[sub_const_matrix] = sat_sub(sat_dim_pret, sat_mask_ecef, id_sat_min_pdop, nr_sat_mask);
sub_const_matrix(:,4) = id_sat_min_pdop(:,4);

% Elevacao e azimute dos sat escolhidos em rela��o ao receptor
 [sat_EE, sat_NN ,sat_UU] = ecef2enu(sub_const_matrix(:,1), sub_const_matrix(:,2), sub_const_matrix(:,3), coord_recep(1), coord_recep(2), coord_recep(3), referenceEllipsoid('wgs84'));
 [sub_const_matrix_elev_azi ] = elev_azi(sat_EE, sat_NN, sat_UU);
 
 
 %% Medi��o das pseudo-dist�ncias
 
 coord_recep(4) = 0; % 4�a vari�vel � o tempo
 
 %Sele��o da vari�ncia do ru�do 
 
    verif=1;
    while(verif == 1)
        verif=0;
        disp('Indique a vari�ncia pretendida do ru�do na medi��o das pseudo-distancias:');
        var_ruido = input('','s');
        var_ruido = str2double(var_ruido);

        if isnan(var_ruido) == 1
             disp ('Valor inv�lido. Insira um valor correcto');
             verif = 1;
        end
    end
    
    verif=1;
    while(verif == 1)
        verif = 0;
        disp('Pretende considerar os erros ionosf�ricos?');
        disp('1 - Sim');
        disp('0 - N�o');

        iono = input('','s');
        iono = str2double(iono);
        disp('');

        if isnan(iono) == 1 || (iono ~= 0 && iono ~=1)
            disp ('Valor inv�lido. Insira um valor correcto');
            verif = 1;
        end
    end

%C�lculo das pseudo-dist�ncias
[pseu_meas] = pseudoranges_calc(sub_const_matrix,coord_recep,sat_dim_pret,iono,var_ruido);


%% Estima��o inicial do recetor

%Estimativa inicial utilizada - Centro geogr�fico de portugal
lat_cg_pt=37;
long_cg_pt=7;
alt_cg_pt=500;
est_ini=lla2ecef([lat_cg_pt long_cg_pt alt_cg_pt]);
est_ini(4) = c*4;

%Estima��o inicial do recetor
est_rec_ini = least_sq_alg(sub_const_matrix,est_ini,coord_recep,sat_dim_pret,var_ruido);
est_rec = est_rec_ini;

%% Percurso percorrido

%C�lculo da posi��o do avi�o ao longo do percurso
percurso = path();

for z=1:100    %100 Simula��es
est_enu(size(percurso,1), 3)=0;   %Matriz para as estimativas
erro_pos(size(percurso,1),3)=0; %Matriz para os erros de posi��o

for k=1:size(percurso,1)

    %Coordenadas EC EF e LLA do percurso realizado
    [p_ecef(1), p_ecef(2), p_ecef(3)] = enu2ecef(percurso(k,1), percurso(k,2), percurso(k,3), latitude, longitude, altitude, referenceEllipsoid('wgs84'));
    p_lla = ecef2lla([p_ecef(1) p_ecef(2) p_ecef(3)]);
    
    %Posi��o dos sat�lites em ECEF para cada instante de simula��o
    t_sim=k;
    sat_matrix = calc_sat_pos(square_A,mean_anom,e0,arg_perigeu,RAAN,rate_of_right_ascend,alfa,t_sim,t_st1, numero_sat,id);
    
    %Eleva��o e Azimute dos sat�lites para cada posi��o do recetor
    [sat_E , sat_N , sat_U] = ecef2enu(sat_matrix(:,1),sat_matrix(:,2),sat_matrix(:,3), p_lla(1), p_lla(2) , p_lla(3) , referenceEllipsoid('wgs84'));
    sat_elev_azi = elev_azi(sat_E , sat_N , sat_U);
   %
    %C�lculo dos sat�lites vis�veis em cada instante com a restri��o da
    %m�scara
    [nr_sat_mask, sat_mask_elev_azi, sat_mask_ecef] = sat_view(sat_matrix, sat_elev_azi, percurso(k,6));
    
    %Check number of satelites available
    if nr_sat_mask >= sat_dim_pret
        [min_pdop , id_sat_min_pdop] = pdop_min( sat_mask_ecef, nr_sat_mask, p_ecef, sat_dim_pret);
    
        %Coordenadas ECEF dos sat que minimizam o PDOP
        [sub_const_matrix] = sat_sub(sat_dim_pret, sat_mask_ecef, id_sat_min_pdop, nr_sat_mask);
    else
        disp('N� de sat�lites visiveis insuficientes. O programa terminou.')
        break
    end
    
    %Introduz o ID dos sat�lites na 4� coluna da matriz
    sub_const_matrix(:,4) = id_sat_min_pdop(:,4);
    
    p_ecef(4)=est_rec_ini(4);
    %Mede as pseudo-dist�ncias entre o recetor e os sat�lites escolhidos
    [pseu_meas] = pseudoranges_calc(sub_const_matrix, p_ecef, sat_dim_pret, iono, var_ruido);
    
    %C�lculo da posi��o do recetor atrav�s do alg. de Minimos Quadrados
    %Garante que se for o percurso DE que n� de sat = 4
    if k>201
        est_rec = least_sq_alg(sub_const_matrix, est_rec, p_ecef, 4, var_ruido);
    else
        est_rec = least_sq_alg(sub_const_matrix, est_rec, p_ecef, sat_dim_pret, var_ruido);
    end
    
    %Erro absoluto e quadr�tico cometido com a estimativa
    erro_pos(k)=sqrt((est_rec(1)-p_ecef(1))^2+(est_rec(2)-p_ecef(2))^2+(est_rec(3)-p_ecef(3))^2);
    erro_quadratico(z,k) = erro_pos(k)^2;
    
    mse_AB(z)=mean(erro_quadratico(z,1:101));
    mse_BC(z)=mean(erro_quadratico(z,101:151));
    mse_CD(z)=mean(erro_quadratico(z,151:201));
    mse_DE(z)=mean(erro_quadratico(z,201:251));
    
    %Calculo das coord ENU das posi��es estimadas
    [est_E, est_N, est_U]=ecef2enu(est_rec(1), est_rec(2), est_rec(3), latitude, longitude, altitude, referenceEllipsoid('wgs84'));
    est_enu(k,:)=[est_E est_N est_U];   
    est_enu(252,:)=[0 0 0];
    
    vel_est(k)=sqrt((est_enu(k,1)-est_enu(k+1,1))^2+(est_enu(k,2)-est_enu(k+1,2))^2);
    erro_vel(k)=sqrt((vel_est(k)-percurso(k,4))^2);
    erro_quad_vel(z,k)=erro_vel(k)^2;
    
    mse_vel(z)=mean(erro_quad_vel(z,:));
end

 est_enu=est_enu(1:end-1,:);
 vel_est=vel_est(1:end-1);
 vel_est(101)=vel_est(100);
 erro_vel=erro_vel(1:end-1);
 erro_vel(101)=erro_vel(100);

end
