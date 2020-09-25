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

