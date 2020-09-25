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

