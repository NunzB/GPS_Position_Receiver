%% CONVERTE DE ENU PARA AZIMUTE E ELEVAÇÃO
% enu em qualquer unidade (desde que sempre a mesma)
% Devolve em GRAUS

function sat_elev_azi=elev_azi(e,n,u)

for i =1:size(e(:,1))

    azim(i,1)=rad2deg(atan2(e(i,1),n(i,1)));
    elev(i,1)=rad2deg(asin(u(i,1)/sqrt(e(i,1)^2+n(i,1)^2+u(i,1)^2)));
    
end
%el=rad2deg(atan2(u,sqrt(n^2+e^2)));
sat_elev_azi = [azim, elev];
end