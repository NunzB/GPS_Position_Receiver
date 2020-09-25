%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             PLOTS               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure
polar(sat_elev_azii(:,2),sat_elev_azii(:,1),'o')
title('Azimute and Elevation of all the Satellites')
figure
polar(sat_mask_elev_azii(:,2),sat_mask_elev_azii(:,1),'o')
title('Azimute and Elevation of the Visible Satellites')
figure
polar(sub_const_matrix_elev_azi(:,2),sub_const_matrix_elev_azi(:,1),'o')
title('Azimute and Elevation of the Visible Satellites with PDOP Minimization')

% % Percurso Total
figure
subplot(2,2,[1 3])
plot(est_enu(:,1),est_enu(:,2),percurso(:,1),percurso(:,2))
title('Posi��es do receptor e da estimativa')
xlabel('East')
ylabel('North')
legend('Estimativas','Posi��o Exacta')
hold on
subplot(2,2,2)
plot(erro_pos(1:size(est_enu)))
title('Erro de Posi��o')
xlabel('Tempo de Simula��o')
ylabel('Erro (metro)')
hold on
subplot(2,2,4)
plot(erro_vel)
title('Erro de Velocidade')
xlabel('Tempo de Simula��o')
ylabel('Erro (metro/s)')

% Percurso AB
figure
subplot(2,1,1)
plot(mse_AB)
title('Erros m�dios quadr�ticos para o percurso AB')
xlabel('N�mero de simula��es')
ylabel('MSE (m^2)')
hold on
subplot(2,1,2)
plot(sqrt(mse_AB))
title('Raiz quadrada erros m�dios quadr�ticos para o percurso AB')
xlabel('N�mero de simula��es')
ylabel('RMSE (m)')


% Percurso BC
figure
subplot(2,1,1)
plot(mse_BC)
title('Erros m�dios quadr�ticos para o percurso BC')
xlabel('N�mero de simula��es')
ylabel('MSE (m^2)')
hold on
subplot(2,1,2)
plot(sqrt(mse_BC))
title('Raiz quadrada erros m�dios quadr�ticos para o percurso BC')
xlabel('N�mero de simula��es')
ylabel('RMSE (m)')


% Percurso CD
figure
subplot(2,1,1)
plot(mse_CD)
title('Erros m�dios quadr�ticos para o percurso CD')
xlabel('N�mero de simula��es')
ylabel('MSE (m^2)')
hold on
subplot(2,1,2)
plot(sqrt(mse_CD))
title('Raiz quadrada erros m�dios quadr�ticos para o percurso CD')
xlabel('N�mero de simula��es')
ylabel('RMSE (m)')

% Percurso DE
figure
subplot(2,1,1)
plot(mse_DE)
title('Erros m�dios quadr�ticos para o percurso DE')
xlabel('N�mero de simula��es')
ylabel('MSE (m^2)')
hold on
subplot(2,1,2)
plot(sqrt(mse_DE))
title('Raiz quadrada erros m�dios quadr�ticos para o percurso DE')
xlabel('N�mero de simula��es')
ylabel('RMSE (m)')


% Velocidade exacta e velocidade estimada
figure
plot(mean(erro_quadratico_v))
hold on
plot(percurso(:,4))
title('Velocidade do receptor e estimativa')
ylabel('Velocidade (m/s)')
xlabel('Tempo de Simula��o (s)')
legend('Estimativas','Velocidade Exacta')
