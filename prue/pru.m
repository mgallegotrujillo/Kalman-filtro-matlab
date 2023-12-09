clear all
clc
duracion=10;
dt=0.1;
% kalman(duracion, dt)

% function kalman(duracion, dt)
%
% Simulación del filtro de Kalman para un vehículo que circula por una carretera.
% INPUTS
% duracion = duración de la simulación (seg)
% dt = step size (seconds)

%Sistema
a = [1 dt; 0 1]; % Matriz de transicion
b = [dt^2/2; dt]; % Matriz de entrada
c = [1 0]; % Matriz de sensores
x = [0; 0]; % Estado inicial
xposhat=x(1);
xvelhat=x(2);

%Filtro de kalman
xhat = x; % Estado inicial estimacion
DesviacionRuidomedicion = 10; % ruido de medición de posición (m)
DesviacionRuidoproceso = 0.2; % Ruido aceleracion (m/seg^2)
Sz = DesviacionRuidomedicion^2; % covarianza de error de medicion
Sw = DesviacionRuidoproceso^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]; %covarianza error del proceso
P = Sw; % Cobarianza estimacion inicial

P1_1=(DesviacionRuidoproceso^2*dt^4)/4;
P1_2=(DesviacionRuidoproceso^2*dt^3)/2;
P2_1=(DesviacionRuidoproceso^2*dt^3)/2;
P2_2=DesviacionRuidoproceso^2*dt^2;

% Inicializacion vectores para graficacion

pos = []; % Vector Posicion real 
poshat = []; % Vector Posicion estimada 
posmedido = []; % Vector Posicion medida 

vel = []; % Vector Velocidad real 
velhat = []; % Vector Velocidad estimada 


w=0;
for t=0:dt:duracion
% Entrada Aceleracion m/seg^2.
w=w+dt;
u=5*sin(w*5);
% Simulacion 

%Simulacion ruido del proceso
Ruidoproceso=DesviacionRuidoproceso*[(dt^2/2)*randn;dt*randn];

% Simulacion sistema
x=a*x+b*u+Ruidoproceso;

%Simulacion ruido de la medicion
RuidoMedicion=DesviacionRuidomedicion*randn;

% Simulacion medida
y=c*x+RuidoMedicion;

% 
%Kalman

if (0)%calculos matlab
    
    
% Prediccion

%P1
%Prediccion Estado
xhat=a*xhat+b*u;


%P2

%Covarianza Inovacion.
s=c*P*c'+Sz;

% Proyeccion del error
P=a*P*a'-a*P*c'*inv(s)*c*P*a'+Sw;


%Actualizacion

%A1
%Ganancia de Kalman
K=a*P*c'*inv(s);

%A2
% Inovacion.
Inn=y-c*xhat;
%Actualizacion Estimado
xhat=xhat+K*Inn;

else %calculos para programar
UA=u;

ymedida=y;

P1_1=P1_1 + P2_1*dt + dt*(P1_2 + P2_2*dt) + (DesviacionRuidoproceso^2*dt^4)/4 - (P1_1*(P1_1 + P2_1*dt))/(DesviacionRuidomedicion^2 + P1_1) - (P1_2*dt*(P1_1 + P2_1*dt))/(DesviacionRuidomedicion^2 + P1_1);
P1_2=P1_2 + P2_2*dt + (DesviacionRuidoproceso^2*dt^3)/2 - (P1_2*(P1_1 + P2_1*dt))/(DesviacionRuidomedicion^2 + P1_1);
P2_1=P2_1 + P2_2*dt + (DesviacionRuidoproceso^2*dt^3)/2 - (P1_1*P2_1)/(DesviacionRuidomedicion^2 + P1_1) - (P1_2*P2_1*dt)/(DesviacionRuidomedicion^2 + P1_1);
P2_2=DesviacionRuidoproceso^2*dt^2 + P2_2 - (P1_2*P2_1)/(DesviacionRuidomedicion^2 + P1_1);


xposhat=xposhat + dt*xvelhat + (UA*dt^2)/2 - (((UA*dt^2)/2 + xvelhat*dt + xposhat - ymedida)*(P1_1 + P2_1*dt + dt*(P2_1 + P2_2*dt + (DesviacionRuidoproceso^2*dt^3)/2 - (P1_1*P2_1)/(DesviacionRuidomedicion^2 + P1_1) - (P1_2*P2_1*dt)/(DesviacionRuidomedicion^2 + P1_1)) + dt*(P1_2 + P2_2*dt) + (DesviacionRuidoproceso^2*dt^4)/4 - (P1_1*(P1_1 + P2_1*dt))/(DesviacionRuidomedicion^2 + P1_1) - (P1_2*dt*(P1_1 + P2_1*dt))/(DesviacionRuidomedicion^2 + P1_1)))/(DesviacionRuidomedicion^2 + P1_1);
xvelhat=xvelhat + UA*dt - (((UA*dt^2)/2 + xvelhat*dt + xposhat - ymedida)*(P2_1 + P2_2*dt + (DesviacionRuidoproceso^2*dt^3)/2 - (P1_1*P2_1)/(DesviacionRuidomedicion^2 + P1_1) - (P1_2*P2_1*dt)/(DesviacionRuidomedicion^2 + P1_1)))/(DesviacionRuidomedicion^2 + P1_1);

end

% Crear vectores de datos para graficacion
pos = [pos; x(1)];
posmedido = [posmedido; y];
poshat = [poshat; xposhat];
vel = [vel; x(2)];
velhat = [velhat; xvelhat];
end

% Graficacion
close all;
t=0:dt:duracion;

figure()
subplot(2,2,1)
plot(t,pos,'b-',t,posmedido,'g-.', t,poshat,'r--','LineWidth',2);
grid on
xlabel('Tiempo (seg)');
ylabel('Posicion (m)');
title('Posicon Vehiculo (True, Measured, and Estimated)')
legend('Real','Medida','Estimada')

subplot(2,2,2)
plot(t,pos-posmedido,'g-.',t,pos-poshat,'r--','LineWidth',2);
grid on 
xlabel('Tiempo (seg)');
ylabel('Position Error (feet)');
title('Error Posicion');
legend('Medida','Estimada')


subplot(2,2,3)
plot(t,vel,'b-', t,velhat,'r--','LineWidth',2);
grid on
xlabel('Tiempo (seg)');
ylabel('Velocidad (m/seg)');
title('Velocidad');
legend('Real','Estimada')

subplot(2,2,4)
plot(t,vel-velhat,'r--','LineWidth',2);
grid on
xlabel('Tiempo (seg)');
ylabel('Error Velocidad (m/seg)');
title('Error Velocidad');
legend('Estimada')