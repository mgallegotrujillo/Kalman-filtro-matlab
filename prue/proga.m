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

dt=sym('dt','real');
xpos=sym('xpos','real');
xvel=sym('xvel','real');

a = [1 dt; 0 1]; % Matriz de transicion
b = [dt^2/2; dt]; % Matriz de entrada
c = [1 0]; % Matriz de sensores
x = [xpos; xvel]; % Estado inicial




DesviacionRuidomedicion=sym('DesviacionRuidomedicion','real');
DesviacionRuidoproceso=sym('DesviacionRuidoproceso','real');

%Filtro de kalman

xhat = x % Estado inicial estimacion
% DesviacionRuidomedicion = 10; % ruido de medición de posición (m)
% DesviacionRuidoproceso = 0.2; % Ruido aceleracion (m/seg^2)

Sz = DesviacionRuidomedicion^2; % covarianza de error de medicion
Sw = DesviacionRuidoproceso^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]; %covarianza error del proceso
% P = Sw; % Cobarianza estimacion inicial
P=sym('P', [2 2])


% Inicializacion vectores para graficacion

pos = []; % Vector Posicion real 
poshat = []; % Vector Posicion estimada 
posmedido = []; % Vector Posicion medida 

vel = []; % Vector Velocidad real 
velhat = []; % Vector Velocidad estimada 


syms UA
% Entrada Aceleracion m/seg^2.
u=UA;

% Simulacion 

%Simulacion ruido del proceso
Ruidoproceso=DesviacionRuidoproceso*[(dt^2/2)*randn;dt*randn];

% Simulacion sistema
x=a*x+b*u+Ruidoproceso;

%Simulacion ruido de la medicion
RuidoMedicion=DesviacionRuidomedicion*randn;

% Simulacion medida
syms ymedida
y=c*x+RuidoMedicion;



%Kalman


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
Inn=ymedida-c*xhat;
%Actualizacion Estimado
xhat=xhat+K*Inn;
