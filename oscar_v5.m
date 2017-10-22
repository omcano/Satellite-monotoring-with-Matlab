%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Trabalho SCDS. Óscar Martínez Cano 45604616D     05/2014
%   Nº alumno: 63994
%
% É um software para simular a posição na superfície da Terra dum satélite
% qualquer a partir de seus elementos keplerianos obtidos de Internet e 
% também calcular a distancia desde um observador ate o satélite, o ângulo
% de elevação e o azimute.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all;close all; %Limpar tudo

%% Constantes
rad2deg = 180/pi;        % radianes para graus
deg2rad = pi/180;        % graus para radianes
Npontos = 1000;          % Numero de pontos por orbita (Resoloção)

%% Entradas
disp('-------------------------------');
disp('Satelites disponiveis');  
disp('     (1) ISS (ZARYA)');  
disp('     (2) Iridium  33');
disp('     (3) GPS BIIA-10 (PRN 32)');
disp('     (4) NOAA 15');
disp('-------------------------------');  
satSelect = input('Escolha umos dos satelite disponiveis: ');
disp('-------------------------------');  
Norb = input('Escolha o numero de orbitas para simular: ');
disp('-------------------------------');
disp('Posiçoes disponiveis'); 
disp('     (1) Lisboa');  
disp('     (2) Madrid');
disp('     (3) Río de Janeiro');
disp('     (4) Pekin');
disp('-------------------------------');  
posOb = input('Escolha posição do observador: ');
disp('-------------------------------'); 
%% Modelo da Terra
mu = 3.986005e14;                   % G*MT [m^3*s^-2]
radEc = 6378.137;                   % Radio ecuatorial [Km] 
Tterra = 86164.091;                 % dia sideral [s]
Wt = 2*pi/Tterra;                   % Velocidade angular da terra [rad/s]

%% Dados keplerianos das órbitas
%%de estaçoes, sat de comunicaçoes, meteorologicos e GPS obtidos de Internet 
TXT = fopen('satellites.txt','w');                                                 % cria uma TXT para guardar os dados
strURL = urlread('http://www.celestrak.com/NORAD/elements/stations.txt');          % ler os dados da URL
fprintf(TXT,strURL);                                                               % escreve os dados da url no TXT
strURL = urlread('http://www.celestrak.com/NORAD/elements/iridium-33-debris.txt');
fprintf(TXT,strURL);
strURL = urlread('http://www.celestrak.com/NORAD/elements/gps-ops.txt');
fprintf(TXT,strURL);
strURL = urlread('http://www.celestrak.com/NORAD/elements/weather.txt');
fprintf(TXT,strURL);
fclose(TXT);                                                                       % fecha o documento

%% Seleccionar o satelite da lista descargada de internet
if satSelect == 1                        % ISS (ZARYA): estaçoes   
    salto = 1;    nome = 'ISS (ZARYA)';
elseif satSelect ==2                     % Iridium  33 : comunicaçoes
    salto = 103;  nome = 'Iridium  33';
elseif satSelect ==3                     % GPS BIIA-10 (PRN 32)  : navigations
    salto = 1228; nome = 'GPS BIIA-10 (PRN 32)';
elseif satSelect == 4                    % NOAA 15 : weather
    salto = 1324; nome = 'NOAA 15';
end
%% Seleccionar posição do observador
if posOb == 1                         
    latObs = 38.717;       longObs = -9.167;
elseif posOb ==2                     
    latObs = 40.418889;    longObs = -3.691944;
elseif posOb ==3                     
    latObs = -22.9;        longObs = -43.233056;
elseif posOb == 4                    
    latObs = 39.905;       longObs = 116.391389;
end
%% Guardar as dados keplerianos do TXT em variaveis de MatLab
[Nlinha1,Nsat,Nsat2,epoch,decay,cero,desco,cero2,Set] = textread('satellites.txt','%d %s %s %f %f %s %s %d %d',1,'headerlines',salto);
[Nlinha2,Ncata,Ainq,RAAN,excent,argPer,meanA,meanMot] = textread('satellites.txt','%d %d %f %f %d %f %f %f',1,'headerlines',salto+1);
excent = excent*1e-7;             %Correcçao da excentricidade

%% Dados keplerianos

Torb = (1/meanMot)*24*3600;       % Periodo orbital [s]
Worb = 2*pi/Torb;                 % Velocidade angular [rad/s]

    omegaM = RAAN    *  deg2rad;  % [0-360] Anngulo equinocio-nodulo ascendente [rad]
    inq    = Ainq    *  deg2rad;  % [0 90] Angulo de inclinaçao da órbita [rad]
    omegam = argPer  *  deg2rad;  % [0-360]Angulo nodulo ascendente-periastro  [rad]

    sEixeM = nthroot((mu*(Torb/(2*pi))^2),3) *1e-3;  %Semi eixe maior [Km]
    excOrb = excent;              % excentricidade [0 1]
    t0 = epoch;                   % Tempo de pasagem no perigeu
%% Vector de tempos
t = linspace(0,Norb*Torb,Norb*Npontos);   % vector de tempos linearmente espçados Npontos/Norbitas

%% Calculo das anomalias
n = 2*pi/Torb;            % movimiento médio [rad/s]
M0 = meanA * deg2rad;     % anomalia media incial [rad]
M = M0 + n*(t-t0);        % anmalia média nem cualquier instante [rad]
    
for j=1:Norb*Npontos
    E(j) = kepler_E(excOrb,M(j));                               % eccentric anomaly [rad]
    theta(j) = 2*atan(sqrt((1+excOrb)/(1-excOrb))*tan(E(j)/2)); % verdadera anomalia [rad]
end

%% Calculo do trazo da trayectoria na terra.
p = sEixeM*(1-excOrb^2);                 %raio da elipse quando se cruza o eixo dos yy, i.e [Km]
r = p./(1 + excOrb*cos(theta + omegam)); % raio orbital [Km]
ra = sEixeM*(1+excOrb);                  % distança apogeu-foco [Km]
rp = sEixeM*(1-excOrb);                  % distança perigeu-foco [Km]

xs = r.*cos(theta).*cos(RAAN)-r.*sin(theta).*cos(inq).*sin(RAAN); % ECI x-coordinate SAT [km]
ys = r.*cos(theta).*sin(RAAN)+r.*sin(theta).*cos(inq).*cos(RAAN); % ECI y-coordinate SAT [km]
rs = p./(1+excOrb.*cos(theta-omegam));                            % norm of radius SAT   [km]        
    
latTraRad = asin(sin(inq).*sin(theta));
lonTraRad = (atan2(ys./rs,xs./rs)-Wt*t);

latiTrazo = latTraRad *rad2deg; % Latitude [rad-> deg]
longTrazo = lonTraRad *rad2deg;% Longitude[rad]

%% Calculo do distacia ao setélite, elevaçao e azimute
Pobservador = [latObs longObs];   %lat-long

latOb = Pobservador(1)* deg2rad;
lonOb = Pobservador(2)* deg2rad;

L = lonTraRad-lonOb;
posNul = zeros(1,length(L));
for i = 1:length(L)
    if (L(i)<0)          %L não pode ser menor do que 0
        L(i) = L(i) + 2*pi;
        %posNul(i)=1;
    end
end

  fiM = acos(cos(latTraRad).*cos(L).*cos(latOb)+ sin(latOb).*sin(latTraRad));

dis = sqrt(radEc.^2 + r.^2 - 2.*radEc.*r.*cos(fiM));   %[Km]  

elev = atan((cos(fiM)-radEc./r)./sqrt(1-cos(fiM).^2)) *rad2deg;

azimute = asin(sin(L).*cos(fiM)./sin(fiM)) *rad2deg;
%% Representaçao dos dados

ah = axes('unit', 'normalized', 'position', [0 0 1 1]);  % create an axes that spans the whole gui
bg = imread('cielo.jpg'); imagesc(bg);                   % import the background image and show it on the axes
set(ah,'handlevisibility','off','visible','off')         % prevent plotting over the background and turn the axis off
uistack(ah, 'bottom');                                   % making sure the background is behind all the other uicontrols
hold on
planisfer = worldmap('world');                           % carrega planisferio
setm(planisfer, 'FFaceColor','blue')                     % cor de fondo do planisferio
geoshow('landareas.shp', 'FaceColor', [0.4 0.9 0.2])     % carrega os continentes
geoshow('worldlakes.shp', 'FaceColor', 'cyan')           % carrega os lagos
geoshow('worldrivers.shp', 'Color', 'cyan')              % carrega os rios
 
title(nome,'color','red','fontSize',18)                  % titulo
plotm(latiTrazo,longTrazo,'Color','white','linewidth',1,'LineStyle','--')   % trazo da orbita
plotm(Pobservador,'color','red','LineWidth',2,'Marker','o');

for i=1:length(t)
    plotm(latiTrazo(i),longTrazo(i),'Marker','o',...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','r',...
                'MarkerSize',4);
    drawnow;
        disp(['Distaça: ', num2str(dis(i)),' Km. ',...
              'Elevaçao: ', num2str(elev(i)),' Graus. ',...
              'Azimute: ', num2str(azimute(i)),' Graus. ']);
    pause(0.00001); % pause n secons
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%OMC