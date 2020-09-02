%% satellite Orbital Elements

satellite = struct('M',1,'SOI',0,'m_anom0',0, ...
               'ra',1,'rp',1,'a',1,'e',1,'h',1,'p',1,'u_sol',1,'pos',eye(N,3));

% i parametri non necessari sono commentati
           
% satellite.M = ;                               %mass [kg]
% satellite.u_sol = G*(M_sun+satellite.M);      %[m^3/s^2]
% satellite.ra = ;                              %distanza afelio [km]
% satellite.rp = ;                              %distanza perielio [km]

W = 4.57297;                                    %Longitude of ascending node [rad]            
w = 1.88201;                                    %argomento del periasse [rad] (longitude_of_perihelio - W = w, vedi wiki)          
i = 2.93393;                                    %Orbital inclination  [rad]

satellite.a = 3.98959e+08;                                                  %semiasse maggiore [km]
satellite.e =  0.292987;                                                    %modulo di eccentricità [km]
satellite.p = satellite.a*(1-satellite.e^2);                                %valore di semibraccio medio [km]
% satellite.h = (satellite.u_sol*satellite.a*1e3*(1-satellite.e^2))^0.5;    %modulo h [m^2/s]

satellite.Meulero = eulerT(W,w,i);        %matrice di rotazione

%calcolo dei vettori posizione
for i = 1:N
    r = satellite.p/(1+satellite.e*cos(Tet(i)));       %modulo del raggio posizione

    X = r*cos(Tet(i));                                  %componenti nel piano orbitale
    Y = r*sin(Tet(i));

    satellite.pos(i,1) = satellite.Meulero(1,1)*X + satellite.Meulero(1,2)*Y;     %componenti nello spazio
    satellite.pos(i,2) = satellite.Meulero(2,1)*X + satellite.Meulero(2,2)*Y;
    satellite.pos(i,3) = satellite.Meulero(3,1)*X + satellite.Meulero(3,2)*Y;
end

    