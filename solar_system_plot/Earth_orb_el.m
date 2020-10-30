%% Earth Orbital Elements

terra = struct('M',1,'SOI',0.00621*AU,'m_anom0',0, ...
               'ra',1,'rp',1,'a',1,'e',1,'h',1,'p',1,'u_sol',1,'pos',eye(N,3));

terra.M = 59.736e23;                  %mass [kg]
terra.u_sol = G*(M_sun+terra.M);      %[m^3/s^2]
terra.ra = 152.097701e6;              %distanza afelio [km]
terra.rp = 147.097074e6;              %distanza perielio [km]
W = -11.26064*(pi/180);               %Longitude of ascending node [rad]            
w = (102.94719*(pi/180))-W;           %argomento del periasse [rad] (longitude_of_perihelio - W = w, vedi wiki)          
i = 0.00005*(pi/180);                 %Orbital inclination  

terra.a = (terra.ra+terra.rp)/2;                       %semiasse maggiore [km]
terra.e = (terra.ra-terra.rp)/(terra.ra+terra.rp);     %modulo di eccentricità [km]
terra.p = terra.a*(1-terra.e^2);                       %valore di semibraccio medio [km]
terra.h = (terra.u_sol*terra.a*1e3*(1-terra.e^2))^0.5; %modulo h [m^2/s]

terra.Meulero = eulerT(W,w,i);        %matrice di rotazione

%calcolo dei vettori posizione
for i = 1:N
    r = terra.p/(1+terra.e*cos(Tet(i)));       %modulo del raggio posizione

    X = r*cos(Tet(i));                         %componenti nel piano orbitale
    Y = r*sin(Tet(i));

    terra.pos(i,1) = terra.Meulero(1,1)*X + terra.Meulero(1,2)*Y;     %componenti nello spazio 3D
    terra.pos(i,2) = terra.Meulero(2,1)*X + terra.Meulero(2,2)*Y;
    terra.pos(i,3) = terra.Meulero(3,1)*X + terra.Meulero(3,2)*Y;
end

    