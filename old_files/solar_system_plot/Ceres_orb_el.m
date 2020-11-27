%% Ceres Orbital Elements

ceres = struct('M',1,'SOI',0.00050486*AU,'m_anom0',0, ...
               'ra',1,'rp',1,'a',1,'e',1,'h',1,'p',1,'u_sol',1,'pos',eye(N,3));

ceres.M = 8.958e20;                   %mass [kg]
ceres.u_sol = G*(M_sun+ceres.M);      %[m^3/s^2]
ceres.ra = 3.80951e+08;               %distanza afelio [km]
ceres.rp = 4.46428e+08;               %distanza perielio [km]
W = 80.30553090445737*(pi/180);       %Longitude of ascending node [rad]            
w = 73.59769469844186*(pi/180);       %argomento del periasse [rad] (longitude_of_perihelio - W = w, vedi wiki)          
i = 10.59406719506626*(pi/180);       %Orbital inclination



ceres.a = 2.769165148633284*AU;                        %semiasse maggiore [km]
ceres.e = .0760090265983052;                           %modulo di eccentricità [km]
ceres.p = ceres.a*(1-ceres.e^2);                       %valore di semibraccio medio [km]
ceres.h = (ceres.u_sol*ceres.a*1e3*(1-ceres.e^2))^0.5; %modulo h [m^2/s]

ceres.Meulero = eulerT(W,w,i);        %matrice di rotazione

%calcolo dei vettori posizione
for i = 1:N
    r = ceres.p/(1+ceres.e*cos(Tet(i)));       %modulo del raggio posizione

    X = r*cos(Tet(i));                         %componenti nel piano orbitale
    Y = r*sin(Tet(i));

    ceres.pos(i,1) = ceres.Meulero(1,1)*X + ceres.Meulero(1,2)*Y;     %componenti nello spazio
    ceres.pos(i,2) = ceres.Meulero(2,1)*X + ceres.Meulero(2,2)*Y;
    ceres.pos(i,3) = ceres.Meulero(3,1)*X + ceres.Meulero(3,2)*Y;
end
