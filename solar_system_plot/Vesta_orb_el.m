%% Vesta Orbital Elements

vesta = struct('M',1,'SOI',0.00026244*AU,'m_anom0',0, ...
               'ra',1,'rp',1,'a',1,'e',1,'h',1,'p',1,'u_sol',1,'pos',eye(N,3));

vesta.M = 2.589e20;                   %mass [kg]
vesta.u_sol = G*(M_sun+vesta.M);      %[m^3/s^2]
vesta.ra = 3.8467e+08;                %distanza afelio [km]
vesta.rp = 3.2197e+08;                %distanza perielio [km]
W = 103.810803*(pi/180);              %Longitude of ascending node [rad]            
w = 150.72854*(pi/180);               %argomento del periasse [rad] (longitude_of_perihelio - W = w, vedi wiki)                  
i = 7.141771*(pi/180);                %Orbital inclination  

vesta.a = (vesta.ra+vesta.rp)/2;                       %semiasse maggiore [km]
vesta.e = (vesta.ra-vesta.rp)/(vesta.ra+vesta.rp);     %modulo di eccentricità [km]
vesta.p = vesta.a*(1-vesta.e^2);                       %valore di semibraccio medio [km]
vesta.h = (vesta.u_sol*vesta.a*1e3*(1-vesta.e^2))^0.5; %modulo h [m^2/s]

vesta.Meulero = eulerT(W,w,i);        %matrice di rotazione

%calcolo dei vettori posizione
for i = 1:N
    r = vesta.p/(1+vesta.e*cos(Tet(i)));       %modulo del raggio posizione

    X = r*cos(Tet(i));                         %componenti nel piano orbitale
    Y = r*sin(Tet(i));

    vesta.pos(i,1) = vesta.Meulero(1,1)*X + vesta.Meulero(1,2)*Y;     %componenti nello spazio
    vesta.pos(i,2) = vesta.Meulero(2,1)*X + vesta.Meulero(2,2)*Y;
    vesta.pos(i,3) = vesta.Meulero(3,1)*X + vesta.Meulero(3,2)*Y;
end

    