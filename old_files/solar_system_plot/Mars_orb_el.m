%% Mars Orbital Elements

marte =    struct('M',1,'SOI',0.00385*AU,'m_anom0',0, ...
                      'ra',1,'rp',1,'a',1,'e',1,'h',1,'p',1,'u_sol',1,'pos',eye(N,3));

marte.M = 6.4184e23;                  %mass [kg]
marte.u_sol = G*(M_sun+marte.M);      %[m^3/s^2]
marte.ra = 249.23e6;                  %distanza afelio [km]
marte.rp = 206.62e6;                  %distanza perielio [km]
W = 49.57854*(pi/180);                %Longitude of ascending node [rad]            
w = (336.04084*(pi/180))-W;           %argomento del periasse [rad] (longitude_of_perihelio - W = w, vedi wiki)          
i = 1.85061*(pi/180);                 %Orbital inclination  

marte.a = (marte.ra+marte.rp)/2;                       %semiasse maggiore [km]
marte.e = (marte.ra-marte.rp)/(marte.ra+marte.rp);     %modulo di eccentricità [km]
marte.p = marte.a*(1-marte.e^2);                       %valore di semibraccio medio [km]
marte.h = (marte.u_sol*marte.a*1e3*(1-marte.e^2))^0.5; %modulo h [m^2/s]

marte.Meulero = eulerT(W,w,i);        %matrice di rotazione

%calcolo dei vettori posizione
for i = 1:N
    r = marte.p/(1+marte.e*cos(Tet(i)));       %modulo del raggio posizione

    X = r*cos(Tet(i));                         %componenti nel piano orbitale
    Y = r*sin(Tet(i));

    marte.pos(i,1) = marte.Meulero(1,1)*X + marte.Meulero(1,2)*Y;     %componenti nello spazio
    marte.pos(i,2) = marte.Meulero(2,1)*X + marte.Meulero(2,2)*Y;
    marte.pos(i,3) = marte.Meulero(3,1)*X + marte.Meulero(3,2)*Y;
end

    