function [mercurio,venere,terra,marte,p] = gener_orb(N,Tet,p_nodes)
    % per generare l'orbita 2D! devo ricevere i parametri
    % ra = distanza afelio e rp = distanza perielio
    %le unità di misura della distanza sono in km

    AU = 149597870.7;               %in km
    
    mercurio = struct('M',1,'SOI',0.00075*AU,'m_anom0',0, ...
                      'ra',1,'rp',1,'a',1,'e',1,'h',1,'p',1,'u_sol',1,'pos',eye(N,3));
                  
    venere =   struct('M',1,'SOI',0.00411*AU,'m_anom0',0, ...
                      'ra',1,'rp',1,'a',1,'e',1,'h',1,'p',1,'u_sol',1,'pos',eye(N,3),...
                      'Meulero' ,zeros(3,3));
                  
    terra =    struct('M',1,'SOI',0.00621*AU,'m_anom0',0, ...
                      'ra',1,'rp',1,'a',1,'e',1,'h',1,'p',1,'u_sol',1,'pos',eye(N,3));
                  
    marte =    struct('M',1,'SOI',0.00385*AU,'m_anom0',0, ...
                      'ra',1,'rp',1,'a',1,'e',1,'h',1,'p',1,'u_sol',1,'pos',eye(N,3));
    
    %definizione costanti universali
    G = 6.67e-11;         %in forma SI
    M_sun = 1.9891e30;              %kg
    
    %traccio orbita terrestre
    M_earth = 59.736e23;
    u_sol = G*(M_sun+M_earth);         %in m^3/s^2
    ra = 152.097701e6;                  %km
    rp = 147.097074e6;                  %km
    W = -11.26064*(pi/180);                     %rad
    w = (102.94719*(pi/180))-W;                 %rad
    i = 0.00005*(pi/180);                       %rad
    
    a = (ra+rp)/2;                         %semiasse maggiore (Km)
    ee = (ra-rp)/(ra+rp);       %modulo di eccentricità
    p = a*(1-ee^2);                         %valore di semibraccio medio (km)
    h = (u_sol*a*1e3*(1-ee^2))^0.5;   %modulo h espresso in m^2/s
    
    M = eulerT(W,w,i);
    %calcolo dei vettori posizione
    for i = 1:N
        r = p/(1+ee*cos(Tet(i)));               %modulo del raggio posizione
    
        X = r*cos(Tet(i));
        Y = r*sin(Tet(i));
        
        terra.pos(i,1) = M(1,1)*X + M(1,2)*Y;
        terra.pos(i,2) = M(2,1)*X + M(2,2)*Y;
        terra.pos(i,3) = M(3,1)*X + M(3,2)*Y;
    end
    
    terra.ra = ra;
    terra.rp = rp;
    terra.a = a;
    terra.e = ee;
    terra.p = p;
    terra.u_sol = u_sol;
    terra.h = h;
    terra.M = M_earth;
    
    
    %traccio orbita venere
    M_venus = 48.685e23;                        %kg
    venere.M = M_venus;
    venere.u_sol = G*(M_sun+M_venus);        %m^3/s^2
    venere.ra = 108.94e6;                    %km
    venere.rp = 107.48e6;                    %km
    W = 76.68069*(pi/180);                      %rad
    w = (131.53298*(pi/180))-W;                 %rad
    i = 3.39471*(pi/180);                       %rad
    
    venere.a = (venere.ra+venere.rp)/2;                      %semiasse maggiore (km)
    venere.e = (venere.ra-venere.rp)/(venere.ra+venere.rp);        %modulo di eccentricità
    venere.p = venere.a*(1-venere.e^2);                              %valore di semibraccio medio (km)
    venere.h = (venere.u_sol*venere.a*1e3*(1-venere.e^2))^0.5;            %modulo h (m^2/s)
    
    M = eulerT(W,w,i);
    venere.Meulero = M;
    %calcolo posizione venere
    for i = 1:N
        r = venere.p/(1+venere.e*cos(Tet(i)));               %modulo del raggio posizione
    
        X = r*cos(Tet(i));
        Y = r*sin(Tet(i));
        
        venere.pos(i,1) = M(1,1)*X + M(1,2)*Y;
        venere.pos(i,2) = M(2,1)*X + M(2,2)*Y;
        venere.pos(i,3) = M(3,1)*X + M(3,2)*Y;

    end
    
    %traccio orbita marte
    M_marte = 6.4184e23;                        %kg
    marte.M = M_marte;
    marte.u_sol = G*(M_sun+M_marte);         %m^3/s^2
    marte.ra = 249.23e6;                     %km
    marte.rp = 206.62e6;                     %km
    W = 49.57854*(pi/180);                      %rad
    w = (336.04084*(pi/180))-W;                 %rad
    i = 1.85061*(pi/180);                       %rad
    
    marte.a = (marte.ra+marte.rp)/2;                          %semiasse maggiore (km)
    marte.e = (marte.ra-marte.rp)/(marte.ra+marte.rp);        %modulo di eccentricità
    marte.p = marte.a*(1-marte.e^2);                          %valore di semibraccio medio (km)
    marte.h = (marte.u_sol*marte.a*1e3*(1-marte.e^2))^0.5;    %modulo h (m^2/s)
    
    M = eulerT(W,w,i);
    %calcolo posizione marte
    for i = 1:N
        r = marte.p/(1+marte.e*cos(Tet(i)));               %modulo del raggio posizione
        
        X = r*cos(Tet(i));
        Y = r*sin(Tet(i));
        
        marte.pos(i,1) = M(1,1)*X + M(1,2)*Y;
        marte.pos(i,2) = M(2,1)*X + M(2,2)*Y;
        marte.pos(i,3) = M(3,1)*X + M(3,2)*Y;

    end
    
    %traccio orbita mercurio
    M_merc = 3.3011e23;                         %kg
    mercurio.M = M_merc;
    mercurio.u_sol = G*(M_sun+M_merc);      %m^3/s^2
    mercurio.ra = 69.82e6;                   %km
    mercurio.rp = 46e6;                      %km
    W = 48.34*(pi/180);                         %rad
    w = (77.46*(pi/180))-W;                     %rad
    i = 7.006*(pi/180);                         %rad
    
    mercurio.a = (mercurio.ra+mercurio.rp)/2;                            %semiasse maggiore (km)
    mercurio.e = (mercurio.ra-mercurio.rp)/(mercurio.ra+mercurio.rp);    %modulo di eccentricità
    mercurio.p = mercurio.a*(1-mercurio.e^2);                            %valore di semibraccio medio (km)
    mercurio.h = (mercurio.u_sol*mercurio.a*1e3*(1-mercurio.e^2))^0.5;   %modulo h (m^2/s)
    
    M = eulerT(W,w,i);
    %calcolo posizione mercurio
    for i = 1:N
        r = mercurio.p/(1+mercurio.e*cos(Tet(i)));               %modulo del raggio posizione
    
        X = r*cos(Tet(i));
        Y = r*sin(Tet(i));
        
        mercurio.pos(i,1) = M(1,1)*X + M(1,2)*Y;
        mercurio.pos(i,2) = M(2,1)*X + M(2,2)*Y;
        mercurio.pos(i,3) = M(3,1)*X + M(3,2)*Y;

 
    end
    
    %codice riconoscimento nodi pianeti
    if p_nodes == 1
        for i = 1:N-1
            %mercurio
            if (mercurio.pos(i,3) >= 0 && mercurio.pos(i+1,3) <0)
                plot3(mercurio.pos(i,1),mercurio.pos(i,2),mercurio.pos(i,3),'black x');hold on;
            end
            if (mercurio.pos(i+1,3) >= 0 && mercurio.pos(i,3) <0)
                plot3(mercurio.pos(i,1),mercurio.pos(i,2),mercurio.pos(i,3),'black x');hold on;
            end
            
            %venere
            if (venere.pos(i,3) >= 0 && venere.pos(i+1,3) <0)
                plot3(venere.pos(i,1),venere.pos(i,2),venere.pos(i,3),'black x');hold on;
            end
            if (venere.pos(i+1,3) >= 0 && venere.pos(i,3) <0)
                plot3(venere.pos(i,1),venere.pos(i,2),venere.pos(i,3),'black x');hold on;
            end
            
            %terra
            if (terra.pos(i,3) >= 0 && terra.pos(i+1,3) <0)
                plot3(terra.pos(i,1),terra.pos(i,2),terra.pos(i,3),'black x');hold on;
            end
            if (terra.pos(i+1,3) >= 0 && terra.pos(i,3) <0)
                plot3(terra.pos(i,1),terra.pos(i,2),terra.pos(i,3),'black x');hold on;
            end
            
            %marte
            if (marte.pos(i,3) >= 0 && marte.pos(i+1,3) <0)
                plot3(marte.pos(i,1),marte.pos(i,2),marte.pos(i,3),'black x');hold on;
            end
            if (marte.pos(i+1,3) >= 0 && marte.pos(i,3) <0)
                plot3(marte.pos(i,1),marte.pos(i,2),marte.pos(i,3),'black x');hold on;
            end

        end
    end
end