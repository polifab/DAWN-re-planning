%conosco le posizioni di venere e navicella nel momento in cui ho finito il
%flyby, da queste ne conosco già le velocità, e, in base ad r0 (distanza in
%modulo del periasse con venere del moto iperbolico, costruisco la rotta
navi2 = struct;
pos_venere2 = [2.2639e7     1.0537e8    1.3316e5]';                %posizione venere uscita sole
pos_navi2 =     [2.4826e7     1.0489e8    350.9542]';               %posizione navicella uscita sole
vel_navi2 =      [-32.44104   12.22933    2.081501]';               %velocità uscita navicella rispetto sole
vel_venere2 =  [-34.3527     7.1799        2.0810    ]';                %velocità venere uscita rispetto il sole

Ceulers_orb2 = [venere.Meulero    [0 0 0]';
                         0   0   0                 1                  ];                              %matrice cambio variabili venere->sole
navi2.f = 5.34677035598;
Ceulers_orb1 = [cos(navi2.f)    -sin(navi2.f)   0   sqrt(dot(pos_venere2,pos_venere2))*cos(navi2.f);...
                           sin(navi2.f)     cos(navi2.f)    0   sqrt(dot(pos_venere2,pos_venere2))*sin(navi2.f);...
                           0                     0                      1   0;...
                           0                     0                      0   1];
navi2.Ceuler = Ceulers_orb2*Ceulers_orb1;
                     
pos_navi2v = (inv(navi2.Ceuler))*[pos_navi2; 1];                            %vettore posizione navicella uscita rispetto venere
vel_navi2v = (inv(navi2.Ceuler))*[vel_navi2-vel_venere2; 0];         %velocità navicella uscita rispetto venere

v_mod2 = sqrt(dot(vel_navi2v,vel_navi2v));                     %modulo velocità vinf+

%navi2.R0 = r0+6051.84;
%navi2.V0 = sqrt(v_mod2^2+(2*6.67e-20*venere.M/navi2.R0));           %velocità al periasse in modulo
%navi2.a = 6.67e-20*venere.M/(v_mod2^2);                                             %semiasse maggiore
%navi2.e = 1+navi2.R0/navi2.a;                                                                          %ellitticità
%navi2.p = -navi2.a*(1-navi2.e^2);                                                           

%creo il vettore di angoli per la simulazione del satellite
%navi2.Tet_max = acos((1/navi2.e)*(navi2.p/venere.SOI - 1));                 %angolo di simulazione
%navi2.Tet = -navi2.Tet_max:ris:navi2.Tet_max;                                                           %vettore degli angoli per la simulazione grafica

 
navi2.posF = [pos_navi2v(1,1); pos_navi2v(2,1); pos_navi2v(3,1)];
navi2.velF = [vel_navi2v(1,1); vel_navi2v(2,1); vel_navi2v(3,1)];
[navi2.a,navi2.e,navi2.i_orb,navi2.W,navi2.w,navi2.t_anom,navi2.p,navi2.ra,navi2.rp] = rv_to_orb(navi2.posF, navi2.velF, venere.M*6.67e-11*1e-9);

navi2.Tet = -navi2.t_anom:ris:navi2.t_anom;
%navi2.w = pi+atan(pos_navi2v(2,1)/pos_navi2v(1,1))-navi2.Tet_max;             %argomento del periasse              
%navi2.W = 0;
%navi2.i_orb = 0;
navi2.Meulero = eulerT(navi2.W,navi2.w,navi2.i_orb);                            %trasformazione di eulero, necessaria per lo studio

%creazione orbita navicella
[navi2.dimr, navi2.dimc] = size(navi2.Tet);
navi2.pos = zeros(navi2.dimc,3);
for ind = 1:navi2.dimc
    navi2.r = navi2.p/(1+navi2.e*cos(navi2.Tet(ind)));
    navi2.x = navi2.r*cos(navi2.Tet(ind));
    navi2.y = navi2.r*sin(navi2.Tet(ind));
    
    navi2.pos(ind,1) = navi2.Meulero(1,:)*[navi2.x ; navi2.y; 0];
    navi2.pos(ind,2) = navi2.Meulero(2,:)*[navi2.x ; navi2.y; 0];
    navi2.pos(ind,3) = navi2.Meulero(3,:)*[navi2.x ; navi2.y; 0];
end

ff = figure;
plot3(navi2.pos(:,1),navi2.pos(:,2),navi2.pos(:,3),'magenta');hold on;
plot3(0,0,0,'green .','MarkerSize',10);
axis equal;
xlabel('x(km)');
ylabel('y(km)');
zlabel('z(km)');

%plot moto della navicella
%navi2.EE = (navi2.e+cos(navi2.Tet(1)))/(1+navi2.e*cos(navi2.Tet(1)));
%navi2.ecc_anom0 = log(navi2.EE+sqrt(navi2.EE^2-1));                                              %in rad anomalia eccentrica
navi2.ecc_anom0 = 2*atanh(sqrt((navi2.e-1)/(navi2.e+1))*tan(0.5*navi2.Tet(1)));          %in radianti
navi2.mean_anom0 = -navi2.ecc_anom0 + navi2.e*sinh(navi2.ecc_anom0);                   %anomalia media iniziale
navi2.n = sqrt(6.67e-20*venere.M/(navi2.a^3));                                                              %rad/s

%uso il tempo in secondi simulando la posizione ogni mezz'ora
navi2.tempo_volo = (2*abs(navi2.mean_anom0)/navi2.n)/3600;                                %tempo di volo in ore
navi2.tempo_sim = 0:30*60:navi2.tempo_volo*3600;                                                                    %vettore temporale di simulazione
[navi2.dimr, navi2.dimc] = size(navi2.tempo_sim);
navi2.p = plot3(0,0,0,'green .','MarkerSize',10);hold on;
for ind = 1:navi2.dimc
    %prima calcolo l'anomalia vera da quella media
    navi2.manom = navi2.mean_anom0+navi2.n*navi2.tempo_sim(ind);
    [navi2.ecc_anom_ven, navi2.teta_nav] = kepler1(navi2.manom,navi2.e);
    
    if  real(navi2.teta_nav) > pi
        navi2.teta_nav = -2*pi+navi2.teta_nav;
    end
        
    
    navi2.pos_lettura = fix(real(navi2.teta_nav/ris))+0.5*42570+1;
    
    delete(navi2.p);
    navi2.p = plot3(navi2.pos(navi2.pos_lettura,1),navi2.pos(navi2.pos_lettura,2),navi2.pos(navi2.pos_lettura,3),...
                              'magenta .', 'MarkerSize',7); hold on;
    pause(0.1);
end
close ff