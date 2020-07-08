%% PRIMA PARTE DEL CODICE

vol3.datiSoi = datiSoi;                            %struttura dati pos/vel venere
vol3.pos_init = vol2.pnav;                      %posizione iniziale ingresso SOI navicella
 
vol3.day_init = i;                                     %salvo il giorno di ingresso del satellite nell'orbita del pianeta

vol3.Vpmin = [vol3.datiSoi.dati_venere(1,2); vol3.datiSoi.dati_venere(2,2); vol3.datiSoi.dati_venere(3,2)];        %salvo la velocità di venere all'inizio del volo
u_sol = venere.u_sol*1e-9;                                                                                          %mi rispetto sole
pippo = sqrt(u_sol/vol1.p);
vol3.Vsmin = vol2.Meulero*[-pippo*sin(vol2.teta_nav); pippo*(vol2.e+cos(vol2.teta_nav)); 0];         %velocità navicella inizio volo rispetto sole
vol3.Vinfmin = vol3.Vsmin-vol3.Vpmin;                                                                                                   %velocità ingresso Soi (coordinate rispetto sole)

vol3.Ct_s1 = [eye(3,3), vol3.datiSoi.dati_venere(:,1); 0 0 0 1]*[eulerT(0,vol3.datiSoi.tanom_ven,0), zeros(3,1); 0 0 0 1];
vol3.Ct_s3 = [venere.Meulero, zeros(3,1); 0 0 0 1];
vol3.Ct_s = vol3.Ct_s1*vol3.Ct_s3;                                                                                                  %trasformazione coordinate frame terra->sole
vol3.Vinfmin = inv(vol3.Ct_s)*[vol3.Vinfmin; 0];
vol3.Vinfmin = [vol3.Vinfmin(1,1); vol3.Vinfmin(2,1); vol3.Vinfmin(3,1)];                                     %vettore velocità navicella ingresso soi (rispetto frame venere)

vol3.posmin = inv(vol3.Ct_s)*[vol3.pos_init; 1];
vol3.posmin = [vol3.posmin(1,1); vol3.posmin(2,1); vol3.posmin(3,1)];                                           %posizione ingresso soi (coordinate rispetto frame venere)
vol3.u_ven = 6.67e-20*venere.M;                                                                                     

[vol3.a,vol3.e,vol3.i_orb,vol3.W,vol3.w,vol3.t_anom0,vol3.p,vol3.ra,vol3.rp] = rv_to_orb(vol3.posmin,vol3.Vinfmin,vol3.u_ven);   %ho i parametri orbitali

if vol3.t_anom0 > pi
    vol3.t_anom0 = -(2*pi-vol3.t_anom0);
end

vol3.Tet = vol3.t_anom0:riso:-vol3.t_anom0;

vol3.Meuler = eulerT(vol3.W,vol3.w,vol3.i_orb);
[vol3.sizem,vol3.sizec] = size(vol3.Tet);
vol3.pos_nav = zeros(vol3.sizec,3);
for z = 1:vol3.sizec
    vol3.r = vol3.p/(1+vol3.e*cos(vol3.Tet(z)));
    vol3.pos_fly = vol3.r*cos(vol3.Tet(z));
    vol3.pos_fly = [vol3.pos_fly; vol3.r*sin(vol3.Tet(z)); 0];              %posizione navicella sistema di riferimento iperbole
    
    vol3.pos_navA(z,:) = vol3.Meuler*vol3.pos_fly;                              %posizione navicella sistema di riferimento venere
end

%plot prima orbita, quella di ingresso
vol3.F2 = figure('Name','Simulazione Venere', 'Units','pixel','Position',[0 0 800 600]);
vol3.p_orbA = plot3(vol3.pos_navA(:,1),vol3.pos_navA(:,2),vol3.pos_navA(:,3),'magenta -.'); hold on;
vol3.p_orbV = plot3(0,0,0,'black .','MarkerSize',18);
axis equal;
xlabel('x(km)');
ylabel('y(km)');
zlabel('z(km)');

%%
%calcolo e plot della seconda orbita, scelgo di "agganciarmi" prima al
%pianeta e successivamente per passi arrivare all'orbita finale
vol3.r = vol3.p/(1+vol3.e*cos(0));
vol3.pos_fly = vol3.r*cos(0);
vol3.pos_fly = [vol3.pos_fly; vol3.r*sin(0); 0];              %posizione navicella nel periasse dell'iperbole = apoasse ellisse attorno venere
    
vol3.pos_startB = vol3.Meuler*vol3.pos_fly;           %posizione navicella sistema di riferimento venere, cambio in seconda rotta
vol3.rpB_target = 6051.84+500;                                          %periasse target
vol3.pippo = sqrt(venere.M*6.67e-20/vol3.p);
vol3.vpA = vol3.Meuler*[-vol3.pippo*sin(0); vol3.pippo*(vol3.e+cos(0));0];
vol3.vel_test = 1.5;
vol3.vaB = vol3.vpA-vol3.Meuler*[0;vol3.vel_test;0];          %velocità dopo accensione motori, serve per il calcolo dei parametri orbitali

[vol3.aB,vol3.eB,vol3.i_orbB,vol3.WB,vol3.wB,vol3.t_anom0B,vol3.pB,vol3.raB,vol3.rpB] = rv_to_orb(vol3.pos_startB,vol3.vaB,vol3.u_ven);
%loop di autocalcolo della velocità necessaria per arrivare al target
while(vol3.rpB >= vol3.rpB_target)
    vol3.vel_test = vol3.vel_test + 0.0001;
    vol3.vaB = vol3.vpA-vol3.Meuler*[0;vol3.vel_test;0];          %velocità dopo accensione motori, serve per il calcolo dei parametri orbitali
    [vol3.aB,vol3.eB,vol3.i_orbB,vol3.WB,vol3.wB,vol3.t_anom0B,vol3.pB,vol3.raB,vol3.rpB] = rv_to_orb(vol3.pos_startB,vol3.vaB,vol3.u_ven);
end

vol3.TetB = 0:riso:2*pi;
vol3.MeulerB = eulerT(vol3.WB,vol3.wB,vol3.i_orbB);
[vol3.sizemB,vol3.sizecB] = size(vol3.TetB);
vol3.pos_navB = zeros(vol3.sizecB,3);
for z = 1:vol3.sizecB
    vol3.r = vol3.pB/(1+vol3.eB*cos(vol3.TetB(z)));
    vol3.pos_fly = vol3.r*cos(vol3.TetB(z));
    vol3.pos_fly = [vol3.pos_fly; vol3.r*sin(vol3.TetB(z)); 0];              %posizione navicella sistema di riferimento iperbole
    
    vol3.pos_navB(z,:) = vol3.MeulerB*vol3.pos_fly;                              %posizione navicella sistema di riferimento venere
end

vol3.p_orbB = plot3(vol3.pos_navB(:,1),vol3.pos_navB(:,2),vol3.pos_navB(:,3),'magenta --'); hold on;

%%
%calcolo dell'orbita finale
vol3.r = vol3.pB/(1+vol3.eB*cos(0));
vol3.pos_fly = vol3.r*cos(0);
vol3.pos_fly = [vol3.pos_fly; vol3.r*sin(0); 0];                    %posizione navicella nel periasse dell'ellisse = punto orbita cicolare attorno venere
    
vol3.pos_startC = vol3.MeulerB*vol3.pos_fly;                    %posizione navicella sistema di riferimento venere, cambio in terza rotta
vol3.vel_des3 = sqrt(vol3.u_ven/norm(vol3.pos_startC));   %modulo della velocità necessaria per passare da orbita ellittica in quella circolare
vol3.pippo = sqrt(venere.M*6.67e-20/vol3.pB);
vol3.vpB = vol3.MeulerB*[-vol3.pippo*sin(0); vol3.pippo*(vol3.eB+cos(0));0];
vol3.delta_V2 = norm(vol3.vpB)-vol3.vel_des3;
vol3.vC = vol3.vpB-vol3.MeulerB*[0;vol3.delta_V2;0];     %velocità finale

[vol3.aC,vol3.eC,vol3.i_orbC,vol3.WC,vol3.wC,vol3.t_anom0C,vol3.pC,vol3.raC,vol3.rpC] = rv_to_orb(vol3.pos_startC,vol3.vC,vol3.u_ven);

vol3.TetC = 0:riso:2*pi;
vol3.MeulerC = eulerT(vol3.WC,vol3.wC,vol3.i_orbC);
[vol3.sizemC,vol3.sizecC] = size(vol3.TetC);
vol3.pos_navC = zeros(vol3.sizecC,3);
for z = 1:vol3.sizecC
    vol3.r = vol3.pC/(1+vol3.eC*cos(vol3.TetC(z)));
    vol3.pos_fly = vol3.r*cos(vol3.TetC(z));
    vol3.pos_fly = [vol3.pos_fly; vol3.r*sin(vol3.TetC(z)); 0];              %posizione navicella sistema di riferimento iperbole
    
    vol3.pos_navC(z,:) = vol3.MeulerC*vol3.pos_fly;                              %posizione navicella sistema di riferimento venere
end

vol3.p_orbC = plot3(vol3.pos_navC(:,1),vol3.pos_navC(:,2),vol3.pos_navC(:,3),'magenta'); hold on;

%% calcolo dell'orbita target polare, venere è inclinata rispetto il suo piano orbitale di 177.40°
vol3.MeulerT = eulerT(0,0,(-90+177.40)*pi/180);                           %matrice trasformazione coordinate orbitali
vol3.RT = norm(vol3.pos_startC);
vol3.TetT = 0:riso:2*pi;

[~,vol3.sizecT]= size(vol3.TetT);
vol3.pos_navT = zeros(vol3.sizecT,3);
for z = 1:vol3.sizecT
    vol3.pos_fly = vol3.RT*cos(vol3.TetT(z));
    vol3.pos_fly = [vol3.pos_fly; vol3.RT*sin(vol3.TetT(z)); 0];
    
    vol3.pos_navT(z,:) = vol3.MeulerT*vol3.pos_fly;
end
vol3.p_orbT = plot3(vol3.pos_navT(:,1),vol3.pos_navT(:,2),vol3.pos_navT(:,3),'magenta'); hold on;


%% simulazione di volo
if simula == 1
    simulazione_venere
end

%chiudo tutto
close(vol3.F2)

%cancello un pò di cose, inizia la parte di codice che deve progettare il
%rientro sull terra dopo un mese di volo
flag_volo = 3;                                                                                                  %fine di simulazione volo andata
delete(vol2.plot_posN);
delete(vol2.plot_orbN);


%calcolo giorno rientro
day_burn