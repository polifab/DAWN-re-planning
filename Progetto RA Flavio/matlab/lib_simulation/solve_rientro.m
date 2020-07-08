%% codice che calcola il rientro da venere alla terra
%idea è di calcolarsi giorno a giorno la possibilità di una finestra di
%rientro andando a confrontare la posizione della terra, il rientro è
%un'Hohmann e per questo devo trovare una traiettoria ellittica per il
%rientro, considerando anche però la fuga e il rientro su orbite
%iperboliche

%j = 236;                                                                                                                            %giorno di partenza relativo
%day_start3 = 587+j;                                                                                              %giorno di partenza assoluto
%r = struct;                                                                                                                         %struttura volo rientro
%r.teta_exit_venere = 89*pi/180;                                                                                        %angolo uscita piano orbitale venere

%[~,r.ven_t,r.ter_t,~,r.dati_start] = pianeti_pos(J2018+587+j,mercurio,venere,terra,marte,riso);

r.Ceulero_venere = [venere.Meulero, r.ven_t; 0 0 0 1];                                                 %trasformazioni coordinate venere->sole nel giorno di uscita soi

r.pos_start = r.ven_t;                                                                                                       %salvo posizione di partenza

r.angolo = atan2(r.ven_t(2,1),r.ven_t(1,1));                                                                    %si commette un errore sicuro nel calcolo, lo dovrò correggere con un offset
if r.angolo < 0
    r.angolo = r.angolo + 2*pi;
end

%ora cerco la possibile posizione della terra, per conoscere la distanza
%dell'afelio

r.index = search_posT(r.angolo,terra.pos);                                       %calcolo l'indice del vettore posizione da leggere

%calcolo parametri orbitali
r.rp1 = norm(r.ven_t);
r.ra1 = norm(terra.pos(r.index,:)) + r.correzione_pos;                      %la correzione della posizione è dovuto all'approssimazione del calcolo dell'angolo di posizione di venere
r.a1 = 0.5*(r.rp1+r.ra1);

r.semiPeriod = sqrt(pi*pi*(r.a1^3)/(6.67e-20*1.9891e30))/(24*3600);

%[~,~,r.ter_t,~,~] = pianeti_pos(J2018+587+j+fix(r.semiPeriod),mercurio,venere,terra,marte,riso);

r.pos_fin = terra.pos(r.index,:)';                                                        %salvo il punto di arrivo, uso questo punto, è più preciso
r.u_sol = 1.9891e30*6.67e-20;


%% codice calcolo uscita venere
r.u_ven = 6.67e-20*venere.M;
r.v_peri_req = sqrt((r.u_sol/r.a1)*(r.ra1/r.rp1));                                     %modulo velocità richiesta per raggiungere la terra
r.vel_ven = r.dati_start.dati_venere(:,2);                                                 %velocità del pianeta venere
r.v_inf_venere_norm = r.v_peri_req-sqrt(r.u_sol/r.rp1);                            %modulo velocità uscita soi venere richiesta
r.a_ip_venere = -r.u_ven/(r.v_inf_venere_norm)^2;
%sono ad energia minima, parto dall'orbita di parcheggio di 500km
r.rp_ip_venere = 6051.8+500;
r.e_ip_venere = 1-r.rp_ip_venere/r.a_ip_venere;                             %calcolo ellitticità uscita
r.p_ip_venere = r.a_ip_venere*(1-r.e_ip_venere^2);

%codice soluzione 3D, idea è di uscire dalla soi sullo stesso piano di
%venere
r.pos_inf_venere = venere.SOI*[cos(r.teta_exit_venere); sin(r.teta_exit_venere);0];                             %posizione uscita soi rispetto frame venere
r.vel_inf_venere = r.v_inf_venere_norm*[cos(r.teta_exit_venere); sin(r.teta_exit_venere);0];               %velocità uscita soi rispetto frame venere

%riferisco tutto rispetto il sole
r.pos_inf_venere = r.Ceulero_venere*[r.pos_inf_venere; 1];
r.pos_inf_venere = [r.pos_inf_venere(1,1); r.pos_inf_venere(2,1); r.pos_inf_venere(3,1)];
r.vel_inf_venere = r.Ceulero_venere*[r.vel_inf_venere; 0];
r.vel_inf_venere = [r.vel_inf_venere(1,1); r.vel_inf_venere(2,1); r.vel_inf_venere(3,1)];
r.vel_nav_start = r.vel_inf_venere + r.vel_ven;

%codice generazione orbita interna
r.teta_inf_venere = acos((1/r.e_ip_venere)*(r.p_ip_venere/venere.SOI -1));
%correzione angolo negativo
if r.teta_inf_venere < 0
    r.teta_inf_venere = r.teta_inf_venere+2*pi;
end

r.teta_off_venere = r.teta_exit_venere - r.teta_inf_venere;
r.Meuler_ip_venere = eulerT(0,r.teta_off_venere,0);

r.tet_ip_v = 0:riso:r.teta_inf_venere;
[~,r.sizec] = size(r.tet_ip_v);
r.pos_ip_ven = zeros(r.sizec,3);

%creo orbita uscita iperbolica
for j = 1:r.sizec
    r.r = r.p_ip_venere/(1+r.e_ip_venere*cos(r.tet_ip_v(j)));
    
    X = r.r*cos(r.tet_ip_v(j));
    Y = r.r*sin(r.tet_ip_v(j));
    
    r.pos_ip_ven(j,:) = r.Meuler_ip_venere*[X; Y; 0];
end

%creo orbita parcheggio circolare
r.posC_v = zeros(Nl,3);
for j = 1:Nl
    r.r = r.rp_ip_venere;               %modulo del raggio posizione
    
    X = r.r*cos(teta(j));
    Y = r.r*sin(teta(j));
       
    r.posC_v(j,1) = X;
    r.posC_v(j,2) = Y;
    r.posC_v(j,3) = 0;
end

r.Fig = figure('Name','Uscita da Venere', 'Units','pixel','Position',[0 0 1024 768])

r.plot_parc = plot3(r.posC_v(:,1),r.posC_v(:,2),r.posC_v(:,3),'magenta');hold on;
r.plot_in = plot3(r.pos_ip_ven(:,1),r.pos_ip_ven(:,2),r.pos_ip_ven(:,3),'magenta');hold on;
plot3(0,0,0,'green .','MarkerSize',15);hold on;
xlabel('x(km)');
ylabel('y(km)');
zlabel('z(km)');
axis equal;

%simulazione volo uscita soi
if simula == 1
    simulazione_exit_venere
end


close (r.Fig)



%% codice generazione orbita interstellare
[r.aA,r.eA,r.i_orbA,r.WA,r.wA,r.t_anom0A,r.pA,r.raA,r.rpA] = rv_to_orb(r.pos_inf_venere, r.vel_nav_start, r.u_sol);

r.MeuleroA = eulerT(r.WA,r.wA,r.i_orbA);
r.posA = zeros(Nl,3);

for j = 1:Nl
    r.r = r.pA/(1+r.eA*cos(teta(j)));               %modulo del raggio posizione
    
    X = r.r*cos(teta(j));
    Y = r.r*sin(teta(j));
       
    r.posA(j,1) = r.MeuleroA(1,1)*X + r.MeuleroA(1,2)*Y;
    r.posA(j,2) = r.MeuleroA(2,1)*X + r.MeuleroA(2,2)*Y;
    r.posA(j,3) = r.MeuleroA(3,1)*X + r.MeuleroA(3,2)*Y;
end

r.plotOrb = plot3(r.posA(:,1),r.posA(:,2),r.posA(:,3),'magenta');

vol4 = struct;
%calcolo dell'anomalia zero
 [vol4.p_start1, vol4.ind_start1,...
  vol4.p_start2, vol4.ind_start2] = punto_cont(r.pos_inf_venere',r.posA,1,1,1,62832);

vol4.appoggio = vol4.ind_start2*riso;
vol4.m_anom0 = vol4.appoggio - r.eA*sin(vol4.appoggio);         %non parto preciso dal periasse quando sono uscito dalla soi

vol4.a = r.aA;                                                                                %salvo valore aperiasse
vol4.e = r.eA;                                                                                %salvo valore ellitticità
%simulazione volo, in questa simulazione ci devo inserire il riconoscimento
%della posizione, idealmente so che l'orbita è inclinata quindi devo
%attuare una correzione orbitale in modo da poter intercettare la terra
%uso la semplice conoscenza della i (inclinazione
%orbitale di venere), quindi applico un delta_vel di correzione)
%per ora creo solo le orbite, la simulazione del volo è lasciata risolvere
%alla funzione simulate_solar_system

r.MeuleroB = eulerT(r.WA,r.wA,0);                                           %ecco la correzione orbitale

r.posB = zeros(Nl,3);

for j = 1:Nl
    r.r = r.pA/(1+r.eA*cos(teta(j)));               %modulo del raggio posizione
    
    X = r.r*cos(teta(j));
    Y = r.r*sin(teta(j));
       
    r.posB(j,1) = r.MeuleroB(1,1)*X + r.MeuleroB(1,2)*Y;
    r.posB(j,2) = r.MeuleroB(2,1)*X + r.MeuleroB(2,2)*Y;
    r.posB(j,3) = r.MeuleroB(3,1)*X + r.MeuleroB(3,2)*Y;
    
end
acca = 1;
