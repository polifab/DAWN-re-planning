%simulazione di volo prima attorno alla terra, sull'orbita di parcheggio,
%poi faccio una prima manovra di cambio traiettoria circolare per poi
%fuggire dalla terra                                        FUNZIONA


%% volo 1             ORBITA PARCHEGGIO EQUATORIALE
%prima parte di volo, prima cerco il punto di contatto tra le due orbite
%circolari, per la prima manovra
%spazzolo il vettore por_navC da 5524 a 5807
%spazzolo il vettore pos_navB da 26969 a 28502

[vol0.s.p_start1, vol0.s.ind_start1,...
 vol0.s.p_start2, vol0.s.ind_start2] = punto_cont(vol0.pos_navC,vol0.pos_navB,5524,5667,54350,57505);           %ho il punto di partenza e l'indice da cui ho letto la posizione nel vettore pos_navC

%simulazione semplice, anomalia eccentrica e anomalia vera coincidono,
%conoscere quindi Teta mi da conoscenza anche di E e di M allo stesso tempo
%da notare che l'angolo di inizio simulazione lo derivo facilmente dalla
%posizione di lettura
vol0.s.t_anom0 = vol0.s.ind_start1*riso;
vol0.s.n = sqrt(6.67e-20*terra.M/(vol0.rC)^3);                                  %moto medio orbita attorno alla terra
vol0.s.tempo_sim = time_sim;                                                                 %velocità simulazione temporale
vol0.s.numero_giri = 4;
vol0.s.t_anom = vol0.s.t_anom0;
vol0.s.temp = 0;

    while vol0.s.t_anom < vol0.s.t_anom0+2*vol0.s.numero_giri*pi
        vol0.s.pos_fly = vol0.rC*cos(vol0.s.t_anom);
        vol0.s.pos_fly = [vol0.s.pos_fly; vol0.rC*sin(vol0.s.t_anom); 0];              %posizione navicella sistema di riferimento iperbole
    
        vol0.s.pos = vol0.MeulerC*vol0.s.pos_fly;
        vol0.s.plot = plot3(vol0.s.pos(1,1),vol0.s.pos(2,1),vol0.s.pos(3,1),'magenta .','MarkerSize',6);hold on;
        pause(0.05);
        delete(vol0.s.plot)
        vol0.s.temp = vol0.s.temp + vol0.s.tempo_sim;
        vol0.s.t_anom = vol0.s.t_anom0 + vol0.s.n*(vol0.s.temp);
    end
    
    
%% volo 2                CAMBIO ORBITA
%so già la posizione da cui parto, lo faccio girare per un pò e poi parte
%per lo spazio
vol0.s.t_anom0 = vol0.s.ind_start2*riso;
vol0.s.t_anom = vol0.s.t_anom0;
[vol0.s.p_start1, vol0.s.ind_start1,~,~] = punto_cont(vol0.pos_navB,vol0.pos_startB,1,62832,1,1);
vol0.s.t_fin = vol0.s.ind_start1*riso;
vol0.s.dteta = vol0.s.t_fin-vol0.s.t_anom0;
vol0.s.temp_in = vol0.s.temp;

while vol0.s.t_anom < vol0.s.t_anom0+2*pi+vol0.s.dteta
    vol0.s.pos_fly = vol0.rC*cos(vol0.s.t_anom);
    vol0.s.pos_fly = [vol0.s.pos_fly; vol0.rC*sin(vol0.s.t_anom); 0];              %posizione navicella sistema di riferimento iperbole
    
    vol0.s.pos = vol0.MeulerB*vol0.s.pos_fly;
    vol0.s.plot = plot3(vol0.s.pos(1,1),vol0.s.pos(2,1),vol0.s.pos(3,1),'magenta .','MarkerSize',6);hold on;
    pause(0.05);
   delete(vol0.s.plot)
   vol0.s.temp = vol0.s.temp + vol0.s.tempo_sim;
   vol0.s.t_anom = vol0.s.t_anom0 + vol0.s.n*(vol0.s.temp-vol0.s.temp_in);
end

%% volo 3               USCITA SOI TERRA
vol0.s.temp_in2 = vol0.s.temp;
vol0.s.m_anom = 0;                                                  %anomalia media uscita
vol0.s.n2 = sqrt(6.67e-20*terra.M/(-vol0.a^3));      %moto medio uscita hyper
vol0.s.tempo_sim2 = time_sim;                                     %nuovo tempo di simulazione, per simulazione più veloce
vol0.s.tet = 0;
vol0.s.M = vol0.s.m_anom;

while vol0.s.tet <= vol0.t_anom0
    [~,vol0.s.tet] = kepler2(vol0.s.M,vol0.e);
    if vol0.s.tet > pi
        vol0.s.tet = vol0.s.tet-2*pi;
    end
    
    vol0.s.r = vol0.p/(1+vol0.e*cos(vol0.s.tet));
    vol0.s.pos_fly = vol0.s.r*cos(vol0.s.tet);
    vol0.s.pos_fly = [vol0.s.pos_fly; vol0.s.r*sin(vol0.s.tet); 0];              %posizione navicella sistema di riferimento iperbole
    
    vol0.s.pos = vol0.Meuler*vol0.s.pos_fly;
    vol0.s.plot = plot3(vol0.s.pos(1,1),vol0.s.pos(2,1),vol0.s.pos(3,1),'magenta .','MarkerSize',6);hold on;
    pause(0.05);
    delete(vol0.s.plot)
    vol0.s.temp = vol0.s.temp+vol0.s.tempo_sim2;
    vol0.s.M = (vol0.s.temp-vol0.s.temp_in2)*vol0.s.n2;
end