function [mer,ven,ter,mar,dati] = pianeti_pos(temp,mercurio,venere,terra,marte,ris)
    tmp = temp;
    risoluzione = ris;
    mean_anom_mer0 = mercurio.m_anom0;
    mean_anom_ven0 = venere.m_anom0;
    mean_anom_ter0 = terra.m_anom0;
    mean_anom_mar0 = marte.m_anom0;
    dati = struct;
    
    %calcolo r(t) di mercurio
    mer = eye(3,1);
    manom = mean_anom_t(mean_anom_mer0,mercurio.u_sol,mercurio.a,tmp);   %mean anomaly of mercury from J2018
    [ecc_anom_mer, teta_mer] = kepler1(manom,mercurio.e);
    pos_lettura = fix(teta_mer/risoluzione)+1;
    mer(1,1) = mercurio.pos(pos_lettura,1);
    mer(2,1) = mercurio.pos(pos_lettura,2);
    mer(3,1) = mercurio.pos(pos_lettura,3);
    
    %calcolo r(t) di venere
    ven = eye(3,1);
    manom = mean_anom_t(mean_anom_ven0,venere.u_sol,venere.a,tmp);       %mean anomaly of venus from J2018
    [ecc_anom_ven, teta_ven] = kepler1(manom,venere.e);
    pos_lettura = fix(teta_ven/risoluzione)+1;
    ven(1,1) = venere.pos(pos_lettura,1);
    ven(2,1) = venere.pos(pos_lettura,2);
    ven(3,1) = venere.pos(pos_lettura,3);
    
    %serve per me
    dati.dati_venere = zeros(3,2);       %posizione e velocità vettori
    dati.dati_venere(1,1) = venere.pos(pos_lettura,1);
    dati.dati_venere(2,1) = venere.pos(pos_lettura,2);
    dati.dati_venere(3,1) = venere.pos(pos_lettura,3);
    acca = sqrt(6.67e-20*1.9891e30/venere.p);
    dati.dati_venere(1,2) = venere.Meulero(1,1)*(-acca)*sin(teta_ven)+venere.Meulero(1,2)*acca*(venere.e+cos(teta_ven));
    dati.dati_venere(2,2) = venere.Meulero(2,1)*(-acca)*sin(teta_ven)+venere.Meulero(2,2)*acca*(venere.e+cos(teta_ven));
    dati.dati_venere(3,2) = venere.Meulero(3,2)*(-acca)*sin(teta_ven)+venere.Meulero(3,2)*acca*(venere.e+cos(teta_ven));
    dati.tanom_ven = teta_ven;
    
    %calcolo r(t) della terra
    ter = eye(3,1);
    manom = mean_anom_t(mean_anom_ter0,terra.u_sol,terra.a,tmp);         %mean anomaly of earth from J2018
    [ecc_anom_ter,teta_ter] = kepler1(manom,terra.e);
    pos_lettura = fix(teta_ter/risoluzione)+1;
    ter(1,1) = terra.pos(pos_lettura,1);
    ter(2,1) = terra.pos(pos_lettura,2);
    ter(3,1) = terra.pos(pos_lettura,3);
    
    %serve per me
    dati.dati_terra = zeros(3,2);       %posizione e velocità vettori
    dati.dati_terra(1,1) = terra.pos(pos_lettura,1);
    dati.dati_terra(2,1) = terra.pos(pos_lettura,2);
    dati.dati_terra(3,1) = terra.pos(pos_lettura,3);
    acca = sqrt(6.67e-20*1.9891e30/terra.p);
    dati.dati_terra(1,2) = terra.Meulero(1,1)*(-acca)*sin(teta_ter)+terra.Meulero(1,2)*acca*(terra.e+cos(teta_ter));
    dati.dati_terra(2,2) = terra.Meulero(2,1)*(-acca)*sin(teta_ter)+terra.Meulero(2,2)*acca*(terra.e+cos(teta_ter));
    dati.dati_terra(3,2) = terra.Meulero(3,2)*(-acca)*sin(teta_ter)+terra.Meulero(3,2)*acca*(terra.e+cos(teta_ter));
    dati.tanom_ter = teta_ter;
    
    %calcolo r(t) di marte
    mar = eye(3,1);
    manom = mean_anom_t(mean_anom_mar0,marte.u_sol,marte.a,tmp);         %mean anomaly of mars from J2018
    [ecc_anom_mar, teta_mar] = kepler1(manom,marte.e);
    pos_lettura = fix(teta_mar/risoluzione)+1;
    mar(1,1) = marte.pos(pos_lettura,1);
    mar(2,1) = marte.pos(pos_lettura,2);
    mar(3,1) = marte.pos(pos_lettura,3);


end