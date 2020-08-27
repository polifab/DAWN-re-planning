    %trovo le posizioni della Terra e di Venere, il lancio lo voglio fare
    %quando la terra passa per il proprio punto nodale il 12 Marzo 2018
    %il giorno corrispondente dal primo gennaio 2018 è i = 71
    
    %    tet_lancio = 115*pi/180;             %angolo di sgancio dalla soi della terra
    %   tet_end = 130*pi/180;                   %angolo ingresso soi di venere nel suo sistema di riferimento
    %   teta_long_terra = 3.35*pi/180;
    %   teta_long_ven = -6*pi/180;
    %   vol1.tet_end = tet_end;
    
    vol1 = struct;
    
    tet_lancio = 216*pi/180;             %angolo di sgancio dalla soi della terra
    tet_end = 339.185*pi/180;                   %angolo ingresso soi di venere nel suo sistema di riferimento
    teta_long_terra = 3.0916*pi/180;
    teta_long_ven = 15.1*pi/180;          %day_start = 250 end 329
    vol1.tet_end = tet_end;
   
    %devo risolvere la posizione di venere quando passa per il nodo target
    %e vi passerà: tra il 12 e 13 Aprile 2018 la prima volta (i=104);
    %              tra il 24 e 25 Novembre 2018 la seconda volta(i=328/329)
    %              il 7 Luglio 2019 per la terza volta (i=553)
    %              il 17 Febbraio 2020 per la quarta volta (i = 778)
    t1 = temp(day_start);
    t_target = temp(329);             %ho scelto come data di arrivo il 
    
    
    %parte risoluzione volo
    [a,b,ter,c,dati] = pianeti_pos(J2018+t1,mercurio,venere,terra,marte,riso);
    Ct_s1 = [eye(3,3), ter; 0 0 0 1];
    Ct_s2 = [eulerT(0,dati.tanom_ter,0), zeros(3,1); 0 0 0 1];
    Ct_s3 = [terra.Meulero, zeros(3,1); 0 0 0 1];
    r1_s = [terra.SOI*cos(tet_lancio)*cos(teta_long_terra); terra.SOI*sin(tet_lancio)*cos(teta_long_terra); terra.SOI*sin(teta_long_terra)];           %posizione riferito rispetto frame terra piano lancio
    r1_s = Ct_s1*Ct_s2*Ct_s3*[r1_s;1];                                                      %posizione riferito rispetto frame sole                                                                             
    r1 = [r1_s(1,1); r1_s(2,1); r1_s(3,1)];                                                    %posizione rispetto il sole del lancio!
            %plot3(r1(1,1),r1(2,1),r1(3,1),'magenta .', 'MarkerSize',5)
            %plot3(ter(1,1),ter(2,1),ter(3,1),'black .','MarkerSize',8)

            
    
    [a,b,ter,c,dati] = pianeti_pos(J2018+t_target,mercurio,venere,terra,marte,riso);
    posi_ven = [dati.dati_venere(1,1); dati.dati_venere(2,1); dati.dati_venere(3,1)];
    Ct_s1 = [eye(3,3), posi_ven; 0 0 0 1]*[eulerT(0,dati.tanom_ven,0), zeros(3,1); 0 0 0 1];
    %Ct_s2 = [eulerT(0,dati.tanom_ven,0), zeros(3,1); 0 0 0 1];
    Ct_s3 = [venere.Meulero, zeros(3,1); 0 0 0 1];
    vol1.Ct_ven_sol = Ct_s1*Ct_s3;
    r2_f = [venere.SOI*cos(tet_end)*cos(teta_long_ven); venere.SOI*sin(tet_end)*cos(teta_long_ven); venere.SOI*sin(teta_long_ven)];
    vol1.pos_start = r2_f;
    r2_f = Ct_s1*Ct_s3*[r2_f;1];                                                         %posizione vettore navicella rispetto sole finale SE(3)
    r2 = [r2_f(1,1); r2_f(2,1); r2_f(3,1)];                                                       %vettore posizione navicella rispetto SOLE
    vol1.r2_f_sole = r2;
        %plot3(r2(1,1),r2(2,1),r2(3,1),'magenta .', 'MarkerSize',5)
        %plot3(dati.dati_venere(1,1), dati.dati_venere(2,1), dati.dati_venere(3,1),'green .','MarkerSize',9)
    
    %il codice che ho scelto, crea la traiettoria del satellite giorno dopo
    %giorno, in modo da poter scegliere l'ottima finestra di lancio
    
    u = (6.67e-11*1.9891e30)*1e-9;
    
%     [v1, v2, extremal_distances, exitflag] = lambert(r1', r2', t_target-t1,% 0, u);
    [v1, v2] = lambert(r1', r2', t_target-t1,'pro');
    v1 = v1';
    v2 = v2';
    %[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(ra,va,terra.u_sol);
    [vol1.a,vol1.e,vol1.i_orb,vol1.W,vol1.w,vol1.t_anom0,vol1.p,vol1.ra,vol1.rp] = rv_to_orb(r1,v1,u);
    
    
    vol1.Meulero = eulerT(vol1.W,vol1.w,vol1.i_orb);
    
    vol1.pos = zeros(Nl,3);
    for j = 1:Nl
        vol1.r = vol1.p/(1+vol1.e*cos(teta(j)));               %modulo del raggio posizione
    
        X = vol1.r*cos(teta(j));
        Y = vol1.r*sin(teta(j));
        
        vol1.pos(j,1) = vol1.Meulero(1,1)*X + vol1.Meulero(1,2)*Y;
        vol1.pos(j,2) = vol1.Meulero(2,1)*X + vol1.Meulero(2,2)*Y;
        vol1.pos(j,3) = vol1.Meulero(3,1)*X + vol1.Meulero(3,2)*Y;
    end
    E0 = atan(sqrt((1-vol1.e)/(1+vol1.e))*tan(vol1.t_anom0*0.5))*2;
    vol1.mean_anom0 = E0-vol1.e*sin(E0);
    if vol1.mean_anom0 < 0
        vol1.mean_anom0 = 2*pi + vol1.mean_anom0;
    end
    
    vol1.v_init = v1;
    vol1.pos_init = r1;
    clear tet_lancio t1 a b ter c Ct_s r1_s tet_emd posi_ven Ct_f r2_f v1 v2 r1 r2 X Y E0
    
    %v_orbitN = zeros(400,6);
    %v_orbitN(i,1) = v1(1,1);v_orbitN(i,2) = v1(2,1);v_orbitN(i,3) = v1(3,1);
    %v_orbitN(i,4) = v2(1,1);v_orbitN(i,5) = v2(2,1);v_orbitN(i,6) = v2(3,1);
    %v_orbitN(i,7) = sqrt(dot(v1,v1));
    %v_orbitN(i,8) = sqrt(dot(v2,v2));
    