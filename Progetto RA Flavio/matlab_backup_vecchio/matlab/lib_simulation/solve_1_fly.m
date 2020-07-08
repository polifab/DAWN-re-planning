    %trovo le posizioni della Terra e di Venere, il lancio lo voglio fare
    %quando la terra passa per il proprio punto nodale il 12 Marzo 2018
    %il giorno corrispondente dal primo gennaio 2018 è i = 71
    t1 = temp(i);
    [a,b,ter,c] = pianeti_pos(J2018+t1,mercurio,venere,terra,marte,ris);
    r1 = [ter(1,1); ter(2,1); ter(3,1)];
    
    %la posizione target è quando venere è nel suo nodo
    r2 = [venere.pos(53260,1);venere.pos(53260,2);0];
    
    %devo risolvere la posizione di venere quando passa per il nodo target
    %e vi passerà: tra il 12 e 13 Aprile 2018 la prima volta (i=104);
    %              tra il 24 e 25 Novembre 2018 la seconda volta(i=328/329)
    %              il 7 Luglio 2019 per la terza volta (i=553)
    %              il 17 Febbraio 2020 per la quarta volta (i = 778)
    
    t_target = 553;             %ho scelto come data di arrivo il 7 Luglio 2019
    
    %il codice che ho scelto, crea la traiettoria del satellite giorno dopo
    %giorno, in modo da poter scegliere l'ottima finestra di lancio
    
    u = (6.67e-11*1.9891e30)*1e-9;
    
    [v1, v2, extremal_distances, exitflag] = lambert(r1', r2', t_target-t1, 1, u);
    v1 = v1';
    v2 = v2';
    %[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(ra,va,terra.u_sol);
    [a,e,i_orb,W,w,t_anom,p,ra,rp] = rv_to_orb(r1,v1,u);
    
    
    M = eulerT(W,w,i_orb);
    
    pos = zeros(Nl,3);
    for j = 1:Nl
        r = p/(1+e*cos(teta(j)));               %modulo del raggio posizione
    
        X = r*cos(teta(j));
        Y = r*sin(teta(j));
        
        pos(j,1) = M(1,1)*X + M(1,2)*Y;
        pos(j,2) = M(2,1)*X + M(2,2)*Y;
        pos(j,3) = M(3,1)*X + M(3,2)*Y;
    end
    plot_orbN = plot3(pos(:,1),pos(:,2),pos(:,3),'magenta');hold on;
    
    %v_orbitN = zeros(400,6);
    v_orbitN(i,1) = v1(1,1);v_orbitN(i,2) = v1(2,1);v_orbitN(i,3) = v1(3,1);
    v_orbitN(i,4) = v2(1,1);v_orbitN(i,5) = v2(2,1);v_orbitN(i,6) = v2(3,1);
    v_orbitN(i,7) = sqrt(dot(v1,v1));
    v_orbitN(i,8) = sqrt(dot(v2,v2));
    