    %in questa funzione voglio studiare la seconda parte della rotta,
    %l'idea usata è quella della navicella che è passata per venere con un
    %fly-by e voglio rincontrarla nello stesso nodo dopo un periodo di
    %volo
    
    t1 = 553;
    [a,b,ter,c] = pianeti_pos(J2018+t1,mercurio,venere,terra,marte,ris);
    r1 = [ter(1,1); ter(2,1); ter(3,1)];
    
    %la posizione target è quando venere è nel suo nodo
    r2 = [venere.pos(21843,1);venere.pos(21843,2);0];
    
    %devo risolvere la posizione di venere quando passa per il nodo target
    %e vi passerà: tra il 12 e 13 Aprile 2018 la prima volta (i=104);
    %              tra il 24 e 25 Novembre 2018 la seconda volta(i=308/309)
    %              il 7 Luglio 2019 per la terza volta (i=553)
    
    t_target = 553;             %ho scelto come data di arrivo il 7 Luglio 2019
    
    %il codice che ho scelto, crea la traiettoria del satellite giorno dopo
    %giorno, in modo da poter scegliere l'ottima finestra di lancio
    
    u = (6.67e-11*1.9891e30)*1e-9;
    
    [v1, v2, extremal_distances, exitflag] = lambert(r1', r2', t_target-t1, 2, u);
    v1 = v1';
    v2 = v2';
    %[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(ra,va,terra.u_sol);
    [a,e,i_orb,W,w,t_anom,p,ra,rp] = rv_to_orb(r1,v1,u);
    
    
    M = eulerT(W,w,i_orb);
    
    pos = zeros(N,3);
    for i = 1:N
        r = p/(1+e*cos(teta(i)));               %modulo del raggio posizione
    
        X = r*cos(teta(i));
        Y = r*sin(teta(i));
        
        pos(i,1) = M(1,1)*X + M(1,2)*Y;
        pos(i,2) = M(2,1)*X + M(2,2)*Y;
        pos(i,3) = M(3,1)*X + M(3,2)*Y;
    end
    plo = plot3(pos(:,1),pos(:,2),pos(:,3),'magenta');hold on;
