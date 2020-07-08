%script che crea le orbite di uscite dalla prima orbita di parcheggio
%equatoriale fino alla traiettoria iperbolica

[~,~,~,~,vol0.datiSoi] = pianeti_pos(J2018+day_start-1,mercurio,venere,terra,marte,riso);                     %serve per velocità terra giorno uscita soi(velocità rispetto il sole)
vol0.VpP = [vol0.datiSoi.dati_terra(1,2); vol0.datiSoi.dati_terra(2,2); vol0.datiSoi.dati_terra(3,2)];   %salvo la velocità della terra alla fine dell'uscita (rispetto sole)
vol0.VsP = vol1.v_init;                                                                                                                      %salvo la velocità finale della navicella uscita soi (rispetto sole)
vol0.VinfP = vol0.VsP-vol0.VpP;                                                                                                     %velocità uscita naviella soi (coordinate rispetto frame sole)

vol0.Ct_s1 = [eye(3,3), vol0.datiSoi.dati_terra(:,1); 0 0 0 1]*[eulerT(0,vol0.datiSoi.tanom_ter,0), zeros(3,1); 0 0 0 1];
vol0.Ct_s3 = [terra.Meulero, zeros(3,1); 0 0 0 1];
vol0.Ct_s = vol0.Ct_s1*vol0.Ct_s3;                                                                                                  %trasformazione coordinate frame terra->sole
vol0.VinfP = inv(vol0.Ct_s)*[vol0.VinfP; 0];
vol0.VinfP = [vol0.VinfP(1,1); vol0.VinfP(2,1); vol0.VinfP(3,1)];                                                   %velocità uscita soi rispetto frame terra

vol0.pos_fin = vol1.pos_init;                                                                                                               %salvo la posizione finale uscita soi (coordinate rispetto sole)
vol0.pos_fin = inv(vol0.Ct_s)*[vol0.pos_fin; 1];                                                                                %coordinate posizione uscita soi (coordinate rispetto terra!)
vol0.pos_fin = [vol0.pos_fin(1,1); vol0.pos_fin(2,1); vol0.pos_fin(3,1)];
vol0.u_ter = 6.67e-20*terra.M;

[vol0.a,vol0.e,vol0.i_orb,vol0.W,vol0.w,vol0.t_anom0,vol0.p,vol0.ra,vol0.rp] = rv_to_orb(vol0.pos_fin,vol0.VinfP,vol0.u_ter);   %ho i parametri orbitali

vol0.Tet = 0:riso:vol0.t_anom0;

vol0.Meuler = eulerT(vol0.W,vol0.w,vol0.i_orb);
[vol0.sizem,vol0.sizec] = size(vol0.Tet);
for z = 1:vol0.sizec
    vol0.r = vol0.p/(1+vol0.e*cos(vol0.Tet(z)));
    vol0.pos_fly = vol0.r*cos(vol0.Tet(z));
    vol0.pos_fly = [vol0.pos_fly; vol0.r*sin(vol0.Tet(z)); 0];              %posizione navicella sistema di riferimento iperbole
    
    vol0.pos_ter(z,:) = vol0.Meuler*vol0.pos_fly;                              %posizione navicella sistema di riferimento terra
end

vol0.F2 = figure('Name','Orbita iniziale Terra', 'Units','pixel','Position',[0 0 800 600]);
vol0.p_orb = plot3(vol0.pos_ter(:,1),vol0.pos_ter(:,2),vol0.pos_ter(:,3),'magenta'); hold on;
vol0.p_orbV = plot3(0,0,0,'black .','MarkerSize',18);
axis equal;
xlabel('x(km)');
ylabel('y(km)');
zlabel('z(km)');

%codice creazione traiettoria precedente all'uscita
vol0.u_ter = 6.67e-20*terra.M;
vol0.r = vol0.p/(1+vol0.e*cos(0));
vol0.pos_fly = vol0.r*cos(0);
vol0.pos_fly = [vol0.pos_fly; vol0.r*sin(0); 0];                    %posizione navicella nel periasse dell'iperbole = punto orbita circolare terra sgancio    
vol0.pos_startB = vol0.Meuler*vol0.pos_fly;                    %posizione navicella sistema di riferimento terra, cambio in rotta uscita
vol0.vel_des3 = sqrt(vol0.u_ter/norm(vol0.pos_startB));   %modulo della velocità necessaria per passare da orbita ellittica in quella circolare
vol0.pippo = sqrt(vol0.u_ter/vol0.p);
vol0.vpA = vol0.Meuler*[-vol0.pippo*sin(0); vol0.pippo*(vol0.e+cos(0));0];
vol0.delta_V2 = norm(vol0.vpA)-vol0.vel_des3;
vol0.vB = vol0.vpA-vol0.Meuler*[0;vol0.delta_V2;0];     %velocità finale

[vol0.aB,vol0.eB,vol0.i_orbB,vol0.WB,vol0.wB,vol0.t_anom0B,vol0.pB,vol0.raB,vol0.rpB] = rv_to_orb(vol0.pos_startB,vol0.vB,vol0.u_ter);

vol0.TetB = 0:riso:2*pi;
vol0.MeulerB = eulerT(vol0.WB,vol0.wB,vol0.i_orbB);
[vol0.sizemB,vol0.sizecB] = size(vol0.TetB);
vol0.pos_navB = zeros(vol0.sizecB,3);
for z = 1:vol0.sizecB
    vol0.r = vol0.pB/(1+vol0.eB*cos(vol0.TetB(z)));
    vol0.pos_fly = vol0.r*cos(vol0.TetB(z));
    vol0.pos_fly = [vol0.pos_fly; vol0.r*sin(vol0.TetB(z)); 0];              %posizione navicella sistema di riferimento iperbole
    
    vol0.pos_navB(z,:) = vol0.MeulerB*vol0.pos_fly;                              %posizione navicella sistema di riferimento terra
end

vol0.p_orbB = plot3(vol0.pos_navB(:,1),vol0.pos_navB(:,2),vol0.pos_navB(:,3),'magenta'); hold on;

%plotto orbita circolare equatoriale
vol0.TetC = vol0.TetB;
vol0.MeulerC = eulerT(0,0,23.45*pi/180);
[vol0.sizemC,vol0.sizecC] = size(vol0.TetC);
vol0.pos_navC = zeros(vol0.sizecC,3);
vol0.rC = 6378.136+200;                                                                          %distanza di parcheggio orbita equatoriale
for z = 1:vol0.sizecC
    vol0.pos_fly = vol0.rC*cos(vol0.TetC(z));
    vol0.pos_fly = [vol0.pos_fly; vol0.rC*sin(vol0.TetC(z)); 0];              %posizione navicella sistema di riferimento iperbole
    
    vol0.pos_navC(z,:) = vol0.MeulerC*vol0.pos_fly;                              %posizione navicella sistema di riferimento terra
end
vol0.p_orbC = plot3(vol0.pos_navC(:,1),vol0.pos_navC(:,2),vol0.pos_navC(:,3),'magenta'); hold on;

%simulazione volo
if simula == 1
    simulazione_exit_T
end

close(vol0.F2)

vol1.plot_orbN = plot3(vol1.pos(:,1),vol1.pos(:,2),vol1.pos(:,3),'magenta');hold on