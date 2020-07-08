%% CODICE RIENTRO SOI TERRA
% Ritorno orbita di parcheggio di 200km (rp dell'iperbole deve essere di
% 6571 km!

%% Codice creazione plot orbita iperboilica di rientro
final = struct;
final.pos_inf_min = vol4.pnav;                                                                                  %salvo la posizione della navicella rif_SOLE
final.pos_terra = ter_t;                                                                                               %salvo la posizione della terra rif_SOLE
final.vel_terra = datiSoi.dati_terra(:,2);                                                                    %salvo la velocità della terra rif_SOLE
final.ans = sqrt(terra.u_sol*1e-9/r.pA);
final.vel_nav = r.MeuleroB*[-final.ans*sin(vol4.t_anom); ...
                                                  final.ans*(vol4.e+cos(vol4.t_anom)); ...
                                                  0];                                                                            %velocità navicella rif_SOLE
                                              
final.C1 = [eye(3), final.pos_terra; 0 0 0 1];
final.C2 = [eulerT(0,datiSoi.tanom_ter,0), zeros(3,1); 0 0 0 1];
final.C3 = [terra.Meulero, zeros(3,1); 0 0 0 1];
final.Ceulero = final.C1*final.C2*final.C3;                                                              %matrice trasformazione coordinate terra->sole

%calcolo posizione navicella rispetto frame terra
final.ans = inv(final.Ceulero)*[final.pos_inf_min; 1];
final.pos_inf_min = [final.ans(1,1); final.ans(2,1); final.ans(3,1)];                           %vettore posizione navicella rispetto frame terra

%calcolo velocità navicella rispetto frame terra
final.vel_inf_min = final.vel_nav - final.vel_terra;
final.ans = inv(final.Ceulero)*[final.vel_inf_min; 0];
final.vel_inf_min = [final.ans(1,1); final.ans(2,1); final.ans(3,1)];                           %vettore velocità navicella rispetto frame terra

%ora che ho il necessario, posso creare l'orbita interna alla soi della
%terra
final.u_ter = 6.67e-20*terra.M;

[final.a,final.e, final.i_orb, final.W,final.w, ...
 final.t_anom0, final.p, final.ra, final.rp] = rv_to_orb(final.pos_inf_min,final.vel_inf_min,final.u_ter);
if (final.t_anom0 > pi)
    final.t_anom0 = (2*pi-final.t_anom0);
end

final.Meulero = eulerT(final.W, final.w, final.i_orb);                                                      %calcolo matrice di rotazione

final.Tet = -final.t_anom0:riso:final.t_anom0;
[~,final.sizeC] = size(final.Tet);
final.orb = zeros(final.sizeC,3);

for j = 1:final.sizeC
    final.r = final.p/(1+final.e*cos(final.Tet(j)));
    
    final.orb(j,:) = (final.Meulero*[final.r*cos(final.Tet(j)); final.r*sin(final.Tet(j)); 0])';               %calcolo della posizione vettoriale    
end

%apro nuova finestra di plot e plotto l'orbita
final.F3 = figure('Name','Rientro Terra', 'Units','pixel','Position',[0 0 800 600]);
final.plotOrb = plot3(final.orb(:,1), final.orb(:,2), final.orb(:,3), 'magenta');hold on                               %plot orbita
final.plotTer = plot3(0,0,0, 'black .', 'MarkerSize', 18);hold on                                                                %plot terra

%% Seconda parte di codice
final.RorbC = final.rp;                                                         %salvo il raggio dell'orbita circolare
final.MeuleroC = eulerT(0,0,final.i_orb);                           %creo la matrice di rotazione
final.TetC = 0:riso:2*pi;

final.sizeC = size(final.TetC,2);
final.orbC = zeros(final.sizeC,3);

for j = 1:final.sizeC
    final.r = final.RorbC;
    
    final.orbC(j,:) = (final.Meulero*[final.r*cos(final.TetC(j)); final.r*sin(final.TetC(j)); 0])';
end

final.plotOrbC = plot3(final.orbC(:,1), final.orbC(:,2), final.orbC(:,3), 'magenta');hold on

%% Simulazione volo

if simula == 1
    simulazione_rientroTerra
end

flag_volo = 5;