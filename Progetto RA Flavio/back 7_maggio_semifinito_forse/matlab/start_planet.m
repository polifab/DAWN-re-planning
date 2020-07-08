%variabili globali
close all
clear all

%flag usati per plottare i nodi e le simulazioni della navicella
p_nod = 0;
simula = 0;

risoluzione = 0.0001;                                                                                           %meglio non toccare questo valore
Fs =figure('Name','Sistema solare', 'Units','pixel','Position',[0 0 800 600]);

Tet = 0:risoluzione:2*pi;                    %vettore angoli teta
N = fix(2*pi/risoluzione)+1;                 %indice


%generazione orbite, posizione pianeti e moto di essi da J2000
[Mercurio,Venere,Terra,Marte,p] = gener_orb(N,Tet,p_nod);

[Mercurio, Venere, Terra, Marte] = simulate_solar_system(Mercurio,Venere,Terra,Marte,N,Tet,risoluzione,simula);