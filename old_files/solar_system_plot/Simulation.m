%% Simulation 3D of orbits and planets

close all
clear all

%definizione risoluzione della simulazione
risoluzione = 0.0001;               %risoluzione dell'orbita     
Tet = 0:risoluzione:2*pi;           %vettore angoli teta da 0 a 2pi
N = fix(2*pi/risoluzione)+1;      	%indice

%generazione della geometria delle orbite dei pianeti e plot 3d
Fs =figure('Name','Sistema solare', 'Units','pixel','Position',[0 0 800 600]);
orbit_generator;

% generazione dell'orbita del satellite da Vesta a Ceres
% (gli elementi orbitali sono stati calcolati con Lambert)

% qui va inserito il codice del calcolo della traiettoria

%i parametri necessari sono [a e W w i] da decidere come passarli
satellite_orb_el;

%generazionee plot delle posizioni dei pianeti nel tempo
plot_planets;

