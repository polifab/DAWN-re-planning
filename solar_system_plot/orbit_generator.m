%definizione costanti universali
G = 6.67e-11;            %in forma SI
M_sun = 1.9891e30;       %[kg]
AU = 149597870.7;        %[km]

% i seguenti script restituiscono delle struct per ogni pianeta, nelle
% quali sono contenuti gli elementi orbitali di ciascun pianeta e un
% vettore delle posizioni lungo l'ellissi
% (l'anomalia vera è definita più avanti)

Earth_orb_el;
Mars_orb_el;
Vesta_orb_el;
Ceres_orb_el;