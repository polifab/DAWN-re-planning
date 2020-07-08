%in vol1 ho il vettore pnav che mi dice la posizione spaziale della
%navicella
%in datiSoi ho il vettore dati_venere dove la prima colonna è la posizione
%di venere e la seconda è la velocità vettore di venere
%in più ho il vettore ter_t che è la posizione della terra nel giorno t
disp([num2str(norm(vol1.pnav-ven_t))])

if norm(vol1.pnav-ter_t)-terra.SOI <= 10000              %il 5000 è l'errore accettabile di misura
    disp(["sono nella SOI della terra"])
end
if norm(vol1.pnav-ven_t)-venere.SOI <= 10000
    disp(["sono nella SOI di venere"])
    %avvio il codice di calcolo del flyby
    fly_by
end