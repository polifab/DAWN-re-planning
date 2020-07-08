j = 249;                                                                                                                            %giorno di partenza relativo
day_start3 = 587+j;                                                                                              %giorno di partenza assoluto
r = struct;                                                                                                                         %struttura volo rientro
r.teta_exit_venere = 89*pi/180;                                                                                        %angolo uscita piano orbitale venere
r.correzione_pos =  13.9245e6;                                                                                         %valore ottenuto da simulazioni per correggere il moto della navicella

[~,r.ven_t,r.ter_t,~,r.dati_start] = pianeti_pos(J2018+587+j,mercurio,venere,terra,marte,riso);