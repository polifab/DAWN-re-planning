
if norm(vol2.pnav-ven_t)-venere.SOI <= 50000 && i >= day_start2+2
    disp(["sono nella SOI di venere"])
    %avvio il codice di calcolo del flyby, con controllo immissione orbita
    %attorno venere
    solve_sim_venus
end