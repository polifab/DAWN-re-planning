%controllo se entro nella SOI della terra nel rientro

if norm(vol4.pnav-ter_t) <= terra.SOI
    flag_volo = 5;
    solve_SOI_fin
end