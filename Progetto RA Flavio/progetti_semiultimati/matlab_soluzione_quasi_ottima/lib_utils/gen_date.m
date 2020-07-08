function [p] = gen_date(t)

    %conversione data da numero di giorni in data del calendario gregoriano
    data = t+1;
    i_anno = 1;       %settato a 1 vuol dire 2018
    mese = ["Gennaio","Febbraio","Marzo","Aprile","Maggio","Giugno","Luglio", ...
            "Agosto","Settembre","Ottobre","Novembre","Dicembre"];
    giorni =  [31,28,31,30,31,30,31,31,30,31,30,31];
    giorniB = [31,29,31,30,31,30,31,31,30,31,30,31];
    anni = [2018,2019,2020,2021,2022,2023,2024,2025,2026,2027,2028,2029,2030];
    
    i_mese = 1;
    data = data - giorni(i_mese);
    while (data > 0)
        if (i_mese == 12)
            i_anno = i_anno + 1;
            i_mese = 0;
        end
        i_mese = i_mese + 1;
        if ((anni(i_anno)/4 - fix(anni(i_anno)/4)) > 0)         %riconoscimento se siamo in anno bisestile
            data = data - giorni(i_mese);
        else
            data = data - giorniB(i_mese);
        end
    end

    if ((anni(i_anno)/4 - fix(anni(i_anno)/4)) > 0)         %riconoscimento se siamo in anno bisestile
        data = data + giorni(i_mese);
    else
        data = data + giorniB(i_mese);
    end
    
    p = title(data +" "+ mese(i_mese) +" "+ anni(i_anno));hold on;

end