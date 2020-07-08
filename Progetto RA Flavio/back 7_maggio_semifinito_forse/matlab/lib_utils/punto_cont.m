%funzione usata per la ricerca del punto di contatto tra due orbite, do
%anche gli indici di inizio e fine di ricerca

function [p_contatto, i_start, p2_contatto, i_start2] = punto_cont(vect1,vect2,in1,fin1,in2,fin2)
    pos_lettura1 = 0;
    pos_lettura2 = 0;
    dist = 1e15;                                                                                         %assegnamento random per garantire l'inizio dell'algoritmo
    for i = in1:fin1
        for j = in2:fin2
            norma = norm(vect1(i)-vect2(j));
            if norma < dist
                dist = norma;
                pos_lettura1 = i;
                pos_lettura2 = j;
            end
        end
    end
    
    p_contatto = vect1(pos_lettura1,:)';
    i_start = pos_lettura1;
    p2_contatto = vect2(pos_lettura2,:)';
    i_start2 = pos_lettura2;
end