

%varizione vincoli
%OPZIONALE
%idee variazioni vincoli  
% {+:facili da cambiare
% -:noiose}
% #myGraph
% +stocasticità: doppia, riga
% +tipo grafo: ciclico(agN,stoc), albero(agN,stoc), binomiale(agN,stoc,prob)
% #progeg
% +range valori (P_,W
% -creazione matrice P,c,d,D
% +cambiare UB(ones), LB(zeroes), H (eyes), b(ones)
% #dualsub
% (+ mettere b1=b, bk=0, k!=1)
% --altri vincoli non utilizzati
% variazione numero agenti
    for agN=2:5
        pause(0)
        %variazione dimensione agenti
        for agni=agN+1:agN+4
            pause(0)
        end
    end
    
    
    DATI DA SALVARE
    
    costi (primalRA,dualRa, (dual e primal del più caratteristico)
    errore convergenza( uno for:agN, uno for:agni)=errore(lambda)
    
    salvare dati: save(filename,variables)
    filename=strcat(partevariabile,'dcs22')