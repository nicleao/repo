function PrintParada(x,crit,Peso,Pcg,S,AR,CD_corrida,CL_corrida,MTOW)
%Critérios de Parada (Convencional - Nícolas 04/06/21):
%0 : Terminou de rodar
%1 : Problema com função cordas
%2 : Restrição Geométrica
%3 : Seção muito pequena
%4 : Razão de Aspecto Baixa
%5 : Razão de Aspecto da EH
%6 : Nervura com altura superior a 10cm
%7 : Possível tail-strike
%8 : Volume de cauda horizontal
%9 : Margem Estática
%10 : Cm0 não positivo
%11 : Instável Longitudinal
%12 : Instável Látero-Direcional
%13 : CL na Corrida muito baixo
%14 : Estol de ponta
%15 : Incidência Estola
%16 : Sem Rotate
%17 : Gradiente de Subida Negativo
%18 : Gradiente De Subita Lamentável =(
%19 : Deflexão do Profundor
%20 : Erro Aleatório do AVL (Try/Catch)

cargapaga=MTOW-Peso;
fileIDBOM = fopen('backup.txt','a');
fprintf(fileIDBOM,'\n%s\n\nx = [ %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f]',datetime,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11));
fprintf(fileIDBOM,'\n\nMTOW = %.4f <<<\n\nPeso = %.4f\n\nCarga Paga = %.4f\n\nPosiçao do CG [raiz] = %.4f',MTOW,Peso,cargapaga,Pcg);
fprintf(fileIDBOM,'\n\n Área: %.4f | AR: %.4f | CDcorrida: %.4f | CLcorrida: %.4f | Parada: %f',S,AR,CD_corrida,CL_corrida,crit);
fprintf(fileIDBOM,'\n------------------------------------------------------------------------------------------------------------------------------');
fclose all;
end