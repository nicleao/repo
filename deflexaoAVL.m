function [de] = deflexaoAVL(alpha,beta,file_aviao)

ID1 = ('temp_commands.dat');
ID2 = ('Ftotal_deflexao.dat');

fileID = fopen(ID1,'w');

fprintf(fileID,  'load %s \n',file_aviao);
fprintf(fileID,  'oper \n');
fprintf(fileID,  'a \na \n');
fprintf(fileID,  '%.4f\n',alpha); 
fprintf(fileID,  'b \nb \n');
fprintf(fileID,  '%.4f\n',beta); 
fprintf(fileID,  'D2 \nPM \n0 \n');
fprintf(fileID,  'X \n');
fprintf(fileID,  'FT \n%s \n',ID2);
fprintf(fileID,  'O\n');
fclose('all');

exec = strcat('avl.exe < ',ID1);         % Execução de comandos no AVL
[~,~] = dos(exec);
delete(ID1)

fileID2 = fopen(ID2,'r');
if fileID2 < 3
    fprintf('Aqui o miseravel: %s, %s\n', file_aviao, ID2);
    de=NaN;
    return
end

frewind(fileID2)
de = textscan(fileID2,'   elevator        =   %f','Headerlines',30);
de = de{1};
fclose(fileID2);

delete(ID2)
end