%% Observações Gerais
% Função para cálculo dos coeficientes aerodinâmicos dos componentes de um convencional
% a partir da asa e da EH em .avl e um range de ângulos.
% Obs1.: Só funciona pra EH e EV retangulares
% Obs2.: Só funciona pra asas iguais

%%
function [CL,CD,Cm] = coefsAVL(alpha1,alpha2,beta,file_asa,file_eh)
tic

% Coeficientes da Asa
CLw  = NaN(1,10);
CDw  = NaN(1,10);
Cmw  = NaN(1,10);

% Coeficientes da EH
CLeh = NaN(1,10);
CDeh = NaN(1,10);
Cmeh = NaN(1,10);

i = 1;
%% Pré-AVL

for alpha = linspace(alpha1,alpha2,alpha2-alpha1+1) % Considera 10 ângulos no intervalo de alpha1 a alpha2
    
    %Asa
    TC1 = strcat('temp_commands_',int2str(alpha),'.dat');
    IDw = strcat('polarAVL_',int2str(alpha),'.dat');
    
    fileID = fopen(TC1,'w');

    fprintf(fileID,  'load %s \n',file_asa);
    fprintf(fileID,  'oper \n');
    fprintf(fileID,  'a \na \n');
    fprintf(fileID,  '%.4f\n',alpha); 
    fprintf(fileID,  'b \nb \n');
    fprintf(fileID,  '%.4f\n',beta); 
    fprintf(fileID,  'X \nFN \n');
    fprintf(fileID,  '%s \n',IDw);
    fprintf(fileID,  'o\n\n\n\n\n\n\n\n\n\n');
    fclose('all');
    
    %EH
    TC2 = strcat('temp_commands_EH_',int2str(alpha),'.dat');
    IDh = strcat('polarAVL_EH_',int2str(alpha),'.dat');
    
    fileID = fopen(TC2,'w');
    
    fprintf(fileID,  'load %s \n',file_eh);
    fprintf(fileID,  'oper \n');
    fprintf(fileID,  'a \na \n');
    fprintf(fileID,  '%.4f\n',(alpha)); 
    fprintf(fileID,  'b \nb \n');
    fprintf(fileID,  '%.4f\n',beta); 
    fprintf(fileID,  'X \nFN \n');
    fprintf(fileID,  '%s \n',IDh);
    fprintf(fileID,  'o\n\n\n\n\n\n\n\n\n\n');
    fclose('all');
end

%% AVL

parfor alpha = linspace(alpha1,alpha2,alpha2-alpha1+1)
    exec = strcat('avl.exe < ',strcat('temp_commands_',int2str(alpha),'.dat'));
    [~,~] = dos(exec);
end
parfor alpha = linspace(alpha1,alpha2,alpha2-alpha1+1)
    exec = strcat('avl.exe < ',strcat('temp_commands_EH_',int2str(alpha),'.dat'));
    [~,~] = dos(exec);
end

%% Pós-AVL

for alpha = linspace(alpha1,alpha2,alpha2-alpha1+1)
    %Asa
    IDw = strcat('polarAVL_',int2str(alpha),'.dat');
    fileIDw = fopen(IDw,'r');
    dataw1 = textscan(fileIDw,'Sref =   %f       Cref =    %f   Bref =    %f','Headerlines',4);

    dataw2 = textscan(fileIDw,[
        '%f %f %f %f %f %f %f %f %f %f   Wing\r\n'...
        '%f %f %f %f %f %f %f %f %f %f   Wing (YDUP)\r\n'],'Headerlines',3);

    %EH
    IDh = strcat('polarAVL_EH_',int2str(alpha),'.dat');
    fileIDh = fopen(IDh,'r');
    datah1 = textscan(fileIDh,'Sref =   %f       Cref =    %f   Bref =    %f','Headerlines',4);

    datah2 = textscan(fileIDh,[
        '%f %f %f %f %f %f %f %f %f %f   EH\r\n'...
        '%f %f %f %f %f %f %f %f %f %f   EH (YDUP)\r\n'],'Headerlines',3);
    



%% CL, CD e Cm das polares
    
    % Asa
    CLw(i) = (dataw2{3} + dataw2{13});
    CDw(i) = (dataw2{4} + dataw2{14});
    Cmw(i) = (dataw2{5} + dataw2{15});
      
    % EH
    CLeh(i) = datah2{3} + datah2{13};
    CDeh(i) = datah2{4} + datah2{14};
    Cmeh(i) = datah2{5} + datah2{15};
    
    i = i+1;
end

%% Ajeitando os outputs

ang = linspace(alpha1,alpha2,alpha2-alpha1+1);

CL.w = @(alpha)interp1(ang,CLw,alpha);
CD.w = @(alpha)interp1(ang,CDw,alpha);
Cm.w = @(alpha)interp1(ang,Cmw,alpha);

CL.h = @(alpha)interp1(ang,CLeh,alpha);
CD.h = @(alpha)interp1(ang,CDeh,alpha);
Cm.h = @(alpha)interp1(ang,Cmeh,alpha);

%% Limpando a casa

fclose('all');
delete(file_asa)
delete(file_eh)
parfor alpha = linspace(alpha1,alpha2,alpha2-alpha1+1)
    delete(strcat('temp_commands_',int2str(alpha),'.dat'))
    delete(strcat('polarAVL_',int2str(alpha),'.dat'))
    delete(strcat('temp_commands_EH_',int2str(alpha),'.dat'))
    delete(strcat('polarAVL_EH_',int2str(alpha),'.dat'))
end
fprintf('tempo coefs : %f \n',toc)
end