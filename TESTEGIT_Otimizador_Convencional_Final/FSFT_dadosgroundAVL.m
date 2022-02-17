function[Yc,FsCL0,Areac,ai,cl_norm,cl,CDi_corrida,CL_corrida] = FSFT_dadosgroundAVL(alpha,beta,ro,V,file_aviao_ground)
%% Comandos
ID1 = strcat('comandosAVLground.dat');
ID2 = strcat('Fsurface_ground.dat');
ID3 = strcat('Ftotal_ground.dat');
ID4 = strcat('FsCL0.dat');

if exist(ID1,'file') == 2      % Evita arquivos já existentes
    delete(ID1)
    delete(ID2)
    delete(ID3)
    delete(ID4)
end

fileID = fopen(ID1,'w');
fprintf(fileID,  'load %s \n',file_aviao_ground);
fprintf(fileID,  'oper \n');
fprintf(fileID,  'a \na \n');
fprintf(fileID,  '%.4f\n',0); 
fprintf(fileID,  'D1 \nD1 \n 0\n');
fprintf(fileID,  'D2 \nD2 \n 0\n');
fprintf(fileID,  'X \n');
fprintf(fileID,  'FS \n%s \n',ID2);
% fprintf(fileID,  'O\n'); %tirar para paralelo
fprintf(fileID,  'FT \n%s \n',ID3);
fprintf(fileID,  'O\n\n\n\n\n\n\n\n\n\n\n');

%% Cl0asaAVLcomandos
fprintf(fileID,  'load %s \n',file_aviao_ground);
fprintf(fileID,  'oper \n');
fprintf(fileID,  'M \n');
fprintf(fileID,  'V \n%.4f \n',V);
fprintf(fileID,  'D \n%.4f \n\n',ro);
fprintf(fileID,  'oper \n');
fprintf(fileID,  'a \na \n');
fprintf(fileID,  '%.4f\n',alpha); 
fprintf(fileID,  'b \nb \n');
fprintf(fileID,  '%.4f\n',beta); 
fprintf(fileID,  'D1 \nD1\n 0\n');
fprintf(fileID,  'D2 \nPM\n 0\n');
fprintf(fileID,  'X \n');
fprintf(fileID,  'FS \n%s \n',ID4);
fprintf(fileID,  'O\n\n\n\n\n\n\n\n\n\n\n');
fclose('all');

%% Análises
exec = strcat('avl.exe < ',ID1);
[~,~] = dos(exec);

delete(ID1);

fileID2 = fopen(ID2,'r');
data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',15);
CDi_asa = data{2};

clear data
formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %[^\n\r]';
data = textscan(fileID2,formatSpec,'Headerlines',4);

Yc = data{:,2};
Cc = data{:,3};
Areac = data{:,4};
ai = rad2deg(data{:,6});
cl_norm = data{:,7};
cl = data{:,8};

% clear data
% data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',38);
% CDi_sup = data{2};
% 
% CDi_asa = (CDi_inf+CDi_sup)/2;

clear data
data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',45);
CDi_eh = data{2};

clear data
data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',69);
CDi_ev = data{2};

fclose all;

%% Sustentação para o ângulo de corrida
fileID3 = fopen(ID3,'r');
CL_corrida = textscan(fileID3,'  CLtot = %f','Headerlines',23);
CL_corrida = CL_corrida{1};

%% Arrasto induzido para o ângulo de corrida
CDi_corrida = textscan(fileID3,'  CDtot = %f','Headerlines',0);
CDi_corrida = CDi_corrida{1};

if CDi_corrida < 0
    disp('Warning: CD negativo!')
    CDi_corrida = abs(CDi_corrida);
end

if CDi_eh < 0
    CDi_eh = abs(CDi_eh);
    fact = abs(CDi_corrida-CDi_asa)/CDi_eh;
    CDi_corrida = CDi_corrida + 2*fact*CDi_eh;
end

if CDi_ev < 0
    CDi_ev = abs(CDi_ev);
    fact = abs(CDi_corrida-CDi_asa)/CDi_ev;
    CDi_corrida = CDi_corrida + 2*fact*CDi_ev;
end


FsCL0 = ID4;
nada = 0;
fclose all;
delete(ID2)
delete(ID3)
end