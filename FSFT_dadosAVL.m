% 2021 - TA
function [Yc,Cc,Areac,ai,cl_norm,cl,CDi_polar,CL_polar] = FSFT_dadosAVL(alpha,beta,fileASA)

ID1 = ('temp_commands.dat');
ID2 = ('Fsurface_polar.dat');
ID3 = ('Ftotal_polar.dat');

fileID = fopen(ID1,'w');

fprintf(fileID,  'load %s \n',fileASA);
fprintf(fileID,  'oper \n');
fprintf(fileID,  'a \na \n');
fprintf(fileID,  '%.4f\n',alpha); 
fprintf(fileID,  'b \nb \n');
fprintf(fileID,  '%.4f\n',beta); 
fprintf(fileID,  ' D1\nD1\n 0\n');
fprintf(fileID,  'D2 \nD2\n 0\n');
fprintf(fileID,  'X \n');
fprintf(fileID,  'FS \n%s \n',ID2);
fprintf(fileID,  'FT \n%s \n',ID3);
fclose('all');

exec = strcat('avl.exe < ',ID1);  % Executa comandos no XFOIL
[~,~] = dos(exec);

delete(ID1)

%% Forças por Strip, Arrasto Induzido Polar e Sustentação para o alpha dado
fileID2 = fopen(ID2,'r');

clear data
data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',15);
CDi_asa = data{2};

%% Forças por Strip
formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %[^\n\r]';
data = textscan(fileID2,formatSpec,'Headerlines',4);

Yc = data{:,2};
Cc = data{:,3};
Areac = data{:,4};
ai = rad2deg(data{:,6});
cl_norm = data{:,7};
cl = data{:,8};
%% FS para EH e EV

clear data
data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',45);
CDi_eh = data{2};

clear data
data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',69);
CDi_ev = data{2};

%% Sustentação para o alpha dado
fileID3 = fopen(ID3,'r');
CL_polar = textscan(fileID3,'  CLtot = %f','Headerlines',23);
CL_polar = CL_polar{1};
%%

CDi_polar = textscan(fileID3,'  CDtot = %f','Headerlines',0);
CDi_polar = CDi_polar{1};
CDi_polar = abs(CDi_polar);

if CDi_eh < 0
    CDi_eh = abs(CDi_eh);
    fact = abs(CDi_polar-CDi_asa)/CDi_eh;
    CDi_polar = CDi_polar + 2*fact*CDi_eh;
end

if CDi_ev < 0
    CDi_ev = abs(CDi_ev);
    fact = abs(CDi_polar-CDi_asa)/CDi_ev;
    CDi_polar = CDi_polar + 2*fact*CDi_ev;    
end

fclose all;

delete(ID1)
delete(ID2)
delete(ID3)
end