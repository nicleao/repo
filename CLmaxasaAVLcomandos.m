function [CLmaxrun,cl,Yc] = CLmaxasaAVLcomandos(alpha,beta,ro,V,fileASA)

ID1 = strcat('comandosAVL.dat');
ID2 = strcat('FsCLmax.dat');

if exist(ID1,'file') == 2      % Evita arquivos já existentes
    delete(ID1)
    delete(ID2)
end

fileID = fopen(ID1,'w');

fprintf(fileID,  'load %s \n',fileASA);
fprintf(fileID,  'oper \n');
fprintf(fileID,  'M \n');
fprintf(fileID,  'V \n%.4f \n',V);
fprintf(fileID,  'D \n%.4f \n\n',ro);
fprintf(fileID,  'a \na \n');
fprintf(fileID,  '%.4f\n',alpha); 
fprintf(fileID,  'b \nb \n');
fprintf(fileID,  '%.4f\n',beta); 
fprintf(fileID,  'D1 \nD1\n 0\n');
fprintf(fileID,  'D2 \nD2\n 0\n');
fprintf(fileID,  'X \n');
fprintf(fileID,  'FS \n%s \n',ID2);
fprintf(fileID,  'O\n');
fclose('all');

exec = strcat('avl.exe < ',ID1);
[~,~] = dos(exec);

% antigo (para caso de algum problema)
% fprintf(fileID,  'load %s ',fileASA);
% fprintf(fileID,  char([13,10]));
% fprintf(fileID,  'oper \n');
% fprintf(fileID,  char([13,10]));
% fprintf(fileID,  'M \n');
% fprintf(fileID,  char([13,10]));
% fprintf(fileID,  'V') ;
% fprintf(fileID,  char([13,10]));
% fprintf(fileID, '%.4f' ,V);
% fprintf(fileID,  char([13,10]));
% fprintf(fileID,  'D ' );
% fprintf(fileID,  char([13,10]));
% fprintf(fileID,  '%.4f ' ,ro);
% fprintf(fileID,  char([13,10]));
% fprintf(fileID,  ['oper' char([13,10])]);
% fprintf(fileID,  ['a' char([13,10]) 'a ' char([13,10])]);
% fprintf(fileID,  '%.4f',alpha); 
% fprintf(fileID,  char([13,10]));
% fprintf(fileID, [ 'b ' char([13,10]) 'b' char([13,10])]);
% fprintf(fileID,  '%.4f',beta);
% fprintf(fileID,  char([13,10]));
% fprintf(fileID,  ['D1 ' char([13,10]) 'D1' char([13,10]) ' 0' char([13,10])]);
% fprintf(fileID,  ['D2 ' char([13,10]) 'PM' char([13,10]) ' 0' char([13,10])]);
% fprintf(fileID,  ['X ' char([13,10])]);
% fprintf(fileID,  ['FS ' char([13,10]) ]);
% fprintf(fileID, 'FsCLmax%u.dat ',i);
% fprintf(fileID,  char([13,10]));
% fprintf(fileID,  ['O' char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) char([13,10]) ]);

filename = ID2;
fileID2 = fopen(filename,'r');
data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',15);
CLmaxrun = data{1};
% CDrunasa = data{2};

formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %[^\n\r]';
data = textscan(fileID2,formatSpec,'Headerlines',4);
cl = data{:,7};
Yc = data{:,2};

% clear data
% data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',38);
% CLmaxrunsup = data{1};
% % CDrunsup = data{2};
% 
% clear data
% data = textscan(fileID2,formatSpec,'Headerlines',4);
% clsup = data{:,7};

% CLmaxrun = (CLmaxruninf+CLmaxrunsup)/2;

%CDrun = (CDruninf+CDrunsup)/2;

%nada = 0;
fclose all;

delete(ID1)
delete(ID2)

end