function [Interpolado]=InterpolaXfoil(perfil0,perfil1,Peso)

if exist('PerfilInt.dat') == 2
    delete('PerfilInt.dat')
end

fID=fopen('ComandosXfoil.dat','w');

fprintf(fID, 'inte\n');
fprintf(fID, 'f\n');
fprintf(fID, '%s\n',perfil0);
fprintf(fID, 'f\n');
fprintf(fID, '%s\n',perfil1);
fprintf(fID, '%f\n',Peso);
fprintf(fID, 'PerfilInt\n');
fprintf(fID, 'pcop\n');
fprintf(fID, 'save\n');
fprintf(fID, 'PerfilInt.dat\n');
fclose(fID); 

[~,~] = dos('xfoil.exe < ComandosXfoil.dat');

delete('ComandosXfoil.dat')

Interpolado='PerfilInt.dat';
% fID2=fopen('PerfilTemporarioInterpolado.dat','r');
% 
% Pontos = textscan(fID2,'     %f    %f','Headerlines',1);
% Pontos=cell2mat(Pontos);
% 
% fclose all;
% delete('PerfilTemporarioInterpolado.dat')
end