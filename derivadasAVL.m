function [CLa,Clb,Cma,Cnb,Clp,Clr,Cmq,Cnp,Cnr,CLde,Clda,Cmde,Cnda,CYda,CYdl,Cldl,Cndl,H] = derivadasAVL(alpha,beta,mac,file_aviao)

ID1 = strcat('temp_commands.dat');
ID2 = strcat('derivadasAVL.dat');

if exist(ID1,'file') == 2      % Evita arquivos já existentes
    delete(ID1)
    delete(ID2)
end

fileID = fopen(ID1,'w');

fprintf(fileID,  'load %s \n',file_aviao);
fprintf(fileID,  'oper \n');
fprintf(fileID,  'a \na \n');
fprintf(fileID,  '%.4f\n',alpha); 
fprintf(fileID,  'b \nb \n');
fprintf(fileID,  '%.4f\n',beta); 
fprintf(fileID,  'X \nST \n');
fprintf(fileID,  '%s \n',ID2);
fprintf(fileID,  'o\n\n\n\n\n\n\n\n\n\n');
fclose('all');

exec = strcat('avl.exe < ',ID1); 
[~,~] = dos(exec);

fileID2 = fopen(ID2,'r');
data=textscan(fileID2,[
    ' z'' force CL |    CLa = %f    CLb = %f\r\n' ...
    ' y  force CY |    CYa = %f    CYb = %f\r\n' ...
    ' x'' mom.  Cl''|    Cla = %f    Clb = %f\r\n' ...
    ' y  mom.  Cm |    Cma = %f    Cmb = %f\r\n' ...
    ' z'' mom.  Cn''|    Cna = %f    Cnb = %f\r\n'],'Headerlines',39);

CLa = data{1};
Clb = data{6};
Cma = data{7};
Cnb = data{10};

data=textscan(fileID2,[
    ' z'' force CL |    CLp = %f    CLq = %f    CLr = %f\r\n'...
    ' y  force CY |    CYp = %f    CYq = %f    CYr = %f\r\n'...
    ' x'' mom.  Cl''|    Clp = %f    Clq = %f    Clr = %f\r\n'...
    ' y  mom.  Cm |    Cmp = %f    Cmq = %f    Cmr = %f\r\n'...
    ' z'' mom.  Cn''|    Cnp = %f    Cnq = %f    Cnr = %f\r\n'],'Headerlines',2);

Clp = data{7};
Clr = data{9};
Cmq = data{11};
Cnp = data{13};
Cnr = data{15};

data = textscan(fileID2,[
    ' z'' force CL |   CLd1 = %f    CLd2 = %f   CLd3 = %f\r\n'...
    ' y  force CY |   CYd1 = %f     CYd2 = %f  CYd3 = %f\r\n'...
    ' x'' mom.  Cl''|   Cld1 = %f   Cld2 = %f   Cld3 = %f\r\n'...
    ' y  mom.  Cm |   Cmd1 = %f     Cmd2 = %f  Cmd3 = %f\r\n'...
    ' z'' mom.  Cn''|   Cnd1 = %f   Cnd2 = %f   Cnd3 = %f\r\n'],'Headerlines',2);

CLde = data{2};
CYda = data{4};
CYdl = data{6};
Clda = data{7};
Cldl = data{9};
Cmde = data{11};
Cnda = data{13};
Cndl = data{15};

H = textscan(fileID2,' Neutral point  Xnp = %f','Headerlines',4);
H = H{1}/mac;
fclose('all');

delete(ID1)
delete(ID2)
end