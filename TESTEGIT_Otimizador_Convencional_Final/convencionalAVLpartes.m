%% Observações gerais
% Função que cria dois arquivos .avl da configuração covencional, um da asa e outro da EH retangular. 
% Arquivo de comandos para rodar dentro do AVL e criar o avião, apenas.
% IMPORTANTE: Os arquivos estão referenciados em seus respectivos Centros Aerodinâmicos, na altura de seu bordo de ataque.

%% Inputs
% S: área da asa [m^2]
% mac: corda média aerodinâmica da asa [m]
% b: envergadura da asa [m]
% posb: vetor com posição das cordas da asa [m]
% offset1: offset da segunda seção da asa [m]
% offset2: offset da ponta da asa [m]
% died: ângulo de diedro [°]
% died_emp: ângulo de diedro da empenagem [°]
% ail:  distância, na semi-envergadura, onde se inicia o aileron
% c_ail: corda da asa quando inicia o aileron
% chord: vetor com cordas da asa [m]
% Xhinge_eh: [%] da corda da hinge da EH
% Xhinge_a: [%] da corda da hinge do aileron
% SgnDup_aileron: sinal de deflexão para DUPLICATE do AVL (-1 para convencional, -2 para diferencial)
% perfilr: perfil da raiz
% perfilm: perfil do meio
% perfilp: perfil da ponta
% perfilh: perfil da EH
% beh: envergadura da EH [m]
% offseteh: offset na EH [m]
% creh: corda da raiz da EH [m]
% cteh: corda da ponta da EH [m]
% iw: incidência da asa [°]
% ih: incidência da EH [°]
% dist_CA: offset do CA da MAC [m] (distância longitudinal entre BA da raiz e CA da mac)
% Seh: área da EH [m^2]
% maceh: corda média aerodinâmica da EH [m]

%% Função para escrita das partes da aeronave no AVL

function[file_asa,file_eh] = convencionalAVLpartes(S,mac,b,posb,offset1,offset2,died,chord,perfilr,perfilm,perfilp,perfilh,beh,offseteh,creh,cteh,iw,ih,dist_CA,Seh,maceh)

file_asa = 'ConvencionalAVLasa.avl';
file_eh = 'ConvencionalAVLeh.avl';

%% ARQUIVO DA ASA

fileID = fopen(file_asa,'w');

fprintf(fileID,  'Asa Convencional TA 2021\n');
fprintf(fileID,  '#Mach \n');
fprintf(fileID,  '0.0 \n');

fprintf(fileID,  '#IYsym IZsym Zsym \n');
fprintf(fileID, '%-8d%-8d%-8d\n\n',0,0,0);

fprintf(fileID,  '#Sref Cref Bref \n');
fprintf(fileID, '%.4f %.4f %.4f \n\n',S,mac,b);

fprintf(fileID, '#Xref Yref Zref \n');
fprintf(fileID, '%-8d%-8d%-8d\n\n',0,0,0);

fprintf(fileID, '#CDp \n');
fprintf(fileID,'%-8.4f\n\n',0);

%--------------------------------------------------------------------------------------

fprintf(fileID,  '#------------------------------------------------------------- \n');
fprintf(fileID,  'SURFACE \n');
fprintf(fileID,  'Wing \n\n');

fprintf(fileID,  '#Nchordwise  Ccpace  Nspanwise  Sspace \n');
fprintf(fileID,  '20           1       20         0 \n\n');

fprintf(fileID,'COMPONENT\n');
fprintf(fileID,  '1 \n\n');

fprintf(fileID,  'YDUPLICATE \n');
fprintf(fileID,  '0.0 \n\n');

fprintf(fileID,'TRANSLATE\n');
fprintf(fileID,'%-8.4f%-8.4f%-8.4f\n\n',-dist_CA,0,0);

fprintf(fileID,  'ANGLE \n');
fprintf(fileID, '%.4f \n\n',iw);

% Seção da raiz
Xle = 0;
Yle = 0;
Zle = 0;

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,chord(1));

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilr);

% Seção do meio 
Xle = offset1;
Yle = posb(2);
Zle = posb(2)*tand(died);

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,chord(2));

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilm);

% Seção da ponta
Xle = offset2;
Yle = b/2;
Zle = (b/2)*tand(died);

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,chord(3));

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilp);

%% ARQUIVO EH

fileID = fopen(file_eh,'w');

fprintf(fileID,  'EH Convencional TA 2021\n');
fprintf(fileID,  '#Mach \n');
fprintf(fileID,  '0.0 \n');

fprintf(fileID,  '#IYsym IZsym Zsym \n');
fprintf(fileID, '%-8d%-8d%-8d\n\n',0,0,0);

fprintf(fileID,  '#Sref Cref Bref \n');
fprintf(fileID, '%.4f %.4f %.4f \n\n',Seh,maceh,beh);

fprintf(fileID, '#Xref Yref Zref \n');
fprintf(fileID, '%-8d%-8d%-8d\n\n',0,0,0);

fprintf(fileID, '#CDp \n');
fprintf(fileID,'%-8.4f\n\n',0);

%--------------------------------------------------------------------------------------

fprintf(fileID,  '#------------------------------------------------------------- \n');
fprintf(fileID,  'SURFACE \n');
fprintf(fileID,  'EH \n\n');

fprintf(fileID,  '#Nchordwise  Ccpace  Nspanwise  Sspace \n');
fprintf(fileID,  '10             2       20         -2 \n\n');

fprintf(fileID,'COMPONENT\n');
fprintf(fileID,  '2 \n\n');

fprintf(fileID,  'YDUPLICATE \n');
fprintf(fileID,  '0.0 \n\n');

fprintf(fileID,'TRANSLATE\n');
fprintf(fileID,'%-8.4f%-8.4f%-8.4f\n\n',-0.25*creh,0,0);

fprintf(fileID,  'ANGLE \n');
fprintf(fileID, '%.4f \n\n',ih);

%--------------------------------------------------------------------------------------
% Seção da raiz

Xle = 0;
Yle = 0;
Zle = 0;

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,creh);

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilh);

% Seção da ponta

Xle = offseteh;
Yle = beh/2;
Zle = (beh/2)*tand(0);

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.3f    0  \n\n',Xle,Yle,Zle,cteh);

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilh);

%%
fclose('all');

end