%% Observa��es gerais
% Fun��o que cria um arquivo .avl com a configura��o convencional.
% Arquivo de comandos para rodar dentro do AVL e somente criar o avi�o.
% IMPORTANTE: Referenciado no CG

% Entradas:
% S: �rea de refer�ncia [m^2]
% mac: Corda de refer�ncia [m]
% b: Envergadura de refer�ncia [m]
% posb: Vetor com posi��o das cordas da asa [m]
% offset: Vetor offset da asa [m]
% died: Diedro [�]
% ail: Dist�ncia onde se inicia o aileron na semi-envergadura [m]
% c_ail: Corda da asa quando inicia o aileron
% Pcg: porcentagem do CG em rela��o a ra�z [%]
% chord: vetor com cordas da asa [m]
% Xhinge_eh: [%] da corda da hinge da EH
% Xhinge_a: [%] da corda da hinge do aileron
% SgnDup_aileron: sinal de deflex�o para superf�cies duplicadas [-1,0,1]
% perfilr: Perfil da raiz da asa
% perfilm: Perfil o meio da asa
% perfilp: Perfil da ponta da asa
% perfilh: Perfil da EH  
% perfilv: Perfil da EV
% beh: Envergadura na EH
% offseteh: Offset da ponta da EH
% creh: Corda da raiz da EH
% cteh: Corda da ponta da EH
% bev: Envergadura na EV
% crev: Corda da raiz da EV
% ctev: Corda da ponta da EV
% offsetev: Offset na EV
% x_emp: Dist�ncia entre bordo de ataque da asa e bordo de ataque da empenagem (raiz)
% h_emp: Altura entra asa e EH (Bordo de ataque da asa at� bordo de ataque da EH)
% iw: �ngulo de incid�ncia da asa
% ih: �ngulo de incid�ncia da EH
% nEV: Quantidade de EVs
% H_CG: Altura CG a partir da altura da asa
% ef_solo:
% h_asa: Altura da asa em rela��o ao solo (intradorso da raiz ate o solo)
% Zcg: altura (estimada) do CG em rela��o ao intradorso da corda da raiz [m]
% i_corrida_w: �ngulo de corrida da asa [�]
% X: se��o reto-trapezoidal [m]
% Xhinge_elevon: posi��o da hinge [% da corda]

%% Fun��o para escrita da aeronave no AVL
function[file_aviao] = convencionalAVL(S, mac, b, posb, offset, died, ail, c_ail, Pcg, chord, Xhinge_eh, Xhinge_a, SgnDup_aileron, perfilr, perfilm, perfilp, perfilh, perfilv, beh, offseteh, creh, cteh, bev, crev, ctev, offsetev, x_emp, h_emp, iw, ih, nEV, H_CG, ef_solo, h_asa)

% Valores de entrada e condi��es
if ail <= posb(2)
      porcentagem = ail/posb(2);
      [Interpolado] = InterpolaXfoil(perfilr,perfilm,porcentagem);
else
      porcentagem = (ail-posb(2))/(posb(3)-posb(2));
      [Interpolado] = InterpolaXfoil(perfilm,perfilp,porcentagem); 
end

if ef_solo == 0
    file_aviao = 'convencionalAVL.avl';
else
    file_aviao = 'convencionalAVL_ground.avl';
end

% ARQUIVO DO AVI�O
fileID = fopen(file_aviao,'w');

fprintf(fileID,  'Avi�o Convencional AVL 2021\n');
fprintf(fileID,  '#Mach \n');
fprintf(fileID,  '0.0 \n');

fprintf(fileID,  '#IYsym IZsym Zsym \n');
if ef_solo == 0
    fprintf(fileID, '%-8d%-8d%-8d\n\n',0,0,0);
else
    fprintf(fileID, '%-8d%-8d%-8d\n\n',0,1,-(H_CG+h_asa));
end

fprintf(fileID, '#Sref Cref Bref \n');
fprintf(fileID, '%.4f %.4f %.4f \n\n',S,mac,b);

fprintf(fileID, '#Xref Yref Zref \n');
fprintf(fileID, '%-8d%-8d%-8d\n\n',0,0,0);

fprintf(fileID, '#CDp \n');
fprintf(fileID,'%-8.4f\n\n',0);

%% Asa

fprintf(fileID,  '#------------------------------------------------------------- \n');
fprintf(fileID,  'SURFACE \n');
fprintf(fileID,  'Wing \n\n');

fprintf(fileID,  '#Nchordwise  Ccpace  Nspanwise  Sspace \n');
fprintf(fileID,  '20             1       20         0 \n\n');

fprintf(fileID,'COMPONENT\n');
fprintf(fileID,  '1 \n\n');

fprintf(fileID,  'YDUPLICATE \n');
fprintf(fileID,  '0.0 \n\n');

fprintf(fileID,'TRANSLATE\n');
fprintf(fileID,'%-8.4f%-8.4f%-8.4f\n\n',-Pcg*chord(1),0,-H_CG);

fprintf(fileID,  'ANGLE \n');
fprintf(fileID, '%.4f \n\n',iw);

% Se��o da raiz
Xle = offset(1);
Yle = 0;
Zle = 0;

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,chord(1));

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilr);

% Se��o do meio e come�o do aileron
if ail < posb(2)
        % Aileron
        Xle = offset(2)*ail/posb(2);
        Yle = ail;
        Zle = ail*tand(died);
      

        fprintf(fileID, '#------------------------------------------------------------- \n');
        fprintf(fileID, 'SECTION \n');
        fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
        fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,c_ail);

        fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
        fprintf(fileID, 'CONTROL \n');
        fprintf(fileID, 'aileron  -1   %.4f    0    0   0    %.4f \n\n', (1-Xhinge_a),SgnDup_aileron);
           

        fprintf(fileID, 'AFILE \n');     
        fprintf(fileID, '%s \n\n',Interpolado);
        
        % Meio
        Xle = offset(2);
        Yle = posb(2);
        Zle = posb(2)*tand(died);

        fprintf(fileID, '#------------------------------------------------------------- \n');
        fprintf(fileID, 'SECTION \n');
        fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
        fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,chord(2));

        fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
        fprintf(fileID, 'CONTROL \n');
        fprintf(fileID, 'aileron  -1   %.4f    0    0   0    %.4f \n\n', (1-Xhinge_a),SgnDup_aileron);

        fprintf(fileID, 'AFILE \n');
        fprintf(fileID, '%s \n\n',perfilm);
        
        
else
        % Meio
        Xle = offset(2);
        Yle = posb(2);
        Zle = posb(2)*tand(died);

        fprintf(fileID, '#------------------------------------------------------------- \n');
        fprintf(fileID, 'SECTION \n');
        fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
        fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,chord(2));

        fprintf(fileID, 'AFILE \n');
        fprintf(fileID, '%s \n\n',perfilm);
        
        % Aileron
        Xle = (offset(3)*(ail-posb(2)))/((b/2)-posb(2)) + offset(2);
        Yle = ail;
        Zle = ail*tand(died);

        fprintf(fileID, '#------------------------------------------------------------- \n');
        fprintf(fileID, 'SECTION \n');
        fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
        fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,c_ail);

        fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
        fprintf(fileID, 'CONTROL \n');
        fprintf(fileID, 'aileron  -1   %.4f    0    0   0    %.4f \n\n', (1-Xhinge_a),SgnDup_aileron);

        fprintf(fileID, 'AFILE \n');        
        fprintf(fileID, '%s \n\n',Interpolado);
end

% Se��o da ponta
Xle = offset(2) + offset(3);
Yle = b/2;
Zle = (b/2)*tand(died);

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,chord(3));

fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
fprintf(fileID, 'CONTROL \n');
fprintf(fileID, 'aileron  -1   %.4f    0    0   0    %.4f \n\n', (1-Xhinge_a),SgnDup_aileron);

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilp);

%% EH

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
fprintf(fileID,'%-8.4f%-8.4f%-8.4f\n\n',(-Pcg*chord(1) + x_emp),0,h_emp-H_CG);

fprintf(fileID,  'ANGLE \n');
fprintf(fileID, '%.4f \n\n',ih);

% Se��o da raiz
Xle = 0;
Yle = 0;
Zle = 0;

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,creh);

fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
fprintf(fileID, 'CONTROL \n');
fprintf(fileID, 'elevator  1   %.4f    0    0   0    1 \n\n', (1-Xhinge_eh));

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilh);

% Se��o da ponta
Xle = offseteh;
Yle = beh/2;
Zle = 0;

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.3f    0  \n\n',Xle,Yle,Zle,cteh);

fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
fprintf(fileID, 'CONTROL \n');
fprintf(fileID, 'elevator  1   %.4f    0    0   0    1 \n\n', (1-Xhinge_eh));

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilh);

%% EV

fprintf(fileID,  '#------------------------------------------------------------- \n');
fprintf(fileID,  'SURFACE \n');
fprintf(fileID,  'EV \n\n');

fprintf(fileID,  '#Nchordwise  Ccpace  Nspanwise  Sspace \n');
fprintf(fileID,  '8             2       10         -2 \n\n');

fprintf(fileID,'COMPONENT\n');
fprintf(fileID,  '3 \n\n');

if nEV ~= 1
    fprintf(fileID,  'YDUPLICATE \n');
    fprintf(fileID,  '0.0 \n\n');
end

fprintf(fileID,'TRANSLATE\n');
fprintf(fileID,'%-8.4f%-8.4f%-8.4f\n\n',(-Pcg*chord(1) + x_emp),0,h_emp-H_CG);

fprintf(fileID,  'ANGLE \n');
fprintf(fileID, '%.4f \n\n',0);

% Se��o da raiz
Xle = 0;
if nEV == 1
    Yle = 0;
else 
    Yle = bh/2;
end
Zle = 0;

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,crev);

fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
fprintf(fileID, 'CONTROL \n');
fprintf(fileID, 'rudder  1   %.4f    0    0   0    1 \n\n', 0);

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilv);

% Ponta

Xle = offsetev;
if nEV == 1
    Yle = 0;
else 
    Yle = bh/2;
end
Zle = bev;

fprintf(fileID, '#------------------------------------------------------------- \n');
fprintf(fileID, 'SECTION \n');
fprintf(fileID, '#Xle      Yle       Zle      Chord    Ainc \n');
fprintf(fileID, '%.4f   %.4f    %.4f    %.4f    0  \n\n',Xle,Yle,Zle,ctev);

fprintf(fileID, '#Cname   Cgain  Xhinge  HingeVec     SgnDup \n');
fprintf(fileID, 'CONTROL \n');
fprintf(fileID, 'rudder  1   %.4f    0    0   0    1 \n\n', 0);

fprintf(fileID, 'AFILE \n');
fprintf(fileID, '%s \n\n',perfilv);

fclose('all');
end