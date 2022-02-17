%% Observações Gerais
% Essa função retorna o valor da corda em cada posição da envergadura para
% diversas seções da asa. Por motivos estruturais e construtivos, tentar
% manter ao máximo em 5 variações de corda.
% Alterem com cuidado e anotem a data e alterações muito grandes abaixo em
% "Versões"

% Entradas:
%| chord: vetor com os valores das cordas nas seções desejadas, deve ser necessariamente em ordem decrescente (da raiz para a ponta) [m]
%| posb: posição de cada variação de corda na envergadura, deve ser necessariamente em ordem crescente (de 0 até b/2) [m]
%| X: seção reto-trapezoidal (em [m] para uma semi-envergadura)
%| b: envegadura [m]
%| offset: offset de cada variação de corda na envergadura, deve ser necessariamente da raíz para a ponta

% Saidas:
%| y: posição equivalente na envergadura
%| x: offset do BA em cada posição da envergadura
%| c: valor da corda em cada posição da envergadura

% Versões:
%(04/02/17) Versão 1.0 - Programa base inicial
%(06/02/17) Versão 1.1 - Variáveis y,c antes das condições para evitar erro de output de função (Gustavo)
%(01/03/18) Versão 1.2 - Inserido input de offsets na envergadura para obter a distribuição de offset na envergadura --> obter output 'offset_mac' (Bal)

%% Função corda

function [y,c,x] = Corda(chord,posb,X,b,offset)

    y = (0:0.01:1)*b/2;        
    c = NaN(1,length(y));
    x = NaN(1,length(y));
    
    if posb(1) ~= 0
        display('Erro: Primeira posição de variação de corda deve ser a da raiz (valor 0)')
        return
    end

%     if posb(size(posb,2)) ~= b/2
%         display('Erro: Última posição de variação de corda deve ser a da ponta (valor b/2)')
%         return
%     end

%     if posb(2) <= X
%         display('Erro: Primeira variação de corda deve ocorrer depois da seção reto-trapezoidal')
%         return
%     end

    if issorted(posb) ~= 1 
        display('Erro: Posição das cordas na envergadura deve ser crescente')
        return
    end
    
    if issorted(fliplr(chord)) ~= 1
        display('Erro: Valores das cordas devem ser decrescentes')
        return
    end

    if size(posb,2) > 5
        display('Atenção: Não aconselhável realizar tantas variações de corda')
    end

    if length(posb) ~= length(chord)
        display('Erro: Vetor de posição e de cordas deve ser do mesmo tamanho')
        return
    end

%%
    j = 2;
    
    for i=1:length(y)        

        if y(i) > posb(j)
            j = j+1;
        end

        if y(i) <= X 
            c(i) = chord(1);
            x(i) = 0;
            
        elseif j == 2 
            c(i) = polyval(polyfit([X posb(j)],[chord(1) chord(j)],1),y(i));    
            x(i) = polyval(polyfit([X posb(j)],[0 offset(j)],1),y(i));
        
        else
            c(i) = polyval(polyfit([posb(j-1) posb(j)],[chord(j-1) chord(j)],1),y(i));
            x(i) = polyval(polyfit([posb(j-1) posb(j)],[offset(j-1) offset(j)],1),y(i));
        end
        
    end