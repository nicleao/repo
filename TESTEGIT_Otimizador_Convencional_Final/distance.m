function [distance] = distance(funcao,inicio,fim,n)

switch n
% n = 1: distance btw 2 points (y axis)
% n = 2: curve length
    
    case 1
        if(funcao(inicio,2)*funcao(fim,2) > 0)
            distance = abs((norm(funcao(inicio,2)) - norm(funcao(fim,2))));
        else
            distance = abs((norm(funcao(inicio,2)) + norm(funcao(fim,2))));
        end
        
    case 2
        
        distance = 0;
        for i=inicio:(fim-1)
    
            if(funcao(i,2)*funcao(i,2) > 0)
                distance = distance + abs((norm(funcao(i,:)) - norm(funcao(i+1,:))));
            else
                distance = distance + abs((norm(funcao(i,:)) + norm(funcao(i+1,:))));
            end
    
        end
        
end

[distance];
