function [A] = area(pontos,n)

A = 0;

for i=2:(n-1)
    
    A = A + ( pontos(i,1)*pontos(i+1,2) - pontos(i,2)*pontos(i+1,1) );

end

A = A + ( pontos(n,1)*pontos(2,2) - pontos(n,2)*pontos(2,1) );

[A] = A/2;

