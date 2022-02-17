
function limpa_lixo
    clc
    
    parfor i=-50:50
       warning ('off','all');
       delete(strcat('FsCLmax_',int2str(i),'.dat'))
       delete(strcat('comandosAVL_',int2str(i),'.dat'))
       delete(strcat('temp_commands_',int2str(i),'.dat'))
       delete(strcat('polarAVL_',int2str(i),'.dat'))
       delete(strcat('temp_commands_EH_',int2str(i),'.dat'))
       delete(strcat('polarAVL_EH_',int2str(i),'.dat'))
    end
    warning ('on','all');
end