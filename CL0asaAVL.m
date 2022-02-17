function [CLmaxrun,CDrun,cl,Yc] = CL0asaAVL(FsCL0)


ID1 = FsCL0;

fileID2 = fopen(ID1,'r');
data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',15);
CLmaxrun = data{1};
CDrun = data{2};

formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %[^\n\r]';
data = textscan(fileID2,formatSpec,'Headerlines',4);
cl = data{:,7};
Yc = data{:,2};

% clear data
% data = textscan(fileID2,'CLsurf  =   %f     CDsurf  =   %f','Headerlines',38);
% CLmaxrunsup = data{1};
% CDrunsup = data{2};
% 
% clear data
% data = textscan(fileID2,formatSpec,'Headerlines',4);
% clsup = data{:,7};
% cl = (clinf+clsup)/2;
% 
% CLmaxrun = (CLmaxruninf+CLmaxrunsup)/2;
% CDrun = (CDruninf+CDrunsup)/2;

fclose all;
% delete(FsCL0)

end