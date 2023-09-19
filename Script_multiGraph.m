clear
clc

[fname,chemin]=uigetfile('*.*','MultiSelect','on');%loads file

if isequal([fname chemin],[0,0])
    return
else
 datname=[chemin fname];
    d=importdata(datname,',',1);
    data=d.data;
    

figure
% hold on
% k=1; factor=k*data(1,1);
% plot(data(:,1))
% for i=2:size(data,2)
%    
%     
%    
%    k= k+1; factor=k;%*data(1,i)*max(data(:,i));
%    plot(factor*data(:,i))
% end
% hold off
sr=2;dt=1/sr;
xt=1:size(data,2);
time=(0:size(data,1)-1).*dt;
hold on
data1=data';
for i = 1:numel(xt)
      plot3(time,i*ones(size(time)),data1(i,:))
end
hold off
grid on

xlabel('time (sec)')
ylabel('Numero Cellule')  
zlabel('Fluo')

% t = 0:0.01:5*pi;
% y = sin(t);
% z = 0:pi/12:2*pi;
% clf()
% axes()
% hold on
% for i = 1:numel(z)
%     plot3(t,z(i)*ones(size(t)),y);
% end
% grid on
% xlabel('t')
% ylabel('o.5 pi')  % Replace with ylabel(['0.5',char(960)]) for pi symbol
% zlabel('sin(t)')
% view(12.6, 27.6)
end