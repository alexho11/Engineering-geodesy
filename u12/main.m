% Inggeo Uebung 12
% Hsin-Feng Ho 3378849
clc
clear all
close all
%% Import Data
load data.mat 

%% Aufgabe a
zeta_1 = data(1:20,4) - data(1:20,3); % Höhenanomalie 1 - 20
sigma_HN = 0.001; % Standardabweichung Normalhöhen
sigma_e = 0.005;  % Standardabweichung Ellipslid Höhe
sigma_zeta = sqrt(sigma_HN^2 + sigma_e^2);  % Fehlerfortpflanzung
figure
scatter3(data(1:20,2),data(1:20,1),zeta_1)

s_x = mean(data(1:20,2));
s_y = mean(data(1:20,1));
% hold on 
% scatter3(s_x,s_y,mean(zeta_1))
%% Aufgabe b
A = [ones(20,1), data(1:20,1) - s_y, data(1:20,2) - s_x, (data(1:20,1) - s_y).* (data(1:20,2) - s_x), (data(1:20,1) - s_y).^2, (data(1:20,2) - s_x) .^2]; 
x_d = inv(A'*A)*A'*zeta_1 ;


%% Aufgabe c
r = 20 - length(x_d);

Sigma_a = sigma_zeta^2 \ (A' * A); % test


sigma_a = sqrt(diag(Sigma_a));
T = abs(x_d - 0) ./ sigma_a;
Q = tinv(1 - 0.025 / length(x_d), r);   % Quantil
idx = find(T < Q);


a_list = cell(6,1);
T_list = cell(6,1);



%%
i = 1;
id = zeros(6,1) * NaN;

check = zeros(6,1) * NaN;
check_list = 1:6;

while ~isempty(idx)
    a_list{i} = x_d;
    T_list{i} = T;
    id(i) = find(T == min(T));
    check(i) = check_list(id(i));
    check_list(id(i)) = [];
    A(:,id(i)) = [];
    x_d = inv(A'*A)*A' * zeta_1;

        
    r = 20 - length(x_d);  
    Sigma_a = sigma_zeta^2 * inv(A' * A); 
   
    sigma_a = sqrt(diag(Sigma_a));
    T = abs(x_d - 0) ./ sigma_a;
    Q = tinv(1 - 0.025 / length(x_d), r);
    idx = find(T < Q);
    i = i + 1;
end
T_list{i} = T;
a_list{i} = x_d;
xq = min(data(1:20,2)):50:max(data(1:20,2));
yq = min(data(1:20,1)):50:max(data(1:20,1));
[xq,yq] = meshgrid(xq,yq);
vq = griddata(data(1:20,2),data(1:20,1),zeta_1,xq,yq);
figure,hold on 
mesh(xq,yq,vq)
xlabel("x")
ylabel("y")
zlabel("Höhenanomalie")


%% Aufgabe d
zeta_2 = x_d(1)+x_d(2).*(data(21:30,1)-s_y)+x_d(3).*(data(21:30,1)-s_y).*(data(21:30,2)-s_x);
NH_under =  data(21:30,4) - zeta_2; % Normalhöhen 21 - 30
data(21:30,3) = NH_under;


%% Aufgabe e
F = [eye(10),-ones(10,1),-(data(21:30,1)-s_y),-(data(21:30,1)-s_y).*(data(21:30,2)-s_x)];
[~,l] = size(F);
Sigma_big = zeros(l,l);
Sigma_big(1:10,1:10) = 0.005^2 * eye(10);
Sigma_big(11:l,11:l) = Sigma_a;

Sigma_nh = F * Sigma_big * F';

sigma_nh = sqrt(diag(Sigma_nh));