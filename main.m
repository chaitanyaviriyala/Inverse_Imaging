% Make sure you close all before executing 

% Author : Chaitanya Krishna V
% TP IMINV
% Date : 07 - 12 - 2020
% Time : 13h45 - 18h00
% Language : English

clc
close all

%% Load the data and visualise 

load('TP1_Exemple.mat');  
% Once we load, we see that there are two variables gth and Z
% gth = 453 * 1 the *gain* of the detector
% and Z = 2530 * 453  the *true pixel intensity*
W = transpose(gth) .* Z; 
fth = log(gth);
fth_centre = fth-mean(fth);
colormap gray; 
imagesc(W);
figure()
colormap gray
imagesc(Z);


%% Function Phi 

abs_x = -2:0.01:2;
phi1 = abs(abs_x);
phi2 = abs_x.^2;
plot(abs_x,phi1,abs_x,phi2);
legend('Absolute Value', 'Squared')

%% Method 0 

g_meth0 = ones(size(gth));
f_meth0 = log(g_meth0);
err_0 = sqrt((sum((f_meth0-fth_centre).^2)));



%% Method 1 

V=log(W);
Y=log(Z);
f_meth1=mean(V,1);
f_meth1_centre = f_meth1-mean(f_meth1);
%g_meth1=exp(f_meth1)';
err_1 = sqrt((sum((f_meth1_centre'-fth_centre).^2)));


%% Differences between method 0 and 1

plot(fth_centre)
hold on
plot(f_meth0)
hold on
plot(f_meth1_centre)
legend('Theoritical Gain centered','Method 0 Gain centred','Method 1 Gain centred')


%% Meth 2

f_meth2=zeros(size(gth)); 
for i =2:453
    deltaV = V(:,i)-V(:,i-1);
    delta_f = median(deltaV); f_meth2(i) = f_meth2(i-1)+delta_f;
end
g_meth2=exp(f_meth2);
f_meth2_centre = f_meth2-mean(f_meth2);
err_2 = sqrt((sum((f_meth2_centre-fth_centre).^2)));
close all
plot(fth_centre)
hold on
plot(f_meth0)
hold on
plot(f_meth1_centre)
hold on
plot(f_meth2_centre)
legend('Theoritical Gain','Method 0 Gain','Method 1 Gain','Method 2 Gain')
%% Meth3

D = diag(ones(1,453))+diag(-ones(1,452),-1);
D(1,1)=0;
%mu = mean(Y,'all');
s=mean(V,1)';
errs3=zeros(20,1);
mus = logspace(-4,6,20);
i=0;
mu_min=mus(1);
min=100000;
for mu=mus
    i=i+1;
    f_meth3=s-inv(D'*D/mu+eye(453))*s;
    %g_meth3=exp(f_meth3);
    f_meth3_centre = f_meth3-mean(f_meth3);
    err_3 = sqrt((sum((f_meth3_centre-fth_centre).^2))); errs3(i)=err_3;
    if err_3<min
        min=err_3;
        mu_min=mu;
    end
end
semilogx(mus,errs3);
xlabel("mu'");
ylabel('Error');
title('Method 3 (phi(x)=x^2)');
f_meth3=s-inv(D'*D/mu_min+eye(453))*s; 
%g_meth3=exp(f_meth3);
f_meth3_centre = f_meth3-mean(f_meth3);
err_3 = sqrt((sum((f_meth3_centre-fth_centre).^2)));


%% Meth4
close all
errs=zeros(10,1);
mus = logspace(-4,6,10);
i=0;
min=10^5;
mu_min=mus(1);

for mu=mus
    i=i+1;
    D=diff(eye(453));
    Vt=V';
    lbd = mu;
    %mean(Y,'all');
    f_meth4 = MAPL1(Vt,D,lbd);
    f_meth4_centre= f_meth4-mean(f_meth4);
    err_4 = sqrt((sum((f_meth4_centre-fth_centre).^2))); errs(i,1) = err_4;
    if err_4<min
        min=err_4;
        mu_min=mu;
    end
end
 
semilogx(mus,errs); 
xlabel("mu'");
ylabel('Error'); 
title('Method 4 (phi(x)=|x|)');
ylim([-1 6])
lbd = mu_min;
%mean(Y,'all');
f_meth4 = MAPL1(Vt,D,lbd);
f_meth4_centre= f_meth4-mean(f_meth4);
err_4 = sqrt((sum((f_meth4_centre-fth_centre).^2)));


%% Comparaison methods 3,4,5
close all
%close all
%plot(fth_centre)
%hold on
%plot(f_meth2_centre)
%hold on
plot(abs(f_meth3_centre-fth_centre)/max(abs(fth_centre)))
hold on
plot(abs(f_meth4_centre-fth_centre)/max(abs(fth_centre)))
legend('Method 3 Gain centered','Method 4 Gain centered')
%'Gain méthode 2 centré','Gain méthode 3 centré','Gain méthode 4 centré')% 

%% Last Graph TP
close all
plot(fth_centre)
hold on
plot(f_meth0)
hold on
plot(f_meth3_centre)
hold on
plot(f_meth4_centre)
legend('Theorical Gain ','Method 0 Gain','Method 3 Gain','Method 4 Gain')
