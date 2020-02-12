%%%%%% This code generates the excitatory postsynaptic currents of
%%%%%% facilitating, depressing and pseudo-linear excitatory syanpses
%%%%%% based on Tsodyks-Markram synaptic model
%%%%%% Written by AmirAli Farokhniaee, PhD (aafarokh@gmail.com)
%%%%%% Initiated 2016, editions 2017
%%% This code produces the plots in figures 2 and 3 of the following paper:
%%%  Farokhniaee, A., & McIntyre, C. C. (2019). 
%%%  Theoretical principles of deep brain stimulation induced
%%%  synaptic suppression. Brain Stimulation, 12(6), 1402-1409. 
%%%  https://doi.org/10.1016/j.brs.2019.07.005
%%~ It also computes the area under the synaptic currents and compares with
%%~ theoretical values.
clearvars

dt=.1; ti=dt; tf=10000;
t=ti:dt:tf;

% DBS input
fdbs=1:200;
T=ones(1,length(fdbs));
dbsi=500/dt; dbsf=10000/dt;

% I Kernel time constant
taus=1.75; %3

%input spike train
sp=zeros(length(fdbs),length(t));

% Synapse parameters % Each column represents E1, E2 and E3 respectively
tauf=[670,17,326];
taud=[138,671,329];
U=[.09,.5,.29];
% A=[0.0025,0.0025,0.0025];
A=[1,1,1];
n=1;
A=n*A;

% Compute EPSC
u=zeros(length(fdbs),length(t));
x=ones(length(fdbs),length(t));
I=zeros(length(fdbs),length(t));
% It=zeros(length(fdbs),length(t));
EPSC=zeros(length(A),length(fdbs),length(t));
select_time=dbsf-50000:dbsf;
It=zeros(length(fdbs),length(select_time));
M_I=ones(length(A),length(fdbs));   
mi=zeros(length(A),1);
M_Iall=ones(length(A),length(fdbs));   
area=zeros(length(A),length(fdbs));
area1=zeros(length(A),length(fdbs));
% Sc=zeros(length(A),length(fdbs));

for p=1:3
    for j=1:length(fdbs)
        T(j)=round(1000/fdbs(j)/dt);
        ts=dbsi:T(j):dbsf;
        sp(j,ts)=1/dt;
            for i=1:length(t)-1 
                u(j,(i+1)) = u(j,i)+dt*(-(u(j,i)/tauf(p))+U(p)*(1-u(j,i))*sp(j,i)); 
                x(j,(i+1)) = x(j,i) + dt*((1/taud(p))*(1-x(j,i)) - u(j,i+1)*x(j,i)*sp(j,i));
                I(j,(i+1)) = I(j,i) + dt*((-1/taus)*I(j,i) + A(p)*u(j,i+1)*x(j,i)*sp(j,i));
            end
        EPSC(p,j,:)= I(j,:);     
%         M_Iall(p,j)=max(I(j,:));
        It(j,:)=I(j,select_time);
        M_I(p,j)=max(It(j,:));
%         mi(p)=max(M_Iall(p,j));
%         M_I(p,j)=M_I(p,j)./mi(p);
t1{p,j}=ts(end-1)+1:ts(end);   %last period EPSC curve
It1{p,j}=I(j,t1{p,j});
area1(p,j)=trapz(t1{p,j},It1{p,j})/10;       %area under the EPSC curve for 1 EPSC
area(p,j)=area1(p,j)*j;                      %area under the EPSC curve in 1 second
    end  
end

%gain peak frequency
theta=1000/sqrt(tauf(1)*taud(1)*U(1)); %valid only for facilitating synapse

%Make figure
freq1=20; freq2=130;
figure;
ax1=subplot(4,3,1);
plot(t,squeeze(EPSC(1,freq1,:)),'k','LineWidth',1); ylabel({'EPSC (nA)';'Facilitating'},'FontWeight','bold')
xlim([450 1000]); ylim([0 .6]);

ax2=subplot(4,3,2);
plot(t,squeeze(EPSC(1,freq2,:)),'k','LineWidth',1); %ylabel({'I_{syn} (nA)'; 'EPSC'},'FontWeight','bold')
xlim([450 1000]); ylim([0 .6]);

ax3=subplot(4,3,3);
scatter(fdbs,squeeze(M_I(1,:,1)),'k','.'); ylabel({'Facilitating synapse';'EPSC_{st} amplitude (nA)'},'FontWeight','bold')
hold on
% plot((1./fdbs)+.008,'--','LineWidth',1); zoom xon; %ylim([0 .14])
ylim([0 .6]);

ax4=subplot(4,3,4);
plot(t,squeeze(EPSC(2,freq1,:)),'k','LineWidth',1); ylabel({'EPSC (nA)';'Depressing'},'FontWeight','bold')
xlim([450 1000]); ylim([0 .6]);

ax5=subplot(4,3,5);
plot(t,squeeze(EPSC(2,freq2,:)),'k','LineWidth',1); %ylabel({'I_{syn} (nA)'; 'EPSC'},'FontWeight','bold')
xlim([450 1000]); ylim([0 .6]);

ax6=subplot(4,3,6);
scatter(fdbs,squeeze(M_I(2,:,1)),'k','.'); ylabel({'Depressing synapse';'EPSC_{st} amplitude (nA)'},'FontWeight','bold')
ylim([0 .6]);

ax7=subplot(4,3,7);
plot(t,squeeze(EPSC(3,freq1,:)),'k','LineWidth',1); ylabel({'EPSC (nA)';'Pseudo-linear'},'FontWeight','bold')
xlim([450 1000]); ylim([0 .6]);

ax8=subplot(4,3,8);
plot(t,squeeze(EPSC(3,freq2,:)),'k','LineWidth',1); %ylabel({'I_{syn} (nA)'; 'EPSC'},'FontWeight','bold')
xlim([450 1000]); ylim([0 .6]);

ax9=subplot(4,3,9);
scatter(fdbs,squeeze(M_I(3,:,1)),'k','.'); ylabel({'Pseudo-linear synapse';'EPSC_{st} amplitude (nA)'},'FontWeight','bold')
ylim([0 .6]);

ax10=subplot(4,3,10);
plot(t,sp(freq1,:),'k','LineWidth',1); zoom xon; ylabel(['Input ',num2str(freq1),' Hz'],'FontWeight','bold')
xlim([450 1000]); 
xlabel('Time (ms)','FontWeight','bold'); 

ax11=subplot(4,3,11);
plot(t,sp(freq2,:),'k','LineWidth',1); zoom xon; ylabel(['Input ',num2str(freq2),' Hz'],'FontWeight','bold')
xlim([450 1000]); 
xlabel('Time (ms)','FontWeight','bold'); 

ax12=subplot(4,3,12);
scatter(fdbs,squeeze(M_I(1,:,1)),'.'); zoom xon; hold on;
scatter(fdbs,squeeze(M_I(2,:,1)),'.'); hold on
scatter(fdbs,squeeze(M_I(3,:,1)),'.'); hold on
xlabel('Frequency (Hz)','FontWeight','bold'); 
ylabel({'All synapses';'EPSC_{st} amplitude (nA)'},'FontWeight','bold')
ylim([0 .6]);

figure
scatter(fdbs,squeeze(M_I(1,:,1)),'filled'); zoom xon; hold on;
scatter(fdbs,squeeze(M_I(2,:,1)),'filled'); hold on
scatter(fdbs,squeeze(M_I(3,:,1)),'filled'); hold on
xlabel('Frequency (Hz)','FontWeight','bold'); 
ylabel({'EPSC_{st} amplitude (nA)'},'FontWeight','bold')
ylim([0 .5])
set(gca,'FontSize',12,'FontWeight','bold')

%% Integrals
S1=zeros(length(A),length(fdbs));
S=zeros(length(A),length(fdbs));
for j=1:3
for i=1:length(fdbs)
S1(j,i) = -M_I(j,i)*taus*(exp(-T(i)/taus)-1); %The integral of one EPSC at the steady state
S(j,i) = S1(j,i)*i;
end
end

area1f=area1(1,:);
area1d=area1(2,:);
area1p=area1(3,:);
areaf=area(1,:);
aread=area(2,:);
areap=area(3,:);

f_weight=.45; d_weight=.38; p_weight=.18; 
area_tot=f_weight*areaf+d_weight*aread+p_weight*areap;

figure; title('Area under 1 EPSC'); hold on
for p=1:3
scatter(fdbs,area1(p,:));
hold on
end
legend('F','D','P')
for p=1:3
plot(fdbs,S1(p,:),'Linewidth',1);
hold on
end
xlabel('DBS frequency (Hz)')
ylabel('S_1')
set(gca,'FontSize',12,'FontWeight','bold')

figure; title('Area under EPSCs in 1 second of stimulation'); hold on
for p=1:3
scatter(fdbs,area(p,:),'filled');
hold on
end
scatter(fdbs,area_tot,'filled','k');

legend('F','D','P','Total')
% for p=1:3
% plot(fdbs,S(p,:),'Linewidth',1);
% hold on
% end
% plot(ff,sf,'LineWidth',1)
xlabel('DBS frequency (Hz)')
ylabel('S')
set(gca,'FontSize',12,'FontWeight','bold')


