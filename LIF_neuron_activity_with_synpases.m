%%%%%% This code computes an LIF neuron activity before and during DBS.
%%%%%% Written by AmirAli Farokhniaee, PhD (aafarokh@gmail.com)
%%%%%% Initiated and editions 2017-2018.
%%% This code particularly produces the plots in figs. 4 and 5
%%% of the following paper:
%%%  Farokhniaee, A., & McIntyre, C. C. (2019). 
%%%  Theoretical principles of deep brain stimulation induced
%%%  synaptic suppression. Brain Stimulation, 12(6), 1402-1409. 
%%%  https://doi.org/10.1016/j.brs.2019.07.005
tic 
clearvars

% transmission + synaptic delay: td 
td=2; %2 ms for trasmission and .5 ms for synaptic delay

dt=.1; ti=dt; tf=61000+td;%tf=1500+td;%tf=61000+td; %in mili seconds 
t=ti:dt:tf;

% DBS input
fdbs=130;
T=ones(1,length(fdbs));
% dbsi=(100)/dt; dbsf=1100/dt; %in mili seconds
dbsi=(dt)/dt; dbsf=61000/dt;

%Poissonian input
fr=10; %for fr Hz baseline poissonian firing from other cells 
[spikes,tsp]=poissonSpikeGen(fr,tf/1000,1,dt/1000);
tp=find(spikes==1); 
% ssp=zeros(1,length(t));
% ssp(tp)=1; %uncomment for stochastic model (adding noise to the system)

%noise term
% wght=0;   %no noise
wght=.5; %default noise
% wght=5;   %high noise
kisi=wght*randn(1,length(t));

% I Kernel time constant
taus=.75;  %For excitatory synapse

% transmission + synaptic delay: td 
td=td/dt; %convert to simulation step scale

%input spike train
sp=zeros(length(fdbs),length(t));

% Synapse parameters % Each column 1,2,3 means F,D,P respectively and each row means
% Excitatory and inhibitory synapse (1: excitatory, 2: inhibitory)
% In this study we just used the first row, excitstory synapses.
tauf=[670,17,326; 376,21,62];
taud=[138,671,329; 45,706,144];
U=[.09,.5,.29; .016,.25,.32];
A=[.0025,.0025,.0025; .0025,.0025,.0025];
% n=10; A=n*A;  % change the strength of A (order of magnitude of totall number of synapses)
ie=ones(1,2);
w=1;

fid=2.5; %synaptic fidelity 
we=fid*200; wi=0; 
% Percentage of excitatory and inhibitory synapses:
ne=we*[45,38,17]; %original: 45,38,17
% ne=zeros(1,3);
% for 1 synapse n1=1 and so forth (approximately giving 2 pA exc. current)
% ne=10;  % for 10 synapses (approximately giving 20 pA exc. current)
% ne=100; % for 100 synapses (approximately giving 200 pA exc. current)
% ne=1000;% for 1000 synapses (approximately giving 2 nA exc. current)
% ni=wi*[13,10,6];     % for 1 synapse (approximately giving 10 pA inhibitory current)
ni=wi*[8,76,16]; %ne=ni;
% ni=zeros(1,3);
% ni=10;  % for 10 synapses (approximately giving 100 pA inh. current)
% ni=100; % for 100 synapses (approximately giving 1 nA inh. current)
% ni=1000;% for 1000 synapses (approximately giving 10 nA inh. current)
A=[ne.*A(1,:);ni.*A(2,:)];

% Compute EPSC
u=zeros(length(fdbs),length(t));
x=ones(length(fdbs),length(t));
I=zeros(length(fdbs),length(t));
Iwo=zeros(length(fdbs),length(t));
% It=zeros(length(fdbs),length(t));
PSC=zeros(length(ie),length(A),length(fdbs),length(t));
% IPSC=zeros(length(A),length(fdbs),length(t));

% Compute EPSP (passive mechanism, membrane potential)
tau_memb=40;
r=10^2;   %M Ohm
v=zeros(length(fdbs),length(t));
PSP=zeros(length(ie),length(A),length(fdbs),length(t));
% IPSP=zeros(length(A),length(fdbs),length(t));

% Neuron parameters: (for ~20 Hz base firing .56 and for ~8-10 Hz choose .26)
Cm= 1; Rm=100; Ie=.26; %(for deterministic model)
% Ie=.16; %subthreshold firing (for noise purpose, stochastic model)
El=-70; Vth=-54; 
Vreset=-80;

% % Neuron parameters: (for 62.5 Hz base firing)
% Cm= 1; Rm=100; Ie=1.52; %(for deterministic model)
% % Ie=.18; %subthreshold firing (for noise purpose, stochastic model)
% El=-70; Vth=-54; 
% Vreset=-80;

% Compute neuron firing pattern with and without synaptic input:
V=zeros(length(fdbs),length(t));
Vn=zeros(length(ie),length(A),length(fdbs),length(t));
V_all=zeros(length(fdbs),length(t));
% Vn_all=zeros(length(ie),length(A),length(fdbs),length(t));
Vin=zeros(1,length(t));

wk=10; %Poissonian weight
poiss=wk*rand(1,length(sp)).*sp(1,:);
for i=1:length(t)-1
Vin(i+1) = Vin(i) + (dt/Cm)*(((El-Vin(i))/Rm) + Ie + poiss(i) + kisi(i));
if Vin(i+1)>= Vth+kisi(i)
    Vin(i)=0+kisi(i);
    Vin(i+1)=Vreset+kisi(i);
end
end

for q=1:length(ie)
    if q==1 
        w=1;
    else
        w=-1;
    end
for p=1:length(A)
    for j=1:length(fdbs)
        T(j)=round((1000/fdbs(j))/dt);
        dbs=dbsi:T(j):dbsf;
        ts=[tp,dbs]; %uncomment for Poissonian+DBS
%         ts=dbs;      %uncomment for DBS only
        sp(j,ts)=1/dt;
            for i=td+1:length(t)-1  
                u(j,(i+1)) = u(j,i) + dt*(-(u(j,i)/tauf(q,p))+U(q,p)*(1-u(j,i))*sp(j,i-td)); 
                x(j,(i+1)) = x(j,i) + dt*((1/taud(q,p))*(1-x(j,i)) - u(j,i+1)*x(j,i)*sp(j,i-td));
                I(j,(i+1)) = I(j,i) + dt*((-1/taus)*I(j,i) + A(q,p)*u(j,i+1)*x(j,i)*sp(j,i-td));
                Iwo(j,(i+1)) = Iwo(j,i) + dt*((-1/taus)*Iwo(j,i) + A(q,p)*sp(j,i-td));
                v(j,(i+1)) = v(j,i) + dt*(((-v(j,i)+r*I(j,i)))/tau_memb);
                %Replace I with Iwo for no depletion of synaptic conduction
                V(j,(i+1)) = V(j,i) + (dt/Cm)*(((El-V(j,i))/Rm) +Ie + w*I(j,i) + poiss(i)+ kisi(i));
                    if  V(j,i+1)>= Vth+kisi(i)
                        V(j,i)=0+kisi(i);
                        V(j,i+1)=Vreset+kisi(i);
                    end
            end
            %replace I with Iwo for no depletion
        PSC(q,p,j,:)= w*I(j,:); %IPSC(p,j,:)= -I(j,:);
        PSP(q,p,j,:)= w*v(j,:); %IPSP(p,j,:)= -v(j,:);
        Vn(q,p,j,:)= V(j,:);

    end
end
end

PSC_exc=sum(PSC(1,:,:,:),2);
PSC_inh=sum(PSC(2,:,:,:),2);
PSC_all=PSC_exc+PSC_inh;

PSP_exc=sum(PSP(1,:,:,:),2);
PSP_inh=sum(PSP(2,:,:,:),2);
PSP_all=PSP_exc+PSP_inh;

for j=1:length(fdbs)
for i=1:length(t)-1  
                V_all(j,(i+1)) = V_all(j,i) + (dt/Cm)*(((El-V_all(j,i))/Rm) + PSC_all(1,1,j,i) +Ie + poiss(i) + kisi(i));
                    if  V_all(j,i+1)>= Vth+kisi(i)
                        V_all(j,i)=0+kisi(i);
                        V_all(j,i+1)=Vreset +kisi(i);
                    end
end
end
            
%% Make figure with arbitrary selection of synapse and DBS frequency (Figure 4 in the paper)
EI=1;    % Choose 1 for excitatory and 2 for inhibitory
syn=1;   % Choose 1 for F, 2 for D and 3 for P synaptic types
freq=1;  % The dseired DBS frequency to be illustrated
figure;
ax1=subplot(4,1,1);
hold on
title(['Neuron activity pattern from excitatory facilitating synapses at DBS ',num2str(freq),' (Hz)'],'FontWeight','bold')
plot(t,sp(freq,:),'LineWidth',1); zoom xon; ylabel('Input','FontWeight','bold')
ax2=subplot(4,1,2);
plot(t,squeeze(PSC(EI,syn,freq,:)),'LineWidth',1); ylabel('EPSC (nA)','FontWeight','bold')
hold on
plot(t,squeeze(PSC(EI,syn+1,freq,:)),'LineWidth',1); ylabel('EPSC (nA)','FontWeight','bold')
hold on
plot(t,squeeze(PSC(EI,syn+2,freq,:)),'LineWidth',1); ylabel('EPSC (nA)','FontWeight','bold')
ylim([0 15])
ax3=subplot(4,1,3);
plot(t,squeeze(PSC_all(1,1,freq,:)),'LineWidth',1);  
ylabel('Total EPSC (nA)','FontWeight','bold')
ylim([0 15])
ax4=subplot(4,1,4);
plot(t,squeeze(Vn(EI,syn,freq,:)),'LineWidth',1);
hold on
plot(t,Vin,'--','LineWidth',1); 
ylabel({'Model neurn potential';' with and without synapse';' (mv)'},'FontWeight','bold')
xlabel('Time (ms)','FontWeight','bold'); 
ylim([-100 5])
linkaxes([ax1,ax2,ax3,ax4],'x')
%% Compute firing rate of the LIF rneuron without synaptic input:
r_isi_without_syn=(1000/dt)*length(find(Vin(dbsi:dbsf)>=Vth))/((dbsf-dbsi));
disp(['LIF rate without any synaptic connection = ',num2str(r_isi_without_syn),' (Hz)'])
%% Compute firing rate of the LIF neuron with synaptic input:
r_isi_with_syn=(1000/dt)*length(find(Vn(EI,syn,freq,dbsi:dbsf)>=Vth))/((dbsf-dbsi));
disp(['LIF rate with a fraction of synapses during DBS',num2str(freq*10),'Hz = ',num2str(r_isi_with_syn),' (Hz)'])
%% Compute firing rate of the LIF neuron with all synaptic inputs:
r_isi_with_all_syn=(1000/dt)*length(find(V_all(freq,dbsi:dbsf)>=Vth))/((dbsf-dbsi));
disp(['LIF rate with all synapses during DBS', num2str(freq*10),'Hz = ',num2str(r_isi_with_all_syn),' Hz'])
%% Raster plot (Figure 5 in the paper) and PSTH for 130 Hz:
for q=1 
    for sq=1:2
dbsT=round((1000/fdbs(q)/dt));
width=1;
edges=0:width:dbsT;
psth=zeros(1,round(dbsT/width)+1);
figure;
% title(['LIF neuron raster plot at DBS ',num2str(fdbs(q)),' Hz, fidelity= ',num2str(fid*100),'%'],'FontWeight','bold'); 
hold on

if sq==1
for i=dbsi:dbsT:dbsf-dbsT 
    [xx,zz]=find(V_all(q,(i:i+dbsT))>=Vth);
    hh=hist(zz,edges);
    psth=psth+hh;   
scat=scatter(zz*dt,(i*dt/1000)*ones(1,length(xx)),16,'k','filled'); hold on
end
axis ij
axis off
xlim([0 dbsT*dt])
ylim([0 60])

figure;
plot(edges*dt,psth,'LineWidth',1);
xlabel('DBS pulse time period (ms)')
ylabel('Number of spikes')
xlim([0 dbsT*dt]) 
hold on
set(gca,'FontSize',12,'FontWeight','bold')
% saveas(fig,['DBS_',num2str(fdbs(q)),num2str(fid),'fidelity'],'jpg')
end

if sq==2
for i=dbsi:dbsT:dbsf-dbsT 
    [xx,zz]=find(V_all(q,(i:i+dbsT))>=Vth);
    hh=hist(zz,edges);
    psth=psth+hh;   
scat=scatter(zz*dt,(i*dt/1000)*ones(1,length(xx)),121,'k','square','MarkerFaceColor','k'); hold on
end
axis ij
% axis off
xlim([0 dbsT*dt]) 
ylim([42 42.2])
% set(gca,'FontSize',12,'FontWeight','bold')
end
    end
end
toc