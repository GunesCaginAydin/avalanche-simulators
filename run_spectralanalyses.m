clear
close
clc
addpath(genpath(pwd));

conntype = 3;
lambda = 200;
[W,Ampli,NE,NI,typ,xyz] = Get_Connectivity(conntype,lambda);
N = NE + NI;

response_fn = @(x) tanh(x).*(x>0);
beta_param = 1;
alpha_param = 0.1;

t_min = 0;
t_max = 1800;
n_batch = 3;

%% Base Case: Constant Input

I = 1e-3*ones(N,t_max*2);
args = {};

spike_times = [];
spike_ids   = [];
time = [];
Fds = [];
sh = 0;
for n = 1:n_batch
    if n==1
    init_state = zeros(2,N);
    init_state(1,:) = rand(1,N)<.3;
    init_state(2,:) = ~init_state(1,:);
    else
    init_state = network_state;
    end
    
    [sp_times,sp_ids,network_state,Irc,tc] = ...
        Gillespie_EImodel(W,response_fn,beta_param,alpha_param,I,t_min,t_max,init_state,args);

    shift = (n-1)*t_max;
    spike_times = [spike_times (sp_times + shift)];
    spike_ids   = [spike_ids sp_ids];

    tauR = .01;
    tauD = .5;
    K = .6; 
    q = 5; 
    Fm = 10; 
    
    resol_F = 1/15;
    [t,fds] = SpikesToFluoresence(sp_times,sp_ids,N,resol_F,tauR,tauD,K,q,Fm);
    time = [time (t+sh)];
    sh = time(end);
    Fds = [Fds;fds];
    clear t fds
end

Tmax = t_max+shift;

E_spike_times = spike_times(spike_ids<=NE);
E_spike_ids = spike_ids(spike_ids<=NE);

I_spike_times = spike_times(spike_ids>NE);
I_spike_ids = spike_ids(spike_ids>NE);

Iintrs = trapz(Irc,2)/size(Irc,2);
Iint = trapz(I,2)/size(I,2);
rcerr = mean(abs(Iintrs - Iint));
tr = linspace(0,size(I,2),size(I,2));

Rates = zeros(1,N);
for i=1:N
Rates(i) = sum(spike_times(spike_ids==i))/(Tmax);
end
MeanRate = mean(Rates);
firing_rates = Rates;

dT = 1;
bins = 0:dT:Tmax;
nbins =length(bins)-1;
Pop = zeros(1,nbins);
PopE = zeros(1,nbins);
PopI = zeros(1,nbins);
for i = 1:nbins
    Pop(i) = sum( spike_times>bins(i) & spike_times<bins(i+1) )/N;
    PopE(i) = sum( E_spike_times>bins(i) & E_spike_times<bins(i+1) )/NE;
    PopI(i) = sum( I_spike_times>bins(i) & I_spike_times<bins(i+1) )/NI;
end
EIratio = PopE./(PopE + PopI);

SD = std(Fds);
F = (Fds - repmat(mean(Fds),[size(Fds,1),1]))./repmat(SD,[size(Fds,1),1]);
Raster_F = F > 3;
Rates_F = sum(Raster_F)/(Tmax);
MeanRate_F = mean(Rates_F);

%% Spectral Analysis
Pop_fft = fft(Pop);
PopE_fft = fft(PopE);
PopI_fft = fft(PopI);
Tns = size(Pop,2);
freqns = 1;
freqrns = (-Tns/2:Tns/2-1)*(freqns/Tns);

Raster_fft = zeros(size(Raster_F'));
for neuron=1:size(Raster_fft,1)
    Raster_fft(neuron,:) = squeeze(fft(Raster_F(:,neuron)'));
end
PopRaster_fft = fft(mean(Raster_F,2)');

Ts = size(Raster_F,1);
freqs = 1;
freqrs = (-Ts/2:Ts/2-1)*(freqs/Ts);

figure(1)
f1 = tiledlayout(2,1);
nexttile
plot(freqrns,fftshift(PopE_fft),'LineWidth',1)
hold on
plot(freqrns,fftshift(PopI_fft),'LineWidth',1)
legend('Exhitatory Activity','Inhibitory Activity')
ylim([-2,4])

nexttile
plot(freqrns,fftshift(Pop_fft),'LineWidth',1)
legend('Whole Population Activity')
ylim([-2,4])

xlabel(f1,'Frequency')
ylabel(f1,'Frequency Spectra')
title(f1,'Frequency Spectra of Non-Spatial Population Activity')

figure(2)
f2 = tiledlayout(2,1);
nexttile
plot(freqrs,fftshift(PopRaster_fft(1,:)),'LineWidth',1)
legend('Raster')
ylim([-5,20])

nexttile
plot(freqrns,fftshift(Pop_fft),'LineWidth',1)
legend('Whole Population Activity')
ylim([-2,4])

xlabel(f1,'Frequency')
ylabel(f1,'Frequency Spectra')
title(f1,'Frequency Spectra of Non-Spatial Population Activity')
%% Raster Plot of Population Activity and Spectral Analysis
figure(3)

% Raster: spiking activty:
subplot(3,1,1)
ii = find(spike_times<Tmax);
spT = spike_times(ii);
spID = spike_ids(ii);
col = lines(2);
hold on
plot(spT(spID<=NE),spID(spID<=NE),'.','markersize',3,'color',col(1,:))
plot(spT(spID>NE),spID(spID>NE),'.','markersize',3,'color',col(2,:))
set(gca,'ylim',[0 N],'xtick',[],'xcolor','w','fontsize',9)
box off
xlabel('Time [s]','fontsize',9)
ylabel('neuron ID','fontsize',9)
text(.32,.98,sprintf(' \\lambda = %g \\mum \n \\phi = %2.2f',lambda,Ampli\i),'units','normalized','EdgeColor','k','BackgroundColor','w')
box off
text(-.2,1.02,'B','fontsize',12,'units','normalized','fontweight','bold')

subplot(3,1,2)
% Population activity:
hold on
plot(bins(1:end-1),PopE,'color',col(1,:))
plot(bins(1:end-1),PopI,'color',col(2,:))
set(gca,'fontsize',9)
box off
xlabel('Time [s]','fontsize',9)
ylabel('Pop. activity','fontsize',9)

subplot(3,1,3)
% Population activity fluorescence:
Eact = sum(Raster_F(:,1:NE),2)/NE;
Iact = sum(Raster_F(:,NE+1:end),2)/NI;
hold on
plot(time,Eact,'color',col(1,:))
plot(time,Iact,'color',col(2,:))
set(gca,'fontsize',9)
box off
xlabel('Time [s]','fontsize',9)
ylabel('Pop. activity','fontsize',9)

%% Discovering Arnold Tongues: Sinusoidal Inputs on a w-A Grid
f = 0.01:0.001:0.05;
A = 0.01:0.001:0.05;
Pop_grid = cell(length(f),length(A));
Raster_grid = cell(length(f),length(A));

for ai=1:length(A)
    for fi=1:length(f)
        t = t_min:t_max;
        I = A + A*sin(f*t);
        args = {};
        
        spike_times = [];
        spike_ids   = [];
        time = [];
        Fds = [];
        sh = 0;
        for n = 1:n_batch
            if n==1
            init_state = zeros(2,N);
            init_state(1,:) = rand(1,N)<.3;
            init_state(2,:) = ~init_state(1,:);
            else
            init_state = network_state;
            end
            
            [sp_times,sp_ids,network_state,Irc,tc] = ...
                Gillespie_EImodel(W,response_fn,beta_param,alpha_param,I,t_min,t_max,init_state,args);
        
            shift = (n-1)*t_max;
            spike_times = [spike_times (sp_times + shift)];
            spike_ids   = [spike_ids sp_ids];
        
            tauR = .01;
            tauD = .5;
            K = .6; 
            q = 5; 
            Fm = 10; 
            
            resol_F = 1/15;
            [t,fds] = SpikesToFluoresence(sp_times,sp_ids,N,resol_F,tauR,tauD,K,q,Fm);
            time = [time (t+sh)];
            sh = time(end);
            Fds = [Fds;fds];
            clear t fds
        end
        
        Tmax = t_max+shift;
        
        E_spike_times = spike_times(spike_ids<=NE);
        E_spike_ids = spike_ids(spike_ids<=NE);
        
        I_spike_times = spike_times(spike_ids>NE);
        I_spike_ids = spike_ids(spike_ids>NE);
        
        Iintrs = trapz(Irc,2)/size(Irc,2);
        Iint = trapz(I,2)/size(I,2);
        rcerr = mean(abs(Iintrs - Iint));
        tr = linspace(0,size(I,2),size(I,2));
        
        Rates = zeros(1,N);
        for i=1:N
        Rates(i) = sum(spike_times(spike_ids==i))/(Tmax);
        end
        MeanRate = mean(Rates);
        firing_rates = Rates;
        
        dT = 1;
        bins = 0:dT:Tmax;
        nbins =length(bins)-1;
        Pop = zeros(1,nbins);
        PopE = zeros(1,nbins);
        PopI = zeros(1,nbins);
        for i = 1:nbins
            Pop(i) = sum( spike_times>bins(i) & spike_times<bins(i+1) )/N;
            PopE(i) = sum( E_spike_times>bins(i) & E_spike_times<bins(i+1) )/NE;
            PopI(i) = sum( I_spike_times>bins(i) & I_spike_times<bins(i+1) )/NI;
        end
        EIratio = PopE./(PopE + PopI);
        
        SD = std(Fds);
        F = (Fds - repmat(mean(Fds),[size(Fds,1),1]))./repmat(SD,[size(Fds,1),1]);
        Raster_F = F > 3;
        Rates_F = sum(Raster_F)/(Tmax);
        MeanRate_F = mean(Rates_F);
    
        Pop_grid{fi,ai} = mean(EIratio);
        Raster_grid{fi,ai} = Raster_F;
    end
end

%% Arnold Tongues on a w-A Grid