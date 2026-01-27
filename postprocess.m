% Analyzed multiple data from different simulations and plots them on the 
% same axis for comparative studies.

%% Import Data

d1 = load("./Data/gillespie_sinusoidal_amp0001_freq001");
d2 = load("./Data/gillespie_sinusoidal_amp0001_freq001");
d3 = load("./Data/gillespie_sinusoidal_amp0001_freq001");
d = [d1 d2 d3];
names = {"A:0.001","A:0.01","A:0.1"};
t_max = 1800;

lambda = {};
spike_time = {};
spike_id = {};
PopE = {};
PopI = {};
Raster_F = {};

Size = {};
Duration = {};
tau = {};
Alpha = {};
Ts = {};
S = {};

Size_F = {};
Duration_F = {};
tau_F = {};
Alpha_F = {};
Ts_F = {};
S_F = {};

Size_Space = {};
Duration_Space = {};
tau_Space = {};
Alpha_Space = {};
Ts_Space = {};
S_Space = {};

for i=1:length(d)
    lambda = {lambda; d(i).lambda};
    spike_time = {spike_time; d(i).spike_time};
    spike_id = {spike_id; d(i).spike_id};
    PopE = {PopE; d(i).PopE};
    PopI = {PopI; d(i).PopI};
    Raster_F = {Raster_F d(i).Raster_F};
    
    Size = {Size d(i).Size};
    Duration = {Duration d(i).Duration};
    tau = {tau d(i).tau};
    Alpha = {Alpha d(i).Alpha};
    Ts = {Ts d(i).Ts};
    S = {S d(i).S};
    
    Size_F = {Size_F d(i).Size_F};
    Duration_F = {Duration_F d(i).Duration_F};
    tau_F = {tau_F d(i).tau_F};
    Alpha_F = {Alpha_F d(i).Alpha_F};
    Ts_F = {Ts_F d(i).Ts_F};
    S_F = {S_F d(i).S_F};
    
    Size_Space = {Size_Space d(i).Size_Space};
    Duration_Space = {Duration_Space d(i).Duration_Space};
    tau_Space = {tau_Space d(i).tau_Space};
    Alpha_Space = {Alpha_Space d(i).Alpha_Space};
    Ts_Space = {Ts_Space d(i).Ts_Space};
    S_Space = {S_Space d(i).S_Space};    
end


%% Raster

figure(1)
xSize = 17; ySize = 12.5;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')

axes('Position',[.065 .8 .22 .18])
x = 0:1000;
for i=1:length(d)
    Pr = exp(-x/lambda{i});
    area(x,Pr,'FaceColor',[.3 .7 .88],'EdgeColor',[.2 .4 .6],'LineWidth',1)
end
legend(names)
xlabel('Distance \itr\rm [\mum]','fontsize',9)
ylabel('Connection prob.','fontsize',9)
text(.29,.6,'P(\itr\rm) = exp(-\itr\rm/\lambda)','units','normalized','fontsize',10)
box off
text(-.2,1.01,'A','fontsize',12,'units','normalized','fontweight','bold')

axes('Position',[.065 .4 .22 .3])
% Raster: spiking activty:
t_lim = 10*60; % 10 min
for i=1:length(d)
    ii = find(spike_times{i}<t_max);
    spT = spike_times{i}(ii);
    spID = spike_ids{i}(ii);
    col = lines(2);
    hold on
    plot(spT(spID<=NE),spID(spID<=NE),'.','markersize',3,'color',col(1,:))
    plot(spT(spID>NE),spID(spID>NE),'.','markersize',3,'color',col(2,:))
    text(.32,.98,sprintf(' \\lambda = %g \\mum \n \\phi = %2.2f',lambda{i},Ampli\i),'units','normalized','EdgeColor','k','BackgroundColor','w')
end
legend(names)
set(gca,'xlim',[t_lim 2*t_lim],'YLim',[0 N],'xtick',[],'xcolor','w','fontsize',9)
box off
xlabel('Time [s]','fontsize',9)
ylabel('neuron ID','fontsize',9)
box off
text(-.2,1.02,'B','fontsize',12,'units','normalized','fontweight','bold')

axes('Position',[.065 .26 .22 .11])
% Population activity:
for i=1:length(d)
    hold on
    plot(bins(1:end-1),PopE,'color',col(1,:))
    plot(bins(1:end-1),PopI,'color',col(2,:))
end
set(gca,'xlim',[t_lim 2*t_lim],'fontsize',9)
box off
xlabel('Time [s]','fontsize',9)
ylabel('Pop. activity','fontsize',9)

axes('Position',[.065 .08 .22 .11])
% Population activity fluorescence:
%Eact = sum(Raster_F(:,1:NE),2)/NE;
%Iact = sum(Raster_F(:,NE+1:end),2)/NI;
%hold on
%plot(time(time<t_lim),Eact(time<t_lim),'color',col(1,:))
%plot(time(time<t_lim),Iact(time<t_lim),'color',col(2,:))
%set(gca,'xlim',[t_lim 2*t_lim],'fontsize',9)
%box off
%xlabel('Time [s]','fontsize',9)
%ylabel('Pop. activity','fontsize',9)

%% Nonspatial Spiking Activity

axes('Position',[.38 .745 .16 .195])
% avalanches spiking activity
% Size distribution:
for i=1:length(d)
    hold on
    Bins = logspace(log10(min(Size{i})),1.1*log10(max(Size{i})),12);
    [x,n]=get_pdfbins(Size{i},Bins);
    plot(x,n,'o-','color',lines(1),'markerfacecolor',lines(1),'markersize',3) 
end
set(gca,'xscale','log','yscale','log','fontsize',9)
xlabel('Avalanche size S','fontsize',9)
ylabel('Probability density','fontsize',9)

xlim = [.8*x(find(~isnan(x),1,'first')) 1.2*x(find(~isnan(x),1,'last'))];

xo = xlim;
for i=1:length(d)
    hold on
    y = xo.^(-tau{i});
    plot(xo,y,'--','color','k','linewidth',1)
    text(.1,.1,sprintf('\\tau = %2.2f',tau{i}),'units','normalized','fontsize',8)
end
set(gca,'xlim',xlim)
text(.1,.26,'P(S) \sim S^{-\tau}','units','normalized','fontsize',9)

box off
text(-.27,1.05,'C','fontsize',12,'units','normalized','fontweight','bold')


axes('Position',[.60 .745 .16 .195])
% avalanches spiking activity
% Duration distribution:
for i=1:length(d)
    hold on
    Bins = logspace(log10(min(Duration{i})),log10(max(Duration{i})),9);
    [x,n]=get_pdfbins(Duration{i},Bins);
    plot(x,n,'o-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)   
end
set(gca,'xscale','log','yscale','log','fontsize',9)
xlabel('Avalanche duration T','fontsize',9)
ylabel('Probability density','fontsize',9)


xlim = [.8*x(find(~isnan(x),1,'first')) 1.2*x(find(~isnan(x),1,'last'))];
xo = xlim;
for i=1:length(d)
    y = xo.^(-Alpha{i});
    plot(xo,y,'--','color','k','linewidth',1)
    text(.1,.1,sprintf('\\alpha = %2.2f',Alpha{i}),'units','normalized','fontsize',8)
end
set(gca,'xlim',xlim)
text(.1,.26,'P(T) \sim T^{-\alpha}','units','normalized','fontsize',9)

box off
text(-.27,1.05,'D','fontsize',12,'units','normalized','fontweight','bold')



axes('Position',[.82 .745 .16 .195])
% avalanches spiking activity
% <S>(T):
for i=1:length(d)
    plot(Ts{i},S{i},'s-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)
end
set(gca,'xscale','log','yscale','log')

xlabel('Avalanche duration T','fontsize',9)
ylabel('<S>(T)','fontsize',9)

text(.4,.26,'<S>(T) \sim T^{1/\sigma\nuz}','units','normalized','fontsize',9)
text(.4,.1,sprintf('\\sigma\\nuz = %2.2f',sigmaNuZ),'units','normalized','fontsize',8)

box off

text(-.27,1.05,'E','fontsize',12,'units','normalized','fontweight','bold')

%% Nonspatial Calcium Events

axes('Position',[.38 .41 .16 .195])
% avalanches fluorescence
% Size distribution:
for i=1:length(d)
    hold on
    Bins = logspace(log10(min(Size_F{i})),1.1*log10(max(Size_F{i})),12);
    [x,n]=get_pdfbins(Size_F{i},Bins);
    plot(x,n,'o-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)
end
set(gca,'xscale','log','yscale','log','fontsize',9)
xlabel('Avalanche size S','fontsize',9)
ylabel('Probability density','fontsize',9)


xlim = [.8*x(find(~isnan(x),1,'first')) 1.2*x(find(~isnan(x),1,'last'))];

xo = xlim;
for i=1:length(d)
    y = xo.^(-tau_F{i});
    plot(xo,y,'--','color','k','linewidth',1)
    text(.1,.1,sprintf('\\tau = %2.2f',tau_F{i}),'units','normalized','fontsize',8)
end
set(gca,'xlim',xlim)
text(.1,.26,'P(S) \sim S^{-\tau}','units','normalized','fontsize',9)

box off
text(-.27,1.05,'F','fontsize',12,'units','normalized','fontweight','bold')


axes('Position',[.60 .41 .16 .195])
% avalanches fluorescence
% Duration distribution:
for i=1:length(d)
    hold on
    Bins = logspace(log10(min(Duration_F{i})),log10(max(Duration_F{i})),12);
    [x,n]=get_pdfbins(Duration_F{i},Bins);
    plot(x,n,'o-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)   
end
set(gca,'xscale','log','yscale','log','fontsize',9)
xlabel('Avalanche duration T','fontsize',9)
ylabel('Probability density','fontsize',9)


xlim = [.8*x(find(~isnan(x),1,'first')) 1.2*x(find(~isnan(x),1,'last'))];
xo = xlim;
for i=1:length(d)
    y = xo.^(-Alpha_F{i});
    plot(xo,y,'--','color','k','linewidth',1)
    text(.1,.1,sprintf('\\alpha = %2.2f',Alpha_F{i}),'units','normalized','fontsize',8)
end
set(gca,'xlim',xlim)
text(.1,.26,'P(T) \sim T^{-\alpha}','units','normalized','fontsize',9)

box off
text(-.27,1.05,'G','fontsize',12,'units','normalized','fontweight','bold')


axes('Position',[.82 .41 .16 .195])
% avalanches fluorescence
% <S>(T):
for i=1:length(d)
    plot(Ts_F{i},S_F{i},'s-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)
end
set(gca,'xscale','log','yscale','log')

xlabel('Avalanche duration T','fontsize',9)
ylabel('<S>(T)','fontsize',9)

text(.4,.26,'<S>(T) \sim T^{1/\sigma\nuz}','units','normalized','fontsize',9)
text(.4,.1,sprintf('\\sigma\\nuz = %2.2f',sigmaNuZ_F),'units','normalized','fontsize',8)

box off
text(-.27,1.05,'H','fontsize',12,'units','normalized','fontweight','bold')

%% Spatial Calcium Events

axes('Position',[.38 .07 .16 .195])
% avalanches fluorescence
% Spatial definition
% Size distribution:
for i=1:length(d)
    hold on
    Bins = logspace(log10(min(Size_Space{i})),1.1*log10(max(Size_Space{i})),10);
    [x,n]=get_pdfbins(Size_Space{i},Bins);
    plot(x,n,'o-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)   
end
set(gca,'xscale','log','yscale','log','fontsize',9)
xlabel('Avalanche size S','fontsize',9)
ylabel('Probability density','fontsize',9)


xlim = [.8*x(find(~isnan(x),1,'first')) 1.2*x(find(~isnan(x),1,'last'))];

xo = xlim;
for i=1:length(d)
    y = xo.^(-tau_Space{i});
    plot(xo,y,'--','color','k','linewidth',1)
end
set(gca,'xlim',xlim)
text(.1,.26,'P(S) \sim S^{-\tau}','units','normalized','fontsize',9)
text(.1,.1,sprintf('\\tau = %2.2f',tau_Space{i}),'units','normalized','fontsize',8)

box off
text(-.27,1.05,'I','fontsize',12,'units','normalized','fontweight','bold')


axes('Position',[.60 .07 .16 .195])
% avalanches fluorescence
% Spatial definition
% Duration distribution:
for i=1:length(d)
    hold on
    Bins = logspace(log10(min(Duration_Space{i})),log10(max(Duration_Space{i})),13);
    [x,n]=get_pdfbins(Duration_Space{i},Bins);
    plot(x,n,'o-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)  
end
set(gca,'xscale','log','yscale','log','fontsize',9)
xlabel('Avalanche duration T','fontsize',9)
ylabel('Probability density','fontsize',9)


xlim = [.8*x(find(~isnan(x),1,'first')) 1.2*x(find(~isnan(x),1,'last'))];
xo = xlim;
for i=1:length(d)
    y = xo.^(-Alpha_Space{i});
    plot(xo,y,'--','color','k','linewidth',1)
    text(.1,.1,sprintf('\\alpha = %2.2f',Alpha_Space{i}),'units','normalized','fontsize',8)
end
set(gca,'xlim',xlim)
text(.1,.26,'P(T) \sim T^{-\alpha}','units','normalized','fontsize',9)

box off
text(-.27,1.05,'J','fontsize',12,'units','normalized','fontweight','bold')


axes('Position',[.82 .07 .16 .195])
% avalanches fluorescence
% Spatial definition
% <S>(T):
for i=1:length(d)
    plot(Ts_Space{i},S_Space{i},'s-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)
end
set(gca,'xscale','log','yscale','log','xlim',[3 50])

xlabel('Avalanche duration T','fontsize',9)
ylabel('<S>(T)','fontsize',9)

text(.4,.26,'<S>(T) \sim T^{1/\sigma\nuz}','units','normalized','fontsize',9)
text(.4,.1,sprintf('\\sigma\\nuz = %2.2f',sigmaNuZ_Space),'units','normalized','fontsize',8)

box off
text(-.27,1.05,'K','fontsize',12,'units','normalized','fontweight','bold')


annotation('textbox',[0.5 .999 0.4 .01],'string','Non-spatial avalanches: spiking activity','fontsize',12,'edgecolor','none','fontweight','bold')
annotation('textbox',[0.5 .669 0.4 .01],'string','Non-spatial avalanches: calcium events','fontsize',12,'edgecolor','none','fontweight','bold')
annotation('textbox',[0.52 .320 0.4 .01],'string','Spatial avalanches: calcium events','fontsize',12,'edgecolor','none','fontweight','bold')

