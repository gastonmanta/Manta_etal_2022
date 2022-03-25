
%% TSG data


opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = [3, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["fecha", "longitud", "latitud", "salinidad", "temperatura", "fluor", "conductividad", "sigmat", "fecha_instrumento"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "datetime"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "fecha", "InputFormat", "yyyy-MM-dd HH:mm:ss.SSS");
opts = setvaropts(opts, "fecha_instrumento", "InputFormat", "yyyy-MM-dd HH:mm:ss.SSS");

% Import the data
tsg = readtable("/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/ancap_gamboa_termosal.csv", opts);

clear opts

tsg.dates=datenum(tsg.fecha) 

tsg_dates=datenum(tsg.fecha) 

cloro53 = [0.66 0.54 0.12]
z53= [5 20 70]

%10/apr/2016 to 10/may/2016
%% cruise track and instruments location plot, feel happy to improve it
% 
% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/gamboa_cnvs
% load ancap_ctd_stations;


load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/Perfiles_CTD/Excel/ancap_ctd_stations_old.mat')

load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/gamboa_cnvs/ancap_ctd_stations.mat')


%round to 0 the hour min seg
dates=datevec(dates);dates(:,4:end)=0;dates=datenum(dates);


cd /Users/gaston/Documents/Phd/satelite_gamboa
ncdisp('cloro_gamboa.nc')

cloro=ncread('cloro_gamboa.nc','CHL');lcloro=log10(cloro);
lon_cloro=ncread('cloro_gamboa.nc','lon');
time_cloro=ncread('cloro_gamboa.nc','time');time_cloro=double(time_cloro+ datenum('1900-01-01 00:00:00'));
lat_cloro=ncread('cloro_gamboa.nc','lat');

% ncdisp('sst_gamboa.nc')
% 
% lon_sst=ncread('sst_gamboa.nc','lon');
% time_sst=ncread('sst_gamboa.nc','time');time_sst=double(time_sst/86400)++ datenum('1981-01-01 00:00:00');
% lat_sst=ncread('sst_gamboa.nc','lat');
% sst=ncread('sst_gamboa.nc','analysed_sst');sst=sst-273.15;



%% plot
% cd /Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016
% vw = VideoWriter('eddies_gamboa_cloa3','MPEG-4');
% vw.FrameRate = 2; 


% open(vw);

time_step=min(dates):max(dates)

cd /Users/gaston/Documents/Phd/daily_eddies

input_dir='/Users/gaston/Documents/Phd/daily_eddies'

dia_a_plotear=datevec(time_step(1)+13);i=16

% %12
% dia_a_plotear=datevec(time_step(1)+1);i=16-12


% %12
% dia_a_plotear=datevec(time_step(1)+1);i=16-23

%for i=1:length(time_step)



date_of_interest=time_step(1)+i-1;

figure(i);hold on

load (([input_dir,'/adt_',datestr(date_of_interest,'yyyymmdd'),'.mat']))



% yellow el cuadrado amarillo % m_plot(xb, yb, 'y-', 'LineWidth', 2);
% x1=-56;x2=-50;y1=-38;y2=-33.5;
% xb = [x1, x2, x2, x1, x1];yb = [y1, y1, y2, y2, y1];

WL_plot=-60;EL_plot=-44;SL_plot=-41;NL_plot=-29;



figure
%para el 24/4
m_proj('mercator','lon',[-56.8 -48.5], 'lat',[-38 -32.5]);
% 
% %para el 12/4
% m_proj('mercator','lon',[-56.2 -49.1], 'lat',[-38.4 -33.85]);

hold on

m_pcolor(lon_cloro,lat_cloro,lcloro(:,:,i)');shading interp
cmocean('delta');caxis([-1 1])


% m_pcolor(lon_sst,lat_sst,sst(:,:,i)');shading interp
% % colormap(flipud(brewermap(49,'Spectral')));
% colormap(jet);caxis([11 26])

%[ax,~]=m_contfbar([.18 .42],.865,[10 25],[10:1:25],'endpiece','no','axfrac',.05);
% cmocean('delta');caxis([-1 1])


title(datestr(time_step(i),'dd/mm/yy'))

% [~,hc]=m_contourf(lon_sst,lat_sst,sst(:,:,sx)',21,'linestyle','none');
% set(hc,'EdgeColor','none'); hold on
% caxis([16 26])
% t=colorbar; colormap(jet(21))
% caxis_tmp=caxis;
% set(get(t,'ylabel'),'String', 'SST(ï¿½C)','rotation',-270,'Position',[3 (caxis_tmp(1)+caxis_tmp(2))/2]);
% title(datestr(date_of_interest,'dd-mmmm-yyyy'));hold on  

[CS,CH]=m_etopo2('contour',[-2000 -500],'color','g','linewidth',1);

%[~,hc]=m_contour(X,Y,ADT'*100,[-40:5:100],'color',[.6 .6 .6]);
m_contour(X,Y,ADT'*100,[-40:5:100],'color',[.7 .7 .7]);

m_contour(X,Y,ADT'*100,[0 0],'c','linewidth',1);


m_contour(X,Y,ADT'*100,[25 25],'y','linewidth',1);


m_contour(X,Y,ADT'*100,[50 50],'m','linewidth',1);


%[CS,CH]=m_etopo2('contour',[-1000 -200],'color',[.1 .1 .1]);%caxis([-5500 7000]);colormap(flip(gray(7)));  

% m_contour(X,Y,ADT'*100,[-40:5:100],'w');

% yellow el cuadrado amarillo % m_plot(xb, yb, 'y-', 'LineWidth', 2);
% x1=-56;x2=-50;y1=-38;y2=-33.5;
% xb = [x1, x2, x2, x1, x1];yb = [y1, y1, y2, y2, y1];

[lonx,latx]=meshgrid(X,Y);
t=1;
 
m_quiver(lonx(1:t:end,1:t:end),latx(1:t:end,1:t:end),U(1:t:end,1:t:end)',V(1:t:end,1:t:end)','k','autoscalefactor',2);

% [ax,~]=m_contfbar([.18 .42],.865,[10 25],[10:1:25]);
% title(ax,'SST (ºC)')
% 
% [ax,~]=m_contfbar([.1 .33],.71,[-1 1],[-1:.1:1]);
% 
% title(ax,'Log_1_0 Chlorophyll (mg.m^-^3)')


 Cont_CEs=cat (3,Cyclonic_Cell(:,5), Cyclonic_Cell(:,6));
 Cont_AEs=cat (3,Anticyclonic_Cell(:,5), Anticyclonic_Cell(:,6));

% 

id_time=1;
date_num=date_of_interest;

id_not_empty=cellfun(@isempty,Cont_CEs(:,id_time,1))==0;
id_in_cyclo=find(cellfun(@max,Cont_CEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,1))<EL_plot & ...
    cellfun(@max,Cont_CEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,2))<NL_plot);


id_not_empty=cellfun(@isempty,Cont_AEs(:,id_time,1))==0;

id_in_anti=find(cellfun(@max,Cont_AEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,1))<EL_plot & ...
    cellfun(@max,Cont_AEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,2))<NL_plot);

% All_eddy_AR=[AR_Traj{:,2}];

for loop_cyclo_id=id_in_anti'
    m_plot(Cont_AEs{loop_cyclo_id,id_time,1},Cont_AEs{loop_cyclo_id,id_time,2},'r','LineWidth',1)
end



for loop_anti_id=id_in_cyclo'
    m_plot(Cont_CEs{loop_anti_id,id_time,1},Cont_CEs{loop_anti_id,id_time,2},'b','LineWidth',1)
end

% ind=find(dates==dates(1)+i-1);
% 
%  ind=find(dates==date_of_interest)
% 
% m_plot(lon(ind),lat(ind),'og','markerfacecolor','g')

%  ind=find(dates==date_of_interest+2)
% 
%  m_plot(lon(39),lat(39),'dy','markersize',15,'markerfacecolor',[.1 .3 .6],'linewidth',2)
% 
%  m_plot(lon(53),lat(53),'dy','markersize',15,'markerfacecolor',[.9 .8 .5],'linewidth',2)

%m_plot(lon,lat,'*k')

S1 = shaperead('ZEEU.shp');
m_plot([S1.X],[S1.Y], 'k');

% m_scatter(tsg.longitud,tsg.latitud,15,tsg.salinidad)
% %m_scatter(ancapgamboatermosal.longitud,ancapgamboatermosal.latitud,15,ancapgamboatermosal.fluor)
% m_scatter(tsg.longitud,tsg.latitud,15,tsg.temperatura)
% caxis([32 37])

m_gshhs_h('patch',[0.81818 0.77647 0.70909]);

%m_quiver(lonn(1:tt:end,1:tt:end),latt(1:tt:end,1:tt:end),u(1:tt:end,1:tt:end)',v(1:tt:end,1:tt:end)','color',[.3 .3 .3],'autoscalefactor',ii);

 m_grid('linestyle','none');
 

set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(gcf,'color','w');

cd /Users/gaston/Documents/Phd/satelite_gamboa

% savefig('cloro_ancap_snapsho2t')
% print(gcf,'-painters','-depsc2','-r600','cloro_ancap_snapshot2')





% 
 %% 4 panel TSG plot
% figure
% subplot(2,2,3)
% 
% m_proj('mercator','lon',[-56 -50], 'lat',[-38 -33.8]);
% hold on
% [CS,CH]=m_etopo2('contour',[ -5000 -4000 -1000 -200],'color',[.5 .5 .5]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% cmocean('delta');%caxis([0 4])
% hold on
% m_scatter(tsg.longitud,tsg.latitud,15,tsg.fluor)
% 
% caxis([0 .8])
% m_plot([S1.X],[S1.Y], 'k');
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% m_grid
% 
% title('Fluorescence')
% colorbar
% 
% subplot(2,2,1)
% 
% m_proj('mercator','lon',[-56 -50], 'lat',[-38 -33.8]);
% hold on
% [CS,CH]=m_etopo2('contour',[ -5000 -4000 -1000 -200],'color',[.5 .5 .5]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% hold on
% m_scatter(tsg.longitud,tsg.latitud,15,tsg.salinidad)
% %caxis([-1 .2])
% 
% m_plot([S1.X],[S1.Y], 'k');
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% m_grid
% 
% title('Salinity')
% colorbar
% 
% caxis([32 37])
% cmocean('haline');
% 
% 
% 
% subplot(2,2,2)
% m_proj('mercator','lon',[-56 -50], 'lat',[-38 -33.8]);
% hold on
% [CS,CH]=m_etopo2('contour',[ -5000 -4000 -1000 -200],'color',[.5 .5 .5]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% hold on
% m_scatter(tsg.longitud,tsg.latitud,15,tsg.temperatura)
% %caxis([-1 .2])
% colorbar
% cmocean('thermal');%caxis([0 4])
% m_plot([S1.X],[S1.Y], 'k');
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% m_grid
% 
% caxis([10 25])
% title('Temperature')
% 
% set(gcf,'color','w');
% 
% 
% 
% subplot(2,2,4)
% 
% m_proj('mercator','lon',[-56 -50], 'lat',[-38 -33.8]);
% hold on
% [CS,CH]=m_etopo2('contour',[ -5000 -4000 -1000 -200],'color',[.5 .5 .5]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% hold on
%  m_scatter(tsg.longitud,tsg.latitud,15,datenum(tsg.fecha))
% 
% caxis([min(datenum(tsg.fecha)) max(datenum(tsg.fecha))])
% colorbar
% cbdate
% m_plot([S1.X],[S1.Y], 'k');
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% m_grid
% 
% title('Date')
% 
% set(gcf,'color','w');
% 
% 
% % savefig('TSG_gamboa')
% % 
% % print(gcf,'-painters','-depsc2','-r600','TSG_gamboa')
% 



%%

ff=datenum(tsg.fecha)


ind=find(ff>(datenum('04-23-16')) & ff<(datenum('04-25-16')))

tsg_filament=(tsg(ind,:));

aa=(tsg.salinidad./tsg.salinidad);

%aa=(tsg_filament.salinidad./tsg_filament.salinidad);
% sall=tsg_filament.salinidad
% llon=tsg_filament.longitud
% llat=tsg_filament.latitud
% [SA, in_ocean] = gsw_SA_from_SP(sall,aa,llon,llat);


[SA, in_ocean] = gsw_SA_from_SP(tsg.salinidad,aa,tsg.longitud,tsg.latitud);

CT = gsw_CT_from_t(SA,tsg.temperatura,aa);



figure
% m_proj('mercator','lon',[-57 -48], 'lat',[-38 -32.5]);
m_proj('mercator','lon',[-56.8 -48.5], 'lat',[-38 -32.5]);

hold on

  m_plot(lon(39),lat(39),'dy','markersize',15,'markerfacecolor',[.1 .3 .6],'linewidth',2)
  m_plot(lon(53),lat(53),'dy','markersize',15,'markerfacecolor',[.9 .8 .5],'linewidth',2)



%m_scatter(tsg_filament.longitud,tsg_filament.latitud,25,SA,'filled')


m_scatter(tsg.longitud(19680:21400),tsg.latitud(19680:21400),25,SA(19680:21400),'filled')


% cmocean('haline');

caxis([32 37.1])

[ax,~]=m_contfbar([.08 .3],.763,[32 37],[32:0.2: 37]);
title(ax,'S_A (g kg^-^1)')



 m_plot(lon(39),lat(39),'dy','markersize',15,'markerfacecolor','none','linewidth',3)
 m_plot(lon(53),lat(53),'dy','markersize',15,'markerfacecolor','none','linewidth',3)



m_plot([S1.X],[S1.Y], 'k');

m_grid('linestyle','none');


set(findall(gcf,'-property','Fontweight'),'Fontweight','normal')
set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(gcf,'color','w');


colormap(jet)

% 
% print(gcf,'-painters','-depsc2','-r600','sal_snapshot')
% 
% savefig('sal_snapshot')
% 



%% 


for i=1:length(tsg.longitud)-1
dist(i)=sw_dist([tsg.longitud(i+1) tsg.longitud(i)],[tsg.latitud(i+1) tsg.latitud(i)],'km');
end




transectt=cumsum(dist(19680:21400));



%% adcp data
load adcp_clean_gamboa

i=-40;
acros= (u.*cosd(i))+(v*sind(i));
along= -(u*sind(i))+(v*cosd(i));

%%
[ind1 ind11]=find_close_value(timeadcp,tsg_dates(19680))


[ind2 ind22]=find_close_value(timeadcp,tsg_dates(21400))



%acrossup=nmean(acros(1:3,ind11:ind22),1);
acrossup=nmean(acros(1:3,:),1);



figure


plot(transectt,SA(19680:21400))
hold on
yyaxis right
plot(transectt,CT(19680:21400))




figure
%plot(transectt,tsg.fluor(19680:21400))

plot(latadcp(ind11:ind22),acrossup(ind11:ind22))

hold on

plot(latadcp(ind11:ind22),acrossup(ind11:ind22))


hold on
yyaxis right
plot(tsg.latitud(19680:21400),SA(19680:21400))


figure
%plot(transectt,tsg.fluor(19680:21400))


ssa2=1;ssa=1;

plot(latadcp(ind11:ind22),acrossup(ind11:ind22))
hold on

plot(latadcp(ind11+ssa2:ind11+ssa),acrossup(ind11+ssa2:ind11+ssa),'r')




nmean(acrossup(ind11+ssa2:ind11+ssa))


nmean(acrossup(ind11+ssa2:ind11+ssa))



alongsup=nmean(along(1:3,:),1);


speed=sqrt(v.*v+u.*u);



speedsup=nmean(speed(1:3,:),1);



nmean(speedsup(ind11+ssa2:ind11+ssa))










jq=1;kk=1;


x1 = (1:length(tsg.longitud(28015-jq*kk:30758)));


for i=1:length(x1)-1

distance(i)=sw_dist([tsg.latitud(i) tsg.latitud(i+1)],[tsg.longitud(i) tsg.longitud(i+1)],'km');

end



distanc=ncumsum(distance);x1=[distanc distanc(end)];


x1=[distanc distanc(end) distanc(end)]';

acrossup(ind11+ssa2:ind11+ssa)

x1=tsg.latitud(19680:21400);


x2 = x1;

x3=x2;

y1 = CT(19680:21400);
y2 = SA(19680:21400);

%y2 = tsg.salinidad(28015-jq*kk:30758+jq);
y3 = movmean(tsg.fluor (19680:21400),10);

x4 = latadcp((ind11:ind22));


y4 = movmean(acrossup(ind11:ind22),10);
y6 = movmean(alongsup(ind11:ind22),10);
y5 = movmean(speedsup(ind11:ind22),10);

% ylat = tsg.latitud(28015-jq*kk:30758+jq);
% ylon = tsg.longitud(28015-jq*kk:30758+jq);

% x1=ylat;
% 
% zz=10;
% 
% % y1=movmedian(y1,zz,'omitnan');
% % y2=movmedian(y2,zz,'omitnan');
% % y3=movmedian(y3,zz,'omitnan');
% 
% y4=movmedian(y4,zz,'omitnan');
% y5=movmedian(y5,zz,'omitnan');
% 
% y4=movmedian(y4,zz,'omitnan');
% y5=movmedian(y5,zz,'omitnan');



figure
[ax,hlines,fh] = jzplotys({x1,y1,x1,y2,x1,y3,x4,y4,x4,y5},[1 1 1 1 1],[x1(end) x1(1)],100);



legend(cat(1,hlines{:}),'Temperature','Salinity','Fluorescene','Across-slope vel.','Along-slope vel.','orientation','horizontal','location',[.2 .05 .01 .01],'box','off')

ylabel(ax(1),'\Theta (ºC)');

ylabel(ax(3),'S_A (g kg^-^1)');

ylabel(ax(5),'Fluorescene (mg m^-^3)');

ylabel(ax(7),'Across slope vel. (m s^-^1)');

ylabel(ax(9),'Speed (m s^-^1)');

set(gcf,'color','w');set(findall(gcf,'-property','FontSize'),'FontSize',12)

set(gcf,'color','w');set(findall(gcf,'-property','FontSize'),'FontSize',12.5)

set(gcf,'color','w');set(findall(gcf,'-property','linewidth'),'linewidth',0.8)































sw_dist([-36.287 -35.866],[-52.352 -51.775],'km')

latadcp(ind11:ind22)

[dmin,imin]=find_close_valuexy(-52.352,-36.287,lonadcp(ind11:ind22),latadcp(ind11:ind22))

lonadcp(ind11:ind22),latadcp(ind11:ind22)

[dmin,imin]=pdist2([-52.352 -36.287],[lonadcp(ind11:ind22) latadcp(ind11:ind22)],'euclidean','smallest',1)

[ss, ssa]=min(dmin)


[dmin,imin]=pdist2([-51.775 -35.866],[lonadcp(ind11:ind22) latadcp(ind11:ind22)],'euclidean','smallest',1)

[ss2, ssa2]=min(dmin)



fi

ssa


ff=datenum(tsg.fecha)


ind=find(ff>(datenum('04-23-16')) & ff<(datenum('04-25-16')))

tsg_filamentt=(tsg(ind,:));



aa=(tsg_filamentt.salinidad./tsg_filamentt.salinidad);

sall=tsg_filamentt.salinidad
llon=tsg_filamentt.longitud
llat=tsg_filamentt.latitud

 [SA, in_ocean] = gsw_SA_from_SP(sall,aa,llon,llat);




% 
% 
% figure;
% scatter(tsg.salinidad,tsg.temperatura,5,tsg.fluor,'filled')
% colorbar
% 
% xlim([34 37])
% 
% 



figure
m_proj('mercator','lon',[-57 -48], 'lat',[-38 -32.5]);
hold on

% 
% m_scatter(tsg_filament.longitud,tsg_filament.latitud,35,log10(tsg_filament.fluor),'filled')
% 
% cmocean('delta');
% 
% caxis([-1 -.5])



m_scatter(tsg_filamentt.longitud,tsg_filamentt.latitud,25,SA,'filled')

cmocean('haline');

caxis([32 37.1])

[ax,~]=m_contfbar([.05 .34],.763,[32 37],[32:0.2: 37]);
title(ax,'S_A (kg.m^-^3)')

m_plot([S1.X],[S1.Y], 'k');

m_grid('linestyle','none');




figure;
plot(tsg_filamentt.fecha,tsg_filamentt.salinidad); hold on
yyaxis right
plot(tsg_filamentt.fecha,tsg_filamentt.fluor); hold on
% hold on
% plot(tsg_filamentt.fecha,tsg_filamentt.temperatura); hold on

hold on
plot(tsg_filamentt.fecha,tsg_filamentt.temperatura); hold on
ylim([19 25])



for i=1:length(tsg_filamentt.fecha)-1
distance(i)=sw_dist([tsg_filamentt.latitud(i) tsg_filamentt.latitud(i+1)],[tsg_filamentt.longitud(i) tsg_filamentt.longitud(i+1)],'km'); 
end


cum_dist=cumsum(distance);


% for i=1:length(tsg_filamentt.fecha)-1
% gradfluo(i)=(tsg_filamentt.fluor(i+1)-tsg_filamentt.fluor(i))/distance(i); 
% end
% 

figure;
% plot(tsg_filamentt.fecha,tsg_filamentt.salinidad); hold on
% yyaxis right
%  plot(tsg_filamentt.fecha,tsg_filamentt.fluor); hold on


 figure('Renderer', 'painters', 'Position', [150 150 570 200])

 plot(cum_dist,tsg_filamentt.salinidad(1:end-1)); hold on
 yyaxis right
  plot(cum_dist,tsg_filamentt.fluor(1:end-1)); hold on

% longitud del filamento 1720 2370
sw_dist([tsg_filamentt.latitud(1720) tsg_filamentt.latitud(2370)],[tsg_filamentt.longitud(1720) tsg_filamentt.longitud(2370)],'km')

% hold on
% plot(tsg_filamentt.fecha,tsg_filamentt.temperatura); hold on



% 
% hold on
% plot(tsg_filamentt.fecha,tsg_filamentt.temperatura); hold on
% ylim([19 25])




figure('Renderer', 'painters', 'Position', [150 150 570 200])
plot(cum_dist,tsg_filamentt.salinidad(1:end-1),'b'); hold on

ylabel ('S_A (g.kg^-^1)')

 yyaxis right
  plot(cum_dist,tsg_filamentt.fluor(1:end-1)); hold on
 xlim([120 267])

ylabel ('Fluorescence (mg.m^-^3)')


set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',13)
xlabel ('Distance (km)')

cd /Users/gaston/Documents/Phd/satelite_gamboa

grid on
% 
% print(gcf,'-painters','-depsc2','-r600','tsg_filament_timeseries')
% savefig('tsg_filament_timeseries')



figure('Renderer', 'painters', 'Position', [150 150 570 200])
plot(cum_dist,tsg_filamentt.temperatura(1:end-1),'b'); hold on

ylabel (' \Theta (ºC)')

 yyaxis right
  plot(cum_dist,tsg_filamentt.fluor(1:end-1)); hold on
 xlim([120 267])

ylabel ('Fluorescence (mg.m^-^3)')


set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',13)
xlabel ('Distance (km)')

% savefig('tsg_filament_timeseries2')
%  print(gcf,'-painters','-depsc2','-r600','tsg_filament_timeseries2')



%% CTD profiles in the filament
load('/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa/gamboa_cnvs/ancap_ctd_stations.mat')

% creates a matrix with each variable
stations=1:length(depth); ppres=0:max(cellfun(@max, depth));

max_depth=cellfun(@max, depth);

%create a empty matrix
temp_new=NaN(length(ppres),length(stations));
sal_new=NaN(length(ppres),length(stations));
ox_new=NaN(length(ppres),length(stations));
flu2_new=NaN(length(ppres),length(stations));
turb_new=NaN(length(ppres),length(stations));
pres_new=NaN(length(ppres),length(stations));

%put values at thier corresponding depth

    for kk=1:length(stations)
  % for pp=press(1):length(press) 

    press=round(depth{kk});
    tempp=(temp{kk});
     sall=(sal{kk});
       oxx=(oxy{kk});
         flu2=(fluo{kk});
          turbb=(turb{kk});
%      
for pp=press(1):length(press) 
% temp_new(pp,kk)=tempp(pp);%clear tempp;clear press;
%ppp=pp;
 ppp=pp-press(1)+1;
       temp_new(pp,kk)=tempp(ppp);
        sal_new(pp,kk)=sall(ppp); 
        ox_new(pp,kk)=oxx(ppp);
           flu2_new(pp,kk)=flu2(ppp);
          turb_new(pp,kk)=turbb(ppp);
         pres_new(pp,kk)=press(ppp);
           
    end
    end


temp_new(1400:1800,81)=movmedian(temp_new(1400:1800,81),20);
temp_new(1300:1400,58)=movmedian(temp_new(1300:1400,58),20);

sal_new(1400:1800,81)=movmedian(sal_new(1400:1800,81),20);
sal_new(1300:1400,58)=movmedian(sal_new(1300:1400,58),20);
sal_new(760:920,35)=sal_new(758,35);

ox_new(1400:1800,81)=movmedian(ox_new(1400:1800,81),20);
ox_new(1300:1400,58)=movmedian(ox_new(1300:1400,58),20);
ox_new(760:920,35)=ox_new(758,35);

flu2_new(:,35)=NaN;ox_new(:,35)=NaN;

% 
% sal_new(sal_new<0.5)=NaN;sal_new(sal_new>39)=NaN;
% ox_new(ox_new<0.01)=NaN;ox_new(ox_new>10)=NaN;
% temp_new(temp_new<-3)=NaN;temp_new(temp_new>32)=NaN;

theta=sw_ptmp(sal_new,temp_new,pres_new,0);

% theta(theta<-3)=NaN;theta(theta>32)=NaN;

lon3=repmat(lon,4219,1);
lat3=repmat(lat,4219,1);

[SA, in_ocean] = gsw_SA_from_SP(sal_new,pres_new,lon,lat);

CT = gsw_CT_from_t(SA,temp_new,pres_new);

%gamma_GP = gamma_GP_from_SP_pt(SA,CT,pres_new,lon,lat);


theta=sw_ptmp(sal_new,temp_new,pres_new,1);


gamma= gamma_GP_from_SP_pt(sal_new,theta,pres_new,lon,lat);





S1 = shaperead('ZEEU.shp');
S2 = shaperead('eez_v11.shp');
S3 = shaperead('sudamerica.shp');

%%


ind=find(dates<(datenum('4-25-2016')) & dates>(datenum('4-17-2016')))

figure;plot(temp_new(:,ind),-pres_new(:,ind));ylim([-1000 0])


%filamento: 22(38-54-63) y 24(39-53)



filament=[39 53];


zz=-200;

figure('Renderer', 'painters', 'Position', [150 150 280 400])
ha = tight_subplot(1,3,[.01 .01],[.1 .01],[.13 .01])

axes(ha(1));box on; grid on
hold on
 
plot(SA(7,39),0,'dy','markersize',8,'markerfacecolor',[.1 .3 .6],'linewidth',1)
plot(SA(7,53),0,'dy','markersize',8,'markerfacecolor',[.6 .5 .2],'linewidth',1)

plot(SA(:,39),-pres_new(:,39),'color',[.1 .3 .6],'linewidth',.7);
plot(SA(:,53),-pres_new(:,53),'color',[.6 .5 .2],'linewidth',.7);

ylim([zz 0])
xlabel ('S_A (g.kg^-^1)')

yticks([-200:20:0])

yticklabels({'200','180','160','140','120','100','80','60','40','20','0'})

xticks([33:1:37])
xticklabels([33:1:37])

xlim([32.5 36.9])
ylabel('Depth (m)')



axes(ha(2));box on; grid on
hold on

plot(temp_new(7,39),0,'dy','markersize',8,'markerfacecolor',[.1 .3 .6],'linewidth',1)
plot(temp_new(7,53),0,'dy','markersize',8,'markerfacecolor',[.6 .5 .2],'linewidth',1)

plot(temp_new(:,39),-pres_new(:,39),'color',[.1 .3 .6],'linewidth',.7);
plot(temp_new(:,53),-pres_new(:,53),'color',[.6 .5 .2],'linewidth',.7);

ylim([zz 0])

yticklabels('')

xticks([10:4:22])

xticklabels([10:4:22])

xlabel (' \Theta (ºC)')


axes(ha(3));box on; grid on
hold on

flu2_neww=flu2_new-(nmean(flu2_new(200:end,:)));

plot(flu2_neww(7,39),0,'dy','markersize',8,'markerfacecolor',[.1 .3 .6],'linewidth',1)
plot(flu2_neww(7,53),0,'dy','markersize',8,'markerfacecolor',[.6 .5 .2],'linewidth',1)

plot(flu2_neww(:,39),-pres_new(:,39),'color',[.1 .3 .6],'linewidth',.7);
plot(flu2_neww(:,53),-pres_new(:,53),'color',[.6 .5 .2],'linewidth',.7);



plot(cloro53,-z53,'*','color',[.6 .5 .2],'markersize',10)


ylim([zz 0])
yticklabels('')

xticks([0.2:.2:.6])
xticklabels([0.2:.2:.6])

xlabel ('Fluo (mg.m^-^3)')

xlim([0 0.72])

% axes(ha(4))
% plot(turb_new(:,filament),-pres_new(:,filament));ylim([-300 0])


set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',12)
% 
% savefig('cloro_snapshot_profiles')
% 
% print(gcf,'-painters','-depsc2','-r600','cloro_snapshot_profiles')
% 








%% video plot figure S2

WL_plot=-60;EL_plot=-44;SL_plot=-41;NL_plot=-29;

cd /Users/gaston/Documents/Phd/daily_eddies

input_dir='/Users/gaston/Documents/Phd/daily_eddies'

dia_a_plotear=datevec(time_step(1)+13);i=16

% %12
% dia_a_plotear=datevec(time_step(1)+1);i=16-12


% %12
% dia_a_plotear=datevec(time_step(1)+1);i=16-23

%for i=1:length(time_step)




input_dir='/Users/gaston/Documents/adts_eddies'

output_dir='/Users/gaston/Documents/eddies_zeeu/eddies_cloro/'; % outputfigure

close all

 %2016           4           9           0           0           0

ttsg=datevec(tsg_dates);ttsg=datenum(ttsg(:,1:3));

for i=1:32

date_of_interest=datenum('2016-04-09')+(i)-1;

figure(i);hold on

load (([input_dir,'/adt_',datestr(date_of_interest,'yyyymmdd'),'.mat']))

m_pcolor(lon_cloro,lat_cloro,lcloro(:,:,i)');shading interp
cmocean('delta');caxis([-1 1])




m_proj('mercator','lon',[-56.8 -48.5], 'lat',[-38 -32.5]);
hold on

title(datestr(date_num,'dd/mm/yy'))

m_contour(X,Y,ADT'*100,[-40:5:100],'color',[.7 .7 .7]);

%[CS,CH]=m_etopo2('contour',[-1000 -200],'color',[.1 .1 .1]);%caxis([-5500 7000]);colormap(flip(gray(7)));  

[CS,CH]=m_etopo2('contour',[-2000 -500],'color','g','linewidth',1);

[lonx,latx]=meshgrid(X,Y);
t=1;
 
m_quiver(lonx(1:t:end,1:t:end),latx(1:t:end,1:t:end),U(1:t:end,1:t:end)',V(1:t:end,1:t:end)','k','autoscalefactor',2);

 Cont_CEs=cat (3,Cyclonic_Cell(:,5), Cyclonic_Cell(:,6));
 Cont_AEs=cat (3,Anticyclonic_Cell(:,5), Anticyclonic_Cell(:,6));


id_time=1;
date_num=date_of_interest;

id_not_empty=cellfun(@isempty,Cont_CEs(:,id_time,1))==0;
id_in_cyclo=find(cellfun(@max,Cont_CEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,1))<EL_plot & ...
    cellfun(@max,Cont_CEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_CEs(id_not_empty,id_time,2))<NL_plot);


id_not_empty=cellfun(@isempty,Cont_AEs(:,id_time,1))==0;

id_in_anti=find(cellfun(@max,Cont_AEs(id_not_empty,id_time,1))>WL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,1))<EL_plot & ...
    cellfun(@max,Cont_AEs(id_not_empty,id_time,2))>SL_plot & cellfun(@min,Cont_AEs(id_not_empty,id_time,2))<NL_plot);

% All_eddy_AR=[AR_Traj{:,2}];

for loop_cyclo_id=id_in_anti'
    m_plot(Cont_AEs{loop_cyclo_id,id_time,1},Cont_AEs{loop_cyclo_id,id_time,2},'r','LineWidth',1)
end



for loop_anti_id=id_in_cyclo'
    m_plot(Cont_CEs{loop_anti_id,id_time,1},Cont_CEs{loop_anti_id,id_time,2},'b','LineWidth',1)
end

S1 = shaperead('ZEEU.shp');
m_plot([S1.X],[S1.Y], 'k');

m_gshhs_h('patch',[0.81818 0.77647 0.70909]);

%m_quiver(lonn(1:tt:end,1:tt:end),latt(1:tt:end,1:tt:end),u(1:tt:end,1:tt:end)',v(1:tt:end,1:tt:end)','color',[.3 .3 .3],'autoscalefactor',ii);



ind=find(ttsg==date_of_interest);

% m_scatter(tsg.longitud(ind),tsg.latitud(ind),5,tsg.salinidad(ind));
% caxis([34 37])

 m_plot(tsg.longitud(ind),tsg.latitud(ind),'.m','markersize',10);
% caxis([34 37])


 m_grid('linestyle','none');
 

set(findall(gcf,'-property','FontSize'),'FontSize',13)
set(gcf,'color','w');


print( gcf, '-dpng','-r150',[output_dir,datestr(date_num,'yyyymmdd')])

end



% print(gcf,'-dpng','-r600','cloro_ancap_snapshot')







