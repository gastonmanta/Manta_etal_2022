%% ctd data

 load('ancap_ctd_stations.mat')
load('ancap_ctd_stations2.mat')

% cd '/Users/gaston/Desktop'


sal_new= sgolayfilt(sal_new, 4, 9);sal_new=movmean(sal_new,5,1);
temp_new = sgolayfilt(temp_new, 4, 9);temp_new=movmean(temp_new,5,1);
ox_new = sgolayfilt(ox_new, 4, 9);ox_new=movmean(ox_new,5,1);


[SA, in_ocean] = gsw_SA_from_SP(sal_new,pres_new,lon,lat);

CT = gsw_CT_from_t(SA,temp_new,pres_new);

%gamma_GP = gamma_GP_from_SP_pt(SA,CT,pres_new,lon,lat);


theta=sw_ptmp(sal_new,temp_new,pres_new,1);


gamma= gamma_GP_from_SP_pt(sal_new,theta,pres_new,lon,lat);


lon=lon';lat=lat';
lon3=repmat(lon,4219,1);
lat3=repmat(lat,4219,1);


S1 = shaperead('ZEEU.shp');

theta=sw_ptmp(sal_new,temp_new,pres_new,0);


rho = gsw_rho(SA,CT,pres_new);

rr=rho/1000;

% ox is in mg/l *0.7 to ml/l, *44.661 to μmol/l and*rr to μmol/kg
oxy_new=ox_new*44.661.*rr*0.7;


for j=1:length(rho)-1

for i=1:82

drho(j,i)=rho(j,i)-rho(j+1,i);

end
end


spicy=gsw_spiciness0(SA,CT);
lon3=repmat(lon,4219,1);
lat3=repmat(lat,4219,1);

 
% figure
% scatter3(lon3(:),lat3(:),-pres_new(:),25,spicy(:),'filled');%hold on
% 
% zlim([-500 0])
% 
% caxis([0 5])
% 
% view(18,40)
% 
% 
% figure;
% 
% plot(SA,CT)
% hold on
% 
% scatter(SA(:),CT(:),5,spicy(:),'filled')
% 
% 
% caxis([0 5])
% 
% 
% xlim([33.5 37])

%% adcp data
load adcp_clean_gamboa

i=-40;
acros= (u.*cosd(i))+(v*sind(i));
along= -(u*sind(i))+(v*cosd(i));

%% TSG data


opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = [2, Inf];
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

tsg_dates=datenum(tsg.fecha);

[SAtsg, in_ocean] = gsw_SA_from_SP(tsg.salinidad,1,tsg.longitud ,tsg.latitud);

CTtsg = gsw_CT_from_t(SAtsg,tsg.temperatura,1);


tsg_wms=tsg.salinidad

ind=find(tsg_wms<=33.5);tsg_wms(ind)=1;


ind=find(tsg_wms>33.5 & tsg_wms<=34);tsg_wms(ind)=2;



figure
m_proj('mercator','lon',[-56 -50], 'lat',[-38 -33.8]);
m_scatter(tsg.longitud,tsg.latitud,5,tsg_wms,'filled')
colorbar
caxis([1 3])
hold on
[CS,CH]=m_etopo2('contour',[-2000 -500],'color','g','linewidth',1);
m_grid
% 
% figure
% m_proj('mercator','lon',[-56 -50], 'lat',[-38 -33.8]);
% m_scatter(tsg.longitud,tsg.latitud,5,tsg.salinidad,'filled')
% colorbar
% caxis([33.5 34])
% hold on
% [CS,CH]=m_etopo2('contour',[-2000 -500],'color','g','linewidth',1);
% m_grid

%% meteo

opts = delimitedTextImportOptions("NumVariables", 9);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
%opts.VariableNames = ["YYYYMMDDThhmmsssss", "Longitudedegrees_east", "Latitudedegrees_north", "WindSpeedMetresPerSecond", "WindDirectionDegrees", "AirTemperatureDegreesCelsius", "HumidityPercent", "SolarRadiationWattsPerSquareMetre", "AirPressureHectopascals"];
opts.VariableNames = ["YYYYMMDDThhmmsssss", "lon", "lat", "spd", "dir", "temp", "hum", "rad", "pres"];

opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "YYYYMMDDThhmmsssss", "InputFormat", "yyyy-MM-dd HH:mm:ss.SSS");

% Import the data
meteo = readtable("/Users/gaston/Documents/Phd/Datos_CTD_phd/Gamboa_ANCAP_May_2016/gamboa_meteo_ancap.csv", opts);
clear opts

[meteo.uw,meteo.vw]=uvposta(meteo.spd,meteo.dir);

TT = table2timetable(meteo);

tair = retime(TT,'hourly','mean');tair(1:9,:)=[];tair(764:end,:)=[];

% 2016-04-08 11:00:00.000

meteoh = retime(TT,'hourly','mean');

TT = table2timetable(meteo);

meteom = retime(TT,'minutely','mean');tair(1:526,:)=[];tair(45722:end,:)=[];
% % 2016-04-08 11:00:00.000
% TT2 = table2timetable(tsg);
% 
% TT2 = table2timetable(tsg);
% 
% tsst = retime(TT2,'minutely','mean');
% 
% deltat=tair.temp-tsst.temperatura;

wspdspd=movmean(tair.spd,30);
%% legs

deep_leg=[11 34 35 58 59 82];
deep_leg=[82 59 58 35 34 11];
deep_leg2=[81 60 57 36 33 12];
deep_leg12=[80 61 56 37 32 13];

deep_leg13=[79 62 55 38 31 14 10];

shallowest_deep_leg=flip([8 16 28 41 52 65 76]);

shallow_deep_leg=flip([9 15 29 40 53 64 77]);

mid_deep_leg=flip([10 14 30 39 54 63 78]);


deep_leg=[82 59 58 35 34 11];

deep_leg2=[81 60 57 36 33 12];

deep_leg12=[80 61 56 37 32 13];

deep_leg13=[79 62 55 38 31 14 10];

shallowest_deep_leg=flip([8 16 28 41 52 65 76]);

shallow_deep_leg=([9 15 29 40 52 65 76]);

mid_deep_leg=flip([10 14 30 39 53 64 77]);


%% stations plot

figure
plot(tsg.longitud,tsg.latitud,'.','markersize',10); hold on
plot(lonadcp,latadcp,'.y','markersize',5); hold on

stations=1:82;
for k=1:length(stations)
plot(lon(k),lat(k),'r+')
text(lon(k),lat(k),num2str(stations(k)),'fontsize',13,'color','r')
end

legend('TSG','ADCP','CTD')

title('The data')
set(gcf,'color','w');set(findall(gcf,'-property','FontSize'),'FontSize',14)

%% sections

figure; hold on
plot(lonadcp,latadcp,'.k','markersize',4); hold on


[minVal,idx]=find_close_value(lonadcp,-51.6314)


[minVal,idx]=find_close_value(lonadcp,-53.7876)


hold on; plot(lonadcp(8497:9584),latadcp(8497:9584),'*r')


[minVal,idx]=find_close_value(lonadcp,-51.8779)


[minVal,idx]=find_close_value(lonadcp,-54.1981)


hold on; plot(lonadcp(9693:10763),latadcp(9693:10763),'*g')


[minVal,idx]=find_close_value(lonadcp,-53.3978)


[minVal,idx]=find_close_value(lonadcp,-51.2319)

hold on; plot(lonadcp(6938:8020),latadcp(6938:8020),'*y')

hold on; plot(lonadcp(8020:8125),latadcp(8020:8125),'*b')


[minVal,idx]=find_close_value(lonadcp,-51.61692)


[minVal,idx]=find_close_value(lonadcp,-53.5428)

[minVal,idx]=find_close_value(latadcp,-37.2062)

hold on; plot(lonadcp(823:1583),latadcp(823:1583),'*m')


% % acros north shelf
% [minVal,idx]=find_close_value(latadcp,-34.1995)
% 
% [minVal,idx]=find_close_value(latadcp,-34.7589)
% 
% hold on; plot(lonadcp(12882:13210),latadcp(12882:13210),'*y')



[minVal,idx]=find_close_value(latadcp,-36.4664)


[minVal,idx]=find_close_value(lonadcp,-54.0383)

[minVal,idx]=find_close_value(lonadcp,-53.431)

hold on; plot(lonadcp(8030:8334),latadcp(8030:8334),'*c')


[minVal,idx]=find_close_value(lonadcp,-54.7978)


[minVal,idx]=find_close_value(lonadcp,-54.7978)

hold on; plot(lonadcp(16508:17523),latadcp(16508:17523),'*g')



[minVal,idx]=find_close_value(latadcp,-36.6421)

[minVal,idx]=find_close_value(latadcp,-36.2076)

hold on; plot(lonadcp(10981:11493),latadcp(10981:11493),'*k')


[minVal,idx]=find_close_value(latadcp,-37.3572)

[minVal,idx]=find_close_value(latadcp,-35.9888)

[minVal,idx]=find_close_value(latadcp,-34.3541)

[minVal,idx]=find_close_value(latadcp,-34.542)

[minVal,idx]=find_close_value(latadcp,-36.2004)

[minVal,idx]=find_close_value(latadcp,-36.1955)




[minVal,idx]=find_close_value(latadcp,-34.5431)


[minVal,idx]=find_close_value(latadcp,-36.4605)


[minVal,idx]=find_close_value(latadcp,-36.4816)


[minVal,idx]=find_close_value(latadcp,-36.4816)

[minVal,idx]=find_close_value(latadcp,-35.2059)


[minVal,idx]=find_close_value(latadcp,-36.2017)


[minVal,idx]=find_close_value(lonadcp,-53.9697)

[minVal,idx]=find_close_value(lonadcp,-53.5548)


[minVal,idx]=find_close_value(lonadcp,-54.0396)

[minVal,idx]=find_close_value(lonadcp,-54.0396)


[minVal,idx]=find_close_value(lonadcp,-54.5548)


[minVal,idx]=find_close_value(lonadcp,-53.3978)

[minVal,idx]=find_close_value(lonadcp,-53.7438)

[minVal,idx]=find_close_value(latadcp,-36.9546)


[minVal,idx]=find_close_value(latadcp,-36.2149)


[minVal,idx]=find_close_value(lonadcp,-53.9848)

[minVal,idx]=find_close_value(lonadcp,-53.5448)

[minVal,idx]=find_close_value(lonadcp,-53.4229)


% 
% uu16=squeeze(nmean(u16(1:5,:)));vv16=squeeze(nmean(v16(1:5,:)));
% 
% uu8=squeeze(nmean(u8(1:10,:)));vv8=squeeze(nmean(v8(1:10,:)));
% 

%lon8(1:1105)=[];lat8(1:1105)=[];u8(1:1105)=[];u8(1:1105)=[];

%lon8(1:1105)=[];lat8(1:1105)=[];u8(1:1105)=[];u8(1:1105)=[];

% lon8(1:1105)=[];lat8(1:1105)=[];u8(1:1105)=[];u8(1:1105)=[];



[minVal,idx]=find_close_value(lonadcp,-53.0974)

[minVal,idx]=find_close_value(lonadcp,-52.9884)






[minVal,idx]=find_close_value(tsg.longitud,-51.6399)


[minVal,idx]=find_close_value(tsg.longitud,-53.7711)









%%
acros_sup=nmean(acros(1:3,:));
along_sup=nmean(along(1:3,:));

v_sup=nmean(v(1:3,:));
u_sup=nmean(u(1:3,:));




acros_100=nmean(acros(4:6,:));
along_100=nmean(along(4:6,:));

v_100=nmean(v(1:3,:));
u_100=nmean(u(1:3,:));




acros_sup_redi=griddata(lonadcp(8497:9584),latadcp(8497:9584),acros_sup(8497:9584),tsg.longitud(28015:30758),tsg.latitud(28015:30758));
along_sup_redi=griddata(lonadcp(8497:9584),latadcp(8497:9584),along_sup(8497:9584),tsg.longitud(28015:30758),tsg.latitud(28015:30758));

v_sup_redi=griddata(lonadcp(8497:9584),latadcp(8497:9584),v_sup(8497:9584),tsg.longitud(28015:30758),tsg.latitud(28015:30758));
u_sup_redi=griddata(lonadcp(8497:9584),latadcp(8497:9584),u_sup(8497:9584),tsg.longitud(28015:30758),tsg.latitud(28015:30758));



%% section

lat3=repmat(lat,4219,1);
lon3=repmat(lon,4219,1);
tt=6;
ylim([0 500])

lonadcp(8497:9584);

latadcp(8497:9584);



latred=-37.3:0.005:-35.0;


latred=-37.3:0.005:-35.0;
zx=1:500;


[llatredx,zxx]=meshgrid(latred,zx);


alongg=griddata(latadcp(8497:9584),-z,along(:,8497:9584),latred',zx);


across=griddata(latadcp(8497:9584),-z,acros(:,8497:9584),latred',zx);

shallow_deep_leg=[9 15 29 40 53 64 77];



acrossm=nmean(nmean(across(70:130,201:241)))


saa=griddata(lat3(:,shallow_deep_leg),pres_new(:,shallow_deep_leg),SA(:,shallow_deep_leg),latred',zx);


saam=movmean(saa,10,2);saam=movmean(saam,5,1);saam=movmean(saam,10,2);


saam=fillmissing(saam,'nearest',2);


acrosf=fillmissing(across,'nearest');


alongf=fillmissing(alongg,'nearest');


theta=movmean(CT,5,1);


gamma= gamma_GP_from_SP_pt(sal_new,theta,pres_new,lon,lat);



%% THE FIGURE

ccolor=[1 1 1];

figure('Renderer', 'painters', 'Position', [150 150 750 550])
subplot(4,1,1)
hold on

contourf(llatredx,zxx,alongf,[-.9:.1:.9],'linestyle','none'); shading interp

caxis([-.9 1]);cmocean('balance',18)
axis ij
hold on


% contour(llatredx,zxx,saam,[33.4:0.1:34.1],'k'); shading interp
% cmocean('balance');colorbar;caxis([-1 1])

for i=1:length(shallow_deep_leg)
xline(lat(shallow_deep_leg(i)),'linestyle','--')
end
xlim([lat(shallow_deep_leg(1)) lat(shallow_deep_leg(end))])

colorbar

[C,hContour] = contour(lat(shallow_deep_leg),1:4219,gamma(:,shallow_deep_leg),[24 24.5 25 25.5 26 26.35 26.75 27.1] ,'color',ccolor, 'showtext','on','linewidth',.75,'Labelspacing',500);
set(hContour,'showtext','on');clabel(C,hContour,'FontSize',11,'FontWeight','bold','Color',ccolor);
contour(llatredx,zxx,saa,[34.05 34.05],'k','linewidth',1);% shading interp


ylim([0 500])
xlim([lat(shallow_deep_leg(1))-0.005 lat(shallow_deep_leg(end))+0.005])

xticks([-37.2:0.2:-35.2]);xticklabels({});


aa = colorbar;
aa.Label.String = 'Along (m s^-^1 )';
box on;%grid on;

ylabel('(dbar)')


subplot(4,1,2)
hold on

contourf(llatredx,zxx,acrosf,[-.9:.1:.9],'linestyle','none'); shading interp

%pcolor(llatredx,zxx,acrosf); shading interp

[C,hContour] = contour(lat(shallow_deep_leg),1:4219,gamma(:,shallow_deep_leg),[24 24.5 25 25.5 26 26.35 26.75 27.1] ,'color',ccolor, 'showtext','on','linewidth',.75,'Labelspacing',500);
set(hContour,'showtext','on');clabel(C,hContour,'FontSize',11,'FontWeight','bold','Color',ccolor);
contour(llatredx,zxx,saa,[34.05 34.05],'k','linewidth',1);% shading interp

%contourf(llatredx,zxx,acrosf,[-.9:0.1:.9],'linestyle','none'); shading interp

caxis([-.9 .9]);cmocean('balance',19)
axis ij
hold on
xlim([lat(shallow_deep_leg(1)) lat(shallow_deep_leg(end))])
%contour(llatredx,zxx,saam,[33.4:0.1:34.1],'k'); shading interp
%cmocean('balance');colorbar;caxis([-1 1])
% 
% 

for i=1:length(shallow_deep_leg)
xline(lat(shallow_deep_leg(i)),'linestyle','--')
end


xticks([-37.2:0.2:-35.2]);xticklabels({});


aa = colorbar;
aa.Label.String = 'Across (m s^-^1 )';


ylim([0 500])
xlim([lat(shallow_deep_leg(1))-0.005 lat(shallow_deep_leg(end))+0.005])
box on;%grid on;

ylabel('(dbar)')

subplot(4,1,3)
hold on
SAa=movmean(SA,5,1);
pcolor(lat(shallow_deep_leg),pres_new(:,shallow_deep_leg),SAa(:,shallow_deep_leg)); shading interp

%contour(lat(shallow_deep_leg),pres_new(:,shallow_deep_leg),SA(:,shallow_deep_leg),[34 34],'k');


yy=1;pp=1;ee=2;
scatter(tsg.latitud(27870-yy:ee:30761+pp),tsg.latitud(27870-yy:ee:30761+pp)*0,25,SAtsg((27870-yy:ee:30761+pp)),'filled');


%contour(lat(shallow_deep_leg),1:4219,gamma(:,shallow_deep_leg),'color',ccolor,'showtext','on');% shading interp

[C,hContour] = contour(lat(shallow_deep_leg),1:4219,gamma(:,shallow_deep_leg),[24 24.5 25 25.5 26 26.35 26.75 27.1] ,'color',ccolor, 'showtext','on','linewidth',.75,'Labelspacing',500);
set(hContour,'showtext','on');clabel(C,hContour,'FontSize',11,'FontWeight','bold','Color',ccolor);
contour(llatredx,zxx,saa,[34.05 34.05],'k','linewidth',1);% shading interp


for i=1:length(shallow_deep_leg)
xline(lat(shallow_deep_leg(i)),'linestyle','--','color','k')
end

caxis([33.5 37.05]);

cmocean('haline',10)
colorbar

aa = colorbar;
aa.Label.String = 'S_A (g kg^-^1 )';


ylim([0 500])
xlim([lat(shallow_deep_leg(1))-0.005 lat(shallow_deep_leg(end))+0.005])

axis ij
hold on
box on;%grid on;
ylabel('(dbar)')


xticks([-37.2:0.2:-35.2]);xticklabels({});



subplot(4,1,4)

%subplot(5,2,9:10)

hold on
CTa=movmean(CT,5,1);CTa=movmean(CTa,3,1);
pcolor(lat(shallow_deep_leg),pres_new(:,shallow_deep_leg),CTa(:,shallow_deep_leg)); shading interp
colorbar
%pcolor(lat(shallow_deep_leg),pres_new(:,shallow_deep_leg),CTa(:,shallow_deep_leg)); shading interp
hold on

[C,hContour] = contour(lat(shallow_deep_leg),1:4219,gamma(:,shallow_deep_leg),[24 24.5 25 25.5 26 26.35 26.75 27.1] ,'color',ccolor, 'showtext','on','linewidth',.75,'Labelspacing',500);
set(hContour,'showtext','on');clabel(C,hContour,'FontSize',11,'FontWeight','bold','Color',ccolor);
contour(llatredx,zxx,saa,[34.05 34.05],'k','linewidth',1);% shading interp


axis ij
hold on

% yy=1;pp=1;ee=20;
scatter(tsg.latitud(27870-yy:ee:30761+pp),tsg.latitud(27870-yy:ee:30761+pp)*0,25,CTtsg((27870-yy:ee:30761+pp)),'filled');


for i=1:length(shallow_deep_leg)
xline(lat(shallow_deep_leg(i)),'linestyle','--','color','k')
end

caxis([4 24])
cmocean('thermal',10);



% aa = colorbar;
% aa.Label.String = '\Theta (ºC)';


%xlabel('Latitude')



ylabel('(dbar)')


set(gcf,'color','w');box on;%grid on;%set(findall(gcf,'-property','FontSize'),'FontSize',14)

ylim([0 500])
xlim([lat(shallow_deep_leg(1))-0.005 lat(shallow_deep_leg(end))+0.005])

xticks([-37.2:0.2:-35.2]);

xticklabels({'37.2ºS','37ºS','36.8ºS','36.6ºS','36.4ºS','36.2ºS','36ºS','35.8ºS','35.6ºS','35.4ºS','35.2ºS'});



% 
%  savefig('red_transect')
% 
%  print(gcf,'-painters','-depsc2','-r600','red_transect')
% 





%% RED TRANSECT

load('/Users/gaston/Documents/Phd/daily_eddies/adt_20160503.mat')

 WL_plot=-60;EL_plot=-44;SL_plot=-41;NL_plot=-29;
m_proj('mercator','lon',[-54.5 -51.2], 'lat',[-37.5 -34.8]);

date_of_interest=datenum('2016-05-03')%+(i)-1;

% m_proj('mercator','lon',[-56 -50], 'lat',[-38.3 -33.7]);hold on
% load (([input_dir,'/adt_',datestr(date_of_interest,'yyyymmdd'),'.mat']))
% ii=find(time_sst==date_of_interest);
% m_pcolor(lon_sst,lat_sst,sst(:,:,ii)');shading interp
% % [CS,CH]=m_etopo2('contour',[ -5000 -4000 -3000 -2000 -1000 -200],'color',[.65 .65 .65]);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% [CS,CH]=m_etopo2('contour',[-1000 -200],'color','g','linewidth',1);%caxis([-5500 7000]);colormap(flip(gray(7)));  
% caxis([10 25])


figure%('Renderer', 'painters', 'Position', [150 150 500 400])
%m_proj('mercator','lon',[-56 -50], 'lat',[-38 -33.8]);
m_proj('mercator','lon',[-54.5 -51.2], 'lat',[-37.5 -34.8]);


[CS,CH]=m_etopo2('contour',[-2500 -500],'color','g','linewidth',1);
hold on
m_plot([S1.X],[S1.Y], 'k');

yy=000;pp=000

%m_scatter(tsg.longitud(27870-yy:30761+pp),tsg.latitud(27870-yy:30761+pp),25,SAtsg((27870-yy:30761+pp)));
caxis([33.7 37])
cmocean('haline',10)

% 
jj=100;
jj=0;

hold on; m_plot(lonadcp(8497-jj:9584+jj),latadcp(8497-jj:9584+jj),'.k','markersize',4)

% ax = gca;
% salinity_rdlp = colormap(ax);
% save('salinity_rdlp','salinity_rdlp')

uu=nmean(u(1:3,:));vv=nmean(v(1:3,:));

%jj=1500;t=20;kq=0;
 jj=1;t=20;kq=1;



m_contour(X,Y,ADT'*100,[-40:5:100],'color',[.5 .5 .5]);
hold on

m_contour(X,Y,ADT'*100,[25 25],'y','linewidth',1);
m_contour(X,Y,ADT'*100,[0 0],'c','linewidth',1);
m_contour(X,Y,ADT'*100,[50 50],'m','linewidth',1);





%m_quiver(lonxx(1:t:end,1:t:end),latxx(1:t:end,1:t:end),uaug(1:t:end,1:t:end)',vaug(1:t:end,1:t:end)','k','autoscalefactor',2);

[lonx,latx]=meshgrid(X,Y);
t=1;
%m_quiver(lonx(1:t:end,1:t:end),latx(1:t:end,1:t:end),U(1:t:end,1:t:end)',V(1:t:end,1:t:end)','color',[.35 .35 .35],'autoscalefactor',1.1);

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





load('termals')

for i=1:length(shallow_deep_leg)
m_plot(lon(shallow_deep_leg(i)),lat(shallow_deep_leg(i)),'+','Markersize',25,'linewidth',2,'color',termals(i,:))
hold on
end


t=12;
m_quiver(lonadcp(8497-jj*kq:t:9584+jj),latadcp(8497-jj*kq:t:9584+jj),uu(8497-jj*kq:t:9584+jj)',vv(8497-jj*kq:t:9584+jj)','k','autoscalefactor',0.8);
hold on


m_plot([S1.X],[S1.Y], 'k');
m_grid('linestyle','none')

set(gcf,'color','w');


set(findall(gcf,'-property','FontSize'),'FontSize',12)



print(gcf,'-painters','-depsc2','-r600','red_transect_location')


savefig('red_transect_location')




load('ancap_grillado.mat');nlon=[-56:0.25:-50];nlat=[-38:0.25:-33.75];
% [X,Y,Z]=meshgrid([-56:0.25:-50],[-38:0.25:-33.75],1:4219);

%figure('Renderer', 'painters', 'Position', [150 150 500 400])


figure
m_proj('mercator','lon',[-54.5 -51.2], 'lat',[-37.5 -34.8]);

%m_proj('mercator','lon',[-56 -50], 'lat',[-38 -33.8]);
hold on
m_pcolor(nlon,nlat,squeeze(nmean(ssg(:,:,40:120),3)));shading interp;
shading interp;cmocean('haline',8);caxis([33 37])

m_contour(nlon,nlat,squeeze(nmean(ttg(:,:,40:120),3)),[10:2:20],'k','showtext','on');
[CS,CH]=m_etopo2('contour',[-2000 -500],'color','g','linewidth',1);
hold on
%m_plot(lon,lat,'+k')

m_plot([S1.X],[S1.Y], 'k');
hold on

colorbar
m_grid('linestyle','none');


aa = colorbar;
aa.Label.String = 'S_A (g kg^-^1 )';

set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gcf,'color','w');

savefig('red_transect_sal40120')
 print(gcf,'-painters','-depsc2','-r600','red_transect_sal40120')

figure;
theta_sdiag_background
for i=1:length(shallow_deep_leg)
%plot(SA(:,shallow_deep_leg(i)),-pres_new(:,shallow_deep_leg(i)),'linewidth',1,'color',termals(i,:))
hold on
% plot(SAa(:,shallow_deep_leg(i)),CTa(:,shallow_deep_leg(i)),'linewidth',1,'color',termals(i,:))

   SAa=movmean(SA,5,1);SAa=movmean(SAa,3,1);
  CTa=movmean(CT,5,1);CTa=movmean(CTa,3,1);
% 
% 
 plot(SAa(:,shallow_deep_leg(i)),CTa(:,shallow_deep_leg(i)),'color',termals(i,:),'linewidth',1.5)
% plot(SAa(:,shallow_deep_leg(i)),CTa(:,shallow_deep_leg(i)),'linewidth',1,'color',termals(i,:))
% 
% hold on
end

xlim([33.5 37]);
ylim([2 24]);


set(gcf,'color','w');set(findall(gcf,'-property','FontSize'),'FontSize',12)

grid on

set(findall(gcf,'-property','FontSize'),'FontSize',12)


print(gcf,'-painters','-depsc2','-r600','red_transect_TS')

savefig('red_transect_TS')







figure;
%theta_sdiag_background
for i=1:length(shallow_deep_leg)
plot(SA(:,shallow_deep_leg(i)),-pres_new(:,shallow_deep_leg(i)),'linewidth',1,'color',termals(i,:))
hold on
% plot(SAa(:,shallow_deep_leg(i)),CTa(:,shallow_deep_leg(i)),'linewidth',1,'color',termals(i,:))

%   SAa=movmean(SA,5,1);SAa=movmean(SAa,3,1);
%   CTa=movmean(CT,5,1);CTa=movmean(CTa,3,1);
% 
% 
% plot(SAa(:,shallow_deep_leg(i)),CTa(:,shallow_deep_leg(i)),'.','color',termals(i,:))
% plot(SAa(:,shallow_deep_leg(i)),CTa(:,shallow_deep_leg(i)),'linewidth',1,'color',termals(i,:))
% 
% hold on
end

%xlim([33.5 37]);
%ylim([2 24]);



%%
















load('/Users/gaston/Documents/Phd/daily_eddies/adt_20160503.mat')

% valor=colocalize_cruise(lonadcp,latadcp,timeadcp,v_sup',lon,lat,time);


for i=1:length(lonadcp)
valorv(i)=colocalize_cruise(X,Y,1,V,lonadcp(i),latadcp(i),1);
end


for i=1:length(lonadcp)
valoru(i)=colocalize_cruise(X,Y,1,U,lonadcp(i),latadcp(i),1);
end

vcoloc=valorv(8497:9584);


ucoloc=valoru(8497:9584);







figure;

subplot(2,1,1)
plot(vcoloc)
hold on
plot(v_sup(8497:9584))

plot(v(1,8497:9584))

plot(v(7,8497:9584))

plot(v(30,8497:9584))




i=-40;
acroscoloc= (ucoloc.*cosd(i))+(vcoloc*sind(i));
alongcoloc= -(ucoloc*sind(i))+(vcoloc*cosd(i));



for i=1:31

[ro(i),pval(i)]=corr(vcoloc',v(i,8497:9584)','rows','complete')

[rmsee(i)]=rmse(vcoloc',v(i,8497:9584)')

end



for i=1:31

[rou(i),pvalu(i)]=corr(ucoloc',u(i,8497:9584)','rows','complete')

[rmseeu(i)]=rmse(ucoloc',u(i,8497:9584)')

end




for i=1:31

[roa(i),pvala(i)]=corr(acroscoloc',acros(i,8497:9584)','rows','complete')

[rmseea(i)]=rmse(acroscoloc',acros(i,8497:9584)')

end



for i=1:31

[roal(i),pvalal(i)]=corr(alongcoloc',along(i,8497:9584)','rows','complete')

[rmseeal(i)]=rmse(alongcoloc',along(i,8497:9584)')

end



figure;
subplot(2,1,1)
plot(rmseea)
hold on
plot(rmseeal)
subplot(2,1,2)
plot(roa)
hold on
plot(roal)



for i=1:31

[roa(i),pvala(i)]=corr(acroscoloc',acros(i,8497:9584)','rows','complete')

[rmseea(i)]=rmse(acroscoloc',acros(i,8497:9584)')

end





for i=1:31

[rouc(i),pvaluc(i)]=corr(ucoloc',u(i,8497:9584)','rows','complete')

[rmseeu(i)]=rmse(ucoloc',u(i,8497:9584)')

end




for i=1:31

[ro(i),pval(i)]=corr(vcoloc',v(i,8497:9584)','rows','complete')

[rmsee(i)]=rmse(vcoloc',v(i,8497:9584)')

end



for i=1:31

[rou(i),pvalu(i)]=corr(ucoloc',u(i,8497:9584)','rows','complete')

[rmseeu(i)]=rmse(ucoloc',u(i,8497:9584)')

end







distancia=sw_dist([lonadcp(8497) lonadcp(9584)],[latadcp(8497) latadcp(9584)],'km')


xc=linspace(1,280,1088);

figure('Renderer', 'painters', 'Position', [150 150 750 600])

subplot(4,1,1)
plot(xc,along(1,8497:9584))
hold on
plot(xc,along(6,8497:9584))
plot(xc,along(18,8497:9584))
plot(xc,alongcoloc)

axis tight

legend('28m','109m','300m','Altimetry')
ylabel('Along vel (m s ^-^1)')

xlabel('Distance (km)')
subplot(4,1,2)
plot(xc,acros(1,8497:9584))
hold on
plot(xc,acros(5,8497:9584))
plot(xc,acros(10,8497:9584))
plot(xc,acroscoloc)

axis tight

ylabel('Across vel (m s ^-^1)')
legend('28m','109m','300m','Altimetry','location','best')

xlabel('Distance (km)')

subplot(4,1,3)
plot(z(1:31),roa)
hold on
plot(z(1:31),roal)
axis tight
ylabel('CORR'); xlabel('Depth level (m)')
legend('Across','Along','location','best')

xlim([0.2 1])


subplot(4,1,4)
plot(z(1:31),rmseea)
hold on
plot(z(1:31),rmseeal)
legend('Across','Along','location','best')
axis tight
ylabel('RMSE'); xlabel('Depth level (m)')


set(gcf,'color','w');




savefig('corr_altimetry')


[p,q]=corr(vcoloc',v_sup(8497:9584)')

[pu,qu]=corr(ucoloc',u_sup(8497:9584)')





% 
% ind=find(saam<34)
% 
% 
% 
% sw_dist([llatredx(1) llatredx(2)],[-55 -55],'km')
% 
% 
% 
% 
% t7=griddata(aa,ppres,temp_new(:,71:82),d7,ppres);
% 
% 
% figure
% 
% subplot(tt,1,1)
% 
% pcolor(along(:,8497:9584));shading interp
% 
% 
% cmocean('balance');colorbar;caxis([-1 1])
% 
% 
% subplot(tt,1,2)
% 
% 
% pcolor(lat(shallow_deep_leg),pres_new(:,1),SA(:,shallow_deep_leg));shading interp
% 
% 
% shallow_deep_leg(1)
% 
% shallow_deep_leg(end)
% 
% 
% griddata(lat(shallow_deep_leg),lon(shallow_deep_leg)
% 
% 
% 
% 
% 
% axis ij
% 
% ylim([0 500])
% caxis([33.5 37])
% 
% cmocean('haline')
% colorbar
% % 
% % jq=250;
% % m_scatter(tsg.longitud(28015-jq:30758+jq),tsg.latitud(28015-jq:30758+jq),5,tsg.salinidad(28015-jq:30758+jq),'filled')
% % 
% % caxis([33.5 36])
% % 
% % figure
% % m_proj('mercator','lon',[-56 -50], 'lat',[-38 -33.8]);
% % jq=550;kk=2;
% % [CS,CH]=m_etopo2('contour',[-2000 -500],'color','g','linewidth',1);
% % 
% % %m_scatter(tsg.longitud(28015-jq*kk:30758+jq),tsg.latitud(28015-jq*kk:30758+jq),5,tsg.temperatura(28015-jq*kk:30758+jq),'filled')
% % %caxis([10 25])
% % 
% % m_scatter(tsg.longitud(28015-jq*kk:30758+jq),tsg.latitud(28015-jq*kk:30758+jq),5,tsg.salinidad(28015-jq*kk:30758+jq),'filled')
% % caxis([33.6 37])
% % 
% % m_plot([S1.X],[S1.Y], 'k');
% % 
% % m_grid
% 
% 
% 
% 
% 
% colormap(jet)
% 
% 
% figure
% jq=250;kk=5;
% jq=0;kk=0;
% 
% subplot(2,1,1)
% plot(tsg.temperatura(28015-jq*kk:30758+jq)); hold on
% yyaxis right
% plot(tsg.salinidad(28015-jq*kk:30758+jq)); hold on
% yyaxis right
% subplot(2,1,2)
% plot(tsg.fluor(28015-jq*kk:30758+jq)); hold on
% 
% 
% 
% figure
% jq=250;kk=5;
% 
% jq=250;kk=5;
% scatter(tsg.longitud(28015-jq*kk:28015),tsg.latitud(28015-jq*kk:28015),5,tsg.temperatura(28015-jq*kk:28015),'filled')
% 
% figure
% jq=250;kk=5;
% scatter(tsg.longitud(28015-jq*kk:28015),tsg.latitud(28015-jq*kk:28015),5,tsg.salinidad(28015-jq*kk:28015),'filled')
% 
% 
% jq=1;kk=1;
% 
% 
% x1 = (1:length(tsg.longitud(28015-jq*kk:30758)));
% 
% 
% for i=1:length(x1)-1
% 
% distance(i)=sw_dist([tsg.latitud(i) tsg.latitud(i+1)],[tsg.longitud(i) tsg.longitud(i+1)],'km');
% 
% end
% 
% 
% 
% distanc=ncumsum(distance);x1=[distanc distanc(end)];
% 
% 
% x1=[distanc distanc(end) distanc(end)]';
% 
% 
% 
% x2 = x1;
% x3 = x1;
% y1 = tsg.temperatura(28015-jq*kk:30758+jq);
% y2 = SAtsg(28015-jq*kk:30758+jq);
% 
% %y2 = tsg.salinidad(28015-jq*kk:30758+jq);
% y3 = tsg.fluor(28015-jq*kk:30758+jq);
% 
% y4 = [acros_sup_redi ;acros_sup_redi(end); acros_sup_redi(end)];
% y5 = [along_sup_redi; along_sup_redi(end) ;along_sup_redi(end)];
% 
% ylat = tsg.latitud(28015-jq*kk:30758+jq);
% ylon = tsg.longitud(28015-jq*kk:30758+jq);
% 
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
% 
% 
% 
% 
% 
% 
% 
% 
% figure;plot(ylat,y2)
% hold on;yyaxis right;plot(ylat,y1)
% 
% 
% 
% 
% 
% % y4 = u_sup_redi
% % y5 = v_sup_redi
% 
% 
% %x3 = x1
% 
% 
% % y4 = sin(x1);
% % y5 = fliplr(2*x1.^2);
% % y6 = 7*cos(x1);
% % y7 = 7*log(x1+1.2);
% 
% 
% 
% 
% figure
% [ax,hlines,fh] = jzplotys({x1,y1,x1,y2,x1,y3,x1,y4,x1,y5},[1 1 1 1 1],[x1(1) x1(end)],100);
% 
% legend(cat(1,hlines{:}),'Temperature','Salinity','Fluorescene','Across-slope vel.','Along-slope vel.','orientation','horizontal','location',[.2 .05 .01 .01],'box','off')
% 
% ylabel(ax(1),'\Theta (ºC)');
% 
% ylabel(ax(3),'S_A (g kg^-^1)');
% 
% ylabel(ax(5),'Fluorescene (mg m^-^3)');
% 
% ylabel(ax(7),'Across slope vel. (m s^-^1)');
% 
% ylabel(ax(9),'Along slope vel. (m s^-^1)');
% 
% set(gcf,'color','w');set(findall(gcf,'-property','FontSize'),'FontSize',12)
% 
% 
% set(gcf,'color','w');set(findall(gcf,'-property','FontSize'),'FontSize',12.5)
% 
% set(gcf,'color','w');set(findall(gcf,'-property','linewidth'),'linewidth',1.1)
% 
% 
% 
% 
% 
% %savefig('TSG_red_transect')
% %print(gcf,'-painters','-depsc2','-r600','TSG_red_transect')
% 
% 
% 
% %xlabel('Latitude')
% 
% 
% colormap(flipud(brewermap('Spectral')));
% 
% colormap(flipud(brewermap([],'Spectral')));
% 
% colormap(flipud(brewermap(7,'Spectral')));
% 
% 
% hold on
% 
% figure
% 
% m_proj('mercator','lon',[-56 -50], 'lat',[-38 -33.8]);
% 
% % m_proj('mercator','lon',[-54 -51.5], 'lat',[-37.5 -34.9]);
% 
% m_plot(lonadcp,latadcp,'.','markersize',4,'color',[.9 .9 .9]); hold on
% 
% [CS,CH]=m_etopo2('contour',[-2000 -500],'color','g','linewidth',1);
% hold on
% 
% hold on; m_plot(lonadcp(8497:9584),latadcp(8497:9584),'.r','markersize',2)
% 
% 
% % cmocean('thermal',7);ax = gca;
% % termals = colormap(ax);
% load('termals')
% 
% for i=1:length(shallow_deep_leg)
% m_plot(lon(shallow_deep_leg(i)),lat(shallow_deep_leg(i)),'+','Markersize',15,'linewidth',5,'color',termals(i,:))
% hold on
% end
% 
% 
% m_plot([S1.X],[S1.Y], 'k');
% m_grid('linestyle','none')
% 
% set(gcf,'color','w');set(findall(gcf,'-property','FontSize'),'FontSize',12)
% 
% % (gcf,'-painters','-depsc2','-r600','red_transect_location')
% 
% 
% 
% figure;
% %theta_sdiag_background
% for i=1:length(shallow_deep_leg)
% 
% plot(SA(:,shallow_deep_leg(i)),-pres_new(:,shallow_deep_leg(i)),'linewidth',1,'color',termals(i,:))
% hold on
% 
% % plot(SAa(:,shallow_deep_leg(i)),CTa(:,shallow_deep_leg(i)),'linewidth',1,'color',termals(i,:))
% 
% %   SAa=movmean(SA,5,1);SAa=movmean(SAa,3,1);
% %   CTa=movmean(CT,5,1);CTa=movmean(CTa,3,1);
% % 
% % 
% % plot(SAa(:,shallow_deep_leg(i)),CTa(:,shallow_deep_leg(i)),'.','color',termals(i,:))
% % plot(SAa(:,shallow_deep_leg(i)),CTa(:,shallow_deep_leg(i)),'linewidth',1,'color',termals(i,:))
% % 
% % hold on
% end
% 
% xlim([33.5 37]);
% ylim([2 24]);
% set(gcf,'color','w');set(findall(gcf,'-property','FontSize'),'FontSize',15)
% 
% 
%  %print(gcf,'-painters','-depsc2','-r600','red_transect_TS')
% 
% 
% 
% % 
% % ylabel(ax(1),'Y1');
% % ylabel(ax(3),'Y2');
% % ylabel(ax(5),'Y3');
% % ylabel(ax(7),'Y4');
% 
% 
% 
% 
% 
% figure
% m_scatter(tsg.longitud,tsg.latitud,5,tsg.salinidad,'filled')
% 
% %"landscape" or "portrait"
% m_vec(1,lonadcp(8497-jj:t:9584+jj),latadcp(8497-jj:t:9584+jj),uu(8497-jj:t:9584+jj),vv(8497-jj:t:9584+jj),'k');
% 
% 
% [minVal,idx]=find_close_value(tsg.longitud,-51.6314)
% 
% 
% [dmin,imin]=pdist2([tsg.longitud tsg.latitud],[-53.7876 -37.1887],'euclidean','smallest',1);
% 
% [dmin,imin]=pdist2([tsg.longitud tsg.latitud],[-51.6451 -35.0454],'euclidean','smallest',1);
% 
% hold on; m_plot(tsg.longitud(28015:30758),tsg.latitud(28015:30758),'.b','markersize',6)
% 
% 
% 
% 
% [minVal,idx]=find_close_value(tsg.longitud,-53.7876)
% 
% [minVal,idx]=find_close_value(tsg.latitud,-37.1887)
% 
% 
% 
% shallow_deep_leg=([9 15 29 40 52 65 76]);
% 
% for k=1:length(shallow_deep_leg)
% m_plot(lon(shallow_deep_leg(k)),lat(shallow_deep_leg(k)),'+k')
% m_text(lon(shallow_deep_leg(k)),lat(shallow_deep_leg(k)),num2str(stations(shallow_deep_leg(k))),'fontsize',13,'color','k')
% end
% 
% m_grid('linestyle','none')
% 
% set(gcf,'color','w');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% confluence_leg=[9 15 29 40 52 65 76]
% 
% figure;
% subplot(4,3,1)
% 
% pcolor(lat16(7791:8873),z16,vv_chuy(:,7791:8873));shading interp;% axis ij
% 
% 
% ylim([-500 0])
% 
% cmocean('balance');colorbar;caxis([-1 1])
% 
% %xlim([-37.2 -35.2])
% 
% subplot(4,3,2)
% 
% pcolor(lat16(7791:8873),z16,uu_chuy(:,7791:8873));shading interp;% axis ij
% 
% ylim([-500 0])
% 
% cmocean('balance');colorbar;caxis([-1 1])
% %xlim([-37.2 -35.2])
% 
% 
% subplot(4,3,3)
% 
% pcolor(lat3(:,shallowest_deep_leg),z_new(:,shallow_deep_leg),sal_new(:,shallow_deep_leg)); shading interp; axis ij;
% 
% ylim([0 500]); colorbar
% 
% 
% %xlim([-37.2 -35.2])
% 
% 
% subplot(4,1,4)
% 
% %contourf(lat3(:,shallow_deep_leg),z_new(:,shallow_deep_leg),temp_new(:,shallow_deep_leg)); %shading interp;
% %axis ij;
% 
% pcolor(lat3(:,shallow_deep_leg),z_new(:,shallow_deep_leg),temp_new(:,shallow_deep_leg));
%  
% shading interp; axis ij;
% 
% ylim([0 800]); colorbar
% 
% %xlim([-37.2 -35.2])
% 
% 
% subplot(4,1,3)
% 
% pcolor(lat3(:,shallowest_deep_leg),z_new(:,mid_deep_leg),sal_new(:,mid_deep_leg)); shading interp; axis ij;
% 
% ylim([0 500]); colorbar
% %xlim([-37.2 -35.2])
% 
% 
% subplot(4,1,4)
% 
% %contourf(lat3(:,shallow_deep_leg),z_new(:,shallow_deep_leg),temp_new(:,shallow_deep_leg)); %shading interp;
% %axis ij;
% 
% pcolor(lat3(:,shallow_deep_leg),z_new(:,mid_deep_leg),temp_new(:,mid_deep_leg));
%  
% shading interp; axis ij;
% 
% ylim([0 800]); colorbar
% 
% 
% 
% 
% 
% %%
% 
% 
% %% Plot ADT and FSL for Gaston's paper
% clc
% clearvars
% 
% addpath(genpath('/Users/rlaxenaire/OneDrive - Florida State University/On_Work/Missions/21_07_Tara_Microbiomes/Scripts'))
% 
% Main_Dir='/Users/rlaxenaire/OneDrive - Florida State University/On_Work/Papiers/Co_Auth/2021_Manta';
% input_dir_FSLE='/Users/rlaxenaire/OneDrive - Florida State University/Boulot/Data/FSLE/delayed-time';
% input_dir_ADT='/Users/rlaxenaire/OneDrive - Florida State University/Boulot/Data/AVISO/Global/Delayed_Time/DT_2018';
% 
% 
% WL=-60;
% EL=-45;
% SL=-43;
% NL=-32;
% 
% WL_zoom=-58;
% EL_zoom=-50;
% SL_zoom=-38;
% NL_zoom=-33;
% 
% plot_visible='off';
% 
% % Format Figures
% res_fig='-r200';
% format_fig='-dpng';
% 
% %% Bathy 
% dx_topo=.25;dy_topo=dx_topo;
% SL_topo=SL-dx_topo;NL_topo=NL+dx_topo;WL_topo=WL-dy_topo;EL_topo=EL+dy_topo;
% 
% [Bathy_full, ~] = etopo(fullfile('/Users/rlaxenaire/OneDrive - Florida State University/Boulot/Data','Etopo','etopo1_ice_c.flt'),1 ...
%     ,[SL_topo NL_topo] ...
%     ,[WL_topo EL_topo]);
% 
% Lon_topo=linspace(WL_topo,EL_topo,size(Bathy_full,2))';
% Lat_topo=linspace(SL_topo,NL_topo,size(Bathy_full,1))';
% 
% %% FSLE
% 
% cd /Users/gaston/Downloads
% ncdisp('dt_global_allsat_madt_fsle_20160505_20180704.nc')
% 
% 
% fsle=ncread('dt_global_allsat_madt_fsle_20160505_20180704.nc','fsle_max');
% 
% lonf=ncread('dt_global_allsat_madt_fsle_20160505_20180704.nc','lon');
% lonf=lonf-360;
% 
% latf=ncread('dt_global_allsat_madt_fsle_20160505_20180704.nc','lat');
% 
% load '/Users/gaston/Documents/adts_eddies/adt_20160505.mat'
% 
% % all_files=dir(fullfile(input_dir_FSLE,'dt_global_allsat_madt_fsle_*.nc'));
% % 
% % Dir_fig_FSLE=fullfile(Main_Dir,'Figures','FSLE');if exist(Dir_fig_FSLE,'dir')==0;mkdir(Dir_fig_FSLE);end;
% % 
% % for loop_all_files=1:size(all_files,1)
% %     
% %     filename_fsle=fullfile(all_files(loop_all_files).folder,all_files(loop_all_files).name);
% %         time_FSLE=double(ncread(filename_fsle,'time')+datenum(1950,01,01,00,00,00));
% %     filename_fig_fsle=fullfile(Dir_fig_FSLE,['FSLE_',datestr(time_FSLE,'yyyymmdd'),'.',format_fig(3:end)]);
% %             filename_fig_fsle_Zoom=fullfile(Dir_fig_FSLE,['Zoom_FSLE_',datestr(time_FSLE,'yyyymmdd'),'.',format_fig(3:end)]);
% % 
% %     if exist(filename_fig_fsle,'file');continue;end;
% %     
% %     all_ADT=dir(fullfile(input_dir_ADT,['dt_global_allsat_phy_l4_',datestr(time_FSLE,'yyyymmdd'),'_*.nc']));
% %     filename_adt=fullfile(all_ADT(1).folder,all_ADT(1).name);
% % 
% %     
% %     
% %     if ~exist('lon_FSLE','var')
% %         lon_FSLE=ncread(filename_fsle,'lon');
% %         lon_FSLE(lon_FSLE>180)=lon_FSLE(lon_FSLE>180)-360;
% %         
% %         id_out_lon=find(lon_FSLE<WL-1 | lon_FSLE>EL+1);
% %         lon_FSLE(id_out_lon)=[];
% %         
% %         lat_FSLE=ncread(filename_fsle,'lat');
% %         id_out_lat=find(lat_FSLE<SL-1 | lat_FSLE>NL+1);
% %         lat_FSLE(id_out_lat)=[];
% %         
% %         RGB=flipud(gray(16));
% %         
% %          lon_ADT=ncread(filename_adt,'longitude');
% %         lon_ADT(lon_ADT>180)=lon_ADT(lon_ADT>180)-360;
% %          lat_ADT=ncread(filename_adt,'latitude');
% %     end
% %     
% %     FSLE=ncread(filename_fsle,'fsle_max');
% %     FSLE(id_out_lon,:)=[];
% %     FSLE(:,id_out_lat)=[];
% %     
% %      ADT=ncread(filename_adt,'adt');
%     
%     % Plot map
% 
% %    clf=figure('visible',plot_visible);
% 
% S1 = shaperead('ZEEU.shp');
%  figure
%    m_proj('lambert','lon',[WL EL], 'lat',[SL  -30]);hold on
%    %  m_proj('lambert','lon',[WL EL], 'lat',[-46  -30]);hold on
% 
% hold on
%     m_pcolor(lonf,latf,fsle');shading interp;
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);    m_grid;
% [CS,CH]=m_etopo2('contour',[-2000 -500],'color','g','linewidth',1);    
%     m_contour(X,Y,ADT',[0 0],'color','c','linewidth',1);
%     m_contour(X,Y,ADT',[.25 .25],'color','y','linewidth',1);
%     m_contour(X,Y,ADT',[.5 .5],'color','m','linewidth',1);
%  m_grid('linestyle','none');
% 
% colormap(gray)
%     caxis([-.6 0])
% hold on
% m_plot([S1.X],[S1.Y], 'k');
% set(gcf,'color','w');
% 
% print(gcf,'-painters','-depsc2','-r600','flse')
% 
% 
% 
% 
% 
% 
% 
% 
% 
%     title(['Aviso''s FSLEs the ',datestr(double(time_FSLE),'dd/mm/yy')]);
%     h=colorbar;ylabel(h,'Backward FSLEs [1/day]')
%     print(clf,format_fig,res_fig,filename_fig_fsle)
% 
% %     close(clf)
% %     
% %     % Plot Zoom
% %     clf=figure('visible',plot_visible);
% %     m_proj('lambert','lon',[WL_zoom EL_zoom], 'lat',[SL_zoom  NL_zoom]);hold on
% %     colormap(RGB);
% %     caxis([-.6 .15])
% %     m_pcm(lon_FSLE,lat_FSLE,FSLE');shading interp;
% %     m_gshhs_l('patch',[.7 .7 .7]);
% %     m_grid;
% %         m_contour(Lon_topo,Lat_topo,Bathy_full,[-2000 -500],'color','g','linewidth',1);
% %     
% %     m_contour(lon_ADT,lat_ADT,ADT',[0 0],'color','c','linewidth',1);
% %     m_contour(lon_ADT,lat_ADT,ADT',[.25 .25],'color','y','linewidth',1);
% %     m_contour(lon_ADT,lat_ADT,ADT',[.5 .5],'color','m','linewidth',1);
% % 
% %     title(['Aviso''s FSLEs the ',datestr(double(time_FSLE),'dd/mm/yy')]);
% %     h=colorbar;ylabel(h,'Backward FSLEs [1/day]')
% %     print(clf,format_fig,res_fig,filename_fig_fsle_Zoom)
% % 
% %     close(clf)
% %     
% % end
% 
% 
