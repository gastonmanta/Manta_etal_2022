
cd /Users/gaston/Documents/Phd/Database_South_Atl/Nc_filess/

adt=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1554459059166_swa_adt.nc','adt');
u=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1554459393439_swa_ugos.nc','ugos');
v=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1554459685393__swa_vgos.nc','vgos');
lonx=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1554459685393__swa_vgos.nc','longitude');lonx=lonx-360;
latx=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1554459685393__swa_vgos.nc','latitude');

time_adt=ncread('dataset-duacs-rep-global-merged-allsat-phy-l4_1554459059166_swa_adt.nc','time');time_adt=double(time_adt+datenum('1950-01-01 00:00:00'))


% load ('/Users/gaston/Documents/Phd/satelite_gamboa/percentage_presence_eddies.mat');cpo=cpo(:,:,end-3650:end);apo=apo(:,:,end-3650:end);
% cpom=nsum(cpo,3).*100./length(cpo);
% apom=nsum(apo,3).*100./length(cpo);



S1 = shaperead('ZEEU.shp');


% ind=find(apom==0);apom(ind)=NaN;
% ind=find(cpom==0);cpom(ind)=NaN;


 [lonxx,latxx]=meshgrid(lonx,latx);%lonxx=lonxx-360;
% 
% um=nmean(u(:,:,end-3650:end),3);vm=nmean(v(:,:,end-3650:end),3);am=nmean(adt(:,:,end-3650:end),3);
% 
% uu=u(:,:,end-3650:end);
% vv=v(:,:,end-3650:end);

adtm=nanmean(adt,3);
um=nanmean(u,3);
vm=nanmean(v,3);
ke=(um.^2+vm.^2)./2;

%adtrms=rms(permute(adt,[3 1 2]));adtrms=squeeze(adtrms);

kem=nanmean(ke,3);
ekem=nanmean(((u.^2+v.^2)./2-kem),3);

ekekem=ekem./kem;
ekekem=ekem-kem;


% cpom=nsum(cpo,3).*100./length(cpo);
% apom=nsum(apo,3).*100./length(cpo);

% um=nmean(u(:,:,end-3650:end),3);vm=nmean(v(:,:,end-3650:end),3);

am=nmean(adt(:,:,end-3650:end),3);

t=2;


load ('sst_paper_final.mat')


for i=1:length(sst)
[FX(:,:,i),FY(:,:,i)] = gradient(sst(:,:,i));
end

for i=1:length(sst)
grads(:,:,i)=sqrt(FX(:,:,i).^2+FY(:,:,i).^2);
end


gg=nmean(grads(:,:,end-3650:end),3);




[FXm,FYm] = gradient(nmean(sst,3));

gradsm=sqrt(FXm.^2+FYm.^2);


gradsm=nmean(grads,3);


% %for ADT
% 
% for i=1:length(adt)
% [FXa(:,:,i),FYa(:,:,i)] = gradient(adt(:,:,i));
% end
% 
% for i=1:length(adt)
% gradsa(:,:,i)=sqrt(FXa(:,:,i).^2+FYa(:,:,i).^2);
% end
% 
% 
% 
% 
% [FXm,FYm] = gradient(nmean(sst,3));
% 
% gradsm=sqrt(FXm.^2+FYm.^2);




%%

t=2;
figure%('Renderer', 'painters', 'Position', [100 100 500 550])
%m_proj('lambert','long',[-60 -44], 'lat',[-45 -30]); hold on
m_proj('lambert','long',[-60.125 -44.625], 'lat',[-44.12 -31.5]); hold on

%m_proj('lambert','long',[-60.125 -44.625], 'lat',[-46 -30]); hold on


m_contourf(lonx,latx,ekekem',[-.2025:0.025: .2],'linestyle','none');

%[CS,CH]=m_etopo2('contour',[-2500 -1000 -200],'color','g','linewidth',1);
[CS,CH]=m_etopo2('contour',[-2500 -500],'color','g','linewidth',1);

m_quiver(lonxx(1:t:end,1:t:end),latxx(1:t:end,1:t:end),um(1:t:end,1:t:end)',vm(1:t:end,1:t:end)','k','autoscalefactor',2);


hold on
m_contour(lonx,latx,am',[0.5 .5],'m','linewidth',1);
m_contour(lonx,latx,am',[0.0 .0],'c','linewidth',1);
m_contour(lonx,latx,am',[0.25 0.25],'y','linewidth',1);
m_plot([S1.X],[S1.Y], 'k');


colormap(rednblue(17));caxis([-.2 .2])


[ax,h]=m_contfbar([.22 .42],.855,[-.2:0.025: .2],[-.2:0.025: .2]);
title(ax, ' EKE - MKE (m^2 s^-^2)','Fontweight','normal')


m_grid('linestyle','none')
m_gshhs_f('patch',[0.81818 0.77647 0.70909]);
set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',13)


x1=-58.5;x2=-54;y1=-34;y2=-31.5;
xb = [x1, x2, x2, x1, x1];yb = [y1, y1, y2, y2, y1];
m_patch(xb, yb,[0.81818 0.77647 0.70909], 'Linestyle','none') 

m_text(-60,-31.9,'(a)','FontSize',18)


% x1=-56;x2=-50;y1=-38;y2=-33.5;
% xb = [x1, x2, x2, x1, x1];yb = [y1, y1, y2, y2, y1];
%  m_plot(xb, yb, '-', 'LineWidth', 1.5,'color', [1.0000    0.5         0]);


cd /Users/gaston/Desktop
print(gcf,'-painters','-depsc2','-r600','EKE_MKE')

% print(gcf,'-painters','-depsc2','-r600','figure1a_gamboa_may')




%%

t=2;
figure%('Renderer', 'painters', 'Position', [100 100 500 550])
%m_proj('lambert','long',[-60 -44], 'lat',[-45 -30]); hold on
m_proj('lambert','long',[-60.125 -44.625], 'lat',[-44.12 -31.5]); hold on

m_contourf(lon_sst,lat_sst,gradsm',[0.2:0.1:1.4],'linestyle','none');

%[CS,CH]=m_etopo2('contour',[-2500 -1000 -200],'color','g','linewidth',1);
[CS,CH]=m_etopo2('contour',[-2500 -500],'color','g','linewidth',1);

m_quiver(lonxx(1:t:end,1:t:end),latxx(1:t:end,1:t:end),um(1:t:end,1:t:end)',vm(1:t:end,1:t:end)','k','autoscalefactor',2);
hold on


m_contour(lonx,latx,am',[0.5 .5],'m','linewidth',1);
m_contour(lonx,latx,am',[0.0 .0],'c','linewidth',1);
m_contour(lonx,latx,am',[0.25 0.25],'y','linewidth',1);
m_plot([S1.X],[S1.Y], 'k');


 cmocean('amp',11);caxis([0.2 1.4])


[ax,h]=m_contfbar([.22 .42],.855,[0.2:0.1:1.4],[0.2:0.1:1.41]);

title(ax,'\nablaSST (ºC 25km^-^1) ','Fontweight','normal')


m_grid('linestyle','none')
m_gshhs_f('patch',[0.81818 0.77647 0.70909]);
set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',13)

% x1=-58.5;x2=-54;y1=-34;y2=-31.5;
% xb = [x1, x2, x2, x1, x1];yb = [y1, y1, y2, y2, y1];
m_patch(xb, yb,[0.81818 0.77647 0.70909], 'Linestyle','none') 

m_text(-60,-31.9,'(b)','FontSize',18);



cd /Users/gaston/Desktop
print(gcf,'-painters','-depsc2','-r600','grads')


%  m_plot(xb, yb, '-', 'LineWidth', 1.5,'color', [1.0000    0.5         0]);





%%

%% AEs
% 
% figure%('Renderer', 'painters', 'Position', [100 100 500 550])
% %m_proj('lambert','long',[-60 -44], 'lat',[-45 -30]); hold on
% m_proj('lambert','long',[-60.125 -44.625], 'lat',[-44.12 -31.5]); hold on
% m_contourf(X,Y,(apom)',[0:5:60],'linestyle','none');
% %[CS,CH]=m_etopo2('contour',[-1000 -200],'color','g','linewidth',1);
% [CS,CH]=m_etopo2('contour',[-2000 -500],'color','g','linewidth',1);
% 
% caxis([0 60])
% m_quiver(lonxx(1:t:end,1:t:end),latxx(1:t:end,1:t:end),um(1:t:end,1:t:end)',vm(1:t:end,1:t:end)','k','autoscalefactor',2);
% hold on
% m_contour(lonx,latx,am',[0.5 .5],'m','linewidth',1);
% m_contour(lonx,latx,am',[0.0 .0],'c','linewidth',1);
% m_contour(lonx,latx,am',[0.25 0.25],'y','linewidth',1);
% m_plot([S1.X],[S1.Y], 'k');
% 
% cmocean('amp',11)
% 
% m_grid('linestyle','none')
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% set(gcf,'color','w');
% set(findall(gcf,'-property','FontSize'),'FontSize',13)
% 
% 
% x1=-59;x2=-53.7;y1=-33.3;y2=-31.55;
% xb = [x1, x2, x2, x1, x1];yb = [y1, y1, y2, y2, y1];
% m_patch(xb, yb,[0.81818 0.77647 0.70909], 'Linestyle','none') 
% 
% 
% %m_plot(xb, yb, '-', 'LineWidth', 1.5,'color', [1.0000    0.5         0]);
% 
% 
% [ax,h]=m_contfbar([.22 .42],.855,[0:5:50],[0:5:60],'FontSize',12);
% title(ax,'Time spent by AEs (%)','Fontweight','normal')
% 
% 
% m_text(-60,-31.9,'(c)','FontSize',18)
% 
% cd /Users/gaston/Desktop
% 
% print(gcf,'-painters','-depsc2','-r600','time_AEs')
% print(gcf,'-dpng','-r600','time_AEs')
% 
% %% CEs
% 
% figure%('Renderer', 'painters', 'Position', [100 100 500 550])
% %m_proj('lambert','long',[-60 -44], 'lat',[-45 -30]); hold on
% m_proj('lambert','long',[-60.125 -44.625], 'lat',[-44.12 -31.5]); hold on
% m_contourf(X,Y,(cpom)',[0:5:60],'linestyle','none');
% [CS,CH]=m_etopo2('contour',[-2000 -500],'color','g','linewidth',1);
% caxis([0 60])
% m_quiver(lonxx(1:t:end,1:t:end),latxx(1:t:end,1:t:end),um(1:t:end,1:t:end)',vm(1:t:end,1:t:end)','k','autoscalefactor',2);
% hold on
% m_contour(lonx,latx,am',[0.5 .5],'m','linewidth',1);
% m_contour(lonx,latx,am',[0.0 .0],'c','linewidth',1);
% m_contour(lonx,latx,am',[0.25 0.25],'y','linewidth',1);
% m_plot([S1.X],[S1.Y], 'k');
% 
% cmocean('amp',11)
% 
% m_grid('linestyle','none')
% m_gshhs_h('patch',[0.81818 0.77647 0.70909]);
% set(gcf,'color','w');
% set(findall(gcf,'-property','FontSize'),'FontSize',13)
% 
% 
% x1=-59;x2=-53.7;y1=-33.3;y2=-31.55;
% xb = [x1, x2, x2, x1, x1];yb = [y1, y1, y2, y2, y1];
% m_patch(xb, yb,[0.81818 0.77647 0.70909], 'Linestyle','none') 
% 
% 
% %m_plot(xb, yb, '-', 'LineWidth', 1.5,'color', [1.0000    0.5         0]);
% 
% 
% [ax,h]=m_contfbar([.22 .42],.855,[0:5:50],[0:5:60],'FontSize',12);
% title(ax,'Time spent by CEs (%)','Fontweight','normal')
% 
% 
% m_text(-60,-31.9,'(d)','FontSize',18)
% 
% cd /Users/gaston/Desktop
% 
% print(gcf,'-painters','-depsc2','-r600','time_CEs')
% print(gcf,'-dpng','-r600','time_CEs')
% 
