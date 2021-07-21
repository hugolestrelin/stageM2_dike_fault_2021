%clear all
%close all

load parametres.mat

year=3600*24*365;

repres='./results/';

fid=fopen([repres,'t.data'],'r');
t=fread(fid,'real*8');
fclose(fid);

fid=fopen([repres,'v.data'],'r');
v=fread(fid,'real*8');
fclose(fid);

fid=fopen([repres,'th.data'],'r');
th=fread(fid,'real*8');
fclose(fid);

fid=fopen([repres,'u.data'],'r');
u=fread(fid,'real*8');
fclose(fid);

fid=fopen([repres,'s.data'],'r');
s=fread(fid,'real*8');
fclose(fid);

fid=fopen([repres,'tau.data'],'r');
tau=fread(fid,'real*8');
fclose(fid);

fid=fopen([repres,'fco.data'],'r');
fco=fread(fid,'real*8');
fclose(fid);

fid=fopen([repres,'dbc.data'],'r');
dbc=fread(fid,'real*8');
fclose(fid);

fid=fopen([repres,'dpc.data'],'r');
dpc=fread(fid,'real*8');
fclose(fid);

fid=fopen([repres,'tdw.data'],'r');
tdw=fread(fid,'real*8');
fclose(fid);

fid=fopen([repres,'sop.data'],'r');
sop=fread(fid,'real*8');
fclose(fid);

fid=fopen([repres,'pop.data'],'r');
pop=fread(fid,'real*8');
fclose(fid);

iend=min([length(t) length(v) length(th) length(u) length(s) length(tau)]);
t=t(1:iend);v=v(1:iend);th=th(1:iend);u=u(1:iend);s=s(1:iend);tau=tau(1:iend);dbc=dbc(1:iend);dpc=dpc(1:iend);

sop=reshape(sop,iend,[]);

figure(1);clf;
subplot(3,2,1,'align');
%semilogy(t*dc/(year*vp),v*vp,'-k');
semilogy(t,v,'-k');
%xlim([0 5]);
ylim([1e-30 1e10]);
xlabel('time');ylabel('Slip rate');

subplot(3,2,2,'align');
%semilogy(t*dc/(year*vp),th*dc/vp,'-k');
plot(t,-par.k*(par.c*u-t-par.u0),'-k');%plot(t,-par.k*(par.c*u-t-par.u0+dbc),'-k');

%ylim([0 100]);
%xlim([0 5]);
xlabel('time ');ylabel('Tension')

subplot(3,2,3,'align');
plot(t,par.c*u,'-k');
hold on
plot(t,(u*par.c),'--b')
hold on
plot(t,tdw,'--c')
plot(t+1000,t,'--r');
%xlim([0 5]);
xlabel('time');ylabel('Extension')
leg=legend('Fault','Total(Fault)','cumuldike','Location','NorthWest');
set(leg,'Interpreter','Latex');

subplot(3,2,4,'align');
plot(t,s,'-b');
hold on
plot(t,tau*par.b,'-k')
%ylim([0 100]);xlim([0 5]);
xlabel('time');ylabel('stress')
leg=legend('Normal $\sigma$','Shear $\tau$','Location','SouthWest');
set(leg,'Interpreter','Latex');

cumulsigop=[];
for i=1:tpsi2-1
    load(['results/' num2str(i)])
    cumulsigop=[cumulsigop sigOP];
end
cumulsigop=reshape(cumulsigop,[],tpsi2-1);
subplot(3,2,5,'align');
ssigOP=size(sigOP);
h=imagesc(cumulsigop','XData',[X(1,ixOP(1),1),X(1,ixOP(ssigOP(2)),1)],'YData',[0,par.tfinal]);
colormap((redblue))%hot;gray
%h.FaceColor = 'interp';
caxis([(param.rho*param.g*z0/(param.rho*param.g*param.z))*2,0])
colorbar
hold on
plot(catdike.where,catdike.t,'p','MarkerEdgeColor','yellow','MarkerFaceColor','yellow','MarkerSize',15)
%plot(200,catseisme.t,'>','MarkerEdgeColor','white','MarkerFaceColor','white','MarkerSize',5)
xlabel('position along sill');ylabel('time')
set(gca,'Ydir','normal')
ylim([1100 t]);
%xlim([200 280]);
%leg=legend('Normal $\sigma$','Shear $\tau$','Location','SouthWest');
%set(leg,'Interpreter','Latex');

cumulsigop_size=size(cumulsigop);
subplot(3,2,6,'align');
r=plot(1:tpsi2-1,cumulsigop(1,:),'color','k');%find(sigOP==max(sigOP)) ; find(xOP==catdike.where(1))
hold on
r1=plot(1:tpsi2-1,cumulsigop(round(cumulsigop_size(1)/2),:),'g--');%find(sigOP==max(sigOP)) ; find(xOP==catdike.where(1))
hold on
r2=plot(1:tpsi2-1,cumulsigop(cumulsigop_size(1),:),'co');%find(sigOP==max(sigOP)) ; find(xOP==catdike.where(1))
set(r,'XData',round((0:par.tfinal/(tpsi2-2):par.tfinal)))
set(r1,'XData',round((0:par.tfinal/(tpsi2-2):par.tfinal)))
set(r2,'XData',round((0:par.tfinal/(tpsi2-2):par.tfinal)))
xlabel('time');ylabel('sig_xx value around the dike position')
set(gca,'Ydir','normal')
%leg=legend('Normal $\sigma$','Shear $\tau$','Location','SouthWest');
%set(leg,'Interpreter','Latex');

%h=pcolor([X(1,ixOP(1),1):round(X(1,ixOP(ssigOP(2)),1)/ssigOP(2)):X(1,ixOP(ssigOP(2)),1)],round([0:par.tfinal/(isave-2):par.tfinal]),cumulsigop');

%figure(1);
%print(1,'-depsc','fig_qdnormalf_v2.eps');

% figure(2)
% [Xsop,Zsop]=meshgrid(1:55,1:iend);
% h=imagesc(1:55,1:101,cumulsigop); %[1:55],[1:iend], Xsop,Zsop
% colormap(redblue)
% h.FaceColor = 'interp';
% caxis([-3,1])
% colorbar
