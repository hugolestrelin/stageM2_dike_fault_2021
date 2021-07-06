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

iend=min([length(t) length(v) length(th) length(u) length(s) length(tau)]);
t=t(1:iend);v=v(1:iend);th=th(1:iend);u=u(1:iend);s=s(1:iend);tau=tau(1:iend);dbc=dbc(1:iend);dpc=dpc(1:iend);

sop=reshape(sop,iend,[]);

figure(1);clf;
subplot(2,2,1,'align');
%semilogy(t*dc/(year*vp),v*vp,'-k');
semilogy(t,v,'-k');
%xlim([0 5]);
ylim([1e-5 1e10]);
xlabel('time');ylabel('Slip rate');

subplot(2,2,2,'align');
%semilogy(t*dc/(year*vp),th*dc/vp,'-k');
plot(t,-par.k*(par.c*u-t-par.u0),'-k');%plot(t,-par.k*(par.c*u-t-par.u0+dbc),'-k');

%ylim([0 100]);
%xlim([0 5]);
xlabel('time ');ylabel('Tension')

subplot(2,2,3,'align');
plot(t,par.c*u,'-k');
hold on
plot(t,(u*par.c),'--b')
hold on
plot(t,tdw,'--c')
plot(t,t,'--r');
%xlim([0 5]);
xlabel('time');ylabel('Extension')
leg=legend('Fault','Total(Fault)','Tectonic','Location','NorthWest');
set(leg,'Interpreter','Latex');

subplot(2,2,4,'align');
plot(t,s,'-b');
hold on
plot(t,tau*par.b,'-k')
%ylim([0 100]);xlim([0 5]);
xlabel('time');ylabel('stress')
leg=legend('Normal $\sigma$','Shear $\tau$','Location','SouthWest');
set(leg,'Interpreter','Latex');

figure(1);
print(1,'-depsc','fig_qdnormalf_v2.eps');

% figure(2)
% [Xsop,Zsop]=meshgrid(1:55,1:iend);
% h=imagesc(1:55,1:iend,sop); %[1:55],[1:iend], Xsop,Zsop
% %colormap(redblue)
% %h.FaceColor = 'interp';
% caxis([-3,1])
% colorbar
