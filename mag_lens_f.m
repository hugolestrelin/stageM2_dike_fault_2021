function [sig] = mag_lens_f(rx,ry,P,L,nu,wd,x0) % rx and ry matrix n*m*o ,tau,alpha]
    r=((rx-x0).^2+(ry-wd).^2).^(1/2);
    r1=((rx-x0+L/2).^2+(ry-wd).^2).^(1/2); %gauche <0
    r2=((L/2-rx-x0).^2+(ry-wd).^2).^(1/2); %droite >0
    theta=zeros(size(rx));
    theta1=zeros(size(rx));
    theta2=zeros(size(rx));
    theta(ry>wd & rx-x0==0)=rad2deg(pi/2);
    theta1(ry>wd & rx-x0==0)=atand(abs(ry(ry>wd & rx-x0==0)-wd)./abs(rx(ry>wd & rx-x0==0)-x0+L/2));
    theta2(ry>wd & rx-x0==0)=rad2deg(pi)-(atand(abs(ry(ry>wd & rx-x0==0)-wd)./abs(L/2-rx(ry>wd & rx-x0==0)-x0)));
    theta(ry>wd & (0<rx-x0)&(rx-x0<L/2))=(atand(abs(ry(ry>wd & (0<rx-x0)&(rx-x0<L/2))-wd)./abs(rx(ry>wd & (0<rx-x0)&(rx-x0<L/2))-x0)));
    theta1(ry>wd & (0<rx-x0)&(rx-x0<L/2))=(atand(abs(ry(ry>wd & (0<rx-x0)&(rx-x0<L/2))-wd)./abs(rx(ry>wd & (0<rx-x0)&(rx-x0<L/2))-x0+L/2)));
    theta2(ry>wd & (0<rx-x0)&(rx-x0<L/2))=rad2deg(pi)-(atand(abs(ry(ry>wd & (0<rx-x0)&(rx-x0<L/2))-wd)./abs(L/2-rx(ry>wd & (0<rx-x0)&(rx-x0<L/2))-x0)));
    theta(ry>wd & rx-x0==L/2)=(atand(abs(ry(ry>wd & rx-x0==L/2)-wd)./abs(rx(ry>wd & rx-x0==L/2)-x0)));
    theta1(ry>wd & rx-x0==L/2)=(atand(abs(ry(ry>wd & rx-x0==L/2)-wd)./abs(rx(ry>wd & rx-x0==L/2)-x0+L/2)));
    theta2(ry>wd & rx-x0==L/2)=rad2deg((pi/2));
    theta(ry>wd & rx-x0>L/2)=(atand(abs(ry(ry>wd & rx-x0>L/2)-wd)./abs(rx(ry>wd & rx-x0>L/2)-x0)));
    theta1(ry>wd & rx-x0>L/2)=(atand(abs(ry(ry>wd & rx-x0>L/2)-wd)./abs(rx(ry>wd & rx-x0>L/2)-x0+L/2)));
    theta2(ry>wd & rx-x0>L/2)=(atand(abs(ry(ry>wd & rx-x0>L/2)-wd)./abs(L/2-rx(ry>wd & rx-x0>L/2)-x0)));
    theta(ry>wd & (-L/2<rx-x0)&(rx-x0<0))=rad2deg(pi)-(atand(abs(ry(ry>wd & (-L/2<rx-x0)&(rx-x0<0))-wd)./abs(rx(ry>wd & (-L/2<rx-x0)&(rx-x0<0))-x0)));
    theta1(ry>wd & (-L/2<rx-x0)&(rx-x0<0))=(atand(abs(ry(ry>wd & (-L/2<rx-x0)&(rx-x0<0))-wd)./abs(rx(ry>wd & (-L/2<rx-x0)&(rx-x0<0))+L/2-x0)));
    theta2(ry>wd & (-L/2<rx-x0)&(rx-x0<0))=rad2deg(pi)-(atand(abs(ry(ry>wd & (-L/2<rx-x0)&(rx-x0<0))-wd)./abs(L/2-(rx(ry>wd & (-L/2<rx-x0)&(rx-x0<0))-x0))));
    theta(ry>wd & rx-x0==-L/2)=rad2deg(pi)-(atand(abs(ry(ry>wd & rx-x0==-L/2)-wd)./abs((rx(ry>wd & rx-x0==-L/2)-x0))));
    theta1(ry>wd & rx-x0==-L/2)=rad2deg(pi/2);
    theta2(ry>wd & rx-x0==-L/2)=rad2deg(pi)-(atand(abs(ry(ry>wd & rx-x0==-L/2)-wd)./abs(L/2-(rx(ry>wd & rx-x0==-L/2)-x0))));
    theta(ry>wd & rx-x0<-L/2)=rad2deg(pi)-(atand(abs(ry(ry>wd & rx-x0<-L/2)-wd)./abs(rx(ry>wd & rx-x0<-L/2)-x0)));
    theta1(ry>wd & rx-x0<-L/2)=rad2deg(pi)-(atand(abs(ry(ry>wd & rx-x0<-L/2)-wd)./abs(rx(ry>wd & rx-x0<-L/2)-x0+L/2)));
    theta2(ry>wd & rx-x0<-L/2)=rad2deg(pi)-(atand(abs(ry(ry>wd & rx-x0<-L/2)-wd)./abs(L/2-(rx(ry>wd & rx-x0<-L/2)-x0))));
    theta(ry==wd & rx-x0==0)=rad2deg(pi/2);
    theta2(ry==wd & rx-x0==0)=rad2deg(pi);
    theta2(ry==wd &(0<rx-x0)&(rx-x0<L/2))=rad2deg(pi);
    theta2(ry==wd & rx-x0==L/2)=rad2deg(pi/2);
    theta(ry==wd & (-L/2<rx-x0)&(rx-x0<0))=pi;
    theta2(ry==wd & (-L/2<rx-x0)&(rx-x0<0))=rad2deg(pi);
    theta(ry==wd & rx-x0==-L/2)=rad2deg(pi);
    theta1(ry==wd & rx-x0==-L/2)=rad2deg(pi/2);
    theta2(ry==wd & rx-x0==-L/2)=rad2deg(pi);
    theta(ry==wd & rx-x0<-L/2)=rad2deg(pi);
    theta1(ry==wd & rx-x0<-L/2)=rad2deg(pi);
    theta2(ry==wd & rx-x0<-L/2)=rad2deg(pi);
    theta(ry<wd & rx-x0==0)=rad2deg(3*pi/2);
    theta1(ry<wd & rx-x0==0)=rad2deg(2*pi)-atand(abs(wd-ry(ry<wd & rx-x0==0))./abs(rx(ry<wd & rx-x0==0)-x0+L/2));
    theta2(ry<wd & rx-x0==0)=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & rx-x0==0))./abs(L/2-rx(ry<wd & rx-x0==0)-x0)));
    theta(ry<wd & (0<rx-x0)&(rx-x0<L/2))=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & (0<rx-x0)&(rx-x0<L/2)))./abs(rx(ry<wd & (0<rx-x0)&(rx-x0<L/2))-x0)));
    theta1(ry<wd & (0<rx-x0)&(rx-x0<L/2))=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & (0<rx-x0)&(rx-x0<L/2)))./abs(rx(ry<wd & (0<rx-x0)&(rx-x0<L/2))-x0+L/2)));
    theta2(ry<wd & (0<rx-x0)&(rx-x0<L/2))=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & (0<rx-x0)&(rx-x0<L/2)))./abs(L/2-rx(ry<wd & (0<rx-x0)&(rx-x0<L/2))-x0)));
    theta(ry<wd & rx-x0==L/2)=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & rx-x0==L/2))./abs(rx(ry<wd & rx-x0==L/2)-x0)));
    theta1(ry<wd & rx-x0==L/2)=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & rx-x0==L/2))./abs(rx(ry<wd & rx-x0==L/2)-x0+L/2)));
    theta2(ry<wd & rx-x0==L/2)=rad2deg(3*pi/2);
    theta(ry<wd & rx-x0>L/2)=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & rx-x0>L/2))./abs(rx(ry<wd & rx-x0>L/2)-x0)));
    theta1(ry<wd & rx-x0>L/2)=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & rx-x0>L/2))./abs(rx(ry<wd & rx-x0>L/2)-x0+L/2)));
    theta2(ry<wd & rx-x0>L/2)=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & rx-x0>L/2))./abs(L/2-rx(ry<wd & rx-x0>L/2)-x0)));
    theta(ry<wd & (-L/2<rx-x0)&(rx-x0<0))=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & (-L/2<rx-x0)&(rx-x0<0)))./abs(rx(ry<wd & (-L/2<rx-x0)&(rx-x0<0))-x0)));
    theta1(ry<wd & (-L/2<rx-x0)&(rx-x0<0))=rad2deg(2*pi)-(atand(abs(wd-ry(ry<wd & (-L/2<rx-x0)&(rx-x0<0)))./abs(rx(ry<wd & (-L/2<rx-x0)&(rx-x0<0))-x0+L/2)));
    theta2(ry<wd & (-L/2<rx-x0)&(rx-x0<0))=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & (-L/2<rx-x0)&(rx-x0<0)))./abs(L/2-(rx(ry<wd & (-L/2<rx-x0)&(rx-x0<0))-x0))));
    theta(ry<wd & rx-x0==-L/2)=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & rx-x0==-L/2))./abs(rx(ry<wd & rx-x0==-L/2)-x0)));
    theta1(ry<wd & rx-x0==-L/2)=rad2deg(3*pi/2);
    theta2(ry<wd & rx-x0==-L/2)=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & rx-x0==-L/2))./abs(L/2-(rx(ry<wd & rx-x0==-L/2)-x0))));
    theta(ry<wd & rx-x0<-L/2)=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & rx-x0<-L/2))./abs(rx(ry<wd & rx-x0<-L/2)-x0)));
    theta1(ry<wd & rx-x0<-L/2)=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & rx-x0<-L/2))./abs(rx(ry<wd & rx-x0<-L/2)-x0+L/2)));
    theta2(ry<wd & rx-x0<-L/2)=rad2deg(pi)+(atand(abs(wd-ry(ry<wd & rx-x0<-L/2))./abs(L/2-(rx(ry<wd & rx-x0<-L/2)-x0))));
    a=L/2;
    sig.xx=P*(r./(r1.*r2).^(1/2).*(cosd(theta-(theta1+theta2)/2))-1-((a^2)*r./((r1.*r2).^(3/2))).*(sind(theta)).*(sind(3/2*(theta1+theta2))));
    sig.xy=P*(a^2)*r./(r1.*r2).^(3/2).*(sind(theta)).*(cosd(3/2*(theta1+theta2)));
    sig.yy=P*(r./(r1.*r2).^(1/2).*(cosd(theta-(theta1+theta2)/2))-1+((a^2)*r./((r1.*r2).^(3/2))).*(sind(theta)).*(sind(3/2*(theta1+theta2))));
    sig.zz=nu*(sig.xx+sig.yy);
end