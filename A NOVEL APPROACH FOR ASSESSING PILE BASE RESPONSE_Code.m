%A NOVEL APPROACH FOR ASSESSING PILE BASE RESPONSE
clear all
close all
%pile diameter [m]
d = 31.75/1000;
rp = d/2;
Ap = pi*d^2/4;
%vertical effective stress [MPa]
sv0 = (490/1000*17+50)/1000;
%relative density
Dr0 = 0.87;
%Maximum and minimum void ratio
emax = 0.985;
emin = 0.611;
e0 = -(emax-emin)*Dr0+emax;
%critical state friction angle [°]
ficv = 31.6;
%coefficient of lateral earth pressure at rest
%k0 = 1-sin(degtorad(ficv));
k0 = 0.4;
%mean effective pressure
p0 = sv0*(1+2*k0)/3;
%initial ocrahedral shear stress [MPa]
toct0 = sqrt(2)/3*(sv0-k0*sv0);
%Hardin & Black (1966) coefficient
Cg = 900;
eg = 2.17;
ng = 0.40;
%soil stiffness [MPa]
G0 = Cg*(eg-e0)^2/(1+e0)*(p0/0.101)^ng*0.101;
%G0 = 50;
%Poisson's coefficient
ni = 0.15;
%Bulk modulus [MPa]
K0 = 2*G0*(1+ni)/(3*(1-2*ni));
%Lee & Salgado parameter
f = 1-0.07*Dr0;
g = 0.5*Dr0+0.0095;
%g = 0.20;
%Papadimitriou et al. parameters
%alfa1 = 0.5;
%g1 = 0.00055*((1-alfa1)/(2*alfa1^2));
%Bolton parameters
A = 3;
Qr = 10;
R = 1;
%initial angle of dilatancy [°]
psi0 = (A*(Dr0*(Qr-log(p0*1000))-R))/0.8;
if psi0 < 0
   psi0 = 0;
end
if psi0 > 20
   psi0 = 20;
end
%initial peak friction angle [°]
fip0 = ficv+0.8*psi0;

%Maximum load [MN]
Qmax = 0.00255;

%numbers of load steps
nstep = 100;
dQ = Qmax/nstep;
dQb = dQ/Ap;
for i = 1 : nstep+1
    Q(i,1)=0+dQ*(i-1);
    Q_kN(i,1)=Q(i,1)*1000;
    Qb(i,1)=Q(i,1)/Ap;
end

%Maximum relative depth
Zmax_d = 10;
Zmax = Zmax_d*d;
%Maximum relative radial distance
Rmax_d = 5;
Rmax = Rmax_d*d;
%relative depth discretization
dz_d = 0.1/5;
%relative radial discretization
dr_d = 0.0625;


dz=dz_d*d;
dr=dr_d*d;
%Coordinates of each node
for i = 1 : (Rmax_d/dr_d)+1
    for j = 1 : (Zmax_d/dz_d)+1
        r(j+(i-1)*((Zmax_d/dz_d)+1),1)=0+dr*(i-1);
        z(j+(i-1)*((Zmax_d/dz_d)+1),1)=0+dz*(j-1);
    end
end

for i = 1 : length(r)
    if r(i,1) == 0
    Area(i,1)=0;
    else
    Area(i,1)=pi*r(i,1)^2-pi*r(i-(Zmax_d/dz_d)-1,1)^2;
    end
    if z(i,1) == 0
    Area(i,1)=0;
    end
    V(i,1)=Area(i,1)*dz;
end

%number of nodes
n = length(z);

pilex=[0;d/2;d/2];
piley=[0;0;-2*d];
pilex2=[0;d/2;d/2;0];
piley2=[0;0;-2*d;-2*d];

figure (1)
hold on
grid on
title('Mesh','fontsize',14,'fontname','Georgia','color','k')
xlabel('horizontal distance from the pile center [m]','fontsize',12,'fontname','Georgia','color','k')
ylabel('vertical distance from the pile base [m]','fontsize',12,'fontname','Georgia','color','k')
xlim([0,(Rmax_d*d)+1*dr])
ylim([-2*d,(Zmax_d*d)+1*dz])
set (gca, 'ydir', 'reverse')
fill(pilex2,piley2,[0.9 0.9 0.9])
plot(pilex,piley,'linewidth',5,'color',[0.5 0.5 0.5])
plot([0 0],[0 (Zmax_d*d)+1*max(dz:dz)],':','linewidth',1.5,'color','b')
plot([0 (Rmax_d*d)+1*max(dz:dz)],[0 0],':','linewidth',1.5,'color','b')
scatter(r,z,2,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',1,'MarkerFaceColor',[0 0 0])

Gav(1:nstep+1,1)=zeros;
Kav(1:nstep+1,1)=zeros;
Ub_rel_c = 0;
for j = 1 : nstep+1
    sum_dDGt = 0;
    sum_dD = 0;
    sum_dUKt = 0;
    sum_dU = 0;
for i = 1 : n
    Gmax(i,1)=G0;
    Gs(i,1)=G0;
    Ks(i,1)=K0;
    Gav(1,1)=G0;
    Kav(1,1)=K0;
    psi(i,1)=psi0;
    fi(i,1)=fip0;
    e(i,1)=e0;
    Dr(i,1)=Dr0;
    p(i,1)=p0;
    zn=z(i,1);
    rn=r(i,1);
    S0(i,1) = asin(2/(sqrt((zn/rp)^2+(1+rn/rp)^2)+sqrt((zn/rp)^2+(1-rn/rp)^2)));
    N0(i,1) = ((zn^2/rp^2+rn^2/rp^2-1)^2+4*(zn/rp)^2)^(1/4);
    T0(i,1) = atan((2*zn/rp)/(zn^2/rp^2+rn^2/rp^2-1));
    if T0(i,1)<0
        T0(i,1)= T0(i,1)+pi;
    end
    if zn == 0 
    if rn == 0
    T0(i,1) = pi;
    end
    end
    if zn == 0 
    if rn == rp
    T0(i,1) = 0;
    end
    end
    I101(i,1)=sqrt(2/pi)*S0(i,1);
    I103(i,1)=sqrt(2/pi)*(N0(i,1)^(-1))*sin(T0(i,1)/2);
    if zn == 0 
    if rn == rp
    I103(i,1) = 0;
    end
    end
    I121(i,1)=sqrt(2/pi)*(rp/rn)*(1-N0(i,1)*sin(T0(i,1)/2));
    if rn == 0
    I121(i,1) = 0;
    end
    if zn == 0
    I121(i,1) = 0;
    end
    I123(i,1)=sqrt(2/pi)*(rp/rn)*N0(i,1)^-1*(cos(T0(i,1)/2)-zn/rp*sin(T0(i,1)/2));
    if rn == 0
    I123(i,1) = 0;
    end
    if zn == 0
    I123(i,1) = 0;
    end
    I105(i,1)=sqrt(2/pi)*N0(i,1)^-3*(zn/rp*sin(3*T0(i,1)/2)-cos(3*T0(i,1)/2));
    if zn == 0
    I105(i,1) = 0;
    end
    I125(i,1)=sqrt(2/pi)*(rn/rp)*N0(i,1)^-3*sin(3*T0(i,1)/2);
    if zn == 0
    I125(i,1) = 0;
    end
    if j == 1
        uz(i,j)=0;
        ur(i,j)=0;
        dsz(i,j)=0;
        sz(i,j)=0;
        dsr(i,j)=0;
        sr(i,j)=0;
        dsh(i,j)=0;
        sh(i,j)=0;
        dtrz(i,j)=0;
        trz(i,j)=0;
        dez(i,j)=0;
        ez(i,j)=0;
        der(i,j)=0;
        er(i,j)=0;
        deh(i,j)=0;
        eh(i,j)=0;
        dgrz(i,j)=0;
        grz(i,j)=0;
        devol(i,j)=0;
        evol(i,j)=0;
        Gav_cal=G0;
        Kav_cal=K0;
    end
    if j > 1
    Gav_cal=Gav(j-1,1);
    Kav_cal=Kav(j-1,1);
    uz(i,j)=(sqrt(2/pi)*dQ/(4*rp)*((3*Kav_cal+4*Gav_cal)/(Gav_cal*(6*Kav_cal+2*Gav_cal))*I101(i,1)+1/(2*Gav_cal)*(zn/rp)*I103(i,1)))+uz(i,j-1);
    ur(i,j)=(sqrt(2/pi)*dQ/(4*rp)*(-3/(6*Kav_cal+2*Gav_cal)*I121(i,1)+1/(2*Gav_cal)*(zn/rp)*I123(i,1)))+ur(i,j-1);
    Ub(j,1)=uz(1,j);
    Ub_mm(j,1)=Ub(j,1)*1000;
    Ub_rel(j,1)=Ub(j,1)/d*100;
    Ub_rel_c=Ub_rel(j,1);
    %Elastic solution
    Ub_el(j,1)=Q(j,1)/(2*d*G0/(1-ni));
    dsz(i,j)=(sqrt(2/pi)*dQ/(4*rp^2)*(I103(i,1)+(zn/rp)*I105(i,1)));
    sz(i,j)=dsz(i,j)+sz(i,j-1);
    if rn<rp || rn == rp
        if zn == 0
            dsz(i,j) = dQb;
            sz(i,j) = Qb(j,1);
        end
    end
    dsr(i,j)=(sqrt(2/pi)*dQ/(4*rp^2)*(I103(i,1)-(zn/rp)*I105(i,1)-2*Gav_cal*(rp/rn)*(3/(6*Kav_cal+2*Gav_cal)*I121(i,1)-1/(2*Gav_cal)*(zn/rp)*I123(i,1))));
    sr(i,j)=dsr(i,j)+sr(i,j-1);
    if rn == 0
        dsr(i,j) = 0;
        sr(i,j) = 0;
    end
    dsh(i,j)=(sqrt(2/pi)*dQ/(4*rp^2)*((3*Kav_cal-2*Gav_cal)/(3*Kav_cal+Gav_cal)*I103(i,1)+2*Gav_cal*(rp/rn)*(3/(6*Kav_cal+2*Gav_cal)*I121(i,1)-1/(2*Gav_cal)*(zn/rp)*I123(i,1))));
    sh(i,j)=dsh(i,j)+sh(i,j-1);
    if rn == 0
        dsh(i,j)=0;
        sh(i,j) = 0;
    end
    dtrz(i,j)=(sqrt(2/pi)*dQ/(4*rp^2)*(zn/rp)*I125(i,1));
    trz(i,j)=dtrz(i,j)+trz(i,j-1);
    Matrix(1,1)=Kav_cal+4*Gav_cal/3;
    Matrix(2,1)=Kav_cal-2*Gav_cal/3;
    Matrix(3,1)=Kav_cal-2*Gav_cal/3;
    Matrix(4,1)=0;
    Matrix(1,2)=Kav_cal-2*Gav_cal/3;
    Matrix(2,2)=Kav_cal+4*Gav_cal/3;
    Matrix(3,2)=Kav_cal-2*Gav_cal/3;
    Matrix(4,2)=0;
    Matrix(1,3)=Kav_cal-2*Gav_cal/3;
    Matrix(2,3)=Kav_cal-2*Gav_cal/3;
    Matrix(3,3)=Kav_cal+4*Gav_cal/3;
    Matrix(4,3)=0;
    Matrix(1,4)=0;
    Matrix(2,4)=0;
    Matrix(3,4)=0;
    Matrix(4,4)=Gav_cal;
    Vett_tens(1,1)=dsr(i,j);
    Vett_tens(2,1)=dsh(i,j);
    Vett_tens(3,1)=dsz(i,j);
    Vett_tens(4,1)=dtrz(i,j);
    Vett_def=Vett_tens'*inv(Matrix);
    der(i,j)=Vett_def(1,1);
    deh(i,j)=Vett_def(1,2);
    dez(i,j)=Vett_def(1,3);
    dgrz(i,j)=Vett_def(1,4);
    %dez(i,j)=(sqrt(2/pi)*dQ/(4*rp^2)*(3/(6*Kav_cal+2*Gav_cal)*I103(i,1)+1/(2*Gav_cal)*(zn/rp)*I105(i,1)));
    ez(i,j)=dez(i,j)+ez(i,j-1);
    %der(i,j)=(sqrt(2/pi)*dQ/(4*rp^2)*(3/(6*Kav_cal+2*Gav_cal)*(I103(i,1)-(rp/rn)*I121(i,1))-1/(2*Gav_cal)*(zn/rp)*(I105(i,1)-(rp/rn)*I123(i,1))));
    er(i,j)=der(i,j)+er(i,j-1);
    %deh(i,j)=(sqrt(2/pi)*dQ/(4*rp^2)*(3/(6*Kav_cal+2*Gav_cal)*(rp/rn)*I121(i,1)-1/(2*Gav_cal)*(zn/rp)*I123(i,1)));
    eh(i,j)=deh(i,j)+eh(i,j-1);
    %dgrz(i,j)=(sqrt(2/pi)*dQ/(4*rp^2)*(1/Gav_cal)*(zn/rp)*I125(i,1));
    grz(i,j)=dgrz(i,j)+grz(i,j-1);
    if rn == 0
        der(i,j)=0;
        er(i,j) = 0;
        deh(i,j)=0;
        eh(i,j) = 0;
    end
    end
    devol(i,j)=dez(i,j)+der(i,j)+deh(i,j);
    evol(i,j) = ez(i,j)+er(i,j)+eh(i,j);
    %keep in mind that we are computing the (cumulative) increment of octaedral shear
    %stress respect the initial one
    toct(i,j)=1/3*sqrt((sr(i,j)-sh(i,j))^2+(sz(i,j)-sr(i,j))^2+(sz(i,j)-sh(i,j))^2+6*trz(i,j)^2);
    e(i,j)=-(evol(i,j)*(1+e0)-e0);
    if e(i,j) < emin
       e(i,j) = emin;
    end
    if e(i,j) > emax
       e(i,j) = emax;
    end
    p(i,j)=p0+(sz(i,j)+sr(i,j)+sh(i,j))/3;
    dp(i,j)=(dsz(i,j)+dsr(i,j)+dsh(i,j))/3;
    dU(i,j)=1/2*dp(i,j)*devol(i,j);
    dD(i,j)=1/2*(dsz(i,j)*dez(i,j)+dsr(i,j)*der(i,j)+dsh(i,j)*deh(i,j)+dtrz(i,j)*dgrz(i,j)-dp(i,j)*devol(i,j));
    Dr(i,j)=(emax-e(i,j))/(emax-emin);
    psi(i,j)=(A*(Dr(i,j)*(Qr-log(p(i,j)*1000))-R))/0.8;
    if psi(i,j)<0
       psi(i,j)=0;
    end
    fi(i,j)=ficv+0.8*psi(i,j);
    toct_max(i,j)=2*sqrt(2)*sin(degtorad(fi(i,j)))/(3-sin(degtorad(fi(i,j))))*p(i,j);
    if toct(i,j)>toct_max(i,j)-toct0
       toct(i,j)=toct_max(i,j)-toct0;
       p(i,j)=(toct_max(i,j)*sqrt(3/2)/(2*sin(degtorad(fi(i,j)))/(sqrt(3)*(3-sin(degtorad(fi(i,j)))))))/3;
    end
    %Gmax(i,j)=Cg*(eg-e(i,j))^2/(1+e(i,j))*(p(i,j)/0.101)^ng*0.101;
    %toct1(i,j)=sqrt(2)/3*2*alfa1*G0*g1*(p(i,j)/p0);
    %toct1(i,j)=sqrt(2)/3*2*alfa1*G0*g1;
    %toct1(i,j)=sqrt(2)/3*2*alfa1*Gmax(i,j)*g1;
    %Csi1(i,j)=(toct1(i,j))/(toct_max(i,j)-toct0);
    %g(i,j)=log(1/f*(1-alfa1/(p(i,j)/p0)^ng))/log(Csi1(i,j));
    %g(i,j)=0.25;
    %Gmax(i,j)=Cg*(eg-e(i,j))^2/(1+e(i,j))*(p(i,j)/0.101)^ng*0.101;
    Gs(i,j)=((1-f*(toct(i,j)/(toct_max(i,j)-toct0))^g)*(p(i,j)/p0)^ng)*G0;
    Gt(i,j)=Gs(i,j)/(1+(1-f*(toct(i,j)/(toct_max(i,j)-toct0))^g)^-1*f*g*(toct(i,j)/(toct_max(i,j)-toct0))^g);
    Kt(i,j)=(p(i,j)/p0)^0.5*K0;
    Ks(i,j)=K0*(1-0.5)*(p(i,j)-p0)/(p0^0.5*p(i,j)^(1-0.5)-p0);
    nis(i,j)=(3*Ks(i,j)-2*Gs(i,j))/(2*(3*Ks(i,j)+Gs(i,j)));
    if p(i,j)==p0
       Ks(i,j)=K0;
    end
    sum_dDGt=sum_dDGt+dD(i,j)*Gt(i,j)*V(i,1);
    sum_dD=sum_dD+dD(i,j)*V(i,1);
    Gav(j,1)=sum_dDGt/sum_dD;
    sum_dUKt=sum_dUKt+dU(i,j)*Kt(i,j)*V(i,1);
    sum_dU=sum_dU+dU(i,j)*V(i,1);
    Kav(j,1)=sum_dUKt/sum_dU;
    GtG0(j,1) = Gav(j,1)/G0;
    KtK0(j,1) = Kav(j,1)/K0;
end
end

nstep_2=length(Ub_rel);
for i = 1 : (Rmax_d/dr_d)+1
    j(i,1)= (Zmax_d/dz_d+1)*(i-1)+1;
    uz_z0(i,1)=uz(j(i),nstep_2);
    uz_z05(i,1)=uz(j(i)+0.5/dz_d,nstep_2);
    ur_z05(i,1)=ur(j(i)+0.5/dz_d,nstep_2);
    uz_z1(i,1)=uz(j(i)+1/dz_d,nstep_2);
    uz_z2(i,1)=uz(j(i)+2/dz_d,nstep_2);
    r_z0(i,1)=r(j(i),1);
end

%zv=[0:0.2:4];
%for zvz = 1 : length(zv)
%for i = 1 : (Rmax_d/dr_d)+1
%    j(i,1)= (Zmax_d/dz_d+1)*(i-1)+1;
%    zvvar=zv(zvz);
%    uz_zv(i,zvz)=uz(round(j(i)+zvvar/dz_d),nstep_2);
%    ur_zv(i,zvz)=ur(round(j(i)+zvvar/dz_d),nstep_2);
%    zdis_zv(i,zvz)=uz_zv(i,zvz)+zv(zvz)*d;
%    zdis_zv_d(i,zvz)=zdis_zv(i,zvz)/d;
%    rdis_zv(i,zvz)=ur_zv(i,zvz)+r(j(i),1);
%    rdis_zv_d(i,zvz)=rdis_zv(i,zvz)/d;
%    z0_zv_d(i,zvz)=zv(zvz);
%    r_z0_d(i,1)=r(j(i),1)/d;
%end
%end

figure (2)
hold on
grid on
title('load-settlement curve','fontsize',14,'fontname','Georgia','color','k')
xlabel('load, Q [MN]','fontsize',12,'fontname','Georgia','color','k')
ylabel('relative settlement, w/d [%]','fontsize',12,'fontname','Georgia','color','k')
%xlim([0,(Rmax_d*d)+1*dr])
%ylim([-2*d,(Zmax_d*d)+1*dz])
set (gca, 'ydir', 'reverse')
plot(Q,Ub/d*100,'linewidth',2,'color','k')
plot(Q,Ub_el/d*100,'linewidth',2,'color','r')
scatter(Q,Ub/d*100,20,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',1,'MarkerFaceColor',[0 0 0])
legend('proposed model','elastic solution')

figure (3)
hold on
grid on
%title('load-settlement','fontsize',20,'fontname','Georgia','color','k')
xlabel('horizontal distance from pile center [m]','fontsize',12,'fontname','Georgia','color','k')
ylabel('vertical displacement [mm]','fontsize',12,'fontname','Georgia','color','k')
%xlim([0,(Rmax_d*d)+1*dr])
ylim([0,1.2*max(Ub*1000)])
set (gca, 'ydir', 'reverse')
plot(r_z0,uz_z0*1000,'linewidth',2,'color','k')
plot(r_z0,uz_z05*1000,'-.','linewidth',1,'color','k')
plot(r_z0,uz_z1*1000,'--','linewidth',1,'color','k')
plot(r_z0,uz_z2*1000,'-','linewidth',1,'color','k')
legend('z/d=0','z/d=0.5','z/d=1','z/d=2')
%scatter(r_z0,uz_z0*1000,20,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',1,'MarkerFaceColor',[0 0 0])

figure (4)
hold on
grid on
xlabel('horizontal distance from the pile center [m]','fontsize',14,'fontname','Georgia','color','k')
ylabel('vertical distance from the pile base [m]','fontsize',14,'fontname','Georgia','color','k')
xlim([0,3*d])
ylim([-d,3*d])
set (gca, 'ydir', 'reverse')
fill(pilex2,piley2,[0.9 0.9 0.9])
plot(pilex,piley,'linewidth',5,'color',[0.5 0.5 0.5])
plot([0 0],[0 (Zmax_d*d)+1*max(dz:dz)],':','linewidth',1.5,'color','b')
plot([0 (Rmax_d*d)+1*max(dz:dz)],[0 0],':','linewidth',1.5,'color','b')
%scatter(r,z,2,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',1,'MarkerFaceColor',[0 0 0])
quiver(r,z,ur(:,nstep+1),uz(:,nstep+1),2,'k','linewidth',1)
axis equal

[rq,zq] = meshgrid(0:dr/5:Rmax, 0:dz/5:Zmax);
sz_graph = griddata(r,z,sz(:,nstep+1),rq,zq);
p_graph = griddata(r,z,p(:,nstep+1),rq,zq);
toct_graph = griddata(r,z,toct(:,nstep+1),rq,zq);
trz_graph = griddata(r,z,trz(:,nstep+1),rq,zq);
GsG0_graph = griddata(r,z,Gs(:,nstep+1)/G0,rq,zq);
KtG0_graph = griddata(r,z,Ks(:,nstep+1)/K0,rq,zq);
psi_graph = griddata(r,z,psi(:,nstep+1),rq,zq);
%nis_graph = griddata(r,z,nis(:,nstep+1),rq,zq);
%g_graph = griddata(r,z,g(:,nstep+1),rq,zq);

figure (5)
hold on
grid on
title('vertical stress [MPa]','fontsize',14,'fontname','Georgia','color','k')
xlabel('horizontal distance from the pile center [m]','fontsize',12,'fontname','Georgia','color','k')
ylabel('vertical distance from the pile base [m]','fontsize',12,'fontname','Georgia','color','k')
xlim([-d,5*d])
ylim([-d,5*d])
set (gca, 'ydir', 'reverse')
fill(pilex2,piley2,[0.9 0.9 0.9])
plot(pilex,piley,'linewidth',5,'color',[0.5 0.5 0.5])
plot([0 0],[0 (Zmax_d*d)+1*max(dz:dz)],':','linewidth',1.5,'color','b')
plot([0 (Rmax_d*d)+1*max(dz:dz)],[0 0],':','linewidth',1.5,'color','b')
surf(rq,zq,sz_graph)
shading flat
colormap jet
colorbar
%view(2)
%scatter(r,z,1,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',1,'MarkerFaceColor',[0 0 0])
%axis equal

figure (6)
hold on
grid on
title('mean pressure [MPa]','fontsize',14,'fontname','Georgia','color','k')
xlabel('horizontal distance from the pile center [m]','fontsize',12,'fontname','Georgia','color','k')
ylabel('vertical distance from the pile base [m]','fontsize',12,'fontname','Georgia','color','k')
xlim([-d,5*d])
ylim([-d,5*d])
set (gca, 'ydir', 'reverse')
fill(pilex2,piley2,[0.9 0.9 0.9])
plot(pilex,piley,'linewidth',5,'color',[0.5 0.5 0.5])
plot([0 0],[0 (Zmax_d*d)+1*max(dz:dz)],':','linewidth',1.5,'color','b')
plot([0 (Rmax_d*d)+1*max(dz:dz)],[0 0],':','linewidth',1.5,'color','b')
surf(rq,zq,p_graph)
shading flat
colormap jet
colorbar
%view(2)

figure (7)
hold on
grid on
title('octaedral stress [MPa]','fontsize',14,'fontname','Georgia','color','k')
xlabel('horizontal distance from the pile center [m]','fontsize',12,'fontname','Georgia','color','k')
ylabel('vertical distance from the pile base [m]','fontsize',12,'fontname','Georgia','color','k')
xlim([-d,5*d])
ylim([-d,5*d])
set (gca, 'ydir', 'reverse')
fill(pilex2,piley2,[0.9 0.9 0.9])
plot(pilex,piley,'linewidth',5,'color',[0.5 0.5 0.5])
plot([0 0],[0 (Zmax_d*d)+1*max(dz:dz)],':','linewidth',1.5,'color','b')
plot([0 (Rmax_d*d)+1*max(dz:dz)],[0 0],':','linewidth',1.5,'color','b')
surf(rq,zq,toct_graph)
shading flat
colormap jet
colorbar
%view(2)

figure (8)
hold on
grid on
title('shear stress in r-z plane [MPa]','fontsize',14,'fontname','Georgia','color','k')
xlabel('horizontal distance from the pile center [m]','fontsize',12,'fontname','Georgia','color','k')
ylabel('vertical distance from the pile base [m]','fontsize',12,'fontname','Georgia','color','k')
xlim([-d,5*d])
ylim([-d,5*d])
set (gca, 'ydir', 'reverse')
fill(pilex2,piley2,[0.9 0.9 0.9])
plot(pilex,piley,'linewidth',5,'color',[0.5 0.5 0.5])
plot([0 0],[0 (Zmax_d*d)+1*max(dz:dz)],':','linewidth',1.5,'color','b')
plot([0 (Rmax_d*d)+1*max(dz:dz)],[0 0],':','linewidth',1.5,'color','b')
surf(rq,zq,trz_graph)
shading flat
colormap jet
colorbar

figure (9)
hold on
grid on
title('Gs/G0','fontsize',14,'fontname','Georgia','color','k')
xlabel('horizontal distance from the pile center [m]','fontsize',12,'fontname','Georgia','color','k')
ylabel('vertical distance from the pile base [m]','fontsize',12,'fontname','Georgia','color','k')
xlim([-d,10*d])
ylim([-d,10*d])
set (gca, 'ydir', 'reverse')
fill(pilex2,piley2,[0.9 0.9 0.9])
plot(pilex,piley,'linewidth',5,'color',[0.5 0.5 0.5])
plot([0 0],[0 (Zmax_d*d)+1*max(dz:dz)],':','linewidth',1.5,'color','b')
plot([0 (Rmax_d*d)+1*max(dz:dz)],[0 0],':','linewidth',1.5,'color','b')
surf(rq,zq,GsG0_graph)
shading flat
colormap jet
colorbar

figure (10)
hold on
grid on
title('Ks/K0','fontsize',14,'fontname','Georgia','color','k')
xlabel('horizontal distance from the pile center [m]','fontsize',12,'fontname','Georgia','color','k')
ylabel('vertical distance from the pile base [m]','fontsize',12,'fontname','Georgia','color','k')
xlim([-d,10*d])
ylim([-d,10*d])
set (gca, 'ydir', 'reverse')
fill(pilex2,piley2,[0.9 0.9 0.9])
plot(pilex,piley,'linewidth',5,'color',[0.5 0.5 0.5])
plot([0 0],[0 (Zmax_d*d)+1*max(dz:dz)],':','linewidth',1.5,'color','b')
plot([0 (Rmax_d*d)+1*max(dz:dz)],[0 0],':','linewidth',1.5,'color','b')
surf(rq,zq,KtG0_graph)
shading flat
colormap jet
colorbar

figure (11)
hold on
grid on
title('angle of dilatancy [°]','fontsize',14,'fontname','Georgia','color','k')
xlabel('horizontal distance from the pile center [m]','fontsize',12,'fontname','Georgia','color','k')
ylabel('vertical distance from the pile base [m]','fontsize',12,'fontname','Georgia','color','k')
xlim([-d,5*d])
ylim([-d,5*d])
set (gca, 'ydir', 'reverse')
fill(pilex2,piley2,[0.9 0.9 0.9])
plot(pilex,piley,'linewidth',5,'color',[0.5 0.5 0.5])
plot([0 0],[0 (Zmax_d*d)+1*max(dz:dz)],':','linewidth',1.5,'color','b')
plot([0 (Rmax_d*d)+1*max(dz:dz)],[0 0],':','linewidth',1.5,'color','b')
surf(rq,zq,psi_graph)
shading flat
colormap jet
colorbar

% figure (12)
% hold on
% grid on
% title('coefficient g','fontsize',14,'fontname','Georgia','color','k')
% xlabel('horizontal distance from the pile center [m]','fontsize',12,'fontname','Georgia','color','k')
% ylabel('vertical distance from the pile base [m]','fontsize',12,'fontname','Georgia','color','k')
% xlim([-d,5*d])
% ylim([-d,5*d])
% set (gca, 'ydir', 'reverse')
% fill(pilex2,piley2,[0.9 0.9 0.9])
% plot(pilex,piley,'linewidth',5,'color',[0.5 0.5 0.5])
% plot([0 0],[0 (Zmax_d*d)+1*max(dz:dz)],':','linewidth',1.5,'color','b')
% plot([0 (Rmax_d*d)+1*max(dz:dz)],[0 0],':','linewidth',1.5,'color','b')
% surf(rq,zq,g_graph)
% shading flat
% colormap jet
% colorbar