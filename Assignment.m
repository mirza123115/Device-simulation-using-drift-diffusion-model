clear all
close all
clc

%%
e=1.6*10e-19;
mun = 0.14;
mup = 0.04;
Na = 10^15;
Nd = 10^15;
Dn = 0.003618;
Dp = 0.001034;
Vt = 0.025855;
T = 300;
Xn = 6.1204e-7;
Xp = 6.1204e-7;
npo = 2.25e5;
pno = 2.25e5;
Ln = (Dn*5e-7)^0.5;
Lp = (Dp*5e-7)^0.5;
ni = 1.5e10;
pl = 5e-6;
nl = 5e-6;
nxp = 2.58e15;
pxn = 2.58e15;
ep0=8.854*10e-12;

%%
syms x
rop(x)=-e*Na+0*x;
a=linspace(0,-Xn,10);
subplot(2,3,1)
plot(a,rop(a))
hold on
b=linspace(0,Xp,10);
ron(x)=e*Nd+0*x;
plot(b,ron(b))
title('charge density vs x')
xlabel('device length ---> x') 
ylabel('charge density ron,rop')

exp=int(rop)/ep0;
exp=exp-e*Na*Xp/ep0;
subplot(232)
plot(a,exp(a))
hold on
exn=int(ron)/ep0;
exn=exn-e*Na*Xn/ep0;
plot(b,exn(b))
title('electric field vs x')
xlabel('device length ---> x') 
ylabel('Electric field Exp,Exn')

vxp=-int(exp);
vxp=vxp+(e*Na*Xp^2)/ep0;
subplot(233)
plot(a,vxp(a))
hold on
vxn=-int(exn);
vxn=vxn+(e*Na*Xp^2)/ep0;
plot(b,vxn(b))
title('voltage vs x')
xlabel('device length ---> x') 
ylabel('Voltage distribuation Vx')

%%
syms dpnx(x)
ode = diff(dpnx,x,2) -dpnx/(Lp^2) == 0;
cond1 = dpnx(-Xp) == nxp;
cond2 = dpnx(-pl) == npo;
conds = [cond1 cond2];
ySol(x) = dsolve(ode,conds);
ySol = simplify(ySol);
x1= linspace(-pl,-Xp,100);
z1=ySol(x1);
subplot(234)
plot(x1,z1)
hold on

syms dnpx(x)
ode2 = diff(dnpx,x,2) -dnpx/(Lp^2) == 0;
cond1 = dnpx(Xp) == nxp;
cond2 = dnpx(pl) == npo;

conds = [cond1 cond2];
ySol2(x) = dsolve(ode2,conds);
ySol2 = simplify(ySol2);

x2= linspace(Xn,pl,100);
z2=ySol2(x2);
plot(x2,z2)
title('monority carrier concentration vs x')
xlabel('device length ---> x') 
ylabel('Minoritty carrier dpnx,dnpx')



%%
jpn=e*Dn*diff(ySol);
subplot(235)
plot(x1,jpn(x1))
hold on
jnp=-e*Dn*diff(ySol2);
plot(x2,jnp(x2))
hold on

j=vpa(jpn(Xp)+jnp(Xn));
jnn=j-vpa(jpn(x1));
plot(x1,jnn)
hold on
jpp=j-vpa(jnp(x2));
plot(x2,jpp)
title('curreent density jpn,jpp,jnp,jnn vs x')
xlabel('device length ---> x') 
ylabel('curreent density jpn,jpp,jnp,jnn')


%%

syms vab
jd(vab)=(e*Dp*pno/Lp)+(e*Dn*npo/Ln)*((2.7182.^(vab/Vt))-1);

vd=0:.01:1;
subplot(236)
plot(vd,jd(vd))
title('current density vs v')
xlabel('Applied voltage') 
ylabel('device current density Jd')