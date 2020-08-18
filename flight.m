%% 直升机飞行性能计算
%%%%%%%%%%%%%%%%数据输入%%%%%%%%%%%%%%
clear
% clc
%% 直升机参数
M=385;
g=9.8;
G=M*g;%直升机重量
%%%%此处可以更改%%%%%%%%%%
delta=0.2;%旋翼轴前倾角
%%%%%%%%%%%%%%%%%%%%%%%%%
rou=1.225;%空气密度
kapa=0.92;%叶端损失系数
%% 主旋翼参数
R=3.2;%主旋翼半径
A=pi*R^2;%桨盘面积
k=2;%主旋翼桨叶片数
omiga=60;%主旋翼转速
b=0.22;%主旋翼弦长
theta1=-8*pi/180;%主旋翼扭转角
sigma=k*b/pi/R;%主旋翼实度
a8=0.1337/0.1*5.731;%主旋翼升力线斜率
Cx=0.008;%主旋翼阻力系数
Ib= 38.52;%质量惯矩
Mb=18.8;%桨叶质量矩
gamab=a8*b*rou*R^4/Ib;%桨叶的洛克数
e= 0;% 挥舞铰
p=G/pi/R/R;%桨盘载荷
Kchui=p/1000;%K 垂直，一个系数
%% 尾桨参数
Rt=0.64;%尾桨半径
At=pi*Rt^2;
kt=2;%尾桨桨叶片数
bt=0.12;%尾桨弦长
omigat=312;%尾桨转速
theta1t=-4*pi/180;%尾桨扭转角
sigmat=kt*bt/pi/Rt;%尾桨实度
a8t=0.1226/0.1*5.731;%尾桨升力线斜率
Cxt=0.034;%尾桨阻力系数
gamat=0.0001;
%% 重心位置
lambda1=-0.1;lambda1t=-0.1;vx=0;vy=0;vz=0;
% syms vx vy vz lambda1 lambda1t
%旋翼、尾桨减去直升机重心位置的相对坐标
%位置 1 重心(3.12,0.6,0)
xm=0.05;ym=1.0;zm=0;     %桨毂中心坐标(3.17,1.6,0)
xt=4.08;yt=0.24;zt=0;    %尾桨中心坐标(7.2,0.84,0)
%位置 2 重心(3.17,0.6,0)
% xm=0;ym=1.0;zm=0;
% xt=4.03;yt=0.24;zt=0;
%%位置 3 重心(3.22,0.6,0)
% xm=-0.05;ym=1.0;zm=0;
% xt=3.98;yt=0.24;zt=0;
%% 旋翼气动力
syms fai7 fait A1 B1 faiT theta gama
x=[fai7 A1 B1 fait theta gama];%配平量
theta0=fai7-0.7*theta1;
vxm = vx*cos(delta)-vy*sin(delta);
vym = vx*sin(delta)+vy*cos(delta);
vzm = vz;
V = sqrt(vxm^2+vzm^2);% 水平速度[m/s]
v1 = sqrt(-0.5*V^2+sqrt((0.5*V^2)^2+(vym/(2*rou*A))^2)); % 入流速度[m/s]
mu = V/(omiga*R);% 前进比 [less]
muX = vxm/(omiga*R);% x 轴前进比 [less]
muY = vym/(omiga*R);% y 轴前进比 [less]
muZ = vzm/(omiga*R);
lambda0 = muY;
%旋翼挥舞运动系数
a0 =0.5*kapa*gamab*(1/4*(1+mu^2)*theta0+(1/5+1/6*mu^2)*theta1+1/3*lambda1-1/3*mu*B1)-g*Mb/(omiga^2*Ib);%旋翼锥度角a0
a1s =(8/3*theta0+2*theta1+2*lambda1)*mu/(1-1/2*mu^2)-B1*(1+2*mu^2/(1-1/2*mu^2));%后倒角a1s
b1s = 4/3*mu*a0/(1+1/2*mu^2)+A1;                                                %侧倒角b1s
CT =kapa*a8*sigma*((1/3+1/2*mu^2)*theta0+1/4*(1+mu^2)*theta1+1/2*lambda1-1/2*mu*B1);%拉力系数
CH =1/2*sigma*Cx*mu+1/3*kapa*a8*sigma*(a1s*(theta0+3/4*theta1+9/4*lambda1)...       %后向力系数
-1/2*a0*b1s+1/4*mu*(a0^2+a1s^2-lambda1*(2*theta0+theta1)))...
-1/3*kapa*a8*sigma*(3/4*(mu*a1s-lambda1)*B1-1/2*a0*A1);

CS =1/3*kapa*a8*sigma*(b1s*(theta0+3/4*theta1)*(1+3/2*mu^2)+b1s*9/4*(1/3*mu*a1s+lambda1)...%侧向力系数
-9/4*a0*mu*((theta0+3/4*theta1)-2*lambda1)+1/2*a0*a1s*(1-6*mu^2))...
-kapa*a8*sigma*(1/4*(mu*a1s+lambda1)*A1-1/6*(a0*(1+3*mu^2)-3*mu*b1s)*B1);

mk = 1/4*sigma*Cx*(1+mu^2)+kapa*a8*sigma*(theta0*(-1/3*lambda1)+theta1*(-1/4*lambda1)...  %反扭矩系数
-A1*(1/6*muX*a0-1/8*b1s-1/4*muZ*lambda1-1/8*muX*muZ*a1s+1/16*(muZ^2-muX^2)*b1s)...
-B1*(-1/6*muZ*a0-1/8*(a1s)-1/4*muX*lambda1+1/8*muX*muZ*b1s+1/16*(muZ^2-muX^2)*a1s)...
-1/4*(mu^2)*a0^2-1/8*(1+3/2*muX^2+1/2*muZ^2)*a1s^2-1/8*b1s^2*(1+1/2*muX^2+3/2*muZ^2)...
+1/3*(muZ*a1s+muX*b1s)*a0+1/4*muX*muZ*a1s*b1s-1/2*lambda1*(muX*a1s-muZ*b1s+lambda1));

T = 1/2*CT*rou*A*(omiga*R)^2;% 拉力 [N]
H = 1/2*CH*rou*A*(omiga*R)^2;% 后向力
S = 1/2*CS*rou*A*(omiga*R)^2;% 侧向力
Mk = 1/2*mk*rou*A*omiga^2*R^3;%反扭矩[Nm]
%桨毂力矩
MGx = k/2*Mb*omiga^2*e*b1s;%滚转力矩
MGz = k/2*Mb*omiga^2*e*a1s;%俯仰力矩
%%旋翼空气动力在体轴系的力和力矩分量
Fxm = T*sin(delta)-H*cos(delta);
Fym = T*cos(delta)+H*sin(delta);
Fzm = S;
Mxm = MGx*cos(delta)-Mk*sin(delta)+Fzm*ym-Fym*zm;
Mym =-MGx*sin(delta)-Mk*cos(delta)+Fxm*zm+Fzm*xm;
Mzm = MGz-Fym*xm-Fxm*ym;
%% 尾桨空气动力
theta0t=fait-0.7*theta1t;
Vt = sqrt(vy^2+vx^2);% 水平速度[m/s]
v1t = sqrt(-0.5*Vt^2+sqrt((0.5*Vt^2)^2+(vy/(2*rou*At))^2));% 入流速度[m/s]
lambda_pr =(v1t+vx)/(omigat*Rt);% 入流比 [less]
mut = Vt/(omigat*Rt);% 前进比 [less]
muXt = vx/(omigat*Rt);% x 轴前进比 [less]
muYt = vy/(omigat*Rt);
muZt = vz/(omigat*Rt);% y 轴前进比 [less]
lambda0t = muZt;
%挥舞运动系数
a0t =0.5*kapa*gamat*(1/4*(1+mut^2)*theta0t+(1/5+1/6*mut^2)*theta1t+1/3*lambda1t-1/3*mut*0);%挥舞锥度角a0t
a1st =(8/3*theta0t+2*theta1t+2*lambda1t)*mut/(1-1/2*mut^2)-0*(1+2*mut^2/(1-1/2*mut^2));    %后倒角a1st
b1st = 4/3*mut*a0t/(1+1/2*mut^2)+0;                                                        %侧倒角b1st
CTt =kapa*a8t*sigmat*((1/3+1/2*mut^2)*theta0t+1/4*(1+mut^2)*theta1t+1/2*lambda1t-1/2*mut*0);%拉力系数
CHt =1/2*sigmat*Cx*mut+1/3*kapa*a8t*sigmat*(a1st*(theta0t+3/4*theta1t+9/4*lambda1t)...     %后向力系数
-1/2*a0t*b1st+1/4*mut*(a0t^2+a1st^2-lambda1t*(2*theta0t+theta1t)))...
-1/3*kapa*a8t*sigmat*(3/4*(mut*a1st-lambda1t)*0-1/2*a0t*0);
CSt =1/3*kapa*a8t*sigmat*(b1st*(theta0t+3/4*theta1t)*(1+3/2*mut^2)+b1st*9/4*(1/3*mut*a1st+lambda1t)...%侧向力系数
-9/4*a0t*mut*((theta0t+3/4*theta1t)-2*lambda1t)+1/2*a0t*a1st*(1-6*mut^2))...
-kapa*a8t*sigmat*(1/4*(mut*a1st+lambda1t)*0-1/6*(a0t*(1+3*mut^2)-3*mut*b1st)*0);
A1t = 0;B1t = 0;
mkt = 1/4*sigmat*Cx*(1+mut^2)-kapa*a8t*sigmat*(theta0t*(1/3*lambda1t)+theta1t*(-1/4*lambda1t)... %反扭矩系数
-A1t*(1/6*muXt*a0t-1/8*b1st-1/4*muZt*lambda1t-1/8*muXt*muZt*a1st+1/16*(muZt^2-muXt^2)*b1st)...
-B1t*(-1/6*muZt*a0t-1/8*(a1st)-1/4*muXt*lambda1t+1/8*muXt*muZt*b1st+1/16*(muZt^2-muXt^2)*a1st)...
-1/4*(mut^2)*a0t^2-1/8*(1+3/2*muXt^2+1/2*muZt^2)*a1st^2-1/8*b1st^2*(1+1/2*muXt^2+3/2*muZt^2)...
+1/3*(muZt*a1st+muXt*b1st)*a0t+1/4*muXt*muZt*a1st*b1st-1/2*lambda1t*(muXt*a1st-muZt*b1st+lambda1t));

Tt = 1/2*CTt*rou*At*(omigat*Rt)^2;% 拉力 [N]
Ht = 1/2*CHt*rou*At*(omigat*Rt)^2;% 后向力
St = 1/2*CSt*rou*At*(omigat*Rt)^2;% 侧向力
Mkt = 1/2*mkt*rou*At*omigat^2*Rt^3;%反扭矩[Nm]
%%尾桨拉力在体轴系的力与力矩分量
Fxt =-Ht;
Fyt =-St;
Fzt = Tt;
Mzt =-Mkt-Fxt*yt-Fyt*xt;
Mxt = Fzt*yt;
Myt = Fzt*xt;

syms theta gama
Fx = Fxm+Fxt-G*sin(-theta);
Fz =-(Fym+Fyt-G*cos(-theta)*cos(gama));
Fy = Fzm+Fzt+G*cos(theta)*sin(gama);
Mx = Mxm+Mxt;
Mz =-(Mym+Myt);
My = Mzm+Mzt;
%% 配平计算
%直升机悬停时的简化平衡方程
f1=T*(delta-a1s)-G*theta;
f2=T-G*(1+Kchui);
f3=T*b1s+Tt+G*gama;
f4=T*b1s*ym+Tt*yt;
f5=-Mk+Tt*(-xt)+T*b1s*xm;
f6=T*(xm-ym*delta)+T*a1s*ym;
%% 牛顿法求解
f =[f1; f2; f3; f4; f5; f6];
df = jacobian(f, x);%x=[fai7 A1 B1 fait theta gama];
X0=[(7/180)*pi (5/180)*pi (-1/180)*pi (16/180)*pi (1.7/180)*pi (9/180)*pi];
% X0=[(fai7)(A1)(B1)(fait)(theta)(gama)];
F=double(subs(f,{fai7, A1, B1, fait, theta, gama},{X0}));%subs()是符号计算函数
DF=double(subs(df,{fai7, A1, B1, fait, theta, gama},{X0}));%df
X1=transpose(X0)-DF^(-1)*F;%X0-f/df
for j=1:500
    n= norm(X1-transpose(X0));%二范数 （返回向量X1-X0的p范数）
    if n< 1e-60
    break;
    end
    X0=transpose(X1);
    F=double(subs(f,{fai7, A1, B1, fait, theta, gama},{X0}));
    DF=double(subs(df,{fai7, A1, B1, fait, theta, gama},{X0}));
end
disp(['配平量:' '旋翼总距、' '横向周期变距、' '纵向周期变距、' '尾桨总距、' '俯仰角、' '侧倾角γ'])
disp(['配平结果:']) % x =[fai7 A1 B1 fait theta gama]
X=(X1*180/pi)'
disp(['迭代次数为:' num2str(j) '次'])
dfx=diff(Fxm,vxm);
% %% 稳定性分析
% % matrix A
% element = [m m m I_x I_y I_z 1 1 1];
% matrixA = diag(element);
% matrixA(5,4) = -I_xy;matrixA(4,5) = -I_xy;
% matrixA(5,4)
% % matrix B


