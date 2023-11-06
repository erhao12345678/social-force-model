clear
clc
room = load('room.txt');
ex = load('exit.txt');
num = 100; %定义疏散人员数目。
r = 0.21; %定义个体半径。
vexp = normrnd(1.34, 0.37, num, 1); %定义个体期望速度。
tau = 0.5; %定义松弛时间。
t = 0.002; %定义时间步长。
T = 0.05; %定义时距阈值。
k = 1500;
kappa = 3000;
c = 25; %定义物理力中弹性力、摩擦力、粘滞力系数，按照个体质量为80kg换算。
repeat = 1; %定义重复计算次数。
t_total = zeros(repeat + 1, 1); %记录总疏散时间,前repeat项为每次用时，最后一项为平均用时。
TT = zeros(repeat, 3); %三列分别为第一个个体成功疏散的时间、最后一个个体成功疏散的时间、流量。

%=================================开始重复计算=================================
for re =10:repeat
    crowd_in_room = Occ(room, num, r); 
    
%=========================设置参数?0?3、β、δ，设定q初值==========================
    q0 = 0.5*ones(num,1);           %设定q的初值。
    epsilon = 0.99*ones(num,1);      %?0?3，修改数据：将0.01修改为所需数值。
    delta = 0.99*ones(num,1);        %δ，修改数据：将0.01修改为所需数值。
    beta_i = 0.1*ones(num,1);      %β，修改数据：将0.01修改为所需数值。
    
%=================================开始计算=====================================
    NT = [];
    time = 0;
    evac = num; %未逃出人数。
    Ncalculation = 0; %重复计算次数。
    Nrecord = 0; %绘图数据记录。
    while evac > 0
        x = crowd_in_room(:, 2);
        y = crowd_in_room(:, 3);
        vx = crowd_in_room(:, 4);
        vy = crowd_in_room(:, 5);
        v = sqrt(vx.^2 + vy.^ 2); %每个个体的坐标、分速度、总速度。
        exx = ex(4); %出口x位置。
        exy = ex(5) + ex(6) / 2; %出口中心y位置。
        d_crowd_exit = sqrt((x - exx).^2 + (y - exy).^2); %所有个体距出口距离。
        for i = 1:num
            if crowd_in_room(i,6) ~= 0
                exy1 = ex(5);
                exy2 = ex(5) + ex(6);
                dx = x(i) - exx;
                dy = y(i) - exy;
                e = zeros(num, 1); %初始化每个个体的期望速度方向角度。
                dy1 = y(i) - exy1;
                dy2 = y(i) - exy2;
                alpha = zeros(num,1);
                di = sqrt((x - x(i)).^2 + (y - y(i)).^2);
                mm = find(di < 5);
                alpha(mm) = (5 - di(mm)) / 5; %不考虑距离大于5m的人的影响。
                alpha(i) = 0; %自己无影响。
                if (abs(dy) > ex(6) / 2 - r) && (dx < 0) %此时个体不在出口对应的条带内，即不可能直接向右走离开房间。
                    e(i) = acos((-dx) / sqrt(dx^2 + dy^2));
                    if dy > 0
                        e(i) = 2 * pi - e(i);
                    end
                end
%初始化加速度。
                ax_D = (vexp(i) * cos(e(i)) - vx(i)) / tau;
                ay_D = (vexp(i) * sin(e(i)) - vy(i)) / tau; %自驱动力加速度。
                ax_S = 0;
                ay_S = 0; %初始化心理力加速度。
                ax_G = 0;
                ay_G = 0; %初始化物理力加速度。
                CandidateNoTTC = 0; %初始化预期最可能碰撞的个体编号。
                CandidateTTC = 100; %初始化一个很大的预期碰撞时间。
                ax_S_H = 0;
                ay_S_H = 0; %初始化时距心理力加速度，H代表时距心理力。
                d1 = sqrt(dx^2 + dy1^2); %个体距离上角点距离。
                d2 = sqrt(dx^2 + dy2^2); %个体距离下角点距离。
                
%=================================个体与上下墙壁相互作用==================================
                if (x(i) > 0) && (x(i) < room(4) - room(2)) %房间内部。
                    if y(i) - r < room(3) %与下墙壁发生挤压。
                        ysl = - y(i) + r + room(3);
                        ax_G = ax_G - kappa * ysl * vx(i); %摩擦力部分。
                        ay_G = ay_G + k * ysl; %弹力部分，方向向上。
                    end
                    if y(i) + r > room(5) - room(3) %与上墙壁发生挤压。
                        ysl = y(i) + r - (room(5) - room(3));
                        ax_G = ax_G - kappa * ysl * vx(i); %摩擦力部分。
                        ay_G = ay_G - k * ysl; %弹力部分，方向向下。
                    end
                    if (y(i) + r + vy(i) * T > room(5) - room(3)) || (y(i) - r + vy(i) * T < room(3)) %预期与上下墙壁发生碰撞。
                        ay_S = ay_S - vy(i) / tau;
                    end
                end
                if (x(i) > room(4) - room(2)) && (x(i) < room(4) - room(2) + ex(7)) %出口处，仅考虑物理力即可。
                    if y(i) - r < ex(5) %与出口下墙壁发生挤压。
                        ysl = -y(i) + r + ex(5);
                        ax_G = ax_G - kappa * ysl * vx(i); %摩擦力部分。
                        ay_G = ay_G - k * ysl; %弹力部分，方向向上。
                    end
                    if y(i) + r > ex(5) + ex(6) %与出口上墙壁发生挤压。
                        ysl = y(i) + r - (ex(5) + ex(6));
                        ax_G = ax_G - kappa * ysl * vx(i); %摩擦力部分。
                        ay_G = ay_G - k * ysl; %弹力部分，方向向下。
                    end
                end
                                
%=================================个体与右墙壁相互作用==================================
                if (dx + r > 0) && (abs(dy) >= ex(6) / 2) %与右墙壁发生挤压。
                    ax_G = ax_G - k * (dx + r); %弹力部分。
                    ay_G = ay_G - kappa * (dx + r) * vy(i); %摩擦力部分。
                elseif (d1 < r) && (dy1 > 0) %与出口下角点发生碰撞，不仅要求个体中心与下角点距离小于r，同时如果个体中心位于下角点之下仍然属于与墙壁碰撞。
                    ax_G = ax_G + k * (r - d1) * dx / d1;
                    ay_G = ay_G + k * (r - d1) * dy1 / d1; %弹力部分。
                    vt = (vx(i) * dy1 - vy(i) * dx) / d1;
                    vn = (vx(i) * dx + vy(i) * dy1) / d1; %速度沿个体中心与下角点连线正交分解。vn为法向，vt为切向。
                    ax_G = ax_G - kappa * vt * (r - d1) * dy1 / d1 - c * vn * dx / d1;
                    ay_G = ay_G + kappa * vt * (r - d1) * dx / d1 - c * vn * dy1 / d1; %摩擦力部分。
                    if vn < 0 %如果速度方向指向角点，则产生一个阻止继续碰撞的心理力加速度。
                        ax_S = ax_S - vn / tau * dx / d1; 
                        ay_S = ay_S - vn / tau * dy1 / d1;
                    end
                elseif (d2 < r) && (dy2 < 0) %与出口上角点发生碰撞，不仅要求个体中心与上角点距离小于r，同时如果个体中心位于上角点之上仍然属于与墙壁碰撞。
                    ax_G = ax_G + k * (r - d2) * dx / d2;
                    ay_G = ay_G + k * (r - d2) * dy2 / d2;
                    vt = (vx(i) * dy2 - vy(i) * dx) / d2;
                    vn = (vx(i) * dx + vy(i) * dy2) / d2;
                    ax_G = ax_G - kappa * vt * (r - d2) * dy2 / d2 - c * vn * dx / d2;
                    ay_G = ay_G + kappa * vt * (r - d2) * dx / d2 - c * vn * dy2 / d2;
                    if vn < 0
                        ax_S = ax_S - vn / tau * dx / d2;
                        ay_S = ay_S - vn / tau * dy2 / d2;
                    end
                elseif (vx(i) > 0) && (-(dx + r) / vx(i) < T) && ((vy(i) * ((dx + r)) / vx(i) < dy2) || (vy(i) * ((dx + r)) / vx(i) > dy1)) %预计与右墙壁发生碰撞。
                    ax_S = ax_S - vx(i) / tau;
                elseif v(i) > 1e-3 %预计与出口角点发生碰撞。
                    d1p = abs(vy(i) * dx - vx(i) * dy1) / v(i);
                    d2p = abs(vy(i) * dx - vx(i) * dy2) / v(i); %移动路线上距离上下角点的最近距离：点到直线的最短距离为垂线。
                    l_1 = sqrt(d1^2 - d1p^2) - sqrt(r^2 - d1p^2);
                    l_2 = sqrt(d2^2 - d2p^2) - sqrt(r^2 - d2p^2); %移动到碰撞点的距离。
                    if (d1p < r) && (l_1 < v(i) * T) %预计会与出口下角点发生碰撞。
                        xc = x(i) + l_1 * vx(i) / v(i);
                        yc = y(i) + l_1 * vy(i) / v(i); %碰撞点坐标。
                        xcw = xc - exx;
                        ycw = yc - exy1;
                        dcw = sqrt(xcw^2 + ycw^2); %出口下角点与碰撞点之间的距离。
                        vn = (vx(i) * xcw + vy(i) * ycw) / dcw;
                        if (vn < 0) && (xcw < 0)
                            ax_S = ax_S - vn / tau * xcw / dcw;
                            ay_S = ay_S - vn / tau * ycw / dcw;
                        end
                    end
                    if (d2p < r) && (l_2 < v(i) * T) %预计会与出口上角点发生碰撞。
                        xc = x(i) + l_2 * vx(i) / v(i);
                        yc = y(i) + l_2 * vy(i) / v(i);
                        xcw = xc - exx;
                        ycw = yc - exy2;
                        dcw = sqrt(xcw^2 + ycw^2);
                        vn = (vx(i) * xcw + vy(i) * ycw) / dcw;
                        if (vn < 0) && (xcw < 0)
                            ax_S = ax_S - vn / tau * xcw / dcw;
                            ay_S = ay_S - vn / tau * ycw / dcw;
                        end
                    end
                end
                
%=================================个体之间相互作用==================================
                xi_j = x(i) - x;
                yi_j = y(i) - y; %i、j个体间坐标差。
                vxi_j = vx(i) - vx;
                vyi_j = vy(i) - vy;
                vi_j = sqrt(vxi_j.^2 + vyi_j.^2); %i、j个体间相对速度。
                di_j = sqrt(xi_j.^2 + yi_j.^2);
                for j = 1:num
                    if (j ~= i) && (crowd_in_room(j,6) ~= 0) && (di_j(j) < 5) %考虑不包含i本身在内的、尚未成功疏散的、范围2m以内的个体。
                        if r + r - di_j(j) > 0 %i、j已经发生挤压。
                            ax_G = ax_G + k * (r + r - di_j(j)) * xi_j(j) / di_j(j);
                            ay_G = ay_G + k * (r + r - di_j(j)) * yi_j(j) / di_j(j); %弹力部分。
                            vt = (vxi_j(j) * yi_j(j) - vyi_j(j) * xi_j(j)) / di_j(j);
                            vn = (vxi_j(j) * xi_j(j) + vyi_j(j) * yi_j(j)) / di_j(j);
                            ax_G = ax_G - kappa * vt * (r + r - di_j(j)) * yi_j(j) / di_j(j) - c * vn * xi_j(j) / di_j(j);
                            ay_G = ay_G + kappa * vt * (r + r - di_j(j)) * xi_j(j) / di_j(j) - c * vn * yi_j(j) / di_j(j); %摩擦力部分。
                            if cos(e(i)) * xi_j(j) + sin(e(i)) * yi_j(j) < 0 %i的期望速度在i、j连线上有分量，使得二者有进一步挤压的趋势，此时有心理力作用。
                                if vx(i) * xi_j(j) + vy(i) * yi_j(j) < 0 %在i、j连线方向上，i的速度指向j，认为j为静止状态，触发时域加速度，使此方向上期望速度为0。
                                    ax_S_H = vexp(i) * cos(e(i)) / tau;
                                    ay_S_H = vexp(i) * sin(e(i)) / tau;
                                end
                                if vxi_j(j) * xi_j(j) + vyi_j(j) * yi_j(j) < 0 %在i、j连线方向上，i、j相对速度指向j，认为j在移动，触发TTC加速度，使得在以j为参考系中i的期望速度为0。
                                    ax_S = ax_S - ((vexp(i) * cos(e(i)) - vx(j)) * xi_j(j) + (vexp(i) * sin(e(i)) - vy(j)) * yi_j(j)) / (di_j(j) * tau) * xi_j(j) / di_j(j);
                                    ay_S = ay_S - ((vexp(i) * cos(e(i)) - vx(j)) * xi_j(j) + (vexp(i) * sin(e(i)) - vy(j)) * yi_j(j)) / (di_j(j) * tau) * yi_j(j) / di_j(j);
                                end
                            end
                        else
                            if v(i) > 1e-3 %条件1。
                                dp1 = abs(vy(i) * xi_j(j) - vx(i) * yi_j(j)) / v(i);
                                if (vx(i) * xi_j(j) + vy(i) * yi_j(j) < 0) && (dp1 < r + r) && (sqrt(di_j(j)^2 - dp1^2) - sqrt((r + r)^2 - dp1^2) < v(i) * T) %i、j尚未发生碰撞，但i的速度指向j，预计有可能碰撞且碰撞时间小于阈值，触发时域加速度。
                                    ax_S_H = vexp(i) * cos(e(i)) / tau;
                                    ay_S_H = vexp(i) * sin(e(i)) / tau;
                                end
                            end
                            if vi_j(j) > 1e-3 %条件2。
                                dp1 = abs(vyi_j(j) * xi_j(j) - vxi_j(j) * yi_j(j)) / vi_j(j);
                                if (vxi_j(j) * xi_j(j) + vyi_j(j) * yi_j(j) < 0) && (dp1 < r + r) %i、j尚未发生碰撞，但i、j相对速度指向j，预计有可能碰撞。
                                    ttc = (sqrt(di_j(j)^2 - dp1^2) - sqrt((r + r)^2 - dp1^2)) / vi_j(j);
                                    if (ttc < T) && (ttc < CandidateTTC) %筛选最可能发生碰撞的一个个体，触发TTC加速度。
                                        CandidateTTC = ttc;
                                        CandidateNoTTC = j;
                                    end
                                end
                            end
                        end
                    end
                end
                if CandidateNoTTC ~= 0
                    j = CandidateNoTTC;
                    dp = abs(vyi_j(j) * (xi_j(j)) - vxi_j(j) * (yi_j(j))) / vi_j(j);
                    l0 = sqrt(di_j(j)^2 - dp^2) - sqrt((r + r)^2 - dp^2);
                    xc = x(i) + l0 * vx(i) / v(i);
                    yc = y(i) + l0 * vy(i) / v(i);
                    xcj = xc - x(j);
                    ycj = yc - y(j);
                    dcj = sqrt(xcj^2 + ycj^2);
                    if cos(e(i)) * xcj + sin(e(i)) * ycj < 0 %因为i、j尚未发生碰撞，所以TTC加速度使碰撞时的相对法向速度为0。
                        ax_S = ax_S - (vxi_j(j) * xcj + vyi_j(j) * ycj) / (dcj * tau) * xcj / dcj;
                        ay_S = ay_S - (vxi_j(j) * xcj + vyi_j(j) * ycj) / (dcj * tau) * ycj / dcj;
                    end
                end
                %如果j满足条件1，则会产生时域加速度；如果j满足条件2（是与i预计碰撞时间最小的即最容易与i发生碰撞的），则会产生TTC加速度。j有可能全部满足，则会产生两种加速度。
                %只要其他个体满足条件1，就会产生时域加速度；只有最满足条件2的个体，才会产生TTC加速度。
                
%=================================更新个体速度、位置==================================
                ax_S = ax_S - ax_S_H;
                ay_S = ay_S - ay_S_H;
                tt = find(crowd_in_room(:, 6) ~= 0); 
                ax_S = (1 - q0(i)) * ax_S;
                ay_S = (1 - q0(i)) * ay_S;
                ax = ax_S + ax_G + ax_D;
                ay = ay_S + ay_G + ay_D;
                if d_crowd_exit(i) == min(d_crowd_exit(tt)) %距离出口最近的个体认为其按照期望速度移动。
                   ax = ax_D;
                   ay = ay_D;
                end
                vx_new = vx(i) + ax * t;
                vy_new = vy(i) + ay * t;
                tmp = sqrt(vx_new^2 + vy_new^2);
                if tmp >= vexp(i)
                    vx_new = vx_new * vexp(i) / tmp;
                    vy_new = vy_new * vexp(i) / tmp;
                end
                crowd_in_room(i,2) = x(i) + (vx_new + vx(i)) * t / 2;
                crowd_in_room(i,3) = y(i) + (vy_new + vy(i)) * t / 2;
                crowd_in_room(i,4) = vx_new;
                crowd_in_room(i,5) = vy_new;
                boun = 2; %认为个体在出口通道内及以外2m内都会对其他个体的疏散产生影响。
                if (crowd_in_room(i,2) > ex(4) + ex(7) + boun)
                    crowd_in_room(i,6) = 0;
                    crowd_in_room(i, 2) = 100;
					crowd_in_room(i, 3) = 100; %将成功疏散的个体坐标放到远处。
                    evac = evac - 1;
                    NT = [NT;time num - evac];
                end
                if (ex(4) <= crowd_in_room(i,2)) && (crowd_in_room(i,2) <= ex(4) + ex(7) + boun) && ((crowd_in_room(i,3) <= ex(5)) || (crowd_in_room(i,3) >= ex(5) + ex(6)))
                    crowd_in_room(i,6) = 0;
                    crowd_in_room(i, 2) = 100;
					crowd_in_room(i, 3) = 100; 
                    evac = evac - 1;
                    NT = [NT;time num - evac];
                end
                if sum(crowd_in_room(:, 6)) == num - 1 %记录第一个个体疏散成功的时间。
					TT(re, 1) = NT(1,1);
                end
                if sum(crowd_in_room(:, 6)) == 0 %记录最后一个个体疏散成功的时间。
					TT(re, 2) = NT(num,1);
                end
            end
            cal_i = epsilon .* alpha;
            if sum(cal_i) ~= 0 %若是等于0，则该个体不受其他个体的影响。
                gamma_i = sum(cal_i .* delta);
                q_star = sum(cal_i .* q0) / sum(cal_i);
                nervous_i = 1 - (1 - q_star) * (1 - q0(i));
                calm_i = q0(i) * q_star;
                dq = gamma_i * (beta_i(i) * nervous_i + (1 - beta_i(i)) * calm_i) * t;
                q0(i) = q0(i) + dq;
                if q0(i) >= 1 %q不能超过1。
                    q0(i) = 0.99;
                end
            end
        end
        
%===================采集绘图数据，记录0s、25s、50s……时刻计算数据================        
        if mod(Ncalculation,1250) == 0
            Nrecord = Nrecord + 1;
            MT(1,:,Nrecord) = x;
            MT(2,:,Nrecord) = y;
        end
        time = time + t;
        Ncalculation = Ncalculation + 1;
    end
    TT(re, 3) = (num - 1) / (TT(re, 2) - TT(re, 1));
    t_total(re,1) = TT(re, 2);
end
t_total(repeat + 1,1) = mean(t_total(1:repeat,1));
disp('Evacuation Finished!')

%=============================运动轨迹绘图=================================
for picture = 1:Nrecord
    figure(picture)
%墙壁轮廓。
    X = 0:0.001:room(4);
    sizeX = size(X);
    Y = zeros(sizeX(2),1);
    plot(X,Y,'-k','linewidth',2)
    hold on
    X = 0:0.001:room(4);
    sizeX = size(X);
    Y = room(5) * ones(sizeX(2),1);
    plot(X,Y,'-k','linewidth',2)
    hold on
    Y = 0:0.001:ex(5);
    sizeY = size(Y);
    X = room(4) * ones(sizeY(2),1);
    plot(X,Y,'-k','linewidth',2)
    hold on
    Y = (ex(5) + ex(6)):0.001:room(5);
    sizeY = size(Y);
    X = room(4) * ones(sizeY(2),1);
    plot(X,Y,'-k','linewidth',2)
    hold on
    Y = 0:0.001:room(5);
    sizeY = size(Y);
    X = zeros(sizeY(2),1);
    plot(X,Y,'-k','linewidth',2)
    hold on
    axis([-1 (ex(4) + ex(7) + boun) -1 (room(5) + 1)])
    hold on
%个体位置。
    circleangle = 0:(pi / 200):(2 * pi);
    for people = 1:num
        plot((MT(1,people,picture) + r * cos(circleangle)), (MT(2,people,picture) + r * sin(circleangle)), '-r')
        hold on
    end
end