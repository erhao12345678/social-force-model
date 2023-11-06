clear
clc
room = load('room.txt');
ex = load('exit.txt');
num = 100; %������ɢ��Ա��Ŀ��
r = 0.21; %�������뾶��
vexp = normrnd(1.34, 0.37, num, 1); %������������ٶȡ�
tau = 0.5; %�����ɳ�ʱ�䡣
t = 0.002; %����ʱ�䲽����
T = 0.05; %����ʱ����ֵ��
k = 1500;
kappa = 3000;
c = 25; %�����������е�������Ħ������ճ����ϵ�������ո�������Ϊ80kg���㡣
repeat = 1; %�����ظ����������
t_total = zeros(repeat + 1, 1); %��¼����ɢʱ��,ǰrepeat��Ϊÿ����ʱ�����һ��Ϊƽ����ʱ��
TT = zeros(repeat, 3); %���зֱ�Ϊ��һ������ɹ���ɢ��ʱ�䡢���һ������ɹ���ɢ��ʱ�䡢������

%=================================��ʼ�ظ�����=================================
for re =10:repeat
    crowd_in_room = Occ(room, num, r); 
    
%=========================���ò���?0?3���¡��ģ��趨q��ֵ==========================
    q0 = 0.5*ones(num,1);           %�趨q�ĳ�ֵ��
    epsilon = 0.99*ones(num,1);      %?0?3���޸����ݣ���0.01�޸�Ϊ������ֵ��
    delta = 0.99*ones(num,1);        %�ģ��޸����ݣ���0.01�޸�Ϊ������ֵ��
    beta_i = 0.1*ones(num,1);      %�£��޸����ݣ���0.01�޸�Ϊ������ֵ��
    
%=================================��ʼ����=====================================
    NT = [];
    time = 0;
    evac = num; %δ�ӳ�������
    Ncalculation = 0; %�ظ����������
    Nrecord = 0; %��ͼ���ݼ�¼��
    while evac > 0
        x = crowd_in_room(:, 2);
        y = crowd_in_room(:, 3);
        vx = crowd_in_room(:, 4);
        vy = crowd_in_room(:, 5);
        v = sqrt(vx.^2 + vy.^ 2); %ÿ����������ꡢ���ٶȡ����ٶȡ�
        exx = ex(4); %����xλ�á�
        exy = ex(5) + ex(6) / 2; %��������yλ�á�
        d_crowd_exit = sqrt((x - exx).^2 + (y - exy).^2); %���и������ھ��롣
        for i = 1:num
            if crowd_in_room(i,6) ~= 0
                exy1 = ex(5);
                exy2 = ex(5) + ex(6);
                dx = x(i) - exx;
                dy = y(i) - exy;
                e = zeros(num, 1); %��ʼ��ÿ������������ٶȷ���Ƕȡ�
                dy1 = y(i) - exy1;
                dy2 = y(i) - exy2;
                alpha = zeros(num,1);
                di = sqrt((x - x(i)).^2 + (y - y(i)).^2);
                mm = find(di < 5);
                alpha(mm) = (5 - di(mm)) / 5; %�����Ǿ������5m���˵�Ӱ�졣
                alpha(i) = 0; %�Լ���Ӱ�졣
                if (abs(dy) > ex(6) / 2 - r) && (dx < 0) %��ʱ���岻�ڳ��ڶ�Ӧ�������ڣ���������ֱ���������뿪���䡣
                    e(i) = acos((-dx) / sqrt(dx^2 + dy^2));
                    if dy > 0
                        e(i) = 2 * pi - e(i);
                    end
                end
%��ʼ�����ٶȡ�
                ax_D = (vexp(i) * cos(e(i)) - vx(i)) / tau;
                ay_D = (vexp(i) * sin(e(i)) - vy(i)) / tau; %�����������ٶȡ�
                ax_S = 0;
                ay_S = 0; %��ʼ�����������ٶȡ�
                ax_G = 0;
                ay_G = 0; %��ʼ�����������ٶȡ�
                CandidateNoTTC = 0; %��ʼ��Ԥ���������ײ�ĸ����š�
                CandidateTTC = 100; %��ʼ��һ���ܴ��Ԥ����ײʱ�䡣
                ax_S_H = 0;
                ay_S_H = 0; %��ʼ��ʱ�����������ٶȣ�H����ʱ����������
                d1 = sqrt(dx^2 + dy1^2); %��������Ͻǵ���롣
                d2 = sqrt(dx^2 + dy2^2); %��������½ǵ���롣
                
%=================================����������ǽ���໥����==================================
                if (x(i) > 0) && (x(i) < room(4) - room(2)) %�����ڲ���
                    if y(i) - r < room(3) %����ǽ�ڷ�����ѹ��
                        ysl = - y(i) + r + room(3);
                        ax_G = ax_G - kappa * ysl * vx(i); %Ħ�������֡�
                        ay_G = ay_G + k * ysl; %�������֣��������ϡ�
                    end
                    if y(i) + r > room(5) - room(3) %����ǽ�ڷ�����ѹ��
                        ysl = y(i) + r - (room(5) - room(3));
                        ax_G = ax_G - kappa * ysl * vx(i); %Ħ�������֡�
                        ay_G = ay_G - k * ysl; %�������֣��������¡�
                    end
                    if (y(i) + r + vy(i) * T > room(5) - room(3)) || (y(i) - r + vy(i) * T < room(3)) %Ԥ��������ǽ�ڷ�����ײ��
                        ay_S = ay_S - vy(i) / tau;
                    end
                end
                if (x(i) > room(4) - room(2)) && (x(i) < room(4) - room(2) + ex(7)) %���ڴ������������������ɡ�
                    if y(i) - r < ex(5) %�������ǽ�ڷ�����ѹ��
                        ysl = -y(i) + r + ex(5);
                        ax_G = ax_G - kappa * ysl * vx(i); %Ħ�������֡�
                        ay_G = ay_G - k * ysl; %�������֣��������ϡ�
                    end
                    if y(i) + r > ex(5) + ex(6) %�������ǽ�ڷ�����ѹ��
                        ysl = y(i) + r - (ex(5) + ex(6));
                        ax_G = ax_G - kappa * ysl * vx(i); %Ħ�������֡�
                        ay_G = ay_G - k * ysl; %�������֣��������¡�
                    end
                end
                                
%=================================��������ǽ���໥����==================================
                if (dx + r > 0) && (abs(dy) >= ex(6) / 2) %����ǽ�ڷ�����ѹ��
                    ax_G = ax_G - k * (dx + r); %�������֡�
                    ay_G = ay_G - kappa * (dx + r) * vy(i); %Ħ�������֡�
                elseif (d1 < r) && (dy1 > 0) %������½ǵ㷢����ײ������Ҫ������������½ǵ����С��r��ͬʱ�����������λ���½ǵ�֮����Ȼ������ǽ����ײ��
                    ax_G = ax_G + k * (r - d1) * dx / d1;
                    ay_G = ay_G + k * (r - d1) * dy1 / d1; %�������֡�
                    vt = (vx(i) * dy1 - vy(i) * dx) / d1;
                    vn = (vx(i) * dx + vy(i) * dy1) / d1; %�ٶ��ظ����������½ǵ����������ֽ⡣vnΪ����vtΪ����
                    ax_G = ax_G - kappa * vt * (r - d1) * dy1 / d1 - c * vn * dx / d1;
                    ay_G = ay_G + kappa * vt * (r - d1) * dx / d1 - c * vn * dy1 / d1; %Ħ�������֡�
                    if vn < 0 %����ٶȷ���ָ��ǵ㣬�����һ����ֹ������ײ�����������ٶȡ�
                        ax_S = ax_S - vn / tau * dx / d1; 
                        ay_S = ay_S - vn / tau * dy1 / d1;
                    end
                elseif (d2 < r) && (dy2 < 0) %������Ͻǵ㷢����ײ������Ҫ������������Ͻǵ����С��r��ͬʱ�����������λ���Ͻǵ�֮����Ȼ������ǽ����ײ��
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
                elseif (vx(i) > 0) && (-(dx + r) / vx(i) < T) && ((vy(i) * ((dx + r)) / vx(i) < dy2) || (vy(i) * ((dx + r)) / vx(i) > dy1)) %Ԥ������ǽ�ڷ�����ײ��
                    ax_S = ax_S - vx(i) / tau;
                elseif v(i) > 1e-3 %Ԥ������ڽǵ㷢����ײ��
                    d1p = abs(vy(i) * dx - vx(i) * dy1) / v(i);
                    d2p = abs(vy(i) * dx - vx(i) * dy2) / v(i); %�ƶ�·���Ͼ������½ǵ��������룺�㵽ֱ�ߵ���̾���Ϊ���ߡ�
                    l_1 = sqrt(d1^2 - d1p^2) - sqrt(r^2 - d1p^2);
                    l_2 = sqrt(d2^2 - d2p^2) - sqrt(r^2 - d2p^2); %�ƶ�����ײ��ľ��롣
                    if (d1p < r) && (l_1 < v(i) * T) %Ԥ�ƻ�������½ǵ㷢����ײ��
                        xc = x(i) + l_1 * vx(i) / v(i);
                        yc = y(i) + l_1 * vy(i) / v(i); %��ײ�����ꡣ
                        xcw = xc - exx;
                        ycw = yc - exy1;
                        dcw = sqrt(xcw^2 + ycw^2); %�����½ǵ�����ײ��֮��ľ��롣
                        vn = (vx(i) * xcw + vy(i) * ycw) / dcw;
                        if (vn < 0) && (xcw < 0)
                            ax_S = ax_S - vn / tau * xcw / dcw;
                            ay_S = ay_S - vn / tau * ycw / dcw;
                        end
                    end
                    if (d2p < r) && (l_2 < v(i) * T) %Ԥ�ƻ�������Ͻǵ㷢����ײ��
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
                
%=================================����֮���໥����==================================
                xi_j = x(i) - x;
                yi_j = y(i) - y; %i��j���������
                vxi_j = vx(i) - vx;
                vyi_j = vy(i) - vy;
                vi_j = sqrt(vxi_j.^2 + vyi_j.^2); %i��j���������ٶȡ�
                di_j = sqrt(xi_j.^2 + yi_j.^2);
                for j = 1:num
                    if (j ~= i) && (crowd_in_room(j,6) ~= 0) && (di_j(j) < 5) %���ǲ�����i�������ڵġ���δ�ɹ���ɢ�ġ���Χ2m���ڵĸ��塣
                        if r + r - di_j(j) > 0 %i��j�Ѿ�������ѹ��
                            ax_G = ax_G + k * (r + r - di_j(j)) * xi_j(j) / di_j(j);
                            ay_G = ay_G + k * (r + r - di_j(j)) * yi_j(j) / di_j(j); %�������֡�
                            vt = (vxi_j(j) * yi_j(j) - vyi_j(j) * xi_j(j)) / di_j(j);
                            vn = (vxi_j(j) * xi_j(j) + vyi_j(j) * yi_j(j)) / di_j(j);
                            ax_G = ax_G - kappa * vt * (r + r - di_j(j)) * yi_j(j) / di_j(j) - c * vn * xi_j(j) / di_j(j);
                            ay_G = ay_G + kappa * vt * (r + r - di_j(j)) * xi_j(j) / di_j(j) - c * vn * yi_j(j) / di_j(j); %Ħ�������֡�
                            if cos(e(i)) * xi_j(j) + sin(e(i)) * yi_j(j) < 0 %i�������ٶ���i��j�������з�����ʹ�ö����н�һ����ѹ�����ƣ���ʱ�����������á�
                                if vx(i) * xi_j(j) + vy(i) * yi_j(j) < 0 %��i��j���߷����ϣ�i���ٶ�ָ��j����ΪjΪ��ֹ״̬������ʱ����ٶȣ�ʹ�˷����������ٶ�Ϊ0��
                                    ax_S_H = vexp(i) * cos(e(i)) / tau;
                                    ay_S_H = vexp(i) * sin(e(i)) / tau;
                                end
                                if vxi_j(j) * xi_j(j) + vyi_j(j) * yi_j(j) < 0 %��i��j���߷����ϣ�i��j����ٶ�ָ��j����Ϊj���ƶ�������TTC���ٶȣ�ʹ������jΪ�ο�ϵ��i�������ٶ�Ϊ0��
                                    ax_S = ax_S - ((vexp(i) * cos(e(i)) - vx(j)) * xi_j(j) + (vexp(i) * sin(e(i)) - vy(j)) * yi_j(j)) / (di_j(j) * tau) * xi_j(j) / di_j(j);
                                    ay_S = ay_S - ((vexp(i) * cos(e(i)) - vx(j)) * xi_j(j) + (vexp(i) * sin(e(i)) - vy(j)) * yi_j(j)) / (di_j(j) * tau) * yi_j(j) / di_j(j);
                                end
                            end
                        else
                            if v(i) > 1e-3 %����1��
                                dp1 = abs(vy(i) * xi_j(j) - vx(i) * yi_j(j)) / v(i);
                                if (vx(i) * xi_j(j) + vy(i) * yi_j(j) < 0) && (dp1 < r + r) && (sqrt(di_j(j)^2 - dp1^2) - sqrt((r + r)^2 - dp1^2) < v(i) * T) %i��j��δ������ײ����i���ٶ�ָ��j��Ԥ���п�����ײ����ײʱ��С����ֵ������ʱ����ٶȡ�
                                    ax_S_H = vexp(i) * cos(e(i)) / tau;
                                    ay_S_H = vexp(i) * sin(e(i)) / tau;
                                end
                            end
                            if vi_j(j) > 1e-3 %����2��
                                dp1 = abs(vyi_j(j) * xi_j(j) - vxi_j(j) * yi_j(j)) / vi_j(j);
                                if (vxi_j(j) * xi_j(j) + vyi_j(j) * yi_j(j) < 0) && (dp1 < r + r) %i��j��δ������ײ����i��j����ٶ�ָ��j��Ԥ���п�����ײ��
                                    ttc = (sqrt(di_j(j)^2 - dp1^2) - sqrt((r + r)^2 - dp1^2)) / vi_j(j);
                                    if (ttc < T) && (ttc < CandidateTTC) %ɸѡ����ܷ�����ײ��һ�����壬����TTC���ٶȡ�
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
                    if cos(e(i)) * xcj + sin(e(i)) * ycj < 0 %��Ϊi��j��δ������ײ������TTC���ٶ�ʹ��ײʱ����Է����ٶ�Ϊ0��
                        ax_S = ax_S - (vxi_j(j) * xcj + vyi_j(j) * ycj) / (dcj * tau) * xcj / dcj;
                        ay_S = ay_S - (vxi_j(j) * xcj + vyi_j(j) * ycj) / (dcj * tau) * ycj / dcj;
                    end
                end
                %���j��������1��������ʱ����ٶȣ����j��������2������iԤ����ײʱ����С�ļ���������i������ײ�ģ���������TTC���ٶȡ�j�п���ȫ�����㣬���������ּ��ٶȡ�
                %ֻҪ����������������1���ͻ����ʱ����ٶȣ�ֻ������������2�ĸ��壬�Ż����TTC���ٶȡ�
                
%=================================���¸����ٶȡ�λ��==================================
                ax_S = ax_S - ax_S_H;
                ay_S = ay_S - ay_S_H;
                tt = find(crowd_in_room(:, 6) ~= 0); 
                ax_S = (1 - q0(i)) * ax_S;
                ay_S = (1 - q0(i)) * ay_S;
                ax = ax_S + ax_G + ax_D;
                ay = ay_S + ay_G + ay_D;
                if d_crowd_exit(i) == min(d_crowd_exit(tt)) %�����������ĸ�����Ϊ�䰴�������ٶ��ƶ���
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
                boun = 2; %��Ϊ�����ڳ���ͨ���ڼ�����2m�ڶ���������������ɢ����Ӱ�졣
                if (crowd_in_room(i,2) > ex(4) + ex(7) + boun)
                    crowd_in_room(i,6) = 0;
                    crowd_in_room(i, 2) = 100;
					crowd_in_room(i, 3) = 100; %���ɹ���ɢ�ĸ�������ŵ�Զ����
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
                if sum(crowd_in_room(:, 6)) == num - 1 %��¼��һ��������ɢ�ɹ���ʱ�䡣
					TT(re, 1) = NT(1,1);
                end
                if sum(crowd_in_room(:, 6)) == 0 %��¼���һ��������ɢ�ɹ���ʱ�䡣
					TT(re, 2) = NT(num,1);
                end
            end
            cal_i = epsilon .* alpha;
            if sum(cal_i) ~= 0 %���ǵ���0����ø��岻�����������Ӱ�졣
                gamma_i = sum(cal_i .* delta);
                q_star = sum(cal_i .* q0) / sum(cal_i);
                nervous_i = 1 - (1 - q_star) * (1 - q0(i));
                calm_i = q0(i) * q_star;
                dq = gamma_i * (beta_i(i) * nervous_i + (1 - beta_i(i)) * calm_i) * t;
                q0(i) = q0(i) + dq;
                if q0(i) >= 1 %q���ܳ���1��
                    q0(i) = 0.99;
                end
            end
        end
        
%===================�ɼ���ͼ���ݣ���¼0s��25s��50s����ʱ�̼�������================        
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

%=============================�˶��켣��ͼ=================================
for picture = 1:Nrecord
    figure(picture)
%ǽ��������
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
%����λ�á�
    circleangle = 0:(pi / 200):(2 * pi);
    for people = 1:num
        plot((MT(1,people,picture) + r * cos(circleangle)), (MT(2,people,picture) + r * sin(circleangle)), '-r')
        hold on
    end
end