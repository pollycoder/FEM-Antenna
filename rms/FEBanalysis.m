clear,clc, close all
load pos.mat
pos = [x,y,z];
pos = [pos;[0 0 0]]; % 单位：m
num = [1:1:(size(pos,1))]'; %符号记载
points = [6, 12, 24, 36, 36];
IEN = IEN_all(num, points); %可以在迭代开始的时候只生成一次，循环使用
precious_z = @(pos_xy) (pos_xy(:,1).^2+pos_xy(:,2).^2)./(4*2.17); %计算准确的抛物面句柄
rms = loss_cal(IEN, pos, precious_z) %计算rms的时候需要用到IEN和pos


% E = 2.1e11; % Pa
% rho = 7850; % kg/m^3


% E = 3.9e11; % Pa
% rho = 1840; % kg/m^3

% A = pi*(5e-3)^2; % m^2
% m_node = rho*1e-2*pi*(2e-2)^2; % kg
% K = assembely_K(LM, pos, E, A);
% M = assembely_M(LM, pos, rho, A, m_node);

% %analysis frequency
% [V, D] = eig(M^(-1)*K);
% f = real(sqrt(diag(D))/(2*pi));
% f = sort(f);
% f(1:50)
% plot(f, 'o')


function plot_IEN(pos, IEN)
    % 基于IEN绘制单元，检验IEN的正确性
    for i = 1:size(IEN, 2)
        terminal = [pos(IEN(1,i),:); pos(IEN(2,i),:); pos(IEN(3,i),:)];
        center = sum(terminal)/3; 
        plot3(terminal(1:2,1), terminal(1:2,2), terminal(1:2,3), 'b')
        hold on
        plot3(terminal(2:3,1), terminal(2:3,2), terminal(2:3,3), 'b')
        plot3([terminal(1,1), terminal(3,1)], [terminal(1,2), terminal(3,2)], [terminal(1,3), terminal(3,3)], 'b')
        %plot3(terminal(:,1), terminal(:,2), terminal(:,3), 'b')
        scatter3(center(1), center(2), center(3), 'r')
        hold on
    end
    view([0,0,90])
end


function plot_all(pos, LM)
    % 基于LM绘制单元，检验LM的正确性
    for i = 1:size(LM, 2)
        terminal = [pos(LM(1,i),:); pos(LM(2,i),:)];
        plot3(terminal(:,1), terminal(:,2), terminal(:,3), 'k')
        hold on
    end
    view([0,0,90])
end

function K = assembely_K(LM, pos, E, A)
    n = size(LM, 2);
    K = zeros(3*size(pos,1), 3*size(pos,1));
    for i = 1:n
        terminal = [pos(LM(1,i),:); pos(LM(2,i),:)];
        l = terminal(2,:) - terminal(1,:); % 行向量
        L = norm(l);
        l = l/L;
        k = l'*l;
        k = E*A/L*[k, -k; -k, k];
        K(3*LM(1,i)-2:3*LM(1,i), 3*LM(1,i)-2:3*LM(1,i)) = K(3*LM(1,i)-2:3*LM(1,i), 3*LM(1,i)-2:3*LM(1,i)) + k(1:3, 1:3);
        K(3*LM(1,i)-2:3*LM(1,i), 3*LM(2,i)-2:3*LM(2,i)) = K(3*LM(1,i)-2:3*LM(1,i), 3*LM(2,i)-2:3*LM(2,i)) + k(1:3, 4:6);
        K(3*LM(2,i)-2:3*LM(2,i), 3*LM(1,i)-2:3*LM(1,i)) = K(3*LM(2,i)-2:3*LM(2,i), 3*LM(1,i)-2:3*LM(1,i)) + k(4:6, 1:3);
        K(3*LM(2,i)-2:3*LM(2,i), 3*LM(2,i)-2:3*LM(2,i)) = K(3*LM(2,i)-2:3*LM(2,i), 3*LM(2,i)-2:3*LM(2,i)) + k(4:6, 4:6);
    end
end

function M = assembely_M(LM, pos, rho, A, m_node)
    n = size(LM, 2);
    %M = zeros(3*size(pos,1), 3*size(pos,1));
    M  = m_node*diag(ones(3*size(pos,1),1));
    for i = 1:n
        terminal = [pos(LM(1,i),:); pos(LM(2,i),:)];
        l = terminal(2,:) - terminal(1,:); % 行向量
        L = norm(l);
        m = rho*A*L/2*diag(ones(6,1));
        M(3*LM(1,i)-2:3*LM(1,i), 3*LM(1,i)-2:3*LM(1,i)) = M(3*LM(1,i)-2:3*LM(1,i), 3*LM(1,i)-2:3*LM(1,i)) + m(1:3, 1:3);
        M(3*LM(1,i)-2:3*LM(1,i), 3*LM(2,i)-2:3*LM(2,i)) = M(3*LM(1,i)-2:3*LM(1,i), 3*LM(2,i)-2:3*LM(2,i)) + m(1:3, 4:6);
        M(3*LM(2,i)-2:3*LM(2,i), 3*LM(1,i)-2:3*LM(1,i)) = M(3*LM(2,i)-2:3*LM(2,i), 3*LM(1,i)-2:3*LM(1,i)) + m(4:6, 1:3);
        M(3*LM(2,i)-2:3*LM(2,i), 3*LM(2,i)-2:3*LM(2,i)) = M(3*LM(2,i)-2:3*LM(2,i), 3*LM(2,i)-2:3*LM(2,i)) + m(4:6, 4:6);
    end
end