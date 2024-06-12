function rms = loss_cal(IEN, pos, func)
    % 计算整体的loss-rms
    loss_all = 0;
    area_all = 0;
    for i = 1:size(IEN,2)
        IEN_part = IEN(:,i);
        [loss, area] = loss_part(IEN_part, pos, func);
        loss_all = loss_all + loss;
        area_all = area_all + area;
    end
    rms = sqrt(loss_all/area_all);
end


function [loss0, area] = loss_part(IEN_part, pos, func)
    pos_part = pos(IEN_part',:);
    part = @(eta1, eta2, eta3, data) eta1*data(1,:) + eta2*data(2,:) + eta3*data(3,:);
    diffz = @(eta1, eta2, eta3, data) func(part(eta1, eta2, eta3, data(:,1:2)))-part(eta1, eta2, eta3, data(:,3));
    area = triangle_area(pos_part);
    loss = 0;
    % loss = loss+diffz(1/3, 1/3, 1/3, pos_part).^2*(-27/96);
    % loss = loss+diffz(0.6, 0.2, 0.2, pos_part).^2*(25/96);
    % loss = loss+diffz(0.2, 0.6, 0.2, pos_part).^2*(25/96);
    % loss = loss+diffz(0.2, 0.2, 0.6, pos_part).^2*(25/96);
    loss = loss + diffz(0.81685, 0.091576, 0.091576, pos_part).^2*0.10995./2;
    loss = loss + diffz(0.091576, 0.81685, 0.091576, pos_part).^2*0.10995./2;
    loss = loss + diffz(0.091576, 0.091576, 0.81685,  pos_part).^2*0.10995./2;
    loss = loss + diffz(0.108103, 0.445948, 0.445948, pos_part).^2*0.22338./2;
    loss = loss + diffz(0.445948, 0.108103, 0.445948, pos_part).^2*0.22338./2;
    loss = loss + diffz(0.445948, 0.445948, 0.108103, pos_part).^2*0.22338./2;
    loss0 = loss*area*2;
end

function area = triangle_area(data)
    A = data(1,:);
    B = data(2,:);
    C = data(3,:);
    AB = B - A;
    AC = C - A;
    cross_product = cross(AB, AC);
    area = 0.5 * norm(cross_product);
end