function IEN = IEN_all(num, points)
    % IEN_all: generate IEN matrix based on num and points
    n = max(size(points));
    start = 1;
    IEN = [];
    for i = 1:n
        if i == 1
            %IEN = [IEN, IEN_coplanar(num(start:start+points(i)-1))];
            IEN = [IEN, IEN_bificial(num(end), num(start:start+points(i)-1))];
            %disp(IEN)
            start = start + points(i);
        else
            %IEN = [IEN, IEN_coplanar(num(start:start+points(i)-1))];
            start1 = start - points(i-1);
            IEN = [IEN, IEN_bificial(num(start1:start1+points(i-1)-1), num(start:start+points(i)-1))];
            %Draw_bificial(figure, pos(start1:start1+points(i-1)-1,:), pos(start:start+points(i)-1,:))
            start = start + points(i);
        end
    end
end

function IEN = IEN_coplanar(num)
    n = size(num, 1);
    IEN = [];
    for i = 1:n-1
        IEN = [IEN, [num(i); num(i+1)]];
    end
    IEN = [IEN, [num(1) ; num(n)]];
end

function IEN = IEN_bificial( num1, num2)
    n1 = size(num1, 1);
    n2 = size(num2, 1);
    IEN = [];
    diff = n2-n1;
    if diff == 0
        % 两者相等
        for i = 1:n2
            if i ~= n2
                IEN = [IEN, [num1(i); num1(i+1); num2(i)]];
                IEN = [IEN, [num2(i); num2(i+1); num1(i+1)]];
            else
                IEN = [IEN, [num1(1);num1(i); num2(i)]];
                IEN = [IEN, [num2(1);num2(i); num1(1)]];
            end
        end
    else
        % 两者不等
        if n1 == 1 %初始情况
            for i = 1:n2
                if i~= n2
                    IEN = [IEN, [num1(1); num2(i); num2(i+1)]];
                else
                    IEN = [IEN, [num1(1); num2(i); num2(1)]];
                end
            end
        else % 非初始情况
            gap = n1/diff; % 间隔多少跳变一次
            skip = 0;
            for i = 1:n1
                if i == 1 % 初始情况第一个点
                    IEN = [IEN, [num1(i); num2(i+skip); num2(end)]];

                    IEN = [IEN, [num1(i); num1(end); num2(end)]]; % the down floor
                    %IEN = [IEN, [num1(i); num1(end); num2(1)] ];
                    % find if skip
                    if mod(i, gap) == 0
                        skip = skip + 1;
                        IEN = [IEN, [num1(i); num2(i+skip); num2(i+skip-1)] ];
                    end

                else
                    IEN = [IEN, [num1(i-1); num1(i); num2(i+skip-1)] ]; % the down floor
                    IEN = [IEN, [num1(i); num2(i+skip); num2(i+skip-1)]];
                    % find if skip
                    if mod(i, gap) == 0
                        skip = skip + 1;
                        IEN = [IEN, [num1(i); num2(i+skip); num2(i+skip-1)]];
                    end
                end

            end
        end
    end
end


