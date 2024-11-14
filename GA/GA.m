clear, clc, close all
tic
Pc = 0.8;           %交叉概率
Pm = 0.1;           %变异概率
Pn = 0.2;           %逆转概率
max_iter = 3000;    %进化次数
pop_size = 500;     %种群个数
power = [0.6, 2.4, 1.6, 0.4, 0.4, 0.4, 1.5, 1.5, 2, 3.5, 3, 1.5, 3.2];
kfzd = [0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1];
t = [8, 8, 8, 8, 8, 8, 4, 4, 4, 16, 20, 6, 10];
lower = [1, 1, 29, 1, 33, 65, 21, 65, 69, 1, 29, 85, 29];
upper = [96, 96, 88, 24, 56, 88, 32, 84, 96, 96, 60, 96, 92];
tline = 96;
pr = zeros(1, tline);
for j = 1:1:tline
    if j<=24
        pr(j) = 0.3;
    elseif j>=69 && j<=88
        pr(j) = 0.9;
    elseif (j>=37 && j<=48) || (j>=61 && j<=68)
        pr(j) = 0.7;
    elseif (j>=25 && j<=36) || (j>=49 && j<=60) || (j>=89 && j<=96)
        pr(j) = 0.5;
    end
end
recordz = zeros(length(t), max_iter);
xz = zeros(length(t), tline);
parfor ip = 1:length(t)
    [recordz(ip,:), xz(ip,:)] = ...
        main_procedure(Pc, Pm, Pn, max_iter, pop_size, power, kfzd, t, lower, upper, pr, ip);
end
plot(sum(recordz,1))
xlabel('迭代次数')
ylabel('Loss')
disp('最优调度')
disp(num2str(xz))
toc