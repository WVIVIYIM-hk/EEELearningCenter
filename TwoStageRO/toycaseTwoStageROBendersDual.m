%% Clean up
clear; clc; close all;

%% 1. 参数定义 (Parameter Definition)
num_nodes = 3; %3个节点
f = [400; 414; 326];    % 固定建设成本 Fixed cost
a = [18; 25; 20];       % 单位容量建设成本 Unit capacity cost
C_trans = [             % 单位运输成本 Transportation cost
    22, 33, 24;
    33, 23, 30;
    20, 25, 27
];
K = [800; 800; 800];% 最大允许建设容量

% 不确定集参数
d_bar = [206; 274; 220]; % 基准需求
d_hat = [40; 40; 40]; % 最大需求偏差
Gamma = 1.8;             

%  惩罚系数 (Penalty Cost): 代表切负荷的高昂成本。
%  它的作用是保证子问题永远可行。
%  只要这个值比正常的运输成本(20-30)大很多即可。
PENALTY_COST = 5000; 
% [附录 A-4.3] 计算最大可能总需求，用于主问题预约束
% 简单估算：所有基准值 + Gamma * 最大偏差 (最保守估计)
Max_Total_Demand_Est = sum(d_bar) + Gamma * max(d_hat);

%% 2. 初始化
LB = -inf;
UB = inf;
epsilon = 1e-3;
Max_Iter = 30;

Cut_Pi = []; % = pi_sp 由子问题对偶问题求解得到pi_sp具体数值，用于主问题的割平面
Cut_Constant = []; % = sum(lambda_sp .* d) 由子问题对偶问题求解得到lambda_sp、g具体数值，用于主问题的割平面

% 主问题变量
y = binvar(num_nodes, 1);% 0,1变量(3,1)
z = sdpvar(num_nodes, 1);% 连续决策变量(3,1)
eta = sdpvar(1, 1);% 子问题辅助变量(1,1)

fprintf('Starting Robust Benders-Dual (Penalized Method)...\n');
fprintf('%-5s | %-12s | %-12s | %-10s\n', 'Iter', 'LB', 'UB', 'Gap');

for iter = 1:Max_Iter
    
    %% 3. 求解主问题 (Master Problem)
    % 目标函数
    Obj_MP = f'*y + a'*z + eta;
    % 主问题约束
    Constraints_MP = [];
    % 容量建设约束
    Constraints_MP = [Constraints_MP, z <= K .* y];
    % 设施容量非负
    Constraints_MP = [Constraints_MP, z >= 0];
    % 第二阶段运输成本非负
    Constraints_MP = [Constraints_MP, eta >= 0];
    % [附录 A-4.3] 添加总容量下限约束 
    % 这能避免第一轮 z=0 导致的无界，加速收敛
    Constraints_MP = [Constraints_MP, sum(z) >= Max_Total_Demand_Est];

    % 添加 Benders Cuts
    % eta >= pi_k * z + sum(lambda_k * d_k)
    for k = 1:length(Cut_Constant)
        pi_vec = Cut_Pi(:, k);
        const_val = Cut_Constant(k);
        %
        Constraints_MP = [Constraints_MP, eta >= pi_vec' * z + const_val];
    end
    
    % 求解
    ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
    sol_MP = optimize(Constraints_MP, Obj_MP, ops);
    
    if sol_MP.problem == 0 
        % 求解成功，更新y、z值传递给子问题
        y_star = value(y);
        z_star = value(z);
        eta_star = value(eta);
        % 更新主问题的下界
        Current_LB = value(Obj_MP);
    else
        % 第一代若无界(unbounded)，手动给初值
        if iter == 1
             y_star = [1;1;1]; z_star = K; eta_star = 0; Current_LB = -inf;
        else
             error('Master Problem Failed. Error Code: %d. Msg: %s', sol_MP.problem, sol_MP.info);
        end
    end
    LB = Current_LB;
    
    %% 4. 求解子问题 (Subproblem)

    % 对偶变量
    % pi_sp: 对应容量约束 z - sum(x) >= 0
    % lambda_sp: 对应需求约束 sum(x) - d >= 0
    pi_sp = sdpvar(num_nodes, 1);      
    lambda_sp = sdpvar(num_nodes, 1);  
    
    % 不确定变量g(3,1)
    g = sdpvar(num_nodes, 1);
    % 实际需求
    d = d_bar + d_hat .* g; 
    % 子问题对偶问题的目标函数
    Obj_SPDP=sum(d.*lambda_sp)-sum(z_star.*pi_sp);
    
    % 约束
    Constraints_SPDP=[];
    % 对偶约束: C_ij + pi_i - lambda_j >= 0
    for i = 1:num_nodes
        for j = 1:num_nodes 
            Constraints_SPDP = [Constraints_SPDP, C_trans(i,j) + pi_sp(i) - lambda_sp(j) >= 0];
        end
    end
    % 对偶变量非负
    Constraints_SPDP = [Constraints_SPDP,pi_sp>=0];
    Constraints_SPDP = [Constraints_SPDP,lambda_sp>=0];
    % 【关键修正】 对偶变量上界 (防止 Unbounded)
    % 只要系统出现不可行（容量不足），lambda 和 pi 就会顶到这个上界
    % 这相当于告诉主问题：由于容量不足，产生了巨大的惩罚成本
    Constraints_SPDP = [Constraints_SPDP, lambda_sp <= PENALTY_COST];
    Constraints_SPDP = [Constraints_SPDP, pi_sp <= PENALTY_COST];

    % 不确定集约束
    Constraints_SPDP = [Constraints_SPDP,g>=0,g<=1];% 范围
    Constraints_SPDP = [Constraints_SPDP,sum(g)<=Gamma];% 总预算
    Constraints_SPDP = [Constraints_SPDP,g(1)+g(2)<=1.2];% 局部约束
    % 求解
    ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
    ops.gurobi.NonConvex = 2;% 解决 g*lambda 的非凸性
    sol_SP = optimize(Constraints_SPDP, -Obj_SPDP, ops); % Maximize

    if sol_SP.problem ~= 0
        % 若报错，显示原因
        error('Subproblem Failed. Status: %d. Info: %s. Please check Solver!', sol_SP.problem, sol_SP.info);
    end
    
    % 提取结果
    eta_val = value(Obj_SPDP);

    % 注意：我们在Cut里需要的 pi 是对应 "sum(x) <= z" 的对偶变量（通常 <= 0）。
    % 我们的 KKT 里定义的是 "z - sum(x) >= 0" 的对偶变量 pi_sp >= 0。
    % 它们的关系是：Value(Cut_pi) = -Value(pi_sp)
    % 物理意义：增加容量 z，会降低成本，所以 z 的系数应为负。
    pi_val = -value(pi_sp); 
    lambda_val = value(lambda_sp);
    d_val = value(d);

    % 更新 UB
    Current_UB = f'*y_star + a'*z_star + eta_val;
    UB = min(UB, Current_UB);

    % 主问题上下界是否收敛
    Gap = abs(UB - LB)/abs(LB);
    fprintf('%-5d | %-12.2f | %-12.2f | %-10.4f\n', iter, LB, UB, Gap);

    if Gap <= epsilon
        % 主问题上下界收敛，完成求解
        fprintf('Converged successfully!\n');
        break;
    end

    % 更新 Cut
    Cut_Pi = [Cut_Pi, pi_val];
    Cut_Constant = [Cut_Constant, lambda_val' * d_val];

end
% 
if iter == Max_Iter
    warning('Reached Max Iterations without full convergence.');
end

%% 输出最终结果
fprintf('\n=== Final Results ===\n');
fprintf('Facilities Built (y): %s\n', mat2str(round(y_star')));
fprintf('Capacities Built (z): %s\n', mat2str(round(z_star')));
fprintf('Worst-Case Total Cost: %.2f\n', UB);
fprintf('Solver Used: %s\n', sol_MP.solver);