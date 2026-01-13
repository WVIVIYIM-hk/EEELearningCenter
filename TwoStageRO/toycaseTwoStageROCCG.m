%% Clean up
clear; clc; close all;

%% 1. 参数定义 (Parameter Definition)
num_facilities = 3; %3个源点
num_customers = 3; %3个终点
% 成本参数
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
PENALTY_COST = 5000; 
% [附录 A-4.3] 计算最大可能总需求，用于主问题预约束
% 简单估算：所有基准值 + Gamma * 最大偏差 (最保守估计)
Max_Total_Demand_Est = sum(d_bar) + Gamma * max(d_hat);

%% 2. 初始化
LB = -inf;
UB = inf;
epsilon = 1e-3;
Max_Iter = 30;

% 用于存储识别出的"最坏场景" (Worst-case Scenarios)
% CCG 的核心：保存每一次子问题找到的 d 值
Scenario_Set = {}; 

% 初始给一个基准场景 (Nominal Scenario) 启动算法
Scenario_Set{1} = d_bar; 

fprintf('Starting CCG Algorithm (Replication of Table 1)...\n');
fprintf('%-5s | %-12s | %-12s | %-10s\n', 'Iter', 'LB', 'UB', 'Gap');
for iter = 1:Max_Iter
%% ==========================================
    %% 3. 求解主问题 (Master Problem, MP)
    %% ==========================================
    % CCG 主问题：包含所有已知场景的变量 x^{(l)}
    
    % 第一阶段变量
    y = binvar(num_facilities, 1); % 选址
    z = sdpvar(num_facilities, 1); % 容量
    eta = sdpvar(1, 1);            % 辅助变量 (代表第二阶段的最坏成本)
    
    % 目标函数
    Obj_MP = f'*y + a'*z + eta;
    % 主问题约束
    Constraints_MP = [];
    % 容量建设约束
    Constraints_MP = [Constraints_MP, z <= K .* y]; 
    % 设施容量非负
    Constraints_MP = [Constraints_MP, z >= 0];
    % [附录 A-4.3] 添加总容量下限约束 
    % 这能避免第一轮 z=0 导致的无界，加速收敛
    Constraints_MP = [Constraints_MP, sum(z) >= Max_Total_Demand_Est];

    % --- CCG 核心：为每个已知场景 l 添加一套约束 ---
    for l = 1:length(Scenario_Set)
        d_current = Scenario_Set{l}; % 取出第 l 个场景的需求数值
        
        % 创建针对场景 l 的第二阶段变量 x^(l) (num_facilities x num_customers)
        % 变量命名技巧：在 YALMIP 中由于循环变长，建议临时定义，放入约束
        x_l = sdpvar(num_facilities, num_customers, 'full'); 
    
        % 需求限制
        Constraints_MP = [Constraints_MP, sum(x_l, 1)' >= d_current];
        % % 主问题增加slack 和 PENALTY_COST的步骤
        % % 允许切负荷(slack)以保证可行性，但施加巨额惩罚
        % slack_l = sdpvar(num_customers, 1); 
        % 
        % Constraints_MP = [Constraints_MP, sum(x_l, 1)' + slack_l >= d_current];
        % Constraints_MP = [Constraints_MP, slack_l >= 0];
    
        % 容量限制
        Constraints_MP = [Constraints_MP, sum(x_l, 2) <= z];
        % 运输量非负
        Constraints_MP = [Constraints_MP, x_l >= 0];
        
        % 第二阶段成本约束（Epigraph 约束）eta 必须大于所有已知场景的成本
        Cost_Scenario = sum(sum(C_trans .* x_l));
        % Cost_Scenario = sum(sum(C_trans .* x_l)) + PENALTY_COST *
        % sum(slack_l); % 加入松弛后的第二阶段成本约束
        Constraints_MP = [Constraints_MP, eta >= Cost_Scenario(:)];
    end
    
    % 求解 MP
    ops = sdpsettings('solver', 'gurobi', 'verbose', 0);
    sol_MP = optimize(Constraints_MP, Obj_MP, ops);
    
    if sol_MP.problem ~= 0
        error('Master Problem Failed! Info: %s', sol_MP.info);
    end
    
    % 更新下界 LB
    % LB = c*y* + eta*
    y_star = value(y);
    z_star = value(z);
    eta_star = value(eta);
    Current_LB = value(Obj_MP);
    LB = Current_LB; % CCG 的 MP 随着约束增加，LB 单调递增
    
    %% ==========================================
    %% 4. 求解子问题 (Subproblem, SP)
    %% ==========================================
    % 给定 y*, z*，寻找使成本最大的最坏需求 d
    % 这是一个 Max-Min 问题，转化为 Max-Dual
    
    % 对偶变量和不确定性变量
    pi_sp = sdpvar(num_facilities, 1);      % 对应容量约束 (dual variable)
    lambda_sp = sdpvar(num_customers, 1);   % 对应需求约束 (dual variable)
    g = sdpvar(num_customers, 1);           % 不确定性变量
    
    d_real = d_bar + d_hat .* g;            % 实际需求表达式
    
    % --- 目标函数 (Maximize Dual Objective) ---
    % Original Primal: Min cx s.t. Ax >= d, Bx <= z
    % Dual Obj: lambda*d - pi*z
    % Gurobi NonConvex 处理双线性项: lambda * d_real
    Obj_SPDP = sum(lambda_sp .* d_real) - sum(pi_sp .* z_star);
    
    % 子问题约束
    Constraints_SPDP = [];
    
    % 对偶约束 (Dual Feasibility): lambda_j - pi_i <= c_ij
    for i = 1:num_facilities
        for j = 1:num_customers
            Constraints_SPDP = [Constraints_SPDP, lambda_sp(j) - pi_sp(i) <= C_trans(i,j)];
        end
    end
    
    % 对偶变量非负 & 有界 (防止 Unbounded)
    Constraints_SPDP = [Constraints_SPDP, 0 <= pi_sp <= PENALTY_COST];
    Constraints_SPDP = [Constraints_SPDP, 0 <= lambda_sp <= PENALTY_COST];
    
    % 不确定集约束 (Uncertainty Set)
    Constraints_SPDP = [Constraints_SPDP, 0 <= g <= 1];
    Constraints_SPDP = [Constraints_SPDP, sum(g) <= Gamma];
    Constraints_SPDP = [Constraints_SPDP, g(1) + g(2) <= 1.2];
    
    % 求解子问题
    ops_sp = sdpsettings('solver', 'gurobi', 'verbose', 0);
    ops_sp.gurobi.NonConvex = 2; % 开启非凸支持
    
    sol_SP = optimize(Constraints_SPDP, -Obj_SPDP, ops_sp);
    
    if sol_SP.problem ~= 0
        error('Subproblem Failed! Info: %s', sol_SP.info);
    end
    
    % 提取最坏场景
    d_worst = value(d_real);
    Actual_Subproblem_Cost = value(Obj_SPDP);
    
    % 更新上界 UB
    % UB = min(UB, FirstStageCost + MaxRecourseCost)
    Current_UB = f'*y_star + a'*z_star + Actual_Subproblem_Cost;
    UB = min(UB, Current_UB);
    
    % 计算 Gap
    Gap = abs(UB - LB) / (abs(LB) + 1e-5);
    
    fprintf('%-5d | %-12.2f | %-12.2f | %-10.4f\n', iter, LB, UB, Gap);
    
    %% ==========================================
    %% 5. 收敛检查与列生成
    %% ==========================================
    if Gap <= epsilon
        fprintf('\nCCG Converged!\n');
        break;
    else
        % --- CCG 核心步骤：添加新列 (Column Generation) ---
        % 将找到的最坏需求 d_worst 加入场景集合
        % 下一次 MP 会自动为这个新场景创建对应的变量 x 和约束
        Scenario_Set{end+1} = d_worst;
    end
    
end

%% 结果输出
fprintf('\n=== Final Results (Comparison with Table 1) ===\n');
fprintf('Facilities Built (y): %s\n', mat2str(round(value(y)')));
fprintf('Capacities Built (z): %s\n', mat2str(round(value(z)')));
fprintf('Optimal Cost: %.2f\n', UB);
fprintf('Total Iterations: %d\n', iter);