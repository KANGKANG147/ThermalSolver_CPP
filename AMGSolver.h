#pragma once
#include "LinearAlgebra.h"
#include <list>

// AMG 的每一层
struct AMGLevel {
    SparseMatrix A; // 当前层的系统矩阵
    SparseMatrix P; // 插值算子 (Coarse -> Fine)
    SparseMatrix R; // 限制算子 (Fine -> Coarse, 通常 R = P^T)

    Vector x;       // 当前层的解 (guess)
    Vector b;       // 当前层的右端项

    int level_id;
};

class AMGSolver {
public:
    std::vector<AMGLevel> levels;
    int max_levels = 10;
    int coarse_size_limit = 50; // 当矩阵小于这个大小时停止粗化

    // =====================================
    // Setup Phase: 构建 AMG 层级结构
    // =====================================
    void setup(const SparseMatrix& A_fine) {
        levels.clear();

        // Level 0 是最细的网格
        AMGLevel fine_level;
        fine_level.A = A_fine;
        fine_level.level_id = 0;
        levels.push_back(fine_level);

        std::cout << "[AMG Setup] Level 0: " << A_fine.rows << " vars." << std::endl;

        // 循环粗化
        while (levels.size() < max_levels) {
            AMGLevel& current = levels.back();

            // 1. 如果足够小，停止
            if (current.A.rows <= coarse_size_limit) break;

            // 2. 生成插值算子 P
            SparseMatrix P = build_interpolation(current.A);

            // 如果没法粗化了 (比如全是独立点)
            if (P.cols == 0 || P.cols == current.A.rows) break;

            SparseMatrix R = P.transpose();

            // 3. Galerkin 投影生成粗网格矩阵: A_coarse = R * A * P
            // 这里用到了我们刚才实现的矩阵乘法
            SparseMatrix RA = SparseMatrix::multiply(R, current.A);
            SparseMatrix Ac = SparseMatrix::multiply(RA, P);

            // 4. 保存当前层的算子
            current.P = P;
            current.R = R;

            // 5. 创建下一层
            AMGLevel coarse_level;
            coarse_level.A = Ac;
            coarse_level.level_id = levels.size();
            levels.push_back(coarse_level);

            std::cout << "[AMG Setup] Level " << coarse_level.level_id
                << ": " << Ac.rows << " vars (Coarsening Ratio: "
                << (double)current.A.rows / Ac.rows << ")" << std::endl;
        }
    }

    // =====================================
    // Solve Phase: V-Cycle 求解
    // =====================================
    void solve(const Vector& b, Vector& x) {
        // 将 b 和 x 放入 level 0
        levels[0].b = b;
        levels[0].x = x;

        // 开始递归 V-Cycle
        v_cycle(0);

        // 取回结果
        x = levels[0].x;
    }

private:
    // 递归 V-Cycle
    void v_cycle(int level_idx) {
        AMGLevel& lvl = levels[level_idx];

        // 1. 底部求解 (Coarsest Grid Solve)
        if (level_idx == levels.size() - 1) {
            // 最粗网格直接用 Gauss-Seidel 迭代多次求精确解
            // 或者用高斯消元(这里简化为多步迭代)
            lvl.A.smooth_gauss_seidel(lvl.b, lvl.x, 50);
            return;
        }

        // 2. 前平滑 (Pre-smoothing)
        lvl.A.smooth_gauss_seidel(lvl.b, lvl.x, 2); // 2 次迭代消除高频误差

        // 3. 计算残差 r = b - Ax
        Vector r(lvl.A.rows);
        lvl.A.residual(lvl.b, lvl.x, r);

        // 4. 限制残差到下一层 r_coarse = R * r_fine
        AMGLevel& next_lvl = levels[level_idx + 1];
        next_lvl.b.resize(next_lvl.A.rows);
        lvl.R.multiply(r, next_lvl.b);

        // 初始化下一层的猜想 x = 0 (我们求解的是误差方程 Ae = r)
        next_lvl.x.assign(next_lvl.A.rows, 0.0);

        // 5. 递归调用
        v_cycle(level_idx + 1);

        // 6. 延拓修正量 e_fine = P * e_coarse
        Vector e_fine(lvl.A.rows);
        lvl.P.multiply(next_lvl.x, e_fine);

        // 7. 修正当前层解 x = x + e_fine
        for (size_t i = 0; i < lvl.x.size(); ++i) lvl.x[i] += e_fine[i];

        // 8. 后平滑 (Post-smoothing)
        lvl.A.smooth_gauss_seidel(lvl.b, lvl.x, 2);
    }

    // =====================================
    // AMG 核心: 经典 Ruge-Stuben 插值构建
    // =====================================
    SparseMatrix build_interpolation(const SparseMatrix& A) {
        int N = A.rows;
        const double theta = 0.25; // 强连接阈值

        // --- Step 1: 建立强连接图 ---
        std::vector<std::vector<int>> strong_con(N); // i 强依赖于 j
        std::vector<std::vector<int>> strong_depend(N); // j 被 i 强依赖 (transpose)

        for (int i = 0; i < N; ++i) {
            double max_neg = 0.0;
            // 找对角线和行最大负值
            for (int k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
                if (A.col_indices[k] != i && A.values[k] < 0)
                    max_neg = std::max(max_neg, -A.values[k]);
            }

            // 判定强连接
            for (int k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
                int j = A.col_indices[k];
                if (i != j && -A.values[k] >= theta * max_neg) {
                    strong_con[i].push_back(j);
                    strong_depend[j].push_back(i);
                }
            }
        }

        // --- Step 2: C/F 分裂 (经典算法) ---
        // 简单实现：贪心算法
        // 权重 = 该点被多少个未决定的点强依赖
        std::vector<int> weights(N);
        for (int i = 0; i < N; ++i) weights[i] = strong_depend[i].size();

        std::vector<int> cf_status(N, 0); // 0=Undecided, 1=C, 2=F
        int num_c_points = 0;

        // 简单策略：优先选权重大的为 C 点
        // 真实 Ruge-Stuben 会动态更新权重，这里用静态简化版加速
        for (int i = 0; i < N; ++i) {
            if (cf_status[i] != 0) continue;

            // 选 i 为 C 点
            cf_status[i] = 1;
            num_c_points++;

            // 它的所有强邻居变成 F 点
            for (int j : strong_depend[i]) {
                if (cf_status[j] == 0) cf_status[j] = 2;
            }
        }
        // 剩下的全部设为 C (防止孤立)
        for (int i = 0; i < N; ++i) {
            if (cf_status[i] == 0) { cf_status[i] = 1; num_c_points++; }
        }

        // 建立 C 点索引映射 (Coarse Grid Index)
        std::vector<int> c_idx_map(N, -1);
        int c_count = 0;
        for (int i = 0; i < N; ++i) {
            if (cf_status[i] == 1) c_idx_map[i] = c_count++;
        }

        // --- Step 3: 构建插值矩阵 P ---
        // P 的维度是 N x num_c_points
        MatrixBuilder pb(N); // 用于构建 P

        for (int i = 0; i < N; ++i) {
            if (cf_status[i] == 1) {
                // 如果是 C 点，直接注入 (值=1.0)
                pb.add(i, c_idx_map[i], 1.0);
            }
            else {
                // 如果是 F 点，根据强连接的 C 点插值
                // 公式： w_ij = - a_ij / ( a_ii + sum_weak(a_ik) )  (简化版直接插值)

                double sum_strong_c = 0.0;
                double diag = A.get_diag(i);

                // 收集分母和分子
                for (int k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
                    int j = A.col_indices[k];
                    if (cf_status[j] == 1) { // 是 C 点邻居
                        // 检查是否强连接，如果是，加入插值权重
                        // 这里简化为所有 C 点邻居都插值
                        sum_strong_c += A.values[k];
                    }
                }

                // 避免除零
                if (std::abs(sum_strong_c) < 1e-10) {
                    // 孤立 F 点？虽然罕见，但给个 fallback
                    continue;
                }

                // 简单的等权分配 (Direct Interpolation 简化)
                // 真实公式很复杂，这里假设 F 点的值是周围 C 点的加权平均
                // w_j = a_ij / sum(a_ik for k in C)
                for (int k = A.row_ptr[i]; k < A.row_ptr[i + 1]; ++k) {
                    int j = A.col_indices[k];
                    if (cf_status[j] == 1) {
                        double weight = A.values[k] / sum_strong_c;
                        pb.add(i, c_idx_map[j], weight);
                    }
                }
            }
        }

        SparseMatrix P = pb.build();
        P.cols = num_c_points; // 修正列数
        return P;
    }
};