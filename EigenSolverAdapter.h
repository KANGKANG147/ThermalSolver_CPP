#pragma once

// 引入你原有的线性代数定义
#include "LinearAlgebra.h" 

// 引入 Eigen 求解器模块
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <iostream>

class EigenSolverAdapter {
public:
    // 定义映射类型：让 Eigen 直接读你的 std::vector，不复制内存
    // RowMajor 对应 CSR 格式
    using EigenSpMatMap = Eigen::Map<const Eigen::SparseMatrix<double, Eigen::RowMajor, int>>;

    // =================================================================================
    // 求解器选择
    // 热传导矩阵通常是对称正定(SPD)的，使用 ConjugateGradient (CG) 比 BiCGSTAB 快且省内存。
    // 预处理使用 IncompleteCholesky (IC) 替代 IncompleteLUT。
    // =================================================================================

    // 方案 A (推荐): 针对对称正定矩阵 (速度最快)
    using SolverType = Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::Lower | Eigen::Upper, Eigen::IncompleteCholesky<double>>;

    // 方案 B (备用): 针对非对称矩阵 (稳定双共轭梯度法 (BiCGSTAB) + IncompleteLUT 预处理,如果以后引入了流体平流项导致不对称，可切回这个)
    // using SolverType = Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::IncompleteLUT<double>>;

    SolverType solver;

    // 1. 设置矩阵 (Setup)
    // 返回 bool: true 表示分解成功，false 表示矩阵奇异或失败
    bool setup(const SparseMatrix& A) {
        // 【关键】构建 Map
        // 直接使用 A.rows, A.cols 以及 val/ind/ptr 的数据指针
        EigenSpMatMap mat_map(
            A.rows,
            A.cols,
            static_cast<int>(A.values.size()),
            A.row_ptr.data(), // 行偏移指针
            A.col_indices.data(), // 列索引指针
            A.values.data()      // 数值指针
        );

        // 配置求解器参数 (根据需要调整)
        solver.setTolerance(1e-6); // 收敛残差

        // 【注意】CG 不需要 droptol/fillfactor，但 IncompleteCholesky 可能需要 shift 防止非正定
        // 如果是 BiCGSTAB + ILUT，则保留下面的设置：
        // solver.preconditioner().setDroptol(0.001);
        // solver.preconditioner().setFillfactor(5);
        
        // 分析矩阵结构并计算预处理
        solver.compute(mat_map);

        if (solver.info() != Eigen::Success) {
            std::cerr << "[Eigen] Setup failed! Matrix might be singular or not SPD." << std::endl;
            return false;
        }
        return true;
    }

    // 2. 求解 (Solve)
    // 对应原来的 amg_solver.solve(b, x)
    bool solve(const Vector& b, Vector& x) {
        // 将 std::vector 映射为 Eigen::Vector，零拷贝
        Eigen::Map<const Eigen::VectorXd> b_map(b.data(), b.size());
        Eigen::Map<Eigen::VectorXd> x_map(x.data(), x.size());

        // 使用初猜值 (Warm Start)
        x_map = solver.solveWithGuess(b_map, x_map);

        if (solver.info() != Eigen::Success) {
            std::cerr << "[Eigen] Solve failed! Error code: " << solver.info() << " (NoConvergence=1, NumericalIssue=2)" << std::endl;
            return false;
        }
        return true;
    }
};
