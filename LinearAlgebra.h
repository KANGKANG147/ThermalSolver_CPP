#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <iomanip>

// 简单的向量定义
using Vector = std::vector<double>;

// CSR 格式稀疏矩阵 (Compressed Sparse Row)
// 用于存储方程组 Ax = b 中的 A
class SparseMatrix {
public:
    int rows = 0;
    int cols = 0;

    // CSR 三大数组
    std::vector<double> values;     // 非零元素的值
    std::vector<int> col_indices;   // 对应的列索引
    std::vector<int> row_ptr;       // 每一行的起始位置

    SparseMatrix() {}
    SparseMatrix(int r, int c) : rows(r), cols(c) {
        row_ptr.resize(rows + 1, 0);
    }

    // 矩阵-向量乘法: y = A * x
    // 这是所有迭代法(包括AMG)中最频繁的操作
    void multiply(const Vector& x, Vector& y) const {
        assert(x.size() == cols);
        if (y.size() != rows) y.resize(rows);
        for (int i = 0; i < rows; ++i) {
            double sum = 0.0;
            for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                sum += values[k] * x[col_indices[k]];
            }
            y[i] = sum;
        }
    }

    // 2. 矩阵转置 (用于构建限制算子 R = P^T)
    SparseMatrix transpose() const {
        SparseMatrix T(cols, rows);
        std::vector<int> count(cols, 0);

        // 统计每一列有多少非零元素
        for (int idx : col_indices) count[idx]++;

        // 构建转置矩阵的 row_ptr
        T.row_ptr[0] = 0;
        for (int i = 0; i < cols; ++i) {
            T.row_ptr[i + 1] = T.row_ptr[i] + count[i];
        }

        // 填充数据
        std::vector<int> current_pos = T.row_ptr;
        for (int i = 0; i < rows; ++i) {
            for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                int col = col_indices[k];
                double val = values[k];

                int dest = current_pos[col];
                T.values.resize(T.row_ptr.back()); // 确保大小够
                T.col_indices.resize(T.row_ptr.back());

                T.values[dest] = val;
                T.col_indices[dest] = i; // 原来的行号变成列号
                current_pos[col]++;
            }
        }
        return T;
    }

    // 3. 矩阵-矩阵乘法 C = A * B (SpGEMM) - AMG 的核心难点
    // 用于计算粗网格算子 A_coarse = R * A * P
    static SparseMatrix multiply(const SparseMatrix& A, const SparseMatrix& B) {
        assert(A.cols == B.rows);
        SparseMatrix C(A.rows, B.cols);

        // 临时变量用于累加一行
        std::vector<double> dense_row(B.cols, 0.0);
        std::vector<int> nnz_cols;
        std::vector<bool> is_set(B.cols, false);

        C.row_ptr[0] = 0;
        for (int i = 0; i < A.rows; ++i) {
            // 计算 C 的第 i 行 = A 的第 i 行 * B
            nnz_cols.clear();

            for (int ka = A.row_ptr[i]; ka < A.row_ptr[i + 1]; ++ka) {
                int col_a = A.col_indices[ka]; // B 的行
                double val_a = A.values[ka];

                for (int kb = B.row_ptr[col_a]; kb < B.row_ptr[col_a + 1]; ++kb) {
                    int col_b = B.col_indices[kb];
                    double val_b = B.values[kb];

                    if (!is_set[col_b]) {
                        is_set[col_b] = true;
                        nnz_cols.push_back(col_b);
                    }
                    dense_row[col_b] += val_a * val_b;
                }
            }

            // 排序并存入 CSR
            std::sort(nnz_cols.begin(), nnz_cols.end());
            for (int col : nnz_cols) {
                if (std::abs(dense_row[col]) > 1e-20) { // 过滤数值零
                    C.values.push_back(dense_row[col]);
                    C.col_indices.push_back(col);
                }
                // 清理临时状态
                dense_row[col] = 0.0;
                is_set[col] = false;
            }
            C.row_ptr[i + 1] = C.values.size();
        }
        return C;
    }

    // 4. Gauss-Seidel 平滑器 (Smoother)
    void smooth_gauss_seidel(const Vector& b, Vector& x, int iterations) const {
        for (int iter = 0; iter < iterations; ++iter) {
            for (int i = 0; i < rows; ++i) {
                double sigma = 0.0;
                double diag = 0.0;

                for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                    int col = col_indices[k];
                    if (col == i) diag = values[k];
                    else sigma += values[k] * x[col];
                }

                if (std::abs(diag) > 1e-20)
                    x[i] = (b[i] - sigma) / diag;
            }
        }
    }

    // 5. 计算残差 r = b - Ax
    void residual(const Vector& b, const Vector& x, Vector& r) const {
        if (r.size() != rows) r.resize(rows);
        for (int i = 0; i < rows; ++i) {
            double ax = 0.0;
            for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                ax += values[k] * x[col_indices[k]];
            }
            r[i] = b[i] - ax;
        }
    }

    // 6. 获取对角线元素
    double get_diag(int row) const {
        for (int k = row_ptr[row]; k < row_ptr[row + 1]; ++k) {
            if (col_indices[k] == row) return values[k];
        }
        return 1.0; // Fallback
    }
};

class MatrixBuilder {
private:
    struct Entry { int col; double val; };
    std::vector<std::vector<Entry>> rows_data;
    int n_rows;

public:
    MatrixBuilder(int n) : n_rows(n) { rows_data.resize(n); }

    void add(int r, int c, double v) {
        if (r >= 0 && r < n_rows) {
            rows_data[r].push_back({ c, v });
        }
        else {
            std::cerr << "[Error] MatrixBuilder add out of bounds: row " << r << std::endl;
        }
    }

    // 构建 CSR 矩阵
    SparseMatrix build() {
        SparseMatrix mat(n_rows, n_rows);
        mat.cols = 0; // 重新计算最大列

        for (int i = 0; i < n_rows; ++i) {
            // 1. 获取当前行的数据引用
            auto& row_vec = rows_data[i];
            // 2. 按列号排序 (CSR 要求行内有序)
            std::sort(row_vec.begin(), row_vec.end(), [](const Entry& a, const Entry& b) {
                return a.col < b.col;
            });

            // 3. 填充 CSR 数组并合并重复项
            for (size_t k = 0; k < row_vec.size(); ++k) {
                // 如果是当前行的最后一个元素，或者下一个元素的列号不同 -> 写入
                // 如果下一个元素列号相同 -> 累加到当前元素，稍后写入
                if (k + 1 < row_vec.size() && row_vec[k].col == row_vec[k + 1].col) {
                    row_vec[k + 1].val += row_vec[k].val;
                    continue; // 跳过当前，由下一个写入
                }
                // 写入真实数据
                mat.values.push_back(row_vec[k].val);
                mat.col_indices.push_back(row_vec[k].col);

                // 记录最大列号
                if (row_vec[k].col >= mat.cols) {
                    mat.cols = row_vec[k].col + 1;
                }
            }
            // 更新 row_ptr
            mat.row_ptr[i + 1] = mat.values.size();
        }
        return mat;
    }
};
