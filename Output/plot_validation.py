import os
import sys
import numpy as np
try:
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib
except ImportError:
    print("缺库，请运行: pip install pandas matplotlib numpy")
    sys.exit(1)

# ==================== 物理参数 (必须一致) ====================
CONFIG = {
    "T_init_K": 1000.0,  # 初始 1000 K (726.85 C)
    "m":        10.0,    # 质量 (1x1x10)
    "Cp":       100.0,   # 比热
    "A":        1.0,     # 面积 (单面辐射因为另一面绝热? 不，通常是两面)
                         # 你的 config 里前后都 INSULATED，
                         # 但代码里辐射是算 Front 的，所以面积是 1.0
    "epsilon":  1.0,     # 发射率
    "sigma":    5.67e-8, # 斯蒂芬-玻尔兹曼常数
    "dt":       1.0
}

def main():
    # 设置中文
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False

    # 1. 读取 CSV (自动找)
    csv_path = 'D:\\毕业设计\\code_cpp\\ThermalSolver_CPP\\Output\\results.csv'
    if not os.path.exists(csv_path):
        print("找不到 results.csv")
        return
    df = pd.read_csv(csv_path)

    # 2. 找温度列 (最后一列或者 TestGroup)
    temp_col = df.columns[-1]
    for col in df.columns:
        if "TestGroup" in col: temp_col = col; break
    
    print(f"分析列: {temp_col}")
    sim_temp_c = df[temp_col].values
    sim_temp_k = sim_temp_c + 273.15 # 换算成开尔文
    
    time_sec = np.arange(len(sim_temp_k)) * CONFIG['dt']

    # 3. 计算理论解 (深空冷却公式)
    # T(t) = [ T0^-3 + (3 * eps * sigma * A / m * Cp) * t ] ^ (-1/3)
    
    K = (3 * CONFIG['epsilon'] * CONFIG['sigma'] * CONFIG['A']) / (CONFIG['m'] * CONFIG['Cp'])
    print(f"辐射时间系数 K = {K:.4e}")
    
    # 理论计算
    term = np.power(CONFIG['T_init_K'], -3) + K * time_sec
    theory_temp_k = np.power(term, -1/3)
    theory_temp_c = theory_temp_k - 273.15

    # 4. 截取前 600 秒
    limit = 600
    t_plot = time_sec[:limit]
    sim_plot = sim_temp_c[:limit]
    theory_plot = theory_temp_c[:limit]

    # 5. 计算误差
    mae = np.mean(np.abs(sim_plot - theory_plot))
    print(f"MAE = {mae:.4f} °C")

    # 6. 画图
    plt.figure(figsize=(10, 6), dpi=120)
    plt.plot(t_plot, theory_plot, 'r--', linewidth=2, label='理论解析解 (Stefan-Boltzmann)')
    plt.plot(t_plot, sim_plot, 'b-', linewidth=2, alpha=0.7, label='仿真结果')
    plt.title(f'辐射冷却验证 (初始1000K)\n平均误差 MAE = {mae:.4f} °C')
    plt.xlabel('时间 (s)')
    plt.ylabel('温度 (°C)')
    plt.grid(True)
    plt.legend()
    plt.savefig('Validation_Radiation.png')
    plt.show()

if __name__ == "__main__":
    main()