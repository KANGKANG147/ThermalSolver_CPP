import os
import sys
import numpy as np

try:
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib
except ImportError:
    print("缺库")
    sys.exit(1)

# ==================== 论文级物理参数 ====================
PARAMS = {
    "Material": "Iron (铁)",
    "Density":  7800.0,   # kg/m3
    "Cp":       500.0,    # J/kgK
    "Thickness": 0.01,    # m (1cm)
    "Area":     1.0,      # m2 (单面散热)
    "h":        10.0,     # W/m2K
    "T_init":   100.0,    # C
    "T_env":    20.0,     # C
    "dt":       1.0       # s
}

# 自动计算质量
PARAMS["Mass"] = PARAMS["Area"] * PARAMS["Thickness"] * PARAMS["Density"]

def main():
    # 设置中文字体
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False

    csv_path = 'D:\\毕业设计\\code_cpp\\ThermalSolver_CPP\\Output\\results.csv'
    if not os.path.exists(csv_path):
        print("未找到 results.csv")
        return

    df = pd.read_csv(csv_path)
    
    # 自动找温度列
    temp_col = df.columns[-1]
    for col in df.columns:
        if "TestGroup" in col or "Iron" in col: temp_col = col; break
    
    print(f"当前分析列: {temp_col}")
    print(f"计算质量: {PARAMS['Mass']:.2f} kg")
    
    # 时间常数 Tau
    tau = (PARAMS['Mass'] * PARAMS['Cp']) / (PARAMS['h'] * PARAMS['Area'])
    print(f"理论时间常数 (Tau): {tau:.2f} 秒")

    sim_temp = df[temp_col].values
    time_sec = np.arange(len(sim_temp)) * PARAMS['dt']

    # ==================== 理论解析解 (牛顿冷却) ====================
    # T(t) = T_env + (T0 - T_env) * exp(-t / tau)
    theory_temp = PARAMS['T_env'] + (PARAMS['T_init'] - PARAMS['T_env']) * np.exp(-time_sec / tau)

    # 截取前 5000 秒 (看前 1.5 小时即可)
    limit = 3600
    if len(time_sec) > limit:
        t_plot = time_sec[:limit]
        sim_plot = sim_temp[:limit]
        theory_plot = theory_temp[:limit]
    else:
        t_plot = time_sec
        sim_plot = sim_temp
        theory_plot = theory_temp

    # 计算误差
    mae = np.mean(np.abs(sim_plot - theory_plot))
    print(f"========================================")
    print(f"平均绝对误差 (MAE): {mae:.4f} °C")
    print(f"========================================")

    # 画图
    plt.figure(figsize=(10, 6), dpi=150)
    plt.style.use('bmh')
    
    plt.plot(t_plot, theory_plot, 'r--', linewidth=2.5, label='理论解析解')
    plt.plot(t_plot, sim_plot, 'b-', linewidth=2.0, alpha=0.7, label='本文求解器仿真结果')
    
    plt.title(f'对流验证 (铁板 1cm, 无辐射)\nMAE = {mae:.4f} °C', fontsize=14)
    plt.xlabel('时间 (s)', fontsize=12)
    plt.ylabel('温度 (°C)', fontsize=12)
    plt.legend(fontsize=12)
    plt.grid(True)
    
 #   plt.savefig('Output/Validation_Convection_Thesis.png')
    plt.show()

if __name__ == "__main__":
    main()