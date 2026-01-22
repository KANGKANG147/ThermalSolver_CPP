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

# ==================== 物理参数 ====================
PARAMS = {
    "Density":  7800.0,   # kg/m3 (Iron)
    "Cp":       500.0,    # J/kgK
    "k":        40.0,     # W/mK (导热率)
    "Thickness": 0.01,    # m (1cm)
    "Area":     1.0,      # m2 (板子面积)
    "L_dist":   1.0,      # 传热距离 (两个质心的距离: 0.5 + 0.5 = 1.0m)
    "L_contact": 1.0,     # 接触边长 (1m)
    "dt":       10.0      # 仿真步长 (s)
}

# 计算质量 (单块板)
Mass = PARAMS["Area"] * PARAMS["Thickness"] * PARAMS["Density"]
# 计算热容 C
C_thermal = Mass * PARAMS["Cp"]

# 计算导热系数 K (Conductance)
# 接触面积 A_c = 边长 * 厚度
A_contact = PARAMS["L_contact"] * PARAMS["Thickness"]
# K = k * A_c / L_dist
K_cond = (PARAMS["k"] * A_contact) / PARAMS["L_dist"]

print(f"质量 m={Mass:.2f} kg, 热容 C={C_thermal:.1f} J/K")
print(f"接触面积={A_contact} m2, 导热系数 K={K_cond:.4f} W/K")

def main():
    plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'Arial Unicode MS']
    plt.rcParams['axes.unicode_minus'] = False

    csv_path = 'D:\\毕业设计\\code_cpp\\ThermalSolver_CPP\\Output\\results.csv'
    if not os.path.exists(csv_path):
        print("未找到 results.csv")
        return

    df = pd.read_csv(csv_path)
    
    # 自动找两块板的温度列
    # 假设列名包含 "HotPart" 和 "ColdPart"
    col_hot = None
    col_cold = None
    for col in df.columns:
        if "HotPart" in col: col_hot = col
        if "ColdPart" in col: col_cold = col
    
    if not col_hot or not col_cold:
        print("错误：在CSV里找不到 HotPart 或 ColdPart 的列")
        print("现有列名:", df.columns)
        return

    sim_hot = df[col_hot].values
    sim_cold = df[col_cold].values
    time_sec = np.arange(len(sim_hot)) * PARAMS['dt']

    # ==================== 理论解析解 ====================
    # 两个物体互传，温差指数衰减
    # T_avg = (100 + 0) / 2 = 50
    # 时间常数 tau = C / (2 * K)  <-- 注意这里是 2K，因为两边都在变
    tau = C_thermal / (2 * K_cond)
    print(f"理论时间常数 tau = {tau:.1f} 秒")

    # T_hot(t) = 50 + 50 * exp(-t/tau)
    # T_cold(t) = 50 - 50 * exp(-t/tau)
    theory_hot = 50.0 + 50.0 * np.exp(-time_sec / tau)
    theory_cold = 50.0 - 50.0 * np.exp(-time_sec / tau)

    # 计算误差
    mae_hot = np.mean(np.abs(sim_hot - theory_hot))
    
    # 画图
    plt.figure(figsize=(10, 6), dpi=150)
    plt.style.use('bmh')
    
    plt.plot(time_sec, theory_hot, 'r--', linewidth=2, label='理论值 (Hot)')
    plt.plot(time_sec, theory_cold, 'b--', linewidth=2, label='理论值 (Cold)')
    plt.plot(time_sec, sim_hot, 'r-', alpha=0.6, linewidth=2, label='仿真值 (Hot)')
    plt.plot(time_sec, sim_cold, 'b-', alpha=0.6, linewidth=2, label='仿真值 (Cold)')
    
    plt.title(f'热传导验证 (绝热接触)\nMAE = {mae_hot:.4f} °C', fontsize=14)
    plt.xlabel('时间 (s)', fontsize=12)
    plt.ylabel('温度 (°C)', fontsize=12)
    plt.legend()
    plt.grid(True)
    
 #   plt.savefig('Output/Validation_Conduction.png')
    plt.show()

if __name__ == "__main__":
    main()