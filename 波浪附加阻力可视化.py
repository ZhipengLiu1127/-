import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import PchipInterpolator

# matplotlib字体设置，解决中文显示问题
plt.rcParams['font.family'] = 'Microsoft YaHei'  # 使用微软雅黑字体
plt.rcParams['axes.unicode_minus'] = False  # 解决坐标轴负数的负号显示问题

file_path = r"corrected_speeds_and_stall_factors_3.csv"  # 请修改为实际的文件路径
data = pd.read_csv(file_path)
# 1. 对数据按 L_over_lambda 排序，并去除空值
data_sorted = data.sort_values(by='L_over_lambda').dropna(subset=['L_over_lambda', 'wave_res'])
x = data_sorted['L_over_lambda'].values
y = data_sorted['wave_res'].values

# 2. 若存在重复 x 值，需要进行分组聚合（例如对同一个 x 取 y 的平均）
df_interp = pd.DataFrame({'x': x, 'y': y})
df_grouped = df_interp.groupby('x', as_index=False).mean()
x_unique = df_grouped['x'].values
y_unique = df_grouped['y'].values

# 3. 判断有效数据量，若足够则进行PCHIP插值
if len(x_unique) >= 2:
    # 在 x 的范围内生成平滑插值节点
    x_smooth = np.linspace(x_unique.min(), x_unique.max(), 300)

    # 使用 PCHIP 进行保形插值
    pchip = PchipInterpolator(x_unique, y_unique)
    y_smooth = pchip(x_smooth)

    # 4. 绘图
    plt.figure(figsize=(8, 5))

    # 原始数据散点
    plt.scatter(x_unique, y_unique, color='blue', alpha=0.5, label='Data points')
    # PCHIP平滑曲线
    plt.plot(x_smooth, y_smooth, color='red', label='PCHIP fit')

    plt.xlabel('L/λ')
    plt.ylabel('$C_{RW}$')
    plt.title('波浪附加阻力与 L/λ 关系')
    plt.legend()
    plt.grid(True)
    plt.show()
else:
    print("数据量过少，无法进行插值。请检查数据中是否包含足够的有效点。")