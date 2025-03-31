import pandas as pd
import numpy as np
from calme_water import calculate_clame_resistance
from DPM import propeller
from wind import wind_resistance
from irregular_wave import wave_resistance_spectrum
import matplotlib.pyplot as plt
#from scipy.interpolate import PchipInterpolator

# 读取数据文件
file_path = r"C:\Users\15585\Desktop\小论文写作素材\船舶失速航线优化\matched_data_3.csv"  # 请修改为实际的文件路径
data = pd.read_csv(file_path)


# matplotlib字体设置，解决中文显示问题
plt.rcParams['font.family'] = 'Microsoft YaHei'  # 使用微软雅黑字体
plt.rcParams['axes.unicode_minus'] = False  # 解决坐标轴负数的负号显示问题

# 假设已经有输入参数
D_P = 5.25  # 螺旋桨直径
mode_input = "Heavy Ballast, Cranes"  # 载重状态
DWT = 28280  # 载重吨，单位为t
B = 27.2  # 型宽
T = 9.82  # 吃水深度
Lpp = 160.4  # 两柱间长
Cp = 0.78  # 菱形系数
C_B = 0.77 #方形系数
g = 9.8  # 重力加速度 (m/s^2)

def process_data(row, DWT, B, T, Lpp, Cp, D_P, mode_input,C_B):
    # 检查输入的参数是否有空值
    if row['speed'] is None or row['wind_u'] is None or row['wind_v'] is None or row['heading'] is None:
        raise ValueError("输入的CSV数据缺少必要的列（如speed, wind_u, wind_v, heading）")

    # 1. 将船速从节 (knots) 转换为米/秒 (m/s)
    speed_mps = row['speed'] * 0.514444  # 转换公式：1 节 = 0.514444 米/秒

    # 2. 计算静水阻力
    R_calm = calculate_clame_resistance(Lpp, B, T, speed_mps,C_B)  # 使用转换后的船速

    # 3. 计算风阻力
    wind_res = wind_resistance(
        row['wind_u'],  # wind_u 列对应 u
        row['wind_v'],  # wind_v 列对应 v
        row['heading'],  # 传递 heading
        speed_mps,  # 使用转换后的船速
        mode_input,  # 传递载重状态
        B,  # 型宽
        T  # 吃水深度
    )

    # 4. 计算波浪阻力
    if row['wave_period'] is None or row['wave_height'] is None or row['wave_direction'] is None:
        raise ValueError("输入的CSV数据缺少波浪相关列（如wave_period, wave_height, wave_direction）")

    T_F = T  # 使用船舶吃水深度T作为前垂线的吃水
    T_A = T  # 使用船舶吃水深度T作为后垂线的吃水
    T_min = data['wave_period'].min() #获取最小周期
    T_max = data['wave_period'].max() #获取最大波浪周期

    #print(
        #f"T_F: {T_F}, T_A: {T_A}, wave_period: {row['wave_period']}, wave_height: {row['wave_height']}, wave_direction: {row['wave_direction']},T_min: {T_min}")

    wave_res = wave_resistance_spectrum(
        Lpp,
        B,
        speed_mps,  # 使用转换后的船速
        T,  # D（深度）
        T_F,
        T_A,
        C_B,
        row['wave_period'],
        T_min,
        T_max,
        DWT,  # 载重吨
        row['wave_height'],  # 传递 Hs（显著波高）
        row['wave_direction'],  # 传递波浪方向
        row['heading'],  # 传递航向 heading

    )

    # 5. 计算总阻力
    R_total = R_calm + wind_res + (wave_res)
    #print(f"R_total: {R_total}")

    # 6. 使用直接功率法计算修正后的船速
    n = row['rpm']  # 从CSV文件的rpm列获取主机转速
    Vs, Vsc = propeller(R_total, Cp, speed_mps, Lpp, B, T, DWT, D_P, n, R_calm, C_B)

    # 7. 计算失速系数，假设失速系数为修正速度与最大理论速度的比值
    stall_factor = Vsc / Vs

    # 8. 计算 L/λ（波长计算）
    wavelength = (g * row['wave_period'] ** 2) / (2 * np.pi)  # 计算波长 λ
    L_over_lambda = Lpp / wavelength  # 计算 L/λ

    return pd.Series([Vs, Vsc, stall_factor, L_over_lambda, wave_res],
                     index=['original_speed', 'corrected_speed', 'stall_factor', 'L_over_lambda', 'wave_res'])


# 对每一行数据应用计算过程，并将结果作为新列加入data中
data[['original_speed', 'corrected_speed', 'stall_factor', 'L_over_lambda', 'wave_res']] = data.apply(
    lambda row: process_data(row, DWT, B, T, Lpp, Cp, D_P, mode_input, C_B), axis=1)

# 输出部分结果
print(data[['timestamp', 'Lat', 'Long', 'heading', 'speed', 'corrected_speed', 'stall_factor']])

# 保存计算结果到新CSV文件
data.to_csv('corrected_speeds_and_stall_factors_3.csv', index=False)




