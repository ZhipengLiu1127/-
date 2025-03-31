import numpy as np
#import math
from wave import wave_resistance
from wave_spectrum import calculate_combined_spectrum
 
def wave_resistance_spectrum(Lpp, B, Vs, D, T_F, T_A, C_B,
                                T,T_min,T_max,DWT,
                             Hs, wave_direction,heading):
    """
    使用 ITTC 标准计算不规则波的附加阻力：
    R_total = ∫∫ [ R_wave(ω, α) * S(ω, α) / zeta_A^2 ] dα dω。
    输入数据：
    heading:船舶航向
    wave_direction: 波浪方向
    Lpp: 两柱间长
    B: 型宽
    D: 吃水深度 米
    DWT: 载重吨
    T: 波浪周期
    T_min: 最小波浪周期 
    T_max: 最大波浪周期
    Vs: 船舶对水速度
    Vg: 入射波群速度
    omega: 波浪频率（单位：弧度/秒）
    omega_p:波谱峰值频率
    omega_l： 波谱最低频率
    C_B: 方块系数
    T_F: 前垂线的吃水，单位：米
    T_A: 后垂线的吃水，单位：米
    alpha: 船首与波向相对方向夹角

    Hs:显著波高 米
    zeta_A:波浪振幅
    V_displacement: 排水量 （用载重吨代替）
    WP: 波浪谱（叠加方向扩散函数版）
    d_omega： 频率步长
    d_alpha:方向步长
    

    输出结果：
    

    """
    g = 9.8  # 重力加速度 (m/s^2)
    rho_s = 1050  # 海水密度 (kg/m^3)
    def calculate_alpha(heading, wave_direction):
        # 将波浪传播方向转换为波浪来向
        wave_arrival_direction = (wave_direction + 180) % 360

        # 计算船舶航向与波浪来向的角度差
        alpha = heading - wave_arrival_direction

        # 确保角度在0到360度之间
        alpha = (alpha + 360) % 360  # 如果使用度数

        # 如果使用弧度制，可以转换为弧度并确保在0到2pi之间
        alpha_rad = np.deg2rad(alpha)  # 转换为弧度
        alpha_rad_1 = (alpha_rad + 2 * np.pi) % (2 * np.pi)  # 确保在0到2pi之间
        # 如果 α 超过 π，则取其补角，确保返回值在 0 到 π 之间
        if alpha_rad_1 > np.pi:
            alpha_rad_2 = 2 * np.pi - alpha_rad_1
        else:
            alpha_rad_2 = alpha_rad_1

        return np.deg2rad(alpha)

    alpha = calculate_alpha(heading, wave_direction)
   

    omega = 2 * np.pi / T  # 波浪频率（单位：弧度/秒）
    omega_p = 2 * np.pi / T_min  # 波谱峰值频率
    omega_l = 2 * np.pi / T_max   #波谱最低频率
    alpha_min = 0
    alpha_max = 2*np.pi
    # 手动设置方向步长
    d_alpha = np.pi / 18  # 10° 一步
    alpha_array = np.arange(alpha_min, alpha_max + d_alpha, d_alpha)

    # 手动设置频率步长
    d_omega = omega_p / 10  #
    omega_array = np.arange(omega_l, omega_p + d_omega, d_omega)
    R_sum = 0.0

    for omega in omega_array:
        for alpha in alpha_array:
            S_val, zeta_A = calculate_combined_spectrum(
                omega=omega,
                omega_p=omega_p,
                T = T,
                Hs=Hs,
                alpha=alpha,
                theta_m=wave_direction,
                d_omega=d_omega
            )

            R_single = wave_resistance(
                Lpp=Lpp, B=B, Vs=Vs, T=T,
               heading=heading,wave_direction=wave_direction,
                T_max=T_max,T_min=T_min,D=D, T_F=T_F, T_A=T_A,
                C_B=C_B,zeta_A=zeta_A,Hs=Hs,DWT=DWT
            )

            R_sum += (R_single * S_val / zeta_A**2) * d_alpha

    R_irregular =2* R_sum * d_omega
    #print(f"R_irregular: {R_irregular}")

    # 无量纲化处理
    R3 = R_irregular / (rho_s * g * (zeta_A ** 2) *((B ** 2) / Lpp))
    #x = Lpp/(g*T**2/(2*np.pi))
    #print(f"R3: {R3}")
    return R3
