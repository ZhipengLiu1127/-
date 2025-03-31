import numpy as np
import math
import scipy.integrate as integrate
from wave_spectrum import calculate_combined_spectrum
from scipy.integrate import dblquad

def wave_resistance(Lpp, B, T, Vs, D, T_F, T_A, heading,wave_direction,Hs,DWT,C_B,T_min):
    """
    计算船舶波浪阻力
    
    参数说明:
    heading:船舶航向
    wave_direction: 波浪方向
    Lpp: 两柱间长
    B: 型宽
    D: 吃水深度 米
    DWT: 载重吨
    T: 波浪周期
    T_min: 最小波浪周期
    Fr: 弗汝德数 Fr = Vs / (g * Lpp)**0.5
    Fr_rel: 修正后的弗汝德数 Fr_rel = (Vs - Vg) / (g * Lpp)**0.5
    Vs: 船舶对水速度
    Vg: 入射波群速度
    omega: 波浪频率（单位：弧度/秒）
    omega_p:波谱峰值频率
    C_B: 方块系数
    kyy: 螺距回转半径与Lpp比值
    T_F: 前垂线的吃水，单位：米
    T_A: 后垂线的吃水，单位：米
    alpha: 船首与波向相对方向夹角
    Hs:显著波高 米
    zeta_A:波浪振幅
    V_displacement: 排水量 （用载重吨代替）
    WP: 波浪谱（叠加方向扩散函数版）

    返回:
    R1:长波增阻
    R2:短波增阻
    R: 总波浪阻力（规则波下的总阻力）

    
    """
    
    # 常数定义
    g = 9.8  # 重力加速度 (m/s^2)
    rho_s = 1.05  # 海水密度 (kg/m^3)
    omega = 2*np.pi/T #波浪频率（单位：弧度/秒）
    omega_p = 2*np.pi/T_min  #波谱峰值频率
    zeta_A = Hs/2
    #V_displacement = DWT
   #C_B = V_displacement/(Lpp*B*D)
    kyy = 0.25*Lpp
    Vg = g/(2*omega)
    # 计算弗汝德数
    Fr = Vs / (g * Lpp)**0.5
    Fr_rel = (Vs - Vg) / (g * Lpp)**0.5
    # 计算各项 E1,E2
    E1 = np.arctan(0.495 * B / Lpp)  # 船首的角度
    E2 = np.arctan(0.495 * B / Lpp)  # 船尾的角度
    Vg = (g*T)/(2*np.pi)
    print(f"Fr: {Fr}")
    (f"Fr_rel: {Fr_rel}")
    # 取最大吃水
    T_deep = max(T_F, T_A)
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

            return min(alpha_rad_2, abs(np.pi-alpha_rad_2))
    alpha = calculate_alpha(heading,wave_direction)
    print(f"alpha: {alpha}")
    # 计算长波附加阻力 R1
    def R1(Lpp, Vs, Vg, C_B, omega, B, zeta_A, kyy,T_A,T_F, alpha):
        # 计算omega bar
        omegabar = 2.142 * (kyy ** (1./3.)) * np.sqrt(Lpp / (2 * np.pi * g)) * ((C_B / 0.65) ** 0.17) * \
               (1 - (0.111 / C_B) * (np.log(B / T_deep) - np.log(2.75))) * \
               ((-1.377 * Fr**2 + 1.157 * Fr) * np.abs(np.cos(alpha)) + (0.618 * (13 + np.cos(2 * alpha))) / 14) * omega
        print(f"omegabar: {omegabar}")
        # 计算辅助项
        def a1(alpha, Fr, C_B, B, T_deep):
            if 0 <= alpha <= np.pi / 2:
                # 打印 B / T_deep 的计算结果
                print(f"B: {B}, T_deep: {T_deep}, B / T_deep: {B / T_deep}")

                # 确保 B / T_deep > 0
                if B / T_deep <= 0:
                    print(f"Warning: B / T_deep is less than or equal to 0. B: {B}, T_deep: {T_deep}")
                    return None  # 如果 B / T_deep 小于等于 0，返回 None

                # 打印 np.log(B / T_deep) 的结果
                log_value = np.log(B / T_deep)
                print(f"np.log(B / T_deep): {log_value}")

                # 正常计算 a1
                return (0.87 / C_B) ** ((1 + Fr) * np.cos(alpha)) * log_value ** -1 * ((1 + 2 * np.cos(alpha)) / 3)
            elif alpha == np.pi:
                if Vs > Vg and Fr_rel >= 0.12:
                    return (0.87 / C_B) ** (1 + Fr_rel) * (np.log(B / T_deep)) ** -1
                else:
                    return (0.87 / C_B) * (np.log(B / T_deep)) ** -1

        def a2(alpha, Fr, Vs, Vg, Fr_rel):
            if 0 <= alpha <= np.pi / 2:
                if Fr < 0.12:
                    return 0.0072 + 0.1676 * Fr
                else:
                    return Fr ** 1.5 * np.exp(-3.5 * Fr)
            elif alpha == np.pi:
                if Vs <= Vg:
                    return 0.0072 * (2 * Vs / Vg - 1)
            elif Vs > Vg and Fr_rel < 0.12:
                return 0.0072 + 0.1676 * Fr_rel
            else:
                return Fr_rel ** 1.5 * np.exp(-3.5 * Fr_rel)

        a3 = 1.0 + 28.7 * np.arctan(np.abs(T_A - T_F) / Lpp)
        print(f"a3: {a3}")
        
        b1 = 11.0 if omegabar < 1 else -8.5
        print(f"b1: {b1}")
        
        def d1(Lpp, B, C_B, T_A, T_F):
            if omegabar < 1:
                return 566 * ((Lpp * C_B) / B) ** (-2.66)
            else:
                return (-1)*566 * (Lpp / B) ** (-2.66) * (4 - 125 * np.arctan(np.abs(T_A - T_F) / Lpp))
            #print(f"d1: {d1(Lpp, B, C_B, T_A, T_F)}")
            # 调试打印 a1, a2 和 d1 的结果
                # 计算长波附加阻力 R1
        a1_value = a1(alpha, Fr, C_B, B, T_deep)
        a2_value = a2(alpha, Fr, Vs, Vg, Fr_rel)
        d1_value = d1(Lpp, B, C_B, T_A, T_F)

                # 检查a1, a2, d1的值是否有效
        print(f"a1_value: {a1_value}, a2_value: {a2_value}, d1_value: {d1_value}")

                # 如果任何参数无效，打印错误信息并跳过计算
        if not np.isfinite(a1_value) or not np.isfinite(a2_value) or not np.isfinite(d1_value):
            print("Error: Invalid values encountered in a1, a2, or d1.")
            return None
        # 计算 R1 长波附加阻力
        #R1_value = 3859.2 * rho_s * g * zeta_A**2 * (B**2 / Lpp) * C_B**1.34 * (kyy**2) * \
                   #a1(alpha, Fr, C_B, B, T_deep) * a2(alpha, Fr, Vs, Vg, Fr_rel) * a3 *\
                    #omegabar**b1 * np.exp(b1 / d1(Lpp, B, C_B, T_A, T_F) * (1 - omegabar ** d1(Lpp, B, C_B, T_A, T_F)))
        R1_value = 3859.2 * rho_s * g * zeta_A**2 * (B**2 / Lpp) * C_B**1.34 * (kyy**2) *\
                    a1_value*a2_value*a3*(omegabar**b1) * np.exp((b1 / d1_value)*(1 - omegabar ** d1_value))
        print(f"R1_value: {R1_value}")
        return R1_value

    # 计算短波附加阻力 R2
    def R2(Lpp, B, zeta_A, omega, Vs, C_B, Fr, alpha, E1, E2, T_deep, g):

        # Calculate the draft coefficient α_T*
        def alpha_T_star(T_star, Lpp, omega, g):
            term = (2 * np.pi * g) / (omega**2 * Lpp)
            if term <= 2.5:
                return 1 - np.exp(-4 * np.pi * (T_star / (2 * np.pi * g / omega**2) - T_star / (2.5 * Lpp)))
            else:
                return 0
        #计算f_a
        def f_a(alpha,E1):
            if 0<=alpha<=E1:
                f_a = np.cos(alpha)
                return f_a
            else:
                return 0
        # 计算各项 R_AWR_1, R_AWR_2, R_AWR_3, R_AWR_4
        def R_AWR_1(rho_s, g, B, zeta_A, alpha_T_star, omega, V_s, C_B, Fr, alpha, E1):
            
            term_1 = np.sin(E1 + alpha)**2 + (2 * omega * V_s / g) * (np.cos(alpha) - np.cos(E1) * np.cos(E1 + alpha))
            term_2 = (0.87 / C_B) ** ((1 + 4 * np.sqrt(Fr)) *  f_a(alpha, E1))
            return (2.25 / 4) * rho_s * g * B * zeta_A**2 * alpha_T_star * term_1 * term_2 

        def R_AWR_2(rho_s, g, B, zeta_A, alpha_T_star, omega, V_s, C_B, Fr, alpha, E1):
            
            term_1 = np.sin(E1 - alpha)**2 + (2 * omega * V_s / g) * (np.cos(alpha) - np.cos(E1) * np.cos(E1 - alpha))
            term_2 = (0.87 / C_B) ** ((1 + 4 * np.sqrt(Fr)) *  f_a(alpha, E1))
            return (2.25 / 4) * rho_s * g * B * zeta_A**2 * alpha_T_star * term_1 * term_2
        
        def R_AWR_3(rho_s, g, B, zeta_A, alpha_T_star, omega, V_s, alpha, E2):
            term_1 = np.sin(E2 - alpha)**2 + (2 * omega * V_s / g) * (np.cos(alpha) - np.cos(E2) * np.cos(E2 - alpha))
            return -(2.25 / 4) * rho_s * g * B * zeta_A**2 * alpha_T_star * term_1

        def R_AWR_4(rho_s, g, B, zeta_A, alpha_T_star, omega, V_s, alpha, E2):
            term_1 = np.sin(E2 + alpha)**2 + (2 * omega * V_s / g) * (np.cos(alpha) - np.cos(E2) * np.cos(E2 + alpha))
            return -(2.25 / 4) * rho_s * g * B * zeta_A**2 * alpha_T_star * term_1

        
        # 计算阻力
        alpha_T_star_1 = alpha_T_star(T_deep, Lpp, omega, g)
        R_AWR_1_val = R_AWR_1(rho_s, g, B, zeta_A, alpha_T_star_1, omega, Vs, C_B, Fr, alpha, E1)
        R_AWR_2_val = R_AWR_2(rho_s, g, B, zeta_A, alpha_T_star_1, omega, Vs, C_B, Fr, alpha, E1)


        
        if C_B <= 0.75:
            T_star_3_4 = T_deep * (4 + np.sqrt(np.abs(np.cos(alpha)))) / 5
        else:
            T_star_3_4 = T_deep * (2 + np.sqrt(np.abs(np.cos(alpha)))) / 3
        alpha_T_star_3_4 = alpha_T_star(T_star_3_4, Lpp, omega, g)
        
        R_AWR_3_val = R_AWR_3(rho_s, g, B, zeta_A, alpha_T_star_3_4, omega, Vs, alpha, E2)
        R_AWR_4_val = R_AWR_4(rho_s, g, B, zeta_A, alpha_T_star_3_4, omega, Vs, alpha, E2)
        
        R2_value = R_AWR_1_val + R_AWR_2_val + R_AWR_3_val + R_AWR_4_val
        print(f"R2_value: {R2_value}")
        return R2_value
    
    # 计算 R1 和 R2 的总阻力
    print(f"T_F: {T_F}, T_A: {T_A}, alpha: {alpha}, omega: {omega}, C_B: {C_B}, kyy: {kyy}")

    R1_value = R1(Lpp, Vs, Vg, C_B, omega, B,zeta_A, kyy, T_A, T_F, alpha)
    R2_value = R2(Lpp, B,zeta_A, omega, Vs, C_B, Fr, alpha, E1, E2, T_deep,g)
    
    # 总波浪阻力
    total_R = (R1_value + R2_value)
    print(f"total_R: {total_R}")
    
    test_spectrum = calculate_combined_spectrum(omega, omega_p, T, Hs, alpha, wave_direction) #叠加海浪谱（有方向扩散函数版）
    print(f"test_spectrum: {test_spectrum}")

    #total_R_1 =


    #无量纲化处理
    R3 = total_R/(rho_s*g*(zeta_A**2)*(B**2)/Lpp)
    print(f"R3: {R3}")
    
    return total_R