import numpy as np

def calculate_combined_spectrum(omega, omega_p, T, Hs, alpha, theta_m,d_omega):
    """
    计算综合波浪谱，结合频率谱与方向分布函数。

    参数定义：
    omega: 波浪频率
    omega_p: 波谱峰值频率
    T: 波浪特征周期 (平均周期 T)
    Hs: 三一波高
    alpha: 方向角度（通常是0到2π之间）
    theta_m: 波浪主方向(ERA5波向值替代)
    s: 方向扩展参数，通常为1（风浪）或75（长周期波） #此处取s=1
    d_omega： 频率步长

    返回：
    综合波浪谱值
    zeta_A
    """
    s = 1
    
    # JONSWAP波浪谱计算，返回谱幅值
    def Jonswap(omega, omega_p, T, Hs):
        """
        JONSWAP波浪谱计算，返回谱幅值。

        参数定义：
        omega: 波浪频率
        omega_p: 波谱峰值频率
        T: 波浪特征周期 (平均周期 T)
        Hs: 三一波高

        返回：
        波浪谱的幅值
        """
        if omega <= omega_p:
            delta_p = 0.07
        else:
            delta_p = 0.09

        # 波浪谱部分（频率谱）
        g = -((0.206 * omega * T - 1) / delta_p * np.sqrt(2))**2
        spectrum = ((114 * Hs**2) / ((T**4) * (omega**5))) * np.exp((-691) / (T**4) * (omega**4)) * 3.3**np.exp(g)

        return spectrum
    #ISSC波浪谱计算
    def ISSC(omega,T,Hs):
        """
        ISSC波谱计算
        参数定义：
        omega: 波浪频率
        T: 波浪特征周期 (平均周期 T)
        Hs: 三一波高

        """
        A = 173*Hs**2
        B = (T**4)*(omega**5)
        C = (-691/(T**4)*(omega**4))

        spectrum = A/B*np.exp(C)
        return spectrum

    def direction_distribution(alpha_deg, theta_m_deg, s):
        """
        安全归一化方向谱 G(α)，确保返回 float 值，不为 None。
        """
        # 单位转换：角度 -> 弧度
        alpha = np.radians(alpha_deg)
        theta_m = np.radians(theta_m_deg)

        # 方向差（弧度，限制在 -π ~ π）
        angle_diff = (alpha - theta_m + np.pi) % (2 * np.pi) - np.pi

        # 安全计算 cos(angle_diff)
        cos_val = np.cos(angle_diff)
        cos_val = np.clip(cos_val, -1.0, 1.0)  # 限制在 [-1, 1]

        # 原始方向谱值（强制为非负）
        G_raw = max(cos_val ** (2 * s), 0.0)

        # 归一化方向谱（积分值为1）
        angles = np.linspace(0, 2 * np.pi, 360)
        cos_vals = np.cos(angles - np.radians(theta_m_deg))
        cos_vals = np.clip(cos_vals, -1.0, 1.0)
        G_all = np.power(cos_vals, 2 * s)
        G_all[G_all < 0] = 0
        normalization = np.trapz(G_all, angles)

        # 防止除零
        if normalization == 0:
            return 0.0

        G_alpha = G_raw / normalization
        return float(G_alpha)  # 确保类型是 float，不是 None

    # 计算频率谱
    #spectrum = Jonswap(omega, omega_p, T, Hs)
    spectrum = ISSC(omega,T,Hs)
    #print(f"spectrum: {spectrum}")

    # 计算方向分布函数
    G_alpha = direction_distribution(alpha, theta_m, s)
    #print(f"G_alpha: {G_alpha}")

    if G_alpha is None:
        print("警告：G_alpha 计算失败，alpha=", alpha, "theta_m_deg=", theta_m)
    # 计算综合谱值，频率谱和方向分布函数的乘积
    combined_spectrum = spectrum * G_alpha
    #print(f"combined_spectrum: { combined_spectrum}")

    # **计算 zeta_A**
    zeta_A = np.sqrt(2.0 * combined_spectrum * d_omega)
    #print(f"zeta_A : {zeta_A }")
    return combined_spectrum,zeta_A
