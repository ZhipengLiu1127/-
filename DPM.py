import numpy as np

def propeller(R_total,Cp,Vs,Lpp,B,D,DWT,D_P,n,R_calm,C_B):
    """
    Lpp: 两柱间长
    B: 型宽
    D: 吃水深度 米
    DWT: 载重吨
    Cp :菱形系数
    t:  推力减额系数
    T:  主机推力
    Kt: 推进系数
    kq: 扭矩系数
    J: 进速系数
    D_P: 螺旋桨直径
    n: 螺旋桨转速
    ua: 螺旋桨入流速度
    Vs:船舶对水速度
    omega: 船舶伴流系数
    V_displacement: 排水量 （用载重吨代替）
    rho_s  : 海水密度 (kg/m^3)
    D_P: 螺旋桨直径
    eta_r : 传输效率  尾机型取0.98 中机型取0.97
    eta_o : 螺旋桨效率
    eta_d : 推进效率
    tao_pp : 螺旋桨载荷系数
    Vsc : 修正后的船速
    P_id: 理想状态下的船舶功率（静水中的船舶功率）
    R_total : 总阻力（静水阻力加风浪附加阻力）
    R_calm : 静水阻力
    

    K_T_0 = a_t*J**2+b_t*J+c_t
    K_Q_0 = a_q*J**2+b_q*J+c_q
    """
    #.................................常数计算..............................#
    rho_s = 1.05  # 海水密度 (kg/m^3)
    #V_displacement = DWT
    #C_B = V_displacement/(Lpp*B*D)
    omega = 0.5*C_B-0.05  #盖勒公式
    ua = Vs*(1-omega)
    t = 0.5*Cp-0.12   #汉克谢尔推力减额系数
    J_0 = ua/(n*D_P)
    eta_r = 0.98 
    #................................KT 二次拟合系数.........................#
    a_t = -0.149
    b_t = -0.26
    c_t = 0.2937
    
    #................................KQ 二次拟合系数..........................#
    a_q = -0.1917
    b_q = -0.1788
    c_q = 0.3022
    
    #....................船舶理想状态功率计算..........................................#
    T = R_total/(1-t)
    #print(f"R_total :{R_total}")
    #print(f"T: {T}")
    K_T_0 = a_t*J_0**2+b_t*J_0+c_t  #由二项式拟合得到系数
    K_Q_0 = a_q*J_0**2+b_q*J_0+c_q  #由二项式拟合得到系数
    tao_pp = T/((1-omega)**2*rho_s*Vs**2*D**2)
    P_id = R_calm*Vs
    #print(f"R_calm: {R_calm}")
    #J_1 = (-b_t-np.sqrt(b_t**2-4*(a_t-tao_pp)*c_t))/(2*(a_t-tao_pp))

    #..............................船舶实际功率计算..................................#
    J_1 = (-b_t-np.sqrt(b_t**2-4*(a_t-tao_pp)*c_t))/(2*(a_t-tao_pp))
    K_T_1 = a_t*J_1**2+b_t*J_1+c_t
    K_Q_1 = a_q*J_1**2+b_q*J_1+c_q
    eta_o = (J_1/2*np.pi)*(K_T_1/K_Q_1)
    eta_d = eta_o*eta_r*((1-t)/(1-omega))
    Vsc = (n*J_1*D_P)/(1-omega)  #船舶修正船速

    #................................附加功率计算..................................#
    P_delta = ((R_total-R_calm)*Vsc)/eta_d
    P_d = P_id + P_delta
    return Vs,Vsc
    






