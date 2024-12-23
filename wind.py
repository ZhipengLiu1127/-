import numpy as np
import math 
#from math import sqrt 
#from math import cos
import scipy.integrate as integrate
def wind_resistance(u,v,heading,V,mode_input,B,D):
    """
    计算风附加阻力
    u, v:ERA5 风u,v数据
    Cw: 风阻系数
    p: 空气密度
    B:船舶型宽
    D: 船舶吃水
    At: 船体水线上部分，横剖面上的投影面积
    Vw: 风与船的相对速度
    heading:船舶航向
    V:船舶速度
    mode_input:载重状态
    """
     # 风阻系数表
    wind_resistance_data = {
        "Angle of attack [°]": [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180],
        "Heavy Ballast, Cranes": [-0.75, -0.91, -0.94, -0.87, -0.80, -0.59, -0.34, -0.18, 0.02, 0.11, 0.10, 0.22, 0.44, 0.65, 0.78, 0.81, 0.85, 0.82, 0.70],
        "Ballast, Cranes": [-0.71, -0.84, -0.89, -0.82, -0.70, -0.49, -0.23, -0.12, 0.01, -0.02, 0.02, 0.12, 0.36, 0.55, 0.71, 0.76, 0.79, 0.80, 0.68]
    }

    p = 1.226 #空气密度单位 kg/m**3

    At=0.6*B*D #经验公式

    ship_speed = V
    # 计算风速和风向
    true_wind_speed = math.sqrt(u**2 + v**2)
    true_wind_direction = (180 + math.degrees(math.atan2(-u, -v))) % 360  # 调整到 0–360 范围
    # 计算相对风向角
    relative_wind_angle = true_wind_direction - heading
    if relative_wind_angle < 0:
        relative_wind_angle += 360  # 确保角度在0到360之间
    elif relative_wind_angle >= 360:
        relative_wind_angle -= 360

    # 将相对风向角调整到 0–180（对称性）
    relative_wind_angle = min(relative_wind_angle, 360 - relative_wind_angle)

    # 计算相对风速（使用矢量合成法）
    Vw = math.sqrt(true_wind_speed**2 + ship_speed**2 - 2 * true_wind_speed * ship_speed * math.cos(math.radians(relative_wind_angle)))
    
    # 查找风阻系数
    angle_list = wind_resistance_data["Angle of attack [°]"]
    coefficients = wind_resistance_data[mode_input]
    if relative_wind_angle in angle_list:
        Cw = coefficients[angle_list.index(relative_wind_angle)]
    else:
        # 线性插值计算
        for i in range(len(angle_list) - 1):
            if angle_list[i] < relative_wind_angle < angle_list[i + 1]:
                x0, x1 = angle_list[i], angle_list[i + 1]
                y0, y1 = coefficients[i], coefficients[i + 1]
                Cw = y0 + (y1 - y0) * (relative_wind_angle - x0) / (x1 - x0)
                break

    Rw = 0.5*Cw*At*p*Vw**2
    return Rw