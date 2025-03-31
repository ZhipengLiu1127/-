import numpy as np
import matplotlib.pyplot as plt

# 读取提取的曲线数据
data_kt = np.loadtxt("curve_1_actual_data.csv", delimiter=",", skiprows=1)  # K_T 数据
data_kq = np.loadtxt("curve_2_actual_data.csv", delimiter=",", skiprows=1)  # K_Q 数据

# 分离 J 和 K_T / K_Q 数据
J_kt = data_kt[:, 0]
K_T = data_kt[:, 1]
J_kq = data_kq[:, 0]
K_Q = data_kq[:, 1]

# 数据过滤（确保 J 和 K_T 在合理范围内）
valid_kt = (J_kt >= 0.0) & (J_kt <= 1.0) & (K_T >= 0.0) & (K_T <= 0.3)
valid_kq = (J_kq >= 0.0) & (J_kq <= 1.0) & (K_Q >= 0.0) & (K_Q <= 0.3)

J_kt = J_kt[valid_kt]
K_T = K_T[valid_kt]
J_kq = J_kq[valid_kq]
K_Q = K_Q[valid_kq]

# 拟合多项式
degree = 2
coeff_kt = np.polyfit(J_kt, K_T, degree)
coeff_kq = np.polyfit(J_kq, K_Q, degree)

# 创建多项式函数
fit_kt = np.poly1d(coeff_kt)
fit_kq = np.poly1d(coeff_kq)

# 可视化拟合结果
J_fit = np.linspace(0, 1, 100)  # 在 J 范围内生成更多点

plt.figure(figsize=(12, 6))

# 绘制 K_T 曲线
plt.subplot(1, 2, 1)
plt.scatter(J_kt, K_T, label="Extracted Data (K_T)", color="blue")
plt.plot(J_fit, fit_kt(J_fit), label="Fitted Curve (K_T)", color="red")
plt.xlabel("J (Advance Ratio)")
plt.ylabel("K_T (Thrust Coefficient)")
plt.legend()
plt.title("K_T vs J")
plt.grid()

# 绘制 K_Q 曲线
plt.subplot(1, 2, 2)
plt.scatter(J_kq, K_Q, label="Extracted Data (K_Q)", color="blue")
plt.plot(J_fit, fit_kq(J_fit), label="Fitted Curve (K_Q)", color="red")
plt.xlabel("J (Advance Ratio)")
plt.ylabel("K_Q (Torque Coefficient)")
plt.legend()
plt.title("K_Q vs J")
plt.grid()

plt.tight_layout()
plt.show()

# 打印拟合多项式
print("Fitted Polynomial for K_T:")
print(fit_kt)
print("\nFitted Polynomial for K_Q:")
print(fit_kq)

# 保存拟合系数
np.savetxt("kt_coefficients_corrected.csv", coeff_kt, delimiter=",", header="a0,a1,a2", comments="")
np.savetxt("kq_coefficients_corrected.csv", coeff_kq, delimiter=",", header="b0,b1,b2", comments="")
print("拟合系数已保存到 'kt_coefficients_corrected.csv' 和 'kq_coefficients_corrected.csv'")
