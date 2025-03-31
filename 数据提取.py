import matplotlib.pyplot as plt
import numpy as np
from matplotlib.image import imread

# 加载图片文件
image_path = "静水功率曲线.png"  # 替换为您的图片路径
img = imread(image_path)

# 显示图片并让用户标定坐标轴范围
fig, ax = plt.subplots(figsize=(8, 6))
ax.imshow(img)
ax.set_title("请点击坐标轴上的四个标定点（J_min, J_max, K_min, K_max）")
plt.xlabel("J (Advance Ratio)")  # 根据图片轴标签修改
plt.ylabel("K_T or K_Q")  # 根据图片轴标签修改

# 提示用户点击四个标定点
print("请在图中依次点击坐标轴上的四个标定点：")
print("1. J_min 对应的像素点 (左下角)")
print("2. J_max 对应的像素点 (右下角)")
print("3. K_min 对应的像素点 (左下角)")
print("4. K_max 对应的像素点 (左上角)")

# 获取用户点击的标定点
calibration_points = plt.ginput(4, timeout=0, show_clicks=True)
calibration_points = np.array(calibration_points)
plt.close()

# 提取标定点像素位置
J_min_pixel, KT_min_pixel = calibration_points[0]  # J_min 对应的像素位置
J_max_pixel, _ = calibration_points[1]             # J_max 对应的像素位置
_, KT_max_pixel = calibration_points[3]            # K_max 对应的像素位置

# 自动推断物理范围
# 提示用户输入坐标轴上的实际值（程序仅需输入一次）
print("请输入 J_min 和 J_max 的物理值（例如 0.0 和 1.0）：")
J_min, J_max = map(float, input().split())
print("请输入 K_T_min 和 K_T_max 的物理值（例如 0.0 和 0.3）：")
KT_min, KT_max = map(float, input().split())

# 转换函数：像素坐标到实际物理坐标
def pixel_to_actual(pixel_data, x_min_pix, x_max_pix, y_min_pix, y_max_pix, x_min, x_max, y_min, y_max):
    x_pixel, y_pixel = pixel_data[:, 0], pixel_data[:, 1]
    
    # 转换 J（横坐标）
    x_actual = (x_pixel - x_min_pix) / (x_max_pix - x_min_pix) * (x_max - x_min)
    
    # 转换 K_T（纵坐标），并修正 y 轴方向
    y_actual = (y_pixel - y_min_pix) / (y_max_pix - y_min_pix) * (y_max - y_min)
    #y_actual = y_max - y_actual  # 反转方向

    return np.column_stack((x_actual, y_actual))

# 显示图片并开始选取曲线点
fig, ax = plt.subplots(figsize=(8, 6))
ax.imshow(img)
ax.set_title("Click to extract data points for each curve. Press Enter to finish each curve.")
plt.xlabel("J (Advance Ratio)")  # 根据图片轴标签修改
plt.ylabel("K_T or K_Q")  # 根据图片轴标签修改

# 存储多条曲线的数据
all_curves = []

while True:
    print("请在图中点击选取点，完成一条曲线后按 Enter 键。如果提取完成，直接关闭窗口。")

    # 手动选取点
    points = plt.ginput(n=-1, timeout=0, show_clicks=True)  # n=-1 表示不限制点数
    if not points:  # 如果没有选择点，则跳出循环
        print("未选取数据点，跳过当前曲线。")
        break

    # 将提取的像素点数据存储
    points = np.array(points)
    
    # 将像素坐标转换为实际坐标
    actual_points = pixel_to_actual(points, J_min_pixel, J_max_pixel, KT_min_pixel, KT_max_pixel, J_min, J_max, KT_min, KT_max)
    all_curves.append(actual_points)
    print(f"已提取并转换 {len(points)} 个数据点，继续提取下一条曲线或关闭窗口结束。")

    # 是否继续
    continue_input = input("是否继续提取下一条曲线？(y/n): ").strip().lower()
    if continue_input == 'n':
        print("数据提取已完成，程序即将退出。")
        break

# 关闭图片显示
plt.close()

# 处理多条曲线数据
for i, curve in enumerate(all_curves):
    x_data = curve[:, 0]
    y_data = curve[:, 1]
    print(f"曲线 {i + 1} 的实际 X 数据 (J):", x_data)
    print(f"曲线 {i + 1} 的实际 Y 数据 (K_T 或 K_Q):", y_data)

    # 保存每条曲线的数据到单独的文件
    np.savetxt(f"curve_{i + 1}_actual_data.csv", curve, delimiter=",", header="J,K", comments="")
    print(f"曲线 {i + 1} 的实际数据已保存到 curve_{i + 1}_actual_data.csv 文件中。")

print("所有曲线数据提取完成！")
