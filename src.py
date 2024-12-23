import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.interpolate import griddata
from calme_water import calculate_clame_resistance
from wind import wind_resistance
from wave import wave_resistance
from datetime import datetime

# 1. 创建WGS84网格
def create_wgs84_grid(lon_min, lon_max, lat_min, lat_max, resolution):
    """
    创建WGS84网格，指定经纬度范围和分辨率
    """
    lon_grid = np.arange(lon_min, lon_max, resolution)
    lat_grid = np.arange(lat_min, lat_max, resolution)
    lon_grid, lat_grid = np.meshgrid(lon_grid, lat_grid)
    return lon_grid, lat_grid

# 读取风浪数据
def read_wind_wave_data(wind_file, wave_file, lat, lon, time_range):
    """
    从NC文件中读取风和波浪数据
    wind_file: ERA5风速数据文件
    wave_file: 波浪数据文件
    lat: 船舶的纬度
    lon: 船舶的经度
    time_range: 时间范围，格式为 [start_time, end_time]，时间格式为字符串
    """
    # 打开wind数据和wave数据
    wind_data = Dataset(wind_file, 'r')
    wave_data = Dataset(wave_file, 'r')
    
    # 读取时间信息
    times = pd.to_datetime(wind_data.variables['time'][:], unit='hours', origin=pd.Timestamp('1970-01-01'))
    
    # 筛选指定时间范围的数据
    time_mask = (times >= pd.to_datetime(time_range[0])) & (times <= pd.to_datetime(time_range[1]))
    
    # 纬度、经度坐标和风速
    lats = wind_data.variables['latitude'][:]
    lons = wind_data.variables['longitude'][:]
    
    # 找到最近的经纬度点
    lat_idx = np.argmin(np.abs(lats - lat))
    lon_idx = np.argmin(np.abs(lons - lon))
    
    # 提取指定经纬度和时间范围的风速和波浪数据
    u_wind = wind_data.variables['u10'][time_mask, lat_idx, lon_idx]
    v_wind = wind_data.variables['v10'][time_mask, lat_idx, lon_idx]
    wave_height = wave_data.variables['wave_height'][time_mask, lat_idx, lon_idx]
    wave_period = wave_data.variables['wave_period'][time_mask, lat_idx, lon_idx]
    
    return u_wind, v_wind, wave_height, wave_period, times[time_mask]

# 3. 读取船舶航迹数据
def read_ship_route(route_file):
    """
    读取船舶的航迹数据
    """
    route_data = pd.read_csv(route_file)
    route_data['datetime'] = pd.to_datetime(route_data['datetime'])
    return route_data

# 4. 根据时间和经纬度提取风浪数据
def extract_wind_wave_for_time(lon, lat, wind_speed, wave_height, wave_period, times, query_time, lon_grid, lat_grid):
    """
    根据指定时间和经纬度提取风浪数据，并插值到网格位置
    """
    # 找到对应时间的风浪数据索引
    time_idx = np.abs(times - query_time).argmin()  # 找到与查询时间最接近的索引

    # 提取对应时间的风速、波高和波周期
    wind_speed_time = wind_speed[time_idx]
    wave_height_time = wave_height[time_idx]
    wave_period_time = wave_period[time_idx]
    
    # 插值到目标网格位置
    wind_speed_grid = griddata((lon, lat), wind_speed_time, (lon_grid, lat_grid), method='linear')
    wave_height_grid = griddata((lon, lat), wave_height_time, (lon_grid, lat_grid), method='linear')
    wave_period_grid = griddata((lon, lat), wave_period_time, (lon_grid, lat_grid), method='linear')
    
    return wind_speed_grid, wave_height_grid, wave_period_grid

# 5. 计算船舶失速率
def calculate_stall_rate(route_data, lon_grid, lat_grid, wind_speed_grid, wave_height_grid, wave_period_grid):
    stall_rates = []
    
    # 遍历航迹数据，计算每个位置的失速率
    for index, row in route_data.iterrows():
        latitude = row['latitude']
        longitude = row['longitude']
        heading = row['heading']
        speed = row['speed']
        timestamp = row['datetime']
        
        # 计算船舶在该位置对应的网格单元
        grid_x = np.abs(lon_grid[0, :] - longitude).argmin()
        grid_y = np.abs(lat_grid[:, 0] - latitude).argmin()
        
        # 获取该网格单元的风浪数据
        wind_speed = wind_speed_grid[grid_y, grid_x]
        wave_height = wave_height_grid[grid_y, grid_x]
        wave_period = wave_period_grid[grid_y, grid_x]
        
        # 静水阻力计算
        static_resistance = calculate_clame_resistance(speed, ...)  # 填充船舶参数
        
        # 风浪附加阻力
        wind_res = wind_resistance(u=wind_speed, v=0, heading=heading, ship_speed=speed, At=..., mode_input='Heavy Ballast, Cranes')
        wave_res = wave_resistance(Lpp=..., B=..., T=wave_period, Vs=speed, Vg=0, C_B=..., kyy=..., T_F=..., T_A=..., heading=heading, wave_direction=0, Hs=wave_height)
        
        # 计算总阻力和失速率
        total_resistance = static_resistance + wind_res + wave_res
        stall_rate = total_resistance / speed
        stall_rates.append(stall_rate)
    
    return stall_rates

# 6. 可视化失速率
def visualize_stall_rate(route_data, stall_rates):
    plt.figure(figsize=(10, 6))
    plt.plot(route_data['datetime'], stall_rates, label='Stall Rate')
    plt.xlabel('Time')
    plt.ylabel('Stall Rate')
    plt.title('Stall Rate along the Route')
    plt.legend()
    plt.grid(True)
    plt.show()

# 7. 主程序
def main(wind_file, wave_file, route_file, lon_min, lon_max, lat_min, lat_max, resolution=0.1):
    # 创建WGS84网格
    lon_grid, lat_grid = create_wgs84_grid(lon_min, lon_max, lat_min, lat_max, resolution)
    
    # 读取风浪数据
    times, lon, lat, wind_speed, wave_height, wave_period = read_wind_wave_data(wind_file, wave_file)
    
    # 读取船舶航迹数据
    route_data = read_ship_route(route_file)
    
    stall_rates = []
    
    # 逐条船舶航迹数据，提取对应的风浪数据，并计算失速率
    for timestamp in route_data['datetime']:
        # 根据时间和经纬度提取风浪数据
        wind_speed_grid, wave_height_grid, wave_period_grid = extract_wind_wave_for_time(
            lon, lat, wind_speed, wave_height, wave_period, times, timestamp, lon_grid, lat_grid)
        
        # 计算失速率
        stall_rate = calculate_stall_rate(route_data, lon_grid, lat_grid, wind_speed_grid, wave_height_grid, wave_period_grid)
        stall_rates.append(stall_rate)
    
    # 可视化失速率
    visualize_stall_rate(route_data, stall_rates)

# 执行主程序
if __name__ == "__main__":
    wind_file = 'wind_data.nc'   # 风速数据文件
    wave_file = 'wave_data.nc'   # 波浪数据文件
    route_file = 'ship_route.csv' # 船舶航迹数据文件
    lon_min, lon_max, lat_min, lat_max = 120, 130, 30, 40  # 经度和纬度范围
    main(wind_file, wave_file, route_file, lon_min, lon_max, lat_min, lat_max)