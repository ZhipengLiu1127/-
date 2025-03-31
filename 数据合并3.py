import pandas as pd
import xarray as xr


# 读取船舶航行的CSV文件并转换时间
def load_ship_data(file_path):
    ship_data = pd.read_csv(file_path)
    # 将year, month, day, hour, minute列合并成一个datetime对象
    ship_data['timestamp'] = pd.to_datetime(
        ship_data[['Year', 'Month', 'Day', 'Hour', 'Minute']])
    
    # 将日本时间（JST）转换为标准时（UTC），JST比UTC快9小时
    ship_data['timestamp'] = ship_data['timestamp'] - pd.Timedelta(hours=9)
    
    return ship_data


# 从风数据的nc文件中提取给定时间和经纬度的风速数据
def extract_wind_data(file_path, time, lat, lon):
    """
    从风数据的NetCDF文件中提取给定时间和经纬度的风速数据。
    """
    # 加载风数据
    wind_data = xr.open_dataset(file_path)

    # 确保时间维度正确，如果是 'valid_time' 而不是 'time'
    if 'valid_time' in wind_data.coords:
        time_coord = 'valid_time'
    else:
        time_coord = 'time'

    # 使用时间、纬度和经纬度来获取最近的风速数据
    closest_time = wind_data.sel({time_coord: time}, method='nearest')  # 最近时间匹配

    # 使用正确的经纬度维度名：'latitude' 和 'longitude'
    closest_wind = closest_time.sel(latitude=lat, longitude=lon, method='nearest')  # 最近经纬度匹配

    # 提取风速数据（例如u10和v10）
    wind_u = closest_wind['u10'].values.item()  # 东西风速
    wind_v = closest_wind['v10'].values.item()  # 南北风速

    return wind_u, wind_v


# 从浪数据的nc文件中提取给定时间和经纬度的波高、波浪主方向和波浪主周期数据
def extract_wave_data(file_path, time, lat, lon):
    """
    从浪数据的NetCDF文件中提取给定时间和经纬度的波高、波浪主方向和波浪主周期数据。
    """
    # 加载浪数据
    wave_data = xr.open_dataset(file_path)

    # 确保时间维度正确，如果是 'valid_time' 而不是 'time'
    if 'valid_time' in wave_data.coords:
        time_coord = 'valid_time'
    else:
        time_coord = 'time'

    # 使用时间、纬度和经纬度来获取最近的波高、主方向和主周期数据
    closest_time = wave_data.sel({time_coord: time}, method='nearest')  # 最近时间匹配

    # 使用正确的经纬度维度名：'latitude' 和 'longitude'
    closest_wave = closest_time.sel(latitude=lat, longitude=lon, method='nearest')  # 最近经纬度匹配

    # 提取波高、主方向和主周期数据
    wave_height = closest_wave['swh'].values.item()  # 波高
    wave_direction = closest_wave['mwd'].values.item()  # 波浪主方向
    wave_period = closest_wave['mwp'].values.item()  # 波浪主周期

    return wave_height, wave_direction, wave_period

# 匹配船舶数据和风浪数据，返回包含匹配结果的新DataFrame
def match_ship_and_wind_wave(ship_data, wind_file_path, wave_file_path):
    """
    匹配船舶数据和风浪数据，返回包含匹配结果的新DataFrame。
    """
    matched_data = []

    # 遍历每条船舶记录
    for index, row in ship_data.iterrows():
        time = row['timestamp']
        lat = row['Lat']
        lon = row['Long']

        # 提取该时间和经纬度对应的风速数据
        wind_u, wind_v = extract_wind_data(wind_file_path, time, lat, lon)

        # 提取该时间和经纬度对应的波高数据
        wave_height,wave_direction,wave_period = extract_wave_data(wave_file_path, time, lat, lon)

        # 提取船舶的必要信息（如heading, speed, rpm）
        ship_info = {
            'timestamp': time,
            'heading': row['GPS_bearing'],  # 假设'Heading'列存在
            'speed': row['OG_speed'],  # 假设'Speed'列存在
            'rpm': row['Eng_rpm'],  # 假设'RPM'列存在
            'wind_u': wind_u,  # 风速东西分量
            'wind_v': wind_v,  # 风速南北分量
            'wave_height': wave_height,  # 波高
            'wave_direction': wave_direction,  # 波浪主方向
            'wave_period': wave_period  # 波浪主周期
        }

        # 将匹配的数据添加到结果列表中
        matched_data.append(ship_info)

    # 将结果转为DataFrame
    matched_df = pd.DataFrame(matched_data)
    return matched_df


# 主函数：读取船舶数据、匹配风浪数据并保存为新的CSV
def main(ship_file_path, wind_file_path, wave_file_path, output_file_path):
    # 1. 读取船舶航行数据
    ship_data = load_ship_data(ship_file_path)

    # 2. 匹配船舶数据和风浪数据
    matched_df = match_ship_and_wind_wave(ship_data, wind_file_path, wave_file_path)

    # 3. 将匹配的结果保存为新的CSV文件
    matched_df.to_csv(output_file_path, index=False)

    # 打印匹配的结果（前5行）
    print(matched_df.head())


# 示例文件路径
ship_file_path = r"C:\Users\15585\Desktop\shipdata\case8-6.csv"  # 替换为船舶航行数据的CSV文件路径
wind_file_path = r"D:\data\wind_2013.nc"  # 替换为风速数据的NetCDF文件路径
wave_file_path = r"D:\data\wave_2013.nc"  # 替换为波高数据的NetCDF文件路径
output_file_path = 'matched_ship_wind_wave_data_2.csv'  # 输出文件路径

# 执行主函数
main(ship_file_path, wind_file_path, wave_file_path, output_file_path)
