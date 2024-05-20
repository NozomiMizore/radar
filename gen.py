import csv
import random

def generate_random_radar_data(num_radars, lon_range, lat_range, elev_range, r_range, output_file):
    """
    随机生成雷达数据并写入CSV文件

    :param num_radars: 生成的雷达数量
    :param lon_range: 经度的范围 (min_lon, max_lon)
    :param lat_range: 纬度的范围 (min_lat, max_lat)
    :param elev_range: 高程的范围 (min_elev, max_elev)
    :param r_range: 探测半径的范围 (min_r, max_r)
    :param output_file: 输出的 CSV 文件名
    """
    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        # 写入表头
        writer.writerow(['x', 'y', 'z', 'r'])

        for _ in range(num_radars):
            x = random.uniform(*lon_range)
            y = random.uniform(*lat_range)
            z = random.uniform(elev_range[0] - 100, elev_range[1] + 100)  # 比高程范围大一点
            r = random.randint(*r_range)
            writer.writerow([x, y, z, r])

# 参数配置
num_radars = 256  # 生成256个雷达数据
lon_range = (109.99958381761098, 114.99958381761098)  # 经度范围
lat_range = (25.000417247801863, 30.000417247801863)  # 纬度范围
elev_range = (-26, 2124)  # 高程范围
r_range = (1, 50)  # 探测半径范围，1千米到50千米
output_file = './data/radar_data.csv'  # 输出文件名

# 生成并写入数据
generate_random_radar_data(num_radars, lon_range, lat_range, elev_range, r_range, output_file)
