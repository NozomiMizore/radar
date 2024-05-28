import numpy as np
import csv
import random
import pandas as pd

def generate_random_radar_data(num_radars, lon_range, lat_range, elev_range, r_range, output_file):

    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        # 写入表头
        writer.writerow(['x', 'y', 'z', 'r'])

        for _ in range(num_radars):
            x = random.uniform(*lon_range)
            y = random.uniform(*lat_range)
            z = random.uniform(elev_range[0] , elev_range[1])
            r = random.randint(*r_range)
            writer.writerow([x, y, z, r])


if __name__ == "__main__":

    lon = pd.read_csv('./data/lon.csv', header=None).squeeze()
    lat = pd.read_csv('./data/lat.csv', header=None).squeeze()
    elevation = pd.read_csv('./data/elevation.csv', header=None).values
    # 计算经度、纬度和高程的最大值和最小值
    lon_min, lon_max = lon.min(), lon.max()
    lat_min, lat_max = lat.min(), lat.max()
    elevation_min, elevation_max = np.min(elevation), np.max(elevation)

    # 输出最大值和最小值
    print(f"Longitude - Min: {lon_min}, Max: {lon_max}")
    print(f"Latitude - Min: {lat_min}, Max: {lat_max}")
    print(f"Elevation - Min: {elevation_min}, Max: {elevation_max}")

    num_radars = 8
    lon_range = (lon_min, lon_max)
    lat_range = (lat_min, lat_max)
    elev_range = (0, elevation_max/8)
    r_range = (1, 10)  # 探测半径范围为 1 到 10 千米
    output_file = './data/radar_data.csv'
    
    generate_random_radar_data(num_radars, lon_range, lat_range, elev_range, r_range, output_file)