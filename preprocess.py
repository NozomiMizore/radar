from osgeo import gdal
import numpy as np


def read_tif(path_tif, tif_file):
    dataset = gdal.Open(path_tif + tif_file)
    im_width = dataset.RasterXSize
    im_height = dataset.RasterYSize
    im_bands = dataset.RasterCount

    if im_bands == 1:
        band = dataset.GetRasterBand(1)
        data = band.ReadAsArray()
    
    im_geotrans = dataset.GetGeoTransform()
    x_range = range(0, im_width)
    y_range = range(0, im_height)
    x, y = np.meshgrid(x_range, y_range)
    lon = im_geotrans[0] + x * im_geotrans[1] + y * im_geotrans[2]
    lat = im_geotrans[3] + x * im_geotrans[4] + y * im_geotrans[5]
    
    return lon, lat, data


if __name__ == "__main__":

    # 示例
    lon, lat, elevation = read_tif('./data/', '59_07.tif')

    # 保存数据以供C程序读取
    np.savetxt('./data/lon.csv', lon, delimiter=',')
    np.savetxt('./data/lat.csv', lat, delimiter=',')
    np.savetxt('./data/elevation.csv', elevation, delimiter=',')