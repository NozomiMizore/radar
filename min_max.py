from osgeo import gdal
import numpy as np
def preproceess_read_tif(path_tif,tif_file):
        dataset = gdal.Open(path_tif + tif_file)
        print(dataset.ReadAsArray().shape)
        im_width = dataset.RasterXSize #栅格矩阵的列数
        im_height = dataset.RasterYSize #栅格矩阵的行数
        im_bands = dataset.RasterCount #波段数
        
        if im_bands == 1:
            band = dataset.GetRasterBand(im_bands)
            data = band.ReadAsArray()  
    
        im_geotrans = dataset.GetGeoTransform()#获取仿射矩阵信息
        x_range = range(0, im_width)
        y_range = range(0, im_height)
        x, y = np.meshgrid(x_range, y_range)
        lon = im_geotrans[0] + x * im_geotrans[1] + y * im_geotrans[2]
        lat = im_geotrans[3] + x * im_geotrans[4] + y * im_geotrans[5]
        
        lon_t = lon[:].flatten()
        lat_t = lat[:].flatten()
        data_t = data[:].flatten()        
        return lon_t, lat_t ,data_t

if __name__ == "__main__":
    path_out = './data/'
    # tif_file = '59_08.tif'
    tif_file = '59_07.tif'
    lon, lat, elevation = preproceess_read_tif(path_out,tif_file)
    
    # 计算经度、纬度和高程的最大值和最小值
    lon_min, lon_max = np.min(lon), np.max(lon)
    lat_min, lat_max = np.min(lat), np.max(lat)
    elevation_min, elevation_max = np.min(elevation), np.max(elevation)
    # 输出最大值和最小值
    print(tif_file)
    print(f"Longitude - Min: {lon_min}, Max: {lon_max}")
    print(f"Latitude - Min: {lat_min}, Max: {lat_max}")
    print(f"Elevation - Min: {elevation_min}, Max: {elevation_max}")