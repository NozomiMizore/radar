#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#define R 6371e3  // 地球半径，单位：米
#define DEG2RAD(deg) ((deg) * M_PI / 180.0)
#define LINE_SIZE 65535*4
typedef struct {
    double x;  // 经度
    double y;  // 纬度
    double z;  // 高程（米）
    int r;     // 探测半径（千米）
} Radar;

void read_csv_data(const char* filename, double* data, int size) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Cannot open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
    char line[256];
    int i = 0;
    while (i < size) {
        fgets(line, sizeof(line), file);
        data[i] = atof(line);
        i++;
    }

    fclose(file);
}

void read_ele_data(const char* filename, double*** data, int width, int height) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Cannot open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    *data = (double**)malloc(height * sizeof(double*));
    if (!(*data)) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < width; i++) {
        (*data)[i] = (double*)malloc(width * sizeof(double));
        if (!(*data)[i]) {
            fprintf(stderr, "Memory allocation failed\n");
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }

    char line[LINE_SIZE];
    int i = 0;
    int j = 0;
    while (fgets(line, sizeof(line), file) && i < height) {
        char* token = strtok(line, ",");
        while(token != NULL && j < width){
            (*data)[i][j] = atof(token);
            token = strtok(NULL, ",");
            j++;
        }
        i++;
        j = 0;
    }
    fclose(file);
}

void read_radar_data(const char* filename, Radar** radars, int* num_radars) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Cannot open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    char line[256];
    int count = 0;
    while (fgets(line, sizeof(line), file)) {
        count++;
    }
    rewind(file);

    // Skip header
    fgets(line, sizeof(line), file);

    *num_radars = count - 1;
    *radars = (Radar*)malloc(*num_radars * sizeof(Radar));
    if (!(*radars)) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    int i = 0;
    while (fgets(line, sizeof(line), file)) {
        sscanf(line, "%lf,%lf,%lf,%d", &(*radars)[i].x, &(*radars)[i].y, &(*radars)[i].z, &(*radars)[i].r);
        i++;
    }

    fclose(file);
}

double haversine(double lon1, double lat1, double lon2, double lat2) {
    double phi1 = DEG2RAD(lat1);
    double phi2 = DEG2RAD(lat2);
    double delta_phi = DEG2RAD(lat2 - lat1);
    double delta_lambda = DEG2RAD(lon2 - lon1);

    double a = sin(delta_phi / 2.0) * sin(delta_phi / 2.0) +
               cos(phi1) * cos(phi2) * sin(delta_lambda / 2.0) * sin(delta_lambda / 2.0);
    double c = 2.0 * atan2(sqrt(a), sqrt(1 - a));

    double distance = R * c;
    return distance;
}

int is_blocked(double rx, double ry, double rz, int t_x_index, int t_y_index, double* lons, double* lats, double** elevs, double interval, int width, int height, int num_samples) {
    double t_z = elevs[t_x_index][t_y_index];
    double t_x = lons[t_x_index];
    double t_y = lats[t_y_index];
    int flag_x = 0;
    int flag_y = 0;
    if(fabs(t_x - rx) < 2*interval && fabs(t_y - ry) < 2*interval)
        return 0;
    
    for(int i = 0; i  <= num_samples; i++){
        double t = (double) i / num_samples;
        double lon = rx + t * (t_x - rx);
        double lat = ry + t * (t_y - ry);
        double ele = rz + t * (t_z - rz);

        int lon_index = (int)((lon - lons[0]) / interval);
        int lat_index = (int)((lats[0] - lat) / interval);

        if(elevs[lon_index][lat_index] > ele){
            return 1;
        }
    }
    return 0;  // 没有阻挡
}

void calculate_radar_coverage(Radar* radars, int num_radars, double* lons, double* lats, double** elevs, int width, int height, double interval, int num_samples) {
    FILE* output_file = fopen("res/result_07.csv", "w");
    if (!output_file) {
        fprintf(stderr, "result.csv for writing\n");
        exit(EXIT_FAILURE);
    }

    #pragma omp parallel for num_threads(4)
    for (int i = 0; i < num_radars; i++) {
        double min_lon = radars[i].x - radars[i].r * 0.0091;
        double max_lon = radars[i].x + radars[i].r * 0.0091;
        double min_lat = radars[i].y - radars[i].r * 0.01;
        double max_lat = radars[i].y + radars[i].r * 0.01;
        int min_lon_index = 0;
        int max_lon_index = width-1;
        int min_lat_index = 0;
        int max_lat_index = height-1;
        for(; min_lon_index < width; min_lon_index++){
            if(lons[min_lon_index] >= min_lon)
                break;
        }
        for(;max_lon_index > min_lon_index; max_lon_index--){
            if(lons[max_lon_index] <= max_lon)
                break;
        }
        for(; min_lat_index < height; min_lat_index++){
            if(lats[min_lat_index] <= max_lat)
                break;
        }
        for(; max_lat_index > min_lat_index; max_lat_index--){
            if(lats[max_lat_index] >= min_lat)
                break;
        }
        for(int j = min_lon_index; j <= max_lon_index; j++){
            for(int k = min_lat_index; k <= max_lat_index; k++){
                double vertical_distance = fabs(radars[i].z - elevs[j][k]);
                if(vertical_distance > radars[i].r * 1000)
                    continue;
                double horizontal_distance = haversine(radars[i].x, radars[i].y, lons[j], lats[k]);
                double distance = sqrt(horizontal_distance * horizontal_distance + vertical_distance * vertical_distance);
                if (distance <= radars[i].r * 1000) {  // 距离在探测半径内
                    if (!is_blocked(radars[i].x, radars[i].y, radars[i].z, j, k, lons, lats, elevs, interval, width, height, num_samples)) {
                        #pragma omp critical
                        {
                            fprintf(output_file, "Radar %d can detect point at (%lf, %lf, %lf)\n", i, lons[j], lats[k], elevs[j][k]);
                        }
                    }
                }
            }
        }
    }

    fclose(output_file);
}

int main() {
    int width = 6001;
    int height = 6001;
    int num_samples = 100;
    double* lons = (double*)malloc(width * sizeof(double));
    double* lats = (double*)malloc(height * sizeof(double));
    double** elevs;

    read_csv_data("./data/lon_07.csv", lons, width);
    read_csv_data("./data/lat_07.csv", lats, height);
    read_ele_data("./data/elevation_07.csv", &elevs, width, height);
    double interval = fabs(lons[1] - lons[0]);
    
    int num_radars;
    Radar* radars;
    read_radar_data("./data/radar_data_07.csv", &radars, &num_radars);

    time_t start, end;
    start = time(NULL);
    calculate_radar_coverage(radars, num_radars, lons, lats, elevs, width, height, interval, num_samples);
    end = time(NULL);
    printf("time=%lf\n",difftime(end,start));
    
    free(lons);
    free(lats);
    for (int i = 0; i < width; i++) {
        free(elevs[i]);
    }
    free(elevs);
    free(radars);

    return 0;
}
