#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define R 6371e3  // 地球半径，单位：米

typedef struct {
    double x;  // 经度
    double y;  // 纬度
    double z;  // 高程（米）
    int r;     // 探测半径（千米）
} Radar;

void read_csv(const char* filename, double** data, int* size) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Cannot open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // 获取行数
    char line[256];
    int count = 0;
    while (fgets(line, sizeof(line), file)) {
        count++;
    }
    rewind(file);

    *size = count;
    *data = (double*)malloc(count * sizeof(double));
    if (!(*data)) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    int i = 0;
    while (fgets(line, sizeof(line), file)) {
        (*data)[i++] = atof(line);
    }

    fclose(file);
}

Radar* read_radar_data(const char* filename, int* num_radars) {
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
    Radar* radars = (Radar*)malloc(*num_radars * sizeof(Radar));
    if (!radars) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    int i = 0;
    while (fgets(line, sizeof(line), file)) {
        sscanf(line, "%lf,%lf,%lf,%d", &radars[i].x, &radars[i].y, &radars[i].z, &radars[i].r);
        i++;
    }

    fclose(file);
    return radars;
}

double haversine(double lon1, double lat1, double lon2, double lat2) {
    double phi1 = lat1 * M_PI / 180.0;
    double phi2 = lat2 * M_PI / 180.0;
    double delta_phi = (lat2 - lat1) * M_PI / 180.0;
    double delta_lambda = (lon2 - lon1) * M_PI / 180.0;

    double a = sin(delta_phi / 2.0) * sin(delta_phi / 2.0) +
               cos(phi1) * cos(phi2) * sin(delta_lambda / 2.0) * sin(delta_lambda / 2.0);
    double c = 2.0 * atan2(sqrt(a), sqrt(1 - a));

    double distance = R * c;
    return distance;
}

void calculate_radar_coverage(Radar* radars, int num_radars, double* lons, double* lats, double* elevs, int num_points) {
    FILE* output_file = fopen("./res/result.txt", "w+");
    if (!output_file) {
        fprintf(stderr, "Cannot open file radar_coverage_results.csv for writing\n");
        exit(EXIT_FAILURE);
    }
    #pragma omp parallel for num_threads(4)
    for (int i = 0; i < num_radars; i++) {
        for (int j = 0; j < num_points; j++) {
            double horizontal_distance = haversine(radars[i].x, radars[i].y, lons[j], lats[j]);
            double vertical_distance = fabs(radars[i].z - elevs[j]);
            double distance = sqrt(horizontal_distance * horizontal_distance + vertical_distance * vertical_distance);

            if (distance <= radars[i].r * 1000) {  // 距离在探测半径内
                #pragma omp critical
                {
                    fprintf(output_file, "Radar %d can detect point at (%lf, %lf, %lf)\n", i, lons[j], lats[j], elevs[j]);
                }
            }
        }
    }
    fclose(output_file);
}

int main() {
    double *lons, *lats, *elevs;
    int num_points;

    read_csv("./data/lon.csv", &lons, &num_points);
    read_csv("./data/lat.csv", &lats, &num_points);
    read_csv("./data/elevation.csv", &elevs, &num_points);

    int num_radars;
    Radar* radars = read_radar_data("./data/radar_data.csv", &num_radars);

    calculate_radar_coverage(radars, num_radars, lons, lats, elevs, num_points);

    free(lons);
    free(lats);
    free(elevs);
    free(radars);

    return 0;
}

