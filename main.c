#include <stdio.h>
#include <stdlib.h>

#include "dip.h"

#define max_buffer_size 256

const int stdout_precision = 10;

void usage(char** argv) {
	printf("[USAGE]: %s field_inclination field_declination magnetization_inclination magnetization_declination x_dipole y_dipole z_dipole x_observation y_observation z_observation\n", argv[0]);
	printf("[USAGE]: *_inclination in degrees\n");
	printf("[USAGE]: *_declination in degrees\n");
	fflush(stdout);
}

int main(int argc, char** argv) {

	if (argc == 1 || argc == 2) {
		usage(argv);
		return 0;
	} else if (argc == 10) {
		double field_inclination;
		double field_declination;
		double magnetization_inclination;
		double magnetization_declination;
		dip_Vector source = {0};
		int number_observations = 0;

		sscanf(argv[1], "%lf", &field_inclination);
		sscanf(argv[2], "%lf", &field_declination);
		sscanf(argv[3], "%lf", &magnetization_inclination);
		sscanf(argv[4], "%lf", &magnetization_declination);
		sscanf(argv[5], "%lf", &source.x);
		sscanf(argv[6], "%lf", &source.y);
		sscanf(argv[7], "%lf", &source.z);
		sscanf(argv[9], "%d", &number_observations);

		dip_Vector field         = dip_incdectounitvec(field_inclination, field_declination);
		dip_Vector magnetization = dip_incdectounitvec(magnetization_inclination, magnetization_declination);

		FILE* observation_file = fopen(argv[8], "r");
		if (observation_file == NULL) {
			fprintf(stderr, "[ERROR]: Could not open %s", argv[8]);
			return 1;
		}

		char buffer[max_buffer_size] = {0}; 
		for (int i = 0; i < number_observations; i++) {
			dip_Vector observation = {0};
			char* status = fgets(buffer, max_buffer_size, observation_file);
			if (status == NULL) {
				fclose(observation_file);
				fprintf(stderr, "[ERROR]: Could not read file contents into buffer\n");
				return 1;
			}
			sscanf(buffer, "%lf %lf %lf", &observation.x, &observation.y, &observation.z);

			double total_field = 0;
			if (dip_dipole(magnetization, field, observation, source, &total_field)) {
				fclose(observation_file);
				fprintf(stderr, "[ERROR]: Failed to compute dipole\n");
			}
			fprintf(stdout, "%.*e\n", stdout_precision, total_field);
		}

		fclose(observation_file);

	} else if (argc == 11) {
		double field_inclination;
		double field_declination;
		double magnetization_inclination;
		double magnetization_declination;
		dip_Vector source = {0};
		dip_Vector observation = {0};

		sscanf(argv[1], "%lf", &field_inclination);
		sscanf(argv[2], "%lf", &field_declination);
		sscanf(argv[3], "%lf", &magnetization_inclination);
		sscanf(argv[4], "%lf", &magnetization_declination);
		sscanf(argv[5], "%lf", &source.x);
		sscanf(argv[6], "%lf", &source.y);
		sscanf(argv[7], "%lf", &source.z);
		sscanf(argv[8], "%lf", &observation.x);
		sscanf(argv[9], "%lf", &observation.y);
		sscanf(argv[10], "%lf", &observation.z);

		dip_Vector field         = dip_incdectounitvec(field_inclination, field_declination);
		dip_Vector magnetization = dip_incdectounitvec(magnetization_inclination, magnetization_declination);

		double total_field = 0;
		if (dip_dipole(magnetization, field, observation, source, &total_field)) {
			fprintf(stderr, "[ERROR]: Failed to compute dipole\n");
		}

		fprintf(stdout, "%.*e\n", stdout_precision, total_field);
	} else {
		usage(argv);
		fprintf(stderr, "[ERROR]: Got %d arguments\n", argc - 1);
		return 1;
	}

	return 0;
}
