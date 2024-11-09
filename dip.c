#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816
#endif

#define max_buffer_size 256

const int stdout_precision = 10;
const double magnetic_permeability = 4 * M_PI * 1e-7;

typedef struct {
	double x;
	double y;
	double z;
} Vector;

double degree_to_radians(double degree) {
	return M_PI * degree / 180.0;
}

Vector incdec_to_unitvec(double inclination, double declination) {
	double inclination_rad = degree_to_radians(inclination);
	double declination_rad = degree_to_radians(declination);

	double x = cos(inclination_rad) * sin(declination_rad);
	double y = cos(inclination_rad) * cos(declination_rad);
	double z = sin(inclination_rad);

	return (Vector){.x = x, .y = y, .z = z};
}

void usage(char** argv) {
	printf("[USAGE]: %s field_inclination field_declination magnetization_inclination magnetization_declination x_dipole y_dipole z_dipole x_observation y_observation z_observation\n", argv[0]);
	printf("[USAGE]: *_inclination in degrees\n");
	printf("[USAGE]: *_declination in degrees\n");
	fflush(stdout);
}

double dot(Vector a, Vector b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

double dipole(Vector magnetization, Vector field, Vector observation, Vector source) {
	if (observation.x == source.x && observation.y == source.y && observation.z == source.z) {
		fprintf(stderr, "[ERROR]: Observation and source locations are the same\n");
		return 0.0;
	}

	Vector displacement = {
		.x = observation.x - source.x,
		.y = observation.y - source.y,
		.z = observation.z - source.z
	};
	double distance = sqrt(displacement.x * displacement.x + 
			displacement.y * displacement.y + 
			displacement.z * displacement.z);

	double m_dot_d = dot(magnetization, displacement);
	Vector bfield = {0};
	bfield.x = 3.0 * m_dot_d * displacement.x / distance / distance / distance / distance / distance -
		magnetization.x / distance / distance / distance;
	bfield.y = 3.0 * m_dot_d * displacement.y / distance / distance / distance / distance / distance -
		magnetization.y / distance / distance / distance;
	bfield.z = 3.0 * m_dot_d * displacement.z / distance / distance / distance / distance / distance -
		magnetization.z / distance / distance / distance;

	// bfield.x *= magnetic_permeability / 3.0;
	// bfield.y *= magnetic_permeability / 3.0;
	// bfield.z *= magnetic_permeability / 3.0;

	// printf("[INFO]: B-field x %lf nT\n", bfield.x * 1e9);
	// printf("[INFO]: B-field y %lf nT\n", bfield.y * 1e9);
	// printf("[INFO]: B-field z %lf nT\n", bfield.z * 1e9);
	double total_field = dot(field, bfield);
	total_field *= magnetic_permeability / 3.0;
	return total_field;
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
		Vector source = {0};
		int number_observations = 0;

		sscanf(argv[1], "%lf", &field_inclination);
		sscanf(argv[2], "%lf", &field_declination);
		sscanf(argv[3], "%lf", &magnetization_inclination);
		sscanf(argv[4], "%lf", &magnetization_declination);
		sscanf(argv[5], "%lf", &source.x);
		sscanf(argv[6], "%lf", &source.y);
		sscanf(argv[7], "%lf", &source.z);
		sscanf(argv[9], "%d", &number_observations);

		Vector field         = incdec_to_unitvec(field_inclination, field_declination);
		Vector magnetization = incdec_to_unitvec(magnetization_inclination, magnetization_declination);

		FILE* observation_file = fopen(argv[8], "r");
		if (observation_file == NULL) {
			fprintf(stderr, "[ERROR]: Could not open %s", argv[8]);
			return 1;
		}

		char buffer[max_buffer_size] = {0}; 
		for (int i = 0; i < number_observations; i++) {
			Vector observation = {0};
			char* status = fgets(buffer, max_buffer_size, observation_file);
			if (status == NULL) {
				fclose(observation_file);
				fprintf(stderr, "[ERROR]: Could not read file contents into buffer\n");
				return 1;
			}
			sscanf(buffer, "%lf %lf %lf", &observation.x, &observation.y, &observation.z);

			double total_field = dipole(magnetization, field, observation, source);
			fprintf(stdout, "%.*e\n", stdout_precision, total_field);
		}

		fclose(observation_file);

	} else if (argc == 11) {
		double field_inclination;
		double field_declination;
		double magnetization_inclination;
		double magnetization_declination;
		Vector source = {0};
		Vector observation = {0};

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

		Vector field         = incdec_to_unitvec(field_inclination, field_declination);
		Vector magnetization = incdec_to_unitvec(magnetization_inclination, magnetization_declination);

		double total_field = dipole(magnetization, field, observation, source);

		fprintf(stdout, "%.*e\n", stdout_precision, total_field);
	} else {
		usage(argv);
		fprintf(stderr, "[ERROR]: Got %d arguments\n", argc - 1);
		return 1;
	}

	// printf("[INFO]: Program \"%s\"\n", argv[0]);

	// printf("[INFO]: field inclination %lf degree\n", field_inclination);
	// printf("[INFO]: field declination %lf degree\n", field_declination);
	// printf("[INFO]: magnetization inclination %lf degree\n", magnetization_inclination);
	// printf("[INFO]: magnetization declination %lf degree\n", magnetization_declination);
	// printf("[INFO]: dipole easting  %lf meter\n", dipole.x);
	// printf("[INFO]: dipole northing %lf meter\n", dipole.y);
	// printf("[INFO]: dipole vertical %lf meter\n", dipole.z);
	// printf("[INFO]: observation easting %lf  meter\n", observation.x);
	// printf("[INFO]: observation northing %lf meter\n", observation.y);
	// printf("[INFO]: observation vertical %lf meter\n", observation.z);


	// printf("[INFO]: Total-field anomaly %lf nT\n", total_field * 1e9);
	return 0;
}
