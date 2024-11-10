#ifndef DIP_H
#define DIP_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592307816
#endif

const double dip_magnetic_permeability = 4 * M_PI * 1e-7;

typedef struct {
	double x;
	double y;
	double z;
} dip_Vector;

double dip_degtorad(double degree) {
	return M_PI * degree / 180.0;
}

dip_Vector dip_incdectounitvec(double inclination, double declination) {
	double inclination_rad = dip_degtorad(inclination);
	double declination_rad = dip_degtorad(declination);

	double x = cos(inclination_rad) * sin(declination_rad);
	double y = cos(inclination_rad) * cos(declination_rad);
	double z = sin(inclination_rad);

	return (dip_Vector){.x = x, .y = y, .z = z};
}

double dip_dot(dip_Vector a, dip_Vector b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

int dip_dipole(dip_Vector magnetization, dip_Vector field, dip_Vector observation, dip_Vector source, double* total_field_anomaly) {
	if (observation.x == source.x && observation.y == source.y && observation.z == source.z) {
		fprintf(stderr, "[ERROR]: Observation and source locations are the same\n");
		return 1;
	}

	dip_Vector displacement = {
		.x = observation.x - source.x,
		.y = observation.y - source.y,
		.z = observation.z - source.z
	};
	double distance = sqrt(displacement.x * displacement.x + 
			displacement.y * displacement.y + 
			displacement.z * displacement.z);

	double m_dot_d = dip_dot(magnetization, displacement);
	dip_Vector bfield = {0};
	bfield.x = 3.0 * m_dot_d * displacement.x / distance / distance / distance / distance / distance -
		magnetization.x / distance / distance / distance;
	bfield.y = 3.0 * m_dot_d * displacement.y / distance / distance / distance / distance / distance -
		magnetization.y / distance / distance / distance;
	bfield.z = 3.0 * m_dot_d * displacement.z / distance / distance / distance / distance / distance -
		magnetization.z / distance / distance / distance;

	*total_field_anomaly = dip_magnetic_permeability / 3.0 * dip_dot(field, bfield);

	return 0;
}

#endif // DIP_H
