

#include <atmosphere.h>

atmosphere_t get_atmospheric_conditions(double altitude) {
    double rho, temp, p;

    if (altitude < 11000) {
        // Troposphere
        temp = 15.04 - 0.00649 * altitude;
        p = 101.29 * pow((temp + 273.15) / 288.08, 5.256);
    } else if (altitude < 25000) {
        // Lower stratosphere
        temp = -56.46;
        p = 22.65 * exp(1.73 - 0.000157 * altitude);
    } else {
        // Upper stratosphere
        temp = -131.21 + 0.00299 * altitude;
        p = 2.488 / pow((temp + 273.15) / 216.6, 11.388);
    }
    rho = p / (0.2869 * (temp + 273.15));
    atmosphere_t atm = {
        .rho = rho,
        .temp = temp,
        .p = p * 1000.0
    };

    return atm;
}
