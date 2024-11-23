/** BELLA experiment track reconstruction framework
 *
 * (c) 2024 Lawrence Berkeley National Laboratory
 *
 * Mozilla Public License Version 2.0
 */

// System include(s).
#include <iostream>
#include <fstream>

// The main routine
//

bool is_in_magnet(const double x, const double y, const double z)
{
    bool is_x_in_magnet = (x >= 40.f && x <= 50.f) || (x >= 210.f && x <= 220.f);
    bool is_y_in_magnet = (std::abs(y) <= 10.f);
    bool is_z_in_magnet = (std::abs(z) <= 10.f);

    return is_x_in_magnet && is_y_in_magnet && is_z_in_magnet;
}

int main(int argc, char *argv[])
{

    // Residual file
    std::ofstream bfield_file;
    bfield_file.open("bfield.txt");

    // Cell size = 10 mm
    double spacing = 10.f;

    const double start_x = -100.f;
    const double end_x = 1000.f;
    const double start_z = -500.f;
    const double end_z = 500.f;
    const double start_y = -500.f;
    const double end_y = 500.f;

    double x{start_x};
    double y{start_y};
    double z{start_z};
    while (x < end_x)
    {
        while (y < end_y)
        {
            while (z < end_z)
            {
                // @NOTE: What is the b direction?
                const double bx = 0.f;
                const double bz = 0.f;
                const double by = is_in_magnet(x, y, z) ? 0.5f : 0.f; // in Tesla
                bfield_file << x << " " << y << " " << z << " "
                            << bx << " " << by << " " << bz << "\n";

                z += spacing;
            }
            y += spacing;
            z = start_z;
        }
        x += spacing;
        y = start_y;
        z = start_z;
    }

    bfield_file.close();
    return 1;
}