// Copyright (C) 2026 Martí Raya Moreno
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
// implied. See the License for the specific language governing
// permissions and limitations under the License.
//
// @file
// Reads the elphbolt outputs and generates using a generalized 
// Vermeersch model (see https://doi.org/10.1063/1.4948968) 
// the cross-plane thermal conductivity. 
// Note that the model was originally created fitting to RTA MC results
// therefore it is a first approach to beyond-RTA. In other words
// For beyond-RTA such a solution does not exist, however, as first approximation 
// one can modify the formula replacing the mean free paths by its generalized counterpart. 
// This assumes that beyond-RTA contributions are suppressed by the boundaries 
// similarly to the RTA-predicted contributions.

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <iomanip>

using namespace std;

constexpr int naxis = 3;
constexpr double kB = 8.617333262e-5; // eV/K

// ------------------------------------------------------------
// Generate T directory names compatible with elphbolt
// ------------------------------------------------------------
inline std::string get_T_folder_name(double T) {

    int exp = static_cast<int>(std::floor(std::log10(std::abs(T)))) + 1;
    double mant = T / std::pow(10.0, exp);

    std::ostringstream oss;
    oss << std::uppercase << std::fixed << std::setprecision(3) << mant;

    // Build exponent part
    oss << "E"
        << (exp >= 0 ? "+" : "-")
        << std::setw(2) << std::setfill('0') << std::abs(exp);

    return "T"+oss.str();

}

// ------------------------------------------------------------
// Read phonon velocities
// ------------------------------------------------------------
void read_vels(const string& fname,
               vector<vector<vector<double>>>& vels,
               int& nbands)
{
    ifstream fin(fname);
    string line;

    while (getline(fin, line)) {
        stringstream ss(line);
        vector<double> data;
        double x;

        while (ss >> x) data.push_back(x);

        nbands = data.size() / naxis;

        vector<vector<double>> v_at_q(nbands, vector<double>(naxis));
        for (int ib = 0; ib < nbands; ++ib)
            for (int ax = 0; ax < naxis; ++ax)
                v_at_q[ib][ax] = data[nbands * ax + ib];

        vels.push_back(v_at_q);
    }
}

// ------------------------------------------------------------
// Read energies
// ------------------------------------------------------------
void read_energies(const string& fname,
                   vector<vector<double>>& ens)
{
    ifstream fin(fname);
    string line;

    while (getline(fin, line)) {
        stringstream ss(line);
        vector<double> row;
        double x;
        while (ss >> x) row.push_back(x);
        ens.push_back(row);
    }
}

// ------------------------------------------------------------
// Read response function F
// ------------------------------------------------------------
void read_response(const int nbands, const double T,
                   vector<vector<vector<double>>>& F)
{
    
    const auto T_folder = get_T_folder_name(T);

    for (int ib = 0; ib < nbands; ++ib) {
        string fname = T_folder + "/nodrag_F0_" + to_string(ib + 1);
        ifstream fin(fname);

        string line;
        int iq = 0;

        while (getline(fin, line)) {
            if (ib == 0)
                F.push_back(vector<vector<double>>());

            stringstream ss(line);
            vector<double> v(3);
            ss >> v[0] >> v[1] >> v[2];

            if (ib == 0)
                F[iq].resize(nbands);

            F[iq][ib] = v;
            iq++;
        }
    }
}

// ------------------------------------------------------------
// Compute Bose Einstein distribution
// ------------------------------------------------------------
inline double Bose(const double e, const double T)
{
    return 1.0 / std::expm1(e / (kB * T));
}

vector<vector<double>> compute_bose(const vector<vector<double>>& ens,
                                    const double T)
{
    int nq = ens.size();
    int nb = ens[0].size();

    vector<vector<double>> bose(nq, vector<double>(nb));

    for (int iq = 0; iq < nq; ++iq) {
        for (int ib = 0; ib < nb; ++ib) {
            double e = ens[iq][ib];

            // Avoid acoustic modes at Gamma
            if (fabs(e) < 1e-12)
                bose[iq][ib] = 0.0;
            else
                bose[iq][ib] = Bose(e, T);
        }
    }
    return bose;
}


// ------------------------------------------------------------
// Compute GMFP out-of-plane
// ------------------------------------------------------------
vector<vector<double>> compute_gmfp(const vector<vector<vector<double>>>& F,
                                    const vector<double>& normal,
                                    double T,
                                    const vector<vector<double>>& ens)
{
    int nq = ens.size();
    int nb = ens[0].size();

    vector<vector<double>> gmfp(nq, vector<double>(nb, 0.0));

    for (int iq = 0; iq < nq; ++iq) {
        for (int ib = 0; ib < nb; ++ib) {
            if (ib < 3 && iq == 0) {
                gmfp[iq][ib] = 0.0;
            } else {
                double dot = 0.0;
                for (int ax = 0; ax < 3; ++ax)
                    dot += F[iq][ib][ax] * normal[ax];

                gmfp[iq][ib] = T * dot / ens[iq][ib];
            }
        }
    }
    return gmfp;
}

// ------------------------------------------------------------
// Get from the std out from elphbolt the volume
// ------------------------------------------------------------
double getPrimitiveCellVolume(const string& filename) {
    ifstream file(filename);
    string line;
    const string target = "Primitive cell volume =";

    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }

    while (getline(file, line)) {
        size_t foundPos = line.find(target);
        if (foundPos != string::npos) {
            size_t equalSignPos = line.find('=', foundPos);
            string dataPart = line.substr(equalSignPos + 1);

            stringstream ss(dataPart);
            double volume;
            if (ss >> volume) {
                return volume;
            }
        }
    }

    throw runtime_error("The volume was not found in the file");
    return -1.0;
}

// ------------------------------------------------------------
// Get from the std out from elphbolt the temperature
// ------------------------------------------------------------
double getT(const string& filename) {
    ifstream file(filename);
    string line;
    const string target = "Crystal temperature = ";

    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }

    while (getline(file, line)) {
        size_t foundPos = line.find(target);
        if (foundPos != string::npos) {
            size_t equalSignPos = line.find('=', foundPos);
            string dataPart = line.substr(equalSignPos + 1);

            stringstream ss(dataPart);
            double temperature;
            if (ss >> temperature) {
                return temperature;
            }
        }
    }

    throw runtime_error("The temperature was not found in the file");
    return -1.0;
}

// ------------------------------------------------------------
// Generate a thickness using a geometric progression
// ------------------------------------------------------------
vector<double> generate_thickness(double t_init, double t_final, int n_step) {
    vector<double> thickness(n_step);
    
    iota(thickness.begin(), thickness.end(), 0.0);
    double ratio = pow(t_final / t_init, 1.0 / (n_step - 1));
    transform(thickness.begin(), thickness.end(), thickness.begin(),
        [t_init, ratio](double i) {
            return t_init * pow(ratio, i);
        });

    return thickness;
}

// ------------------------------------------------------------
// Main
// ------------------------------------------------------------
int main(int argc, char* argv[]){

    if (argc != 8) {
        cout << "Description: " << endl;
	cout << "This program reads the elphbolt output and generates using the Vermeersch model\nthe cross-plane thermal conductivity" << endl; 
        cout << "Usage: " << endl;
	cout << "cross_plane ELPHBOLT_STDOUT n[0] n[1] n[2] thickness_begin[nm] thickness_end[nm] nsteps" << endl;
        return EXIT_SUCCESS;
    }

    auto T = getT(argv[1]);
    auto volume = getPrimitiveCellVolume(argv[1]);
    vector<double> normal = {atof(argv[2]), atof(argv[3]), atof(argv[4])};
    auto thickness = generate_thickness(atof(argv[5]),atof(argv[6]),atoi(argv[7]));

    const auto normalization_constant = sqrt(inner_product(normal.begin(), normal.end(), normal.begin(), 0.0));
    transform(normal.begin(), normal.end(), normal.begin(),
            [normalization_constant](double d) { return d / normalization_constant; });


    vector<vector<vector<double>>> vels;
    vector<vector<double>> ens;
    vector<vector<vector<double>>> F;

    int nbands = 0;

    read_vels("ph.vels_fbz", vels, nbands);
    read_energies("ph.ens_fbz", ens);
    read_response(nbands, T, F);

    int nq = ens.size();

    auto occ   = compute_bose(ens, T);
    auto gmfp = compute_gmfp(F, normal, T, ens);

    cout << "# INPUT => T [K] : " << T << endl;
    cout << "# INPUT => norm : " << normal[0] << '\t' << normal[1] << '\t' << normal[2] << endl;
    cout << "# INPUT => File : " << argv[1] << endl;
    cout << "# INPUT => Crystal volume [nm**3] : " << volume << endl;
    cout << "# INPUT => Smaller thickness [nm] : " << atof(argv[5]) << endl;
    cout << "# INPUT => Larger thickness [nm] : " << atof(argv[6]) << endl;
    cout << "# INPUT => N steps : " << atoi(argv[7]) << endl;
    cout << "################################# " << endl;
    cout << "# thickness [nm]  kappa_cross-plane [W/(K·m)]" << endl; 
    for (double t : thickness) {

        double kappa = 0.0;

        for (int iq = 0; iq < nq; ++iq) {
            for (int ib = 0; ib < nbands; ++ib) {
                double Kn = fabs(gmfp[iq][ib]) / t;
                double S  = 1.0 / (1.0 + 2.0 * Kn);
		auto df = occ[iq][ib] * ( 1.0 + occ[iq][ib]);

		auto F_proj = 0.0;
		for (int ax = 0; ax < 3; ++ax) F_proj = F[iq][ib][ax] * normal[ax];

                double vproj =
                    (vels[iq][ib][0] * normal[0] +
                     vels[iq][ib][1] * normal[1] +
                     vels[iq][ib][2] * normal[2]);
                kappa += S * ens[iq][ib] * df * vproj * F_proj;
            }
        }

        kappa = kappa * 1.602176634e-19 * 1.0e21 / kB / T/ nq / volume;
        cout << t << '\t' << kappa << std::endl;
    }

    return 0;
}


