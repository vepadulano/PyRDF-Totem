#ifndef _parameters_global_h_
#define _parameters_global_h_

#include <string>
#include <vector>
#include <map>
#include <cmath>

double timestamp0 = 1444860000;

string storageDir;

vector<string> distilledNtuples;

vector<AlignmentSource> alignmentSources;
Analysis anal;
Environment env;

string unsmearing_file;
string unsmearing_object;

string luminosity_data_file;

//----------------------------------------------------------------------------------------------------

void Init_global()
{
    // environment settings
    env.InitNominal();

    // binning
    // TODO
    anal.t_min = 0.02; anal.t_max = 3.5;
    anal.t_min_full = 0.; anal.t_max_full = 4.0;

    // approximate (time independent) resolutions
    // TODO
    anal.si_th_y_1arm = 3.1E-6 / sqrt(2.);
    anal.si_th_y_1arm_unc = 0.E-6 / sqrt(2.);

    anal.si_th_y_2arm = anal.si_th_y_1arm / sqrt(2.);
    anal.si_th_y_2arm_unc = 0E-6;

    anal.si_th_x_1arm_L = 0E-6;
    anal.si_th_x_1arm_R = 0E-6;
    anal.si_th_x_1arm_unc = 0E-6;

    anal.si_th_x_2arm = 0E-6;
    anal.si_th_x_2arm_unc = 0E-6;

    // analysis settings
    anal.th_x_lcut = -1.;
    anal.th_x_hcut = +1.;
}

//----------------------------------------------------------------------------------------------------

void Init_global_45b_56t()
{
    anal.th_y_lcut_L = 30E-6; anal.th_y_lcut_R = 33.5E-6; anal.th_y_lcut = 34.5E-6;
    anal.th_y_hcut_L = 102E-6; anal.th_y_hcut_R = 102E-6; anal.th_y_hcut = 100E-6;
}

//----------------------------------------------------------------------------------------------------

void Init_global_45t_56b()
{
    anal.th_y_lcut_L = 27E-6; anal.th_y_lcut_R = 27.5E-6; anal.th_y_lcut = 28.5E-6;
    anal.th_y_hcut_L = 102E-6; anal.th_y_hcut_R = 102E-6; anal.th_y_hcut = 100E-6;
}

#endif
