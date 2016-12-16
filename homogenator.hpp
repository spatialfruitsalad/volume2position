//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
#ifndef HOMOGENATOR_H_123456
#define HOMOGENATOR_H_123456
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <limits>
#include "tomogram.hpp"
#include "cSummerLinear.h"

class homogenator
{
    static constexpr unsigned int NumberOfBins = 40;
public:
    template <class VALUE>
    void Process (tomogram3d<VALUE>& matrix, tomogram3d<VALUE>& out, std::string foldername, double outside = 0.05)
    {
        std::cout << "radial homogenator " << std::endl;
        std::vector<double> vec_r;  // vector for radial bins
        std::vector<double> vec_mean;   // vector for mean values (same length as vec_r)
        std::vector<int> vec_counts;    // vector for voxel counts (same length as vec_r)

        unsigned nx = matrix.get_sx();  
        unsigned ny = matrix.get_sy();
        unsigned nz = matrix.get_sz();

        // create empty tomogram
        out.clear();
        out.initialise(nx,ny,nz,0);

        // get center for radial calculations
        unsigned xh = nx/2;
        unsigned yh = ny/2;

        // some geometry
        double rwidth = sqrt(xh*xh + yh*yh);    // half diagonal
        double rdelta = rwidth/(NumberOfBins+1);    // number of bins


        /// initialize vectors with empty values
        for (unsigned int i = 0 ; i != NumberOfBins; ++i)
        {
            double rmid = rdelta * static_cast<double>(i);  // radial bin position
            vec_r.push_back(rmid);                          
            vec_mean.push_back(0);
            vec_counts.push_back(0);
        }


        std::cout << "calculating mean values in " << NumberOfBins << " bins" << std::endl;
        // fill vectors
        for (unsigned int k = 0; k != nz; ++k)
        for (unsigned int j = 0; j != ny; ++j)
        for (unsigned int i = 0; i != nx; ++i)
        {
            // get bin 
            double r = sqrt( (i-nx/2)*(i-nx/2) + (j- ny/2)*(j-ny/2));
            unsigned idx = floor(r/rdelta);   

            // get value from tomogram
            VALUE v = matrix.get(i,j,k);
            
            // increase counters
            vec_mean[idx]+=v;
            vec_counts[idx]+=1;
        }

        std::cout << "normalize mean values " << std::endl;
        // correct for different voxel count numbers
        for (unsigned i = 0; i != vec_mean.size(); ++i)
        {
            vec_mean[i] /=static_cast<double>(vec_counts[i]);   
        }

        std::cout << "save homogenator.dat" << std::endl;
        // save output file
        {
        std::ofstream meanFile;
        std::string meanFileName = foldername + "/homogenator.dat";
        meanFile.open (meanFileName); 
        meanFile << "#rmid\tmean\tcounts" << std::endl;
        for (unsigned i = 0; i != vec_mean.size(); ++i)
        {
            meanFile << vec_r[i] << "\t" << vec_mean[i] << "\t" << vec_counts[i] << std::endl;
        }
        }

        // histogram class
        cSummerLinear summer(32);
        // maximum value: we want to be between 0 and vMAX/2 (not using the last bit to prevent intersample clipping in further calculations) 
        VALUE vMAX = std::numeric_limits<VALUE>::max();
        // adapt homogenatorization
        // and calculate histogram
        std::cout << "homogenize Tomogram" << std::endl;
        for (unsigned int k = 0; k != nz; ++k)
        for (unsigned int j = 0; j != ny; ++j)
        for (unsigned int i = 0; i != nx; ++i)
        {
            // omit outside voxels radially
            double r = sqrt( (i-nx/2)*(i-nx/2) + (j- ny/2)*(j-ny/2));
            if (r > 0.5*nx*(1.0 - outside))
            {
                out.set(i,j,k, 0);
                continue;
            }


            // homogenize tomogram
            unsigned idx = floor(r/rdelta);   
            double factor = 0;
            // no interpolation in the last bins (not interesting anyway)
            if( idx >= vec_mean.size() -2)
            {
                factor = vec_mean[idx];
            }
            else
            {
                // linear interpolation
                double inter = r/rdelta - idx;  // between 0 and 1
                factor = vec_mean[idx] * (1.0 - inter) + vec_mean[idx+1] * inter; 
            }

            // calculate new value
            VALUE v = matrix.get(i,j,k);
            VALUE newValue = static_cast<VALUE>(static_cast<double>(v) / factor *static_cast<double>(vMAX/2) );

            out.set(i,j,k,newValue);
            summer.Add(newValue);
        }
        {
        std::ofstream sumFile;
        std::string sumFileName = foldername + "/histogram_hom.dat";
        sumFile.open (sumFileName); 

        sumFile << summer;
        }

    };

};

#endif
