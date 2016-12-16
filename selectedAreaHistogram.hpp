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
#ifndef SELECTEDAREAHISTOGRAM_HPP_GUARD_12345
#define SELECTEDAREAHISTOGRAM_HPP_GUARD_12345

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <numeric>

#include "tomogram.hpp"

struct area
{
    unsigned int z;
    unsigned int x1, x2;
    unsigned int y1, y2;

    area()
    {
        z = 0;
        x1 = x2 = y1 = y2 = 0;
    }

    unsigned int calculateArea()
    {
        return ((x2-x1)*(y2-y1));
    }
    inline bool contains(unsigned int i, unsigned int j, unsigned int k)
    {
        if (k != z) return false;
        if (i > x2 || i < x1) return false;
        if (j > y2 || j < y1) return false;
        return true;
    }
};

class selectedAreaHistogram
{
    std::vector<area> m_areas;

    

public:

    selectedAreaHistogram (unsigned long nx, unsigned long ny, unsigned long nz)
    {
        std::cout << "selectedAreaHistogram Constructor" << std::endl;
        for (unsigned long i = 0; i != nz; ++i)
        {
            area a;
            a.z = i;
            a.x1 = 0;
            a.y1 = 0;
            a.x2 = nx;
            a.y2 = ny;
            m_areas.push_back(a);
        }
    }

    template <typename VALUE>
    void Process (tomogram3d<VALUE>& t, std::vector<double>& binvec, std::vector<double>& histvec, std::string foldername, unsigned int bins = 150)
    {
        std::cout << "start processing tomogram" << std::endl;

        unsigned int nx, ny, nz;
        nx = t.get_sx();
        ny = t.get_sy();
        nz = t.get_sz();

        std::vector<VALUE> greyvalues;
        long double mean = 0.0;
        unsigned int pixelsCounted = 0;
        // get all selected areas and create a vector of all grey values
        for (unsigned long k = 0; k != nz; ++k)
        for (unsigned long j = 0; j != ny; ++j)
        for (unsigned long i = 0; i != nx; ++i)
        {
            // check if the current voxel is in one of the areas
            for(auto it = m_areas.begin(); it != m_areas.end(); ++it)
            {
                if((*it).contains(i,j,k))
                {
                    pixelsCounted++;           
                    VALUE v = t.get(i,j,k);
                    greyvalues.push_back(v);
                    mean += v;
                }
            }
        }
        std::cout << "voxels counted: " << pixelsCounted << " which is a fraction of " << static_cast<double>(pixelsCounted)/static_cast<double>((nx*ny*nz)) << " of all voxels" << std::endl;
        mean /= static_cast<double>(pixelsCounted);
        std::cout << "mean value of: " << mean << std::endl;

        auto min =  (*std::min_element(greyvalues.begin(), greyvalues.end()));
        auto max =  (*std::max_element(greyvalues.begin(), greyvalues.end()));
        std::cout << "min: " << static_cast<int>(min) << std::endl << "max: " << static_cast<int>(max) << std::endl;

        auto diff = max - min;
        // create a histogram of all selected grey values
        long double binsize = static_cast<long double>(diff)/static_cast<long double>(bins);

        binvec.resize(bins + 1, 0);
        for (size_t i = 0; i != binvec.size(); ++i)
        {
           double gv = min + i * binsize;
           binvec.at(i) = gv;
        }
        histvec.resize(bins + 1, 0);
        for(auto n : greyvalues)
        {
            long double u = static_cast<long double>(n - min);
            long double v = static_cast<long double>(binsize);
            size_t idx = static_cast<size_t>(u/v);
            histvec.at(idx) += 1;
        }

        {
        std::cout << "write selected Area Histogram" << std::endl;
        std::ofstream out(foldername + "/selectedAreaHist.dat");
        for (size_t i = 0; i != histvec.size(); ++i)
        {
            long double gv = binvec.at(i);
            auto f = histvec.at(i);
            out << gv << " " << f << std::endl;      
        }
        out.close();
        }
    };
};

#endif
