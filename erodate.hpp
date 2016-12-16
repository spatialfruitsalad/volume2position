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
#ifndef ERODATE
#define ERODATE

#include <iostream>
#include <limits>

#include "tomogram.hpp"

#include "edm.hpp"

class Erodate
{
public:
    // erosion with a filter mask
    static void erodate2(tomogram3d<uint8_t>& matrix_in, int masksize, tomogram3d<uint8_t>& matrix_out)
    {
        std::cout << "erodate: masksize " << masksize << "  ... " << std::flush;

        int sx = matrix_in.get_sx();
        int sy = matrix_in.get_sy();
        int sz = matrix_in.get_sz();
        matrix_out.clear();
        matrix_out.initialise (sx, sy, sz, 0);

        int tol = (masksize-1)/2;

        for (int k = 0; k != sz; ++k)
        for (int j = 0; j != sy; ++j)
        for (int i = 0; i != sx; ++i)
        {
            if (matrix_in.get(i,j,k) != 0) 
            {
                int imin = std::max (i-tol, 0);
                int imax = std::min (i+tol, int (sx-1));
                int jmin = std::max (j-tol, 0);
                int jmax = std::min (j+tol, int (sy-1));
                int kmin = std::max (k-tol, 0);
                int kmax = std::min (k+tol, int (sz-1));

                uint8_t test_min = std::numeric_limits<uint8_t>::max();
                for (int ii = imin; ii <= imax; ++ii)
                for (int jj = jmin; jj <= jmax; ++jj)
                for (int kk = kmin; kk <= kmax; ++kk)
                {
                    test_min = std::min(test_min, matrix_in.get(ii,jj,kk));
                }

                matrix_out.set(i,j,k,test_min);
            }
        }
        std::cout << " done" << std::endl;
    };

    // uses an already calculated edm 
    static void erodate_binary_sphere(tomogram3d<uint8_t>& matrix_in, tomogram3d<float>& edm, float depth)
    {
        std::cout << "erodate sphere: erosion depth " << depth << "  ... " << std::flush;

        int sx = matrix_in.get_sx();
        int sy = matrix_in.get_sy();
        int sz = matrix_in.get_sz();

        for (int k = 0; k != sz; ++k)
        for (int j = 0; j != sy; ++j)
        for (int i = 0; i != sx; ++i)
            if (edm.get(i,j,k) < depth + 0.001 )
                matrix_in.set(i,j,k,0);

        std::cout << " done" << std::endl;
    };

    // just an overload if no edm has been calculated.
    static void erodate_binary_sphere(tomogram3d<uint8_t>& matrix_in, float depth)
    {
        std::cout << "erodate sphere: erosion depth " << depth << "  ... " << std::flush;
        tomogram3d<float> edm;
        edm::create_edm(matrix_in, edm);
        erodate_binary_sphere(matrix_in, edm, depth);
    };
};

//void dilate_binary(tomogram3d<uint8_t>* matrix_in, int masksize);

#endif //ERODATE
