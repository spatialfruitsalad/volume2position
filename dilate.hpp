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
#ifndef DILATE_GUARD_H_123456
#define DILATE_GUARD_H_123456

class Dilate
{
public:
    static void Process(tomogram3d<uint8_t>& matrix_in, int masksize)
    {
        std::cout << "dilate: masksize " << masksize << "  ... " << std::flush;

        tomogram3d<uint8_t> dilate_map;

        int sx = matrix_in.get_sx();
        int sy = matrix_in.get_sy();
        int sz = matrix_in.get_sz();

        std::cout << "initialize dilate_map " << "  ... " << std::flush;
        dilate_map.initialise (sx, sy, sz, 1);
        std::cout << " done ... " << std::flush;

        int tol = (masksize-1)/2;

        for (int k = 0; k != sz; ++k)
        {
            for (int j = 0; j != sy; ++j)
            {
                for (int i = 0; i != sx; ++i)
                {
                    if (matrix_in.get(i,j,k) != 1)
                    {
                        int imin = std::max (i-tol, 0);
                        int imax = std::min (i+tol, int (sx-1));
                        int jmin = std::max (j-tol, 0);
                        int jmax = std::min (j+tol, int (sy-1));
                        int kmin = std::max (k-tol, 0);
                        int kmax = std::min (k+tol, int (sz-1));

                        uint8_t test_max = 0;
                        for (int ii = imin; ii <= imax; ++ii)
                            for (int jj = jmin; jj <= jmax; ++jj)
                                for (int kk = kmin; kk <= kmax; ++kk)
                                {
                                    test_max = std::max(test_max, matrix_in.get(ii,jj,kk));
                                    if (test_max == 1) goto end_voxel;
                                }
end_voxel:
                        dilate_map.set(i,j,k,test_max);
                    }
                }
            }
        }
        std::swap(matrix_in, dilate_map);
        std::cout << " done" << std::endl;
    }
};
#endif
