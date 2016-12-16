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
#ifndef HOSHEN_KOPELMAN_HPP
#define HOSHEN_KOPELMAN_HPP

#include "tomogram.hpp"

class HoshenKopelman
{
public:
    static int lookup (int x, std::vector<int> & rename) {
        assert (x >= 0);
        int ret = rename.at(x);
        if (ret < 0)
            return lookup (-ret, rename);
        else
            return ret;
    }

    static int merge_clusters (int a, int b, std::vector<int> & rename) {
        assert (a != b);
        if (a > b)
            return merge_clusters (b, a, rename);
        rename.at (b) = -a;
        return a;
    }

    static bool voxset (int i, int j, int k, tomogram3d <uint8_t>& matrix_in) {
        return matrix_in.get(i,j,k);
    }

    static void Process (tomogram3d <uint8_t>& matrix_in, tomogram3d<int>& out) {
        // code by fabian.
        std::cout << "start hoshen_kopelman ... " << std::flush;

        std::vector <int> rename;

        int nx = matrix_in.get_sx();
        int ny = matrix_in.get_sy();
        int nz = matrix_in.get_sz();

        // initialise
        rename.clear ();
        rename.push_back (0);
        out.clear();
        out.initialise (nx, ny, nz, 0);

        for (int k = 0; k != nz; ++k){
            for (int j = 0; j != ny; ++j)
                for (int i = 0; i != nx; ++i){
                    if (voxset (i,j,k,matrix_in)) {
                        int imin = std::max (i-1, 0);
                        int imax = std::min (i+1, int (nx-1));
                        int jmin = std::max (j-1, 0);
                        int jmax = std::min (j+1, int (ny-1));
                        int kmin = std::max (k-1, 0);
                        int kmax = std::min (k+1, int (nz-1));
                        int cluster_id = rename.size (); // possibly start a new cluster.

                        // this is all very suboptimal, but coding it the
                        // _right_ way generates a 100-line if-then-else block.
                        for (int ii = imin; ii <= imax; ++ii)
                            for (int jj = jmin; jj <= jmax; ++jj)
                                for (int kk = kmin; kk <= kmax; ++kk) {
                                    if (voxset (ii,jj,kk,matrix_in) && out.get(ii,jj,kk) != 0) {
                                        int oclid = lookup (out.get(ii,jj,kk),rename);
                                        if (cluster_id == (int)rename.size ())
                                            cluster_id = oclid;
                                        else if (cluster_id != oclid)
                                            cluster_id = merge_clusters (oclid, cluster_id, rename);
                                    }
                                }

                        // start a new cluster
                        if (cluster_id == (int)rename.size ())
                            rename.push_back (cluster_id);

                        out.set(i,j,k,cluster_id);
                    }
                }
        }

        int newindex = 0;
        for (int i = 0; i != (int)rename.size (); ++i)
            if (rename.at(i) >= 0)
                rename.at(i) = newindex++;

        // resolve all indirections in the indirection map
        for (int i = 0; i != (int)rename.size (); ++i)
            rename[i] = lookup (i, rename);

        // resolve the indirections on the facing planes
        // in the cluster map (other data is not needed)
    //    rename.insert (rename.begin (), -1);
        for (int i = 0; i != nx; ++i)
        for (int j = 0; j != ny; ++j)
        for (int k = 0; k != nz; ++k)
            out.set(i,j,k, rename[out.get(i,j,k)]);

        std::cout << " done" << std::endl;
    }
};
#endif //HOSHEN_KOPELMAN_HPP
