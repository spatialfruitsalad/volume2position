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
#ifndef REMOVE_INSIDE_HOLES_2_HPP
#define REMOVE_INSIDE_HOLES_2_HPP

#include "tomogram.hpp"

class removeInsideHoles
{
public:
    static int lookup_in_2 (int x, std::vector<int> & rename)
    {
        int ret = rename.at(x);
        if (ret < 0)
            return lookup_in_2 (-ret, rename);
        else
            return ret;
    };
    static int merge_clusters_in_2 (int a, int b, std::vector<int> & rename)
    {
        if (a > b)
            return merge_clusters_in_2 (b, a, rename);
        rename.at (b) = -a;
        return a;
    };

    static bool voxset_in_2 (int i, int j, int k, tomogram3d <uint8_t>& matrix_in)
    {
        return !matrix_in.get(i,j,k);
    };

    // arguments:
    // matrix_in  - binary image to fill holes
    static void Process (tomogram3d <uint8_t> & matrix_in)
    {
        std::cout << "removing inside holes ... " << std::flush;

        tomogram3d<int> clustermap;

        std::vector <int> rename;

        int nx = matrix_in.get_sx();
        int ny = matrix_in.get_sy();
        int nz = matrix_in.get_sz();

        // initialise
        rename.clear ();

        clustermap.initialise (nx, ny, nz, -1);

        // this is basically hoshen kopelman for the black clusters
        for (int k = 0; k != nz; ++k)
        {
            for (int j = 0; j != ny; ++j)
                for (int i = 0; i != nx; ++i)
                {
                    if (voxset_in_2 (i,j,k,matrix_in))
                    {
                        int imin = std::max (i-1, 0);
                        int imax = std::min (i+1, int (nx-1));
                        int jmin = std::max (j-1, 0);
                        int jmax = std::min (j+1, int (ny-1));
                        int kmin = std::max (k-1, 0);
                        int kmax = std::min (k+1, int (nz-1));
                        int cluster_id = rename.size (); // possibly start a new cluster.


                        int pos[6][3] = {{imin,j,k},{imax,j,k},{i,jmin,k},{i,jmax,k},{i,j,kmin},{i,j,kmax}};
                        for (int a = 0; a < 6; ++a)
                        {
//                        std::cout << "a = " << a << std::endl;
                            int ii = pos[a][0];
                            int jj = pos[a][1];
                            int kk = pos[a][2];
//                        std::cout << "ii = " << ii << "  jj = " << jj << "  kk = " << kk << std::endl;
                            if (voxset_in_2 (ii,jj,kk,matrix_in) && clustermap.get(ii,jj,kk) != -1)
                            {
                                int oclid = lookup_in_2 (clustermap.get(ii,jj,kk),rename);
                                if (cluster_id == (int)rename.size ())
                                    cluster_id = oclid;
                                else if (cluster_id != oclid)
                                    cluster_id = merge_clusters_in_2 (oclid, cluster_id, rename);
                            }
                        }



                        // start a new cluster
                        if (cluster_id == (int)rename.size ())
                            rename.push_back (cluster_id);

                        clustermap.set(i,j,k,cluster_id);
                    }
                }
        }



        // resolve all indirections in the indirection map
        for (int i = 0; i != (int)rename.size (); ++i)
            rename[i] = lookup_in_2 (i, rename);
        int max_clustercount = 0;
        // resolve the indirections on the facing planes
        // in the cluster map (other data is not needed)
        rename.insert (rename.begin (), -1);
        for (int i = 0; i != nx; ++i)
            for (int j = 0; j != ny; ++j)
                for (int k = 0; k != nz; ++k)
                {
                    if(rename[clustermap.get(i,j,k)+1] > max_clustercount) max_clustercount = rename[clustermap.get(i,j,k)+1];
                    clustermap.set(i,j,k, rename[clustermap.get(i,j,k)+1]);
                }
        rename.erase (rename.begin ());

        std::vector<int> clustersizes (max_clustercount + 1);
        for (int k = 0; k != nz; ++k)
            for (int j = 0; j != ny; ++j)
                for (int i = 0; i != nx; ++i)
                {
                    if (clustermap.get(i,j,k) >= 0)
                    {
                        clustersizes.at(clustermap.get(i,j,k))++;
                    }
                }


        tomogram3d <uint8_t> matrix_out;
        matrix_out.initialise(nx,ny,nz,0);
        for (int k = 0; k != nz; ++k)
        for (int j = 0; j != ny; ++j)
        for (int i = 0; i != nx; ++i)
        {
            bool fill = false;
            if (clustermap.get(i,j,k) == -1) 
            {
                fill = true;
            }

            else if (clustermap.get(i,j,k) > 0 && clustersizes.at(clustermap.get(i,j,k)) < 100000) 
            {
                fill = true;
            }

            if(fill == true)
            {
                matrix_out.set(i,j,k,1);
            }
        }


        std::cout << " done" << std::endl;

        std::swap (matrix_in, matrix_out);
    }
};

#endif //REMOVE_INSIDE_HOLES_2_HPP
