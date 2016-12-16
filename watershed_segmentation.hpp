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
#ifndef WATERSHEDSEGMENTATION_HPP
#define WATERSHEDSEGMENTATION_HPP

#include "tomogram.hpp"
#include "gauss_filter.hpp"
#include "erodate.hpp"
#include "dilate.hpp"
#include "edm.hpp"
#include "hk.hpp"
#include "common.hpp"


class watershedSegmentation
{
    // helper function 1
    static int lookup_2 (int x, std::vector<int> & rename)
    {
        assert (x >= 0);
        int ret = rename.at(x);
        if (ret < 0)
            return lookup_2 (-ret, rename);
        else
            return ret;
    }

    // helper function 2
    static int find_label (unsigned long nx, unsigned long ny, unsigned long nz, unsigned long i, unsigned long j, unsigned long k, tomogram3d<int>& clustermap, tomogram3d<float>& edm )
    {
        unsigned long imin = std::max (long(i -1), long(0));
        unsigned long imax = std::min (long(i +1), long (nx)-1);
        unsigned long jmin = std::max (long(j -1), long(0));
        unsigned long jmax = std::min (long(j +1), long (ny)-1);
        unsigned long kmin = std::max (long(k -1), long(0));
        unsigned long kmax = std::min (long(k +1), long (nz)-1);


        float grad = 0.;
        long grad_label = 0;
        unsigned long ii_s = 0, jj_s = 0, kk_s = 0;
        for (unsigned long ii = imin; ii <= imax; ++ii)
            for (unsigned long jj = jmin; jj <= jmax; ++jj)
                for (unsigned long kk = kmin; kk <= kmax; ++kk)
                {
                    float diff = edm.get(ii,jj,kk) - edm.get(i,j,k);
                    //float distance = sqrt(pow(float(ii-i),2)+pow(float(jj-j),2)+pow(float(kk-k),2));
                    float distance = sqrt(float((ii-i)*(ii-i) + (jj-j)*(jj-j) + (kk-k)*(kk-k)));
                    if (distance == 0) continue;
                    float slope = diff/distance;
                    if (slope > grad)
                    {
                        grad = slope;
                        grad_label = clustermap.get(ii,jj,kk);
                        ii_s = ii;
                        jj_s = jj;
                        kk_s = kk;
                    }
                }
        if (grad < 0.0001)
            return 0;
        else if (grad_label == 0)
            return find_label(nx,ny,nz,ii_s,jj_s,kk_s, clustermap, edm);
        else
            return grad_label;

    }
public:



    // pass this the binary image, put out a labelled image
    static void Process (tomogram3d<uint8_t>& binary,  tomogram3d<int>& out, float erosion_depth, std::string imagefoldername, int minimumCluster= 50000)
    {
        std::cout << "start segmentation ... " << std::endl;
        int nx = binary.get_sx();
        int ny = binary.get_sy();
        int nz = binary.get_sz();

        std::cout << "calculate edm" << std::endl;

        tomogram3d<float> edm;
        edm::create_edm(binary, edm);

        std::cout << "sqrt edm ..." << std::flush;
        for (int z = 0; z<nz; z++)
        for (int y = 0; y<ny; y++)
        for (int x = 0; x<nx; x++){
            float value = sqrt( edm.get(x,y,z) );
            edm.set(x,y,z,value);
        }
        std::cout << " done" << std::endl;

        // erodate binary image so particles don't overlap anymore. erosiondepth has to be chosen respectively
        tomogram3d<uint8_t> erodated_binary = binary;
        Erodate::erodate_binary_sphere(erodated_binary, edm, erosion_depth);
        WriteImages(erodated_binary, imagefoldername, "5_erosion", nz);
        
        // filter edm for smoother look
        gaussFilter::gauss_filter_mask_5(edm , 5, 0.8);
        
        WriteImages(edm, imagefoldername, "6_edm", nz);

        // identify individual clusters
        HoshenKopelman::Process(erodated_binary, out);
        erodated_binary.clear();
        // rename vector for HoshenKopelman
        std::vector <int> rename;

        rename.clear ();
        rename.push_back(0);

        // find max
        int ellipsoid_number = out.find_max();
        std::cout << "ellipsoid_number = " << ellipsoid_number << std::endl;
        for (int i = 1; i<= ellipsoid_number; i++)
            rename.push_back(i);

        // label previously eroded voxels
        int init = ellipsoid_number;
        for (int k = 0; k != nz; ++k)
            for (int j = 0; j != ny; ++j)
                for (int i = 0; i != nx; ++i)
                    if (binary.get(i,j,k) == 1 && out.get(i,j,k) == 0)
                    {
                        out.set (i, j, k, init);
                        rename.push_back(init);
                        init++;
                    }
    
        // check for slope in filtered edm
        for (int k = 0; k != nz; ++k)
        {
            for (int j = 0; j != ny; ++j)
            for (int i = 0; i != nx; ++i)
            {
                if (binary.get(i,j,k) == 1 && out.get(i,j,k) > ellipsoid_number)
                {
                    int imin = std::max (i-1, 0);
                    int imax = std::min (i+1, int (nx-1));
                    int jmin = std::max (j-1, 0);
                    int jmax = std::min (j+1, int (ny-1));
                    int kmin = std::max (k-1, 0);
                    int kmax = std::min (k+1, int (nz-1));

                    float grad = 0.;
                    int label = out.get(i,j,k);
                    int grad_label = out.get(i,j,k);

                    for (int ii = imin; ii <= imax; ++ii)
                    for (int jj = jmin; jj <= jmax; ++jj)
                    for (int kk = kmin; kk <= kmax; ++kk)
                    {
                        float diff = edm.get(ii,jj,kk) - edm.get(i,j,k);
                        //sweis @Fabian: Why??? pow is so slow 
                        //float distance = sqrt(pow(float(ii-i),2)+pow(float(jj-j),2)+pow(float(kk-k),2));
                        float distance = sqrt( (ii-i)*(ii-i) + (jj-j)*(jj-j) + (kk-k)*(kk-k));
                        if (distance == 0) continue;
                        float slope = diff/distance;
                        if (slope > grad)
                        {
                            grad = slope;
                            grad_label = out.get(ii,jj,kk);
                        }
                    }
                    if (grad_label != label) 
                    {
                        rename.at(label) = - grad_label;
                    }
                    else
                    {
                        rename.at(label) = grad_label;
                    }
                }
            }
        }


        int clustercount = 1;
        for (size_t i = 0; i != rename.size (); ++i)
        {
            if (rename.at(i) > 0)
            {
                rename.at(i) = clustercount;
                clustercount++;
            }
        }

        std::cout << "start checking rename vektor" << std::endl;
        int not_added = 0;
        for (unsigned int i = 1; i != rename.size (); ++i)
        {
            if (rename.at(i) == 0)
            {
                not_added++;
            }
        }
        std::cout << "not_added = " << not_added << std::endl;


        std::cout << "ellipsoid_number = " << ellipsoid_number << std::endl;
        std::cout << "clustercount = " << clustercount << std::endl;


        // resolve all indirections in the rename map
        for (int i = 0; i != (int)rename.size (); ++i)
        {
            if (watershedSegmentation::lookup_2 (i,rename) > ellipsoid_number)
                rename.at(i) = 0;
            else
                rename.at(i) = watershedSegmentation::lookup_2 (i, rename);
        }

        for (int k = 0; k != nz; ++k)
        for (int j = 0; j != ny; ++j)
        for (int i = 0; i != nx; ++i)
        {
            int newvalue = rename.at(out.get(i,j,k));
            out.set(i,j,k, newvalue);
        }

        std::cout << "                             done" << std::endl;

        // count the number of voxels in each cluster
        std::vector<int> clustersize (ellipsoid_number + 1);
        for (int k = 0; k != nz; ++k)
        for (int j = 0; j != ny; ++j)
        for (int i = 0; i != nx; ++i)
            clustersize.at(out.get(i,j,k))++;

        // remove clusters smaller than a certain threshold
        for (int k = 0; k != nz; ++k)
        for (int j = 0; j != ny; ++j)
        for (int i = 0; i != nx; ++i)
        {
            if (clustersize.at(out.get(i,j,k)) < minimumCluster) out.set(i,j,k,0);
        }

        erodated_binary = binary;
        Dilate::Process(erodated_binary, 3);
        edm::create_edm(erodated_binary, edm);

        int added_later = 0;

        for (int k = 0; k != nz; ++k)
        for (int j = 0; j != ny; ++j)
        for (int i = 0; i != nx; ++i)
        {
            if (binary.get(i,j,k) == 1 && out.get(i,j,k) == 0)
            {
                out.set(i,j,k,find_label(nx, ny, nz, i, j, k, out, edm));
                added_later++;
            }
        }
        std::cout << "added_later = " << added_later << std::endl;
        edm.clear();
    };
};

#endif //SEGMENTATION_HPP
