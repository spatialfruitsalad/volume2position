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
#ifndef GAUSS_FILTER
#define GAUSS_FILTER

#include <iostream>
#include <cmath>
#include "tomogram.hpp"


class gaussFilter
{
public:
    gaussFilter(unsigned int masksize, double sigma) : m_masksize(masksize)
    {
        // TODO initial checks
        if(masksize <= 1)
        {
            throw "gauss filter: masksize must be larger than 1!";
        }
        if (sigma <= 0)
        {
            throw "guass filter: sigma must be larger than 0";
        }

        std::cout << "creating gauss Filter with parameters:" << std::endl;
        std::cout << "\tmasksize = " << masksize << "\n" << "\tsigma = " << sigma << std::endl;

        // create mask
        tol = (m_masksize-1)/2;
        mask1D.resize(masksize);
        double sum = 0;

        for (unsigned int i = 0; i != masksize; ++i)
        {
            double distance_square = double(i-tol) * double(i-tol);
            double val = exp(-1.0*distance_square/(2.0*sigma*sigma));
            mask1D.at(i) = val;
            sum += val;
        }
        // normalize mask
        for (unsigned int i = 0; i != masksize; ++i)
        {
            mask1D.at(i) = mask1D.at(i)/sum;
        }

        std::cout << "filter created" << std::endl;
    };


    template <typename VALUE>
    void Process(tomogram3d<VALUE>& matrix_in)
    {
        std::cout << "applying gauss filter" << std::endl;
        
        int sx = matrix_in.get_sx();
        int sy = matrix_in.get_sy();
        int sz = matrix_in.get_sz();
        
        std::cout << "\tPerforming Filter in X" << std::endl;
        std::vector<VALUE> row_x (sx);
        for (int k = 0; k != sz; ++k)
            for (int j = 0; j != sy; ++j)
            {
                for (int i = 0; i != sx; ++i)
                {
                    row_x[i] = matrix_in.get(i,j,k);
                }
                for (int i = 0; i != sx; ++i)
                {
                    int imin = std::max (i-tol, 0);
                    int imax = std::min (i+tol, int (sx-1));
                    double new_val = 0;
                    for (int ii = imin; ii <= imax; ++ii)
                    {
                        new_val += row_x.at(ii) * mask1D.at(ii-i+tol);
                    }
                    matrix_in.set(i,j,k,new_val);
                }
            }

        std::cout << "\tPerforming Filter in Y" << std::endl;
        std::vector<VALUE> row_y (sy);
        for (int k = 0; k != sz; ++k)
            for (int i = 0; i != sx; ++i)
            {
                for (int j = 0; j != sy; ++j)
                {
                    row_y[j] = matrix_in.get(i,j,k);
                }
                for (int j = 2; j != sy - 2; ++j)
                {
                    int jmin = std::max (j-tol, 0);
                    int jmax = std::min (j+tol, int (sy-1));
                    double new_val = 0;
                    for (int jj = jmin; jj <= jmax; ++jj)
                    {
                        new_val += row_y.at(jj)*mask1D.at(jj-j+tol);
                    }
                    matrix_in.set(i,j,k,new_val);
                }
            }

        std::cout << "\tPerforming Filter in Z" << std::endl;
        std::vector<VALUE> row_z (sz);
        for (int i = 0; i != sx; ++i)
            for (int j = 0; j != sy; ++j)
            {
                for (int k = 0; k != sz; ++k)
                {
                    row_z[k] = matrix_in.get(i,j,k);
                }
                for (int k = 2; k != sz - 2; ++k)
                {
                    int kmin = std::max (k-tol, 0);
                    int kmax = std::min (k+tol, int (sz-1));
                    double new_val = 0;
                    for (int kk = kmin; kk <= kmax; ++kk)
                    {
                        new_val += row_z.at(kk)*mask1D.at(kk-k+tol);
                    }
                    matrix_in.set(i,j,k,new_val);
                }
            }

        std::cout << "done" << std::endl;
    };





    // fast version of gaussian filter, which works on the three axis separately
    template <typename VALUE>
    static void gauss_filter_mask_5 (tomogram3d<VALUE>& matrix_in, int masksize, float sigma){

    std::cout << "gauss_filter_neu: masksize " << masksize << "  ... " << std::flush;

    int sx = matrix_in.get_sx();
    int sy = matrix_in.get_sy();
    int sz = matrix_in.get_sz();

    int tol = (masksize-1)/2;

    std::vector <float> mask1D (masksize);
    float sum = 0;
    
    for (int i = 0; i != masksize; ++i){
        float distance_square = pow(float(i-tol),2);
        float val = exp(-1*distance_square/(2*sigma*sigma));
        mask1D.at(i) = val;
        sum += val;
    }

    for (int i = 0; i != masksize; ++i){
        mask1D.at(i) = mask1D.at(i)/sum;
    }

    std::vector<VALUE> row_x (sx);
    for (int k = 0; k != sz; ++k)
    for (int j = 0; j != sy; ++j){
        for (int i = 0; i != sx; ++i){
            row_x.at(i) = matrix_in.get(i,j,k);
        }
        for (int i = 0; i != sx; ++i){
            int imin = std::max (i-tol, 0);
            int imax = std::min (i+tol, int (sx-1));
            float new_val = 0;
            for (int ii = imin; ii <= imax; ++ii){
                new_val += row_x.at(ii)*mask1D.at(ii-i+tol);
            }
            matrix_in.set(i,j,k,new_val);
        }
    }

    std::vector<VALUE> row_y (sy);
    for (int k = 0; k != sz; ++k)
    for (int i = 0; i != sx; ++i){
        for (int j = 0; j != sy; ++j){
            row_y.at(j) = matrix_in.get(i,j,k);
        }
        for (int j = 2; j != sy - 2; ++j){
            int jmin = std::max (j-tol, 0);
            int jmax = std::min (j+tol, int (sy-1));
            float new_val = 0;
            for (int jj = jmin; jj <= jmax; ++jj){
                new_val += row_y.at(jj)*mask1D.at(jj-j+tol);
            }
            matrix_in.set(i,j,k,new_val);
        }
    }

    std::vector<VALUE> row_z (sz);
    for (int i = 0; i != sx; ++i)
    for (int j = 0; j != sy; ++j){
        for (int k = 0; k != sz; ++k){
            row_z.at(k) = matrix_in.get(i,j,k);
        }
        for (int k = 2; k != sz - 2; ++k){
            int kmin = std::max (k-tol, 0);
            int kmax = std::min (k+tol, int (sz-1));
            float new_val = 0;
            for (int kk = kmin; kk <= kmax; ++kk){
                new_val += row_z.at(kk)*mask1D.at(kk-k+tol);
            }
            matrix_in.set(i,j,k,new_val);
        }
    }
    std::cout << " done" << std::endl;
}
private:
    unsigned int m_masksize;
    int tol;
    //double m_sigma;
    std::vector <double> mask1D;
};

#endif //GAUSS_FILTER
