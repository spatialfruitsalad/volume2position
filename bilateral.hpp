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
#ifndef BILATERAL_H_GUARD_123456
#define BILATERAL_H_GUARD_123456
#include "tomogram.hpp"

////////////////////////////
// bilateral filter implementation as described in: 1998 by C. Tomasi et al. DOI: 10.1109/ICCV.1998.710815
////////////////////////////
class bilateralFilter
{
public:
    bilateralFilter(unsigned int masksize, double sigma, double similarity) : m_masksize(masksize), m_sigma(sigma), m_sigmaSimilarity(similarity)
    {
        std::cout << "creating bilateral Filter with parameters:" << std::endl;
        std::cout << "\tmasksize = " << masksize << "\n\tsigma = " << sigma << "\n\tsimilarity = " << similarity;
        // initial checks
        if(masksize <= 1)
        {
            throw std::string("bilateral filter: masksize must be larger than 1!");
        }
        if (sigma <= 0 || similarity <= 0)
        {
            throw std::string("bilateral filter: sigma/similarity must be larger than 0");
        }


        // create mask
        tol = (m_masksize-1)/2;
        std::cout << "\n\ttol = " << tol << std::endl;
        mask.initialise (masksize, masksize, masksize, 0);
        double sum = 0;
        
        // precalculate gaussian mask (this speeds up the process alot
        for (unsigned int k = 0; k != masksize; ++k)
        {
            for (unsigned int j = 0; j != masksize; ++j)
            {
                for (unsigned int i = 0; i != masksize; ++i)
                {
                    double distance_square = (i-tol)*(i-tol) + (j-tol)*(j-tol) + (k-tol)*(k-tol);
                    double val = exp(-1*distance_square/(2*m_sigma*m_sigma));
                    mask.set(i,j,k, val);
                    sum += val;
                }
            }
        }

        // creating gauss table (also a big speedup)
        double x0 = 0;
        double x1 = 4 * m_sigmaSimilarity; // dumb cutoff

        // prerendered version of the gauss
        unsigned int steps = 500;
        double dx = (x1 - x0)/static_cast<double>(steps + 1);
        for (double x = x0; x <= x1; x += dx)
        {
            gaussTable.push_back(gauss(m_sigmaSimilarity, x));
        }

    };

    // helper function to get gauss function values from prerendered gaussTable
    double getGauss (double v)
    {
        v = fabs (v);
        if (v >= 4.0*m_sigmaSimilarity)
            return 0;
        double dx = (4*m_sigmaSimilarity)/static_cast<double>(501);
        unsigned int idx = floor(v/dx);
        return gaussTable[idx];
    }

    // function to calculate the gauss
    static double gauss(double sigma, double v)
    {
        v = fabs(v);
        if (v > 4.0*sigma) return 0;
        double f = v/sigma;
        // i don't care about the prefactor right here, since the computation will have to calculate the weight anyway.
        return exp (-0.5 * (f * f)  );
    }

    // This is the function to do the work
    template <typename VALUE>
    void Process(tomogram3d<VALUE>& matrix_in, tomogram3d<VALUE>& matrix_out)
    {
        std::cout << "applying bilateral filter" << std::endl;

        int sx = matrix_in.get_sx();
        int sy = matrix_in.get_sy();
        int sz = matrix_in.get_sz();

        matrix_out.clear();
        matrix_out.initialise (sx, sy, sz, 0);

        // outer loop over the whole image
        for (int k = 0; k != sz; ++k)
        {
            std::cout << k  << " " << std::flush;
            for (int j = 0; j != sy; ++j)
            {
                for (int i = 0; i != sx; ++i)
                {
                    // original value from the old image
                    VALUE orig = matrix_in.get(i,j,k);

                    // calculate boundaries so I don't read from otside the image
                    int imin = std::max (i-tol, 0);
                    int imax = std::min (i+tol, int (sx-1));
                    int jmin = std::max (j-tol, 0);
                    int jmax = std::min (j+tol, int (sy-1));
                    int kmin = std::max (k-tol, 0);
                    int kmax = std::min (k+tol, int (sz-1));

                    double new_val = 0;
                    double sum = 0;
                    // inner loop over the filter kernel
                    for (int ii = imin; ii <= imax; ++ii)
                    {
                        for (int jj = jmin; jj <= jmax; ++jj)
                        {
                            for (int kk = kmin; kk <= kmax; ++kk)
                            {
                                VALUE other = matrix_in.get(ii,jj,kk);  
                                double similarityFactor = getGauss( other - orig);
                                double closenessFactor = mask.get(ii-i+tol, jj-j+tol, kk-k+tol);
                                double currentWeight = closenessFactor * similarityFactor;
                                sum += currentWeight;
                                new_val += other * currentWeight;
                            }
                        }
                    }
                    matrix_out.set(i,j,k,new_val/sum);
                }
            }
        }
        std::cout << "\ndone" << std::endl;
    };
private:
    unsigned int m_masksize;
    double m_sigma;
    double m_sigmaSimilarity;
    int tol;
    tomogram3d<double> mask;
    std::vector<double> gaussTable;
};
#endif
