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
#ifndef CSEGMENTATIONBYOTS_H_12345
#define CSEGMENTATIONBYOTS_H_12345

#include <vector>

class cSegmentationByOtsu 
{
    public:
    static void NormalizeHisto (std::vector<double> & vals)
    {
        double sum = 0;
        for (double v : vals)
        {
            sum += v;
        }
        for(size_t i = 0; i != vals.size(); ++i)
        {
            vals[i] /= sum;
        }
    };
        
    // pass this function a histogram of grey values
    static double CalculateSegmentationThreshold (std::vector<double> const & _vals)
    {
        // main code by fabian
        long GRAYLEVEL = _vals.size();
        double prob[GRAYLEVEL], omega[GRAYLEVEL]; /* prob of graylevels */
        double myu[GRAYLEVEL];   /* mean value for separation */
        double max_sigma, sigma[GRAYLEVEL]; /* inter-class variance */
        long i; /* Loop variable */
        int threshold; /* threshold for binarization */

        std::cout << "Otsu's binarization process starts now." << std::endl;
        std::cout << "Working on a histogram with " << GRAYLEVEL << " bins" << std::endl;


        /* calculation of probability density */
        for ( i = 0; i < GRAYLEVEL; i ++ ) {
            prob[i] = _vals.at(i);
        }

        /* omega & myu generation */
        omega[0] = prob[0];
        myu[0] = 0.0;       /* 0.0 times prob[0] equals zero */
        for (i = 1; i < GRAYLEVEL; i++) {
            omega[i] = omega[i-1] + prob[i];
            myu[i] = myu[i-1] + i*prob[i];
        }

        /* sigma maximization
           sigma stands for inter-class variance 
           and determines optimal threshold value */
        threshold = 0;
        max_sigma = 0.0;
        for (i = 0; i < GRAYLEVEL-1; i++) 
        {
            if (omega[i] != 0.0 && omega[i] != 1.0)
            {
                sigma[i] = pow(myu[GRAYLEVEL-1]*omega[i] - myu[i], 2) / (omega[i]*(1.0 - omega[i]));
            }
            else
            {
                sigma[i] = 0.0;
            }
            if (sigma[i] > max_sigma) 
            {
                max_sigma = sigma[i];
                threshold = i;
            }
        }
        std::cout << std::endl << "threshold value ="  << threshold << std::endl;
        return threshold;  
    }

    private:
    // don't call this, just use the static method
    cSegmentationByOtsu() {};


};

#endif
