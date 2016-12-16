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
#ifndef COMMON_H_GUARD_123456
#define COMMON_H_GUARD_123456

#include "tomogram.hpp"

// helper function to print automatically three slices of a tomogram, at 10%,50/% and 90% of the total height of the tomogram stack
template <typename VALUE>
inline void WriteImages(tomogram3d<VALUE>& t, std::string path, std::string name, unsigned long nz, bool randomcolors = false)
{
    unsigned long sliceTenPercent = nz/10;
    unsigned long sliceFiftyPercent = nz/2;
    unsigned long sliceNinetyPercent = nz - sliceTenPercent;

    std::string Tenstring = path +    std::to_string(sliceTenPercent)   + "_" + name  + ((randomcolors)? ".ppm" : ".pgm");
    std::string Fiftystring = path +  std::to_string(sliceFiftyPercent) + "_" + name  + ((randomcolors)? ".ppm" : ".pgm");
    std::string Ninetystring = path + std::to_string(sliceNinetyPercent)+ "_" + name  + ((randomcolors)? ".ppm" : ".pgm");

    if(randomcolors)
    {
        t.write_ppm(Tenstring, sliceTenPercent);
        t.write_ppm(Fiftystring, sliceFiftyPercent);
        t.write_ppm(Ninetystring, sliceNinetyPercent);
    }
    else
    {
        t.write_pgm(Tenstring, sliceTenPercent);
        t.write_pgm(Fiftystring, sliceFiftyPercent);
        t.write_pgm(Ninetystring, sliceNinetyPercent);
    }
};

#endif
