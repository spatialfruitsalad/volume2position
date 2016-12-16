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
#include <iostream>
#include <sys/stat.h>   // for mkdir
#include <string>
#include <vector>
#include <map>
#include "tomogram.hpp"
#include "watershed_segmentation.hpp"
#include "bilateral.hpp"
#include "homogenator.hpp"
#include "cSegmentationByOtsu.h"
#include "selectedAreaHistogram.hpp"
#include "remove_inside_holes_2.hpp"
#include "common.hpp"

int main( int argc, char *argv[] )
{

///////////////////////////////
// Parse command line arguments
///////////////////////////////
    if (argc != 5)
    {
        std::cerr << "wrong number of Command line arguments" << std::endl;
        std::cerr << "label2color [inFileName] [nx] [ny] [nz]";
        return 0;
    }

    std::string inFileName = argv[1];

    // the size of the tomogram
    unsigned long  nx = atoi(argv[2]);
    unsigned long  ny = atoi(argv[3]);
    unsigned long  nz = atoi(argv[4]);

///////////////////////////////
// Load Tomogram
///////////////////////////////
    std::cout << "loading tomogram from " << inFileName << std::endl;

    tomogram3d<uint16_t> original;
    original.read_raw(inFileName, "uint16_t", nx, ny, nz);

    std::cout << "tomogram successfully loaded" << std::endl;
    
    original.write_ppm(inFileName + "_51.ppm", 51);    
    original.write_ppm(inFileName + "_100.ppm", 100);    
    original.write_ppm(inFileName + "_255.ppm", 255);    
    original.write_ppm(inFileName + "_400.ppm", 400);    

    
    return 0;
}
