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

// a small helper function to create folders
void CreateFolders(std::string outfoldername)
{
    if(outfoldername.empty())
    {
        throw std::string("outfilepath is empty");
    }
    char lastCharOfFolder = *outfoldername.rbegin();  
    if (lastCharOfFolder != '/')
        outfoldername += '/';

    std::string imagefoldername= outfoldername + "images/";
    std::string thresfoldername= outfoldername + "threshold_histogram/";
   
    mkdir(outfoldername.c_str(),0755);    
    mkdir(imagefoldername.c_str(),0755);    
    mkdir(thresfoldername.c_str(),0755);    

}


int main( int argc, char *argv[] )
{

///////////////////////////////
// Parse command line arguments
///////////////////////////////
    if (argc != 9)
    {
        std::cerr << "wrong number of Command line arguments" << std::endl;
        std::cerr << "./ip_analyzing_gm\n[inFileName] \n[outFolderName] \n[nx] [ny] [nz] [erosiondepth] [sigma_g] [sigma_p]";
        return 0;
    }

    std::string inFileName = argv[1];
    std::string outFolderName = argv[2];

    CreateFolders(outFolderName);

    // the size of the tomogram
    unsigned long  nx = atoi(argv[3]);
    unsigned long  ny = atoi(argv[4]);
    unsigned long  nz = atoi(argv[5]);

    // erosion depth for watershed segmentation of the particles
    int erosiondepth = atoi(argv[6]);


    // bilateral filter parameters
    double bf_sigma_Pos = atof(argv[7]);
    double bf_sigma_Col = atof(argv[8]);
    //double bf_sigma_Pos = 1.75;
    //double bf_sigma_Col = 2000;
    // cutoff for spatial gauss
    unsigned bf_size = static_cast<unsigned>(bf_sigma_Pos * 2.5);

///////////////////////////////
// Load Tomogram
///////////////////////////////
    std::cout << "loading tomogram from " << inFileName << std::endl;

    tomogram3d<uint16_t> original;
    original.read_raw(inFileName, "uint16_t", nx, ny, nz);

    std::cout << "tomogram successfully loaded" << std::endl;
    
    // write some images of the original tomogram
    WriteImages(original, outFolderName + "/images/", "0_original", nz);
    
///////////////////////////////
// Apply bilateral filter
///////////////////////////////
    {
    std::cout << "applying filter" << std::endl;
    // create filter object
    bilateralFilter B (bf_size, bf_sigma_Pos, bf_sigma_Col);
    tomogram3d<uint16_t> filtered;
    // apply filter
    B.Process(original, filtered);
    // original data is no longer needed, so swap it with filtered and dismiss filtered.
    std::swap(original, filtered);
    filtered.clear();
    original.write_raw(outFolderName + "/bilateral.raw", "uint16_t");
    // write some slices of the bilateral tomogram
    WriteImages(original, outFolderName + "/images/", "1_bilateral", nz);
    }

///////////////////////////////
// Create Homogenized Tomogram
///////////////////////////////
    tomogram3d<uint16_t> hom;
    // class for homogenization
    homogenator H;
    H.Process(original, hom, outFolderName);    
    
    //original is no longer needed
    original.clear();
    // write some slices of the homogenized tomogram
    hom.write_raw(outFolderName + "/homogenized.raw", "uint16_t");
    WriteImages(hom, outFolderName + "/images/", "2_homogenized", nz);

///////////////////////////////
// Apply selected Area Histogram
///////////////////////////////
    std::cout << "calculate histogram" << std::endl;
    // this will calculate a histogram of the whole tomogram 
    selectedAreaHistogram sah(nx,ny, nz);

    // selected area histogram will fill two vectors: bins will contain the bin information (list of the bins, identified by the first grey value in each bin.), 
    std::vector<double> bins;
    // hist will include the corresponding greyvalue counts
    std::vector<double> hist;
    sah.Process(hom, bins, hist, outFolderName, 150);

    // the black outside region should not get counted so just remove that from the histogram
    hist[0] = 0;
///////////////////////////////
// Calculate threshold with Otsu's method
///////////////////////////////
    std::cout << "calculate Threshold" << std::endl;
    cSegmentationByOtsu::NormalizeHisto(hist);
    double T = cSegmentationByOtsu::CalculateSegmentationThreshold(hist);
    std::cout << "binthresh from otsu: " << T << std::endl; 
    // since otsu doesn't need to know anything about the bins it is acting on, it will just return the bin ID (this is T). 
    // To get the greyvalue threshold (t1_gv) the next line is needed.
    uint16_t t1_gv = static_cast<uint16_t>(bins.at(T));

    std::cout << "threshol from otsu: " << static_cast<int>(t1_gv) << std::endl;
    
    // this writes a gnuplot script for visualizing otsu thresholding and histogram calculation
    {
    std::cout << "write gnuplot script for histogram" << std::endl;
    std::ofstream out(outFolderName + "/plot_selectedAreaHistogram.gpl");
    out << "#!/usr/bin/gnuplot -persist"<< std::endl;
    out << "t_low = " << static_cast<int>(t1_gv)<< std::endl;
    out << "set xlabel 'greyvalue'" << std::endl;
    out << "set ylabel 'occurrence'" << std::endl;
    out << "set xr [1:]" << std::endl;
    out << "set arrow from t_low, 0 to t_low, graph(0, 1) nohead" << std::endl;
    out << "pl 'selectedAreaHist.dat' w l" << std::endl;
    out.close();
    }
///////////////////////////////
// Create Particle binray image
///////////////////////////////
    std::cout << std::endl << "make binarized versions" << std::endl;
    tomogram3d<uint8_t> binary_particles;
    binary_particles.initialise(nx,ny,nz,0);

    // loop over all voxels and compare with calculated threshold
    for (unsigned long k = 0; k != nz; ++k)
    for (unsigned long j = 0; j != ny; ++j)
    for (unsigned long i = 0; i != nx; ++i)
    {
        uint16_t v = hom.get(i,j,k);
        if (v >= t1_gv)
        { 
            binary_particles.set(i,j,k,1); 
        }
    }
   
    WriteImages(binary_particles,  outFolderName+"images/", "3_binary", nz);

///////////////////////////////
// Remove inside holes from particle image
///////////////////////////////
    tomogram3d<uint8_t> binary_water;
    binary_water.initialise(nx,ny,nz,0);
    removeInsideHoles::Process(binary_particles);
    WriteImages(binary_particles, outFolderName+"images/", "4_binary_woHolesInside", nz); 
    binary_particles.write_raw(outFolderName + "/particles_binary.raw", "uint8_t");

///////////////////////////////
// Create labeled image
///////////////////////////////
    tomogram3d<int> out;
    try
    {
        // this calculates the segmented tomogram
        watershedSegmentation::Process(binary_particles, out, erosiondepth, outFolderName + "/images/", 100);

        // remove all clusters that touch the border
        out.cleanBorder(0.075);

        // write some images 
        WriteImages(out,outFolderName+"images/","7_labeled" , nz);
        WriteImages(out,outFolderName+"images/","8_labeled_color" , nz, true);

    }
    catch (std::string& e)
    {
        std::cerr << e << std::endl;
        return 1;
    }

///////////////////////////////
// Calculate centroids and save tomogram (converted to save disk space)
///////////////////////////////
    // calculate centroids and volumes, also
    // convert to uint16_t for writing
    std::cout << "convert and calculate" << std::endl;
    std::map<unsigned long, unsigned int> volumes;
    std::map<unsigned long, double> xPos;
    std::map<unsigned long, double> yPos;
    std::map<unsigned long, double> zPos;
    
    tomogram3d<uint16_t> safe;
    safe.initialise(nx, ny, nz, 0);
    std::cout << "counting voxels" << std::endl;

    
    for (unsigned int k = 0; k != nz; ++k)
    {
        std::cout << k << " ";
    for (unsigned int j = 0; j != ny; ++j)
    for (unsigned int i = 0; i != nx; ++i)
    {
        uint16_t v = static_cast<uint16_t>(out.get(i,j,k));
        safe.set(i,j,k, v);
        volumes[v] += 1;
        xPos[v] += i;    
        yPos[v] += j;    
        zPos[v] += k;    
    }
    }
    std::cout << "finished, writing centroids" << std::endl;
    // write centroids
    std::string filename = outFolderName +  "/particles_centroids_ed" + std::to_string(static_cast<int>(erosiondepth)) +  ".dat";
    std::cout << "filename: " << filename << std::endl;
    std::ofstream myOut(filename);
    if(!myOut.good())
    {
        std::cerr << "error: cannot write" << std::endl;
        return 1;
    }
    myOut << "#1_x #2_y #3_z #4_l #5_vol" << std::endl;
    myOut << "#boxsz=" << nz << ", boundary_condition= none " << std::endl;
    for(size_t i = 0; i != volumes.size(); ++i)
    {
        double v = static_cast<double>(volumes[i]);
        if(v == 0) continue;
        myOut << " " << xPos[i]/v << " " << yPos[i]/v << " " << zPos[i]/v << " " << i << " "  << volumes[i] << std::endl;
    }
    myOut.close();


    // save images, write tomogram
    WriteImages(safe, outFolderName+"images/", "5_particles_woInsideHoles", nz); 

    filename = outFolderName +  "particles_labelled_ed" + std::to_string(static_cast<int>(erosiondepth)) +  ".raw";
    safe.write_raw(filename, "uint16_t");

    return 0;
}
