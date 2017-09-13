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
#ifndef TOMOGRAM
#define TOMOGRAM

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <stdint.h>
#include <set>
#include <random>
#include "colors.hpp"

template <typename VALUE>
class tomogram3d
{
public:
    void initialise(unsigned long a,unsigned long b,unsigned long c, VALUE value)
    {
        sx = a;
        sy = b;
        sz = c;
        std::cout << "resizing tomogram to " << a << " "  << b << " "  << c << std::endl;
        data.resize(sx*sy*sz);
        for (unsigned long i = 0; i < sx*sy*sz; ++i)
        {
            data.at(i) = value;
        }
        std::cout << "resize done" << std::endl;
    }

    void clear()
    {
        data.clear();
        sx = 0;
        sy = 0;
        sz = 0;
    }

    unsigned long get_sx() const
    {
        return sx;
    }
    unsigned long get_sy() const
    {
        return sy;
    }
    unsigned long get_sz() const
    {
        return sz;
    }


    VALUE get (unsigned long x, unsigned long y, unsigned long z) const
    {
        unsigned long pos = (x + y*sx + z*sx*sy);
        if (pos > data.size()) throw std::string("GET: index exceeds image bounds " + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + " " + std::to_string(pos));
        return data[pos];
    }

    void set (unsigned long x, unsigned long y, unsigned long z, VALUE value)
    {
        unsigned long pos = (x + y*sx + z*sx*sy);
        if (pos > data.size()) throw std::string(" SET: index exceeds image bounds"+ std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + " " + std::to_string(pos));
        data[pos] = value;
    }


    VALUE find_min() const
    {
        VALUE minvalue = *std::min_element (data.begin(),data.end());
        return minvalue;
    }

    VALUE find_max() const
    {
        VALUE maxvalue = *std::max_element (data.begin(),data.end());
        return maxvalue;
    }

    void write_raw (std::string outfilename, std::string type) const
    {
        std::cout << "write " << outfilename << " ... " << std::flush;
        std::ofstream outfile;
        outfile.open(outfilename.c_str());

        if (type == "uint8_t")
        {
            for (unsigned long i = 0; i<sx*sy*sz; i++)
            {
                uint8_t number;
                if (data.at(i) < 0 || data.at(i) > std::numeric_limits<uint8_t>::max())
                    throw std::runtime_error ("Cant't write file: Numbers are to big for selected class type");
                number = uint8_t(data.at(i));
                outfile.write ((char *)&number,sizeof (number));
            }
        }

        else if (type == "uint16_t")
        {
            for (unsigned long i = 0; i<sx*sy*sz; i++)
            {
                uint16_t number;
                if (data.at(i) < 0 || data.at(i) > std::numeric_limits<uint16_t>::max())
                    throw std::runtime_error ("Cant't write file: Numbers are to big for selected class type");
                number = uint16_t(data.at(i));
                outfile.write ((char *)&number,sizeof (number));
            }
        }

        else if (type == "uint32_t")
        {
            for (unsigned long i = 0; i<sx*sy*sz; i++)
            {
                uint32_t number;
                if (data.at(i) < 0 || data.at(i) > static_cast<uint32_t>(std::numeric_limits<uint32_t>::max()))
                    throw std::runtime_error ("Cant't write file: Numbers are to big for selected class type");
                number = uint32_t(data.at(i));
                outfile.write ((char *)&number,sizeof (number));
            }
        }

        else if (type == "int")
        {
            for (unsigned long i = 0; i<sx*sy*sz; i++)
            {
                int number;
                if (data.at(i) > std::numeric_limits<int>::max())
                    throw std::runtime_error ("Cant't write file: Numbers are to big for selected class type");
                number = int(data.at(i));
                outfile.write ((char *)&number,sizeof (number));
            }
        }
        else
        {
            throw std::runtime_error ("can't write raw type");
        }

        outfile.close();
        std::cout << " done" << std::endl;
    }

    void write_pgm (std::string outfilename, unsigned long picnumber) const
    {
        std::cout << "write " << outfilename << " ... " << std::flush;
        std::ofstream myfile;
        myfile.open (outfilename.c_str());
        myfile << "P2" << std::endl;
        myfile << sx << " " << sy << std::endl;

        VALUE maxvalue = *std::max_element (data.begin()+(picnumber*sx*sy),data.begin()+((picnumber+1)*sx*sy));

        myfile << int(maxvalue) << std::endl;
        for (unsigned long j = 0; j != sy; ++j)
            for (unsigned long i = 0; i != sx; ++i)
            {
                myfile << int(get(i,j,picnumber)) << std::endl;
            }
        myfile.close();
        std::cout << " done" << std::endl;
    }

    void write_ppm(std::string outfilename, unsigned long picnumber) const
    {
        std::cout << "write " << outfilename << " ... " << std::flush;
        std::ofstream myfile;
        myfile.open (outfilename.c_str());
        myfile << "P3" << std::endl;
        myfile << sx << " " << sy << "\n255"  << std::endl;
        VALUE maxvalue = *std::max_element (data.begin()+(picnumber*sx*sy),data.begin()+((picnumber+1)*sx*sy));
        std::vector<int> cvals((maxvalue+1)*3, 0);

        std::default_random_engine rng;
        std::uniform_real_distribution<double> d(0,1);

        unsigned int N = 15;
        // create color list
        std::cout << "create " << maxvalue << " colors" << std::endl;


        unsigned int idx = 0;
        for (unsigned long k = 1; k <= maxvalue; ++k)
        {
            hsv ch;
            ch.h = 360.0/static_cast<double>(N) * static_cast<double>(idx);
            ch.s = 0.5 + 0.5 * d(rng);
            ch.v = 0.3 + 0.7 * d(rng);
            rgb c = hsv2rgb(ch);
            idx++;
            if (idx >= N) idx = 0;
            //std::cout << "r= " << c.r*255.0 << std::endl;
            //std::cout << "g= " << c.g*255.0 << std::endl;
            //std::cout << "b= " << c.b << std::endl;
            cvals[k*3 + 0] = (c.r*255.0);
            cvals[k*3 + 1] = (c.g*255.0);
            cvals[k*3 + 2] = (c.b*255.0);
        }

        for (unsigned long j = 0; j != sy; ++j)
        {
            for (unsigned long i = 0; i != sx; ++i)
            {
                int val = int(get(i,j,picnumber));

                int r = cvals.at(val*3 + 0);
                int g = cvals.at(val*3 + 1);
                int b = cvals.at(val*3 + 2);
                myfile << r << " " << g << " " << b << " ";
            }
            myfile << std::endl;
        }
        myfile.close();
        std::cout << " done" << std::endl;
    }

    void read_raw (std::string infilename, std::string type, unsigned long a, unsigned long b, unsigned long c)
    {
        sx = a;
        sy = b;
        sz = c;
        data.resize(sx*sy*sz);
        std::cout << "size: " << a*b*c << "; arraysize: " << data.size() << std::endl;
        std::ifstream infile;
        infile.open(infilename.c_str());
        if(infile.fail())
        {
            throw std::string("cannot open raw file");
        }
        if (type == "uint8_t")
        {
            for (unsigned long i = 0; i<sx*sy*sz; i++)
            {
                uint8_t number;
                infile.read ((char *)&number,sizeof (number));
                if (number > std::numeric_limits<VALUE>::max())
                    throw std::runtime_error ("Can't read file: Numbers are to big for selected class type");
                data.at(i) = VALUE(number);
            }
        }
        else if (type == "uint16_t")
        {
            for (unsigned long i = 0; i<sx*sy*sz; i++)
            {
                uint16_t number;
                infile.read ((char *)&number,sizeof (number));
                if (number > std::numeric_limits<VALUE>::max())
                    throw std::runtime_error ("Can't read file: Numbers are to big for selected class type");
                data.at(i) = VALUE(number);
            }
        }
        else if (type == "int32_t")
        {
            for (unsigned long i = 0; i<sx*sy*sz; i++)
            {
                int32_t number;
                infile.read ((char *)&number,sizeof (number));
                if (number > std::numeric_limits<VALUE>::max())
                    throw std::runtime_error ("Can't read file: Numbers are to big for selected class type");
                data.at(i) = VALUE(number);
            }
        }
        else if (type == "int")
        {
            for (unsigned long i = 0; i<sx*sy*sz; i++)
            {
                int number;
                infile.read ((char *)&number,sizeof (number));
                if (number > std::numeric_limits<VALUE>::max())
                    throw std::runtime_error ("Can't read file: Numbers are to big for selected class type");
                data.at(i) = VALUE(number);
            }
        }
        else
        {
            throw std::runtime_error ("can't read raw type");
        }

        infile.close();
    }

    void read_bin (std::string infilename)
    {
        std::cout << "reading " << infilename << " ... " << std::flush;
        std::ifstream ifile;
        ifile.open(infilename.c_str());
        ifile >> sx >> sy >> sz;
        data.resize(sx*sy*sz);

        for (unsigned int i = 0; i < sx*sy*sz; i++)
        {
            for(;;)
            {
                int c = ifile.get();
                switch (c)
                {
                case EOF:
                    std::cerr << "end of file while reading structure";
                    abort ();
                case '\n':
                case ' ':
                case '\t':
                    // ignore
                    continue;
                default:
                    int x = c - '0';
                    if (x < 0 || x > 10)
                    {
                        std::cerr << "invalid character in file: "
                                  << (char)(c) << "\n";
                        abort ();
                    }
                    data.at(i)= (VALUE) x;
                }
                break;
            }
        }
        ifile.close();
        std::cout << " done" << std::endl;
    }
    
    void cleanBorder(double outside = 0.005)
    {
        unsigned nx = get_sx();
        unsigned ny = get_sy();
        unsigned nz = get_sz();
        
        std::set<VALUE> markList; // all particle labels contained in here will be set to 0
        double zExcludeFraction = 0.05;
        unsigned zTop = zExcludeFraction * nz;
        unsigned zBot = nz - zTop;
    
        for (unsigned long k = 0; k != nz; ++k)
        for (unsigned long j = 0; j != ny; ++j)
        for (unsigned long i = 0; i != nx; ++i)
        {
            // exclude all particles that touch the boundary
            
            auto v = get(i,j,k);
            
            if( v != 0)
            {
                // TOP && bottom
                if(k < zTop || k > zBot)
                {
                    markList.insert(v);
                }
                double r = sqrt( (i-nx/2)*(i-nx/2) + (j- ny/2)*(j-ny/2));
                // radial
                if (r > 0.5*nx*(1.0 - outside))
                {
                    markList.insert(v);
                }
            }
        }
        std::cout << "CleanUp: removing " << markList.size() << " clusters" << std::endl;
// remove all clusters that touch the boundary
        for (unsigned int k = 0; k != nz; ++k)
        for (unsigned int j = 0; j != ny; ++j)
        for (unsigned int i = 0; i != nx; ++i)
        {
            auto v = get(i,j,k);
            if( markList.count(v) == 1)
            {
                set(i,j,k,0);
            }
        }
    }


private:
    std::vector <VALUE> data;
    unsigned long sx,sy,sz;
};



#endif // TOMOGRAM
