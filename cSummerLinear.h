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
#ifndef CSUMMER_H_INCLUDED
#define CSUMMER_H_INCLUDED

#include<map>
#include <iostream>
#include <fstream>
#include <string>
class cSummerLinear 
{
public:
    cSummerLinear (double BinSize):m_dBinSize(BinSize) 
    {
        if(m_dBinSize < 0)
        {
            throw std::string ("negative BinSize is not supported in ISummer.");
        }
    };
    
    ~cSummerLinear()
    {
        m_mapHistogram.clear();
    };

    void Add(double Radius)
    {
        if(Radius < 0.0)
        {
            throw std::string ("negative Keys are not supported by cSummerLinear.");
        }
        int idx = static_cast<int>(Radius/m_dBinSize);
        //std::cout << idx << " ";
        if ( !m_mapHistogram.count(idx) )
        {
            m_mapHistogram.insert(std::pair<unsigned int,double>(idx,0.0));
        }
        m_mapHistogram.find(idx)->second +=1.0;
    };


    void Add (double Key, double value)
    {
        if(Key < 0.0)
        {
            throw std::string ("negative Keys are not supported by cSummerLinear.");
        }
        unsigned int idx = static_cast<unsigned int>(Key/m_dBinSize);

        if ( !m_mapHistogram.count(idx) )
        {
            m_mapHistogram.insert(std::pair<unsigned int,double>(idx,0.0));
        }
        m_mapHistogram.find(idx)->second += value;

    };

    friend std::ostream& operator << (std::ostream &f, const cSummerLinear& S)
    {
        if (!S.m_mapHistogram.empty())
        {
            std::cout << "saving histogram up to bin " << (--S.m_mapHistogram.end())->first << std::endl;
            for ( unsigned long i = 0; i < (--S.m_mapHistogram.end())->first; ++i)
            {
                f << static_cast<double>(S.m_dBinSize * i);
                f << "    ";
                if ( S.m_mapHistogram.count(i) )
                {
                    f << S.m_mapHistogram.find(i)->second;
                }
                else
                {
                    f << 0;
                }

                f << std::endl;
            }

        }
        return f;
    };

    // normalizes all Values to a max value of 1
    void NormalizePeak ()
    {
        // first: determine max value in Summer
        double t_dMax = 0;
        std::map<unsigned long,double>::const_iterator cit;
        for (cit = m_mapHistogram.begin(); cit != m_mapHistogram.end(); cit++)
        {
            if ( cit->second > t_dMax )
            {
                t_dMax = cit->second;
            }
        }


        std::map<unsigned long,double>::iterator it;
        for ( it = m_mapHistogram.begin(); it != m_mapHistogram.end(); ++it)
        {
            it->second /= t_dMax;
        }
        std::cout << m_mapHistogram.size() << " Values normalized with max of " << t_dMax << std::endl;

    };

    unsigned int GetSize()
    {
     return m_mapHistogram.size();
    } ;

    std::map<unsigned long, double>& GetMap() {return m_mapHistogram;};

    void clear ()
    {
        m_mapHistogram.clear();
    }


    void clean ()
    {

        std::map<unsigned long,double>::iterator it;
        for ( it = m_mapHistogram.begin(); it != m_mapHistogram.end(); ++it)
        {
            if( (*it).second <= 0)
            {
                m_mapHistogram.erase(it);
                it = m_mapHistogram.begin();
            }
        }
    }

private:
    // this map stores pairs of double and int where the first int is the indexed square radius via m_dBinSize
    // and the latter int is the number of counted incidents
    std::map<unsigned long,double> m_mapHistogram;

    double m_dBinSize;

};


#endif // CSUMMER_H_INCLUDED
