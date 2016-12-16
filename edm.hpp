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
#ifndef EDM_H_INCLUDEGUARD_123456
#define EDM_H_INCLUDEGUARD_123456

#include <assert.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>

#include "tomogram.hpp"


class edm
{
public:
    template <class T>
    inline static T square(const T &x)
    {
        return x*x;
    };

// fabian stuff -,-
#define INF 1E20

    static float *dt(float *f, int n)
    {
        float *d = new float[n];
        int *v = new int[n];
        float *z = new float[n+1];
        int k = 0;
        v[0] = 0;
        z[0] = -INF;
        z[1] = +INF;
        for (int q = 1; q <= n-1; q++)
        {
            float s  = ((f[q]+square(q))-(f[v[k]]+square(v[k])))/(2*q-2*v[k]);
            while (s <= z[k])
            {
                k--;
                s  = ((f[q]+square(q))-(f[v[k]]+square(v[k])))/(2*q-2*v[k]);
            }
            k++;
            v[k] = q;
            z[k] = s;
            z[k+1] = +INF;
        }

        k = 0;
        for (int q = 0; q <= n-1; q++)
        {
            while (z[k+1] < q)
                k++;
            d[q] = square(q-v[k]) + f[v[k]];
        }

        delete [] v;
        delete [] z;
        return d;
    }
    template <typename VALUE>
    static void create_edm(tomogram3d<VALUE>& matrix_in, tomogram3d<float>& edm)
    {

        std::cout << "calculating edm ... " << std::flush;

        int sx = matrix_in.get_sx();
        int sy = matrix_in.get_sy();
        int sz = matrix_in.get_sz();

        edm.initialise(sx,sy,sz,-1);

        for (int z = 0; z < sz; z++)
            for (int y = 0; y < sy; y++)
                for (int x = 0; x < sx; x++)
                    edm.set(x,y,z,float(matrix_in.get(x,y,z))*INF);


        float *f = new float[std::max(sx,sz)]; //fixme sy, not quadratic case

        for (int z = 0; z < sz; z++)
        {
            for (int y = 0; y < sy; y++)
            {
                for (int x = 0; x < sx; x++)
                {
                    f[x] = edm.get(x, y, z);
                }
                float *d = dt(f, sx);
                for (int x = 0; x < sx; x++)
                {
                    edm.set (x, y, z, d[x]);
                }
                delete [] d;
            }
        }

        for (int y = 0; y < sy; y++)
        {
            for (int x = 0; x < sx; x++)
            {
                for (int z = 0; z < sz; z++)
                {
                    f[z] = edm.get( x, y, z);
                }
                float *d = dt(f, sz);
                for (int z = 0; z < sz; z++)
                {
                    edm.set (x, y, z, d[z]);
                }
                delete [] d;
            }

        }

        for (int x = 0; x < sx; x++)
        {
            for (int z = 0; z < sz; z++)
            {
                for (int y = 0; y < sy; y++)
                {
                    f[y] = edm.get( x, y, z);
                }
                float *d = dt(f, sy);
                for (int y = 0; y < sy; y++)
                {
                    edm.set (x, y, z, d[y]);
                }
                delete [] d;
            }

        }
        delete []f;
        std::cout << " done" << std::endl;
    };
};

#endif // EDM
