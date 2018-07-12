#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <ctime>
#include "Twig3D.h"

#define MAX_CORES 396

inline int Above(int i) { return i - 1; }
inline int Right(int j) { return j + 1; }

void Recurse(Twig3D& T, D3* alg)
{
    if(T.nopencells == 0 || T.ndeadcells == alg->i)
    {
        const int x = T.ncells - 1;
        const int y = T.ndeadcells;
        int tid=0;
        //const int tid = omp_get_thread_num();
        alg->twigs_num[tid]++;
        alg->coefficients[tid][x][y]++;
        return;
    }
    Retract3D r;
    for(int i = 0; i < alg->basic_twigs.size(); i++) {
        if(operation_star(*alg, &T, alg->basic_twigs[i], &r))
            Recurse(T, alg);
        undo_operation_star(&T, &r);
    }
}

void ParseResults(vector<vector<long long> >& coefficients, int i, long long run) {
    FILE* fp;
    string filename = "run" + std::to_string(run) + ".bin";
    fp = fopen (filename.c_str(), "wb");
    
    cout << "f:=(x,y)->x/(1";
    for(int x = 0; x < 4*i+1; x++)
    {
        for(int y = 0; y < i+1; y++)
        {
            fwrite(&x, sizeof(int), 1, fp);
            fwrite(&y, sizeof(int), 1, fp);
            fwrite(&coefficients[x][y], sizeof(long long), 1, fp);
            
            if(coefficients[x][y] == 0)
                continue;
            cout << " -";
            
            switch(coefficients[x][y])
            {
                case 1:
                    break;
                default:
                    cout << coefficients[x][y] << "*";
            }
            
            switch(x)
            {
                case 0:
                    break;
                case 1:
                    cout << "x*";
                    break;
                default:
                    cout << "x^" << x << "*";
            }
            
            switch(y)
            {
                case 0:
                    break;
                case 1:
                    cout << "y";
                    break;
                default:
                    cout << "y^" << y;
            }
        }
    }
    cout << ");\n";
    cout << "func:=(x,y)->f(x,y)/x;" << endl;
    cout << "sol := func(s,z/s);" << endl;
    cout << "d := denom(sol);" << endl;
    cout << "dis := discrim(dd,s);" << endl;
    cout << "ro := fsolve(dis=0,z,maxsols=10000);" << endl;
    cout << "r:=evalf(1/ro[1]);" << endl;
    cout << "save r, ro, \"ub" << unsigned(i) << ".m\"" << endl;
    fclose(fp);
}

int main(int argc, char *argv[])
{
    int i = stoi(argv[1]);
    // Default is no division and a single run
    const int divide = argc > 2 ? atoi(argv[2]) : 1;
    const int run = argc > 3 ? atoi(argv[3]) : 0;
    
    cout << "generating 3D twigs with " << i << " dead cells." << endl;
    
    vector <vector< vector< long long > > > coefficients3;
    vector<long long> twig_num3D;
    
    // Preprocessing phase
    D3 compute3DUpperBound(i, -1, 1);
    compute3DUpperBound.findAll3DContexts();
    compute3DUpperBound.createBasicTwigs();
    
    
        fillLayer(1, 0, run, &compute3DUpperBound);
        cout << "layer = " << compute3DUpperBound.layer.size() << endl;
        if(compute3DUpperBound.layer.size() > 0) {
            Recurse(compute3DUpperBound.layer[0], &compute3DUpperBound);
        }
    
    coefficients3 = compute3DUpperBound.coefficients;
    twig_num3D = compute3DUpperBound.twigs_num;
    
    cout << "done with twig generation" << endl;
    
    vector<vector<long long> > summed_coefficients(4*i+1);
    for(int j = 0; j < 4*i+1; j++) {
        summed_coefficients[j].resize(i+1);
        for(int k = 0; k < i+1; k++)
            summed_coefficients[j][k] = 0;
    }
    
    long long twigs_num = 0;
    
    for(int n = 0; n < 1; n++)
    {
        twigs_num += twig_num3D[n];
        for(int x = 0; x < 4*i+1; x++)
        {
            for(int y = 0; y < i+1; y++) {
                summed_coefficients[x][y] += coefficients3[n][x][y];
            }
        }
    }
    
    return 0;
}
