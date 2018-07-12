#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unordered_map>
#include <omp.h>
#include <algorithm>
#include <ctime>
#include "Twig.h"

// Each node has 12 cores and we may use up to 33 nodes.
#define MAX_CORES 396

inline int Above(int i) { return i - 1; }
inline int Right(int j) { return j + 1; }

void RecurseV2(Twig& T, D2* alg)
{
    if(T.nopencells == 0 || T.ndeadcells == alg->i)
    {
        const int x = T.ncells - 1;
        const int y = T.ndeadcells;
        const int tid = omp_get_thread_num();
        alg->twigs_num[tid]++;
        alg->coefficients[tid][x][y]++;
        return;
    }
    
    Retract r;
    pair<int,int> open_cell(T.open_cells.front().first, T.open_cells.front().second);
    if(T.operation_star(1, 1, &r))
    {
        RecurseV2(T, alg);
        T.undo_operation_star(1, 1, open_cell, r);
    }
    else
        T.undo_operation_star(1, 1, open_cell, r, false);
    
    if(T.operation_star(2, 2, &r))
    {
        RecurseV2(T, alg);
        T.undo_operation_star(2, 2, open_cell, r);
    }
    else
        T.undo_operation_star(2, 2, open_cell, r, false);
    
    if(T.operation_star(3, 3, &r))
    {
        RecurseV2(T, alg);
        T.undo_operation_star(3, 3, open_cell, r);
    }
    else
        T.undo_operation_star(3, 3, open_cell, r, false);
    
    if(T.operation_star(4, 2, &r))
    {
        RecurseV2(T, alg);
        T.undo_operation_star(4, 2, open_cell, r);
    }
    else
        T.undo_operation_star(4, 2, open_cell, r, false);
    
    if(T.operation_star(5, 3, &r))
    {
        RecurseV2(T, alg);
        T.undo_operation_star(5, 3, open_cell, r);
    }
    else
        T.undo_operation_star(5, 3, open_cell, r, false);
}

void ParseResults(vector<vector<long long> >& coefficients, int i, long long run) {
    FILE* fp;
    string filename = "run" + std::to_string(run) + ".bin";
    fp = fopen (filename.c_str(), "wb");
    
    cout << "f:=(x,y)->x/(1";
    for(int x = 0; x < 2*i+1; x++)
    {
        for(int y = 0; y < i+1; y++)
        {
            // Write the coefficients from the current run to disk
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
    long long prev = -1;
    int divide;
    long long run = 0;
    if(argc > 2)
        divide = atoi(argv[2]);
    if(argc > 3)
        run = atoll(argv[3]);
    
    const int max_range = MAX_CORES / divide;
    const int layer_size = MAX_CORES - 2;
    
    int tid,nthreads = 1;
    char *cpu_name;
    double time_initial,time_current,time;
    
    // MPI code
    MPI_Init(&argc,&argv);
    time_initial  = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &tid);
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);
    cpu_name = (char *)calloc(80,sizeof(char));
    gethostname(cpu_name,80);
    time_current  = MPI_Wtime();
    time  = time_current - time_initial;
    
    int num_cores;
    
#pragma omp parallel
    {
#pragma omp critical
        num_cores++;
    }
    
    vector <vector< vector< long long > > > coefficients2;
    vector<long long> twig_num2D;
    D2* compute2DUpperBound = new D2(i, prev, num_cores);
    if(tid + run * max_range < layer_size && tid + run * max_range < (run+1) * max_range) {
        compute2DUpperBound->fillLayer(layer_size+2, tid, run);
        if(compute2DUpperBound->layer.size() > 0) {
            RecurseV2(compute2DUpperBound->layer[tid + run * max_range], compute2DUpperBound);
        }
    }
    coefficients2 = compute2DUpperBound->coefficients;
    twig_num2D = compute2DUpperBound->twigs_num;
    
    vector<vector<long long> > summed_coefficients(2*i+1);
    for(int j = 0; j < 2*i+1; j++) {
        summed_coefficients[j].resize(i+1);
        for(int k = 0; k < i+1; k++)
            summed_coefficients[j][k] = 0;
    }
    
    // Gather information from all cores
    
    long long twigs_num = 0;
    
    for(int n = 0; n < num_cores; n++)
    {
        twigs_num += twig_num2D[n];
        for(int x = 0; x < 2*i+1; x++)
        {
            for(int y = 0; y < i+1; y++) {
                summed_coefficients[x][y] += coefficients2[n][x][y];
            }
        }
    }
    
    long long total_twigs = 0;
    MPI_Reduce(&twigs_num, &total_twigs, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    vector<vector<long long> > total_coefficients(2*i+1);
    for(int j = 0; j < 2*i+1; j++) {
        total_coefficients[j].resize(i+1);
        for(int k = 0; k < i+1; k++)
            total_coefficients[j][k] = 0;
    }
    
    for(int x = 0; x < 2*i+1; x++)
    {
        for(int y = 0; y < i+1; y++)
        {
            MPI_Reduce(&summed_coefficients[x][y], &total_coefficients[x][y], 1,
                       MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    }
    
    if(tid == 0) {
        cout << "Generated " << total_twigs << " twigs." << endl;
        ParseResults(total_coefficients, i, run);
    }
    
    MPI_Finalize();
    return 0;
}
