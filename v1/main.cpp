#include "Twig.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unordered_map>
#include "Twig3D.h"
#include <omp.h>
#include <ctime>

inline int Above(int i) { return i - 1; }
inline int Right(int j) { return j + 1; }

void Recurse(Twig& T, D2* alg)
{
    if(T.nopencells == 0 || T.ndeadcells == alg->i)
    {
        const int x = T.ncells - 1;
        const int y = T.ndeadcells;
        const int tid = omp_get_thread_num();
        alg->twigs_num[tid]++;
        if(alg->estimate != -1 && static_cast<double>(alg->twigs_num[tid] * alg->twigs_num.size()) >= alg->perc * alg->estimate) {
            cout << "Estimated at " << alg->perc << endl;
            alg->perc += 0.1;
        }
        alg->coefficients[tid][x][y]++;
        return;
    }
    
    Twig newT(T);
    #pragma omp task firstprivate(newT)
    {
        Retract r;
        pair<int,int> open_cell(newT.open_cells.front().first, newT.open_cells.front().second);
        if(newT.operation_star(1, 1, &r))
        {
            Recurse(newT, alg);
            newT.undo_operation_star(1, 1, open_cell, r);
        }
    }
    #pragma omp task firstprivate(newT)
    {
        Retract r;
        pair<int,int> open_cell(newT.open_cells.front().first, newT.open_cells.front().second);
        if(newT.operation_star(2, 2, &r))
        {
            Recurse(newT, alg);
            newT.undo_operation_star(2, 2, open_cell, r);
        }
    }
    #pragma omp task firstprivate(newT)
    {
        Retract r;
        pair<int,int> open_cell(newT.open_cells.front().first, newT.open_cells.front().second);
        if(newT.operation_star(3, 3, &r))
        {
            Recurse(newT, alg);
            newT.undo_operation_star(3, 3, open_cell, r);
        }
    }
    #pragma omp task firstprivate(newT)
    {
        Retract r;
        pair<int,int> open_cell(newT.open_cells.front().first, newT.open_cells.front().second);
        if(newT.operation_star(4, 2, &r))
        {
            Recurse(newT, alg);
            newT.undo_operation_star(4, 2, open_cell, r);
        }
    }
    #pragma omp task firstprivate(newT)
    {
        Retract r;
        pair<int,int> open_cell(newT.open_cells.front().first, newT.open_cells.front().second);
        if(newT.operation_star(5, 3, &r))
        {
            Recurse(newT, alg);
            newT.undo_operation_star(4, 2, open_cell, r);
        }
    }
}

void RecurseV2(Twig& T, D2* alg)
{
    if(T.nopencells == 0 || T.ndeadcells == alg->i)
    {
        const int x = T.ncells - 1;
        const int y = T.ndeadcells;
        const int tid = omp_get_thread_num();
        alg->twigs_num[tid]++;
        if(alg->estimate != -1 && static_cast<double>(alg->twigs_num[tid] * alg->twigs_num.size()) >= alg->perc * alg->estimate) {
            cout << "Estimated at " << alg->perc << endl;
            alg->perc += 0.1;
        }
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

void ParseResults(int dim, unordered_map<int, unordered_map<int, long long> >& coefficients, int i) {
    cout << "f:=(x,y)->x/(1";
    for(unordered_map<int, unordered_map<int, long long> >::iterator itx = coefficients.begin() ;
        itx != coefficients.end() ; itx++)
    {
        for(unordered_map<int, long long>::iterator ity = itx->second.begin() ; ity != itx->second.end() ; ity++)
        {
            cout << " -";
            
            switch(ity->second)
            {
                case 1:
                    break;
                default:
                    cout << ity->second << "*";
            }
            
            switch(itx->first)
            {
                case 0:
                    break;
                case 1:
                    cout << "x*";
                    break;
                default:
                    cout << "x^" << itx->first << "*";
            }
            
            switch(ity->first)
            {
                case 0:
                    break;
                case 1:
                    cout << "y";
                    break;
                default:
                    cout << "y^" << ity->first;
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
}

int main(int argc, char *argv[])
{
    int dimension = stoi(argv[1]);
    int i = stoi(argv[2]);
    long long prev = -1;
    if(argc > 3)
        prev = atoll(argv[3]);
    
    /*int tid,nthreads;
    char *cpu_name;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &tid);
    MPI_Comm_size(MPI_COMM_WORLD, &nthreads);*/
    
    int num_threads = 396;
    
    cout << "#CPU: " << num_threads << endl;
    
    clock_t begin = clock();

    vector<unordered_map<int, unordered_map<int, long long> > > coefficients2, coefficients3;
    vector<long long> twig_num2D, twig_num3D;
    // 2D
    if(dimension == 2) {
        D2 compute2DUpperBound(i, prev, num_threads);
        // the single twig
        Twig T(1,0,1);
        T.cells[ROOTY][ROOTX].type = OPEN;
        T.cells[ROOTY][ROOTX].context = LC1;
        T.open_cells.push_back(pair<int,int> (ROOTY, ROOTX));
        #pragma omp parallel shared(T, compute2DUpperBound)
        {
            #pragma omp single
            {
                Recurse(T, &compute2DUpperBound);
            }
        }
        /*compute2DUpperBound.fillLayer(num_threads);
#pragma omp parallel for shared (compute2DUpperBound)
        for(int i = 0; i < compute2DUpperBound.layer.size(); i++) {
            {
                RecurseV2(compute2DUpperBound.layer[i], &compute2DUpperBound);
            }
        }*/

        coefficients2 = compute2DUpperBound.coefficients;
        twig_num2D = compute2DUpperBound.twigs_num;
    }
    else { // 3D
        cout << "generating 3D twigs with " << i << " dead cells." << endl;
        D3 compute3DUpperBound(i, prev, num_threads);
        compute3DUpperBound.findAll3DContexts();
        compute3DUpperBound.createBasicTwigs();
        vector<openCell> open_cell;
        vector<Point3> Xs;
        open_cell.push_back(pair<Point3,int>(Point3(0,0,0),0));
        Twig3D T3(1, open_cell, Xs, 0);
        #pragma omp parallel shared(T3, compute3DUpperBound)
        {
            #pragma omp single
            {
                Recurse3D(T3, &compute3DUpperBound);
            }
        }
        coefficients3 = compute3DUpperBound.coefficients;
        twig_num3D = compute3DUpperBound.twigs_num;
    }
    
    vector<unordered_map<int, unordered_map<int, long long> > >& coefficients(dimension == 2 ? coefficients2 : coefficients3);
    vector<long long>& twig_num = (dimension == 2 ? twig_num2D : twig_num3D);
    
    unordered_map<int, unordered_map<int, long long> > summed_coefficients;
    long long twigs_num = 0;
    
    for(int n = 0; n < num_threads; n++)
    {
        twigs_num += twig_num[n];
        for(unordered_map<int, unordered_map<int, long long> >::iterator itx = coefficients[n].begin() ; itx != coefficients[n].end() ; itx++)
        {
            for(unordered_map<int, long long>::iterator ity = itx->second.begin() ; ity != itx->second.end() ; ity++) {
                summed_coefficients[itx->first][ity->first] += ity->second;
            }
        }
    }
    
    cout << "Generated " << twigs_num << " twigs." << endl;
    ParseResults(dimension, summed_coefficients, i);
    
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    
    cout << "Total time: " << elapsed_secs/60.0 << " minutes (" << elapsed_secs/(60.0*60.0) << " hours)" << endl;
    
    return 0;
}
