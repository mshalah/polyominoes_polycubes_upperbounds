#ifndef mira_twigs_Twig_h
#define mira_twigs_Twig_h

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <fstream>
#include <cstring>
#include <unordered_map>
#include <deque>
#include <omp.h>

#define ROOTX 24
#define ROOTY 24
#define MAX_SIZE 2*ROOTX

using namespace std;

// All L contexts of a cell (page 594)
typedef enum {
    LC1, LC2, LC3, LC4, LC5, LC6, LC7, LC8
} Lcontext;

// types of a cell
typedef enum {
    FREE, DEAD, OPEN, X
} Type;

struct cell
{
    cell() : type(FREE) {};
    Type type;
    Lcontext context;
};

struct Retract {
    pair<int,int> Rcell, Ucell, Ccell;
    bool RwasX, UwasX, CwasX;
    Retract() : RwasX(true), UwasX(true), CwasX(true) {};
};

Lcontext rotate90(Lcontext L)
{
    switch(L)
    {
        case LC1:
            return LC3;
            break;
        case LC2:
            return LC4;
            break;
        case LC3:
            return LC2;
            break;
        case LC4:
            return LC1;
            break;
        case LC5:
            return LC7;
            break;
        case LC6:
            return LC8;
            break;
        case LC7:
            return LC6;
            break;
        case LC8:
            return LC5;
            break;
    }
    return L;
}

Lcontext rotate180(Lcontext L)
{
    switch(L)
    {
        case LC1:
            return LC2;
            break;
        case LC2:
            return LC1;
            break;
        case LC3:
            return LC4;
            break;
        case LC4:
            return LC3;
            break;
        case LC5:
            return LC6;
            break;
        case LC6:
            return LC5;
            break;
        case LC7:
            return LC8;
            break;
        case LC8:
            return LC7;
            break;
    }
    return L;
}

Lcontext rotate270(Lcontext L)
{
    switch(L)
    {
        case LC1:
            return LC4;
            break;
        case LC2:
            return LC3;
            break;
        case LC3:
            return LC1;
            break;
        case LC4:
            return LC2;
            break;
        case LC5:
            return LC8;
            break;
        case LC6:
            return LC7;
            break;
        case LC7:
            return LC5;
            break;
        case LC8:
            return LC6;
            break;
    }
    return L;
}


Lcontext reflect(Lcontext L, int reflection)
{
    switch(reflection)
    {
        case 1:
        case 2: // reflect around x
            switch(L)
        {
            case LC1:
                return LC6;
                break;
            case LC2:
                return LC5;
                break;
            case LC3:
                return LC7;
                break;
            case LC4:
                return LC8;
                break;
            case LC5:
                return LC2;
                break;
            case LC6:
                return LC1;
            case LC7:
                return LC3;
                break;
            case LC8:
                return LC4;
                break;
        }
            break;
        case 3:
        case 4: // reflect around y
            switch(L)
        {
            case LC1:
                return LC5;
                break;
            case LC2:
                return LC6;
                break;
            case LC3:
                return LC8;
                break;
            case LC4:
                return LC7;
                break;
            case LC5:
                return LC1;
                break;
            case LC6:
                return LC2;
            case LC7:
                return LC4;
                break;
            case LC8:
                return LC3;
                break;
        }
            break;
    }
    return L;
}

// rotate and flip an L context
Lcontext rotate_and_flip(Lcontext L, int rot, int ref)
{
    Lcontext ret(L);
    switch(rot)
    {
        case 0:
            ret = L;
            break;
        case 90:
            ret = rotate90(L);
            break;
        case 180:
            ret = rotate180(L);
            break;
        case 270:
            ret = rotate270(L);
            break;
    }
    ret = reflect(ret,ref);
    return ret;
}

class Twig{
    
public:
    
    // number of cells
    int ncells;
    // number of dead cells
    int ndeadcells;
    // number of open cells
    int nopencells;
    // weight of a    twig
    //std::string weight;
    // cords of a root cell
    // pair<int,int> root;
    // table
    vector<vector<cell> > cells;
    // linearly ordered open cells
    deque<pair<int,int> >  open_cells;
    // Twig number (1-5) (if the twing is in set L(Figure 7 page 594))
    //int Ltype;
    //vector<int> seq;
    
    Twig(){};
    
    // constructor
    Twig(int n, int d, int o) : /*root(pair<int,int>(ROOTX,ROOTY)),*/ cells(MAX_SIZE)
    {
        ncells = n;
        ndeadcells = d;
        nopencells = o;
        cell free_cell;
        for(int i = 0; i < MAX_SIZE; i++) {
            for(int j = 0; j < MAX_SIZE; j++)
                cells[i].push_back(free_cell);
        }
    }
    
    // copy constructor
    Twig(const Twig& T) : ncells(T.ncells), ndeadcells(T.ndeadcells), nopencells(T.nopencells) // ,
    // root(T.root)
    {
        cells = T.cells;
        open_cells = T.open_cells;
    }
    
    bool Done(int i) {
        if(nopencells == 0 || ndeadcells == i)
            return true;
        else
            return false;
    }
    
    bool operation_star(int Ltype, int Lsize, Retract* r)
    {
        // get first open cell according to the linear order
        int opencell1i = open_cells.front().first,
        opencell1j = open_cells.front().second;
        // open cell is now dead
        open_cells.pop_front();
        cells[opencell1i][opencell1j].type = DEAD;
        nopencells--;
        ndeadcells++;
        ncells += (Lsize -1);
        
        int rightcelli = 0, rightcellj = 0, uppercelli = 0, uppercellj = 0, corneri = 0, cornerj = 0, rot = 0, ref = 0;
        
        // L context of the open cell
        switch(cells[opencell1i][opencell1j].context)
        {
                // coordinates of the cell to the right, above and corner according to the L context
            case LC1:
                rightcelli = 0;
                rightcellj = 1;
                uppercelli = -1;
                uppercellj = 0;
                corneri = -1;
                cornerj = 1;
                break;
            case LC2:
                rightcelli = 0;
                rightcellj = -1;
                uppercelli = 1;
                uppercellj = 0;
                corneri = 1;
                cornerj = -1;
                rot = 180;
                break;
            case LC3:
                rightcelli = 1;
                rightcellj = 0;
                uppercelli = 0;
                uppercellj = 1;
                corneri = 1;
                cornerj = 1;
                rot = 90;
                break;
            case LC4:
                rightcelli = -1;
                rightcellj = 0;
                uppercelli = 0;
                uppercellj = -1;
                corneri = -1;
                cornerj = -1;
                rot = 270;
                break;
            case LC5:
                rightcelli = 0;
                rightcellj = -1;
                uppercelli = -1;
                uppercellj = 0;
                corneri = -1;
                cornerj = -1;
                ref = 3;
                break;
            case LC6:
                rightcelli = 0;
                rightcellj = 1;
                uppercelli = 1;
                uppercellj = 0;
                corneri = 1;
                cornerj = 1;
                ref = 1;
                break;
            case LC7:
                rightcelli = -1;
                rightcellj = 0;
                uppercelli = 0;
                uppercellj = 1;
                corneri = -1;
                cornerj = 1;
                rot = 90;
                ref = 2;
                break;
            case LC8:
                rightcelli = 1;
                rightcellj = 0;
                uppercelli = 0;
                uppercellj = -1;
                corneri = 1;
                cornerj = -1;
                rot = 90;
                ref = 3;
                break;
        }
        
        int RCi = opencell1i+rightcelli, RCj =  opencell1j+rightcellj, UCi = opencell1i+uppercelli, UCj = opencell1j+uppercellj,
        CCi = opencell1i+corneri, CCj = opencell1j+cornerj;
        
        r->Rcell = pair<int,int> (RCi, RCj);
        r->Ucell = pair<int,int> (UCi, UCj);
        r->Ccell = pair<int,int> (CCi, CCj);
        
        if(RCi >= MAX_SIZE || RCj >= MAX_SIZE || UCi >= MAX_SIZE || UCj >= MAX_SIZE || CCi >= MAX_SIZE || CCj >= MAX_SIZE
           || RCi < 0 || RCj < 0 || UCi < 0 || UCj < 0 || CCi < 0 || CCj < 0) {
            printf("leak!\n");
            exit(1);
        }
        
        switch(Ltype)
        {
            case 1:
                // condition (*)
                if(!(cells[RCi][RCj].type == DEAD || cells[RCi][RCj].type == OPEN)) {
                    if(cells[RCi][RCj].type == FREE)
                        r->RwasX = false;
                    cells[RCi][RCj].type = X;
                }
                if(!(cells[UCi][UCj].type == DEAD || cells[UCi][UCj].type == OPEN)) {
                    if(cells[UCi][UCj].type == FREE)
                        r->UwasX = false;
                    cells[UCi][UCj].type = X;
                }
                break;
                
            case 2:
                if(cells[UCi][UCj].type == DEAD || cells[UCi][UCj].type == OPEN || cells[UCi][UCj].type == X)
                    return false;
                
                if(!(cells[RCi][RCj].type == DEAD || cells[RCi][RCj].type == OPEN)) {
                    if(cells[RCi][RCj].type == FREE)
                        r->RwasX = false;
                    cells[RCi][RCj].type = X;
                }
                if(!(cells[CCi][CCj].type == DEAD || cells[CCi][CCj].type == OPEN)) {
                    if(cells[CCi][CCj].type == FREE)
                        r->CwasX = false;
                    cells[CCi][CCj].type = X;
                }
                // new open cell
                cells[UCi][UCj].type = OPEN;
                cells[UCi][UCj].context = rotate_and_flip(LC5,rot,ref);
                open_cells.push_back(pair<int,int> (UCi, UCj));
                nopencells++;
                break;
                
            case 3:
                if(cells[UCi][UCj].type == DEAD || cells[UCi][UCj].type == OPEN || cells[UCi][UCj].type == X
                   || cells[CCi][CCj].type == DEAD || cells[CCi][CCj].type == OPEN || cells[CCi][CCj].type == X)
                    return false;
                if(!(cells[RCi][RCj].type == DEAD || cells[RCi][RCj].type == OPEN)) {
                    if(cells[RCi][RCj].type == FREE)
                        r->RwasX = false;
                    cells[RCi][RCj].type = X;
                }
                
                cells[CCi][CCj].type = OPEN;
                cells[UCi][UCj].type = OPEN;
                cells[UCi][UCj].context = rotate_and_flip(LC5,rot,ref);
                cells[CCi][CCj].context = rotate_and_flip(LC7,rot,ref);
                open_cells.push_back(pair<int,int> (UCi, UCj));
                open_cells.push_back(pair<int,int> (CCi, CCj));
                nopencells += 2;
                break;
                
            case 4:
                if(cells[RCi][RCj].type == DEAD || cells[RCi][RCj].type == OPEN || cells[RCi][RCj].type == X)
                    return false;
                if(!(cells[UCi][UCj].type == DEAD || cells[UCi][UCj].type == OPEN)) {
                    if(cells[UCi][UCj].type == FREE)
                        r->UwasX = false;
                    cells[UCi][UCj].type = X;
                }
                cells[RCi][RCj].type = OPEN;
                cells[RCi][RCj].context = rotate_and_flip(LC7,rot,ref);
                open_cells.push_back(pair<int,int> (RCi, RCj));
                nopencells++;
                break;
                
            case 5:
                if(cells[UCi][UCj].type == DEAD || cells[UCi][UCj].type == OPEN || cells[UCi][UCj].type == X
                   || cells[RCi][RCj].type == DEAD || cells[RCi][RCj].type == OPEN || cells[RCi][RCj].type == X)
                    return false;
                cells[RCi][RCj].type = OPEN;
                cells[UCi][UCj].type = OPEN;
                cells[RCi][RCj].context = rotate_and_flip(LC7,rot,ref);
                cells[UCi][UCj].context = rotate_and_flip(LC5,rot,ref);
                open_cells.push_back(pair<int,int> (RCi, RCj));
                open_cells.push_back(pair<int,int> (UCi, UCj));
                nopencells += 2;
                break;
        }
        return true;
    }
    
    void undo_operation_star(int Ltype, int Lsize, pair<int,int> open_cell, Retract& r, bool fix = true) {
        if(fix) {
            switch(Ltype)
            {
                case 1:
                    if(!r.RwasX)
                        cells[r.Rcell.first][r.Rcell.second].type = FREE;
                    if(!r.UwasX)
                        cells[r.Ucell.first][r.Ucell.second].type = FREE;
                    break;
                    
                case 2:
                    if(!r.RwasX)
                        cells[r.Rcell.first][r.Rcell.second].type = FREE;
                    if(!r.CwasX)
                        cells[r.Ccell.first][r.Ccell.second].type = FREE;
                    
                    // This cell must have been free
                    cells[r.Ucell.first][r.Ucell.second].type = FREE;
                    open_cells.pop_back();
                    nopencells--;
                    break;
                    
                case 3:
                    if(!r.RwasX)
                        cells[r.Rcell.first][r.Rcell.second].type = FREE;
                    
                    cells[r.Ccell.first][r.Ccell.second].type = FREE;
                    cells[r.Ucell.first][r.Ucell.second].type = FREE;
                    open_cells.pop_back();
                    open_cells.pop_back();
                    nopencells -= 2;
                    break;
                    
                case 4:
                    if(!r.UwasX)
                        cells[r.Ucell.first][r.Ucell.second].type = FREE;
                    
                    cells[r.Rcell.first][r.Rcell.second].type = FREE;
                    open_cells.pop_back();
                    nopencells--;
                    break;
                    
                case 5:
                    cells[r.Rcell.first][r.Rcell.second].type = FREE;
                    cells[r.Ucell.first][r.Ucell.second].type = FREE;
                    open_cells.pop_back();
                    open_cells.pop_back();
                    nopencells -= 2;
                    break;
            }
        }
        open_cells.push_front(open_cell);
        cells[open_cell.first][open_cell.second].type = OPEN;
        nopencells++;
        ndeadcells--;
        ncells -= (Lsize -1);
    }
};

class D2 {
public:
    int i;
    vector<long long> twigs_num;
    deque<Twig> layer;
    long long estimate;
    // double perc;
    vector <vector< vector< long long > > > coefficients;
    D2(int i_ = 1, long long e_ = -1, int cpu_ = 1) : i(i_), twigs_num(cpu_), estimate(e_ != -1 ? 4 * e_ : -1), coefficients(cpu_) {
        for(int j = 0; j < coefficients.size(); j++) {
            coefficients[j].resize(2*i+1);
            for(int k = 0; k < 2*i+1; k++) {
                coefficients[j][k].resize(i+1);
                for(int y = 0; y < i+1; y++)
                    coefficients[j][k][y] = 0;
            }
        }
        
    }
    void Aggregate(const Twig& T) {
        const int x = T.ncells - 1;
        const int y = T.ndeadcells;
        twigs_num[0]++;
        coefficients[0][x][y]++;
    }
    // This function performs a BFS traversal over the tree in Figure 7 until it fills the BFS queue
    // with as many nodes as there computing cores.
    void fillLayer(int ncpu, int tid, int run) {
        Twig T(1,0,1);
        T.cells[ROOTY][ROOTX].type = OPEN;
        T.cells[ROOTY][ROOTX].context = LC1;
        T.open_cells.push_front(pair<int,int> (ROOTY, ROOTX));
        layer.push_front(T);
        // dummies
        Retract r;
        pair<int,int> open_cell(T.open_cells.front().first, T.open_cells.front().second);
        
        // The second condition is meant for small values of i
        while(layer.size() + 5 < ncpu && layer.size() > 0) {
            Twig curr = layer.front();
            layer.pop_front();
            Twig t1(curr);
            if(t1.operation_star(1, 1, &r))
            {
                if(!t1.Done(i))
                    layer.push_back(t1);
                else if(tid == 0 && run == 0)
                    Aggregate(t1);
            }
            Twig t2(curr);
            if(t2.operation_star(2, 2, &r))
            {
                if(!t2.Done(i))
                    layer.push_back(t2);
                // The first condition is to make sure that we don't double count the twig
                // namely, only core 0 counts it (tid=0) and that we only count it in the
                // 1st run, again to prevent double counting
                else if(tid == 0 && run == 0)
                    Aggregate(t2);
            }
            Twig t3(curr);
            if(t3.operation_star(3, 3, &r))
            {
                if(!t3.Done(i))
                    layer.push_back(t3);
                else if(tid == 0 && run == 0)
                    Aggregate(t3);
            }
            Twig t4(curr);
            if(t4.operation_star(4, 2, &r))
            {
                if(!t4.Done(i))
                    layer.push_back(t4);
                else if(tid == 0 && run == 0)
                    Aggregate(t4);
            }
            Twig t5(curr);
            if(t5.operation_star(5, 3, &r))
            {
                if(!t5.Done(i))
                    layer.push_back(t5);
                else if(tid == 0 && run == 0)
                    Aggregate(t5);
            }
        }
    }
};
#endif
