//
//  Twig3D.h
//  mira_twigs
//
//  Created by Mira Shalah on 9/28/17.
//  Copyright (c) 2017 Mira Shalah. All rights reserved.
//

#ifndef mira_twigs_Twig3D_h
#define mira_twigs_Twig3D_h

#include <unordered_map>
#include <vector>
#include <deque>
#include <tuple>
#include "Point3.h"
#include <assert.h>

#define ROOTX3D 10
#define ROOTY3D 10
#define ROOTZ3D 10
#define MAXSIZE3D 2*ROOTX3D

using namespace std;

typedef pair<Point3, int> openCell;

// types of a cell
typedef enum {
    FREE, DEAD, OPEN, X
} Type;

struct Retract3D {
    Point3 open_cell;
    vector<tuple<int,int,int>> was_free;
    vector<tuple<int,int,int>> cells_opened;
};

// 3-Dimensional L-context
struct L3 {
public:
    vector<Point3> points;
    L3() : points(4) {}
    L3(Point3 p1, Point3 p2, Point3 p3, Point3 p4) : points({p1,p2,p3,p4}) {}
    L3 Transform(const Matrix<int> m) {
        return L3(m*points[0], m*points[1], m*points[2], m*points[3]);
    }
    pair<int,int> sum(char* axis1 = NULL, char* axis2 = NULL) {
        Point3 s(points[0]+points[1]+points[2]+points[3]);
        if(s.X() == 0) {
            if(axis1 != NULL) {
                *axis1 = 'Y';
                *axis2 = 'Z';
            }
            return pair<int,int> (s.Y(), s.Z());
        }
        if(s.Y() == 0) {
            if(axis1 != NULL) {
                *axis1 = 'X';
                *axis2 = 'Z';
            }
            return pair<int,int> (s.X(), s.Z());
        }
        if(axis1 != NULL) {
            *axis1 = 'X';
            *axis2 = 'Y';
        }
        return pair<int,int> (s.X(), s.Y());
    }
};

std::ostream& operator<<(std::ostream& os, L3& L)
{
    os << L.points[0] << " " << L.points[1] << " " << L.points[2] << " " << L.points[3];
    return os;
};

struct cell3D
{
    Type type;
    int Lcontext;
    cell3D() : type(FREE), Lcontext(0) {}
    cell3D(Type t, int lc) : type(t), Lcontext(lc) {}
    bool operator==(const cell3D& c) const {
        return type == c.type;
    }
};

// 3-Dimensional twig
class Twig3D {
public:
    // number of cells
    int ncells;
    // number of dead cells
    int ndeadcells;
    // number of open cells
    int nopencells;
    Point3 root;
    // table
    vector<vector<vector<cell3D> > > cells;
    vector<Point3>  forbidden_cells;
    // linearly ordered open cells
    deque<Point3> open_cells;
    
    bool operator==(const Twig3D& t) const {
        if(t.ncells != ncells)
            return false;
        if(t.ndeadcells != ndeadcells)
            return false;
        if(t.nopencells != nopencells)
            return false;
        if(forbidden_cells != t.forbidden_cells)
            return false;
        if(open_cells != t.open_cells)
            return false;
        for(int i = 0; i < cells.size(); i++) {
            for(int j = 0; j < cells[i].size(); j++) {
                for(int k = 0; k < cells[i][j].size(); k++)
                    if(!(cells[i][j][k] == t.cells[i][j][k])) {
                        cout << cells[i][j][k].Lcontext << " " << cells[i][j][k].type << endl;
                        cout << t.cells[i][j][k].Lcontext << " " << t.cells[i][j][k].type << endl;
                        return false;
                    }
            }
        }
        return true;
    }
    
    // constructor
    Twig3D(int n, vector<openCell> opencells, vector<Point3> Xs, int d = 1) :
           root(Point3(ROOTX3D, ROOTY3D, ROOTZ3D)), cells(MAXSIZE3D)
    {
        ncells = n;
        ndeadcells = d;
        nopencells = opencells.size();
        for(int i = 0; i < MAXSIZE3D; i++) {
            cells[i].resize(MAXSIZE3D);
            for(int j = 0; j < MAXSIZE3D; j++)
                cells[i][j].resize(MAXSIZE3D);
        }
        cells[root.X()][root.Y()][root.Z()] = cell3D(DEAD, 0);
        for (int ind = 0; ind < opencells.size(); ind++) {
            const openCell t(opencells[ind]);
            int i = root.X() + t.first.X(),
                j = root.Y() + t.first.Y(),
                k = root.Z() + t.first.Z();
            // direction
            open_cells.push_back(Point3(t.first.X(), t.first.Y(), t.first.Z()));
            cells[i][j][k].type = OPEN;
            cells[i][j][k].Lcontext = t.second;
        }
        for(int i = 0; i < Xs.size(); i++) {
            const Point3 t(Xs[i]);
            cells[root.X() + t.X()][root.Y() + t.Y()][root.Z() + t.Z()].type = X;
            forbidden_cells.push_back(t);
        }
    }
    
    bool Done(int i) {
        if(nopencells == 0 || ndeadcells == i)
            return true;
        else
            return false;
    }
};

class D3 {
public:
    vector<L3> D3LC;
    vector<Twig3D> basic_twigs;
    int i;
    vector<long long> twigs_num;
    long long estimate;
    double perc;
    vector<Matrix<int> > maps;
    unordered_map< int , unordered_map<int, int> > LtoNYZ;
    unordered_map< int , unordered_map<int, int> > LtoNXZ;
    unordered_map< int , unordered_map<int, int> > LtoNXY;
    vector <vector< vector< long long > > > coefficients;
    deque<Twig3D> layer;
    // 24 3-dimensional L-contexts overall
    D3(int i_ = 1, long long e_ = -1, int cpu_ = 1) : D3LC(24), i(i_), twigs_num(cpu_), maps(24), estimate(e_ != -1 ? 8 * e_ : -1), perc(0.1), coefficients(cpu_) {
        for(int j = 0; j < coefficients.size(); j++) {
            coefficients[j].resize(4*i+1);
            for(int k = 0; k < 4*i+1; k++) {
                coefficients[j][k].resize(i+1);
                for(int y = 0; y < i+1; y++)
                    coefficients[j][k][y] = 0;
            }
        }
    };
    void findAll3DContexts() {
        Matrix<int> I(1,0,0,
                      0,1,0,
                      0,0,1);
        Matrix<int> rotate90X(1,0,0,
                              0,0,-1,
                              0,1,0);
        Matrix<int> rotate90Y(0,0,1,
                              0,1,0,
                              -1,0,0);
        Matrix<int> rotate90Z(0,-1,0,
                              1,0,0,
                              0,0,1);
        Matrix<int> reflectXY(1,0,0,
                              0,1,0,
                              0,0,-1);
        Matrix<int> reflectXZ(1,0,0,
                              0,-1,0,
                              0,0,1);
        Matrix<int> reflectYZ(-1,0,0,
                              0,1,0,
                              0,0,1);
        
        // L-contexts on the (initial) YZ plane
        int pi = 0;
        // the initial Lcontext
        D3LC[pi] = L3({Point3(0,-1,0), Point3(0,-1,-1), Point3(0,0,-1), Point3(0,1,-1)});
        
        Matrix<int> tmp(I);
        for(int theta = 0; theta <= 270; theta += 90) {
            D3LC[pi] = D3LC[0].Transform(tmp);
            LtoNYZ[D3LC[pi].sum().first][D3LC[pi].sum().second] = pi;
            maps[pi] = tmp;
            // Rotate by 90
            tmp = rotate90X*tmp;
            pi++;
        }
        // Reflect around XZ
        tmp = reflectXZ;
        for(int theta = 0; theta <= 270; theta += 90) {
            D3LC[pi] = D3LC[0].Transform(tmp);
            LtoNYZ[D3LC[pi].sum().first][D3LC[pi].sum().second] = pi;
            maps[pi] = tmp;
            // Rotate by 90
            tmp = rotate90X*tmp;
            pi++;
        }
    
        // L-contexts on the XZ plane
        tmp = rotate90Z;
        // Rotate by 90
        for(int theta = 0; theta <= 270; theta += 90) {
            D3LC[pi] = D3LC[0].Transform(tmp);
            LtoNXZ[D3LC[pi].sum().first][D3LC[pi].sum().second] = pi;
            maps[pi] = tmp;
            tmp = rotate90Y*tmp;
            pi++;
        }
        // Reflect
        tmp = reflectYZ*rotate90Z;
        for(int theta = 0; theta <= 270; theta += 90) {
            D3LC[pi] = D3LC[0].Transform(tmp);
            LtoNXZ[D3LC[pi].sum().first][D3LC[pi].sum().second] = pi;
            maps[pi] = tmp;
            tmp = rotate90Y*tmp;
            pi++;
        }
        
        // L-contexts on the XY plane
        tmp = rotate90X*rotate90Z;
        for(int theta = 0; theta <= 270; theta += 90) {
            D3LC[pi] = D3LC[0].Transform(tmp);
            LtoNXY[D3LC[pi].sum().first][D3LC[pi].sum().second] = pi;
            maps[pi] = tmp;
            tmp = rotate90Z*tmp;
            pi++;
        }
        tmp = reflectXZ*rotate90X*rotate90Z;
        for(int theta = 0; theta <= 270; theta += 90) {
            D3LC[pi] = D3LC[0].Transform(tmp);
            LtoNXY[D3LC[pi].sum().first][D3LC[pi].sum().second] = pi;
            maps[pi] = tmp;
            tmp = rotate90Z*tmp;
            pi++;
        }
        
        // for(L3& l3 : D3LC)
        //    cout << l3 << endl;
    }
    
    void createBasicTwigs() {
        Point3 x(1,0,0), y(0,1,0), z(0,0,1), x_(-1,0,0), c(0,1,1);
        vector<openCell> open_cells;
        vector<Point3> Xcells;
        
        /************************************* Category 1 ************************************/
        // Twig 1
        Xcells.push_back(x);
        Xcells.push_back(y);
        Xcells.push_back(z);
        Xcells.push_back(x_);
        basic_twigs.push_back(Twig3D(1, open_cells, Xcells));
        
        // Twig 2
        Xcells.clear();
        
        Xcells.push_back(y);
        Xcells.push_back(z);
        Xcells.push_back(x_);
        open_cells.push_back(openCell(Point3(1,0,0),9));
        basic_twigs.push_back(Twig3D(2, open_cells, Xcells));
        
        // Twig 3
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(x);
        Xcells.push_back(y);
        Xcells.push_back(z);
        open_cells.push_back(openCell(Point3(-1,0,0),15));
        basic_twigs.push_back(Twig3D(2, open_cells, Xcells));
        
        // Twig 4
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(y);
        Xcells.push_back(z);
        open_cells.push_back(openCell(Point3(1,0,0),9));
        open_cells.push_back(openCell(Point3(-1,0,0),15));
        basic_twigs.push_back(Twig3D(3, open_cells, Xcells));
        
        /************************************* Category 2 ************************************/
        
        // Twig 5
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(x);
        Xcells.push_back(x_);
        Xcells.push_back(z);
        open_cells.push_back(openCell(Point3(0,1,0),7));
        basic_twigs.push_back(Twig3D(2, open_cells, Xcells));
        
        // Twig 6
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(z);
        Xcells.push_back(x_);
        open_cells.push_back(openCell(Point3(1,0,0),9));
        open_cells.push_back(openCell(Point3(0,1,0),7));
        basic_twigs.push_back(Twig3D(3, open_cells, Xcells));
        
        // Twig 7
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(x);
        Xcells.push_back(z);
        open_cells.push_back(openCell(Point3(0,1,0),7));
        open_cells.push_back(openCell(Point3(-1,0,0),15));
        basic_twigs.push_back(Twig3D(3, open_cells, Xcells));
        
        // Twig 8
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(z);
        open_cells.push_back(openCell(Point3(1,0,0),9));
        open_cells.push_back(openCell(Point3(0,1,0),7)),
        open_cells.push_back(openCell(Point3(-1,0,0),15));
        basic_twigs.push_back(Twig3D(4, open_cells, Xcells));
        
        /************************************* Category 3 ************************************/
        
        // Twig 9
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(x);
        Xcells.push_back(x_);
        open_cells.push_back(openCell(Point3(0,1,0),7));
        open_cells.push_back(openCell(Point3(0,0,1),4));
        basic_twigs.push_back(Twig3D(3, open_cells, Xcells));
        
        // Twig 10
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(x_);
        open_cells.push_back(openCell(Point3(1,0,0),9));
        open_cells.push_back(openCell(Point3(0,1,0),7));
        open_cells.push_back(openCell(Point3(0,0,1),4));
        basic_twigs.push_back(Twig3D(4, open_cells, Xcells));
        
        // Twig 11
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(x);
        open_cells.push_back(openCell(Point3(0,1,0),7));
        open_cells.push_back(openCell(Point3(-1,0,0),15)),
        open_cells.push_back(openCell(Point3(0,0,1),4));
        basic_twigs.push_back(Twig3D(4, open_cells, Xcells));
        
        // Twig 12
        Xcells.clear();
        open_cells.clear();
        
        open_cells.push_back(openCell(Point3(1,0,0),9));
        open_cells.push_back(openCell(Point3(0,1,0),7));
        open_cells.push_back(openCell(Point3(-1,0,0),15));
        open_cells.push_back(openCell(Point3(0,0,1),4));
        basic_twigs.push_back(Twig3D(5, open_cells, Xcells));
        
        /************************************* Category 4 ************************************/

        // Twig 13
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(x);
        Xcells.push_back(y);
        Xcells.push_back(x_);
        Xcells.push_back(c);
        open_cells.push_back(openCell(Point3(0,0,1),4));
        basic_twigs.push_back(Twig3D(2, open_cells, Xcells));
        
        // Twig 14
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(x);
        Xcells.push_back(y);
        Xcells.push_back(x_);
        open_cells.push_back(openCell(Point3(0,0,1),4));
        open_cells.push_back(openCell(Point3(0,1,1),7));
        basic_twigs.push_back(Twig3D(3, open_cells, Xcells));
        
        // Twig 15
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(y);
        Xcells.push_back(x_);
        open_cells.push_back(openCell(Point3(1,0,0),9));
        open_cells.push_back(openCell(Point3(0,0,1),8));
        basic_twigs.push_back(Twig3D(3, open_cells, Xcells));
        
        // Twig 16
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(x);
        Xcells.push_back(y);
        open_cells.push_back(openCell(Point3(-1,0,0),15));
        open_cells.push_back(openCell(Point3(0,0,1),12));
        basic_twigs.push_back(Twig3D(3, open_cells, Xcells));
        
        // Twig 17
        Xcells.clear();
        open_cells.clear();
        
        Xcells.push_back(y);
        open_cells.push_back(openCell(Point3(1,0,0),9));
        open_cells.push_back(openCell(Point3(-1,0,0),15));
        open_cells.push_back(openCell(Point3(0,0,1),8));
        basic_twigs.push_back(Twig3D(4, open_cells, Xcells));
    }
    
    void Aggregate(const Twig3D& T) {
        const int x = T.ncells - 1;
        const int y = T.ndeadcells;
        twigs_num[0]++;
        coefficients[0][x][y]++;
    }
};

// TODO: make alg and basic_twig const
bool operation_star(D3& alg, Twig3D* conf, const Twig3D& basic_twig, Retract3D* r)
{
    //assert(r->cells_opened.size()==0);
    // get the direction of the first open cell according to the linear order
    Point3 conf_open_cell = conf->open_cells.front();
    r->open_cell = conf_open_cell;
    
    int opencell1i = conf_open_cell.X(),
        opencell1j = conf_open_cell.Y(),
        opencell1k = conf_open_cell.Z();
    // open cell is now dead
    conf->open_cells.pop_front();
    conf->nopencells--;
    conf->ndeadcells++;
    
    int grid_cordi = opencell1i + ROOTX3D,
        grid_cordj = opencell1j + ROOTY3D,
        grid_cordk = opencell1k + ROOTZ3D;
    conf->cells[grid_cordi][grid_cordj][grid_cordk].type = DEAD;
    
    // index of the L conetxt of the open cell
    const int context = conf->cells[grid_cordi][grid_cordj][grid_cordk].Lcontext;
    // transformation to be applied to the basic twig
    Matrix<int> mapping = alg.maps[context];
    
    // make the forbidden cells of the basic twig forbidden cells of the twig
    for(int j = 0; j < basic_twig.forbidden_cells.size(); j++) {
        const Point3 p(basic_twig.forbidden_cells[j]);
        Point3 np = mapping * p;
        // no need: conf->forbidden_cells.push_back(p);
        Point3 forbidden_cell(grid_cordi + np.X(),
                              grid_cordj + np.Y(),
                              grid_cordk + np.Z());
        int fi = forbidden_cell.X(), fj = forbidden_cell.Y(), fk = forbidden_cell.Z();
        // set to forbidden
        if(conf->cells[fi][fj][fk].type == FREE) {
            conf->cells[fi][fj][fk].type = X;
            r->was_free.push_back(tuple<int,int,int>(fi, fj, fk));
        }
    }
    
    // make the open cells of the basic twig open cells of the configuration
    deque<Point3> basic_open_cells = basic_twig.open_cells;
    while(!basic_open_cells.empty()) {
        // direction of the open cell
        Point3 dir = basic_open_cells.front();
        int original_context = basic_twig.cells[dir.X()+ROOTX3D]
                                               [dir.Y()+ROOTY3D]
                                               [dir.Z()+ROOTZ3D].Lcontext;
        basic_open_cells.pop_front();
        Point3 p(dir);
        // transform direction
        p = mapping * p;
        Point3 open_cell(grid_cordi + p.X(),
                         grid_cordj + p.Y(),
                         grid_cordk + p.Z());
        int oci = open_cell.X(), ocj = open_cell.Y(), ock = open_cell.Z();
        if(conf->cells[oci][ocj][ock].type == DEAD ||
           conf->cells[oci][ocj][ock].type == OPEN ||
           conf->cells[oci][ocj][ock].type == X) {
            return false;
        }
        L3 c = alg.D3LC[original_context],
        newc = c.Transform(mapping);
        char axis1, axis2;
        pair<int,int> s = newc.sum(&axis1, &axis2);
        int new_context;
        if(axis1 == 'X') {
            if(axis2 == 'Y')
                new_context = alg.LtoNXY[s.first][s.second];
            else
                new_context = alg.LtoNXZ[s.first][s.second];
        } else {
            new_context = alg.LtoNYZ[s.first][s.second];
        }
        conf->cells[oci][ocj][ock].type = OPEN;
        conf->cells[oci][ocj][ock].Lcontext = new_context;
        conf->open_cells.push_back(conf_open_cell + p);
        r->cells_opened.push_back(tuple<int,int,int>(oci, ocj, ock));
        conf->nopencells++;
        conf->ncells++;
    }
    return true;
}

void undo_operation_star(Twig3D* conf, Retract3D* r) {
    for(int i = 0; i < r->was_free.size(); i++) {
        int x = get<0>(r->was_free[i]),
            y = get<1>(r->was_free[i]),
            z = get<2>(r->was_free[i]);
        conf->cells[x][y][z].type = FREE;
    }
    // These cells must have been free
    for(int i = 0; i < r->cells_opened.size(); i++) {
        conf->open_cells.pop_back();
        int x = get<0>(r->cells_opened[i]),
            y = get<1>(r->cells_opened[i]),
            z = get<2>(r->cells_opened[i]);
        conf->cells[x][y][z].type = FREE;
    }
    conf->open_cells.push_front(r->open_cell);
    int grid_cordi = r->open_cell.X() + ROOTX3D,
        grid_cordj = r->open_cell.Y() + ROOTY3D,
        grid_cordk = r->open_cell.Z() + ROOTZ3D;
    conf->cells[grid_cordi][grid_cordj][grid_cordk].type = OPEN;
    conf->nopencells++;
    conf->ndeadcells--;
    conf->ncells -= r->cells_opened.size();
    conf->nopencells -= r->cells_opened.size();
    r->cells_opened.clear();
    r->was_free.clear();
}

// This function performs a BFS traversal over the tree in Figure 7 (but for 3D) until it fills the BFS
// queue with as many nodes as there computing cores.
void fillLayer(int ncpu, int tid, int run, D3* alg) {
    vector<openCell> open_cell;
    vector<Point3> Xs;
    open_cell.push_back(pair<Point3,int>(Point3(0,0,0),0));
    Twig3D T3(1, open_cell, Xs, 0);
    alg->layer.push_front(T3);
    // The second condition is meant for small values of i
    while(alg->layer.size() + 17 < ncpu && alg->layer.size() > 0) {
        Twig3D curr = alg->layer.front();
        alg->layer.pop_front();
        for(int i = 0; i < alg->basic_twigs.size(); i++) {
            Twig3D t_i(curr);
            Retract3D r;
            if(operation_star(*alg, &t_i, alg->basic_twigs[i], &r)) {
                if(!t_i.Done(alg->i))
                    alg->layer.push_back(t_i);
                // The first condition is to make sure that we don't double count the twig
                // namely, only core 0 counts it (tid=0) and that we only count it in the
                // 1st run, again to prevent double counting
                else if(tid == 0 && run == 0)
                    alg->Aggregate(t_i);
            }
        }
    }
};

#endif
