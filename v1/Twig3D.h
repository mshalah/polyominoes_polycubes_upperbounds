//
//  Twig3D.h
//  mira_twigs
//
//  Created by Mira Shalah on 9/28/17.
//

#ifndef mira_twigs_Twig3D_h
#define mira_twigs_Twig3D_h

#include <unordered_map>
#include "utils.h"
#include "Twig.h"
#include <omp.h>

#define ROOTX3D 10
#define ROOTY3D 10
#define ROOTZ3D 10
#define MAXSIZE3D 2*ROOTX3D

using namespace std;

typedef pair<Point3, int> openCell;

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
}

struct cell3D
{
    Type type;
    int Lcontext;
    cell3D() : type(FREE), Lcontext(0) {}
    cell3D(Type t, int lc) : type(t), Lcontext(lc) {}
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
    queue<Point3>  open_cells;
    
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
            open_cells.push(Point3(t.first.X(), t.first.Y(), t.first.Z()));
            cells[i][j][k].type = OPEN;
            cells[i][j][k].Lcontext = t.second;
        }
        for(int i = 0; i < Xs.size(); i++) {
            const Point3 t(Xs[i]);
            cells[root.X() + t.X()][root.Y() + t.Y()][root.Z() + t.Z()].type = X;
            forbidden_cells.push_back(t);
        }
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
    vector<unordered_map<int, unordered_map<int, long long> > > coefficients;
    // 24 3-dimensional L-contexts overall
    D3(int i_ = 1, long long e_ = -1, int cpu_ = 1) : D3LC(24), i(i_), twigs_num(cpu_), maps(24), estimate(e_ != -1 ? 8 * e_ : -1), perc(0.1), coefficients(cpu_) {};
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
        open_cells.push_back(openCell(Point3(-1,0,0),15)),
        open_cells.push_back(openCell(Point3(0,0,1),8));
        basic_twigs.push_back(Twig3D(4, open_cells, Xcells));
    }
};

bool operation_star(D3& alg, Twig3D* conf, const Twig3D& basic_twig)
{
    // get the direction of the first open cell according to the linear order
    Point3 conf_open_cell = conf->open_cells.front();
    int opencell1i = conf_open_cell.X(),
    opencell1j = conf_open_cell.Y(),
    opencell1k = conf_open_cell.Z();
    // open cell is now dead
    conf->open_cells.pop();
    
    int grid_cordi = opencell1i + ROOTX3D,
    grid_cordj = opencell1j + ROOTY3D,
    grid_cordk = opencell1k + ROOTZ3D;
    conf->cells[grid_cordi][grid_cordj][grid_cordk].type = DEAD;
    
    conf->nopencells--;
    conf->ndeadcells++;
    conf->ncells += (basic_twig.ncells - 1);
    
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
        if(!(conf->cells[fi][fj][fk].type == DEAD || conf->cells[fi][fj][fk].type == OPEN))
            conf->cells[fi][fj][fk].type = X;
    }
    
    // make the open cells of the basic twig open cells of the configuration
    queue<Point3>  basic_open_cells = basic_twig.open_cells;
    while(!basic_open_cells.empty()) {
        // direction of the open cell
        Point3 dir = basic_open_cells.front();
        int original_context = basic_twig.cells[dir.X()+ROOTX3D]
        [dir.Y()+ROOTY3D]
        [dir.Z()+ROOTZ3D].Lcontext;
        basic_open_cells.pop();
        Point3 p(dir);
        // transform direction
        p = mapping * p;
        Point3 open_cell(grid_cordi + p.X(),
                         grid_cordj + p.Y(),
                         grid_cordk + p.Z());
        int oci = open_cell.X(), ocj = open_cell.Y(), ock = open_cell.Z();
        if(conf->cells[oci][ocj][ock].type == DEAD ||
           conf->cells[oci][ocj][ock].type == OPEN ||
           conf->cells[oci][ocj][ock].type == X)
            return false;
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
        conf->open_cells.push(conf_open_cell + p);
        conf->nopencells++;
    }
    return true;
}

void Recurse3D(const Twig3D& T, D3* alg)
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
    
    Twig3D newT(T);
#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[0]))
        {
            Recurse3D(newT, alg);
        }
    }
    
#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[1]))
        {
            Recurse3D(newT, alg);
        }
    }

#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[2]))
        {
            Recurse3D(newT, alg);
        }
    }
    
#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[3]))
        {
            Recurse3D(newT, alg);
        }
    }

#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[4]))
        {
            Recurse3D(newT, alg);
        }
    }

#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[5]))
        {
            Recurse3D(newT, alg);
        }
    }

#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[6]))
        {
            Recurse3D(newT, alg);
        }
    }

#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[7]))
        {
            Recurse3D(newT, alg);
        }
    }
    
#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[8]))
        {
            Recurse3D(newT, alg);
        }
    }

#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[9]))
        {
            Recurse3D(newT, alg);
        }
    }

#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[10]))
        {
            Recurse3D(newT, alg);
        }
    }
    
#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[11]))
        {
            Recurse3D(newT, alg);
        }
    }

#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[12]))
        {
            Recurse3D(newT, alg);
        }
    }

#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[13]))
        {
            Recurse3D(newT, alg);
        }
    }

#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[14]))
        {
            Recurse3D(newT, alg);
        }
    }

#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[15]))
        {
            Recurse3D(newT, alg);
        }
    }

#pragma omp task firstprivate(newT)
    {
        if(operation_star(*alg, &newT, alg->basic_twigs[16]))
        {
            Recurse3D(newT, alg);
        }
    }
};

#endif
