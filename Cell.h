#pragma once
#include<vector>
//#include"NSFD_solver.h"   //仅仅需要声明就行，不用把头文件包含进来
//#include"Rules_set.h"
class NSFD_solver;   //for 友元
class Rules_set;       //for 友元
class Cell
{
public:         
    friend NSFD_solver;    //暂时不考虑严格封装，先用友元取得便利 
    friend Rules_set;                    
    Cell(const double& elevation,const double &h0,const double& area,const double& roughness);    

private:
    //const int               _index;                              //the index of cell   0 1 2 3...
    const int               _type= GROUND;               //the type of the cell      
    const double        _ele;                                    //元胞相对于某基准面的形心处高度。(m)
    const double        _area;                                  //the area of the cell(m2)    
    const double         _n;            //the manning's coefficient              
    double                 _ht = 0.;                //在t时刻的深度。也用来表示初始深度。(m)    
    ////double                 _hnt = 0.;         //在t+1时刻的深度
    //double                 _vt = 0.;             //t时刻的水量，也表示初始水量
    //double                 _vnt = 0.;              //t+1时刻的水量
    double                 _sum_g= 0.;               // (G-G^T)I 项 
    //const  double        _dstorage;        //depression storage  (m)   算法里面暂时没有包括
    

public:
    enum CELLTYPE
    {        
        GROUND,
        BUILDING,
        GRASS,
        WATER,
        ROAD,
        INDUSTRIAL,
        JUNCTION,
        BOUNDARY
    };
};







