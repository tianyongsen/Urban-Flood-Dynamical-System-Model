#pragma once
#include<vector>
//#include"NSFD_solver.h"   //������Ҫ�������У����ð�ͷ�ļ���������
//#include"Rules_set.h"
class NSFD_solver;   //for ��Ԫ
class Rules_set;       //for ��Ԫ
class Cell
{
public:         
    friend NSFD_solver;    //��ʱ�������ϸ��װ��������Ԫȡ�ñ��� 
    friend Rules_set;                    
    Cell(const double& elevation,const double &h0,const double& area,const double& roughness);    

private:
    //const int               _index;                              //the index of cell   0 1 2 3...
    const int               _type= GROUND;               //the type of the cell      
    const double        _ele;                                    //Ԫ�������ĳ��׼������Ĵ��߶ȡ�(m)
    const double        _area;                                  //the area of the cell(m2)    
    const double         _n;            //the manning's coefficient              
    double                 _ht = 0.;                //��tʱ�̵���ȡ�Ҳ������ʾ��ʼ��ȡ�(m)    
    ////double                 _hnt = 0.;         //��t+1ʱ�̵����
    //double                 _vt = 0.;             //tʱ�̵�ˮ����Ҳ��ʾ��ʼˮ��
    //double                 _vnt = 0.;              //t+1ʱ�̵�ˮ��
    double                 _sum_g= 0.;               // (G-G^T)I �� 
    //const  double        _dstorage;        //depression storage  (m)   �㷨������ʱû�а���
    

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







