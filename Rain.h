#pragma once
#include <map>
#include<string>
#include<utility>
#include<vector>
class Rain
{
private:        
    std::map<double/*ʱ��*/, double/*����ǿ��*/> _rain;      //����ǿ����ʱ��仯     ʱ�䣨min),����ǿ�ȣ�mm/min)       

    double  _sumRain = 0;//�ۼƽ�����(mm)
    double _presentRain = 0;//�ò�ʱ������  (mm)   

    inline  void update_sumRain(const double& pp) { _sumRain += pp; };            // �����ۼƽ�����

  //��������
    double intergrationRain(std::map<double, double> ::iterator iter1, std::map<double, double>::iterator iter2, double t1, double t2);
public:
    std::vector<std::pair<int, int> > _aimcell;                //�ý�������Ӧ��Ŀ��Ԫ����ż� ������Ϊ��-1,-1����˵����Ӧ����Ԫ��
    //method
    Rain(std::vector<std::pair<int, int> >&aimcell, std::map<double, double>& rain);//���캯��
    ~Rain();
    double cal_presentRain(const double& tb, const double& tstep);//����ò�ʱ������  mm
    inline double get_sumRain() { return _sumRain/1000.; };                 //��ȡ�ۼƽ�����    
    inline double get_presentRain() { return _presentRain/1000.; }      //��ȡ��ǰ������ ������Ϊ��
    inline  bool is_forAllCells();
    //inline bool is_forAllCells();  /1000. 
};
bool Rain::is_forAllCells()
{
     if (_aimcell[0].first == -1 && _aimcell[0].second == -1) 
         return true;
     else return false; 
}

