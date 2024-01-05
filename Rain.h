#pragma once
#include <map>
#include<string>
#include<utility>
#include<vector>
class Rain
{
private:        
    std::map<double/*时间*/, double/*降雨强度*/> _rain;      //降雨强度随时间变化     时间（min),降雨强度（mm/min)       

    double  _sumRain = 0;//累计降雨量(mm)
    double _presentRain = 0;//该步时降雨量  (mm)   

    inline  void update_sumRain(const double& pp) { _sumRain += pp; };            // 更新累计降雨量

  //辅助函数
    double intergrationRain(std::map<double, double> ::iterator iter1, std::map<double, double>::iterator iter2, double t1, double t2);
public:
    std::vector<std::pair<int, int> > _aimcell;                //该降雨所对应的目标元胞编号集 ，若是为（-1,-1）则说明对应所有元胞
    //method
    Rain(std::vector<std::pair<int, int> >&aimcell, std::map<double, double>& rain);//构造函数
    ~Rain();
    double cal_presentRain(const double& tb, const double& tstep);//计算该步时降雨量  mm
    inline double get_sumRain() { return _sumRain/1000.; };                 //获取累计降雨量    
    inline double get_presentRain() { return _presentRain/1000.; }      //获取当前降雨量 并换算为米
    inline  bool is_forAllCells();
    //inline bool is_forAllCells();  /1000. 
};
bool Rain::is_forAllCells()
{
     if (_aimcell[0].first == -1 && _aimcell[0].second == -1) 
         return true;
     else return false; 
}

