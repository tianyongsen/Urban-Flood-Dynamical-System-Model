#include"rain.h"
#include<map>
#include<utility>

Rain::Rain(std::vector<std::pair<int, int> >& aimcell, std::map<double, double>& rain):
    _aimcell(aimcell),
    _rain(rain)      //map的复制构造
{
};
Rain::~Rain()
{
    _rain.clear();
}



double Rain::cal_presentRain(const double& tb, const double& tstep)//计算该步时降雨量
{
    double t = tb / 60.;          //转换为分钟
    double te = (tb + tstep) / 60.;  //     转换为分钟
    double tend = _rain.rbegin()->first;

    if (te <= 0 || t >= tend) _presentRain = 0;      //降雨未开始或者降雨已结束
    else  if (tend > te && te > 0 && t < 0)
    {
        auto upper = _rain.lower_bound(te);
        _presentRain = intergrationRain(_rain.begin(), upper, 0, upper->second - te);
    }
    else if (te > tend && t < 0)
    {
        _presentRain = intergrationRain(_rain.begin(), --_rain.end(), 0, 0);
    }
    else if (0 < te && te < tend && t >= 0)
    {
        auto lower = _rain.upper_bound(t);
        auto upper = _rain.upper_bound(te);
        double dt1 = t - (--lower)->first;
        double dt2 = upper->first - te;
        _presentRain = intergrationRain(lower, upper, dt1, dt2);
    }
    else if (te > tend && 0 < t && t < tend)
    {
        auto lower = _rain.upper_bound(t);
        double dt1 = t - (--lower)->first;
        _presentRain = intergrationRain(lower, --_rain.end(), dt1, 0);

    }
    update_sumRain(_presentRain);    //更新累计水深；  为什么搞这么麻烦哪，给以后丰富留个形式    
    return _presentRain;
}

//线性多段函数积分运算
double Rain::intergrationRain(std::map<double, double> ::iterator iter1, std::map<double, double>::iterator iter2, double t1, double t2)  //迭代器的方式能否这样传递
{
    //iter1是第一个小于等于当前运行时间t的降雨时间区间点迭代器
    //iter2是第一个大于等于t+tstep的降雨时间区间点迭代器
    //t1=t-iter1->first;
    //t2=iter2-(t+tsetp);

    double sum = 0, dx = 0, dy = 0;
    std::map<double, double>::iterator temp = iter1;
    if (++temp == iter2)
    {
        dx = iter2->first - iter1->first;
        dy = iter2->second - iter1->second;
        sum = (dy / dx * (t1 - t2) + iter2->second + iter1->second) * (dx - t1 - t2) / 2.;
    }
    else
    {
        //计算后区间
        temp = iter2;
        --temp;
        dx = iter2->first - temp->first;
        dy = iter2->second - temp->second;
        sum = (iter2->second + temp->second - dy / dx * t2) * (dx - t2) / 2;
        --iter2;
        //计算前区间
        temp = iter1;
        ++temp;
        dx = temp->first - iter1->first;
        dy = temp->second - iter1->second;
        sum += (iter1->second + dy / dx * t1 + temp->second) * (dx - t1) / 2.;
        //计算中间区间     

        for (temp = ++iter1; temp != iter2; ++temp)
        {
            ++iter1;
            sum += (temp->second + (iter1)->second) * (iter1->first - temp->first) / 2.;
        }
    }
    return sum;
}





