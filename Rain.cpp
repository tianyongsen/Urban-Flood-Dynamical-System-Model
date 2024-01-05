#include"rain.h"
#include<map>
#include<utility>

Rain::Rain(std::vector<std::pair<int, int> >& aimcell, std::map<double, double>& rain):
    _aimcell(aimcell),
    _rain(rain)      //map�ĸ��ƹ���
{
};
Rain::~Rain()
{
    _rain.clear();
}



double Rain::cal_presentRain(const double& tb, const double& tstep)//����ò�ʱ������
{
    double t = tb / 60.;          //ת��Ϊ����
    double te = (tb + tstep) / 60.;  //     ת��Ϊ����
    double tend = _rain.rbegin()->first;

    if (te <= 0 || t >= tend) _presentRain = 0;      //����δ��ʼ���߽����ѽ���
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
    update_sumRain(_presentRain);    //�����ۼ�ˮ�  Ϊʲô����ô�鷳�ģ����Ժ�ḻ������ʽ    
    return _presentRain;
}

//���Զ�κ�����������
double Rain::intergrationRain(std::map<double, double> ::iterator iter1, std::map<double, double>::iterator iter2, double t1, double t2)  //�������ķ�ʽ�ܷ���������
{
    //iter1�ǵ�һ��С�ڵ��ڵ�ǰ����ʱ��t�Ľ���ʱ������������
    //iter2�ǵ�һ�����ڵ���t+tstep�Ľ���ʱ������������
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
        //���������
        temp = iter2;
        --temp;
        dx = iter2->first - temp->first;
        dy = iter2->second - temp->second;
        sum = (iter2->second + temp->second - dy / dx * t2) * (dx - t2) / 2;
        --iter2;
        //����ǰ����
        temp = iter1;
        ++temp;
        dx = temp->first - iter1->first;
        dy = temp->second - iter1->second;
        sum += (iter1->second + dy / dx * t1 + temp->second) * (dx - t1) / 2.;
        //�����м�����     

        for (temp = ++iter1; temp != iter2; ++temp)
        {
            ++iter1;
            sum += (temp->second + (iter1)->second) * (iter1->first - temp->first) / 2.;
        }
    }
    return sum;
}





