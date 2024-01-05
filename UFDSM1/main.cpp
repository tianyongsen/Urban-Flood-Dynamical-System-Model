// CA_NSFD_simple_v.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include<string>
#include <vector>
#include <iomanip>
#include<sstream>
#include<stdexcept>
#include"Cellfield.h"
#include"NSFD_solver.h"
#include"Default_Parameters.h"
#include"ReportOutput.h"
inline void cal_runtime(clock_t& start_t, clock_t& end_t, std::ostringstream& msg);
int main()
{    
   
    Default_Parameters  paras;  
    ReportOutput    report(paras.report_file,paras.report_inteval, paras.report_index);  //初始化报告对象
    clock_t start_t = clock();  // 初始化开始时间
    try {
        //主程序
        //======================================================
        Cellfield   cellfield;
        //创建cellfiled 
        cellfield.construct_Cellfiled(paras.inputfile_ele, paras.inputfile_roughness, paras.inputfile_rain);
        //创建求解器
        NSFD_solver nsfd_slover(cellfield, paras.time_step, paras.end_time, paras.q);
        //求解并将结束输出到报告文件
        nsfd_slover.run_NSFD_solver(report);
    }
    catch (std::exception& e)
        {
            //Project::writerErrorsMessage();     
            //Project::writeErrorsMessage();       //怎么向外输出信息哪? 
            //这里暂时没有向报告文件中输出信息，之后补充吧
            std::cout << "\nThere were errors";
            std::cout << '\n' << e.what() << std::endl;     //窗口输出
        }
    clock_t end_t = clock();   //结束时间
    std::ostringstream msg;
    cal_runtime(start_t, end_t, msg);  //
    //msg<<  可添加其他输出信息。 
    std::cout << msg.str() << std::endl; 
}

inline void cal_runtime(clock_t& start_t, clock_t& end_t, std::ostringstream& msg)
{
    double cpu_t = (double(end_t) - double(start_t)) / CLOCKS_PER_SEC;
    msg << "\n  Simulation completed in ";
    if (cpu_t < 0.001) msg << "< 0.001 sec.";
    else  msg << std::setprecision(3) << cpu_t << " sec.";
}

