#pragma once
#include<sstream>
#include<iomanip>
#include<cmath>
#include<string>
#include<iomanip>
#include"Default_Parameters.h"  


#include"ReportOutput.h"

    static const int width = 12;
    static const int precis = 3;
    ReportOutput::ReportOutput(const std::string& report_file,const int& report_inteval, std::vector<std::pair<int, int>>& report_index):
        report_file(report_file),report_inteval(report_inteval),report_index(report_index) 
    {    
            
    };
    ReportOutput::~ReportOutput()
    {}

    

    void ReportOutput::open_file(std::string& filename)
        //功能：打开指定文件名的文件
    {      
        close_fout();     //首先关闭已绑定文件        
        fout.open(filename.c_str(), std::ios::out | std::ios::trunc);  //ios::trunc   即truncate 打开一个文件，如果已经存在则删除原来的内容。其实为默认选择，可以不写。
        if (!fout.is_open()) throw  "FILE ERROR: can not reportfile";
        openstate = true;     
    }
    void ReportOutput::close_fout()
        //功能：关闭文件输入流
    {
        fout.close();     //关闭文件输出流
        openstate = false;
    }
    void ReportOutput::writeRunSummary(const double& runtime)   //输出运行时间等等 
    {
        fout << std::left;     //左对齐
        fout << "the  all runtime is :" << std::setprecision(3) <<runtime << " sec." << std::endl;;
    }
    void  ReportOutput::writeResult2file(const std::string& filename, std::vector<int>& report_time,
        std::vector<std::vector<double>>& report_h, std::vector<std::pair<int, int>>& report_index, 
        bool horizontal /*= false*/)
        //report_h记录着指定元胞在指定时刻的水深信息 其中：行为time, 列为指定元胞        
        //
    {
        report_file= filename;
        _horizontal = horizontal;
        open_file(report_file);    //打开文件
        if (horizontal)
            write_horizontal(report_time, report_h, report_index);
        else
            write_vertical(report_time, report_h, report_index);
    }

    bool ReportOutput::if_report(const int  ite_num)
    {
        //可以写为内联函数
        if (ite_num == report_num)   // 当前迭代计数下是否报告。注：0时刻是报告的，          
            return true;    
        return false;           
    }

    void ReportOutput::writeResult2file_every()  // 将每个报告时刻都输出为相应时刻的txt文件，适用于输出全部元胞。
    //功能：拼接report_file和report_num生成新的文件名，并将report_h输入到该文件
    {
        //1.进入这个函数默认 report_h已经承接了数据  
        if (report_h.empty())
            throw "report_h is empty, there is nothing to report";                   //为空直接抛出异常 
        //2.不为空时
            //2.1拼接report_file和report_num生成新的文件名，并打开 
        std::string current_file = report_file +std::to_string(report_num)+".txt";    //新文件名
        open_file(current_file);   //  打开文件
            //2.2 写入表头内容
        writeHeading();
            //2.3写入报告内容
        write_report2file();      
            //写入完成后，关闭文件输出流，并保存文件
        close_fout();
            
    //3 更新报告计数  
        report_num += report_inteval;   //更新到下一个报告计数
    }
    void ReportOutput::writeResult2file_single()  //将所有报告时刻输出到一个文件，行为时间，列为元胞编号，适用于输出特定元胞
    {

    }

    void ReportOutput::writeHeading()  //写入标头内容
        //功能：根据report_h 写入 标头内容
    {
        //进入该函数，默认fout，已经打开
        if (!openstate)
            throw "in writeHeading(),the file is not openning ";
        //ncols         483
        //    nrows         201
        int    xllcorner = 263976;
        int    yllcorner = 664410;
        //    cellsize      2
        int    NODATA_value = -9999;
        fout << std::left;     //左对齐
        fout << "ncols" << "      " << report_h[0].size() << std::endl;   //写入列数
        fout << "nrows" << "      " << report_h.size() << std::endl;   //写入行数
         //写入坐标 ，需要根据不同算例更改。之所以要输出这个是为了方便从ASCII转为Raster
        fout << "xllcorner" << "      " << xllcorner << std::endl;
        fout << "yllcorner" << "      " << yllcorner<< std::endl; 
        fout << "cellsize" << "      " << Default_Parameters::Default_length << std::endl;   //写入元胞尺寸
        fout << "NODATA_value" << "      " << NODATA_value << std::endl;   //写入无数据标志
    }

    void ReportOutput::write_report2file()
        //功能：向文件写入报告内容
    {        
        //循环输出
        for (int i = 0; i < report_h.size(); i++)
        {
            for (int j = 0; j < report_h[0].size(); j++)
            {
                fout << std::left;    //左对齐
                double& x = report_h[i][j];
                double absX = abs(x);
                //处理绝对值太小的情况，直接置为0
                if (absX < 1.0e-4)   x = 0.0;   //0.1mm以下为0;为了不阻拦表达负的异常值，这里取为[-0.1,0.1] mm 
                //处理绝对值太大的情况，科学计数
                if (absX > 1.0e5)
                {
                    fout << std::scientific;
                }
                else
                {
                    fout << std::fixed;        //设置下一个输出数据的格式
                }     //一般情况固定长度输出  
                fout << std::setprecision(4) << x;    //设置输出精度，  <<fiexd<<std::setprecision(4)<<x 为保留小数点后4位。
                fout << " ";   //输出空格
            }
            if(i != report_h.size()-1)   //i不是最后一行，则换行到下一行，i是最后一行，不换行，结束输出。
                fout << std::endl;  
        }
    }

    void ReportOutput::write_horizontal(std::vector<int>& report_time, std::vector<std::vector<double>>& report_h,
        std::vector<std::pair<int, int>>& report_index)
    {

    }
    void ReportOutput::write_vertical(std::vector<int>& report_time, std::vector<std::vector<double>>& report_h,
        std::vector<std::pair<int, int>>& report_index)
    {
        
        fout << std::left;     //左对齐   
        fout << "time |                                 reported cells                           | ";
        fout << std::endl;
        fout << "----------------------------------------------------------------------";
        fout << std::endl;
        //================
        int time_num = report_time.size();
        for (int i = 0; i < time_num; i++)
        {
            fout << 1;
        }

        
    }

    



