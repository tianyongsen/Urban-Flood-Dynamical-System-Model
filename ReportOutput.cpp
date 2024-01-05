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
        //���ܣ���ָ���ļ������ļ�
    {      
        close_fout();     //���ȹر��Ѱ��ļ�        
        fout.open(filename.c_str(), std::ios::out | std::ios::trunc);  //ios::trunc   ��truncate ��һ���ļ�������Ѿ�������ɾ��ԭ�������ݡ���ʵΪĬ��ѡ�񣬿��Բ�д��
        if (!fout.is_open()) throw  "FILE ERROR: can not reportfile";
        openstate = true;     
    }
    void ReportOutput::close_fout()
        //���ܣ��ر��ļ�������
    {
        fout.close();     //�ر��ļ������
        openstate = false;
    }
    void ReportOutput::writeRunSummary(const double& runtime)   //�������ʱ��ȵ� 
    {
        fout << std::left;     //�����
        fout << "the  all runtime is :" << std::setprecision(3) <<runtime << " sec." << std::endl;;
    }
    void  ReportOutput::writeResult2file(const std::string& filename, std::vector<int>& report_time,
        std::vector<std::vector<double>>& report_h, std::vector<std::pair<int, int>>& report_index, 
        bool horizontal /*= false*/)
        //report_h��¼��ָ��Ԫ����ָ��ʱ�̵�ˮ����Ϣ ���У���Ϊtime, ��Ϊָ��Ԫ��        
        //
    {
        report_file= filename;
        _horizontal = horizontal;
        open_file(report_file);    //���ļ�
        if (horizontal)
            write_horizontal(report_time, report_h, report_index);
        else
            write_vertical(report_time, report_h, report_index);
    }

    bool ReportOutput::if_report(const int  ite_num)
    {
        //����дΪ��������
        if (ite_num == report_num)   // ��ǰ�����������Ƿ񱨸档ע��0ʱ���Ǳ���ģ�          
            return true;    
        return false;           
    }

    void ReportOutput::writeResult2file_every()  // ��ÿ������ʱ�̶����Ϊ��Ӧʱ�̵�txt�ļ������������ȫ��Ԫ����
    //���ܣ�ƴ��report_file��report_num�����µ��ļ���������report_h���뵽���ļ�
    {
        //1.�����������Ĭ�� report_h�Ѿ��н�������  
        if (report_h.empty())
            throw "report_h is empty, there is nothing to report";                   //Ϊ��ֱ���׳��쳣 
        //2.��Ϊ��ʱ
            //2.1ƴ��report_file��report_num�����µ��ļ��������� 
        std::string current_file = report_file +std::to_string(report_num)+".txt";    //���ļ���
        open_file(current_file);   //  ���ļ�
            //2.2 д���ͷ����
        writeHeading();
            //2.3д�뱨������
        write_report2file();      
            //д����ɺ󣬹ر��ļ���������������ļ�
        close_fout();
            
    //3 ���±������  
        report_num += report_inteval;   //���µ���һ���������
    }
    void ReportOutput::writeResult2file_single()  //�����б���ʱ�������һ���ļ�����Ϊʱ�䣬��ΪԪ����ţ�����������ض�Ԫ��
    {

    }

    void ReportOutput::writeHeading()  //д���ͷ����
        //���ܣ�����report_h д�� ��ͷ����
    {
        //����ú�����Ĭ��fout���Ѿ���
        if (!openstate)
            throw "in writeHeading(),the file is not openning ";
        //ncols         483
        //    nrows         201
        int    xllcorner = 263976;
        int    yllcorner = 664410;
        //    cellsize      2
        int    NODATA_value = -9999;
        fout << std::left;     //�����
        fout << "ncols" << "      " << report_h[0].size() << std::endl;   //д������
        fout << "nrows" << "      " << report_h.size() << std::endl;   //д������
         //д������ ����Ҫ���ݲ�ͬ�������ġ�֮����Ҫ��������Ϊ�˷����ASCIIתΪRaster
        fout << "xllcorner" << "      " << xllcorner << std::endl;
        fout << "yllcorner" << "      " << yllcorner<< std::endl; 
        fout << "cellsize" << "      " << Default_Parameters::Default_length << std::endl;   //д��Ԫ���ߴ�
        fout << "NODATA_value" << "      " << NODATA_value << std::endl;   //д�������ݱ�־
    }

    void ReportOutput::write_report2file()
        //���ܣ����ļ�д�뱨������
    {        
        //ѭ�����
        for (int i = 0; i < report_h.size(); i++)
        {
            for (int j = 0; j < report_h[0].size(); j++)
            {
                fout << std::left;    //�����
                double& x = report_h[i][j];
                double absX = abs(x);
                //�������ֵ̫С�������ֱ����Ϊ0
                if (absX < 1.0e-4)   x = 0.0;   //0.1mm����Ϊ0;Ϊ�˲�������︺���쳣ֵ������ȡΪ[-0.1,0.1] mm 
                //�������ֵ̫����������ѧ����
                if (absX > 1.0e5)
                {
                    fout << std::scientific;
                }
                else
                {
                    fout << std::fixed;        //������һ��������ݵĸ�ʽ
                }     //һ������̶��������  
                fout << std::setprecision(4) << x;    //����������ȣ�  <<fiexd<<std::setprecision(4)<<x Ϊ����С�����4λ��
                fout << " ";   //����ո�
            }
            if(i != report_h.size()-1)   //i�������һ�У����е���һ�У�i�����һ�У������У����������
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
        
        fout << std::left;     //�����   
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

    



