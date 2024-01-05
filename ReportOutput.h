#pragma once
#include<string>
#include<fstream>
#include<vector>
class ReportOutput
{
public:
    ReportOutput(const std::string& report_file, const int& report_inteval, std::vector<std::pair<int, int>>& report_index);
    ~ReportOutput();

    bool is_open() { return openstate; }
    void writeRunSummary(const double & runtime);   //�������ʱ��ȵ� ֮����ӹ���
    void writeResult2file(const std::string& filename, std::vector<int> &report_time,
        std::vector<std::vector<double>>& report_h, std::vector<std::pair<int, int>>& report_index,bool horizontal=false);
    bool if_report(const int  ite_num);
    bool if_Allcells_report() { return report_index.empty(); };
    const std::vector<std::pair<int, int>>& get_report_index() { return report_index; };    //���س����ã���֪���᲻��������
    void writeResult2file_every();  // ��ÿ������ʱ�̶����Ϊ��Ӧʱ�̵�txt�ļ������������ȫ��Ԫ����
    void writeResult2file_single();  //�����б���ʱ�������һ���ļ�����Ϊʱ�䣬��ΪԪ����ţ�����������ض�Ԫ����


    std::vector<std::vector<double>> report_h = {};     //�н飬�нӽ����֮�������txt�ļ�
private:
    std::string report_file= "";
    std::ofstream fout;       
    int report_inteval=0 ;                                      //����ʱ����
    int report_num = 0;
    std::vector<std::pair<int, int>> report_index = {};   //��Ϊ��ʱ����ʾ�������Ԫ����Ϣ����Ϊ��ʱԪ��Ϊ {100,20}
    bool openstate = false;
    bool _horizontal = false;    //false ʱ������Ԫ������true ʱ�����Ԫ������  
    void write_horizontal(std::vector<int>& report_time, std::vector<std::vector<double>>& report_h,
        std::vector<std::pair<int, int>>& report_index);
    void write_vertical(std::vector<int>& report_time, std::vector<std::vector<double>>& report_h,
        std::vector<std::pair<int, int>>& report_index);

    void open_file(std::string& filename);
    void close_fout();
    void writeHeading();  //д���ͷ����
    void write_report2file();  //���ļ�д�뱨������ 

 
};
