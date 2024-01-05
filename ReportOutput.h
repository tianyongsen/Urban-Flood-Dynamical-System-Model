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
    void writeRunSummary(const double & runtime);   //输出运行时间等等 之后添加功能
    void writeResult2file(const std::string& filename, std::vector<int> &report_time,
        std::vector<std::vector<double>>& report_h, std::vector<std::pair<int, int>>& report_index,bool horizontal=false);
    bool if_report(const int  ite_num);
    bool if_Allcells_report() { return report_index.empty(); };
    const std::vector<std::pair<int, int>>& get_report_index() { return report_index; };    //返回常引用，不知道会不会有问题
    void writeResult2file_every();  // 将每个报告时刻都输出为相应时刻的txt文件，适用于输出全部元胞。
    void writeResult2file_single();  //将所有报告时刻输出到一个文件，行为时间，列为元胞编号，适用于输出特定元胞。


    std::vector<std::vector<double>> report_h = {};     //中介，承接结果，之后输出到txt文件
private:
    std::string report_file= "";
    std::ofstream fout;       
    int report_inteval=0 ;                                      //报告时间间隔
    int report_num = 0;
    std::vector<std::pair<int, int>> report_index = {};   //当为空时，表示输出所有元胞信息。不为空时元素为 {100,20}
    bool openstate = false;
    bool _horizontal = false;    //false 时间纵向，元胞横向。true 时间横向，元胞纵向  
    void write_horizontal(std::vector<int>& report_time, std::vector<std::vector<double>>& report_h,
        std::vector<std::pair<int, int>>& report_index);
    void write_vertical(std::vector<int>& report_time, std::vector<std::vector<double>>& report_h,
        std::vector<std::pair<int, int>>& report_index);

    void open_file(std::string& filename);
    void close_fout();
    void writeHeading();  //写入标头内容
    void write_report2file();  //向文件写入报告内容 

 
};
