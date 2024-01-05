#pragma once
#include<string>
#include<vector>
#include<utility>
//the programme is preliminary version, so just consider the simple cases --
//-- some parameters are from the input files, and the others are set default values there
class Default_Parameters
{
public:
    Default_Parameters();
	static constexpr double Default_area=4;   //m^2
	static constexpr double Default_h0=0.;      //the initial depth  
	static constexpr double Default_length=2;  //the common length of cell
    //static constexpr double Default_general_n=0.05;  //the general manning's coef  which is signed as 0 in the roughness map
	//static constexpr double Default_road_n=0.02;      //the manning's coef for road which is signed as 1 in the roughness map

    static constexpr double Default_general_n =0.05;  //the general manning's coef  which is signed as 0 in the roughness map
    static constexpr double Default_road_n = 0.05;      //the manning's coef for road which is signed as 1 in the roughness map
	
    //====================================
	//�����ļ�    
    //====================================
    std::string  inputfile_ele = "../inp_and_result/input_file/paper/paper_dem.txt";
    std::string  inputfile_roughness = "../inp_and_result/input_file/paper/paper_n.txt";
    std::string  inputfile_rain = "../inp_and_result/input_file/paper/paper_rain.txt";   

   

    //===================================
    // �������
    //===================================
    
    std::string  report_file= "../inp_and_result/result_file/paper_result_2/result";     //�ļ�ǰ׺�����ս�������ʱ���� 
    //const double report_inteval= 40;         //Ҫ���������ʱ�䲽����Ҫ��Ȼ�����ʱ�̲�����Ӧ�ļ���ʱ��  �籨����Ϊ10s������Ϊ3s
    const int report_inteval = 3000;          //�õ���������������ʱ�̡���ΪNSFDһ�㲻����ĳ������ʱ�̱������ʱ��Ϊreport_inteval*phi(h)
    std::vector<int> report_time = { };
    //std::vector<std::pair<int, int>> report_index = { {114,352},{113,352},{115,352},{114,353},{114,354}, {72,280},{ 161,189 },{128,112} ,{81,165},
    //    {127,298},{55,366},{81,165}, {99,122} };   //��ʽ�ο�
    std::vector<std::pair<int, int>> report_index = {};   //��Ϊ��ʱ����ʾ�������Ԫ����Ϣ 


    //=====================================
    // ������
    //=====================================
    const double time_step = 0.1;      //ԭ��0.01
    const double end_time = 5400; 
    const double q =0.1;	//ԭ��0.01
};

