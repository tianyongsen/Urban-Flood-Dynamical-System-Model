#pragma once
#include<vector>
#include<string>
#include<unordered_map>
#include<map>
#include <time.h>
#include"Cell.h"
#include"Input.h"
#include"Rain.h"
//#include"Four_tuple.h"
//#include"NSFD_solver.h"
typedef  std::vector<std::vector<double>>  Matrix;
struct Four_tuple;
size_t Four_tuple_hash(const Four_tuple& ft);	
class NSFD_solver;     //for  ��Ԫ
//class Cell;
class Cellfield
{
private:
	std::vector<std::vector<Cell*>>          _cells;   //ֻ������򵥵ľ���Ԫ����������ڷ��Ӿ������й���	
	int							            		         _rows=0;  //����
	int														_cols = 0;  //����
	std::vector<Rain*>                               _rains;    //����  
	Input												     _input;	

	////Ԫ���߽��� �ͱ��ϵ�ͨ�� �� ������Լ�����߽����ɵ�ֵָ���ֵ...
	//... ����Ԫ�����򣬾�����Ȼ����ṹ�����ֽṹ�����������󣩽ṹ��Ӧ�����Գ���������ֶ�Ӧ��ϵ...
	//...���㷨�߼��ͼ��㸴�Ӷȡ���������ṹ�Ͷ�Ӧ��ϵҲ�ᱻNSFD������̳С�
	std::vector<std::vector<double> >      _links ;     //�������̺��˶Ծ���߽������ӳ�� 
	//std::unordered_map<const Four_tuple, double , hash_name> _emap;   //��������ó���̫��������ֱ���þ���

	void construct_cells();
	void construct_rains();
	void construct_links();
	//void construct_emap(); 


public:
	friend NSFD_solver;    //��ʱ�������ϸ��װ��������Ԫȡ�ñ���
	void construct_Cellfiled(std::string& ele_file, std::string& roughness_file, std::string& rain_file);	
	~Cellfield();	
};







