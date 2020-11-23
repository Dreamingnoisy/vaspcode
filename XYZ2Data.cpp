/*XYZ2Data命令
* 描述：将.xyz格式转换成Lammps Data格式
* 语法：XYZ2Data -i 输入文件;-o 输出文件;-d data文件的头描述
* 版本：Prototype
* 程序范式：面向过程
* 作者：DingChangjie [Ding2020(at)mail.ustc.edu.cn]
* TODO:修正命令行参数不够robust的缺陷
*/
# include <iostream>
using std::cin; using std::cout; using std::endl;
# include <string>
using std::string;
# include <sstream>
using std::stringstream;
# include <fstream>
using std::ifstream; using std::ofstream;
#include <set>
using std::set;
# include <vector>
using std::vector;

int main(int argc, char **argv)
{
    string argument[3]; //存放命令行参数
    //命令行参数解析：将参数按顺序重新排列
    for (int i=1;i<argc;i++)
    {
        string tmp = argv[i];//强制类型转换，便于用if判断
        if (tmp == "-i")
        {
            argument[0] = argv[i+1]; // -i
        }
        else if (tmp == "-o")
        {
            argument[1] = argv[i+1]; // -o
        }
        else if (tmp == "-d")
        {
            argument[2] = argv[i+1]; //-d
        }        
    }
    
    //通过流打开.xyz文件
    ifstream datain(argument[0]);
    //扫描输入文件，获取必要信息
    string read_line; //read_line用于保存当前读取的行
    int num_atoms; //保存.xyz文件给出的原子数
    if (getline(datain,read_line))
    {
        stringstream tmp;
        tmp << read_line;
        tmp >> num_atoms;
    }
    getline(datain,read_line); //跳过.xyz的注释行
    //自动获取atomtype信息
    set<string> atom_types; //将atom type定义为一个集合
    while(getline(datain,read_line))
    {
        string tmp_type = read_line.substr(0,2);//获得当前行的type信息。注意元素符号至多两个字符
        atom_types.insert(tmp_type);//再将其写入集合  
    }
    int num_types = atom_types.size();//保存.xyz文件给出的type数
    vector<string> types_map;//用于集合映射
    set<string> ::iterator itr = atom_types.begin();
    for (itr=atom_types.begin(); itr != atom_types.end(); itr++)//执行集合映射，将元素符号映射为type编号
    {
        types_map.push_back(*itr);
    }
    //流指针送回数据部分
    datain.clear();
    datain.seekg(0,std::ios::beg);
    getline(datain,read_line);
    getline(datain,read_line);

    //写入data格式的文件
    ofstream dataout(argument[1]);
    //文件头
    dataout << argument[2] <<endl;//写入description行
    dataout << "\n";
    dataout << num_atoms << " " << "atoms" <<endl;
    dataout << "\n";
    dataout << num_types <<" " << "atom types" <<endl;
    dataout << "\n";
    //插入#Masses#标记
    dataout << "#Masses#" <<endl;
    dataout << "\n";
    dataout << "Atoms" <<endl;
    dataout << "\n";
    //文件体
    double q = 0.0000; //默认电荷量
    int molecule = 1; //默认分子标记
    int tmp_loop = 1;
    while(getline(datain,read_line))
    {
        string tmp_type = read_line.substr(0,2); //获得当前行的type信息
        string tmp_position = read_line.substr(3); //获得当前行的原子坐标信息
        for (int i=0;i<num_types;i++) //集合解映射后，格式化写入
        {
            if (types_map[i] == tmp_type)
            {
                dataout << tmp_loop << " " << molecule << " " << i+1 << " " << q << " " << tmp_position <<endl;
                tmp_loop++;
                break;
            }
        }

    }

    //关闭流对象
    datain.close();
    dataout.close();
}