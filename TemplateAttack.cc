#include "TemplateAttack.hh"
#include "matio.h"
#include "eigen3/Eigen/Dense"
#include <iostream>
#include <string.h>
#include <math.h>
#include "boost/multi_array.hpp"
#include <ctime>



using namespace std;
using namespace Eigen;


double Exp(double x){
  return exp(x);
}


VectorXd mean(MatrixXd m, bool ColBased = true){
  if(ColBased){
    return m.colwise().mean();
  }
  else{
    return m.rowwise().mean();
  }  
}


MatrixXd cov(MatrixXd m, bool ColBased = true){
  return (m.rowwise() - m.colwise().mean()).transpose()*(m.rowwise()-m.colwise().mean()) / (m.rows()-1);  
}


VectorXd mvnpdf(MatrixXd data, VectorXd m, MatrixXd c){
  double denom = sqrt(c.determinant() * pow(2*M_PI, m.size()));
  VectorXd result(m.size());
  for (int i = 0; i < data.rows(); i++){
    
    result(i) = exp(-0.5*(data.row(i)-m.transpose()) * c.inverse() * (data.row(i)-m.transpose()).transpose()) / denom;
  }
  return result;
}


void test0(){
  MatrixXd m(5,2);
  m << 1,3,
    2,5,
    1.1, 3.2,
    1.2, 3.7,
    1.7, 4.2;
  MatrixXd data(2,2);
  data << 1.5 , 4,
    1.2, 3.3;
  cout << m << endl;
  cout << mean(m) << endl;
  cout << cov(m) << endl;
  cout << mvnpdf(data, mean(m), cov(m));
}



void display(matvar_t * element){
  cout << "Name: ";
  if (element->name != NULL){
    cout << element->name;
  }
  cout << " | ";
  cout << "Rank: " << element->rank << " | ";
  cout << "Dims: ";

  size_t* dims = element->dims;
  for (int i = 0; i < element->rank; i++){
    cout << dims[i] << " ";
  }
  cout << " | ";

  cout << "Class type: ";
  switch(element->class_type){
  case 0:
    cout << "empty";
    break;
  case 1:
    cout << "cell";
    break;
  case 6:
    cout << "double";
    break;
  default:
    cout << element->class_type;
    break;
  }
  cout << " | ";

  cout << "Data type: ";
  switch(element->data_type){
  case 0:
    cout << "unkown" ;
    break;
  case 9:
    cout << "double";
    break;
  case 21:
    cout << "cell";
    break;
  case 23:
    cout << "array";
    break;
  }
  cout << endl;  
}


void timedMessage(const char* mes, time_t * now){
  time(now);
  cout << mes << ctime(now);
}


matvar_t* readCellInfo(mat_t * mat, const char* name){
  matvar_t * varinfo = Mat_VarReadInfo(mat, name);
  if (varinfo == NULL){
    cerr << "Error reading " << name << " variable" << endl;    
  }
  if (varinfo->class_type != MAT_C_CELL){
    cerr << "Error reading " << name << ": not cell" << endl;
  }  
  return varinfo;
}


MatrixXd readCellData(mat_t* mat, matvar_t * varinfo, int row, int col){
  matvar_t * cell = Mat_VarGetCell(varinfo, col * (varinfo->dims)[0] + row);
  Mat_VarReadDataAll(mat, cell);
  MatrixXd data((cell->dims)[0], (cell->dims)[1]);
  double * cell_data = static_cast<double*>(cell->data);
  for(int i = 0; i < (cell->dims)[0]; i++){ // rows
    for(int j = 0; j < (cell->dims)[1]; j++){  // cols
      data(i,j) = cell_data[i+j*(cell->dims)[0]];
    }
  }
  return data;
}


matvar_t* readInfo(mat_t* mat, const char* name){
  matvar_t * varinfo = Mat_VarReadInfo(mat, name);
  if (varinfo == NULL){
    cerr << "Error reading " << name << " variable" << endl;
  }
  if (varinfo->class_type != MAT_C_DOUBLE){
    cerr << "Error reading " << name << ": not double" << endl;
  }
  return varinfo;  
}


MatrixXd readData(mat_t * mat, matvar_t* varinfo){
  matvar_t* cell = varinfo;
  Mat_VarReadDataAll(mat, cell);
  MatrixXd data((cell->dims)[0], (cell->dims)[1]);
  double * cell_data = static_cast<double*>(cell->data);
  for(int i = 0; i < (cell->dims)[0]; i++){ // rows
    for(int j = 0; j < (cell->dims)[1]; j++){  // cols
      data(i,j) = cell_data[i+j*(cell->dims)[0]];
    }
  }
  return data;
}


MatrixXd indexSlice(MatrixXd data, MatrixXd select){
  MatrixXd index;
  if(select.cols() == 1){
    index = select.transpose();
  }
  else{
    index = select.row(1);
  }
  MatrixXd result(data.rows(), index.cols());
  for(int i = 0; i<index.cols();i++){
    result.col(i) = data.col(index(0,i));
  }
  return result;
}



int main(int argc, char** argv){
  char * set = set6;
  time_t now = time(NULL);
  TrainRange(round){
    cout << "Round " << round << endl;
    char filename [100];
    memset(filename, 0,100);
    sprintf(filename, "%stemplate7_%s_28000_Set%d_average.mat",template_prefix, set, round+1);
    cout << "Opening template file " << filename ;
    mat_t* mat = Mat_Open(filename, MAT_ACC_RDONLY);
    if (mat == NULL){
      cerr << "Error opening mat file: " << filename << endl;
      return 1;
    }
    cout << " Done" << endl;

    cout << "Setting variable info and matrices..." << endl;
    
    timedMessage("Setting Attack ... ", &now);
    matvar_t * Attack_varinfo = readCellInfo(mat, "Attack");
    timedMessage("Creating array Attack ...", &now);
    boost::multi_array<MatrixXd, 1> Attack_matrix{boost::extents[_TemplateRange]};

    timedMessage("Setting SysMean ... ", &now);
    matvar_t * SysMean_varinfo = readCellInfo(mat, "SysMean");
    timedMessage("Setting SysCov ... ", &now);
    matvar_t * SysCov_varinfo = readCellInfo(mat, "SysCov");
    
  
    timedMessage("Creating array SysMean SysCov", &now);
    boost::multi_array<MatrixXd, 2> SysMean_matrix{boost::extents[OneRun][_TemplateRange]};
    boost::multi_array<MatrixXd, 2> SysCov_matrix{boost::extents[OneRun][_TemplateRange]};

    timedMessage("Setting CPUMean ... ", &now);    
    matvar_t * CPUMean_varinfo = readCellInfo(mat, "CPUMean");
    timedMessage("Setting CPUCov ... ", &now);
    matvar_t * CPUCov_varinfo = readCellInfo(mat, "CPUCov");
    boost::multi_array<MatrixXd, 2> CPUMean_matrix{boost::extents[OneRun][_TemplateRange]};
    boost::multi_array<MatrixXd, 2> CPUCov_matrix{boost::extents[OneRun][_TemplateRange]};

    
    timedMessage("Setting IcacheMean ... ", &now);
    matvar_t * IcacheMean_varinfo = readCellInfo(mat, "IcacheMean");
    timedMessage("Setting IcacheCov ... ", &now);
    matvar_t * IcacheCov_varinfo = readCellInfo(mat, "IcacheCov");
    boost::multi_array<MatrixXd, 2> IcacheMean_matrix{boost::extents[OneRun][_TemplateRange]};
    boost::multi_array<MatrixXd, 2> IcacheCov_matrix{boost::extents[OneRun][_TemplateRange]};

   
    timedMessage("Setting DcacheMean ... ", &now);    
    matvar_t * DcacheMean_varinfo = readCellInfo(mat, "DcacheMean");
    timedMessage("Setting DcacheCov ... ", &now);
    matvar_t * DcacheCov_varinfo = readCellInfo(mat, "DcacheCov");
    boost::multi_array<MatrixXd, 2> DcacheMean_matrix{boost::extents[OneRun][_TemplateRange]};
    boost::multi_array<MatrixXd, 2> DcacheCov_matrix{boost::extents[OneRun][_TemplateRange]};

    
    timedMessage("Setting L2cacheMean ... ", &now);
    matvar_t * L2cacheMean_varinfo;
    timedMessage("Setting L2cacheCov ... ", &now);
    matvar_t * L2cacheCov_varinfo;
    boost::multi_array<MatrixXd, 2> L2cacheMean_matrix{boost::extents[OneRun][_TemplateRange]};
    boost::multi_array<MatrixXd, 2> L2cacheCov_matrix{boost::extents[OneRun][_TemplateRange]};    
    if(strstr(set, "l2") != NULL){
      L2cacheMean_varinfo = readCellInfo(mat, "L2cacheMean");
      L2cacheCov_varinfo = readCellInfo(mat, "L2cacheCov");      
    }

		 
    timedMessage("Setting DRAMMean DRAMCov ... ", &now);
    matvar_t * DRAMMean_varinfo;
    matvar_t * DRAMCov_varinfo;
    boost::multi_array<MatrixXd, 2> DRAMMean_matrix{boost::extents[OneRun][_TemplateRange]};
    boost::multi_array<MatrixXd, 2> DRAMCov_matrix{boost::extents[OneRun][_TemplateRange]};    
    if(strstr(set, "dram") != NULL){
      DRAMMean_varinfo = readCellInfo(mat, "DRAMMean");
      DRAMCov_varinfo = readCellInfo(mat, "DRAMCov");  
    }

    
    timedMessage("Setting PS ... ", &now);
    matvar_t * PS_varinfo = readCellInfo(mat, "PS");
    boost::multi_array<MatrixXd, 2> PS_matrix{boost::extents[OneRun][_TemplateRange]};

    timedMessage("Setting Col ... ", &now);
    matvar_t * Col_varinfo = readCellInfo(mat, "Col");
    boost::multi_array<MatrixXd, 2> Col_matrix{boost::extents[OneRun][_TemplateRange]};

    timedMessage("Start reading templates...", &now);
    TemplateRange(key){
      timedMessage("Reading key template " , &now);
      Attack_matrix[key] = readCellData(mat, Attack_varinfo, 0, key);
      for(int ite = 0; ite < OneRun; ite++){	
	Col_matrix[ite][key] = readCellData(mat, Col_varinfo, ite, key);
	
	SysMean_matrix[ite][key] = readCellData(mat, SysMean_varinfo, ite, key);
	SysCov_matrix[ite][key] = readCellData(mat, SysCov_varinfo, ite, key);
	  
	CPUMean_matrix[ite][key] = readCellData(mat, CPUMean_varinfo, ite, key);
	IcacheMean_matrix[ite][key] = readCellData(mat, IcacheMean_varinfo, ite, key);
	DcacheMean_matrix[ite][key] = readCellData(mat, DcacheMean_varinfo, ite, key);

	CPUCov_matrix[ite][key] = readCellData(mat, CPUCov_varinfo, ite, key);
	IcacheCov_matrix[ite][key] = readCellData(mat, IcacheCov_varinfo, ite, key);
	DcacheCov_matrix[ite][key] = readCellData(mat, DcacheCov_varinfo, ite, key);

    	if(strstr(set, "l2") != NULL){
	  L2cacheMean_matrix[ite][key] = readCellData(mat, L2cacheMean_varinfo, ite, key);
	  L2cacheCov_matrix[ite][key] = readCellData(mat, L2cacheCov_varinfo, ite, key);
	}
	if(strstr(set, "dram") != NULL){
	  DRAMMean_matrix[ite][key] = readCellData(mat, DRAMMean_varinfo, ite, key);
	  DRAMCov_matrix[ite][key] = readCellData(mat, DRAMCov_varinfo, ite, key);  
	}
	
	matvar_t* CurrPS_cell = Mat_VarGetCell(PS_varinfo, key*OneRun+ite);
	matvar_t* CurrPS_keep = Mat_VarGetStructFieldByName(CurrPS_cell, "keep", 0);
	PS_matrix[ite][key] = readData(mat, CurrPS_keep);
      }
    }



    for (int ite = 0; ite < OneRun; ite++){      
      TemplateRange(testkey){
	MatrixXd CurrTest = indexSlice(indexSlice(Attack_matrix[testkey], PS_matrix[ite][testkey]), Col_matrix[ite][testkey]);	
	TemplateRange(key){
	MatrixXd CurrSysMean = SysMean_matrix[ite][key];
	MatrixXd CurrSysCov = SysCov_matrix[ite][key];

	MatrixXd CurrSynthMean = (CPUMean_matrix[ite][key] + IcacheMean_matrix[ite][key] + DcacheMean_matrix[ite][key]);
	MatrixXd CurrSynthCov = (CPUCov_matrix[ite][key] + IcacheCov_matrix[ite][key] + DcacheCov_matrix[ite][key]);
	if(strstr(set, "l2") != NULL){
	  CurrSynthMean = (CurrSynthMean + L2cacheMean_matrix[ite][key]);
	  CurrSynthCov = (CurrSynthCov + L2cacheCov_matrix[ite][key]);  
	}
	if(strstr(set, "dram") != NULL){
	  CurrSynthMean = (CurrSynthMean + DRAMMean_matrix[ite][key]);
	  CurrSynthCov = (CurrSynthCov + DRAMCov_matrix[ite][key]);  
	}
	

	}
      }
    }
  }
  return 0;

}


