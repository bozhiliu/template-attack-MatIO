#include "TemplateAttack.hh"
#include "matio.h"
#include "eigen3/Eigen/Dense"
#include <iostream>
#include <math.h>



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


int main(int argc, char** argv){
  TrainRange(round){
    cout << "Round " << round << endl;
    const char * filename = "../test.mat";
    mat_t* mat = Mat_Open(filename, MAT_ACC_RDONLY);
    if (mat == NULL){
      cerr << "Error opening mat file" << endl;
      return 1;
    }

    matvar_t * Attack_varinfo = readCellInfo(mat, "z");
    display(Attack_varinfo);
    
    MatrixXd d1 = readCellData(mat, Attack_varinfo, 1,2);

    matvar_t * ainfo = readInfo(mat, "p");
    MatrixXd d2 = readData(mat, ainfo);

    cout << d1 << endl;
    cout << d2 << endl;
    
    for (int ite = 0; ite < OneRun; ite++){
      TemplateRange(testkey){
	



	
      }
    }
  }
  return 0;

}


