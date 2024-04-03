#include <iostream>
#include "Vector_scratch.h"

template <typename T, int Rows, int Cols>
class Matrix {
private:
    int rows;
    int cols;
    Vector_scratch<T> mat;

public:
    Matrix() : rows(Rows), cols(Cols), mat(Rows * Cols) {

    }

    T& operator()(int i, int j) {
        if (i >= rows || j >= cols) {
            throw std::out_of_range("Index out of bounds");
        }
        return mat[i * cols + j];
    }



    Vector_scratch<T> operator[](int index){
        if(index>=rows){
            throw std::out_of_range("Index out of bounds");
        }
        int start=index*cols;
        int end=start+cols;

        Vector_scratch<T> res(cols);

        for(int i=start;i<end;i++){
            res.push_back(mat[i]);
        }
        return res;

    }

    int numRows() const {
        return rows;
    }

    int numCols() const {
        return cols;
    }

    int size(){
        return rows*cols;
    }

    void put(const T& value, int i, int j) {
        if (i >= rows || j >= cols || i < 0 || j < 0) {
            throw std::out_of_range("Index out of bounds");
        }
        mat[i * cols + j] = value;
    }

    void print(){
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                cout<<(*this)(i,j)<<" ";
            }
            cout<<"\n";
        }
    }

    
        
};

