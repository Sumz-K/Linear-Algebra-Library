#include<iostream>
#include"Matrix.hpp"

using namespace std;

int main(){
    int ele=4;
    Matrix<int,2,2> m1;
    for(int i=0;i<2;i++){
        for(int  j=0;j<2;j++){
            m1[i][j]=ele--;
        }

    }

    m1.print();

    Matrix<int,2,2> m2={{0,1},{1,0}};
    

    //m2=m2+2;
    int a=1;
    m2=m2+a;


    m2=m2*2;
    m2.print();

    cout<<m2.accumulate();

}
