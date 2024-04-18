
#include"Matrix.hpp"


int main()
{
    std::cout<<"\t\t########   Vector #######\n";
    Vector<int> vector1 = {1, 2, 3, 4, 5};
    cout<<"accumulate\n"<<vector1.accumulate<5>()<<"\n";
    
    std::cout << "Vector initialized with size 5: ";
    vector1.print(); 
    std::cout<<"Normalisation:\n";
    auto r=normalise(vector1);
    r.print();
    std::cout<<"\n";
    Vector<int> vector2(vector1);
    std::cout << "Copied vector v2: ";
    vector2.print();
    std::cout<<"\n";
    vector1.push_back(6);
    vector1.push_back(7);
    std::cout << "Vector after pushing back elements 6,7: ";
    vector1.print();
    std::cout<<"\n\n";
    vector1.pop_back();
    std::cout << "Vector after popping back: ";
    vector1.print(); 
    std::cout<<"\n\n";
    vector1.push_back_many(7,8,9,10,11);
    std::cout << "Vector after pushing back many elements: ";
    vector1.print(); 
    std::cout<<"\n\n";
    std::cout << "Element at index 3 in vector1: " << vector1[3] << "\n\n";
    std::cout << "Size of vector2: " << vector1.size() <<"\n\n"; 
    std::cout << "Max size of vector2: " << vector1.maxsize() <<"\n\n"; 
    vector1.del();
    std::cout << "Vector1 after deletion: ";
    vector1.print(); 
    std::cout<<"\n";
    std::cout << "Reversed vector1: ";

    Vector<int> v1 = {1,2,3,4,5,6,7,8,9,10,11};

    v1.reverse().print(); 
    std::cout<<"\n\n";
    std::cout<<"\n\n";
    std::cout << "Magnitude of vector1: " << v1.magnitude() <<"\n\n";
    Vector<int> vector3 = {1, 2, 3};
    //cout<<vector3[0];
    Vector<int> vector4 = {4, 5, 6};
    std::cout<<"Vector<int> vector3 = {1, 2, 3};\nVector<int> vector4 = {4, 5, 6};\n\n";
    std::cout << "Dot product of vector3 and vector4: " << dot(vector3, vector4) <<"\n\n";
    std::cout << "Euclidean distance between vector3 and vector4: " << euclid_distance(vector3, vector4) <<"\n\n";
    std::cout << "Manhattan distance between vector3 and vector4: " << manhattan_distance(vector3, vector4) <<"\n\n";
    vector3.sort();
    std::cout<<"sorting vector3 {1, 2, 3}\n";
    std::cout << "Sorted vector3: ";
    vector3.print();
    std::cout<<"\n\n";
    std::cout<<"\t\t########   Matrix #######\n";
    Matrix<int, 3, 3> matrix1;
    std::cout << "Matrix initialized with zeros:\n";
    matrix1.print();
    Matrix<int, 3, 3> matrix2 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    std::cout << "Matrix initialized using initializer list:\n";
    matrix2.print();  
    std::cout<<"\n";
    std::cout << "Element at row 1, column 2: " << matrix2[1][2] << std::endl; // Should print 6
    std::cout<<"\n";
    // Add two matrices
    Matrix<int, 3, 3> matrix3 = matrix1 + 2;
    std::cout << "Matrix after addition:\n";
    matrix3.print(); 
    Matrix<int, 3, 3> matrix5 = matrix2 * 2;
    std::cout << "Matrix after multiplication by scalar:\n";
    matrix5.print();
    Matrix<int, 3, 3> matrix6 = matrix2.transpose();
    std::cout << "Transposed matrix:\n";
    matrix6.print();
    
    Matrix<int, 2, 2> matrix7 = matrix2.slice<2, 2>(0, 2, 0, 2);
    std::cout << "Sliced matrix:0-0 to 2-2\n";
    matrix7.print();
    matrix2.sort();
    std::cout << "Sorted matrix:\n";
    matrix2.print();
    Matrix<int, 2, 3> mat1 = {{1, 2, 3}, {4, 5, 6}};
    Matrix<int, 3, 2> mat2 = {{7, 8}, {9, 10}, {11, 12}};
    Matrix<int, 2, 2> mat_result = multiply(mat1, mat2);
    std::cout<<"Multipling matrices :\nMatrix<int, 2, 3> mat1 = {{1, 2, 3}, {4, 5, 6}}\nMatrix<int, 3, 2> mat2 = {{7, 8}, {9, 10}, {11, 12}}\n\n";
    std::cout << "Result of matrix multiplication:\n";
    mat_result.print();
    Matrix<int, 2, 2> matdet = {{1, 2}, {3, 4}};
    std::cout << "Determinant of the matrix: \n{1, 2}\n{3, 4}\n" << determinant(matdet) << std::endl;
    Matrix<int, 2, 2> mat_det = {{2, 0}, {1, 3}};
    std::cout<<"eigen values of this matrix\n{2, 0}\n{1, 3}\n";
    int a = 1;
    int b = -(mat_det[0][0] + mat_det[1][1]);
    int c = mat_det[0][0] * mat_det[1][1] - mat_det[0][1] * mat_det[1][0];

    auto eigen2=solveQuadraticEquation(a,b,c);
    eigen2.print();
    Matrix<int, 3, 3> myMatrix;
    int values[9]={6,2,3,3,1,1,10,3,4};
        for (int i = 0; i < myMatrix.numRows(); ++i) {
            for (int j = 0; j < myMatrix.numCols(); ++j) {
            
                myMatrix[i][j] = values[i*3+j];
        }
    }
    std::cout<<"3*3 matrix\n";
    myMatrix.print();
    Matrix<int, 3,3> adj=adjoint(myMatrix);
    std::cout<<"adjoint of matrix is :\n";
    adj.print();
    std::cout << "Inverse of the matrix:\n";
    auto inv=inverse(myMatrix);


    cout<<inv<<"\n";
    inv.print();


    Matrix<int,4,4> m99={{1,4,1,6},{-1,1,2,10},{3,6,1,8},{4,1,1,1}};
    cout<<determinant(m99);

}
