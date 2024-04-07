#include "Vector.hpp"
#include <vector>

int main()
{
    Vector<int> v1 = {10, 5, 6};
    Vector<int> v2 = {10, 5, 6};
    // Vector<int>v2;
    // v2=v1*9;
    // //v2.pop_back();
    // v2.push_back_many(4,5,6);

    // v2=v2.reverse();
    // v2.pop_back();
    // v2.print();
    // cout<<v2.accumulate();
    

    auto res=dot(v1,v2);
    cout<<res<<" "<<typeid(res).name();

}
