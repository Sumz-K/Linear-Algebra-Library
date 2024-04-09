#include <iostream>
#include <type_traits>
#define INIT 1
#define MAX 100
using namespace std;

template <typename T>
class Vector_scratch
{
public:
    T *container;
    int capacity;
    int current;


    T* getContainer(){
        return this->container;
    }

    int getCapacity(){
        return this->capacity;
    }
    int getCurrent(){
        return this->current;
    }

    Vector_scratch(int size)
    {
        capacity = size;
        current = 0;
        container = new T[capacity];

        for (int i = 0; i < capacity; i++)
        {
            container[i] = 0;
        }
    }

    Vector_scratch()
    {
        capacity = INIT;
        current = 0;
        container = new T[capacity];
        for (int i = 0; i < capacity; i++)
        {
            container[i] = 0;
        }
    }

    Vector_scratch(initializer_list<T> list)
    {
        capacity = list.size();
        container = new T[capacity];
        int i = 0;
        for (auto &elem : list)
        {
            container[i++] = elem;
        }
        current = capacity;
    }

    Vector_scratch(const Vector_scratch<T>& other) : capacity(other.capacity), current(other.current) {
        container = new T[capacity];
        for (int i = 0; i < current; ++i) {
            container[i] = other.container[i];
        }
    }


    void print()
    {
        if(this->capacity==0){
            cout<<"No memory allocated to the vector\n";
            return ;
        }
        for (int i = 0; i < current; i++)
        {
            cout << container[i] << " ";
        }
        cout << "\n";
    }

    void push_back(T data)
    {
        if (current < capacity)
        {
            container[current] = data;
            current += 1;
        }
        else
        {
            T *expanded = new T[2 * capacity];
            for (int i = 0; i < capacity; i++)
            {
                expanded[i] = container[i];
            }
            delete[] container;
            capacity = capacity * 2;
            container = expanded;
            container[current] = data;
            current += 1;
        }
    }

    void pop_back()
    {
        current -= 1;
    }

    void pop_front()
    {
        for (int i = 0; i < current; i++)
        {
            container[i] = container[i + 1];
        }
        current -= 1;
    }

    T &operator[](int index)
    {
        if (index >= capacity || index < 0)
        {
            int newSize = (index >= capacity) ? index + 1 : capacity * 2;
            T *expanded = new T[newSize];
            for (int i = 0; i < newSize; i++)
            {
                expanded[i] = (i < capacity) ? container[i] : 0;
            }
            delete[] container;
            capacity = newSize;
            container = expanded;
        }
        if (index >= current)
        {
            current = index + 1;
        }
        return container[index];
    }
    
    const T& operator[](int index) const{
        if (index < 0 || index >= current) {
            throw std::out_of_range("Index out of range");
        }
        return container[index];
    }
    int size()
    {
        return current;
    }

    int maxsize()
    {
        return capacity;
    }

    template <typename... Args>
    void push_back_many(Args &&...args)
    {
        (this->push_back(std::forward<Args>(args)), ...);
    }

    T accumulate() const
    {
        return accumulate_helper(std::make_index_sequence<MAX>{});
    }

    template <size_t... Is>
    T accumulate_helper(std::index_sequence<Is...>) const
    {
        return (container[Is] + ...);
    }
    void del()
    {
        

        if constexpr (std::is_class_v<T>)
        {
            for (int i = 0; i < current; ++i)
            {
                container[i].~T();
            }
        }
        this->capacity=0;
        delete[] container; 
    }


};
