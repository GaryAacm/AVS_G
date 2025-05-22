#ifndef AVS_API_MEMORY_H
#define AVS_API_MEMORY_H

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include<cstring>
#include <inttypes.h>
template<typename T>
class Memory
{
public:
    Memory()
    {
    }

    ~Memory()
    {
        mem_free();
    }

    void mem_free()
    {
        if(m_data)
        {
            free(m_data);
            m_data = nullptr;
            m_size = 0;
            m_capacity = 0;
        }
    }

    T* mem_malloc(uint32_t size)
    {
        m_capacity = size;
        m_data = (T*)malloc(m_capacity*sizeof(T));
        if(m_data == NULL)
        {
            printf("%s %d err\n", __FUNCTION__, __LINE__);
            exit(1);
        }
        return m_data;
    }

    T* mem_calloc(uint32_t size)
    {
        m_capacity = size;
        m_data = (T*)calloc(m_capacity, sizeof(T));
        if(m_data == NULL)
        {
            printf("%s %d err\n", __FUNCTION__, __LINE__);
            exit(1);
        }
        return m_data;
    }

    T* mem_realloc(uint32_t size)
    {
        m_capacity = size;
        m_data = (T*)realloc(m_data, m_capacity*sizeof(T));
        if(m_data == NULL)
        {
            printf("%s %d err\n", __FUNCTION__, __LINE__);
            exit(1);
        }
        return m_data;
    }

    void clear()
    {
        m_size = 0;
    }

    void push_back(T val)
    {
        if(m_capacity==m_size)
        {   
            m_capacity *= 2; 
            m_data = (T*)realloc(m_data, m_capacity*sizeof(T));
        }
        m_data[m_size++] = val;
    }

    void append(const T *ptr, uint32_t num)
    {
        if(m_capacity-m_size<=num)
        {   
            m_capacity += 2*num; 
            m_data = (T*)realloc(m_data, m_capacity*sizeof(T));
        }
        memcpy(m_data+m_size, ptr, num*sizeof(T));
        m_size += num;
    }

    void copy(T *ptr, uint32_t num)
    {
        if(m_capacity<num)
        {   
            m_capacity += 2*num; 
            m_data = (T*)realloc(m_data, m_capacity*sizeof(T));
        }
        memcpy(m_data, ptr, num*sizeof(T));
        m_size = num;
    }

    T &operator[](uint32_t idx)
    {
        return m_data[idx];
    }

    void swap(Memory<T> &arry)
    {
        T *tmp = this->m_data;
        this->m_data = arry.getdata();
        arry.setdata(tmp);
        uint32_t sizet = this->m_size;
        this->m_size = arry.size();
        arry.setsize(sizet);
        uint32_t capacityt = this->m_capacity;
        this->m_capacity = arry.capacity();
        arry.setcapacity(capacityt);
    }


    T *getdata(){return m_data;}
    T *getCurPtr()
    {
        return &m_data[m_size];
    }
    void updatesize(uint32_t size)
    {
        m_size += size;
    }
    uint32_t size(){return m_size;}
    uint32_t capacity(){return m_capacity;}
    void setdata(T *ptr)
    {
        m_data = ptr;
    }
    void setsize(uint32_t size)
    {
        m_size = size;
    }
    void setcapacity(uint32_t capacity)
    {
        m_capacity = capacity;
    }
private:
    T *m_data = nullptr;
    uint32_t m_size = 0; //实际个数,不是字节数
    uint32_t m_capacity = 0; //最大个数
};

#endif