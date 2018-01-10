#include "sorting.h"
#include <iostream>
#include <vector>
#include <cmath>

template<class T>
void printArray(T* arr, size_t n)
{
    for(size_t i = 0; i < n; ++i)
    {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
}

template<class T>
void printVector(const std::vector<T> &arr)
{
    for(const auto &e: arr)
        std::cout << e << " ";
    std::cout << std::endl;
}

int main()
{
    std::vector<double> a = {10, -99.4, 8, 7, 6, 5, 4, 3, 2, 100, 400, -900, 3.14, 15, 60};
    bubble_sort<double>(a, 0, a.size() - 1);
    std::cout << "Bubble sorting" << std::endl;
    printVector(a);

    //Insertion sorting

    std::vector<double> b = {2018.0, -9000.0, 180.4, -170.5, 6.4, 5000, 4.5, 3.7, 222.4, -12.3};
    insertion_sort<double>(b, 0, b.size() - 1);
    std::cout << "Insertion sorting 1" << std::endl;
    printVector(b);

    std::vector<double> b2 = {2017.0, -9700.0, 190.4, -160.5, 65.4, 500, 40.5, -3.7, 232.4, -1212.3};
    insertion_sort2<double>(b2, 0, b2.size() - 1);
    std::cout << "Insertion sorting 2" << std::endl;
    printVector(b2);

    //Shell sorting

    double c[] = {8,70,-17.4,6,5,4,100,3,2,1};

    shell_sort<double>(c, 10, [](const auto &a, const auto &b){ return a > b; });
    std::cout << "Shell sorting" << std::endl;
    printArray(c, 10);

    //Selection sorting

    double d[] = {1,2,3,4,5,6,7,8,9,10};

    selection_sort(d, 10, false);
    std::cout << "Selection sorting" << std::endl;
    printArray(d, 10);

    //Shaker sorting

    double e[] = {10, -99.4, 8, 7, 6, 5, 4, 3, 2, 100};

    shaker_sort(e, 10);
    std::cout << "Shaker sorting" << std::endl;
    printArray(e, 10);

    //Gnome sorting

    double g[] = {500, 287, 194, 125, 600, 400, -99, 185, 237, 537};
    auto compare = [](const auto &a, const auto &b){ return a > b; };
    gnome_sort<double>(g, 10, compare);
    std::cout << "Gnome sorting" << std::endl;
    printArray(g, 10);

    //Merge sorting

    std::cout << "Merge sorting" << std::endl;
    std::vector<int> vct;
    for(size_t i = 1; i <= 40; ++i)
        vct.push_back(i);

    srand(time(0));
    for(size_t i = 0; i < vct.size(); ++i)
        std::swap(vct[i], vct[rand() % (vct.size() - i) + i]);

    printVector(vct);

    merge_sort<int>(vct, 0, vct.size());

    printVector(vct);

    //Quick sorting
    std::cout << "Quick sorting" << std::endl;
    std::vector<double> vct2;
    for(size_t i = 1; i <= 19; ++i)
        vct2.push_back(pow(-1, i) * (i * i + 1));

    printVector(vct2);

    quick_sort_stack<double>(vct2, 0, vct2.size());

    printVector(vct2);

    std::cout << "Quick sorting recursively" << std::endl;
    std::vector<double> vct5;
    for(size_t i = 1; i <= 20; ++i)
        vct5.push_back(pow(-1, i) * (i + 1));

    printVector(vct5);

    quick_sort_recursively<double>(vct5, 0, vct5.size() - 1);

    printVector(vct5);

    std::cout << "Quick sorting recursively with 3 elements median" << std::endl;
    std::vector<double> vct8;
    for(size_t i = 1; i <= 100; ++i)
        vct8.push_back(pow(-1, i) * (i + 1));

    printVector(vct8);

    quick_sort_median3<double>(vct8, 0, vct8.size() - 1);

    printVector(vct8);

    //Heap sorting

    std::cout << "Heap sorting" << std::endl;
    std::vector<double> vct3;
    for(size_t i = 1, j = 7; i <= 8 && j > 0; ++i, --j)
        vct3.push_back(pow(-1, i + j) * (i * j * i) / 5 * 7);

    printVector(vct2);

    heap_sort<double>(vct3);

    printVector(vct3);

    std::cout << "******Order statistic**********" << std::endl;

    std::vector<double> vct4 {10,9,-8,700,6,5,40,3,2,7,200};

    printVector(vct4);

    for(int i = 1; i <= 10; ++i)
    {
        std::cout << "k" << i << " = " << find_order_statistic(vct4, 0, vct4.size() - 1, i) << std::endl;
    }

    heap_sort<double>(vct4);

    printVector(vct4);

    std::cout << "&vct[2] = " << &vct[2] << std::endl;

    return 0;
}
