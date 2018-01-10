#ifndef SORTING_H
#define SORTING_H

#include <vector>
#include <stack>
#include <algorithm>
#include <functional>
#include <iostream>

/*
 * 1.Быстрая сортировка с поиском медианы трех элементов
 * 2.Быстрая сортировка с разделением на три части
*/

#define MAXSTACK 2048
#define THRESHOLD 11

template<class T>
long partition(std::vector<T> &arr, long left_bound, long right_bound)
{
    T pivot = arr[(left_bound + right_bound) >> 1];
    auto i = left_bound;
    auto j = right_bound;

    while(i < j)
    {
        while(arr[i] < pivot) ++i;
        while(arr[j] > pivot) --j;
        if(i < j) std::swap(arr[i++], arr[j--]);
    }
    return j;
}

template<class T>
T find_order_statistic(std::vector<T> &arr, long begin, long end, long order)
{
    long left = begin;
    long right = end;
    for(;;)
    {
        auto mid = partition(arr, left, right);
        if(mid == order - 1)
            return arr[mid];
        else if(order - 1 < mid)
            right = mid;
        else
            left = mid + 1;
    }
}

template<class T>
void bubble_sort(std::vector<T> &arr, long begin, long end,
                 std::function<int(const T &a, const T &b)> comparator =
                 [](const auto &a, const auto &b){
                    return a == b ? 0 : (a > b ? 1 : -1);
                 })
{
    bool isSorted = false;
    while (!isSorted) {
        isSorted = true;
        int j = begin;
        for(long i = begin; i < end - j; ++i)
        {
            if(comparator(arr[i], arr[i + 1]) > 0)
            {
                std::swap(arr[i], arr[i + 1]);
                isSorted = false;
            }
        }
        ++j;
    }
}

/*
 * Алгоритм работает одинаково в лучшем и худшем случаях Teta(N^2)
*/

template<class T>
void insertion_sort(std::vector<T> &arr, long begin, long end,
                    std::function<int(const T &a, const T &b)> comparator =
                    [](const auto &a, const auto &b){
                       return a == b ? 0 : (a > b ? 1 : -1);
                    })
{
    for(long i = end; i > begin; --i)
    {
        if(comparator(arr[i], arr[i - 1]) < 0)
            std::swap(arr[i], arr[i - 1]);
    }

    for(long i = begin + 2; i <= end; ++i)
    {
        int j = i;
        T tmp = arr[j];
        while(j > 0 && comparator(tmp, arr[j - 1]) < 0)
        {
            arr[j] = arr[j - 1];
            --j;
        }
        arr[j] = tmp;
    }
}

template<class T>
void insertion_sort2(std::vector<T> &arr, long begin, long end,
                     std::function<int(const T &a, const T &b)> comparator =
                     [](const auto &a, const auto &b){
                        return a == b ? 0 : (a > b ? 1 : -1);
                     })
{
    for(long i = begin + 1; i <= end; ++i)
    {
        T tmp = arr[i];
        size_t j = i;
        for(;j > 0 && comparator(tmp, arr[j - 1]) < 0; --j)
        {
            arr[j] = arr[j - 1];
        }
        arr[j] = tmp;
    }
}

template<class T>
void shell_sort(T* arr, size_t size,
                std::function<bool(const T &a, const T &b)> compare =
                [](const auto &a, const auto &b){ return a < b; } )
{
    if(!arr) return;
    size_t h = 1;
    while(h <= size / 9) h = 3*h + 1;
    std::cout << "h = " << h << std::endl;
    for(; h > 0; h /= 3)
    {
        for(size_t i = h; i < size; ++i)
        {
            size_t j = i;
            T tmp = arr[i];
            while(j >= h && compare(tmp, arr[j - h]))
            {
                arr[j] = arr[j - h];
                j -= h;
            }
            arr[j] = tmp;
        }
    }
}

/*
 * 1.Сложность во всех случаях O(N^2)
 * 2.Алгоритм неустойчив
 * 3.Минимальное количество обменов
*/

template<class T>
void selection_sort(T* arr, size_t size, bool isAscendingOrder = true)
{
    if(!arr) return;
    for(size_t i = 0; i < size - 1; ++i)
    {
        size_t desiredIndex;
        if(isAscendingOrder)
        {
            size_t min_index = i;
            for(size_t j = min_index + 1; j < size; ++j)
            {
                if(arr[j] < arr[min_index]) min_index = j;
            }
            desiredIndex = min_index;
        }
        else
        {
            int max_index = i;
            for(size_t j = max_index + 1; j < size; ++j)
            {
                if(arr[j] > arr[max_index]) max_index = j;
            }
            desiredIndex = max_index;
        }

        if(desiredIndex != i) std::swap(arr[i], arr[desiredIndex]);
    }
}

template<class T>
void shaker_sort(T arr[], size_t size, bool isAscendingOrder = true)
{
    for(size_t left = 0, right = size - 1; left < right;)
    {
        for(size_t i = left; i < right - 1; ++i)
        {
            if(isAscendingOrder ? arr[i] > arr[i + 1] : arr[i] < arr[i + 1])
                std::swap(arr[i], arr[i + 1]);
        }
        if(isAscendingOrder)
            --right;
        else
            ++left;
        for(size_t i = right; i > left; --i)
        {
            if(isAscendingOrder ? arr[i] < arr[i - 1] : arr[i] > arr[i - 1])
                std::swap(arr[i], arr[i - 1]);
        }
        if(isAscendingOrder)
            ++left;
        else
            --right;
    }
}

template <class T>
void gnome_sort(T arr[], size_t size,
                std::function<bool(const T &a, const T &b)> compare =
                [](const auto &a, const auto &b){ return a < b; })
{
    size_t i {0};
    while (i < size) {
        if(i == 0 || compare(arr[i - 1], arr[i]))
        {
            ++i;
        }
        else
        {
            std::swap(arr[i], arr[i - 1]);
            --i;
        }
    }
}

template<class T>
void merge_sort(std::vector<T> &a, long begin, long end,
                std::function<bool(const T &a, const T &b)> compare =
                [](const auto &a, const auto &b){ return a < b; })
{
    if(end - begin < 2) return;
    if(end - begin == 2)
    {
        if(!compare(a[begin], a[begin + 1]))
            std::swap(a[begin], a[begin + 1]);
        return;
    }
    auto mid = (begin + end) >> 1;

    merge_sort(a, begin, mid, compare);
    merge_sort(a, mid, end, compare);

    std::vector<T> res;
    long pos_left = begin;
    long pos_right = mid;

    while (pos_left < mid && pos_right < end) {
        if(compare(a[pos_left], a[pos_right]))
            res.push_back(a[pos_left++]);
        else
            res.push_back(a[pos_right++]);
    }

    while(pos_right < end) res.push_back(a[pos_right++]);
    while(pos_left < mid) res.push_back(a[pos_left++]);

    for(auto i = begin; i < end; ++i)
        a[i] = res[i - begin];
}

template<class T>
long find_median_index(const std::vector<T> &arr, long left, long mid, long right)
{
    std::vector<T> a = {arr[left], arr[mid], arr[left]};
    bubble_sort(a, 0, a.size() - 1);
    if(a[1] == arr[left])
        return left;
    else if (a[1] == arr[mid])
        return mid;
    else
        return right;
}

template<class T>
void quick_sort_median3(std::vector<T> &arr, long begin, long end,
                        std::function<int(const T &a, const T &b)> compare =
                        [](const auto &a, const auto &b){
                                return a == b ? 0 : (a > b ? +1 : -1);
                        })
{
    if(end - begin <= THRESHOLD)
    {
        insertion_sort2(arr,begin,end,compare);
        return;
    }
    auto middle_pos = (begin + end) >> 1;
    auto median_pos = find_median_index(arr,begin,middle_pos, end);
    std::swap(arr[median_pos],arr[middle_pos]);
    //auto mid = partition(arr, left, right);
    auto mid = partition(arr,begin,end);
    quick_sort_median3(arr,begin,mid,compare);
    quick_sort_median3(arr,mid + 1,end,compare);
}

template<class T>
void quick_sort_recursively(std::vector<T> &arr, long begin, long end,
                            std::function<int(const T &a, const T &b)> compare =
                            [](const auto &a, const auto &b){
                                    return a == b ? 0 : (a > b ? +1 : -1);
                            })
{
    long i = begin;
    long j = end;
    auto pivot_pos = (begin + end) >> 1;
    T pivot = arr[pivot_pos];
    while(i <= j)
    {
        while(compare(arr[i], pivot) < 0) ++i;
        while(compare(arr[j], pivot) > 0) --j;
        if(i <= j) std::swap(arr[i++], arr[j--]);
    }
    if(j > begin) quick_sort_recursively(arr, begin, j);
    if(i < end) quick_sort_recursively(arr, i, end);
}

template<class T>
void quick_sort_stack(std::vector<T> &arr, size_t begin, size_t end,
                std::function<int(const T &a, const T &b)> compare =
                [](const auto &a, const auto &b){
                        return a == b ? 0 : (a > b ? +1 : -1);
                })
{
    if(begin < 0 || begin >= arr.size() || end <= 0 || end > arr.size())
        return;
    long i, j, left_bound = begin, right_bound = end - 1;
    std::stack<std::pair<long, long>> st;
    T pivot;
    long pivot_pos;
    st.push(std::make_pair(left_bound, right_bound));
    while(!st.empty())
    {
        left_bound = st.top().first;
        right_bound = st.top().second;
        st.pop();
        while(left_bound < right_bound)
        {
            pivot_pos = (left_bound + right_bound) >> 1;
            i = left_bound;
            j = right_bound;
            pivot = arr[pivot_pos];

            while(i <= j)
            {
                 while (compare(arr[i], pivot) < 0) ++i;
                 while (compare(arr[j], pivot) > 0) --j;
                 if(i <= j) std::swap(arr[i++], arr[j--]);
            }

            if(i < pivot_pos) //Right-hand side is greater
            {
                if(i < right_bound)
                    st.push(std::make_pair(i, right_bound));
                right_bound = j;
            }
            else //Left-hand side is greater
            {
                if(j > left_bound)
                    st.push(std::make_pair(left_bound, j));
                left_bound = i;
            }
        }
    }
}

template<class T>
void quick_sort_arraystack(std::vector<T> &arr, size_t begin, size_t end,
                std::function<int(const T &a, const T &b)> compare =
                [](const auto &a, const auto &b){
                        return a == b ? 0 : (a < b ? +1 : -1);
                })
{
    if(begin < 0 || begin >= arr.size() || end <= 0 || end > arr.size())
        return;
    long i, j;
    long left_bound, right_bound;
    long left_bound_stack[MAXSTACK], right_bound_stack[MAXSTACK];
    long stack_pos {1};
    T pivot;
    long pivot_pos;

    left_bound_stack[1] = begin;
    right_bound_stack[1] = end - 1;
    while(stack_pos != 0)
    {
        left_bound = left_bound_stack[stack_pos];
        right_bound = right_bound_stack[stack_pos];
        --stack_pos;

        while(left_bound < right_bound)
        {
             pivot_pos = (left_bound + right_bound) >> 1;
             i = left_bound;
             j = right_bound;
             pivot = arr[pivot_pos];

             while(i <= j)
             {
                  while (compare(arr[i], pivot) > 0) ++i;
                  while (compare(arr[j], pivot) < 0) --j;
                  if(i <= j) std::swap(arr[i++], arr[j--]);
             }

             if(i < pivot_pos) //Right-hand side is greater
             {
                 if(i < right_bound)
                 {
                     ++stack_pos;
                     left_bound_stack[stack_pos] = i;
                     right_bound_stack[stack_pos] = right_bound;
                 }
                 right_bound = j;
             }
             else //Left-hand side is greater
             {
                 if(j > left_bound)
                 {
                     ++stack_pos;
                     left_bound_stack[stack_pos] = left_bound;
                     right_bound_stack[stack_pos] = j;
                 }
                 left_bound = i;
             }
        }
    }
}

template<class T>
void heapify(std::vector<T> &arr, long k, long n)
{
    if(k < 0 || k >= long(arr.size()) || n < 0 || n >= long(arr.size()) || k >= n) return;
    T next = arr[k];
    long child_index;
    while(k <= n >> 1)
    {
        child_index = 2 * k;
        // We choose greater child
        if(child_index < n && arr[child_index] < arr[child_index + 1])
            ++child_index;
        if(next < arr[child_index])
        {
            std::swap(arr[k], arr[child_index]);
            k = child_index;
        }
        else
            break;
    }
    arr[k] = next;
}

template<class T>
void heap_sort(std::vector<T> &arr)
{
    //Here we create binary heap
    for(long i = (arr.size() >> 1) - 1; i >=0; --i)
        heapify(arr, i, arr.size() - 1);

    for(long i = arr.size() - 1; i > 0; --i)
    {
        std::swap(arr[0], arr[i]);
        heapify(arr, 0, i - 1);
    }
}



#endif // SORTING_H
