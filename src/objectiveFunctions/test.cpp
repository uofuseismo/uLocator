#include <iostream>
int main()
{
    int n = 6;
    int k = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            std::cout << k << " " << ( i*(2*n - i - 3) )/2 + j - 1 << std::endl;
            k = k + 1;
        }
getchar();
    }
}

