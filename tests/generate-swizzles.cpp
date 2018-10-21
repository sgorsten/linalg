#include <iostream>
#include <vector>
#include <algorithm>

bool unique(std::vector<int> n)
{
    std::sort(begin(n), end(n));
    return std::unique(begin(n), end(n)) == end(n);
}

int main()
{
    const char * const sets[] {"xyzw", "rgba", "stpq"};
    for(int n : {1,2,3,4})
    {
        std::cout << "        union\n        {\n";
        std::cout << "            detail::element_storage<T[" << n << "]> elems;\n";

        // Provide scalar accessors
        std::cout << "           ";
        for(int i=0; i<n; ++i)
        {
            std::cout << " _scalar<T,T[" << n << "]," << i << "> ";
            std::cout << sets[0][i] << ",";
            std::cout << sets[1][i] << ",";
            std::cout << sets[2][i] << ";";
        }
        std::cout << std::endl;

        // Provide two-element swizzles
        for(int i=0; i<n; ++i)
        {
            std::cout << "           ";
            for(int j=0; j<n; ++j)
            {
                std::cout << (unique({i,j}) ? " _lswizzle" : " _rswizzle");
                std::cout << "<T,T[" << n << "]," << i << "," << j << "> ";
                std::cout << sets[0][i] << sets[0][j] << ",";
                std::cout << sets[1][i] << sets[1][j] << ",";
                std::cout << sets[2][i] << sets[2][j] << ";";
            }
            std::cout << std::endl;
        }

        // Provide three-element swizzles
        for(int i=0; i<n; ++i)
        {
            for(int j=0; j<n; ++j)
            {
                std::cout << "           ";
                for(int k=0; k<n; ++k)
                {
                    std::cout << (unique({i,j,k}) ? " _lswizzle" : " _rswizzle");
                    std::cout << "<T,T[" << n << "]," << i << "," << j << "," << k << "> ";
                    std::cout << sets[0][i] << sets[0][j] << sets[0][k] << ",";
                    std::cout << sets[1][i] << sets[1][j] << sets[1][k] << ",";
                    std::cout << sets[2][i] << sets[2][j] << sets[2][k] << ";";
                }
                std::cout << std::endl;
            }
        }

        // Provide four-element swizzles
        for(int i=0; i<n; ++i)
        {
            for(int j=0; j<n; ++j)
            {
                for(int k=0; k<n; ++k)
                {
                    std::cout << "           ";
                    for(int l=0; l<n; ++l)
                    {
                        std::cout << (unique({i,j,k,l}) ? " _lswizzle" : " _rswizzle");
                        std::cout << "<T,T[" << n << "]," << i << "," << j << "," << k << "," << l << "> ";
                        std::cout << sets[0][i] << sets[0][j] << sets[0][k] << sets[0][l] << ",";
                        std::cout << sets[1][i] << sets[1][j] << sets[1][k] << sets[1][l] << ",";
                        std::cout << sets[2][i] << sets[2][j] << sets[2][k] << sets[2][l] << ";";
                    }
                    std::cout << std::endl;
                }
            }
        }

        std::cout << "        };\n\n";
    }
}