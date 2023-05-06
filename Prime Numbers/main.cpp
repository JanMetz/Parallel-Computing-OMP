#include <iostream>
#include <vector>
#include <string>
#include <omp.h>

#include "erastotenes_sieve.hpp"

void printIf(const std::vector<int> &result, const short dispOpt)
{
    switch(dispOpt)
    {
        case 0:                                                                     break;
        case 1: std::cout << "Found " << result.size() << " primes: " << std::endl; break;
        case 2: ErastotenesSieve::printPrimes(result);                              break;
        default: std::cout << "Error:  Invalid display option selected!" << std::endl;
    }

    std::cout << "Process finished!" << std::endl; 
}

int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        std::cout << "Wrong number of arguments provided! Invoke template: "
                  << ".\\app.exe minPrimeNumber maxPrimeNumber calculationsMode displayMode\n"
                  << "Possible calculation modes: \"seq_a\" for sequential additive, \"seq_d\" for sequential divisive, "
                  << "dom\" for domain, \"fun\" for functional, \"div\" for parallel divisive.\n"
                  << "Possible display modes: \"0\" for not displaying anything, \"1\" for displaying found primes quantity, "
                  << "2\" for displaying all found prime numbers." << std::endl;

        return 1;
    }

    const long long minNum = std::stoll(argv[1]);
    const long long maxNum = std::stoll(argv[2]);
    const short dispOpt = std::stoi(argv[3]);
    const std::string calcOpt = argv[4];
    
    if (maxNum < minNum)
    {
        std::cout << "Error: Minimal number is bigger than maximal!" << std::endl;
        return 1;
    }
    
    if (minNum < 2)
    {
        std::cout << "Error: Minimal number is smaller than 2!" << std::endl;
        return 1;
    }

    if (maxNum - minNum + 1 < omp_get_max_threads())
    {
        std::cout << "Error: Span is smaller than number of threads!" << std::endl;
        return 1;
    }

    ErastotenesSieve es(minNum, maxNum, omp_get_max_threads());

    if (calcOpt == "seq_a")
    {
        std::vector<int> result{es.findPrimesSequential_add(maxNum)};
        printIf(result, dispOpt);
    }

    else if (calcOpt == "seq_d")
    {
        std::vector<int> result{es.findPrimesSequential_div()};
        printIf(result, dispOpt);
    }

    else if (calcOpt == "div")
    {
        std::vector<int> result{es.findPrimesDiv()};
        printIf(result, dispOpt);
    }

    else if (calcOpt == "fun")
    {
        std::vector<int> result{es.findPrimesFunctional()}; 
        printIf(result, dispOpt);
    }

    else if (calcOpt == "dom")
    {
        std::vector<int> result{es.findPrimesDomain()}; 
        printIf(result, dispOpt);
    }
    
    else
    {
        std::cout << "Error: Invalid calculations option selected!" << std::endl;
        return 1;
    }

    return 0;
}
