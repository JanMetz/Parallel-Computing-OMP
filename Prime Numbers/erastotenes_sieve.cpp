#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <omp.h>

#include "erastotenes_sieve.hpp"

ErastotenesSieve::ErastotenesSieve(const int minNum, const int maxNum, const int threadsNum) : 
                                   mMinNum(minNum), mMaxNum(maxNum), mThreadsNum(threadsNum)
{

}

void ErastotenesSieve::printPrimes(const std::vector<int> &primes)
{
    std::cout << "Found " << primes.size() << " primes: " << std::endl;

    int i = 0;
    for (const auto &prime : primes)
    {
        std::cout << prime << '\t';
        if (i % 10 == 9) 
            std::cout << std::endl;
        ++i;
    }
}

std::vector<int> ErastotenesSieve::findPrimesSequential_add(const int max)  const
{
    std::vector<bool> isPrime(max - 1, true);

    for (int divider = 2; divider <= max; ++divider)
    {
        if (isPrime[divider - 2] == false)
            continue;

        removeMultiples(divider, isPrime, 2);
    }

    std::vector<int> startingPrimes;
    fillPrimesList(isPrime, startingPrimes, 2);

    return startingPrimes;
}

bool ErastotenesSieve::isPrime_div(const int number)
{
    const int limit = floor(sqrt(number));
    int potential_div = 2;

    while (potential_div <= limit) 
    {
        if (number % potential_div == 0) 
            return false;

        ++potential_div;
    }

    return true;
}

std::vector<int> ErastotenesSieve::findPrimesSequential_div() const
{
    std::vector<int> primes;
    for (int i = mMinNum; i < mMaxNum; ++i)
    {
        if (ErastotenesSieve::isPrime_div(i))
            primes.emplace_back(i);
    }

    return primes;
}

std::vector<int> ErastotenesSieve::findPrimesDomain() const
{
    const std::vector<std::vector<int>> ranges{getRanges()};
    const std::vector<int> startingPrimes{findPrimesSequential_add(static_cast<int>(sqrt(mMaxNum)))};

    std::vector<std::vector<int>> primesMulti(mThreadsNum);

    #pragma omp parallel firstprivate(ranges, startingPrimes) shared(primesMulti) 
    {
        const int threadID = omp_get_thread_num();

        const int lowerLimit = ranges[threadID][0];
        const int upperLimit = ranges[threadID][1];
        
        std::vector<bool> isPrime((upperLimit - lowerLimit + 1), true);
        for (const auto &prime : startingPrimes)
        {
            int i = 0;
            while (prime * i < lowerLimit)
                ++i;
        
            const int offset = prime * i - lowerLimit;
            for (int multiple = 0; multiple + offset < isPrime.size(); multiple += prime)
                isPrime[offset + multiple] = false;
        }

        std::vector<int> primes;
        fillPrimesList(isPrime, primes, lowerLimit);
	
        primesMulti[threadID] = primes;
    }

    std::vector<int> primes;
    for (const auto &prime : startingPrimes)
        if (prime >= mMinNum)
            primes.emplace_back(prime);

    for (int i = 0; i < mThreadsNum; ++i)
	    primes.insert(primes.end(), primesMulti[i].begin(), primesMulti[i].end());

    return primes;
}

std::vector<int> ErastotenesSieve::findPrimesDiv() const
{
    std::vector<std::vector<bool>> isPrimeMulti(mThreadsNum);
    for (auto &row : isPrimeMulti)
        row = std::vector<bool>(mMaxNum - mMinNum + 1, true);

    #pragma omp parallel for schedule(guided) shared(isPrimeMulti)
    for (int num = mMinNum; (num - mMinNum) < isPrimeMulti[omp_get_thread_num()].size(); ++num)
        isPrimeMulti[omp_get_thread_num()][num - mMinNum] = ErastotenesSieve::isPrime_div(num);

    std::vector<bool> isPrime;
    combinePrimesLists(isPrimeMulti, isPrime);

    std::vector<int> primes;
    fillPrimesList(isPrime, primes, mMinNum);

    return primes;
}

std::vector<std::vector<int>> ErastotenesSieve::getRanges() const
{
    const int step = static_cast<int>((mMaxNum - mMinNum) / mThreadsNum);
    std::vector<std::vector<int>> ranges;

    int nextNumber = mMinNum;
    for (int i = 0; i < mThreadsNum - 1; ++i)
    {
        ranges.push_back({ nextNumber, nextNumber + step - 1});
        nextNumber += step;
    }

    ranges.push_back({ nextNumber, mMaxNum });

    return ranges;
}

std::vector<int> ErastotenesSieve::findPrimesFunctional() const
{
    const std::vector<int> startingPrimes{findPrimesSequential_add(static_cast<int>(sqrt(mMaxNum)))};

    std::vector<std::vector<bool>> isPrimeMulti(mThreadsNum);
    for (auto &row : isPrimeMulti)
        row = std::vector<bool>(mMaxNum - mMinNum + 1, true);

    #pragma omp parallel for schedule(guided) firstprivate(startingPrimes) shared(isPrimeMulti)
    for (int i = 0; i < startingPrimes.size(); ++i)
        removeMultiples(startingPrimes[i], isPrimeMulti[omp_get_thread_num()], mMinNum);

    std::vector<bool> isPrime;
    combinePrimesLists(isPrimeMulti, isPrime);

    std::vector<int> primes;
    fillPrimesList(isPrime, primes, mMinNum);

    return primes;
}

void ErastotenesSieve::fillPrimesList(const std::vector<bool> &isPrime, std::vector<int> &primes, const int add) const
{
    for (int i = 0; i < isPrime.size(); ++i)
    {
        if (isPrime[i])
            primes.emplace_back(i + add);
    }
}

void ErastotenesSieve::removeMultiples(const int prime, std::vector<bool> &isPrime, const int sub) const
{
    for (int multiple = prime * 2; multiple - sub < isPrime.size(); multiple += prime)
        isPrime[multiple - sub] = false;
}

void ErastotenesSieve::combinePrimesLists(const std::vector<std::vector<bool>> &threadsPrimesLists, std::vector<bool> &isPrime) const
{
    int i = 0;
    #pragma omp parallel for schedule(guided) shared(isPrime) firstprivate(mThreadsNum, range1, threadsPrimesList)
    for (; i < threadsPrimesLists[0].size(); ++i)
    {
        bool combined = true;
        for (int j = 0; j < mThreadsNum; ++j)
            combined *= threadsPrimesLists[j][i];

        isPrime.emplace_back(combined);
    }
    
    for (; i < threadsPrimesLists.back().size(); ++i)
        isPrime.emplace_back(threadsPrimesLists.back()[i]);
}