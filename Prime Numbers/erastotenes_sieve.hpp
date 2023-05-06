#include <vector>

class ErastotenesSieve
{
public:
    ErastotenesSieve(const int minNum, const int maxNum, const int threadsNum);

    static void printPrimes(const std::vector<int> &primes);
    std::vector<int> findPrimesSequential_add(const int max) const;
    std::vector<int> findPrimesSequential_div() const;
    std::vector<int> findPrimesFunctional() const;
    std::vector<int> findPrimesDomain() const;
    std::vector<int> findPrimesDiv() const;

private:
    static bool isPrime_div(const int number);
    std::vector<std::vector<int>> getRanges() const;
    void combinePrimesLists(const std::vector<std::vector<bool>> &threadsPrimesLists, std::vector<bool> &isPrime) const;
    void removeMultiples(const int prime, std::vector<bool> &isPrime, const int sub) const;
    void fillPrimesList(const std::vector<bool> &isPrime, std::vector<int> &primes, const int add) const;

    int mThreadsNum;
    int mMinNum;
    int mMaxNum;
};
