#include <iostream>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <cmath>

using std::vector;
using std::pair;

using dmatrix = vector<vector<double>>;

void printMatrix(const dmatrix &matrix) {
    // фиксированное количество знаков после запятой
    std::cout << std::fixed;
    std::cout << std::setprecision(4);

    // вывод коэффициентов
    std::cout << "id  ";
    for (size_t i = 0; i < matrix.size(); ++i) {
        std::cout << i << "        ";
    }
    std::cout << std::endl;

    // сам вывод матрицы
    for (size_t i = 0; i < matrix.size(); ++i) {
        // вывод индекса строки
        std::cout << i << "  ";
        for (size_t j = 0; j < matrix.front().size(); ++j) {
            if (matrix[i][j] >= 0) {
                std::cout << " ";
            }
            std::cout << matrix[i][j] << "  ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void printResult(const vector<double> &vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        std::cout << "x" << i + 1 << " = " << vec[i] << std::endl;
    }
}

// принимает расширенную матрицу, возвращает вектор коэффициентов
vector<double> solve(dmatrix matrix) {
    // проверка на пустую матрицу
    if (matrix.empty()) {
        throw std::runtime_error("Empty matrix");
    }

    // матрица должна иметь размер n * (n + 1)
    if (matrix.size() != (matrix.front().size() - 1)) {
        throw std::runtime_error("Incorrect matrix dimensions");
    }

    for (size_t i = 0; i < matrix.size(); ++i) {
        // поиск максимального элемента
        size_t jmax = i;

        for (size_t j = i + 1; j < matrix.size(); ++j) {
            if (std::abs(matrix[j][i]) > std::abs(matrix[jmax][i])) {
                jmax = j;
            }
        }

        // максимальный элемент не должен быть нулевым
        if (matrix[jmax][i] == 0) {
            throw std::runtime_error("Unsolvable matrix");
        }

        if (jmax != i) {
            // перемещение строки с максимальным элементов наверх
            std::swap(matrix[i], matrix[jmax]);
        }

        // нормализация строки
        for (size_t j = i + 1; j < matrix.front().size(); ++j) {
            matrix[i][j] /= matrix[i][i];
        }
        matrix[i][i] = 1;
        // вычитание текущей строки из остальных
        for (size_t j = 0; j < matrix.size(); ++j) {
            if (i != j) {
                for (size_t k = i + 1; k < matrix.front().size(); ++k) {
                    matrix[j][k] -= matrix[j][i] * matrix[i][k];
                }
                matrix[j][i] = 0;
            }
        }
    }

    vector<double> x(matrix.size());

    for (size_t i = 0; i < matrix.size(); ++i) {
        x[i] = matrix[i].back();
    }

    return x;
}

// числа с плавающей точкой должны сравниваться с учетом погрешности
bool doubleEqual(double d1, double d2) {
    return std::abs(d1 - d2) < 0.00001;
}

void runTest(const dmatrix &input,
             const vector<double> &expected) {
    // номер теста
    static int testNumber = 0;

    ++testNumber;

    auto result = solve(input);

    // флаг правильности теста
    bool isOk = true;

    if (expected.size() == result.size()) {
        for (size_t i = 0; i < expected.size(); ++i) {
            if (!doubleEqual(expected[i], result[i])) {
                isOk = false;
                break;
            }
        }
    } else {
        isOk = false;
    }

    std::cout << "[TEST " << testNumber << "] ";

    if (isOk) {
        std::cout << "OK";
    } else {
        std::cout << "FAILED";
    }

    std::cout << std::endl;
}

int main() {
    vector<pair<dmatrix, vector<double>>> testCases;

    testCases.emplace_back(dmatrix{{5, -13, 13, -5, -10, -14},
                                   {5, -7, 12, 6, 6, 57},
                                   {-8, 11, 1, -8, -1, 101},
                                   {-1, 5, -9, -7, 13, -43},
                                   {8, -1, -2, -10, -5, -82}},
                           vector<double>{-4, 7, 9, 2, 1});

    testCases.emplace_back(dmatrix{{2, -2, -4, -2, 3, -22},
                                   {8, 8, -6, -7, -4, -72},
                                   {4, -3, -1, -8, 2, -15},
                                   {10, -5, 4, 1, 7, 69},
                                   {-7, -5, 7, 3, -6, 70}},
                           vector<double>{3, -5, 6, 4, -2});

    testCases.emplace_back(dmatrix{{3, 1, 5, -2, 3, 35},
                                   {4, 3, -7, 5, 6, 54},
                                   {-7, 5, 4, 1, -1, -96},
                                   {1, 4, 1, -3, -10, -71},
                                   {6, -9, -8, -8, -2, 59}},
                           vector<double>{10, -7, 2, 5, 4});

    testCases.emplace_back(dmatrix{{3, 1, 5, -2, 3, 35},
                                   {4, 3, -7, 5, 6, 54},
                                   {-7, 5, 4, 1, -1, -96},
                                   {1, 4, 1, -3, -10, -71},
                                   {6, -9, -8, -8, -2, 59}},
                           vector<double>{0, 0, 0, 0, 0});


    for (size_t i = 0; i < testCases.size(); ++i) {
        runTest(testCases[i].first, testCases[i].second);
    }

    return 0;
}
