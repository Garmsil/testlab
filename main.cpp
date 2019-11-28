#include <iostream>
#include <vector>
#include <stdexcept>
#include <iomanip>

using std::vector;

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

    std::cout << "Source matrix:" << std::endl;
    printMatrix(matrix);

    for (size_t i = 0; i < matrix.size(); ++i) {
        std::cout << "[Iteration " << i + 1 << "]" << std::endl;
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

        std::cout << "Maximum module of " << i << " column is " << matrix[jmax][i]
                  << " in row " << jmax << std::endl;

        if (jmax != i) {
            std::cout << "Swapping rows " << i << " and " << jmax << std::endl;

            // перемещение строки с максимальным элементов наверх
            std::swap(matrix[i], matrix[jmax]);

            printMatrix(matrix);
        }


        std::cout << "Normalizing row " << i << " with division by " << matrix[i][i] << std::endl;
        // нормализация строки
        for (size_t j = i + 1; j < matrix.front().size(); ++j) {
            matrix[i][j] /= matrix[i][i];
        }
        matrix[i][i] = 1;

        printMatrix(matrix);

        std::cout << "Subtracting row " << i << " from others" << std::endl;
        // вычитание текущей строки из остальных
        for (size_t j = 0; j < matrix.size(); ++j) {
            if (i != j) {
                for (size_t k = i + 1; k < matrix.front().size(); ++k) {
                    matrix[j][k] -= matrix[j][i] * matrix[i][k];
                }
                matrix[j][i] = 0;
            }
        }

        printMatrix(matrix);
    }

    vector<double> x(matrix.size());

    for (size_t i = 0; i < matrix.size(); ++i) {
        x[i] = matrix[i].back();
    }

    return x;
}

int main() {
    dmatrix matrix{
        {5, -13, 13, -5, -10, -14},
        {5, -7, 12, 6, 6, 57},
        {-8, 11, 1, -8, -1, 101},
        {-1, 5, -9, -7, 13, -43},
        {8, -1, -2, -10, -5, -82},
        };

    auto result = solve(matrix);

    std::cout << "Equation system solved!" << std::endl;
    printResult(result);

    return 0;
}
