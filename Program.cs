using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace numAnalysis2
{
    //Решить указанную систему линейных алгебраических уравнений Ax = b, используя методы:
    //-метод Гаусса
    //-модификация метода Гаусса
    //-LU-алгоритм
    //Оценить погрешность полученного решения по правой части.
    //Вычислить определитель матрицы A, используя указанные выше методы.
    //Построить обратную матрицу A^-1, используя метод Гаусса.
    //
    //Вариант 5
    // 76x1 + 21x2 + 6x3 - 34x4 = 142
    // 12x1 - 114x2 + 8x3 + 9x4 = 83
    // 16x1 + 24x2 - 100x3 + 35x4 = 121
    // 23x1 - 8x2 + 5x3 - 75x4 = 85

    //Класс для описания методов вызываемых от экземпляров матрицы double[,]
    public static class DoubleMatrixExtension
    {
        public static double[,] ColumnSwap(this double[,] matrix, int num1, int num2)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                matrix[i, num1] += matrix[i, num2];
                matrix[i, num2] = matrix[i, num1] - matrix[i, num2];
                matrix[i, num1] -= matrix[i, num2];
            }
            return matrix;
        }

        public static double[,] RowSwap(this double[,] matrix, int num1, int num2)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                matrix[num1, i] += matrix[num2, i];
                matrix[num2, i] = matrix[num1, i] - matrix[num2, i];
                matrix[num1, i] -= matrix[num2, i];
            }
            return matrix;
        }

        public static double DecisionError(this double[,] matrix, double[] freeMatrixVector, double[] resultVector)
        {
            double[] deltaB = new double[matrix.GetLength(0)];
            for(int i = 0; i < matrix.GetLength(0); i++)
            {
                for(int j = 0; j < matrix.GetLength(1); j++)
                {
                    deltaB[i] += matrix[i, j] * resultVector[j];
                }
                deltaB[i] = Math.Abs(deltaB[i] - freeMatrixVector[i]);
            }
            double maxFreeVector = 0;
            foreach(double d in freeMatrixVector)
            {
                if(Math.Abs(d) > maxFreeVector)
                {
                    maxFreeVector = Math.Abs(d);
                }
            }
            return deltaB.Max() / maxFreeVector;
        }
    }

    class Program
    { 
        public static void MatrixOutput(double[,] matrix, double[] resultVector)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    Console.Write($"{Math.Round(matrix[i, j], 4),9}");
                }
                Console.WriteLine($" \td = {resultVector[i]}\n");
            }

        }
        public static double[] Gauss(double[,] matrixS, double[] freeMatrixVectorS)
        {
            double[,] matrix = new double[matrixS.GetLength(0), matrixS.GetLength(0)];
            for (int i = 0; i < matrixS.GetLength(0); i++)
            {
                for(int j = 0; j < matrixS.GetLength(1); j++)
                {
                    matrix[i, j] = matrixS[i, j];
                }
            }
            double[] freeMatrixVector = new double[freeMatrixVectorS.Length];
            for(int i = 0; i < freeMatrixVectorS.Length; i++)
            {
                freeMatrixVector[i] = freeMatrixVectorS[i];
            }

            //Прямой ход
            for(int i = 0; i < matrix.GetLength(0) - 1; i++)
            {
                double divider = matrix[i, i];
                for(int j = i; j < matrix.GetLength(1); j++)
                {
                    matrix[i, j] /= divider;
                }
                freeMatrixVector[i] /= divider;
                for(int j = i + 1; j < matrix.GetLength(1); j++)
                {
                    double factor = matrix[j, i];
                    for (int k = 0; k < matrix.GetLength(0); k++)
                    {
                        matrix[j, k] -= matrix[i, k] * factor;
                    }
                    freeMatrixVector[j] -= freeMatrixVector[i] * factor;
                }
            }

            freeMatrixVector[freeMatrixVector.Length - 1] /= matrix[matrix.GetLength(0) - 1, matrix.GetLength(1) - 1];
            matrix[matrix.GetLength(0) - 1, matrix.GetLength(1) - 1] = 1;

            Console.WriteLine("Gauss method:\n\tupper triangular matrix");
            MatrixOutput(matrix, freeMatrixVector);

            //Обратный ход
            for (int i = matrix.GetLength(0) -2; i >= 0; i--)
            {
                double num = 0;
                for (int j = matrix.GetLength(1) - 1; j > i; j--)
                {
                    num += matrix[i, j] * freeMatrixVector[j];
                }
                freeMatrixVector[i] -= num;
            }

            Console.WriteLine("Result :");
            for (int i = 0; i < freeMatrixVector.Length; i++)
            {
                Console.WriteLine($"x[{i + 1}] = {freeMatrixVector[i]}");
            }
            return freeMatrixVector;
        }

        public static double[] GaussModifyColumnMax(double[,] matrixS, double[] freeMatrixVectorS)
        {
            double[,] matrix = new double[matrixS.GetLength(0), matrixS.GetLength(0)];
            for (int i = 0; i < matrixS.GetLength(0); i++)
            {
                for (int j = 0; j < matrixS.GetLength(1); j++)
                {
                    matrix[i, j] = matrixS[i, j];
                }
            }
            double[] freeMatrixVector = new double[freeMatrixVectorS.Length];
            for (int i = 0; i < freeMatrixVectorS.Length; i++)
            {
                freeMatrixVector[i] = freeMatrixVectorS[i];
            }

            //Прямой ход
            for (int i = 0; i < matrix.GetLength(0) - 1; i++)
            {
                //Нахождение ведущего элемента по столбцу
                #region maxInColumn
                int swapNum = 0;
                double max = 0;
                for (int j = i; j < matrix.GetLength(0); j++)
                {
                    if (Math.Abs(matrix[j, i]) > Math.Abs(max))
                    {
                        max = matrix[j, i];
                        swapNum = j;
                    }
                }
                if (swapNum != i)
                {
                    matrix = matrix.RowSwap(i, swapNum);
                    freeMatrixVector[i] += freeMatrixVector[swapNum];
                    freeMatrixVector[swapNum] = freeMatrixVector[i] - freeMatrixVector[swapNum];
                    freeMatrixVector[i] -= freeMatrixVector[swapNum];
                }
                #endregion

                double divider = matrix[i, i];
                for (int j = i; j < matrix.GetLength(1); j++)
                {
                    matrix[i, j] /= divider;
                }
                freeMatrixVector[i] /= divider;
                for (int j = i + 1; j < matrix.GetLength(1); j++)
                {
                    double factor = matrix[j, i];
                    for (int k = 0; k < matrix.GetLength(0); k++)
                    {
                        matrix[j, k] -= matrix[i, k] * factor;
                    }
                    freeMatrixVector[j] -= freeMatrixVector[i] * factor;
                }
            }

            freeMatrixVector[freeMatrixVector.Length - 1] /= matrix[matrix.GetLength(0) - 1, matrix.GetLength(1) - 1];
            matrix[matrix.GetLength(0) - 1, matrix.GetLength(1) - 1] = 1;

            Console.WriteLine("\nGauss method (max element in column):\n\tupper triangular matrix");
            MatrixOutput(matrix, freeMatrixVector);

            //Обратный ход
            for (int i = matrix.GetLength(0) - 2; i >= 0; i--)
            {
                double num = 0;
                for (int j = matrix.GetLength(1) - 1; j > i; j--)
                {
                    num += matrix[i, j] * freeMatrixVector[j];
                }
                freeMatrixVector[i] -= num;
            }

            Console.WriteLine("Result :");
            for (int i = 0; i < freeMatrixVector.Length; i++)
            {
                Console.WriteLine($"x[{i + 1}] = {freeMatrixVector[i]}");
            }
            return freeMatrixVector;
        }

        public static double[] GaussModifyRowMax(double[,] matrixS, double[] freeMatrixVectorS)
        {
            double[,] matrix = new double[matrixS.GetLength(0), matrixS.GetLength(0)];
            for (int i = 0; i < matrixS.GetLength(0); i++)
            {
                for (int j = 0; j < matrixS.GetLength(1); j++)
                {
                    matrix[i, j] = matrixS[i, j];
                }
            }
            double[] freeMatrixVector = new double[freeMatrixVectorS.Length];
            for (int i = 0; i < freeMatrixVectorS.Length; i++)
            {
                freeMatrixVector[i] = freeMatrixVectorS[i];
            }

            //Прямой ход
            int[] vectorX = { 1, 2, 3, 4 };
            for (int i = 0; i < matrix.GetLength(0) - 1; i++)
            {
                //Нахождение ведущего элемента по строке
                #region maxInRow
                int swapNum = 0;
                double max = 0;
                for (int j = i; j < matrix.GetLength(0); j++)
                {
                    if (Math.Abs(matrix[i, j]) > Math.Abs(max))
                    {
                        max = matrix[i, j];
                        swapNum = j;
                    }
                }
                if (swapNum != i)
                {
                    matrix = matrix.ColumnSwap(i, swapNum);
                    vectorX[i] += vectorX[swapNum];
                    vectorX[swapNum] = vectorX[i] - vectorX[swapNum];
                    vectorX[i] -= vectorX[swapNum];
                }
                #endregion

                double divider = matrix[i, i];
                for (int j = i; j < matrix.GetLength(1); j++)
                {
                    matrix[i, j] /= divider;
                }
                freeMatrixVector[i] /= divider;
                for (int j = i + 1; j < matrix.GetLength(1); j++)
                {
                    double factor = matrix[j, i];
                    for (int k = 0; k < matrix.GetLength(0); k++)
                    {
                        matrix[j, k] -= matrix[i, k] * factor;
                    }
                    freeMatrixVector[j] -= freeMatrixVector[i] * factor;
                }
            }

            freeMatrixVector[freeMatrixVector.Length - 1] /= matrix[matrix.GetLength(0) - 1, matrix.GetLength(1) - 1];
            matrix[matrix.GetLength(0) - 1, matrix.GetLength(1) - 1] = 1;

            Console.WriteLine("\nModify Gauss method (max element in row):\n\tupper triangular matrix");
            MatrixOutput(matrix, freeMatrixVector);

            //Обратный ход
            for (int i = matrix.GetLength(0) - 2; i >= 0; i--)
            {
                double num = 0;
                for (int j = matrix.GetLength(1) - 1; j > i; j--)
                {
                    num += matrix[i, j] * freeMatrixVector[j];
                }
                freeMatrixVector[i] -= num;
            }

            Console.WriteLine("Result :");
            for (int i = 0; i < freeMatrixVector.Length; i++)
            {
                Console.WriteLine($"x[{vectorX[i]}] = {freeMatrixVector[i]}");
            }

            //Сортировка результирующего вектора
            bool IsNotSortArray = true;
            while (IsNotSortArray)
            {
                IsNotSortArray = false;
                for (int i = 0; i < vectorX.Length; i++)
                {
                    if (vectorX[i] != i + 1)
                    {
                        IsNotSortArray = true;
                        double elem = freeMatrixVector[i];
                        freeMatrixVector[i] = freeMatrixVector[vectorX[i] - 1];
                        freeMatrixVector[vectorX[i] - 1] = elem;

                        elem = vectorX[i];
                        vectorX[i] = vectorX[vectorX[i] - 1];
                        vectorX[(int)elem - 1] = (int)elem;
                    }
                }
            }
            return freeMatrixVector;
        }

        public static double[] GaussModify(double[,] matrixS, double[] freeMatrixVectorS)
        {
            double[,] matrix = new double[matrixS.GetLength(0), matrixS.GetLength(0)];
            for (int i = 0; i < matrixS.GetLength(0); i++)
            {
                for (int j = 0; j < matrixS.GetLength(1); j++)
                {
                    matrix[i, j] = matrixS[i, j];
                }
            }
            double[] freeMatrixVector = new double[freeMatrixVectorS.Length];
            for (int i = 0; i < freeMatrixVectorS.Length; i++)
            {
                freeMatrixVector[i] = freeMatrixVectorS[i];
            }

            //Прямой ход
            double[] vectorX = { 1, 2, 3, 4 };
            for (int i = 0; i < matrix.GetLength(0) - 1; i++)
            {
                //Нахождение ведущего элемента по матрице
                #region maxInMatrix
                int swapCol = 0, swapRow = 0;
                double max = 0;
                for (int j = i; j < matrix.GetLength(0); j++)
                {
                    for(int k = i; k < matrix.GetLength(1); k++)
                    {
                        if (Math.Abs(matrix[j, k]) > Math.Abs(max))
                        {
                            max = matrix[j, k];
                            swapCol = k;
                            swapRow = j;
                        }
                    }
                }
                if (swapRow != i)
                {
                    matrix = matrix.RowSwap(i, swapRow);
                    freeMatrixVector[i] += freeMatrixVector[swapRow];
                    freeMatrixVector[swapRow] = freeMatrixVector[i] - freeMatrixVector[swapRow];
                    freeMatrixVector[i] -= freeMatrixVector[swapRow];
                }
                if (swapCol != i)
                {
                    matrix = matrix.ColumnSwap(i, swapCol);
                    vectorX[i] += vectorX[swapCol];
                    vectorX[swapCol] = vectorX[i] - vectorX[swapCol];
                    vectorX[i] -= vectorX[swapCol];
                }
                #endregion

                double divider = matrix[i, i];
                for (int j = i; j < matrix.GetLength(1); j++)
                {
                    matrix[i, j] /= divider;
                }
                freeMatrixVector[i] /= divider;
                for (int j = i + 1; j < matrix.GetLength(1); j++)
                {
                    double factor = matrix[j, i];
                    for (int k = 0; k < matrix.GetLength(0); k++)
                    {
                        matrix[j, k] -= matrix[i, k] * factor;
                    }
                    freeMatrixVector[j] -= freeMatrixVector[i] * factor;
                }
            }

            freeMatrixVector[freeMatrixVector.Length - 1] /= matrix[matrix.GetLength(0) - 1, matrix.GetLength(1) - 1];
            matrix[matrix.GetLength(0) - 1, matrix.GetLength(1) - 1] = 1;

            Console.WriteLine("\nGauss method (modify):\n\tupper triangular matrix");
            MatrixOutput(matrix, freeMatrixVector);

            //Обратный ход
            for (int i = matrix.GetLength(0) - 2; i >= 0; i--)
            {
                double num = 0;
                for (int j = matrix.GetLength(1) - 1; j > i; j--)
                {
                    num += matrix[i, j] * freeMatrixVector[j];
                }
                freeMatrixVector[i] -= num;
            }

            Console.WriteLine("Result :");
            for (int i = 0; i < freeMatrixVector.Length; i++)
            {
                Console.WriteLine($"x[{vectorX[i]}] = {freeMatrixVector[i]}");
            }

            //Сортировка результирующего вектора
            bool IsNotSortArray = true;
            while (IsNotSortArray)
            {
                IsNotSortArray = false;
                for (int i = 0; i < vectorX.Length; i++)
                {
                    if (vectorX[i] != i + 1)
                    {
                        IsNotSortArray = true;
                        double elem = freeMatrixVector[i];
                        freeMatrixVector[i] = freeMatrixVector[(int)vectorX[i] - 1];
                        freeMatrixVector[(int)vectorX[i] - 1] = elem;

                        elem = vectorX[i];
                        vectorX[i] = vectorX[(int)vectorX[i] - 1];
                        vectorX[(int)elem - 1] = (int)elem;
                    }
                }
            }
            return freeMatrixVector;
        }

        public static double[] LUAlg(double[,] matrixS, double[] freeMatrixVectorS)
        {
            double[,] matrix = new double[matrixS.GetLength(0), matrixS.GetLength(0)];
            for (int i = 0; i < matrixS.GetLength(0); i++)
            {
                for (int j = 0; j < matrixS.GetLength(1); j++)
                {
                    matrix[i, j] = matrixS[i, j];
                }
            }
            double[] freeMatrixVector = new double[freeMatrixVectorS.Length];
            for (int i = 0; i < freeMatrixVectorS.Length; i++)
            {
                freeMatrixVector[i] = freeMatrixVectorS[i];
            }

            //Нахождение матриц L и U
            double[,] U = new double[matrix.GetLength(0), matrix.GetLength(1)];
            double[,] L = new double[matrix.GetLength(0), matrix.GetLength(1)];
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                U[0, i] = matrix[0, i];
                L[i, 0] = matrix[i, 0] / U[0, 0];
            }
            for(int i = 1; i < matrix.GetLength(0); i++)
            {
                for(int j = i; j < matrix.GetLength(1); j++)
                {
                    double sumU = 0, sumL = 0;
                    for(int k = 0; k < i; k++)
                    {
                        sumU += U[k, j] * L[i, k];
                        sumL += U[k, i] * L[j, k];
                    }
                    U[i, j] = matrix[i, j] - sumU;
                    L[j, i] = (matrix[j, i] - sumL) / U[i, i];
                }
            }

            Console.WriteLine("\nLU method:\n\tupper and lower triangular matrices");
            Console.WriteLine("\nmatrix U:");
            MatrixOutput(U, freeMatrixVector);

            Console.WriteLine("\nmatrix L:");
            MatrixOutput(L, freeMatrixVector);

            //Нахождение вектора y
            double[] yVector = new double[matrix.GetLength(0)];
            yVector[0] = freeMatrixVector[0];
            for(int i = 1; i < matrix.GetLength(0); i++)
            {
                double num = 0;
                for(int j = 0; j < i; j++)
                {
                    num += yVector[j] * L[i, j];
                }
                yVector[i] = freeMatrixVector[i] - num;
            }
            Console.WriteLine("Vector y:");
            for(int i = 0; i < yVector.Length; i++)
            {
                Console.WriteLine($"y[{i + 1}] = {yVector[i]}");
            }

            //Нахожденеие результирующего вектора
            freeMatrixVector[freeMatrixVector.Length - 1] = yVector[yVector.Length - 1] / U[U.GetLength(0) - 1, U.GetLength(1) - 1];
            for(int i = matrix.GetLength(1) - 2; i >= 0; i--)
            {
                double num = 0;
                for(int j = matrix.GetLength(1) - 1; j > i; j--)
                {
                    num += freeMatrixVector[j] * U[i, j];
                }
                freeMatrixVector[i] = (yVector[i] - num) / U[i, i];
            }
            Console.WriteLine("\nResult :");
            for (int i = 0; i < freeMatrixVector.Length; i++)
            {
                Console.WriteLine($"x[{i + 1}] = {freeMatrixVector[i]}");
            }

            return freeMatrixVector;
        }

        static void Main(string[] args)
        {
            //double[,] matrix = { {76, 21, 6, -34 }, {12, -114, 8, 9 }, {16, 24, -100, 35 }, {23, -8, 5, -75 } };
            //double[] freeMatrixVector = { 142, 83, 121, 85 };

            double[,] matrix = { { 2, 1, -1, 1 }, { 0.4, 0.5, 4, -8.5 }, { 0.3, -1, 1, 5.2 }, { 1, 0.2, 2.5, -1 } };
            double[] freeMatrixVector = { 2.7, 21.9, -3.9, 9.9 };
            double[] resultVector = Gauss(matrix, freeMatrixVector);
            Console.WriteLine($"Gauss error :\t{matrix.DecisionError(freeMatrixVector, resultVector)}");
            resultVector = GaussModifyRowMax(matrix, freeMatrixVector);
            Console.WriteLine($"Gauss row modify error :\t{matrix.DecisionError(freeMatrixVector, resultVector)}");
            resultVector = GaussModifyColumnMax(matrix, freeMatrixVector);
            Console.WriteLine($"Gauss column modify error :\t{matrix.DecisionError(freeMatrixVector, resultVector)}");
            resultVector = GaussModify(matrix, freeMatrixVector);
            Console.WriteLine($"Gauss modify error :\t{matrix.DecisionError(freeMatrixVector, resultVector)}");
            resultVector = LUAlg(matrix, freeMatrixVector);
            Console.WriteLine($"LU algoritm error :\t{matrix.DecisionError(freeMatrixVector, resultVector)}");

            ////Console.WriteLine("Start matrix");
            ////for (int i = 0; i < matrix.GetLength(0); i++)
            ////{
            ////    for (int j = 0; j < matrix.GetLength(1); j++)
            ////    {
            ////        Console.Write($"{matrix[i, j],6}");
            ////    }
            ////    Console.WriteLine($"\t|{freeMatrixVector[i]}");
            ////}
        }
    }
}
