namespace FDM;
public struct Matrix /// Структура матрицы
{
    public int N;
    public int maxIter;
    public double EPS;
    public int shift1, dshift, shift2; 
    public double[] di, du1, du2, dl1, dl2, pr, x, absolut_x;
}

public static class Helper
{
    //* Вычисление нормы вектора
    public static double Norm(double[] array) {
        double norm = 0;
        for (int i = 0; i < array.Count(); i++)
            norm += array[i] * array[i];
        return Sqrt(norm);
    }
}
