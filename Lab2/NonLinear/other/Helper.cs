namespace NonLinear.helper;

public struct Elem     /// Структура конечного элемента
{
    public double x1;  /// Координата начала КЭ
    public double x2;  /// Координата конца КЭ
    public double h;   /// Длина КЭ

    public Elem(double _x1, double _x2, double _h) {
        x1 = _x1; x2 = _x2; h = _h;
    }

    public override string ToString() => $"{x1,20} {x2,24} {h,24}";
}

public struct SLAU                 /// Структура СЛАУ
{
    public double[] di, dl, du;    /// Матрица
    public double[] q, q1;         /// Вектора решений
    public double[] f;             /// Правая часть
    public int MaxIter;            /// Максимальное количество итераций
    public double EPS;             /// Точность

    //* Обнуление СЛАУ (левой и правой частей)
    public void Clear() {
        Array.Clear(di, 0, di.Length);
        Array.Clear(dl, 0, dl.Length);
        Array.Clear(du, 0, du.Length);
        Array.Clear(f , 0, f.Length);
    }

    //* Перемножение матрицы на вектор
    public void Mult(double[] vec, out double[] res_vec) {
        res_vec = new double[vec.Length];
        res_vec[0] = di[0] * q[0] + du[0] * q[1];
        for (int i = 1; i < vec.Length - 1; i++) {
            res_vec[i]  = di[i] * q[i];
            res_vec[i] += du[i] * q[i + 1]; 
            res_vec[i] += dl[i] * q[i - 1];
        }
        res_vec[vec.Length - 1] = di[vec.Length - 1] * q[vec.Length - 1] + dl[vec.Length - 1] * q[vec.Length - 2];
    }

    //* Метод прогонки
    public double[] Progonka() {
        var alpha = new double[q.Length];
        var betta = new double[q.Length];
        
        // Подсчет альф и бетт
        alpha[0] = -du[0] / di[0]; 
        betta[0] = f[0] / di[0];
        for (int i = 1; i < q.Length; i++) {
            alpha[i] = -du[i] / (di[i] + dl[i] * alpha[i - 1]);
            betta[i] = (f[i] - dl[i] * betta[i - 1]) / (di[i] + dl[i] * alpha[i - 1]);
        }

        q[q.Length - 1] = betta[q.Length - 1];
        for (int i = q.Length - 2; i >= 0; i--)
            q[i] = alpha[i] * q[i + 1] + betta[i];
        return q;
    }
}

public static class Helper
{
    //: ***************** Перечисления ***************** :\\
    public enum Method {
        Iteration, 
        Newton
    }
    //: ***************** Перечисления ***************** :\\

    //* Вычисление нормы вектора
    public static double Norm(double[] array) {
        double norm = 0;
        for (int i = 0; i < array.Count(); i++)
            norm += array[i] * array[i];
        return Sqrt(norm);
    }

    //* Окно помощи при запуске (если нет аргументов или по команде)
    public static void ShowHelp() {
        WriteLine("----Команды----                        \n" + 
        "-help             - показать справку             \n" + 
        "-i                - входной файл                 \n" + 
        "-m iteration      - метод простой итерации       \n" + 
        "-m newton         - метод Ньютона                \n" + 
        "-func 'number'    - номер функции                \n");
    }
}
