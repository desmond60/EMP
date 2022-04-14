namespace Harmonic.other;

public struct Issue    /// Структура сетки
{
    public int N_X { get; }                               /// Количество узлов по Оси X
    public int N_Y { get; }                               /// Количество узлов по Оси Y
    public int N_Z { get; }                               /// Количество узлов по Оси Z
    public int Count_Node => N_X * N_Y * N_Z;             /// Общее количество узлов
    public int Count_Elem => (N_X - 1)*(N_Y - 1)*(N_Z-1); /// Общее количество КЭ

    public Node[]  Nodes;      /// Узлы
    public Elem[]  Elems;      /// КЭ
    public Kraev[] Kraevs;     /// Краевые           

    public Issue(int n_x, int n_y, int n_z, Node[] nodes, Elem[] elem, Kraev[] kraevs) {
        this.N_X    = n_x;
        this.N_Y    = n_y;
        this.N_Z    = n_z;
        this.Nodes  = nodes;
        this.Elems  = elem;
        this.Kraevs = kraevs;
    }
}

public struct Node     /// Структура Узла
{
    public double x;  /// Координата X 
    public double y;  /// Координата Y 
    public double z;  /// Координата Z 

    public Node(double _x, double _y, double _z) {
        x = _x; y = _y; z = _z;
    }

    public override string ToString() => $"{x,20} {y,24} {z,24}";
}

public struct Elem     /// Структура КЭ
{
    public double[] Node;  /// Узлы КЭ

    public Elem(params double[] node) { Node = node; }

    public override string ToString() {
        StringBuilder str_elem = new StringBuilder();
        str_elem.Append($"{Node[0],5}");
        for (int i = 1; i < Node.Count(); i++)
            str_elem.Append($"{Node[i],8}");
        return str_elem.ToString();
    }
}

public struct Kraev    /// Структура краевого
{
    public double[] Node;  /// Узлы КЭ

    public Kraev(params double[] node) { Node = node; }

    public override string ToString() {
        StringBuilder str_elem = new StringBuilder();
        str_elem.Append($"{Node[0],5}");
        for (int i = 1; i < Node.Count(); i++)
            str_elem.Append($"{Node[i],8}");
        return str_elem.ToString();
    }
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
