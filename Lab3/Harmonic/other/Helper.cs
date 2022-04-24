namespace Harmonic.other;

public struct Grid    /// Структура сетки
{
    public int Count_Node  { get; set; }    /// Общее количество узлов
    public int Count_Elem  { get; set; }    /// Общее количество КЭ
    public int Count_Kraev { get; set; }    /// Количество краевых

    public Node[]  Nodes;      /// Узлы
    public Elem[]  Elems;      /// КЭ
    public Kraev[] Kraevs;     /// Краевые           

    public Grid(int count_node, int count_elem, int count_kraev, Node[] nodes, Elem[] elem, Kraev[] kraevs) {
        this.Count_Node  = count_node;
        this.Count_Elem  = count_elem;
        this.Count_Kraev = count_kraev;
        this.Nodes       = nodes;
        this.Elems       = elem;
        this.Kraevs      = kraevs;
    }
}

public struct Node     /// Структура Узла
{
    public double x { get; set; }  /// Координата X 
    public double y { get; set; }  /// Координата Y 
    public double z { get; set; }  /// Координата Z 

    public Node(double _x, double _y, double _z) {
        x = _x; y = _y; z = _z;
    }

    public void Deconstructor(out double x, 
                              out double y,
                              out double z) 
    {
        x = this.x;
        y = this.y;
        z = this.z;
    }

    public override string ToString() => $"{x,20} {y,24} {z,24}";
}

public struct Elem     /// Структура КЭ
{
    public int[] Node;  /// Узлы КЭ

    public Elem(params int[] node) { Node = node; }

    public void Deconstructor(out int[] nodes) { nodes = this.Node; }

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
    public int[] Node;   /// Узлы КЭ
    public double   Value;  /// Значение краевого (на )

    public Kraev(params int[] node) { Node = node; Value = 0; }

    public Kraev(int node, double value) {
        this.Node  = new int[1] { node };
        this.Value = value;
    }

    public void Deconstructor(out int[] nodes) { nodes = this.Node; }

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
    public Vector di, gl, gu;            /// Матрица
    public int[] ig, jg; 
    public Vector f;                     /// Правая часть
    public int N;                        /// Размерность
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
