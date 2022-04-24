using System.Threading.Tasks.Dataflow;
namespace Harmonic;
public class Generate
{
    //* Данные 
    protected Vector start  { get; set; }      /// Начальная точка
    protected Vector end    { get; set; }      /// Конечная точка
    protected double hx     { get; set; }      /// Шаг по Оси X
    protected double hy     { get; set; }      /// Шаг по Оси Y
    protected double hz     { get; set; }      /// Шаг по Оси Z 
    protected double kx     { get; set; }      /// Коэффициент деления по Оси X  
    protected double ky     { get; set; }      /// Коэффициент деления по Оси Y  
    protected double kz     { get; set; }      /// Коэффициент деления по Оси Z  

    protected string Path   { get; set; }      /// Путь к папке с задачей 
    
    private int N_X;                                                                      /// Количество узлов по Оси X
    private int N_Y;                                                                      /// Количество узлов по Оси Y
    private int N_Z;                                                                      /// Количество узлов по Оси Z
    private int Count_Node  => N_X * N_Y * N_Z;                                           /// Общее количество узлов
    private int Count_Elem  => (N_X-1)*(N_Y-1)*(N_Z-1);                                   /// Общее количество КЭ
    private int Count_Kraev => 2*(N_X-1)*(N_Y-1) + 2*(N_X-1)*(N_Z-1) + 2*(N_Y-1)*(N_Z-1); /// Количество краевых

    //* Конструктор
    public Generate(Data data, string Path) {
        (start, end, hx, hy, hz, kx, ky, kz) = data;
        this.Path = Path;

        // Подсчет количества узлов на Осях 
        N_X = kx != 1
            ? (int)(Log(1 - (end[0] - start[0])*(kx - 1) / hx) / Log(kx) + 2)
            : (int)((end[0] - start[0]) / hx + 1);
        N_Y = ky != 1
            ? (int)(Log(1 - (end[1] - start[1])*(ky - 1) / hy) / Log(ky) + 2)
            : (int)((end[1] - start[1]) / hy + 1);
        N_Z = kz != 1
            ? (int)(Log(1 - (end[2] - start[2])*(kz - 1) / hz) / Log(kz) + 2)
            : (int)((end[2] - start[2]) / hz + 1);
    }

    //* Функция генерации
    public Grid generate() {
        Path += "/grid";
        Directory.CreateDirectory(Path);
        Node[]  nodes  = generate_coords(); //? Генерация координат
        Elem[]  elems  = generate_elems();  //? Генерация КЭ
        Kraev[] kraevs = generate_kraevs(); //? Генерация краевых

        return new Grid(Count_Node, Count_Elem, Count_Kraev, nodes, elems, kraevs);
    }

    //* Функция трансформация сетки под синус и косинус
    public Grid transformation(Grid grid) {

        var newNodes = new Node[2*grid.Count_Node];
        var newElems = new Elem[grid.Count_Elem];
        var newKraev = new Kraev[grid.Count_Kraev];

        // Дублирую координаты для расчетов функции синуса и косинуса (так удобнее)
        for (int i = 0, id = 0; i < grid.Count_Node; i++) {
            double x, y, z;
            grid.Nodes[i].Deconstructor(out x, out y, out z);
            newNodes[id++] = new Node(x, y, z);
            newNodes[id++] = new Node(x, y, z);
        }

        // Проходим по элементам и делаем узлы (2i и 2i+1)
        for (int i = 0; i < grid.Count_Elem; i++) {
            int[] node = new int[8];
            grid.Elems[i].Deconstructor(out node);
            
            var newNode = new int[16];
            for (int j = 0, id = 0; j < node.Length; j++) {
                newNode[id++] = 2*node[j];
                newNode[id++] = 2*node[j] + 1;
            }
            newElems[i] = new Elem(newNode);
        }

        // Проходим по краевым и делаем узлы (2i и 2i+1)
        for (int i = 0; i < grid.Count_Kraev; i++) {
            int[] node = new int[4];
            grid.Kraevs[i].Deconstructor(out node);
            
            var newNode = new int[8];
            for (int j = 0, id = 0; j < node.Length; j++) {
                newNode[id++] = 2*node[j];
                newNode[id++] = 2*node[j] + 1;
            }
            newKraev[i] = new Kraev(newNode);
        }

        // Вычислим значения для каждого узла краевого 
        List<int> lNodes = new List<int>();
        for (int i = 0; i < newKraev.Length; i++)
            for (int j = 0; j < 8; j++)
                lNodes.Add(newKraev[i].Node[j]);
        lNodes = lNodes.Distinct().OrderBy(x => x).ToList();
                
        var newnewKraev = new Kraev[lNodes.Count()];
        for (int i = 0; i < lNodes.Count(); i++) {
            int id = lNodes[i];
            Vector node = new Vector(new double[] { newNodes[id].x, newNodes[id].y, newNodes[id].z });
            if (id % 2 == 0)
                newnewKraev[i] = new Kraev(id, Us(node));
            else 
                newnewKraev[i] = new Kraev(id, Uc(node)); 
        }

        return new Grid(2*grid.Count_Node, grid.Count_Elem, grid.Count_Kraev, newNodes, newElems, newnewKraev);
    }

    //* Генерация координат
    private Node[] generate_coords() {
        Vector X_vec = generate_array(start[0], end[0], hx, kx, N_X);
        Vector Y_vec = generate_array(start[1], end[1], hy, ky, N_Y);
        Vector Z_vec = generate_array(start[2], end[2], hz, kz, N_Z);
    
        Node[] nodes = new Node[Count_Node];

        for (int i = 0; i < N_X; i++)
            for (int j = 0; j < N_Y; j++)
                for (int k = 0; k < N_Z; k++)
                    nodes[k*N_X*N_Y + j*N_X + i] = new Node(X_vec[i], Y_vec[j], Z_vec[k]);

        File.WriteAllText(Path + "/coords.txt", String.Join("\n", nodes));
        return nodes;
    }

    //* Генерация массива по Оси (с шагом и коэффицентом)
    private Vector generate_array(double start, double end, double h, double k, int n) {
        var coords = new Vector(n);
        coords[0]     = start;
        coords[n - 1] = end;
        for (int i = 1; i < n; i++, h *= k) 
            coords[i] = coords[i - 1] + h;
        
        return coords;
    }

    //* Генерация КЭ
    private Elem[] generate_elems() {
        Elem[] elems = new Elem[Count_Elem];

        for (int i = 0, id = 0; i < N_X - 1; i++)
            for (int j = 0; j < N_Y - 1; j++)
                for (int k = 0; k < N_Z - 1; k++, id++)
                    elems[id] = new Elem(
                         k   *N_X*N_Y +  j   *N_Y + i,
                         k   *N_X*N_Y +  j   *N_Y + i + 1,
                         k   *N_X*N_Y + (j+1)*N_X + i,
                         k   *N_X*N_Y + (j+1)*N_X + i + 1,
                        (k+1)*N_X*N_Y +  j   *N_Y + i,
                        (k+1)*N_X*N_Y +  j   *N_Y + i + 1,
                        (k+1)*N_X*N_Y + (j+1)*N_X + i,
                        (k+1)*N_X*N_Y + (j+1)*N_X + i + 1
                    );

        File.WriteAllText(Path + "/elems.txt", String.Join("\n", elems));
        return elems;
    }

    //* Генерация краевых
    private Kraev[] generate_kraevs() {
        Kraev[] kraevs = new Kraev[Count_Kraev];
        int id = 0;

        // Нижний слой
        for (int i = 0; i < N_X - 1; i++)
            for (int j = 0; j < N_Y - 1; j++, id++)
                kraevs[id] = new Kraev(
                    j    *N_X + i,
                    j    *N_X + i + 1,
                    (j+1)*N_X + i,
                    (j+1)*N_X + i + 1
                );

        // Верхний слой
        for (int i = 0; i < N_X - 1; i++)
            for (int j = 0; j < N_Y - 1; j++, id++)
                kraevs[id] = new Kraev(
                    (N_Z-1)*N_X*N_Y + j    *N_X + i,
                    (N_Z-1)*N_X*N_Y + j    *N_X + i + 1,
                    (N_Z-1)*N_X*N_Y + (j+1)*N_X + i,
                    (N_Z-1)*N_X*N_Y + (j+1)*N_X + i + 1
                );

        // Правый слой 
        for (int i = 0; i < N_X - 1; i++)
            for (int j = 0; j < N_Z - 1; j++, id++)
                kraevs[id] = new Kraev(
                    j    *N_X*N_Y + i,
                    j    *N_X*N_Y + i + 1,
                    (j+1)*N_X*N_Y + i,
                    (j+1)*N_X*N_Y + i + 1
                );

        // Левый слой 
        for (int i = 0; i < N_X - 1; i++)
            for (int j = 0; j < N_Z - 1; j++, id++)
                kraevs[id] = new Kraev(
                    j    *N_X*N_Y + (N_Y-1)*N_X + i,
                    j    *N_X*N_Y + (N_Y-1)*N_X + i + 1,
                    (j+1)*N_X*N_Y + (N_Y-1)*N_X + i,
                    (j+1)*N_X*N_Y + (N_Y-1)*N_X + i + 1
                );

        // Лицевой слой
        for (int i = 0; i < N_Y - 1; i++)
            for (int j = 0; j < N_Z - 1; j++, id++)
                kraevs[id] = new Kraev(
                    j    *N_X*N_Y + i    *N_X,
                    j    *N_X*N_Y + (i+1)*N_X,
                    (j+1)*N_X*N_Y + i    *N_X,
                    (j+1)*N_X*N_Y + (i+1)*N_X
                );

        // Задний слой
        for (int i = 0; i < N_Y - 1; i++)
            for (int j = 0; j < N_Z - 1; j++, id++)
                kraevs[id] = new Kraev(
                    j    *N_X*N_Y + i    *N_X + N_X-1,
                    j    *N_X*N_Y + (i+1)*N_X + N_X-1,
                    (j+1)*N_X*N_Y + i    *N_X + N_X-1,
                    (j+1)*N_X*N_Y + (i+1)*N_X + N_X-1
                );

        File.WriteAllText(Path + "/kraevs.txt", String.Join("\n", kraevs));
        return kraevs;
    }

}
